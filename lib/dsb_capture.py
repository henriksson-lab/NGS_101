"""Bounded constructions for paired-end break capture libraries."""
from dataclasses import dataclass

import seqprimers as sp
from chemdraw import Construct, Segment, feature
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, SequencingPlan, StepKind
from protocol_model import role


@dataclass(frozen=True)
class BreakCaptureDesign:
    name: str
    proximal_adapter: str
    distal_adapter: str
    immobilization: str
    direct_break_read: bool

    def __post_init__(self):
        if not all((self.name, self.proximal_adapter, self.distal_adapter,
                    self.immobilization)):
            raise ValueError("break capture requires named proximal/distal ends and immobilization")
        if self.proximal_adapter == self.distal_adapter:
            raise ValueError("the native-break and fragmentation ends must remain distinguishable")


def _illumina_break_library(name: str, *, inferred: bool, tagged: bool,
                            p5_arm: str, junction: str, index_arm: str, p7_end: str,
                            index_length: int = 6,
                            index_inferred: bool = False,
                            proximal_tag: str | None = None,
                            distal_tag: str | None = None) -> tuple[Construct, tuple]:
    """Build a two-end single-index library from caller-supplied end sequences.

    The P5/Read 1 arm, the ligation junction base(s), the Index 1/Read 2 arm and
    the P7 end are protocol facts: pass the published adapter or primer bases (or
    canonical bases with ``inferred=True``).  The constructor only checks that the
    TruSeq Read 1, Index 1 and Read 2 primers land on the assembled molecule.
    """
    if tagged and (not proximal_tag or not distal_tag):
        raise ValueError("tagged break library needs both end tags")
    if not all((p5_arm, index_arm, p7_end)):
        raise ValueError("break library needs explicit P5, index-read and P7 end sequences")
    index_inferred = inferred or index_inferred
    left = Segment("P5 + Read 1 arm", p5_arm, "p5", inferred=inferred)
    payload = ([Segment("proximal read tag", proximal_tag, "cbc",
                        feature=feature("proximal_break_end", "feature_barcode", "fixed")),
                role("break-to-random genomic DNA"),
                Segment("distal read tag", distal_tag, "cbc",
                        feature=feature("distal_fragment_end", "feature_barcode", "fixed"))]
               if tagged else [role("DNA beginning at repaired DSB")])
    junction_parts = [Segment("ligation junction", junction, inferred=inferred)] if junction else []
    lib = Construct([left, *payload, *junction_parts,
                     Segment("Index 1 + Read 2 arm", index_arm, "r2", inferred=inferred),
                     Segment(f"{index_length}-base i7 reverse complement", "I" * index_length, "cbc", placeholder=True,
                             feature=feature("sample_index_i7", "sample_index", "unknown"),
                             inferred=index_inferred),
                     Segment("P7 reverse complement", p7_end, "p7", inferred=inferred)], name=name)
    primers = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["R2"])
    problems = sp.verify(lib, primers, required_roles=tuple(p.role for p in primers))
    if problems:
        raise ValueError(f"{name}: " + "; ".join(problems))
    return lib, primers


def _finish_pcr(graph, prefix, template, final, primers, endpoint, reviewed, *, label,
                forward_state=None, reverse_state=None, step_label=None,
                forward_label=None, reverse_label=None, polymerase_label=None,
                primer_source=None, primers_inferred=False, geometry=None,
                unavailable_reason=None):
    """Endpoint PCR with role-level default labels (no index or kit claim)."""
    forward = graph.input(prefix + "_forward", forward_label or label + " PCR forward primer",
                          MaterialKind.OLIGO, state=forward_state, source=primer_source,
                          inferred=primers_inferred)
    reverse = graph.input(prefix + "_reverse", reverse_label or label + " PCR reverse primer",
                          MaterialKind.OLIGO, state=reverse_state, source=primer_source,
                          inferred=primers_inferred)
    polymerase = graph.input(prefix + "_polymerase", polymerase_label or label + " PCR polymerase",
                             MaterialKind.CATALYST)
    library = graph.product(prefix + "_library", label + " sequencing library",
                            MaterialKind.LIBRARY, state=final)
    graph.react(prefix + "_pcr",
                step_label or ("Amplify the " + label + " library"),
                StepKind.PCR,
                ((template, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (library,), geometry=geometry or PCRGeometry.unavailable(
                    final, reason=unavailable_reason or "The final sequencing construct is source-backed, while a single universal PCR-binding geometry is not available for every kit revision."))
    graph.sequence(library, SequencingPlan(library=final, primers=tuple(primers)))
    reviewed.extend((prefix + "_forward", prefix + "_reverse", prefix + "_polymerase",
                     prefix + "_pcr", prefix + "_library"))
    return library
