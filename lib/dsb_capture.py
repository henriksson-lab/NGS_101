"""Bounded constructions for paired-end break capture libraries."""
from dataclasses import dataclass

import illumina as il
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
                            index_length: int = 6,
                            index_inferred: bool = False,
                            proximal_tag: str | None = None,
                            distal_tag: str | None = None) -> tuple[Construct, tuple]:
    """Build a two-end library while preserving inferred index and adapter provenance."""
    if tagged and (not proximal_tag or not distal_tag):
        raise ValueError("tagged break library needs both end tags")
    index_inferred = inferred or index_inferred
    left = Segment("P5 + Read 1 arm", il.TRUSEQ_P5_FULL, "p5", inferred=inferred)
    payload = ([Segment("proximal read tag", proximal_tag, "cbc",
                        feature=feature("proximal_break_end", "feature_barcode", "fixed")),
                role("break-to-random genomic DNA"),
                Segment("distal read tag", distal_tag, "cbc",
                        feature=feature("distal_fragment_end", "feature_barcode", "fixed"))]
               if tagged else [role("DNA beginning at repaired DSB")])
    lib = Construct([left, *payload, Segment("dA junction", "A", inferred=inferred),
                     Segment("Index 1 + Read 2 arm", il.INDEX1_PRIMER, "r2", inferred=inferred),
                     Segment(f"{index_length}-base i7 reverse complement", "I" * index_length, "cbc", placeholder=True,
                             feature=feature("sample_index_i7", "sample_index", "unknown"),
                             inferred=index_inferred),
                     Segment("P7 reverse complement", il.P7_RC, "p7", inferred=inferred)], name=name)
    primers = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["R2"])
    problems = sp.verify(lib, primers, required_roles=tuple(p.role for p in primers))
    if problems:
        raise ValueError(f"{name}: " + "; ".join(problems))
    return lib, primers


def _finish_pcr(graph, prefix, template, final, primers, endpoint, reviewed, *, label,
                forward_state=None, reverse_state=None, step_label=None):
    forward = graph.input(prefix + "_forward", label + " PCR forward primer",
                          MaterialKind.OLIGO, state=forward_state)
    reverse = graph.input(prefix + "_reverse", label + " PCR reverse primer",
                          MaterialKind.OLIGO, state=reverse_state)
    polymerase = graph.input(prefix + "_polymerase", label + " PCR polymerase",
                             MaterialKind.CATALYST)
    library = graph.product(prefix + "_library", label + " sequencing library",
                            MaterialKind.LIBRARY, state=final)
    graph.react(prefix + "_pcr",
                step_label or ("Amplify and index the " + label + " library"),
                StepKind.PCR,
                ((template, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (library,), geometry=PCRGeometry.unavailable(
                    final, reason="The final sequencing construct is source-backed, while a single universal PCR-binding geometry is not available for every kit revision."))
    graph.sequence(library, SequencingPlan(library=final, primers=tuple(primers)))
    reviewed.extend((prefix + "_forward", prefix + "_reverse", prefix + "_polymerase",
                     prefix + "_pcr", prefix + "_library"))
    return library
