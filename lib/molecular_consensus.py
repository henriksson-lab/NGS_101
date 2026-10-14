"""Reusable molecular states for UMI-family and strand-aware consensus methods.

The constructors keep two meanings separate: a single-strand UID identifies one
amplified template strand, whereas a duplex UID is copied into both strands before
amplification and can therefore reconnect Watson and Crick families.
"""
from __future__ import annotations

from dataclasses import dataclass

from batch_ngs import seg
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Scene, Segment, feature
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from protocol_model import molecule, role


def uid_segment(name: str, length: int, identity: str, *, group: str | None = None,
                part: str | None = None) -> Segment:
    """Construct a random molecular identifier with validated positive length."""
    if length < 1:
        raise ValueError("a molecular UID must contain at least one base")
    return seg(name, "N" * length, "umi", placeholder=True,
               feature=feature(identity, "umi", "random", group=group, part=part))


@dataclass(frozen=True)
class SingleStrandUID:
    """A primer-borne UID assigned before redundant target amplification."""

    length: int = 14

    def assignment(self) -> Scene:
        target = [seg("target body", "X" * 20, placeholder=True),
                  seg("target-primer site", "T" * 12)]
        primer = [seg("universal PCR handle", "H" * 14, placeholder=True),
                  uid_segment("random UID", self.length, "safe_seqs_uid"),
                  seg("target-specific primer", "A" * 12, placeholder=True)]
        sc = Scene(); sc.strand("template", target, label="one target strand")
        sc.anneal("UID primer", primer, to="template",
                  pair=("target-specific primer", "target-primer site"),
                  label="UID-assignment primer",
                  unpaired=("universal PCR handle", "random UID"))
        sc.arrow("UID primer", "limited-cycle UID assignment")
        sc.mark("UID primer", "random UID", "one UID defines one strand family")
        sc.labels("UID primer")
        return sc

    def family(self) -> Scene:
        parts = [seg("universal PCR handle", "H" * 14, placeholder=True),
                 uid_segment("copied UID", self.length, "safe_seqs_uid"),
                 seg("target amplicon", "X" * 32, placeholder=True),
                 seg("opposite universal handle", "K" * 14, placeholder=True)]
        sc = Scene.duplex(parts, label="redundant copies in one UID family")
        sc.mark("top", "copied UID", "group reads before family consensus")
        sc.note("top", "a variant is accepted only when it dominates the UID family")
        sc.labels("top")
        return sc


@dataclass(frozen=True)
class CopiedDuplexUID:
    """A one-ended UID copied across a duplex before library amplification."""

    length: int = 14

    def three_prime_adapter(self) -> Scene:
        parts = [seg("blunt dephosphorylated insert", "X" * 32, placeholder=True),
                 uid_segment("3-prime adapter UID", self.length,
                             "saferseqs_duplex_uid"),
                 seg("3-prime adapter arm", "A" * 14, placeholder=True)]
        sc = Scene(); sc.strand("ligated end", parts, label="one 3′ end of the original duplex")
        sc.junction("ligated end", "blunt dephosphorylated insert",
                    "3-prime adapter UID", "3′ UID-adapter ligation")
        sc.note("ligated end", "the opposite 3′ fragment end receives the same adapter class")
        sc.labels("ligated end")
        return sc

    def copied_barcode(self) -> Scene:
        parts = [seg("5-prime adapter", "H" * 14, placeholder=True),
                 seg("genomic insert", "X" * 32, placeholder=True),
                 uid_segment("copied duplex UID", self.length,
                             "saferseqs_duplex_uid"),
                 seg("3-prime adapter arm", "K" * 14, placeholder=True)]
        sc = Scene.duplex(parts, label="duplex-barcoded original molecule")
        sc.junction("top", "5-prime adapter", "genomic insert", "5′ adapter ligation")
        sc.junction("top", "genomic insert", "copied duplex UID", "3′ UID-adapter ligation")
        sc.mark("top", "copied duplex UID", "same molecular identity on Watson and Crick")
        sc.note("top", "extension copies the UID; 5′-adapter ligation seals the duplex")
        sc.labels("top")
        return sc

    def strand_families(self) -> Scene:
        sc = Scene()
        sc.strand("Watson family", [uid_segment("Watson UID", self.length,
                                                "saferseqs_duplex_uid"),
                                     seg("Watson target", "X" * 30, placeholder=True)],
                  label="Watson-specific hemi-nested PCR")
        sc.strand("Crick family", [uid_segment("Crick UID", self.length,
                                               "saferseqs_duplex_uid"),
                                    seg("Crick target", "Y" * 30, placeholder=True)],
                  label="Crick-specific hemi-nested PCR")
        sc.note("Watson family", "paired by the shared 14-nt UID after family consensus")
        sc.labels("Watson family"); sc.labels("Crick family")
        return sc


def safe_seqs_graph(authority: str, final: Construct, primers) -> ProtocolGraph:
    graph = ProtocolGraph("Safe-SeqS — single-strand UID-family consensus",
                          authority=authority)
    target = graph.input("target_dna", "target DNA containing assay primer sites",
                         MaterialKind.BIOLOGICAL,
                         state=molecule("target DNA", (role("target strands"),)), primary=True)
    uid_primer = graph.input("uid_primer", "target-specific primer carrying a random UID",
                             MaterialKind.OLIGO)
    opposing = graph.input("opposing_primer", "opposing target-specific primer",
                           MaterialKind.OLIGO)
    assignment_pol = graph.input("assignment_polymerase", "UID-assignment PCR polymerase",
                                 MaterialKind.CATALYST)
    assigned = graph.product("uid_assigned", "UID-assigned target strands",
                             MaterialKind.MOLECULAR, state=SingleStrandUID(14).assignment())
    graph.react("assign_uid", "Assign a random UID in two target-specific PCR cycles",
                StepKind.PCR,
                ((target, InputRole.TEMPLATE), (uid_primer, InputRole.STRUCTURAL),
                 (opposing, InputRole.STRUCTURAL),
                 (assignment_pol, InputRole.CATALYST)), (assigned,),
                geometry=PCRGeometry.unavailable(
                    final, reason="The assay-specific target template and opposing site vary by panel."))
    universal_f = graph.input("universal_forward", "universal-handle forward primer",
                              MaterialKind.OLIGO)
    universal_r = graph.input("universal_reverse", "universal-handle reverse primer",
                              MaterialKind.OLIGO)
    family_pol = graph.input("family_polymerase", "UID-family amplification polymerase",
                             MaterialKind.CATALYST)
    library = graph.product("library", "redundantly amplified Safe-SeqS UID families",
                            MaterialKind.LIBRARY, state=final)
    graph.react("amplify_families", "Amplify through the universal handles",
                StepKind.PCR,
                ((assigned, InputRole.TEMPLATE), (universal_f, InputRole.STRUCTURAL),
                 (universal_r, InputRole.STRUCTURAL), (family_pol, InputRole.CATALYST)),
                (library,), geometry=PCRGeometry.unavailable(
                    final, reason="The final TruSeq shell is exact, while target-specific insert coordinates vary."))
    graph.sequence(library, SequencingPlan(library=final, primers=tuple(primers)))
    graph.declare_reviewed_scope([
        "target_dna", "uid_primer", "opposing_primer", "assignment_polymerase",
        "assign_uid", "universal_forward", "universal_reverse", "family_polymerase",
        "amplify_families", "library"], authority)
    return graph.require_publishable()


def saferseqs_graph(authority: str, final: Construct, limitation: str, *,
                    primers: tuple = (), read_uses: dict | None = None) -> ProtocolGraph:
    graph = ProtocolGraph("SaferSeqS — copied duplex UID with strand-specific enrichment",
                          authority=authority)
    dna = graph.input("genomic_dna", "end-repaired genomic duplex DNA",
                      MaterialKind.BIOLOGICAL,
                      state=molecule("end-repaired genomic DNA", (role("original Watson strand"),
                                                                   role("original Crick strand"))),
                      primary=True)
    phosphatase = graph.input("phosphatase", "fragment-end phosphatase",
                              MaterialKind.CATALYST)
    dephosphorylated = graph.product("dephosphorylated", "blunt dephosphorylated duplex DNA",
                                    MaterialKind.MOLECULAR, state=dna.state)
    graph.react("dephosphorylate", "Dephosphorylate repaired fragment ends",
                StepKind.DEPHOSPHORYLATION,
                ((dna, InputRole.TEMPLATE), (phosphatase, InputRole.CATALYST)),
                (dephosphorylated,))
    uid_adapter = graph.input("uid_adapter", "random-UID-bearing 3-prime adapter",
                              MaterialKind.ADAPTER,
                              source=authority + " Supplementary Table 7")
    ligase3 = graph.input("three_prime_ligase", "3-prime adapter ligase",
                          MaterialKind.CATALYST)
    uid_ligated = graph.product("uid_ligated", "3-prime UID-adapter-ligated duplex",
                                MaterialKind.MOLECULAR, state=CopiedDuplexUID(14).three_prime_adapter())
    graph.react("ligate_uid", "Ligate the UID-bearing 3-prime adapter",
                StepKind.LIGATION,
                ((dephosphorylated, InputRole.TEMPLATE), (uid_adapter, InputRole.STRUCTURAL),
                 (ligase3, InputRole.CATALYST)), (uid_ligated,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The published oligo pair is transcribed on the protocol page; its "
                    "two-ended placement on the genomic duplex is not yet encoded as an exact geometry."))
    extension = graph.input("extension_polymerase", "UID-copying DNA polymerase",
                            MaterialKind.CATALYST)
    copied = graph.product("uid_copied", "duplex with the UID copied into both strands",
                           MaterialKind.MOLECULAR, state=CopiedDuplexUID(14).copied_barcode())
    graph.react("copy_uid", "Extend through the UID to copy it into both strands",
                StepKind.EXTENSION,
                ((uid_ligated, InputRole.TEMPLATE), (extension, InputRole.CATALYST)),
                (copied,))
    five_adapter = graph.input("five_prime_adapter", "5-prime library adapter",
                               MaterialKind.ADAPTER,
                               source=authority + " Supplementary Table 7")
    ligase5 = graph.input("five_prime_ligase", "5-prime adapter ligase",
                          MaterialKind.CATALYST)
    sealed = graph.product("sealed", "duplex-barcoded library precursor",
                           MaterialKind.MOLECULAR, state=copied.state)
    graph.react("ligate_five_prime", "Ligate the 5-prime adapter",
                StepKind.LIGATION,
                ((copied, InputRole.TEMPLATE), (five_adapter, InputRole.STRUCTURAL),
                 (ligase5, InputRole.CATALYST)), (sealed,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The published 5-prime adapter is transcribed on the protocol page; its "
                    "two-ended placement on the genomic duplex is not yet encoded as an exact geometry."))
    watson = graph.product("watson_template", "Watson-oriented target molecules",
                           MaterialKind.MOLECULAR, state=sealed.state)
    crick = graph.product("crick_template", "Crick-oriented target molecules",
                          MaterialKind.MOLECULAR, state=sealed.state)
    graph.react("orient", "Partition into Watson- and Crick-specific enrichments",
                StepKind.SPLIT, ((sealed, InputRole.TEMPLATE),), (watson, crick),
                changes_molecule=False)
    polymerase = graph.input("nested_polymerase", "hemi-nested PCR polymerase",
                             MaterialKind.CATALYST)
    reviewed = ["genomic_dna", "phosphatase", "dephosphorylate", "uid_adapter",
                "three_prime_ligase", "ligate_uid", "extension_polymerase", "copy_uid",
                "five_prime_adapter", "five_prime_ligase", "ligate_five_prime",
                "orient", "nested_polymerase"]
    for prefix, template in (("watson", watson), ("crick", crick)):
        forward = graph.input(prefix + "_forward", prefix.title() + " hemi-nested forward primer",
                              MaterialKind.OLIGO, source=authority + " Supplementary Tables 5 and 7")
        reverse = graph.input(prefix + "_reverse", prefix.title() + " hemi-nested reverse primer",
                              MaterialKind.OLIGO, source=authority + " Supplementary Tables 5 and 7")
        library = graph.product(prefix + "_library", prefix.title() + " SaferSeqS library",
                                MaterialKind.LIBRARY, state=final)
        graph.react(prefix + "_pcr", "Hemi-nested PCR of the " + prefix + " strand",
                    StepKind.PCR,
                    ((template, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                     (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                    (library,), geometry=PCRGeometry.unavailable(final, reason=limitation))
        graph.sequence(library, SequencingPlan(
            library=final if primers else None, primers=primers,
            acquired_roles=("Read 1", "Index 1 (i7)", "Index 2 (i5)", "Read 2"),
            read_uses=read_uses or {}, limitation=limitation))
        reviewed.extend((forward.id, reverse.id, prefix + "_pcr", library.id))
    graph.declare_reviewed_scope(reviewed, authority)
    return graph.require_publishable()
