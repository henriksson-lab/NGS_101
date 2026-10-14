"""Source-backed page specifications for targeted 3C and Repli-seq methods."""
from __future__ import annotations

from dataclasses import dataclass
from types import SimpleNamespace

import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry
from batch_ngs import seg, truseq_library
from chemdraw import Construct, MolecularState, Scene, Segment, Workflow, feature
from circular import CircularizationGeometry
from chromosome_capture import (Circular4CTemplate, ThreeCLibrary,
                                inverse_pcr_product, multiway_contact, region)
from replication_labeling import CopyNumberReplication, NascentPulse
from protocol_model import da_tailed
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)


def role(name: str, *, tag=None, molecular_feature=None) -> Segment:
    return Segment(name, f"[{name}]", tag, placeholder=True,
                   feature=molecular_feature)


def molecule(label: str, names: tuple[str, ...], *,
             junctions: tuple[tuple[str, str, str], ...] = ()) -> Scene:
    parts = [region(name) for name in names]
    sc = Scene(); sc.strand("molecule", parts, label=label)
    for left, right, text in junctions:
        sc.junction("molecule", left, right, text)
    sc.labels("molecule")
    return sc


def chain(initial: Scene, steps: tuple[tuple[str, Scene, str], ...]) -> Workflow:
    w = Workflow(MolecularState("Starting material", tuple(initial.rows())))
    for action, state, note in steps:
        w.react(action, state.rows(), name=action, note=note)
    return w


def _three_c(enzyme: str, names: tuple[str, ...]) -> ThreeCLibrary:
    return ThreeCLibrary(enzyme, names)


def workflow_for(key: str) -> Workflow:
    if key == "4c-seq":
        contact = _three_c("DpnII", ("viewpoint fragment", "contacted reporter fragment"))
        circle = Circular4CTemplate(contact, "Csp6I", "viewpoint fragment")
        return chain(molecule("formaldehyde-fixed chromatin", ("viewpoint", "proximal chromatin")), (
            ("Digest chromatin with the primary restriction enzyme and proximity-ligate",
             contact.scene(), "The first junction records proximity in fixed nuclei."),
            ("Reverse crosslinks, digest with a second enzyme and ligate dilute DNA",
             circle.scene(), "The second ligation recircularizes shorter 4C templates."),
            ("Inverse-PCR outward from the chosen viewpoint", inverse_pcr_product(circle),
             "A two-stage PCR adds sample and Illumina sequences to viewpoint-specific products."),
        ))
    if key == "ng-capture-c":
        contact = _three_c("DpnII or NlaIII", ("captured viewpoint fragment", "reporter fragment"))
        adapted = molecule("indexed 3C sequencing library",
                           ("left Illumina adapter", "viewpoint", "contact junction", "reporter", "right Illumina adapter"),
                           junctions=(("left Illumina adapter", "viewpoint", "adapter ligation"),
                                      ("viewpoint", "contact junction", "proximity ligation"),
                                      ("reporter", "right Illumina adapter", "adapter ligation")))
        return chain(molecule("formaldehyde-fixed chromatin", ("viewpoint", "proximal chromatin")), (
            ("Digest and proximity-ligate to make a 3C library", contact.scene(),
             "No biotin fill-in is used to select the contact junction."),
            ("Sonicate to about 200 bp, end-repair and add indexed Illumina adapters", adapted,
             "Random sonication ends distinguish original molecules from PCR duplicates."),
            ("Perform two sequential captures with biotinylated viewpoint oligos", adapted,
             "Streptavidin enrichment retains fragments containing the selected viewpoint."),
        ))
    if key == "umi-4c":
        contact = _three_c("DpnII", ("viewpoint fragment", "reporter fragment"))
        one_end = molecule("single-end-adapter 3C fragment",
                           ("viewpoint", "contact junction", "reporter", "random sonication endpoint",
                            "single-end Illumina adapter"),
                           junctions=(("viewpoint", "contact junction", "proximity ligation"),
                                      ("random sonication endpoint", "single-end Illumina adapter", "adapter ligation")))
        return chain(molecule("formaldehyde-fixed chromatin", ("viewpoint", "proximal chromatin")), (
            ("Digest with DpnII, proximity-ligate and reverse crosslinks", contact.scene(),
             "This is a linear 3C library; UMI-4C does not perform 4C recircularization."),
            ("Sonicate, end-repair and add one 3-prime dA to the ligatable end",
             da_tailed("dA-tailed UMI-4C fragment", "viewpoint-contact fragment with random sonication endpoint"),
             "The random shearing coordinate is already fixed before adapter ligation."),
            ("Ligate the single-end adapter", one_end,
             "The random shearing coordinate acts as the molecular identifier; it is not a synthesized random barcode."),
            ("Amplify with two nested viewpoint PCRs", one_end,
             "The nested viewpoint primer and universal adapter primer create the sequencing library."),
        ))
    if key == "tri-c":
        contact = multiway_contact("NlaIII", "viewpoint fragment",
                                   ("reporter fragment 1", "reporter fragment 2"))
        adapted = molecule("adapter-ligated multiway 3C fragment",
                           ("left TruSeq adapter", "random sonication endpoint", "viewpoint",
                            "reporter 1", "reporter 2", "right TruSeq adapter"),
                           junctions=(("left TruSeq adapter", "random sonication endpoint", "adapter ligation"),
                                      ("viewpoint", "reporter 1", "proximity ligation"),
                                      ("reporter 1", "reporter 2", "proximity ligation"),
                                      ("reporter 2", "right TruSeq adapter", "adapter ligation")))
        return chain(molecule("formaldehyde-fixed chromatin", ("viewpoint", "two proximal regions")), (
            ("Digest with NlaIII and proximity-ligate in nuclei", contact.scene("multiway 3C concatemer"),
             "At least two ligation junctions in one molecule encode a concurrent contact."),
            ("Reverse crosslinks and sonicate gently to about 450 bp", adapted,
             "Longer fragments preserve multiple short NlaIII products and their junctions."),
            ("Add TruSeq adapters, index, then perform double viewpoint capture", adapted,
             "A biotinylated oligo targets the middle of the short viewpoint restriction fragment."),
        ))
    if key == "el-repli-seq":
        pulse = NascentPulse("BrdU", ("early S phase", "late S phase"))
        adapted = molecule("adapter-ligated nascent DNA",
                           ("left Illumina adapter", "BrdU-containing nascent DNA", "right Illumina adapter"),
                           junctions=(("left Illumina adapter", "BrdU-containing nascent DNA", "adapter ligation"),
                                      ("BrdU-containing nascent DNA", "right Illumina adapter", "adapter ligation")))
        return chain(molecule("asynchronously cycling cells", ("replicating chromosomes",)), (
            ("Pulse-label nascent DNA with BrdU", pulse.labelled_scene(),
             "Only DNA synthesized during the two-hour pulse contains BrdU."),
            ("Fix and sort early- and late-S fractions by DNA content", pulse.sorted_scene(),
             "Early and late fractions remain separate libraries."),
            ("Purify, sonicate, end-repair and ligate universal Illumina adapters", adapted,
             "Adapters are installed before immunoprecipitation."),
            ("Denature, immunoprecipitate BrdU DNA, then index and PCR-amplify", adapted,
             "Sequencing coverage is compared as an early/late ratio."),
        ))
    if key == "screpli-seq":
        dosage = CopyNumberReplication("mid-S-phase", "SeqPlex DOP-PCR whole-genome amplification")
        amplified = molecule("single-cell WGA product",
                             ("left SEQXE amplification handle", "cell genome copy",
                              "right SEQXE amplification handle"))
        adapted = molecule("single-cell Illumina library",
                           ("P5 / sample index", "WGA genomic fragment", "sample index / P7"),
                           junctions=(("P5 / sample index", "WGA genomic fragment", "adapter ligation"),
                                      ("WGA genomic fragment", "sample index / P7", "adapter ligation")))
        return chain(molecule("one flow-sorted S-phase cell", ("partly replicated genome",)), (
            ("Use copy number as the replication mark", dosage.states(),
             "Replicated domains have twice the DNA dosage of unreplicated domains in the same cell."),
            ("Whole-genome amplify each cell by SeqPlex DOP-PCR", amplified,
             "G1 cells processed in parallel provide the amplification and mappability reference."),
            ("Fragment WGA DNA, end-repair and add one 3-prime dA to each end",
             da_tailed("dA-tailed single-cell WGA fragment", "WGA genomic fragment"),
             "The amplified genomic fragment is shown before adapter ligation."),
            ("Ligate Illumina adapters", adapted,
             "The protocol permits KAPA LTP or NEBNext Ultra II library construction."),
            ("PCR-index each single-cell library and sequence", adapted,
             "The index assigns reads to one cell; it is not an inline cell barcode."),
        ))
    raise KeyError(key)


@dataclass(frozen=True)
class Definition:
    title: str
    doi: str
    citation: str
    summary: str
    final_parts: tuple[Segment, ...]
    sequencing: str
    standard_truseq: bool
    caveat: str = ""


ENDPOINT = Segment("random sonication endpoint identifier",
                   "[random shearing endpoint used as UMI]", "umi", placeholder=True,
                   feature=feature("sonication_endpoint", "umi", "random",
                                   note="physical shearing coordinate, not a synthesized barcode"))
SAMPLE_INDEX = Segment("PCR sample index", "I" * 8, "cbc", placeholder=True,
                       feature=feature("sample_index", "sample_index", "unknown"))

DEFS = {
    "4c-seq": Definition(
        "4C-seq — circular chromosome conformation capture sequencing",
        "10.1016/j.ymeth.2019.07.014", "Krijger et al., Methods (2020)",
        "Two restriction digests and two ligations convert a viewpoint-centred 3C contact into a circle that can be opened by outward-facing inverse PCR.",
        (role("viewpoint PCR / Illumina end"), SAMPLE_INDEX, role("viewpoint fragment"),
         role("captured reporter"), role("second viewpoint PCR / Illumina end")),
        "A viewpoint-specific reading primer enters the captured reporter beside the primary restriction junction; sample and flow-cell sequences are added by the two PCR stages.", False,
        "The genomic portions of viewpoint primers are locus-specific. Published universal tails and second-round primer sequences are listed below; the chosen genomic portions and index depend on the experiment."),
    "ng-capture-c": Definition(
        "NG Capture-C — next-generation Capture-C", "10.1038/nmeth.3664",
        "Davies et al., Nature Methods (2016)",
        "A sonicated, indexed 3C library undergoes two rounds of capture with biotinylated DNA oligos directed to selected restriction fragments.",
        (role("captured viewpoint"), role("proximity-ligation junction"), role("reporter fragment")),
        "Paired-end Illumina reads recover the viewpoint-containing end and its ligated reporter; random sonication coordinates identify PCR duplicates.", True),
    "umi-4c": Definition(
        "UMI-4C — quantitative targeted chromosome contacts", "10.1038/nmeth.3922",
        "Schwartzman et al., Nature Methods (2016)",
        "Random shearing supplies molecule-specific endpoints before single-end adapter ligation and nested viewpoint PCR, enabling quantitative targeted 3C profiling.",
        (role("nested viewpoint / Illumina end"), SAMPLE_INDEX, role("viewpoint"),
         role("contact reporter"), ENDPOINT, role("universal adapter end")),
        "Paired-end sequencing recovers the viewpoint/contact molecule and the random shearing endpoint used to collapse PCR duplicates.", False,
        "The published ‘UMI’ is the genomic sonication endpoint, not a fixed-length random oligonucleotide; the exact locus-specific bait primers vary by experiment."),
    "tri-c": Definition(
        "Tri-C — targeted multiway chromatin contacts", "10.1038/s41588-018-0253-2",
        "Oudelaar et al., Nature Genetics (2018)",
        "Small NlaIII fragments, gentle sonication and viewpoint capture preserve multiple proximity-ligation junctions in one sequenced molecule.",
        (role("viewpoint"), role("reporter 1"), role("second proximity junction"),
         role("reporter 2"), ENDPOINT),
        "150-bp paired-end Illumina reads identify the viewpoint and two or more reporters within one roughly 450-bp captured 3C fragment.", True),
    "el-repli-seq": Definition(
        "E/L Repli-seq — early/late replication timing", "10.1038/nprot.2017.148",
        "Marchal et al., Nature Protocols (2018)",
        "BrdU-labelled nascent DNA is separated into early- and late-S populations, adapter-ligated, immunoprecipitated and sequenced as two quantitative fractions.",
        (role("BrdU-enriched genomic fragment"),),
        "Ordinary Illumina genomic reads measure coverage in separate early- and late-S libraries; replication timing is their normalized log2 ratio.", True),
    "screpli-seq": Definition(
        "scRepli-seq — single-cell DNA replication sequencing", "10.1038/s41596-020-0378-5",
        "Miura et al., Nature Protocols (2020)",
        "A single S-phase cell is whole-genome amplified and sequenced; twofold copy-number differences distinguish replicated from unreplicated domains.",
        (role("SeqPlex-amplified genomic fragment"),),
        "Each PCR index identifies one single-cell library. Genomic Read 1/Read 2 coverage, normalized to G1 controls, reports replicated and unreplicated domains.", True),
}


def _conformation_graph(key: str, definition, final: Construct,
                        primers: tuple[sp.SeqPrimer, ...]) -> ProtocolGraph:
    """Typed graphs for 4C/Capture-C derivatives with distinct contact ligation."""
    source = f"{definition.citation}; doi:{definition.doi}"
    graph = ProtocolGraph(definition.title, authority=source)
    fixed = graph.input("fixed_chromatin", "formaldehyde-fixed chromatin",
                        MaterialKind.BIOLOGICAL,
                        state=molecule("formaldehyde-fixed chromatin",
                                       ("viewpoint", "proximal chromatin")), primary=True)
    restriction_name = {"4c-seq": "DpnII", "ng-capture-c": "DpnII or NlaIII",
                        "umi-4c": "DpnII", "tri-c": "NlaIII"}[key]
    restriction = graph.input("primary_restriction_enzyme", restriction_name,
                              MaterialKind.CATALYST)
    digested = graph.product(
        "restriction_fragments", "restriction-digested fixed chromatin",
        MaterialKind.MOLECULAR,
        state=molecule("restriction fragments", ("viewpoint fragment",
                                                   "proximal reporter fragments")))
    graph.react("primary_digest", f"Digest fixed chromatin with {restriction_name}",
                StepKind.DIGESTION,
                ((fixed, InputRole.TEMPLATE),
                 (restriction, InputRole.CATALYST)), (digested,))
    ligase = graph.input("proximity_ligase", "proximity-ligation DNA ligase",
                         MaterialKind.CATALYST)
    if key == "tri-c":
        contact_model = multiway_contact(
            "NlaIII", "viewpoint fragment", ("reporter fragment 1", "reporter fragment 2"))
        contact_state = contact_model.scene("multiway 3C concatemer")
        contact_label = "multiway proximity-ligation concatemer"
    else:
        names = (("captured viewpoint fragment", "reporter fragment")
                 if key == "ng-capture-c" else
                 ("viewpoint fragment", "reporter fragment"))
        contact_model = _three_c(restriction_name, names)
        contact_state = contact_model.scene()
        contact_label = "proximity-ligated 3C contact"
    contact = graph.product("contact_product", contact_label,
                            MaterialKind.MOLECULAR, state=contact_state)
    graph.react("proximity_ligation", "Proximity-ligate restriction fragments",
                StepKind.PROXIMITY_LIGATION,
                ((digested, InputRole.TEMPLATE),
                 (ligase, InputRole.CATALYST)), (contact,))
    proteinase = graph.input("crosslink_reversal", "proteinase/cross-link reversal condition",
                             MaterialKind.CONDITION)
    deproteinized = graph.product("deproteinized_3c_dna", "deproteinized 3C DNA",
                                  MaterialKind.MOLECULAR, state=contact_state)
    graph.react("reverse_crosslinks", "Reverse cross-links and recover 3C DNA",
                StepKind.DEPROTEINIZATION,
                ((contact, InputRole.TEMPLATE),
                 (proteinase, InputRole.CONDITION)), (deproteinized,))
    reviewed = [
        "fixed_chromatin", "primary_restriction_enzyme", "primary_digest",
        "proximity_ligase", "proximity_ligation", "crosslink_reversal",
        "reverse_crosslinks",
    ]

    if key == "4c-seq":
        second = graph.input("secondary_restriction_enzyme", "Csp6I",
                             MaterialKind.CATALYST)
        secondary = graph.product("secondary_fragments", "Csp6I-digested 3C products",
                                  MaterialKind.MOLECULAR, state=contact_state)
        graph.react("secondary_digest", "Digest 3C products with Csp6I",
                    StepKind.DIGESTION,
                    ((deproteinized, InputRole.TEMPLATE),
                     (second, InputRole.CATALYST)), (secondary,))
        circle_model = Circular4CTemplate(contact_model, "Csp6I", "viewpoint fragment")
        circular = graph.product("circular_template", "recircularized 4C template",
                                 MaterialKind.MOLECULAR, state=circle_model.scene())
        graph.react("recircularize", "Ligate dilute DNA to form 4C circles",
                    StepKind.CIRCULARIZATION,
                    ((secondary, InputRole.TEMPLATE),
                     (ligase, InputRole.CATALYST)), (circular,),
                    geometry=CircularizationGeometry.known(
                        "double-strand intramolecular ligation",
                        statement="The 4C workflow dilute-ligates the secondary restriction fragments into circles."))
        inverse_forward = graph.input("inverse_forward_primer", "outward viewpoint primer 1",
                                      MaterialKind.OLIGO)
        inverse_reverse = graph.input("inverse_reverse_primer", "outward viewpoint primer 2",
                                      MaterialKind.OLIGO)
        polymerase = graph.input("pcr_polymerase", "PCR polymerase",
                                 MaterialKind.CATALYST)
        first_amplicon = graph.product("inverse_pcr_product", "viewpoint inverse-PCR product",
                                      MaterialKind.MOLECULAR,
                                      state=inverse_pcr_product(circle_model))
        graph.react("inverse_pcr", "Inverse-PCR outward from the viewpoint",
                    StepKind.PCR,
                    ((circular, InputRole.TEMPLATE),
                     (inverse_forward, InputRole.STRUCTURAL),
                     (inverse_reverse, InputRole.STRUCTURAL),
                     (polymerase, InputRole.CATALYST)), (first_amplicon,),
                    geometry=PCRGeometry.unavailable(
                        final, reason="Viewpoint primer sequences are locus-specific."))
        handle_forward = graph.input("handle_forward_primer", "Illumina-handle PCR forward primer",
                                     MaterialKind.OLIGO, inferred=True)
        handle_reverse = graph.input("handle_reverse_primer", "indexed Illumina-handle PCR reverse primer",
                                     MaterialKind.OLIGO, inferred=True)
        library = graph.product("sequencing_library", "4C-seq sequencing library",
                                MaterialKind.LIBRARY, state=final, inferred=True)
        graph.react("handle_pcr", "Add sample and Illumina sequences by PCR",
                    StepKind.PCR,
                    ((first_amplicon, InputRole.TEMPLATE),
                     (handle_forward, InputRole.STRUCTURAL),
                     (handle_reverse, InputRole.STRUCTURAL),
                     (polymerase, InputRole.CATALYST)), (library,),
                    geometry=PCRGeometry.unavailable(
                        final, reason="The exact viewpoint/handle combination varies by experiment."))
        graph.sequence(library, SequencingPlan(limitation=definition.sequencing))
        reviewed.extend((
            "secondary_restriction_enzyme", "secondary_digest", "recircularize",
            "inverse_forward_primer", "inverse_reverse_primer", "pcr_polymerase",
            "inverse_pcr", "handle_forward_primer", "handle_reverse_primer",
            "handle_pcr", "sequencing_library",
        ))
    else:
        sonication = graph.input("sonication", "DNA sonication condition",
                                 MaterialKind.CONDITION)
        sonicated = graph.product("sonicated_3c_dna", "sonicated 3C DNA",
                                  MaterialKind.MOLECULAR, state=deproteinized.state)
        graph.react("sonicate", "Sonicate recovered 3C DNA", StepKind.FRAGMENTATION,
                    ((deproteinized, InputRole.TEMPLATE),
                     (sonication, InputRole.CONDITION)), (sonicated,))
        repair = graph.input("end_repair_enzymes", "end-repair enzyme mix",
                             MaterialKind.CATALYST)
        repaired = graph.product("end_repaired_dna", "end-repaired sonicated 3C DNA",
                                 MaterialKind.MOLECULAR, state=sonicated.state)
        graph.react("end_repair", "End-repair sonicated 3C DNA", StepKind.END_REPAIR,
                    ((sonicated, InputRole.TEMPLATE),
                     (repair, InputRole.CATALYST)), (repaired,))
        tailing = graph.input("a_tailing_enzyme", "dA-tailing enzyme",
                              MaterialKind.CATALYST)
        tailed_state = da_tailed("dA-tailed 3C fragment", "viewpoint-contact fragment")
        tailed = graph.product("da_tailed_dna", "dA-tailed 3C DNA",
                               MaterialKind.MOLECULAR, state=tailed_state)
        graph.react("da_tail", "Add one 3-prime dA to each end", StepKind.TAILING,
                    ((repaired, InputRole.TEMPLATE),
                     (tailing, InputRole.CATALYST)), (tailed,))
        adapters = graph.input("illumina_adapters", "Illumina library adapters",
                               MaterialKind.ADAPTER, inferred=True)
        adapter_ligase = graph.input("adapter_ligase", "adapter DNA ligase",
                                     MaterialKind.CATALYST)
        adapted = graph.product("adapter_ligated_dna", "adapter-ligated 3C DNA",
                                MaterialKind.MOLECULAR, state=final, inferred=True)
        graph.react("adapter_ligation", "Ligate Illumina adapters",
                    StepKind.LIGATION,
                    ((tailed, InputRole.TEMPLATE),
                     (adapters, InputRole.STRUCTURAL),
                     (adapter_ligase, InputRole.CATALYST)), (adapted,),
                    geometry=AdapterLigationGeometry.unavailable(
                        "The defining protocol identifies its library kit or adapter role but does not print one immutable adapter-end reagent."))
        reviewed.extend((
            "sonication", "sonicate", "end_repair_enzymes", "end_repair",
            "a_tailing_enzyme", "da_tail", "illumina_adapters", "adapter_ligase",
            "adapter_ligation",
        ))
        polymerase = graph.input("pcr_polymerase", "library PCR polymerase",
                                 MaterialKind.CATALYST)

        if key == "umi-4c":
            first_forward = graph.input("nested1_viewpoint_primer", "first nested viewpoint primer",
                                        MaterialKind.OLIGO)
            first_reverse = graph.input("nested1_adapter_primer", "first universal-adapter primer",
                                        MaterialKind.OLIGO)
            nested = graph.product("nested_amplicon", "first nested UMI-4C amplicon",
                                   MaterialKind.MOLECULAR, state=adapted.state)
            graph.react("nested_pcr1", "First nested viewpoint PCR", StepKind.PCR,
                        ((adapted, InputRole.TEMPLATE),
                         (first_forward, InputRole.STRUCTURAL),
                         (first_reverse, InputRole.STRUCTURAL),
                         (polymerase, InputRole.CATALYST)), (nested,),
                        geometry=PCRGeometry.unavailable(
                            final, reason="The viewpoint primer is locus-specific."))
            second_forward = graph.input("nested2_viewpoint_primer", "indexed nested viewpoint primer",
                                         MaterialKind.OLIGO)
            second_reverse = graph.input("nested2_adapter_primer", "indexed universal-adapter primer",
                                         MaterialKind.OLIGO)
            library = graph.product("sequencing_library", "UMI-4C sequencing library",
                                    MaterialKind.LIBRARY, state=final, inferred=True)
            graph.react("nested_pcr2", "Second nested viewpoint PCR", StepKind.PCR,
                        ((nested, InputRole.TEMPLATE),
                         (second_forward, InputRole.STRUCTURAL),
                         (second_reverse, InputRole.STRUCTURAL),
                         (polymerase, InputRole.CATALYST)), (library,),
                        geometry=PCRGeometry.unavailable(
                            final, reason="The indexed viewpoint primer remains experiment-specific."))
            graph.sequence(library, SequencingPlan(limitation=definition.sequencing))
            reviewed.extend((
                "pcr_polymerase", "nested1_viewpoint_primer", "nested1_adapter_primer",
                "nested_pcr1", "nested2_viewpoint_primer", "nested2_adapter_primer",
                "nested_pcr2", "sequencing_library",
            ))
        else:
            forward = graph.input("index_forward_primer", "indexed library PCR forward primer",
                                  MaterialKind.OLIGO)
            reverse = graph.input("index_reverse_primer", "indexed library PCR reverse primer",
                                  MaterialKind.OLIGO)
            indexed = graph.product("indexed_library", "indexed 3C library",
                                    MaterialKind.MOLECULAR, state=final)
            graph.react("index_pcr", "Index and amplify the 3C library", StepKind.PCR,
                        ((adapted, InputRole.TEMPLATE),
                         (forward, InputRole.STRUCTURAL),
                         (reverse, InputRole.STRUCTURAL),
                         (polymerase, InputRole.CATALYST)), (indexed,),
                        geometry=PCRGeometry.unavailable(
                            final, reason="Index choice varies while the final TruSeq endpoint is canonical."))
            capture1 = graph.input("capture_oligo_round1", "biotinylated viewpoint capture oligo round 1",
                                   MaterialKind.SELECTION)
            selected1 = graph.product("capture_round1", "first-round viewpoint-enriched library",
                                      MaterialKind.MOLECULAR, state=final)
            reject1 = graph.product("capture_round1_reject", "uncaptured round-1 library",
                                    MaterialKind.WASTE)
            graph.react("capture_round1_step", "First viewpoint capture",
                        StepKind.AFFINITY_CAPTURE,
                        ((indexed, InputRole.TEMPLATE),
                         (capture1, InputRole.SELECTION)), (selected1, reject1))
            graph.discard(reject1, "does not hybridize to the selected viewpoint probe")
            capture2 = graph.input("capture_oligo_round2", "biotinylated viewpoint capture oligo round 2",
                                   MaterialKind.SELECTION)
            library = graph.product("sequencing_library", definition.title + " captured library",
                                    MaterialKind.LIBRARY, state=final)
            reject2 = graph.product("capture_round2_reject", "uncaptured round-2 library",
                                    MaterialKind.WASTE)
            graph.react("capture_round2_step", "Second viewpoint capture",
                        StepKind.AFFINITY_CAPTURE,
                        ((selected1, InputRole.TEMPLATE),
                         (capture2, InputRole.SELECTION)), (library, reject2))
            graph.discard(reject2, "does not survive the second viewpoint capture")
            graph.sequence(library, SequencingPlan(
                final, tuple(primers), tuple(p.role for p in primers)))
            reviewed.extend((
                "pcr_polymerase", "index_forward_primer", "index_reverse_primer",
                "index_pcr", "capture_oligo_round1", "capture_round1_step",
                "capture_oligo_round2", "capture_round2_step", "sequencing_library",
            ))
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _replication_graph(key: str, definition, final: Construct,
                       primers: tuple[sp.SeqPrimer, ...]) -> ProtocolGraph:
    """Typed replication-timing graphs retaining every comparison library."""
    source = f"{definition.citation}; doi:{definition.doi}"
    graph = ProtocolGraph(definition.title, authority=source)
    reviewed: list[str] = []
    plan = SequencingPlan(final, tuple(primers), tuple(p.role for p in primers))
    adapter_geometry = AdapterLigationGeometry.unavailable(
        "The protocol identifies the supported Illumina library preparation but does not print one immutable adapter-end reagent.")

    def finish_short_read_branch(prefix: str, template, label: str, *, wga: bool) -> None:
        current = template
        if wga:
            dop_forward = graph.input(prefix + "_dop_forward", label + " DOP-PCR forward primer",
                                      MaterialKind.OLIGO)
            dop_reverse = graph.input(prefix + "_dop_reverse", label + " DOP-PCR reverse primer",
                                      MaterialKind.OLIGO)
            wga_polymerase = graph.input(prefix + "_wga_polymerase", label + " WGA polymerase",
                                         MaterialKind.CATALYST)
            amplified = graph.product(prefix + "_wga", label + " SeqPlex WGA product",
                                      MaterialKind.MOLECULAR,
                                      state=molecule(label + " WGA product", (
                                          "left SEQXE amplification handle",
                                          "genome copy",
                                          "right SEQXE amplification handle")))
            graph.react(prefix + "_wga_step", "Whole-genome amplify " + label,
                        StepKind.PCR,
                        ((current, InputRole.TEMPLATE),
                         (dop_forward, InputRole.STRUCTURAL),
                         (dop_reverse, InputRole.STRUCTURAL),
                         (wga_polymerase, InputRole.CATALYST)), (amplified,),
                        geometry=PCRGeometry.unavailable(
                            final, reason="SeqPlex uses a degenerate DOP-PCR family rather than one fixed amplicon."))
            current = amplified
            reviewed.extend((prefix + "_dop_forward", prefix + "_dop_reverse",
                             prefix + "_wga_polymerase", prefix + "_wga_step"))
        else:
            extraction = graph.input(prefix + "_dna_extraction", label + " DNA purification condition",
                                     MaterialKind.CONDITION)
            purified = graph.product(prefix + "_purified_dna", label + " purified nascent DNA",
                                     MaterialKind.MOLECULAR, state=current.state)
            graph.react(prefix + "_extract", "Purify " + label + " DNA",
                        StepKind.EXTRACTION,
                        ((current, InputRole.TEMPLATE),
                         (extraction, InputRole.CONDITION)), (purified,))
            current = purified
            reviewed.extend((prefix + "_dna_extraction", prefix + "_extract"))

        sonication = graph.input(prefix + "_sonication", label + " sonication condition",
                                 MaterialKind.CONDITION)
        fragmented = graph.product(prefix + "_fragments", label + " DNA fragments",
                                   MaterialKind.MOLECULAR, state=current.state)
        graph.react(prefix + "_fragment", "Fragment " + label + " DNA",
                    StepKind.FRAGMENTATION,
                    ((current, InputRole.TEMPLATE),
                     (sonication, InputRole.CONDITION)), (fragmented,))
        repair = graph.input(prefix + "_end_repair", label + " end-repair enzymes",
                             MaterialKind.CATALYST)
        repaired = graph.product(prefix + "_repaired", label + " end-repaired DNA",
                                 MaterialKind.MOLECULAR, state=fragmented.state)
        graph.react(prefix + "_repair_step", "End-repair " + label + " DNA",
                    StepKind.END_REPAIR,
                    ((fragmented, InputRole.TEMPLATE),
                     (repair, InputRole.CATALYST)), (repaired,))
        tailing = graph.input(prefix + "_a_tailing", label + " dA-tailing enzyme",
                              MaterialKind.CATALYST)
        tailed = graph.product(prefix + "_tailed", label + " dA-tailed DNA",
                               MaterialKind.MOLECULAR,
                               state=da_tailed(label + " dA-tailed DNA", label + " fragment"))
        graph.react(prefix + "_tail_step", "Add one 3-prime dA to each end",
                    StepKind.TAILING,
                    ((repaired, InputRole.TEMPLATE),
                     (tailing, InputRole.CATALYST)), (tailed,))
        adapters = graph.input(prefix + "_adapters", label + " Illumina adapters",
                               MaterialKind.ADAPTER, inferred=True)
        ligase = graph.input(prefix + "_adapter_ligase", label + " DNA ligase",
                             MaterialKind.CATALYST)
        adapted = graph.product(prefix + "_adapted", label + " adapter-ligated DNA",
                                MaterialKind.MOLECULAR, state=final, inferred=True)
        graph.react(prefix + "_ligation", "Ligate Illumina adapters to " + label,
                    StepKind.LIGATION,
                    ((tailed, InputRole.TEMPLATE),
                     (adapters, InputRole.STRUCTURAL),
                     (ligase, InputRole.CATALYST)), (adapted,),
                    geometry=adapter_geometry)
        reviewed.extend((
            prefix + "_sonication", prefix + "_fragment", prefix + "_end_repair",
            prefix + "_repair_step", prefix + "_a_tailing", prefix + "_tail_step",
            prefix + "_adapters", prefix + "_adapter_ligase", prefix + "_ligation",
        ))

        if not wga:
            denaturing = graph.input(prefix + "_denaturing", label + " DNA denaturation condition",
                                     MaterialKind.CONDITION)
            labelled_strand = graph.product(prefix + "_brdu_strand", label + " BrdU-labelled strand",
                                            MaterialKind.MOLECULAR, state=adapted.state)
            partner = graph.product(prefix + "_unlabelled_partner", label + " unlabelled partner strand",
                                    MaterialKind.WASTE)
            graph.react(prefix + "_denature", "Denature " + label + " adapter-ligated DNA",
                        StepKind.STRAND_SEPARATION,
                        ((adapted, InputRole.TEMPLATE),
                         (denaturing, InputRole.CONDITION)), (labelled_strand, partner))
            graph.discard(partner, "the BrdU-containing nascent strand is the selected analyte")
            antibody = graph.input(prefix + "_brdu_antibody", label + " anti-BrdU capture reagent",
                                   MaterialKind.SELECTION)
            selected = graph.product(prefix + "_selected", label + " BrdU-enriched library template",
                                     MaterialKind.MOLECULAR, state=adapted.state)
            unselected = graph.product(prefix + "_unselected", label + " unselected DNA",
                                       MaterialKind.WASTE)
            graph.react(prefix + "_brdu_capture", "Immunoprecipitate BrdU-labelled " + label + " DNA",
                        StepKind.AFFINITY_CAPTURE,
                        ((labelled_strand, InputRole.TEMPLATE),
                         (antibody, InputRole.SELECTION)), (selected, unselected))
            graph.discard(unselected, "not retained by anti-BrdU immunoprecipitation")
            current = selected
            reviewed.extend((
                prefix + "_denaturing", prefix + "_denature", prefix + "_brdu_antibody",
                prefix + "_brdu_capture",
            ))
        else:
            current = adapted

        forward = graph.input(prefix + "_index_forward", label + " indexed PCR forward primer",
                              MaterialKind.OLIGO)
        reverse = graph.input(prefix + "_index_reverse", label + " indexed PCR reverse primer",
                              MaterialKind.OLIGO)
        pcr_polymerase = graph.input(prefix + "_pcr_polymerase", label + " PCR polymerase",
                                     MaterialKind.CATALYST)
        library = graph.product(prefix + "_library", label + " sequencing library",
                                MaterialKind.LIBRARY, state=final)
        graph.react(prefix + "_index_pcr", "PCR-index the " + label + " library",
                    StepKind.PCR,
                    ((current, InputRole.TEMPLATE),
                     (forward, InputRole.STRUCTURAL),
                     (reverse, InputRole.STRUCTURAL),
                     (pcr_polymerase, InputRole.CATALYST)), (library,),
                    geometry=PCRGeometry.unavailable(
                        final, reason="Sample index choice varies while the final TruSeq endpoint is canonical."))
        graph.sequence(library, plan)
        reviewed.extend((
            prefix + "_index_forward", prefix + "_index_reverse",
            prefix + "_pcr_polymerase", prefix + "_index_pcr", prefix + "_library",
        ))

    if key == "el-repli-seq":
        cells = graph.input("cycling_cells", "asynchronously cycling cells",
                            MaterialKind.BIOLOGICAL,
                            state=molecule("asynchronously cycling cells",
                                           ("replicating chromosomes",)), primary=True)
        brdu = graph.input("brdu", "BrdU pulse-label nucleotide",
                           MaterialKind.MOLECULAR)
        replication_polymerase = graph.input("replication_polymerase",
                                             "active cellular DNA polymerases",
                                             MaterialKind.CATALYST)
        labelled = graph.product("brdu_labelled_cells", "BrdU pulse-labelled S-phase cells",
                                 MaterialKind.MOLECULAR,
                                 state=NascentPulse("BrdU", ("early S phase", "late S phase")).labelled_scene())
        graph.react("pulse_label", "Pulse-label nascent DNA with BrdU",
                    StepKind.EXTENSION,
                    ((cells, InputRole.TEMPLATE), (brdu, InputRole.STRUCTURAL),
                     (replication_polymerase, InputRole.CATALYST)), (labelled,))
        gates = graph.input("facs_gates", "early/late S-phase FACS gates",
                            MaterialKind.SELECTION)
        early = graph.product("early_cells", "early-S BrdU-labelled cells",
                              MaterialKind.MOLECULAR, state=labelled.state)
        late = graph.product("late_cells", "late-S BrdU-labelled cells",
                             MaterialKind.MOLECULAR, state=labelled.state)
        graph.react("sort_fractions", "Sort early- and late-S fractions",
                    StepKind.SELECTION,
                    ((labelled, InputRole.TEMPLATE),
                     (gates, InputRole.SELECTION)), (early, late))
        reviewed.extend(("cycling_cells", "brdu", "replication_polymerase",
                         "pulse_label", "facs_gates", "sort_fractions"))
        finish_short_read_branch("early", early, "early-S", wga=False)
        finish_short_read_branch("late", late, "late-S", wga=False)
    else:
        cells = graph.input("sorted_cells", "S-phase cell plus parallel G1 controls",
                            MaterialKind.BIOLOGICAL,
                            state=molecule("flow-sorted cells", ("S-phase cell", "G1 controls")),
                            primary=True)
        s_cell = graph.product("s_phase_cell", "one S-phase cell",
                               MaterialKind.MOLECULAR, state=cells.state)
        g1_cell = graph.product("g1_control", "parallel G1 control cell",
                                MaterialKind.MOLECULAR, state=cells.state)
        graph.react("separate_cells", "Keep S-phase and G1 cells as separate libraries",
                    StepKind.SPLIT, ((cells, InputRole.TEMPLATE),), (s_cell, g1_cell))
        reviewed.extend(("sorted_cells", "separate_cells"))
        finish_short_read_branch("s_phase", s_cell, "S-phase cell", wga=True)
        finish_short_read_branch("g1", g1_cell, "G1 control", wga=True)
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def get(key: str):
    d = DEFS[key]
    oligos = ()
    if key == "4c-seq":
        oligos = (
            "<b>Reading inverse-PCR primer</b> — 5′-TACACGACGCTCTTCCGATCT-[0–5-nt spacer]-[viewpoint-specific sequence]-3′",
            "<b>Non-reading inverse-PCR primer</b> — 5′-ACTGGAGTTCAGACGTGTGCTCTTCCGATCT-[viewpoint-specific sequence]-3′",
            "<b>Second-round universal forward</b> — 5′-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT-3′",
        )
    if d.standard_truseq:
        # These defining protocols predate routine unique-dual indexing and describe one
        # PCR sample index.  Do not silently modernize their read structure.
        final, primers = truseq_library(list(d.final_parts), d.title + " final library",
                                        dual_index=False)
    else:
        final = Construct([role("Illumina-compatible left end"), *d.final_parts,
                           role("Illumina-compatible right end")], name=d.title + " final library")
        primers = ()
    spec = SimpleNamespace(
        TITLE=d.title, NOTES="01_" + key.replace("-", "_") + ".html",
        SOURCE=f'Defining source: <a href="https://doi.org/{d.doi}">{d.citation}</a>.',
        SUMMARY=d.summary, CAVEAT=d.caveat, FINAL_LIBRARY=final, SEQ_PRIMERS=primers,
        FINAL_CAPTION=d.sequencing, SEQUENCING_INTRO=d.sequencing, READ_LENGTHS=None,
        SEQUENCING_ENDING=d.sequencing,
        SEQUENCING_UNAVAILABLE=("The defining protocol uses locus-specific or custom primer sequence that cannot be represented by one universal binding sequence."
                                if not primers else ""),
        workflow=lambda: workflow_for(key), oligos=lambda: oligos,
    )
    if key in {"4c-seq", "ng-capture-c", "umi-4c", "tri-c"}:
        spec.protocol_graph = lambda: _conformation_graph(key, d, final, primers)
    elif key in {"el-repli-seq", "screpli-seq"}:
        spec.protocol_graph = lambda: _replication_graph(key, d, final, primers)
    return spec
