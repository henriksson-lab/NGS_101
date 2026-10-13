"""Shared molecular models for rare-cell RNA capture, BAGs and CXXC CAP-seq.

The methods collected here share a useful design constraint: the named protocol is
primarily a *selection or partitioning operation*.  Library shells are therefore built
from canonical Illumina constructors where the source names a standard chemistry, and
otherwise remain explicit role tokens rather than guessed bases.
"""
from __future__ import annotations

from dataclasses import dataclass
from types import SimpleNamespace

import illumina as il
import rt
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry
from batch_ngs import nextera_library, seg
from chemdraw import Construct, MolecularState, Scene, Segment, Workflow, feature
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from protocol_model import role, role_library


def _scene(label: str, *parts: Segment) -> Scene:
    sc = Scene()
    sc.strand("molecule", list(parts), label=label)
    sc.labels("molecule")
    return sc


def _workflow(initial: Scene, steps: tuple[tuple[str, Scene, str], ...]) -> Workflow:
    w = Workflow(MolecularState("Starting material", tuple(initial.rows())))
    for action, state, note in steps:
        w.react(action, state.rows(), name=action, note=note)
    return w


TRANSCRIPT = role("transcript-derived insert")
GENOMIC = role("genomic insert")
RAISIN_CST = (("CHAPS", "0.49% w/v"), ("NaCl", "146 mM"),
              ("CaCl2", "1 mM"), ("MgCl2", "21 mM"),
              ("Tris pH 8.0", "10 mM"))
CAP_HUMAN_WASH_MM = 600
CAP_MOUSE_WASH_MM = 560


def raisin_particle() -> Scene:
    sc = Scene()
    sc.strand("nucleus", [role("intact nucleus")], label="RAISIN")
    sc.strand("ribosome", [role("rough ER / attached ribosomes"),
                            role("ribosome-bound poly(A) RNA")], label="retained exterior RNA")
    sc.note("ribosome", "CHAPS/salt extraction retains nucleus-associated ribosomes and mRNA")
    return sc


def smart_cdna_scene() -> Scene:
    sc = Scene()
    sc.strand("cDNA", [role("oligo-dT primer / SMART handle"), TRANSCRIPT,
                        role("copied locked TSO")], label="full-length cDNA")
    sc.note("cDNA", "poly(dT) reverse transcription + template switching; 21-cycle PCR")
    return sc


def raisin_workflow() -> Workflow:
    final, _ = raisin_library()
    return _workflow(_scene("frozen tissue", role("frozen tissue")), (
        ("Extract with CST (RAISIN) buffer", raisin_particle(),
         "CHAPS with salt, Ca²⁺ and Mg²⁺ preserves intact nuclei with attached ribosomes."),
        ("Sort one RAISIN per well, lyse and purify RNA", _scene(
            "single-RAISIN RNA", role("nuclear RNA"), role("ribosome-bound poly(A) RNA")),
         "The plate implementation sorts into TCL/β-mercaptoethanol and uses RNAClean XP."),
        ("Oligo-dT reverse-transcribe with a locked TSO, then amplify", smart_cdna_scene(),
         "Maxima reverse transcriptase is followed by 21 cycles of KAPA HiFi PCR."),
        ("Tagment amplified cDNA with Nextera XT and add sample indexes",
         Scene.duplex(list(final), label="RAISIN RNA-seq library"),
         "The custom barcode-adapter sequences were available only on request."),
    ))


def raisin_library() -> tuple[Construct, tuple]:
    return nextera_library([seg("RAISIN transcript-derived insert", "X" * 36,
                                placeholder=True)], "RAISIN RNA-seq library")


CELL_BARCODE = feature("cell_barcode", "cell_barcode", "whitelist",
                       whitelist="10x GemCode v2")
UMI = feature("umi", "umi", "random")
SAMPLE_INDEX = feature("sample_index_i7", "sample_index", "unknown")


def gemcode_v2_library() -> tuple[Construct, tuple]:
    """The standard 10x v2 endpoint used by droplet RAISIN and MIRACL-seq."""
    lib = Construct([
        Segment("P5", il.P5, "p5"),
        Segment("Read 1 arm", il.TRUSEQ_READ1[4:], "r1"),
        Segment("cell barcode", "B" * 16, "cbc", placeholder=True,
                feature=CELL_BARCODE),
        Segment("UMI", "U" * 10, "umi", placeholder=True, feature=UMI),
        Segment("poly(dT)30", "T" * 30),
        Segment("V", "V", placeholder=True), Segment("N", "N", placeholder=True),
        Segment("RAISIN transcript-derived insert", "X" * 36, placeholder=True),
        Segment("dA junction", "A"),
        Segment("Index 1 / Read 2 arm", il.INDEX1_PRIMER, "r2"),
        Segment("i7 reverse complement", "I" * 8, "cbc", placeholder=True,
                feature=SAMPLE_INDEX),
        Segment("P7 reverse complement", il.P7_RC, "p7"),
    ], name="GemCode v2 RAISIN library")
    primers = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["R2"])
    problems = sp.verify(lib, primers, required_roles=tuple(p.role for p in primers))
    if problems:
        raise ValueError("invalid GemCode v2 endpoint: " + "; ".join(problems))
    return lib, primers


@dataclass(frozen=True)
class MiraclLoading:
    coefficient: int

    def __post_init__(self) -> None:
        if self.coefficient < 2:
            raise ValueError("MIRACL-seq is defined by deliberate droplet overloading")

    def partition(self) -> Scene:
        return _scene(f"{self.coefficient}× overloaded GEM",
                      role("one 10x gel bead"), role("multiple RAISIN nuclei"),
                      role("oil-emulsion partition"))


def miracl_workflow(coefficient: int = 2) -> Workflow:
    loading = MiraclLoading(coefficient)
    final, _ = gemcode_v2_library()
    barcoded = _scene("barcoded first strand", role("10x Read 1 handle"),
                      role("16-nt cell barcode", molecular_feature=CELL_BARCODE),
                      role("10-nt UMI", molecular_feature=UMI), TRANSCRIPT)
    return _workflow(raisin_particle(), (
        (f"Overload the GemCode v2 channel ({coefficient}× shown)", loading.partition(),
         "The defining experiments tested 2×–8× loading; rare-cell homotypic multiplets remain unlikely."),
        ("Lyse in GEMs and perform barcoded reverse transcription", barcoded,
         "One bead barcode labels cDNA from every nucleus sharing that partition."),
        ("Amplify, shear, attach the 5′ adapter and sample index",
         Scene.duplex(list(final), label="MIRACL-seq library"),
         "This is the standard GemCode v2 molecular endpoint; heterotypic multiplets are removed computationally."),
    ))


# Verbatim capture oligos printed for the BAG platform.  The hybrid paper uses both.
BAG_DNA_CAPTURE_WRITTEN = "/5ACryd//iSp18/TGTGTTGGGTGTGTTTGGKKKKKKKGKKKKKKKKNN"
BAG_RNA_CAPTURE_WRITTEN = "/5ACryd//iSp18/AAGCAGTGGTATCAACGCAGAGTNNWNNNSTTTTTTTTTTTTTTTTTTTTTTTTTTTT"

BAG_BARCODE = feature("bag_barcode", "cell_barcode", "combinatorial",
                      group="bag_barcode")
BAG_VARIETAL = feature("varietal_tag", "umi", "combinatorial",
                       group="varietal_tag")


def bag_final_library() -> Construct:
    return role_library("hybrid BAG-seq library", (
        role("Read 1 common sequence"),
        role("BAG barcode round 3", length_bp=6, molecular_feature=BAG_BARCODE),
        role("BAG barcode round 2", length_bp=6, molecular_feature=BAG_BARCODE),
        role("varietal tag round 2", length_bp=4, molecular_feature=BAG_VARIETAL),
        role("BAG barcode round 1", length_bp=5, molecular_feature=BAG_BARCODE),
        role("varietal tag round 1", length_bp=4, molecular_feature=BAG_VARIETAL),
        role("NlaIII-proximal genomic tag", length_bp=4, molecular_feature=BAG_VARIETAL),
        role("CATG NlaIII site", length_bp=4),
        role("DNA- or RNA-derived insert"),
    ))


def bag_workflow() -> Workflow:
    captured = _scene("hybrid BAG", role("polyacrylamide gel matrix"),
                      role("5′-Acrydite random-TG DNA primer"),
                      role("5′-Acrydite poly(dT) RNA primer"),
                      role("single nucleus"))
    copies = _scene("gel-bound copies", role("Acrydite anchor"), GENOMIC,
                    role("Acrydite anchor"), TRANSCRIPT, role("TSO-derived NlaIII site"))
    cut = _scene("NlaIII-cleaved copies", role("Acrydite anchor"),
                 role("DNA- or RNA-derived copy"), role("CATG sticky end"))
    tagged = _scene("three-round tagged molecules", role("UP1"),
                    role("BAG barcode 1 + varietal tag 1"),
                    role("BAG barcode 2 + varietal tag 2"),
                    role("BAG barcode 3"), role("DNA- or RNA-derived copy"))
    final = bag_final_library()
    return _workflow(_scene("single nucleus", role("single sorted nucleus")), (
        ("Co-encapsulate at 85 °C and polymerize the hybrid BAG", captured,
         "Both 5′-Acrydite primers copolymerize into the gel; the lower denaturation temperature preserves RNA."),
        ("Reverse-transcribe RNA and Klenow-extend genomic DNA", copies,
         "Both copy types remain covalently gel-bound; the TSO supplies a double-stranded NlaIII site to cDNA."),
        ("Digest the bound copies with NlaIII", cut,
         "Cleavage exposes the sticky end used by the first split-pool extension."),
        ("Add two ligation/hybridization-extension tags, then a split PCR tag", tagged,
         "The three rounds form a 17-nt BAG barcode and a 12-nt varietal tag."),
        ("Pool and tagment the amplicons", _scene("hybrid BAG-seq library", *tuple(final)),
         "Paired-end 150-bp sequencing starts one read at the tag-bearing end and the other at a random tagmentation end."),
    ))


def cap_library() -> Construct:
    return role_library("CAP-seq library", (role("unmethylated CpG-rich insert"),))


def cap_workflow() -> Workflow:
    adapted = _scene("adapter-ligated genomic library", role("Solexa left adapter"),
                     role("genomic insert with CpG sites"), role("Solexa right adapter"))
    bound = _scene("CXXC-bound library", role("nickel-Sepharose-bound recombinant CXXC"),
                   role("adapter-ligated unmethylated CpG-rich fragment"))
    final = cap_library()
    return _workflow(adapted, (
        ("Bind adapter-ligated input to recombinant CXXC on nickel Sepharose", bound,
         "CXXC affinity selects clustered unmethylated CpGs; the adapters were installed before purification."),
        (f"Wash with {CAP_HUMAN_WASH_MM} mM NaCl (human) or "
         f"{CAP_MOUSE_WASH_MM} mM NaCl (mouse)", bound,
         "The lower-stringency mouse wash retains its relatively CpG-poor islands."),
        ("Elute at 1 M NaCl, pool and sequence",
         _scene("CAP-seq library", *tuple(final)),
         "The defining experiment used single-end Solexa sequencing."),
    ))


def _finish_pcr(graph, prefix, template, final, primers, endpoint, reviewed, *, label):
    forward = graph.input(prefix + "_forward", label + " PCR forward primer",
                          MaterialKind.OLIGO)
    reverse = graph.input(prefix + "_reverse", label + " PCR reverse primer",
                          MaterialKind.OLIGO)
    polymerase = graph.input(prefix + "_polymerase", label + " PCR polymerase",
                             MaterialKind.CATALYST)
    library = graph.product(prefix + "_library", label + " sequencing library",
                            MaterialKind.LIBRARY, state=final)
    graph.react(prefix + "_pcr", "Amplify and index the " + label + " library",
                StepKind.PCR,
                ((template, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (library,), geometry=PCRGeometry.unavailable(
                    final, reason="The source establishes the endpoint but withholds or varies the complete amplification-primer pair."))
    plan = (SequencingPlan(library=final, primers=tuple(primers)) if primers
            else SequencingPlan(limitation=endpoint))
    graph.sequence(library, plan)
    reviewed.extend((prefix + "_forward", prefix + "_reverse", prefix + "_polymerase",
                     prefix + "_pcr", prefix + "_library"))
    return library


def _raisin_graph(spec):
    source = spec.SOURCE
    graph = ProtocolGraph(spec.TITLE, authority=source)
    tissue = graph.input("tissue", "frozen tissue", MaterialKind.BIOLOGICAL,
                         state=_scene("frozen tissue", role("frozen tissue")), primary=True)
    cst = graph.input("cst_buffer", "CST RAISIN extraction buffer",
                      MaterialKind.CONDITION)
    particles = graph.product("raisin_particles", "intact nuclei retaining attached ribosomes",
                              MaterialKind.MOLECULAR, state=raisin_particle())
    graph.react("extract", "Extract RAISIN particles with CST buffer",
                StepKind.EXTRACTION,
                ((tissue, InputRole.TEMPLATE), (cst, InputRole.CONDITION)), (particles,))
    sorter = graph.input("sort_rule", "single-RAISIN fluorescence-sorting rule",
                         MaterialKind.SELECTION)
    selected = graph.product("single_raisin", "one sorted RAISIN per well",
                             MaterialKind.MOLECULAR, state=particles.state)
    excluded = graph.product("excluded_particles", "unsorted or multiplet particles",
                             MaterialKind.WASTE)
    graph.react("sort", "Sort one RAISIN into each well", StepKind.SAMPLING,
                ((particles, InputRole.TEMPLATE), (sorter, InputRole.SELECTION)),
                (selected, excluded), changes_molecule=False)
    graph.discard(excluded, "does not meet the one-particle-per-well selection")
    rt_primer = graph.input("rt_primer", "oligo(dT) SMART reverse-transcription primer",
                            MaterialKind.OLIGO, inferred=True,
                            source=source + " Exact custom primer bases were available only on request.")
    rtase = graph.input("reverse_transcriptase", "Maxima reverse transcriptase",
                        MaterialKind.CATALYST)
    first = graph.product("first_strand", "RAISIN transcript first-strand cDNA",
                          MaterialKind.MOLECULAR, state=smart_cdna_scene(), inferred=True,
                          source=source + " Locked TSO sequence is not printed.")
    graph.react("reverse_transcribe", "Oligo-dT reverse-transcribe RAISIN RNA",
                StepKind.REVERSE_TRANSCRIPTION,
                ((selected, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
                 (rtase, InputRole.CATALYST)), (first,))
    tso = graph.input("locked_tso", "locked SMART template-switch oligo",
                      MaterialKind.OLIGO, inferred=True,
                      source=source + " Exact locked TSO sequence is not printed.")
    switched = graph.product("switched", "template-switched full-length RAISIN cDNA",
                             MaterialKind.MOLECULAR, state=smart_cdna_scene(), inferred=True,
                             source=source + " Exact locked TSO sequence is not printed.")
    graph.react("template_switch", "Template-switch at the transcript 5-prime end",
                StepKind.TEMPLATE_SWITCHING,
                ((first, InputRole.TEMPLATE), (tso, InputRole.STRUCTURAL),
                 (rtase, InputRole.CATALYST)), (switched,))
    pre_f = graph.input("preamp_forward", "SMART cDNA preamplification forward primer",
                        MaterialKind.OLIGO, inferred=True,
                        source=source + " Custom primer bases were available only on request.")
    pre_r = graph.input("preamp_reverse", "SMART cDNA preamplification reverse primer",
                        MaterialKind.OLIGO, inferred=True,
                        source=source + " Custom primer bases were available only on request.")
    pre_poly = graph.input("preamp_polymerase", "KAPA HiFi polymerase",
                           MaterialKind.CATALYST)
    amplified = graph.product("amplified", "21-cycle amplified full-length cDNA",
                              MaterialKind.MOLECULAR, state=switched.state)
    graph.react("preamp", "Amplify full-length cDNA for 21 cycles", StepKind.PCR,
                ((switched, InputRole.TEMPLATE), (pre_f, InputRole.STRUCTURAL),
                 (pre_r, InputRole.STRUCTURAL), (pre_poly, InputRole.CATALYST)),
                (amplified,), geometry=PCRGeometry.unavailable(
                    spec.FINAL_LIBRARY, reason="Custom barcode-adapter primer sequences were available only on request."))
    transposome = graph.input("transposome", "Nextera XT loaded transposome",
                              MaterialKind.ADAPTER)
    tn5 = graph.input("tn5", "Nextera XT transposase", MaterialKind.CATALYST)
    tagged = graph.product("tagged", "tagmented amplified RAISIN cDNA",
                           MaterialKind.MOLECULAR, state=amplified.state)
    waste = graph.product("tn5_waste", "same-end or non-library tagmentation products",
                          MaterialKind.WASTE)
    graph.react("tagment", "Tagment amplified RAISIN cDNA", StepKind.TAGMENTATION,
                ((amplified, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged, waste))
    graph.discard(waste, "does not form a heterologous amplifiable Nextera library")
    reviewed = ["tissue", "cst_buffer", "extract", "sort_rule", "sort", "rt_primer",
                "reverse_transcriptase", "reverse_transcribe", "locked_tso", "template_switch",
                "preamp_forward", "preamp_reverse", "preamp_polymerase", "preamp",
                "transposome", "tn5", "tagment"]
    _finish_pcr(graph, "library", tagged, spec.FINAL_LIBRARY, spec.SEQ_PRIMERS,
                spec.SEQUENCING_INTRO, reviewed, label="RAISIN RNA-seq")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _miracl_graph(spec):
    source = spec.SOURCE
    graph = ProtocolGraph(spec.TITLE, authority=source)
    particles = graph.input("raisin_particles", "unsorted RAISIN nuclei",
                            MaterialKind.BIOLOGICAL, state=raisin_particle(), primary=True)
    droplets = graph.input("overload_condition", "deliberately overloaded GemCode-v2 partition condition",
                           MaterialKind.CONDITION)
    occupied = graph.product("occupied_gems", "overloaded GEMs containing one bead and multiple RAISIN nuclei",
                             MaterialKind.MOLECULAR, state=MiraclLoading(2).partition())
    unresolved = graph.product("unresolved_gems", "empty or uninterpretable GEMs",
                               MaterialKind.WASTE)
    graph.react("partition", "Overload RAISIN nuclei into GemCode-v2 partitions",
                StepKind.PARTITION,
                ((particles, InputRole.TEMPLATE), (droplets, InputRole.CONDITION)),
                (occupied, unresolved), changes_molecule=False)
    graph.discard(unresolved, "does not yield a usable gel-bead barcode partition")
    bead_primer = graph.input("bead_primer", "GemCode-v2 barcode/UMI oligo(dT) bead primer",
                              MaterialKind.OLIGO)
    rtase = graph.input("reverse_transcriptase", "droplet reverse transcriptase",
                        MaterialKind.CATALYST)
    cdna = graph.product("barcoded_cdna", "partition-barcoded RAISIN transcript cDNA",
                         MaterialKind.MOLECULAR,
                         state=_scene("barcoded first strand", role("16-nt cell barcode",
                             molecular_feature=CELL_BARCODE), role("10-nt UMI",
                             molecular_feature=UMI), TRANSCRIPT))
    graph.react("reverse_transcribe", "Lyse in GEMs and perform barcoded reverse transcription",
                StepKind.REVERSE_TRANSCRIPTION,
                ((occupied, InputRole.TEMPLATE), (bead_primer, InputRole.STRUCTURAL),
                 (rtase, InputRole.CATALYST)), (cdna,))
    reviewed = ["raisin_particles", "overload_condition", "partition", "bead_primer",
                "reverse_transcriptase", "reverse_transcribe"]
    _finish_pcr(graph, "library", cdna, spec.FINAL_LIBRARY, spec.SEQ_PRIMERS,
                spec.SEQUENCING_INTRO, reviewed, label="MIRACL-seq GemCode-v2")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _bag_graph(spec):
    source = spec.SOURCE
    graph = ProtocolGraph(spec.TITLE, authority=source)
    nucleus = graph.input("nucleus", "single sorted nucleus containing genomic DNA and RNA",
                          MaterialKind.BIOLOGICAL,
                          state=_scene("single nucleus", role("genomic DNA"), role("RNA")),
                          primary=True)
    gel = graph.input("gel_condition", "hybrid-BAG polymerization and co-encapsulation condition",
                      MaterialKind.CONDITION)
    bag = graph.product("bag", "one nucleus co-encapsulated in a hybrid BAG",
                        MaterialKind.MOLECULAR,
                        state=_scene("hybrid BAG", role("polyacrylamide gel matrix"),
                                     role("Acrydite random-TG primer"),
                                     role("Acrydite poly(dT) primer"), role("single nucleus")))
    excluded = graph.product("excluded", "uncaptured nuclei or empty gels", MaterialKind.WASTE)
    graph.react("encapsulate", "Co-encapsulate and polymerize the hybrid BAG",
                StepKind.PARTITION,
                ((nucleus, InputRole.TEMPLATE), (gel, InputRole.CONDITION)),
                (bag, excluded), changes_molecule=False)
    graph.discard(excluded, "does not retain a one-nucleus hybrid gel particle")
    rt_primer = graph.input("rna_capture_primer", "Acrydite poly(dT) RNA capture primer",
                            MaterialKind.OLIGO)
    rtase = graph.input("reverse_transcriptase", "reverse transcriptase",
                        MaterialKind.CATALYST)
    rna_copy = graph.product("rna_copy", "gel-bound transcript cDNA",
                             MaterialKind.MOLECULAR, state=bag.state)
    graph.react("reverse_transcribe", "Reverse-transcribe gel-captured RNA",
                StepKind.REVERSE_TRANSCRIPTION,
                ((bag, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
                 (rtase, InputRole.CATALYST)), (rna_copy,))
    dna_primer = graph.input("dna_capture_primer", "Acrydite random-TG DNA capture primer",
                             MaterialKind.OLIGO)
    klenow = graph.input("klenow", "Klenow DNA polymerase", MaterialKind.CATALYST)
    copies = graph.product("copies", "gel-bound genomic and transcript-derived copies",
                           MaterialKind.MOLECULAR, state=rna_copy.state)
    graph.react("copy_dna", "Klenow-extend genomic DNA in the same BAG",
                StepKind.EXTENSION,
                ((rna_copy, InputRole.TEMPLATE), (dna_primer, InputRole.STRUCTURAL),
                 (klenow, InputRole.CATALYST)), (copies,))
    nla = graph.input("nlaiii", "NlaIII restriction endonuclease", MaterialKind.CATALYST)
    cut = graph.product("cut", "NlaIII-cleaved gel-bound DNA and cDNA copies",
                        MaterialKind.MOLECULAR, state=copies.state)
    graph.react("digest", "Digest gel-bound copies with NlaIII", StepKind.DIGESTION,
                ((copies, InputRole.TEMPLATE), (nla, InputRole.CATALYST)), (cut,))
    current = cut
    reviewed = ["nucleus", "gel_condition", "encapsulate", "rna_capture_primer",
                "reverse_transcriptase", "reverse_transcribe", "dna_capture_primer",
                "klenow", "copy_dna", "nlaiii", "digest"]
    for round_number in (1, 2):
        prefix = f"round{round_number}"
        adapter = graph.input(prefix + "_adapter", f"BAG barcode/varietal-tag adapter round {round_number}",
                              MaterialKind.ADAPTER)
        ligase = graph.input(prefix + "_ligase", f"BAG tagging ligase round {round_number}",
                             MaterialKind.CATALYST)
        tagged = graph.product(prefix + "_product", f"BAG-tagged molecules after round {round_number}",
                               MaterialKind.MOLECULAR, state=current.state)
        graph.react(prefix + "_ligate", f"Add BAG barcode/tag round {round_number}",
                    StepKind.LIGATION,
                    ((current, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
                     (ligase, InputRole.CATALYST)), (tagged,),
                    geometry=AdapterLigationGeometry.unavailable(
                        "The papers establish barcode lengths and order but do not expose one complete base-level tagging construct."))
        reviewed.extend((prefix + "_adapter", prefix + "_ligase", prefix + "_ligate"))
        current = tagged
    round3_f = graph.input("round3_forward", "third-round split-PCR forward primer",
                           MaterialKind.OLIGO)
    round3_r = graph.input("round3_reverse", "third-round split-PCR reverse primer",
                           MaterialKind.OLIGO)
    polymerase = graph.input("round3_polymerase", "third-round PCR polymerase",
                             MaterialKind.CATALYST)
    amplified = graph.product("amplified", "three-round BAG-tagged amplicons",
                              MaterialKind.MOLECULAR, state=current.state)
    graph.react("round3_pcr", "Add the third BAG barcode by split PCR", StepKind.PCR,
                ((current, InputRole.TEMPLATE), (round3_f, InputRole.STRUCTURAL),
                 (round3_r, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (amplified,), geometry=PCRGeometry.unavailable(
                    spec.FINAL_LIBRARY, reason="The source defines barcode geometry but not one complete universal amplicon."))
    transposome = graph.input("transposome", "library transposome", MaterialKind.ADAPTER)
    tn5 = graph.input("tn5", "Tn5 transposase", MaterialKind.CATALYST)
    library = graph.product("library", "hybrid BAG-seq library", MaterialKind.LIBRARY,
                            state=spec.FINAL_LIBRARY)
    waste = graph.product("tn5_waste", "same-end or non-library tagmentation products",
                          MaterialKind.WASTE)
    graph.react("tagment", "Pool and tagment the BAG-tagged amplicons",
                StepKind.TAGMENTATION,
                ((amplified, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (library, waste))
    graph.discard(waste, "does not retain the tag-bearing end in a sequenceable library")
    graph.sequence(library, SequencingPlan(limitation=spec.SEQUENCING_ENDING))
    reviewed.extend(("round3_forward", "round3_reverse", "round3_polymerase", "round3_pcr",
                     "transposome", "tn5", "tagment", "library"))
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _cap_graph(spec):
    source = spec.SOURCE
    graph = ProtocolGraph(spec.TITLE, authority=source)
    adapted = graph.input("adapted_dna", "adapter-ligated genomic DNA library",
                          MaterialKind.LIBRARY,
                          state=_scene("adapter-ligated genomic library",
                                       role("Solexa left adapter"),
                                       role("genomic insert with CpG sites"),
                                       role("Solexa right adapter")), primary=True)
    cxxc = graph.input("cxxc", "immobilized recombinant CXXC on nickel Sepharose",
                       MaterialKind.SELECTION)
    bound = graph.product("bound", "CXXC-bound unmethylated CpG-rich fragments",
                          MaterialKind.MOLECULAR, state=adapted.state)
    depleted = graph.product("depleted", "CXXC-unbound genomic library", MaterialKind.WASTE)
    graph.react("bind", "Bind adapter-ligated input to immobilized CXXC",
                StepKind.AFFINITY_CAPTURE,
                ((adapted, InputRole.TEMPLATE), (cxxc, InputRole.SELECTION)),
                (bound, depleted))
    graph.discard(depleted, "lacks sufficient clustered unmethylated CpG affinity")
    wash = graph.input("wash_condition", "species-adjusted high-salt wash condition",
                       MaterialKind.CONDITION)
    washed = graph.product("washed", "washed CXXC-bound CpG-rich fragments",
                           MaterialKind.MOLECULAR, state=bound.state)
    wash_waste = graph.product("wash_waste", "fragments removed by high-salt washing",
                               MaterialKind.WASTE)
    graph.react("wash", "Wash the CXXC-bound library", StepKind.CLEANUP,
                ((bound, InputRole.TEMPLATE), (wash, InputRole.CONDITION)),
                (washed, wash_waste))
    graph.discard(wash_waste, "does not remain bound at the chosen species-specific stringency")
    elution = graph.input("elution_condition", "1-M NaCl elution condition",
                          MaterialKind.CONDITION)
    library = graph.product("library", "CXXC-enriched CAP-seq library",
                            MaterialKind.LIBRARY, state=spec.FINAL_LIBRARY)
    graph.react("elute", "Elute and pool CXXC-enriched fragments", StepKind.EXTRACTION,
                ((washed, InputRole.TEMPLATE), (elution, InputRole.CONDITION)), (library,))
    graph.sequence(library, SequencingPlan(limitation=spec.SEQUENCING_ENDING))
    reviewed = ["adapted_dna", "cxxc", "bind", "wash_condition", "wash",
                "elution_condition", "elute", "library"]
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _attach_graph(key, spec):
    factories = {
        "raisin-rna-seq": _raisin_graph,
        "miracl-seq": _miracl_graph,
        "hybrid-bag-seq": _bag_graph,
        "cap-seq": _cap_graph,
    }
    spec.protocol_graph = lambda: factories[key](spec)
    return spec


def _spec(key, **fields):
    return _attach_graph(key, SimpleNamespace(**fields))


def get(key: str):
    if key == "raisin-rna-seq":
        final, primers = raisin_library()
        return _spec(key,
            TITLE="RAISIN RNA-seq — intact nuclei with ribosome-bound RNA",
            NOTES="01_raisin_rna_seq.html",
            SOURCE='Defining source: <a href="https://doi.org/10.1016/j.cell.2020.08.003">Drokhlyansky et al., Cell (2020)</a>.',
            SUMMARY="CST extraction preserves an intact nucleus together with attached ribosomes and their poly(A) RNA; the plate implementation uses SMART-seq2 followed by Nextera XT.",
            CAVEAT="The paper does not print the locked TSO or custom barcode-adapter sequences; those workflow states are labelled by role rather than assigned bases.",
            FINAL_LIBRARY=final, SEQ_PRIMERS=primers, READ_LENGTHS=None,
            FINAL_CAPTION="A Nextera XT library made from plate-amplified RAISIN cDNA.",
            SEQUENCING_INTRO="Standard Nextera primers bind the completed library as computed below.",
            SEQUENCING_ENDING="", SEQUENCING_UNAVAILABLE="",
            workflow=raisin_workflow, oligos=lambda: ())
    if key == "miracl-seq":
        final, primers = gemcode_v2_library()
        return _spec(key,
            TITLE="MIRACL-seq — droplet overloading for rare-cell RNA-seq",
            NOTES="01_miracl_seq.html",
            SOURCE='Defining source: <a href="https://doi.org/10.1016/j.cell.2020.08.003">Drokhlyansky et al., Cell (2020)</a>.',
            SUMMARY="MIRACL-seq deliberately overloads GemCode v2 droplets with unsorted RAISIN nuclei, increasing rare-cell recovery while computationally removing heterotypic multiplets.",
            CAVEAT="MIRACL-seq changes partition loading, not the GemCode v2 oligo chemistry; the final construct therefore uses the standard v2 molecular endpoint.",
            FINAL_LIBRARY=final, SEQ_PRIMERS=primers, READ_LENGTHS={"Read 1": 26},
            FINAL_CAPTION="GemCode v2 3′ gene-expression library from overloaded RAISIN partitions.",
            SEQUENCING_INTRO="Read 1 records the 16-nt bead barcode and 10-nt UMI; Read 2 interrogates the transcript-derived insert.",
            SEQUENCING_ENDING="", SEQUENCING_UNAVAILABLE="",
            workflow=miracl_workflow, oligos=lambda: ())
    if key == "hybrid-bag-seq":
        return _spec(key,
            TITLE="Hybrid BAG-seq — joint single-nucleus DNA and RNA sequencing",
            NOTES="01_hybrid_bag_seq.html",
            SOURCE='Defining source: <a href="https://doi.org/10.1186/s13059-025-03790-5">Li et al., Genome Biology (2025)</a>; capture-primer sequences were checked against the original <a href="https://doi.org/10.1101/gr.253047.119">BAG-seq paper</a>.',
            SUMMARY="Acrydite random-TG and poly(dT) primers copolymerized into one gel capture genomic DNA and RNA-derived cDNA from the same nucleus before three rounds of split-pool tagging.",
            CAVEAT="The papers establish the barcode geometry and print the capture primers, but do not expose one complete base-level final adapter construct; unreleased regions are shown by role.",
            FINAL_LIBRARY=bag_final_library(), SEQ_PRIMERS=(), READ_LENGTHS=None,
            FINAL_CAPTION="Role-level paired-end library; Read 1 begins with the split-pool tags and Read 2 begins at a random tagmentation end.",
            SEQUENCING_INTRO="", SEQUENCING_ENDING="NovaSeq 6000 paired-end 150-bp reads were used. Read 1 begins with an 80-base tag/adapter block containing the 17-nt BAG barcode and 12-nt varietal tag.",
            SEQUENCING_UNAVAILABLE="The custom final-library sequencing-primer sites are not printed in the authoritative sources.",
            workflow=bag_workflow,
            oligos=lambda: (f"5′-{BAG_DNA_CAPTURE_WRITTEN}-3′ — Acrydite random-TG DNA capture primer",
                             f"5′-{BAG_RNA_CAPTURE_WRITTEN}-3′ — Acrydite poly(dT) RNA capture primer"))
    if key == "cap-seq":
        return _spec(key,
            TITLE="CAP-seq — CXXC affinity purification sequencing",
            NOTES="01_cap_seq.html",
            SOURCE='Defining source: <a href="https://doi.org/10.1371/journal.pgen.1001134">Illingworth et al., PLOS Genetics (2010)</a>.',
            SUMMARY="Adapter-ligated genomic DNA is selected on immobilized recombinant CXXC for clusters of unmethylated CpGs, then eluted and sequenced.",
            CAVEAT="The defining paper refers to earlier Solexa library methods without printing the adapter or sequencing-primer bases; the library ends are therefore shown by role.",
            FINAL_LIBRARY=cap_library(), SEQ_PRIMERS=(), READ_LENGTHS=None,
            FINAL_CAPTION="Solexa library enriched for fragments carrying clustered unmethylated CpGs.",
            SEQUENCING_INTRO="", SEQUENCING_ENDING="The defining study used single-end Solexa reads to map CXXC-enriched fragments.",
            SEQUENCING_UNAVAILABLE="Exact sequencing-primer sites cannot be drawn because the cited library protocol does not print the adapter sequences.",
            workflow=cap_workflow, oligos=lambda: ())
    raise KeyError(key)
