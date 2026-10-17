"""Source-owned protocol definition and reaction graph.

Authority: Cho et al., Cell 184:3559 (2021), PMC8238917 JATS (published author
manuscript), STAR Methods lines 211-285 (seed library, 1st-Seq, HDMI-array
processing, tissue, RT, second strand, library PCR) and analysis lines 341-343.
"""
from __future__ import annotations

from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            SourceClaim, SourceDiscrepancy, StepKind)
from protocol_specs import Definition, get, workflow
from rna_library_operations import molecule
import seqscope_chem as sc
from seqscope_chem import oligos as seqscope_oligos

SOURCE_LABEL = "Cho et al., Cell (2021)"
DOI = "10.1016/j.cell.2021.05.010"
AUTHORITY = f"{SOURCE_LABEL}, PMC8238917 STAR Methods; doi:{DOI}"

SEQUENCING = (
    "1st-Seq: the HDMI-DraI seed is clustered on a MiSeq v3 flow cell and read single-end "
    "with the custom Read1-DraI primer (25 cycles for HDMI-DraI, 37 for HDMI32-DraI); the run "
    "stops after Read 1 and the flow cell is retrieved. 2nd-Seq: the PCR library is read "
    "paired-end (100–150 bp; Illumina and BGI platforms) with the TruSeq primers. Read 1 "
    "gives the HDMI (the first 20 bases for HDMI-DraI, 30 for HDMI32-DraI, matched to the "
    "reverse complement of the 1st-Seq read); the first 9 bases of Read 2 are the random-"
    "primer UMI, followed by transcript sequence. The i7 read length is not stated.")

FINAL_LIBRARY, SEQ_PRIMERS = sc.final_library()

DEFINITION = Definition(
    "Seq-Scope — flow-cell spatial transcriptomics", DOI, SOURCE_LABEL,
    "An Illumina flow cell is first decoded as a dense HDMI barcode array and then reused "
    "to capture spatially resolved tissue transcripts.",
    "truseq", ("spatial HDMI", "HR1 / oligo-dT", "transcript cDNA", "random-primer UMI"),
    SEQUENCING, "RNA",
    "The open Methods print all primary and alternative seed oligos, both custom first-run "
    "primers, the random second-strand primer and both PCR rounds with phosphorothioates. "
    "The final library is assembled from those printed oligos for the HDMI-DraI (20-nt HDMI) "
    "seed used for the liver and colon data. Seed bridge amplification and MiSeq read "
    "chemistry are on-instrument and are not drawn base by base. The HDMI-DraI seed is "
    "labelled with \"P5\" at its 5′ end and \"P7\" at its 3′ end in the paper, although "
    "those bases are the Illumina P7 adapter and the reverse complement of P5; segment "
    "names here follow the bases.")


def _source_workflow():
    return workflow(molecule("Illumina flow cell", ("random HDMI oligos",)), (
        ("Cluster and sequence the HDMI array in situ",
         molecule("decoded capture surface", ("mapped HDMI cluster", "poly(dT) capture end")),
         "The first sequencing run maps each spatial barcode to an imaging coordinate."),
    ))


SPEC = get('seq-scope', definition=DEFINITION, workflow_fn=_source_workflow,
           final_override=(FINAL_LIBRARY, SEQ_PRIMERS), oligos_fn=seqscope_oligos)


def _strip(seq):
    return seq.replace("*", "")


def _cleanup(graph, reviewed, ident, label, template, selection_label, kept_label,
             waste_label, waste_reason, state=None, kind=StepKind.CLEANUP,
             material=MaterialKind.MOLECULAR):
    matrix = graph.input(ident + "_matrix", selection_label, MaterialKind.SELECTION)
    kept = graph.product(ident + "_kept", kept_label, material,
                         state=state if state is not None else template.state)
    waste = graph.product(ident + "_waste", waste_label, MaterialKind.WASTE)
    graph.react(ident, label, kind,
                ((template, InputRole.TEMPLATE), (matrix, InputRole.SELECTION)),
                (kept, waste))
    graph.discard(waste, waste_reason)
    reviewed.extend((ident + "_matrix", ident))
    return kept


def _local_graph():
    graph = ProtocolGraph(DEFINITION.title, authority=AUTHORITY)
    r = []
    exact = sc.graph_states()
    seed = graph.input("seed", "HDMI-DraI ssDNA seed library (IDT Ultramer)",
                       MaterialKind.OLIGO, primary=True)
    flow_cell = graph.input("flow_cell", "MiSeq v3 flow cell and on-instrument cluster reagents",
                            MaterialKind.CONDITION)
    clusters = graph.product("clusters", "HDMI clusters on the MiSeq flow cell",
                             MaterialKind.MOLECULAR,
                             state=molecule("HDMI cluster", ("surface-anchored HDMI-DraI copies",)))
    graph.react("cluster", "Load 100 pM HDMI-DraI and generate clusters on the MiSeq",
                StepKind.MAGIC, ((seed, InputRole.TEMPLATE), (flow_cell, InputRole.CONDITION)),
                (clusters,),
                reason="On-instrument bridge amplification and strand linearisation are not described in the paper.")
    read1_drai = graph.input("read1_drai", "Read1-DraI custom Read 1 primer", MaterialKind.OLIGO,
                             state=exact["read1_drai"])
    sbs = graph.input("first_seq_sbs", "MiSeq v3 SBS polymerase and reagents (manual mode, 25 cycles single end, no re-synthesis)",
                      MaterialKind.CATALYST)
    decoded = graph.product("decoded", "coordinate-decoded clusters with Read1-DraI duplex",
                            MaterialKind.MOLECULAR,
                            state=molecule("decoded cluster", (
                                "P7-anchored HDMI-DraI strand", "Read1-DraI extension product")))
    graph.react("first_seq", "1st-Seq: read the HDMI with Read1-DraI and retrieve the flow cell",
                StepKind.EXTENSION,
                ((clusters, InputRole.TEMPLATE), (read1_drai, InputRole.STRUCTURAL),
                 (sbs, InputRole.CATALYST)), (decoded,))
    drai = graph.input("drai", "DraI (NEB R0129), 1× CutSmart, 37 °C overnight", MaterialKind.CATALYST)
    cut = graph.product("drai_cut", "clusters with exposed 3′ oligo-dT", MaterialKind.MOLECULAR,
                        state=molecule("exposed capture strand", (
                            "P7 adapter", "TruSeq Read 1", "HDMI", "HR1", "oligo-dT (30 T)")))
    released = graph.product("drai_released", "released DraI-adapter / P5-complement ends",
                             MaterialKind.WASTE)
    graph.react("drai_digest", "DraI cuts TTT|AAA to remove the P5 end and expose oligo-dT",
                StepKind.DIGESTION,
                ((decoded, InputRole.TEMPLATE), (drai, InputRole.CATALYST)), (cut, released))
    graph.discard(released, "washed off the flow cell")
    exo1 = graph.input("exo1_array", "Exonuclease I (NEB M2903), 37 °C 45 min", MaterialKind.CATALYST)
    trimmed = graph.product("exo1_array_kept", "Read1-DraI-protected HDMI strands",
                            MaterialKind.MOLECULAR, state=cut.state)
    nonspecific = graph.product("exo1_array_waste", "non-specific single-stranded DNA",
                                MaterialKind.WASTE)
    graph.react("exo1_array_digest", "Exo I removes non-specific ssDNA", StepKind.DIGESTION,
                ((cut, InputRole.TEMPLATE), (exo1, InputRole.CATALYST)), (trimmed, nonspecific))
    graph.discard(nonspecific, "degraded by Exo I and washed away")
    naoh = graph.input("naoh_array", "water, 3× 0.1 N NaOH 5 min, 0.1 M Tris pH 7.5, water washes",
                       MaterialKind.CONDITION)
    array = graph.product("hdmi_array", "single-stranded HDMI-array (oligo-dT 3′ ends)",
                          MaterialKind.MOLECULAR, state=cut.state)
    primer_waste = graph.product("read1_drai_removed", "denatured Read1-DraI strand",
                                 MaterialKind.WASTE)
    graph.react("denature_array", "Denature and wash away the Read1-DraI strand",
                StepKind.STRAND_SEPARATION,
                ((trimmed, InputRole.TEMPLATE), (naoh, InputRole.CONDITION)), (array, primer_waste))
    graph.discard(primer_waste, "removed by NaOH washes")
    tissue = graph.input("tissue", "10 µm fresh-frozen section on the disassembled HDMI-array, "
                         "4% formaldehyde, H&E imaging", MaterialKind.BIOLOGICAL)
    release = graph.input("rna_release", "collagenase I 37 °C 20 min, then 1 mg/mL pepsin in 0.1 M HCl 37 °C 10 min",
                          MaterialKind.CONDITION)
    rna = graph.product("released_rna", "poly(A) RNA released from the fixed section",
                        MaterialKind.MOLECULAR, state=molecule("tissue RNA", ("poly(A) RNA",)))
    graph.react("release_rna", "Release RNA from the fixed tissue over the array",
                StepKind.PERMEABILIZATION,
                ((tissue, InputRole.TEMPLATE), (release, InputRole.CONDITION)), (rna,))
    captured = graph.product("captured", "tissue poly(A) RNA hybridized to oligo-dT",
                             MaterialKind.MOLECULAR,
                             state=molecule("capture duplex", ("HDMI capture strand", "poly(A) RNA")))
    graph.react("capture", "Capture released RNA on the HDMI-array oligo-dT",
                StepKind.HYBRIDIZATION,
                ((array, InputRole.TEMPLATE), (rna, InputRole.STRUCTURAL)), (captured,))
    rt = graph.input("maxima", "Maxima H− RT, Ficoll PM-400, actinomycin D, 42 °C overnight",
                     MaterialKind.CATALYST)
    hybrid = graph.product("rt_hybrid", "surface cDNA:RNA hybrid", MaterialKind.MOLECULAR,
                           state=molecule("surface cDNA", (
                               "P7 adapter", "TruSeq Read 1", "HDMI", "HR1", "oligo-dT",
                               "first-strand cDNA")))
    graph.react("reverse_transcribe", "Extend the array oligo-dT on captured RNA",
                StepKind.REVERSE_TRANSCRIPTION,
                ((captured, InputRole.TEMPLATE), (captured, InputRole.SELF_PRIMER),
                 (rt, InputRole.CATALYST)), (hybrid,))
    exo1b = graph.input("exo1_tissue", "Exonuclease I, 37 °C 45 min", MaterialKind.CATALYST)
    extended = graph.product("exo1_tissue_kept", "RNA-hybridized cDNA strands",
                             MaterialKind.MOLECULAR, state=hybrid.state)
    unextended = graph.product("exo1_tissue_waste", "capture oligos not hybridized to mRNA",
                               MaterialKind.WASTE)
    graph.react("exo1_tissue_digest", "Exo I removes DNA not hybridized with mRNA",
                StepKind.DIGESTION,
                ((hybrid, InputRole.TEMPLATE), (exo1b, InputRole.CATALYST)), (extended, unextended))
    graph.discard(unextended, "degraded by Exo I")
    pk = graph.input("proteinase_k", "tissue digestion buffer with proteinase K, 37 °C 40 min",
                     MaterialKind.CATALYST)
    cleared = graph.product("tissue_cleared", "array cDNA:RNA after tissue removal",
                            MaterialKind.MOLECULAR, state=hybrid.state)
    debris = graph.product("tissue_debris", "digested tissue", MaterialKind.WASTE)
    graph.react("digest_tissue", "Digest the tissue from the array", StepKind.DEPROTEINIZATION,
                ((extended, InputRole.TEMPLATE), (pk, InputRole.CATALYST)), (cleared, debris))
    graph.discard(debris, "washed off the array")
    naoh2 = graph.input("naoh_rna", "water, 3× 0.1 N NaOH 5 min, Tris, water washes",
                        MaterialKind.CONDITION)
    first = graph.product("first_strand", "single-stranded surface cDNA", MaterialKind.MOLECULAR,
                          state=hybrid.state)
    rna_waste = graph.product("rna_removed", "template mRNA", MaterialKind.WASTE)
    graph.react("strip_rna", "Remove all mRNA with NaOH", StepKind.STRAND_SEPARATION,
                ((cleared, InputRole.TEMPLATE), (naoh2, InputRole.CONDITION)), (first, rna_waste))
    graph.discard(rna_waste, "removed by NaOH washes")
    random_primer = graph.input("random_primer", "TruSeq Read2-conjugated N9 random primer (IDT)",
                                MaterialKind.OLIGO)
    klenow = graph.input("klenow", "Klenow fragment (3′→5′ exo−), NEBuffer 2, 37 °C 2 h",
                         MaterialKind.CATALYST)
    duplex = graph.product("second_strand", "surface cDNA paired with random-primed second strand",
                           MaterialKind.MOLECULAR, state=sc.second_strand_duplex())
    graph.react("second_strand_synthesis", "Random-primed second-strand synthesis on the array",
                StepKind.EXTENSION,
                ((first, InputRole.TEMPLATE), (random_primer, InputRole.STRUCTURAL),
                 (klenow, InputRole.CATALYST)), (duplex,))
    washed = _cleanup(graph, r, "wash_unbound", "Wash away DNA not bound to the array",
                      duplex, "three water washes", "array-bound second-strand duplexes",
                      "unbound primers and unanchored products", "washed off before elution")
    naoh3 = graph.input("naoh_elute", "2× 30 µL 0.1 N NaOH 5 min, neutralized with 3 M KOAc",
                        MaterialKind.CONDITION)
    eluted = graph.product("eluted", "eluted second strand", MaterialKind.MOLECULAR,
                           state=duplex.state)
    surface = graph.product("surface_left", "array-anchored first-strand cDNA", MaterialKind.WASTE)
    graph.react("elute_second_strand", "Elute the second strand with NaOH",
                StepKind.STRAND_SEPARATION,
                ((washed, InputRole.TEMPLATE), (naoh3, InputRole.CONDITION)), (eluted, surface))
    graph.discard(surface, "remains attached to the HDMI-array, which is not used further")
    purified = _cleanup(graph, r, "ampure_18", "AMPure XP 1.8× cleanup", eluted,
                        "AMPure XP 1.8×, 40 µL elution", "purified second strand",
                        "AMPure supernatant", "not bound by the beads")
    kapa = graph.input("kapa_pcr1", "KAPA HiFi HotStart ReadyMix, 13–15 cycles", MaterialKind.CATALYST)
    f1 = graph.input("pcr1_forward", "first-round forward primer (3′ phosphorothioates)",
                     MaterialKind.OLIGO, state=exact["pcr1_forward"])
    r1 = graph.input("pcr1_reverse", "first-round reverse primer (3′ phosphorothioates)",
                     MaterialKind.OLIGO, state=exact["pcr1_reverse"])
    pcr1 = graph.product("pcr1", "first-round PCR product", MaterialKind.MOLECULAR,
                         state=sc.pcr1_product())
    graph.react("pcr1_amplify", "First-round library PCR", StepKind.PCR,
                ((purified, InputRole.TEMPLATE), (f1, InputRole.STRUCTURAL),
                 (r1, InputRole.STRUCTURAL), (kapa, InputRole.CATALYST)), (pcr1,),
                geometry=PCRGeometry.exact(
                    sc.second_strand_duplex(), PCRPrimer("first-round forward", sc.PCR1_FORWARD),
                    PCRPrimer("first-round reverse", sc.PCR1_REVERSE),
                    statement="Both primers are printed in Methods (JATS line 279)."))
    pcr1_clean = _cleanup(graph, r, "ampure_12", "AMPure XP 1.2× cleanup", pcr1,
                          "AMPure XP 1.2×", "purified first-round product",
                          "AMPure supernatant", "not bound by the beads")
    kapa2 = graph.input("kapa_pcr2", "KAPA HiFi HotStart ReadyMix, 8–9 cycles, 10 µL of 2 nM template",
                        MaterialKind.CATALYST)
    f2 = graph.input("pcr2_forward", "second-round P5 forward primer", MaterialKind.OLIGO,
                     state=exact["pcr2_forward"])
    r2 = graph.input("pcr2_reverse", "second-round P7 / 8-mer index reverse primer",
                     MaterialKind.OLIGO)
    library = graph.product("pcr2", "indexed Seq-Scope library", MaterialKind.LIBRARY,
                            state=FINAL_LIBRARY)
    graph.react("pcr2_amplify", "Second-round indexing PCR", StepKind.PCR,
                ((pcr1_clean, InputRole.TEMPLATE), (f2, InputRole.STRUCTURAL),
                 (r2, InputRole.STRUCTURAL), (kapa2, InputRole.CATALYST)), (library,),
                geometry=PCRGeometry.exact(
                    sc.pcr1_product(),
                    PCRPrimer("second-round forward", sc.PCR2_FORWARD, anneal_nt=23),
                    PCRPrimer("second-round reverse",
                              sc.PCR2_REVERSE_5P + "N" * 8 + sc.PCR2_REVERSE_3P, anneal_nt=20),
                    statement="Both primers are printed in Methods (JATS line 281); the index is an unassigned 8-mer."))
    gel = _cleanup(graph, r, "gel_select", "Gel-excise 400–850 bp (Zymoclean)", library,
                   "agarose gel, Zymoclean Gel DNA Recovery", "400–850 bp library",
                   "out-of-range gel products", "outside the excised 400–850 bp window",
                   kind=StepKind.SELECTION, material=MaterialKind.LIBRARY)
    final = _cleanup(graph, r, "ampure_final", "AMPure XP 0.6–0.7× cleanup", gel,
                     "AMPure XP 0.6–0.7×", "final Seq-Scope library",
                     "AMPure supernatant", "not bound by the beads",
                     material=MaterialKind.LIBRARY)
    graph.sequence(final, SequencingPlan(library=FINAL_LIBRARY, primers=SEQ_PRIMERS,
                                         read_uses={"Read 1": "extract",
                                                    "Index 1 (i7)": "demultiplex",
                                                    "Read 2": "extract_then_include"}))
    graph.annotate(SourceDiscrepancy(
        "seed", "Seed end labels differ from the printed bases",
        (SourceClaim("STAR Methods, Generation of Seed HDMI-oligo Library (PMC8238917 JATS line 217)",
                     "backbone: (P5 sequence) (TR1) (HDMI) (HR1) (Oligo-dT) (DraI) (DraI-adapter) (P7 sequence)"),
         SourceClaim("Same section, HDMI-DraI and HDMI32-DraI sequences (JATS lines 221–225)",
                     "5′ CAAGCAGAAGACGGCATACGAGAT … GTGTAGATCTCGGTGGTCGCCGTATCATT 3′, i.e. the "
                     "Illumina P7 adapter and the reverse complement of P5")),
        "Segments are named for their bases (P7 adapter at the 5′ end, P5 complement at the 3′ "
        "end); the printed bases are used unchanged.",
        "The label/base identity is computed against the canonical Illumina sequences.",
        "The paper may label the ends by their side in the final library (TR1 side = P5 after "
        "PCR); it does not say so."))
    r.extend(("seed", "flow_cell", "cluster", "read1_drai", "first_seq_sbs", "first_seq",
              "drai", "drai_digest", "exo1_array", "exo1_array_digest", "naoh_array",
              "denature_array", "tissue", "rna_release", "release_rna", "capture", "maxima",
              "reverse_transcribe", "exo1_tissue", "exo1_tissue_digest", "proteinase_k",
              "digest_tissue", "naoh_rna", "strip_rna", "random_primer", "klenow",
              "second_strand_synthesis", "naoh_elute", "elute_second_strand", "kapa_pcr1",
              "pcr1_forward", "pcr1_reverse", "pcr1_amplify", "kapa_pcr2", "pcr2_forward",
              "pcr2_reverse", "pcr2_amplify", "ampure_final_kept"))
    graph.declare_reviewed_scope(r, AUTHORITY)
    return graph.require_publishable()


def protocol_graph():
    return _local_graph()


SPEC.protocol_graph = protocol_graph
