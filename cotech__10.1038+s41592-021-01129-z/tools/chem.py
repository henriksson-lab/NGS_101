"""CoTECH (Xiong et al., Nature Methods 2021) — role-level declaration.

Source scope (round 2026-10-17-breadth-first-08):
* Published Supplementary Table 1 (Springer MOESM3 XLSX): custom primers, PA-Tn5 T5/T7
  barcode adaptors, barcoded oligo-dT RT primers, Nextera i5 and TruSeq i5/i7 index
  primers.  Oligos below are transcribed from it verbatim.
* The published article text is paywalled and was not read.  Step order (pool and
  redistribute, preamplify, purify and split, cDNA Tn5 tagmentation, two-round DNA
  PCR, second-round index added during library PCR) comes from the bioRxiv preprint
  10.1101/2020.10.15.339226v1 Results lines 150-172 and Figure 1 legend lines 347-358.
  Those steps are labelled provisional/inferred until the published Methods are read.
"""
from types import SimpleNamespace
from pcr import PCRGeometry
from protocol_model import chain, role_library
from multimodal_protocols import Definition, InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind, ident, molecule, role

DEFINITION = Definition("CoTECH — chromatin occupancy and transcriptome", "10.1038/s41592-021-01129-z",
        "Xiong et al., Nature Methods (2021)", "Antibody-targeted pA–Tn5 carrying a well-specific T5/T7 adaptor combination and a well-barcoded poly(dT) RT primer mark chromatin and RNA in the same first-round well; after pooling and redistribution, each secondary well's products receive a second-round index during library PCR.",
        "Separate chromatin and RNA libraries; a cell is the first-round well barcode(s) plus the second-round PCR index of its secondary well.", "01_cotech.html")

PREPRINT = "bioRxiv 10.1101/2020.10.15.339226v1 (preprint; published Methods not read)"
TABLE_S1 = "Supplementary Table 1 (MOESM3), published"

# Verbatim from published Supplementary Table 1 (modification notation as printed).
OLIGOS = (
    ("Mosaic-A_primer (MEA)", "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"),
    ("common annealing primer", "5Phos/CTGTCTCTTATACACATCT"),
    ("TSO", "AAGCAGTGGTATCAACGCAGAGTACATrGrG+G"),
    ("IS primer", "AAGCAGTGGTATCAACGCAGAGT"),
    ("P2 primer", "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC"),
    ("connector A", "TCGTCGGCAGCGTCTCCACGC"),
    ("connector B", "GTCTCGTGGGCTCGGCTGTCC"),
    ("Truseq 1st round forward primer", "ACACTCTTTCCCTACACGACGCTCTTCCGATCTTCGTCGGCAGCGTCTCCACGC"),
    ("Truseq 1st round reverse primer", "GACTGGAGTTCAGACGTGTGCTCTTCCGATCTGTCTCGTGGGCTCGGCTGTCCC"),
    ("T5-1 (one of 24 PA-T5 barcode primers)", "TCGTCGGCAGCGTCTCCACGCTATAGCCTGCGATCGAGGACGGCAGATGTGTATAAGAGACAG"),
    ("T7-1 (one of 25 PA-T7 barcode primers)", "GTCTCGTGGGCTCGGCTGTCCCTGTCCCGAGTAATCACCGTCTCCGCCTCAGATGTGTATAAGAGACAG"),
    ("RT-barcode_01 (one of 96 barcoded oligo dT)", "TCAGACGTGTGCTCTTCCGATCTAACGTGATDDDDDDDDttttttttttttttttttttttttt"),
)

R1 = ident("first-round RT well barcode", "cotech_round1_rt", "cell_barcode")
T5 = ident("T5 barcode (first-round well, with T7)", "cotech_round1_t5", "cell_barcode")
T7 = ident("T7 barcode (first-round well, with T5)", "cotech_round1_t7", "cell_barcode")
R2 = ident("second-round PCR index of the secondary well", "cotech_round2", "cell_barcode")
UMI = ident("8-nt UMI (D window)", "umi", "umi", length=8, tag="umi")


def _cotech_graph(d):
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    cells = graph.input("cells", "antibody-bound permeabilized cells containing poly(A) RNA",
                        MaterialKind.BIOLOGICAL, state=_cotech().initial, primary=True)
    transposome = graph.input("barcoded_patn5",
                              "pA–Tn5 assembled with a well-specific T5/T7 barcode-primer combination, each annealed to the 5Phos common annealing primer",
                              MaterialKind.ADAPTER, source=TABLE_S1 + "; " + PREPRINT)
    tn5 = graph.input("tn5", "protein-A–Tn5 transposase", MaterialKind.CATALYST)
    tagged = graph.product("tagged_cells", "cells with target-proximal T5/T7-barcoded chromatin",
                           MaterialKind.MOLECULAR, state=molecule("tagged chromatin branch", (
                               role("connector A"), T5, role("mosaic end (T5 side)"),
                               role("epitope-proximal DNA"), role("mosaic end (T7 side)"),
                               T7, role("connector B"))))
    graph.react("tagment", "Tagment antibody-targeted chromatin in the first-round well",
                StepKind.TAGMENTATION,
                ((cells, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged,))
    rt_primer = graph.input("rt_primer", "first-round well-barcoded oligo dT with 8-nt UMI (RT-barcode_01–96)",
                            MaterialKind.OLIGO, source=TABLE_S1)
    tso = graph.input("tso", "template-switch oligo AAGCAGTGGTATCAACGCAGAGTACATrGrG+G",
                      MaterialKind.OLIGO, source=TABLE_S1)
    rt = graph.input("reverse_transcriptase", "reverse transcriptase", MaterialKind.CATALYST)
    paired = graph.product("paired_branches", "cells containing first-round-barcoded chromatin and full-length cDNA",
                           MaterialKind.MOLECULAR,
                           state=molecule("full-length cDNA branch", (
                               role("TSO / IS site"), role("transcript cDNA"), role("poly(dT)"),
                               UMI, R1, role("TruSeq Read 2 handle"))))
    graph.react("reverse_transcribe", "Reverse-transcribe poly(A) RNA with the well's barcoded oligo dT and template switch",
                StepKind.REVERSE_TRANSCRIPTION,
                ((tagged, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
                 (tso, InputRole.STRUCTURAL), (rt, InputRole.CATALYST)), (paired,))
    second_well = graph.input("second_well", "pooling and redistribution into secondary 96-/384-well plates",
                              MaterialKind.CONDITION, source=PREPRINT)
    redistributed = graph.product("redistributed", "redistributed cells, one secondary well each",
                                  MaterialKind.MOLECULAR, state=paired.state, source=PREPRINT)
    graph.react("redistribute", "Pool and redistribute cells into secondary wells (each well is followed separately)",
                StepKind.GENERIC,
                ((paired, InputRole.TEMPLATE), (second_well, InputRole.CONDITION)),
                (redistributed,), changes_molecule=False, source=PREPRINT)
    pre_f = graph.input("preamp_forward", "INFERRED — preamplification forward primers (connector A / IS primer are Table S1 candidates; assignment not read)",
                        MaterialKind.OLIGO, source=PREPRINT, inferred=True)
    pre_r = graph.input("preamp_reverse", "INFERRED — preamplification reverse primers (connector B / P2 primer are Table S1 candidates; assignment not read)",
                        MaterialKind.OLIGO, source=PREPRINT, inferred=True)
    pre_pol = graph.input("preamp_polymerase", "preamplification polymerase", MaterialKind.CATALYST,
                          source=PREPRINT)
    preamped = graph.product("preamplified", "preamplified cDNA and tagmented chromatin from one secondary well",
                             MaterialKind.MOLECULAR, state=paired.state, source=PREPRINT)
    graph.react("preamplify", "Preamplify cDNA and tagmented chromatin in each secondary well",
                StepKind.PCR,
                ((redistributed, InputRole.TEMPLATE), (pre_f, InputRole.STRUCTURAL), (pre_r, InputRole.STRUCTURAL),
                 (pre_pol, InputRole.CATALYST)), (preamped,), source=PREPRINT,
                geometry=PCRGeometry.unavailable(
                    role_library("CoTECH preamplification product", (role("chromatin or cDNA"),)),
                    reason="The preprint names preamplification but its primers and cycling were not read; published Methods unread."))
    purified = graph.product("purified", "AMPure XP-purified preamplified DNA", MaterialKind.MOLECULAR,
                             state=paired.state, source=PREPRINT)
    cleanup_waste = graph.product("cleanup_waste", "bead supernatant (primers, short products)",
                                  MaterialKind.WASTE, source=PREPRINT)
    graph.react("purify", "Purify with AMPure XP beads", StepKind.CLEANUP,
                ((preamped, InputRole.TEMPLATE),), (purified, cleanup_waste),
                changes_molecule=False, source=PREPRINT)
    graph.discard(cleanup_waste, "removed with the bead supernatant")
    chromatin = graph.product("chromatin_branch", "chromatin half", MaterialKind.MOLECULAR,
                              state=tagged.state, source=PREPRINT)
    cdna_half = graph.product("rna_branch", "cDNA half", MaterialKind.MOLECULAR,
                              state=paired.state, source=PREPRINT)
    graph.react("split_libraries", "Split into halves for DNA and RNA library preparation",
                StepKind.SPLIT, ((purified, InputRole.TEMPLATE),), (chromatin, cdna_half),
                changes_molecule=False, source=PREPRINT)
    reviewed = ["cells", "barcoded_patn5", "tn5", "tagment", "rt_primer", "tso",
                "reverse_transcriptase", "reverse_transcribe", "second_well",
                "redistribute", "preamp_forward", "preamp_reverse", "preamp_polymerase", "preamplify",
                "purify", "split_libraries"]

    # DNA: common TruSeq first-round primers then well-specific TruSeq i5/i7 second round.
    dna1_f = graph.input("dna_round1_forward", "Truseq 1st round forward primer (Read 1 + connector A)",
                         MaterialKind.OLIGO, source=TABLE_S1)
    dna1_r = graph.input("dna_round1_reverse", "Truseq 1st round reverse primer (Read 2 + connector B)",
                         MaterialKind.OLIGO, source=TABLE_S1)
    dna_pol = graph.input("dna_polymerase", "DNA library PCR polymerase", MaterialKind.CATALYST)
    dna_mid = molecule("chromatin first-round PCR product", (
        role("TruSeq Read 1 site"), role("connector A"), T5, role("mosaic end (T5 side)"),
        role("chromatin insert"), role("mosaic end (T7 side)"), T7, role("connector B"),
        role("TruSeq Read 2 site")))
    dna_round1 = graph.product("dna_round1", "chromatin first-round PCR product",
                               MaterialKind.MOLECULAR, state=dna_mid, source=TABLE_S1 + "; " + PREPRINT)
    graph.react("dna_pcr1", "First-round chromatin PCR with common TruSeq adaptor primers",
                StepKind.PCR, ((chromatin, InputRole.TEMPLATE), (dna1_f, InputRole.STRUCTURAL), (dna1_r, InputRole.STRUCTURAL),
                               (dna_pol, InputRole.CATALYST)), (dna_round1,),
                source=TABLE_S1 + "; " + PREPRINT,
                geometry=PCRGeometry.unavailable(
                    role_library("chromatin first-round product", (role("chromatin insert"),)),
                    reason="Primer bases are printed but the T5/T7 adaptor boundaries and cycling were not modelled base by base."))
    dna_idx = graph.input("dna_i5_primer", "well-specific TruSeq i5 index primer (Truseq-501–516)",
                          MaterialKind.OLIGO, source=TABLE_S1)
    dna_i7 = graph.input("dna_i7_primer", "well-specific TruSeq i7 index primer (Truseq-701–708)",
                         MaterialKind.OLIGO, source=TABLE_S1)
    dna_final = role_library("CoTECH chromatin library roles", (
        R2, role("TruSeq Read 1 site"), role("connector A"), T5, role("mosaic end (T5 side)"),
        role("target-proximal chromatin insert"), role("mosaic end (T7 side)"),
        T7, role("connector B"), role("TruSeq Read 2 site"), role("i7 index")))
    dna_library = graph.product("chromatin_library", "CoTECH chromatin sequencing library",
                                MaterialKind.LIBRARY, state=dna_final, source=TABLE_S1 + "; " + PREPRINT)
    graph.react("dna_pcr2", "Second-round PCR adds the secondary-well TruSeq i5/i7 index pair",
                StepKind.PCR, ((dna_round1, InputRole.TEMPLATE), (dna_idx, InputRole.STRUCTURAL), (dna_i7, InputRole.STRUCTURAL),
                               (dna_pol, InputRole.CATALYST)), (dna_library,),
                source=TABLE_S1 + "; " + PREPRINT,
                geometry=PCRGeometry.unavailable(
                    dna_final, reason="Index primer family printed; the index pair per secondary well is not assigned."))
    graph.sequence(dna_library, SequencingPlan(limitation=
        "Read lengths and which reads cover the T5/T7 barcodes were not read from the published Methods."))

    # RNA: Tn5 tagmentation of cDNA; i5/i7 PCR amplifies the 3' (RT-primer) end.
    rna_tn5 = graph.input("cdna_tn5", "INFERRED — Tn5 loaded with MEA/common annealing primer (Table S1 lists them; loading not read)",
                          MaterialKind.ADAPTER, source=TABLE_S1, inferred=True)
    rna_tn5_enzyme = graph.input("cdna_tn5_enzyme", "Tn5 transposase", MaterialKind.CATALYST, source=PREPRINT)
    tagmented = graph.product("cdna_tagmented", "tagmented cDNA fragments", MaterialKind.MOLECULAR,
                              state=molecule("tagmented 3' cDNA fragment", (
                                  role("mosaic end"), role("cDNA fragment"), role("poly(dT)"),
                                  UMI, R1, role("TruSeq Read 2 handle"))), source=PREPRINT)
    graph.react("cdna_tagment", "Tagment preamplified full-length cDNA with Tn5",
                StepKind.TAGMENTATION, ((cdna_half, InputRole.TEMPLATE), (rna_tn5, InputRole.STRUCTURAL),
                                        (rna_tn5_enzyme, InputRole.CATALYST)), (tagmented,), source=PREPRINT)
    rna_idx = graph.input("rna_i5_primer", "Nextera i5 index primer (i5-N501–N526)",
                          MaterialKind.OLIGO, source=TABLE_S1)
    rna_i7 = graph.input("rna_i7_primer", "TruSeq i7 index primer (Truseq-701–708)",
                         MaterialKind.OLIGO, source=TABLE_S1)
    rna_pol = graph.input("rna_polymerase", "RNA library PCR polymerase", MaterialKind.CATALYST)
    rna_final = role_library("CoTECH transcriptome library roles", (
        R2, role("Nextera Read 1 / mosaic end"), role("cDNA fragment"), role("poly(dT)"),
        UMI, R1, role("TruSeq Read 2 site"), role("i7 index")))
    rna_library = graph.product("rna_library", "CoTECH 3' transcriptome sequencing library",
                                MaterialKind.LIBRARY, state=rna_final, source=TABLE_S1 + "; " + PREPRINT)
    graph.react("rna_pcr", "Amplify the 3' cDNA end and add the second-round RNA index",
                StepKind.PCR, ((tagmented, InputRole.TEMPLATE), (rna_idx, InputRole.STRUCTURAL), (rna_i7, InputRole.STRUCTURAL),
                               (rna_pol, InputRole.CATALYST)), (rna_library,), source=PREPRINT,
                geometry=PCRGeometry.unavailable(
                    rna_final, reason="Which of i5/i7 carries the secondary-well index, and the selected pair, were not read."))
    graph.sequence(rna_library, SequencingPlan(limitation=
        "Read lengths and the read covering RT barcode/UMI were not read from the published Methods."))
    reviewed += ["dna_round1_forward", "dna_round1_reverse", "dna_polymerase", "dna_pcr1", "dna_round1",
                 "dna_i5_primer", "dna_i7_primer", "dna_pcr2", "chromatin_library", "cdna_tn5",
                 "cdna_tn5_enzyme", "cdna_tagment", "rna_i5_primer", "rna_i7_primer", "rna_polymerase",
                 "rna_pcr", "rna_library"]
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def protocol_graph():
    return _cotech_graph(DEFINITION)


def _cotech():
    return chain(molecule("antibody-bound permeabilized cells", (
        role("target chromatin epitope"), role("poly(A) RNA"))), (
        ("Tagment with pA–Tn5 carrying the well's T5/T7 barcode combination",
         molecule("tagged chromatin branch", (T5, role("epitope-proximal DNA"), T7)),
         "The antibody localizes pA–Tn5; the T5/T7 combination identifies the first-round well."),
        ("Reverse-transcribe with the well's barcoded oligo dT and template switch",
         molecule("full-length cDNA", (role("TSO site"), role("transcript cDNA"), UMI, R1)),
         "Chromatin and RNA receive first-round barcodes of the same well before pooling."),
        ("Pool, redistribute, preamplify, purify and split (preprint order)",
         molecule("preamplified products", (role("chromatin half"), role("cDNA half"))),
         "Provisional order from the bioRxiv preprint."),
        ("Library PCR adds the secondary-well index to each branch",
         molecule("matched CoTECH libraries", (T5, T7, R1, R2, role("chromatin or 3' cDNA insert"))),
         "The first-round barcode plus the second-round PCR index identify a cell."),
    ))


d = DEFINITION
FINAL_LIBRARY = role_library(d.title + " final library roles", (R1, role("chromatin or cDNA insert")), index=R2)
SPEC = SimpleNamespace(
    TITLE=d.title,
    NOTES=d.notes,
    SOURCE=f'Defining source: <a href="https://doi.org/{d.doi}">{d.citation}</a>.',
    SUMMARY=d.summary,
    CAVEAT=("The published article text was not available to this review; oligos come from its published "
            "Supplementary Table 1, while the order pool → redistribute → preamplify → purify/split → "
            "library-PCR indexing is taken provisionally from the bioRxiv preprint (10.1101/2020.10.15.339226). "
            "Preamplification primers and Tn5 loading for the cDNA branch are INFERRED from the table's primer list."),
    FINAL_LIBRARY=FINAL_LIBRARY,
    SEQ_PRIMERS=(),
    FINAL_CAPTION=d.endpoint,
    SEQUENCING_INTRO=d.endpoint,
    READ_LENGTHS=None,
    SEQUENCING_ENDING=d.endpoint,
    SEQUENCING_UNAVAILABLE="Read lengths and the reads covering the first-round barcodes and UMI were not read from the published Methods; no universal read layout is drawn.",
    workflow=_cotech,
    oligos=lambda: tuple(f"{name}: 5′-{seq}-3′" for name, seq in OLIGOS),
)
SPEC.protocol_graph = protocol_graph
globals().update(vars(SPEC))
