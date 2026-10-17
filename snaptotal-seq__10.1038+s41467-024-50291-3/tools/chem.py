"""snapTotal-seq: MALBAC/template-switch total-RNA amplification and tag library.

Authority: Niu, Luo & Zong, Nat. Commun. 15:5941 (2024), PMC11251146 JATS
Methods (cell lysis/RT/amplification and "Sequencing library construction")
and Supplementary Data 12 (MOESM15 workbook, "Sequence of DNA oligos used in
snapTotal-seq").  The page draws the paper's standard-Illumina-primer option.
"""
import illumina as il
import seqprimers as sp
from batch_ngs import seg
from chemdraw import Construct, feature, revcomp
from endprep import (ThreePrimeOverhangAdapter, cohesive_ligation_geometry,
                     dA_tailed_scene, repair_and_dA_tail)
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            SourceClaim, SourceDiscrepancy, StepKind)

TITLE = "snapTotal-seq — single-cell total-RNA sequencing"
NOTES = "01_snaptotal-seq.html"
DOI = "https://doi.org/10.1038/s41467-024-50291-3"
SOURCE = (f'Defining source: <a href="{DOI}">Niu, Luo &amp; Zong, <i>Nature Communications</i> '
          '15:5941 (2024)</a>, Methods and Supplementary Data 12.')
SUMMARY = ("Low-temperature MALBAC priming and template switching amplify total RNA; a later "
           "double-strand conversion adds an 8-nt cell barcode to only one cDNA end. After "
           "pooling and sonication, a second conversion adds a handle to tagged fragments, a "
           "Y adapter is ligated and PCR keeps the tagged end opposite P5.")

# Supplementary Data 12 (MOESM15), transcribed verbatim (spaces removed).
GAT27 = "GTGAGTGATGGTTGAGGATGTGTGGAG"
RT_PRIMERS = {"GAT27-5N-3G": GAT27 + "NNNNN" + "GGG",
              "GAT27-5N-3T": GAT27 + "NNNNN" + "TTT",
              "GAT27-5N-dT": GAT27 + "NNNNN" + "T" * 12}
CELL_BARCODE_PRIMER = "GTGTGCTCTTCCGATCT" + "NNNNNNNN" + "AGGAG" + "AGTGTGAGTGATGGTTGAGGATGTGTGGAG"
ADAPTER_DSC = "GCACGACATCTGCTAACGCAGTA" + "GTGTGCTCTTCCGATCT"
TRUSEQ_P5_INDEX = ("AATGATACGGCGACCACCGAGATCTACAC" + "NNNNNNNN"
                   + "ACACTCTTTCCCTACACGACGCTCTTCCGATCT")
Y_COMPLEMENT = "GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCTCGTATGCCGTCTTCTGCTTG"  # printed P*
ADAPTER_PCR = "GCACGACATCTGCTAACGCAGTA"
ILLUMINA_P5 = "AATGATACGGCGACCACCGAGA"
ADAPTER_P7 = "CAAGCAGAAGACGGCATACGAGATGCACGACATCTGCTAACGCAGTA"
TRUSEQ_P7_ARM = "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"
TRUSEQ_P7_EXAMPLE = "CAAGCAGAAGACGGCATACGAGAT" + "CGAGTAAT" + TRUSEQ_P7_ARM  # example index
CUSTOM_READ2 = "AGAGGTGAGTGAGTGATGGTTGAGGATGTGTGGAG"
CUSTOM_I7 = "CTCCACACATCCTCAACCATCACTCACTCACCTCT"

# The spacer between the cell barcode and the UMI is the cell-barcode primer's
# 3' 35 nt (AGGAG + 30 nt), read in Read 2 as bases 9-43 (Methods: barcode 1-8, UMI 44-48).
_SPACER = CELL_BARCODE_PRIMER[-35:]
if not _SPACER.endswith(GAT27):
    raise AssertionError("cell-barcode primer must end in GAT27")
if revcomp(CUSTOM_READ2) != CUSTOM_I7:
    raise AssertionError("custom i7 primer is printed as the reverse complement of Read 2")

FINAL_LIBRARY = Construct([
    seg("P5", il.P5, "p5"),
    seg("i5 library barcode", "J" * 8, "cbc", placeholder=True,
        feature=feature("library_barcode_i5", "sample_index", "unknown",
                        note="Supplementary Data 12 'cell barcode 2 (library barcode)' NNNNNNNN")),
    seg("Read 1 arm", il.TRUSEQ_READ1, "r1"),
    seg("captured RNA sequence (with RT-primer 3′ anchor)", "X" * 42, placeholder=True),
    seg("UMI", "U" * 5, "umi", placeholder=True, feature=feature("umi", "umi", "random")),
    seg("GAT27/AGGAG spacer", revcomp(_SPACER)),
    seg("cell barcode", "B" * 8, "cbc", placeholder=True,
        feature=feature("cell_barcode", "cell_barcode", "unknown",
                        note="printed only as NNNNNNNN in Supplementary Data 12")),
    seg("Read 2 / Index 1 arm", revcomp(TRUSEQ_P7_ARM), "r2"),
    seg("i7 reverse complement", "I" * 8, "cbc", placeholder=True,
        feature=feature("sample_index_i7", "sample_index", "unknown",
                        note="CGAGTAAT is printed only as an example index")),
    seg("P7 reverse complement", il.P7_RC, "p7"),
], name="snapTotal-seq standard-primer library")
SEQ_PRIMERS = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["I2"], sp.TRUSEQ["R2"])
_problems = sp.verify(FINAL_LIBRARY, SEQ_PRIMERS,
                      required_roles=tuple(p.role for p in SEQ_PRIMERS))
if _problems:
    raise ValueError("; ".join(_problems))

FINAL_CAPTION = ("Standard-Illumina option (TruSeq-P7 primer in the final PCR): Read 1 reads the "
                 "captured RNA from the ligated Y-adapter end; Read 2 reads the 8-nt cell barcode "
                 "(bases 1–8), the 35-nt primer-derived spacer and the 5-nt UMI (bases 44–48).")
SEQUENCING_INTRO = ("Drawn for the standard-primer alternative. The original NextSeq 500 runs used "
                    "Adapter-P7 and custom Read 2 / i7 primers, whose site is discussed in the "
                    "discrepancy notice above.")
Y_ADAPTER = ThreePrimeOverhangAdapter("snapTotal-seq Y-adapter ligating stem",
                                      "GCTCTTCCGATCT", "GATCGGAAGAGC")


def oligos():
    rows = [f"<b>{name}</b> — 5′-{seq}-3′ (NNNNN: UMI)" for name, seq in RT_PRIMERS.items()]
    rows += [f"<b>GAT27</b> (cDNA PCR) — 5′-{GAT27}-3′",
             f"<b>Cell barcode primer</b> — 5′-{CELL_BARCODE_PRIMER}-3′ (NNNNNNNN: cell barcode)",
             f"<b>Adapter-DSC</b> — 5′-{ADAPTER_DSC}-3′",
             f"<b>Truseq-P5-index primer</b> (Y-adapter strand) — 5′-{TRUSEQ_P5_INDEX}-3′",
             f"<b>Y-shape complementary sequence</b> — 5′-P*{Y_COMPLEMENT}-3′ (P*: as printed, 5′ phosphate)",
             f"<b>Adapter-PCR</b> — 5′-{ADAPTER_PCR}-3′; <b>Illumina-P5</b> — 5′-{ILLUMINA_P5}-3′",
             f"<b>Adapter-P7</b> (custom-sequencing final PCR) — 5′-{ADAPTER_P7}-3′",
             f"<b>Truseq-P7</b> (standard final PCR; CGAGTAAT is an example index) — 5′-{TRUSEQ_P7_EXAMPLE}-3′",
             f"<b>Custom Read 2</b> — 5′-{CUSTOM_READ2}-3′; <b>custom Index i7</b> — 5′-{CUSTOM_I7}-3′"]
    return rows


def protocol_graph():
    src = f"Niu et al. 2024 ({DOI}) Methods"
    graph = ProtocolGraph("snapTotal-seq", authority=f"Defining source: {DOI}")
    rna = graph.input("total_rna", "single-cell lysate total RNA (42 °C, 3.5 min lysis)",
                      MaterialKind.BIOLOGICAL, primary=True)
    malbac = graph.input("malbac_primer", "GAT27-5N-3G / -3T / -dT primer mix (in lysis buffer)",
                         MaterialKind.OLIGO, state=tuple(RT_PRIMERS.values()))
    rt = graph.input("reverse_transcriptase", "Maxima H Minus reverse transcriptase",
                     MaterialKind.CATALYST)
    cdna = graph.product("first_strand_cdna", "multiply primed total-RNA cDNA", MaterialKind.MOLECULAR)
    graph.react("reverse_transcription",
                "10 cycles 8→15→20→30→42→50 °C multiple annealing and extension",
                StepKind.REVERSE_TRANSCRIPTION,
                ((rna, InputRole.TEMPLATE), (malbac, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (cdna,))
    switched = graph.product("switched_cdna", "template-switched cDNA with GAT27 at both ends",
                             MaterialKind.MOLECULAR)
    graph.react("template_switch",
                "Template-switch onto the same MALBAC primers (no separate TSO); 74 °C inactivation",
                StepKind.TEMPLATE_SWITCHING,
                ((cdna, InputRole.TEMPLATE), (malbac, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (switched,))
    gat27 = graph.input("gat27_primer", "GAT27 cDNA amplification primer", MaterialKind.OLIGO,
                        state=GAT27)
    gotaq = graph.input("gotaq", "GoTaq Flexi DNA polymerase", MaterialKind.CATALYST)
    amplified = graph.product("amplified_cdna", "GAT27-amplified cDNA", MaterialKind.MOLECULAR)
    graph.react("malbac_pcr", "Amplify cDNA with GAT27, 23–26 cycles (63 °C anneal)", StepKind.PCR,
                ((switched, InputRole.TEMPLATE), (gat27, InputRole.STRUCTURAL),
                 (gotaq, InputRole.CATALYST)), (amplified,),
                geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY, reason="A single GAT27 primer amplifies cDNA through identical "
                    "ends; the intermediate cDNA state is not drawn base by base.",
                    symmetric_species=True))
    barcode_primer = graph.input("barcode_dsc_primer", "cell barcode primer (8-nt cell barcode)",
                                 MaterialKind.OLIGO, state=CELL_BARCODE_PRIMER)
    deep_vent = graph.input("polymerase", "Deep Vent (exo-) DNA polymerase", MaterialKind.CATALYST)
    barcoded = graph.product("barcoded_cdna", "one-end cell-barcoded cDNA", MaterialKind.MOLECULAR)
    graph.react("barcode_conversion",
                "Double-strand conversion: 15 cycles 63/72 °C, EDTA stop",
                StepKind.EXTENSION,
                ((amplified, InputRole.TEMPLATE), (barcode_primer, InputRole.STRUCTURAL),
                 (deep_vent, InputRole.CATALYST)), (barcoded,))
    shear = graph.input("shearing_condition", "pool 4 µL/cell, 1.1X AMPure, Covaris S220 to "
                        "150–250 bp, 1.8X AMPure", MaterialKind.CONDITION)
    sheared = graph.product("sheared_fragments", "pooled sheared cDNA fragments",
                            MaterialKind.MOLECULAR)
    graph.react("shear", "Pool, purify and sonicate the barcoded cDNA", StepKind.FRAGMENTATION,
                ((barcoded, InputRole.TEMPLATE), (shear, InputRole.CONDITION)), (sheared,))
    adapter_dsc = graph.input("adapter_dsc", "Adapter-DSC primer", MaterialKind.OLIGO,
                              state=ADAPTER_DSC)
    handled = graph.product("handle_fragments",
                            "fragments; tagged ones gain the Adapter-DSC handle (untagged carried along)",
                            MaterialKind.MOLECULAR)
    graph.react("enrich", "Second double-strand conversion: 6 cycles 56/72 °C, 1.4X AMPure",
                StepKind.EXTENSION,
                ((sheared, InputRole.TEMPLATE), (adapter_dsc, InputRole.STRUCTURAL),
                 (deep_vent, InputRole.CATALYST)), (handled,))
    klenow = graph.input("tailing_polymerase", "Klenow fragment (3′→5′ exo-), dATP",
                         MaterialKind.CATALYST)
    tailed_construct = repair_and_dA_tail(Construct([seg("tagged cDNA fragment", "X" * 36,
                                                         placeholder=True)]))
    tailed = graph.product("tailed_fragments", "dA-tailed fragments", MaterialKind.MOLECULAR,
                           state=dA_tailed_scene(tailed_construct))
    graph.react("da_tail", "Add 3′ dA (room temperature, 30 min; no separate end repair is "
                "described), 1.4X AMPure", StepKind.TAILING,
                ((handled, InputRole.TEMPLATE), (klenow, InputRole.CATALYST)), (tailed,))
    adapter = graph.input("y_adapter", "Y-shape adapter (Truseq-P5-index primer + P*-complement)",
                          MaterialKind.ADAPTER, state=(TRUSEQ_P5_INDEX, Y_COMPLEMENT))
    ligase = graph.input("ligase", "Quick Ligase", MaterialKind.CATALYST)
    precursor = graph.product("adapter_ligated", "Y-adapter-ligated fragments",
                              MaterialKind.MOLECULAR)
    graph.react("adapter_ligation", "Ligate the Y adapter (RT, 20 min), EDTA, 1.2X AMPure",
                StepKind.LIGATION,
                ((tailed, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (precursor,),
                geometry=cohesive_ligation_geometry(
                    tailed_construct, Y_ADAPTER,
                    statement="Supplementary Data 12 prints both Y-adapter strands; the "
                    "Truseq-P5-index strand ends in 3′-T and the complement carries P*."))
    adapter_pcr = graph.input("adapter_pcr_primer", "Adapter-PCR primer", MaterialKind.OLIGO,
                              state=ADAPTER_PCR)
    p5 = graph.input("p5_primer", "Illumina-P5 primer", MaterialKind.OLIGO,
                     state=ILLUMINA_P5)
    amplified_lib = graph.product("ligation_pcr_product", "handle-to-P5 amplified fragments",
                                  MaterialKind.MOLECULAR)
    graph.react("ligation_pcr", "Amplify ligation products: Adapter-PCR + Illumina-P5, 10 cycles",
                StepKind.PCR,
                ((precursor, InputRole.TEMPLATE), (adapter_pcr, InputRole.STRUCTURAL),
                 (p5, InputRole.STRUCTURAL), (deep_vent, InputRole.CATALYST)), (amplified_lib,),
                geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY, reason="The intermediate handle-to-P5 amplicon is not drawn; "
                    "Adapter-PCR anneals inside the Adapter-DSC handle."))
    dsn = graph.input("dsn", "duplex-specific nuclease (Evrogen): 100 ng, 95 °C denature, "
                      "80 °C 3 h re-annealing, 80 °C 15 min", MaterialKind.CATALYST)
    normalized = graph.product("dsn_treated", "DSN-normalized fragments", MaterialKind.MOLECULAR)
    digested = graph.product("dsn_digested", "re-annealed abundant (ribosomal) duplexes digested by DSN",
                             MaterialKind.WASTE)
    graph.react("dsn_treatment", "Remove ribosomal reads with DSN, EDTA quench, 1.2X AMPure",
                StepKind.DIGESTION,
                ((amplified_lib, InputRole.TEMPLATE), (dsn, InputRole.CATALYST)),
                (normalized, digested))
    graph.discard(digested, "DSN digests re-annealed duplexes; the products are purified with "
                  "1.2X AMPure XP beads (Methods)")
    enriched = graph.product("enriched_fragments", "undigested fragments re-amplified",
                             MaterialKind.MOLECULAR)
    graph.react("enrichment_pcr", "Enrichment PCR: Adapter-DSC + Illumina-P5, 5 cycles, 1.2X AMPure",
                StepKind.PCR,
                ((normalized, InputRole.TEMPLATE), (adapter_dsc, InputRole.STRUCTURAL),
                 (p5, InputRole.STRUCTURAL), (deep_vent, InputRole.CATALYST)), (enriched,),
                geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY, reason="The intermediate Adapter-DSC-to-P5 amplicon is not drawn."))
    truseq_p7 = graph.input("p7_primer", "Truseq-P7 primer (8-nt i7 index)", MaterialKind.OLIGO,
                             state=TRUSEQ_P7_EXAMPLE)
    library = graph.product("library", FINAL_LIBRARY.name, MaterialKind.LIBRARY, state=FINAL_LIBRARY)
    graph.react("library_pcr", "Final PCR, standard option: Truseq-P7 + Illumina-P5, 5 cycles",
                StepKind.PCR,
                ((enriched, InputRole.TEMPLATE), (p5, InputRole.STRUCTURAL),
                 (truseq_p7, InputRole.STRUCTURAL), (deep_vent, InputRole.CATALYST)), (library,),
                geometry=PCRGeometry.exact(
                    FINAL_LIBRARY, PCRPrimer("Illumina-P5", ILLUMINA_P5),
                    PCRPrimer("Truseq-P7", TRUSEQ_P7_EXAMPLE),
                    statement="Methods and Supplementary Data 12: Truseq-P7 replaces Adapter-P7 "
                    "and anneals to the cell-barcode primer's 5′ GTGTGCTCTTCCGATCT."))
    graph.sequence(library, SequencingPlan(
        library=FINAL_LIBRARY, primers=SEQ_PRIMERS,
        acquired_roles=tuple(p.role for p in SEQ_PRIMERS),
        read_uses={p.role: ("demultiplex" if "Index" in p.role else "include")
                   for p in SEQ_PRIMERS}))
    graph.annotate(SourceDiscrepancy(
        "library", "Custom Read 2 / Index i7 primers versus the cell-barcode primer",
        (SourceClaim("Methods (sequencing paragraph) and Supplementary Data 12, custom primers",
                     f"Read 2 5′-{CUSTOM_READ2}-3′; Index i7 5′-{CUSTOM_I7}-3′ (its reverse complement)"),
         SourceClaim("Methods (DSC paragraph) and Supplementary Data 12 step 1, cell barcode primer",
                     f"…NNNNNNNN-{_SPACER}-3′: the 8 nt before GAT27 are AGGAGAGT, "
                     "where the custom Read 2 primer has AGAGGTGA")),
        "The page draws only the paper's standard-primer option (Truseq-P7, standard "
        "Illumina primers), whose Read 2 layout (barcode 1–8, UMI 44–48) the Methods state.",
        "Computed from the printed sequences: the custom primers agree with the template "
        "only over GAT27; the custom i7 primer's 3′ 8 nt therefore do not pair with the "
        "drawn barcode-side spacer.",
        "The source does not explain how the custom index read primes on this template or "
        "which sequence is intended; the custom-primer run layout is not drawn."))
    graph.declare_reviewed_scope((
        "total_rna", "malbac_primer", "reverse_transcriptase", "reverse_transcription",
        "template_switch", "gat27_primer", "gotaq", "malbac_pcr",
        "barcode_dsc_primer", "polymerase", "barcode_conversion", "shearing_condition", "shear",
        "adapter_dsc", "enrich", "tailing_polymerase", "da_tail", "y_adapter", "ligase",
        "adapter_ligation", "adapter_pcr_primer", "p5_primer", "ligation_pcr", "dsn",
        "dsn_treatment", "enrichment_pcr", "p7_primer", "library_pcr", "library"), src)
    return graph.require_publishable()
