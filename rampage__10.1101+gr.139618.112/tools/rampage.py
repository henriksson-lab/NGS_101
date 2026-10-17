"""RAMPAGE promoter-profiling workflow.

Source scope (bounded pass 2026-10-17-breadth-first-16):
- Batut et al. 2013 Genome Research (PMC3530677 JATS) Methods "Library preparation
  and sequencing", "Sequencing reads alignment" and duplicate removal.
- Batut & Gingeras 2013 Current Protocols 25B.11 (PMC4372803 author manuscript JATS):
  Synthetic oligonucleotides, Sequence barcodes and Basic Protocol steps 1-42.
- Typeset Current Protocols Supplement 104 pages 25B.11.3-25B.11.5 as reproduced in
  Philippe Batut's CSHL thesis (Appendix 1, PDF pp. 122-124): materials list and
  Table 25B.11.1 (the 40 rampage_TS_01-40 TSOs).
"""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
import illumina as il
import seqprimers as sp
from batch_ngs import seg
from chemdraw import Construct, feature, revcomp
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)

TITLE = "RAMPAGE — paired-end promoter profiling"
NOTES = "01_rampage.html"
SOURCE = 'Defining paper: <a href="https://doi.org/10.1101/gr.139618.112">Batut et al. (2013)</a>; detailed protocol: <a href="https://doi.org/10.1002/0471142727.mb25b11s104">Batut &amp; Gingeras (2013)</a>.'
SUMMARY = "Template switching tags 5′-complete cDNA, cap trapping enriches genuine capped molecules, and paired-end sequencing connects a base-resolution transcription start site to downstream transcript sequence."
CAVEAT = ("The 40 six-base TSO barcodes are read from Table 25B.11.1 of the typeset "
          "Current Protocols article as reproduced in Philippe Batut's CSHL thesis "
          "(the publisher's Supplementary File 1 remains behind PMC's download gate). "
          "Which barcode each sample of the defining study used is printed in that "
          "paper's Supplemental Methods/Table S1, which was not read, so the inline "
          "barcode is drawn as a placeholder. The library PCR is shown as an endpoint "
          "rather than reconstructed duplex geometry.")

# Batut & Gingeras 2013, Synthetic oligonucleotides (PMC4372803 lines 188-228).
TSO_HANDLE = "TAGTCGAACTGAAGGTCTCCAGCA"          # rampage_TS_** 5' part = rampage_r1
RT_HANDLE = "TAGTCGAACGAAGGTCTCCGAACCGCTCTTCCGATCT"   # rampage_RT 5' of (N)15
PCR_F_TAIL_3P = "TAGTCGAACTGAAGG"                # rampage_F 3' part after P5
SBS8 = il.PE_READ2_PRIMER   # Illumina SBS8 Read 2 primer
R1_PRIMER = TSO_HANDLE
R2_PRIMER = SBS8
RAMPAGE_F = il.P5 + PCR_F_TAIL_3P
RAMPAGE_R = il.P7 + SBS8
if not TSO_HANDLE.startswith(PCR_F_TAIL_3P):
    raise ValueError("rampage_F 3' end must copy the TSO handle")


def common_suffix(a: str, b: str) -> int:
    n = 0
    while n < min(len(a), len(b)) and a[-1 - n] == b[-1 - n]:
        n += 1
    return n


# rampage_R carries the complete SBS8 site, but only its 3' end repeats the RT-primer
# handle: the first RT_HANDLE bases are replaced by the PCR primer's 5' tail.
RT_R_OVERLAP = common_suffix(RT_HANDLE, SBS8)

# Table 25B.11.1, typeset Current Protocols Suppl. 104 p. 25B.11.5 (thesis PDF p. 124),
# bold library-identification barcodes of rampage_TS_01 ... rampage_TS_40 in order.
TSO_BARCODES = (
    "AGGTAA", "CACTAC", "TTGGTC", "GTGTCA", "GCCGAA", "TGATCA", "CTGTAT", "CGACTG",
    "TTCCAG", "ACTCTT", "GGATAC", "TTAACG", "GAGTGC", "AAGGAC", "CGCGTT", "ATGCGT",
    "ATAAGC", "ATCTCC", "TAACTC", "ATAGAG", "AGCCTA", "TGTAGT", "AAACGG", "CCTACG",
    "ACTAGA", "CCCTCT", "GGTATA", "GATCCC", "CAATGT", "GCGTTG", "AACTGA", "CCAATA",
    "GCGACT", "GGGGAT", "TCTTCC", "GAACAT", "TCGAAG", "TGCTGC", "CTGCTA", "GACGTG")


def min_pairwise_hamming(codes=TSO_BARCODES) -> int:
    return min(sum(x != y for x, y in zip(a, b))
               for i, a in enumerate(codes) for b in codes[i + 1:])


BARCODE = feature("sample_inline", "sample_index", "whitelist",
                  whitelist="Batut & Gingeras 2013 Table 25B.11.1 (rampage_TS_01-40)")
RT_N15 = feature("rt_random_15mer", "umi", "degenerate",
                 note="RT-primer (N)15; used with alignment coordinates as a "
                      "pseudo-random barcode for duplicate collapsing")


def tso_oligo():
    return [seg("Read 1 handle", TSO_HANDLE, "r1"),
            seg("library barcode", "B"*6, "cbc", placeholder=True, feature=BARCODE),
            seg("rGrGrG (riboguanosines)", "GGG")]


def rt_primer_oligo():
    return [seg("RT handle", RT_HANDLE, "r2"),
            seg("random 15-mer", "N"*15, placeholder=True, feature=RT_N15)]


FINAL_LIBRARY = Construct([
    seg("P5", il.P5, "p5"), seg("custom Read 1 site", R1_PRIMER, "r1"),
    seg("inline library barcode", "B"*6, "cbc", placeholder=True, feature=BARCODE),
    seg("TSO-derived GGG", "GGG"),
    seg("5-prime-complete cDNA", "X"*42, placeholder=True),
    seg("random-primed 15 nt", "N"*15, placeholder=True, feature=RT_N15),
    seg("SBS8 Read 2 site reverse complement", revcomp(R2_PRIMER), "r2"),
    seg("P7 reverse complement", il.P7_RC, "p7")], name="RAMPAGE library")
SEQ_PRIMERS = (
    sp.custom("Read 1", "RAMPAGE custom Read 1 (rampage_r1)", R1_PRIMER,
              "Batut & Gingeras 2013, Synthetic oligonucleotides"),
    sp.custom("Read 2", "Illumina SBS8 Read 2", R2_PRIMER,
              "Batut & Gingeras 2013, Synthetic oligonucleotides"),
)
# Batut et al. 2013 Methods: HiSeq 2x101 for embryo, larva/pupa libraries;
# the adult library was GAIIx 2x76.
READ_LENGTHS = {"Read 1": 101, "Read 2": 101}
if sp.verify(FINAL_LIBRARY, SEQ_PRIMERS, required_roles=("Read 1", "Read 2")):
    raise ValueError("RAMPAGE custom sequencing primers do not land on the final library")
FINAL_CAPTION = ("Paired-end RAMPAGE library (embryo and larva/pupa HiSeq runs, 2×101; the "
                 "adult GAIIx run was 2×76). Read 1 starts with the 6-base library barcode "
                 "and the TSO-derived G stretch, then the transcription start; Read 2 starts "
                 "in the random-primed 15 nt.")
SEQUENCING_INTRO = ("Read 1 (custom rampage_r1) identifies the library barcode and the "
                    "transcription start boundary; Read 2 (standard SBS8) supplies the "
                    "RT-primer 15-mer and downstream transcript connectivity.")


def oligos():
    return [
        "rampage_TS_** (TSO): 5′-TAGTCGAACTGAAGGTCTCCAGCA<b>NNNNNN</b>rGrGrG-3′ "
        f"— six-base barcode, one of {len(TSO_BARCODES)} (Table 25B.11.1); last three residues riboguanosines",
        f"rampage_RT: 5′-{RT_HANDLE}(N)15-3′",
        f"rampage_F: 5′-{RAMPAGE_F}-3′",
        f"rampage_R: 5′-{RAMPAGE_R}-3′ — only its 3′ {RT_R_OVERLAP} bases repeat the rampage_RT handle",
        f"rampage_r1 (custom Read 1): 5′-{R1_PRIMER}-3′; Read 2: Illumina SBS8 5′-{R2_PRIMER}-3′",
    ]


def protocol_graph():
    graph = ProtocolGraph(TITLE, authority=SOURCE)
    rna = graph.input("rna", "DNase I-treated total RNA (5 µg; 2–3 µg when pooling)",
                      MaterialKind.MOLECULAR, primary=True)
    tex = graph.input("tex", "Terminator 5′-phosphate-dependent exonuclease, 30 °C 90 min",
                      MaterialKind.CATALYST)
    tex_rna = graph.product("tex_rna", "RNA depleted of 5′-monophosphate species (rRNA)",
                            MaterialKind.MOLECULAR)
    tex_digest = graph.product("tex_digest", "digested 5′-monophosphate RNA",
                               MaterialKind.WASTE)
    graph.react("tex_digest_step", "Degrade 5′-monophosphate RNA, then RNAClean XP cleanup",
                StepKind.DIGESTION,
                ((rna, InputRole.TEMPLATE), (tex, InputRole.CATALYST)),
                (tex_rna, tex_digest))
    graph.discard(tex_digest, "digestion products are not recovered by the RNAClean XP cleanup")
    rt_primer = graph.input("rt_primer", "rampage_RT: RT handle plus (N)15",
                            MaterialKind.OLIGO,
                            state=Construct(rt_primer_oligo(), name="rampage_RT"))
    tso = graph.input("tso", "rampage_TS_**: barcoded rGrGrG-ended TSO",
                      MaterialKind.OLIGO, state=Construct(tso_oligo(), name="rampage_TS_**"))
    rt = graph.input("reverse_transcriptase",
                     "SuperScript III with sorbitol/trehalose and betaine; 22 °C 1 min, 42 °C 30 min, 75 °C 15 min",
                     MaterialKind.CATALYST)
    first = graph.product("first_strand", "random-primed cDNA reaching RNA 5′ ends",
                          MaterialKind.MOLECULAR)
    graph.react("reverse_transcribe", "Random-primed RT in the same reaction as the TSO",
                StepKind.REVERSE_TRANSCRIPTION,
                ((tex_rna, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (first,))
    complete = graph.product("complete", "template-switched cDNA:RNA hybrids carrying the barcoded Read 1 arm",
                             MaterialKind.MOLECULAR)
    graph.react("template_switch", "Template switching at the RNA 5′ end installs the TSO",
                StepKind.TEMPLATE_SWITCHING,
                ((first, InputRole.TEMPLATE), (tso, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (complete,))
    pooling = graph.input("pooling", "RNAClean XP; optionally after qPCR (rampage_F/R) quantification and equimolar pooling of barcoded RT products",
                          MaterialKind.SELECTION)
    pooled = graph.product("pooled", "purified (optionally pooled) barcoded cDNA:RNA hybrids",
                           MaterialKind.MOLECULAR)
    rt_waste = graph.product("rt_cleanup_waste", "RNAClean XP supernatant and washes",
                             MaterialKind.WASTE)
    graph.react("pool", "Purify RT products; barcoded samples may be pooled here",
                StepKind.CLEANUP,
                ((complete, InputRole.TEMPLATE), (pooling, InputRole.SELECTION)),
                (pooled, rt_waste), changes_molecule=False)
    graph.discard(rt_waste, "bead supernatant and ethanol washes are discarded")
    periodate = graph.input("periodate", "NaIO4 in NaOAc pH 4.5, 45 min on ice; glycerol/Tris stop; RNAClean XP",
                            MaterialKind.CONDITION)
    oxidized = graph.product("oxidized", "hybrids with oxidized cap and 3′-terminal riboses",
                             MaterialKind.MOLECULAR)
    graph.react("oxidize", "Oxidize diol-bearing riboses (cap and RNA 3′ ends)",
                StepKind.CHEMICAL_CONVERSION,
                ((pooled, InputRole.TEMPLATE), (periodate, InputRole.CONDITION)), (oxidized,))
    biotin = graph.input("biotin_hydrazide", "biotin hydrazide (long arm), Na citrate pH 6.0, 14–15 h RT",
                         MaterialKind.CONDITION)
    biotinylated = graph.product("biotinylated", "hybrids biotinylated at the cap and RNA 3′ ends",
                                 MaterialKind.MOLECULAR)
    graph.react("biotinylate", "Biotinylate the oxidized riboses",
                StepKind.CHEMICAL_CONVERSION,
                ((oxidized, InputRole.TEMPLATE), (biotin, InputRole.CONDITION)), (biotinylated,))
    rnase = graph.input("rnase_i", "RNase I, 37 °C 30–60 min, 65 °C 5 min; RNAClean XP",
                        MaterialKind.CATALYST)
    trimmed = graph.product("trimmed", "hybrids with single-stranded RNA (and its 3′ biotin) removed",
                            MaterialKind.MOLECULAR)
    rnase_frag = graph.product("rnase_fragments", "digested single-stranded RNA",
                               MaterialKind.WASTE)
    graph.react("rnase_i_digest", "RNase I leaves biotin only on caps of 5′-complete hybrids",
                StepKind.DIGESTION,
                ((biotinylated, InputRole.TEMPLATE), (rnase, InputRole.CATALYST)),
                (trimmed, rnase_frag))
    graph.discard(rnase_frag, "digested RNA is not recovered by the RNAClean XP cleanup")
    captrap = graph.input("streptavidin", "tRNA-blocked MPG streptavidin beads; washes with buffers 1–4",
                          MaterialKind.SELECTION)
    captured = graph.product("captured", "bead-bound cap-trapped 5′-complete hybrids",
                             MaterialKind.MOLECULAR)
    unbound = graph.product("unbound", "unbiotinylated hybrids and wash supernatants",
                            MaterialKind.WASTE)
    graph.react("cap_trap", "Streptavidin pulldown (cap trapping)", StepKind.AFFINITY_CAPTURE,
                ((trimmed, InputRole.TEMPLATE), (captrap, InputRole.SELECTION)),
                (captured, unbound))
    graph.discard(unbound, "the supernatant and washes are discarded")
    naoh = graph.input("naoh", "50 mM NaOH 10 min RT, neutralised with Tris-HCl pH 7.0; AMPure XP",
                       MaterialKind.CONDITION)
    released = graph.product("released", "eluted 5′-complete cDNA", MaterialKind.MOLECULAR)
    beads = graph.product("spent_beads", "streptavidin beads", MaterialKind.WASTE)
    graph.react("release_cdna", "Alkaline elution of cDNA from the beads", StepKind.EXTRACTION,
                ((captured, InputRole.TEMPLATE), (naoh, InputRole.CONDITION)), (released, beads))
    graph.discard(beads, "beads are left on the magnet when the eluate is transferred")
    forward = graph.input("pcr_forward", "rampage_F (P5 + TSO-handle 5′ 15 nt)", MaterialKind.OLIGO,
                          state=Construct([seg("P5", il.P5, "p5"),
                                           seg("TSO handle 5′ part", PCR_F_TAIL_3P, "r1")],
                                          name="rampage_F"))
    reverse = graph.input("pcr_reverse", "rampage_R (P7 + complete SBS8 site)", MaterialKind.OLIGO,
                          state=Construct([seg("P7", il.P7, "p7"), seg("SBS8 site", SBS8, "r2")],
                                          name="rampage_R"))
    polymerase = graph.input("pcr_polymerase", "Ex Taq HS, 16 cycles", MaterialKind.CATALYST)
    amplified = graph.product("amplified", "amplified RAMPAGE library before size selection",
                              MaterialKind.MOLECULAR, state=FINAL_LIBRARY)
    graph.react("library_pcr", "PCR adds P5/P7; rampage_R completes the SBS8 Read 2 site",
                StepKind.PCR,
                ((released, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (amplified,), geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY,
                    reason=(f"rampage_R anneals through only the 3′ {RT_R_OVERLAP} bases of the "
                            "rampage_RT handle and replaces the rest with its own tail; this "
                            "compact model does not reconstruct the single-stranded cDNA as a PCR duplex.")))
    ampure = graph.input("size_selection", "AMPure XP 0.52× (remove large) then 0.70× (keep medium)",
                         MaterialKind.SELECTION)
    library = graph.product("library", FINAL_LIBRARY.name, MaterialKind.LIBRARY,
                            state=FINAL_LIBRARY)
    rejected = graph.product("rejected_sizes", "large-insert beads and short-fragment supernatant",
                             MaterialKind.WASTE)
    graph.react("size_select", "Double-sided AMPure XP size selection (~300–1,000 bp)",
                StepKind.SELECTION,
                ((amplified, InputRole.TEMPLATE), (ampure, InputRole.SELECTION)),
                (library, rejected), changes_molecule=False)
    graph.discard(rejected, "first-round beads and second-round supernatant are discarded")
    graph.sequence(library, SequencingPlan(library=FINAL_LIBRARY,
                                           primers=tuple(SEQ_PRIMERS),
                                           read_lengths=READ_LENGTHS))
    graph.declare_reviewed_scope((
        "rna", "tex", "tex_digest_step", "rt_primer", "tso", "reverse_transcriptase",
        "reverse_transcribe", "template_switch", "pooling", "pool", "periodate", "oxidize",
        "biotin_hydrazide", "biotinylate", "rnase_i", "rnase_i_digest", "streptavidin",
        "cap_trap", "naoh", "release_cdna", "pcr_forward", "pcr_reverse", "pcr_polymerase",
        "library_pcr", "size_selection", "size_select", "library"), SOURCE)
    return graph.require_publishable()
