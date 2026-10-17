"""SPATAC-seq protocol page model."""
from pathlib import Path
import sys

sys.path[:0] = [str(Path(__file__).resolve().parent),
                str(Path(__file__).resolve().parents[2] / "lib")]

import spatac_model as spa
from protocol_graph import (InputRole, MaterialKind, SourceClaim, SourceDiscrepancy,
                            StepKind)


TITLE = "SPATAC-seq — split-pool ligation-based single-cell ATAC-seq"
NOTES = "01_spatac_seq.html"
SOURCE = ('<a href="https://doi.org/10.1038/s41556-024-01449-0">Sun et al., '
          '<i>Nature Cell Biology</i> (2024)</a>, Figure 1 and Supplementary Table 1; '
          '<a href="https://github.com/Lan-lab/SPATAC-seq">authors\' processing code</a>; '
          'published Methods (nuclei isolation, transposome assembly, transposition, '
          'hybridization–ligation, sub-library PCR); GEO GSE215901 and GSE243256 sample records.')
SUMMARY = ("A two-ended barcode is installed during custom-Tn5 transposition, then "
           "successive split-pool ligations append a second barcode to the Read-1 end "
           "and a third barcode to the Read-2 end. Four 8-nt pieces identify one nucleus.")
CAVEAT = ("SPATAC-seq is a split-pool single-cell method, not Spatial-ATAC-seq. The "
          "representative sequences below are one member of each published oligo panel; "
          "the linked research note identifies the complete tables and whitelist.")
FINAL_LIBRARY, SEQ_PRIMERS = spa.library("accessible genomic DNA")
# Acquired cycles: GEO GSE243256 series overall design ("SPATAC-seq ... libraries were
# sequenced with paired-end 150-bp reads on Illumina NovaSeq 6000 platform or
# DNBSEQ-T7") and the article's Reporting Summary (PE150).  The authors' processing
# spans (Read 1 to 115, Read 2 to 100) are described in SEQUENCING_INTRO.
READ_LENGTHS = {
    "Read 1": 150,
    "Index 1 (i7)": 8,
    "Index 2 (i5)": 8,
    "Read 2": 150,
}
# Spans the authors' processing code reads (Lan-lab/SPATAC-seq); not acquired cycles.
PROCESSING_SPANS = {"Read 1": 115, "Index 1 (i7)": 8, "Index 2 (i5)": 8, "Read 2": 100}
SHOW_READ_CONTENT = True
FINAL_CAPTION = ("The hybrid final library uses a TruSeq Read 1 site on the left and a "
                 "Nextera S7/Read 2 site on the right. ** marks the two split-pool "
                 "ligation bonds. Barcode pieces are structured parts of one 32-nt "
                 "cell identifier; i5 and i7 identify the PCR sublibrary.")
SEQUENCING_INTRO = ("The Reporting Summary states that all libraries were sequenced "
                    "paired-end 150 bp on Illumina NovaSeq 6000 or HiSeq X Ten (or on MGI "
                    "T7/MGISEQ-2000), and GEO GSE243256 likewise states paired-end 150-bp "
                    "reads; the cycle counts shown are therefore 150/8/8/150. The authors' "
                    "processing code uses only Read 1 up to cycle 115 and Read 2 up to "
                    "cycle 100; barcode positions do not change with read length. "
                    "The authors' demonstrated preprocessing retains Read 1 cycles "
                    "66–115 and Read 2 cycles 51–100 as the 50-nt genomic mates. "
                    "The published Methods and the GEO records GSE215901 and GSE243256 "
                    "instead describe the genomic mates as raw Read 1 positions 66–121 and "
                    "Read 2 positions 51–101: "
                    "the start positions agree, the trimming end points differ, and the "
                    "page follows the code. Before "
                    "those positions, Read 1 reports round 2 and the forward half of "
                    "round 1; Read 2 reports round 3 and the reverse half of round 1. "
                    "The two eight-cycle index reads report the PCR sublibrary.")


def oligos():
    return spa.oligos()


def workflow():
    return spa.workflow("accessible genomic DNA")


def final_rows():
    return spa.final_rows(FINAL_LIBRARY)


METHODS_NUCLEI = ("Published Methods, 'Nuclei isolation and fixation from zebrafish embryos' "
                  "and '... from cell lines' (Nature Cell Biology 26:1187-1199)")
GEO_GSE243256 = "GEO GSE243256 sample records (Sample_extract_protocol_ch1, zebrafish embryo nuclei)"
GEO_GSE215901 = "GEO GSE215901 sample records (Sample_extract_protocol_ch1, Hepa/K562 nuclei)"


def prepare_nuclei(graph):
    """Nuclei isolation and fixation from the published Methods.

    The Methods (and GEO GSE243256) describe the zebrafish-embryo route; the cell-line
    route (and GEO GSE215901) lyses Hepa/K562 with Omni-ATAC buffer and then follows the
    embryo fixation.  The formaldehyde volume follows the published Methods; GEO prints a
    different volume (see ``protocol_graph``).
    """
    both = METHODS_NUCLEI + "; " + GEO_GSE243256
    sample = graph.input("sample", "frozen zebrafish embryos (species-mixing runs: Hepa and K562 cells)",
                         MaterialKind.BIOLOGICAL, primary=True, source=both)
    lysis = graph.input("homogenization_buffer",
                        "Homogenization Buffer: 10 mM Tris-HCl pH 7.8, 5 mM CaCl2, 3 mM Mg(Ac)2, "
                        "0.1% NP-40, 0.1 mM EDTA, 320 mM sucrose, 0.015 mM PMSF, 0.15 mM "
                        "β-mercaptoethanol (cell lines: Omni-ATAC lysis buffer, 0.1% NP-40, "
                        "0.1% Tween-20, 0.01% digitonin, 3 min on ice)",
                        MaterialKind.CONDITION, source=both + "; " + GEO_GSE215901)
    nuclei = graph.product("isolated_nuclei", "isolated nuclei washed in NWB",
                           MaterialKind.BIOLOGICAL, source=both)
    graph.react("isolate_nuclei",
                "Homogenize in a 2 ml Dounce (\u2018pipetting ten times each with pestle A "
                "and pestle B\u2019), filter through "
                "40 µm and wash in NWB (10 mM Tris-HCl pH 7.4, 10 mM NaCl, 3 mM MgCl2, "
                "1% BSA, 0.1% Tween-20); 350 g, 5 min, 4 °C",
                StepKind.EXTRACTION,
                ((sample, InputRole.TEMPLATE), (lysis, InputRole.CONDITION)), (nuclei,),
                source=both)
    formaldehyde = graph.input("formaldehyde",
                               "138 µL 37% formaldehyde in 10 mL NWB (final 0.5%)",
                               MaterialKind.CONDITION, source=METHODS_NUCLEI)
    fixed = graph.product("fixed_nuclei", "formaldehyde-fixed nuclei",
                          MaterialKind.BIOLOGICAL, source=both)
    graph.react("fix_nuclei", "Fix nuclei for exactly 10 min on ice", StepKind.CROSS_LINKING,
                ((nuclei, InputRole.TEMPLATE), (formaldehyde, InputRole.CONDITION)), (fixed,),
                source=both)
    glycine = graph.input("glycine", "500 µL 2.5 M glycine", MaterialKind.CONDITION,
                          source=both)
    ready = graph.product("nuclei", "fixed nuclei at 5,000 per µL in NWB",
                          MaterialKind.BIOLOGICAL, source=both)
    graph.react("quench_fixation",
                "Quench with glycine (5 min, room temperature), wash twice in NWB and "
                "dilute to 5,000 nuclei per µL",
                StepKind.QUENCHING,
                ((fixed, InputRole.TEMPLATE), (glycine, InputRole.CONDITION)), (ready,),
                source=both)
    return ready, ("sample", "homogenization_buffer", "isolated_nuclei", "isolate_nuclei",
                   "formaldehyde", "fixed_nuclei", "fix_nuclei", "glycine", "quench_fixation")


def protocol_graph():
    graph = spa.protocol_graph(
        "accessible genomic DNA",
        authority="Defining source: https://doi.org/10.1038/s41556-024-01449-0",
        prepare_nuclei=prepare_nuclei, read_lengths=READ_LENGTHS)
    graph.annotate(SourceDiscrepancy(
        "fix_nuclei", "Formaldehyde volume differs between the Methods and GEO",
        (SourceClaim("Published Methods, 'Nuclei isolation and fixation from zebrafish embryos' "
                     "(article PDF text twin line 3381)",
                     "138 µl 37% formaldehyde (final concentration, 0.5%) in 10 ml NWB"),
         SourceClaim("GEO GSE243256 Sample_extract_protocol_ch1",
                     "140 µl of 37% formaldehyde in 10 ml NWB")),
        "Use the published Methods value (138 µl; 0.5% final).",
        "Both values are printed in the original records; the published article takes "
        "precedence over the deposited sample text.",
        "Whether the GEO value is rounding or a different batch is not stated; both give "
        "about 0.5% formaldehyde, so the difference does not change the modelled chemistry."))
    graph.annotate(SourceDiscrepancy(
        "library_a", "Round-3 barcode length in the drawn final library",
        (SourceClaim("Published article, Extended Data Fig. 1a, row '4th BC-indexed PCR' "
                     "(rendered PDF page 20)",
                     "the Read 2 strand of the final library prints the round-3 barcode as "
                     "[NNNNNNN] (seven N)"),
         SourceClaim("Published article, Extended Data Fig. 1a, row '3rd BC-ligation', and "
                     "Supplementary Table 1 (Tn5_R_L1–L48; Tn5_R_LA-1–48)",
                     "the round-3 barcode is [NNNNNNNN] (eight N); every Tn5_R_L oligo ends "
                     "in an 8-nt barcode and each splint carries its 8-nt reverse complement")),
        "Use an 8-nt round-3 barcode.",
        "Supplementary Table 1 prints every round-3 oligo and splint with 8 barcode bases, "
        "the figure's own round-3 ligation row draws eight, and the figure's 50-nt Read 2 "
        "prefix bracket (barcode + 15-nt leader + 8-nt round-1 half + 19-nt mosaic end) "
        "adds up only with eight.",
        "The seven-N row looks like a drawing slip; the article does not correct it. It "
        "does not change any modelled sequence."))
    return graph
