"""GenIE targeted cDNA readout."""
from batch_ngs import seg, truseq_library
from chemdraw import Row

TITLE = "GenIE — endogenous variant editing with targeted cDNA sequencing"
NOTES = "01_genie.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1093/nar/gkaa960">Cooper et al., <i>Nucleic Acids Research</i> (2020)</a>.'
SUMMARY = "Arrayed CRISPR editing creates endogenous alleles; locus-specific RT-PCR measures each allele's effect on transcript abundance or splicing."
FINAL_LIBRARY, SEQ_PRIMERS = truseq_library(
    [seg("locus-specific cDNA amplicon", "X" * 56, placeholder=True)],
    "GenIE targeted cDNA library")
FINAL_CAPTION = "The insert is the locus-specific cDNA amplicon spanning the edited variant and the transcript feature being measured."
SEQUENCING_INTRO = "Standard paired-end Illumina primers read the indexed targeted amplicon."

def sections():
    return [
        ("Create an allelic series at the endogenous locus", [Row(chunks=[("arrayed guide + HDR donor → edited cell population", "cbc", False)])], "Each well contains one designed edit or control."),
        ("Reverse-transcribe cellular RNA", [Row(chunks=[("RNA → locus-containing cDNA", "r1", False)])], "The molecular readout is RNA, rather than genomic editing frequency alone."),
        ("Amplify the locus from cDNA", [Row(chunks=[("outer locus PCR → nested allele-spanning amplicon", "r2", False)])], "Primer pairs are locus-specific and are tabulated in the paper's supplement."),
        ("Install indexes and sequence", [Row(chunks=[("P5—index—read arm—target cDNA—read arm—index—P7", "cbc", False)])], "Indexed amplicons from the arrayed experiment are pooled for Illumina sequencing."),
    ]
