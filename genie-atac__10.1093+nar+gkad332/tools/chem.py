"""GenIE-ATAC single-adapter targeted accessibility readout."""
from batch_ngs import seg, truseq_library
from chemdraw import Row

TITLE = "GenIE-ATAC — targeted allelic chromatin accessibility"
NOTES = "01_genie-atac.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1093/nar/gkad332">Cooper et al., <i>Nucleic Acids Research</i> (2023)</a>.'
SUMMARY = "Tn5 loaded with one mosaic-end adapter tags open chromatin; locus-specific linear amplification and biotin capture make the edited site sequenceable."
MEDS_A = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
MEDS_REV = "CTGTCTCTTATACACATCTT"
FINAL_LIBRARY, SEQ_PRIMERS = truseq_library(
    [seg("nested locus-to-MEDS-A amplicon", "X" * 62, placeholder=True)],
    "GenIE-ATAC targeted accessibility library")
FINAL_CAPTION = "Each insert joins a locus-specific nested-primer end to a nearby single-adapter Tn5 insertion."
SEQUENCING_INTRO = "The final barcoding PCR installs standard paired-end MiSeq ends around the targeted amplicon."

def oligos():
    return ["<b>MEDS-A</b> — 5′-TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG-3′",
            "<b>MEDS-REV</b> — 5′-/5Phos/CTGTCTCTTATACACATCT/intT/-3′"]

def sections():
    return [
        ("Tagment open chromatin with single-adapter Tn5", [Row(chunks=[("MEDS-A ** accessible genomic fragment ** MEDS-A", "me", False)])], "Using only MEDS-A leaves one universal transposon-primer site."),
        ("Linearly copy the target from a biotinylated primer", [Row(chunks=[("biotin—locus primer → target SNP → MEDS-A", "r1", False)])], "The locus-specific product terminates at a nearby tagmentation site."),
        ("Capture biotinylated products", [Row(chunks=[("biotinylated target strand ● streptavidin bead", "cbc", False)])], "Pulldown removes the overwhelming background of unrelated MEDS-A fragments."),
        ("Semi-nested PCR and index PCR", [Row(chunks=[("nested locus primer → SNP → MEDS-A primer → indexed MiSeq library", "r2", False)])], "A second locus primer and the common MEDS-A primer define the sequenced insert."),
    ]
