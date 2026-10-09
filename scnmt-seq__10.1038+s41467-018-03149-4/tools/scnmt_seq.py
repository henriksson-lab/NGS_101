"""Clark et al. 2018 scNMT-seq chemistry and its two library branches."""
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT/"lib"), str(ROOT/"scbs-seq__10.1038+nmeth.3035"/"tools")]

from base_conversion import bisulfite_path
from batch_ngs import nextera_library, seg
from chemdraw import Row
import scbs

TITLE = "scNMT-seq — RNA, DNA methylation and accessibility"
NOTES = "01_scnmt-seq.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1038/s41467-018-03149-4">Clark et al., <i>Nature Communications</i> (2018)</a>.'
SUMMARY = "M.CviPI marks accessible GpC sites, then bead capture separates RNA for Smart-seq2/Nextera XT from DNA for scBS-seq."


def rna_library():
    return nextera_library([seg("Smart-seq2 cDNA insert", "X"*36, placeholder=True)],
                           "scNMT-seq RNA library")


def dna_library():
    return scbs.final_library(), tuple(scbs.SEQ_PRIMERS)


RNA_LIBRARY, RNA_PRIMERS = rna_library()
DNA_LIBRARY, DNA_PRIMERS = dna_library()


def sections():
    return [
        ("Mark accessible chromatin", [
            Row(chunks=[("nucleosome-free GpC + M.CviPI + SAM → Gp5mC", "w1", False)]),
            Row(chunks=[("nucleosome-protected GpC → unmarked", None, False)]),
        ], "GpC methylation records accessibility before the cell is lysed."),
        ("Separate RNA and DNA", [
            Row(chunks=[("poly(A) RNA → oligo-dT magnetic beads → Smart-seq2 cDNA", "r1", False)]),
            Row(chunks=[("supernatant genomic DNA → bisulfite conversion → scBS-seq", "r2", False)]),
        ], "Repeated bead washes are transferred to the DNA plate to maximize recovery."),
        ("Interpret the DNA branch", [
            Row(chunks=[(bisulfite_path(protected=False).text() + "  (closed unmarked GpC or unmethylated CpG)", None, False)]),
            Row(chunks=[(bisulfite_path(protected=True).text() + "  (accessible M.CviPI-marked GpC or endogenous 5mCpG)", "w1", False)]),
        ], "Sequence context separates exogenous GpC accessibility marks from endogenous CpG methylation."),
    ]
