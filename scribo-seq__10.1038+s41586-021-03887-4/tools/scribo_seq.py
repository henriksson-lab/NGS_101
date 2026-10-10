"""Van Insberghe et al. single-cell ribosome profiling."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
import illumina as il
import seqprimers as sp
from batch_ngs import seg
from chemdraw import Construct, MolecularState, Scene, Workflow, feature, revcomp

TITLE = "scRibo-seq — single-cell ribosome profiling"
NOTES = "01_scribo-seq.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1038/s41586-021-03887-4">Van Insberghe et al. (2021)</a>.'
SUMMARY = "MNase generates ribosome-protected fragments in individual wells. One-pot end repair, 3′ and 5′ RNA-adapter ligation, reverse transcription and indexed PCR preserve a 10-nt UMI and a 10-nt cell index."
CAVEAT = "The supplement names the ordering oligos; the schematic uses the reported read structure and common small-RNA landing sites without presenting the unverified full ordering strings."

SMALL_RNA_R1 = "GTTCAGAGTTCTACAGTCCGACGATC"
SMALL_RNA_I1 = "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC"
SMALL_RNA_I2 = revcomp(SMALL_RNA_R1)

FINAL_LIBRARY = Construct([
    seg("P5", il.P5, "p5"),
    seg("ten-base cell index", "C" * 10, "cbc", placeholder=True,
        feature=feature("cell_index", "cell_barcode", "whitelist", whitelist="published well-specific forward primers")),
    seg("small-RNA Read 1 site", SMALL_RNA_R1, "r1"),
    seg("ten-base UMI", "U" * 10, "umi", placeholder=True,
        feature=feature("umi", "umi", "random")),
    seg("ribosome-protected fragment", "R" * 35, placeholder=True, length_bp=35),
    seg("3-prime adapter / Index 1 arm", SMALL_RNA_I1, "r2"),
    seg("six-base plate index reverse complement", "I" * 6, "cbc", placeholder=True,
        feature=feature("plate_i7", "sample_index", "whitelist", whitelist="published RPI primer set")),
    seg("P7 reverse complement", il.P7_RC, "p7"),
], name="scRibo-seq library")
SEQ_PRIMERS = (
    sp.custom("Read 1", "scRibo small-RNA Read 1", SMALL_RNA_R1, "Van Insberghe et al. Supplementary Table 1"),
    sp.custom("Index 1 (i7)", "scRibo plate-index primer", SMALL_RNA_I1, "Van Insberghe et al. Supplementary Table 1"),
    sp.custom("Index 2 (i5)", "scRibo cell-index primer", SMALL_RNA_I2, "Van Insberghe et al. Supplementary Table 1"),
)
if problems := sp.verify(FINAL_LIBRARY, SEQ_PRIMERS, required_roles=tuple(p.role for p in SEQ_PRIMERS)):
    raise ValueError("scRibo-seq final library: " + "; ".join(problems))
READ_LENGTHS = {"Read 1": 75, "Index 1 (i7)": 6, "Index 2 (i5)": 10}
FINAL_CAPTION = "Read 1 begins with the 10-nt UMI and then enters the 30–45-nt footprint. The six-cycle i7 read is the plate index; the ten-cycle i5 read is the well/cell index."
SEQUENCING_INTRO = "The defining paper specifies 75 cycles Read 1, six i7 cycles and ten i5 cycles on NextSeq 500. Identifier cycle spans are derived from the annotated final construct."

_initial = Scene()
_initial.strand("mRNA", [seg("translated mRNA", "R" * 48, placeholder=True)], label="cycloheximide-stalled translating mRNA")
_initial.mark("mRNA", "translated mRNA", "ribosome shields one footprint")
INITIAL_ROWS = tuple(_initial.rows())
INITIAL_NAME = "Single-cell translating mRNA"


def workflow():
    footprint = Scene(); footprint.strand("RPF", [seg("ribosome-protected RNA", "R" * 35, placeholder=True, length_bp=35)], label="MNase-resistant RPF")
    footprint.mark("RPF", "ribosome-protected RNA", "30–45 nt retained after MNase")
    ligated = Scene(); ligated.strand("RNA", [
        seg("ten-base UMI", "U" * 10, "umi", placeholder=True,
            feature=feature("umi", "umi", "random")),
        seg("5-prime adapter fixed bases", "X" * 12, "r1", placeholder=True),
        seg("RPF", "R" * 35, placeholder=True), seg("3-prime pre-adenylated adapter", "X" * 18, "r2", placeholder=True)], label="doubly ligated footprint RNA")
    ligated.junction("RNA", "5-prime adapter fixed bases", "RPF", "RNA ligation")
    ligated.junction("RNA", "RPF", "3-prime pre-adenylated adapter", "RNA ligation")
    ligated.labels("RNA")
    repaired = Scene(); repaired.strand("RPF", [seg("end-repaired RPF", "R" * 35, placeholder=True)], label="5′-phosphate / 3′-hydroxyl RPF")
    repaired.mark("RPF", "end-repaired RPF", "T4 PNK prepares both ends")
    indexed = Scene.duplex(list(FINAL_LIBRARY), label="indexed single-cell RPF library")
    # Size selection removes molecules, not sequence features from a retained molecule.
    # Keep the complete indexed construct so the final Workflow state still carries the
    # cell index, UMI and plate index installed in the preceding reaction.
    selected = Scene.duplex(list(FINAL_LIBRARY), label="175–185-bp PAGE-selected library")
    selected.mark("top", "ribosome-protected fragment", "30–40-nt footprint-size selection")
    w = Workflow(MolecularState(INITIAL_NAME, INITIAL_ROWS))
    w.react("Digest exposed RNA with MNase", footprint.rows(), name="Ribosome-protected fragment",
            note="The ribosome shields a 30–45-nt footprint in each single-cell well.")
    w.react("Repair footprint ends", repaired.rows(), name="Ligatable RPF",
            note="T4 PNK prepares both ends for the two RNA ligations.")
    w.react("Ligate 3′ and 5′ RNA adapters", ligated.rows(), name="Adapter-ligated RPF",
            note="T4 Rnl2(tr) KQ adds the pre-adenylated 3′ adapter; T4 Rnl1 adds the UMI-bearing 5′ adapter.")
    w.react("Reverse-transcribe and index each well", indexed.rows(), name="Indexed cDNA library",
            note="Maxima H-minus RT copies RNA; well forward and RPI reverse primers install cell and plate indices.")
    w.react("Pool and footprint-size select", selected.rows(), name="Footprint-selected pool",
            note="PAGE retains the 175–185-bp library band, corresponding to a 30–40-nt insert.")
    return w
