"""Original Kaya-Okur et al. CUT&Tag molecular model."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from batch_ngs import duplex_fragment, nextera_library, seg, targeting_rows
from chemdraw import Row, Scene

TITLE = "CUT&Tag — antibody-tethered Tn5 chromatin profiling"
NOTES = "01_cut-tag.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1038/s41467-019-09982-5">Kaya-Okur et al., <i>Nature Communications</i> (2019)</a>.'
SUMMARY = "An antibody recruits Protein A–Tn5 to a chromatin epitope; magnesium activates adapter-loaded Tn5 so cleavage and adapter installation occur at the targeted site in one reaction."
CAVEAT = "The original paper specifies Tn5MEDS-A/B transposomes and indexed PCR. The page uses the repository’s canonical Nextera mosaic-end and primer sequences; protein geometry is schematic rather than a claimed covalent nucleic-acid structure."

FINAL_LIBRARY, SEQ_PRIMERS = nextera_library(
    [seg("target-proximal chromatin", "X" * 36, placeholder=True)], "CUT&Tag library")
FINAL_CAPTION = "A productive heterologous S5/S7 tagmentation product after gap fill and indexed PCR. Tn5 installs the mosaic ends at the antibody-selected chromatin fragment boundaries."
SEQUENCING_INTRO = "Paired reads enter the targeted chromatin from the two Tn5 insertion sites; index reads identify the PCR-added sample indices."

def sections():
    frag = duplex_fragment("epitope-proximal DNA")
    return [
        ("Bind antibody and adapter-loaded pA–Tn5",
         targeting_rows("Tn5 [S5-ME + S7-ME]", "targeted transposome retained after washes"),
         "Primary antibody supplies specificity; Protein A couples the antibody to a Tn5 dimer preloaded with two different mosaic-end adapters."),
        ("Activate targeted tagmentation",
         [Row(chunks=[("S5—ME ** target-proximal DNA ** ME—S7", "me", False)]),
          Row(chunks=[("          Mg²⁺-activated cleavage + adapter transfer", None, False)])],
         "** marks the two coupled cleavage/adapter-transfer junctions. Free pA–Tn5 is washed away before magnesium activation, limiting background."),
        ("Gap-fill and amplify only productive S5/S7 fragments", frag.rows(),
         "Tagmentation leaves the standard 9-bp staggered gap. Extension and indexed PCR complete the two-ended Illumina library."),
    ]
