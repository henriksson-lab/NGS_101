"""SAMOSA-Tag molecular model."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from batch_ngs import joined_scene, seg
from chemdraw import Row
from longread_footprinting import hairpin_tagmented_smrtbell, marked_fiber, pacbio_entry, smrtbell_rows

TITLE = "SAMOSA-Tag — hairpin tagmentation of marked chromatin"
NOTES = "01_samosa-tag.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1038/s41588-024-01748-0">Dennis et al., <i>Nature Genetics</i> (2024)</a>.'
SUMMARY = "EcoGII first records chromatin accessibility as m6dA; hairpin-loaded Tn5 then fragments the marked fibers while installing PacBio-compatible closed ends."
CAVEAT = "Hairpin-Tn5 adapters are published as protocol oligos but are drawn by molecular role here. The final molecule is shown only after gap repair has made a closed SMRTbell."
SMRTBELL = hairpin_tagmented_smrtbell()
FINAL_LIBRARY = None
SEQ_PRIMERS = ()
SEQUENCING_ENDING = pacbio_entry(SMRTBELL)

def sections():
    return [
        ("Mark accessible adenines in nuclei", marked_fiber("EcoGII"),
         "EcoGII and SAM deposit m6dA on exposed DNA while protein-bound intervals remain footprints."),
        ("Tagment with hairpin-loaded Tn5", joined_scene([seg("left hairpin-ME","X"*20,"me",placeholder=True),seg("marked chromatin fragment","X"*36,placeholder=True),seg("right ME-hairpin","X"*20,"me",placeholder=True)],(("left hairpin-ME","marked chromatin fragment","Tn5 transfer"),("marked chromatin fragment","right ME-hairpin","Tn5 transfer")),label="hairpin-tagmented chromatin").rows(),
         "Tn5 couples fragmentation to transfer of a hairpin-bearing mosaic-end adapter at both boundaries."),
        ("Gap-repair the tagmented product", smrtbell_rows(SMRTBELL),
         "Polymerase fills the transposition gaps, leaving one covalently closed PacBio template."),
    ]
