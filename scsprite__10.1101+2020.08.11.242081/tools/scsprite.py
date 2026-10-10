from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from chemdraw import Row
from sprite import final_library, ligation_round_rows, single_cell_barcode

TITLE="scSPRITE — cell and spatial-complex barcoding"
NOTES="01_scsprite.html"
SOURCE='Defining preprint: <a href="https://doi.org/10.1101/2020.08.11.242081">Arrastia et al. (2020)</a>; detailed published method: <a href="https://doi.org/10.1038/s41587-021-00998-1"><i>Nature Biotechnology</i> (2022)</a>.'
SUMMARY="Three split-pool ligations barcode intact nuclei; after selecting individual nuclei and fragmenting chromatin, three more rounds extend the same molecules with a spatial-cluster barcode."
CAVEAT="The first three tag parts identify the cell and participate again in the full six-part spatial-cluster identifier. The data model records both overlapping composite meanings explicitly."
TAGS,IDENTIFIERS=single_cell_barcode()
FINAL_LIBRARY,SEQ_PRIMERS=final_library(TAGS,"scSPRITE library")
for x in IDENTIFIERS: x.bind(FINAL_LIBRARY)
FINAL_CAPTION="Read 1 reports genomic DNA plus DPM; Read 2 traverses the other five tags. The first three tags identify the nucleus and all six identify a spatial cluster."
SEQUENCING_INTRO="Read 1 is at least 120 bp and begins with genomic information plus DPM; Read 2 is at least 95 bp and reads Odd–Even–Odd–Even–Y-even."
READ_LENGTHS={"Read 1":120,"Read 2":95}

def sections(): return [
 ("Crosslink, permeabilize, repair and dA-tail nuclei",[Row(chunks=[("intact nucleus: repaired genomic fragments remain compartmentalized",None,False)])],"Keeping nuclei intact makes the first barcode rounds cell-specific."),
 ("Add three in-nucleus split-pool tags",ligation_round_rows(TAGS[:3]),"DPM, Odd and Even across three 96-well rounds form the cell barcode."),
 ("Select nuclei, lyse and immobilize spatial complexes",[Row(chunks=[("one cell barcode -> many bead-bound chromatin complexes", "cbc",False)])],"About 1,500 filtered nuclei were taken forward; sonication releases crosslinked complexes for bead coupling."),
 ("Extend each complex with three spatial tags",ligation_round_rows(TAGS[3:]),"Odd, Even and Y-even tags distinguish complexes while retaining the first three-part cell identity."),
 ("Reverse crosslinks and PCR-complete the library",[Row(chunks=[("P5 -- DPM -- genomic DNA -- Y-even -- Even -- Odd -- Even -- Odd -- P7",None,False)])],"Only reads with the complete six-tag order are usable for cell and cluster assignment."),
]
