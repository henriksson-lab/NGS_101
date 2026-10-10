from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from chemdraw import Row
from sprite import bulk_barcode, final_library, ligation_round_rows

TITLE="SPRITE — split-pool tagging of spatial complexes"
NOTES="01_sprite.html"
SOURCE='Defining source: <a href="https://doi.org/10.1016/j.cell.2018.05.024">Quinodoz et al., <i>Cell</i> (2018)</a>.'
SUMMARY="Crosslinked multi-molecule complexes are immobilized and repeatedly split across 96 tag wells; molecules that remain in one complex acquire the same five-part barcode without proximity-ligating the genomic fragments to one another."
CAVEAT="Barcode-tag ligations are covalent, but genomic fragments in a cluster are not ligated to each other. Shared tag combinations encode co-membership of a crosslinked spatial complex."
TAGS,CLUSTER_ID=bulk_barcode()
FINAL_LIBRARY,SEQ_PRIMERS=final_library(TAGS,"bulk SPRITE library")
CLUSTER_ID.bind(FINAL_LIBRARY)
FINAL_CAPTION="Read 1 reports the DPM tag and genomic DNA; Read 2 reports the remaining split-pool tag chain."
SEQUENCING_INTRO="Paired-end sequencing reads genomic sequence from the DPM end and reconstructs the cluster identifier from DPM plus the four tag parts read from the opposite end."
READ_LENGTHS={"Read 1":100,"Read 2":75}

def sections(): return [
 ("Crosslink, isolate and fragment chromatin",[Row(chunks=[("crosslinked complex: [DNA fragment] ... [DNA fragment] ... [RNA/protein]",None,False)])],"DSG and formaldehyde preserve multiway complexes; sonication and DNase make sequenceable DNA fragments."),
 ("Immobilize intact complexes and repair DNA ends",[Row(chunks=[("NHS bead -- covalently immobilized spatial complex -- blunt, 5-prime-phosphorylated, dA-tailed DNA",None,False)])],"Low complex loading and detergent washes reduce accidental co-barcoding; repair creates a common ligation substrate."),
 ("Perform five split-pool tag ligations",ligation_round_rows(TAGS),"The published order is DPM, Odd, Even, Odd, Terminal. DPM/Terminal have 9-nt identifiers; Odd/Even tags have 17-nt identifiers."),
 ("Reverse crosslinks and PCR-complete the library",[Row(chunks=[("P5 -- DPM -- genomic DNA -- Terminal -- Odd -- Even -- Odd -- P7",None,False)])],"Terminal and DPM arms provide the two library-PCR entry sites; molecules are grouped after sequencing by the full tag combination."),
]
