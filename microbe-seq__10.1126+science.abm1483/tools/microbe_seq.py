"""Microbe-seq molecular model."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from batch_ngs import nextera_library,seg
from chemdraw import Row,feature
from microbial_singlecell import indexed_mda_rows

TITLE="Microbe-seq — droplet single-microbe genome sequencing"
NOTES="01_microbe-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1126/science.abm1483">Zheng et al., <i>Science</i> (2022)</a>.'
SUMMARY="Individual microbes are lysed and MDA-amplified in droplets, tagmented by reagent fusion, then fused with barcode-bead droplets for indexed PCR."
CAVEAT="The paper establishes droplet operations and library order but not every commercial bead base; those regions are structural placeholders."
CELL=feature("microbeseq_droplet","cell_barcode","whitelist",whitelist="Microbe-seq barcode bead set")
FINAL_LIBRARY,SEQ_PRIMERS=nextera_library([seg("microbe barcode","B"*16,"cbc",placeholder=True,feature=CELL),seg("single-microbe genomic DNA","X"*46,placeholder=True)],"Microbe-seq WGS library")
FINAL_CAPTION="A bead barcode identifies all tagmented MDA fragments derived from one microbial droplet."
SEQUENCING_INTRO="Paired genomic reads are associated with the bead-derived microbe barcode and sample indexes."
READ_LENGTHS={"Read 1":150,"Read 2":150}

def sections(): return [
 ("Encapsulate and lyse individual microbes",[Row(chunks=[("one microbe + lysis droplet",None,False)])],"Lysis remains physically partitioned by cell."),
 ("Fuse with amplification droplets for MDA",indexed_mda_rows()[:1],"Phi29 multiple-displacement amplification creates enough genomic DNA while retaining droplet identity."),
 ("Fuse with Tn5 droplets for tagmentation",indexed_mda_rows()[1:2],"Nextera mosaic ends are transferred to amplified genomic fragments inside droplets."),
 ("Fuse with barcode-bead/PCR droplets",indexed_mda_rows()[2:],"Bead primers install one droplet barcode on the tagmented fragments."),
 ("Break emulsion and complete sequencing adapters",[Row(chunks=[("barcoded genomic amplicons → sample-index PCR → Illumina library",None,False)])],"Bulk PCR completes the standard sequencing arms."),
]
