"""BacDrop molecular model."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from chemdraw import Row
from microbial_singlecell import bacterial_rna_library, droplet_barcode_rows, indexed_random_rt_scene

TITLE="BacDrop — droplet bacterial single-cell RNA-seq"
NOTES="01_bacdrop.html"
SOURCE='Defining source: <a href="https://doi.org/10.1016/j.cell.2023.01.002">Ma et al., <i>Cell</i> (2023)</a>.'
SUMMARY="Fixed bacteria receive a plate-specific RT index, are pooled at several cells per 10x droplet, and acquire a second bead barcode after TdT poly(A)-tailing of cDNA."
CAVEAT="The two-index cell identity is published, while commercial 10x bead sequence is represented by its molecular role."
FINAL_LIBRARY,SEQ_PRIMERS=bacterial_rna_library("bacdrop","BacDrop RNA library")
FINAL_CAPTION="The plate barcode, droplet barcode and UMI remain physically linked to each bacterial cDNA insert."
SEQUENCING_INTRO="Read cycles recover the two-part cell identity and UMI before transcript sequence."
READ_LENGTHS={"Read 1":28,"Read 2":90}

def sections(): return [
 ("Fix, permeabilize and deplete abundant rRNA in cells",[Row(chunks=[("bacterium: RNA — rRNA probes/RNase H — DNase I cleanup",None,False)])],"Universal rRNA probes and RNase H reduce ribosomal molecules before indexing."),
 ("Install the first barcode by random-primed in-situ RT",indexed_random_rt_scene(method="bacdrop").rows(),"Each plate well contributes one RT barcode and each copied molecule receives a UMI."),
 ("Pool cells and partition several bacteria per 10x droplet",[Row(chunks=[("many plate indexes + one droplet bead barcode", "cbc",False)])],"Combinatorial plate and droplet indexes permit super-loading."),
 ("TdT-tail cDNA and transfer the droplet barcode",droplet_barcode_rows("bacdrop"),"A poly(A) tail substitutes for inefficient bacterial template switching; bead-primed second-strand synthesis transfers the second index."),
 ("Enrich cDNA and complete the Illumina library",[Row(chunks=[("two-index cDNA → enrichment PCR → sequencing library",None,False)])],"The cell-associated cDNA is released and amplified after emulsion break."),
]
