"""M3-seq molecular model."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from chemdraw import Row
from microbial_singlecell import bacterial_rna_library,droplet_barcode_rows,indexed_random_rt_scene

TITLE="M3-seq — multiplexed microbial single-cell RNA-seq"
NOTES="01_m3-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/s41564-023-01462-3">McNulty et al., <i>Nature Microbiology</i> (2023)</a>.'
SUMMARY="Random in-situ RT supplies a plate barcode and UMI; 10x ATAC droplets ligate a second barcode, followed by Tn5/T7 amplification and post hoc RNase-H rRNA depletion."
CAVEAT="The page preserves the important RNA→cDNA→T7 RNA→cDNA transitions. Proprietary droplet-barcode bases remain role tokens."
FINAL_LIBRARY,SEQ_PRIMERS=bacterial_rna_library("m3seq","M3-seq RNA library")
FINAL_CAPTION="The final cDNA library retains the round-one barcode, round-two barcode and UMI through T7 amplification and depletion."
SEQUENCING_INTRO="Read cycles recover the combinatorial cell identity, UMI and bacterial transcript insert."
READ_LENGTHS={"Read 1":28,"Read 2":90}

def sections(): return [
 ("Index fixed bacteria by random-primed in-situ RT",indexed_random_rt_scene(method="m3seq").rows(),"The RT primer carries BC1 and a UMI."),
 ("Ligate the second index in 10x ATAC droplets",droplet_barcode_rows("m3seq"),"BC2 is ligated to cell-associated cDNA; BC1+BC2 defines a cell."),
 ("Lyse cells and synthesize the second cDNA strand",[Row(chunks=[("single-stranded indexed cDNA → Klenow/random-primer dsDNA",None,False)])],"Bulk second-strand synthesis makes the library amplifiable."),
 ("Tagment, PCR-add a T7 promoter and transcribe",[Row(chunks=[("dsDNA → Tn5 fragment → T7-promoter PCR → amplified RNA", "me",False)])],"T7 transcription linearly amplifies indexed library molecules."),
 ("Deplete rRNA library molecules",[Row(chunks=[("RNA library + rRNA DNA probes → RNase H cleavage",None,False)])],"Post hoc depletion acts on amplified RNA rather than intact cells."),
 ("Reverse-transcribe depleted RNA and add sequencing adapters",[Row(chunks=[("depleted RNA → cDNA → final adapter PCR",None,False)])],"A final RT and PCR return the selected molecules to an Illumina-compatible cDNA library."),
]
