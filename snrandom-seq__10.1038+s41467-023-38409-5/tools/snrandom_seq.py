"""Xu et al. 2023 snRandom-seq chemistry."""
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[2]; sys.path[:0]=[str(ROOT/"lib")]
from batch_ngs import truseq_library, seg
from chemdraw import MolecularState, Scene, Segment, Workflow, feature
from rna_prep import rna_template, single_strand

TITLE="snRandom-seq — FFPE single-nucleus total RNA"
NOTES="01_snrandom-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/s41467-023-38409-5">Xu et al., <i>Nature Communications</i> (2023)</a>.'
SUMMARY="In-nucleus random/oligo-dT RT and dA tailing are followed by extension from three-round barcoded hydrogel beads; Read 1 carries a 30-nt cell barcode and 8-nt UMI, while Read 2 reports RNA-derived sequence."
CAVEAT="The paper names a Vazyme universal Illumina kit but does not print that vendor adapter's bases. Those standard Illumina-compatible outer arms are therefore inferred; the published inline barcode, UMI and linker architecture is not inferred."

BLOCK_PRIMER="GAGAATGTGAGTGAAGATGTATGGTGANNNNNNN"
RT_HANDLE="GGAGTTGGAGTGAGTGGATGAGTGATGGAAGGAAT"
BEAD_LINKER_12="TGGT"
BEAD_LINKER_23="GAGA"
BEAD_ANCHOR_WRITTEN="/5Acryd/ATTATATATATUGTGAGTGATGGTTGAGGATGTGTGGAGATA"

INITIAL_NAME="fixed, permeabilized nucleus"
INITIAL_ROWS=rna_template("crosslink-reversed nuclear RNA",38)

cbc1=seg("cell barcode 1","B"*10,"cbc",placeholder=True,feature=feature("cell_barcode_1","cell_barcode","combinatorial",group="cell_barcode",part="round 1"))
cbc2=seg("cell barcode 2","C"*10,"cbc",placeholder=True,feature=feature("cell_barcode_2","cell_barcode","combinatorial",group="cell_barcode",part="round 2"))
cbc3=seg("cell barcode 3","D"*10,"cbc",placeholder=True,feature=feature("cell_barcode_3","cell_barcode","combinatorial",group="cell_barcode",part="round 3"))
umi=seg("UMI","U"*8,"umi",placeholder=True,feature=feature("umi","umi","random"))
# A long schematic insert keeps PE150 read-cycle derivation faithful to the documented
# run: Read 2 remains in transcript sequence instead of spuriously reaching identifiers.
insert=[cbc1,seg("TGGT linker",BEAD_LINKER_12),cbc2,seg("GAGA linker",BEAD_LINKER_23),cbc3,umi,seg("priming tract","T"*21),seg("RNA-derived insert","X"*160,placeholder=True)]
FINAL_LIBRARY,SEQ_PRIMERS=truseq_library(insert,"snRandom-seq library",inferred_adapters=True)
FINAL_CAPTION="Published 30-nt split-pool cell barcode, 8-nt UMI and transcript insert inside inferred vendor-standard Illumina outer adapters."
SEQUENCING_INTRO="The paper used paired-end 150-base NovaSeq sequencing and explicitly extracted the 30-base cell barcode plus 8-base UMI from Read 1; Read 2 generated the expression matrix."
READ_LENGTHS={"Read 1":150,"Read 2":150}

def workflow():
    w=Workflow(MolecularState(INITIAL_NAME,tuple(INITIAL_ROWS)))
    blocked_sc=Scene()
    blocked_sc.strand("blocked DNA",[Segment("genomic DNA", "X"*34,placeholder=True),Segment("block-primer extension",BLOCK_PRIMER)],label="blocked genomic DNA")
    blocked_sc.strand("nuclear RNA",[Segment("nuclear RNA","X"*38,placeholder=True)],label="RNA retained in nucleus")
    blocked_sc.labels("blocked DNA"); blocked_sc.labels("nuclear RNA")
    blocked=blocked_sc.rows()
    w.react("Block accessible genomic ssDNA",blocked,name="blocked nucleus",note="A 5′ handle plus N7 primer is extended on genomic DNA before RNA reverse transcription.")
    rt_sc=Scene()
    rt_sc.strand("blocked DNA",[Segment("genomic DNA","X"*34,placeholder=True),Segment("block-primer extension",BLOCK_PRIMER)],label="blocked genomic DNA")
    rt_sc.strand("RNA",[Segment("RNA template","X"*38,placeholder=True)],label="nuclear RNA")
    rt_sc.anneal("cDNA",[Segment("RNA-derived cDNA","x"*38,placeholder=True),Segment("RT handle",RT_HANDLE)],to="RNA",pair=("RNA-derived cDNA","RNA template"),label="random/oligo-dT cDNA",unpaired=("RT handle",))
    rt_sc.labels("blocked DNA"); rt_sc.labels("cDNA")
    rt=rt_sc.rows()
    w.react("Random + oligo-dT reverse transcription",rt,name="in-nucleus cDNA",note="Five random-primer and five oligo(dT)-primer variants undergo twelve annealing ramps from 8 to 42 °C.")
    tail_sc=Scene()
    tail_sc.strand("blocked DNA",[Segment("genomic DNA","X"*34,placeholder=True),Segment("block-primer extension",BLOCK_PRIMER)],label="blocked genomic DNA")
    tail_sc.strand("dA-tailed cDNA",[Segment("RT handle",RT_HANDLE),Segment("RNA-derived cDNA","X"*38,placeholder=True),Segment("3′ dA tail","A"*12)],label="dA-tailed cDNA")
    tail_sc.labels("blocked DNA"); tail_sc.labels("dA-tailed cDNA")
    tailed=tail_sc.rows()
    w.react("TdT dA tailing",tailed,name="dA-tailed cDNA",note="TdT and dATP add the homopolymer used for bead-primer extension.")
    bead_product=[cbc1,Segment("TGGT linker",BEAD_LINKER_12),cbc2,Segment("GAGA linker",BEAD_LINKER_23),cbc3,umi,Segment("poly(dT)","T"*21),Segment("RNA-derived cDNA","X"*38,placeholder=True)]
    barcode_sc=Scene()
    barcode_sc.strand("blocked DNA",[Segment("genomic DNA","X"*34,placeholder=True),Segment("block-primer extension",BLOCK_PRIMER)],label="blocked genomic DNA")
    barcode_sc.strand("barcoded product",bead_product,label="droplet-barcoded cDNA")
    barcode_sc.labels("blocked DNA"); barcode_sc.labels("barcoded product")
    barcoded=barcode_sc.rows()
    w.react("Droplet bead extension",barcoded,name="cell-barcoded cDNA",note="A hydrogel bead assembled by three split-pool ligations contributes three 10-base barcode parts and an 8-base UMI.")
    w.react("PCR-select barcoded cDNA",single_strand("PCR product",bead_product,"amplified barcoded cDNA"),name="purified barcoded amplicons",note="After emulsion breakage, primer PCR amplifies bead-linked cDNA; AMPure purification removes the blocked genomic material.")
    w.react("INFERRED — vendor end prep, ligation and PCR",Scene.duplex(list(FINAL_LIBRARY),label="sequencing library").rows(),name="final library",note="The named Vazyme universal kit fragments, end-repairs/A-tails, ligates an Illumina-compatible adapter and amplifies; its adapter bases are not printed in the paper.")
    return w
