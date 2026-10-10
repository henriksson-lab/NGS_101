"""Slide-DNA-seq molecular model."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
import illumina as il
import nextera as nx
import seqprimers as sp
from chemdraw import Construct,Row,Segment,feature,oligo,revcomp

READ1="GCTTTGCTAACGGTCGAGAGATGTGTATAAGAGACAG"
READ2="CGGATGTTGCACCAGCAGATGTGTATAAGAGACAG"
T5="GAGCTTTGCTAACGGTCGAG"+nx.ME
T7="CTTACGGATGTTGCACCAGC"+nx.ME
FINAL_LIBRARY=Construct([
 Segment("P5",il.P5,"p5"),Segment("i5","J"*8,"cbc",placeholder=True,feature=feature("slide_sample_i5","sample_index","unknown")),Segment("TruSeq Read 1",il.TRUSEQ_READ1,"r1"),
 Segment("bead spatial barcode A","A"*7,"cbc",placeholder=True,feature=feature("slide_spatial_a","spatial_barcode","combinatorial",group="slide_bead",part="block 1",whitelist="sequenced bead array")),
 Segment("bead linker","TCTTCAGCGTTCCCGAGA"),Segment("bead spatial barcode B","B"*7,"cbc",placeholder=True,feature=feature("slide_spatial_b","spatial_barcode","combinatorial",group="slide_bead",part="block 2",whitelist="sequenced bead array")),
 Segment("bridge/ligation handle","GCTCGGACACATGGGCG"),Segment("custom Read 1 + ME",READ1,"r1"),Segment("genomic fragment","X"*42,placeholder=True),Segment("custom Read 2 site",revcomp(READ2),"r2"),
 Segment("i7 reverse complement","I"*8,"cbc",placeholder=True,feature=feature("slide_sample_i7","sample_index","unknown")),Segment("P7 reverse complement",il.P7_RC,"p7")],name="slide-DNA-seq library")
SEQ_PRIMERS=(sp.custom("Read 1","slide-DNA Read 1",READ1,"Supplementary Table 1"),sp.custom("Read 2","slide-DNA Read 2",READ2,"Supplementary Table 1"))
_e=sp.verify(FINAL_LIBRARY,SEQ_PRIMERS,required_roles=("Read 1","Read 2"))
if _e: raise ValueError("slide-DNA library: "+"; ".join(_e))
TITLE="slide-DNA-seq — spatial capture of tissue genomic DNA"
NOTES="01_slide-dna-seq.html"
SOURCE='<a href="https://doi.org/10.1038/s41586-021-04217-4">Zhao et al., <i>Nature</i> (2022)</a>.'
SUMMARY="Histone-depleted tissue DNA is tagmented with custom ends directly on a spatially indexed bead array. Photocleaved bead barcodes ligate to nearby genomic fragments before PCR."
CAVEAT="The two genomic read primers and bead-oligo architecture are source-printed. Sample-index read-primer placement is not asserted because the paper prints PCR primers rather than a complete strand-resolved index-read scheme."
FINAL_CAPTION="Representative captured genomic fragment with the two-block bead barcode and custom paired-end sequencing sites."
SEQUENCING_INTRO="The source-printed custom Read 1 and Read 2 primers enter genomic DNA from opposite custom Tn5 ends."
SEQUENCING_UNAVAILABLE="i5/i7 sample indexes are added in PCR, but the source does not give a complete index-read primer map on the final strand; only the exact genomic read sites are placed."
def oligos():
 yield oligo("Tn5 adapter 1",[Segment("custom handle + ME",T5,"me")],mods="5′ phosphate")
 yield oligo("Tn5 adapter 2",[Segment("custom handle + ME",T7,"me")])
 yield oligo("sequencing primer Read 1",[Segment("primer",READ1,"r1")])
 yield oligo("sequencing primer Read 2",[Segment("primer",READ2,"r2")])
def sections(): return [
 ("Place tissue on a sequenced bead array",[Row(chunks=[("10-µm tissue section over beads with known spatial barcodes",None,False)])],"Each bead's barcode has already been assigned a coordinate by sequencing-by-ligation."),
 ("Remove histones and tagment DNA",[Row(chunks=[("custom Tn5 end ** genomic fragment ** custom Tn5 end", "me",False)])],"HCl exposes DNA and custom-loaded Tn5 fragments it in situ."),
 ("Photocleave and ligate bead barcodes",[Row(chunks=[("bead barcode A — linker — barcode B ** proximal genomic fragment", "cbc",False)])],"The covalent junction preserves local spatial identity after material leaves the slide."),
 ("PCR-complete the library",[Row(chunks=[("P5/i5 — bead barcode — genomic insert — i7/P7",None,False)])],"Custom paired-end primers read the genomic insert."),]
