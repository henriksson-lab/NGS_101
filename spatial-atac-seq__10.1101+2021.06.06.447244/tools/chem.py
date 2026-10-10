"""Spatial-ATAC-seq orthogonal in-tissue barcoding."""
import illumina as il
import nextera as nx
import seqprimers as sp
from batch_ngs import seg
from chemdraw import Construct, Row, feature
TITLE="Spatial-ATAC-seq — spatial chromatin accessibility"
NOTES="01_spatial-atac-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1101/2021.06.06.447244">Deng et al., bioRxiv (2021)</a>; chemistry checked against the published version.'
SUMMARY="Tn5 tags accessible DNA in a tissue section; two perpendicular microfluidic ligations append an A×B coordinate before release and PCR."
CUSTOM_R1="GCCTGTCCGCGGAAGCAGTGGTATCAACGCAGAGTAC"
FINAL_LIBRARY=Construct([seg("P5",il.P5,"p5"),seg("custom Read 1 site",CUSTOM_R1,"r1"),
 seg("spatial barcode B","B"*8,"cbc",placeholder=True,feature=feature("spatial_b","spatial_barcode","whitelist",group="spatial_coordinate",part="B",whitelist="Deng et al. Supplementary Table 3")),
 seg("ligation linker 2","T"*15,placeholder=True),
 seg("spatial barcode A","A"*8,"cbc",placeholder=True,feature=feature("spatial_a","spatial_barcode","whitelist",group="spatial_coordinate",part="A",whitelist="Deng et al. Supplementary Table 3")),
 seg("left mosaic end",nx.ME,"me"),seg("accessible genomic fragment","X"*42,placeholder=True),
 seg("right mosaic end reverse complement",nx.ME_RC,"me"),seg("S7 reverse complement",nx.S7_RC,"s7"),
 seg("i7 reverse complement","I"*8,"cbc",placeholder=True,feature=feature("sample_index_i7","sample_index","unknown")),seg("P7 reverse complement",il.P7_RC,"p7")],name="Spatial-ATAC-seq library")
R1=sp.custom("Read 1","Spatial-ATAC custom Read 1",CUSTOM_R1,"Deng et al. Supplementary Table 3")
SEQ_PRIMERS=(R1,sp.NEXTERA["I1"],sp.NEXTERA["R2"])
FINAL_CAPTION="The barcode A×B pair records the tissue-pixel coordinate adjacent to a tagmented accessible-DNA fragment."
SEQUENCING_INTRO="A custom Read 1 primer reads the spatial-barcode side; the opposite end uses Nextera-compatible index and Read 2 primers."
def oligos(): return ["<b>Tn5ME-A</b> — 5′-/5Phos/CATCGGCGTACGACTAGATGTGTATAAGAGACAG-3′","<b>Tn5ME-B</b> — 5′-GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG-3′","<b>Custom Read 1 primer</b> — 5′-GCCTGTCCGCGGAAGCAGTGGTATCAACGCAGAGTAC-3′"]
def sections(): return [
 ("Tagment accessible chromatin in the fixed section",[Row(chunks=[("tissue chromatin ** linker-bearing Tn5 ends","me",False)])],"Adapters carrying the first ligation linker enter accessible loci."),
 ("Ligate barcode A through parallel channels",[Row(chunks=[("Tn5 end ** barcode A ** linker 2","cbc",False)])],"Each channel assigns one coordinate along the first axis."),
 ("Ligate barcode B through perpendicular channels",[Row(chunks=[("Tn5 end ** barcode A ** barcode B","cbc",False)])],"The intersection of A and B identifies one spatial pixel."),
 ("Reverse crosslinks and amplify",[Row(chunks=[("P5—A×B coordinate—accessible DNA—index—P7","r1",False)])],"Released spatially tagged fragments receive complete sequencing ends by PCR."),]
