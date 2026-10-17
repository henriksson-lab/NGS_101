"""Slide-DNA-seq molecular model."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
import illumina as il
import nextera as nx
import seqprimers as sp
from chemdraw import Construct,Row,Segment,feature,oligo,revcomp
from adapter_ligation import AdapterLigationGeometry, reported_adapter_end
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, StepKind)

# Supplementary Table 1 (Supplementary Information file 3, sheet 1), verbatim.
BEAD_OLIGO_5="GCCGCTACACGACGCTCTTCCGATCT"   # after 5'-TTTT-PCT-
BEAD_LINKER="TCTTCAGCGTTCCCGAGA"
BEAD_OLIGO_3="GCTCGGACACATGGGCG"
BRIDGE="CGTTAGCAAAGCTCCGCCCATGTGTC"
MOSAIC_BOTTOM="CTGTCTCTTATACA"               # 5'-Phos … -ddC-3'
PCR_FWD_BULK="AATGATACGGCGACCACCGAGATCTACAC"+"GACACATGGGCGGAGCTTTGCTAACGGTCGAG"
PCR_FWD_SLIDE_TAIL=il.P5
PCR_FWD_SLIDE_ARM=il.TRUSEQ_READ1               # printed with a *T before the terminal base
PCR_REV_TAIL="CAAGCAGAAGACGGCATACGAGAT"
PCR_REV_ARM="GTTGGCACCAGGCTTACGGATGTTGCACCAGC"
READ1="GCTTTGCTAACGGTCGAGAGATGTGTATAAGAGACAG"
READ2="CGGATGTTGCACCAGCAGATGTGTATAAGAGACAG"
INDEX1_SEQ="GCTGGTGCAACATCCGTAAGCCTGGTGCCAAC"
INDEX2_SEQ="CTCGACCGTTAGCAAAGCTCCGCCCATGTGTC"
SOURCE_OLIGO_EXCEPTIONS={
 "p5_primer":"Supplementary Table 1 prints this primer as three parts around a bracketed '[i5 barcode]' placeholder, so no whole strand occurs verbatim; the flanks are checked by eye against that row.",
 "p7_primer":"Supplementary Table 1 prints this primer as three parts around a bracketed '[i7 barcode]' placeholder, so no whole strand occurs verbatim; the flanks are checked by eye against that row.",
}
T5="GAGCTTTGCTAACGGTCGAG"+nx.ME
T7="CTTACGGATGTTGCACCAGC"+nx.ME

# The amplified strand is assembled from those printed oligos, not retyped: the P5 end
# from the slide-DNA-seq PCR forward primer, the barcode block from the bead oligo, the
# genomic flanks from the two Tn5 adapters, and the P7 end from the PCR reverse primer
# (its 3' 32 bases prime on the Tn5 adapter 2 arm, so only the 12-base remainder is new).
_REV_NEW=revcomp(PCR_REV_ARM[:len(PCR_REV_ARM)-len("CTTACGGATGTTGCACCAGC")])
FINAL_LIBRARY=Construct([
 Segment("P5",il.P5,"p5"),
 Segment("i5","J"*8,"cbc",placeholder=True,feature=feature("slide_sample_i5","sample_index","unknown")),
 Segment("TruSeq Read 1 arm",PCR_FWD_SLIDE_ARM,"r1"),
 Segment("bead spatial barcode A","A"*7,"cbc",placeholder=True,feature=feature("slide_spatial_a","spatial_barcode","combinatorial",group="slide_bead",part="block 1",whitelist="sequenced bead array")),
 Segment("bead linker",BEAD_LINKER),
 Segment("bead spatial barcode B","B"*7,"cbc",placeholder=True,feature=feature("slide_spatial_b","spatial_barcode","combinatorial",group="slide_bead",part="block 2",whitelist="sequenced bead array")),
 Segment("bead 3′ ligation handle",BEAD_OLIGO_3),
 Segment("Tn5 adapter 1 handle + ME",T5,"me"),
 Segment("genomic fragment","X"*42,placeholder=True),
 Segment("Tn5 adapter 2 handle + ME, reverse complement",revcomp(T7),"me"),
 Segment("PCR reverse primer arm",_REV_NEW,"r2"),
 Segment("i7 reverse complement","I"*8,"cbc",placeholder=True,feature=feature("slide_sample_i7","sample_index","unknown")),
 Segment("P7 reverse complement",il.P7_RC,"p7")],name="slide-DNA-seq library")
SEQ_PRIMERS=(
 sp.custom("Read 1","slide-DNA Read 1",READ1,"Supplementary Table 1"),
 sp.custom("Index 1 (i7)","slide-DNA index 1 (sample index)",INDEX1_SEQ,"Supplementary Table 1"),
 sp.custom("Index 2 (i5)","slide-DNA index 2 (spatial barcode read)",INDEX2_SEQ,"Supplementary Table 1",
           "Reads the bead barcode, not an i5 sample index: it primes on the bead 3′ handle/Tn5 adapter 1 junction and runs back across barcode block B, the bead linker and block A."),
 sp.custom("Read 2","slide-DNA Read 2",READ2,"Supplementary Table 1"))
_e=sp.verify(FINAL_LIBRARY,SEQ_PRIMERS,required_roles=tuple(p.role for p in SEQ_PRIMERS))
if _e: raise ValueError("slide-DNA library: "+"; ".join(_e))
TITLE="slide-DNA-seq — spatial capture of tissue genomic DNA"
NOTES="01_slide-dna-seq.html"
SOURCE='<a href="https://doi.org/10.1038/s41586-021-04217-4">Zhao et al., <i>Nature</i> (2022)</a>.'
SUMMARY="Histone-depleted tissue DNA is tagmented with custom ends directly on a spatially indexed bead array. Photocleaved bead barcodes ligate to nearby genomic fragments before PCR."
CAVEAT="All four custom sequencing primers, the bead-oligo architecture and both PCR primers are source-printed; the strand below is assembled from them. The genomic insert is representative, and the i5 and i7 sample-index bases are not published."
FINAL_CAPTION="Representative captured genomic fragment with the two-block bead barcode and custom paired-end sequencing sites."
READ_LENGTHS={"Read 1":57,"Read 2":57,"Index 1 (i7)":8,"Index 2 (i5)":44}
SEQUENCING_INTRO="All four primers are the custom ones printed in Supplementary Table 1. Read 1 and Read 2 enter genomic DNA from opposite custom Tn5 ends; index 1 reads the i7 sample index off the PCR reverse arm, while index 2 is not a sample-index read at all — it primes at the bead-handle/Tn5 junction and runs back across the spatial barcode. Supplementary Methods §1.5 gives the slide-DNA-seq run as Read 1: 57 cycles, Read 2: 57 cycles, Index 1: 8 cycles and Index 2: 44 cycles, and 44 cycles is exactly what is needed to clear barcode block B, the 18-base bead linker and block A. The separate bulk WGS control uses a different 40/40/8/8 recipe and the bulk PCR forward primer."
SEQUENCING_UNAVAILABLE=""
def oligos():
 yield oligo("Tn5 adapter 1",[Segment("custom handle + ME",T5,"me")],mods="5′ phosphate")
 yield oligo("Tn5 adapter 2",[Segment("custom handle + ME",T7,"me")])
 yield oligo("bead oligo",[Segment("5′ spacer + handle",BEAD_OLIGO_5),Segment("spatial barcode block 1","[7-nt split-pool barcode]",placeholder=True),Segment("bead linker",BEAD_LINKER),Segment("spatial barcode block 2","[7-nt split-pool barcode]",placeholder=True),Segment("3′ ligation handle",BEAD_OLIGO_3)],mods="5′ TTTT-PCT photocleavable linker")
 yield oligo("bridge oligo",[Segment("bead-to-Tn5 bridge","CGTTAGCAAAGCTCCGCCCATGTGTC")])
 yield oligo("Tn5 mosaic bottom",[Segment("mosaic bottom before dideoxy block","CTGTCTCTTATACA")],mods="5′ phosphate",three="-ddC-3′")
 yield oligo("slide-DNA PCR forward",[Segment("P5 and i5 flank","AATGATACGGCGACCACCGAGATCTACAC"),Segment("i5 index","[i5 barcode]",placeholder=True),Segment("TruSeq R1 site","ACACTCTTTCCCTACACGACGCTCTTCCGATC"),Segment("protected terminal T","*T",placeholder=True)])
 yield oligo("PCR reverse",[Segment("P7 and i7 flank","CAAGCAGAAGACGGCATACGAGAT"),Segment("i7 index","[i7 barcode]",placeholder=True),Segment("custom genomic-side arm","GTTGGCACCAGGCTTACGGATGTTGCACCAGC")])
 yield oligo("sequencing primer Read 1",[Segment("primer",READ1,"r1")])
 yield oligo("sequencing primer Read 2",[Segment("primer",READ2,"r2")])
 yield oligo("sequencing primer index 1",[Segment("primer",INDEX1_SEQ,"r2")])
 yield oligo("sequencing primer index 2 (spatial barcode read)",[Segment("primer",INDEX2_SEQ,"cbc")])
def sections(): return [
 ("Place tissue on a sequenced bead array",[Row(chunks=[("10-µm tissue section over beads with known spatial barcodes",None,False)])],"Each bead's barcode has already been assigned a coordinate by sequencing-by-ligation."),
 ("Remove histones and tagment DNA",[Row(chunks=[("custom Tn5 end ** genomic fragment ** custom Tn5 end", "me",False)])],"HCl exposes DNA and custom-loaded Tn5 fragments it in situ."),
 ("Photocleave and ligate bead barcodes",[Row(chunks=[("bead barcode A — linker — barcode B ** proximal genomic fragment", "cbc",False)])],"The covalent junction preserves local spatial identity after material leaves the slide."),
 ("PCR-complete the library",[Row(chunks=[("P5/i5 — bead barcode — genomic insert — i7/P7",None,False)])],"Custom paired-end primers read the genomic insert."),]

def protocol_graph():
 source="Defining source: https://doi.org/10.1038/s41586-021-04217-4"
 graph=ProtocolGraph("slide-DNA-seq",authority=source)
 tissue=graph.input("tissue","10-µm tissue section",MaterialKind.BIOLOGICAL,primary=True)
 array=graph.input("bead_array","sequenced spatial barcode bead array",MaterialKind.SELECTION)
 placed=graph.product("placed_tissue","tissue registered to the spatial bead array",
                      MaterialKind.BIOLOGICAL)
 outside=graph.product("outside_array","tissue material outside usable bead positions",
                       MaterialKind.WASTE)
 graph.react("array_registration","Place the tissue on the sequenced bead array",
             StepKind.BINDING,
             ((tissue,InputRole.TEMPLATE),(array,InputRole.SELECTION)),(placed,outside),
             changes_molecule=False)
 graph.discard(outside,"material outside registered bead positions has no spatial barcode")
 acid=graph.input("histone_removal","0.1 N HCl, 5 min (after 0.5% Triton-X permeabilization)",MaterialKind.CONDITION)
 exposed=graph.product("exposed_dna","histone-depleted tissue DNA",MaterialKind.MOLECULAR)
 graph.react("deproteinize","Remove histones to expose tissue DNA",StepKind.DEPROTEINIZATION,
             ((placed,InputRole.TEMPLATE),(acid,InputRole.CONDITION)),(exposed,))
 transposome=graph.input("custom_transposome","Tn5 adapter 1 / adapter 2 with the blocked mosaic bottom strand",
                         MaterialKind.ADAPTER,
                         state=(T5,T7,MOSAIC_BOTTOM))
 tn5=graph.input("tn5","Tn5 transposase",MaterialKind.CATALYST)
 tagged=graph.product("tagged_dna","custom-end-tagmented spatial tissue DNA",
                      MaterialKind.MOLECULAR)
 tag_waste=graph.product("tagmentation_waste","non-library tagmentation products",
                         MaterialKind.WASTE)
 graph.react("tagmentation","Tagment histone-depleted DNA in situ (37 °C, 1 h; repeated 4× in the increased-tagmentation variant)",StepKind.TAGMENTATION,
             ((exposed,InputRole.TEMPLATE),(transposome,InputRole.STRUCTURAL),
              (tn5,InputRole.CATALYST)),(tagged,tag_waste))
 graph.discard(tag_waste,"molecules without usable paired custom insertions do not amplify")
 bead_barcode=graph.input("bead_barcode","photocleaved two-block spatial bead oligo (365 nm, 3 min)",
                          MaterialKind.ADAPTER,
                          state=BEAD_OLIGO_5+"N"*7+BEAD_LINKER+"N"*7+BEAD_OLIGO_3)
 bridge=graph.input("bridge","bridge oligo, 10 µM, hybridized 37 °C 2 h in 4× SSC",
                    MaterialKind.OLIGO,state=BRIDGE)
 ligase=graph.input("ligase","T4 DNA ligase, 12 U/µL, 1 h at room temperature",MaterialKind.CATALYST)
 spatial=graph.product("spatial_dna","spatial-barcode-ligated genomic fragments",
                       MaterialKind.MOLECULAR)
 graph.react("spatial_ligation","Photocleave and ligate bead barcodes to proximal DNA",
             StepKind.LIGATION,
             ((tagged,InputRole.TEMPLATE),(bead_barcode,InputRole.STRUCTURAL),
              (bridge,InputRole.STRUCTURAL),(ligase,InputRole.CATALYST)),(spatial,),
             geometry=reported_adapter_end(
                 "slide-DNA-seq bead oligo",kind="splinted",polarity="3-prime",
                 terminal_sequence=BEAD_OLIGO_3,
                 statement=("The bridge oligo is a splint, derived by pairing the printed strands: its 3′ half "
                            "(CGCCCATGTGTC) is complementary to the bead oligo's 3′ handle and its 5′ half "
                            "(CGTTAGCAAAGCTC) to the 5′-phosphorylated Tn5 adapter 1 handle, so T4 ligase seals the "
                            "bead 3′ OH to that adapter 5′ phosphate. The source states the hybridization, "
                            "photocleavage and ligation but does not draw the junction.")))
 p5=graph.input("p5_primer","slide-DNA-seq indexed P5 PCR forward primer",MaterialKind.OLIGO,
                state=PCR_FWD_SLIDE_TAIL+"N"*8+PCR_FWD_SLIDE_ARM)
 p7=graph.input("p7_primer","indexed P7 PCR reverse primer",MaterialKind.OLIGO,
                state=PCR_REV_TAIL+"N"*8+PCR_REV_ARM)
 polymerase=graph.input("polymerase","library PCR polymerase (Buenrostro 2015 ATAC-seq PCR, annealing changed to 65 °C)",MaterialKind.CATALYST)
 library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=FINAL_LIBRARY)
 graph.react("library_pcr","PCR-complete the spatial library",StepKind.PCR,
             ((spatial,InputRole.TEMPLATE),(p5,InputRole.STRUCTURAL),
              (p7,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(library,),
             geometry=PCRGeometry.unavailable(
                 FINAL_LIBRARY,reason="The paper prints PCR primers but not a complete strand-resolved index-read scheme."))
 graph.sequence(library,SequencingPlan(
     library=FINAL_LIBRARY,primers=SEQ_PRIMERS,
     acquired_roles=("Read 1","Index 1 (i7)","Index 2 (i5)","Read 2"),
     read_lengths=READ_LENGTHS,
     read_uses={"Read 1":"include","Read 2":"include",
                "Index 1 (i7)":"demultiplex","Index 2 (i5)":"extract"}))
 graph.declare_reviewed_scope((
     "tissue","bead_array","array_registration","histone_removal","deproteinize",
     "custom_transposome","tn5","tagmentation","bead_barcode","bridge","ligase",
     "spatial_ligation","p5_primer","p7_primer","polymerase","library_pcr","library"),source)
 return graph.require_publishable()
