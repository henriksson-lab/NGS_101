"""Spear-ATAC accessibility and guide-capture branches."""
import illumina as il
import nextera as nx
import seqprimers as sp
from batch_ngs import seg
from chemdraw import Construct, Row, feature, oligo, revcomp
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
TITLE="Spear-ATAC — single-cell ATAC CRISPR screening"
NOTES="01_spear-atac.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/s41467-021-23213-w">Pierce et al., <i>Nature Communications</i> (2021)</a>.'
SUMMARY="A 10x scATAC-derived workflow captures accessible chromatin and enriches guide identity from the same barcoded nuclei."
CAVEAT=("The paper and its supplementary protocol print the bead oligo, enrichment "
        "primers, custom sequencing primers and both completed library structures. "
        "Variable cell, guide and sample indexes remain typed placeholders, not missing bases.")

CELL=feature("cell_barcode","cell_barcode","whitelist",whitelist="10x scATAC 16-nt whitelist")
GUIDE_BC=feature("guide_barcode","guide_barcode","whitelist",whitelist="Spear-ATAC 10-nt guide-barcode set")
I7=feature("sample_index_i7","sample_index","whitelist",whitelist="selected 8-nt i7 set")

OSP1594="AATGATACGGCGACCACCGAGA"
OSP2053="GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTAAGTATCCCTTGGAGAACCACCTTG"
OMCB1672="GCCACTTTTTCAAGTTGATAACGGACTAGCCTTATTTAAACTTGCTATGCTGTTTCCAGCTTAGCTCTTAAAC"
GUIDE_AFTER_BARCODE="CGAGGCTGAGTGTAGATTCGAGCAAAAAGCACCGACTCGGT"
GUIDE_PRIMER_SITE=OMCB1672
GUIDE_DOWNSTREAM="CCAACAAGGTGGTTCTCCAAGGGATACTT"
P7_TRUSEQ_R2=il.P7+"N"*8+il.TRUSEQ_READ2

ATAC=Construct([
 seg("P5",il.P5,"p5"),seg("cell barcode","N"*16,"cbc",placeholder=True,feature=CELL),
 seg("S5",nx.S5,"s5"),seg("left mosaic end",nx.ME,"me"),
 seg("accessible chromatin fragment","X"*60,placeholder=True),
 seg("right mosaic end reverse complement",nx.ME_RC,"me"),seg("S7 reverse complement",nx.S7_RC,"s7"),
 seg("i7 reverse complement","N"*8,"cbc",placeholder=True,feature=I7),seg("P7 reverse complement",il.P7_RC,"p7")],
 name="Spear-ATAC chromatin library")
AP=(sp.NEXTERA["R1"],sp.NEXTERA["I1"],sp.NEXTERA["I2"],sp.NEXTERA["R2"])

GUIDE=Construct([
 seg("P5",il.P5,"p5"),seg("cell barcode","N"*16,"cbc",placeholder=True,feature=CELL),
 seg("S5",nx.S5,"s5"),seg("left mosaic end",nx.ME,"me"),
 seg("guide barcode","N"*10,"cbc",placeholder=True,feature=GUIDE_BC),
 seg("guide constant region upstream",GUIDE_AFTER_BARCODE),seg("custom Read 1 primer site",GUIDE_PRIMER_SITE,"r1"),
 seg("guide spacer","X"*20,"guide",placeholder=True),seg("guide constant region downstream",GUIDE_DOWNSTREAM),
 seg("TruSeq Read 2 reverse complement",revcomp(il.TRUSEQ_READ2),"r2"),
 seg("i7 reverse complement","N"*8,"cbc",placeholder=True,feature=I7),seg("P7 reverse complement",il.P7_RC,"p7")],
 name="Spear-ATAC guide library")
GUIDE_R1=sp.custom("Read 1","oMCB1672 custom guide Read 1",OMCB1672,"Pierce et al. 2021 Supplementary Note 3")
GUIDE_I1=sp.custom("Index 1 (i7)","guide sample-index primer",revcomp(il.TRUSEQ_READ2),"Pierce et al. 2021 Supplementary Note 3")
GUIDE_I2=sp.NEXTERA["I2"]
GP=(GUIDE_R1,GUIDE_I1,GUIDE_I2)
for library,primers in ((ATAC,AP),(GUIDE,GP)):
 problems=sp.verify(library,primers,required_roles=tuple(p.role for p in primers))
 if problems: raise ValueError(library.name+": "+"; ".join(problems))

FINAL_LIBRARIES=(("Chromatin-accessibility library",ATAC,AP,"Tagmented chromatin fragments carry the droplet cell barcode.","Published Nextera primers read the ATAC branch."),("Guide-identity library",GUIDE,GP,"Guide enrichment preserves the same droplet cell barcode.","Published custom primers read guide spacer, sample index and cell barcode."))

def oligos():
 yield oligo("oSP1594 P5 PCR primer",[seg("P5 PCR",OSP1594,"p5")])
 yield oligo("oSP2053 biotinylated guide linear-PCR primer",[seg("guide enrichment",OSP2053,"r2")],five="5′-/Bio/-",three="-3′")
 yield oligo("oMCB1672 custom guide Read 1 primer",[seg("guide Read 1",OMCB1672,"r1")])
 yield oligo("indexed P7–TruSeq Read 2 PCR primer",[seg("P7",il.P7,"p7"),seg("i7","N"*8,"cbc",placeholder=True,feature=I7),seg("TruSeq Read 2",il.TRUSEQ_READ2,"r2")])
def sections(): return [
 ("Tagment nuclei before partitioning",[Row(chunks=[("accessible chromatin ** Tn5 mosaic ends","me",False)])],"Bulk tagmentation marks open chromatin."),
 ("Barcode nuclei in droplets",[Row(chunks=[("ATAC fragment + gel-bead cell barcode","cbc",False)])],"Gap fill transfers one droplet identity to accessible fragments."),
 ("Amplify the shared barcoded material",[Row(chunks=[("barcoded ATAC products → PCR pool","r1",False)])],"The pool contains genomic ATAC molecules and guide-associated material."),
 ("Complete the two sequencing libraries",[Row(chunks=[("ATAC library | guide-specific linear PCR → guide library","r2",False)])],"The supplementary protocol prints the terminal primers and completed library structures."),]

def protocol_graph():
 source="Defining source: https://doi.org/10.1038/s41467-021-23213-w"
 graph=ProtocolGraph("Spear-ATAC",authority=source)
 nuclei=graph.input("nuclei","CRISPR-screen nuclei",MaterialKind.BIOLOGICAL,primary=True)
 transposome=graph.input("transposome","10x scATAC-compatible transposome",MaterialKind.ADAPTER)
 tn5=graph.input("tn5","Tn5 transposase",MaterialKind.CATALYST)
 tagged=graph.product("tagged_nuclei","nuclei containing tagmented accessible chromatin",
                      MaterialKind.MOLECULAR)
 twaste=graph.product("tagmentation_waste","non-library tagmentation products",MaterialKind.WASTE)
 graph.react("tagment","Tagment nuclei before partitioning",StepKind.TAGMENTATION,
             ((nuclei,InputRole.TEMPLATE),(transposome,InputRole.STRUCTURAL),
              (tn5,InputRole.CATALYST)),(tagged,twaste))
 graph.discard(twaste,"fragments lacking usable library ends do not amplify")
 loading=graph.input("gem_loading","10x droplet loading condition",MaterialKind.CONDITION)
 gems=graph.product("productive_gems","GEMs containing nuclei and barcode beads",MaterialKind.MOLECULAR)
 excluded=graph.product("excluded_gems","empty or multiplet GEMs",MaterialKind.WASTE)
 graph.react("partition","Partition tagmented nuclei with gel-bead barcodes",StepKind.PARTITION,
             ((tagged,InputRole.TEMPLATE),(loading,InputRole.CONDITION)),(gems,excluded))
 graph.discard(excluded,"empty and multiplet GEMs are excluded from accepted cells")
 bead=graph.input("bead_primer","published P5–16-nt-cell-barcode–S5 gel-bead oligo",MaterialKind.OLIGO,
                  source=source)
 opposite=graph.input("opposite_primer","published opposite-end ATAC amplification primer",MaterialKind.OLIGO,
                      source=source)
 polymerase=graph.input("polymerase","gap-fill and PCR polymerase",MaterialKind.CATALYST)
 shared=graph.product("shared_pool","shared cell-barcoded ATAC/guide-associated product pool",
                      MaterialKind.MOLECULAR,source=source)
 graph.react("gem_pcr","Gap-fill and amplify the shared barcoded material",StepKind.PCR,
             ((gems,InputRole.TEMPLATE),(bead,InputRole.STRUCTURAL),
              (opposite,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(shared,),
             geometry=PCRGeometry.unavailable(ATAC,
                 reason="The paper defines the shared barcode and branches but not one contiguous pre-branch molecule."))
 atac_in=graph.product("atac_branch","chromatin-accessibility enrichment branch",MaterialKind.MOLECULAR)
 guide_in=graph.product("guide_branch","guide-identity enrichment branch",MaterialKind.MOLECULAR)
 graph.react("branch","Divide the shared pool into ATAC and guide enrichments",StepKind.SPLIT,
             ((shared,InputRole.TEMPLATE),),(atac_in,guide_in))
 atac_f=graph.input("atac_forward","SI-PCR Primer B / P5-side primer",MaterialKind.OLIGO,
                    source=source)
 atac_r=graph.input("atac_reverse","indexed N7xx ATAC primer",MaterialKind.OLIGO,
                    source=source)
 atac_lib=graph.product("atac_library",ATAC.name,MaterialKind.LIBRARY,state=ATAC,
                        source=source)
 graph.react("atac_pcr","Complete the chromatin-accessibility sequencing library",StepKind.PCR,
             ((atac_in,InputRole.TEMPLATE),(atac_f,InputRole.STRUCTURAL),
              (atac_r,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(atac_lib,),
             geometry=PCRGeometry.exact(
                 ATAC,PCRPrimer("oSP1594 / SI-PCR Primer B",OSP1594),
                 PCRPrimer("indexed N7xx",nx.n7xx_primer("N"*8)),
                 statement="The supplementary protocol prints the P5-side primer and indexed N7xx architecture."))
 guide_linear=graph.input("guide_linear_primer","5′-biotinylated oSP2053 guide-specific linear-extension primer",MaterialKind.OLIGO,source=source)
 guide_product=graph.product("guide_enriched","guide-specific enriched products",MaterialKind.MOLECULAR)
 graph.react("guide_extension","Linearly enrich guide identity from the shared pool",StepKind.EXTENSION,
             ((guide_in,InputRole.TEMPLATE),(guide_linear,InputRole.STRUCTURAL),
              (polymerase,InputRole.CATALYST)),(guide_product,))
 guide_f=graph.input("guide_forward","oSP1594 P5 PCR primer",MaterialKind.OLIGO,source=source)
 guide_r=graph.input("guide_reverse","indexed P7–TruSeq Read 2 primer",MaterialKind.OLIGO,source=source)
 guide_lib=graph.product("guide_library",GUIDE.name,MaterialKind.LIBRARY,state=GUIDE,
                         source=source)
 graph.react("guide_pcr","Complete the guide-identity sequencing library",StepKind.PCR,
             ((guide_product,InputRole.TEMPLATE),(guide_f,InputRole.STRUCTURAL),
              (guide_r,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(guide_lib,),
             geometry=PCRGeometry.exact(
                 GUIDE,PCRPrimer("oSP1594",OSP1594),
                 PCRPrimer("indexed P7–TruSeq Read 2",P7_TRUSEQ_R2),
                 statement="The paper prints oSP1594 and the complete indexed P7–TruSeq Read 2 primer."))
 for material,library,primers in ((atac_lib,ATAC,AP),(guide_lib,GUIDE,GP)):
  graph.sequence(material,SequencingPlan(
      library=library,primers=primers,acquired_roles=tuple(p.role for p in primers),
      read_uses={p.role:("demultiplex" if "Index" in p.role else "include") for p in primers},
      limitation="Variable barcode and index members are represented by their published lengths and roles."))
 graph.declare_reviewed_scope((
     "nuclei","transposome","tn5","tagment","gem_loading","partition","bead_primer",
     "opposite_primer","polymerase","gem_pcr","branch","atac_forward","atac_reverse",
     "atac_pcr","guide_linear_primer","guide_extension","guide_forward","guide_reverse",
     "guide_pcr","atac_library","guide_library"),source)
 return graph.require_publishable()
