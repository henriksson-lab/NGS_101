"""MUSIC molecular model."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
import illumina as il
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct,MolecularState,Row,Scene,Segment,Workflow,feature,oligo
from multimodal_spatial import finish_pcr_branch
from interaction_capture import FixedMolecularComplex
from protocol_graph import InputRole,MaterialKind,ProtocolGraph,StepKind
RNA_LINKER="CGAGGAGCGCTTNNNNN"+"ATAGCATTGC"
DNA_TOP="CTAGACACTGTGCGTATCTNBA"+"A"*30
DNA_BOTTOM="CGAGGAGNNNNNACAACGCACAGTGTCTAGT"
FINAL=Construct([Segment("P5",il.P5,"p5"),Segment("Read 1 arm",il.TRUSEQ_READ1,"r1"),Segment("10x GEM complex barcode","G"*16,"cbc",placeholder=True,feature=feature("music_gem_complex","molecular_complex_barcode","whitelist",group="music_complex",part="GEM",whitelist="10x 3-prime whitelist")),Segment("10x UMI","U"*12,"umi",placeholder=True,feature=feature("music_complex_umi","umi","random")),Segment("cell barcode round 3","C"*14,"cbc",placeholder=True,feature=feature("music_cell_3","cell_barcode","combinatorial",group="music_cell",part="round 3",whitelist="96 barcodes")),Segment("cell barcode round 2","B"*14,"cbc",placeholder=True,feature=feature("music_cell_2","cell_barcode","combinatorial",group="music_cell",part="round 2",whitelist="96 barcodes")),Segment("cell barcode round 1","A"*14,"cbc",placeholder=True,feature=feature("music_cell_1","cell_barcode","combinatorial",group="music_cell",part="round 1",whitelist="96 barcodes")),Segment("RNA or DNA linker","L"*17,placeholder=True),Segment("RNA or DNA insert","X"*48,placeholder=True),Segment("dA junction","A"),Segment("Index 1 / Read 2 arm",il.INDEX1_PRIMER,"r2"),Segment("i7 complex barcode reverse complement","I"*8,"cbc",placeholder=True,feature=feature("music_i7_complex","molecular_complex_barcode","whitelist",group="music_complex",part="i7",whitelist="eight indexes per library")),Segment("P7 reverse complement",il.P7_RC,"p7")],name="MUSIC RNA-or-DNA library molecule")
PRIMERS=(sp.TRUSEQ["R1"],sp.TRUSEQ["I1"],sp.TRUSEQ["R2"])
_e=sp.verify(FINAL,PRIMERS,required_roles=("Read 1","Index 1 (i7)","Read 2"))
if _e:raise ValueError("MUSIC library: "+"; ".join(_e))
TITLE="MUSIC — multinucleic-acid interaction mapping in single cells"
NOTES="01_music.html"
SOURCE='<a href="https://doi.org/10.1038/s41586-024-07239-w">Zhao et al., <i>Nature</i> (2024)</a>.'
SUMMARY="RNA and fragmented DNA receive distinct linkers, three rounds of cell barcodes and a two-part molecular-complex barcode. Reads sharing both identities reveal multiway DNA–DNA and RNA–DNA complexes."
CAVEAT="The DNA linker's single-base sticky end is source-resolved. This page has not yet verified the RNA-linker junction, the strand-resolved three-round cell-barcode overhangs or the exact final PCR-primer pair from the supplementary protocol; those geometries remain explicitly unavailable."
FINAL_LIBRARIES=(("Combined RNA/DNA library",FINAL,PRIMERS,"One molecule contains either an RNA or DNA insert; its linker distinguishes the analyte.","Read 1 gives the 16-bp GEM complex barcode and 12-bp UMI; i7 completes complex identity; Read 2 traverses three cell barcodes, linker and insert."),)
READ_LENGTHS={"Combined RNA/DNA library":{"Read 1":28,"Index 1 (i7)":8,"Read 2":150}}
COMPLEX=FixedMolecularComplex(("RNA molecule","DNA fragment A","DNA fragment B"),"formaldehyde",label="multiway nucleic-acid complex")
def oligos():
 yield oligo("RNA linker",[Segment("DNA portion + 5-nt UMI",RNA_LINKER[:17],"umi",placeholder=True),Segment("RNA portion",RNA_LINKER[17:])],mods="5′-OH; final 10 bases are RNA; 3′-OH")
 yield oligo("DNA linker top",[Segment("phosphorylated top",DNA_TOP,placeholder=True)],mods="5′ phosphate")
 yield oligo("DNA linker bottom",[Segment("bottom with 5-nt UMI",DNA_BOTTOM,"umi",placeholder=True)])
def dna_linker_geometry():
 sc=Scene(); sc.strand("linker",[Segment("Y-linker body","[published Y-linker duplex]",placeholder=True),Segment("single-base dT","T")],label="DNA linker")
 sc.ligation_end("linker","single-base dT",kind="cohesive",polarity="3-prime",owner="adapter",text="published 3′-dT sticky end")
 sub=Scene(); sub.strand("DNA",[Segment("dA-tailed genomic DNA","X"*24,placeholder=True),Segment("genomic 3-prime dA","A")],label="DNA")
 sub.ligation_end("DNA","genomic 3-prime dA",kind="cohesive",polarity="3-prime",owner="substrate",text="3′-dA")
 return AdapterLigationGeometry.known("cohesive",sc,sub,statement="The published MUSIC Y linker has a single-base sticky end complementary to dA-tailed genomic DNA; its full top and bottom oligos remain copyable above.")

def sections():
 wf=Workflow(MolecularState("Crosslinked RNA and fragmented DNA",tuple(COMPLEX.complex().rows())))
 wf.ligate_adapter("Ligate the analyte-specific DNA linker",dna_linker_geometry(),[Row(chunks=[("DNA — dA:T ** Y linker","me",False)])],note="The linker sequence distinguishes DNA and carries a five-base UMI.")
 wf.react("Attach the chimeric RNA linker",[Row(chunks=[("RNA ** chimeric RNA linker",None,False)])],note="RNA ligation is a single-strand end reaction, not a duplex-overhang ligation.")
 wf.react("Ligate three cell barcodes",[Row(chunks=[("round 1 ** round 2 ** round 3 ** linker ** insert","cbc",False)])],note="Three 96-member sets create 96³ cell identities with ordered complementary overhangs.")
 wf.react("Preserve molecular complexes",COMPLEX.complex().rows(),note="Associated RNA and DNA remain separate molecules retained in one fixed complex through nuclear fragmentation.")
 wf.react("Add complex barcodes in 10x GEMs",[Row(chunks=[("16-nt GEM barcode + 12-nt UMI + 8-nt i7 → complex identity","cbc",False)])],note="Molecules sharing GEM and i7 identify the same molecular complex.")
 return [(r.action,r.after.rows,r.note) for r in wf.reactions]


def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 complexes=g.input("complexes","crosslinked molecular complexes containing RNA and fragmented dA-tailed DNA",MaterialKind.BIOLOGICAL,state=sections()[0][1],primary=True)
 dna_linker=g.input("dna_linker","published UMI-bearing Y DNA linker",MaterialKind.ADAPTER,state=dna_linker_geometry().scenes[0])
 dna_ligase=g.input("dna_ligase","DNA-linker ligase",MaterialKind.CATALYST)
 dna_linked=g.product("dna_linked","complexes with the analyte-specific DNA linker attached",MaterialKind.MOLECULAR,state=sections()[0][1])
 g.react("ligate_dna_linker","Ligate the Y linker to dA-tailed genomic DNA",StepKind.LIGATION,((complexes,InputRole.TEMPLATE),(dna_linker,InputRole.STRUCTURAL),(dna_ligase,InputRole.CATALYST)),(dna_linked,),geometry=dna_linker_geometry())

 rna_linker=g.input("rna_linker","published chimeric RNA linker",MaterialKind.ADAPTER,inferred=True,source=SOURCE+" (junction orientation not yet verified in this page model)")
 rna_ligase=g.input("rna_ligase","single-strand RNA ligase",MaterialKind.CATALYST)
 both_linked=g.product("both_linked","complexes with analyte-specific RNA and DNA linkers",MaterialKind.MOLECULAR,state=sections()[1][1],inferred=True,source=SOURCE+" (RNA-linker junction orientation not yet verified)")
 g.react("ligate_rna_linker","Attach the chimeric linker to RNA",StepKind.LIGATION,((dna_linked,InputRole.TEMPLATE),(rna_linker,InputRole.STRUCTURAL),(rna_ligase,InputRole.CATALYST)),(both_linked,),geometry=AdapterLigationGeometry.unavailable("The chimeric linker sequence is published, but its strand-resolved ligating junction has not yet been verified from the supplementary protocol."))

 barcode_ligase=g.input("barcode_ligase","cell-barcode DNA ligase",MaterialKind.CATALYST)
 current=both_linked
 reviewed=["complexes","dna_linker","dna_ligase","ligate_dna_linker","rna_linker","rna_ligase","ligate_rna_linker","barcode_ligase"]
 for number in (1,2,3):
  adapter=g.input(f"cell_barcode_{number}",f"round-{number} cell-barcode adapter",MaterialKind.ADAPTER,inferred=True,source=SOURCE+" (strand-resolved complementary overhang not yet transcribed)")
  indexed=g.product(f"cell_indexed_{number}",f"complex-associated RNA and DNA after cell-barcode round {number}",MaterialKind.MOLECULAR,state=sections()[2][1],inferred=True,source=SOURCE+" (barcode identity known; overhang geometry not yet transcribed)")
  g.react(f"ligate_cell_barcode_{number}",f"Ligate cell-barcode round {number}",StepKind.LIGATION,((current,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),(barcode_ligase,InputRole.CATALYST)),(indexed,),geometry=AdapterLigationGeometry.unavailable(f"Round {number} is source-defined, but its strand-resolved complementary overhang has not yet been verified from the supplementary barcode table."))
  reviewed.extend((adapter.id,f"ligate_cell_barcode_{number}"));current=indexed

 gem_condition=g.input("gem_condition","10x GEM partitioning condition",MaterialKind.CONDITION)
 partitioned=g.product("partitioned","cell-indexed molecular complexes isolated in GEMs",MaterialKind.MOLECULAR,state=sections()[4][1])
 collisions=g.product("collisions","empty or multiply occupied GEM material",MaterialKind.WASTE)
 g.react("partition","Partition molecular complexes into 10x GEMs",StepKind.PARTITION,((current,InputRole.TEMPLATE),(gem_condition,InputRole.CONDITION)),(partitioned,collisions),changes_molecule=False)
 g.discard(collisions,"does not preserve an interpretable one-GEM molecular-complex identity")
 gem_oligo=g.input("gem_oligo","10x GEM complex-barcode and UMI oligo",MaterialKind.OLIGO)
 complex_indexed=g.product("complex_indexed","RNA-or-DNA templates carrying cell and molecular-complex identities",MaterialKind.MOLECULAR,state=FINAL)
 g.react("attach_complex_barcode","Attach the 16-nt GEM barcode and 12-nt UMI",StepKind.ADAPTER_ATTACHMENT,((partitioned,InputRole.TEMPLATE),(gem_oligo,InputRole.STRUCTURAL)),(complex_indexed,))

 reviewed.extend(("gem_condition","partition","gem_oligo","attach_complex_barcode"))
 finish_pcr_branch(g,"combined",complex_indexed,FINAL,PRIMERS,reviewed,label="MUSIC combined RNA/DNA",read_lengths=READ_LENGTHS["Combined RNA/DNA library"],inferred_primers=True,limitation="The exact final PCR-primer pair has not yet been verified from the supplementary protocol.")
 g.declare_reviewed_scope(reviewed,SOURCE)
 return g.require_publishable()
