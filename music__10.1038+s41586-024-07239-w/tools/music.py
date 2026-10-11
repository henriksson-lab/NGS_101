"""MUSIC molecular model."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
import illumina as il
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct,MolecularState,Row,Scene,Segment,Workflow,feature,oligo
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
FINAL_LIBRARIES=(("Combined RNA/DNA library",FINAL,PRIMERS,"One molecule contains either an RNA or DNA insert; its linker distinguishes the analyte.","Read 1 gives the 16-bp GEM complex barcode and 12-bp UMI; i7 completes complex identity; Read 2 traverses three cell barcodes, linker and insert."),)
READ_LENGTHS={"Combined RNA/DNA library":{"Read 1":28,"Index 1 (i7)":8,"Read 2":150}}
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
 wf=Workflow(MolecularState("Crosslinked RNA and fragmented DNA",(Row(chunks=[("RNA … crosslinked complex … dA-tailed DNA",None,False)]),)))
 wf.ligate_adapter("Ligate the analyte-specific DNA linker",dna_linker_geometry(),[Row(chunks=[("DNA — dA:T ** Y linker","me",False)])],note="The linker sequence distinguishes DNA and carries a five-base UMI.")
 wf.react("Attach the chimeric RNA linker",[Row(chunks=[("RNA ** chimeric RNA linker",None,False)])],note="RNA ligation is a single-strand end reaction, not a duplex-overhang ligation.")
 wf.react("Ligate three cell barcodes",[Row(chunks=[("round 1 ** round 2 ** round 3 ** linker ** insert","cbc",False)])],note="Three 96-member sets create 96³ cell identities with ordered complementary overhangs.")
 wf.react("Preserve molecular complexes",[Row(chunks=[("crosslinked complex: RNA + DNA + DNA",None,False)])],note="Associated molecules remain together through nuclear fragmentation.")
 wf.react("Add complex barcodes in 10x GEMs",[Row(chunks=[("16-nt GEM barcode + 12-nt UMI + 8-nt i7 → complex identity","cbc",False)])],note="Molecules sharing GEM and i7 identify the same molecular complex.")
 return [(r.action,r.after.rows,r.note) for r in wf.reactions]
