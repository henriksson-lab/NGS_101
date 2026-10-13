"""sci-Plex nuclear-hashing model."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from batch_ngs import seg,truseq_library
from chemdraw import Row,feature
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, StepKind)

FINAL_LIBRARY,SEQ_PRIMERS=truseq_library([
 seg("sci-RNA cell barcode combination","C"*20,"cbc",placeholder=True,feature=feature("sciplex_cell","cell_barcode","combinatorial",group="cell",part="sci-RNA rounds")),
 seg("UMI","U"*8,"umi",placeholder=True,feature=feature("sciplex_umi","umi","random")),
 seg("hash barcode or cDNA insert","X"*36,placeholder=True),
],"sci-Plex library",dual_index=True)
TITLE="sci-Plex — nuclear hashing before combinatorial RNA indexing"
NOTES="01_sci-plex.html"
SOURCE='<a href="https://doi.org/10.1126/science.aax6234">Srivatsan et al., <i>Science</i> (2020)</a>.'
SUMMARY="Each treated sample is labeled with a well-specific polyadenylated ssDNA hash oligo. Pooled nuclei then undergo sci-RNA-seq, so hash oligos and endogenous transcripts receive the same combinatorial cell indexes."
FINAL_CAPTION="A sci-RNA-seq product representing either endogenous cDNA or a captured nuclear-hash oligo; both retain the same cell barcode combination and UMI."
SEQUENCING_INTRO="The combinatorial index combination identifies a nucleus; hash-derived reads recover treatment/sample identity and transcript-derived reads measure expression."
SHOW_READ_CONTENT=True
def sections(): return [
 ("Treat samples in separate wells",[Row(chunks=[("sample 1 | sample 2 | … | sample 4992",None,False)])],"Perturbation identity initially exists as physical well identity."),
 ("Hash permeabilized nuclei",[Row(chunks=[("well-specific ssDNA hash barcode — poly(A) → nucleus", "cbc",False)])],"Unmodified polyadenylated DNA associates with exposed nuclei and is fixed in place."),
 ("Pool and perform sci-RNA-seq",[Row(chunks=[("hash oligo + endogenous mRNA → same RT/ligation/PCR index combination",None,False)])],"The poly(A) tail makes the hash compatible with the transcript workflow."),]

def protocol_graph():
 source="Defining source: https://doi.org/10.1126/science.aax6234"
 graph=ProtocolGraph("sci-Plex",authority=source)
 nuclei=graph.input("nuclei","permeabilized nuclei from treated samples",
                    MaterialKind.BIOLOGICAL,primary=True)
 hash_oligo=graph.input("hash_oligo","well-specific polyadenylated ssDNA hash",
                        MaterialKind.MOLECULAR)
 hashed=graph.product("hashed_nuclei","hash-labelled nuclei",MaterialKind.BIOLOGICAL)
 unbound=graph.product("unbound_hash","unbound hash oligo",MaterialKind.WASTE)
 graph.react("hash_binding","Associate and fix the well-specific hash to nuclei",
             StepKind.BINDING,
             ((nuclei,InputRole.TEMPLATE),(hash_oligo,InputRole.ANALYTE)),(hashed,unbound))
 graph.discard(unbound,"unbound hash oligo is removed before pooling")
 rt_primer=graph.input("indexed_rt_primer","sci-RNA-seq indexed RT primer",MaterialKind.OLIGO)
 rt=graph.input("reverse_transcriptase","reverse transcriptase",MaterialKind.CATALYST)
 cdna=graph.product("indexed_cdna","hash- and transcript-derived indexed cDNA",
                    MaterialKind.MOLECULAR)
 graph.react("reverse_transcription","Pool nuclei and reverse-transcribe hashes and mRNA",
             StepKind.REVERSE_TRANSCRIPTION,
             ((hashed,InputRole.TEMPLATE),(rt_primer,InputRole.STRUCTURAL),
              (rt,InputRole.CATALYST)),(cdna,))
 ligation_index=graph.input("ligation_index","sci-RNA-seq ligation-round index oligo",
                            MaterialKind.ADAPTER)
 indexed=graph.product("combinatorial_cdna","combinatorially indexed cDNA and hash products",
                       MaterialKind.MOLECULAR)
 graph.react("ligation_indexing","Add the split-pool ligation index",
             StepKind.ADAPTER_ATTACHMENT,
             ((cdna,InputRole.TEMPLATE),(ligation_index,InputRole.STRUCTURAL)),(indexed,))
 p5=graph.input("p5_primer","P5 library PCR primer",MaterialKind.OLIGO)
 p7=graph.input("p7_primer","indexed P7 library PCR primer",MaterialKind.OLIGO)
 polymerase=graph.input("polymerase","library PCR polymerase",MaterialKind.CATALYST)
 library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=FINAL_LIBRARY)
 graph.react("library_pcr","PCR-complete the sci-Plex library",StepKind.PCR,
             ((indexed,InputRole.TEMPLATE),(p5,InputRole.STRUCTURAL),
              (p7,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(library,),
             geometry=PCRGeometry.unavailable(
                 FINAL_LIBRARY,reason="The compact sci-RNA branch model records the completed library rather than separate primer landing objects."))
 graph.sequence(library,SequencingPlan(
     library=FINAL_LIBRARY,primers=SEQ_PRIMERS,
     acquired_roles=tuple(p.role for p in SEQ_PRIMERS),
     read_uses={p.role:("demultiplex" if "Index" in p.role else "extract_then_include") for p in SEQ_PRIMERS}))
 graph.declare_reviewed_scope((
     "nuclei","hash_oligo","hash_binding","indexed_rt_primer","reverse_transcriptase",
     "reverse_transcription","ligation_index","ligation_indexing","p5_primer",
     "p7_primer","polymerase","library_pcr","library"),source)
 return graph.require_publishable()
