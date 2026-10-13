"""Original TAIL-seq molecular model."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from batch_ngs import seg,truseq_library
from chemdraw import Row
from adapter_ligation import AdapterLigationGeometry
from chemdraw import MolecularState, Scene, Workflow
from terminal_rna import three_prime_adapter_scene,terminal_fragment_selection_rows
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from protocol_model import molecule,role

TITLE="TAIL-seq — poly(A)-tail length and terminal-base sequencing"
NOTES="01_tail-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1016/j.molcel.2014.02.007">Chang et al., <i>Molecular Cell</i> (2014)</a>.'
SUMMARY="A biotinylated 3-prime adapter preserves native RNA termini through partial RNase T1 digestion and capture; paired reads identify the transcript and traverse its complete poly(A) tail."
CAVEAT="The defining paper is authoritative for reaction order and paired-read roles. Adapter bases not available in the fetched source are left inferred rather than copied from later TAIL-seq variants."
FINAL_LIBRARY,SEQ_PRIMERS=truseq_library([seg("gene-proximal RNA fragment","X"*34,placeholder=True),seg("native poly(A) tail","A"*24),seg("native terminal additions","X"*6,placeholder=True,length_bp=6),seg("3-prime adapter","D"*18,"r2",placeholder=True,length_bp=18)],"TAIL-seq library",inferred_adapters=True)
FINAL_CAPTION="Read 1 identifies the gene-proximal fragment; the long opposite read crosses terminal additions and the poly(A) tract from the native 3-prime end."
SEQUENCING_INTRO="Paired-end sequencing assigns the transcript with one read and measures the tail with raw fluorescence from the opposite read."
READ_LENGTHS={"Read 1":51,"Read 2":251}
SHOW_READ_CONTENT=True

def workflow():
 w=Workflow(MolecularState("native polyadenylated RNA",(Row(chunks=[("gene RNA—poly(A)—native 3′-OH",None,False)]),)))
 g3=AdapterLigationGeometry.unavailable("The original TAIL-seq source establishes single-strand RNA 3′-end ligation to a biotinylated adapter, but the fetched defining source does not disclose that adapter's terminal bases or activation chemistry; no duplex sticky end is implied.")
 w.ligate_adapter("Ligate the biotinylated 3′ adapter",g3,three_prime_adapter_scene(biotin=True).rows(),name="3′-adapter-ligated RNA",note="The adapter fixes the native terminus and supplies the affinity handle.")
 w.react("Partially digest with RNase T1 and retain terminal fragments",terminal_fragment_selection_rows(),name="captured terminal RNA fragments",note="Only fragments linked to the biotinylated adapter survive purification.")
 g5=AdapterLigationGeometry.unavailable("The original TAIL-seq source establishes single-strand RNA 5′-adapter ligation after end repair, but the fetched defining source does not disclose the adapter terminal bases; this is not a duplex overhang.")
 joined=Scene(); joined.strand("terminal RNA",[seg("5′ RNA adapter","X"*18,placeholder=True),seg("gene fragment","X"*28,placeholder=True),seg("poly(A)","A"*24),seg("3′ adapter","X"*18,placeholder=True)],label="doubly adapter-ligated RNA"); joined.junction("terminal RNA","5′ RNA adapter","gene fragment","RNA ligation")
 w.ligate_adapter("Repair the 5′ end and ligate the second RNA adapter",g5,joined.rows(),name="two-adapter terminal RNA",note="The selected fragment now has priming sites at both ends.")
 w.react("Reverse-transcribe and PCR-complete",[Row(chunks=[("terminal RNA → cDNA → indexed Illumina library",None,False)])],name="TAIL-seq library",note="One read identifies the gene and the other reads the complete tail.")
 return w


def protocol_graph():
 graph=ProtocolGraph(TITLE,authority=SOURCE)
 rna=graph.input("rna","native polyadenylated RNA",MaterialKind.BIOLOGICAL,state=molecule("polyadenylated RNA",(role("gene body"),role("native poly(A) and terminal additions"))),primary=True)
 adapter3=graph.input("three_prime_adapter","biotinylated 3-prime RNA adapter",MaterialKind.ADAPTER,inferred=True,source=SOURCE+" (terminal bases unavailable)")
 ligase3=graph.input("three_prime_ligase","3-prime RNA ligase",MaterialKind.CATALYST)
 tagged=graph.product("tagged","biotin-adapter-tagged native RNA termini",MaterialKind.MOLECULAR,state=three_prime_adapter_scene(biotin=True))
 graph.react("ligate_three_prime","Ligate the biotinylated 3-prime adapter",StepKind.LIGATION,((rna,InputRole.TEMPLATE),(adapter3,InputRole.STRUCTURAL),(ligase3,InputRole.CATALYST)),(tagged,),geometry=AdapterLigationGeometry.unavailable("The fetched defining source does not disclose the single-strand RNA adapter terminal bases or activation chemistry."))
 rnase=graph.input("rnase_t1","RNase T1",MaterialKind.CATALYST)
 fragments=graph.product("fragments","partially digested adapter-bearing terminal RNA fragments",MaterialKind.MOLECULAR,state=tagged.state)
 graph.react("digest","Partially digest RNA with RNase T1",StepKind.DIGESTION,((tagged,InputRole.TEMPLATE),(rnase,InputRole.CATALYST)),(fragments,))
 beads=graph.input("streptavidin","streptavidin capture matrix",MaterialKind.SELECTION)
 captured=graph.product("captured","captured native-terminal RNA fragments",MaterialKind.MOLECULAR,state=fragments.state)
 rejected=graph.product("rejected","non-terminal RNA fragments",MaterialKind.WASTE)
 graph.react("capture","Retain biotin-linked native terminal fragments",StepKind.AFFINITY_CAPTURE,((fragments,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),(captured,rejected))
 graph.discard(rejected,"lacks the biotinylated native 3-prime-end adapter")
 repair=graph.input("end_repair","5-prime RNA-end repair enzyme",MaterialKind.CATALYST)
 repaired=graph.product("repaired","5-prime-repaired captured RNA",MaterialKind.MOLECULAR,state=captured.state)
 graph.react("repair","Repair the selected RNA 5-prime end",StepKind.END_REPAIR,((captured,InputRole.TEMPLATE),(repair,InputRole.CATALYST)),(repaired,))
 adapter5=graph.input("five_prime_adapter","5-prime RNA adapter",MaterialKind.ADAPTER,inferred=True,source=SOURCE+" (terminal bases unavailable)")
 ligase5=graph.input("five_prime_ligase","5-prime RNA ligase",MaterialKind.CATALYST)
 dual=graph.product("dual_adapter","two-adapter terminal RNA",MaterialKind.MOLECULAR,state=repaired.state)
 graph.react("ligate_five_prime","Ligate the second RNA adapter",StepKind.LIGATION,((repaired,InputRole.TEMPLATE),(adapter5,InputRole.STRUCTURAL),(ligase5,InputRole.CATALYST)),(dual,),geometry=AdapterLigationGeometry.unavailable("The fetched defining source does not disclose the single-strand 5-prime adapter terminal bases."))
 rt_primer=graph.input("rt_primer","adapter-complementary RT primer",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (sequence unavailable)")
 rt=graph.input("reverse_transcriptase","reverse transcriptase",MaterialKind.CATALYST)
 cdna=graph.product("cdna","TAIL-seq cDNA",MaterialKind.MOLECULAR,state=dual.state)
 graph.react("reverse_transcribe","Reverse-transcribe the selected terminal RNA",StepKind.REVERSE_TRANSCRIPTION,((dual,InputRole.TEMPLATE),(rt_primer,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(cdna,))
 fwd=graph.input("pcr_forward","library PCR forward primer",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (adapter sequence unavailable)")
 rev=graph.input("pcr_reverse","library PCR reverse primer",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (adapter sequence unavailable)")
 pol=graph.input("pcr_polymerase","library PCR polymerase",MaterialKind.CATALYST)
 library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=FINAL_LIBRARY,inferred=True,source=SOURCE+" (canonical compatible shell inferred)")
 graph.react("library_pcr","Complete the indexed Illumina library",StepKind.PCR,((cdna,InputRole.TEMPLATE),(fwd,InputRole.STRUCTURAL),(rev,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(library,),geometry=PCRGeometry.unavailable(FINAL_LIBRARY,reason="The original adapter bases are not disclosed."))
 graph.sequence(library,SequencingPlan(library=FINAL_LIBRARY,primers=tuple(SEQ_PRIMERS),read_lengths=READ_LENGTHS))
 graph.declare_reviewed_scope(["rna","three_prime_adapter","three_prime_ligase","ligate_three_prime","rnase_t1","digest","streptavidin","capture","end_repair","repair","five_prime_adapter","five_prime_ligase","ligate_five_prime","rt_primer","reverse_transcriptase","reverse_transcribe","pcr_forward","pcr_reverse","pcr_polymerase","library_pcr","library"],SOURCE)
 return graph.require_publishable()
