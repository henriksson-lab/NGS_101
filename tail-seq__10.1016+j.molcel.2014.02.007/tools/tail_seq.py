"""Original TAIL-seq molecular model."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from batch_ngs import seg,truseq_library
from chemdraw import Row
from adapter_ligation import AdapterLigationGeometry
from chemdraw import MolecularState, Scene, Workflow
from terminal_rna import three_prime_adapter_scene,terminal_fragment_selection_rows

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
