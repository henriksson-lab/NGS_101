"""Original TAIL-seq molecular model."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from batch_ngs import seg,truseq_library
from chemdraw import Row
from terminal_rna import three_prime_adapter_scene,terminal_fragment_selection_rows

TITLE="TAIL-seq — poly(A)-tail length and terminal-base sequencing"
NOTES="01_tail-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1016/j.molcel.2014.02.007">Chang et al., <i>Molecular Cell</i> (2014)</a>.'
SUMMARY="A biotinylated 3-prime adapter preserves native RNA termini through partial RNase T1 digestion and capture; paired reads identify the transcript and traverse its complete poly(A) tail."
CAVEAT="The defining paper is authoritative for reaction order and paired-read roles. Adapter bases not available in the fetched source are left inferred rather than copied from later TAIL-seq variants."
FINAL_LIBRARY,SEQ_PRIMERS=truseq_library([seg("gene-proximal RNA fragment","X"*34,placeholder=True),seg("native poly(A) tail","A"*24),seg("native terminal additions","X"*6,placeholder=True),seg("3-prime adapter","D"*18,"r2",placeholder=True)],"TAIL-seq library",inferred_adapters=True)
FINAL_CAPTION="Read 1 identifies the gene-proximal fragment; the long opposite read crosses terminal additions and the poly(A) tract from the native 3-prime end."
SEQUENCING_INTRO="Paired-end sequencing assigns the transcript with one read and measures the tail with raw fluorescence from the opposite read."
READ_LENGTHS={"Read 1":51,"Read 2":251}

def sections(): return [
 ("Ligate a biotinylated adapter to native RNA 3-prime ends",three_prime_adapter_scene(biotin=True).rows(),"The adapter fixes the precise RNA terminus and supplies the affinity handle."),
 ("Partially digest with RNase T1 and retain terminal fragments",terminal_fragment_selection_rows(),"Only fragments linked to the biotinylated 3-prime adapter survive streptavidin purification."),
 ("Repair the fragment 5-prime end and ligate the second RNA adapter",[Row(chunks=[("5-prime RNA adapter — gene fragment — poly(A) — 3-prime adapter",None,False)])],"The selected terminal RNA fragment now has priming sites at both ends."),
 ("Reverse-transcribe and PCR-complete the paired-end library",[Row(chunks=[("terminal RNA → cDNA → indexed Illumina library",None,False)])],"One read identifies the gene and the other reads the full tail and terminal non-A bases."),
]
