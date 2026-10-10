"""snmC2T-seq, the transcriptome+methylome core of snmCAT-seq."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from batch_ngs import seg,truseq_library
from chemdraw import Row,feature

FINAL_LIBRARY,SEQ_PRIMERS=truseq_library([
 seg("6-nt inline cell barcode","B"*6,"cbc",placeholder=True,feature=feature("snmc2t_cell","cell_barcode","whitelist",whitelist="snmC-seq2 random-primer set")),
 seg("random N9","N"*9,placeholder=True),
 seg("RNA-derived fully mC cDNA or genomic DNA","X"*42,placeholder=True),
 seg("Adaptase-added bases","X"*3,placeholder=True,inferred=True),
],"snmC2T-seq mixed DNA/cDNA library",dual_index=True,inferred_adapters=True)
TITLE="snmC2T-seq — transcriptome and methylome in one bisulfite library"
NOTES="01_snmc2t-seq.html"
SOURCE='<a href="https://doi.org/10.1016/j.xgen.2022.100107">Luo et al., <i>Cell Genomics</i> (2022)</a> (published snmCAT-seq study).'
SUMMARY="Smart-seq reverse transcription substitutes 5-methyl-dCTP, making RNA-derived cDNA fully methylated. Genomic DNA and cDNA then share one snmC-seq2 bisulfite library and are separated by read-level non-CG methylation."
CAVEAT="INFERRED — snmC-seq2 adapter families are named by the source, but proprietary Adaptase-added bases and one complete final strand are not printed."
FINAL_CAPTION="One mixed library: highly methylated cDNA reads report RNA, while bisulfite-converted genomic reads report the methylome."
SEQUENCING_INTRO="Read 1 begins with the snmC-seq2 inline cell barcode and N9. Dual indexes identify the well; RNA-versus-DNA identity is computed from non-CG conversion, not from a separate library barcode."
def sections(): return [
 ("Optional NOMe marking",[Row(chunks=[("M.CviPI: accessible GpC → GmC", "w1",False)])],"The broader snmCAT-seq workflow adds accessibility; snmC2T omits this step."),
 ("Full-length cDNA synthesis with 5-methyl-dCTP",[Row(chunks=[("RNA → fully cytosine-methylated cDNA", "tso",False)])],"Smart-seq/Smart-seq2 chemistry copies nuclear RNA while substituting 5-methyl-dCTP for dCTP."),
 ("Bisulfite conversion of mixed cDNA and genomic DNA",[Row(chunks=[("fully mC cDNA remains C-rich | genomic unmethylated C → U/T",None,False)])],"The two molecule classes stay in the same well and library."),
 ("snmC-seq2 random priming and PCR",[Row(chunks=[("inline barcode — N9 — mixed insert — Adaptase tail — dual-index library", "cbc",False)])],"Read-level mCH partitions transcriptome and methylome after sequencing."),]
