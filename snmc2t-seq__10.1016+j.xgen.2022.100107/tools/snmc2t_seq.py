"""snmC2T-seq, the transcriptome+methylome core of snmCAT-seq."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from batch_ngs import seg,truseq_library
from chemdraw import Row,feature
from multimodal_spatial import finish_pcr_branch
from protocol_graph import InputRole,MaterialKind,ProtocolGraph,StepKind

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


def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 nucleus=g.input("nucleus","one nucleus containing genomic DNA and nuclear RNA",MaterialKind.BIOLOGICAL,state=sections()[1][1],primary=True)
 rt_primer=g.input("rt_primer","Smart-seq full-length reverse-transcription primer",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (exact oligo not transcribed here)")
 rt=g.input("rt","template-switching reverse transcriptase",MaterialKind.CATALYST)
 methylated_dctp=g.input("methylated_dctp","5-methyl-dCTP substitution condition",MaterialKind.CONDITION)
 copied=g.product("copied","fully cytosine-methylated RNA-derived cDNA beside genomic DNA",MaterialKind.MOLECULAR,state=sections()[1][1],inferred=True,source=SOURCE+" (exact RT-primer architecture not transcribed here)")
 g.react("reverse_transcribe","Synthesize full-length cDNA while substituting 5-methyl-dCTP",StepKind.REVERSE_TRANSCRIPTION,((nucleus,InputRole.TEMPLATE),(rt_primer,InputRole.STRUCTURAL),(rt,InputRole.CATALYST),(methylated_dctp,InputRole.CONDITION)),(copied,))
 tso=g.input("tso","Smart-seq template-switch oligo",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (exact oligo not transcribed here)")
 switched=g.product("switched","template-switched fully methylated cDNA mixed with genomic DNA",MaterialKind.MOLECULAR,state=copied.state,inferred=True,source=SOURCE+" (exact TSO architecture not transcribed here)")
 g.react("template_switch","Template-switch the nascent methylated cDNA",StepKind.TEMPLATE_SWITCHING,((copied,InputRole.TEMPLATE),(tso,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(switched,))
 bisulfite=g.input("bisulfite","bisulfite-conversion condition",MaterialKind.CONDITION)
 converted=g.product("converted","mixed converted genomic DNA and C-rich methylated cDNA",MaterialKind.MOLECULAR,state=sections()[2][1])
 g.react("bisulfite_convert","Bisulfite-convert the mixed cDNA and genomic DNA",StepKind.CHEMICAL_CONVERSION,((switched,InputRole.TEMPLATE),(bisulfite,InputRole.CONDITION)),(converted,))
 random_primer=g.input("random_primer","snmC-seq2 inline-barcoded N9 random primer",MaterialKind.OLIGO)
 polymerase=g.input("strand_polymerase","random-primer extension polymerase",MaterialKind.CATALYST)
 primed=g.product("primed","inline-barcoded mixed cDNA/genomic-DNA products",MaterialKind.MOLECULAR,state=sections()[3][1])
 g.react("random_prime","Random-prime both molecule classes with the indexed N9 primer",StepKind.EXTENSION,((converted,InputRole.TEMPLATE),(random_primer,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(primed,))
 adaptase=g.input("adaptase","proprietary Adaptase adapter-completion reagent",MaterialKind.ADAPTER,inferred=True,source=SOURCE+" (added bases not disclosed)")
 adaptase_enzyme=g.input("adaptase_enzyme","Adaptase completion enzyme",MaterialKind.CATALYST)
 adapted=g.product("adapted","Adaptase-completed mixed-library PCR template",MaterialKind.MOLECULAR,state=FINAL_LIBRARY,inferred=True,source=SOURCE+" (Adaptase-added bases not disclosed)")
 g.react("adaptase_complete","Add the proprietary Adaptase completion end",StepKind.ADAPTER_ATTACHMENT,((primed,InputRole.TEMPLATE),(adaptase,InputRole.STRUCTURAL),(adaptase_enzyme,InputRole.CATALYST)),(adapted,))
 reviewed=["nucleus","rt_primer","rt","methylated_dctp","reverse_transcribe","tso","template_switch","bisulfite","bisulfite_convert","random_primer","strand_polymerase","random_prime","adaptase","adaptase_enzyme","adaptase_complete"]
 finish_pcr_branch(g,"final",adapted,FINAL_LIBRARY,SEQ_PRIMERS,reviewed,label="snmC2T-seq mixed RNA/DNA",inferred_primers=True,limitation="The source names the snmC-seq2 dual-index PCR family, but its exact primer sequences are not transcribed in this compact module.")
 g.declare_reviewed_scope(reviewed,SOURCE)
 return g.require_publishable()
