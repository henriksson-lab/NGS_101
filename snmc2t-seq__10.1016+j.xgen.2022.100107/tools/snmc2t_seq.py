"""snmC2T-seq, the transcriptome+methylome core of snmCAT-seq."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from chemdraw import Row, Segment, oligo
from pcr import PCRGeometry
from multimodal_spatial import finish_pcr_branch
from protocol_graph import Disposition, Terminal, InputRole,MaterialKind,ProtocolGraph,StepKind
import snmc2_core

FINAL_LIBRARY,SEQ_PRIMERS=snmc2_core.library(
 "RNA-derived fully mC cDNA or genomic DNA",insert_nt=42,
 library_name="snmC2T-seq mixed DNA/cDNA library",feature_id="snmc2t_cell",
 whitelist="snmC-seq2 random-primer set")
TITLE="snmC2T-seq — transcriptome and methylome in one bisulfite library"
NOTES="01_snmc2t-seq.html"
SOURCE='<a href="https://doi.org/10.1016/j.xgen.2022.100107">Luo et al., <i>Cell Genomics</i> (2022)</a> (published snmCAT-seq study).'
SUMMARY="Smart-seq reverse transcription substitutes 5-methyl-dCTP, making RNA-derived cDNA fully methylated. This graph selects the nuclear RNA/methylome branch: genomic DNA and cDNA then share one snmC-seq2 bisulfite library and are separated by read-level non-CG methylation."
CAVEAT=("The paper prints the reverse-transcription, template-switch and cDNA-PCR "
        "oligos. Its methylome finish explicitly reuses the source-exact snmC-seq2 "
        "random-primer and dual-index families. INFERRED — only the short bases "
        "installed by Adaptase are not assigned in this model. Whole cultured cells instead use the earlier snmC-seq finish; their 12-cycle cDNA PCR is not the nuclear branch shown here.")
FINAL_CAPTION="One mixed library: highly methylated cDNA reads report RNA, while bisulfite-converted genomic reads report the methylome."
SEQUENCING_INTRO="Read 1 begins with the snmC-seq2 inline cell barcode and N9. Dual indexes identify the well; RNA-versus-DNA identity is computed from non-CG conversion, not from a separate library barcode."

DT30VN = "AAGCAGUGGUAUCAACGCAGAGUACUTTTTTUTTTTTUTTTTTUTTTTTUTTTTTVN"
N6 = "AAGCAGUGGUAUCAACGCAGAGUACNNNNNN"
TSO = "AAGCAGUGGUAUCAACGCAGAGUGAAUrGrG+G"
ISPCR23 = "AAGCAGUGGUAUCAACGCAGAGU"


def oligos():
 yield "<b>dT30VN_4 (cultured cells)</b> — 5′-" + DT30VN + "-3′"
 yield "<b>dT30VN_4 (human brain nuclei)</b> — 5′-/5SpC3/" + DT30VN + "-3′"
 yield "<b>N6_2 random primer</b> — 5′-/5SpC3/" + N6 + "-3′"
 yield "<b>TSO_3</b> — 5′-/5SpC3/" + TSO + "-3′"
 yield "<b>ISPCR23_2</b> — 5′-/5SpC3/" + ISPCR23 + "-3′"
 yield oligo("snmC-seq2 P5L_AD002_H random primer", [
  Segment("partial Read 1 handle", "TTCCCTACACGACGCTCTTCCGATCT"),
  Segment("representative inline barcode", snmc2_core.INLINE, "cbc"),
  Segment("random N9", "N" * 9, placeholder=True)],
  mods="/5SpC3/ before the printed sequence")
 yield oligo("snmC-seq2 dual-index P5 primer 1", [Segment("P5 primer", snmc2_core.P5_INDEX_1)])
 yield oligo("snmC-seq2 dual-index P7 primer 1", [Segment("P7 primer", snmc2_core.P7_INDEX_1)])
def sections(): return [
 ("Optional NOMe marking",[Row(chunks=[("M.CviPI: accessible GpC → GmC", "w1",False)])],"The broader snmCAT-seq workflow adds accessibility; snmC2T omits this step."),
 ("Full-length cDNA synthesis with 5-methyl-dCTP",[Row(chunks=[("RNA → fully cytosine-methylated cDNA", "tso",False)])],"Smart-seq/Smart-seq2 chemistry copies nuclear RNA while substituting 5-methyl-dCTP for dCTP."),
 ("Bisulfite conversion of mixed cDNA and genomic DNA",[Row(chunks=[("fully mC cDNA remains C-rich | genomic unmethylated C → U/T",None,False)])],"The two molecule classes stay in the same well and library."),
 ("snmC-seq2 random priming and PCR",[Row(chunks=[("inline barcode — N9 — mixed insert — Adaptase tail — dual-index library", "cbc",False)])],"Read-level mCH partitions transcriptome and methylome after sequencing."),]


def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 nucleus=g.input("nucleus","one nucleus containing genomic DNA and nuclear RNA",MaterialKind.BIOLOGICAL,state=sections()[1][1],primary=True)
 rt_primer=g.input("rt_primer","source-printed dT30VN_4 and N6_2 reverse-transcription primers",MaterialKind.OLIGO)
 rt=g.input("rt","SuperScript II reverse transcriptase",MaterialKind.CATALYST)
 methylated_dctp=g.input("methylated_dctp","5-methyl-dCTP substitution condition",MaterialKind.CONDITION)
 copied=g.product("copied","fully cytosine-methylated RNA-derived cDNA beside genomic DNA",MaterialKind.MOLECULAR,state=sections()[1][1])
 g.react("reverse_transcribe","Synthesize full-length cDNA while substituting 5-methyl-dCTP",StepKind.REVERSE_TRANSCRIPTION,((nucleus,InputRole.TEMPLATE),(rt_primer,InputRole.STRUCTURAL),(rt,InputRole.CATALYST),(methylated_dctp,InputRole.CONDITION)),(copied,))
 tso=g.input("tso","source-printed TSO_3 template-switch oligo",MaterialKind.OLIGO)
 switched=g.product("switched","template-switched fully methylated cDNA mixed with genomic DNA",MaterialKind.MOLECULAR,state=copied.state)
 g.react("template_switch","Template-switch the nascent methylated cDNA",StepKind.TEMPLATE_SWITCHING,((copied,InputRole.TEMPLATE),(tso,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(switched,))
 cdna_primer=g.input("cdna_primer","ISPCR23_2 cDNA amplification primer",MaterialKind.OLIGO)
 cdna_reverse=g.input("cdna_reverse","ISPCR23_2 at the opposite cDNA end",MaterialKind.OLIGO)
 cdna_polymerase=g.input("cdna_polymerase","KAPA2G Robust HotStart DNA polymerase",MaterialKind.CATALYST)
 amplified=g.product("amplified_cdna","amplified cDNA beside genomic DNA",MaterialKind.MOLECULAR,state=switched.state)
 g.react("cdna_amplify","Amplify cDNA: 15 cycles for cultured nuclei or 14 for brain nuclei",StepKind.PCR,((switched,InputRole.TEMPLATE),(cdna_primer,InputRole.STRUCTURAL),(cdna_reverse,InputRole.STRUCTURAL),(cdna_polymerase,InputRole.CATALYST)),(amplified,),geometry=PCRGeometry.unavailable(FINAL_LIBRARY,reason="The shared ISPCR23_2 species binds both cDNA ends; the complete pre-bisulfite construct is not yet represented."))
 cleavage=g.input("uracil_cleavage","UDG; Endonuclease VIII additionally used in cultured-cell experiments",MaterialKind.CATALYST)
 cleaved=g.product("cleaved_oligos","unused oligos bearing abasic sites or cleavage products",MaterialKind.MOLECULAR)
 prepared=g.product("prepared_mixture","cDNA/genomic mixture after unused-primer treatment",MaterialKind.MOLECULAR,state=amplified.state)
 g.react("remove_unused_primers","Treat unused DNA primers with UDG (brain); UDG/EndoVIII in cultured-cell experiments",StepKind.DIGESTION,((amplified,InputRole.TEMPLATE),(cleavage,InputRole.CATALYST)),(prepared,cleaved))
 g.mark_terminal(Terminal(cleaved.id,Disposition.UNRESOLVED,reason="UDG creates abasic sites; brain samples rely on subsequent alkaline bisulfite desulfonation for cleavage. Physical recovery of these products is not represented."))
 bisulfite=g.input("bisulfite","bisulfite-conversion condition",MaterialKind.CONDITION)
 converted=g.product("converted","mixed converted genomic DNA and C-rich methylated cDNA",MaterialKind.MOLECULAR,state=sections()[2][1])
 g.react("bisulfite_convert","Bisulfite-convert the mixed cDNA and genomic DNA",StepKind.CHEMICAL_CONVERSION,((prepared,InputRole.TEMPLATE),(bisulfite,InputRole.CONDITION)),(converted,))
 random_primer=g.input("random_primer","snmC-seq2 inline-barcoded N9 random primer",MaterialKind.OLIGO)
 polymerase=g.input("strand_polymerase","random-primer extension polymerase",MaterialKind.CATALYST)
 primed=g.product("primed","inline-barcoded mixed cDNA/genomic-DNA products",MaterialKind.MOLECULAR,state=sections()[3][1])
 g.react("random_prime","Random-prime both molecule classes with the indexed N9 primer",StepKind.EXTENSION,((converted,InputRole.TEMPLATE),(random_primer,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(primed,))
 adaptase=g.input("adaptase","proprietary Adaptase adapter-completion reagent",MaterialKind.ADAPTER,inferred=True,source=SOURCE+" (added bases not assigned in this model)")
 adaptase_enzyme=g.input("adaptase_enzyme","Adaptase completion enzyme",MaterialKind.CATALYST)
 adapted=g.product("adapted","Adaptase-completed mixed-library PCR template",MaterialKind.MOLECULAR,state=FINAL_LIBRARY,inferred=True,source=SOURCE+" (Adaptase-added bases not assigned in this model)")
 g.react("adaptase_complete","Add the proprietary Adaptase completion end",StepKind.ADAPTER_ATTACHMENT,((primed,InputRole.TEMPLATE),(adaptase,InputRole.STRUCTURAL),(adaptase_enzyme,InputRole.CATALYST)),(adapted,))
 reviewed=["cdna_primer","cdna_reverse","cdna_polymerase","cdna_amplify","uracil_cleavage","remove_unused_primers","nucleus","rt_primer","rt","methylated_dctp","reverse_transcribe","tso","template_switch","bisulfite","bisulfite_convert","random_primer","strand_polymerase","random_prime","adaptase","adaptase_enzyme","adaptase_complete"]
 finish_pcr_branch(g,"final",adapted,FINAL_LIBRARY,SEQ_PRIMERS,reviewed,
                   label="snmC2T-seq mixed RNA/DNA", read_lengths={"Read 1":150,"Read 2":150},
                   geometry=snmc2_core.representative_pcr_geometry(FINAL_LIBRARY),
                   limitation="Adaptase-added bases are unassigned; source-printed primer regions are explicit.")
 g.declare_reviewed_scope(reviewed,SOURCE)
 return g.require_publishable()
