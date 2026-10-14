"""INSIGHT barcoded NASBA-to-NGS chemistry."""
import illumina as il
import seqprimers as sp
from batch_ngs import seg
from chemdraw import Construct, Row, feature, revcomp
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, StepKind)

TITLE="INSIGHT — barcoded NASBA followed by pooled sequencing"
NOTES="01_insight.html"
SOURCE='Defining source: <a href="https://doi.org/10.1126/sciadv.abe5054">Wu et al., <i>Science Advances</i> (2021)</a>.'
SUMMARY="Paired 5-nt sample barcodes enter a SARS-CoV-2 amplicon during isothermal NASBA; completed reactions are pooled and receive flow-cell adapters and run indexes in one PCR."
CUSTOM_R1="AGCCAGCTCTGGAGAATTCTAATACGACTCACTATAGGGAGAAGG"
TARGET_F="CCAGCAACTGTTTGTGGACCTA"; TARGET_R="ATGGTTTAACAGGCACAGGTGT"
FINAL_LIBRARY=Construct([
 seg("P5",il.P5,"p5"),seg("i5","J"*8,"cbc",placeholder=True,feature=feature("run_i5","sample_index","random")),
 seg("custom T7-containing read arm",CUSTOM_R1,"r1"),
 seg("right sample barcode","B"*5,"cbc",placeholder=True,feature=feature("nasba_right","sample_index","combinatorial",group="nasba_sample",part="right")),
 seg("viral target amplicon","X"*46,placeholder=True),
 seg("left sample barcode","C"*5,"cbc",placeholder=True,feature=feature("nasba_left","sample_index","combinatorial",group="nasba_sample",part="left")),
 seg("TruSeq Read 2 arm reverse complement",revcomp(il.TRUSEQ_READ2),"r2"),
 seg("i7 reverse complement","I"*8,"cbc",placeholder=True,feature=feature("run_i7","sample_index","random")),
 seg("P7 reverse complement",il.P7_RC,"p7")],name="INSIGHT sequencing amplicon")
R1=sp.custom("Read 1","INSIGHT T7-containing sequencing primer",CUSTOM_R1,"Wu et al. 2021, Table 2")
SEQ_PRIMERS=(R1,sp.TRUSEQ["I1"],sp.TRUSEQ["R2"])
FINAL_CAPTION="The NASBA barcode pair identifies the original sample; the outer PCR supplies run indexes and flow-cell ends."
SEQUENCING_INTRO="INSIGHT uses its T7-containing custom primer from the P5 side and a standard TruSeq primer from the P7 side."
def oligos(): return [
 "<b>Barcoded NASBA forward primer</b> — 5′-TGACTGGAGTTCAGACGTGTGCTCTTCCGATCTNNNNNCCAGCAACTGTTTGTGGACCTA-3′",
 "<b>Barcoded NASBA reverse primer</b> — 5′-AATTCTAATACGACTCACTATAGGGAGAAGGNNNNNACACCTGTGCCTGTTAAACCAT-3′",
 "<b>Indexed P5 library primer</b> — 5′-AATGATACGGCGACCACCGAGATCTACACNNNNNNNNAGCCAGCTCTGGAGAATTCTAATACGACTCACTATAGGGAGAAGG-3′",
 "<b>Indexed P7 library primer</b> — 5′-CAAGCAGAAGACGGCATACGAGATNNNNNNNNGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT-3′",
 "<b>Custom sequencing primer</b> — 5′-AGCCAGCTCTGGAGAATTCTAATACGACTCACTATAGGGAGAAGG-3′"]
def sections(): return [
 ("Reverse transcribe the viral RNA",[Row(chunks=[("viral RNA ← barcoded forward primer","cbc",False)])],"The forward primer contributes one sample barcode and an Illumina handle."),
 ("RNase H and second-strand synthesis create a T7 template",[Row(chunks=[("Illumina arm—left barcode—viral target—right barcode—T7 promoter","cbc",False)])],"The reverse primer supplies the second barcode and T7 promoter."),
 ("T7 transcription and repeated NASBA cycling",[Row(chunks=[("dsDNA template → many ssRNA amplicons → dsDNA templates","r1",False)])],"Reverse transcriptase, RNase H and T7 RNA polymerase amplify isothermally."),
 ("Pool completed reactions and perform universal PCR",[Row(chunks=[("P5—run index—T7 read arm—barcode pair—target—TruSeq arm—run index—P7","cbc",False)])],"Residual NASBA primers are removed before the minimal-cycle pooled PCR.")]

def protocol_graph():
 source="Defining source: https://doi.org/10.1126/sciadv.abe5054"
 graph=ProtocolGraph("INSIGHT",authority=source)
 viral_rna=graph.input("viral_rna","SARS-CoV-2 target RNA",MaterialKind.BIOLOGICAL,primary=True)
 forward=graph.input("forward_primer","barcoded NASBA forward primer",MaterialKind.OLIGO)
 rt=graph.input("reverse_transcriptase","NASBA reverse transcriptase",MaterialKind.CATALYST)
 first=graph.product("first_strand","forward-barcode-bearing first-strand cDNA",
                     MaterialKind.MOLECULAR)
 graph.react("reverse_transcription","Reverse-transcribe viral RNA from the barcoded forward primer",
             StepKind.REVERSE_TRANSCRIPTION,
             ((viral_rna,InputRole.TEMPLATE),(forward,InputRole.STRUCTURAL),
              (rt,InputRole.CATALYST)),(first,))
 rnaseh=graph.input("rnase_h","RNase H",MaterialKind.CATALYST)
 cleared=graph.product("rna_removed_cdna","RNA-depleted first-strand cDNA",MaterialKind.MOLECULAR)
 graph.react("rnase_h_digest","Remove RNA from the RNA–cDNA hybrid",StepKind.DIGESTION,
             ((first,InputRole.TEMPLATE),(rnaseh,InputRole.CATALYST)),(cleared,))
 reverse=graph.input("reverse_primer","barcoded T7-promoter NASBA reverse primer",
                     MaterialKind.OLIGO)
 polymerase=graph.input("dna_polymerase","NASBA DNA polymerase",MaterialKind.CATALYST)
 template=graph.product("t7_template","two-barcode double-stranded T7 template",
                        MaterialKind.MOLECULAR)
 graph.react("second_strand","Synthesize the second strand and install the T7 promoter",
             StepKind.EXTENSION,
             ((cleared,InputRole.TEMPLATE),(reverse,InputRole.STRUCTURAL),
              (polymerase,InputRole.CATALYST)),(template,))
 t7=graph.input("t7_polymerase","T7 RNA polymerase",MaterialKind.CATALYST)
 amplicons=graph.product("nasba_amplicons","NASBA-amplified sample-barcoded products",
                         MaterialKind.MOLECULAR)
 graph.react("nasba_transcription","Transcribe and recycle products through NASBA",
             StepKind.TRANSCRIPTION,
             ((template,InputRole.TEMPLATE),(t7,InputRole.CATALYST)),(amplicons,))
 cleanup=graph.input("primer_cleanup","residual-primer removal reagent",MaterialKind.SELECTION)
 clean=graph.product("clean_amplicons","pooled clean NASBA amplicons",MaterialKind.MOLECULAR)
 primer_waste=graph.product("residual_primers","removed residual NASBA primers",MaterialKind.WASTE)
 graph.react("cleanup","Pool completed reactions and remove residual NASBA primers",
             StepKind.CLEANUP,
             ((amplicons,InputRole.TEMPLATE),(cleanup,InputRole.SELECTION)),(clean,primer_waste))
 graph.discard(primer_waste,"residual NASBA primers are removed before pooled PCR")
 p5=graph.input("p5_primer","universal indexed P5 PCR primer",MaterialKind.OLIGO)
 p7=graph.input("p7_primer","universal indexed P7 PCR primer",MaterialKind.OLIGO)
 pcr_polymerase=graph.input("pcr_polymerase","minimal-cycle PCR polymerase",MaterialKind.CATALYST)
 library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=FINAL_LIBRARY)
 graph.react("library_pcr","Add flow-cell adapters and run indexes by universal PCR",
             StepKind.PCR,
             ((clean,InputRole.TEMPLATE),(p5,InputRole.STRUCTURAL),
              (p7,InputRole.STRUCTURAL),(pcr_polymerase,InputRole.CATALYST)),(library,),
             geometry=PCRGeometry.unavailable(
                 FINAL_LIBRARY,reason="The final construct and custom sequencing primer are exact, while the page does not expose separate universal-PCR landing objects."))
 graph.sequence(library,SequencingPlan(
     library=FINAL_LIBRARY,primers=SEQ_PRIMERS,
     acquired_roles=tuple(p.role for p in SEQ_PRIMERS),
     read_uses={p.role:("demultiplex" if "Index" in p.role else "include") for p in SEQ_PRIMERS}))
 graph.declare_reviewed_scope((
     "viral_rna","forward_primer","reverse_transcriptase","reverse_transcription",
     "rnase_h","rnase_h_digest","reverse_primer","dna_polymerase","second_strand",
     "t7_polymerase","nasba_transcription","primer_cleanup","cleanup","p5_primer",
     "p7_primer","pcr_polymerase","library_pcr","library"),source)
 return graph.require_publishable()
