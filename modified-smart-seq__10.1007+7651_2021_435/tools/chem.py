"""Xi Chen laboratory Modified SMART-Seq chemistry."""
import illumina as il
import nextera as nx
import seqprimers as sp
from batch_ngs import seg
from chemdraw import Construct, Row, feature
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)

TITLE="Modified SMART-Seq — 3′ counting with well barcodes and UMIs"
NOTES="01_modified-smart-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1007/7651_2021_435">Zheng et al., <i>Methods in Molecular Biology</i> (2022)</a>.'
SUMMARY="A barcoded oligo(dT) installs an 8-nt well barcode and 10-nt UMI during reverse transcription; template switching and cDNA PCR precede tagmentation and indexed library PCR."
RT="CTACACGACGCTCTTCCGATCT"
TSO="AAGCAGTGGTATCAACGCAGAGTACATGGG"
FINAL_LIBRARY=Construct([
    seg("P5",il.P5,"p5"),seg("Read 1 arm",RT,"r1"),
    seg("well barcode","B"*8,"cbc",placeholder=True,feature=feature(
        "well_barcode","cell_barcode","whitelist",
        whitelist="Zheng et al. 2022, Table 1 well-barcode set")),
    seg("UMI","U"*10,"umi",placeholder=True,feature=feature("umi","umi","random")),
    seg("transcript cDNA","X"*36,placeholder=True),seg("TSO-derived end",TSO),
    seg("mosaic end reverse complement",nx.ME_RC,"me"),seg("S7 reverse complement",nx.S7_RC,"s7"),
    seg("plate index reverse complement","I"*8,"cbc",placeholder=True,feature=feature(
        "plate_index","sample_index","whitelist",
        whitelist="Zheng et al. 2022, Table 1 plate-barcode set")),
    seg("P7 reverse complement",il.P7_RC,"p7")],name="Modified SMART-Seq library")
SEQ_PRIMERS=(sp.TRUSEQ["R1"],sp.NEXTERA["I1"],sp.NEXTERA["R2"])
FINAL_CAPTION="Read 1 begins with the well barcode and UMI; Read 2 reads transcript sequence from the tagmented end."
SEQUENCING_INTRO="The protocol specifies at least 18 Read 1 cycles, an 8-cycle plate-index read and at least 50 Read 2 cycles."
READ_LENGTHS={"Read 1":18,"Index 1 (i7)":8,"Read 2":75}
SHOW_READ_CONTENT=True
def oligos(): return [
 "<b>Barcoded RT primer</b> — 5′-CTACACGACGCTCTTCCGATCT[8-bp well barcode]NNNNNNNNNNTTTTTTTTTTTTTTTTTTTTTTTTTVN-3′",
 "<b>Template-switching oligo</b> — 5′-AAGCAGTGGTATCAACGCAGAGTACATrGrGrG-3′",
 "<b>cDNA PCR forward / reverse</b> — 5′-CTACACGACGCTCTTCCGATCT-3′ / 5′-AAGCAGTGGTATCAACGCAGAG-3′",
 "<b>Library PCR universal</b> — 5′-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTC-3′",
 "<b>Library PCR indexed</b> — 5′-CAAGCAGAAGACGGCATACGAGAT[8-bp plate barcode]GTCTCGTGGGCTCGG-3′"]
def sections(): return [
 ("Reverse transcription installs well barcode and UMI",[Row(chunks=[("poly(A) RNA ← oligo(dT)—UMI—well barcode—Read 1 arm","umi",False)])],"Each well receives its own RT primer before cells are sorted."),
 ("Template switch at the cDNA 5′ end",[Row(chunks=[("first-strand cDNA ** TSO-derived PCR handle","r1",False)])],"The terminal riboguanosines support template switching."),
 ("Pool wells after cDNA PCR and Exonuclease I cleanup",[Row(chunks=[("well-tagged full-length cDNA pool",None,False)])],"Primer removal precedes pooling and suppresses cross-well mispriming."),
 ("Tagment and add the plate index by library PCR",[Row(chunks=[("P5—well barcode—UMI—cDNA ** mosaic end—plate index—P7","cbc",False)])],"Tn5 fragments amplified cDNA; the Nextera-style index primer identifies the plate.")]

def protocol_graph():
 source="Defining source: https://doi.org/10.1007/7651_2021_435"
 graph=ProtocolGraph("Modified SMART-Seq",authority=source)
 rna=graph.input("rna","polyadenylated single-cell RNA",MaterialKind.BIOLOGICAL,primary=True)
 rt_primer=graph.input("barcoded_rt_primer","barcoded oligo(dT) RT primer",MaterialKind.OLIGO)
 rt=graph.input("reverse_transcriptase","template-switching reverse transcriptase",MaterialKind.CATALYST)
 first=graph.product("first_strand","well-barcode/UMI-tagged first-strand cDNA",MaterialKind.MOLECULAR)
 graph.react("reverse_transcription","Reverse-transcribe from the barcoded oligo(dT)",
             StepKind.REVERSE_TRANSCRIPTION,
             ((rna,InputRole.TEMPLATE),(rt_primer,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(first,))
 tso=graph.input("tso","template-switching oligo",MaterialKind.OLIGO)
 switched=graph.product("switched_cdna","template-switched full-length cDNA",MaterialKind.MOLECULAR)
 graph.react("template_switch","Install the 5-prime PCR handle by template switching",
             StepKind.TEMPLATE_SWITCHING,
             ((first,InputRole.TEMPLATE),(tso,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(switched,))
 cdna_forward=graph.input("cdna_pcr_forward","Table 1 cDNA PCR forward primer",MaterialKind.OLIGO)
 cdna_reverse=graph.input("cdna_pcr_reverse","Table 1 cDNA PCR reverse primer",MaterialKind.OLIGO)
 polymerase=graph.input("pcr_polymerase","cDNA/library PCR polymerase",MaterialKind.CATALYST)
 amplified=graph.product("amplified_cdna","well-tagged amplified full-length cDNA",MaterialKind.MOLECULAR)
 graph.react("cdna_pcr","Amplify the template-switched cDNA",StepKind.PCR,
             ((switched,InputRole.TEMPLATE),(cdna_forward,InputRole.STRUCTURAL),
              (cdna_reverse,InputRole.STRUCTURAL),
              (polymerase,InputRole.CATALYST)),(amplified,),
             geometry=PCRGeometry.unavailable(FINAL_LIBRARY,
                 reason="The cDNA-PCR intermediate is not represented base-for-base.",
                 symmetric_species=True))
 cleanup=graph.input("exonuclease_cleanup","Exonuclease I",MaterialKind.CATALYST)
 cleaned=graph.product("cleaned_cdna","cleaned pooled cDNA",MaterialKind.MOLECULAR)
 waste=graph.product("free_primers","digested residual primers",MaterialKind.WASTE)
 graph.react("cleanup","Remove primers with Exonuclease I and pool wells",StepKind.CLEANUP,
             ((amplified,InputRole.TEMPLATE),(cleanup,InputRole.CATALYST)),(cleaned,waste))
 graph.discard(waste,"residual primers are digested before well pooling")
 transposome=graph.input("transposome","Nextera transposome",MaterialKind.ADAPTER)
 tn5=graph.input("tn5","Tn5 transposase",MaterialKind.CATALYST)
 tagged=graph.product("tagmented_cdna","tagmented well-barcoded cDNA",MaterialKind.MOLECULAR)
 twaste=graph.product("tagmentation_waste","non-library tagmentation products",MaterialKind.WASTE)
 graph.react("tagment","Tagment pooled cDNA",StepKind.TAGMENTATION,
             ((cleaned,InputRole.TEMPLATE),(transposome,InputRole.STRUCTURAL),(tn5,InputRole.CATALYST)),
             (tagged,twaste))
 graph.discard(twaste,"fragments without the selected RT-primer and Nextera ends are not amplified")
 forward=graph.input("library_forward","RT-handle library PCR primer",MaterialKind.OLIGO)
 reverse=graph.input("indexed_reverse","plate-index Nextera PCR primer",MaterialKind.OLIGO)
 library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=FINAL_LIBRARY)
 graph.react("library_pcr","Select RT-handle/Tn5 fragments and add the plate index",StepKind.PCR,
             ((tagged,InputRole.TEMPLATE),(forward,InputRole.STRUCTURAL),
              (reverse,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(library,),
             geometry=PCRGeometry.unavailable(FINAL_LIBRARY,
                 reason="The tagmented pre-PCR fragment is variable; the final primer sites are canonical."))
 graph.sequence(library,SequencingPlan(
     library=FINAL_LIBRARY,primers=SEQ_PRIMERS,
     acquired_roles=tuple(p.role for p in SEQ_PRIMERS),read_lengths=READ_LENGTHS,
     read_uses={"Read 1":"extract","Index 1 (i7)":"demultiplex","Read 2":"include"}))
 graph.declare_reviewed_scope((
     "rna","barcoded_rt_primer","reverse_transcriptase","reverse_transcription","tso",
     "template_switch","cdna_pcr_forward","cdna_pcr_reverse","pcr_polymerase","cdna_pcr",
     "exonuclease_cleanup","cleanup","transposome","tn5","tagment",
     "library_forward","indexed_reverse","library_pcr","library"),source)
 return graph.require_publishable()
