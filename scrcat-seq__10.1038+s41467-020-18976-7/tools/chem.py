"""scRCAT-seq cap-and-tail endpoint libraries."""
from chemdraw import Row
from multibranch import illumina_branch
from multimodal_spatial import finish_pcr_branch
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind
TITLE="scRCAT-seq — single-cell RNA cap and tail sequencing"
NOTES="01_scrcat-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/s41467-020-18976-7">Liu et al., <i>Nature Communications</i> (2020)</a>.'
SUMMARY="Barcoded, UMI-tagged full-length cDNA is tagmented; targeted PCR separately enriches fragments adjacent to the transcript 5′ and 3′ boundaries."
CAVEAT="The RT tag, barcode lengths and selective end PCR are source-exact. The compact completed Nextera shells are marked inferred pending a full base-by-base reconstruction of both selected end products."
TSS,TP=illumina_branch("scRCAT-seq TSS library","5-prime cDNA boundary",architecture="nextera",cell=16,umi=8,inferred_adapters=True)
TES,EP=illumina_branch("scRCAT-seq TES library","3-prime cDNA boundary",architecture="nextera",cell=16,umi=8,inferred_adapters=True)
FINAL_LIBRARIES=(("TSS-end library",TSS,TP,"Tag-adjacent fragments identify the template-switched transcript 5′ end.","Nextera-compatible primers read the TSS branch."),("TES-end library",TES,EP,"Poly(A)-adjacent fragments identify the transcript 3′ end.","Nextera-compatible primers read the TES branch."))
def oligos(): return ["<b>RT primer architecture</b> — 5′-AAGCAGTGGTATCAACGCAGAGTNNNNNNNN[16-bp cell barcode]T₃₀VN-3′"]
def sections(): return [
 ("Tag both boundaries during reverse transcription",[Row(chunks=[("5′ tag—UMI—cell barcode ** full-length cDNA ** poly(A)-side tag","cbc",False)])],"The RT primer and template switch place related tags at the cDNA ends."),
 ("Amplify full-length cDNA",[Row(chunks=[("tagged cDNA → modified Smart-seq2 PCR","tso",False)])],"Plate position and the 16-base cell identifier preserve cell identity."),
 ("Tagment amplified cDNA",[Row(chunks=[("end-tagged cDNA ** internal Tn5 cuts","me",False)])],"Tagmentation creates short fragments while retaining end-adjacent molecules."),
 ("INFERRED — complete the selected end libraries",[Row(chunks=[("tag primer + Tn5 primer → TSS and TES libraries","r2",False)])],"Targeted selection is source-established; the compact outer adapter completion is inferred."),]


def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    cell=g.input("cell","single-cell poly(A) RNA",MaterialKind.BIOLOGICAL,primary=True)
    rt_primer=g.input("rt_primer","source-exact barcoded UMI oligo(dT) RT primer",MaterialKind.OLIGO,state=oligos()[0])
    rt=g.input("reverse_transcriptase","template-switching reverse transcriptase",MaterialKind.CATALYST)
    first_strand=g.product("first_strand","barcoded UMI-tagged first-strand cDNA",MaterialKind.MOLECULAR,state=sections()[0][1])
    g.react("reverse_transcribe","Prime the transcript 3-prime boundary and reverse-transcribe",StepKind.REVERSE_TRANSCRIPTION,((cell,InputRole.TEMPLATE),(rt_primer,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(first_strand,))
    tso=g.input("tso","scRCAT-seq transcript 5-prime template-switch oligo",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (complete TSO not transcribed in compact module)")
    tagged=g.product("tagged","full-length cDNA tagged at both transcript boundaries",MaterialKind.MOLECULAR,state=sections()[0][1],inferred=True,source=SOURCE+" (complete TSO not transcribed in compact module)")
    g.react("template_switch","Tag the transcript 5-prime boundary by template switching",StepKind.TEMPLATE_SWITCHING,((first_strand,InputRole.TEMPLATE),(tso,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(tagged,))
    preamp_f=g.input("preamp_forward","modified Smart-seq2 preamplification forward primer",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (complete primer not transcribed in compact module)")
    preamp_r=g.input("preamp_reverse","modified Smart-seq2 preamplification reverse primer",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (complete primer not transcribed in compact module)")
    preamp_pol=g.input("preamp_polymerase","modified Smart-seq2 PCR polymerase",MaterialKind.CATALYST)
    amplified=g.product("amplified","amplified end-tagged full-length cDNA",MaterialKind.MOLECULAR,state=sections()[1][1],inferred=True,source=SOURCE+" (complete primer pair not transcribed in compact module)")
    g.react("preamp","Amplify full-length end-tagged cDNA",StepKind.PCR,((tagged,InputRole.TEMPLATE),(preamp_f,InputRole.STRUCTURAL),(preamp_r,InputRole.STRUCTURAL),(preamp_pol,InputRole.CATALYST)),(amplified,),geometry=PCRGeometry.unavailable(TSS,reason="The complete modified Smart-seq2 preamplification primer geometry has not been transcribed in this compact module."))
    transposome=g.input("transposome","loaded Tn5 adapters",MaterialKind.ADAPTER,inferred=True,source=SOURCE+" (complete loaded-adapter oligos not transcribed)")
    tn5=g.input("tn5","Tn5 transposase",MaterialKind.CATALYST)
    tagmented=g.product("tagmented","tagmented end-tagged cDNA",MaterialKind.MOLECULAR,state=sections()[2][1],inferred=True,source=SOURCE+" (complete loaded-adapter oligos not transcribed)")
    tagment_waste=g.product("tagment_waste","internal fragments lacking a selected transcript boundary",MaterialKind.WASTE)
    g.react("tagment","Tagment amplified full-length cDNA",StepKind.TAGMENTATION,((amplified,InputRole.TEMPLATE),(transposome,InputRole.STRUCTURAL),(tn5,InputRole.CATALYST)),(tagmented,tagment_waste));g.discard(tagment_waste,"does not retain a transcript-boundary tag for the selected libraries")
    tss_template=g.product("tss_template","TSS-tag-adjacent tagmented molecules",MaterialKind.MOLECULAR,state=TSS,inferred=True,source=SOURCE+" (outer adapter completion inferred)")
    tes_template=g.product("tes_template","TES-tag-adjacent tagmented molecules",MaterialKind.MOLECULAR,state=TES,inferred=True,source=SOURCE+" (outer adapter completion inferred)")
    g.react("select_ends","Divide tag-adjacent molecules for separate TSS and TES enrichment",StepKind.SPLIT,((tagmented,InputRole.TEMPLATE),),(tss_template,tes_template),changes_molecule=False)
    reviewed=["cell","rt_primer","reverse_transcriptase","reverse_transcribe","tso","template_switch","preamp_forward","preamp_reverse","preamp_polymerase","preamp","transposome","tn5","tagment","select_ends"]
    finish_pcr_branch(g,"tss",tss_template,TSS,TP,reviewed,label="scRCAT-seq TSS end",inferred_primers=True,limitation="The complete tag-selective TSS and Nextera index-primer pair has not been transcribed in this compact module.")
    finish_pcr_branch(g,"tes",tes_template,TES,EP,reviewed,label="scRCAT-seq TES end",inferred_primers=True,limitation="The complete tag-selective TES and Nextera index-primer pair has not been transcribed in this compact module.")
    g.declare_reviewed_scope(reviewed,SOURCE)
    return g.require_publishable()
