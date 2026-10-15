"""scCAT-seq physical RNA/nucleus separation."""
import nextera as nx
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "smart-seq2__10.1038+nmeth.2639" / "tools"))
import smartseq2 as ss2
from chemdraw import Construct, Row, Segment
from multibranch import illumina_branch
from multimodal_spatial import finish_pcr_branch
from pcr import PCRGeometry, PCRPrimer
from protocol_extensions import analyte_split
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind
TITLE="scCAT-seq — transcriptome and chromatin accessibility"
NOTES="01_sccat-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/s41467-018-08205-7">Liu et al., <i>Nature Communications</i> (2019)</a>.'
SUMMARY="Gentle lysis separates cytoplasmic RNA from the nucleus of one cell; Smart-seq2 and scATAC-seq then process the matched fractions."
CAVEAT=("The defining paper publishes the complete SMART-seq2 RT, template-switch and "
        "cDNA-preamp oligos and both ATAC preamplification primers. It names, but does "
        "not sequence, the final universal and barcode primers used before BGISEQ-500 "
        "conversion; those endpoint primer geometries therefore remain unavailable.")
RNA,RP=illumina_branch("scCAT-seq RNA library","Smart-seq2 cDNA",architecture="nextera")
ATAC,AP=illumina_branch("scCAT-seq ATAC library","accessible chromatin fragment",architecture="nextera")
FINAL_LIBRARIES=(("Transcriptome library",RNA,RP,"Smart-seq2 cDNA after Nextera library construction.","Nextera-compatible primers read the RNA branch."),("Chromatin-accessibility library",ATAC,AP,"The matched nucleus yields a tagmented accessibility library.","Nextera-compatible primers read the ATAC branch."))
ATAC_PREAMP_FORWARD=nx.READ1_PRIMER
ATAC_PREAMP_REVERSE=nx.READ2_PRIMER

def oligos(): return [
    "<b>Modified oligo-dT primer</b> — 5′-"+ss2.OLIGO_DT+"-3′",
    # The defining paper's Methods print this oligo in full, ACAT linker included, so the
    # shared SMART-seq2 model supplies the bases rather than a second transcription.
    "<b>Template-switching oligo</b> — 5′-"+ss2.TSO+"-3′",
    "<b>cDNA PCR primer</b> — 5′-"+ss2.ISPCR+"-3′",
    "<b>transposase adapter 1 / ATAC preamp forward</b> — 5′-"+ATAC_PREAMP_FORWARD+"-3′",
    "<b>transposase adapter 2 / ATAC preamp reverse</b> — 5′-"+ATAC_PREAMP_REVERSE+"-3′",
]
def sections(): return [
 ("Lyse one cell without disrupting its nucleus",[Row(chunks=[("cell → cytoplasmic RNA + intact nucleus","r1",False)])],"The nucleus is physically retained while cytoplasmic RNA enters the supernatant."),
 ("Separate the matched analytes",analyte_split("cytoplasmic RNA → Smart-seq2","nucleus → scATAC-seq"),"The two branches retain common cell provenance, not a covalent barcode."),
 ("Build the RNA branch",[Row(chunks=[("oligo(dT) RT → template switching → cDNA PCR → tagmentation","tso",False)])],"Smart-seq2 produces the transcriptome library."),
 ("Build the chromatin branch",[Row(chunks=[("nuclear chromatin → Tn5 tagmentation → indexed PCR","me",False)])],"The intact matched nucleus produces the accessibility library."),]


def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    cell=g.input("cell","one intact cell",MaterialKind.BIOLOGICAL,primary=True)
    cytoplasm=g.product("cytoplasm","cytoplasmic poly(A) RNA fraction",MaterialKind.MOLECULAR,state=analyte_split("cytoplasmic RNA → Smart-seq2","nucleus → scATAC-seq")[:1])
    nucleus=g.product("nucleus","matched intact nucleus",MaterialKind.MOLECULAR,state=analyte_split("cytoplasmic RNA → Smart-seq2","nucleus → scATAC-seq")[1:])
    g.react("separate","Gently lyse the cell and separate cytoplasm from its nucleus",StepKind.SPLIT,((cell,InputRole.TEMPLATE),),(cytoplasm,nucleus),changes_molecule=False)

    oligo_dt=g.input("oligo_dt","published Smart-seq2 anchored oligo(dT) primer",MaterialKind.OLIGO,
                     state=Construct((Segment("anchored oligo(dT)",ss2.OLIGO_DT,
                                              placeholder=True),)))
    rt=g.input("reverse_transcriptase","template-switching reverse transcriptase",MaterialKind.CATALYST)
    first_strand=g.product("first_strand","oligo(dT)-primed cytoplasmic first-strand cDNA",MaterialKind.MOLECULAR)
    g.react("reverse_transcribe","Reverse-transcribe cytoplasmic poly(A) RNA",StepKind.REVERSE_TRANSCRIPTION,((cytoplasm,InputRole.TEMPLATE),(oligo_dt,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(first_strand,))
    tso=g.input("tso","published Smart-seq2 template-switch oligo",MaterialKind.OLIGO,
                state=Construct((Segment("scCAT template-switch oligo",
                                         ss2.ISPCR+ss2.TSO_LINKER+"GGG",
                                         placeholder=True),)))
    switched=g.product("switched","template-switched full-length cDNA",MaterialKind.MOLECULAR,state=ss2.preamp_template(linker_inferred=False))
    g.react("template_switch","Template-switch at the cDNA 5-prime end",StepKind.TEMPLATE_SWITCHING,((first_strand,InputRole.TEMPLATE),(tso,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(switched,))
    preamp_f=g.input("preamp_forward","published SMART/ISPCR primer, forward-facing copy",MaterialKind.OLIGO,
                     state=Construct((Segment("ISPCR",ss2.ISPCR),)))
    preamp_r=g.input("preamp_reverse","published SMART/ISPCR primer, reverse-facing copy",MaterialKind.OLIGO,
                     state=Construct((Segment("ISPCR",ss2.ISPCR),)))
    preamp_pol=g.input("preamp_polymerase","Smart-seq2 cDNA PCR polymerase",MaterialKind.CATALYST)
    amplified=g.product("amplified_cdna","amplified full-length transcript cDNA",MaterialKind.MOLECULAR,state=ss2.preamp_template(linker_inferred=False))
    g.react("preamp","Preamplify the Smart-seq2 cDNA",StepKind.PCR,((switched,InputRole.TEMPLATE),(preamp_f,InputRole.STRUCTURAL),(preamp_r,InputRole.STRUCTURAL),(preamp_pol,InputRole.CATALYST)),(amplified,),geometry=ss2.ispcr_geometry(switched.state))

    rna_adapter=g.input("rna_transposome","Nextera transposome for amplified cDNA",MaterialKind.ADAPTER,inferred=True,source=SOURCE+" (loaded-adapter oligos not transcribed on this page)")
    rna_tn5=g.input("rna_tn5","RNA-branch Tn5 transposase",MaterialKind.CATALYST)
    rna_template=g.product("rna_template","tagmented Smart-seq2 cDNA",MaterialKind.MOLECULAR,state=RNA,inferred=True,source=SOURCE+" (loaded-adapter oligos not transcribed on this page)")
    rna_waste=g.product("rna_tagment_waste","non-library cDNA tagmentation products",MaterialKind.WASTE)
    g.react("tagment_rna","Tagment amplified transcript cDNA",StepKind.TAGMENTATION,((amplified,InputRole.TEMPLATE),(rna_adapter,InputRole.STRUCTURAL),(rna_tn5,InputRole.CATALYST)),(rna_template,rna_waste));g.discard(rna_waste,"does not retain a library-eligible heterologous Nextera fragment")

    atac_adapter=g.input("atac_transposome","scATAC-seq transposome carrying the published adapters",MaterialKind.ADAPTER)
    atac_tn5=g.input("atac_tn5","chromatin Tn5 transposase",MaterialKind.CATALYST)
    atac_template=g.product("atac_template","tagmented accessible chromatin from the matched nucleus",MaterialKind.MOLECULAR,state=ATAC,inferred=True,source=SOURCE+" (loaded-adapter oligos not transcribed on this page)")
    atac_waste=g.product("atac_tagment_waste","non-library chromatin tagmentation products",MaterialKind.WASTE)
    g.react("tagment_atac","Tagment accessible chromatin in the matched nucleus",StepKind.TAGMENTATION,((nucleus,InputRole.TEMPLATE),(atac_adapter,InputRole.STRUCTURAL),(atac_tn5,InputRole.CATALYST)),(atac_template,atac_waste));g.discard(atac_waste,"does not retain a library-eligible heterologous scATAC fragment")

    reviewed=["cell","separate","oligo_dt","reverse_transcriptase","reverse_transcribe","tso","template_switch","preamp_forward","preamp_reverse","preamp_polymerase","preamp","rna_transposome","rna_tn5","tagment_rna","atac_transposome","atac_tn5","tagment_atac"]
    finish_pcr_branch(g,"rna",rna_template,RNA,RP,reviewed,label="scCAT-seq transcriptome",inferred_primers=True,limitation="The exact Nextera index-primer pair for the RNA branch has not yet been transcribed on this page.")
    finish_pcr_branch(g,"atac",atac_template,ATAC,AP,reviewed,label="scCAT-seq chromatin accessibility",inferred_primers=True,limitation="The exact scATAC index-primer pair has not yet been transcribed on this page.")
    g.declare_reviewed_scope(reviewed,SOURCE)
    return g.require_publishable()
