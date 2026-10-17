"""scCAT-seq physical RNA/nucleus separation."""
import nextera as nx
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "smart-seq2__10.1038+nmeth.2639" / "tools"))
import smartseq2 as ss2
from chemdraw import Construct, Row, Segment
import seqprimers as sp
from pcr import PCRGeometry, PCRPrimer
from protocol_extensions import analyte_split
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
TITLE="scCAT-seq — transcriptome and chromatin accessibility"
NOTES="01_sccat-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/s41467-018-08205-7">Liu et al., <i>Nature Communications</i> (2019)</a>.'
SUMMARY="Gentle lysis separates cytoplasmic RNA from the nucleus of one cell; Smart-seq2 and scATAC-seq then process the matched fractions."
CAVEAT=("The defining paper publishes the complete SMART-seq2 RT, template-switch and "
        "cDNA-preamp oligos and both ATAC preamplification primers. It names, but does "
        "not sequence, the final ATAC universal and barcode primers, refers the RNA-branch "
        "PCR to Smart-seq2, and then converts all libraries into BGISEQ-500 ssDNA circles "
        "and DNA nanoballs for paired-end 50-bp reads. The endpoint primers, any platform "
        "adapter tails and the circle/DNB oligos are therefore not drawn: each final panel "
        "shows only the Tn5-adapter-flanked core that the published preamp primers (ATAC) "
        "or the shared Tn5 adapter design (RNA, inferred) imply, not an Illumina library.")
SOURCE_METHODS=SOURCE+" Methods, scCAT-seq library preparation and sequencing (PMC6349937 JATS paragraphs 219–229)"

def _core(name, *, inferred):
    """Tn5-adapter-flanked fragment after 72 °C gap fill; insert is a placeholder."""
    return Construct((Segment("transposase adapter 1 (S5)",nx.S5,"s5",inferred=inferred),
                      Segment("mosaic end",nx.ME,"me",inferred=inferred),
                      Segment("accessible or cDNA insert","X"*34,placeholder=True),
                      Segment("mosaic end reverse complement",nx.ME_RC,"me",inferred=inferred),
                      Segment("transposase adapter 2 reverse complement (S7)",nx.S7_RC,"s7",inferred=inferred)),name=name)

ATAC=_core("scCAT-seq ATAC adapter-flanked core (endpoint tails and BGISEQ circle not drawn)",inferred=False)
RNA=_core("scCAT-seq RNA adapter-flanked core (endpoint tails and BGISEQ circle not drawn)",inferred=True)
ENDPOINT_LIMITATION=("All scCAT-seq libraries were converted into BGISEQ-500 ssDNA circles and "
    "DNA nanoballs and sequenced paired-end 50 bp (Methods). The circle/DNB oligos, the "
    "endpoint primers and the BGISEQ sequencing-primer sites are not given in the defining "
    "paper, so no primer landing is computed; Illumina Nextera read primers are not claimed. "
    "This is a model limitation pending the named BGISEQ and Smart-seq2 originals, not an "
    "accepted publication-absence claim.")
_PE50=(sp.read_spec("Sequencing output",50,"BGISEQ-500 read 1 of a 50-bp pair; read entry point on the converted circle not modelled","include"),
       sp.read_spec("Sequencing output",50,"BGISEQ-500 read 2 of a 50-bp pair; read entry point on the converted circle not modelled","include"))
FINAL_LIBRARIES=(("Transcriptome library",RNA,(),"Inferred Tn5-adapter-flanked tagmented cDNA core. The Smart-seq2 PCR primers, any BGISEQ adapter tails and the ssDNA circle that is actually sequenced are not drawn.","BGISEQ-500 paired-end 50 bp."),
                 ("Chromatin-accessibility library",ATAC,(),"Tn5-adapter-flanked accessible-chromatin core implied by the published preamp primers. The unsequenced universal/barcode primer extensions and the BGISEQ ssDNA circle are not drawn.","BGISEQ-500 paired-end 50 bp."))
READ_SPECS={"Transcriptome library":_PE50,"Chromatin-accessibility library":_PE50}
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
    "<b>ATAC universal and barcode primers</b> — named in Methods without bases",
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

    rna_adapter=g.input("rna_transposome","in-house Tn5 transposome for 2 ng amplified cDNA (loaded adapters not stated for this branch)",MaterialKind.ADAPTER,inferred=True,source=SOURCE_METHODS)
    rna_tn5=g.input("rna_tn5","in-house Tn5 transposase, 0.3 µl in 1x TAG buffer",MaterialKind.CATALYST)
    rna_template=g.product("rna_template","tagmented Smart-seq2 cDNA after SDS release of Tn5",MaterialKind.MOLECULAR,state=RNA,inferred=True,source=SOURCE_METHODS+" (adapter identity on the RNA branch inferred)")
    rna_waste=g.product("rna_tagment_waste","non-library cDNA tagmentation products",MaterialKind.WASTE)
    g.react("tagment_rna","Tagment 2 ng cDNA at 55 °C for 10 min, then release Tn5 with 0.1% SDS",StepKind.TAGMENTATION,((amplified,InputRole.TEMPLATE),(rna_adapter,InputRole.STRUCTURAL),(rna_tn5,InputRole.CATALYST)),(rna_template,rna_waste));g.discard(rna_waste,"does not retain a library-eligible heterologous adapter pair")

    atac_adapter=g.input("atac_transposome","in-house Tn5 transposome whose adapters match the published preamp primers",MaterialKind.ADAPTER)
    atac_tn5=g.input("atac_tn5","in-house Tn5 transposase, 0.3 µl in 1x TAG buffer",MaterialKind.CATALYST)
    atac_tagmented=g.product("atac_tagmented","tagmented accessible chromatin from the matched nucleus",MaterialKind.MOLECULAR)
    atac_waste=g.product("atac_tagment_waste","non-library chromatin tagmentation products",MaterialKind.WASTE)
    g.react("tagment_atac","Tagment the matched nucleus 15 min at 37 °C; stop with EDTA at 50 °C",StepKind.TAGMENTATION,((nucleus,InputRole.TEMPLATE),(atac_adapter,InputRole.STRUCTURAL),(atac_tn5,InputRole.CATALYST)),(atac_tagmented,atac_waste));g.discard(atac_waste,"does not retain a library-eligible heterologous adapter pair")
    carrier=g.input("carrier_plasmid","30 ng carrier plasmid DNA added with RLT Plus buffer",MaterialKind.MOLECULAR)
    beads18=g.input("xp_beads_18","AMPure XP beads, 1:1.8 on KingFisher Flex",MaterialKind.SELECTION)
    atac_clean=g.product("atac_clean","purified transposed DNA with carrier plasmid",MaterialKind.MOLECULAR)
    atac_clean_waste=g.product("atac_clean_waste","bead supernatant: protein, buffer and unbound material",MaterialKind.WASTE)
    g.react("purify_atac","Add carrier plasmid and RLT Plus, lyse 15 min at 37 °C, purify with 1.8x XP beads",StepKind.CLEANUP,((atac_tagmented,InputRole.TEMPLATE),(carrier,InputRole.SPIKE_IN),(beads18,InputRole.SELECTION)),(atac_clean,atac_clean_waste));g.discard(atac_clean_waste,"bead supernatant removed")
    ad1=g.input("atac_adapter1","published transposase adapter 1 primer",MaterialKind.OLIGO,state=Construct((Segment("adapter 1",ATAC_PREAMP_FORWARD),)))
    ad2=g.input("atac_adapter2","published transposase adapter 2 primer",MaterialKind.OLIGO,state=Construct((Segment("adapter 2",ATAC_PREAMP_REVERSE),)))
    neb=g.input("atac_preamp_polymerase","NEBNext High-Fidelity 2x PCR Master Mix",MaterialKind.CATALYST)
    atac_pre=g.product("atac_preamplified","8-cycle preamplified adapter-flanked transposed DNA",MaterialKind.MOLECULAR,state=ATAC)
    g.react("atac_preamp","Gap-fill at 72 °C, then preamplify 8 cycles with adapter 1/2 primers and purify (1x XP)",StepKind.PCR,((atac_clean,InputRole.TEMPLATE),(ad1,InputRole.STRUCTURAL),(ad2,InputRole.STRUCTURAL),(neb,InputRole.CATALYST)),(atac_pre,),
            geometry=PCRGeometry.exact(ATAC,PCRPrimer("transposase adapter 1",ATAC_PREAMP_FORWARD,anneal_nt=len(ATAC_PREAMP_FORWARD)),PCRPrimer("transposase adapter 2",ATAC_PREAMP_REVERSE,anneal_nt=len(ATAC_PREAMP_REVERSE)),
            statement="The two published 33/34-nt primers are the full S5+ME and S7+ME adapter strands, so after 72 °C gap fill each anneals across its whole length to one end of a heterologously tagged fragment."))

    reviewed=["cell","separate","oligo_dt","reverse_transcriptase","reverse_transcribe","tso","template_switch","preamp_forward","preamp_reverse","preamp_polymerase","preamp","rna_transposome","rna_tn5","tagment_rna","atac_transposome","atac_tn5","tagment_atac","carrier_plasmid","xp_beads_18","purify_atac","atac_adapter1","atac_adapter2","atac_preamp_polymerase","atac_preamp","atac_preamplified"]
    endpoint_src=SOURCE_METHODS+" (endpoint primer sequence not published in the defining paper)"
    def endpoint(prefix,template,state,label,pcr_label,pcr_reason,polymerase_label):
        f=g.input(prefix+"_forward",label+" forward endpoint primer",MaterialKind.OLIGO,inferred=True,source=endpoint_src)
        r=g.input(prefix+"_reverse",label+" reverse endpoint primer",MaterialKind.OLIGO,inferred=True,source=endpoint_src)
        pol=g.input(prefix+"_polymerase",polymerase_label,MaterialKind.CATALYST)
        linear=g.product(prefix+"_linear",label+" linear double-stranded library before BGISEQ conversion",MaterialKind.MOLECULAR,state=state,inferred=True,source=endpoint_src)
        g.react(prefix+"_pcr",pcr_label,StepKind.PCR,((template,InputRole.TEMPLATE),(f,InputRole.STRUCTURAL),(r,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(linear,),geometry=PCRGeometry.unavailable(state,reason=pcr_reason))
        library=g.product(prefix+"_library",label+" BGISEQ-500 DNA-nanoball library",MaterialKind.LIBRARY,state=state,inferred=True,source=SOURCE_METHODS)
        g.react(prefix+"_bgiseq","Pool 2 pmol, make ssDNA circles and roll them into DNA nanoballs (BGISEQ-500)",StepKind.MAGIC,((linear,InputRole.TEMPLATE),),(library,),
                reason="Methods report ssDNA-circle formation, rolling-circle DNB generation and BGISEQ-500 PE50 sequencing, citing the BGISEQ-500 reference-dataset paper; circle/DNB oligos and the platform adapter junction are not given there and are not modelled. Model gap, not a publication-absence claim.")
        g.sequence(library,SequencingPlan(library=state,read_lengths={},limitation=ENDPOINT_LIMITATION))
        reviewed.extend((f.id,r.id,pol.id,prefix+"_pcr",linear.id,prefix+"_bgiseq",library.id))
    endpoint("rna",rna_template,RNA,"scCAT-seq transcriptome","PCR-amplify tagmented cDNA as in Smart-seq2 (primers not printed by scCAT-seq)",
             "The defining paper defers this PCR to the Smart-seq2 method and does not print its primers or any BGISEQ tails.","RNA-branch library PCR polymerase (not named)")
    sel=g.input("xp_size_select","AMPure XP size selection for 150–700 bp fragments",MaterialKind.SELECTION)
    atac_amp=g.product("atac_amplified","10–16-cycle universal/barcode-primer PCR product",MaterialKind.MOLECULAR,state=ATAC,inferred=True,source=endpoint_src)
    u=g.input("atac_universal","ATAC universal primer (bases not published)",MaterialKind.OLIGO,inferred=True,source=endpoint_src)
    b=g.input("atac_barcode","ATAC barcode primer (bases not published)",MaterialKind.OLIGO,inferred=True,source=endpoint_src)
    neb2=g.input("atac_final_polymerase","NEBNext High-Fidelity 2x PCR Master Mix",MaterialKind.CATALYST)
    g.react("atac_final_pcr","Amplify 10–16 more cycles with 20 µM universal and barcode primers (15 for one cell)",StepKind.PCR,((atac_pre,InputRole.TEMPLATE),(u,InputRole.STRUCTURAL),(b,InputRole.STRUCTURAL),(neb2,InputRole.CATALYST)),(atac_amp,),
            geometry=PCRGeometry.unavailable(ATAC,reason="The universal and barcode primers are named without bases, so their landing sites and added tails cannot be computed."))
    atac_linear=g.product("atac_linear","size-selected ATAC library before BGISEQ conversion",MaterialKind.MOLECULAR,state=ATAC,inferred=True,source=endpoint_src)
    atac_sel_waste=g.product("atac_size_waste","fragments outside 150–700 bp",MaterialKind.WASTE)
    g.react("atac_size_select","Size-select 150–700 bp with XP beads; elute in TE",StepKind.SELECTION,((atac_amp,InputRole.TEMPLATE),(sel,InputRole.SELECTION)),(atac_linear,atac_sel_waste));g.discard(atac_sel_waste,"removed by bead size selection")
    atac_lib=g.product("atac_library","scCAT-seq chromatin accessibility BGISEQ-500 DNA-nanoball library",MaterialKind.LIBRARY,state=ATAC,inferred=True,source=SOURCE_METHODS)
    g.react("atac_bgiseq","Pool 2 pmol, make ssDNA circles and roll them into DNA nanoballs (BGISEQ-500)",StepKind.MAGIC,((atac_linear,InputRole.TEMPLATE),),(atac_lib,),
            reason="Methods report ssDNA-circle formation, rolling-circle DNB generation and BGISEQ-500 PE50 sequencing, citing the BGISEQ-500 reference-dataset paper; circle/DNB oligos and the platform adapter junction are not given there and are not modelled. Model gap, not a publication-absence claim.")
    g.sequence(atac_lib,SequencingPlan(library=ATAC,limitation=ENDPOINT_LIMITATION))
    reviewed.extend(("xp_size_select","atac_amplified","atac_universal","atac_barcode","atac_final_polymerase","atac_final_pcr","atac_linear","atac_size_select","atac_library","atac_bgiseq"))
    g.declare_reviewed_scope(reviewed,SOURCE)
    return g.require_publishable()
