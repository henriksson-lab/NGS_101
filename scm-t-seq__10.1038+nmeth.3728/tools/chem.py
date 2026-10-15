"""scM&T-seq physical RNA/DNA separation."""
from pathlib import Path
from dataclasses import replace
import sys
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"scbs-seq__10.1038+nmeth.3035"/"tools"))
from chemdraw import Construct, Row
from batch_ngs import nextera_library, seg
import illumina as il
from multimodal_spatial import finish_pcr_branch
from pcr import PCRGeometry, PCRPrimer
from protocol_extensions import analyte_split
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind
import rt
import scbs
import smartseq2 as s2
TITLE="scM&T-seq — transcriptome and DNA methylome"
NOTES="01_scm-t-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/nmeth.3728">Angermueller et al., <i>Nature Methods</i> (2016)</a>.'
SUMMARY="Biotinylated oligo(dT) beads remove a cell's mRNA for Smart-seq2 while the genomic-DNA supernatant enters scBS-seq."
CAVEAT=("The defining method establishes bead-based physical separation and explicitly "
        "inherits Smart-seq2/Nextera XT and scBS-seq. Their authoritative oligos are "
        "reused here, so both branch-specific PCR geometries are exact. One inherited "
        "sequence is not exact, though: Angermueller et al. do not print a "
        "template-switch oligo at all &mdash; the RNA half is delegated to the "
        "G&amp;T-seq protocol &mdash; and "
        "Smart-seq2 has two published TSO linkers. Picelli et al. 2013 Supplementary "
        "Table 2 gives the chosen <code>rGrG+G</code> oligo as the 23-nt handle + "
        "<code>AC</code> + <code>rGrG+G</code>, while the <i>Nature Protocols</i> form "
        "carries <code>ACAT</code>. The longer form is modelled, and the whole TSO "
        "linker is therefore drawn dotted as an inherited inference.")


def _tso_linker_inferred(con: Construct) -> Construct:
    """Redraw the inherited Smart-seq2 TSO linker as an inference, not a published base.

    scM&T-seq prints no TSO, so the choice between the two published Smart-seq2 linkers
    (``AC`` and ``ACAT``) is ours.  Marking the segment here, rather than in every caller,
    keeps the uncertainty a property of the state the page renders.
    """
    con.get("TSO linker")      # a lib rename must fail loudly, not silently un-mark it
    return Construct([replace(s,inferred=True) if s.name=="TSO linker" else s
                      for s in con],name=con.name)


RNA,RP=nextera_library([seg("Smart-seq2 cDNA insert","X"*36,placeholder=True)],"scM&T-seq RNA")
DNA=Construct([replace(s,name=s.name.replace("N9","N6"),top="N"*6)
               if s.name in ("oligo1 N9","oligo2 N9") else s
               for s in scbs.final_library()],name="scM&T-seq DNA methylome")
DP=tuple(scbs.SEQ_PRIMERS)
RNA_PREAMP=_tso_linker_inferred(s2.preamp_template())
FINAL_LIBRARIES=(("RNA library",RNA,RP,"Smart-seq2 cDNA after Nextera XT.","Nextera primers declared for the RNA branch."),("DNA methylome library",DNA,DP,"Bisulfite-converted genomic fragments after scBS-seq.","Illumina primers declared for the DNA branch."))

def oligos():
    return [f"Smart-seq2 oligo-dT30VN: 5′-{s2.OLIGO_DT}-3′",
            f"Smart-seq2 TSO (linker inherited, not printed by this paper): "
            f"5′-{rt.SMART_HANDLE}<inf>{s2.TSO_LINKER}</inf>{rt.TSO_G_TAIL_LNA}-3′",
            f"Smart-seq2 ISPCR primer: 5′-{s2.ISPCR}-3′",
            f"scBS-seq PE1: 5′-{scbs.PE1}-3′",
            f"scM&T oligo1 N6: 5′-{scbs.OLIGO1_HANDLE}{'N'*6}-3′",
            f"scM&T oligo2 N6: 5′-{scbs.OLIGO2_HANDLE}{'N'*6}-3′",
            "scBS-seq indexed iPCRTag: P7–i7–GATCGGAAGAGCACACGTCTGAACTCCAGTCAC"]
def sections():
    return [("Capture mRNA on biotinylated oligo(dT) beads",[Row(chunks=[("mRNA poly(A) || biotin—oligo(dT) ● bead", "r1", False)])],"The bead-bound fraction is physically removed from genomic DNA."),("Separate the analytes",analyte_split("bead-bound RNA → Smart-seq2","supernatant DNA → scBS-seq"),"Wash fractions are retained with the DNA to reduce loss."),("Convert and amplify each branch",[Row(chunks=[("RNA: RT → template switch → cDNA PCR → Nextera", "tso", False)]),Row(chunks=[("DNA: bisulfite → random priming → adapter PCR", "r2", False)])],"The branch chemistries remain distinct after separation."),]


def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    lysate=g.input("lysate","one-cell lysate containing genomic DNA and poly(A) RNA",MaterialKind.BIOLOGICAL,state=sections()[0][1],primary=True)
    beads=g.input("oligodt_beads","biotinylated oligo(dT) magnetic beads",MaterialKind.SELECTION)
    rna=g.product("rna_fraction","bead-bound poly(A) RNA",MaterialKind.MOLECULAR,state=sections()[1][1][:1])
    dna=g.product("dna_fraction","genomic-DNA supernatant plus retained washes",MaterialKind.MOLECULAR,state=sections()[1][1][1:])
    g.react("separate","Capture mRNA and separate it from genomic DNA",StepKind.AFFINITY_CAPTURE,((lysate,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),(rna,dna),changes_molecule=False)

    rt_primer=g.input("rt_primer","Smart-seq2 anchored oligo(dT) primer",MaterialKind.OLIGO,state=(seg("SMART handle",rt.SMART_HANDLE,"tso"),seg("oligo-dT linker",s2.DT_LINKER),seg("dT30VN","T"*30+"VN",placeholder=True)))
    rt_enzyme=g.input("rt","template-switching reverse transcriptase",MaterialKind.CATALYST)
    copied=g.product("copied","first-strand Smart-seq2 cDNA",MaterialKind.MOLECULAR,state=RNA_PREAMP)
    g.react("reverse_transcribe","Reverse-transcribe bead-captured RNA",StepKind.REVERSE_TRANSCRIPTION,((rna,InputRole.TEMPLATE),(rt_primer,InputRole.STRUCTURAL),(rt_enzyme,InputRole.CATALYST)),(copied,))
    tso=g.input("tso","Smart-seq2 template-switch oligo",MaterialKind.OLIGO,state=(seg("SMART handle",rt.SMART_HANDLE,"tso"),seg("TSO linker",s2.TSO_LINKER,inferred=True),seg("rGrG+G","GGG","tso")))
    switched=g.product("switched","template-switched full-length cDNA",MaterialKind.MOLECULAR,state=RNA_PREAMP)
    g.react("template_switch","INFERRED — Template-switch the nascent cDNA (the TSO linker is inherited from Smart-seq2; this paper prints no TSO)",StepKind.TEMPLATE_SWITCHING,((copied,InputRole.TEMPLATE),(tso,InputRole.STRUCTURAL),(rt_enzyme,InputRole.CATALYST)),(switched,))
    ispcr_f=g.input("ispcr_forward","Smart-seq2 ISPCR primer at the 5′ handle",MaterialKind.OLIGO,state=(seg("ISPCR",s2.ISPCR,"tso"),))
    ispcr_r=g.input("ispcr_reverse","same Smart-seq2 ISPCR species at the 3′ handle",MaterialKind.OLIGO,state=(seg("ISPCR",s2.ISPCR,"tso"),))
    prepol=g.input("preamp_polymerase","Smart-seq2 cDNA PCR polymerase",MaterialKind.CATALYST)
    preamp=g.product("preamp","amplified full-length Smart-seq2 cDNA",MaterialKind.MOLECULAR,state=RNA_PREAMP)
    g.react("preamp_rna","Amplify full-length cDNA",StepKind.PCR,((switched,InputRole.TEMPLATE),(ispcr_f,InputRole.STRUCTURAL),(ispcr_r,InputRole.STRUCTURAL),(prepol,InputRole.CATALYST)),(preamp,),geometry=s2.ispcr_geometry(RNA_PREAMP))
    transposome=g.input("rna_transposome","Nextera XT transposome",MaterialKind.ADAPTER)
    tn5=g.input("rna_tn5","Tn5 transposase",MaterialKind.CATALYST)
    rna_tagged=g.product("rna_tagged","Nextera-tagged Smart-seq2 cDNA",MaterialKind.MOLECULAR,state=RNA)
    twaste=g.product("rna_tn5_waste","non-library cDNA tagmentation products",MaterialKind.WASTE)
    g.react("tagment_rna","Tagment amplified cDNA with Nextera XT",StepKind.TAGMENTATION,((preamp,InputRole.TEMPLATE),(transposome,InputRole.STRUCTURAL),(tn5,InputRole.CATALYST)),(rna_tagged,twaste));g.discard(twaste,"does not retain an amplifiable unlike-ended cDNA fragment")

    bisulfite=g.input("bisulfite","bisulfite-conversion condition",MaterialKind.CONDITION)
    converted=g.product("converted","bisulfite-converted genomic DNA",MaterialKind.MOLECULAR,state=DNA)
    g.react("bisulfite_convert","Bisulfite-convert the genomic-DNA fraction",StepKind.CHEMICAL_CONVERSION,((dna,InputRole.TEMPLATE),(bisulfite,InputRole.CONDITION)),(converted,))
    random_primer=g.input("random_primer","scM&T oligo1 with a six-base random 3′ end",MaterialKind.OLIGO)
    dna_polymerase=g.input("dna_strand_polymerase","scBS-seq strand-synthesis polymerase",MaterialKind.CATALYST)
    primed=g.product("primed_dna","oligo1-primed bisulfite-converted DNA",MaterialKind.MOLECULAR,state=DNA)
    g.react("prime_dna","Random-prime converted genomic DNA",StepKind.EXTENSION,((converted,InputRole.TEMPLATE),(random_primer,InputRole.STRUCTURAL),(dna_polymerase,InputRole.CATALYST)),(primed,))
    oligo2=g.input("oligo2","scM&T oligo2 with a six-base random 3′ end",MaterialKind.OLIGO)
    dna_adapted=g.product("dna_adapted","oligo1/oligo2-completed scBS-seq template",MaterialKind.MOLECULAR,state=DNA)
    g.react("adapt_dna","Random-prime the complementary strand with scBS-seq oligo2",StepKind.EXTENSION,((primed,InputRole.TEMPLATE),(oligo2,InputRole.STRUCTURAL),(dna_polymerase,InputRole.CATALYST)),(dna_adapted,))

    reviewed=["lysate","oligodt_beads","separate","rt_primer","rt","reverse_transcribe","tso","template_switch","ispcr_forward","ispcr_reverse","preamp_polymerase","preamp_rna","rna_transposome","rna_tn5","tagment_rna","bisulfite","bisulfite_convert","random_primer","dna_strand_polymerase","prime_dna","oligo2","adapt_dna"]
    finish_pcr_branch(g,"rna",rna_tagged,RNA,RP,reviewed,label="scM&T-seq RNA",geometry=s2.nextera_index_geometry(RNA))
    ipcrtag=il.P7+"N"*scbs.IPCRTAG_INDEX_NT+scbs.IPCRTAG_3_ARM
    dna_geometry=PCRGeometry.exact(DNA,PCRPrimer("PE1",scbs.PE1),PCRPrimer("indexed iPCRTag",ipcrtag,anneal_nt=len(scbs.IPCRTAG_3_ARM)),statement="Smallwood et al. print PE1 and the indexed iPCRTag family; their annealing regions face inward on the scBS-seq product.")
    finish_pcr_branch(g,"dna",dna_adapted,DNA,DP,reviewed,label="scM&T-seq DNA methylome",geometry=dna_geometry)
    g.declare_reviewed_scope(reviewed,SOURCE)
    return g.require_publishable()
