"""scM&T-seq physical RNA/DNA separation."""
from pathlib import Path
from dataclasses import replace
import sys
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"scbs-seq__10.1038+nmeth.3035"/"tools"))
sys.path.insert(0,str(ROOT/"smart-seq2__10.1038+nmeth.2639"/"tools"))
from chemdraw import Construct, Row
from batch_ngs import nextera_library, seg
import illumina as il
from multimodal_spatial import finish_pcr_branch
from pcr import PCRGeometry, PCRPrimer
from protocol_extensions import analyte_split
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
import rt
import scbs
import smartseq2 as s2
TITLE="scM&T-seq — transcriptome and DNA methylome"
NOTES="01_scm-t-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/nmeth.3728">Angermueller et al., <i>Nature Methods</i> (2016)</a>.'
SUMMARY="Biotinylated oligo(dT) beads remove a cell's mRNA for Smart-seq2 while the genomic-DNA supernatant enters scBS-seq."
CAVEAT=("The defining method establishes bead-based physical separation and explicitly "
        "inherits Smart-seq2/Nextera XT and scBS-seq. Their authoritative oligos are "
        "reused here. Angermueller et al. do not print a template-switch oligo; their "
        "cited G&amp;T-seq predecessor (Macaulay et al. 2015, Online Methods p. 4) "
        "prints the 23-nt handle + <code>AC</code> + <code>rGrG+G</code>. "
        "Picelli et al. 2013 Supplementary Table 2 also gives <code>AC</code>, "
        "whereas Picelli et al. 2014 <i>Nature Protocols</i> gives <code>ACAT</code>. "
        "The G&amp;T <code>AC</code> form is modelled; its linker remains dotted because "
        "scM&amp;T delegates rather than prints this exact oligo. "
        "Supplementary Fig. 3e reports sequenced DNA fragments made of "
        "first-strand primer concatemers; their exact molecular assembly is not established.")

G_AND_T_TSO_LINKER="AC"
# Angermueller et al. Online Methods "Sample collection & Single-cell sequencing":
# both libraries were sequenced on an Illumina HiSeq2000 with 125-bp paired-end reads.
# Index-read cycles are not stated.
READS_125PE={"Read 1":125,"Read 2":125}

RNA,RP=nextera_library([seg("Smart-seq2 cDNA insert","X"*36,placeholder=True)],"scM&T-seq RNA")
DNA=Construct([replace(s,name=s.name.replace("N9","N6"),top="N"*6)
               if s.name in ("oligo1 N9","oligo2 N9") else s
               for s in scbs.final_library()],name="scM&T-seq DNA methylome")
DP=tuple(scbs.SEQ_PRIMERS)
RNA_PREAMP=Construct([replace(s,top=G_AND_T_TSO_LINKER) if s.name=="TSO linker" else s
                      for s in s2.preamp_template()],name="G&T-inherited full-length cDNA before ISPCR")
FINAL_LIBRARIES=(("RNA library",RNA,RP,"Smart-seq2 cDNA after Nextera XT.","Nextera primers declared for the RNA branch."),("DNA methylome library",DNA,DP,"Bisulfite-converted genomic fragments after scBS-seq.","Illumina primers declared for the DNA branch."))

def oligos():
    return [f"G&T-seq capture oligo-dT30VN (inherited): 5′-biotin-triethyleneglycol-{s2.OLIGO_DT}-3′",
            f"Smart-seq2 TSO (linker inherited, not printed by this paper): "
            f"5′-{rt.SMART_HANDLE}<inf>{G_AND_T_TSO_LINKER}</inf>{rt.TSO_G_TAIL_LNA}-3′",
            f"Smart-seq2 ISPCR primer: 5′-{s2.ISPCR}-3′",
            f"scBS-seq PE1: 5′-{scbs.PE1}-3′",
            f"scM&T oligo1 N6: 5′-{scbs.OLIGO1_HANDLE}{'N'*6}-3′",
            f"scM&T oligo2 N6: 5′-{scbs.OLIGO2_HANDLE}{'N'*6}-3′",
            "scBS-seq indexed iPCRTag: P7–i7–GATCGGAAGAGCACACGTCTGAACTCCAGTCAC"]
def sections():
    return [("Capture mRNA on biotinylated oligo(dT) beads",[Row(chunks=[("mRNA poly(A) || biotin—oligo(dT) ● bead", "r1", False)])],"The bead-bound fraction is physically removed from genomic DNA."),("Separate the analytes",analyte_split("bead-bound RNA → Smart-seq2","supernatant DNA → scBS-seq"),"Wash fractions are retained with the DNA to reduce loss."),("Convert and amplify each branch",[Row(chunks=[("RNA: RT → template switch → cDNA PCR → Nextera", "tso", False)]),Row(chunks=[("DNA: bisulfite → random priming → adapter PCR", "r2", False)])],"The branch chemistries remain distinct after separation."),]


def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    lysate=g.input("lysate","one-cell RLT Plus lysate (1 U/µl SUPERase-In) containing genomic DNA and poly(A) RNA",MaterialKind.BIOLOGICAL,state=sections()[0][1],primary=True)
    beads=g.input("oligodt_beads","biotinylated oligo(dT) magnetic beads",MaterialKind.SELECTION)
    rna=g.product("rna_fraction","bead-bound poly(A) RNA",MaterialKind.MOLECULAR,state=sections()[1][1][:1])
    dna=g.product("dna_fraction","genomic-DNA supernatant plus retained washes",MaterialKind.MOLECULAR,state=sections()[1][1][1:])
    g.react("separate","Capture mRNA and separate it from genomic DNA",StepKind.AFFINITY_CAPTURE,((lysate,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),(rna,dna),changes_molecule=False)

    rt_primer=g.input("rt_primer","G&T-seq 5′-biotin-triethyleneglycol anchored oligo(dT) primer",MaterialKind.OLIGO,
                      state=(seg("SMART handle",rt.SMART_HANDLE,"tso"),seg("oligo-dT linker",s2.DT_LINKER),seg("dT30VN","T"*30+"VN",placeholder=True)),
                      source="Macaulay et al. 2015 G&T-seq Online Methods p. 4, lines 569–576")
    rt_enzyme=g.input("rt","template-switching reverse transcriptase",MaterialKind.CATALYST)
    copied=g.product("copied","first-strand Smart-seq2 cDNA",MaterialKind.MOLECULAR,state=RNA_PREAMP)
    g.react("reverse_transcribe","Reverse-transcribe bead-captured RNA",StepKind.REVERSE_TRANSCRIPTION,((rna,InputRole.TEMPLATE),(rt_primer,InputRole.STRUCTURAL),(rt_enzyme,InputRole.CATALYST)),(copied,))
    tso=g.input("tso","G&T-seq template-switch oligo inherited by scM&T",MaterialKind.OLIGO,state=(seg("SMART handle",rt.SMART_HANDLE,"tso"),seg("TSO linker",G_AND_T_TSO_LINKER,inferred=True),seg("rGrG+G","GGG","tso")))
    switched=g.product("switched","template-switched full-length cDNA",MaterialKind.MOLECULAR,state=RNA_PREAMP)
    g.react("template_switch","INFERRED — Template-switch the nascent cDNA (the AC linker is printed by cited G&T-seq; this paper prints no TSO)",StepKind.TEMPLATE_SWITCHING,((copied,InputRole.TEMPLATE),(tso,InputRole.STRUCTURAL),(rt_enzyme,InputRole.CATALYST)),(switched,))
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

    # Online Methods: DNA eluted into 10 µl H2O after G&T separation; EZ Methylation
    # Direct (Zymo) conversion on purified DNA with AMPure XP beads present;
    # desulphonation/purification on Zymo magnetic beads, eluted into first-strand mix.
    bisulfite=g.input("bisulfite","EZ Methylation Direct bisulfite reagent (Zymo), desulphonated on Zymo magnetic beads",MaterialKind.CONDITION)
    converted=g.product("converted","bisulfite-converted genomic DNA",MaterialKind.MOLECULAR,state=DNA)
    g.react("bisulfite_convert","Bisulfite-convert the bead-purified genomic DNA and desulphonate on magnetic beads",StepKind.CHEMICAL_CONVERSION,((dna,InputRole.TEMPLATE),(bisulfite,InputRole.CONDITION)),(converted,))
    random_primer=g.input("random_primer","scM&T oligo1 with a six-base random 3′ end",MaterialKind.OLIGO,
                          state=(seg("oligo1 inherited handle",scbs.OLIGO1_HANDLE,"r1"),
                                 seg("paper-specified random hexamer","N"*6,placeholder=True)),
                          source="Angermueller et al. Methods: N6 replaces the inherited scBS-seq N9 tail")
    dna_polymerase=g.input("dna_strand_polymerase","scBS-seq strand-synthesis polymerase",MaterialKind.CATALYST)
    primed=g.product("primed_dna","oligo1-primed bisulfite-converted DNA",MaterialKind.MOLECULAR,state=DNA)
    g.react("prime_dna","Random-prime converted genomic DNA with oligo1 (repeated rounds; no first-strand biotin capture)",StepKind.EXTENSION,((converted,InputRole.TEMPLATE),(random_primer,InputRole.STRUCTURAL),(dna_polymerase,InputRole.CATALYST)),(primed,))
    oligo2=g.input("oligo2","scM&T oligo2 with a six-base random 3′ end",MaterialKind.OLIGO,
                   state=(seg("oligo2 inherited handle",scbs.OLIGO2_HANDLE,"r2"),
                          seg("paper-specified random hexamer","N"*6,placeholder=True)),
                   source="Angermueller et al. Methods: N6 replaces the inherited scBS-seq N9 tail")
    dna_adapted=g.product("dna_adapted","oligo1/oligo2-completed scBS-seq template",MaterialKind.MOLECULAR,state=DNA)
    g.react("adapt_dna","Random-prime the complementary strand with scBS-seq oligo2",StepKind.EXTENSION,((primed,InputRole.TEMPLATE),(oligo2,InputRole.STRUCTURAL),(dna_polymerase,InputRole.CATALYST)),(dna_adapted,))

    # Online Methods: an extra 0.8x AMPure XP purification between second-strand
    # synthesis and PCR replaces the omitted streptavidin capture.
    ampure=g.input("dna_ampure","0.8x AMPure XP beads",MaterialKind.SELECTION)
    dna_clean=g.product("dna_clean","AMPure-retained oligo1/oligo2-completed template",MaterialKind.MOLECULAR,state=DNA)
    ampure_waste=g.product("dna_ampure_supernatant","0.8x AMPure XP supernatant (unretained small material, including free primers)",MaterialKind.WASTE)
    g.react("cleanup_dna","Purify second-strand products with 0.8x AMPure XP before PCR",StepKind.CLEANUP,((dna_adapted,InputRole.TEMPLATE),(ampure,InputRole.SELECTION)),(dna_clean,ampure_waste),changes_molecule=False)
    g.discard(ampure_waste,"bead supernatant removed before library PCR")
    reviewed=["lysate","oligodt_beads","separate","rt_primer","rt","reverse_transcribe","tso","template_switch","ispcr_forward","ispcr_reverse","preamp_polymerase","preamp_rna","rna_transposome","rna_tn5","tagment_rna","bisulfite","bisulfite_convert","random_primer","dna_strand_polymerase","prime_dna","oligo2","adapt_dna","dna_ampure","cleanup_dna"]
    finish_pcr_branch(g,"rna",rna_tagged,RNA,RP,reviewed,label="scM&T-seq RNA",read_lengths=READS_125PE,geometry=s2.nextera_index_geometry(RNA))
    ipcrtag=il.P7+"N"*scbs.IPCRTAG_INDEX_NT+scbs.IPCRTAG_3_ARM
    dna_geometry=PCRGeometry.exact(DNA,PCRPrimer("PE1",scbs.PE1),PCRPrimer("indexed iPCRTag",ipcrtag,anneal_nt=len(scbs.IPCRTAG_3_ARM)),statement="Smallwood et al. print PE1 and the indexed iPCRTag family; their annealing regions face inward on the scBS-seq product.")
    dna_forward=g.input("dna_forward","scM&T-seq DNA methylome forward PCR primer",MaterialKind.OLIGO)
    dna_reverse=g.input("dna_reverse","scM&T-seq DNA methylome reverse PCR primer",MaterialKind.OLIGO)
    dna_pcr_polymerase=g.input("dna_polymerase_pcr","scM&T-seq DNA methylome PCR polymerase",MaterialKind.CATALYST)
    dna_library=g.product("dna_library","scM&T-seq DNA methylome sequencing library",MaterialKind.LIBRARY,state=DNA)
    primer_concatemers=g.product("primer_concatemers","first-strand primer concatemers observed among sequenced DNA fragments",MaterialKind.WASTE,
                                 source="Angermueller et al. 2016, Supplementary Fig. 3e caption; exact assembly not established")
    g.react("dna_pcr","Amplify and index the scM&T-seq DNA methylome branch",StepKind.PCR,
            ((dna_clean,InputRole.TEMPLATE),(dna_forward,InputRole.STRUCTURAL),
             (dna_reverse,InputRole.STRUCTURAL),(dna_pcr_polymerase,InputRole.CATALYST)),
            (dna_library,primer_concatemers),geometry=dna_geometry)
    g.sequence(dna_library,SequencingPlan(library=DNA,primers=DP,acquired_roles=tuple(p.role for p in DP),read_lengths=READS_125PE))
    g.discard(primer_concatemers,"observed among sequenced reads and limits genomic alignment; exact formation stage is unresolved")
    reviewed.extend((dna_forward.id,dna_reverse.id,dna_pcr_polymerase.id,"dna_pcr",dna_library.id))
    g.declare_reviewed_scope(reviewed,SOURCE)
    return g.require_publishable()
