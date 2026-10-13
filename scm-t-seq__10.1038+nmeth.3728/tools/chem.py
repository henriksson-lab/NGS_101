"""scM&T-seq physical RNA/DNA separation."""
from chemdraw import Row
from adapter_ligation import AdapterLigationGeometry
from multibranch import illumina_branch
from multimodal_spatial import finish_pcr_branch
from pcr import PCRGeometry
from protocol_extensions import analyte_split
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind
TITLE="scM&T-seq — transcriptome and DNA methylome"
NOTES="01_scm-t-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/nmeth.3728">Angermueller et al., <i>Nature Methods</i> (2016)</a>.'
SUMMARY="Biotinylated oligo(dT) beads remove a cell's mRNA for Smart-seq2 while the genomic-DNA supernatant enters scBS-seq."
CAVEAT="The defining method establishes bead-based physical separation and the Smart-seq2/scBS-seq branch chemistries. Exact branch-specific oligos are not transcribed in this compact page, so their adapter and PCR landing geometries remain explicitly unavailable."
RNA,RP=illumina_branch("scM&T-seq RNA","Smart-seq2 cDNA",architecture="nextera")
DNA,DP=illumina_branch("scM&T-seq DNA","bisulfite-converted genomic DNA")
FINAL_LIBRARIES=(("RNA library",RNA,RP,"Smart-seq2 cDNA after Nextera XT.","Nextera primers declared for the RNA branch."),("DNA methylome library",DNA,DP,"Bisulfite-converted genomic fragments after scBS-seq.","Illumina primers declared for the DNA branch."))
def sections():
    return [("Capture mRNA on biotinylated oligo(dT) beads",[Row(chunks=[("mRNA poly(A) || biotin—oligo(dT) ● bead", "r1", False)])],"The bead-bound fraction is physically removed from genomic DNA."),("Separate the analytes",analyte_split("bead-bound RNA → Smart-seq2","supernatant DNA → scBS-seq"),"Wash fractions are retained with the DNA to reduce loss."),("Convert and amplify each branch",[Row(chunks=[("RNA: RT → template switch → cDNA PCR → Nextera", "tso", False)]),Row(chunks=[("DNA: bisulfite → random priming → adapter PCR", "r2", False)])],"The branch chemistries remain distinct after separation."),]


def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    lysate=g.input("lysate","one-cell lysate containing genomic DNA and poly(A) RNA",MaterialKind.BIOLOGICAL,state=sections()[0][1],primary=True)
    beads=g.input("oligodt_beads","biotinylated oligo(dT) magnetic beads",MaterialKind.SELECTION)
    rna=g.product("rna_fraction","bead-bound poly(A) RNA",MaterialKind.MOLECULAR,state=sections()[1][1][:1])
    dna=g.product("dna_fraction","genomic-DNA supernatant plus retained washes",MaterialKind.MOLECULAR,state=sections()[1][1][1:])
    g.react("separate","Capture mRNA and separate it from genomic DNA",StepKind.AFFINITY_CAPTURE,((lysate,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),(rna,dna),changes_molecule=False)

    rt_primer=g.input("rt_primer","Smart-seq2 anchored oligo(dT) primer",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (exact branch oligo not transcribed here)")
    rt=g.input("rt","template-switching reverse transcriptase",MaterialKind.CATALYST)
    copied=g.product("copied","first-strand Smart-seq2 cDNA",MaterialKind.MOLECULAR,state=RNA,inferred=True,source=SOURCE+" (exact RT intermediate not transcribed here)")
    g.react("reverse_transcribe","Reverse-transcribe bead-captured RNA",StepKind.REVERSE_TRANSCRIPTION,((rna,InputRole.TEMPLATE),(rt_primer,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(copied,))
    tso=g.input("tso","Smart-seq2 template-switch oligo",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (exact branch oligo not transcribed here)")
    switched=g.product("switched","template-switched full-length cDNA",MaterialKind.MOLECULAR,state=RNA,inferred=True,source=SOURCE+" (exact TSO intermediate not transcribed here)")
    g.react("template_switch","Template-switch the nascent cDNA",StepKind.TEMPLATE_SWITCHING,((copied,InputRole.TEMPLATE),(tso,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(switched,))
    ispcr_f=g.input("ispcr_forward","Smart-seq2 cDNA PCR forward primer",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (exact primer not transcribed here)")
    ispcr_r=g.input("ispcr_reverse","Smart-seq2 cDNA PCR reverse primer",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (exact primer not transcribed here)")
    prepol=g.input("preamp_polymerase","Smart-seq2 cDNA PCR polymerase",MaterialKind.CATALYST)
    preamp=g.product("preamp","amplified full-length Smart-seq2 cDNA",MaterialKind.MOLECULAR,state=RNA,inferred=True,source=SOURCE+" (exact preamplification geometry not transcribed here)")
    g.react("preamp_rna","Amplify full-length cDNA",StepKind.PCR,((switched,InputRole.TEMPLATE),(ispcr_f,InputRole.STRUCTURAL),(ispcr_r,InputRole.STRUCTURAL),(prepol,InputRole.CATALYST)),(preamp,),geometry=PCRGeometry.unavailable(RNA,reason="The page has not transcribed the exact Smart-seq2 cDNA-amplification primer pair."))
    transposome=g.input("rna_transposome","Nextera XT transposome",MaterialKind.ADAPTER)
    tn5=g.input("rna_tn5","Tn5 transposase",MaterialKind.CATALYST)
    rna_tagged=g.product("rna_tagged","Nextera-tagged Smart-seq2 cDNA",MaterialKind.MOLECULAR,state=RNA)
    twaste=g.product("rna_tn5_waste","non-library cDNA tagmentation products",MaterialKind.WASTE)
    g.react("tagment_rna","Tagment amplified cDNA with Nextera XT",StepKind.TAGMENTATION,((preamp,InputRole.TEMPLATE),(transposome,InputRole.STRUCTURAL),(tn5,InputRole.CATALYST)),(rna_tagged,twaste));g.discard(twaste,"does not retain an amplifiable unlike-ended cDNA fragment")

    bisulfite=g.input("bisulfite","bisulfite-conversion condition",MaterialKind.CONDITION)
    converted=g.product("converted","bisulfite-converted genomic DNA",MaterialKind.MOLECULAR,state=DNA)
    g.react("bisulfite_convert","Bisulfite-convert the genomic-DNA fraction",StepKind.CHEMICAL_CONVERSION,((dna,InputRole.TEMPLATE),(bisulfite,InputRole.CONDITION)),(converted,))
    random_primer=g.input("random_primer","scBS-seq random primer",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (exact branch oligo not transcribed here)")
    dna_polymerase=g.input("dna_strand_polymerase","scBS-seq strand-synthesis polymerase",MaterialKind.CATALYST)
    primed=g.product("primed_dna","random-primed bisulfite-converted DNA",MaterialKind.MOLECULAR,state=DNA,inferred=True,source=SOURCE+" (exact random-primer architecture not transcribed here)")
    g.react("prime_dna","Random-prime converted genomic DNA",StepKind.EXTENSION,((converted,InputRole.TEMPLATE),(random_primer,InputRole.STRUCTURAL),(dna_polymerase,InputRole.CATALYST)),(primed,))
    dna_adapter=g.input("dna_adapter","scBS-seq terminal adapter",MaterialKind.ADAPTER,inferred=True,source=SOURCE+" (exact adapter not transcribed here)")
    dna_adapted=g.product("dna_adapted","adapter-completed scBS-seq template",MaterialKind.MOLECULAR,state=DNA,inferred=True,source=SOURCE+" (exact adapter geometry not transcribed here)")
    g.react("adapt_dna","Attach the scBS-seq terminal adapter",StepKind.LIGATION,((primed,InputRole.TEMPLATE),(dna_adapter,InputRole.STRUCTURAL)),(dna_adapted,),geometry=AdapterLigationGeometry.unavailable("The compact page does not transcribe the scBS-seq terminal adapter-end geometry."))

    reviewed=["lysate","oligodt_beads","separate","rt_primer","rt","reverse_transcribe","tso","template_switch","ispcr_forward","ispcr_reverse","preamp_polymerase","preamp_rna","rna_transposome","rna_tn5","tagment_rna","bisulfite","bisulfite_convert","random_primer","dna_strand_polymerase","prime_dna","dna_adapter","adapt_dna"]
    finish_pcr_branch(g,"rna",rna_tagged,RNA,RP,reviewed,label="scM&T-seq RNA",inferred_primers=True,limitation="The exact Nextera XT terminal PCR-primer pair has not been transcribed in this page.")
    finish_pcr_branch(g,"dna",dna_adapted,DNA,DP,reviewed,label="scM&T-seq DNA methylome",inferred_primers=True,limitation="The exact scBS-seq terminal PCR-primer pair has not been transcribed in this page.")
    g.declare_reviewed_scope(reviewed,SOURCE)
    return g.require_publishable()
