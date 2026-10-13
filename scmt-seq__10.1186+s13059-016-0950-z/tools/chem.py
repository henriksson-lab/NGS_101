"""scMT-seq joint RNA and methylome workflow."""
from chemdraw import Row
from adapter_ligation import AdapterLigationGeometry
from multibranch import illumina_branch
from multimodal_spatial import finish_pcr_branch
from pcr import PCRGeometry
from protocol_extensions import analyte_split
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind
TITLE="scMT-seq — transcriptome and DNA methylome"
NOTES="01_scmt-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1186/s13059-016-0950-z">Hu et al., <i>Genome Biology</i> (2016)</a>.'
SUMMARY="Cytoplasmic RNA and nuclear genomic DNA from one cell are physically separated for transcriptome sequencing and single-cell methylome sequencing."
CAVEAT="The paper establishes physical cytoplasm/nucleus separation and the two branch workflows. This compact page has not transcribed the exact RNA amplification or methylome adapter/PCR oligos, so those geometries remain explicitly unavailable."
RNA,RP=illumina_branch("scMT-seq RNA library","single-cell cDNA",architecture="nextera")
DNA,DP=illumina_branch("scMT-seq methylome library","bisulfite-converted genomic DNA")
FINAL_LIBRARIES=(("Transcriptome library",RNA,RP,"Amplified single-cell cDNA is converted to an Illumina library.","Declared primers read the RNA branch."),("DNA methylome library",DNA,DP,"Bisulfite-converted nuclear DNA yields the methylome library.","Declared primers read the DNA branch."))
def sections(): return [
 ("Lyse one cell and retain its nucleus",[Row(chunks=[("cell → cytoplasmic RNA + nuclear genomic DNA","r1",False)])],"The fractions are separated before either amplification."),
 ("Separate the matched analytes",analyte_split("cytoplasmic RNA → transcriptome","nuclear DNA → methylome"),"Cell identity is maintained by plate position."),
 ("Convert the RNA fraction",[Row(chunks=[("reverse transcription → cDNA amplification → Illumina library","tso",False)])],"The cytoplasmic branch reports expression."),
 ("Convert the DNA fraction",[Row(chunks=[("bisulfite conversion → random priming → indexed library","r2",False)])],"The nuclear branch reports methylated cytosines."),]


def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 cell=g.input("cell","one intact cell",MaterialKind.BIOLOGICAL,state=sections()[0][1],primary=True)
 cytoplasm=g.product("cytoplasm","cytoplasmic RNA fraction",MaterialKind.MOLECULAR,state=sections()[1][1][:1])
 nucleus=g.product("nucleus","nuclear genomic-DNA fraction",MaterialKind.MOLECULAR,state=sections()[1][1][1:])
 g.react("separate","Lyse the cell while retaining and separating its nucleus",StepKind.SPLIT,((cell,InputRole.TEMPLATE),),(cytoplasm,nucleus),changes_molecule=False)

 rt_primer=g.input("rt_primer","single-cell RNA reverse-transcription primer",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (exact oligo not transcribed here)")
 rt=g.input("rt","reverse transcriptase",MaterialKind.CATALYST)
 cdna=g.product("cdna","cytoplasmic transcript cDNA",MaterialKind.MOLECULAR,state=RNA,inferred=True,source=SOURCE+" (exact RT intermediate not transcribed here)")
 g.react("reverse_transcribe","Reverse-transcribe the cytoplasmic RNA fraction",StepKind.REVERSE_TRANSCRIPTION,((cytoplasm,InputRole.TEMPLATE),(rt_primer,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(cdna,))
 rna_f=g.input("rna_preamp_forward","RNA-branch preamplification forward primer",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (exact primer not transcribed here)")
 rna_r=g.input("rna_preamp_reverse","RNA-branch preamplification reverse primer",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (exact primer not transcribed here)")
 rna_pol=g.input("rna_preamp_polymerase","RNA-branch preamplification polymerase",MaterialKind.CATALYST)
 preamp=g.product("rna_preamp","amplified single-cell cDNA",MaterialKind.MOLECULAR,state=RNA,inferred=True,source=SOURCE+" (exact PCR geometry not transcribed here)")
 g.react("preamp_rna","Amplify the transcript cDNA",StepKind.PCR,((cdna,InputRole.TEMPLATE),(rna_f,InputRole.STRUCTURAL),(rna_r,InputRole.STRUCTURAL),(rna_pol,InputRole.CATALYST)),(preamp,),geometry=PCRGeometry.unavailable(RNA,reason="The exact single-cell RNA preamplification primer pair has not been transcribed in this page."))
 rna_adapter=g.input("rna_adapter","RNA-branch Illumina adapter",MaterialKind.ADAPTER,inferred=True,source=SOURCE+" (exact adapter not transcribed here)")
 rna_adapted=g.product("rna_adapted","adapter-completed transcript cDNA",MaterialKind.MOLECULAR,state=RNA,inferred=True,source=SOURCE+" (exact adapter geometry not transcribed here)")
 g.react("adapt_rna","Attach the RNA-branch Illumina adapter",StepKind.LIGATION,((preamp,InputRole.TEMPLATE),(rna_adapter,InputRole.STRUCTURAL)),(rna_adapted,),geometry=AdapterLigationGeometry.unavailable("The compact page does not transcribe the RNA-branch adapter-end geometry."))

 bisulfite=g.input("bisulfite","bisulfite-conversion condition",MaterialKind.CONDITION)
 converted=g.product("converted","bisulfite-converted nuclear DNA",MaterialKind.MOLECULAR,state=DNA)
 g.react("bisulfite_convert","Bisulfite-convert the nuclear genomic DNA",StepKind.CHEMICAL_CONVERSION,((nucleus,InputRole.TEMPLATE),(bisulfite,InputRole.CONDITION)),(converted,))
 random_primer=g.input("random_primer","methylome random primer",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (exact oligo not transcribed here)")
 strand_pol=g.input("strand_polymerase","methylome strand-synthesis polymerase",MaterialKind.CATALYST)
 primed=g.product("primed","random-primed converted genomic DNA",MaterialKind.MOLECULAR,state=DNA,inferred=True,source=SOURCE+" (exact random-primer architecture not transcribed here)")
 g.react("prime_dna","Random-prime bisulfite-converted DNA",StepKind.EXTENSION,((converted,InputRole.TEMPLATE),(random_primer,InputRole.STRUCTURAL),(strand_pol,InputRole.CATALYST)),(primed,))
 dna_adapter=g.input("dna_adapter","methylome terminal adapter",MaterialKind.ADAPTER,inferred=True,source=SOURCE+" (exact adapter not transcribed here)")
 dna_adapted=g.product("dna_adapted","adapter-completed methylome template",MaterialKind.MOLECULAR,state=DNA,inferred=True,source=SOURCE+" (exact adapter geometry not transcribed here)")
 g.react("adapt_dna","Attach the methylome terminal adapter",StepKind.LIGATION,((primed,InputRole.TEMPLATE),(dna_adapter,InputRole.STRUCTURAL)),(dna_adapted,),geometry=AdapterLigationGeometry.unavailable("The compact page does not transcribe the methylome adapter-end geometry."))

 reviewed=["cell","separate","rt_primer","rt","reverse_transcribe","rna_preamp_forward","rna_preamp_reverse","rna_preamp_polymerase","preamp_rna","rna_adapter","adapt_rna","bisulfite","bisulfite_convert","random_primer","strand_polymerase","prime_dna","dna_adapter","adapt_dna"]
 finish_pcr_branch(g,"rna",rna_adapted,RNA,RP,reviewed,label="scMT-seq transcriptome",inferred_primers=True,limitation="The exact RNA terminal PCR-primer pair has not been transcribed in this page.")
 finish_pcr_branch(g,"dna",dna_adapted,DNA,DP,reviewed,label="scMT-seq methylome",inferred_primers=True,limitation="The exact methylome terminal PCR-primer pair has not been transcribed in this page.")
 g.declare_reviewed_scope(reviewed,SOURCE)
 return g.require_publishable()
