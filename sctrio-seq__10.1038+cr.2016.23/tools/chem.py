"""scTrio-seq cytoplasm/nucleus split and its RNA and DNA branches."""
from chemdraw import Construct, MolecularState, Row, Workflow
from batch_ngs import seg
from adapter_ligation import AdapterLigationGeometry
from endprep import (ThreePrimeOverhangAdapter, dA_tailed_scene,
                     ligation_product_scene, repair_and_dA_tail,
                     cohesive_ligation_geometry)
from multibranch import illumina_branch
from multimodal_spatial import finish_pcr_branch
from protocol_extensions import analyte_split
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind
TITLE="scTrio-seq — genome copy number, methylome and transcriptome"
NOTES="01_sctrio-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/cr.2016.23">Hou et al., <i>Cell Research</i> (2016)</a>.'
SUMMARY="A single cell is mechanically separated into cytoplasm for RNA-seq and nucleus for scRRBS; the DNA library jointly reports copy number and CpG methylation."
CAVEAT="The scRRBS branch has source-defined MspI processing and premethylated dA/T adapter geometry. This page has not yet transcribed the Tang-lineage RNA adapter and terminal PCR oligos, so those reaction geometries remain explicitly unavailable."
RNA,RP=illumina_branch("scTrio-seq RNA","cytoplasmic cDNA")
DNA,DP=illumina_branch("scTrio-seq DNA","MspI reduced-representation insert")
FINAL_LIBRARIES=(("Transcriptome library",RNA,RP,"Cytoplasmic RNA-derived library.","Illumina primers on the RNA branch."),("Genome and methylome library",DNA,DP,"One scRRBS library supports both copy-number and methylation inference.","Illumina primers on the DNA branch."))


def _tailed_fragment():
    return repair_and_dA_tail(Construct([
        seg("MspI reduced-representation insert", "X" * 36, placeholder=True)]))


def _workflow():
    w=Workflow(MolecularState("one cell",(Row(chunks=[("one intact cell",None,False)]),)))
    w.react("Separate cytoplasm from nucleus",analyte_split("cytoplasm: poly(A) RNA","nucleus: genomic DNA"),name="separated analytes",note="The three measurements arise from two physical branches.")
    w.react("Prepare the cytoplasmic transcriptome",[Row(chunks=[("poly(A) RNA → RT → cDNA amplification → Illumina library","r1",False)])],name="RNA library plus nuclear DNA",note="The RNA workflow follows the single-cell Tang lineage.")
    w.react("Digest nuclear DNA with MspI",[Row(chunks=[("C^CGG → fragments with two-base 5′-CG restriction overhangs","me",False)])],name="MspI fragments",note="MspI enriches CpG-containing fragments.")
    tailed=_tailed_fragment()
    w.react("Fill recessed ends and add 3′ dA",dA_tailed_scene(tailed).rows(),name="dA-tailed RRBS fragments",note="End fill removes the original CG overhang before dA-tailing.")
    adapter=ThreePrimeOverhangAdapter("premethylated TruSeq ligating stem","GCTCTTCCGATCT","GATCGGAAGAGC")
    geometry=cohesive_ligation_geometry(tailed,adapter,statement="scTrio-seq follows scRRBS: the premethylated TruSeq fork carries a 3′-T complementary to the filled fragment's 3′-dA.")
    product=ligation_product_scene(DNA,"MspI reduced-representation insert",label="scTrio-seq DNA library")
    w.ligate_adapter("Ligate premethylated adapters",geometry,product.rows(),name="adapter-ligated RRBS library",note="Bisulfite conversion and PCR follow ligation.")
    w.react("Bisulfite convert and amplify",[Row(chunks=[("scRRBS library → methylation calls + copy-number bins","r2",False)])],name="final DNA library",note="Both DNA measurements derive from this one library.")
    return w

def sections():
    return [(r.action,list(r.after.rows),r.note) for r in _workflow().reactions]


def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    cell=g.input("cell","one intact cell",MaterialKind.BIOLOGICAL,state=(Row(chunks=[("one intact cell",None,False)]),),primary=True)
    cytoplasm=g.product("cytoplasm","cytoplasmic RNA fraction",MaterialKind.MOLECULAR,state=analyte_split("cytoplasm: poly(A) RNA","nucleus: genomic DNA")[:1])
    nucleus=g.product("nucleus","nuclear genomic-DNA fraction",MaterialKind.MOLECULAR,state=analyte_split("cytoplasm: poly(A) RNA","nucleus: genomic DNA")[1:])
    g.react("separate","Mechanically separate cytoplasm from nucleus",StepKind.SPLIT,((cell,InputRole.TEMPLATE),),(cytoplasm,nucleus),changes_molecule=False)

    rt_primer=g.input("rt_primer","Tang-lineage poly(A) reverse-transcription primer",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (exact RNA-branch oligo not yet transcribed)")
    rt=g.input("rt","reverse transcriptase",MaterialKind.CATALYST)
    cdna=g.product("cdna","cytoplasmic transcript cDNA",MaterialKind.MOLECULAR,state=RNA,inferred=True,source=SOURCE+" (RNA-branch intermediate oligos not yet transcribed)")
    g.react("reverse_transcribe","Reverse-transcribe the cytoplasmic RNA fraction",StepKind.REVERSE_TRANSCRIPTION,((cytoplasm,InputRole.TEMPLATE),(rt_primer,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(cdna,))
    rna_adapter=g.input("rna_adapter","RNA-branch Illumina adapter",MaterialKind.ADAPTER,inferred=True,source=SOURCE+" (exact adapter not yet transcribed)")
    rna_adapted=g.product("rna_adapted","adapter-completed cytoplasmic cDNA",MaterialKind.MOLECULAR,state=RNA,inferred=True,source=SOURCE+" (exact adapter geometry not yet transcribed)")
    g.react("adapt_rna","Attach the RNA-branch Illumina adapter",StepKind.LIGATION,((cdna,InputRole.TEMPLATE),(rna_adapter,InputRole.STRUCTURAL)),(rna_adapted,),geometry=AdapterLigationGeometry.unavailable("The defining paper follows the Tang single-cell RNA lineage without printing the terminal adapter-end geometry on this page."))

    mspi=g.input("mspi","MspI restriction endonuclease",MaterialKind.CATALYST)
    fragments=g.product("mspi_fragments","MspI reduced-representation genomic fragments",MaterialKind.MOLECULAR,state=Construct([seg("MspI reduced-representation insert","X"*36,placeholder=True)],name="MspI fragments"))
    g.react("digest","Digest nuclear DNA with MspI",StepKind.DIGESTION,((nucleus,InputRole.TEMPLATE),(mspi,InputRole.CATALYST)),(fragments,))
    repair=g.input("repair_enzymes","recessed-end repair enzymes",MaterialKind.CATALYST)
    repaired=g.product("repaired","end-repaired MspI fragments",MaterialKind.MOLECULAR,state=fragments.state)
    g.react("end_repair","Fill recessed MspI fragment ends",StepKind.END_REPAIR,((fragments,InputRole.TEMPLATE),(repair,InputRole.CATALYST)),(repaired,))
    tailing=g.input("tailing_polymerase","dA-tailing polymerase",MaterialKind.CATALYST)
    tailed=_tailed_fragment()
    da=g.product("da_tailed","3-prime-dA-tailed reduced-representation fragments",MaterialKind.MOLECULAR,state=tailed)
    g.react("da_tail","Add one 3-prime dA to repaired fragments",StepKind.TAILING,((repaired,InputRole.TEMPLATE),(tailing,InputRole.CATALYST)),(da,))
    adapter_model=ThreePrimeOverhangAdapter("premethylated TruSeq ligating stem","GCTCTTCCGATCT","GATCGGAAGAGC")
    geometry=cohesive_ligation_geometry(tailed,adapter_model,statement="scTrio-seq follows scRRBS: the premethylated TruSeq fork carries a 3-prime dT complementary to the filled fragment's 3-prime dA.")
    dna_adapter=g.input("dna_adapter","premethylated TruSeq forked adapter",MaterialKind.ADAPTER,state=geometry.scenes[0])
    ligase=g.input("dna_ligase","adapter DNA ligase",MaterialKind.CATALYST)
    dna_adapted=g.product("dna_adapted","premethylated-adapter-ligated RRBS fragments",MaterialKind.MOLECULAR,state=DNA)
    g.react("ligate_dna","Ligate premethylated dT-overhang adapters",StepKind.LIGATION,((da,InputRole.TEMPLATE),(dna_adapter,InputRole.STRUCTURAL),(ligase,InputRole.CATALYST)),(dna_adapted,),geometry=geometry)
    bisulfite=g.input("bisulfite","bisulfite-conversion condition",MaterialKind.CONDITION)
    converted=g.product("converted","bisulfite-converted RRBS library templates",MaterialKind.MOLECULAR,state=DNA)
    g.react("bisulfite_convert","Bisulfite-convert the adapter-ligated DNA",StepKind.CHEMICAL_CONVERSION,((dna_adapted,InputRole.TEMPLATE),(bisulfite,InputRole.CONDITION)),(converted,))

    reviewed=["cell","separate","rt_primer","rt","reverse_transcribe","rna_adapter","adapt_rna","mspi","digest","repair_enzymes","end_repair","tailing_polymerase","da_tail","dna_adapter","dna_ligase","ligate_dna","bisulfite","bisulfite_convert"]
    finish_pcr_branch(g,"rna",rna_adapted,RNA,RP,reviewed,label="scTrio-seq transcriptome",inferred_primers=True,limitation="The exact Tang-lineage terminal PCR-primer pair has not yet been transcribed from the defining source.")
    finish_pcr_branch(g,"dna",converted,DNA,DP,reviewed,label="scTrio-seq genome and methylome",inferred_primers=True,limitation="The exact scRRBS terminal PCR-primer pair has not yet been transcribed from the defining source.")
    g.declare_reviewed_scope(reviewed,SOURCE)
    return g.require_publishable()
