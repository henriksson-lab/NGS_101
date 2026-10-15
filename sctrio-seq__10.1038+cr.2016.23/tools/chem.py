"""scTrio-seq cytoplasm/nucleus split and its RNA and DNA branches."""
from chemdraw import Construct, MolecularState, Row, Segment, Workflow
from batch_ngs import seg
from adapter_ligation import AdapterLigationGeometry
from endprep import (ThreePrimeOverhangAdapter, dA_tailed_scene,
                     ligation_product_scene, repair_and_dA_tail,
                     cohesive_ligation_geometry)
from multibranch import illumina_branch
from multimodal_spatial import finish_pcr_branch
from pcr import PCRGeometry, PCRPrimer
from protocol_extensions import analyte_split
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind
import rrbs
import tang_cdna as tang
TITLE="scTrio-seq — genome copy number, methylome and transcriptome"
NOTES="01_sctrio-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/cr.2016.23">Hou et al., <i>Cell Research</i> (2016)</a>.'
# cr.2016.23 cites Tang et al. for the cDNA amplification but never reprints its oligos:
# UP1/UP2 and the 5'-amine NH2-UP1/NH2-UP2 are published only in the 2009 supplementary table.
TANG_OLIGO_SOURCE='Inherited oligo source: <a href="https://doi.org/10.1038/nmeth.1315">Tang et al., <i>Nature Methods</i> (2009)</a>, supplementary oligo table.'
SUMMARY="A single cell is mechanically separated into cytoplasm for RNA-seq and nucleus for scRRBS; the DNA library jointly reports copy number and CpG methylation."
CAVEAT="The paper exposes the inherited Tang UP1/UP2 cDNA-amplification oligos and the inherited scRRBS chemistry. The later NEBNext Ultra RNA-library kit is named, but its adapter and indexing-primer bases are not printed in the paper; those terminal arms are therefore marked inferred."
RNA,RP=illumina_branch("scTrio-seq RNA","cytoplasmic cDNA",inferred_adapters=True)
DNA,DP=rrbs.final_library(),rrbs.SEQ_PRIMERS
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
    product=ligation_product_scene(DNA,"bisulfite-converted insert",label="scTrio-seq DNA library")
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

    rt_primer=g.input("rt_primer","UP1 oligo-dT reverse-transcription primer",MaterialKind.OLIGO,state=(Segment("UP1",tang.UP1),),source=TANG_OLIGO_SOURCE)
    rt=g.input("rt","reverse transcriptase",MaterialKind.CATALYST)
    first=g.product("first_strand","UP1-primed cytoplasmic cDNA",MaterialKind.MOLECULAR,state=tang.first_strand())
    g.react("reverse_transcribe","Reverse-transcribe poly(A) RNA with UP1",StepKind.REVERSE_TRANSCRIPTION,((cytoplasm,InputRole.TEMPLATE),(rt_primer,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(first,))
    exo=g.input("exonuclease","Exonuclease I",MaterialKind.CATALYST)
    cleaned=g.product("cleaned_first_strand","first-strand cDNA after free-UP1 removal",MaterialKind.MOLECULAR,state=tang.first_strand())
    primer_waste=g.product("primer_waste","digested unused UP1",MaterialKind.WASTE)
    g.react("remove_up1","Remove unused UP1",StepKind.CLEANUP,((first,InputRole.TEMPLATE),(exo,InputRole.CATALYST)),(cleaned,primer_waste))
    g.discard(primer_waste,"unused reverse-transcription primer is digested")
    tdt=g.input("tdt","terminal deoxynucleotidyl transferase",MaterialKind.CATALYST)
    datp=g.input("datp","dATP tailing substrate",MaterialKind.OLIGO)
    tailed_cdna=g.product("tailed_cdna","TdT poly(A)-tailed first-strand cDNA",MaterialKind.MOLECULAR)
    g.react("tail_cdna","Add a second poly(A) tail with TdT",StepKind.TAILING,((cleaned,InputRole.TEMPLATE),(datp,InputRole.STRUCTURAL),(tdt,InputRole.CATALYST)),(tailed_cdna,))
    up2=g.input("up2","UP2 oligo-dT second-strand primer",MaterialKind.OLIGO,state=(Segment("UP2",tang.UP2),),source=TANG_OLIGO_SOURCE)
    second_pol=g.input("second_strand_polymerase","second-strand polymerase",MaterialKind.CATALYST)
    duplex=g.product("double_stranded_cdna",tang.amplified_cdna().name,MaterialKind.MOLECULAR,state=tang.amplified_cdna())
    g.react("second_strand","Synthesize the second strand from UP2",StepKind.EXTENSION,((tailed_cdna,InputRole.TEMPLATE),(up2,InputRole.STRUCTURAL),(second_pol,InputRole.CATALYST)),(duplex,))
    pcr_pol=g.input("cdna_pcr_polymerase","cDNA PCR polymerase",MaterialKind.CATALYST)
    amplified=g.product("amplified_cdna","UP1/UP2-amplified cDNA",MaterialKind.MOLECULAR,state=tang.amplified_cdna())
    g.react("cdna_pcr","Amplify cDNA with UP1 and UP2",StepKind.PCR,((duplex,InputRole.TEMPLATE),(rt_primer,InputRole.STRUCTURAL),(up2,InputRole.STRUCTURAL),(pcr_pol,InputRole.CATALYST)),(amplified,),geometry=PCRGeometry.exact(tang.amplified_cdna(),PCRPrimer("UP2",tang.UP2),PCRPrimer("UP1",tang.UP1),statement="The inherited Tang supplementary oligo table prints UP1 and UP2, and both sites are located on the amplified cDNA model."))
    aup1=g.input("aup1","5-prime-amine-blocked AUP1 (NH2-UP1)",MaterialKind.OLIGO,state=(Segment("AUP1",tang.AUP1),),source=TANG_OLIGO_SOURCE)
    aup2=g.input("aup2","5-prime-amine-blocked AUP2 (NH2-UP2)",MaterialKind.OLIGO,state=(Segment("AUP2",tang.AUP2),),source=TANG_OLIGO_SOURCE)
    blocked=g.product("amine_blocked_cdna","5-prime-amine-blocked amplified cDNA",MaterialKind.MOLECULAR,state=tang.amplified_cdna())
    g.react("blocked_pcr","Re-amplify with 5-prime-amine-blocked AUP1/AUP2",StepKind.PCR,((amplified,InputRole.TEMPLATE),(aup1,InputRole.STRUCTURAL),(aup2,InputRole.STRUCTURAL),(pcr_pol,InputRole.CATALYST)),(blocked,),geometry=PCRGeometry.exact(tang.amplified_cdna(),PCRPrimer("AUP2",tang.AUP2),PCRPrimer("AUP1",tang.AUP1),statement="AUP1/AUP2 have the published UP1/UP2 bases and a source-reported 5-prime amine block."))
    shear=g.input("rna_shearing","Covaris shearing condition",MaterialKind.CONDITION)
    fragments=g.product("rna_fragments","150-250-bp cDNA fragments",MaterialKind.MOLECULAR)
    g.react("fragment_rna_branch","Shear amplified cDNA",StepKind.FRAGMENTATION,((blocked,InputRole.TEMPLATE),(shear,InputRole.CONDITION)),(fragments,))
    rna_adapter=g.input("rna_adapter","NEBNext Ultra DNA Prep adapter/indexing reagents",MaterialKind.ADAPTER,inferred=True,source=SOURCE+" (kit named, bases not disclosed in the paper)")
    rna_adapted=g.product("rna_adapted","adapter-completed cytoplasmic cDNA",MaterialKind.MOLECULAR,state=RNA,inferred=True,source=SOURCE+" (NEBNext Ultra terminal oligo bases not disclosed)")
    g.react("adapt_rna","Build the barcoded Illumina library with NEBNext Ultra",StepKind.LIGATION,((fragments,InputRole.TEMPLATE),(rna_adapter,InputRole.STRUCTURAL)),(rna_adapted,),geometry=AdapterLigationGeometry.unavailable("The paper names NEBNext Ultra DNA Prep E7370 but does not print its adapter-end bases."))

    mspi=g.input("mspi","MspI restriction endonuclease",MaterialKind.CATALYST)
    dna_fragments=g.product("mspi_fragments","MspI reduced-representation genomic fragments",MaterialKind.MOLECULAR,state=rrbs.digested_fragment_scene())
    g.react("digest","Digest nuclear DNA with MspI",StepKind.DIGESTION,((nucleus,InputRole.TEMPLATE),(mspi,InputRole.CATALYST)),(dna_fragments,))
    repair=g.input("repair_enzymes","recessed-end repair enzymes",MaterialKind.CATALYST)
    repaired=g.product("repaired","end-repaired MspI fragments",MaterialKind.MOLECULAR,state=fragments.state)
    g.react("end_repair","Fill recessed MspI fragment ends",StepKind.END_REPAIR,((dna_fragments,InputRole.TEMPLATE),(repair,InputRole.CATALYST)),(repaired,))
    tailing=g.input("tailing_polymerase","dA-tailing polymerase",MaterialKind.CATALYST)
    tailed=_tailed_fragment()
    da=g.product("da_tailed","3-prime-dA-tailed reduced-representation fragments",MaterialKind.MOLECULAR,state=tailed)
    g.react("da_tail","Add one 3-prime dA to repaired fragments",StepKind.TAILING,((repaired,InputRole.TEMPLATE),(tailing,InputRole.CATALYST)),(da,))
    adapter_model=ThreePrimeOverhangAdapter("premethylated TruSeq ligating stem","GCTCTTCCGATCT","GATCGGAAGAGC")
    geometry=cohesive_ligation_geometry(tailed,adapter_model,statement="scTrio-seq follows scRRBS: the premethylated TruSeq fork carries a 3-prime dT complementary to the filled fragment's 3-prime dA.")
    dna_adapter=g.input("dna_adapter","premethylated indexed TruSeq forked adapter",MaterialKind.ADAPTER,state=rrbs.adapter_scene())
    ligase=g.input("dna_ligase","adapter DNA ligase",MaterialKind.CATALYST)
    dna_adapted=g.product("dna_adapted","premethylated-adapter-ligated RRBS fragments",MaterialKind.MOLECULAR,state=DNA)
    g.react("ligate_dna","Ligate premethylated dT-overhang adapters",StepKind.LIGATION,((da,InputRole.TEMPLATE),(dna_adapter,InputRole.STRUCTURAL),(ligase,InputRole.CATALYST)),(dna_adapted,),geometry=geometry)
    bisulfite=g.input("bisulfite","bisulfite-conversion condition",MaterialKind.CONDITION)
    converted=g.product("converted","bisulfite-converted RRBS library templates",MaterialKind.MOLECULAR,state=DNA)
    g.react("bisulfite_convert","Bisulfite-convert the adapter-ligated DNA",StepKind.CHEMICAL_CONVERSION,((dna_adapted,InputRole.TEMPLATE),(bisulfite,InputRole.CONDITION)),(converted,))

    reviewed=["cell","separate","rt_primer","rt","reverse_transcribe","exonuclease","remove_up1","tdt","datp","tail_cdna","up2","second_strand_polymerase","second_strand","cdna_pcr_polymerase","cdna_pcr","aup1","aup2","blocked_pcr","rna_shearing","fragment_rna_branch","rna_adapter","adapt_rna","mspi","digest","repair_enzymes","end_repair","tailing_polymerase","da_tail","dna_adapter","dna_ligase","ligate_dna","bisulfite","bisulfite_convert"]
    finish_pcr_branch(g,"rna",rna_adapted,RNA,RP,reviewed,label="scTrio-seq transcriptome",inferred_primers=True,limitation="The exact Tang-lineage terminal PCR-primer pair has not yet been transcribed from the defining source.")
    finish_pcr_branch(g,"dna",converted,DNA,DP,reviewed,label="scTrio-seq genome and methylome",limitation="The inherited scRRBS paper identifies QP1/QP2 and fixes the final canonical arms, but the converted pre-PCR intermediate is not represented base-for-base.")
    g.declare_reviewed_scope(reviewed,SOURCE)
    return g.require_publishable()
