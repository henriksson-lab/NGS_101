"""scTrio-seq cytoplasm/nucleus split and its RNA and DNA branches."""
from chemdraw import Construct, MolecularState, Row, Workflow
from batch_ngs import seg
from endprep import (ThreePrimeOverhangAdapter, dA_tailed_scene,
                     ligation_product_scene, repair_and_dA_tail,
                     cohesive_ligation_geometry)
from multibranch import illumina_branch
from protocol_extensions import analyte_split
TITLE="scTrio-seq — genome copy number, methylome and transcriptome"
NOTES="01_sctrio-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/cr.2016.23">Hou et al., <i>Cell Research</i> (2016)</a>.'
SUMMARY="A single cell is mechanically separated into cytoplasm for RNA-seq and nucleus for scRRBS; the DNA library jointly reports copy number and CpG methylation."
RNA,RP=illumina_branch("scTrio-seq RNA","cytoplasmic cDNA")
DNA,DP=illumina_branch("scTrio-seq DNA","MspI reduced-representation insert")
FINAL_LIBRARIES=(("Transcriptome library",RNA,RP,"Cytoplasmic RNA-derived library.","Illumina primers on the RNA branch."),("Genome and methylome library",DNA,DP,"One scRRBS library supports both copy-number and methylation inference.","Illumina primers on the DNA branch."))
def _workflow():
    w=Workflow(MolecularState("one cell",(Row(chunks=[("one intact cell",None,False)]),)))
    w.react("Separate cytoplasm from nucleus",analyte_split("cytoplasm: poly(A) RNA","nucleus: genomic DNA"),name="separated analytes",note="The three measurements arise from two physical branches.")
    w.react("Prepare the cytoplasmic transcriptome",[Row(chunks=[("poly(A) RNA → RT → cDNA amplification → Illumina library","r1",False)])],name="RNA library plus nuclear DNA",note="The RNA workflow follows the single-cell Tang lineage.")
    w.react("Digest nuclear DNA with MspI",[Row(chunks=[("C^CGG → fragments with two-base 5′-CG restriction overhangs","me",False)])],name="MspI fragments",note="MspI enriches CpG-containing fragments.")
    tailed=repair_and_dA_tail(Construct([seg("MspI reduced-representation insert","X"*36,placeholder=True)]))
    w.react("Fill recessed ends and add 3′ dA",dA_tailed_scene(tailed).rows(),name="dA-tailed RRBS fragments",note="End fill removes the original CG overhang before dA-tailing.")
    adapter=ThreePrimeOverhangAdapter("premethylated TruSeq ligating stem","GCTCTTCCGATCT","GATCGGAAGAGC")
    geometry=cohesive_ligation_geometry(tailed,adapter,statement="scTrio-seq follows scRRBS: the premethylated TruSeq fork carries a 3′-T complementary to the filled fragment's 3′-dA.")
    product=ligation_product_scene(DNA,"MspI reduced-representation insert",label="scTrio-seq DNA library")
    w.ligate_adapter("Ligate premethylated adapters",geometry,product.rows(),name="adapter-ligated RRBS library",note="Bisulfite conversion and PCR follow ligation.")
    w.react("Bisulfite convert and amplify",[Row(chunks=[("scRRBS library → methylation calls + copy-number bins","r2",False)])],name="final DNA library",note="Both DNA measurements derive from this one library.")
    return w

def sections():
    return [(r.action,list(r.after.rows),r.note) for r in _workflow().reactions]
