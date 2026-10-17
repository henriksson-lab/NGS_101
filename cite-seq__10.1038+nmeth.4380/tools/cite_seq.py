"""Original CITE-seq with 10x 3-prime transcriptome and ADT readout."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from batch_ngs import seg, truseq_library
import illumina as il
import seqprimers as sp
from chemdraw import Construct, Segment, feature, revcomp
from multimodal_spatial import finish_pcr_branch
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind
from protocol_model import molecule, role

TITLE="CITE-seq"
NOTES="01_cite-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/nmeth.4380">Stoeckius et al., <i>Nature Methods</i> (2017)</a>.'
SUMMARY="Antibodies carry poly(A)-tailed DNA tags containing an antibody identity and an amplification handle. The same oligo-dT bead primer that captures mRNA also copies these antibody-derived tags, transferring the droplet cell barcode and UMI before ADT and transcript libraries are separated by size."
CAVEAT=(
    "The original paper demonstrates both Drop-seq and 10x implementations. This page follows its 10x 3′ v2 CBMC implementation. "
    "Its Online Methods, ‘Antibody-oligo sequences’ and ‘10x – CITE-seq’, identify the v2 CBMC branch; "
    "<a href='https://pmc.ncbi.nlm.nih.gov/articles/PMC5669064/bin/NIHMS891423-supplement-1.pdf'>Supplementary Table 2</a> maps its 13 antibodies to six-base identifiers. "
    "Those identifiers crosswalk to the author <a href='https://static-content.springer.com/esm/art%3A10.1038%2Fnmeth.4380/MediaObjects/41592_2017_BFnmeth4380_MOESM3_ESM.pdf'>protocol's v2_BC1–6 and v2_BC8–14</a>; "
    "v2_BC7 appears only in the protocol list and is not assigned to the 13-antibody panel. "
    "The construct below uses a role placeholder for whichever panel antibody is selected. "
    "The complete handle-bearing oligo for each CD8a FACS comparison aliquot is not identified separately in the Methods. "
    "The ADT library is assembled from the printed v2_BC oligo, PE 1.0 P5 and RPI1 primers; the bead poly(dT) length comes from the vendor bead, not this paper, and is drawn as an inferred 30-base placeholder. "
    "Neither the paper nor the author protocol states the acquired read cycles (HiSeq 2500 Rapid Run v2; ADT 10% and cDNA 90% of a lane) or the sequencing primers for the small-RNA arm, so no cycle counts are assigned and only TruSeq Read 1 is declared for the ADT library."
)
CELL=feature("cell_barcode","cell_barcode","whitelist",whitelist="10x-chromium-3prime")
UMI=feature("umi","umi","random")
ADT_FEATURE=feature("antibody_barcode","feature_barcode","unknown")
ADT_I7=feature("adt_sample_index","sample_index","unknown")

# Author protocol (MOESM3, version 2017-07-11) p.6-7, text lines 403-442: printed oligos.
ADT_HANDLE="CCTTGGCACCCGAGAATTCC"   # TruSeq small-RNA Read-2 handle in every v2_BC oligo
V2_BC1="/5AmMC12/CCTTGGCACCCGAGAATTCCAATCACGBAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
PE1_P5="AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGAT*C*T"
RPI1="CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
_RPI_P7, _RPI_INDEX, RPI_ARM = RPI1[:24], RPI1[24:30], RPI1[30:]
if _RPI_P7 != il.P7 or not RPI_ARM.endswith(ADT_HANDLE+"A"):
    raise ValueError("RPI1 transcription no longer splits into P7, six-base index and handle arm")
if PE1_P5.replace("*","") != il.P5[:-4]+il.TRUSEQ_READ1:
    raise ValueError("PE 1.0 P5 transcription differs from the P5 + TruSeq Read 1 arm")

GEX,GEX_P=truseq_library([seg("cell barcode","B"*16,"cbc",placeholder=True,feature=CELL),seg("UMI","U"*10,"umi",placeholder=True,feature=UMI),seg("poly(dT) junction","T"*10,placeholder=True),seg("mRNA cDNA","X"*34,placeholder=True)],"CITE-seq transcript library",dual_index=False)
# ADT library top strand: PE 1.0 P5 end, bead cell barcode/UMI/poly(dT), then the copy of the
# antibody oligo (complement of B, antibody identifier, fixed A, handle) and the RPI1 tail.
ADT=Construct([
    Segment("P5 flank",il.P5[:-4],"p5"),
    Segment("TruSeq Read 1",il.TRUSEQ_READ1,"r1"),
    Segment("cell barcode","B"*16,"cbc",placeholder=True,feature=CELL),
    Segment("UMI","U"*10,"umi",placeholder=True,feature=UMI),
    Segment("bead poly(dT) / ADT poly(A) copy","T"*30,placeholder=True,inferred=True),
    Segment("complement of ordered B","V",placeholder=True),
    Segment("antibody barcode reverse complement","N"*6,"cbc",placeholder=True,feature=ADT_FEATURE),
    Segment("complement of shared A","T"),
    Segment("ADT handle reverse complement",revcomp(ADT_HANDLE),"r2"),
    Segment("RPI arm reverse complement",revcomp(RPI_ARM[:-len(ADT_HANDLE)-1]),"r2"),
    Segment("i7 reverse complement","I"*len(_RPI_INDEX),"cbc",placeholder=True,feature=ADT_I7),
    Segment("P7 reverse complement",il.P7_RC,"p7")],name="CITE-seq ADT library")
ADT_P=(sp.TRUSEQ["R1"],)
_problems=sp.verify(ADT,ADT_P,required_roles=("Read 1",))
if _problems:
    raise ValueError("; ".join(_problems))
FINAL_LIBRARIES=(
    ("Final transcript library",GEX,GEX_P,"The 10x cell barcode and UMI identify each transcript molecule.","The transcript library follows the 10x v2 manufacturer workflow; this page draws its TruSeq shell without asserting acquired cycles."),
    ("Final antibody-derived-tag library",ADT,ADT_P,"The same cell barcode and UMI are followed by the copied antibody tag, its small-RNA handle and the RPI index arm.","Only TruSeq Read 1 is declared: the sequencing primers for the small-RNA arm and the acquired read cycles are not printed in the inspected sources."),)
READ_LENGTHS={
    "Final transcript library":{},
    "Final antibody-derived-tag library":{},
}
READ_USES={
    "Final antibody-derived-tag library":{},
}
SHOW_READ_CONTENT=True

def oligos():
    return [
        "<b>v2_BC1 antibody oligo</b> — 5′-"+V2_BC1+"-3′",
        "<b>10x v2 ADT P5 PCR primer (Illumina PE 1.0)</b> — 5′-"+PE1_P5+"-3′",
        "<b>ADT RPI1 index primer (TruSeq Small RNA)</b> — 5′-"+RPI1+"-3′",
    ]

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    cells=g.input("cells","cells containing poly(A) RNA",MaterialKind.BIOLOGICAL,state=molecule("cells",(role("surface proteins"),role("poly(A) RNA"))),primary=True)
    fc_block=g.input("fc_block","TruStain FcX Fc-receptor block (10 min, 4 °C)",MaterialKind.CONDITION)
    antibodies=g.input("antibodies","streptavidin-antibody conjugates carrying S-S-biotinylated poly(A) DNA tags",MaterialKind.SELECTION)
    stained=g.product("stained","cells bearing antibody-derived DNA tags",MaterialKind.MOLECULAR,state=molecule("CITE-seq stained cells",(role("poly(A) RNA"),role("poly(A)-tailed antibody tag"))))
    unbound=g.product("unbound","unbound antibody tags",MaterialKind.WASTE)
    g.react("stain","Label surface proteins (30 min, 4 °C) and wash 3×",StepKind.BINDING,((cells,InputRole.TEMPLATE),(fc_block,InputRole.CONDITION),(antibodies,InputRole.SELECTION)),(stained,unbound))
    g.discard(unbound,"removed in three staining-buffer washes")
    droplets=g.input("droplets","10x Single Cell 3′ v2 droplet partitioning",MaterialKind.CONDITION)
    occupied=g.product("occupied","occupied GEM products",MaterialKind.MOLECULAR,state=stained.state)
    failed=g.product("failed","empty or unresolved droplet products",MaterialKind.WASTE)
    g.react("partition","Partition stained cells with barcoded beads",StepKind.PARTITION,((stained,InputRole.TEMPLATE),(droplets,InputRole.CONDITION)),(occupied,failed),changes_molecule=False)
    g.discard(failed,"does not retain a resolvable cell barcode")
    bead=g.input("bead_primer","10x v2 cell-barcode/UMI oligo(dT) bead primer",MaterialKind.OLIGO)
    rt=g.input("rt","10x reverse transcription (MMLV-type DNA-dependent polymerase activity)",MaterialKind.CATALYST)
    cdna=g.product("cdna","cell-barcoded, UMI-tagged transcript cDNA",MaterialKind.MOLECULAR,state=molecule("CITE-seq cDNA",(role("cell barcode"),role("UMI"),role("mRNA cDNA"))))
    copied=g.product("copied_adt","cell-barcoded, UMI-tagged antibody-tag copies",MaterialKind.MOLECULAR,state=molecule("CITE-seq ADT copy",(role("cell barcode"),role("UMI"),role("antibody barcode"))))
    g.react("reverse_transcribe","Extend the bead primer on mRNA and on the DNA antibody tags in each GEM",StepKind.REVERSE_TRANSCRIPTION,((occupied,InputRole.TEMPLATE),(bead,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(cdna,copied))
    pool=g.input("cdna_primers","10x v2 cDNA amplification primers",MaterialKind.OLIGO)
    pool_r=g.input("cdna_primers_reverse","10x v2 cDNA amplification primers (second member)",MaterialKind.OLIGO,inferred=True)
    amp_pol=g.input("cdna_polymerase","10x cDNA amplification master mix",MaterialKind.CATALYST)
    amplified=g.product("amplified_pool","amplified cDNA with antibody-derived tags",MaterialKind.MOLECULAR,state=molecule("amplified CITE-seq products",(role("cDNA >300 bp"),role("ADT ~180 bp"))))
    g.react("amplify_cdna","Amplify cDNA for 10 cycles (manufacturer workflow)",StepKind.PCR,((cdna,InputRole.TEMPLATE),(copied,InputRole.TEMPLATE),(pool,InputRole.STRUCTURAL),(pool_r,InputRole.STRUCTURAL),(amp_pol,InputRole.CATALYST)),(amplified,),geometry=PCRGeometry.unavailable(GEX,reason="The 10x v2 cDNA primers are vendor reagents not printed in the paper; how the antibody-tag copies behave in this PCR (no ADT-specific primer is added in the 2017 protocol) is not described."))
    spri=g.input("spri_06x","0.6X SPRI",MaterialKind.CONDITION)
    rna=g.product("rna_branch","bead-bound cDNA >300 bp",MaterialKind.MOLECULAR,state=cdna.state)
    adt=g.product("adt_branch","supernatant ADTs ~180 bp",MaterialKind.MOLECULAR,state=copied.state)
    g.react("split","Separate cDNA (beads) from ADTs (supernatant) with 0.6X SPRI",StepKind.SPLIT,((amplified,InputRole.TEMPLATE),(spri,InputRole.CONDITION)),(rna,adt),changes_molecule=False)
    clean=g.input("spri_2x","two sequential 2X SPRI purifications",MaterialKind.CONDITION)
    adt_clean=g.product("adt_clean","purified ADTs",MaterialKind.MOLECULAR,state=copied.state)
    carry=g.product("primer_carryover","cDNA-amplification primer carryover",MaterialKind.WASTE)
    g.react("purify_adt","Purify the ADT supernatant twice with 2X SPRI",StepKind.CLEANUP,((adt,InputRole.TEMPLATE),(clean,InputRole.CONDITION)),(adt_clean,carry))
    g.discard(carry,"reduced by sequential 2X SPRI so that a ~140 bp TSO-RT-oligo product is not amplified in the ADT PCR (author protocol Figure 3a)")
    reviewed=["cells","fc_block","antibodies","stain","unbound","droplets","partition","occupied","failed","bead_primer","rt","reverse_transcribe","cdna","copied_adt","cdna_primers","cdna_primers_reverse","cdna_polymerase","amplify_cdna","amplified_pool","spri_06x","split","rna_branch","adt_branch","spri_2x","purify_adt","adt_clean","primer_carryover"]
    finish_pcr_branch(g,"gex",rna,GEX,GEX_P,reviewed,label="CITE-seq transcript",inferred_primers=True,read_lengths=READ_LENGTHS["Final transcript library"],limitation="The paper prepares the v2 transcript library according to the manufacturer's instructions; the vendor fragmentation, ligation and sample-index PCR are not redrawn here.")
    finish_pcr_branch(g,"adt",adt_clean,ADT,ADT_P,reviewed,label="CITE-seq antibody-derived tag",read_lengths=READ_LENGTHS["Final antibody-derived-tag library"],read_uses=READ_USES["Final antibody-derived-tag library"],step_label="Amplify the ADT library with PE 1.0 P5 and an i7-indexed RPIx primer (8–12 cycles)",forward_label="Illumina PE 1.0 P5 oligo (3′ phosphorothioates)",reverse_label="TruSeq Small RNA RPIx primer containing the i7 index (RPI1 shown)",polymerase_label="2x KAPA HiFi PCR Master Mix",forward_state=Construct([Segment("Illumina PE 1.0 P5 oligo",PE1_P5.replace("*",""),"p5")],name="Illumina PE 1.0 P5 oligo"),reverse_state=Construct([Segment("Illumina Small RNA RPI1",RPI1,"p7")],name="TruSeq Small RNA RPI1"),primer_source="Stoeckius et al. author protocol (MOESM3, 2017-07-11), 'Amplify ADT sequencing library' and 'Oligos required for ADT library amplification'",limitation="KAPA HiFi, PE 1.0 P5 and TruSeq Small RNA RPIx, 8–12 cycles, then 1.6X SPRI (author protocol p.4–5; Methods: 10 additional cycles); the double-stranded state of the pre-PCR ADT copy is not described, so primer landing is not drawn on it.")
    g.declare_reviewed_scope(reviewed,SOURCE)
    return g.require_publishable()
