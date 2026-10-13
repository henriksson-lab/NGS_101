"""ASAP-seq molecular model."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from antibody_tags import bridge_extension_scene, feature_library
from batch_ngs import nextera_library, seg
from chemdraw import Row, feature
from multimodal_spatial import finish_pcr_branch, modality_split
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind
from protocol_model import molecule, role

TITLE = "ASAP-seq — ATAC and surface-protein profiling"
NOTES = "01_asap-seq.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1038/s41587-021-00927-2">Mimitou et al., <i>Nature Biotechnology</i> (2021)</a>.'
SUMMARY = "Antibody-derived tags are extended with an ATAC-bead capture handle before a 10x single-cell ATAC workflow yields separate chromatin-accessibility and protein-tag libraries."
CAVEAT = "TotalSeq tag and bridge sequences vary by reagent generation. Their published molecular roles are shown without asserting one product's proprietary bases."
CELL=feature("asap_cell","cell_barcode","whitelist",whitelist="10x Chromium ATAC bead barcode")
ATAC, ATAC_P = nextera_library([seg("cell barcode association","C"*16,"cbc",placeholder=True,feature=CELL),seg("accessible genomic DNA", "X" * 42, placeholder=True)], "ASAP-seq ATAC library")
ADT, ADT_P = feature_library("ASAP-seq antibody-tag library", barcode_id="asap_antibody",cell_id="asap_cell")
FINAL_LIBRARIES = (
    ("ATAC library", ATAC, ATAC_P, "10x scATAC bead barcoding and indexed PCR complete the tagmented chromatin library.", "Paired genomic reads use the standard Nextera primer sites."),
    ("Antibody-tag library", ADT, ADT_P, "The extended TotalSeq tag is amplified as a separate Illumina library.", "Read cycles recover antibody barcode, UMI and sample indexes."),
)
READ_LENGTHS = {"ATAC library": {"Read 1": 50, "Read 2": 50}, "Antibody-tag library": {"Read 1": 28, "Read 2": 50}}

def sections():
    return [
        ("Stain fixed, permeabilized cells with DNA-barcoded antibodies", [Row(chunks=[("antibody — feature barcode — UMI — bridge tail", "cbc", False)])], "Surface or intracellular epitopes retain their antibody-derived DNA tags."),
        ("Extend antibody tags across a blocked bridge", bridge_extension_scene().rows(), "The bridge supplies sequence complementary to the scATAC bead oligo; its 3-prime block prevents bridge self-extension."),
        ("Tagment accessible chromatin and partition nuclei", [Row(chunks=[("Tn5 → accessible chromatin fragments; GEM bead → shared cell identity", "me", False)])], "Chromatin and antibody tags from the same cell receive the same droplet barcode."),
        ("Separate and amplify modalities", modality_split("ATAC library", "antibody-tag library"), "Distinct PCRs recover genomic inserts and antibody-derived tags."),
    ]

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    cells=g.input("cells","fixed permeabilized cells containing accessible chromatin",MaterialKind.BIOLOGICAL,state=molecule("ASAP-seq cells",(role("surface and intracellular epitopes"),role("accessible chromatin"))),primary=True)
    antibodies=g.input("antibodies","DNA-barcoded antibodies",MaterialKind.SELECTION)
    stained=g.product("stained","cells retaining antibody-derived tags",MaterialKind.MOLECULAR,state=molecule("stained cells",(role("antibody feature tags"),role("accessible chromatin"))))
    unbound=g.product("unbound","unbound antibody tags",MaterialKind.WASTE)
    g.react("stain","Stain targets with DNA-barcoded antibodies",StepKind.BINDING,((cells,InputRole.TEMPLATE),(antibodies,InputRole.SELECTION)),(stained,unbound))
    g.discard(unbound,"does not report a cell-associated target")
    bridge=g.input("bridge","blocked ATAC-bead capture bridge oligo",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (sequence varies by reagent generation)")
    bridged=g.product("bridged","antibody tags hybridized to blocked bridges",MaterialKind.MOLECULAR,state=stained.state)
    g.react("hybridize_bridge","Hybridize blocked bridges to antibody tags",StepKind.HYBRIDIZATION,((stained,InputRole.TEMPLATE),(bridge,InputRole.STRUCTURAL)),(bridged,))
    polymerase=g.input("bridge_polymerase","bridge-extension polymerase",MaterialKind.CATALYST)
    extended=g.product("extended","antibody tags extended with ATAC-bead capture handles",MaterialKind.MOLECULAR,state=molecule("bridge-extended cell",(role("capturable antibody tag"),role("accessible chromatin"))))
    g.react("extend_bridge","Extend antibody tags across the blocked bridge",StepKind.EXTENSION,((bridged,InputRole.TEMPLATE),(polymerase,InputRole.CATALYST)),(extended,))
    transposome=g.input("transposome","10x scATAC loaded transposome",MaterialKind.ADAPTER)
    tn5=g.input("tn5","Tn5 transposase",MaterialKind.CATALYST)
    tagged=g.product("tagged","cells retaining tagmented chromatin and extended antibody tags",MaterialKind.MOLECULAR,state=extended.state)
    tn5_waste=g.product("tn5_waste","non-library tagmentation products",MaterialKind.WASTE)
    g.react("tagment","Tagment accessible chromatin",StepKind.TAGMENTATION,((extended,InputRole.TEMPLATE),(transposome,InputRole.STRUCTURAL),(tn5,InputRole.CATALYST)),(tagged,tn5_waste))
    g.discard(tn5_waste,"does not retain an amplifiable accessibility fragment")
    droplets=g.input("droplets","10x scATAC GEM partitioning condition",MaterialKind.CONDITION)
    occupied=g.product("occupied","cell-barcoded ATAC and antibody-tag products",MaterialKind.MOLECULAR,state=tagged.state)
    failed=g.product("failed","empty or unresolved GEM products",MaterialKind.WASTE)
    g.react("partition","Partition nuclei and transfer the shared cell identity",StepKind.PARTITION,((tagged,InputRole.TEMPLATE),(droplets,InputRole.CONDITION)),(occupied,failed),changes_molecule=False)
    g.discard(failed,"does not retain a resolvable cell barcode")
    atac=g.product("atac_branch","cell-barcoded ATAC branch",MaterialKind.MOLECULAR,state=occupied.state)
    adt=g.product("adt_branch","cell-barcoded antibody-tag branch",MaterialKind.MOLECULAR,state=occupied.state)
    g.react("split","Separate ATAC and antibody-tag products",StepKind.SPLIT,((occupied,InputRole.TEMPLATE),),(atac,adt),changes_molecule=False)
    reviewed=["cells","antibodies","stain","bridge","hybridize_bridge","bridge_polymerase","extend_bridge","transposome","tn5","tagment","droplets","partition","split"]
    finish_pcr_branch(g,"atac",atac,ATAC,ATAC_P,reviewed,label="ASAP-seq ATAC",read_lengths=READ_LENGTHS["ATAC library"])
    finish_pcr_branch(g,"adt",adt,ADT,ADT_P,reviewed,label="ASAP-seq antibody tag",read_lengths=READ_LENGTHS["Antibody-tag library"])
    g.declare_reviewed_scope(reviewed,SOURCE)
    return g.require_publishable()
