"""Original CITE-seq with 10x 3-prime transcriptome and ADT readout."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from batch_ngs import seg, truseq_library
from chemdraw import Row, feature
from multimodal_spatial import finish_pcr_branch, modality_split
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind
from protocol_model import molecule, role

TITLE="CITE-seq"
NOTES="01_cite-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/nmeth.4380">Stoeckius et al., <i>Nature Methods</i> (2017)</a>.'
SUMMARY="Antibodies carry poly(A)-tailed DNA tags containing an antibody identity and an amplification handle. The same oligo-dT bead primer that captures mRNA also copies these antibody-derived tags, transferring the droplet cell barcode and UMI before ADT and transcript libraries are separated by size."
CAVEAT="The original paper demonstrates both Drop-seq and 10x implementations. This page follows its 10x 3′ implementation; antibody barcode identities are placeholders, not a claim that all panels use one fixed sequence."
CELL=feature("cell_barcode","cell_barcode","whitelist",whitelist="10x-chromium-3prime")
UMI=feature("umi","umi","random")
ADT_FEATURE=feature("antibody_barcode","feature_barcode","unknown")

GEX,GEX_P=truseq_library([seg("cell barcode","B"*16,"cbc",placeholder=True,feature=CELL),seg("UMI","U"*10,"umi",placeholder=True,feature=UMI),seg("poly(dT) junction","T"*10,placeholder=True),seg("mRNA cDNA","X"*34,placeholder=True)],"CITE-seq transcript library",dual_index=False)
ADT,ADT_P=truseq_library([seg("cell barcode","B"*16,"cbc",placeholder=True,feature=CELL),seg("UMI","U"*10,"umi",placeholder=True,feature=UMI),seg("ADT PCR handle","H"*18,"r2",placeholder=True),seg("antibody barcode","A"*15,"cbc",placeholder=True,feature=ADT_FEATURE)],"CITE-seq ADT library",dual_index=False)
FINAL_LIBRARIES=(
    ("Final transcript library",GEX,GEX_P,"The 10x cell barcode and UMI identify each transcript molecule.","The final transcript duplex shows the declared Illumina primer geometry."),
    ("Final antibody-derived-tag library",ADT,ADT_P,"The same cell barcode and UMI are followed by the antibody-specific tag.","The final ADT duplex shows the declared Illumina primer geometry."),)
READ_LENGTHS={
    "Final transcript library":{"Read 1":26,"Index 1 (i7)":8},
    "Final antibody-derived-tag library":{"Read 1":26,"Index 1 (i7)":8},
}
READ_USES={
    "Final antibody-derived-tag library":{"Read 2":"extract"},
}
SHOW_READ_CONTENT=True

def sections():
    return [
        ("Label surface proteins with poly(A)-tailed DNA tags",
         [Row(chunks=[("antibody—[ADT PCR handle][antibody barcode][poly(A)]", "cbc", False)])],
         "Each antibody species is associated with a distinct DNA barcode; the poly(A) tail makes the tag compatible with an oligo-dT single-cell RNA workflow."),
        ("Copy mRNA and antibody tags in the same droplet",
         [Row(chunks=[("bead [cell barcode][UMI]—oligo(dT) || mRNA poly(A)", "umi", False)]),
          Row(chunks=[("bead [cell barcode][UMI]—oligo(dT) || ADT poly(A)", "cbc", False)])],
         "Reverse transcriptase copies RNA and extends on the DNA ADT, transferring the same droplet barcode system to both modalities."),
        ("Separate by size and amplify independently",modality_split("large transcript-derived cDNA","short antibody-derived tags"),
         "SPRI size separation sends mRNA cDNA through the standard gene-expression branch and ADTs through tag-specific PCR."),]

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    cells=g.input("cells","cells containing poly(A) RNA",MaterialKind.BIOLOGICAL,state=molecule("cells",(role("surface proteins"),role("poly(A) RNA"))),primary=True)
    antibodies=g.input("antibodies","poly(A)-tailed DNA-barcoded antibodies",MaterialKind.SELECTION)
    stained=g.product("stained","cells bearing antibody-derived DNA tags",MaterialKind.MOLECULAR,state=molecule("CITE-seq stained cells",(role("poly(A) RNA"),role("poly(A)-tailed antibody tag"))))
    unbound=g.product("unbound","unbound antibody tags",MaterialKind.WASTE)
    g.react("stain","Label surface proteins with poly(A)-tailed DNA tags",StepKind.BINDING,((cells,InputRole.TEMPLATE),(antibodies,InputRole.SELECTION)),(stained,unbound))
    g.discard(unbound,"does not report a cell-associated surface protein")
    droplets=g.input("droplets","10x 3-prime droplet partitioning condition",MaterialKind.CONDITION)
    occupied=g.product("occupied","occupied GEM products",MaterialKind.MOLECULAR,state=stained.state)
    failed=g.product("failed","empty or unresolved droplet products",MaterialKind.WASTE)
    g.react("partition","Partition stained cells with barcoded beads",StepKind.PARTITION,((stained,InputRole.TEMPLATE),(droplets,InputRole.CONDITION)),(occupied,failed),changes_molecule=False)
    g.discard(failed,"does not retain a resolvable cell barcode")
    rna=g.product("rna_branch","poly(A)-RNA capture branch",MaterialKind.MOLECULAR,state=occupied.state)
    adt=g.product("adt_branch","poly(A)-tailed antibody-tag branch",MaterialKind.MOLECULAR,state=occupied.state)
    g.react("split","Separate transcript and antibody-tag products by size",StepKind.SPLIT,((occupied,InputRole.TEMPLATE),),(rna,adt),changes_molecule=False)
    bead=g.input("bead_primer","cell-barcode/UMI oligo(dT) bead primer",MaterialKind.OLIGO)
    rt=g.input("rt","reverse transcriptase",MaterialKind.CATALYST)
    cdna=g.product("cdna","cell-barcoded, UMI-tagged transcript cDNA",MaterialKind.MOLECULAR,state=molecule("CITE-seq cDNA",(role("cell barcode"),role("UMI"),role("mRNA cDNA"))))
    g.react("reverse_transcribe","Reverse-transcribe mRNA on the barcoded bead primer",StepKind.REVERSE_TRANSCRIPTION,((rna,InputRole.TEMPLATE),(bead,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(cdna,))
    pol=g.input("adt_extension_polymerase","antibody-tag barcode-transfer polymerase",MaterialKind.CATALYST)
    copied=g.product("copied_adt","cell-barcoded, UMI-tagged antibody tags",MaterialKind.MOLECULAR,state=molecule("CITE-seq ADT precursor",(role("cell barcode"),role("UMI"),role("antibody barcode"))))
    g.react("copy_adt","Extend the same bead-primer system across antibody tags",StepKind.EXTENSION,((adt,InputRole.TEMPLATE),(bead,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(copied,))
    reviewed=["cells","antibodies","stain","droplets","partition","split","bead_primer","rt","reverse_transcribe","adt_extension_polymerase","copy_adt"]
    finish_pcr_branch(g,"gex",cdna,GEX,GEX_P,reviewed,label="CITE-seq transcript",read_lengths=READ_LENGTHS["Final transcript library"])
    finish_pcr_branch(g,"adt",copied,ADT,ADT_P,reviewed,label="CITE-seq antibody-derived tag",read_lengths=READ_LENGTHS["Final antibody-derived-tag library"],read_uses=READ_USES["Final antibody-derived-tag library"])
    g.declare_reviewed_scope(reviewed,SOURCE)
    return g.require_publishable()
