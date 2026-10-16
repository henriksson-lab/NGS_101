"""ECCITE-seq 5-prime multimodal tag and direct-guide capture."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from batch_ngs import seg,truseq_library
from chemdraw import Row, feature
from multimodal_spatial import finish_pcr_branch, modality_split
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind
from protocol_model import molecule, role

TITLE="ECCITE-seq"
NOTES="01_eccite-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/s41592-019-0392-0">Mimitou et al., <i>Nature Methods</i> (2019)</a>.'
SUMMARY="ECCITE-seq adapts 10x 5′ chemistry so transcript, antibody tag, sample hashtag and directly captured sgRNA products receive a common cell barcode and UMI. Modality-specific handles then support separate libraries, alongside standard immune-receptor enrichment when desired."
CAVEAT=("The schematic separates transcript, antibody tag, sample hashtag and guide-derived libraries. "
        "Supplementary Table 3 in the publisher's MOESM1 prints 49 antibody barcodes "
        "of 15 bases and three of 12 bases (lines 219–387); Tables 4–5 print 12-base "
        "mouse/human hashtags (lines 389–460). Representative examples of all three "
        "sizes are drawn here. "
        "The MOESM3 oligo examples (lines 385–396) instead show eight <code>x</code> "
        "characters despite calling the example barcode 12 nt; the identity of those "
        "eight positions is unresolved. They are not substituted for the Table 3–5 barcodes. "
        "The exact barcode examples are source-transcribed, while the final-library "
        "handle placeholders show branch identity rather than a verified end-to-end "
        "sequence. Immune-receptor enrichment is the inherited 10x V(D)J branch and "
        "is linked in the notes.")
CELL=feature("cell_barcode","cell_barcode","whitelist",whitelist="10x-chromium-5prime")
UMI=feature("umi","umi","random")
ANTIBODY_FEATURE=feature("antibody_barcode","feature_barcode","whitelist",whitelist="ECCITE-seq Supplementary Table 3")
HASHTAG_FEATURE=feature("hashtag_barcode","sample_index","whitelist",whitelist="ECCITE-seq Supplementary Tables 4–5")
GUIDE_FEATURE=feature("guide","guide_barcode","unknown")
ANTIBODY_LONG_BARCODE="TAGAACCGTATCCGT"
ANTIBODY_SHORT_BARCODE="AATAGCGGAGCC"
HASHTAG_A_BARCODE="AGGACCATCCAA"

def library(name,feature):
    return truseq_library([seg("cell barcode","B"*16,"cbc",placeholder=True,feature=CELL),seg("UMI","U"*10,"umi",placeholder=True,feature=UMI),*feature],name,dual_index=False)
GEX,GEX_P=library("ECCITE-seq transcript library",[seg("5-prime cDNA","X"*36,placeholder=True)])
ADT,ADT_P=library("ECCITE-seq 15-base antibody-tag example",[seg("protein-tag handle","H"*18,placeholder=True),seg("antibody barcode, Table 3 first entry",ANTIBODY_LONG_BARCODE,"cbc",feature=ANTIBODY_FEATURE)])
ADT_SHORT,ADT_SHORT_P=library("ECCITE-seq 12-base antibody-tag example",[seg("protein-tag handle","H"*18,placeholder=True),seg("antibody barcode, Table 3 short entry",ANTIBODY_SHORT_BARCODE,"cbc",feature=ANTIBODY_FEATURE)])
HTO,HTO_P=library("ECCITE-seq sample-hashtag example",[seg("hashtag handle","H"*18,placeholder=True),seg("Hashtag_A barcode",HASHTAG_A_BARCODE,"cbc",feature=HASHTAG_FEATURE)])
GDO,GDO_P=library("ECCITE-seq guide-derived oligo library",[seg("sgRNA protospacer","G"*20,"cbc",placeholder=True,feature=GUIDE_FEATURE),seg("guide scaffold capture region","S"*18,placeholder=True)])
FINAL_LIBRARIES=(
    ("Final transcript library",GEX,GEX_P,"5′ transcript sequence follows cell barcode and UMI.","The transcript library's primer geometry is computed on the finished duplex."),
    ("Final antibody-tag library, 15-base example",ADT,ADT_P,"The first antibody barcode in Supplementary Table 3 is TAGAACCGTATCCGT; most entries have 15 bases.","The antibody-tag library's primer geometry is computed on the finished duplex."),
    ("Final antibody-tag library, 12-base example",ADT_SHORT,ADT_SHORT_P,"The last three Supplementary Table 3 antibody barcodes have 12 bases; AATAGCGGAGCC is one example.","This shorter antibody-tag variant uses the same library preparation."),
    ("Final sample-hashtag library",HTO,HTO_P,"Supplementary Table 4 prints mouse Hashtag_A as AGGACCATCCAA; Tables 4–5 use 12-base barcodes.","The hashtag library's primer geometry is computed on the finished duplex."),
    ("Final guide-derived oligo library",GDO,GDO_P,"Direct guide capture reports sgRNA identity with the same cell barcode and UMI.","The guide library's primer geometry is computed on the finished duplex."),)

def sections():
    return [
        ("Prime tags and sgRNAs during 5′ reverse transcription",
         [Row(chunks=[("antibody/hashtag oligo ⇄ bead-TSO-complementary end → cell barcode + UMI", "cbc", False)]),
          Row(chunks=[("sgRNA protospacer—scaffold ⇄ scaffold RT primer → template switch → cell barcode + UMI", "umi", False)])],
         "An antibody tag anneals to the bead-associated template-switch system; a separate RT primer binds the sgRNA scaffold, copies the guide and then template-switches."),
        ("Amplify all barcoded products, then enrich by handle",
         modality_split("5′ gene-expression / optional V(D)J","antibody-derived tags","sample hashtags","guide-derived oligos"),
         "A shared cell identity is installed before modality-specific enrichment. Different handles keep the short tag libraries distinct."),]

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    cells=g.input("cells","cells carrying poly(A) RNA, antibody/hashtag oligos and sgRNA",MaterialKind.BIOLOGICAL,state=molecule("ECCITE-seq cells",(role("poly(A) RNA"),role("antibody and hashtag oligos"),role("sgRNA"))),primary=True)
    droplets=g.input("droplets","10x 5-prime GEM partitioning condition",MaterialKind.CONDITION)
    occupied=g.product("occupied","occupied GEM products",MaterialKind.MOLECULAR,state=cells.state)
    failed=g.product("failed","empty or unresolved GEM products",MaterialKind.WASTE)
    g.react("partition","Partition cells with barcoded gel beads",StepKind.PARTITION,((cells,InputRole.TEMPLATE),(droplets,InputRole.CONDITION)),(occupied,failed),changes_molecule=False)
    g.discard(failed,"does not retain a resolvable cell barcode")
    transcript=g.product("transcript_branch","5-prime transcript branch",MaterialKind.MOLECULAR,state=occupied.state)
    tag=g.product("tag_branch","antibody and hashtag branch",MaterialKind.MOLECULAR,state=occupied.state)
    guide=g.product("guide_branch","direct sgRNA-capture branch",MaterialKind.MOLECULAR,state=occupied.state)
    g.react("split","Separate transcript, tag and guide-derived products",StepKind.SPLIT,((occupied,InputRole.TEMPLATE),),(transcript,tag,guide),changes_molecule=False)
    bead=g.input("bead_oligo","cell-barcode/UMI template-switch bead oligo",MaterialKind.OLIGO)
    rt=g.input("rt","template-switching reverse transcriptase",MaterialKind.CATALYST)
    cdna=g.product("cdna","cell-barcoded 5-prime cDNA",MaterialKind.MOLECULAR,state=molecule("ECCITE transcript cDNA",(role("cell barcode"),role("UMI"),role("5-prime cDNA"))))
    g.react("transcript_rt","Reverse-transcribe and template-switch transcript RNA",StepKind.REVERSE_TRANSCRIPTION,((transcript,InputRole.TEMPLATE),(bead,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(cdna,))
    tag_primer=g.input("tag_primer","tag-to-bead capture primer",MaterialKind.OLIGO)
    tag_pol=g.input("tag_extension_polymerase","tag barcode-transfer polymerase",MaterialKind.CATALYST)
    copied_tag=g.product("copied_tag","cell-barcoded antibody/hashtag products",MaterialKind.MOLECULAR,state=molecule("ECCITE tag precursor",(role("cell barcode"),role("UMI"),role("antibody or hashtag barcode"))))
    g.react("copy_tags","Copy antibody and hashtag oligos with the shared cell identity",StepKind.EXTENSION,((tag,InputRole.TEMPLATE),(tag_primer,InputRole.STRUCTURAL),(tag_pol,InputRole.CATALYST)),(copied_tag,))
    copied_adt=g.product("copied_adt","cell-barcoded antibody-tag products",MaterialKind.MOLECULAR,
                         state=molecule("ECCITE antibody-tag precursor",(role("cell barcode"),role("UMI"),role("antibody barcode"))))
    copied_hto=g.product("copied_hto","cell-barcoded sample-hashtag products",MaterialKind.MOLECULAR,
                         state=molecule("ECCITE hashtag precursor",(role("cell barcode"),role("UMI"),role("12-base hashtag barcode"))))
    g.react("split_tag_classes","Separate antibody tags from sample hashtags by their handles",StepKind.SPLIT,
            ((copied_tag,InputRole.TEMPLATE),),(copied_adt,copied_hto),changes_molecule=False,
            source=SOURCE+" Supplementary Tables 3–5 and Supplementary Protocol tag-library PCR")
    copied_adt_long=g.product("copied_adt_long","antibody-tag products with 15-base barcodes",MaterialKind.MOLECULAR,state=copied_adt.state)
    copied_adt_short=g.product("copied_adt_short","antibody-tag products with 12-base barcodes",MaterialKind.MOLECULAR,state=copied_adt.state)
    g.react("split_antibody_lengths","Resolve published antibody barcode length variants",StepKind.SPLIT,
            ((copied_adt,InputRole.TEMPLATE),),(copied_adt_long,copied_adt_short),changes_molecule=False,
            source=SOURCE+" Supplementary Table 3, final three 12-base barcodes")
    guide_primer=g.input("guide_rt_primer","sgRNA-scaffold RT primer AGCAAGTGAGAAGCATCGTGTCAAAGCACCGACTCGGTGCCAC",MaterialKind.OLIGO,source=SOURCE+" author supplementary protocol p. 5")
    guide_cdna=g.product("guide_cdna","template-switched guide-derived cDNA",MaterialKind.MOLECULAR,state=molecule("ECCITE guide cDNA",(role("cell barcode"),role("UMI"),role("guide protospacer"))))
    g.react("guide_rt","Reverse-transcribe sgRNA from its scaffold primer",StepKind.REVERSE_TRANSCRIPTION,((guide,InputRole.TEMPLATE),(guide_primer,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(guide_cdna,))
    reviewed=["cells","droplets","partition","split","bead_oligo","rt","transcript_rt","tag_primer","tag_extension_polymerase","copy_tags","split_tag_classes","split_antibody_lengths","guide_rt_primer","guide_rt"]
    finish_pcr_branch(g,"gex",cdna,GEX,GEX_P,reviewed,label="ECCITE-seq transcript")
    finish_pcr_branch(g,"adt",copied_adt_long,ADT,ADT_P,reviewed,label="ECCITE-seq 15-base antibody-tag")
    finish_pcr_branch(g,"adt_short",copied_adt_short,ADT_SHORT,ADT_SHORT_P,reviewed,label="ECCITE-seq 12-base antibody-tag")
    finish_pcr_branch(g,"hto",copied_hto,HTO,HTO_P,reviewed,label="ECCITE-seq sample-hashtag")
    finish_pcr_branch(g,"guide",guide_cdna,GDO,GDO_P,reviewed,label="ECCITE-seq guide-derived oligo")
    g.declare_reviewed_scope(reviewed,SOURCE)
    return g.require_publishable()
