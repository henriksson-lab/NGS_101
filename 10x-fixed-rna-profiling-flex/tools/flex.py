"""10x Chromium Fixed RNA Profiling / Flex probe-pair chemistry."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from batch_ngs import seg, truseq_library
from chemdraw import Row, feature
from multimodal_spatial import adjacent_probe_scene, finish_pcr_branch
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind
from protocol_model import molecule, role

TITLE="10x Fixed RNA Profiling (Flex)"
NOTES="01_flex.html"
SOURCE='Commercial protocol: 10x Genomics <a href="https://www.10xgenomics.com/support/cn/flex-gene-expression/documentation/steps/library-prep/chromium-fixed-rna-profiling-reagent-kits-for-singleplexed-samples-with-feature-barcode-technology-for-protein-using-barcode-oligo-capture">Fixed RNA Profiling user guide CG000674 Rev C</a>.'
SUMMARY="Pairs of gene-specific probes hybridize adjacently on formaldehyde-fixed RNA and are ligated. The ligated probe product—not the RNA itself—is then given a GEM cell barcode and UMI and converted into a sequencing library."
CAVEAT="CG000674 Rev C Appendix p. 85 prints the gel-bead oligo and complete gene-expression library layout. The compact schematic still normalizes the platform arms and does not transcribe the printed probe-side constants or Small RNA Read-2 arm. Probe-panel target sequences remain commercial; the drawn generic probe region is inferred."

def oligos():
    return [
        "<b>Flex gel-bead primer</b> — 5′-CTACACGACGCTCTTCCGATCT-[16-base cell barcode]-[12-base UMI]-TTGCTAGGACCG-3′",
        "<b>Published gene-expression reporter junction</b> — 5′-TTGCTAGGACCG-[BC8]-NN-TACGTGCTAACCGCGT-[ligated probe insert]-TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC-3′",
    ]

LIB, PRIMERS = truseq_library([
    seg("cell barcode", "B"*16, "cbc", placeholder=True, inferred=True,
        feature=feature("cell_barcode", "cell_barcode", "whitelist", whitelist="10x Flex whitelist")),
    seg("UMI", "U"*12, "umi", placeholder=True, inferred=True,
        feature=feature("umi", "umi", "random")),
    seg("ligated probe-pair identifier", "X"*30, placeholder=True, inferred=True)],
    "Flex gene-expression library", inferred_adapters=True)
FINAL_LIBRARIES=(("Final probe-derived gene-expression library", LIB, PRIMERS,
                  "INFERRED — This compact final duplex preserves read roles but normalizes platform arms; the source-exact Flex reporter constants are transcribed in Key oligos and research notes.",
                  "The declared Illumina primer sites are shown on the completed probe-derived library."),)

def sections():
    return [
        ("Hybridize adjacent probe halves to fixed RNA", adjacent_probe_scene().rows(),
         "A valid target is represented by two gene-specific probe halves bound next to one another."),
        ("Ligate the adjacent probe pair",
         [Row(chunks=[("[left probe barcode] ** [right probe barcode]", "cbc", False)]),
          Row(chunks=[("                 ** ligation", "me", False)])],
         "Only adjacent, correctly hybridized probe halves form the amplifiable reporter molecule."),
        ("Partition, barcode and amplify",
         [Row(chunks=[("ligated reporter + GEM bead → [cell barcode][UMI][probe-pair identifier]", "cbc", False)]),
          Row(chunks=[("→ pre-amplification → sample-index PCR", None, False)])],
         "INFERRED — The cell barcode and UMI are added after probe ligation; this compact workflow does not draw the published gel-bead sequence base by base."),
    ]

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    rna=g.input("fixed_rna","formaldehyde-fixed target RNA",MaterialKind.BIOLOGICAL,state=molecule("fixed RNA",(role("left probe target"),role("right probe target"))),primary=True)
    left=g.input("left_probe","gene-specific left probe half",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (panel-specific commercial sequence)")
    right=g.input("right_probe","gene-specific right probe half",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (panel-specific commercial sequence)")
    adjacent=g.product("adjacent","adjacent probe halves hybridized on target RNA",MaterialKind.MOLECULAR,state=adjacent_probe_scene())
    g.react("hybridize","Hybridize adjacent probe halves to fixed RNA",StepKind.HYBRIDIZATION,((rna,InputRole.TEMPLATE),(left,InputRole.STRUCTURAL),(right,InputRole.STRUCTURAL)),(adjacent,))
    ligase=g.input("ligase","probe-pair ligase",MaterialKind.CATALYST)
    reporter=g.product("reporter","ligated probe-pair reporter",MaterialKind.MOLECULAR,state=molecule("ligated Flex reporter",(role("left probe barcode"),role("right probe barcode")),junctions=(("left probe barcode","right probe barcode","probe-pair ligation"),)))
    g.react("seal_probe_pair","Ligate the adjacent probe pair",StepKind.NICK_SEALING,((adjacent,InputRole.TEMPLATE),(ligase,InputRole.CATALYST)),(reporter,))
    droplets=g.input("droplets","Flex GEM partitioning condition",MaterialKind.CONDITION)
    occupied=g.product("occupied","partitioned ligated reporters",MaterialKind.MOLECULAR,state=reporter.state)
    failed=g.product("failed","empty or unresolved GEM material",MaterialKind.WASTE)
    g.react("partition","Partition ligated reporters with barcoded gel beads",StepKind.PARTITION,((reporter,InputRole.TEMPLATE),(droplets,InputRole.CONDITION)),(occupied,failed),changes_molecule=False)
    g.discard(failed,"does not retain a resolvable cell-associated reporter")
    bead=g.input("bead_oligo","cell-barcode/UMI gel-bead oligo CTACACGACGCTCTTCCGATCT-N16-N12-TTGCTAGGACCG",MaterialKind.OLIGO,source=SOURCE+" Appendix p. 85")
    pol=g.input("extension_polymerase","barcode-transfer polymerase",MaterialKind.CATALYST)
    barcoded=g.product("barcoded","cell-barcode/UMI-tagged probe reporters",MaterialKind.MOLECULAR,state=molecule("barcoded Flex reporter",(role("cell barcode"),role("UMI"),role("probe-pair identifier"))))
    g.react("barcode","Transfer cell barcode and UMI to the ligated reporter",StepKind.EXTENSION,((occupied,InputRole.TEMPLATE),(bead,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(barcoded,))
    reviewed=["fixed_rna","left_probe","right_probe","hybridize","ligase","seal_probe_pair","droplets","partition","bead_oligo","extension_polymerase","barcode"]
    finish_pcr_branch(g,"flex",barcoded,LIB,PRIMERS,reviewed,label="Flex probe-derived gene-expression",inferred_primers=True)
    g.declare_reviewed_scope(reviewed,SOURCE)
    return g.require_publishable()
