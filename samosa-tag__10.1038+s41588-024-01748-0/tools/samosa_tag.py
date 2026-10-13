"""SAMOSA-Tag molecular model."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from batch_ngs import joined_scene, seg
from chemdraw import Row
from longread_footprinting import hairpin_tagmented_smrtbell, marked_fiber, pacbio_entry, smrtbell_rows
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)

TITLE = "SAMOSA-Tag — hairpin tagmentation of marked chromatin"
NOTES = "01_samosa-tag.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1038/s41588-024-01748-0">Dennis et al., <i>Nature Genetics</i> (2024)</a>.'
SUMMARY = "EcoGII first records chromatin accessibility as m6dA; hairpin-loaded Tn5 then fragments the marked fibers while installing PacBio-compatible closed ends."
CAVEAT = "Hairpin-Tn5 adapters are published as protocol oligos but are drawn by molecular role here. The final molecule is shown only after gap repair has made a closed SMRTbell."
SMRTBELL = hairpin_tagmented_smrtbell()
FINAL_LIBRARY = None
SEQ_PRIMERS = ()
SEQUENCING_ENDING = pacbio_entry(SMRTBELL)
INITIAL_ROWS=tuple(marked_fiber("EcoGII"))
INITIAL_NAME="Native chromatin fiber"

def sections():
    return [
        ("Mark accessible adenines in nuclei", marked_fiber("EcoGII"),
         "EcoGII and SAM deposit m6dA on exposed DNA while protein-bound intervals remain footprints."),
        ("Tagment with hairpin-loaded Tn5", joined_scene([seg("left hairpin-ME","X"*20,"me",placeholder=True),seg("marked chromatin fragment","X"*36,placeholder=True),seg("right ME-hairpin","X"*20,"me",placeholder=True)],(("left hairpin-ME","marked chromatin fragment","Tn5 transfer"),("marked chromatin fragment","right ME-hairpin","Tn5 transfer")),label="hairpin-tagmented chromatin").rows(),
         "Tn5 couples fragmentation to transfer of a hairpin-bearing mosaic-end adapter at both boundaries."),
        ("Gap-repair the tagmented product", smrtbell_rows(SMRTBELL),
         "Polymerase fills the transposition gaps, leaving one covalently closed PacBio template."),
    ]


def protocol_graph():
    graph=ProtocolGraph(TITLE,authority=SOURCE)
    chromatin=graph.input("chromatin","native chromatin fibers",MaterialKind.BIOLOGICAL,
                          primary=True)
    ecogii=graph.input("ecogii","EcoGII adenine methyltransferase",MaterialKind.CATALYST)
    sam=graph.input("sam","S-adenosylmethionine methyl donor",MaterialKind.CONDITION)
    marked=graph.product("marked_fibers","m6dA-marked chromatin fibers",MaterialKind.MOLECULAR,
                         state=INITIAL_ROWS)
    graph.react("mark_accessibility","Mark accessible adenines in nuclei",
                StepKind.CHEMICAL_CONVERSION,
                ((chromatin,InputRole.TEMPLATE),(ecogii,InputRole.CATALYST),
                 (sam,InputRole.CONDITION)),(marked,))
    hairpin_adapter=graph.input("hairpin_transposome","hairpin-loaded mosaic-end adapter",
                                MaterialKind.ADAPTER)
    tn5=graph.input("tn5","Tn5 transposase",MaterialKind.CATALYST)
    tagmented=graph.product("hairpin_tagmented","hairpin-tagmented marked chromatin",
                            MaterialKind.MOLECULAR)
    waste=graph.product("tagmentation_waste","non-SMRTbell tagmentation products",
                        MaterialKind.WASTE)
    graph.react("hairpin_tagment","Fragment and transfer hairpin adapters",
                StepKind.TAGMENTATION,
                ((marked,InputRole.TEMPLATE),(hairpin_adapter,InputRole.STRUCTURAL),
                 (tn5,InputRole.CATALYST)),(tagmented,waste))
    graph.discard(waste,"products lacking the required two hairpin-bearing transfer ends")
    repair=graph.input("gap_repair_polymerase","transposition-gap repair polymerase",
                       MaterialKind.CATALYST)
    library=graph.product("library",SMRTBELL.name,MaterialKind.LIBRARY,state=SMRTBELL)
    graph.react("gap_repair","Repair both transposition gaps to close the SMRTbell",
                StepKind.EXTENSION,
                ((tagmented,InputRole.TEMPLATE),(repair,InputRole.CATALYST)),(library,))
    graph.sequence(library,SequencingPlan(
        limitation="PacBio hairpin adapter bases and primer landing sequence are proprietary; the closed SMRTbell traversal is known."))
    graph.declare_reviewed_scope(("chromatin","ecogii","sam","mark_accessibility",
                                  "hairpin_transposome","tn5","hairpin_tagment",
                                  "gap_repair_polymerase","gap_repair","library"),SOURCE)
    return graph.require_publishable()
