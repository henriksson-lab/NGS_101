"""Illumina Nextera Mate Pair library construction."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from adapter_ligation import reported_adapter_end
from circular import CircularizationGeometry
from chemdraw import Row
from chemdraw import Scene
from mate_pair import (DUPLICATE_JUNCTION, circularization_rows, recovered_truseq_library,
                       recovery_scene, tagmentation_scene)
from endprep import dA_tailed_scene, repair_and_dA_tail
from chemdraw import Construct, Segment
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)

TITLE = "Nextera Mate Pair — circular junction recovery"
NOTES = "01_nextera-mate-pair.html"
SOURCE = ('Commercial defining source: Illumina, <a href="https://support.illumina.com/'
          'content/dam/illumina-support/documents/documentation/chemistry_documentation/'
          'samplepreps_nextera/nexteramatepair/nextera-mate-pair-reference-guide-15035209-02.pdf">'
          '<i>Nextera Mate Pair Library Prep Reference Guide</i>, document 15035209 v02</a>; '
          '<a href="https://support.illumina.com/content/dam/illumina-marketing/documents/'
          'products/technotes/technote_nextera_matepair_data_processing.pdf">junction-sequence technical note</a>.')
SUMMARY = ("Long genomic fragments are Tn5-tagged, gap-filled and circularized. After shearing, "
           "streptavidin selects short fragments spanning the old circle-closing bond; TruSeq "
           "adapters turn those distant genomic ends into an outward-facing read pair.")
CAVEAT = "The drawn duplicated junction is the common product; the vendor also documents single-junction products."

FINAL_LIBRARY, SEQ_PRIMERS = recovered_truseq_library()
FINAL_CAPTION = ("Single-index TruSeq library containing a selected duplicated mate-pair junction. "
                 "The two genomic pieces flanking it were distant ends of one long input fragment.")
SEQUENCING_INTRO = ("Standard TruSeq Read 1, six-base Index 1 and Read 2 primers bind the final "
                    "library. The genomic mates map outward because circularization inverted their relationship.")


def adapter_ligation_scene() -> Scene:
    sc = Scene.duplex(list(FINAL_LIBRARY), label="adapter-ligated captured fragment")
    sc.junction("top", "Read 1 arm", "left distant end", "adapter ligation")
    sc.junction("top", "dA junction", "Index 1 / Read 2 arm", "adapter ligation")
    sc.labels("top")
    return sc


def sections():
    return [
        ("Tagment long genomic DNA", [*tagmentation_scene().rows(),
            Row(chunks=[("strand-displacement fill repairs the nine-base Tn5 gaps", None, False)]),
        ], "Mate Pair Tagment Enzyme fragments and tags long DNA; ** marks the transferred tag boundaries."),
        ("Circularize long fragments", circularization_rows(),
         "Blunt intramolecular ligation closes the two distant ends; exonuclease removes molecules that remain linear."),
        ("Shear and recover the old circle junction", recovery_scene().rows(),
         "Covaris shearing yields short pieces. Streptavidin beads retain the biotin-marked junction-spanning subset."),
        ("End-repair and add one 3-prime dA to each end",
         dA_tailed_scene(repair_and_dA_tail(Construct([Segment("captured mate-pair junction fragment", "X" * 46, placeholder=True)]))).rows(),
         "The captured junction fragment is shown after dA-tailing and before the index adapter is added."),
        ("Ligate TruSeq index adapter", adapter_ligation_scene().rows(),
         "The final library uses the TruSeq DNA LT single-index workflow, not Nextera sequencing primers."),
    ]


def protocol_graph():
    graph=ProtocolGraph(TITLE,authority=SOURCE)
    dna=graph.input("genomic_dna","long genomic DNA",MaterialKind.MOLECULAR,primary=True)
    transposome=graph.input("transposome","Mate Pair Tagment Enzyme loaded adapters",
                            MaterialKind.ADAPTER)
    tn5=graph.input("tn5","mate-pair Tn5 transposase",MaterialKind.CATALYST)
    tagged=graph.product("tagged","long Tn5-tagged genomic fragments",
                         MaterialKind.MOLECULAR,state=tuple(tagmentation_scene().rows()))
    tagment_waste=graph.product("tagment_waste","non-mate-pair tagmentation products",
                                MaterialKind.WASTE)
    graph.react("tagment","Tagment long genomic DNA",StepKind.TAGMENTATION,
                ((dna,InputRole.TEMPLATE),(transposome,InputRole.STRUCTURAL),
                 (tn5,InputRole.CATALYST)),(tagged,tagment_waste))
    graph.discard(tagment_waste,"products outside the long mate-pair fragment range are removed")
    fill_polymerase=graph.input("fill_polymerase","strand-displacement gap-fill polymerase",
                                MaterialKind.CATALYST)
    filled=graph.product("filled","gap-filled long tagged fragments",MaterialKind.MOLECULAR)
    graph.react("gap_fill","Repair the nine-base Tn5 gaps",StepKind.EXTENSION,
                ((tagged,InputRole.TEMPLATE),(fill_polymerase,InputRole.CATALYST)),
                (filled,))
    circle_ligase=graph.input("circle_ligase","blunt-end circularization ligase",
                              MaterialKind.CATALYST)
    circles=graph.product("circles","circularized long mate-pair molecules",
                          MaterialKind.MOLECULAR)
    residual_linear=graph.product("residual_linear","uncircularized long molecules",
                                  MaterialKind.MOLECULAR)
    graph.react("circularize","Close the distant tagged ends by blunt ligation",
                StepKind.CIRCULARIZATION,
                ((filled,InputRole.TEMPLATE),(circle_ligase,InputRole.CATALYST)),
                (circles,residual_linear),geometry=CircularizationGeometry.known(
                    "double-strand intramolecular ligation",
                    statement="The vendor workflow blunt-ligates the two ends of each long tagged fragment."))
    exonuclease=graph.input("exonuclease","linear-DNA exonuclease",MaterialKind.CATALYST)
    linear_waste=graph.product("linear_waste","digested uncircularized DNA",MaterialKind.WASTE)
    graph.react("remove_linear","Digest molecules that remain linear",StepKind.DIGESTION,
                ((residual_linear,InputRole.TEMPLATE),(exonuclease,InputRole.CATALYST)),
                (linear_waste,))
    graph.discard(linear_waste,"uncircularized molecules are destroyed before junction recovery")
    shear=graph.input("shear_condition","Covaris shearing condition",MaterialKind.CONDITION)
    junction_fragments=graph.product("junction_fragments","junction-spanning shear fragments",
                                     MaterialKind.MOLECULAR,state=recovery_scene())
    ordinary_fragments=graph.product("ordinary_fragments","shear fragments lacking the old junction",
                                     MaterialKind.WASTE)
    graph.react("shear","Shear circular molecules into short fragments",StepKind.FRAGMENTATION,
                ((circles,InputRole.TEMPLATE),(shear,InputRole.CONDITION)),
                (junction_fragments,ordinary_fragments))
    streptavidin=graph.input("streptavidin","streptavidin capture beads",
                            MaterialKind.SELECTION)
    captured=graph.product("captured","captured biotin-marked junction fragments",
                           MaterialKind.MOLECULAR,state=recovery_scene())
    uncaptured=graph.product("uncaptured","uncaptured non-junction shear fragments",
                             MaterialKind.WASTE)
    graph.react("capture","Capture the old circle-junction fragments",
                StepKind.AFFINITY_CAPTURE,
                ((junction_fragments,InputRole.TEMPLATE),(streptavidin,InputRole.SELECTION)),
                (captured,uncaptured))
    graph.discard(uncaptured,"fragments not retaining the biotin-marked junction are washed away")
    graph.discard(ordinary_fragments,"ordinary shearing fragments lack the selected mate-pair junction")
    endrepair=graph.input("endrepair","end-repair enzyme mix",MaterialKind.CATALYST)
    repaired=graph.product("repaired","end-repaired captured junction fragments",
                           MaterialKind.MOLECULAR)
    graph.react("end_repair","Repair captured fragment ends",StepKind.END_REPAIR,
                ((captured,InputRole.TEMPLATE),(endrepair,InputRole.CATALYST)),(repaired,))
    tailing=graph.input("tailing_polymerase","dA-tailing polymerase",MaterialKind.CATALYST)
    tailed=graph.product("tailed","3-prime-dA-tailed junction fragments",
                         MaterialKind.MOLECULAR)
    graph.react("da_tail","Add one 3-prime dA to each end",StepKind.TAILING,
                ((repaired,InputRole.TEMPLATE),(tailing,InputRole.CATALYST)),(tailed,))
    adapter=graph.input("truseq_adapter","TruSeq DNA LT index adapter",MaterialKind.ADAPTER)
    ligase=graph.input("adapter_ligase","adapter ligase",MaterialKind.CATALYST)
    adapted=graph.product("adapted","TruSeq-adapted mate-pair junction fragments",
                          MaterialKind.MOLECULAR,state=FINAL_LIBRARY)
    adapter_geometry=reported_adapter_end(
        "TruSeq DNA LT adapter",kind="cohesive",polarity="3-prime",terminal_sequence="T",
        statement="The captured fragment is explicitly dA-tailed before the complementary TruSeq adapter is ligated.")
    graph.react("ligate_adapter","Ligate the TruSeq index adapter",StepKind.LIGATION,
                ((tailed,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),
                 (ligase,InputRole.CATALYST)),(adapted,),geometry=adapter_geometry)
    forward=graph.input("pcr_forward","TruSeq PCR forward primer",MaterialKind.OLIGO)
    reverse=graph.input("pcr_reverse","TruSeq indexed PCR reverse primer",MaterialKind.OLIGO)
    pcr_polymerase=graph.input("pcr_polymerase","library PCR polymerase",
                               MaterialKind.CATALYST)
    library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=FINAL_LIBRARY)
    graph.react("library_pcr","Amplify the captured mate-pair library",StepKind.PCR,
                ((adapted,InputRole.TEMPLATE),(forward,InputRole.STRUCTURAL),
                 (reverse,InputRole.STRUCTURAL),(pcr_polymerase,InputRole.CATALYST)),
                (library,),geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY,reason="The final TruSeq sites are exact; the vendor adapter precursor is represented by its reported ligating end."))
    graph.sequence(library,SequencingPlan(library=FINAL_LIBRARY,
                                          primers=tuple(SEQ_PRIMERS)))
    graph.declare_reviewed_scope((
        "genomic_dna","transposome","tn5","tagment","fill_polymerase","gap_fill",
        "circle_ligase","circularize","exonuclease","remove_linear","shear_condition",
        "shear","streptavidin","capture","endrepair","end_repair",
        "tailing_polymerase","da_tail","truseq_adapter","adapter_ligase",
        "ligate_adapter","pcr_forward","pcr_reverse","pcr_polymerase",
        "library_pcr","library"),SOURCE)
    return graph.require_publishable()
