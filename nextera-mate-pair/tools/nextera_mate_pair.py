"""Illumina Nextera Mate Pair library construction.

Authority: Illumina Nextera Mate Pair Library Prep Reference Guide 15035209 v02 (January 2016)
and the Nextera Mate Pair data-processing technical note (biotinylated junction adapter,
junction sequences, incomplete-enrichment products). Gel-free path drawn; gel-plus size
selection is a recorded branch, not drawn."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from adapter_ligation import reported_adapter_end
from circular import CircularizationGeometry
from chemdraw import Row
from chemdraw import Scene
import illumina as il
import nextera as nx
import seqprimers as sp
from batch_ngs import seg
from chemdraw import circle_rows, feature
from endprep import dA_tailed_scene, repair_and_dA_tail
from chemdraw import Construct, Segment
from pcr import PCRGeometry
from protocol_graph import (Highlight, InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)

DUPLICATE_JUNCTION = nx.ME_RC + nx.ME


def junction_insert(left: str = "left distant end", right: str = "right distant end",
                    *, duplicate: bool = True) -> list:
    """A recovered shearing fragment spanning the original circle-closing bond."""
    junction = ([seg("junction ME reverse complement", nx.ME_RC, "me"),
                 seg("junction ME", nx.ME, "me")]
                if duplicate else [seg("single junction ME", nx.ME_RC, "me")])
    return [seg(left, "L" * 28, placeholder=True), *junction,
            seg(right, "R" * 28, placeholder=True)]


def recovered_truseq_library() -> tuple[Construct, tuple]:
    """Single-index TruSeq library around a selected mate-pair junction."""
    insert = junction_insert()
    lib = Construct([
        # Single-index adapter: P5 and the Read 1 arm are one oligo, so their shared
        # ACAC occurs once (canonical universal adapter, illumina.TRUSEQ_P5_FULL).
        seg("P5", il.TRUSEQ_P5_FULL[:-len(il.TRUSEQ_READ1)], "p5"), seg("Read 1 arm", il.TRUSEQ_READ1, "r1"),
        *insert, seg("dA junction", "A"),
        seg("Index 1 / Read 2 arm", il.INDEX1_PRIMER, "r2"),
        seg("six-base i7 index reverse complement", "I" * 6, "cbc", placeholder=True,
            feature=feature("sample_index_i7", "sample_index", "unknown")),
        seg("P7 reverse complement", il.P7_RC, "p7"),
    ], name="Nextera mate-pair junction library")
    primers = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["R2"])
    problems = sp.verify(lib, primers, required_roles=tuple(p.role for p in primers))
    if problems:
        raise ValueError("mate-pair library: " + "; ".join(problems))
    return lib, primers


def circularization_rows() -> list[Row]:
    """Long tagged molecule closed at the future mate-pair junction."""
    long = Construct([
        seg("left tag", nx.ME_RC, "me"),
        seg("long genomic fragment", "G" * 48, placeholder=True),
        seg("right tag", nx.ME, "me"),
    ], name="long tagmented genomic molecule")
    return circle_rows(long, "blunt intramolecular ligation")


def tagmentation_scene() -> Scene:
    """Long fragment with both transferred junctions anchored to named segments."""
    parts = [seg("left ME tag", nx.ME_RC, "me"),
             seg("long genomic fragment", "G" * 48, placeholder=True),
             seg("right ME tag", nx.ME, "me")]
    sc = Scene.duplex(parts, label="long tagmented molecule")
    sc.junction("top", "left ME tag", "long genomic fragment", "Tn5 transfer")
    sc.junction("top", "long genomic fragment", "right ME tag", "Tn5 transfer")
    sc.labels("top")
    return sc


def recovery_scene() -> Scene:
    """The short, junction-spanning product retained by streptavidin capture."""
    con = Construct(junction_insert(), name="junction-spanning shear product")
    sc = Scene.duplex(list(con), label="captured shear product")
    sc.junction("top", "junction ME reverse complement", "junction ME",
                "circle-closing junction")
    sc.mark("top", "junction ME reverse complement", "biotin / streptavidin capture")
    sc.labels("top")
    return sc


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
            Row(chunks=[("strand displacement repairs the short single-stranded Tn5 gaps", None, False)]),
        ], "Mate Pair Tagment Enzyme fragments long DNA and tags both ends with a biotinylated junction adapter; ** marks the transferred tag boundaries."),
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
    transposome=graph.input("transposome","Mate Pair Tagment Enzyme biotinylated junction adapter",
                            MaterialKind.ADAPTER)
    tn5=graph.input("tn5","mate-pair Tn5 transposase",MaterialKind.CATALYST)
    tagged=graph.product("tagged","long Tn5-tagged genomic fragments",
                         MaterialKind.MOLECULAR,state=tuple(tagmentation_scene().rows()))
    graph.react("tagment","Tagment long genomic DNA (55 °C, 30 min; Zymo column purification)",StepKind.TAGMENTATION,
                ((dna,InputRole.TEMPLATE),(transposome,InputRole.STRUCTURAL),
                 (tn5,InputRole.CATALYST)),(tagged,))
    fill_polymerase=graph.input("fill_polymerase","Strand Displacement Enzyme Mix with dNTPs (20 °C, 30 min)",
                                MaterialKind.CATALYST)
    filled=graph.product("filled","gap-filled long tagged fragments",MaterialKind.MOLECULAR)
    graph.react("gap_fill","Repair the short single-stranded Tn5 gaps",StepKind.EXTENSION,
                ((tagged,InputRole.TEMPLATE),(fill_polymerase,InputRole.CATALYST)),
                (filled,))
    long_beads=graph.input("long_beads","AMPure XP beads, about 0.4x (gel-free 40 µl to 100 µl diluted reaction)",
                           MaterialKind.SELECTION)
    long_fragments=graph.product("long_fragments","purified long gap-filled tagged fragments",
                                 MaterialKind.MOLECULAR)
    short_waste=graph.product("short_waste","short tagmentation fragments (< 1500 bp) and reaction mix",
                              MaterialKind.WASTE)
    graph.react("remove_short","Remove short fragments with AMPure XP",StepKind.SELECTION,
                ((filled,InputRole.TEMPLATE),(long_beads,InputRole.SELECTION)),
                (long_fragments,short_waste),changes_molecule=False)
    graph.discard(short_waste,"the bead supernatant is discarded (guide: removes fragments < 1500 bp)")
    circle_ligase=graph.input("circle_ligase","Circularization Ligase, up to 600 ng DNA in 300 µl, 30 °C overnight",
                              MaterialKind.CATALYST)
    circles=graph.product("circles","circularized long mate-pair molecules",
                          MaterialKind.MOLECULAR,state=tuple(circularization_rows()))
    residual_linear=graph.product("residual_linear","uncircularized long molecules",
                                  MaterialKind.MOLECULAR)
    graph.react("circularize","Close the distant tagged ends by blunt ligation",
                StepKind.CIRCULARIZATION,
                ((long_fragments,InputRole.TEMPLATE),(circle_ligase,InputRole.CATALYST)),
                (circles,residual_linear),geometry=CircularizationGeometry.known(
                    "double-strand intramolecular ligation",
                    statement="The vendor workflow blunt-ligates the two ends of each long tagged fragment."))
    exonuclease=graph.input("exonuclease","Exonuclease added to the circularization reaction (37 °C 30 min, 70 °C 30 min, then Stop Ligation Buffer)",MaterialKind.CATALYST)
    linear_waste=graph.product("linear_waste","digested uncircularized DNA",MaterialKind.WASTE)
    graph.react("remove_linear","Digest molecules that remain linear",StepKind.DIGESTION,
                ((residual_linear,InputRole.TEMPLATE),(exonuclease,InputRole.CATALYST)),
                (linear_waste,))
    graph.discard(linear_waste,"uncircularized molecules are destroyed before junction recovery")
    shear=graph.input("shear_condition","Covaris S2/S220, 40 s at 6 °C (~300–1000 bp, 3′/5′ overhangs)",MaterialKind.CONDITION)
    sheared=graph.product("sheared","sheared circle fragments, junction-spanning and non-junction",
                          MaterialKind.MOLECULAR,state=recovery_scene())
    graph.react("shear","Shear circular molecules into short fragments",StepKind.FRAGMENTATION,
                ((circles,InputRole.TEMPLATE),(shear,InputRole.CONDITION)),(sheared,))
    streptavidin=graph.input("streptavidin","Dynabeads M-280 Streptavidin (20 °C, 15 min; 4 Bead Wash Buffer and 2 RSB washes)",
                            MaterialKind.SELECTION)
    captured=graph.product("captured","captured biotin-marked junction fragments",
                           MaterialKind.MOLECULAR,state=recovery_scene())
    uncaptured=graph.product("uncaptured","non-junction shear fragments lacking the biotinylated adapter",
                             MaterialKind.WASTE)
    graph.react("capture","Capture the old circle-junction fragments",
                StepKind.AFFINITY_CAPTURE,
                ((sheared,InputRole.TEMPLATE),(streptavidin,InputRole.SELECTION)),
                (captured,uncaptured))
    graph.discard(uncaptured,"supernatant and washes are discarded; the technote notes that incomplete enrichment can carry some over")
    endrepair=graph.input("endrepair","End Repair Mix, on bead (30 °C, 30 min)",MaterialKind.CATALYST)
    repaired=graph.product("repaired","end-repaired captured junction fragments",
                           MaterialKind.MOLECULAR)
    graph.react("end_repair","Repair captured fragment ends",StepKind.END_REPAIR,
                ((captured,InputRole.TEMPLATE),(endrepair,InputRole.CATALYST)),(repaired,))
    tailing=graph.input("tailing_polymerase","A-Tailing Mix, on bead (37 °C, 30 min)",MaterialKind.CATALYST)
    tailed=graph.product("tailed","3-prime-dA-tailed junction fragments",
                         MaterialKind.MOLECULAR)
    graph.react("da_tail","Add one 3-prime dA to each end",StepKind.TAILING,
                ((repaired,InputRole.TEMPLATE),(tailing,InputRole.CATALYST)),(tailed,))
    adapter=graph.input("truseq_adapter","TruSeq DNA Adapter Index (one of 12)",MaterialKind.ADAPTER)
    ligase=graph.input("adapter_ligase","Ligation Mix, on bead (30 °C, 10 min; Ligation Stop Buffer, bead washes)",MaterialKind.CATALYST)
    adapted=graph.product("adapted","TruSeq-adapted mate-pair junction fragments",
                          MaterialKind.MOLECULAR,state=FINAL_LIBRARY)
    adapter_geometry=reported_adapter_end(
        "TruSeq DNA LT adapter",kind="cohesive",polarity="3-prime",terminal_sequence="T",
        statement="The captured fragment is explicitly dA-tailed before the complementary TruSeq adapter is ligated.")
    graph.react("ligate_adapter","Ligate the TruSeq index adapter",StepKind.LIGATION,
                ((tailed,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),
                 (ligase,InputRole.CATALYST)),(adapted,),geometry=adapter_geometry)
    forward=graph.input("pcr_forward","PCR Primer Cocktail (P5-end primer; sequence not printed in the guide)",MaterialKind.OLIGO)
    reverse=graph.input("pcr_reverse","PCR Primer Cocktail (P7-end primer; sequence not printed in the guide)",MaterialKind.OLIGO)
    pcr_polymerase=graph.input("pcr_polymerase","Enhanced PCR Mix, on bead; 10 or 15 cycles (98/60/72 °C)",
                               MaterialKind.CATALYST)
    amplified=graph.product("amplified","amplified mate-pair library in PCR supernatant",MaterialKind.LIBRARY,state=FINAL_LIBRARY)
    graph.react("library_pcr","Amplify the captured mate-pair library",StepKind.PCR,
                ((adapted,InputRole.TEMPLATE),(forward,InputRole.STRUCTURAL),
                 (reverse,InputRole.STRUCTURAL),(pcr_polymerase,InputRole.CATALYST)),
                (amplified,),geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY,reason="The final TruSeq sites are exact; the vendor adapter precursor is represented by its reported ligating end."))
    cleanup_beads=graph.input("cleanup_beads","AMPure XP beads on the 45 µl PCR supernatant",MaterialKind.SELECTION)
    library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=FINAL_LIBRARY)
    short_library=graph.product("short_library","short library fragments (< 300 bp) and PCR reagents",MaterialKind.WASTE)
    graph.react("clean_up","Clean up libraries with AMPure XP",StepKind.CLEANUP,
                ((amplified,InputRole.TEMPLATE),(cleanup_beads,InputRole.SELECTION)),
                (library,short_library),changes_molecule=False)
    graph.discard(short_library,"bead supernatant is discarded (guide: removes library fragments < 300 bp)")
    graph.annotate(Highlight("capture","Incomplete enrichment and single-junction products",
        "The data-processing technical note states that non-junction fragments lack the biotinylated adapter and should be removed at enrichment, but incomplete enrichment leaves them as forward-reverse pairs with a < 1,000 bp insert peak. Shearing inside the junction adapter yields a fragment carrying only one original end; single-junction products (one adapter copy) also occur. The drawn duplicated junction is the common product."))
    graph.annotate(Highlight("remove_short","Gel-plus branch not drawn",
        "The gel-free path (1 µg input, 2–15 kb, median 2.5–4 kb) is drawn. The gel-plus path uses 4 µg input, a 100 µl bead addition to the 200 µl reaction, then Pippin Prep or agarose size selection with Zymo recovery before circularization; PCR is 10 or 15 cycles depending on recovered amount and size."))
    graph.sequence(library,SequencingPlan(library=FINAL_LIBRARY,
                                          primers=tuple(SEQ_PRIMERS)))
    graph.declare_reviewed_scope((
        "genomic_dna","transposome","tn5","tagment","fill_polymerase","gap_fill",
        "long_beads","remove_short","cleanup_beads","clean_up",
        "circle_ligase","circularize","exonuclease","remove_linear","shear_condition",
        "shear","streptavidin","capture","endrepair","end_repair",
        "tailing_polymerase","da_tail","truseq_adapter","adapter_ligase",
        "ligate_adapter","pcr_forward","pcr_reverse","pcr_polymerase",
        "library_pcr","library"),SOURCE)
    return graph.require_publishable()
