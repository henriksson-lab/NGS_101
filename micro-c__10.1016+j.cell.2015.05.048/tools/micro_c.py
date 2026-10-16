"""Molecular model for the original Hsieh et al. Micro-C protocol."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from chemdraw import Construct, Scene, Segment, complement_segments
from end_repair import fill_five_prime_overhang
from single_cell_hic import unresolved_illumina_library
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_graph import (Disposition, Terminal, InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)

TITLE = "Micro-C — nucleosome-resolution proximity ligation"
SOURCE = 'Defining source: <a href="https://doi.org/10.1016/j.cell.2015.05.048">Hsieh et al., <i>Cell</i> (2015)</a>.'


def seg(name, top, tag=None, **kw):
    return Segment(name=name, top=top, tag=tag, **kw)


# A representative MNase-generated end, not a recognition sequence.  Both the repair
# strand and positions at which biotin-dA/dC are incorporated are computed from it.
REPAIR = fill_five_prime_overhang("TAGC", substituted_bases=("A", "C"))
SIZE_SELECTION_BP = (250, 350)
PCR_CYCLES = (12, 15)
READ_LENGTHS = {"Read 1": 50, "Read 2": 50}


def mnase_scene():
    top = [seg("left linker", "X" * 10, placeholder=True),
           seg("nucleosomal DNA", "N" * 34, "insert", placeholder=True),
           seg("right linker", "X" * 10, placeholder=True)]
    sc = Scene.duplex(top, label="crosslinked chromatin")
    sc.junction("top", "left linker", "nucleosomal DNA", "MNase", ch="^")
    sc.junction("top", "nucleosomal DNA", "right linker", "MNase", ch="^")
    sc.mark("top", "nucleosomal DNA", "protected mononucleosome")
    return sc


def resected_end_scene():
    top = [seg(f"overhang {i + 1}", b, "me")
           for i, b in enumerate(REPAIR.overhang_5p)]
    body = seg("nucleosomal DNA", "X" * 28, placeholder=True)
    sc = Scene(); sc.strand("top", [*top, body], label="resected end")
    sc.anneal("bottom", complement_segments([body]), to="top",
              pair=("nucleosomal DNA'", "nucleosomal DNA"), label="resected end")
    sc.mark("top", "overhang 1", "5′ single-stranded end", through="overhang 4")
    return sc


def repaired_end_scene():
    top = [seg(f"overhang {i + 1}", b, "me")
           for i, b in enumerate(REPAIR.overhang_5p)]
    body = seg("nucleosomal DNA", "X" * 28, placeholder=True)
    fill = [seg(f"fill {i + 1}", b, "w1" if b in REPAIR.substituted_bases else "me")
            for i, b in enumerate(REPAIR.fill_5p)]
    bottom = [*complement_segments([body]), *fill]
    sc = Scene(); sc.strand("top", [*top, body], label="blunt repaired end")
    sc.anneal("bottom", bottom, to="top", pair=("fill 1", "overhang 4"),
              label="blunt repaired end")
    for i in REPAIR.substituted_positions:
        sc.mark("bottom", f"fill {i + 1}", f"new biotin-d{REPAIR.fill_5p[i]}")
    return sc


def contact_junction():
    """Representative blunt ligation; unlike restriction Hi-C it has no fixed motif."""
    return Construct([
        seg("nucleosome A", "X" * 28, placeholder=True),
        seg("repaired end A", REPAIR.overhang_5p, "w1"),
        seg("repaired end B", "ACGT", "w1"),
        seg("nucleosome B", "X" * 28, placeholder=True),
    ], name="Micro-C proximity-ligation product")


def selected_insert():
    return Construct([
        seg("nucleosome A read end", "X" * 32, placeholder=True),
        *list(contact_junction().slice("repaired end A", "repaired end B")),
        seg("nucleosome B read end", "X" * 32, placeholder=True),
    ], name="biotin-selected Micro-C insert")


def final_library():
    return unresolved_illumina_library(selected_insert(),
        "Micro-C sequencing library", indexed=None)


def protocol_graph():
    graph=ProtocolGraph(TITLE,authority=SOURCE)
    chromatin=graph.input("crosslinked_chromatin","crosslinked chromatin",
                          MaterialKind.BIOLOGICAL,state=tuple(mnase_scene().rows()),primary=True)
    mnase=graph.input("mnase","micrococcal nuclease",MaterialKind.CATALYST)
    mononucleosomes=graph.product("mononucleosomes","mononucleosome-protected chromatin fragments",
                                  MaterialKind.MOLECULAR,state=tuple(mnase_scene().rows()))
    linker_waste=graph.product("linker_waste","MNase linker-digestion products; recovery unassigned",MaterialKind.MOLECULAR)
    graph.react("mnase_digest","Digest linker DNA to mononucleosomes",StepKind.DIGESTION,
                ((chromatin,InputRole.TEMPLATE),(mnase,InputRole.CATALYST)),
                (mononucleosomes,linker_waste))
    graph.mark_terminal(Terminal(linker_waste.id,Disposition.UNRESOLVED,"MNase digestion is source-described; individual digestion products and their recovery through concentration/purification are unassigned."))
    phosphatase=graph.input("phosphatase","Antarctic phosphatase after concentration",MaterialKind.CATALYST)
    dephosphorylated=graph.product("dephosphorylated","dephosphorylated chromatin ends",MaterialKind.MOLECULAR,state=mononucleosomes.state)
    graph.react("dephosphorylate","Dephosphorylate before T4 polymerase end processing",StepKind.DEPHOSPHORYLATION,
                ((mononucleosomes,InputRole.TEMPLATE),(phosphatase,InputRole.CATALYST)),(dephosphorylated,))
    t4pol=graph.input("t4_polymerase","T4 DNA polymerase",MaterialKind.CATALYST)
    resected=graph.product("resected_ends","resected heterogeneous MNase ends",
                           MaterialKind.MOLECULAR,state=tuple(resected_end_scene().rows()))
    graph.react("resect_ends","Create repairable 5′ single-stranded ends",StepKind.END_REPAIR,
                ((dephosphorylated,InputRole.TEMPLATE),(t4pol,InputRole.CATALYST)),(resected,))
    labelled=graph.product("labelled_blunt_ends","biotin-filled blunt nucleosomal ends",
                           MaterialKind.MOLECULAR,state=tuple(repaired_end_scene().rows()))
    graph.react("fill_label","Fill ends with biotin-dATP and biotin-dCTP",StepKind.EXTENSION,
                ((resected,InputRole.TEMPLATE),(t4pol,InputRole.CATALYST)),(labelled,))
    ligase=graph.input("dna_ligase","T4 DNA ligase",MaterialKind.CATALYST)
    contacts=graph.product("contact_products","blunt proximity-ligation products",
                           MaterialKind.MOLECULAR,state=contact_junction())
    graph.react("proximity_ligate","Join nucleosomal ends held in proximity",
                StepKind.PROXIMITY_LIGATION,
                ((labelled,InputRole.TEMPLATE),(ligase,InputRole.CATALYST)),(contacts,))
    exo3=graph.input("exonuclease_iii","Exonuclease III",MaterialKind.CATALYST)
    protected=graph.product("protected_contacts","contacts with internal biotin junctions",
                            MaterialKind.MOLECULAR,state=selected_insert())
    free_ends=graph.product("free_end_waste","exonuclease products from terminally labelled ends",MaterialKind.MOLECULAR)
    graph.react("remove_free_ends","Remove unligated labelled ends",StepKind.DIGESTION,
                ((contacts,InputRole.TEMPLATE),(exo3,InputRole.CATALYST)),(protected,free_ends))
    graph.mark_terminal(Terminal(free_ends.id,Disposition.UNRESOLVED,"The source uses exonucleaseIII to eliminate biotinylated unligated termini, not a demonstrated physical discard of every entire unligated molecule; product recovery remains unassigned."))
    reversal=graph.input("crosslink_reversal","crosslink-reversal and proteinase condition",
                         MaterialKind.CONDITION)
    purified=graph.product("purified_contacts","deproteinized Micro-C contact DNA",
                           MaterialKind.MOLECULAR,state=selected_insert())
    graph.react("reverse_crosslinks","Reverse crosslinks and purify contact DNA",
                StepKind.DEPROTEINIZATION,
                ((protected,InputRole.TEMPLATE),(reversal,InputRole.CONDITION)),(purified,))
    gel=graph.input("gel_selection","250–350-bp gel selection",MaterialKind.SELECTION)
    selected=graph.product("selected_contacts","250–350-bp Micro-C contact inserts",
                           MaterialKind.MOLECULAR,state=selected_insert())
    offsize=graph.product("offsize_contacts","off-size contact molecules",MaterialKind.WASTE)
    graph.react("size_select","Select 250–350-bp contact products",StepKind.SELECTION,
                ((purified,InputRole.TEMPLATE),(gel,InputRole.SELECTION)),(selected,offsize))
    graph.discard(offsize,"molecules outside the reported gel interval are excluded")
    endprep=graph.input("library_endprep","library end-repair and dA-tailing enzymes",
                        MaterialKind.CATALYST)
    tailed=graph.product("tailed_contacts","end-repaired and dA-tailed contact inserts",
                         MaterialKind.MOLECULAR,state=selected_insert())
    graph.react("library_endprep_step","Repair and dA-tail the selected contacts",
                StepKind.END_REPAIR,
                ((selected,InputRole.TEMPLATE),(endprep,InputRole.CATALYST)),(tailed,))
    adapter=graph.input("illumina_adapter","unspecified Illumina adapter",MaterialKind.ADAPTER,
                        inferred=True,source=SOURCE+" (adapter kit not identified)")
    adapted=graph.product("adapted_contacts","Illumina-adapter-ligated contacts",
                          MaterialKind.MOLECULAR,inferred=True,
                          source=SOURCE+" (adapter architecture not identified)")
    graph.react("adapter_ligation","Ligate Illumina adapters",StepKind.LIGATION,
                ((tailed,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),
                 (ligase,InputRole.CATALYST)),(adapted,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The paper does not identify the adapter kit or ligation-ready end geometry."))
    streptavidin=graph.input("streptavidin","streptavidin beads",MaterialKind.SELECTION)
    captured=graph.product("captured_contacts","biotin-junction Micro-C molecules",
                            MaterialKind.MOLECULAR)
    uncaptured=graph.product("uncaptured_contacts","molecules lacking an internal biotin junction",
                             MaterialKind.WASTE)
    graph.react("capture","Capture internal-biotin contact molecules",StepKind.AFFINITY_CAPTURE,
                ((adapted,InputRole.TEMPLATE),(streptavidin,InputRole.SELECTION)),
                (captured,uncaptured))
    graph.discard(uncaptured,"molecules lacking a retained biotinylated contact junction are washed away")
    fwd=graph.input("pcr_forward","unspecified Illumina PCR forward primer",MaterialKind.OLIGO,
                    inferred=True,source=SOURCE+" (primer sequence not reported)")
    rev=graph.input("pcr_reverse","unspecified Illumina PCR reverse primer",MaterialKind.OLIGO,
                    inferred=True,source=SOURCE+" (primer sequence not reported)")
    pol=graph.input("pcr_polymerase","library PCR polymerase",MaterialKind.CATALYST)
    final=final_library()
    library=graph.product("library",final.name,MaterialKind.LIBRARY,state=final,inferred=True,
                          source=SOURCE+" (adapter architecture not identified)")
    graph.react("library_pcr","Amplify the captured library for 12–15 cycles",StepKind.PCR,
                ((captured,InputRole.TEMPLATE),(fwd,InputRole.STRUCTURAL),
                 (rev,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(library,),
                geometry=PCRGeometry.unavailable(
                    final,reason="The source does not report the adapter or primer sequences."))
    graph.sequence(library,SequencingPlan(
        acquired_roles=("Read 1","Read 2"),
        read_uses={"Read 1":"include","Read 2":"include"},
        read_lengths=READ_LENGTHS,
        limitation="Original Sequence analysis and mapping specifies paired50-base reads. Exact ordered primers and their physical sites remain unassembled; original supplements are not checked in this bounded pass."))
    graph.declare_reviewed_scope(("crosslinked_chromatin","mnase","mnase_digest",
      "phosphatase","dephosphorylate","t4_polymerase","resect_ends","fill_label","dna_ligase","proximity_ligate",
      "exonuclease_iii","remove_free_ends","crosslink_reversal","reverse_crosslinks",
      "gel_selection","size_select","library_endprep","library_endprep_step",
      "illumina_adapter","adapter_ligation","streptavidin","capture","pcr_forward",
      "pcr_reverse","pcr_polymerase","library_pcr","library"),SOURCE)
    return graph.require_publishable()
