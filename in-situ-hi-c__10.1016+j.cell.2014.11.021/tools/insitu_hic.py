"""Molecular model for the MboI implementation of Rao et al. in situ Hi-C."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from chemdraw import Construct, Scene, Segment, complement_segments
from restriction import Digest, MBOI
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, MolecularRelation, ProtocolGraph,
                            SequencingPlan, StepKind)
from protocol_model import molecule, role
from single_cell_hic import ContactWorkflow, unresolved_illumina_library


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


DEMO = Digest(MBOI, "X" * 18, "X" * 18)
FILLED = DEMO.fill_in(biotin_base="A")
SHEAR_RANGE_BP = (300, 500)
WORKFLOW = ContactWorkflow("in situ Hi-C", MBOI, "bulk nuclei", "biotin-fill",
                           "PCR", capture=True)
TITLE = "in situ Hi-C — MboI proximity ligation in intact nuclei"
CITATION = 'Defining source: <a href="https://doi.org/10.1016/j.cell.2014.11.021">Rao et al., <i>Cell</i> (2014)</a>.'
LIBRARY_CAVEAT = "Rao et al. name an indexed Illumina adapter but do not identify its kit generation or print its ligation geometry."
SEQUENCING_UNAVAILABLE = "The source does not print or identify the sequencing-primer set."
REQUIRED_ROLES = ("Read 1", "Read 2")


def cleavage_scene() -> Scene:
    top = [seg("left locus", DEMO.left_flank, placeholder=True),
           seg("MboI site", MBOI.site, "me"),
           seg("right locus", DEMO.right_flank, placeholder=True)]
    sc = Scene.duplex(top, label="chromatin DNA")
    sc.junction("top", "left locus", "MboI site", "top cut")
    sc.junction("bottom", "MboI site'", "right locus'", "bottom cut")
    sc.labels("top")
    return sc


def cut_left_scene() -> Scene:
    left = seg("left locus", DEMO.left_flank, placeholder=True)
    over = seg("5-prime GATC overhang", MBOI.overhang, "me")
    sc = Scene(); sc.strand("top", [left], label="left fragment")
    sc.anneal("bottom", [over, *complement_segments([left])], to="top",
              pair=("left locus'", "left locus"), label="left fragment")
    sc.mark("bottom", "5-prime GATC overhang", "cohesive end")
    return sc


def cut_right_scene() -> Scene:
    over = seg("5-prime GATC overhang", MBOI.overhang, "me")
    right = seg("right locus", DEMO.right_flank, placeholder=True)
    sc = Scene(); sc.strand("top", [over, right], label="right fragment")
    sc.anneal("bottom", complement_segments([right]), to="top",
              pair=("right locus'", "right locus"), label="right fragment")
    sc.mark("top", "5-prime GATC overhang", "cohesive end")
    return sc


def filled_left_scene() -> Scene:
    made = FILLED.left_new_top
    bio = FILLED.biotin_left_top[0]
    top = [seg("left locus", DEMO.left_flank, placeholder=True),
           seg("new fill before biotin", made[:bio], "me"),
           seg("new biotin-dA", made[bio:bio + 1], "w1"),
           seg("new fill after biotin", made[bio + 1:], "me")]
    sc = Scene.duplex(top, label="filled left end")
    sc.mark("top", "new fill before biotin", "Klenow fill-in",
            through="new fill after biotin")
    sc.mark("top", "new biotin-dA", "biotin")
    return sc


def filled_right_scene() -> Scene:
    # The new strand is stored 5'->3'; reverse its coordinates for the bottom
    # strand as drawn left-to-right (3'->5').
    overhang = MBOI.overhang
    drawn_bio = len(overhang) - 1 - FILLED.biotin_right_bottom[0]
    top = [seg("existing before opposite biotin", overhang[:drawn_bio], "me"),
           seg("opposite new biotin-dA", overhang[drawn_bio:drawn_bio + 1], "w1"),
           seg("existing after opposite biotin", overhang[drawn_bio + 1:], "me"),
           seg("right locus", DEMO.right_flank, placeholder=True)]
    sc = Scene.duplex(top, label="filled right end")
    sc.mark("bottom", "opposite new biotin-dA'", "new bottom-strand biotin-dA")
    return sc


def contact_junction() -> Construct:
    """One chimeric molecule after two filled MboI ends ligate."""
    junction = FILLED.junction
    top_bio = FILLED.junction_biotin_top[0]
    bottom_bio = FILLED.junction_biotin_bottom[0]
    boundary = MBOI.bottom_cut
    return Construct([
        seg("locus A", "X" * 22, placeholder=True),
        seg("left filled prefix", junction[:top_bio], "me"),
        seg("top biotin-dA", junction[top_bio:top_bio + 1], "w1"),
        seg("left filled suffix", junction[top_bio + 1:boundary], "me"),
        seg("right filled prefix", junction[boundary:bottom_bio], "me"),
        seg("opposite bottom biotin-dA", junction[bottom_bio:bottom_bio + 1], "w1"),
        seg("right filled suffix", junction[bottom_bio + 1:], "me"),
        seg("locus B", "X" * 22, placeholder=True),
    ], name="MboI proximity-ligation product")


contact_product = contact_junction


def sheared_insert() -> Construct:
    return Construct([
        seg("locus A read end", "X" * 28, placeholder=True),
        *list(contact_junction().slice("left filled prefix", "right filled suffix")),
        seg("locus B read end", "X" * 28, placeholder=True),
    ], name="biotin-selected Hi-C insert")


def final_library() -> Construct:
    """Indexed Illumina library; the protocol does not disclose adapter bases."""
    return unresolved_illumina_library(sheared_insert(),
        "in situ Hi-C sequencing library", indexed=True)


def protocol_graph() -> ProtocolGraph:
    """Rao's bulk-nucleus MboI route, including the source-limited library shell."""
    g = ProtocolGraph(TITLE, authority=CITATION)
    sample = g.input("sample", "fixed intact nuclei", MaterialKind.BIOLOGICAL,
                     state=molecule("fixed chromatin", (role("chromosomal DNA"),)),
                     primary=True)
    mboi = g.input("mboi", "MboI restriction enzyme", MaterialKind.CATALYST)
    cut = g.product("cut", "MboI-cut chromatin", MaterialKind.MOLECULAR,
                    state=cleavage_scene(), chemical_state=DEMO)
    g.react("digest", "Digest chromatin in intact nuclei with MboI", StepKind.DIGESTION,
            ((sample, InputRole.TEMPLATE), (mboi, InputRole.CATALYST)), (cut,))

    biotin = g.input("biotin_datp", "biotin-14-dATP", MaterialKind.MOLECULAR)
    klenow = g.input("klenow", "Klenow fill-in polymerase", MaterialKind.CATALYST)
    filled = g.product("filled", "biotin-filled MboI ends", MaterialKind.MOLECULAR,
                       state=filled_left_scene(), chemical_state=FILLED)
    g.react("fill", "Fill MboI ends with biotin-dATP", StepKind.EXTENSION,
            ((cut, InputRole.TEMPLATE), (biotin, InputRole.STRUCTURAL),
             (klenow, InputRole.CATALYST)), (filled,),
            relations=(MolecularRelation.RESTRICTION_FILL,))

    ligase = g.input("contact_ligase", "proximity-ligation DNA ligase",
                     MaterialKind.CATALYST)
    contacts = g.product("contacts", "biotin-marked contact junctions",
                         MaterialKind.MOLECULAR, state=contact_junction())
    g.react("contact_ligation", "Ligate filled chromatin ends in intact nuclei",
            StepKind.PROXIMITY_LIGATION,
            ((filled, InputRole.TEMPLATE), (ligase, InputRole.CATALYST)),
            (contacts,), relations=(MolecularRelation.CONTACT_JUNCTION,))

    shear = g.input("shear", "fragmentation to 300–500 bp", MaterialKind.CONDITION)
    fragments = g.product("fragments", "sheared contact fragments",
                          MaterialKind.MOLECULAR, state=sheared_insert())
    g.react("fragment", "Reverse crosslinks and shear contact DNA",
            StepKind.FRAGMENTATION,
            ((contacts, InputRole.TEMPLATE), (shear, InputRole.CONDITION)),
            (fragments,))
    beads = g.input("streptavidin", "streptavidin beads", MaterialKind.SELECTION)
    selected = g.product("selected", "biotin-enriched contact fragments",
                         MaterialKind.MOLECULAR, state=sheared_insert())
    unbound = g.product("unbound", "fragments without a captured biotin junction",
                        MaterialKind.WASTE)
    g.react("capture", "Capture biotin-bearing contact fragments",
            StepKind.AFFINITY_CAPTURE,
            ((fragments, InputRole.TEMPLATE), (beads, InputRole.SELECTION)),
            (selected, unbound))
    g.discard(unbound, "does not carry a captured biotin-marked contact junction")

    adapter = g.input("adapter", "source-unspecified indexed Illumina adapter",
                      MaterialKind.ADAPTER, inferred=True, source=LIBRARY_CAVEAT)
    adapter_ligase = g.input("adapter_ligase", "adapter ligase", MaterialKind.CATALYST)
    ligated = g.product("ligated", "adapter-ligated contact DNA",
                        MaterialKind.MOLECULAR, state=sheared_insert())
    g.react("adapter_ligation", "Repair, dA-tail and ligate indexed adapters",
            StepKind.LIGATION,
            ((selected, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
             (adapter_ligase, InputRole.CATALYST)), (ligated,),
            geometry=AdapterLigationGeometry.unavailable(LIBRARY_CAVEAT))
    forward = g.input("pcr_forward", "Illumina PCR forward primer", MaterialKind.OLIGO,
                      inferred=True, source=SEQUENCING_UNAVAILABLE)
    reverse = g.input("pcr_reverse", "Illumina PCR reverse primer", MaterialKind.OLIGO,
                      inferred=True, source=SEQUENCING_UNAVAILABLE)
    polymerase = g.input("pcr_polymerase", "library PCR polymerase", MaterialKind.CATALYST)
    final = final_library()
    library = g.product("library", "indexed in situ Hi-C library",
                        MaterialKind.LIBRARY, state=final, inferred=True,
                        source=LIBRARY_CAVEAT)
    g.react("library_pcr", "PCR-amplify the captured contact library", StepKind.PCR,
            ((ligated, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
             (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
            (library,), geometry=PCRGeometry.unavailable(
                final, reason="Rao et al. do not print the PCR primers or adapter bases."))
    g.sequence(library, SequencingPlan(acquired_roles=REQUIRED_ROLES,
                                       limitation=SEQUENCING_UNAVAILABLE))
    g.declare_reviewed_scope(tuple(g.steps) + tuple(g.materials), CITATION)
    return g.require_publishable()
