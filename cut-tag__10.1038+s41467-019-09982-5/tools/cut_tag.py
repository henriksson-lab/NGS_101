"""Original Kaya-Okur et al. CUT&Tag molecular model."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
import nextera as nx
from batch_ngs import duplex_fragment, nextera_library, seg, targeting_rows
from chemdraw import Construct, Row, Scene
from pcr import PCRGeometry
from protocol_graph import (Disposition, Input, InputRole, Material, MaterialKind,
                            ProtocolGraph, SequencingPlan, Step, StepKind, Terminal)

TITLE = "CUT&Tag — antibody-tethered Tn5 chromatin profiling"
NOTES = "01_cut-tag.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1038/s41467-019-09982-5">Kaya-Okur et al., <i>Nature Communications</i> (2019)</a>.'
SUMMARY = "An antibody recruits Protein A–Tn5 to a chromatin epitope; magnesium activates adapter-loaded Tn5 so cleavage and adapter installation occur at the targeted site in one reaction."
CAVEAT = "The original paper specifies Tn5MEDS-A/B transposomes and indexed PCR. The page uses the repository’s canonical Nextera mosaic-end and primer sequences; protein geometry is schematic rather than a claimed covalent nucleic-acid structure."

FINAL_LIBRARY, SEQ_PRIMERS = nextera_library(
    [seg("target-proximal chromatin", "X" * 36, placeholder=True)], "CUT&Tag library")
FINAL_CAPTION = "A productive heterologous S5/S7 tagmentation product after gap fill and indexed PCR. Tn5 installs the mosaic ends at the antibody-selected chromatin fragment boundaries."
SEQUENCING_INTRO = "Paired reads enter the targeted chromatin from the two Tn5 insertion sites; index reads identify the PCR-added sample indices."

PREAMP_TEMPLATE = Construct([
    seg("S5", nx.S5, "s5"), seg("left mosaic end", nx.ME, "me"),
    seg("target-proximal chromatin", "X" * 36, placeholder=True),
    seg("right mosaic end reverse complement", nx.ME_RC, "me"),
    seg("S7 reverse complement", nx.S7_RC, "s7"),
], name="gap-filled CUT&Tag template")

def sections():
    frag = duplex_fragment("epitope-proximal DNA")
    return [
        ("Bind antibody and adapter-loaded pA–Tn5",
         targeting_rows("Tn5 [S5-ME + S7-ME]", "targeted transposome retained after washes"),
         "Primary antibody supplies specificity; Protein A couples the antibody to a Tn5 dimer preloaded with two different mosaic-end adapters."),
        ("Activate targeted tagmentation",
         [Row(chunks=[("S5—ME ** target-proximal DNA ** ME—S7", "me", False)]),
          Row(chunks=[("          Mg²⁺-activated cleavage + adapter transfer", None, False)])],
         "** marks the two coupled cleavage/adapter-transfer junctions. Free pA–Tn5 is washed away before magnesium activation, limiting background."),
        ("Gap-fill and amplify only productive S5/S7 fragments", frag.rows(),
         "Tagmentation leaves the standard 9-bp staggered gap. Extension and indexed PCR complete the two-ended Illumina library."),
    ]


def protocol_graph() -> ProtocolGraph:
    """Source-scoped material graph for targeted transposition and library PCR."""
    source = "Kaya-Okur et al. 2019, doi:10.1038/s41467-019-09982-5"
    graph = ProtocolGraph(TITLE, authority=source)

    def material(key, label, kind, *, state=None, primary=False):
        value = Material(key, label, kind, state=state, primary=primary, source=source)
        graph.add_input(value)
        return value

    nuclei = material("nuclei", "permeabilized nuclei", MaterialKind.BIOLOGICAL,
                      state=targeting_rows("Tn5 [S5-ME + S7-ME]", "targeted transposome"),
                      primary=True)
    antibody = material("primary_antibody", "target-specific primary antibody",
                        MaterialKind.SELECTION)
    loaded = material("loaded_mosaic_end_adapters", "Tn5MEDS-A/B adapter duplexes",
                      MaterialKind.ADAPTER)
    tn5 = material("pa_tn5", "Protein A–Tn5 transposase", MaterialKind.CATALYST)
    targeted = Material("targeted_nuclei", "washed antibody-tethered pA–Tn5 nuclei",
                        MaterialKind.MOLECULAR, state=targeting_rows(
                            "Tn5 [S5-ME + S7-ME]", "targeted transposome retained"),
                        source=source)
    unbound = Material("unbound_targeting_reagents", "unbound targeting reagents",
                       MaterialKind.WASTE, source=source)
    graph.add_step(Step(
        "target_chromatin", "Bind antibody and adapter-loaded pA–Tn5",
        StepKind.AFFINITY_CAPTURE,
        (Input(nuclei.id, InputRole.TEMPLATE),
         Input(antibody.id, InputRole.SELECTION),
         Input(tn5.id, InputRole.CATALYST),
         Input(loaded.id, InputRole.STRUCTURAL)),
        (targeted.id, unbound.id), source=source), (targeted, unbound))
    graph.mark_terminal(Terminal(
        unbound.id, Disposition.DISCARD,
        "removed by the washes before magnesium activation"))

    magnesium = material("magnesium", "magnesium activation condition",
                         MaterialKind.CONDITION)
    tagged = Material(
        "tagmented_chromatin", "productively tagmented target-proximal chromatin",
        MaterialKind.MOLECULAR,
        state=sections()[1][1], source=source)
    non_library = Material(
        "non_library_tn5_products", "non-amplifiable Tn5 products",
        MaterialKind.WASTE, source=source)
    graph.add_step(Step(
        "targeted_tagmentation", "Activate targeted tagmentation",
        StepKind.TAGMENTATION,
        (Input(targeted.id, InputRole.TEMPLATE),
         Input(loaded.id, InputRole.STRUCTURAL),
         Input(tn5.id, InputRole.CATALYST),
         Input(magnesium.id, InputRole.CONDITION)),
        (tagged.id, non_library.id), source=source), (tagged, non_library))
    graph.mark_terminal(Terminal(
        non_library.id, Disposition.DISCARD,
        "symmetric or otherwise non-amplifiable transposition products"))

    polymerase = material("gap_fill_polymerase", "gap-fill polymerase",
                          MaterialKind.CATALYST)
    filled = Material("gap_filled_template", "gap-filled S5/S7 CUT&Tag template",
                      MaterialKind.MOLECULAR, state=PREAMP_TEMPLATE, source=source)
    graph.add_step(Step(
        "gap_fill", "Fill the staggered transposition gaps", StepKind.EXTENSION,
        (Input(tagged.id, InputRole.TEMPLATE),
         Input(polymerase.id, InputRole.CATALYST)),
        (filled.id,), source=source), (filled,))

    n5 = material("n5_index_primer", "indexed N5 PCR primer family", MaterialKind.OLIGO)
    n7 = material("n7_index_primer", "indexed N7 PCR primer family", MaterialKind.OLIGO)
    library = Material("cut_tag_library", "CUT&Tag sequencing library",
                       MaterialKind.LIBRARY, state=FINAL_LIBRARY, source=source)
    graph.add_step(Step(
        "indexed_pcr", "Amplify productive S5/S7 fragments", StepKind.PCR,
        (Input(filled.id, InputRole.TEMPLATE),
         Input(n5.id, InputRole.STRUCTURAL),
         Input(n7.id, InputRole.STRUCTURAL)),
        (library.id,), geometry=PCRGeometry.unavailable(
            PREAMP_TEMPLATE,
            reason="The defining paper names indexed PCR but does not print one immutable selected N5/N7 pair."),
        source=source), (library,))
    graph.mark_terminal(Terminal(library.id, Disposition.SEQUENCING),
                        SequencingPlan(FINAL_LIBRARY, tuple(SEQ_PRIMERS),
                                       tuple(p.role for p in SEQ_PRIMERS)))

    graph.declare_reviewed_scope((
        "nuclei", "primary_antibody", "loaded_mosaic_end_adapters", "pa_tn5",
        "target_chromatin", "magnesium", "targeted_tagmentation",
        "gap_fill_polymerase", "gap_fill", "n5_index_primer", "n7_index_primer",
        "indexed_pcr", "cut_tag_library",
    ), source)
    return graph.require_publishable()
