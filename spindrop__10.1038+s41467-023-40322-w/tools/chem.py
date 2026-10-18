"""Protocol-owned source, molecular workflow and reaction graph."""
from types import SimpleNamespace
from chemdraw import Construct, Scene, Segment, feature
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import molecule, role

DEFINITION = SimpleNamespace(
    title='spinDrop',
    doi='10.1038/s41467-023-40322-w',
    citation='De Jonghe et al., Nature Communications (2023)',
    summary='Fluorescence sorting enriches occupied inDrop-v3 droplets before a second device picoinjects reverse-transcription reagent.',
    endpoint='The inDrop-v3 barcode read identifies cell and molecule; the insert read reports the transcript 3′ end.',
    notes='01_spindrop.html',
)
CELL = role("inDrop-v3 cell-barcode architecture; not yet declared", tag="cbc",
            molecular_feature=feature("spindrop_cell", "cell_barcode", "unknown"))
UMI = role("inDrop-v3 UMI; length not yet declared", tag="umi",
           molecular_feature=feature("spindrop_umi", "umi", "unknown"))
SAMPLE = role("inDrop-v3 library index architecture; not yet declared", tag="cbc",
              molecular_feature=feature("spindrop_library_index", "sample_index", "unknown"))
PAYLOAD = (CELL, UMI, role("3-prime transcript insert"))

TITLE = DEFINITION.title
NOTES = DEFINITION.notes
SOURCE = f'Defining source: <a href="https://doi.org/{DEFINITION.doi}">{DEFINITION.citation}</a>.'
SUMMARY = DEFINITION.summary
CAVEAT = "The downstream inDrop library preparation remains summarized; IVT, adapter construction and recovery intermediates have not been mapped here. The RT primer is UV-released from the bead before picoinjection; it is not supplied by the RT mix."
FINAL_LIBRARY = Construct([role("inDrop-v3 library end; not yet declared"), *PAYLOAD,
                           SAMPLE, role("inDrop-v3 opposite library end; not yet declared")],
                          name="spinDrop inherited library endpoint roles")
SEQ_PRIMERS = ()
FINAL_CAPTION = DEFINITION.endpoint
SEQUENCING_INTRO = DEFINITION.endpoint
READ_LENGTHS = {"Read 1": 61, "Index 1 (i7)": 8, "Index 2 (i5)": 8, "Read 2": 14}
SEQUENCING_ENDING = DEFINITION.endpoint
SEQUENCING_UNAVAILABLE = "The source-stated 61/8/8/14 acquired cycles are retained; the implementation-specific inDrop-v3 primers and identifier geometry remain unassigned in this model."

def oligos():
    return ()


def _droplet_population_state(label, cell_role, primer_label):
    """Cellular material and bead-primer species occupy separate rows."""
    sc = Scene()
    sc.strand("cellular material", [role(cell_role)], label=label+": cellular material")
    sc.strand("barcode primer", [CELL, UMI], label=primer_label)
    return sc


def protocol_graph() -> ProtocolGraph:
    d = DEFINITION
    final = FINAL_LIBRARY
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    inputs = graph.input(
        "cells_and_microgels", "fluorescently labelled cells and inDrop-v3 microgels",
        MaterialKind.BIOLOGICAL,
        state=_droplet_population_state("cell and inDrop-v3 microgel",
            "cell or nucleus", "microgel-bound barcode/UMI primer"), primary=True)
    sorting = graph.input("fads_gate", "occupied/viable/target fluorescence gates",
                          MaterialKind.SELECTION)
    selected = graph.product("selected_droplets", "selected occupied droplets",
                             MaterialKind.MOLECULAR,
                             state=_droplet_population_state("selected droplet",
                                 "selected cell", "microgel-bound barcode/UMI primer"))
    rejected = graph.product("rejected_droplets", "empty, damaged or unwanted droplets",
                             MaterialKind.WASTE)
    graph.react("sort_droplets", "Co-encapsulate and fluorescence-sort droplets",
                StepKind.SELECTION,
                ((inputs, InputRole.TEMPLATE), (sorting, InputRole.SELECTION)),
                (selected, rejected))
    graph.discard(rejected, "fails the inline fluorescence selection gates")
    uv = graph.input("uv_release", "7-minute UV exposure after 25-minute room-temperature lysis", MaterialKind.CONDITION)
    released = graph.product("released_primers", "selected droplets with solubilized bead primers", MaterialKind.MOLECULAR,
        state=_droplet_population_state("lysed selected droplet", "cellular lysate", "UV-released barcode/UMI primer"))
    graph.react("release_primers", "Lyse for 25 min at 23°C, then UV-release barcode primers (7 min)", StepKind.CHEMICAL_CONVERSION,
                ((selected, InputRole.TEMPLATE), (uv, InputRole.CONDITION)), (released,))
    heat = graph.input("heat_treatment", "droplet heat-treatment condition",
                       MaterialKind.CONDITION)
    treated = graph.product("treated_droplets", "heat-treated selected droplets",
                            MaterialKind.MOLECULAR, state=released.state)
    graph.react("heat_treat", "Heat-treat selected droplets (70°C, 10 min), then chill", StepKind.THERMAL_TREATMENT,
                ((released, InputRole.TEMPLATE), (heat, InputRole.CONDITION)),
                (treated,), changes_molecule=False)
    reviewed = ["cells_and_microgels", "fads_gate", "sort_droplets",
                "heat_treatment", "heat_treat", "uv_release", "release_primers"]
    # The released bead primer is already present before RT mix injection.
    rt_mix = graph.input("rt_mix", "picoinjected 1.8X First Strand Buffer, MgCl2, DTT and RNaseOUT", MaterialKind.CONDITION)
    injected = graph.product("rt_injected", "selected droplets after RT-mix picoinjection", MaterialKind.MOLECULAR, state=treated.state)
    graph.react("picoinject", "Picoinject RT solution into lysate droplets at 1:1 flow ratio", StepKind.GENERIC,
                ((treated, InputRole.TEMPLATE), (rt_mix, InputRole.CONDITION)), (injected,), changes_molecule=False)
    rt_primer = graph.input("library_rt_primer", "UV-released inDrop-v3 bead primer already present in the droplets", MaterialKind.OLIGO)
    rt = graph.input("library_rt_enzyme", "SuperScript III, 24 U/µL in injected solution", MaterialKind.CATALYST)
    cdna = graph.product("library_cdna", "barcoded first-strand cDNA after picoinjection", MaterialKind.MOLECULAR,
        state=Construct([CELL, UMI, role("captured transcript cDNA")], name="spinDrop first-strand cDNA roles"))
    graph.react("library_rt", "Reverse-transcribe 2 h at 50 °C, then 70 °C 15 min", StepKind.REVERSE_TRANSCRIPTION,
                ((injected, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL), (rt, InputRole.CATALYST)), (cdna,))
    library = graph.product("library_library", "spinDrop inDrop-v3 sequencing library", MaterialKind.LIBRARY, state=final)
    graph.react("library_unmodelled", "Inherited inDrop library preparation after reverse transcription", StepKind.MAGIC,
                ((cdna, InputRole.TEMPLATE),), (library,),
                reason="Published PMC10409775 Methods, Picoinjection of RT mix (line 293), delegates downstream preparation to inDrop. "
                       "Recovery, second-strand synthesis, IVT, later RNA/cDNA transformations, adapter assembly and final PCR are not yet mapped for this v3 implementation. "
                       "This is a model/source-review boundary; it is not evidence that those methods were unpublished.")
    graph.sequence(library, SequencingPlan(library=final,
        acquired_roles=("Read 1", "Index 1 (i7)", "Index 2 (i5)", "Read 2"), read_lengths=READ_LENGTHS,
        limitation="Published NextSeq profile: R1 61, I1 8, I2 8, R2 14 cycles. The v3 primer and identifier architecture has not been declared here; acquired index cycles do not establish index-region lengths."))
    reviewed.extend(("rt_mix", "picoinject", "library_rt_primer", "library_rt_enzyme", "library_rt",
                     "library_unmodelled", "library_library"))
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()
