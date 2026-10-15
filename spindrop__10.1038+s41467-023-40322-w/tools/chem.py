"""Protocol-owned source, molecular workflow and reaction graph."""
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import chain, molecule, role, role_library
from extended_protocols import (Definition, CELL, UMI, SAMPLE, TARGET, GUIDE, RT_WELL, LIG_WELL, PCR_WELL, _finish_pcr, _rt_then_pcr)

DEFINITION = Definition('spinDrop', '10.1038/s41467-023-40322-w', 'Hochgerner et al., Nature Communications (2023)', 'Fluorescence sorting enriches occupied inDrop-v3 droplets before a second device picoinjects reverse-transcription reagent.', 'The inDrop-v3 barcode read identifies cell and molecule; the insert read reports the transcript 3′ end.', '01_spindrop.html')
PAYLOAD = (CELL, UMI, role("3-prime transcript insert"))

TITLE = DEFINITION.title
NOTES = DEFINITION.notes
SOURCE = f'Defining source: <a href="https://doi.org/{DEFINITION.doi}">{DEFINITION.citation}</a>.'
SUMMARY = DEFINITION.summary
CAVEAT = ""
FINAL_LIBRARY = role_library(TITLE + " final library roles", PAYLOAD, index=SAMPLE)
SEQ_PRIMERS = ()
FINAL_CAPTION = DEFINITION.endpoint
SEQUENCING_INTRO = DEFINITION.endpoint
READ_LENGTHS = None
SEQUENCING_ENDING = DEFINITION.endpoint
SEQUENCING_UNAVAILABLE = "The source defines the read roles, but the primer sequence is panel-, kit- or implementation-dependent; no universal binding oligo is asserted."

def workflow():
    return chain(molecule("fluorescently labelled cell and inDrop-v3 microgel", (role("cell or nucleus"), CELL, UMI)), (
        ("Co-encapsulate and fluorescence-sort occupied, viable or target droplets",
         molecule("selected droplet", (role("selected cell"), CELL, UMI)),
         "Inline FADS removes empty, damaged or unwanted droplets."),
        ("Heat-treat and picoinject an improved reverse-transcription mixture",
         molecule("barcoded first-strand cDNA", (CELL, UMI, role("3-prime transcript cDNA"))),
         "A second device adds RT reagent after sorting, enabling a multistep droplet reaction."),
        ("Recover, amplify and index cDNA",
         molecule("spinDrop library", (CELL, UMI, role("3-prime transcript insert"), SAMPLE)),
         "The final read structure follows the inDrop-v3 barcode scheme."),
    ))


def oligos():
    return ()


def protocol_graph() -> ProtocolGraph:
    d = DEFINITION
    final = FINAL_LIBRARY
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    inputs = graph.input(
        "cells_and_microgels", "fluorescently labelled cells and inDrop-v3 microgels",
        MaterialKind.BIOLOGICAL,
        state=molecule("cell and inDrop-v3 microgel", (
            role("cell or nucleus"), CELL, UMI)), primary=True)
    sorting = graph.input("fads_gate", "occupied/viable/target fluorescence gates",
                          MaterialKind.SELECTION)
    selected = graph.product("selected_droplets", "selected occupied droplets",
                             MaterialKind.MOLECULAR,
                             state=molecule("selected droplet", (
                                 role("selected cell"), CELL, UMI)))
    rejected = graph.product("rejected_droplets", "empty, damaged or unwanted droplets",
                             MaterialKind.WASTE)
    graph.react("sort_droplets", "Co-encapsulate and fluorescence-sort droplets",
                StepKind.SELECTION,
                ((inputs, InputRole.TEMPLATE), (sorting, InputRole.SELECTION)),
                (selected, rejected))
    graph.discard(rejected, "fails the inline fluorescence selection gates")
    heat = graph.input("heat_treatment", "droplet heat-treatment condition",
                       MaterialKind.CONDITION)
    treated = graph.product("treated_droplets", "heat-treated selected droplets",
                            MaterialKind.MOLECULAR, state=selected.state)
    graph.react("heat_treat", "Heat-treat selected droplets", StepKind.THERMAL_TREATMENT,
                ((selected, InputRole.TEMPLATE), (heat, InputRole.CONDITION)),
                (treated,), changes_molecule=False)
    reviewed = ["cells_and_microgels", "fads_gate", "sort_droplets",
                "heat_treatment", "heat_treat"]
    _rt_then_pcr(graph, "library", treated, final, d.endpoint, reviewed,
                 rt_primer_label="picoinjected improved inDrop-v3 RT/bead-primer mixture",
                 label="spinDrop")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()
