"""Protocol-owned source, molecular workflow and reaction graph."""
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import chain, molecule, role, role_library
from extended_protocols import (Definition, CELL, UMI, SAMPLE, TARGET, GUIDE, RT_WELL, LIG_WELL, PCR_WELL, _finish_pcr, _rt_then_pcr)

DEFINITION = Definition('DisCo', '10.1038/s41592-021-01391-1', 'Bues et al., Nature Methods (2022)', 'Machine vision and multilayer valves deterministically co-encapsulate one low-input cell and one mRNA-capture bead for Drop-seq-derived 3′ RNA sequencing.', 'The bead barcode identifies the captured cell and the UMI identifies each transcript molecule.', '01_disco.html')
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
    return chain(molecule("low-input cell suspension and mRNA-capture beads", (role("single cells"), CELL, UMI)), (
        ("Machine-vision position exactly one cell and one bead at separate valves",
         molecule("deterministically paired cell and bead", (role("one cell"), CELL, UMI)),
         "Active positioning replaces Poisson co-loading for samples below about 500 cells."),
        ("Co-eject the pair into one droplet and sort it to the sample outlet",
         molecule("occupied DisCo droplet", (role("lysed cell RNA"), CELL, UMI)),
         "The device rejects empty or mis-paired events before molecular capture."),
        ("Capture poly(A) RNA and prepare a Drop-seq-derived library",
         molecule("DisCo library", (CELL, UMI, role("3-prime transcript insert"), SAMPLE)),
         "The molecular construct is bead-based 3′ RNA-seq; the innovation is deterministic partitioning."),
    ))


def oligos():
    return ()


def protocol_graph() -> ProtocolGraph:
    d = DEFINITION
    final = FINAL_LIBRARY
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    inputs = graph.input(
        "cells_and_beads", "low-input cells and mRNA-capture beads",
        MaterialKind.BIOLOGICAL,
        state=molecule("low-input suspension and capture beads", (
            role("single cells"), CELL, UMI)), primary=True)
    positioning = graph.input("vision_positioning", "machine-vision valve positioning condition",
                              MaterialKind.CONDITION)
    paired = graph.product("paired_cell_bead", "deterministically paired cell and bead",
                           MaterialKind.MOLECULAR,
                           state=molecule("deterministically paired cell and bead", (
                               role("one cell"), CELL, UMI)))
    rejected = graph.product("positioning_rejects", "empty or mis-paired valve events",
                             MaterialKind.WASTE)
    graph.react("position_pair", "Position exactly one cell and one bead",
                StepKind.PARTITION,
                ((inputs, InputRole.TEMPLATE), (positioning, InputRole.CONDITION)),
                (paired, rejected))
    graph.discard(rejected, "fails deterministic one-cell/one-bead positioning")
    droplet_condition = graph.input("droplet_condition", "co-ejection and droplet-sorting condition",
                                    MaterialKind.CONDITION)
    occupied = graph.product("occupied_droplet", "occupied DisCo droplet",
                             MaterialKind.MOLECULAR,
                             state=molecule("occupied DisCo droplet", (
                                 role("lysed cell RNA"), CELL, UMI)))
    failed_droplets = graph.product("failed_droplets", "unsorted or failed droplets",
                                    MaterialKind.WASTE)
    graph.react("make_droplet", "Co-eject the pair and sort the occupied droplet",
                StepKind.PARTITION,
                ((paired, InputRole.TEMPLATE), (droplet_condition, InputRole.CONDITION)),
                (occupied, failed_droplets))
    graph.discard(failed_droplets, "does not reach the selected sample outlet")
    reviewed = ["cells_and_beads", "vision_positioning", "position_pair",
                "droplet_condition", "make_droplet"]
    _rt_then_pcr(graph, "library", occupied, final, d.endpoint, reviewed,
                 rt_primer_label="Drop-seq bead oligo(dT) carrying cell barcode and UMI",
                 label="DisCo")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()
