"""Construct model for the individual dscATAC-seq and dsciATAC-seq pages."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
import illumina as il
import nextera as nx
from chemdraw import Construct, Scene, Segment, feature
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)

C1="TATGCATGAC"; C2="AGTCACTGAG"; BC_NT=7; TN5_BC_NT=6
REP_TN5_BC="AAAGAA"

def _s(n,t,g=None,**kw): return Segment(n,t,g,**kw)

def bead_published_segments():
    return [_s("BC1","N"*BC_NT,"cbc",placeholder=True,feature=feature("cell_bc1","cell_barcode","whitelist",whitelist="published bead barcode set",group="cell_id",part="BC1")),
            _s("phase block","NN",placeholder=True),_s("constant 1",C1),
            _s("BC2","N"*BC_NT,"cbc",placeholder=True,feature=feature("cell_bc2","cell_barcode","whitelist",whitelist="published bead barcode set",group="cell_id",part="BC2")),_s("constant 2",C2),
            _s("BC3","N"*BC_NT,"cbc",placeholder=True,feature=feature("cell_bc3","cell_barcode","whitelist",whitelist="published bead barcode set",group="cell_id",part="BC3")),_s("s5 primer",nx.S5,"s5")]

def tn5_r1(protocol):
    if protocol=="dscATAC-seq": return [_s("s5",nx.S5,"s5"),_s("ME",nx.ME,"me")]
    if protocol=="dsciATAC-seq": return [_s("s5",nx.S5,"s5"),
        _s("Tn5 barcode",REP_TN5_BC,"cbc",feature=feature("cell_tn5","cell_barcode","whitelist",whitelist="published Tn5 barcode set",group="cell_id",part="Tn5")),_s("ME",nx.ME,"me")]
    raise ValueError(protocol)

def bead_priming_scene():
    template=[_s("s5'",nx.S5_RC,"s5"),_s("ME'",nx.ME_RC,"me"),
              _s("insert","XXXXXXXX...",placeholder=True)]
    bead=[_s("unpublished bead 5' arm","XXXXXXXXXXXX",inferred=True,placeholder=True),
          *bead_published_segments()]
    sc=Scene(); sc.strand("tagmented strand",template,label="fragment")
    sc.anneal("released bead oligo",bead,to="tagmented strand",
              pair=("s5 primer","s5'"),label="bead primer",
              unpaired=("unpublished bead 5' arm","BC1","phase block","constant 1",
                        "BC2","constant 2","BC3"))
    sc.arrow("released bead oligo","polymerase")
    return sc

def filled_fragment(protocol):
    return Construct([*tn5_r1(protocol),_s("genomic insert","XXXXXXXX...XXXXXXXX",placeholder=True),
                      _s("ME'",nx.ME_RC,"me"),_s("s7'",nx.S7_RC,"s7")],
                     name=f"{protocol} gap-filled fragment")

def final_library(protocol):
    inner=[*bead_published_segments()]
    if protocol=="dsciATAC-seq": inner.append(_s("Tn5 barcode", "N"*TN5_BC_NT,"cbc",placeholder=True,feature=feature("cell_tn5","cell_barcode","whitelist",whitelist="published Tn5 barcode set",group="cell_id",part="Tn5")))
    elif protocol!="dscATAC-seq": raise ValueError(protocol)
    inner += [_s("ME",nx.ME,"me"),_s("genomic insert","XXXXXXXX...XXXXXXXX",placeholder=True),
              _s("ME'",nx.ME_RC,"me"),_s("s7'",nx.S7_RC,"s7"),
              _s("i7'","N"*8,"cbc",placeholder=True,feature=feature("sample_i7","sample_index","unknown")),_s("P7'",il.P7_RC,"p7")]
    return Construct([_s("unpublished bead/P5 arm","[BEAD RELEASE / P5 / R1 ARM]",inferred=True,
                         placeholder=True),*inner],name=f"{protocol} final library")

def _validate():
    if "".join(s.top for s in tn5_r1("dsciATAC-seq")) != nx.S5+REP_TN5_BC+nx.ME:
        raise ValueError("barcoded Tn5 adaptor no longer matches Supplementary Table 6")
    bead_priming_scene().rows()
_validate()

def protocol_graph(protocol):
    """The defining dscATAC-seq graph; dsciATAC owns its assembly locally."""
    if protocol != "dscATAC-seq":
        raise ValueError("dsciATAC-seq assembles its graph in its own directory")
    source = "Defining source: https://doi.org/10.1038/s41587-019-0147-6"
    graph = ProtocolGraph(protocol, authority=source)
    nuclei = graph.input("nuclei", "intact cells or nuclei", MaterialKind.BIOLOGICAL,
                         primary=True)
    tagged = tagment_cells(graph, nuclei, protocol=protocol)
    filled = fill_transposition_gaps(graph, tagged, protocol=protocol)
    droplets = partition_with_beads(graph, filled)
    copy_bead_barcode(graph, droplets, protocol=protocol, source=source)
    graph.declare_reviewed_scope((
        "nuclei", "transposome", "tn5", "tagment", "gap_fill_polymerase", "gap_fill",
        "droplet_loading", "partition", "bead_oligo", "reverse_primer", "pcr_polymerase",
        "droplet_pcr", "library"), source)
    return graph.require_publishable()


def tagment_cells(graph: ProtocolGraph, nuclei, *, protocol: str):
    transposome = graph.input("transposome",
        "barcoded Tn5 transposome" if protocol == "dsciATAC-seq" else "Tn5 transposome",
        MaterialKind.ADAPTER, state=tuple(tn5_r1(protocol)))
    tn5 = graph.input("tn5", "Tn5 transposase", MaterialKind.CATALYST)
    tagged = graph.product("tagmented_cells", "intact cells containing tagmented chromatin",
                           MaterialKind.MOLECULAR)
    waste = graph.product("tagmentation_waste", "non-library tagmentation products",
                          MaterialKind.WASTE)
    graph.react("tagment", "Tagment accessible chromatin while retaining intact cells",
                StepKind.TAGMENTATION,
                ((nuclei, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged, waste))
    graph.discard(waste, "fragments without usable library ends are not amplified")
    return tagged


def fill_transposition_gaps(graph: ProtocolGraph, tagged, *, protocol: str):
    polymerase = graph.input("gap_fill_polymerase", "gap-fill polymerase", MaterialKind.CATALYST)
    filled = graph.product("gap_filled_cells", "cells containing gap-filled tagmented fragments",
                           MaterialKind.MOLECULAR, state=filled_fragment(protocol))
    graph.react("gap_fill", "Fill transposition gaps", StepKind.EXTENSION,
                ((tagged, InputRole.TEMPLATE), (polymerase, InputRole.CATALYST)), (filled,))
    return filled


def partition_with_beads(graph: ProtocolGraph, filled):
    loading = graph.input("droplet_loading", "cell and bead co-encapsulation condition",
                          MaterialKind.CONDITION)
    droplets = graph.product("productive_droplets", "droplets containing one cell and one barcode bead",
                             MaterialKind.MOLECULAR)
    excluded = graph.product("excluded_droplets", "empty or multiplet droplets", MaterialKind.WASTE)
    graph.react("partition", "Co-encapsulate tagmented cells with super-loaded barcode beads",
                StepKind.PARTITION,
                ((filled, InputRole.TEMPLATE), (loading, InputRole.CONDITION)),
                (droplets, excluded))
    graph.discard(excluded, "empty and multiplet droplets are excluded from single-cell libraries")
    return droplets


def copy_bead_barcode(graph: ProtocolGraph, droplets, *, protocol: str, source: str):
    bead = graph.input("bead_oligo", "published bead-barcode/s5 oligo portion", MaterialKind.OLIGO,
                       state=tuple(bead_published_segments()), inferred=True,
                       source=source + " (5-prime bead/P5 arm is unpublished)")
    reverse = graph.input("reverse_primer", "droplet PCR reverse primer", MaterialKind.OLIGO,
                          inferred=True, source=source + " (sequence unpublished)")
    polymerase = graph.input("pcr_polymerase", "droplet PCR polymerase", MaterialKind.CATALYST)
    final = final_library(protocol)
    library = graph.product("library", final.name, MaterialKind.LIBRARY, inferred=True,
                            source=source + " (bead/P5 arm unpublished)")
    graph.react("droplet_pcr", "Copy the bead barcode and complete the library in droplets",
                StepKind.PCR,
                ((droplets, InputRole.TEMPLATE), (bead, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (library,), geometry=PCRGeometry.unavailable(
                    final, reason="The custom Read-1 and library-primer sequences are unpublished."))
    graph.sequence(library, SequencingPlan(
        limitation="The defining paper reports read roles and cycles but does not publish the custom Read-1 and complete bead/P5 primer sequences.",
        acquired_roles=("Read 1", "Index 1 (i7)", "Read 2"),
        read_lengths={"Read 1": 118, "Index 1 (i7)": 8, "Read 2": 40},
        read_uses={"Read 1": "extract_then_include", "Index 1 (i7)": "demultiplex",
                   "Read 2": "include"}))
    return library
