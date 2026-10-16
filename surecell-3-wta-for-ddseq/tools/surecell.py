"""Black-box molecular model for the commercial SureCell WTA 3' kit.

Illumina's Reference Guide #1000000021452 v01 is authoritative for the reaction flow,
but prints no oligo sequences.  Consequently every reagent-derived molecular region is
created by the inferred-placeholder helper below; no secondary reconstruction enters
the model as established sequence.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from chemdraw import Construct, Scene, Segment, feature
from pcr import PCRGeometry
from protocol_graph import (InputRole,MaterialKind,ProtocolGraph,SequencingPlan,StepKind)


def inferred(name: str, token: str, tag: str | None = None, **kw) -> Segment:
    """Create a visibly uncertain, non-sequence token for an unverified kit region."""
    if not token.startswith("[") or not token.endswith("]"):
        raise ValueError("unverified region tokens must be bracketed descriptions")
    return Segment(name=name, top=token, tag=tag, placeholder=True, inferred=True,
                   note="sequence and exact length not specified in the checked Illumina kit guide", **kw)


def inferred_poly_t(name: str = "oligo-dT", display_bases: int = 12) -> Segment:
    """Draw the documented poly(dT) role without asserting its unverified length."""
    if display_bases < 2:
        raise ValueError("poly(dT) display must contain enough bases to be legible")
    return Segment(name, "T" * display_bases, inferred=True,
                   note="poly(dT) role is documented; exact length is not specified in the checked guide")


def bead_primer() -> Construct:
    return Construct([
        inferred("3' Barcode Mix oligo", "[unverified barcode region]", "cbc",
                 feature=feature("cell_barcode_region", "cell_barcode", "unknown",
                                 note="internal barcode anatomy and length are unresolved")),
        inferred_poly_t(),
    ], name="SureCell bead oligo (schematic, not to scale)")


def first_strand() -> Construct:
    return Construct([
        *bead_primer().segments,
        Segment("antisense cDNA", "XXXXXXXX...XXXXXXXX", placeholder=True),
    ], name="SureCell first-strand cDNA")


def selected_fragment() -> Construct:
    """Bead-end fragment selected by TPP1 and one N7xx adapter during PCR."""
    return Construct([
        inferred("TPP1-selected end", "[TPP1-selected barcode end]", "cbc",
                 feature=feature("cell_barcode_region", "cell_barcode", "unknown",
                                 note="internal barcode anatomy and length are unresolved")),
        inferred_poly_t("poly(T)"),
        Segment("3' cDNA fragment", "XXXXXXXX...XXXXXXXX", placeholder=True),
        inferred("SureCell transposome end", "[transposome-derived end]", "me"),
        inferred("N7xx adapter", "[N7xx adapter, including index]", "p7",
                 feature=feature("sample_i7", "sample_index", "whitelist",
                                 whitelist="Illumina N7xx index set")),
    ], name="SureCell final library (schematic, not to scale)")


def schematic_duplex(con: Construct) -> Scene:
    """Pair known bases and explicitly exempt only non-sequence role placeholders."""
    unpaired = [s.name + "'" for s in con if s.placeholder]
    return Scene.duplex(list(con), unpaired=unpaired)


def read_boundaries() -> tuple[tuple[str, str], ...]:
    """Return only sequencing roles stated by the vendor guide."""
    con = selected_fragment()
    names = [s.name for s in con]
    required = ["TPP1-selected end", "3' cDNA fragment", "N7xx adapter"]
    if any(x not in names for x in required):
        raise ValueError("SureCell construct lacks a required read-layout boundary")
    return (
        ("Read 1", "supplied custom Sequencing Primer SP; molecular start undisclosed"),
        ("Index 1", "the selected N7xx adapter index"),
    )


def protocol_graph() -> ProtocolGraph:
    source=("Illumina SureCell WTA 3' Library Prep Reference Guide, "
            "document 1000000021452 v01 (June 2017)")
    graph=ProtocolGraph("SureCell WTA 3' for ddSEQ",authority=source)
    cells=graph.input("cells","single-cell suspension",MaterialKind.BIOLOGICAL,primary=True)
    beads=graph.input("barcode_beads","3-prime Barcode Mix beads",MaterialKind.ADAPTER,
                      state=bead_primer(),inferred=True,
                      source=source+" (oligo sequences unresolved)")
    droplets=graph.input("droplet_device","ddSEQ droplet generator",MaterialKind.CONDITION)
    captured=graph.product("captured_rna","droplet-confined bead-captured poly(A) RNA",
                           MaterialKind.MOLECULAR,inferred=True,
                           source=source+" (barcode anatomy unresolved)")
    empty=graph.product("empty_droplets","unproductive cell/bead loading outcomes",
                        MaterialKind.WASTE)
    graph.react("partition","Co-encapsulate cells, RT mix and barcode beads",
                StepKind.PARTITION,
                ((cells,InputRole.TEMPLATE),(beads,InputRole.STRUCTURAL),
                 (droplets,InputRole.CONDITION)),(captured,empty))
    graph.discard(empty,"droplets without a productive cell-bead pairing do not yield library")
    rt=graph.input("reverse_transcriptase","reverse transcriptase",MaterialKind.CATALYST)
    first=graph.product("first_strand",first_strand().name,MaterialKind.MOLECULAR,
                        state=first_strand(),inferred=True,
                        source=source+" (barcode anatomy unresolved)")
    graph.react("reverse_transcribe","Reverse-transcribe bead-captured RNA in droplets",
                StepKind.REVERSE_TRANSCRIPTION,
                ((captured,InputRole.TEMPLATE),(beads,InputRole.STRUCTURAL),
                 (rt,InputRole.CATALYST)),(first,))
    purification=graph.input("purification_beads","cDNA purification beads",
                             MaterialKind.SELECTION)
    cleaned=graph.product("cleaned_first_strand","purified bead-barcoded first-strand cDNA",
                          MaterialKind.MOLECULAR,state=first_strand(),inferred=True,
                          source=source+" (barcode anatomy unresolved)")
    free_reagent=graph.product("free_barcode_reagent","unbound barcode reagent",
                               MaterialKind.WASTE)
    graph.react("cleanup","Break droplets and remove unbound barcode reagent",
                StepKind.CLEANUP,
                ((first,InputRole.TEMPLATE),(purification,InputRole.SELECTION)),
                (cleaned,free_reagent))
    graph.discard(free_reagent,"free barcode reagent is removed before second-strand synthesis")
    rna_removal=graph.input("rna_removal","RNA-removal condition",MaterialKind.CONDITION)
    stripped=graph.product("stripped_first_strand","RNA-free first-strand cDNA",
                           MaterialKind.MOLECULAR,state=first_strand(),inferred=True,
                           source=source+" (barcode anatomy unresolved)")
    rna_waste=graph.product("rna_waste","removed RNA template",MaterialKind.WASTE)
    graph.react("remove_rna","Remove the RNA template",StepKind.STRAND_SEPARATION,
                ((cleaned,InputRole.TEMPLATE),(rna_removal,InputRole.CONDITION)),
                (stripped,rna_waste))
    graph.discard(rna_waste,"RNA is removed before replacement-strand synthesis")
    second_primer=graph.input("second_strand_primer","second-strand reagent (sequence unverified)",
                              MaterialKind.OLIGO,inferred=True,
                              source=source+" (sequence unresolved)")
    second_pol=graph.input("second_strand_polymerase","second-strand polymerase",
                           MaterialKind.CATALYST)
    duplex=graph.product("double_stranded_cdna","replacement-strand SureCell cDNA",
                         MaterialKind.MOLECULAR,state=first_strand(),inferred=True,
                         source=source+" (second-strand sequence unresolved)")
    graph.react("second_strand","Synthesize the replacement cDNA strand",StepKind.EXTENSION,
                ((stripped,InputRole.TEMPLATE),(second_primer,InputRole.STRUCTURAL),
                 (second_pol,InputRole.CATALYST)),(duplex,))
    transposome=graph.input("transposome","Nextera SureCell transposome",MaterialKind.ADAPTER,
                            inferred=True,source=source+" (adapter load unresolved)")
    tn5=graph.input("tn5","Tn5 transposase",MaterialKind.CATALYST)
    tagged=graph.product("tagmented_cdna","tagmented SureCell cDNA",MaterialKind.MOLECULAR,
                         inferred=True,source=source+" (adapter sequence unresolved)")
    waste=graph.product("tagment_waste","fragments lacking the bead-derived end",
                        MaterialKind.WASTE)
    graph.react("tagment","Tagment double-stranded cDNA",StepKind.TAGMENTATION,
                ((duplex,InputRole.TEMPLATE),(transposome,InputRole.STRUCTURAL),
                 (tn5,InputRole.CATALYST)),(tagged,waste))
    graph.discard(waste,"fragments lacking the selected bead end are not amplified")
    tpp1=graph.input("tpp1","TPP1 primer (sequence unverified)",MaterialKind.OLIGO,inferred=True,
                     source=source+" (sequence unresolved)")
    n7=graph.input("n7xx","indexed N7xx primer",MaterialKind.OLIGO,inferred=True,
                   source=source+" (sequence unresolved)")
    polymerase=graph.input("pcr_polymerase","library PCR polymerase",MaterialKind.CATALYST)
    final=selected_fragment()
    library=graph.product("library",final.name,MaterialKind.LIBRARY,state=final,
                          inferred=True,source=source+" (all adapter sequences unresolved)")
    graph.react("library_pcr","Select bead-end fragments with TPP1 and N7xx",
                StepKind.PCR,
                ((tagged,InputRole.TEMPLATE),(tpp1,InputRole.STRUCTURAL),
                 (n7,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(library,),
                geometry=PCRGeometry.unavailable(
                    final,reason="The vendor names TPP1 and N7xx but does not publish their sequences or binding sites."))
    graph.sequence(library,SequencingPlan(
        acquired_roles=("Read 1","Index 1 (i7)","Read 2"),
        read_uses={"Read 1":"extract_then_include","Index 1 (i7)":"demultiplex",
                   "Read 2":"include"},
        limitation="The checked kit guide requires custom Sequencing Primer SP but does not specify its bases or run lengths; instrument-specific primer guides remain unavailable."))
    graph.declare_reviewed_scope((
        "cells","barcode_beads","droplet_device","partition","reverse_transcriptase",
        "reverse_transcribe","purification_beads","cleanup","rna_removal","remove_rna",
        "second_strand_primer","second_strand_polymerase","second_strand","transposome",
        "tn5","tagment","tpp1","n7xx","pcr_polymerase","library_pcr","library"),source)
    return graph.require_publishable()
