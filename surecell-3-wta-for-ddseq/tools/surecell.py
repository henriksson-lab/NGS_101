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
from protocol_graph import (Disposition,InputRole,MaterialKind,ProtocolGraph,SequencingPlan,
                            StepKind,Terminal)


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
    guide="Reference Guide 1000000021452 v01"
    graph=ProtocolGraph("SureCell WTA 3' for ddSEQ",authority=source)
    cells=graph.input("cells","single-cell suspension",MaterialKind.BIOLOGICAL,primary=True)
    beads=graph.input("barcode_beads","3-prime Barcode Mix beads",MaterialKind.ADAPTER,
                      state=bead_primer(),inferred=True,
                      source=source+" (oligo sequences unresolved)")
    droplets=graph.input("droplet_device","ddSEQ Single-Cell Isolator cartridge (Cell and Barcode Suspension Mixes)",
                         MaterialKind.CONDITION)
    captured=graph.product("captured_rna","droplet-confined bead-captured poly(A) RNA",
                           MaterialKind.MOLECULAR,inferred=True,
                           source=source+" (barcode anatomy unresolved)")
    empty=graph.product("empty_droplets","droplets without a productive cell-bead pairing",
                        MaterialKind.MOLECULAR)
    graph.react("partition","Co-encapsulate cells, lysis/RT mix and barcode beads",
                StepKind.PARTITION,
                ((cells,InputRole.TEMPLATE),(beads,InputRole.STRUCTURAL),
                 (droplets,InputRole.CONDITION)),(captured,empty),
                source=f"{guide}, Prepare Cell/Barcode Suspension Mixes and Load Cartridge")
    graph.mark_terminal(Terminal(empty.id,Disposition.UNRESOLVED,
        "Empty or bead-free droplets are not physically separated; the emulsion is broken in bulk and only short unbound barcodes are named as removed by the SPB clean-up."))
    rt=graph.input("reverse_transcriptase","reverse transcriptase in the Cell Suspension Mix (enzyme not named)",
                   MaterialKind.CATALYST)
    first=graph.product("first_strand",first_strand().name,MaterialKind.MOLECULAR,
                        state=first_strand(),inferred=True,
                        source=source+" (barcode anatomy unresolved)")
    graph.react("reverse_transcribe","Reverse-transcribe in droplets (37 °C 30 min, 50 °C 60 min, 85 °C 5 min)",
                StepKind.REVERSE_TRANSCRIPTION,
                ((captured,InputRole.TEMPLATE),(beads,InputRole.STRUCTURAL),
                 (rt,InputRole.CATALYST)),(first,),source=f"{guide}, Reverse Transcription (RT program)")
    disruptor=graph.input("droplet_disruptor","Droplet Disruptor",MaterialKind.CONDITION)
    purification=graph.input("purification_beads","Purification Beads (SPB)",
                             MaterialKind.SELECTION)
    cleaned=graph.product("cleaned_first_strand","purified bead-barcoded first-strand cDNA",
                          MaterialKind.MOLECULAR,state=first_strand(),inferred=True,
                          source=source+" (barcode anatomy unresolved)")
    free_reagent=graph.product("free_barcode_reagent","oil layer and SPB supernatant with short fragments (unbound barcodes)",
                               MaterialKind.WASTE)
    graph.react("cleanup","Break emulsion; 90 µl SPB clean-up removes unbound barcodes and combines both outputs",
                StepKind.CLEANUP,
                ((first,InputRole.TEMPLATE),(disruptor,InputRole.CONDITION),
                 (purification,InputRole.SELECTION)),(cleaned,free_reagent),
                source=f"{guide}, Break Emulsion and Clean Up First Strand Synthesis")
    graph.discard(free_reagent,"SPB supernatant (short fragments/unbound barcodes) is removed and discarded")
    ssb=graph.input("second_strand_buffer","Second Strand Buffer (SSB)",MaterialKind.CONDITION)
    sse=graph.input("second_strand_enzyme","Second Strand Enzyme (SSE; enzymes not named)",
                    MaterialKind.CATALYST)
    duplex=graph.product("double_stranded_cdna","replacement-strand SureCell cDNA",
                         MaterialKind.MOLECULAR,state=first_strand(),inferred=True,
                         source=source+" (second-strand sequence unresolved)")
    graph.react("second_strand","Remove RNA template and synthesize the replacement strand (16 °C, 2 h)",
                StepKind.EXTENSION,
                ((cleaned,InputRole.TEMPLATE),(ssb,InputRole.CONDITION),(sse,InputRole.CATALYST)),
                (duplex,),source=f"{guide}, Synthesize Second Strand cDNA")
    clean_cdna=graph.product("clean_cdna","purified double-stranded cDNA",MaterialKind.MOLECULAR,
                             state=first_strand(),inferred=True,
                             source=source+" (second-strand sequence unresolved)")
    cdna_waste=graph.product("cdna_cleanup_waste","SPB supernatant with short fragments",
                             MaterialKind.WASTE)
    graph.react("cdna_cleanup","44 µl SPB clean-up of cDNA (removes short fragments)",StepKind.CLEANUP,
                ((duplex,InputRole.TEMPLATE),(purification,InputRole.SELECTION)),(clean_cdna,cdna_waste),
                source=f"{guide}, Clean Up cDNA")
    graph.discard(cdna_waste,"SPB supernatant removed and discarded")
    transposome=graph.input("transposome","Nextera SureCell transposome (TCE) in Tagment Buffer (TCB)",
                            MaterialKind.ADAPTER,
                            inferred=True,source=source+" (adapter load unresolved)")
    tn5=graph.input("tn5","Tn5 transposase",MaterialKind.CATALYST)
    tagged=graph.product("tagmented_cdna","tagmented SureCell cDNA",MaterialKind.MOLECULAR,
                         inferred=True,source=source+" (adapter sequence unresolved)")
    graph.react("tagment","Tagment double-stranded cDNA (55 °C, 5 min)",StepKind.TAGMENTATION,
                ((clean_cdna,InputRole.TEMPLATE),(transposome,InputRole.STRUCTURAL),
                 (tn5,InputRole.CATALYST)),(tagged,),source=f"{guide}, Tagment cDNA")
    tsb=graph.input("tagment_stop_buffer","Tagment Stop Buffer (TSB)",MaterialKind.CONDITION)
    stopped=graph.product("stopped_tagmentation","stopped tagmentation reaction",MaterialKind.MOLECULAR,
                          inferred=True,source=source+" (adapter sequence unresolved)")
    graph.react("stop_tagment","Add Tagment Stop Buffer, 5 min at room temperature",StepKind.QUENCHING,
                ((tagged,InputRole.TEMPLATE),(tsb,InputRole.CONDITION)),(stopped,),
                changes_molecule=False,source=f"{guide}, Tagment cDNA steps 7-9")
    tpp1=graph.input("tpp1","Tagment PCR Adapter TPP1 (sequence unverified)",MaterialKind.OLIGO,inferred=True,
                     source=source+" (sequence unresolved)")
    n7=graph.input("n7xx","DNA Adapter N7xx, one index per sample",MaterialKind.OLIGO,inferred=True,
                   source=source+" (sequence unresolved)")
    polymerase=graph.input("pcr_polymerase","Tagmentation PCR Mix (TPM)",MaterialKind.CATALYST)
    final=selected_fragment()
    amplified=graph.product("amplified_library","amplified bead-end library",MaterialKind.MOLECULAR,
                            state=final,inferred=True,source=source+" (all adapter sequences unresolved)")
    unselected=graph.product("unselected_fragments","tagmented fragments lacking the TPP1-selected bead end",
                             MaterialKind.MOLECULAR)
    graph.react("library_pcr","Select bead-end fragments with TPP1 and N7xx (15 cycles)",
                StepKind.PCR,
                ((stopped,InputRole.TEMPLATE),(tpp1,InputRole.STRUCTURAL),
                 (n7,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(amplified,unselected),
                source=f"{guide}, Amplify Tagmented cDNA (LA program)",
                geometry=PCRGeometry.unavailable(
                    final,reason="The vendor names TPP1 and N7xx but does not publish their sequences or binding sites."))
    graph.mark_terminal(Terminal(unselected.id,Disposition.UNRESOLVED,
        "Not exponentially amplified by TPP1/N7xx; the guide describes no physical removal, so their fate in the cleaned library is not assigned."))
    once=graph.product("once_cleaned_library","library after the first SPB clean-up",MaterialKind.MOLECULAR,
                       state=final,inferred=True,source=source+" (all adapter sequences unresolved)")
    library=graph.product("library",final.name,MaterialKind.LIBRARY,state=final,
                          inferred=True,source=source+" (all adapter sequences unresolved)")
    for ident,inp,out,label,section in (
            ("library_cleanup",amplified,once,"58 µl SPB bind, two 80% EtOH washes, elute 51 µl RSB (keep 50 µl)","Clean Up Libraries: Bind/Wash/Elute"),
            ("library_second_cleanup",once,library,"Second bind with 30 µl SPB, two washes, elute 22 µl RSB (keep 20 µl)","Clean Up Libraries: Second Bind/Wash/Elute")):
        waste=graph.product(f"{ident}_waste","SPB supernatant with short library fragments",MaterialKind.WASTE)
        graph.react(ident,label,StepKind.CLEANUP,((inp,InputRole.TEMPLATE),(purification,InputRole.SELECTION)),
                    (out,waste),source=f"{guide}, {section}")
        graph.discard(waste,"SPB supernatant removed and discarded")
    graph.sequence(library,SequencingPlan(
        acquired_roles=("Read 1","Index 1 (i7)","Read 2"),
        read_uses={"Read 1":"extract_then_include","Index 1 (i7)":"demultiplex",
                   "Read 2":"include"},
        limitation="The checked kit guide requires custom Sequencing Primer SP but does not specify its bases or run lengths; instrument-specific primer guides remain unavailable."))
    graph.declare_reviewed_scope((
        "cells","barcode_beads","droplet_device","partition","reverse_transcriptase",
        "reverse_transcribe","droplet_disruptor","purification_beads","cleanup",
        "second_strand_buffer","second_strand_enzyme","second_strand","cdna_cleanup",
        "transposome","tn5","tagment","tagment_stop_buffer","stop_tagment","tpp1","n7xx",
        "pcr_polymerase","library_pcr","library_cleanup","library_second_cleanup","library"),source)
    return graph.require_publishable()
