"""Source-backed bulk Multi-CUT&Tag molecular model.

The defining paper also adapts the chemistry to the 10x scATAC workflow. The exact
antibody-coded transposon and bulk sequencing geometry are printed in Table S1, whereas
the vendor cell-barcode oligos are not; this module therefore draws the published bulk
construct exactly and keeps the single-cell branch explicit but separate.
"""
from __future__ import annotations

import illumina as il
import nextera as nx
import seqprimers as sp
from batch_ngs import seg
from chemdraw import Construct, Scene, feature, revcomp
from tagmentation import loaded_tn5_adapter
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from protocol_model import molecule, role

TITLE = "Multi-CUT&Tag — antibody-coded targeted transposition"
NOTES = "01_multi-cut-tag.html"
SOURCE = ('Defining source: <a href="https://doi.org/10.1016/j.molcel.2021.09.019">'
          'Gopalan et al., <i>Molecular Cell</i> (2021)</a>, STAR Methods. Oligo sequences '
          'are from Supplemental Table 1 of the <a href="https://doi.org/10.1101/2021.04.27.441642">'
          'bioRxiv preprint</a>, because the published Table S1 (PMC supplement '
          'NIHMS1743381-supplement-2.pdf) could not be retrieved; the published Methods name the '
          'same A/B oligo series and custom primers.')
SUMMARY = ("Each antibody is coupled to pA–Tn5 loaded with a distinct eight-base adapter "
           "code. Both ends of every chromatin fragment consequently report which "
           "antibody-targeted transposomes made it.")
CAVEAT = ("The exact construct below is the paper’s bulk Multi-CUT&Tag branch. The "
          "single-cell branch passes the same target-coded fragments into the 10x Chromium "
          "Single Cell ATAC workflow. The 10x CG000209 guide prints its invariant gel-bead "
          "arms, and vendor software distributes its barcode whitelist; the exact barcode "
          "acquired by a particular cell varies by experiment. INFERRED: the bulk barcoded "
          "i5/i7 amplification primers are named in the Methods but not printed in the "
          "inspected sources, so the P5–i5 and i7–P7 arms of the bulk library are inferred; "
          "only their positions next to the custom Index2/Index1 primer sites are derived.")

# Table S1. The four entries in each family differ only in the antibody barcode.
READ1 = "TCGTCGGCAGCGTCTCCACGC"
READ2 = "GTCTCGTGGGCTCGGCTGTCCCTGTCC"
INDEX1 = "GGACAGGGACAGCCGAGCCCACGAGAC"
INDEX2 = "GCGTGGAGACGCTGCCGACGA"
CONNECTOR_A = "GCGATCGAGGACGGC"
CONNECTOR_B = "CACCGTCTCCGCCTC"
BARCODES_A = ("TATAGCCT", "ATAGAGGC", "CCTATCCT", "GGCTCTGA")
BARCODES_B = ("CGAGTAAT", "TCTCCGGA", "AATGAGCG", "GGAATCTC")
TN5_ME_REV = nx.ME_RC


def _target_feature(side: str):
    return feature(f"antibody_barcode_{side.lower()}", "feature_barcode", "fixed",
                   group="antibody_identity", part=side,
                   note="one of four Table S1 sequences")


def adapter_a(barcode: str = BARCODES_A[0]):
    if barcode not in BARCODES_A:
        raise ValueError("adapter A needs a published antibody barcode")
    return [seg("Read 1 priming site", READ1, "r1"),
            seg("A-side antibody barcode", barcode, "cbc",
                feature=_target_feature("A")),
            seg("connector A", CONNECTOR_A, "me"),
            seg("A mosaic end", nx.ME, "me")]


def adapter_b(barcode: str = BARCODES_B[1]):
    if barcode not in BARCODES_B:
        raise ValueError("adapter B needs a published antibody barcode")
    return [seg("Read 2 priming site", READ2, "r2"),
            seg("B-side antibody barcode", barcode, "cbc",
                feature=_target_feature("B")),
            seg("connector B", CONNECTOR_B, "me"),
            seg("B mosaic end", nx.ME, "me")]


def loaded_adapter_scene(parts, label: str) -> Scene:
    """A source-listed barcoded top oligo annealed to phosphorylated Tn5MErev."""
    sc = loaded_tn5_adapter(
        parts, mosaic_end=parts[-1].name, label=label,
        bottom_name="Tn5MErev", bottom_label="5′-phosphorylated Tn5MErev",
        bottom_sequence=TN5_ME_REV, mod5="p")
    barcode = next(s.name for s in parts if s.feature)
    sc.mark("Tn5 top", barcode, "8-nt antibody identity")
    return sc


def transferred_fragment() -> Construct:
    """Representative heterologous A/B product after gap fill.

    Different barcode choices at the two ends are deliberately independent: equal codes
    are homogeneous reads, while unequal codes are the paper's mixed reads.
    """
    left = adapter_a(BARCODES_A[0])
    right = adapter_b(BARCODES_B[1])
    right_rc = [seg(s.name + " reverse complement", revcomp(s.top), s.tag,
                    feature=s.feature) for s in reversed(right)]
    return Construct([*left, seg("target-proximal chromatin", "X" * 36,
                                 placeholder=True), *right_rc],
                     name="gap-filled antibody-coded fragment")


TRANSFERRED = transferred_fragment()


def _final_library() -> Construct:
    return Construct([
        seg("P5", il.P5, "p5", inferred=True),
        seg("i5", "J" * 8, "cbc", placeholder=True, inferred=True,
            feature=feature("sample_index_i5", "sample_index", "unknown")),
        *list(TRANSFERRED),
        seg("i7 reverse complement", "I" * 8, "cbc", placeholder=True, inferred=True,
            feature=feature("sample_index_i7", "sample_index", "unknown")),
        seg("P7 reverse complement", il.P7_RC, "p7", inferred=True),
    ], name="bulk Multi-CUT&Tag library")


FINAL_LIBRARY = _final_library()
SEQ_PRIMERS = (
    sp.custom("Read 1", "Multi-CUT&Tag Read1", READ1, "Gopalan et al. Table S1"),
    sp.custom("Index 1 (i7)", "Multi-CUT&Tag Index1", INDEX1,
              "Gopalan et al. Table S1"),
    sp.custom("Index 2 (i5)", "Multi-CUT&Tag Index2", INDEX2,
              "Gopalan et al. Table S1"),
    sp.custom("Read 2", "Multi-CUT&Tag Read2", READ2, "Gopalan et al. Table S1"),
)
_problems = sp.verify(FINAL_LIBRARY, SEQ_PRIMERS)
if _problems:
    raise ValueError("Multi-CUT&Tag final library: " + "; ".join(_problems))

FINAL_CAPTION = ("Representative bulk A/B library. The independently chosen A- and B-side "
                 "codes may identify the same antibody (homogeneous) or two antibodies "
                 "occupying the same fragment (mixed).")
SEQUENCING_INTRO = ("The four custom Table S1 primers land by exact sequence. Each genomic "
                    "read first reports an 8-nt antibody barcode, 15-nt connector and 19-nt "
                    "mosaic end: the 42 bases removed before genomic alignment.")
READ_LENGTHS = {"Read 1": 72, "Index 1 (i7)": 8,
                "Index 2 (i5)": 8, "Read 2": 72}
READ_USES = {"Read 1": "extract_then_include", "Index 1 (i7)": "demultiplex",
             "Index 2 (i5)": "demultiplex", "Read 2": "extract_then_include"}


def oligos():
    rows = []
    for i, bc in enumerate(BARCODES_A, 1):
        rows.append(f"<b>P5_i5_{i}_Universal_Connector_A</b> — 5′-{READ1}{bc}{CONNECTOR_A}{nx.ME}-3′")
    for i, bc in enumerate(BARCODES_B, 1):
        rows.append(f"<b>P7_i7_{i}_Universal_Connector_B</b> — 5′-{READ2}{bc}{CONNECTOR_B}{nx.ME}-3′")
    rows.extend((
        f"<b>Tn5MErev</b> — 5′-/5Phos/{TN5_ME_REV}-3′",
        f"<b>Read1 / Read2</b> — 5′-{READ1}-3′ / 5′-{READ2}-3′",
        f"<b>Index1 / Index2</b> — 5′-{INDEX1}-3′ / 5′-{INDEX2}-3′",
    ))
    return rows


PUBLISHED = ("Gopalan et al. 2021 PMC8604773 STAR Methods (pA-Tn5 loading; bulk "
             "multi-CUT&Tag and library preparation; sequencing; single-cell multi-CUT&Tag)")
TABLE_S1 = ("bioRxiv 2021.04.27.441642 v1 Supplemental Table 1 (published Table S1, "
            "NIHMS1743381-supplement-2.pdf, was not retrievable)")


def protocol_graph():
    graph = ProtocolGraph(TITLE, authority=SOURCE)
    adapters = graph.input(
        "loaded_ab_adapters", "Table S1 A/B barcoded oligos each annealed to /5Phos/Tn5MErev "
        "(45 µM, 95 °C 5 min, 0.1 °C/min to 25 °C)", MaterialKind.ADAPTER,
        source=TABLE_S1,
        state=(loaded_adapter_scene(adapter_a(), "adapter A"),
               loaded_adapter_scene(adapter_b(), "adapter B")))
    pa_tn5 = graph.input("pa_tn5", "6-His-tagged Protein A–Tn5 transposase", MaterialKind.CATALYST)
    antibodies = graph.input("coded_antibodies", "antibody of interest for each barcode pair",
                             MaterialKind.SELECTION)
    coupled = graph.product("ab_pa_tn5_crude", "adapter-loaded Ab·pA–Tn5 with free adapters "
                            "and uncomplexed antibody", MaterialKind.MOLECULAR,
                            state=adapters.state)
    graph.react("load_and_couple", "Load each A/B adapter pair onto pA–Tn5 (glycerol, RT 60 min), "
                "then couple to its antibody (~2:1, 4 °C 4 h/overnight)", StepKind.BINDING,
                ((adapters, InputRole.TEMPLATE), (pa_tn5, InputRole.CATALYST),
                 (antibodies, InputRole.SELECTION)), (coupled,), changes_molecule=False)
    his_beads = graph.input("his_tag_beads", "Dynabeads His-Tag pulldown, 300 mM imidazole "
                            "elution, Amicon buffer exchange", MaterialKind.SELECTION)
    complexes = graph.product("ab_pa_tn5", "purified Ab·pA–Tn5 complexes", MaterialKind.MOLECULAR,
                              state=adapters.state)
    free = graph.product("free_adapters_antibody", "uncomplexed antibody and free adapters",
                         MaterialKind.WASTE)
    graph.react("purify_complexes", "His-tag purify Ab·pA–Tn5 complexes", StepKind.CLEANUP,
                ((coupled, InputRole.TEMPLATE), (his_beads, InputRole.SELECTION)),
                (complexes, free), changes_molecule=False)
    graph.discard(free, "is removed by His-tag pulldown before the complexes meet cells")
    nuclei = graph.input(
        "nuclei", "cells (bulk: 100,000 on concanavalin A beads; single-cell: 200,000 free)",
        MaterialKind.BIOLOGICAL,
        state=molecule("native chromatin", (role("multiple antibody-defined epitopes"),)),
        primary=True)
    bound = graph.product("bound_nuclei", "cells with bound and unbound Ab·pA–Tn5 complexes",
                          MaterialKind.MOLECULAR, state=nuclei.state)
    graph.react("target_chromatin", "Bind complexes overnight at 4 °C in Dig-wash + 2 mM EDTA",
                StepKind.BINDING, ((nuclei, InputRole.TEMPLATE), (complexes, InputRole.STRUCTURAL)),
                (bound,), changes_molecule=False)
    washes = graph.input("target_washes", "three washes (bulk: Dig-med on magnet; single-cell: "
                         "Dig-wash, 600 g pelleting)", MaterialKind.SELECTION)
    targeted = graph.product("targeted_nuclei", "washed cells with bound Ab·pA–Tn5 complexes",
                             MaterialKind.MOLECULAR, state=nuclei.state)
    wash = graph.product("unbound_targeting_reagents", "unbound Ab·pA–Tn5 complexes",
                         MaterialKind.WASTE)
    graph.react("wash_unbound", "Wash away unbound complexes", StepKind.CLEANUP,
                ((bound, InputRole.TEMPLATE), (washes, InputRole.SELECTION)),
                (targeted, wash), changes_molecule=False)
    graph.discard(wash, "is removed before magnesium activation")
    magnesium = graph.input("magnesium", "10 mM MgCl2 in Dig-med buffer, 37 °C 1 h",
                            MaterialKind.CONDITION)
    tagged = graph.product("target_coded_fragments", "tagmented antibody-coded fragments "
                           "(9-nt gaps unfilled)", MaterialKind.MOLECULAR, state=TRANSFERRED)
    graph.react("tagment", "Transfer antibody-coded adapters at target-proximal chromatin",
                StepKind.TAGMENTATION,
                ((targeted, InputRole.TEMPLATE), (complexes, InputRole.STRUCTURAL),
                 (pa_tn5, InputRole.CATALYST), (magnesium, InputRole.CONDITION)),
                (tagged,))
    bulk = graph.product("bulk_branch", "bulk tagmented fragments",
                         MaterialKind.MOLECULAR, state=TRANSFERRED)
    single_cell = graph.product("single_cell_branch", "single-cell tagmentation reaction",
                                MaterialKind.MOLECULAR, state=TRANSFERRED)
    graph.react("split_branches", "Bulk and single-cell experiments are separate preparations",
                StepKind.SPLIT, ((tagged, InputRole.TEMPLATE),), (bulk, single_cell),
                changes_molecule=False)

    extraction = graph.input("stop_extract", "EDTA/SDS/proteinase K 55 °C 1 h; phenol:chloroform, "
                             "ethanol precipitation", MaterialKind.SELECTION)
    extracted = graph.product("extracted_dna", "extracted tagmented DNA", MaterialKind.MOLECULAR,
                              state=TRANSFERRED)
    protein_waste = graph.product("extraction_waste", "digested protein, beads and phases",
                                  MaterialKind.WASTE)
    graph.react("extract_dna", "Stop tagmentation and extract DNA", StepKind.EXTRACTION,
                ((bulk, InputRole.TEMPLATE), (extraction, InputRole.SELECTION)),
                (extracted, protein_waste), changes_molecule=False)
    graph.discard(protein_waste, "phenol phase and supernatants are discarded")
    hifi = graph.input("nebnext_hifi", "NEBNext HiFi 2x PCR Master Mix", MaterialKind.CATALYST)
    filled = graph.product("gap_filled", "gap-filled antibody-coded fragments",
                           MaterialKind.MOLECULAR, state=TRANSFERRED)
    graph.react("gap_fill", "Fill the 9-nt Tn5 gaps: 72 °C 5 min in the PCR mix", StepKind.EXTENSION,
                ((extracted, InputRole.TEMPLATE), (hifi, InputRole.CATALYST)), (filled,))
    bulk_forward = graph.input(
        "bulk_forward", "barcoded i5 PCR primer (sequence not in inspected sources)",
        MaterialKind.OLIGO, inferred=True,
        source="STAR Methods name barcoded i5/i7 primers; preprint Table S1 lists none and the "
               "published Table S1 was not retrievable")
    bulk_reverse = graph.input(
        "bulk_reverse", "barcoded i7 PCR primer (sequence not in inspected sources)",
        MaterialKind.OLIGO, inferred=True,
        source="STAR Methods name barcoded i5/i7 primers; preprint Table S1 lists none and the "
               "published Table S1 was not retrievable")
    amplified = graph.product("bulk_amplified", FINAL_LIBRARY.name + " (in PCR mix)",
                              MaterialKind.MOLECULAR, state=FINAL_LIBRARY)
    graph.react("bulk_pcr", "INFERRED — barcoded i5/i7 primer arms not printed; PCR: 98 °C 30 s; 17 cycles 98 °C 10 s / 63 °C 30 s; 72 °C 1 min",
                StepKind.PCR,
                ((filled, InputRole.TEMPLATE), (bulk_forward, InputRole.STRUCTURAL),
                 (bulk_reverse, InputRole.STRUCTURAL), (hifi, InputRole.CATALYST)), (amplified,),
                geometry=PCRGeometry.unavailable(
                    TRANSFERRED,
                    reason="The barcoded i5/i7 amplification primers are not printed in the inspected "
                           "sources; the P5/i5 and i7/P7 arms are inferred from the custom index primers."))
    ampure = graph.input("ampure_bulk", "1.1x AMPure XP, two 80% ethanol washes", MaterialKind.SELECTION)
    bulk_library = graph.product("bulk_library", FINAL_LIBRARY.name,
                                 MaterialKind.LIBRARY, state=FINAL_LIBRARY)
    pcr_waste = graph.product("bulk_pcr_waste", "PCR mix and unbound primers", MaterialKind.WASTE)
    graph.react("bulk_cleanup", "Clean up the PCR with 1.1x AMPure XP", StepKind.CLEANUP,
                ((amplified, InputRole.TEMPLATE), (ampure, InputRole.SELECTION)),
                (bulk_library, pcr_waste), changes_molecule=False)
    graph.discard(pcr_waste, "bead supernatant and washes are discarded")
    graph.sequence(bulk_library, SequencingPlan(
        library=FINAL_LIBRARY, primers=tuple(SEQ_PRIMERS),
        acquired_roles=tuple(p.role for p in SEQ_PRIMERS), read_lengths=READ_LENGTHS,
        read_uses=READ_USES))

    gel_bead = graph.input(
        "gel_bead_oligo", "10x scATAC gel-bead cell-barcode oligo",
        MaterialKind.OLIGO, inferred=True,
        source="10x CG000209 v1.1 guide and Cell Ranger ATAC barcode whitelist; experiment-specific barcode assignment")
    sc_polymerase = graph.input("single_cell_polymerase", "Chromium barcoding master mix polymerase",
                                MaterialKind.CATALYST)
    cell_barcoded = graph.product(
        "cell_barcoded", "10x-cell-barcoded target-coded fragments",
        MaterialKind.MOLECULAR,
        state=molecule("single-cell Multi-CUT&Tag product", (
            role("10x cell barcode"), role("A/B antibody-coded chromatin fragment"))),
        inferred=True, source="10x CG000209 v1.1 guide; experiment-specific barcode assignment")
    graph.react("cell_barcode", "In GEMs: 72 °C 5 min, then only 1 cycle of 98/59/72 °C (vendor: >10)",
                StepKind.EXTENSION,
                ((single_cell, InputRole.TEMPLATE), (gel_bead, InputRole.STRUCTURAL),
                 (sc_polymerase, InputRole.CATALYST)), (cell_barcoded,))
    silane = graph.input("silane_ampure", "Dynabeads MyOne SILANE and AMPure XP (vendor protocol)",
                         MaterialKind.SELECTION)
    gem_clean = graph.product("gem_cleaned", "post-GEM purified cell-barcoded fragments",
                              MaterialKind.MOLECULAR, state=cell_barcoded.state, inferred=True,
                              source="10x CG000209 v1.1 guide; experiment-specific barcode assignment")
    gem_waste = graph.product("gem_waste", "GEM oil, reagents and bead supernatants", MaterialKind.WASTE)
    graph.react("post_gem_cleanup", "Post-GEM cleanup", StepKind.CLEANUP,
                ((cell_barcoded, InputRole.TEMPLATE), (silane, InputRole.SELECTION)),
                (gem_clean, gem_waste), changes_molecule=False)
    graph.discard(gem_waste, "is discarded")
    sc_forward = graph.input("sc_forward", "10x scATAC SI-PCR forward primer",
                             MaterialKind.OLIGO, inferred=True,
                             source="10x CG000209 v1.1 guide, Appendix A; SI-PCR Primer B forward arm")
    sc_reverse = graph.input("sc_reverse", "10x Single Index Kit reverse primer",
                             MaterialKind.OLIGO, inferred=True,
                             source="10x CG000209 v1.1 guide, Appendix A; variable-index SI-PCR reverse arm")
    sc_pcr_polymerase = graph.input("sc_pcr_polymerase", "single-cell library PCR polymerase",
                                    MaterialKind.CATALYST)
    sc_amplified = graph.product(
        "single_cell_amplified", "amplified single-cell library before size selection",
        MaterialKind.MOLECULAR, inferred=True,
        source="10x CG000209 v1.1 guide, Appendix A; index and cell barcode vary by experiment")
    graph.react("single_cell_pcr", "Sample-index PCR, 14 cycles", StepKind.PCR,
                ((gem_clean, InputRole.TEMPLATE), (sc_forward, InputRole.STRUCTURAL),
                 (sc_reverse, InputRole.STRUCTURAL),
                 (sc_pcr_polymerase, InputRole.CATALYST)), (sc_amplified,),
                geometry=PCRGeometry.unavailable(
                    TRANSFERRED,
                    reason="The paper uses the 10x scATAC workflow. Its CG000209 guide prints the invariant gel-bead and SI-PCR arms; this graph has not yet assembled the kit's variable barcode and index choices."))
    double_select = graph.input("ampure_double", "AMPure XP double size selection (vendor protocol)",
                                MaterialKind.SELECTION)
    sc_library = graph.product(
        "single_cell_library", "single-cell Multi-CUT&Tag sequencing library",
        MaterialKind.LIBRARY, inferred=True,
        source="10x CG000209 v1.1 guide, Appendix A; index and cell barcode vary by experiment")
    size_waste = graph.product("size_selection_waste", "fragments outside the double size window",
                               MaterialKind.WASTE)
    graph.react("single_cell_size_selection", "Double-sided AMPure size selection",
                StepKind.CLEANUP, ((sc_amplified, InputRole.TEMPLATE),
                                   (double_select, InputRole.SELECTION)),
                (sc_library, size_waste), changes_molecule=False)
    graph.discard(size_waste, "is discarded with the bead fractions")
    graph.sequence(sc_library, SequencingPlan(
        acquired_roles=("Read 1", "Index 1 (i7)", "Index 2 (i5)", "Read 2"),
        limitation="The paper establishes the run roles (72/8/16/72 cycles) and the 10x CG000209 guide prints the kit architecture; this model has not yet assembled its variable barcode and index choices."))
    graph.declare_reviewed_scope((
        "loaded_ab_adapters", "pa_tn5", "coded_antibodies", "load_and_couple",
        "his_tag_beads", "purify_complexes", "nuclei", "target_chromatin", "target_washes",
        "wash_unbound", "magnesium",
        "tagment", "split_branches", "stop_extract", "extract_dna", "nebnext_hifi",
        "gap_fill", "bulk_forward", "bulk_reverse", "bulk_pcr", "ampure_bulk", "bulk_cleanup", "bulk_library",
        "gel_bead_oligo", "single_cell_polymerase", "cell_barcode", "silane_ampure",
        "post_gem_cleanup", "sc_forward", "sc_reverse", "sc_pcr_polymerase", "single_cell_pcr", "ampure_double",
        "single_cell_size_selection", "single_cell_library"), PUBLISHED)
    return graph.require_publishable()
