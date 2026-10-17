"""Open-chromatin-anchored chromosome-conformation library chemistry.

HiCAR deliberately joins two unlike sequencing entries: a TruSeq Read 1 side reports the
restriction-fragment partner and a Nextera Read 2 side reports the Tn5-tagged accessible
anchor.  This module owns that hybrid construct and validates all four run primers.
"""
from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import illumina as il
import nextera as nx
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, feature, inline_sequence, revcomp
from circular import CircularizationGeometry
from pcr import PCRGeometry
from protocol_graph import (Disposition, Terminal, InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from protocol_model import chain, closed_circle, molecule, role


HICAR7_I0_WRITTEN = (
    "/5Phos/taGGTCTGAGTTCGCACAAAGGGAAGAGCGTCGTGTGGAGATGTGTATAAGAGACAG"
)
HICAR7_I0 = "TAGGTCTGAGTTCGCACAAAGGGAAGAGCGTCGTGTGG" + nx.ME
TRUSEQ_R1_SPLINT = "ACACGACGCTCTTCCGATCT"
I5 = il.NEBNEXT_I5_SET1["i501"]
I7 = "CAGCCTCG"
I501_PCR = il.nebnext_i5_primer(I5)
NEXTERA_I7_10_L = il.P7 + I7 + nx.S7 + nx.ME

SAMPLE_I5 = feature("hicar.i5", "sample_index", "fixed",
                    note="NEB i501 index carried by the defining library primer")
SAMPLE_I7 = feature("hicar.i7", "sample_index", "fixed",
                    note="Nextera i7-10-L index carried by the defining library primer")


def _seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    # Exact displayed arms come from the later protocol; application to 2022 is inferred.
    kw.setdefault("inferred", True)
    return Segment(name, top, tag, **kw)


def hicar_library() -> tuple[Construct, tuple]:
    """Finished HiCAR amplicon with its asymmetric primer geometry enforced."""
    lib = Construct([
        _seg("P5", il.P5, "p5"),
        _seg("i5 (i501)", I5, "cbc", feature=SAMPLE_I5),
        _seg("TruSeq Read 1 arm", il.TRUSEQ_READ1, "r1"),
        role("CviQI-cut interacting genomic DNA"),
        role("HiCAR7–splint proximity-ligation junction"),
        role("Tn5-tagged accessible genomic DNA"),
        _seg("mosaic end reverse complement", nx.ME_RC, "me"),
        _seg("S7 reverse complement", nx.S7_RC, "s7"),
        _seg("i7 reverse complement", revcomp(I7), "cbc", feature=SAMPLE_I7),
        _seg("P7 reverse complement", il.P7_RC, "p7"),
    ], name="HiCAR open-to-all contact library")
    primers = (sp.TRUSEQ["R1"], sp.NEXTERA["I1"],
               sp.TRUSEQ["I2"], sp.NEXTERA["R2"])
    problems = sp.verify(lib, primers, required_roles=tuple(p.role for p in primers))
    if problems:
        raise ValueError("HiCAR hybrid library: " + "; ".join(problems))
    return lib, primers


def _contact_scene(label: str) -> Scene:
    return molecule(label, (
        role("CviQI-cut interacting genomic DNA"),
        role("HiCAR7 adapter and splint"),
        role("Tn5-tagged accessible genomic DNA"),
    ), junctions=(("CviQI-cut interacting genomic DNA", "HiCAR7 adapter and splint",
                  "proximity ligation"),
                 ("HiCAR7 adapter and splint", "Tn5-tagged accessible genomic DNA",
                  "Tn5 strand transfer")))


def hicar_workflow():
    initial = molecule("formaldehyde-crosslinked nuclei", (
        role("spatially proximal chromatin"), role("accessible regulatory DNA")))
    circle = closed_circle("NlaIII-trimmed HiCAR circle", (
        role("CviQI-cut interacting genomic DNA"),
        role("HiCAR7–splint junction"),
        role("Tn5-tagged accessible genomic DNA"),
    ), closure_label="dilute intramolecular T4 ligation")
    return chain(initial, (
        ("Tagment accessible chromatin with Tn5 loaded with engineered HiCAR7",
         molecule("tagged accessible chromatin", (
             role("chromatin outside the accessible site"),
             role("HiCAR7 transferred strand"),
             role("accessible genomic DNA")),
             junctions=(("HiCAR7 transferred strand", "accessible genomic DNA",
                         "Tn5 strand transfer"),)),
         "HiCAR7 carries the Tn5 mosaic end and a single-stranded 5-prime TA-compatible flank; the defining protocol starts with 30,000–100,000 fixed cells."),
        ("Digest nuclei with CviQI and splint-ligate spatial partners in situ",
         _contact_scene("open-to-all proximity-ligation product"),
         "The TruseqR1 splint bridges the engineered adapter to a CviQI-cut partner. This is the contact-recording covalent junction; Read 1 and Read 2 later report opposite sides."),
        ("Reverse crosslinks and fill the staggered Tn5 gaps with T4 DNA polymerase",
         _contact_scene("repaired linear HiCAR contact"),
         "Gap repair completes the adapter-bearing duplex before the secondary restriction digest."),
        ("Digest purified DNA with NlaIII and circularize at low concentration",
         circle,
         "NlaIII shortens the contact product; dilute T4 ligation closes individual molecules in the 4C-like library-preparation step."),
        ("Treat the circularized material with PmeI before PCR",
         molecule("PmeI-treated HiCAR template pool", (
             role("circularized contact molecules and PmeI cleavage products"),)),
         "The defining Methods, step 5, digest the purified circle pool with PmeI for 1 hour at 37 °C before adding the inverse-PCR mix; the source does not identify every resulting molecular species."),
        ("Inverse-PCR between the splint side and mosaic-end side for 12 cycles",
         molecule("linear dual-index HiCAR amplicon", (
             role("P5, i501 and TruSeq Read 1 arm"),
             role("CviQI-cut interacting genomic DNA"),
             role("HiCAR7–splint contact junction"),
             role("Tn5-tagged accessible genomic DNA"),
             role("mosaic end, S7, i7-10-L and P7"))),
         "Outward primers open the selected circle into the sequencing library. The TruSeq side produces R1 from the interacting partner and the Nextera side produces R2 from accessible chromatin."),
    ))


def oligos():
    return (
        f"Later 2023 HiCAR7-i0 (2022 source-version assignment unresolved; 5-prime phosphate retained) — {HICAR7_I0_WRITTEN}",
        "TruseqR1 splint — exact defining Table S1 oligo remains unresolved in this review",
        f"Later 2023 NEB primer i501 — {inline_sequence(I501_PCR, 'NEB i501 PCR primer')}",
        f"Later 2023 Nextera-pcr-i7-10-L — {inline_sequence(NEXTERA_I7_10_L, 'Nextera i7-10-L PCR primer')}",
    )


def _protocol_graph(final, primers) -> ProtocolGraph:
    source = "Wei et al., Molecular Cell (2022); doi:10.1016/j.molcel.2022.01.023"
    graph = ProtocolGraph("HiCAR", authority=source)
    nuclei = graph.input(
        "fixed_nuclei", "formaldehyde-crosslinked nuclei",
        MaterialKind.BIOLOGICAL,
        state=molecule("formaldehyde-crosslinked nuclei", (
            role("spatially proximal chromatin"), role("accessible regulatory DNA"))),
        primary=True)
    hicar_adapter = graph.input("hicar7_transposome", "engineered HiCAR7-loaded transposome",
                                MaterialKind.ADAPTER, inferred=True,
                                source="Lu et al., 2023 Table 1/2; 2022 Table S1 unavailable",
                                state=Construct([_seg("5-prime phosphorylated HiCAR7-i0 (later protocol)", HICAR7_I0)], name="later HiCAR7-i0 representative"))
    tn5 = graph.input("tn5", "Tn5 transposase", MaterialKind.CATALYST)
    tagged = graph.product(
        "tagged_chromatin", "HiCAR7-tagged accessible chromatin",
        MaterialKind.MOLECULAR,
        state=molecule("tagged accessible chromatin", (
            role("chromatin outside the accessible site"),
            role("HiCAR7 transferred strand"), role("accessible genomic DNA")),
            junctions=(("HiCAR7 transferred strand", "accessible genomic DNA",
                        "Tn5 strand transfer"),)))
    tn5_losses = graph.product("tn5_losses", "candidate Tn5 products outside the represented contact route",
                               MaterialKind.MOLECULAR, inferred=True)
    graph.react("tagment", "INFERRED — tagment chromatin; displayed adapter is later-source i0",
                StepKind.TAGMENTATION,
                ((nuclei, InputRole.TEMPLATE), (hicar_adapter, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged, tn5_losses))
    graph.mark_terminal(Terminal(tn5_losses.id, Disposition.UNRESOLVED,
        "Candidate products, not an experimentally measured population. Lack of the selected anchor configuration does not establish physical removal; downstream propagation is not traced here."))
    quench = graph.input("tagment_quench", "EDTA stop, nuclear wash and SDS/Triton treatment", MaterialKind.CONDITION)
    ready = graph.product("digestion_ready", "washed, permeabilized tagmented nuclei", MaterialKind.MOLECULAR, state=tagged.state)
    graph.react("stop_tagmentation", "Stop Tn5 with EDTA; wash and permeabilize nuclei before CviQI", StepKind.MAGIC,
                ((tagged, InputRole.TEMPLATE), (quench, InputRole.CONDITION)), (ready,),
                reason="Defining Step 1 summarizes EDTA stop, BSA wash, SDS treatment and Triton neutralization; exact Tn5 displacement/adapter-end state is not resolved here.")

    cviqi = graph.input("cviqi", "CviQI restriction endonuclease", MaterialKind.CATALYST)
    digested = graph.product("cviqi_fragments", "CviQI-digested HiCAR chromatin",
                             MaterialKind.MOLECULAR, state=tagged.state)
    graph.react("cviqi_digest", "Digest nuclei with CviQI and TruseqR1 oligo for 2 h at room temperature", StepKind.DIGESTION,
                ((ready, InputRole.TEMPLATE), (cviqi, InputRole.CATALYST)),
                (digested,))
    splint = graph.input("truseq_r1_splint", "published TruseqR1 ligation splint",
                         MaterialKind.OLIGO)
    ligase = graph.input("proximity_ligase", "T4 DNA ligase (in-situ ligation)",
                         MaterialKind.CATALYST)
    contact = graph.product(
        "contact", "open-to-all proximity-ligation product",
        MaterialKind.MOLECULAR,
        state=_contact_scene("open-to-all proximity-ligation product"))
    graph.react("proximity_ligation", "Splint-ligate spatial partners with T4 ligase for 4 h at room temperature",
                StepKind.PROXIMITY_LIGATION,
                ((digested, InputRole.TEMPLATE), (splint, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (contact,))
    reversal = graph.input("crosslink_reversal", "68 °C crosslink reversal, ethanol precipitation and 80% ethanol wash",
                           MaterialKind.CONDITION)
    linear = graph.product("linear_contact", "purified linear HiCAR contact",
                           MaterialKind.MOLECULAR, state=contact.state)
    graph.react("reverse_crosslinks", "Reverse crosslinks at 68 °C and recover DNA by ethanol precipitation",
                StepKind.MAGIC,
                ((contact, InputRole.TEMPLATE), (reversal, InputRole.CONDITION)),
                (linear,), reason="Defining Step 3 states nuclear pelleting, crosslink reversal and ethanol recovery; these operations are summarized without assigning exact chemical termini.")
    fill_polymerase = graph.input("gap_fill_polymerase", "T4 DNA polymerase",
                                  MaterialKind.CATALYST)
    repaired = graph.product("repaired_contact", "Tn5-gap-filled HiCAR contact",
                             MaterialKind.MOLECULAR, state=contact.state)
    graph.react("gap_fill", "Fill Tn5 gaps with T4 DNA polymerase for 30 min at room temperature", StepKind.EXTENSION,
                ((linear, InputRole.TEMPLATE), (fill_polymerase, InputRole.CATALYST)),
                (repaired,))
    heat = graph.input("polymerase_heat_stop", "75 °C for 20 min", MaterialKind.CONDITION)
    inactivated = graph.product("heat_stopped", "gap-repaired DNA after polymerase heat inactivation", MaterialKind.MOLECULAR, state=repaired.state)
    graph.react("inactivate_polymerase", "Heat-inactivate T4 DNA polymerase before NlaIII", StepKind.THERMAL_TREATMENT,
                ((repaired, InputRole.TEMPLATE), (heat, InputRole.CONDITION)), (inactivated,))
    nlaiii = graph.input("nlaiii", "NlaIII restriction endonuclease", MaterialKind.CATALYST)
    shortened = graph.product("nlaiii_fragments", "NlaIII-trimmed HiCAR contact fragments",
                              MaterialKind.MOLECULAR, state=repaired.state)
    graph.react("nlaiii_digest", "Digest repaired contacts with NlaIII",
                StepKind.DIGESTION,
                ((inactivated, InputRole.TEMPLATE), (nlaiii, InputRole.CATALYST)),
                (shortened,))
    spri = graph.input("spri_condition", "0.9X SPRI purification", MaterialKind.CONDITION)
    purified_short = graph.product("purified_short", "SPRI-purified NlaIII fragments", MaterialKind.MOLECULAR, state=shortened.state)
    spri_waste = graph.product("spri_waste", "unrecovered material from NlaIII SPRI purification", MaterialKind.WASTE)
    graph.react("purify_nlaiii", "Purify NlaIII digest with 0.9X SPRI beads", StepKind.CLEANUP,
                ((shortened, InputRole.TEMPLATE), (spri, InputRole.CONDITION)), (purified_short, spri_waste), changes_molecule=False)
    graph.discard(spri_waste, "Removed by the stated 0.9X SPRI cleanup; molecular composition not assigned.")
    circle_ligase = graph.input("circle_ligase", "dilute-circularization T4 DNA ligase",
                                MaterialKind.CATALYST)
    circle_state = closed_circle("NlaIII-trimmed HiCAR circle", (
        role("CviQI-cut interacting genomic DNA"), role("HiCAR7–splint junction"),
        role("Tn5-tagged accessible genomic DNA")),
        closure_label="dilute intramolecular T4 ligation")
    circles = graph.product("circles", "NlaIII-trimmed HiCAR circles",
                            MaterialKind.MOLECULAR, state=circle_state)
    graph.react("circularize", "Dilute to 1 ng/µl and circularize with T4 ligase for at least 2 h",
                StepKind.CIRCULARIZATION,
                ((purified_short, InputRole.TEMPLATE), (circle_ligase, InputRole.CATALYST)),
                (circles,), geometry=CircularizationGeometry.known(
                    "double-strand intramolecular ligation",
                    statement="The HiCAR workflow dilute-ligates NlaIII-trimmed contact molecules into inverse-PCR circles."))
    circle_cleanup = graph.input("circle_cleanup", "Zymo DNA clean & concentrator", MaterialKind.CONDITION)
    purified_circles = graph.product("purified_circles", "purified circularization mixture", MaterialKind.MOLECULAR, state=circle_state)
    circle_waste = graph.product("circle_waste", "unrecovered material from circle purification", MaterialKind.WASTE)
    graph.react("purify_circles", "Purify circularization products before PmeI", StepKind.CLEANUP,
                ((circles, InputRole.TEMPLATE), (circle_cleanup, InputRole.CONDITION)), (purified_circles, circle_waste), changes_molecule=False)
    graph.discard(circle_waste, "Removed by the stated DNA clean & concentrator purification; molecular composition not assigned.")
    pmei = graph.input("pmei", "PmeI restriction endonuclease", MaterialKind.CATALYST)
    pmei_templates = graph.product(
        "pmei_treated_templates", "PmeI-treated HiCAR circle pool",
        MaterialKind.MOLECULAR,
        state=molecule("PmeI-treated HiCAR template pool", (
            role("circularized contact molecules and PmeI cleavage products"),)))
    graph.react("pmei_digest", "Digest circularized material with PmeI before PCR",
                StepKind.DIGESTION,
                ((purified_circles, InputRole.TEMPLATE), (pmei, InputRole.CATALYST)),
                (pmei_templates,))
    i501 = graph.input("i501_primer", "published NEB i501 inverse-PCR primer",
                       MaterialKind.OLIGO, source="Lu et al., 2023 Table 2; application to defining 2022 source unresolved",
                       state=Construct([Segment("NEB i501", I501_PCR)]))
    i7 = graph.input("i7_primer", "published Nextera i7-10-L inverse-PCR primer",
                     MaterialKind.OLIGO, source="Lu et al., 2023 Table 2; application to defining 2022 source unresolved",
                     state=Construct([Segment("Nextera i7-10-L", NEXTERA_I7_10_L)]))
    polymerase = graph.input("pcr_polymerase", "Q5 DNA polymerase",
                             MaterialKind.CATALYST)
    library = graph.product("sequencing_library", "linear dual-index HiCAR amplicon",
                            MaterialKind.LIBRARY, state=final)
    geometry = PCRGeometry.unavailable(final,
        reason="The actual 2022 PmeI-treated template and Table S1 oligos are unresolved. Matching the later protocol primers on the illustrative final amplicon does not establish their sites on that template.")
    graph.react("inverse_pcr", "INFERRED — 12-cycle Q5 PCR; later-source compatible final arms",
                StepKind.PCR,
                ((pmei_templates, InputRole.TEMPLATE), (i501, InputRole.STRUCTURAL),
                 (i7, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (library,), geometry=geometry)
    gel = graph.input("gel_selection", "300–750 bp gel extraction with Zymo DNA recovery", MaterialKind.SELECTION)
    selected = graph.product("selected_library", "300–750 bp recovered HiCAR library", MaterialKind.LIBRARY, state=final)
    excluded = graph.product("gel_excluded", "material outside the recovered gel interval", MaterialKind.WASTE)
    graph.react("gel_extract", "Gel-extract 300–750 bp PCR products", StepKind.SELECTION,
                ((library, InputRole.TEMPLATE), (gel, InputRole.SELECTION)), (selected, excluded), changes_molecule=False)
    graph.discard(excluded, "Excluded by the defining Step 5 gel extraction; composition not assigned.")
    graph.sequence(selected, SequencingPlan(
        final, tuple(primers), tuple(p.role for p in primers),
        {},
        {"Read 1": "include", "Index 1 (i7)": "demultiplex",
         "Index 2 (i5)": "demultiplex", "Read 2": "include"}))
    graph.declare_reviewed_scope((
        "fixed_nuclei", "hicar7_transposome", "tn5", "tagment", "cviqi",
        "cviqi_digest", "truseq_r1_splint", "proximity_ligase", "proximity_ligation",
        "crosslink_reversal", "reverse_crosslinks", "gap_fill_polymerase", "gap_fill",
        "nlaiii", "nlaiii_digest", "circle_ligase", "circularize", "pmei", "pmei_digest", "i501_primer",
        "i7_primer", "pcr_polymerase", "inverse_pcr", "sequencing_library",
        "tagment_quench", "stop_tagmentation", "polymerase_heat_stop", "inactivate_polymerase",
        "spri_condition", "purify_nlaiii", "circle_cleanup", "purify_circles",
        "gel_selection", "gel_extract", "selected_library",
    ), source)
    return graph.require_publishable()


def get():
    final, primers = hicar_library()
    spec = SimpleNamespace(
        TITLE="HiCAR — open-chromatin-anchored chromosome conformation",
        NOTES="01_hicar.html",
        SOURCE='<a href="https://doi.org/10.1016/j.molcel.2022.01.023">Wei et al., <i>Molecular Cell</i> (2022)</a>; defining CviQI/NlaIII implementation.',
        SUMMARY="An engineered Tn5 adapter marks accessible chromatin and is splint-ligated in fixed nuclei to a spatially proximal CviQI fragment. NlaIII digestion, dilute recircularization, PmeI treatment and inverse PCR then recover an asymmetric contact library.",
        CAVEAT="The later authors’ Current Protocols version substitutes MseI at the first digest. This page preserves the defining paper’s CviQI implementation; the displayed HiCAR7-i0 and PCR oligos are verbatim from the later 2023 Tables 1/2. <inf>Defining 2022 Table S1 is gated, so assignment of these exact sequences and final arms to 2022 remains inferred; the TruseqR1 splint sequence and acquired cycle counts are unresolved.</inf>",
        oligos=oligos,
        workflow=hicar_workflow,
        FINAL_LIBRARY=final,
        FINAL_CAPTION="<inf>INFERRED — compatible endpoint using later-source primer arms; exact 2022 oligo assignment remains unresolved.</inf> The hybrid library has a TruSeq entry on the interacting-partner side and a Nextera entry on the accessible-anchor side. The HiCAR7–splint junction remains between those genomic regions; bracketed genomic regions vary in length.",
        SEQ_PRIMERS=primers,
        READ_LENGTHS={},
        SEQUENCING_INTRO="Read 1 enters the CviQI-cut interacting genomic partner from the TruSeq side. Read 2 enters the Tn5-tagged accessible anchor from the Nextera mosaic end. The two index reads identify the library; there is no molecule-level UMI. Exact study acquisition cycles remain unresolved in this bounded review; index lengths in the inferred molecule do not establish acquired cycle counts.",
        SEQUENCING_ENDING="",
        SEQUENCING_UNAVAILABLE="",
    )
    spec.protocol_graph = lambda: _protocol_graph(final, primers)
    return spec


SPEC = get()
globals().update(vars(SPEC))


def protocol_graph():
    return _protocol_graph(SPEC.FINAL_LIBRARY, SPEC.SEQ_PRIMERS)
