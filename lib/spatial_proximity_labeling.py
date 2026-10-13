"""Sequencing readouts of distance-dependent covalent labelling in fixed cells."""
from __future__ import annotations

from types import SimpleNamespace

import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry
from batch_ngs import seg, truseq_library
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from protocol_model import chain, da_tailed, molecule, role


def tsa_library():
    """Compatible Kapa Hyper Prep shell; undisclosed kit oligos stay inferred."""
    lib, all_primers = truseq_library([
        seg("TSA-biotin-labelled genomic fragment", "X" * 46, placeholder=True),
    ], "TSA-seq pulldown library", dual_index=False,
       inferred_adapters=True, i7_nt=6)
    primers = tuple(p for p in all_primers if p.role != "Read 2")
    problems = sp.verify(lib, primers, required_roles=tuple(p.role for p in primers))
    if problems:
        raise ValueError("TSA-seq library: " + "; ".join(problems))
    return lib, primers


def tsa_workflow():
    return chain(molecule("formaldehyde-fixed, permeabilized cells", (
        role("target nuclear compartment"), role("chromosomes at varying distances"))), (
        ("Immunostain the compartment marker and recruit horseradish peroxidase",
         molecule("localized HRP source in an intact nucleus", (
             role("antibody-marked nuclear compartment"), role("HRP at the target"),
             role("surrounding chromosomal DNA"))),
         "The defining demonstration used SON for nuclear speckles and lamin antibodies for the nuclear lamina; a no-primary sample controls nonspecific labelling."),
        ("Generate short-lived tyramide-biotin radicals around the HRP source",
         molecule("distance-decaying DNA labelling field", (
             role("HRP source"), role("strongly biotin-labelled nearby DNA"),
             role("weakly labelled distant DNA"))),
         "Radical concentration and therefore direct DNA biotinylation decay with cytological distance. This is a spatial concentration gradient, not a DNA contact or ligation junction."),
        ("Reverse crosslinks, purify DNA and sonicate to 200–800 bp",
         molecule("sheared labelled genomic DNA", (
             role("random sonication end"), role("genomic DNA carrying sparse biotin"),
             role("random sonication end"))),
         "Most fragments were about 500 bp. Ten percent input was retained before affinity purification."),
        ("Capture biotinylated fragments on streptavidin beads and release DNA",
         molecule("TSA-seq pulldown DNA", (
             role("random sonication end"), role("biotin-enriched genomic fragment"),
             role("random sonication end"))),
         "Stringent washes enrich fragments in proportion to their direct tyramide-biotin labelling; pulldown coverage is compared with the matched input."),
        ("Blunt-end and add one 3-prime dA to each end",
         da_tailed("dA-tailed TSA-seq pulldown DNA", "biotin-enriched genomic insert"),
         "The dA-tailed pulldown fragment is the substrate for Kapa adapter ligation."),
        ("Ligate indexed Kapa adapters and PCR for 8–16 cycles",
         molecule("single-index TSA-seq library", (
             role("P5 and Read 1 arm", inferred=True),
             role("biotin-enriched genomic insert"),
             role("Read 2 arm, six-base sample index and P7", inferred=True))),
         "The defining paper names Kapa Hyper Prep and six-base indexed adapters but does not publish the kit oligo sequences."),
    ))


def tsa_protocol_graph(final, primers) -> ProtocolGraph:
    """Both sequenced TSA pulldown and matched-input branches, with explicit losses."""
    source = "Chen et al., Journal of Cell Biology (2018); doi:10.1083/jcb.201807108"
    graph = ProtocolGraph("TSA-seq", authority=source)
    cells = graph.input(
        "fixed_cells", "formaldehyde-fixed permeabilized cells",
        MaterialKind.BIOLOGICAL,
        state=molecule("formaldehyde-fixed, permeabilized cells", (
            role("target nuclear compartment"), role("chromosomes at varying distances"))),
        primary=True)
    antibody_hrp = graph.input(
        "antibody_hrp", "compartment-specific primary antibody detection complex",
        MaterialKind.SELECTION)
    hrp = graph.input("hrp_catalyst", "antibody-localized horseradish peroxidase",
                      MaterialKind.CATALYST)
    localized = graph.product(
        "localized_hrp", "cells with HRP localized at the target compartment",
        MaterialKind.MOLECULAR,
        state=molecule("localized HRP source in an intact nucleus", (
            role("antibody-marked nuclear compartment"), role("HRP at the target"),
            role("surrounding chromosomal DNA"))))
    unstained = graph.product("unstained_cells", "unbound or unsuccessfully stained cells",
                              MaterialKind.WASTE)
    graph.react(
        "immunostain", "Immunostain the compartment marker and recruit HRP",
        StepKind.BINDING,
        ((cells, InputRole.TEMPLATE), (antibody_hrp, InputRole.SELECTION)),
        (localized, unstained))
    graph.discard(unstained, "does not carry the localized HRP source used for TSA labelling")

    tyramide = graph.input("tyramide_biotin", "tyramide-biotin labelling reagent",
                           MaterialKind.MOLECULAR)
    labelled = graph.product(
        "labelled_chromatin", "distance-dependently biotin-labelled chromatin",
        MaterialKind.MOLECULAR,
        state=molecule("distance-decaying DNA labelling field", (
            role("HRP source"), role("strongly biotin-labelled nearby DNA"),
            role("weakly labelled distant DNA"))))
    graph.react(
        "tyramide_labelling", "Generate short-lived tyramide-biotin radicals",
        StepKind.CHEMICAL_CONVERSION,
        ((localized, InputRole.TEMPLATE), (tyramide, InputRole.STRUCTURAL),
         (hrp, InputRole.CATALYST)), (labelled,))

    reversal = graph.input("crosslink_reversal", "cross-link reversal and DNA-purification condition",
                           MaterialKind.CONDITION)
    purified = graph.product("purified_dna", "purified biotin-labelled genomic DNA",
                             MaterialKind.MOLECULAR, state=labelled.state)
    graph.react(
        "reverse_crosslinks", "Reverse cross-links and purify labelled DNA",
        StepKind.DEPROTEINIZATION,
        ((labelled, InputRole.TEMPLATE), (reversal, InputRole.CONDITION)), (purified,))
    sonication = graph.input("sonication", "200–800-bp sonication condition",
                             MaterialKind.CONDITION)
    sheared = graph.product(
        "sheared_dna", "sheared labelled genomic DNA", MaterialKind.MOLECULAR,
        state=molecule("sheared labelled genomic DNA", (
            role("random sonication end"), role("genomic DNA carrying sparse biotin"),
            role("random sonication end"))))
    graph.react(
        "sonicate", "Sonicate purified DNA to 200–800 bp", StepKind.FRAGMENTATION,
        ((purified, InputRole.TEMPLATE), (sonication, InputRole.CONDITION)), (sheared,))

    input_aliquot = graph.product("input_aliquot", "matched ten-percent input DNA",
                                  MaterialKind.MOLECULAR, state=sheared.state)
    capture_aliquot = graph.product("capture_aliquot", "DNA assigned to streptavidin capture",
                                    MaterialKind.MOLECULAR, state=sheared.state)
    graph.react("split_input", "Retain matched input before affinity purification",
                StepKind.SPLIT, ((sheared, InputRole.TEMPLATE),),
                (input_aliquot, capture_aliquot))
    streptavidin = graph.input("streptavidin_beads", "streptavidin affinity beads",
                              MaterialKind.SELECTION)
    captured = graph.product(
        "captured_dna", "streptavidin-enriched TSA-labelled DNA",
        MaterialKind.MOLECULAR,
        state=molecule("TSA-seq pulldown DNA", (
            role("random sonication end"), role("biotin-enriched genomic fragment"),
            role("random sonication end"))))
    unbound = graph.product("unbound_dna", "unbound weakly labelled DNA", MaterialKind.WASTE)
    graph.react(
        "streptavidin_capture", "Capture biotinylated fragments and release DNA",
        StepKind.AFFINITY_CAPTURE,
        ((capture_aliquot, InputRole.TEMPLATE), (streptavidin, InputRole.SELECTION)),
        (captured, unbound))
    graph.discard(unbound, "is not retained by the TSA-seq affinity pulldown")

    plan = SequencingPlan(
        final, tuple(primers), tuple(p.role for p in primers),
        {"Read 1": 100, "Index 1 (i7)": 6})
    reviewed = [
        "fixed_cells", "antibody_hrp", "hrp_catalyst", "immunostain", "tyramide_biotin",
        "tyramide_labelling", "crosslink_reversal", "reverse_crosslinks",
        "sonication", "sonicate", "split_input", "streptavidin_beads",
        "streptavidin_capture",
    ]

    def prepare(prefix: str, template, label: str):
        repair = graph.input(prefix + "_repair_mix", label + " blunt-end repair mix",
                             MaterialKind.CATALYST)
        repaired = graph.product(prefix + "_repaired", label + " blunt-ended DNA",
                                 MaterialKind.MOLECULAR, state=template.state)
        graph.react(prefix + "_repair", "Blunt-end " + label + " DNA",
                    StepKind.END_REPAIR,
                    ((template, InputRole.TEMPLATE), (repair, InputRole.CATALYST)),
                    (repaired,))
        tailing = graph.input(prefix + "_tailing_enzyme", label + " dA-tailing enzyme",
                              MaterialKind.CATALYST)
        tailed = graph.product(prefix + "_tailed", label + " dA-tailed DNA",
                               MaterialKind.MOLECULAR,
                               state=da_tailed(label + " dA-tailed DNA", label + " genomic insert"))
        graph.react(prefix + "_tail", "Add one 3-prime dA to each end",
                    StepKind.TAILING,
                    ((repaired, InputRole.TEMPLATE), (tailing, InputRole.CATALYST)),
                    (tailed,))
        adapters = graph.input(prefix + "_adapters", label + " indexed Kapa adapters",
                               MaterialKind.ADAPTER, inferred=True)
        ligase = graph.input(prefix + "_ligase", label + " adapter ligase",
                             MaterialKind.CATALYST)
        adapted = graph.product(prefix + "_adapted", label + " adapter-ligated DNA",
                                MaterialKind.MOLECULAR, state=final, inferred=True)
        graph.react(prefix + "_ligation", "Ligate indexed Kapa adapters",
                    StepKind.LIGATION,
                    ((tailed, InputRole.TEMPLATE), (adapters, InputRole.STRUCTURAL),
                     (ligase, InputRole.CATALYST)), (adapted,),
                    geometry=AdapterLigationGeometry.unavailable(
                        "Kapa Hyper Prep is named, but its six-base indexed adapter oligos are proprietary."))
        forward = graph.input(prefix + "_pcr_forward", label + " library PCR forward primer",
                              MaterialKind.OLIGO, inferred=True)
        reverse = graph.input(prefix + "_pcr_reverse", label + " indexed library PCR reverse primer",
                              MaterialKind.OLIGO, inferred=True)
        polymerase = graph.input(prefix + "_pcr_polymerase", label + " library PCR polymerase",
                                 MaterialKind.CATALYST)
        library = graph.product(prefix + "_library", label + " TSA-seq library",
                                MaterialKind.LIBRARY, state=final, inferred=True)
        graph.react(prefix + "_pcr", "PCR-amplify the " + label + " library",
                    StepKind.PCR,
                    ((adapted, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                     (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                    (library,), geometry=PCRGeometry.unavailable(
                        final, reason="The compatible Illumina endpoint is known, but the proprietary Kapa PCR oligos are not published."))
        graph.sequence(library, plan)
        reviewed.extend((
            prefix + "_repair_mix", prefix + "_repair", prefix + "_tailing_enzyme",
            prefix + "_tail", prefix + "_adapters", prefix + "_ligase",
            prefix + "_ligation", prefix + "_pcr_forward", prefix + "_pcr_reverse",
            prefix + "_pcr_polymerase", prefix + "_pcr", prefix + "_library",
        ))

    prepare("pulldown", captured, "pulldown")
    prepare("input", input_aliquot, "matched input")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def get():
    final, primers = tsa_library()
    spec = SimpleNamespace(
        TITLE="TSA-seq — distance mapping from a nuclear compartment",
        NOTES="01_tsa_seq.html",
        SOURCE='<a href="https://doi.org/10.1083/jcb.201807108">Chen et al., <i>Journal of Cell Biology</i> (2018)</a>.',
        SUMMARY="Antibody-localized HRP generates a tyramide-biotin radical gradient in fixed nuclei. Nearby DNA is labelled more strongly than distant DNA, then sonicated, captured with streptavidin and sequenced against matched input.",
        CAVEAT="TSA-seq does not encode distance in a molecule-level barcode or junction: distance is inferred from population enrichment. <inf>The final drawing uses compatible standard Illumina geometry for the Kapa Hyper Prep adapters because their exact six-base indexed oligos were not published.</inf>",
        oligos=lambda: (),
        workflow=tsa_workflow,
        FINAL_LIBRARY=final,
        FINAL_CAPTION="One enriched genomic fragment inside the compatible single-index Kapa/Illumina shell. Biotin enabled pulldown before library construction and is not a sequenced base. <inf>Adapter bases show compatible geometry rather than a source transcription of the proprietary kit oligos.</inf>",
        SEQ_PRIMERS=primers,
        READ_LENGTHS={"Read 1": 100, "Index 1 (i7)": 6},
        SEQUENCING_INTRO="The defining libraries were sequenced for 101 instrument cycles and yielded 100-nt single-end genomic reads plus a six-base sample index. Read density in streptavidin pulldown versus input supplies the distance-dependent signal.",
        SEQUENCING_ENDING="",
        SEQUENCING_UNAVAILABLE="",
    )
    spec.protocol_graph = lambda: tsa_protocol_graph(final, primers)
    return spec
