"""Sequencing readouts of distance-dependent covalent labelling in fixed cells."""
from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry
from batch_ngs import seg, truseq_library
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from protocol_model import chain, da_tailed, molecule, role


def tsa_library():
    """Compatible Illumina shell; the study does not identify its indexed adapters."""
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
         "KAPA HyperPrep performs end repair and A-tailing in one reaction; the dA-tailed fragment is the substrate for indexed-adapter ligation."),
        ("Ligate six-base-indexed adapters and PCR for 8–16 cycles",
         molecule("single-index TSA-seq library", (
             role("P5 and Read 1 arm", inferred=True),
             role("biotin-enriched genomic insert"),
             role("Read 2 arm, six-base sample index and P7", inferred=True))),
         "INFERRED — The defining paper names KAPA HyperPrep chemistry but gives no adapter supplier or oligo bases. The HyperPrep kit did not include adapters."),
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
    peroxide = graph.input(
        "peroxidase_quench", "1.5% H2O2 in CMF-PBS, 1 h RT; endogenous-peroxidase quench",
        MaterialKind.CONDITION)
    quenched = graph.product(
        "quenched_cells", "fixed, permeabilized, peroxidase-quenched cells",
        MaterialKind.MOLECULAR, state=cells.state)
    graph.react(
        "quench_peroxidase", "Quench endogenous peroxidases before immunostaining",
        StepKind.QUENCHING,
        ((cells, InputRole.TEMPLATE), (peroxide, InputRole.CONDITION)), (quenched,))
    antibody_hrp = graph.input(
        "antibody_hrp", "IgG block, primary antibody 20–24 h 4 °C, HRP-conjugated secondary 20–24 h 4 °C",
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
        ((quenched, InputRole.TEMPLATE), (antibody_hrp, InputRole.SELECTION)),
        (localized, unstained))
    graph.discard(unstained, "does not carry the localized HRP source used for TSA labelling")

    tyramide = graph.input("tyramide_biotin", "tyramide-biotin (8 mM stock diluted 1:5,000) with 0.003% H2O2, 10 min",
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

    reversal = graph.input("crosslink_reversal", "SDS/proteinase K lysis 55 °C 8–12 h; 0.2 M NaCl 65 °C 20–24 h; phenol–chloroform, RNase A, ethanol precipitation",
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
    streptavidin = graph.input("streptavidin_beads", "Dynabeads M-270 streptavidin in W&B buffer; stringent washes; proteinase K digestion with free biotin, phenol–chloroform, glycogen-carrier precipitation",
                              MaterialKind.SELECTION)
    captured = graph.product(
        "captured_dna", "streptavidin-enriched TSA-labelled DNA",
        MaterialKind.MOLECULAR,
        state=molecule("TSA-seq pulldown DNA", (
            role("random sonication end"), role("biotin-enriched genomic fragment"),
            role("random sonication end"))))
    unbound = graph.product("unbound_dna", "unbound DNA, predominantly unlabelled (paper estimates ~0.002–0.15% of fragments biotinylated)", MaterialKind.WASTE)
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
        "fixed_cells", "peroxidase_quench", "quench_peroxidase", "antibody_hrp", "hrp_catalyst", "immunostain", "tyramide_biotin",
        "tyramide_labelling", "crosslink_reversal", "reverse_crosslinks",
        "sonication", "sonicate", "split_input", "streptavidin_beads",
        "streptavidin_capture",
    ]

    def prepare(prefix: str, template, label: str):
        erat = graph.input(prefix + "_erat_mix", label + " KAPA HyperPrep end repair and A-tailing enzyme mix",
                           MaterialKind.CATALYST)
        tailed = graph.product(prefix + "_tailed", label + " blunt-ended, 3-prime dA-tailed DNA",
                               MaterialKind.MOLECULAR,
                               state=da_tailed(label + " dA-tailed DNA", label + " genomic insert"))
        graph.react(prefix + "_erat", "Blunt-end and 3-prime dA-tail " + label + " DNA in one ER & AT reaction",
                    StepKind.TAILING,
                    ((template, InputRole.TEMPLATE), (erat, InputRole.CATALYST)),
                    (tailed,))
        adapters = graph.input(prefix + "_adapters", label + " six-base-indexed adapters of unidentified supplier",
                               MaterialKind.ADAPTER, inferred=True)
        ligase = graph.input(prefix + "_ligase", label + " adapter ligase",
                             MaterialKind.CATALYST)
        adapted = graph.product(prefix + "_adapted", label + " adapter-ligated DNA",
                                MaterialKind.MOLECULAR, state=final, inferred=True)
        graph.react(prefix + "_ligation", "Ligate six-base-indexed adapters",
                    StepKind.LIGATION,
                    ((tailed, InputRole.TEMPLATE), (adapters, InputRole.STRUCTURAL),
                     (ligase, InputRole.CATALYST)), (adapted,),
                    geometry=AdapterLigationGeometry.unavailable(
                        "KAPA HyperPrep supplied the ligation chemistry but not adapters; the kit sheet describes 3'-dTMP-overhang adapters generically, but the paper does not identify the adapter supplier, sequences or exact end geometry."))
        forward = graph.input(prefix + "_pcr_forward", label + " library PCR forward primer",
                              MaterialKind.OLIGO, inferred=True)
        reverse = graph.input(prefix + "_pcr_reverse", label + " library PCR reverse primer",
                              MaterialKind.OLIGO, inferred=True)
        polymerase = graph.input(prefix + "_pcr_polymerase", label + " KAPA HiFi polymerase",
                                 MaterialKind.CATALYST)
        library = graph.product(prefix + "_library", label + " TSA-seq library",
                                MaterialKind.LIBRARY, state=final, inferred=True)
        graph.react(prefix + "_pcr", "PCR-amplify the " + label + " library for 8–16 cycles",
                    StepKind.PCR,
                    ((adapted, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                     (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                    (library,), geometry=PCRGeometry.unavailable(
                        final, reason="The paper names KAPA HiFi and 8–16 cycles but not the amplification primers; their sequences and binding sites are unassigned in this model."))
        graph.sequence(library, plan)
        reviewed.extend((
            prefix + "_erat_mix", prefix + "_erat", prefix + "_adapters", prefix + "_ligase",
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
        CAVEAT="TSA-seq does not encode distance in a molecule-level barcode or junction: distance is inferred from population enrichment. <inf>The final drawing uses compatible standard Illumina geometry for six-base-indexed adapters whose supplier and sequences the paper does not identify.</inf>",
        oligos=lambda: (),
        workflow=tsa_workflow,
        FINAL_LIBRARY=final,
        FINAL_CAPTION="One enriched genomic fragment inside a compatible single-index Illumina shell. Biotin enabled pulldown before library construction and is not a sequenced base. <inf>Adapter bases show compatible geometry; the actual ordered adapter oligos remain unidentified.</inf>",
        SEQ_PRIMERS=primers,
        READ_LENGTHS={"Read 1": 100, "Index 1 (i7)": 6},
        SEQUENCING_INTRO="The paper reports 101 sequencing cycles on HiSeq 2000/2500/4000 and 100-nt reads per sample. INFERRED — a single-end genomic read and a 6-cycle index read are model choices from those counts and the 6-nt adaptor barcodes; the paper does not state the read mode or index-read cycles. Read density in streptavidin pulldown versus input supplies the distance-dependent signal.",
        SEQUENCING_ENDING="",
        SEQUENCING_UNAVAILABLE="",
    )
    spec.protocol_graph = lambda: tsa_protocol_graph(final, primers)
    return spec


SPEC = get()
globals().update(vars(SPEC))


def protocol_graph():
    return tsa_protocol_graph(SPEC.FINAL_LIBRARY, SPEC.SEQ_PRIMERS)
