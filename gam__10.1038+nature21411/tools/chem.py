"""GAM preserves spatial association by sampling one nuclear profile.

The group is a collection of independent chromosome pieces selected by one physical
nuclear slice, never one covalent contact chimera.
"""
from __future__ import annotations

from types import SimpleNamespace

import illumina as il
import seqprimers as sp
from batch_ngs import truseq_library
from chemdraw import Construct, Scene, Segment, feature
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from protocol_model import chain, molecule, role


from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "lib"))

def _seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name, top, tag, **kw)


def nuclear_profile_scene(label: str, *, amplified: bool = False) -> Scene:
    """Independent chromosome pieces selected by one physical nuclear slice."""
    sc = Scene()
    for n in range(1, 4):
        parts = []
        if amplified:
            parts.append(_seg(f"5-prime WGA4 adapter {n}", "W" * 30,
                              placeholder=True, inferred=True, length_bp=30))
        parts.append(_seg(f"slice-intersected genomic region {n}", "X" * (24 + 5 * n),
                          placeholder=True, length_bp=30000))
        if amplified:
            parts.append(_seg(f"3-prime WGA4 adapter {n}", "W" * 30,
                              placeholder=True, inferred=True, length_bp=30))
        sc.strand(f"region {n}", parts, label=f"region {n}")
    return sc


def gam_library():
    lib, all_primers = truseq_library([
        _seg("5-prime WGA4 adapter", "W" * 30, placeholder=True,
             inferred=True, length_bp=30),
        role("one genomic region retained in the nuclear profile"),
        _seg("3-prime WGA4 adapter", "W" * 30, placeholder=True,
             inferred=True, length_bp=30),
    ], "GAM nuclear-profile library")
    primers = tuple(p for p in all_primers if p.role != "Read 2")
    problems = sp.verify(lib, primers, required_roles=tuple(p.role for p in primers))
    if problems:
        raise ValueError("GAM library: " + "; ".join(problems))
    return lib, primers


def gam_workflow():
    return chain(molecule("fixed, sucrose-embedded nucleus", (
        role("three-dimensional chromosome arrangement preserved in one nucleus"),)), (
        ("Cut an approximately 220-nm cryosection through a random nuclear plane",
         nuclear_profile_scene("one thin nuclear profile"),
         "A physical slice contains the genomic loci intersecting that plane. Spatially close loci co-occur across profiles more often, but the selected DNA regions are not covalently joined."),
        ("Laser-microdissect one nuclear profile into its own adhesive PCR cap",
         nuclear_profile_scene("one isolated nuclear profile"),
         "Each of the original GAM libraries represents one independently sampled nuclear profile. Empty caps and WGA-amplified bulk genomic DNA are carried as negative and positive controls."),
        ("Lyse, fragment and prepare the profile with the Sigma WGA4 system",
         nuclear_profile_scene("WGA4-prepared profile", amplified=True),
         "Protease extraction is followed by the proprietary library-preparation enzyme. The kit creates about 30 bases of WGA adapter at each amplified genomic end; its oligo sequences are not public."),
        ("Whole-genome amplify the nuclear-profile DNA",
         nuclear_profile_scene("amplified nuclear-profile regions", amplified=True),
         "WGA4 PCR supplies enough material from the subcellular slice for a conventional library. Presence or absence of genomic windows, not molecular counting, is the measurement."),
        ("Purify, size-select and build one indexed TruSeq library per profile",
         molecule("representative indexed GAM library molecule", (
             role("P5, i5 and Read 1 arm"), role("30-nt WGA4 adapter", length_bp=30, inferred=True),
             role("one profile-derived genomic region"), role("30-nt WGA4 adapter", length_bp=30, inferred=True),
             role("Read 2 arm, i7 and P7"))),
         "The defining study used TruSeq DNA HT or TruSeq Nano DNA HT, pooled 96 libraries and size-selected approximately 300–500-nt products."),
    ))


def _finish_pcr(graph, template, spec, prefix, reviewed, label):
    forward = graph.input(prefix + "_forward", label + " forward primer", MaterialKind.OLIGO)
    reverse = graph.input(prefix + "_reverse", label + " reverse primer", MaterialKind.OLIGO)
    polymerase = graph.input(prefix + "_polymerase", label + " polymerase",
                            MaterialKind.CATALYST)
    library = graph.product(prefix + "_library", label + " sequencing library",
                            MaterialKind.LIBRARY, state=spec.FINAL_LIBRARY)
    graph.react(prefix + "_pcr", label, StepKind.PCR,
                ((template, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (library,), geometry=PCRGeometry.unavailable(
                    spec.FINAL_LIBRARY,
                    reason="The commercial intermediate is not disclosed as an exact Construct; the final primer landing sites are checked independently."))
    graph.sequence(library, SequencingPlan(library=spec.FINAL_LIBRARY,
                                           primers=tuple(spec.SEQ_PRIMERS)))
    reviewed.extend((forward.id, reverse.id, polymerase.id, prefix + "_pcr",
                     library.id))


def _gam_graph(spec):
    source = spec.SOURCE
    graph = ProtocolGraph(spec.TITLE, authority=source)
    nucleus = graph.input("nucleus", "fixed sucrose-embedded nucleus",
                          MaterialKind.BIOLOGICAL,
                          state=molecule("fixed nucleus", (role("three-dimensional chromosome arrangement"),)),
                          primary=True)
    sectioning = graph.input("cryosection", "approximately 220-nm cryosection condition",
                             MaterialKind.CONDITION)
    profile = graph.product("profile", "one thin nuclear profile",
                            MaterialKind.MOLECULAR, state=nuclear_profile_scene("profile"))
    excluded = graph.product("excluded_nucleus", "nuclear material outside the sampled plane",
                             MaterialKind.WASTE)
    graph.react("section", "Cut a thin random plane through the nucleus",
                StepKind.SAMPLING,
                ((nucleus, InputRole.TEMPLATE), (sectioning, InputRole.CONDITION)),
                (profile, excluded))
    graph.discard(excluded, "lies outside the sampled nuclear profile")
    cap = graph.input("capture_cap", "adhesive laser-microdissection cap",
                      MaterialKind.SELECTION)
    isolated = graph.product("isolated", "one isolated nuclear profile",
                             MaterialKind.MOLECULAR, state=profile.state)
    rejected = graph.product("uncaptured", "uncaptured section material", MaterialKind.WASTE)
    graph.react("microdissect", "Laser-microdissect one nuclear profile",
                StepKind.SELECTION,
                ((profile, InputRole.TEMPLATE), (cap, InputRole.SELECTION)),
                (isolated, rejected))
    graph.discard(rejected, "is not part of the individually indexed nuclear profile")
    lysis = graph.input("wga4_lysis", "WGA4 proteinase K in single cell lysis and fragmentation buffer",
                        MaterialKind.CATALYST)
    lysed = graph.product("lysed", "deproteinized, fragmented nuclear-profile DNA",
                          MaterialKind.MOLECULAR, state=profile.state)
    graph.react("lyse", "Proteinase K lysis in the cap (50 °C, 4 h), then 99 °C for 4 min",
                StepKind.DEPROTEINIZATION,
                ((isolated, InputRole.TEMPLATE), (lysis, InputRole.CATALYST)), (lysed,))
    wga4 = graph.input("wga4_reagent", "WGA4 library preparation buffer, stabilisation solution and library preparation enzyme",
                       MaterialKind.CONDITION, inferred=True,
                       source=source + " (commercial composition not disclosed)")
    prepared = graph.product("prepared", "WGA4-prepared nuclear-profile DNA",
                             MaterialKind.MOLECULAR,
                             state=nuclear_profile_scene("WGA4-prepared profile", amplified=True))
    graph.react("prepare", "WGA4 library preparation: 95 °C 4 min, then enzyme at 16/24/37 °C (20 min each) and 75 °C 5 min",
                StepKind.MAGIC,
                ((lysed, InputRole.TEMPLATE), (wga4, InputRole.CONDITION)), (prepared,),
                reason="Methods give the reagent names, temperatures and times; the kit's adapter-attachment chemistry and adapter bases are proprietary, and only the 30-bp adaptor length at both ends is stated.")
    wga_f = graph.input("wga_forward", "WGA4 amplification primer", MaterialKind.OLIGO,
                        inferred=True, source=source + " (sequence proprietary)")
    wga_r = graph.input("wga_reverse", "WGA4 opposing amplification primer",
                        MaterialKind.OLIGO, inferred=True,
                        source=source + " (sequence proprietary)")
    wga_polymerase = graph.input("wga_polymerase", "WGA4 polymerase",
                                 MaterialKind.CATALYST)
    amplified = graph.product("amplified", "whole-genome-amplified profile DNA",
                              MaterialKind.MOLECULAR, state=prepared.state)
    graph.react("wga", "Whole-genome amplify the nuclear-profile DNA", StepKind.PCR,
                ((prepared, InputRole.TEMPLATE), (wga_f, InputRole.STRUCTURAL),
                 (wga_r, InputRole.STRUCTURAL),
                 (wga_polymerase, InputRole.CATALYST)), (amplified,),
                geometry=PCRGeometry.unavailable(
                    spec.FINAL_LIBRARY, reason="WGA4 adapter and primer sequences are proprietary."))
    reviewed = ["nucleus", "cryosection", "section", "capture_cap", "microdissect",
                "wga4_lysis", "lyse", "wga4_reagent", "prepare", "wga_forward", "wga_reverse",
                "wga_polymerase", "wga"]
    _finish_pcr(graph, amplified, spec, "library", reviewed,
                "Build the indexed TruSeq GAM library")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def get():
    final, primers = gam_library()
    spec = SimpleNamespace(
        TITLE="GAM — Genome Architecture Mapping by nuclear-profile co-segregation",
        NOTES="01_gam.html",
        SOURCE='<a href="https://doi.org/10.1038/nature21411">Beagrie et al., <i>Nature</i> (2017)</a>; original WGA4 and TruSeq implementation.',
        SUMMARY="A random 220-nm cryosection through one fixed nucleus is laser-microdissected, whole-genome amplified and sequenced as its own indexed library. Loci that repeatedly occur in the same thin profiles are inferred to have been close in three-dimensional space.",
        CAVEAT=("GAM has no contact junction and no molecule-level complex barcode: association is statistical co-segregation across many separately indexed nuclear-profile libraries. <inf>The WGA4 adapter is known to occupy 30 bases at each genomic end, but its proprietary bases are not published.</inf> "
                "Read length: Methods say each pool was sequenced &ldquo;single end 100 bp&rdquo; with the first 30 cycles not imaged (dark cycles); they do not say whether the 30 dark cycles are inside or in addition to the 100. The model uses the printed 100-cycle Read 1 and does not decide how many genomic bases were imaged. "
                "Published discrepancy (positive control only): the Isolation of nuclear profiles paragraph names WGA-amplified E14 mESC genomic DNA, while the Whole genome amplification paragraph adds 1 &micro;l human genomic DNA to an empty lid. The NP libraries modelled here are unaffected; the control's species is not resolved here."),
        oligos=lambda: (),
        workflow=gam_workflow,
        FINAL_LIBRARY=final,
        FINAL_CAPTION="A representative genomic piece from one nuclear profile, flanked by two proprietary 30-base WGA4 adapters inside a dual-index TruSeq library. The sample-index pair identifies the profile; it does not label contacts within a molecule. <inf>WGA4 adapter bases are undisclosed.</inf>",
        SEQ_PRIMERS=primers,
        READ_LENGTHS={"Read 1": 100, "Index 1 (i7)": 8,
                          "Index 2 (i5)": 8},
        SEQUENCING_INTRO="The original study sequenced pools of 96 libraries single-end 100 bp in HiSeq rapid-run mode; the first 30 Read 1 cycles, which traverse the known-length WGA4 adapter, were not imaged (dark cycles; whether they count within the 100 is not stated). Supplementary Table 1 lists 8-nt i7 and i5 indexes for every NP and negative-control library (TruSeq HT or Nano; acquired index cycles are not stated). Read 2 was not used.",
        SEQUENCING_ENDING="",
        SEQUENCING_UNAVAILABLE=""
    )
    spec.protocol_graph = lambda: _gam_graph(spec)
    return spec


SPEC = get()
globals().update(vars(SPEC))


def protocol_graph():
    return _gam_graph(SPEC)
