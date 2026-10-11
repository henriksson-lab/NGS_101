"""Chromatin-conformation assays that preserve groups without proximity ligation.

These methods must not be represented as one covalent contact chimera.  The group is a
collection of independent molecules sharing either a physical compartment or a sampled
nuclear slice; constructors here keep that distinction explicit.
"""
from __future__ import annotations

from types import SimpleNamespace

import illumina as il
import seqprimers as sp
from batch_ngs import truseq_library
from chemdraw import Construct, Scene, Segment, feature
from protocol_model import chain, molecule, role


CHIA_COMPLEX = feature(
    "chia_drop.gem", "molecular_complex_barcode", "whitelist",
    whitelist="10x Chromium Genome v2 GemCode set",
    note="16-nt partition barcode shared by every copied fragment in one GEM")
CHIA_UMI = feature(
    "chia_drop.random_primer", "umi", "random",
    note="8-nt random priming sequence retained beside the GEM barcode")
CHIA_SAMPLE = feature(
    "chia_drop.sample", "sample_index", "unknown",
    note="library sample index; the defining paper does not publish the oligo set")


def _seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name, top, tag, **kw)


def grouped_fragments(label: str, *, copied: bool = False) -> Scene:
    """Three independent DNA molecules in one complex/GEM, never one joined strand."""
    sc = Scene()
    for n in range(1, 4):
        parts = []
        if copied:
            parts.extend((
                _seg(f"GEM barcode {n}", "B" * 16, "cbc", placeholder=True,
                     feature=CHIA_COMPLEX),
                _seg(f"random 8-mer / UMI {n}", "U" * 8, "umi", placeholder=True,
                     feature=CHIA_UMI),
            ))
        parts.append(_seg(f"genomic fragment {n}", "X" * (26 + n * 4),
                          placeholder=True, length_bp=1000 * n))
        sc.strand(f"fragment {n}", parts, label=f"fragment {n}")
        if copied:
            sc.junction(f"fragment {n}", f"GEM barcode {n}",
                        f"random 8-mer / UMI {n}", "barcode-primer synthesis")
            sc.junction(f"fragment {n}", f"random 8-mer / UMI {n}",
                        f"genomic fragment {n}", "random-primer extension")
    return sc


def chia_drop_library():
    """Representative Chromium-v2 product; proprietary arms are marked inferred."""
    lib = Construct([
        _seg("P5", il.P5, "p5", inferred=True),
        _seg("Read 1 arm", il.TRUSEQ_READ1, "r1", inferred=True),
        _seg("16-nt GEM complex barcode", "B" * 16, "cbc", placeholder=True,
             feature=CHIA_COMPLEX),
        _seg("8-nt random primer / UMI", "U" * 8, "umi", placeholder=True,
             feature=CHIA_UMI),
        _seg("copied genomic fragment", "X" * 46, placeholder=True),
        _seg("dA junction", "A", inferred=True),
        _seg("Index 1 / Read 2 arm", il.INDEX1_PRIMER, "r2", inferred=True),
        _seg("i7 sample index", "I" * 8, "cbc", placeholder=True,
             inferred=True, feature=CHIA_SAMPLE),
        _seg("P7 reverse complement", il.P7_RC, "p7", inferred=True),
    ], name="representative ChIA-Drop fragment library")
    primers = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["R2"])
    problems = sp.verify(lib, primers,
                         required_roles=("Read 1", "Index 1 (i7)", "Read 2"))
    if problems:
        raise ValueError("ChIA-Drop library: " + "; ".join(problems))
    return lib, primers


def chia_drop_workflow():
    initial = molecule("formaldehyde-crosslinked nucleus", (
        role("chromatin complex containing spatially associated DNA segments"),))
    return chain(initial, (
        ("Permeabilize nuclei and digest chromatin in situ",
         grouped_fragments("restriction-fragmented crosslinked complex"),
         "The defining implementation compared MboI, HindIII and sonication; longer HindIII/sonicated templates gave more mappable reads. Crosslinks preserve association, and there is no fill-in or proximity ligation."),
        ("Release intact crosslinked chromatin complexes by sonication",
         grouped_fragments("released multi-fragment chromatin complex"),
         "Sonication opens nuclei and brings chromatin to roughly 3–6-kb template lengths without de-crosslinking or purifying free DNA."),
        ("Partition about 0.5 ng into Chromium Genome v2 GEMs",
         grouped_fragments("one physical complex in one GEM"),
         "A gel bead supplies many copies of one 16-nt GemCode. A GEM is intended to contain one chromatin complex, but collisions can place more than one complex in a partition."),
        ("Random-prime overnight isothermal linear amplification at 30 degrees C",
         grouped_fragments("independent copies sharing one GEM barcode", copied=True),
         "Each accessible DNA template anneals to a bead oligo containing P5, a Read 1 site, the common 16-nt barcode and a random 8-mer. Products remain separate molecules: shared barcode, not ligation, records co-partitioning."),
        ("Break GEMs, pool products and complete the Illumina library",
         molecule("representative amplified fragment", (
             role("P5 and Read 1 arm", inferred=True),
             role("16-nt GEM complex barcode", length_bp=16, molecular_feature=CHIA_COMPLEX),
             role("8-nt random-primer UMI", length_bp=8, molecular_feature=CHIA_UMI),
             role("one genomic fragment copy"),
             role("Read 2, sample-index and P7 arm", inferred=True))),
         "The commercial v2 reagents establish the final flow-cell-compatible molecule, but the defining paper publishes adapter roles rather than the proprietary oligo sequences."),
    ))


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


def get(key: str):
    if key == "gam":
        final, primers = gam_library()
        return SimpleNamespace(
            TITLE="GAM — Genome Architecture Mapping by nuclear-profile co-segregation",
            NOTES="01_gam.html",
            SOURCE='<a href="https://doi.org/10.1038/nature21411">Beagrie et al., <i>Nature</i> (2017)</a>; original WGA4 and TruSeq implementation.',
            SUMMARY="A random 220-nm cryosection through one fixed nucleus is laser-microdissected, whole-genome amplified and sequenced as its own indexed library. Loci that repeatedly occur in the same thin profiles are inferred to have been close in three-dimensional space.",
            CAVEAT="GAM has no contact junction and no molecule-level complex barcode: association is statistical co-segregation across many separately indexed nuclear-profile libraries. <inf>The WGA4 adapter is known to occupy 30 bases at each genomic end, but its proprietary bases are not published.</inf>",
            oligos=lambda: (),
            workflow=gam_workflow,
            FINAL_LIBRARY=final,
            FINAL_CAPTION="A representative genomic piece from one nuclear profile, flanked by two proprietary 30-base WGA4 adapters inside a dual-index TruSeq library. The sample-index pair identifies the profile; it does not label contacts within a molecule. <inf>WGA4 adapter bases are undisclosed.</inf>",
            SEQ_PRIMERS=primers,
            READ_LENGTHS={"Read 1": 130, "Index 1 (i7)": 8,
                          "Index 2 (i5)": 8},
            SEQUENCING_INTRO="The original study used single-end sequencing: 30 initial Read 1 cycles traversing the known-length WGA4 adapter were run as dark cycles, followed by 100 imaged genomic cycles. Index reads demultiplex the nuclear-profile libraries; Read 2 was not used.",
            SEQUENCING_ENDING="",
            SEQUENCING_UNAVAILABLE="",
        )
    if key != "chia-drop":
        raise KeyError(key)
    final, primers = chia_drop_library()
    return SimpleNamespace(
        TITLE="ChIA-Drop — ligation-free droplet mapping of multiplex chromatin complexes",
        NOTES="01_chia_drop.html",
        SOURCE='<a href="https://doi.org/10.1038/s41586-019-0949-1">Zheng et al., <i>Nature</i> (2019)</a>; Chromium Genome v2 library and deposited read metadata.',
        SUMMARY="Crosslinked chromatin complexes are fragmented but not proximity-ligated. In a Chromium GEM, independent fragments from one physical complex undergo random-primed linear amplification with the same 16-nt barcode; reads are grouped by that shared partition identity.",
        CAVEAT="A GEM barcode identifies a droplet, not necessarily one biological complex: partition collisions are separated computationally where possible. <inf>The public final molecule uses the compatible standard Illumina arm geometry because the discontinued Chromium Genome v2 oligo sequences were not published in the defining paper.</inf>",
        oligos=lambda: (),
        workflow=chia_drop_workflow,
        FINAL_LIBRARY=final,
        FINAL_CAPTION="One representative copied fragment is shown. The 16-nt barcode is repeated across separate molecules from the same GEM; the random 8-mer both primes copying and supplies a per-copy UMI. <inf>Standard Illumina arm sequences show the compatible v2 geometry, not a claimed transcription of proprietary kit oligos.</inf>",
        SEQ_PRIMERS=primers,
        READ_LENGTHS={"Read 1": 150, "Index 1 (i7)": 8, "Read 2": 150},
        SEQUENCING_INTRO="The defining data were sequenced 2 × 150 bp. Read 1 begins with the 16-bp GemCode and 8-bp random-primer UMI before genomic sequence; Read 2 is genomic. The sample index is read separately. Primer landing positions on the compatible Illumina shell are computed below.",
        SEQUENCING_ENDING="",
        SEQUENCING_UNAVAILABLE="",
    )
