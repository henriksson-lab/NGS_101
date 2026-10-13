"""Reusable molecular states for isothermal amplification chemistries.

These helpers describe topology, not a fictitious single consensus sequence.  LAMP and
multiple-displacement products are populations of branched molecules, so callers supply
the sequence-bearing primers separately and use these states for the reaction workflow.
"""
from __future__ import annotations

from chemdraw import Scene
from protocol_model import molecule, role


def lamp_stem_loop(target: str = "target-derived DNA") -> Scene:
    """First self-priming dumbbell/stem-loop intermediate produced by inner primers."""
    return molecule("LAMP stem-loop intermediate", (
        role("FIP-derived inverted repeat"), role(target),
        role("BIP-derived inverted repeat"), role("self-paired terminal loops")))


def lamp_product(target: str = "target-derived DNA",
                 barcode: str = "FIP sample barcode") -> Scene:
    """Branched concatemer population after loop-primer-assisted displacement."""
    return molecule("branched LAMP concatemer population", (
        role("repeated inverted target units"), role(barcode, tag="cbc"),
        role("stem loops and displaced branches"), role(target)))


def mda_product(label: str = "branched multiple-displacement product") -> Scene:
    """Random-primer/phi29 product without implying one linear final molecule."""
    return molecule(label, (
        role("long copied genomic tracts"), role("displaced single strands"),
        role("new random-primer extension branches")))


def spia_product() -> Scene:
    """Ribo-SPIA cDNA population made through repeated RNase-H priming cycles."""
    return molecule("Ribo-SPIA amplified cDNA population", (
        role("repeatedly displaced transcript-derived cDNA copies"),
        role("complementary DNA products"), role("SPIA handle remnant")))


def proximity_extension_pair() -> Scene:
    """Two antibody probes co-bound to one analyte before polymerase extension."""
    return molecule("antigen-bound PEA probe pair", (
        role("forward antibody + oligo"), role("target protein"),
        role("reverse antibody + oligo"), role("paired hybridization arms")))


def proximity_extension_product() -> Scene:
    """Reporter created when one probe arm templates extension of its neighbour."""
    return molecule("proximity-extension reporter", (
        role("P5 and Read 1 precursor"), role("forward assay barcode", tag="cbc"),
        role("copied hybridization junction"), role("reverse assay barcode", tag="cbc"),
        role("common reverse-PCR handle")))
