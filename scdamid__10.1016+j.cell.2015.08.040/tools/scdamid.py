"""Molecular model for single-cell DamID (Kind et al., Cell 2015)."""
from __future__ import annotations

import illumina as il
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, complement_segments, revcomp

# The accessible defining paper names these oligos but does not print their sequences.
# Exact bases are from the upstream scg_lib_structs reconstruction and are therefore
# created only through secondary(), which makes inferred styling unavoidable.
ADRT = "CTAATACGACTCACTATAGGGCAGCGTGGTCGCGGCCGAGGA"
ADRB = "TCCTCGGCCG"
PREAMP_STUB = "GTGGTCGCGGCCGAGGA"
RANDOM_NT = 4
INDEX_NT = 6
RUN_ROLES = ("Read 1", "Index 1 (i7)")


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


def secondary(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    """A source-bounded segment that can never render as established sequence."""
    return seg(name, top, tag, inferred=True, **kw)


def damid_adaptor_scene() -> Scene:
    paired = ADRT[-len(ADRB):]
    sc = Scene()
    sc.strand("AdRt", [secondary("5' tail", ADRT[:-len(ADRB)]),
                       secondary("paired end", paired)], label="AdRt")
    sc.anneal("AdRb", [secondary("paired end reverse complement", ADRB)],
              to="AdRt", pair=("paired end reverse complement", "paired end"),
              label="AdRb")
    return sc


def dpni_fragment(insert_nt: int = 28) -> Construct:
    """Blunt genomic fragment selected between two Dam-methylated GATC sites."""
    return Construct([
        seg("left GATC", "GATC"),
        seg("Dam-marked genomic DNA", "X" * insert_nt, placeholder=True),
        seg("right GATC", "GATC"),
    ], name="DpnI-selected genomic fragment")


def preamp_product(insert_nt: int = 28) -> Construct:
    """One strand after the single DamID pre-amplification primer acts at both ends."""
    return Construct([
        secondary("N4, left", "N" * RANDOM_NT, placeholder=True),
        secondary("DamID stub, left", PREAMP_STUB),
        secondary("restored GATC, left", "TC"),
        seg("Dam-marked genomic DNA", "X" * insert_nt, placeholder=True),
        secondary("restored GATC, right", "GA"),
        secondary("DamID stub, right", revcomp(PREAMP_STUB)),
        secondary("N4, right", "N" * RANDOM_NT, placeholder=True),
    ], name="scDamID pre-amplification product")


SEQ_PRIMERS = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"])


def final_library(insert_nt: int = 28) -> Construct:
    """Single-index Illumina library; every outer base is secondary-source here."""
    lib = Construct([
        secondary("P5", il.P5, "p5"),
        secondary("Read 1 remainder", il.TRUSEQ_READ1[4:], "r1"),
        *preamp_product(insert_nt),
        secondary("Read 2 reverse complement", revcomp(il.TRUSEQ_READ2), "r2"),
        secondary("i7", "I" * INDEX_NT, "cbc", placeholder=True),
        secondary("P7 reverse complement", il.P7_RC, "p7"),
    ], name="scDamID sequencing library")
    problems = sp.verify(lib, SEQ_PRIMERS, required_roles=RUN_ROLES)
    if problems:
        raise ValueError("invalid scDamID final library: " + "; ".join(problems))
    return lib


def preamp_duplex() -> Scene:
    return Scene.duplex(preamp_product().segments)


def _validate() -> None:
    # This is construction logic: the published adaptor can only be exposed if its short
    # strand actually pairs, and the final run must retain both used primer sites.
    damid_adaptor_scene().rows()
    Scene.duplex(dpni_fragment().segments).rows()
    preamp_duplex().rows()
    final_library()


_validate()
