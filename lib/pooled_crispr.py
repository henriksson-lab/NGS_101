"""Reusable construct models for Zhang-lineage pooled CRISPR screens.

The protocol pages choose a vector and a readout independently.  This module keeps that
separation structural: vector pages cannot silently turn a one-PCR readout into a nested
readout, and both readout builders return finished libraries on which ``seqprimers`` can
locate the primers actually used by the run.
"""
from __future__ import annotations

import crispr
import illumina as il
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Scene, Segment, feature, revcomp
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from protocol_model import molecule, role


def seg(name: str, seq: str, tag: str | None = None, *, placeholder: bool = False,
        **kw) -> Segment:
    return Segment(name, seq, tag, placeholder=placeholder, **kw)


SPACER = "N" * crispr.SPACER_LEN


def guide_oligos() -> tuple[str, str]:
    """BsmBI-compatible oligos for a representative non-G-starting spacer."""
    return crispr.guide_oligos(SPACER)


def cloned_guide_cassette() -> Construct:
    """The sequence restored by guide ligation, with the two ligation seams named."""
    return Construct([
        seg("hU6 promoter", crispr.U6_UPSTREAM + crispr.U6_3PRIME, "s5"),
        seg("Pol III +1 G", crispr.U6_PLUS1, "me"),
        seg("20-nt guide", SPACER, "cbc", placeholder=True,
            feature=feature("guide", "guide_barcode", "unknown")),
        seg("sgRNA scaffold", crispr.SCAFFOLD_V1, "r2"),
    ], name="cloned guide cassette")


def cloning_scene() -> Scene:
    """Duplex cassette with both BsmBI-directed ligation junctions marked."""
    con = cloned_guide_cassette()
    sc = Scene.duplex(list(con), label="cloned cassette")
    sc.junction("top", "hU6 promoter", "Pol III +1 G", "ligation junction 1")
    sc.junction("top", "20-nt guide", "sgRNA scaffold", "ligation junction 2")
    return sc


def vector_map(version: str) -> tuple[tuple[str, str], ...]:
    """Ordered transfer-genome features; cPPT placement is the v1/v2 distinction."""
    common_left = (("5' LTR", "r3"), ("packaging / RRE", "s5"))
    guide = (("hU6-guide-scaffold", "cbc"),)
    expression = (("EFS-Cas9-P2A-Puro", "me"), ("WPRE", "t7"), ("3' LTR", "r2"))
    if version == "v1":
        return common_left + guide + (("cPPT/CTS", "w1"),) + expression
    if version == "v2":
        return common_left + (("cPPT/CTS", "w1"),) + guide + expression
    raise ValueError("version must be 'v1' or 'v2'")


# Joung et al. 2017: a single PCR, with a reverse primer in the invariant scaffold.
JOUNG_STAGGER = "TAAGTAGAG"
JOUNG_U6_SITE = crispr.U6_UPSTREAM[1:] + crispr.U6_3PRIME
JOUNG_SCAFFOLD_SITE_TOP = crispr.SCAFFOLD_V1[50:73]
JOUNG_I7 = "TCGCCTTG"


def one_pcr_library() -> Construct:
    """Finished ~260-bp Joung one-PCR library for a non-G-starting guide."""
    con = Construct([
        seg("P5", il.P5, "p5"),
        seg("Read 1 site", il.TRUSEQ_P5_FULL[len(il.P5):], "s5"),
        seg("stagger", JOUNG_STAGGER, "me"),
        seg("hU6 primer site", JOUNG_U6_SITE, "s5"),
        seg("Pol III +1 G", crispr.U6_PLUS1, "me"),
        seg("20-nt guide", SPACER, "cbc", placeholder=True,
            feature=feature("guide", "guide_barcode", "unknown")),
        seg("scaffold before reverse site", crispr.SCAFFOLD_V1[:50], "r2"),
        seg("reverse-primer site", JOUNG_SCAFFOLD_SITE_TOP, "r2"),
        seg("Read 2 site", revcomp(il.TRUSEQ_READ2), "t7"),
        seg("i7", revcomp(JOUNG_I7), "cbc",
            feature=feature("sample_index_i7", "sample_index", "fixed")),
        seg("P7", revcomp(il.P7), "p7"),
    ], name="Joung one-PCR library")
    problems = sp.verify(con, one_pcr_sequencing_primers(),
                         required_roles=("Read 1", "Index 1 (i7)"))
    if problems:
        raise ValueError("invalid one-PCR library: " + "; ".join(problems))
    return con


def one_pcr_sequencing_primers():
    return [sp.TRUSEQ["R1"], sp.TRUSEQ["I1"]]


# Moffat Lab TKOv3 two-step readout (2019 protocol, LCV2::TKOv3 branch).
MOFFAT_PCR1_F = "GAGGGCCTATTTCCCATGATTC"
MOFFAT_PCR1_R = "GTTGCGAAAAAGAACGTTCACGG"
MOFFAT_PCR2_SITE_F = "TTGTGGAAAGGACGAAACACCG"
MOFFAT_PCR2_SITE_R = "ACTTGCTATTTCTAGCTCTAAAAC"
MOFFAT_I5 = il.NEBNEXT_I5_SET1["i501"]
MOFFAT_I7 = il.NEBNEXT_I7_SET1["i701"]


def two_step_intermediate() -> Construct:
    """PCR1 enrichment product; distal genomic/vector sequence is intentionally schematic."""
    return Construct([
        seg("PCR1 forward site", MOFFAT_PCR1_F, "s5"),
        seg("left vector / genomic flank", "X" * 12, "w1", placeholder=True),
        seg("hU6", crispr.U6_3PRIME + crispr.U6_PLUS1, "s5"),
        seg("20-nt guide", SPACER, "cbc", placeholder=True,
            feature=feature("guide", "guide_barcode", "unknown")),
        seg("scaffold", crispr.SCAFFOLD_V1[:24], "r2"),
        seg("right vector / genomic flank", "X" * 12, "w1", placeholder=True),
        seg("PCR1 reverse site", revcomp(MOFFAT_PCR1_R), "r2"),
    ], name="PCR1 enriched guide locus (~600 bp)")


def two_step_library() -> Construct:
    """Finished dual-index LCV2::TKOv3 PCR2 library from the printed primer sequences."""
    con = Construct([
        seg("P5", il.P5, "p5"),
        seg("i5", MOFFAT_I5, "cbc",
            feature=feature("sample_index_i5", "sample_index", "fixed")),
        seg("Read 1 site", il.TRUSEQ_READ1, "s5"),
        seg("PCR2 hU6 site", MOFFAT_PCR2_SITE_F, "s5"),
        seg("20-nt guide", SPACER, "cbc", placeholder=True,
            feature=feature("guide", "guide_barcode", "unknown")),
        seg("PCR2 scaffold site", revcomp(MOFFAT_PCR2_SITE_R), "r2"),
        seg("Read 2 site", revcomp(il.TRUSEQ_READ2), "t7"),
        seg("i7", revcomp(MOFFAT_I7), "cbc",
            feature=feature("sample_index_i7", "sample_index", "fixed")),
        seg("P7", revcomp(il.P7), "p7"),
    ], name="Moffat two-step PCR library")
    problems = sp.verify(con, two_step_sequencing_primers(),
                         required_roles=("Read 1", "Index 1 (i7)", "Index 2 (i5)"))
    if problems:
        raise ValueError("invalid two-step library: " + "; ".join(problems))
    return con


def two_step_sequencing_primers():
    return [sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["I2"]]


def _pcr(graph, prefix, template, final, primers, label, reviewed):
    forward = graph.input(prefix + "_forward", label + " forward primer", MaterialKind.OLIGO)
    reverse = graph.input(prefix + "_reverse", label + " reverse primer", MaterialKind.OLIGO)
    polymerase = graph.input(prefix + "_polymerase", label + " polymerase",
                            MaterialKind.CATALYST)
    library = graph.product(prefix + "_library", label + " sequencing library",
                            MaterialKind.LIBRARY, state=final)
    graph.react(prefix + "_pcr", label, StepKind.PCR,
                ((template, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (library,), geometry=PCRGeometry.unavailable(
                    final, reason="The genomic/vector template is larger than the exact final amplicon construct."))
    graph.sequence(library, SequencingPlan(library=final, primers=tuple(primers)))
    reviewed.extend((forward.id, reverse.id, polymerase.id, prefix + "_pcr", library.id))
    return library


def one_pcr_graph(authority: str, *, name: str = "Joung one-PCR CRISPR screen readout"):
    graph = ProtocolGraph(name, authority=authority)
    gdna = graph.input("genomic_dna", "genomic DNA from the pooled screen",
                       MaterialKind.BIOLOGICAL,
                       state=molecule("screen genomic DNA", (role("integrated guide cassette"),)),
                       primary=True)
    reviewed = ["genomic_dna"]
    _pcr(graph, "readout", gdna, one_pcr_library(), one_pcr_sequencing_primers(),
         "One-PCR guide-cassette readout", reviewed)
    graph.declare_reviewed_scope(reviewed, authority)
    return graph.require_publishable()


def two_step_graph(authority: str, *, name: str = "Moffat two-step CRISPR screen readout"):
    graph = ProtocolGraph(name, authority=authority)
    gdna = graph.input("genomic_dna", "genomic DNA from the LCV2::TKOv3 screen",
                       MaterialKind.BIOLOGICAL,
                       state=molecule("screen genomic DNA", (role("integrated guide cassette"),)),
                       primary=True)
    pcr1_f = graph.input("pcr1_forward", "PCR1 locus forward primer", MaterialKind.OLIGO)
    pcr1_r = graph.input("pcr1_reverse", "PCR1 locus reverse primer", MaterialKind.OLIGO)
    pcr1_pol = graph.input("pcr1_polymerase", "PCR1 polymerase", MaterialKind.CATALYST)
    intermediate = graph.product("pcr1_product", "approximately 600-bp guide-locus enrichment product",
                                 MaterialKind.MOLECULAR, state=two_step_intermediate())
    graph.react("pcr1", "PCR1 enriches the integrated guide locus", StepKind.PCR,
                ((gdna, InputRole.TEMPLATE), (pcr1_f, InputRole.STRUCTURAL),
                 (pcr1_r, InputRole.STRUCTURAL), (pcr1_pol, InputRole.CATALYST)),
                (intermediate,), geometry=PCRGeometry.unavailable(
                    two_step_intermediate(), reason="The surrounding integrated-vector/genomic template is schematic."))
    reviewed = ["genomic_dna", "pcr1_forward", "pcr1_reverse", "pcr1_polymerase", "pcr1"]
    _pcr(graph, "pcr2", intermediate, two_step_library(), two_step_sequencing_primers(),
         "PCR2 adds the dual-index sequencing arms", reviewed)
    graph.declare_reviewed_scope(reviewed, authority)
    return graph.require_publishable()


def vector_screen_graph(version: str, authority: str, *, name: str):
    graph = ProtocolGraph(name, authority=authority)
    vector = graph.input("vector", f"lentiCRISPR {version} transfer vector",
                         MaterialKind.MOLECULAR,
                         state=molecule(f"lentiCRISPR {version} plasmid", (role("BsmBI cloning cassette"),)),
                         primary=True)
    guide = graph.input("guide_duplex", "annealed BsmBI-compatible guide oligos",
                        MaterialKind.OLIGO)
    ligase = graph.input("cloning_ligase", "guide-cloning DNA ligase",
                         MaterialKind.CATALYST)
    cloned = graph.product("cloned_vector", "guide-bearing transfer vector",
                           MaterialKind.MOLECULAR, state=cloned_guide_cassette())
    graph.react("clone_guide", "Ligate the guide duplex into the BsmBI-opened cassette",
                StepKind.LIGATION,
                ((vector, InputRole.TEMPLATE), (guide, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (cloned,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The two BsmBI-directed seams are checked in cloning_scene; a single adapter-end object does not represent a two-ended plasmid insertion."))
    packaging = graph.input("packaging_system", "lentiviral packaging system",
                            MaterialKind.CONDITION)
    virus = graph.product("virus", "guide-library lentiviral particles",
                          MaterialKind.MOLECULAR, state=cloned.state)
    graph.react("package", "Package the pooled transfer-vector library",
                StepKind.GENERIC,
                ((cloned, InputRole.TEMPLATE), (packaging, InputRole.CONDITION)),
                (virus,), changes_molecule=False)
    cells = graph.input("cells", "screen target-cell population", MaterialKind.BIOLOGICAL,
                        state=molecule("screen cells", (role("cellular genomes"),)), primary=True)
    integration = graph.input("integration_catalyst", "lentiviral integration machinery",
                              MaterialKind.CATALYST)
    screened = graph.product("screened_cells", "selected cells carrying integrated guides",
                             MaterialKind.MOLECULAR, state=cells.state)
    graph.react("integrate", "Transduce cells and integrate one guide cassette",
                StepKind.INTEGRATION,
                ((cells, InputRole.TEMPLATE), (virus, InputRole.STRUCTURAL),
                 (integration, InputRole.CATALYST)), (screened,))
    extraction = graph.input("extraction_condition", "genomic-DNA extraction condition",
                             MaterialKind.CONDITION)
    gdna = graph.product("genomic_dna", "screen genomic DNA containing guide cassettes",
                         MaterialKind.MOLECULAR, state=screened.state)
    graph.react("extract", "Extract genomic DNA after pooled selection",
                StepKind.EXTRACTION,
                ((screened, InputRole.TEMPLATE), (extraction, InputRole.CONDITION)), (gdna,))
    reviewed = ["vector", "guide_duplex", "cloning_ligase", "clone_guide",
                "packaging_system", "package", "cells", "integration_catalyst",
                "integrate", "extraction_condition", "extract"]
    _pcr(graph, "readout", gdna, one_pcr_library(), one_pcr_sequencing_primers(),
         "Compatible one-PCR guide-cassette readout", reviewed)
    graph.declare_reviewed_scope(reviewed, authority)
    return graph.require_publishable()
