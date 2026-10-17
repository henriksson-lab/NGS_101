"""ChIA-Drop preserves chromatin groups without proximity ligation.

The group is a collection of independent molecules sharing a GEM barcode, never one
covalent contact chimera.
"""
from __future__ import annotations

from types import SimpleNamespace

import illumina as il
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, feature
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_graph import (Disposition, Terminal, InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from protocol_model import chain, molecule, role


from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "lib"))

CHIA_COMPLEX = feature(
    "chia_drop.gem", "molecular_complex_barcode", "whitelist",
    whitelist="10x Chromium Genome v2 GemCode set",
    note="16-nt partition barcode shared by every copied fragment in one GEM")
CHIA_UMI = feature(
    "chia_drop.random_primer", "umi", "random",
    note="8-nt random priming sequence retained beside the GEM barcode")
CHIA_SAMPLE = feature(
    "chia_drop.sample", "sample_index", "unknown",
    note="library sample index; exact assigned oligo remains unresolved in this bounded review")



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
                          placeholder=True))
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
        ("Permeabilize nuclei and digest chromatin with HindIII in situ",
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
         "The source states end repair, A-tailing, PN-120262 adapter ligation, Sample Index PCR and size selection around 600 bp. Exact commercial oligos remain unresolved here; final arms are inferred."),
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
                    reason="The commercial intermediate and PCR primers are unresolved in this bounded review; the compatible inferred final primer sites are checked independently."))
    reviewed.extend((forward.id, reverse.id, polymerase.id, prefix + "_pcr",
                     library.id))
    return library


def _chia_graph(spec):
    source = spec.SOURCE
    graph = ProtocolGraph(spec.TITLE, authority=source)
    chromatin = graph.input("chromatin", "formaldehyde-cross-linked nucleus",
                            MaterialKind.BIOLOGICAL,
                            state=molecule("cross-linked nucleus", (role("multiplex chromatin complexes"),)),
                            primary=True)
    nuclease = graph.input("restriction_enzyme", "HindIII (selected published restriction branch)",
                           MaterialKind.CATALYST)
    digested = graph.product("digested", "restriction-fragmented cross-linked complexes",
                             MaterialKind.MOLECULAR, state=grouped_fragments("digested complexes"))
    graph.react("digest", "Digest chromatin with HindIII at 37 °C overnight; no proximity ligation",
                StepKind.DIGESTION,
                ((chromatin, InputRole.TEMPLATE), (nuclease, InputRole.CATALYST)),
                (digested,))
    sonication = graph.input("sonication", "complex-release sonication condition",
                             MaterialKind.CONDITION)
    released = graph.product("released", "released multi-fragment chromatin complexes",
                             MaterialKind.MOLECULAR, state=grouped_fragments("released complexes"))
    graph.react("release", "Release intact complexes by sonication",
                StepKind.FRAGMENTATION,
                ((digested, InputRole.TEMPLATE), (sonication, InputRole.CONDITION)),
                (released,))
    gem = graph.input("gem_condition", "Chromium Genome v2 GEM partition condition",
                      MaterialKind.CONDITION)
    partitioned = graph.product("partitioned", "complex-bearing GEMs, including possible collisions",
                                MaterialKind.MOLECULAR, state=released.state)
    no_template = graph.product("no_template_gems", "candidate GEMs lacking accessible chromatin template",
                                MaterialKind.MOLECULAR, inferred=True)
    graph.react("partition", "Partition chromatin complexes into GEMs",
                StepKind.PARTITION,
                ((released, InputRole.TEMPLATE), (gem, InputRole.CONDITION)),
                (partitioned, no_template))
    graph.mark_terminal(Terminal(no_template.id, Disposition.UNRESOLVED,
        "Supplementary Methods p. 5 associates short barcode-only reads with empty or inaccessible-template GEMs. Their detailed downstream products are not traced here; this is not a physical discard. Colliding complexes remain in the main amplification path."))
    bead_oligo = graph.input("bead_oligo", "GemCode barcode/random-primer bead oligo",
                             MaterialKind.OLIGO, inferred=True,
                             source=source + " (exact commercial oligo not resolved in this review)")
    polymerase = graph.input("isothermal_polymerase", "isothermal extension polymerase",
                             MaterialKind.CATALYST)
    copied = graph.product("copied", "independent copies sharing one GEM barcode",
                           MaterialKind.MOLECULAR,
                           state=grouped_fragments("barcoded copies", copied=True))
    graph.react("amplify", "Random-prime linear amplification at 30 °C overnight",
                StepKind.EXTENSION,
                ((partitioned, InputRole.TEMPLATE), (bead_oligo, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (copied,))
    reviewed = ["chromatin", "restriction_enzyme", "digest", "sonication", "release",
                "gem_condition", "partition", "bead_oligo", "isothermal_polymerase",
                "amplify"]
    recovered = graph.product("recovered", "pooled, purified barcoded DNA copies",
                              MaterialKind.MOLECULAR,
                              state=grouped_fragments("recovered independent copies", copied=True))
    recovery = graph.input("gem_recovery", "GEM breakup and DNA purification reagents",
                           MaterialKind.CONDITION)
    graph.react("recover", "Release barcoded amplicons from all GEMs and purify DNA",
                StepKind.MAGIC,
                ((copied, InputRole.TEMPLATE), (recovery, InputRole.CONDITION)), (recovered,),
                reason="Supplementary Methods p. 6 states release and purification but the kit-era recovery mechanism and exact termini are not resolved in this bounded review. Copies remain separate molecules; collision filtering occurs after sequencing.")
    end_reagents = graph.input("end_repair_reagents", "commercial end-repair reagents",
                               MaterialKind.CATALYST)
    repaired = graph.product("repaired", "end-repaired barcoded DNA", MaterialKind.MOLECULAR,
                             state=molecule("end-repaired copies", (role("independent barcoded DNA; exact ends unresolved"),)))
    graph.react("end_repair", "End-repair purified barcoded DNA", StepKind.END_REPAIR,
                ((recovered, InputRole.TEMPLATE), (end_reagents, InputRole.CATALYST)), (repaired,))
    tail_reagents = graph.input("tail_reagents", "commercial A-tailing reagents", MaterialKind.CATALYST)
    tailed = graph.product("tailed", "A-tailed barcoded DNA", MaterialKind.MOLECULAR,
                           state=molecule("A-tailed copies", (role("barcoded DNA after source-stated A-tailing; exact ends unresolved"),)))
    graph.react("a_tail", "A-tail end-repaired DNA", StepKind.TAILING,
                ((repaired, InputRole.TEMPLATE), (tail_reagents, InputRole.CATALYST)), (tailed,))
    adapter = graph.input("multiplex_adapter", "Chromium Multiplex PN-120262 adapter", MaterialKind.ADAPTER)
    ligase = graph.input("multiplex_ligase", "commercial adapter-ligation reagents", MaterialKind.CATALYST)
    ligated = graph.product("ligated", "adapter-ligated barcoded DNA", MaterialKind.MOLECULAR,
                            state=molecule("ligated copies", (role("one genomic copy bearing barcode and commercial adapter; exact junction unresolved"),)))
    graph.react("adapter_ligation", "Ligate library adapter to independent DNA copies", StepKind.LIGATION,
                ((tailed, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL), (ligase, InputRole.CATALYST)),
                (ligated,), geometry=AdapterLigationGeometry.unavailable(
                    "PN-120262 ligation is stated in Supplementary Methods p. 6. Exact adapter termini and junction are unresolved here; the historical vendor original is a deferred follow-up. This is library-adapter ligation, not proximity ligation."))
    amplified = _finish_pcr(graph, ligated, spec, "library", reviewed,
                            "INFERRED — Sample Index PCR; compatible final arms shown")
    size_condition = graph.input("size_condition", "post-PCR size selection around 600 bp", MaterialKind.SELECTION)
    selected = graph.product("selected", "size-selected ChIA-Drop library around 600 bp",
                             MaterialKind.LIBRARY, state=spec.FINAL_LIBRARY)
    excluded = graph.product("size_excluded", "material outside the recovered size fraction", MaterialKind.WASTE)
    graph.react("size_select", "Size-select PCR products around 600 bp", StepKind.SELECTION,
                ((amplified, InputRole.TEMPLATE), (size_condition, InputRole.SELECTION)),
                (selected, excluded), changes_molecule=False)
    graph.discard(excluded, "Excluded by the stated post-PCR size selection around 600 bp; precise cutoffs and composition remain unresolved.")
    graph.sequence(selected, SequencingPlan(library=spec.FINAL_LIBRARY,
                                            primers=tuple(spec.SEQ_PRIMERS),
                                            read_lengths=spec.READ_LENGTHS,
                                            read_uses={"Read 1": "extract_then_include", "Read 2": "undetermined"}))
    reviewed.extend(("gem_recovery", "recover", "end_repair_reagents", "end_repair",
                     "tail_reagents", "a_tail", "multiplex_adapter", "multiplex_ligase",
                     "adapter_ligation", "size_condition", "size_select", "selected"))
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def get():
    final, primers = chia_drop_library()
    spec = SimpleNamespace(
        TITLE="ChIA-Drop — ligation-free droplet mapping of multiplex chromatin complexes",
        NOTES="01_chia_drop.html",
        SOURCE='<a href="https://doi.org/10.1038/s41586-019-0949-1">Zheng et al., <i>Nature</i> (2019)</a>; Supplementary Methods pp. 4–6 and Supplementary Note 1 pp. 15–17.',
        SUMMARY="Crosslinked chromatin complexes are fragmented but not proximity-ligated. In a Chromium GEM, independent fragments from one physical complex undergo random-primed linear amplification with the same 16-nt barcode; reads are grouped by that shared partition identity.",
        CAVEAT="A GEM barcode identifies a droplet, not necessarily one biological complex: partition collisions are separated computationally where possible. <inf>The public final molecule uses the compatible standard Illumina arm geometry as an inferred illustration; exact commercial oligos and the historical kit guides remain unresolved in this bounded review.</inf> The graph selects the published HindIII branch; MboI and sonication-only preparations are alternatives.",
        oligos=lambda: (),
        workflow=chia_drop_workflow,
        FINAL_LIBRARY=final,
        FINAL_CAPTION="One representative copied fragment is shown. The 16-nt barcode is repeated across separate molecules from the same GEM; the random 8-mer both primes copying and supplies a per-copy UMI. <inf>Standard Illumina arm sequences show the compatible v2 geometry, not a claimed transcription of proprietary kit oligos.</inf>",
        SEQ_PRIMERS=primers,
        READ_LENGTHS={"Read 1": 150, "Read 2": 150},
        SEQUENCING_INTRO="The defining data were sequenced 2 × 150 bp. Read 1 begins with the 16-bp GemCode and 8-bp random-primer UMI before genomic sequence; Read 2 is genomic. The published analysis used R1 only; R2 was acquired but reserved for future analysis. The sample-index cycle count remains unresolved; the eight-base index drawn here is inferred compatibility geometry. Primer landing positions on the compatible Illumina shell are computed below.",
        SEQUENCING_ENDING="",
        SEQUENCING_UNAVAILABLE=""
    )
    spec.protocol_graph = lambda: _chia_graph(spec)
    return spec


SPEC = get()
globals().update(vars(SPEC))


def protocol_graph():
    return _chia_graph(SPEC)
