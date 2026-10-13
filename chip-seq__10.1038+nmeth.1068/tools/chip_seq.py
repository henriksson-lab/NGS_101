"""Robertson et al. 2007 ChIP-seq and original single-read Solexa library."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import illumina as il
import seqprimers as sp
from chemdraw import Construct, Row, Scene, Segment, revcomp
from chromatin_epigenetics import antibody_enzyme_rows
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from protocol_model import molecule, role


def seg(name, top, tag=None, **kw):
    return Segment(name, top, tag, **kw)


TITLE = "ChIP-seq — conventional immunoprecipitation library"
NOTES = "01_chip-seq.html"
SOURCE = 'Defining protocol: <a href="https://doi.org/10.1038/nmeth.1068">Robertson et al., <i>Nature Methods</i> (2007)</a>.'
SUMMARY = "Crosslink and sonicate chromatin, immunoprecipitate target-associated DNA, reverse crosslinks, and prepare the recovered fragments for the original single-read Solexa Genome Analyzer."
CAVEAT = "This 2007 protocol predates TruSeq and sample indexing. It used Illumina’s Genomic DNA Sample Prep adapter and PCR primers 1.1/2.1, followed by one genomic-DNA sequencing read—not a modern paired-end dual-index endpoint."

GENOMIC_ADAPTER_BOTTOM = il.TRUSEQ_READ1
GENOMIC_PCR2 = "CAAGCAGAAGACGGCATACGAGCTCTTCCGATCT"
# PCR primer 2.1 has one 3'-terminal T paired to the insert's dA junction; the rest
# determines the phosphorylated adapter oligo rather than duplicating it independently.
GENOMIC_ADAPTER_TOP = revcomp(GENOMIC_PCR2)[1:]


def adapter_scene() -> Scene:
    sc = Scene()
    sc.strand("bottom", [
        seg("Read 1 fork arm", GENOMIC_ADAPTER_BOTTOM[:-13], "r1"),
        seg("genomic stem complement", GENOMIC_ADAPTER_BOTTOM[-13:-1], "r1"),
        seg("3′ T", GENOMIC_ADAPTER_BOTTOM[-1], "r1"),
    ], label="Genomic Adapter Oligo")
    sc.anneal("top", [
        seg("genomic stem", GENOMIC_ADAPTER_TOP[:12], "p7"),
        seg("genomic fork arm", GENOMIC_ADAPTER_TOP[12:], "p7"),
    ], to="bottom", pair=("genomic stem", "genomic stem complement"),
              label="Genomic Adapter Oligo", mod5="p",
              unpaired=("genomic fork arm",))
    sc.ligation_end("bottom", "3′ T", kind="cohesive", polarity="3-prime",
                    owner="adapter", text="3′-T ligation overhang")
    return sc


SEQ_PRIMERS = (sp.custom(
    "Read 1", "Genomic DNA sequencing primer", il.TRUSEQ_READ1,
    'Illumina "Genomic DNA Sample Prep Kit" / adapter-sequence guide'),)


def final_library() -> Construct:
    lib = Construct([
        seg("P5 / PCR primer 1.1", il.P5, "p5"),
        seg("Read 1 remainder", il.TRUSEQ_READ1[4:], "r1"),
        seg("ChIP-enriched DNA", "X" * 36, placeholder=True),
        seg("dA junction", "A"),
        seg("Genomic adapter / PCR primer 2.1 arm", GENOMIC_ADAPTER_TOP, "p7"),
    ], name="2007 ChIP-seq library")
    problems = sp.verify(lib, SEQ_PRIMERS, required_roles=("Read 1",))
    if problems:
        raise ValueError("invalid original ChIP-seq library: " + "; ".join(problems))
    return lib


FINAL_LIBRARY = final_library()
FINAL_CAPTION = "Original unindexed single-read Solexa library made with Genomic Adapter Oligo Mix and PCR primers 1.1/2.1."
SEQUENCING_INTRO = "The defining experiment generated single-end approximately 27-nt reads. It had no sample index and no Read 2."


def sections():
    return [
        ("Crosslink and fragment chromatin",
         [Row(chunks=[("protein—DNA chromatin → formaldehyde crosslink → sonicated fragments", None, False)])],
         "Mechanical shearing defines fragment ends independently of the target protein."),
        ("Immunoprecipitate the target",
         antibody_enzyme_rows("bead capture", "retain target-associated fragments"),
         "A target-specific antibody enriches crosslinked protein–DNA complexes; crosslinks are then reversed."),
        ("Repair and add one 3-prime dA to each end",
         Scene.duplex([seg("3-prime dA", "A"), seg("ChIP-enriched DNA", "X" * 36, placeholder=True), seg("3-prime dA", "A")], label="dA-tailed ChIP DNA").rows(),
         "The 3-prime dA overhangs are shown before the forked adapter is introduced."),
        ("Ligate the Genomic Adapter", adapter_scene().rows(),
         "The original Illumina single-read fork has a 3′-T ligation overhang. PCR primers 1.1 and 2.1 complete the flow-cell arms."),
    ]


def protocol_graph():
    graph = ProtocolGraph(TITLE, authority=SOURCE)
    cells = graph.input("cells", "input cells", MaterialKind.BIOLOGICAL,
                        state=molecule("cells", (role("native protein-DNA chromatin"),)),
                        primary=True)
    formaldehyde = graph.input("formaldehyde", "formaldehyde cross-linking condition",
                               MaterialKind.CONDITION)
    fixed = graph.product("fixed_chromatin", "formaldehyde-cross-linked chromatin",
                          MaterialKind.MOLECULAR, state=cells.state)
    graph.react("crosslink", "Cross-link protein-DNA chromatin",
                StepKind.CROSS_LINKING,
                ((cells, InputRole.TEMPLATE), (formaldehyde, InputRole.CONDITION)),
                (fixed,))
    sonication = graph.input("sonication", "chromatin-sonication condition",
                             MaterialKind.CONDITION)
    fragments = graph.product("chromatin_fragments", "sonicated cross-linked chromatin",
                              MaterialKind.MOLECULAR, state=fixed.state)
    graph.react("fragment", "Sonicate cross-linked chromatin", StepKind.FRAGMENTATION,
                ((fixed, InputRole.TEMPLATE), (sonication, InputRole.CONDITION)),
                (fragments,))
    antibody = graph.input("target_antibody", "target-specific antibody and capture beads",
                           MaterialKind.SELECTION)
    selected = graph.product("selected_chromatin", "target-associated chromatin fragments",
                             MaterialKind.MOLECULAR, state=fragments.state)
    flowthrough = graph.product("ip_flowthrough", "unbound chromatin fragments",
                                MaterialKind.WASTE)
    graph.react("immunoprecipitate", "Immunoprecipitate target-associated chromatin",
                StepKind.AFFINITY_CAPTURE,
                ((fragments, InputRole.TEMPLATE), (antibody, InputRole.SELECTION)),
                (selected, flowthrough))
    graph.discard(flowthrough, "is not retained by the target-directed immunoprecipitation")
    reversal = graph.input("crosslink_reversal", "cross-link reversal condition",
                           MaterialKind.CONDITION)
    dna = graph.product("recovered_dna", "de-cross-linked ChIP-enriched DNA",
                        MaterialKind.MOLECULAR, state=selected.state)
    graph.react("reverse_crosslinks", "Reverse cross-links and recover ChIP DNA",
                StepKind.DEPROTEINIZATION,
                ((selected, InputRole.TEMPLATE), (reversal, InputRole.CONDITION)),
                (dna,))
    repair = graph.input("end_repair_mix", "end-repair and phosphorylation enzymes",
                         MaterialKind.CATALYST)
    repaired = graph.product("end_repaired", "end-repaired phosphorylated ChIP DNA",
                             MaterialKind.MOLECULAR, state=dna.state)
    graph.react("end_repair", "Repair and phosphorylate recovered DNA ends",
                StepKind.END_REPAIR,
                ((dna, InputRole.TEMPLATE), (repair, InputRole.CATALYST)), (repaired,))
    tailing = graph.input("tailing_polymerase", "dA-tailing polymerase",
                          MaterialKind.CATALYST)
    datp = graph.input("datp", "dATP tailing nucleotide", MaterialKind.MOLECULAR)
    tailed = graph.product(
        "da_tailed", "3′-dA-tailed ChIP DNA", MaterialKind.MOLECULAR,
        state=Scene.duplex([seg("3-prime dA", "A"),
                            seg("ChIP-enriched DNA", "X" * 36, placeholder=True),
                            seg("opposite 3-prime dA", "A")], label="dA-tailed ChIP DNA"))
    graph.react("da_tail", "Add one 3′ dA to each repaired fragment end",
                StepKind.TAILING,
                ((repaired, InputRole.TEMPLATE), (datp, InputRole.STRUCTURAL),
                 (tailing, InputRole.CATALYST)), (tailed,))
    adapter = graph.input("genomic_adapter", "Illumina Genomic Adapter Oligo Mix",
                          MaterialKind.ADAPTER, state=adapter_scene())
    ligase = graph.input("adapter_ligase", "adapter DNA ligase", MaterialKind.CATALYST)
    ligated = graph.product("adapter_ligated", "Genomic-Adapter-ligated ChIP DNA",
                            MaterialKind.MOLECULAR, state=FINAL_LIBRARY)
    graph.react("adapter_ligation", "Ligate the original Illumina Genomic Adapter",
                StepKind.LIGATION,
                ((tailed, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (ligated,),
                geometry=AdapterLigationGeometry.known(
                    "cohesive", adapter_scene(),
                    statement="The archived Genomic Adapter oligos expose the 3′ dT that pairs with the library fragment's 3′ dA."))
    size_matrix = graph.input("size_selection", "gel size-selection matrix",
                              MaterialKind.SELECTION)
    sized = graph.product("size_selected", "size-selected adapter-ligated ChIP DNA",
                          MaterialKind.MOLECULAR, state=FINAL_LIBRARY)
    rejected = graph.product("rejected_sizes", "out-of-range library fragments",
                             MaterialKind.WASTE)
    graph.react("size_select", "Select the adapter-ligated library size range",
                StepKind.SELECTION,
                ((ligated, InputRole.TEMPLATE), (size_matrix, InputRole.SELECTION)),
                (sized, rejected))
    graph.discard(rejected, "falls outside the retained library-size range")
    primer1 = graph.input("pcr_primer_1_1", "Genomic PCR primer 1.1",
                          MaterialKind.OLIGO)
    primer2 = graph.input("pcr_primer_2_1", "Genomic PCR primer 2.1",
                          MaterialKind.OLIGO)
    polymerase = graph.input("pcr_polymerase", "library PCR polymerase",
                             MaterialKind.CATALYST)
    library = graph.product("sequencing_library", FINAL_LIBRARY.name,
                            MaterialKind.LIBRARY, state=FINAL_LIBRARY)
    graph.react("library_pcr", "Amplify with Genomic PCR primers 1.1 and 2.1",
                StepKind.PCR,
                ((sized, InputRole.TEMPLATE), (primer1, InputRole.STRUCTURAL),
                 (primer2, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (library,), geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY,
                    reason="The adapter-ligated pre-PCR fork is represented through its published end geometry rather than as one exact PCR template."))
    graph.sequence(library, SequencingPlan(
        library=FINAL_LIBRARY, primers=tuple(SEQ_PRIMERS),
        acquired_roles=("Read 1",), read_lengths={"Read 1": 27}))
    graph.declare_reviewed_scope((
        "cells", "formaldehyde", "crosslink", "sonication", "fragment",
        "target_antibody", "immunoprecipitate", "crosslink_reversal",
        "reverse_crosslinks", "end_repair_mix", "end_repair", "tailing_polymerase",
        "datp", "da_tail", "genomic_adapter", "adapter_ligase",
        "adapter_ligation", "size_selection", "size_select", "pcr_primer_1_1",
        "pcr_primer_2_1", "pcr_polymerase", "library_pcr", "sequencing_library"),
        SOURCE)
    return graph.require_publishable()
