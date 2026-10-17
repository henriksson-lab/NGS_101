"""Source-owned protocol definition and reaction graph."""
from __future__ import annotations

import illumina as il
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, feature, oligo, revcomp
from pcr import PCRGeometry
from protocol_specs import Definition, get
from adapter_ligation import AdapterLigationGeometry
from adapter_ligation import reported_adapter_end
from protocol_graph import InputRole
from protocol_graph import MaterialKind
from protocol_graph import ProtocolGraph
from protocol_graph import StepKind, SequencingPlan, Highlight, AnnotationSlot
from rna_library_operations import molecule

# Bae et al. 2023, Supplementary Table 5.  The paper publishes eight matched
# sample-index sets; set 1 is the representative molecule drawn here.  The three Ns
# are the per-molecule UMI, not an undisclosed sequence.
CODEC_SID5_1 = "CTTGAACGGACTGTCCAC"
CODEC_SID7_1 = "CACCGAGCGTTAGACTAC"
CODEC_ADAP5_1 = il.P5 + CODEC_SID5_1 + "T"
CODEC_ADAP7_1 = il.P7 + CODEC_SID7_1 + "T"
CODEC_LINKER_CORE_5 = "GCTAGACTCTGACGTGTTGATCCTCGAAGC"
CODEC_LINKER_CORE_7 = revcomp(CODEC_LINKER_CORE_5)
CODEC_BRID5_1 = (revcomp(CODEC_SID5_1) + "NNN" + revcomp(il.TRUSEQ_READ1)
                 + "TTACATAGTTATCC" + CODEC_LINKER_CORE_5)
CODEC_BRID7_1 = (revcomp(CODEC_SID7_1) + "NNN" + revcomp(il.TRUSEQ_READ2)
                 + "CAATCTATAAGTT" + CODEC_LINKER_CORE_7)

def codec_oligos():
    """Copyable representative CODEC adapter set, preserving ordering modifications."""
    return (
        oligo("LD4-adap5-1", [Segment("published ordering sequence",
                                      CODEC_ADAP5_1[:-1] + "*T", placeholder=True)],
              three="-3′"),
        oligo("LD4-adap7-1", [Segment("published ordering sequence",
                                      CODEC_ADAP7_1[:-1] + "*T", placeholder=True)],
              three="-3′"),
        oligo("LD4-brid5-1", [Segment("5-prime bridge", CODEC_BRID5_1)],
              five="5′-/5Phos/-", three="-3′"),
        oligo("LD4-brid7-1", [Segment("7-prime bridge", CODEC_BRID7_1)],
              five="5′-/5Phos/-", three="-3′"),
    )


def codec_adapter_scene() -> Scene:
    """Source Table 5 set-1 bridge pair; the central duplex is checked."""
    sc = Scene()
    sc.strand("bridge5", [Segment("SID5 complement", revcomp(CODEC_SID5_1)),
        Segment("UMI5", "NNN", "umi", feature=feature("adapter_umi5", "umi", "random")),
        Segment("Read 1 site complement", revcomp(il.TRUSEQ_READ1), "r1"),
        Segment("bridge5 spacer", "TTACATAGTTATCC"),
        Segment("central duplex5", CODEC_LINKER_CORE_5)], label="LD4-brid5-1", mod5="P-")
    sc.anneal("bridge7", [Segment("SID7 complement", revcomp(CODEC_SID7_1)),
        Segment("UMI7", "NNN", "umi", feature=feature("adapter_umi7", "umi", "random")),
        Segment("Read 2 site complement", revcomp(il.TRUSEQ_READ2), "r2"),
        Segment("bridge7 spacer", "CAATCTATAAGTT"),
        Segment("central duplex7", CODEC_LINKER_CORE_7)],
        to="bridge5", pair=("central duplex7", "central duplex5"), label="LD4-brid7-1", mod5="P-")
    return sc


def codec_final_library() -> Construct:
    """Primer-addressable linear view of the linked product in Fig. 1b.

    The source establishes the outward read order on both sides of the central bridge:
    primer, 3-nt UMI, 18-nt matched sample index, terminal T, original-strand copy.
    The bridge itself stays a role token because the figure does not print the extended
    product as one contiguous strand; the four input oligos above remain fully exact.
    """
    sid_set = "CODEC Supplementary Table 5 matched LD4 sample-index sets"
    return Construct([
        Segment("P7", il.P7, "p7"),
        Segment("Crick-strand copy", "x" * 28, placeholder=True),
        Segment("adap7 terminal T", "T"),
        Segment("SID7", CODEC_SID7_1, "cbc",
                feature=feature("sample_index_sid7", "sample_index", "whitelist",
                                whitelist=sid_set)),
        Segment("Read 2 UMI", "N" * 3, "umi",
                feature=feature("umi_r2", "umi", "random")),
        Segment("Read 2 primer site reverse complement", revcomp(il.TRUSEQ_READ2), "r2"),
        Segment("CODEC bridge linker", "[CODEC bridge linker]", "me",
                placeholder=True),
        Segment("Read 1 primer site", il.TRUSEQ_READ1, "r1"),
        Segment("Read 1 UMI", "N" * 3, "umi",
                feature=feature("umi_r1", "umi", "random")),
        Segment("SID5", CODEC_SID5_1, "cbc",
                feature=feature("sample_index_sid5", "sample_index", "whitelist",
                                whitelist=sid_set)),
        Segment("adap5 terminal T", "T"),
        Segment("Watson-strand copy", "X" * 28, placeholder=True),
        Segment("P5 reverse complement", il.P5_RC, "p5"),
    ], name="CODEC linked-duplex sequencing library")


def _p(name: str, base: str, n: int, tag: str | None = None) -> Segment:
    return Segment(name, base * n, tag, placeholder=True)


def codec_linked_duplex() -> Scene:
    """One CODEC product containing readable copies of both original strands."""
    parts = [_p("Watson copy", "X", 28), _p("CODEC linkage", "L", 14, "me"),
             _p("Crick copy", "x", 28)]
    sc = Scene(); sc.strand("linked duplex", parts, label="single CODEC cluster molecule")
    sc.junction("linked duplex", "Watson copy", "CODEC linkage", "strand linkage")
    sc.junction("linked duplex", "CODEC linkage", "Crick copy", "strand linkage")
    sc.note("linked duplex", "one read pair carries both complements of the original duplex")
    sc.labels("linked duplex")
    return sc



DEFINITION = Definition("CODEC — single-duplex error correction", "10.1038/s41588-023-01376-0", "Bae et al., Nature Genetics (2023)", "An adapter quadruplex and strand-displacing extension place readable copies of both strands of one original DNA duplex into a single Illumina cluster molecule.", "custom", ("Watson copy", "CODEC linkage", "Crick copy"), "Standard TruSeq Read 1 and Read 2 primers bind within the two published bridge oligos and read outward through a 3-nt UMI, an 18-nt inline sample index and one original strand copy.", "DNA", "Supplementary Table 5 publishes all four adapter oligos in eight matched sample-index sets. Set 1 is shown as a representative; the linked post-extension product remains role-scaled where the source figure does not print one contiguous strand.")


from protocol_specs import workflow


def _source_workflow():
    return workflow(molecule("original genomic duplex", ("Watson strand", "Crick strand")), (
        ("End-repair and ligate the preassembled CODEC adapter quadruplex",
         molecule("adapter-bound duplex", ("left CODEC adapter", "original duplex", "right CODEC adapter"), joined=True,
                  junctions=(("left CODEC adapter", "original duplex", "adapter ligation"),
                             ("original duplex", "right CODEC adapter", "adapter ligation"))),
         "The four adapter oligos, modifications and matched index sets are printed in Supplementary Table 5."),
        ("Strand-displacing extension copies and physically links both complements",
         codec_linked_duplex(), "Both original strands become recoverable in one cluster molecule."),
        ("PCR and outward Illumina sequencing", codec_linked_duplex(),
         "Reads recover insert, UMI and index information for within-pair consensus."),
    ))

SPEC = get('codec', definition=DEFINITION,
           workflow_fn=_source_workflow,
           final_override=codec_final_library(),
           primers_override=(sp.TRUSEQ["R1"], sp.TRUSEQ["R2"]),
           oligos_fn=codec_oligos)


def _local_graph(d, final, primers):
    source = f"{d.source_label}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    duplex = graph.input("duplex", "original genomic Watson/Crick duplex",
                         MaterialKind.BIOLOGICAL,
                         state=molecule("original genomic duplex", ("Watson strand", "Crick strand")),
                         primary=True)
    methods = source + "; Methods 'DNA samples' and 'CODEC' (PMC10181940 JATS lines 202, 210)"
    covaris = graph.input("covaris", "Covaris ultrasonication to a mean of 150 bp (not cfDNA)",
                          MaterialKind.CONDITION, source=methods)
    sheared = graph.product("sheared", "sheared genomic duplex", MaterialKind.MOLECULAR,
                            state=duplex.state)
    graph.react("shear", "Shear genomic DNA (2.5-20 ng input to CODEC)", StepKind.FRAGMENTATION,
                ((duplex, InputRole.TEMPLATE), (covaris, InputRole.CONDITION)), (sheared,),
                source=methods)
    repair = graph.input("end_repair_enzyme", "NEBNext Ultra II End Prep (E7103L)", MaterialKind.CATALYST,
                         source=methods)
    end_prepped = graph.product("end_prepped", "end-repaired, dA-tailed genomic duplex",
                                MaterialKind.MOLECULAR, state=duplex.state)
    graph.react("end_repair", "End-repair and dA-tail the original duplex", StepKind.END_REPAIR,
                ((sheared, InputRole.TEMPLATE), (repair, InputRole.CATALYST)), (end_prepped,))
    pre_beads = graph.input("pre_ligation_ampure", "1.8x AMPure XP, elute in 40 µl low TE",
                            MaterialKind.SELECTION, source=methods + " exception (1)")
    repaired = graph.product("repaired", "purified end-prepped duplex", MaterialKind.MOLECULAR,
                             state=duplex.state)
    end_prep_waste = graph.product("end_prep_waste", "End Prep reaction components",
                                   MaterialKind.WASTE)
    graph.react("pre_ligation_cleanup", "AMPure cleanup before ligation (CODEC exception 1)",
                StepKind.CLEANUP,
                ((end_prepped, InputRole.TEMPLATE), (pre_beads, InputRole.SELECTION)),
                (repaired, end_prep_waste), source=methods)
    graph.discard(end_prep_waste, "removed in the 1.8x bead supernatant")
    adapter = graph.input("codec_adapter", "four-oligo CODEC adapter quadruplex (annealed 85 °C 3 min, -1 °C/min to 20 °C, 12 h RT; 500 nM in ligation)",
                          MaterialKind.ADAPTER, state=codec_adapter_scene(),
                          source=source + "; Supplementary Table 5")
    ligase = graph.input("ligase", "NEBNext Ultra II ligation mix (ligation extended to 1 h)", MaterialKind.CATALYST)
    deadenylase = graph.input("deadenylase", "5′-deadenylase added to ligation", MaterialKind.CATALYST,
                             source="Bae et al. 2023 CODEC Methods; NEB M0331S")
    adapted = graph.product("adapted", "CODEC-adapter-bound duplex", MaterialKind.MOLECULAR,
                            state=molecule("adapter-bound duplex", (
                                "left CODEC adapter", "original duplex", "right CODEC adapter"),
                                junctions=(("left CODEC adapter", "original duplex", "adapter ligation"),
                                           ("original duplex", "right CODEC adapter", "adapter ligation"))))
    blank = graph.product("blank_ligation", "Putative blank-ligation byproduct (no insert)", MaterialKind.WASTE)
    double = graph.product("double_ligation", "Putative double-ligation byproduct", MaterialKind.WASTE)
    graph.react("ligate", "Ligate the preassembled CODEC adapter quadruplex",
                StepKind.LIGATION,
                ((repaired, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST), (deadenylase, InputRole.CATALYST)), (adapted, blank, double),
                geometry=AdapterLigationGeometry.known(
                    "cohesive",
                    *reported_adapter_end(
                        "LD4-adap5-1 arm", kind="cohesive", polarity="3-prime",
                        terminal_sequence="T", body_sequence=CODEC_ADAP5_1[:-1],
                        statement="Supplementary Table 5 prints the phosphorothioate-linked terminal T."
                    ).scenes,
                    *reported_adapter_end(
                        "LD4-adap7-1 arm", kind="cohesive", polarity="3-prime",
                        terminal_sequence="T", body_sequence=CODEC_ADAP7_1[:-1],
                        statement="Supplementary Table 5 prints the phosphorothioate-linked terminal T."
                    ).scenes,
                    statement=("The paper dA-tails the insert and Supplementary Table 5 prints a "
                               "terminal phosphorothioate-linked T on each ligating adapter arm.")))
    graph.discard(blank, "Extended Data Fig. 2b proposes insert-free T/T mismatched ligation; this is a likely mechanism, not a proved exact product.")
    graph.discard(double, "Extended Data Fig. 2b proposes double ligation involving T/T mismatch; the product is shown by role only.")
    strand_polymerase = graph.input("strand_polymerase", "phi29 DNA polymerase (M0269L), 0.2 mM dNTP, 30 °C 20 min",
                                    MaterialKind.CATALYST, source=methods + " exception (4)")
    linked = graph.product("linked", "physically linked copies of both original complements",
                           MaterialKind.MOLECULAR, state=codec_linked_duplex())
    intermolecular = graph.product("intermolecular", "Putative intermolecular-extension byproduct", MaterialKind.WASTE)
    graph.react("link_duplex", "Copy and physically link both complements by strand displacement",
                StepKind.EXTENSION,
                ((adapted, InputRole.TEMPLATE), (strand_polymerase, InputRole.CATALYST)),
                (linked, intermolecular))
    graph.discard(intermolecular, "Extended Data Fig. 2b proposes extension on another ligation product rather than the opposite strand; no exact sequence is assigned.")
    ext_beads = graph.input("extension_ampure", "0.75x AMPure XP", MaterialKind.SELECTION,
                            source=methods + " exception (4)")
    linked_clean = graph.product("linked_clean", "purified linked-duplex products", MaterialKind.MOLECULAR,
                                 state=codec_linked_duplex())
    ext_waste = graph.product("extension_waste", "extension reaction components and short material",
                              MaterialKind.WASTE)
    graph.react("extension_cleanup", "0.75x AMPure cleanup after strand-displacing extension",
                StepKind.CLEANUP, ((linked, InputRole.TEMPLATE), (ext_beads, InputRole.SELECTION)),
                (linked_clean, ext_waste), source=methods)
    graph.discard(ext_waste, "removed in the 0.75x bead supernatant")
    reviewed = ["duplex", "covaris", "shear", "end_repair_enzyme", "end_repair",
                "pre_ligation_ampure", "pre_ligation_cleanup", "codec_adapter", "ligase",
                "ligate", "strand_polymerase", "link_duplex", "extension_ampure", "extension_cleanup"]
    amplification_primers = graph.input("amplification_primers", "xGen Library Amplification Primer Mix (IDT 1077677)", MaterialKind.OLIGO,
                                        source="Bae et al. 2023 CODEC Methods")
    amplification_partner = graph.input("amplification_partner", "xGen Library Amplification Primer Mix partner (IDT 1077677)", MaterialKind.OLIGO,
                                        source="Bae et al. 2023 CODEC Methods; ordered strands unassigned")
    pcr_enzyme = graph.input("pcr_enzyme", "KAPA HiFi HotStart ReadyMix (KAPA Library Amplification manual, 2 min extension)", MaterialKind.CATALYST,
                             source="Bae et al. 2023 CODEC Methods; Roche 07958935001")
    amplified = graph.product("amplified", "amplified CODEC library before cleanup", MaterialKind.MOLECULAR, state=final)
    graph.react("library_pcr", "Amplify the concatenated library with xGen primers", StepKind.PCR,
                ((linked_clean, InputRole.TEMPLATE), (amplification_primers, InputRole.STRUCTURAL),
                 (amplification_partner, InputRole.STRUCTURAL), (pcr_enzyme, InputRole.CATALYST)), (amplified,),
                geometry=PCRGeometry.unavailable(final,
                    reason="The model has not assigned the ordered xGen 1077677 primer strands or a contiguous post-extension amplicon; the original CODEC Methods names the commercial mix. This is undeclared model data, not accepted source absence."))
    post_beads = graph.input("post_pcr_ampure", "two 0.65x AMPure XP cleanups (cfDNA: second as 0.5x/0.65x double-sided)",
                             MaterialKind.SELECTION, source=methods + " exception (6)")
    library = graph.product("library", "CODEC linked-duplex sequencing library", MaterialKind.LIBRARY, state=final)
    pcr_waste = graph.product("pcr_waste", "PCR reaction components and short products", MaterialKind.WASTE)
    graph.react("post_pcr_cleanup", "Two 0.65x AMPure cleanups after PCR", StepKind.CLEANUP,
                ((amplified, InputRole.TEMPLATE), (post_beads, InputRole.SELECTION)),
                (library, pcr_waste), source=methods)
    graph.discard(pcr_waste, "removed in the 0.65x bead supernatants")
    graph.sequence(library, SequencingPlan(library=final, primers=tuple(primers),
                    acquired_roles=("Read 1", "Read 2"),
                    read_lengths={"Read 1":166,"Read 2":166},
                    read_uses={"Read 1":"extract_then_include", "Read 2":"extract_then_include"}))
    reviewed.extend(("deadenylase", "amplification_primers", "amplification_partner", "pcr_enzyme", "library_pcr", "amplified",
                     "post_pcr_ampure", "post_pcr_cleanup", "library"))
    graph.annotate(Highlight(target="library", title="Representative adapter set and sequencing run",
        explanation="Supplementary Table 5 prints eight matched four-oligo sets; this page displays set 1. Extended Data Fig. 1f and Supplementary Text also describe NNNA/NNNT read-start variants; the fixed three-position UMI model applies to set 1. The original Methods reports NovaSeq 6000 paired 166- or 260-cycle CODEC runs; the read spans shown select the 166-cycle option. Sample indices are inline within Read 1/Read 2, rather than separately acquired index reads. The concatenated endpoint remains a region-level model, not a source-printed contiguous strand.",
        slot=AnnotationSlot.AFTER_WORKFLOW))
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def protocol_graph():
    return _local_graph(DEFINITION, SPEC.FINAL_LIBRARY, SPEC.SEQ_PRIMERS)

SPEC.protocol_graph = protocol_graph
