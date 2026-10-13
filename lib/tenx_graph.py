"""Shared typed provenance for Chromium single-cell RNA library construction."""
from __future__ import annotations

from adapter_ligation import reported_adapter_end
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import molecule, role


READ2_ADAPTER_END = reported_adapter_end(
    "10x Read-2 adapter", kind="cohesive", polarity="3-prime", terminal_sequence="T",
    statement="The Chromium RNA workflows ligate a 3-prime-dT adapter to dA-tailed cDNA fragments.")


def partition_material(graph, material, reviewed, *, label="cells or nuclei"):
    """Partition one Chromium input and account for non-resolvable GEM products."""
    gems = graph.input("gems", "Chromium GEM partitioning condition",
                       MaterialKind.CONDITION)
    occupied = graph.product("occupied", f"occupied GEMs containing {label}",
                             MaterialKind.MOLECULAR, state=material.state)
    failed = graph.product("failed", "empty, multiplet or unresolved GEM products",
                           MaterialKind.WASTE)
    graph.react("partition", f"Partition {label} with barcoded gel beads",
                StepKind.PARTITION,
                ((material, InputRole.TEMPLATE), (gems, InputRole.CONDITION)),
                (occupied, failed), changes_molecule=False)
    graph.discard(failed, "does not retain a resolvable single-cell barcode")
    reviewed.extend(("gems", "partition"))
    return occupied


def chromium_pcr(graph, prefix, template, product_state, reviewed, *, label,
                 library=False, geometry=None, inferred_primers=False,
                 limitation="The page authority does not print both PCR-primer sequences."):
    """Add one explicitly primed Chromium PCR and return its product material."""
    source = graph.authority + (" (commercial primer bases not asserted here)"
                                if inferred_primers else "")
    forward = graph.input(prefix + "_forward", label + " forward primer",
                          MaterialKind.OLIGO, inferred=inferred_primers, source=source)
    reverse = graph.input(prefix + "_reverse", label + " reverse primer",
                          MaterialKind.OLIGO, inferred=inferred_primers, source=source)
    polymerase = graph.input(prefix + "_polymerase", label + " polymerase",
                             MaterialKind.CATALYST)
    product = graph.product(prefix + "_product", label + " product",
                            MaterialKind.LIBRARY if library else MaterialKind.MOLECULAR,
                            state=product_state)
    graph.react(prefix, label, StepKind.PCR,
                ((template, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (product,), geometry=(geometry or PCRGeometry.unavailable(
                    product_state, reason=limitation)))
    reviewed.extend((forward.id, reverse.id, polymerase.id, prefix))
    return product


def chromium_rna_graph(*, title, authority, final, primers, first_strand,
                       amplified_cdna, read_lengths, five_prime=False,
                       inferred_intermediates=False):
    """Build the common Chromium 3′/5′ gene-expression provenance graph.

    Version-specific modules still own every sequence and molecular state.  This function
    owns only operations common to the published workflows and therefore makes partition
    loss, GEM recovery, both PCRs, fragmentation, end preparation, adapter ligation and
    the sequencing endpoint consistent across versions.
    """
    g = ProtocolGraph(title, authority=authority)
    cells = g.input("cells", "cells containing poly(A) RNA", MaterialKind.BIOLOGICAL,
                    state=molecule("cells", (role("poly(A) RNA"),)), primary=True)
    reviewed = ["cells"]
    occupied = partition_material(g, cells, reviewed, label="cells")

    primer_label = ("poly(dT) RT primer" if five_prime
                    else "cell-barcode/UMI oligo(dT) gel-bead primer")
    primer = g.input("rt_primer", primer_label, MaterialKind.OLIGO,
                     inferred=inferred_intermediates,
                     source=authority + (" (intermediate bases inferred from the final construct)"
                                         if inferred_intermediates else ""))
    rt = g.input("rt", "template-switching reverse transcriptase", MaterialKind.CATALYST)
    first = g.product("first_strand", "first-strand transcript cDNA",
                      MaterialKind.MOLECULAR, state=first_strand,
                      inferred=inferred_intermediates,
                      source=authority + (" (intermediate structure inferred)"
                                          if inferred_intermediates else ""))
    g.react("reverse_transcribe", "Capture and reverse-transcribe poly(A) RNA",
            StepKind.REVERSE_TRANSCRIPTION,
            ((occupied, InputRole.TEMPLATE), (primer, InputRole.STRUCTURAL),
             (rt, InputRole.CATALYST)), (first,))

    tso_label = ("cell-barcode/UMI gel-bead template-switch oligo" if five_prime
                 else "template-switch oligo")
    tso = g.input("tso", tso_label, MaterialKind.OLIGO,
                  inferred=inferred_intermediates,
                  source=authority + (" (intermediate bases inferred from the final construct)"
                                      if inferred_intermediates else ""))
    switched = g.product("switched", "template-switched cell-barcoded cDNA",
                         MaterialKind.MOLECULAR, state=first_strand,
                         inferred=inferred_intermediates,
                         source=authority + (" (intermediate structure inferred)"
                                             if inferred_intermediates else ""))
    g.react("template_switch", "Install the opposite cDNA amplification handle",
            StepKind.TEMPLATE_SWITCHING,
            ((first, InputRole.TEMPLATE), (tso, InputRole.STRUCTURAL),
             (rt, InputRole.CATALYST)), (switched,))

    recovery = g.input("gem_recovery", "GEM break and cDNA recovery condition",
                       MaterialKind.CONDITION)
    recovered = g.product("recovered", "recovered cell-barcoded cDNA",
                          MaterialKind.MOLECULAR, state=first_strand)
    g.react("recover", "Break GEMs and recover cDNA", StepKind.EXTRACTION,
            ((switched, InputRole.TEMPLATE), (recovery, InputRole.CONDITION)),
            (recovered,), changes_molecule=False)

    pre_f = g.input("preamp_forward", "cDNA preamplification forward primer",
                    MaterialKind.OLIGO, inferred=True,
                    source=authority + " (commercial primer bases not asserted here)")
    pre_r = g.input("preamp_reverse", "cDNA preamplification reverse primer",
                    MaterialKind.OLIGO, inferred=True,
                    source=authority + " (commercial primer bases not asserted here)")
    pre_pol = g.input("preamp_polymerase", "cDNA preamplification polymerase",
                      MaterialKind.CATALYST)
    amplified = g.product("amplified", "amplified cell-barcoded cDNA",
                          MaterialKind.MOLECULAR, state=amplified_cdna,
                          inferred=inferred_intermediates,
                          source=authority + (" (intermediate structure inferred)"
                                              if inferred_intermediates else ""))
    g.react("preamp", "Amplify recovered full-length cDNA", StepKind.PCR,
            ((recovered, InputRole.TEMPLATE), (pre_f, InputRole.STRUCTURAL),
             (pre_r, InputRole.STRUCTURAL), (pre_pol, InputRole.CATALYST)),
            (amplified,), geometry=PCRGeometry.unavailable(
                first_strand, reason="The page authority does not print both cDNA preamplification primers."))

    fragmentase = g.input("fragmentase", "cDNA fragmentation enzyme",
                          MaterialKind.CATALYST)
    fragments = g.product("fragments", "fragmented barcode-bearing cDNA",
                          MaterialKind.MOLECULAR)
    g.react("fragment", "Fragment amplified cDNA", StepKind.FRAGMENTATION,
            ((amplified, InputRole.TEMPLATE), (fragmentase, InputRole.CATALYST)),
            (fragments,))
    end_mix = g.input("end_mix", "end-repair and dA-tailing enzymes",
                      MaterialKind.CATALYST)
    tailed = g.product("tailed", "end-repaired dA-tailed cDNA fragments",
                       MaterialKind.MOLECULAR)
    g.react("end_prep", "End-repair and dA-tail cDNA fragments", StepKind.END_REPAIR,
            ((fragments, InputRole.TEMPLATE), (end_mix, InputRole.CATALYST)), (tailed,))
    adapter = g.input("read2_adapter", "3-prime-dT Read-2 adapter",
                      MaterialKind.ADAPTER, inferred=inferred_intermediates,
                      source=authority + (" (adapter reconstructed from final library)"
                                          if inferred_intermediates else ""))
    ligase = g.input("ligase", "adapter ligase", MaterialKind.CATALYST)
    adapted = g.product("adapted", "Read-2-adapted barcode-bearing fragments",
                        MaterialKind.MOLECULAR)
    g.react("ligate", "Ligate the Read-2 adapter", StepKind.LIGATION,
            ((tailed, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
             (ligase, InputRole.CATALYST)), (adapted,), geometry=READ2_ADAPTER_END)

    index_f = g.input("index_forward", "sample-index PCR forward primer",
                      MaterialKind.OLIGO, inferred=True,
                      source=authority + " (commercial primer bases not asserted here)")
    index_r = g.input("index_reverse", "sample-index PCR reverse primer",
                      MaterialKind.OLIGO, inferred=True,
                      source=authority + " (commercial primer bases not asserted here)")
    index_pol = g.input("index_polymerase", "sample-index PCR polymerase",
                        MaterialKind.CATALYST)
    library = g.product("library", "indexed Chromium gene-expression library",
                        MaterialKind.LIBRARY, state=final)
    g.react("index_pcr", "Select the barcode-bearing end and add sample indexes",
            StepKind.PCR,
            ((adapted, InputRole.TEMPLATE), (index_f, InputRole.STRUCTURAL),
             (index_r, InputRole.STRUCTURAL), (index_pol, InputRole.CATALYST)),
            (library,), geometry=PCRGeometry.unavailable(
                final, reason="The page authority does not print both sample-index PCR primers."))
    g.sequence(library, SequencingPlan(
        library=final, primers=tuple(primers),
        acquired_roles=tuple(p.role for p in primers), read_lengths=read_lengths))
    reviewed.extend(("rt_primer", "rt", "reverse_transcribe", "tso",
                     "template_switch", "gem_recovery", "recover", "preamp_forward",
                     "preamp_reverse", "preamp_polymerase", "preamp", "fragmentase",
                     "fragment", "end_mix", "end_prep", "read2_adapter", "ligase",
                     "ligate", "index_forward", "index_reverse", "index_polymerase",
                     "index_pcr", "library"))
    g.declare_reviewed_scope(reviewed, authority)
    return g.require_publishable()
