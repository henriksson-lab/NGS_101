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
                 limitation="Exact sequences for both PCR primers have not been declared in this model."):
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


class ChromiumRNAOps:
    """Bounded Chromium reactions; each protocol module chooses their order."""

    def __init__(self, graph: ProtocolGraph, *, inferred_intermediates: bool = False):
        self.graph = graph
        self.authority = graph.authority
        self.inferred = inferred_intermediates
        self.reviewed: list[str] = []

    def cells(self):
        cells = self.graph.input("cells", "cells containing poly(A) RNA",
                                 MaterialKind.BIOLOGICAL,
                                 state=molecule("cells", (role("poly(A) RNA"),)),
                                 primary=True)
        self.reviewed.append("cells")
        return cells

    def partition(self, cells):
        return partition_material(self.graph, cells, self.reviewed, label="cells")

    def reverse_transcribe(self, occupied, first_strand, *, primer_label):
        g, authority = self.graph, self.authority
        primer = g.input("rt_primer", primer_label, MaterialKind.OLIGO,
                         inferred=self.inferred,
                         source=authority + (" (intermediate bases inferred from the final construct)"
                                             if self.inferred else ""))
        rt = g.input("rt", "template-switching reverse transcriptase",
                     MaterialKind.CATALYST)
        first = g.product("first_strand", "first-strand transcript cDNA",
                          MaterialKind.MOLECULAR, state=first_strand,
                          inferred=self.inferred,
                          source=authority + (" (intermediate structure inferred)"
                                              if self.inferred else ""))
        g.react("reverse_transcribe", "Capture and reverse-transcribe poly(A) RNA",
                StepKind.REVERSE_TRANSCRIPTION,
                ((occupied, InputRole.TEMPLATE), (primer, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (first,))
        self.reviewed.extend(("rt_primer", "rt", "reverse_transcribe"))
        return first, rt

    def template_switch(self, first, rt, first_strand, *, tso_label):
        g, authority = self.graph, self.authority
        tso = g.input("tso", tso_label, MaterialKind.OLIGO,
                      inferred=self.inferred,
                      source=authority + (" (intermediate bases inferred from the final construct)"
                                          if self.inferred else ""))
        switched = g.product("switched", "template-switched cell-barcoded cDNA",
                             MaterialKind.MOLECULAR, state=first_strand,
                             inferred=self.inferred,
                             source=authority + (" (intermediate structure inferred)"
                                                 if self.inferred else ""))
        g.react("template_switch", "Install the opposite cDNA amplification handle",
                StepKind.TEMPLATE_SWITCHING,
                ((first, InputRole.TEMPLATE), (tso, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (switched,))
        self.reviewed.extend(("tso", "template_switch"))
        return switched

    def recover(self, switched, first_strand):
        g = self.graph
        recovery = g.input("gem_recovery", "GEM break and cDNA recovery condition",
                           MaterialKind.CONDITION)
        recovered = g.product("recovered", "recovered cell-barcoded cDNA",
                              MaterialKind.MOLECULAR, state=first_strand)
        g.react("recover", "Break GEMs and recover cDNA", StepKind.EXTRACTION,
                ((switched, InputRole.TEMPLATE), (recovery, InputRole.CONDITION)),
                (recovered,), changes_molecule=False)
        self.reviewed.extend(("gem_recovery", "recover"))
        return recovered

    def preamplify(self, recovered, first_strand, amplified_cdna):
        g, authority = self.graph, self.authority
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
                              inferred=self.inferred,
                              source=authority + (" (intermediate structure inferred)"
                                                  if self.inferred else ""))
        g.react("preamp", "Amplify recovered full-length cDNA", StepKind.PCR,
                ((recovered, InputRole.TEMPLATE), (pre_f, InputRole.STRUCTURAL),
                 (pre_r, InputRole.STRUCTURAL), (pre_pol, InputRole.CATALYST)),
                (amplified,), geometry=PCRGeometry.unavailable(
                    first_strand, reason="Exact cDNA preamplification primer sequences have not been declared in this model."))
        self.reviewed.extend(("preamp_forward", "preamp_reverse", "preamp_polymerase", "preamp"))
        return amplified

    def fragment(self, amplified):
        g = self.graph
        fragmentase = g.input("fragmentase", "cDNA fragmentation enzyme",
                              MaterialKind.CATALYST)
        fragments = g.product("fragments", "fragmented barcode-bearing cDNA",
                              MaterialKind.MOLECULAR)
        g.react("fragment", "Fragment amplified cDNA", StepKind.FRAGMENTATION,
                ((amplified, InputRole.TEMPLATE), (fragmentase, InputRole.CATALYST)),
                (fragments,))
        self.reviewed.extend(("fragmentase", "fragment"))
        return fragments

    def end_prepare(self, fragments):
        g = self.graph
        end_mix = g.input("end_mix", "end-repair and dA-tailing enzymes",
                          MaterialKind.CATALYST)
        tailed = g.product("tailed", "end-repaired dA-tailed cDNA fragments",
                           MaterialKind.MOLECULAR)
        g.react("end_prep", "End-repair and dA-tail cDNA fragments", StepKind.END_REPAIR,
                ((fragments, InputRole.TEMPLATE), (end_mix, InputRole.CATALYST)), (tailed,))
        self.reviewed.extend(("end_mix", "end_prep"))
        return tailed

    def ligate_read2(self, tailed):
        g, authority = self.graph, self.authority
        adapter = g.input("read2_adapter", "3-prime-dT Read-2 adapter",
                          MaterialKind.ADAPTER, inferred=self.inferred,
                          source=authority + (" (adapter reconstructed from final library)"
                                              if self.inferred else ""))
        ligase = g.input("ligase", "adapter ligase", MaterialKind.CATALYST)
        adapted = g.product("adapted", "Read-2-adapted barcode-bearing fragments",
                            MaterialKind.MOLECULAR)
        g.react("ligate", "Ligate the Read-2 adapter", StepKind.LIGATION,
                ((tailed, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (adapted,), geometry=READ2_ADAPTER_END)
        self.reviewed.extend(("read2_adapter", "ligase", "ligate"))
        return adapted

    def index_pcr(self, adapted, final):
        g, authority = self.graph, self.authority
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
                    final, reason="Exact sample-index PCR primer sequences have not been declared in this model."))
        self.reviewed.extend(("index_forward", "index_reverse", "index_polymerase",
                              "index_pcr", "library"))
        return library

    def sequence(self, library, final, primers, read_lengths):
        self.graph.sequence(library, SequencingPlan(
            library=final, primers=tuple(primers),
            acquired_roles=tuple(p.role for p in primers), read_lengths=read_lengths))
        self.graph.declare_reviewed_scope(self.reviewed, self.authority)
        return self.graph.require_publishable()
