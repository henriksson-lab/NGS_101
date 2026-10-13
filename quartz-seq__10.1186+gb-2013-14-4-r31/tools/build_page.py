#!/usr/bin/env python3
"""Build the diagram-centric Quartz-Seq chemistry page."""
from __future__ import annotations

import html
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import quartz_seq as Q
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Segment, oligo, protocol_graph_panels
from endprep import dA_tailed_scene, repair_and_dA_tail
from page import head, info, table
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from rna_amplification import poly_a_tagging_amplification

OUT = HERE.parent / "quartz-seq.html"


def seg(name, top, tag=None, **kw):
    return Segment(name=name, top=top, tag=tag, **kw)


def oligos() -> str:
    rows = [
        oligo("RT primer", Q.rt_primer()),
        oligo("Tagging primer", Q.tagging_primer()),
        oligo("Suppression primer", [seg("G", "G"), seg("M handle", Q.M, "tso")],
              mods="/5Amine/"),
        oligo("TRSU adapter strand", [seg("P5 + Read 1", Q.TRSU, "p5+r1")],
              three="-3' (terminal phosphorothioate)"),
        oligo("TRSI-2 adapter strand", [seg("Index 1 primer", Q.TRSI2[:33], "r2"),
                                         seg("i7", Q.TRSI2_INDEX, "cbc"),
                                         seg("P7 reverse complement", Q.TRSI2[-24:], "p7")],
              mods="/5Phos/", three="-3' (terminal phosphorothioate)"),
    ]
    return "<h2>Key oligos</h2>\n<seq>\n" + "\n".join(rows) + "\n</seq>"


def protocol_graph() -> ProtocolGraph:
    source = "https://doi.org/10.1186/gb-2013-14-4-r31"
    graph = ProtocolGraph("Quartz-Seq", authority=source)
    mrna = graph.input("mrna", "polyadenylated RNA from one cell",
                       MaterialKind.BIOLOGICAL, state=Q.mrna_rt_scene(), primary=True)
    rt_primer = graph.input("rt_primer", "Quartz-Seq RT primer",
                            MaterialKind.OLIGO, state=Construct(Q.rt_primer()))
    reverse_transcriptase = graph.input("reverse_transcriptase",
                                       "SuperScript III reverse transcriptase",
                                       MaterialKind.CATALYST)
    first = graph.product("first_strand", "first-strand cDNA",
                          MaterialKind.MOLECULAR, state=Q.first_strand())
    graph.react("reverse_transcription", "Reverse-transcribe polyadenylated RNA",
                StepKind.REVERSE_TRANSCRIPTION,
                ((mrna, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
                 (reverse_transcriptase, InputRole.CATALYST)), (first,))

    exonuclease = graph.input("exonuclease_i", "exonuclease I",
                              MaterialKind.CATALYST)
    cleaned = graph.product("primer_cleaned_cdna", "first strand after primer removal",
                            MaterialKind.MOLECULAR, state=Q.first_strand())
    primer_waste = graph.product("primer_waste", "digested unused RT primer",
                                 MaterialKind.WASTE)
    graph.react("remove_unused_primer", "Digest unused RT primer",
                StepKind.CLEANUP,
                ((first, InputRole.TEMPLATE),
                 (exonuclease, InputRole.CATALYST)), (cleaned, primer_waste))
    graph.discard(primer_waste, "unused RT primer is digested before tailing")

    rnase_h = graph.input("rnase_h", "RNase H", MaterialKind.CATALYST)
    rna_removed = graph.product("rna_removed_cdna", "RNA-depleted first-strand cDNA",
                                MaterialKind.MOLECULAR, state=Q.first_strand())
    graph.react("remove_rna", "Remove the RNA template with RNase H",
                StepKind.DIGESTION,
                ((cleaned, InputRole.TEMPLATE),
                 (rnase_h, InputRole.CATALYST)), (rna_removed,))
    amplified, shared_reviewed = poly_a_tagging_amplification(
        graph, rna_removed, tailed_state=Q.tailed_first_strand(),
        tagging_primer_label="Quartz-Seq tagging primer",
        tagging_primer_state=Construct(Q.tagging_primer()),
        duplex_state=Q.wta_amplicon(),
        suppression_primer_label="Quartz-Seq suppression primer",
        suppression_primer_sequence=Q.SUPPRESSION_PRIMER,
        amplified_state=Q.wta_amplicon())

    shear = graph.input("shearing_condition", "physical shearing condition",
                        MaterialKind.CONDITION)
    fragment = Construct([seg("sheared cDNA fragment", "X" * 36,
                              placeholder=True)])
    sheared = graph.product("sheared_cdna", "sheared amplified cDNA",
                            MaterialKind.MOLECULAR, state=fragment)
    graph.react("shear", "Shear amplified cDNA", StepKind.FRAGMENTATION,
                ((amplified, InputRole.TEMPLATE),
                 (shear, InputRole.CONDITION)), (sheared,))
    end_repair = graph.input("end_repair_enzymes", "end-repair enzyme mix",
                             MaterialKind.CATALYST)
    repaired = graph.product("end_repaired_cdna", "end-repaired cDNA fragments",
                             MaterialKind.MOLECULAR, state=fragment)
    graph.react("end_repair", "End-repair the sheared fragments",
                StepKind.END_REPAIR,
                ((sheared, InputRole.TEMPLATE),
                 (end_repair, InputRole.CATALYST)), (repaired,))
    da_enzyme = graph.input("da_tailing_enzyme", "dA-tailing polymerase",
                            MaterialKind.CATALYST)
    da_product = repair_and_dA_tail(fragment)
    da_tailed = graph.product("da_tailed_cdna", "3′-dA-tailed cDNA fragments",
                              MaterialKind.MOLECULAR,
                              state=dA_tailed_scene(da_product))
    graph.react("da_tail", "Add 3′ dA overhangs", StepKind.TAILING,
                ((repaired, InputRole.TEMPLATE),
                 (da_enzyme, InputRole.CATALYST)), (da_tailed,))

    adapter = graph.input("truseq_adapter", "TRSU/TRSI-2 forked adapter",
                          MaterialKind.ADAPTER, state=Q.adapter_scene())
    ligase = graph.input("dna_ligase", "DNA ligase", MaterialKind.CATALYST)
    ligated = graph.product("ligated_library", "adapter-ligated library molecules",
                            MaterialKind.MOLECULAR, state=Q.final_library())
    ligation_geometry = AdapterLigationGeometry.known(
        "cohesive", dA_tailed_scene(da_product), Q.adapter_scene(),
        statement="The source-defined adapter 3′-T pairs with each substrate 3′-dA.")
    graph.react("ligate_adapter", "Ligate the single-index TruSeq adapter",
                StepKind.LIGATION,
                ((da_tailed, InputRole.TEMPLATE),
                 (adapter, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (ligated,),
                geometry=ligation_geometry)

    beads = graph.input("size_selection_beads", "AMPure XP size-selection beads",
                        MaterialKind.SELECTION)
    selected = graph.product("selected_ligation_products",
                             "size-selected adapter-ligated molecules",
                             MaterialKind.MOLECULAR, state=Q.final_library())
    adapter_waste = graph.product("adapter_dimer_waste",
                                  "adapter dimers and rejected fragments",
                                  MaterialKind.WASTE)
    graph.react("size_select", "Remove adapter dimers by two bead-selection rounds",
                StepKind.SELECTION,
                ((ligated, InputRole.TEMPLATE),
                 (beads, InputRole.SELECTION)), (selected, adapter_waste))
    graph.discard(adapter_waste, "rejected by the high-salt AMPure size selection")

    tpc1 = PCRPrimer("TPC1", Q.TPC1)
    tpc2 = PCRPrimer("TPC2", Q.TPC2)
    forward = graph.input("tpc1", "TPC1 library-PCR primer", MaterialKind.OLIGO,
                          state=Construct([seg("TPC1", Q.TPC1)]))
    reverse = graph.input("tpc2", "TPC2 library-PCR primer", MaterialKind.OLIGO,
                          state=Construct([seg("TPC2", Q.TPC2)]))
    library_polymerase = graph.input("library_polymerase", "KAPA HiFi polymerase",
                                     MaterialKind.CATALYST)
    final = Q.final_library()
    library = graph.product("library", "Quartz-Seq sequencing library",
                            MaterialKind.LIBRARY, state=final)
    graph.react("library_pcr", "Amplify the adapter-ligated library with TPC1/TPC2",
                StepKind.PCR,
                ((selected, InputRole.TEMPLATE),
                 (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL),
                 (library_polymerase, InputRole.CATALYST)), (library,),
                geometry=PCRGeometry.exact(
                    final, tpc1, tpc2,
                    statement="Table 1 prints TPC1 and TPC2."))
    graph.sequence(library, SequencingPlan(
        library=final, primers=Q.SEQ_PRIMERS,
        acquired_roles=tuple(p.role for p in Q.SEQ_PRIMERS),
        read_lengths={"Read 1": 50, "Index 1 (i7)": 6, "Read 2": 50}))
    reviewed = [
        "mrna", "rt_primer", "reverse_transcriptase", "reverse_transcription",
        "exonuclease_i", "remove_unused_primer", "rnase_h", "remove_rna",
        *shared_reviewed, "shearing_condition", "shear", "end_repair_enzymes",
        "end_repair", "da_tailing_enzyme", "da_tail", "truseq_adapter",
        "dna_ligase", "ligate_adapter", "size_selection_beads", "size_select",
        "tpc1", "tpc2", "library_polymerase", "library_pcr", "library",
    ]
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def sequencing() -> str:
    lib = Q.final_library()
    rows = []
    for primer, hit in Q.primer_landings():
        rows.append((html.escape(primer.role), html.escape(primer.name),
                     html.escape(", ".join(hit.covers)),
                     f"<code>{html.escape(hit.reads)}</code>&hellip;"))
    return f'''<h2>Sequencing</h2>
{info('Paired-end 50-base reads with a 6-base Index 1. Cell identity is the library tube/index; there is no molecular cell barcode or UMI.')}
{sp.diagram(lib, Q.SEQ_PRIMERS)}
{table(("Read", "Primer", "Primer site", "First bases"), rows)}
{sp.read_content_table(lib, Q.SEQ_PRIMERS,
    {"Read 1":50,"Index 1 (i7)":6,"Read 2":50},
    acquired_roles=tuple(p.role for p in Q.SEQ_PRIMERS))}
'''


def main() -> None:
    graph = protocol_graph()
    construction = "<h2>Library construction</h2>\n" + protocol_graph_panels(graph)
    page = "\n".join([head("Quartz-Seq library chemistry"),
        '<div class="wrap" data-protocol-complete="true">',
        '<h1>Quartz-Seq &mdash; poly(A)-tagging whole-transcript amplification</h1>',
        info('Defining source: <a href="https://doi.org/10.1186/gb-2013-14-4-r31">Sasagawa et al., <i>Genome Biology</i> (2013)</a>.'),
        oligos(), construction, sequencing(), '</div>'])
    OUT.write_text(page, encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
