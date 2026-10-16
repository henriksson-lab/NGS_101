#!/usr/bin/env python3
"""Build the diagram-centric Quartz-Seq2 chemistry page."""
from __future__ import annotations

import html
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import quartz_seq2 as Q
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Segment, oligo, protocol_graph_panels
from endprep import dA_tailed_scene, repair_and_dA_tail
from page import head, info, table
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from rna_amplification import poly_a_tagging_amplification

OUT = HERE.parent / "quartz-seq2.html"


def seg(name, top, tag=None, **kw):
    return Segment(name=name, top=top, tag=tag, **kw)


def oligos() -> str:
    rows = [
        oligo("1536-well RT primer", Q.rt_primer()),
        oligo("Tagging primer", Q.tagging_primer()),
        oligo("gM suppression primer", [seg("G", "G"), seg("M handle", Q.M, "tso")]),
        oligo("rYshapeP5 short strand", [seg("stem", Q.RYSHAPE_P5[:12], "r1"),
                                          seg("unpaired arm", Q.RYSHAPE_P5[12:], "r1")]),
        oligo("rYshapeP7LT06 long strand",
              [seg("P7", Q.RYSHAPE_P7_LT06[:24], "p7"),
               seg("i7", Q.LT06_OLIGO_INDEX, "cbc"),
               seg("TruSeq Read 2", Q.RYSHAPE_P7_LT06[-34:], "r2")]),
        oligo("P5-gMac selective PCR primer",
              [seg("truncated P5", Q.P5_GMAC[:28], "p5"), seg("TT", "TT"),
               seg("gM + AC", Q.P5_GMAC[30:], "tso")],
              mods="two 3'-terminal phosphorothioates"),
        oligo("Read1DropQuartz sequencing primer", [seg("custom R1", Q.READ1_DROPQUARTZ, "r1")]),
    ]
    return "<h2>Key oligos</h2>\n<seq>\n" + "\n".join(rows) + "\n</seq>"


def protocol_graph() -> ProtocolGraph:
    source = "https://doi.org/10.1186/s13059-018-1407-3"
    graph = ProtocolGraph("Quartz-Seq2", authority=source)
    mrna = graph.input("mrna", "polyadenylated RNA from one well",
                       MaterialKind.BIOLOGICAL, state=Q.mrna_rt_scene(), primary=True)
    rt_primer = graph.input("rt_primer", "well-barcoded Quartz-Seq2 RT primer",
                            MaterialKind.OLIGO, state=Construct(Q.rt_primer()))
    reverse_transcriptase = graph.input("reverse_transcriptase",
                                       "reverse transcriptase", MaterialKind.CATALYST)
    first = graph.product("first_strand", "well-barcoded first-strand cDNA",
                          MaterialKind.MOLECULAR, state=Q.first_strand())
    graph.react("reverse_transcription",
                "Reverse-transcribe with the well-barcoded RT primer",
                StepKind.REVERSE_TRANSCRIPTION,
                ((mrna, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
                 (reverse_transcriptase, InputRole.CATALYST)), (first,))

    peers = graph.input("peer_wells", "barcoded first-strand cDNA from other wells",
                        MaterialKind.MOLECULAR, state=Q.first_strand())
    pooled = graph.product("pooled_first_strand", "pooled barcoded first-strand cDNA",
                           MaterialKind.MOLECULAR, state=Q.first_strand())
    graph.react("pool", "Pool wells after molecular cell indexing", StepKind.POOL,
                ((first, InputRole.TEMPLATE), (peers, InputRole.TEMPLATE)),
                (pooled,), changes_molecule=False)
    column = graph.input("purification_column", "silica purification column",
                         MaterialKind.SELECTION)
    cleaned = graph.product("cleaned_first_strand", "column-purified first-strand cDNA",
                            MaterialKind.MOLECULAR, state=Q.first_strand())
    cleanup_waste = graph.product("cleanup_waste", "unretained primers and small products",
                                  MaterialKind.WASTE)
    graph.react("column_cleanup", "Column-purify the pooled cDNA",
                StepKind.CLEANUP,
                ((pooled, InputRole.TEMPLATE), (column, InputRole.SELECTION)),
                (cleaned, cleanup_waste))
    graph.discard(cleanup_waste, "not retained by the post-pooling column cleanup")

    rnase_h = graph.input("rnase_h", "RNase H", MaterialKind.CATALYST)
    rna_removed = graph.product("rna_removed_cdna", "RNA-depleted first-strand cDNA",
                                MaterialKind.MOLECULAR, state=Q.first_strand())
    graph.react("remove_rna", "Remove the RNA template with RNase H",
                StepKind.DIGESTION,
                ((cleaned, InputRole.TEMPLATE),
                 (rnase_h, InputRole.CATALYST)), (rna_removed,))
    amplified, shared_reviewed = poly_a_tagging_amplification(
        graph, rna_removed, tailed_state=Q.tailed_first_strand(),
        tagging_primer_label="Quartz-Seq2 tagging primer",
        tagging_primer_state=Construct(Q.tagging_primer()),
        duplex_state=Q.wta_amplicon(), suppression_primer_label="gM primer",
        suppression_primer_sequence=Q.GM_PRIMER,
        amplified_state=Q.wta_amplicon())

    shear = graph.input("shearing_condition", "physical shearing condition",
                        MaterialKind.CONDITION)
    fragment = Construct([seg("3-prime cDNA fragment", "X" * 36,
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

    adapter = graph.input("truncated_y_adapter", "truncated rYshape adapter",
                          MaterialKind.ADAPTER, state=Q.truncated_adapter_scene())
    ligase = graph.input("dna_ligase", "DNA ligase", MaterialKind.CATALYST)
    selective_template = Q.selective_pcr_template()
    ligated = graph.product("ligated_template", "adapter-ligated 3′ cDNA fragment",
                            MaterialKind.MOLECULAR, state=selective_template)
    geometry = AdapterLigationGeometry.known(
        "cohesive", dA_tailed_scene(da_product), Q.truncated_adapter_scene(),
        statement="The long adapter strand exposes the source-defined 3′-T overhang.")
    graph.react("ligate_adapter", "Ligate the truncated Y adapter",
                StepKind.LIGATION,
                ((da_tailed, InputRole.TEMPLATE),
                 (adapter, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (ligated,), geometry=geometry)

    p5_primer = PCRPrimer("P5-gMac", Q.P5_GMAC)
    tpc2_primer = PCRPrimer("TPC2", Q.TPC2)
    p5 = graph.input("p5_gmac", "P5-gMac selective PCR primer",
                     MaterialKind.OLIGO,
                     state=Construct([seg("P5-gMac", Q.P5_GMAC)]))
    tpc2 = graph.input("tpc2", "TPC2 library-PCR primer", MaterialKind.OLIGO,
                       state=Construct([seg("TPC2", Q.TPC2)]))
    polymerase = graph.input("library_polymerase", "KAPA HiFi polymerase",
                             MaterialKind.CATALYST)
    final = Q.final_library()
    library = graph.product("library", "Quartz-Seq2 sequencing library",
                            MaterialKind.LIBRARY, state=final)
    graph.react("selective_pcr",
                "Select the RT-primer end and complete the sequencing library",
                StepKind.PCR,
                ((ligated, InputRole.TEMPLATE), (p5, InputRole.STRUCTURAL),
                 (tpc2, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (library,),
                geometry=PCRGeometry.exact(
                    final, p5_primer, tpc2_primer,
                    statement="Supplementary Table S4 prints P5-gMac and TPC2."))
    graph.sequence(library, SequencingPlan(
        library=final, primers=Q.SEQ_PRIMERS,
        acquired_roles=tuple(p.role for p in Q.SEQ_PRIMERS),
        read_lengths={"Read 1": 23, "Index 1 (i7)": 6, "Read 2": 63}))
    reviewed = [
        "mrna", "rt_primer", "reverse_transcriptase", "reverse_transcription",
        "peer_wells", "pool", "purification_column", "column_cleanup", "rnase_h",
        "remove_rna", *shared_reviewed, "shearing_condition", "shear",
        "end_repair_enzymes", "end_repair", "da_tailing_enzyme", "da_tail",
        "truncated_y_adapter", "dna_ligase", "ligate_adapter", "p5_gmac", "tpc2",
        "library_polymerase", "selective_pcr", "library",
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
{info('Quartz-Seq2 v3.2 uses Read 1: 23 cycles, Index 1: 6, and Read 2: 63. Read 1 reports the 15-base cell barcode followed by the 8-base UMI; Read 2 reports cDNA.')}
{sp.diagram(lib, Q.SEQ_PRIMERS)}
{table(("Read", "Primer", "Primer site", "First bases"), rows)}
<h3>Read 1 barcode layout</h3>
{table(("Cycles", "Content"), Q.read1_layout(), scroll=False)}
{sp.read_content_table(lib, Q.SEQ_PRIMERS,
    {"Read 1":23,"Index 1 (i7)":6,"Read 2":63},
    acquired_roles=tuple(p.role for p in Q.SEQ_PRIMERS))}
'''


def main() -> None:
    graph = protocol_graph()
    construction = "<h2>Library construction</h2>\n" + protocol_graph_panels(graph)
    page = "\n".join([head("Quartz-Seq2 library chemistry"),
        '<div class="wrap" data-protocol-complete="true">',
        '<h1>Quartz-Seq2 &mdash; barcoded 3-prime single-cell RNA sequencing</h1>',
        info('Defining source: <a href="https://doi.org/10.1186/s13059-018-1407-3">Sasagawa et al., <i>Genome Biology</i> (2018)</a>.'),
        info('Published-source discrepancy: Additional file 5, Table S4, v3.2 sheet '
             '(MOESM5 extracted lines 389–393) labels its heading “SeqLv15mer” but the '
             'next line calls the same v3.2 barcodes “14 mer.” The printed eMDRT0001 '
             'barcode <code>ATCAATCCTATCTGC</code> has 15 bases, as do the v3.2 '
             'barcodes inspected in that sheet; the article Methods (cell-barcode design) '
             'also specifies 15 nt for v3.2. This model uses 15 nt. The conflicting '
             '“14 mer” line appears to be a table heading error; no correction notice '
             'was found in the inspected sources.'),
        oligos(), construction, sequencing(), '</div>'])
    OUT.write_text(page, encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
