"""Construct model for 10x Chromium Single Cell 3' Gene Expression v2 (CG000108)."""
from __future__ import annotations

import illumina as il
import rt
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, complement_segments, feature, revcomp
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2]))
from protocol_families.tenx_graph import ChromiumRNAOps, READ2_ADAPTER_END
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import InputRole, MaterialKind, StepKind
from protocol_graph import ProtocolGraph

PARTIAL_R1 = il.TRUSEQ_READ1[11:]
TSO_DNA = rt.SMART_HANDLE + "ACAT"
CELL_BARCODE = feature("cell_barcode", "cell_barcode", "whitelist",
                       whitelist="10x-chromium-3prime-v2")
UMI = feature("umi", "umi", "random")
SAMPLE_INDEX = feature("sample_index_i7", "sample_index", "unknown")
FEATURES = {"cell barcode": CELL_BARCODE, "UMI": UMI,
            "i7 reverse complement": SAMPLE_INDEX}


def seg(name, top, tag=None, **kw):
    kw.setdefault("feature", FEATURES.get(name))
    return Segment(name=name, top=top, tag=tag, **kw)


def bead_oligo():
    return [seg("Partial Read 1", PARTIAL_R1, "r1"),
            seg("cell barcode", "B" * 16, "cbc", placeholder=True),
            seg("UMI", "U" * 10, "umi", placeholder=True),
            seg("poly(dT)30", "T" * 30), seg("V", "V", placeholder=True),
            seg("N", "N", placeholder=True)]


def tso():
    return [seg("SMART handle", rt.SMART_HANDLE, "tso"), seg("ACAT", "ACAT", "tso"),
            seg("rGrGrG", "GGG", "tso")]


def capture_scene():
    mrna = [seg("mRNA", "X" * 32, placeholder=True), seg("poly(A)", "A" * 30)]
    sc = Scene(); sc.strand("mRNA", mrna, label="mRNA")
    sc.anneal("bead primer", bead_oligo(), to="mRNA", pair=("poly(dT)30", "poly(A)"),
              label="released gel-bead primer", unpaired=("V", "N"))
    sc.arrow("bead primer", "reverse transcriptase")
    return sc


def first_strand():
    return Construct([*bead_oligo(), seg("antisense cDNA", "X" * 34, placeholder=True),
                      seg("CCC", rt.UNTEMPLATED_TAIL, "tso")], name="v2 first strand")


def switch_scene():
    first = first_strand(); sc = Scene(); sc.strand("first", list(first), label="first strand")
    sc.anneal("TSO", tso(), to="first", pair=("rGrGrG", "CCC"), label="TSO", above=True)
    sc.arrow("first", "RT switches template")
    return sc


def amplified_cdna():
    return Construct([*bead_oligo(), seg("antisense cDNA", "X" * 34, placeholder=True),
                      seg("CCC", rt.UNTEMPLATED_TAIL, "tso"),
                      seg("copied TSO", revcomp(TSO_DNA), "tso")],
                     name="v2 amplified cDNA")


def adapter_top():
    return [seg("Read 2 adapter", il.INDEX1_PRIMER, "r2")]


def adapter_bottom():
    return [seg("stem complement", il.STEM_COMPLEMENT, "r2"), seg("3' T", "T", "r2")]


def adapter_scene():
    sc = Scene(); sc.strand("top", adapter_top(), label="Read 2 adapter")
    sc.anneal("bottom", adapter_bottom(), to="top", pair=("stem complement", "Read 2 adapter"),
              label="short strand")
    sc.mark("bottom", "3' T", "ligation overhang")
    return sc


SEQ_PRIMERS = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["R2"])


def final_library(insert_nt=36):
    lib = Construct([seg("P5", il.P5, "p5"), seg("Read 1 arm", il.TRUSEQ_READ1[4:], "r1"),
        seg("cell barcode", "B" * 16, "cbc", placeholder=True),
        seg("UMI", "U" * 10, "umi", placeholder=True), seg("poly(dT)30", "T" * 30),
        seg("V", "V", placeholder=True), seg("N", "N", placeholder=True),
        seg("cDNA insert", "X" * insert_nt, placeholder=True), seg("dA junction", "A"),
        seg("Index 1 / Read 2 arm", il.INDEX1_PRIMER, "r2"),
        seg("i7 reverse complement", "I" * 8, "cbc", placeholder=True),
        seg("P7 reverse complement", il.P7_RC, "p7")], name="10x 3' v2 library")
    problems = sp.verify(lib, SEQ_PRIMERS,
                         required_roles=("Read 1", "Index 1 (i7)", "Read 2"))
    if problems: raise ValueError("invalid v2 library: " + "; ".join(problems))
    return lib


def read1_layout():
    return [(x.cycles, x.feature.label) for x in
            sp.feature_spans(final_library(), SEQ_PRIMERS, {"Read 1": 26})]



def switched_cdna():
    return Construct(list(amplified_cdna()), name="v2 template-switched cDNA")


def cdna_primers():
    return (PCRPrimer("cDNA forward primer", PARTIAL_R1),
            PCRPrimer("cDNA reverse primer", TSO_DNA))


def index_primers():
    # Source prints an N8 sample-index window; it represents one oligo in a well mix.
    return (PCRPrimer("sample-index forward primer", il.NEBNEXT_UNIVERSAL_PRIMER[:49],
                      anneal_nt=13),
            PCRPrimer("sample-index reverse primer",
                      il.P7+"N"*8+il.TRUSEQ_READ2[:21]))


def pre_index_library(insert_nt=36):
    return Construct([*bead_oligo(),
                      seg("cDNA insert", "X"*insert_nt, placeholder=True),
                      seg("dA junction", "A"), *adapter_top()],
                     name="v2 adapter-ligated cDNA fragment")


def _source_pcr(ops, prefix, label, template, product_state, primers, *, library=False):
    """The source prints both primers; locate their sites rather than claiming absence."""
    g=ops.graph
    forward, reverse=primers
    f=g.input(prefix+"_forward", forward.name, MaterialKind.OLIGO,
              state=Construct([Segment(forward.name, forward.sequence)]))
    r=g.input(prefix+"_reverse", reverse.name, MaterialKind.OLIGO,
              state=Construct([Segment(reverse.name, reverse.sequence)]))
    polymerase=g.input(prefix+"_polymerase", label+" polymerase", MaterialKind.CATALYST)
    product=g.product("library" if library else "amplified", label+" product",
                      MaterialKind.LIBRARY if library else MaterialKind.MOLECULAR,
                      state=product_state)
    g.react(prefix, label, StepKind.PCR,
            ((template,InputRole.TEMPLATE),(f,InputRole.STRUCTURAL),
             (r,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)), (product,),
            geometry=PCRGeometry.exact(template.state,forward,reverse,
                statement="Both primer architectures are printed in CG000108 Rev A Fig. 4."))
    ops.reviewed.extend((f.id,r.id,polymerase.id,prefix,product.id))
    return product


def _source_ligation(ops, tailed):
    g=ops.graph
    adapter=g.input("read2_adapter", "3-prime-dT Read-2 adapter", MaterialKind.ADAPTER,
                    state=adapter_scene())
    ligase=g.input("ligase", "adapter ligase", MaterialKind.CATALYST)
    adapted=g.product("adapted", "Read-2-adapted barcode-bearing fragment",
                      MaterialKind.MOLECULAR,state=pre_index_library())
    g.react("ligate", "Ligate the Read-2 adapter", StepKind.LIGATION,
            ((tailed,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),
             (ligase,InputRole.CATALYST)),(adapted,),geometry=READ2_ADAPTER_END)
    ops.reviewed.extend((adapter.id,ligase.id,"ligate"))
    return adapted

def protocol_graph():
    authority = ('Primary chemistry source: <a href="https://www.10xgenomics.com/support">'
                 '10x Genomics Technical Note CG000108 Rev A</a>.')
    graph = ProtocolGraph("10x Chromium Single Cell 3-prime Gene Expression v2", authority=authority)
    ops = ChromiumRNAOps(graph, inferred_intermediates=False)
    first_state = first_strand()
    final = final_library()
    cells = ops.cells()
    occupied = ops.partition(cells)
    first, rt_enzyme = ops.reverse_transcribe(
        occupied, first_state, primer_label="cell-barcode/UMI oligo(dT) gel-bead primer")
    switched = ops.template_switch(
        first, rt_enzyme, switched_cdna(), tso_label="template-switch oligo")
    recovered = ops.recover(switched, switched_cdna())
    amplified = _source_pcr(ops, "preamp", "Amplify recovered full-length cDNA",
                            recovered, amplified_cdna(), cdna_primers())
    fragments = ops.fragment(amplified)
    tailed = ops.end_prepare(fragments)
    adapted = _source_ligation(ops, tailed)
    library = _source_pcr(ops, "index_pcr", "Select barcode-bearing end and add sample index",
                          adapted, final, index_primers(), library=True)
    ops.sequence(library, final, SEQ_PRIMERS, {"Read 1": 26, "Index 1 (i7)": 8, "Read 2": 98})
    return graph
