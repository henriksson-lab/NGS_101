"""Shared construct builder for vendor-defined 10x 5' V(D)J chemistries."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import illumina as il
import rt
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, feature, revcomp
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import molecule, role
from tenx_graph import READ2_ADAPTER_END, chromium_pcr, partition_material

PARTIAL_R1 = il.TRUSEQ_READ1[11:]
TSO = "TTTCTTATATGGG"
VDJ_FORWARD = "GATCTACACTCTTTCCCTACACGACGC"
CDNA_FORWARD = "CTACACGACGCTCTTCCGATCT"
CDNA_REVERSE = "AAGCAGTGGTATCAACGCAGAG"
POLY_DT_NT = 30
INDEX_NT = 10
CELL_BARCODE = feature("cell_barcode", "cell_barcode", "whitelist",
                       whitelist="10x-5prime-vdj")
UMI_FEATURE = feature("umi", "umi", "random")
I5_FEATURE = feature("sample_index_i5", "sample_index", "unknown")
I7_FEATURE = feature("sample_index_i7", "sample_index", "unknown")
FEATURES = {"cell barcode": CELL_BARCODE, "UMI": UMI_FEATURE,
            "i5": I5_FEATURE, "i7": I7_FEATURE, "i7 oligo": I7_FEATURE}


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    kw.setdefault("feature", FEATURES.get(name))
    return Segment(name=name, top=top, tag=tag, **kw)


class VDJChemistry:
    """A complete versioned chemistry; construction validates every molecular interlock."""

    def __init__(self, *, version: str, umi_nt: int, guide: str) -> None:
        if umi_nt not in (10, 12):
            raise ValueError("supported 10x 5' V(D)J UMI lengths are 10 (v2) and 12 (v3)")
        self.version = version
        self.umi_nt = umi_nt
        self.guide = guide
        self._validate()

    def gel_bead_segments(self) -> list[Segment]:
        return [seg("partial Read 1", PARTIAL_R1, "r1"),
                seg("cell barcode", "N" * 16, "cbc", placeholder=True),
                seg("UMI", "N" * self.umi_nt, "umi", placeholder=True),
                seg("TSO", TSO[:-3], "tso"), seg("rGrGrG", TSO[-3:], "tso")]

    def polydt_segments(self) -> list[Segment]:
        return [seg("SMART handle", rt.SMART_HANDLE, "tso"), seg("AC", "AC", "tso"),
                seg("poly(dT)30", "T" * POLY_DT_NT),
                seg("VN", "VN", placeholder=True)]

    def index_p5_segments(self) -> list[Segment]:
        return [seg("P5", il.P5, "p5"),
                seg("i5", "N" * INDEX_NT, "cbc", placeholder=True),
                seg("partial Read 1", il.TRUSEQ_READ1[:24], "r1")]

    def index_p7_segments(self) -> list[Segment]:
        return [seg("P7", il.P7, "p7"),
                seg("i7 oligo", "N" * INDEX_NT, "cbc", placeholder=True),
                seg("partial Read 2", il.TRUSEQ_READ2[:21], "r2")]

    def capture_scene(self) -> Scene:
        mrna = [seg("mRNA", "XXXXXXXX...XXXXXXXX", placeholder=True),
                seg("poly(A)", "A" * POLY_DT_NT)]
        sc = Scene()
        sc.strand("mRNA", mrna, label="mRNA")
        sc.anneal("RT primer", self.polydt_segments(), to="mRNA",
                  pair=("poly(dT)30", "poly(A)"), label="poly(dT) RT primer",
                  unpaired=("VN",))
        sc.arrow("RT primer", "reverse transcriptase")
        return sc

    def template_switch_scene(self) -> Scene:
        first = [*self.polydt_segments(),
                 seg("antisense cDNA", "XXXXXXXX...XXXXXXXX", placeholder=True),
                 seg("CCC", rt.UNTEMPLATED_TAIL, "tso")]
        sc = Scene()
        sc.strand("first strand", first, label="first strand cDNA")
        sc.anneal("gel-bead TSO", self.gel_bead_segments(), to="first strand",
                  pair=("rGrGrG", "CCC"), label="gel-bead TSO", above=True)
        sc.arrow("first strand", "RT switches template")
        return sc

    def vdj_amplicon(self) -> Construct:
        return Construct([
            seg("GATCT extension", "GATCT", "p5"),
            seg("Read 1", il.TRUSEQ_READ1, "r1"),
            seg("cell barcode", "N" * 16, "cbc", placeholder=True),
            seg("UMI", "N" * self.umi_nt, "umi", placeholder=True),
            seg("TSO", TSO, "tso"),
            seg("5' UTR–V(D)J–C", "XXXXXXXX...XXXXXXXX", placeholder=True),
            seg("inner constant-region primer", "NNNNNNNNNNNNNNNNNNNN",
                placeholder=True),
        ], name=f"10x 5' V(D)J {self.version} nested amplicon")

    def cdna_template(self) -> Construct:
        """Representative full-length cDNA carrying both published amplification sites."""
        return Construct([
            seg("partial Read 1", CDNA_FORWARD, "r1"),
            seg("cell barcode", "N" * 16, "cbc", placeholder=True),
            seg("UMI", "N" * self.umi_nt, "umi", placeholder=True),
            seg("TSO", TSO, "tso"),
            seg("5' UTR–transcript–poly(A)", "XXXXXXXX...XXXXXXXX", placeholder=True),
            seg("cDNA reverse-primer site", revcomp(CDNA_REVERSE), "r2"),
        ], name=f"10x 5' {self.version} full-length cDNA")

    def cdna_pcr(self) -> PCRGeometry:
        return PCRGeometry.exact(
            self.cdna_template(), PCRPrimer("cDNA Forward", CDNA_FORWARD),
            PCRPrimer("cDNA Reverse", CDNA_REVERSE),
            statement=f"{self.guide} prints both cDNA-primer sequences.")

    def vdj_template(self, *, inner: bool = False) -> Construct:
        return Construct([
            seg("partial Read 1", CDNA_FORWARD, "r1"),
            seg("cell barcode", "N" * 16, "cbc", placeholder=True),
            seg("UMI", "N" * self.umi_nt, "umi", placeholder=True),
            seg("TSO", TSO, "tso"),
            seg("5' UTR–V(D)J", "XXXXXXXX...XXXXXXXX", placeholder=True),
            seg("constant region", "XXXXXXXXXXXXXXXXXXXX", placeholder=True,
                note="outer or inner pool target; exact transcript coordinate not reported"),
        ], name=("inner" if inner else "outer") + " V(D)J PCR template")

    def vdj_pcr(self, *, inner: bool) -> PCRGeometry:
        stage = "inner Mix 2" if inner else "outer Mix 1"
        return PCRGeometry.regional(
            self.vdj_template(inner=inner),
            PCRPrimer("V(D)J universal forward", VDJ_FORWARD, anneal_nt=11),
            reverse_pool=f"{stage} constant-region primer pool",
            reverse_region="constant region",
            statement=(f"{self.guide} prints every {stage} reverse-primer sequence and "
                       "states that the pool binds constant regions, but does not map "
                       "each oligo to a transcript coordinate."))

    def pre_index_library(self) -> Construct:
        return Construct([
            seg("GATCT extension", "GATCT", "p5"),
            seg("Read 1", il.TRUSEQ_READ1, "r1"),
            seg("cell barcode", "N" * 16, "cbc", placeholder=True),
            seg("UMI", "N" * self.umi_nt, "umi", placeholder=True),
            seg("TSO", TSO, "tso"),
            seg("V(D)J fragment", "XXXXXXXX...XXXXXXXX", placeholder=True),
            seg("dA junction", "A"),
            seg("Index 1 / Read 2 arm", il.INDEX1_PRIMER, "r2"),
        ], name=f"10x 5' V(D)J {self.version} adapter-ligated fragment")

    def index_pcr(self) -> PCRGeometry:
        forward = "".join(s.top for s in self.index_p5_segments())
        reverse = "".join(s.top for s in self.index_p7_segments())
        return PCRGeometry.exact(
            self.pre_index_library(), PCRPrimer("TT i5 / P5 primer", forward),
            PCRPrimer("TT i7 / P7 primer", reverse),
            statement=f"{self.guide} prints both Dual Index TT Set A primer structures.")

    def adapter_top(self) -> list[Segment]:
        return [seg("Read 2 adapter", il.INDEX1_PRIMER, "r2")]

    def adapter_bottom(self) -> list[Segment]:
        return [seg("12-bp stem complement", il.STEM_COMPLEMENT, "r2"),
                seg("3' T", "T", "r2")]

    def adapter_scene(self) -> Scene:
        sc = Scene()
        sc.strand("long oligo", self.adapter_top(), label="long adaptor oligo")
        sc.anneal("short oligo", self.adapter_bottom(), to="long oligo",
                  pair=("12-bp stem complement", "Read 2 adapter"),
                  label="short adaptor oligo")
        sc.mark("short oligo", "3' T", "3'-T ligation overhang")
        return sc

    def final_library(self) -> Construct:
        return Construct([
            seg("P5", il.P5, "p5"),
            seg("i5", "N" * INDEX_NT, "cbc", placeholder=True),
            seg("Read 1", il.TRUSEQ_READ1, "r1"),
            seg("cell barcode", "N" * 16, "cbc", placeholder=True),
            seg("UMI", "N" * self.umi_nt, "umi", placeholder=True),
            seg("TSO", TSO, "tso"),
            seg("V(D)J fragment", "XXXXXXXX...XXXXXXXX", placeholder=True),
            seg("dA junction", "A"),
            seg("Index 1 / Read 2 arm", il.INDEX1_PRIMER, "r2"),
            seg("i7", "N" * INDEX_NT, "cbc", placeholder=True),
            seg("P7'", il.P7_RC, "p7"),
        ], name=f"10x 5' V(D)J {self.version} library")

    def protocol_graph(self) -> ProtocolGraph:
        authority=(f'Primary source: <a href="https://www.10xgenomics.com/support">'
                   f'10x Genomics User Guide {self.guide}</a>.')
        g=ProtocolGraph(f"10x Chromium Single Cell 5-prime V(D)J {self.version}",
                        authority=authority)
        cells=g.input("cells","cells containing receptor and poly(A) transcripts",
                      MaterialKind.BIOLOGICAL,
                      state=molecule("cells",(role("V(D)J transcripts"),role("poly(A) RNA"))),
                      primary=True)
        reviewed=["cells"]
        occupied=partition_material(g,cells,reviewed,label="cells")
        rt_primer=g.input("rt_primer","poly(dT) RT primer",MaterialKind.OLIGO)
        rt_enzyme=g.input("rt","template-switching reverse transcriptase",MaterialKind.CATALYST)
        first=g.product("first","first-strand full-length cDNA",MaterialKind.MOLECULAR,
                        state=self.cdna_template())
        g.react("reverse_transcribe","Reverse-transcribe poly(A) transcripts",
                StepKind.REVERSE_TRANSCRIPTION,
                ((occupied,InputRole.TEMPLATE),(rt_primer,InputRole.STRUCTURAL),
                 (rt_enzyme,InputRole.CATALYST)),(first,))
        tso=g.input("tso","barcoded gel-bead template-switch oligo",MaterialKind.OLIGO)
        switched=g.product("switched","cell-barcoded full-length cDNA",MaterialKind.MOLECULAR,
                           state=self.cdna_template())
        g.react("template_switch","Template-switch onto the barcoded gel-bead oligo",
                StepKind.TEMPLATE_SWITCHING,
                ((first,InputRole.TEMPLATE),(tso,InputRole.STRUCTURAL),
                 (rt_enzyme,InputRole.CATALYST)),(switched,))
        reviewed.extend(("rt_primer","rt","reverse_transcribe","tso","template_switch"))
        amplified=chromium_pcr(g,"cdna_pcr",switched,self.cdna_template(),reviewed,
                               label="Amplify full-length cDNA",geometry=self.cdna_pcr())
        outer=chromium_pcr(g,"outer_vdj_pcr",amplified,self.vdj_template(inner=True),reviewed,
                           label="Outer V(D)J PCR with Mix 1",geometry=self.vdj_pcr(inner=False))
        inner=chromium_pcr(g,"inner_vdj_pcr",outer,self.vdj_amplicon(),reviewed,
                           label="Inner V(D)J PCR with Mix 2",geometry=self.vdj_pcr(inner=True))
        fragmentase=g.input("fragmentase","V(D)J amplicon fragmentation enzyme",
                            MaterialKind.CATALYST)
        fragments=g.product("fragments","fragmented barcode-bearing V(D)J amplicons",
                            MaterialKind.MOLECULAR)
        g.react("fragment","Fragment nested V(D)J amplicons",StepKind.FRAGMENTATION,
                ((inner,InputRole.TEMPLATE),(fragmentase,InputRole.CATALYST)),(fragments,))
        end_mix=g.input("end_mix","end-repair and dA-tailing enzymes",MaterialKind.CATALYST)
        tailed=g.product("tailed","end-repaired dA-tailed V(D)J fragments",
                         MaterialKind.MOLECULAR)
        g.react("end_prep","End-repair and dA-tail V(D)J fragments",StepKind.END_REPAIR,
                ((fragments,InputRole.TEMPLATE),(end_mix,InputRole.CATALYST)),(tailed,))
        adapter=g.input("read2_adapter","3-prime-dT Read-2 adapter",MaterialKind.ADAPTER)
        ligase=g.input("ligase","adapter ligase",MaterialKind.CATALYST)
        adapted=g.product("adapted","adapter-ligated barcode-bearing V(D)J fragments",
                          MaterialKind.MOLECULAR,state=self.pre_index_library())
        g.react("ligate","Ligate the Read-2 adapter",StepKind.LIGATION,
                ((tailed,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),
                 (ligase,InputRole.CATALYST)),(adapted,),geometry=READ2_ADAPTER_END)
        reviewed.extend(("fragmentase","fragment","end_mix","end_prep",
                         "read2_adapter","ligase","ligate"))
        library=chromium_pcr(g,"index_pcr",adapted,self.final_library(),reviewed,
                             label="Select the barcode-bearing fragment by dual-index PCR",
                             library=True,geometry=self.index_pcr())
        lengths={"Read 1":26 if self.umi_nt==10 else 28,
                 "Index 1 (i7)":10,"Index 2 (i5)":10,"Read 2":90}
        g.sequence(library,SequencingPlan(library=self.final_library(),
            primers=self.seq_primers,acquired_roles=tuple(p.role for p in self.seq_primers),
            read_lengths=lengths))
        reviewed.append(library.id)
        g.declare_reviewed_scope(reviewed,authority)
        return g.require_publishable()

    @property
    def seq_primers(self):
        return (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["I2"],
                sp.TRUSEQ["R2"])

    def _validate(self) -> None:
        bead = "".join(s.top for s in self.gel_bead_segments())
        if bead != PARTIAL_R1 + "N" * 16 + "N" * self.umi_nt + TSO:
            raise ValueError(f"{self.version} gel-bead primer does not match its guide")
        if VDJ_FORWARD != "GATCT" + il.TRUSEQ_READ1[:22]:
            raise ValueError("V(D)J forward primer no longer reconstructs the Read 1 end")
        if "A" + il.INDEX1_PRIMER != revcomp(il.TRUSEQ_READ2):
            raise ValueError("dA plus adaptor no longer reconstructs the Read 2 site")
        self.capture_scene().rows()
        self.template_switch_scene().rows()
        self.adapter_scene().rows()
        self.cdna_pcr().rows()
        self.vdj_pcr(inner=False).rows()
        self.vdj_pcr(inner=True).rows()
        self.index_pcr().rows()
        Scene.duplex(list(self.vdj_amplicon())).rows()
        Scene.duplex(list(self.final_library())).rows()
        problems = sp.verify(self.final_library(), self.seq_primers)
        if problems:
            raise ValueError(f"invalid {self.version} read layout: " + "; ".join(problems))
