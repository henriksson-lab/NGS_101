"""Shared construct builder for vendor-defined 10x 5' V(D)J chemistries."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import illumina as il
import rt
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, feature, revcomp

PARTIAL_R1 = il.TRUSEQ_READ1[11:]
TSO = "TTTCTTATATGGG"
VDJ_FORWARD = "GATCTACACTCTTTCCCTACACGACGC"
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
        Scene.duplex(list(self.vdj_amplicon())).rows()
        Scene.duplex(list(self.final_library())).rows()
        problems = sp.verify(self.final_library(), self.seq_primers)
        if problems:
            raise ValueError(f"invalid {self.version} read layout: " + "; ".join(problems))
