"""Reusable constructs for the protocol-class gap batch.

The helpers here make final-library completeness a construction invariant: a short-read
library is returned only after its declared sequencing primers have been located.
"""
from __future__ import annotations

from dataclasses import dataclass

import illumina as il
import nextera as nx
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry
from chemdraw import (Construct, Row, Scene, Segment, bases_pair, complement_segments,
                      feature, revcomp)
from endprep import EndPreparedDuplex, LigationJunction, dA_tailed_scene, ligation_product_scene


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name, top, tag, **kw)


# Pass as ``i5_nt``/``i7_nt`` when the authoritative source does not state an index
# length: the window is drawn as a bracketed role token, which carries no molecular
# length, instead of an assumed run of index letters.
UNASSIGNED_LENGTH = None


def _check_index_lengths(*lengths: int | None, what: str) -> None:
    if any(n is not UNASSIGNED_LENGTH and n < 1 for n in lengths):
        raise ValueError(f"{what} needs positive (or explicitly unassigned) index lengths")


def _p5_segment(next_arm: str, dual_index: bool, inferred: bool) -> Segment:
    """The P5 end of a PCR-installed left arm.

    With an i5 window the primer is P5 + i5 + arm.  Without one, P5 and the following
    arm are a single primer sequence, so any bases they share (P5's ACAC and TruSeq
    Read 1's leading ACAC) occur once, not twice.  The joined TruSeq arm is checked
    against the canonical 58-nt universal primer.
    """
    p5 = il.P5 if dual_index else il.P5[:len(il.P5) - _overlap(il.P5, next_arm)]
    if (not dual_index and next_arm == il.TRUSEQ_READ1
            and p5 + next_arm != il.TRUSEQ_P5_FULL):
        raise AssertionError("single-index TruSeq left arm differs from TRUSEQ_P5_FULL")
    return seg("P5", p5, "p5", inferred=inferred)


def index_window(name: str, letter: str, nt: int | None, identifier: str, *,
                 inferred: bool = False, sample_feature=None) -> Segment:
    """A sample-index window of ``nt`` placeholder bases, or of unassigned length."""
    if nt is UNASSIGNED_LENGTH:
        return seg(name, f"[{name.split()[0]} index; length not stated]", "cbc",
                   placeholder=True,
                   inferred=inferred,
                   feature=sample_feature or feature(identifier, "sample_index", "unknown",
                                                     note="index length not stated by the source"))
    return seg(name, letter * nt, "cbc", placeholder=True, inferred=inferred,
               feature=sample_feature or feature(identifier, "sample_index", "unknown"))


def truseq_library(insert: list[Segment], name: str, *, dual_index: bool = True,
                    inferred_adapters: bool = False,
                    inferred_junction: bool = False,
                    i5_nt: int | None = 8, i7_nt: int | None = 8,
                    i5_feature=None, i5_name: str = "i5") -> tuple[Construct, tuple]:
    """Canonical TruSeq library with validated R1/I1/(I2)/R2 sites.

    ``inferred_junction`` marks only the single dA junction base as inferred, for
    callers whose source establishes the adapters but not the insert tailing.
    """
    _check_index_lengths(i5_nt, i7_nt, what="a TruSeq indexed library")
    left = [_p5_segment(il.TRUSEQ_READ1, dual_index, inferred_adapters)]
    if dual_index:
        left.append(index_window(i5_name, "J", i5_nt, "sample_index_i5",
                                 inferred=inferred_adapters, sample_feature=i5_feature))
    left.append(seg("Read 1 arm", il.TRUSEQ_READ1, "r1", inferred=inferred_adapters))
    lib = Construct([*left, *insert,
        seg("dA junction", "A", inferred=inferred_adapters or inferred_junction),
        seg("Index 1 / Read 2 arm", il.INDEX1_PRIMER, "r2", inferred=inferred_adapters),
        index_window("i7 reverse complement", "I", i7_nt, "sample_index_i7",
                     inferred=inferred_adapters),
        seg("P7 reverse complement", il.P7_RC, "p7", inferred=inferred_adapters)], name=name)
    primers = ((sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["I2"], sp.TRUSEQ["R2"])
               if dual_index else (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["R2"]))
    roles = tuple(p.role for p in primers)
    problems = sp.verify(lib, primers, required_roles=roles)
    if problems:
        raise ValueError(f"{name}: " + "; ".join(problems))
    return lib, primers


def truseq_amplicon_library(insert: list[Segment], name: str, *,
                             dual_index: bool = True,
                             inferred_adapters: bool = False,
                             i5_nt: int | None = 8,
                             i7_nt: int | None = 8) -> tuple[Construct, tuple]:
    """TruSeq-compatible PCR amplicon, without inventing a dA-ligation junction.

    Amplicon workflows install the terminal sequencing shell through tailed primers.
    Keeping this constructor distinct from :func:`truseq_library` prevents those
    protocols from inheriting the single dA base that belongs only to adapter ligation.
    """
    _check_index_lengths(i5_nt, i7_nt, what="a TruSeq indexed amplicon")
    left = [_p5_segment(il.TRUSEQ_READ1, dual_index, inferred_adapters)]
    if dual_index:
        left.append(index_window("i5", "J", i5_nt, "sample_index_i5",
                                 inferred=inferred_adapters))
    left.append(seg("Read 1 arm", il.TRUSEQ_READ1, "r1",
                    inferred=inferred_adapters))
    lib = Construct([*left, *insert,
        # PCR installs the complete reverse complement of the Read 2 arm.  The ligation
        # constructor omits its leading A because that base is supplied by the insert's
        # dA junction; an amplicon has no such junction and must carry all 34 bases here.
        seg("Index 1 / Read 2 arm", revcomp(il.TRUSEQ_READ2), "r2",
            inferred=inferred_adapters),
        index_window("i7 reverse complement", "I", i7_nt, "sample_index_i7",
                     inferred=inferred_adapters),
        seg("P7 reverse complement", il.P7_RC, "p7",
            inferred=inferred_adapters)], name=name)
    primers = ((sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["I2"], sp.TRUSEQ["R2"])
               if dual_index else (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["R2"]))
    problems = sp.verify(lib, primers, required_roles=tuple(p.role for p in primers))
    if problems:
        raise ValueError(f"{name}: " + "; ".join(problems))
    return lib, primers


# ------------------------------------------------------- PCR-free, full-length adaptor
# A PCR-free library is finished the moment the adaptor is sealed: the adaptor itself
# carries P5, the sample index(es) and P7, so no PCR primer installs a platform sequence
# and the insert is flanked by two ligation junctions rather than two primer sites.
#   Illumina Adapter Sequences #1000000002694 v22, "TruSeq Single Indexes" / "TruSeq DNA
#   and RNA CD Indexes" / "IDT for Illumina–TruSeq DNA and RNA UD Indexes":
#     "A-tailing is performed before adapter ligation. For example, the additional A base
#      is in parentheses in the i7 adapter, as follows.
#      Index 1 (i7) Adapters
#      (A)GATCGGAAGAGCACACGTCTGAACTCCAGTCAC[i7]ATCTCGTATGCCGTCTTCTGCTTG"
#     "Index 2 (i5) Adapters
#      AATGATACGGCGACCACCGAGATCTACAC[i5]ACACTCTTTCCCTACACGACGCTCTTCCGATCT"
#     "TruSeq Universal Adapter
#      5′ AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"
DT_JUNCTION = "dT junction"
DA_JUNCTION = "dA junction"
TRUSEQ_STEM_NT = len(il.STEM)


def _overlap(left: str, right: str) -> int:
    """Length of the longest suffix of ``left`` that is also a prefix of ``right``."""
    n = min(len(left), len(right))
    while n and left[-n:] != right[:n]:
        n -= 1
    return n


def _slice(s: Segment, a: int, b: int, name: str | None = None) -> Segment:
    """A sub-region of one segment, keeping its presentation and identifier semantics."""
    return Segment(name or s.name, s.top[a:b], s.tag, s.placeholder, s.inferred,
                   None if s.bottom is None else s.bottom[a:b], s.note, s.feature)


@dataclass(frozen=True)
class FullLengthAdapter:
    """A pre-indexed adaptor that is already flow-cell complete before it is ligated.

    ``universal`` and ``indexed`` are the two synthesized oligos, 5′→3′, in the
    granularity the sealed library is drawn with: the universal arm ends in the ligating
    3′ dT, and the indexed arm begins with the stem it shares with that arm.  Because both
    platform arms live here, a library assembled from this object receives P5, P7 and every
    sample index by ligation; there is nothing left for a PCR primer to install.

    The fork is a property of the reagent, not of the sealed product, so it is derived
    here rather than drawn by a caller: :meth:`scene` splits the shared stem out of both
    arms, pairs only that stem, leaves both arms unpaired and declares the 3′ dT as a
    cohesive adapter end.  The constructor refuses arms whose stems do not base-pair, so a
    Y-fork that could not anneal cannot be rendered.
    """

    name: str
    universal: tuple[Segment, ...]
    indexed: tuple[Segment, ...]
    stem_nt: int = TRUSEQ_STEM_NT
    inferred: bool = False

    def __post_init__(self) -> None:
        if len(self.universal) < 2 or len(self.indexed) < 2:
            raise ValueError(f"{self.name}: a forked adaptor needs a body on both arms")
        if self.stem_nt < 1:
            raise ValueError(f"{self.name}: the paired stem must be at least one base")
        dt = self.universal[-1]
        if len(dt) != 1 or dt.placeholder:
            raise ValueError(f"{self.name}: the universal arm must end in its ligating base")
        arm, idx = self.universal[-2], self.indexed[0]
        if len(arm) <= self.stem_nt or len(idx) <= self.stem_nt:
            raise ValueError(f"{self.name}: both arms must extend beyond the shared stem")
        if revcomp(idx.top[:self.stem_nt]) != arm.top[-self.stem_nt:]:
            raise ValueError(f"{self.name}: the two arms do not share a paired stem")

    @property
    def overhang_3(self) -> str:
        """The single 3′ base the adaptor presents to a tailed insert."""
        return self.universal[-1].top

    def universal_sequence(self) -> str:
        return "".join(s.top for s in self.universal)

    def indexed_sequence(self) -> str:
        return "".join(s.top for s in self.indexed)

    def scene(self) -> Scene:
        """The ligation-ready reagent: paired stem, two unpaired arms, 3′ dT overhang."""
        n, (arm, dt) = self.stem_nt, (self.universal[-2], self.universal[-1])
        top = [*self.universal[:-2],
               _slice(arm, 0, len(arm) - n),
               _slice(arm, len(arm) - n, len(arm), name="universal stem"), dt]
        idx = self.indexed[0]
        bottom = [_slice(idx, 0, n, name="indexed stem"),
                  _slice(idx, n, len(idx)), *self.indexed[1:]]
        sc = Scene()
        sc.strand("universal arm", top, label="universal arm")
        sc.anneal("indexed arm", bottom, to="universal arm",
                  pair=("indexed stem", "universal stem"), label="indexed arm",
                  unpaired=tuple(s.name for s in bottom[1:]))
        sc.mark("universal arm", "universal stem", f"{n}-bp paired stem; both arms fork")
        sc.overhang("universal arm", dt.name, polarity="3-prime", owner="adapter",
                    text=f"3′ d{self.overhang_3} ligating end")
        sc.labels("universal arm")
        return sc

    def geometry(self, substrate: EndPreparedDuplex, *,
                 statement: str) -> AdapterLigationGeometry:
        """The two compatible reagent ends, with the fork that carries P5/index/P7."""
        if not substrate.five_prime_phosphate:
            raise ValueError(f"{self.name}: ligation needs a 5′-phosphorylated substrate")
        if not bases_pair(self.overhang_3, substrate.three_prime_overhang):
            raise ValueError(f"{self.name}: adaptor and substrate overhangs are incompatible")
        return AdapterLigationGeometry.known(
            "cohesive", dA_tailed_scene(substrate), self.scene(), statement=statement)

    def ligate(self, insert: list[Segment], *,
               insert_overhang: str = "A") -> list[Segment]:
        """The sealed segment list: one adaptor copy at each end of ``insert``.

        Both junctions are validated, not just the one that needs a drawn base.  The left
        junction's ``A`` is the complement of the adaptor's own 3′ dT and therefore comes
        from the drawing convention; the right junction's dA is the insert's own tailed
        base and is drawn as a segment.
        """
        if not insert:
            raise ValueError(f"{self.name}: a ligated library needs an insert")
        left = LigationJunction(self.overhang_3, insert_overhang)
        right = LigationJunction(self.overhang_3, insert_overhang)
        return [*self.universal, *insert,
                seg(DA_JUNCTION, right.insert_overhang, inferred=self.inferred,
                    note=f"the insert's own 3′ d{right.insert_overhang}, sealed to the "
                         f"adaptor's 3′ d{left.adapter_overhang}"),
                *self.indexed]


def truseq_full_length_adapter(*, dual_index: bool = True, i5_nt: int | None = 8,
                                i7_nt: int | None = 8,
                                i5_feature=None, i5_name: str = "i5",
                                inferred: bool = False,
                                name: str = "full-length TruSeq adaptor") -> FullLengthAdapter:
    """The indexed TruSeq adaptor as ordered: P5/i5/Read 1 arm plus Read 2 arm/i7/P7."""
    _check_index_lengths(i7_nt, *((i5_nt,) if dual_index else ()),
                         what="a full-length indexed adaptor")
    carried = "carried by the ligated full-length adaptor, not by a PCR primer"
    universal = [seg("P5", il.P5, "p5", inferred=inferred)]
    if dual_index:
        read1 = il.TRUSEQ_READ1
        universal.append(index_window(
            i5_name, "J", i5_nt, "sample_index_i5", inferred=inferred,
            sample_feature=i5_feature or feature("sample_index_i5", "sample_index",
                                                 "unknown", note=carried)))
    else:
        # One oligo, so P5 and the Read 1 arm share their ACAC overlap instead of
        # repeating it the way a two-primer PCR shell does.
        read1 = il.TRUSEQ_READ1[_overlap(il.P5, il.TRUSEQ_READ1):]
        if il.P5 + read1 != il.TRUSEQ_P5_FULL:
            raise ValueError("derived universal adaptor differs from the canonical oligo")
    universal += [seg("Read 1 arm", read1[:-1], "r1", inferred=inferred),
                  seg(DT_JUNCTION, read1[-1], "r1", inferred=inferred,
                      note="the adaptor's ligating 3′ dT; its partner is the insert's 3′ dA")]
    indexed = [seg("Index 1 / Read 2 arm", il.INDEX1_PRIMER, "r2", inferred=inferred),
               index_window("i7 reverse complement", "I", i7_nt, "sample_index_i7",
                            inferred=inferred,
                            sample_feature=feature("sample_index_i7", "sample_index",
                                                   "unknown", note=carried)),
               seg("P7 reverse complement", il.P7_RC, "p7", inferred=inferred)]
    return FullLengthAdapter(name, tuple(universal), tuple(indexed), TRUSEQ_STEM_NT, inferred)


def truseq_pcr_free_library(insert: list[Segment], name: str, *, dual_index: bool = True,
                             inferred_adapters: bool = False,
                             i5_nt: int | None = 8, i7_nt: int | None = 8,
                             i5_feature=None, i5_name: str = "i5",
                             insert_overhang: str = "A",
                             adapter: FullLengthAdapter | None = None,
                             ) -> tuple[Construct, tuple]:
    """A PCR-free TruSeq library: every platform sequence arrives by ligation.

    The sequence is the architecture a PCR-free kit produces and is kept distinct from
    :func:`truseq_library` for provenance, not for bases: here both flanks are the ligated
    adaptor, both junctions are sealed nicks, and the run has no amplification step that
    could install or copy an arm.  Keeping the constructors separate stops an amplified
    protocol from inheriting a two-ended ligation claim, and stops a PCR-free protocol from
    silently presenting its indexes as PCR-added.

    The returned construct is linear in the project's one-top-strand convention, as every
    other ligation product is; the reagent's Y-fork lives in
    :meth:`FullLengthAdapter.geometry`, which is where a ligation's end geometry belongs.
    """
    adapter = adapter or truseq_full_length_adapter(
        dual_index=dual_index, i5_nt=i5_nt, i7_nt=i7_nt, i5_feature=i5_feature,
        i5_name=i5_name, inferred=inferred_adapters)
    lib = Construct(adapter.ligate(insert, insert_overhang=insert_overhang), name=name)
    primers = ((sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["I2"], sp.TRUSEQ["R2"])
               if dual_index else (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["R2"]))
    problems = sp.verify(lib, primers, required_roles=tuple(p.role for p in primers))
    if problems:
        raise ValueError(f"{name}: " + "; ".join(problems))
    return lib, primers


def truseq_pcr_free_rows(library: Construct, insert_first: str, *, label: str = "",
                          adapter_group: str = "full-length TruSeq adaptor",
                          highlight_components: bool = True) -> list[Row]:
    """Draw the sealed product so the insert keeps its own dA and both nicks are marked."""
    return ligation_product_scene(
        library, insert_first, DA_JUNCTION, label=label or library.name,
        adapter_group=adapter_group,
        highlight_components=highlight_components).rows()


def nextera_library(insert: list[Segment], name: str, *,
                    inferred_adapters: bool = False, dual_index: bool = True,
                    i7_nt: int | None = 8,
                    i5_nt: int | None = 8) -> tuple[Construct, tuple]:
    """Heterologous S5/S7 Tn5 product after indexed PCR, with validated run sites."""
    _check_index_lengths(i7_nt, *((i5_nt,) if dual_index else ()),
                         what="a Nextera indexed library")
    left = [_p5_segment(nx.S5, dual_index, inferred_adapters)]
    if dual_index:
        left.append(index_window("i5", "J", i5_nt, "sample_index_i5",
                                 inferred=inferred_adapters))
    left.append(seg("S5", nx.S5, "s5", inferred=inferred_adapters))
    lib = Construct([*left,
        seg("left mosaic end", nx.ME, "me", inferred=inferred_adapters), *insert,
        seg("right mosaic end reverse complement", nx.ME_RC, "me",
            inferred=inferred_adapters),
        seg("S7 reverse complement", nx.S7_RC, "s7", inferred=inferred_adapters),
        index_window("i7 reverse complement", "I", i7_nt, "sample_index_i7",
                     inferred=inferred_adapters),
        seg("P7 reverse complement", il.P7_RC, "p7", inferred=inferred_adapters)], name=name)
    primers = ((sp.NEXTERA["R1"], sp.NEXTERA["I1"], sp.NEXTERA["I2"], sp.NEXTERA["R2"])
               if dual_index else (sp.NEXTERA["R1"], sp.NEXTERA["I1"], sp.NEXTERA["R2"]))
    problems = sp.verify(lib, primers, required_roles=tuple(p.role for p in primers))
    if problems:
        raise ValueError(f"{name}: " + "; ".join(problems))
    return lib, primers


def duplex_fragment(label: str = "genomic fragment", length: int = 34) -> Scene:
    return Scene.duplex([seg(label, "X" * length, placeholder=True)], label=label)


def joined_scene(parts: list[Segment], junctions: tuple[tuple[str, str, str], ...],
                 *, label: str, duplex: bool = True) -> Scene:
    """Draw covalently joined parts; ``Scene.junction`` validates every boundary."""
    if duplex:
        sc = Scene.duplex(parts, label=label)
        strand = "top"
    else:
        sc = Scene()
        strand = "molecule"
        sc.strand(strand, parts, label=label)
    for left, right, kind in junctions:
        sc.junction(strand, left, right, kind)
    sc.labels(strand)
    return sc


def spatial_rt_scene(*, barcode_parts: tuple[tuple[str, int], ...], umi: int,
                     surface: str) -> Scene:
    """Poly(A)-primed spatial RT; barcode provenance is encoded in named segments."""
    mrna = [seg("transcript", "X" * 28, placeholder=True), seg("poly(A)", "A" * 18)]
    primer = [seg(surface, "X" * 6, placeholder=True)]
    primer.extend(
        seg(name, "B" * n, "cbc", placeholder=True,
            feature=feature(f"spatial_barcode_{i}", "spatial_barcode", "combinatorial",
                            group="spatial_barcode", part=name))
        for i, (name, n) in enumerate(barcode_parts, 1)
    )
    if umi:
        primer.append(seg("UMI", "U" * umi, "umi", placeholder=True,
                          feature=feature("umi", "umi", "random")))
    primer.append(seg("poly(dT)", "T" * 18))
    sc = Scene(); sc.strand("mRNA", mrna, label="tissue mRNA")
    sc.anneal("capture primer", primer, to="mRNA", pair=("poly(dT)", "poly(A)"),
              label="spatial capture oligo")
    sc.arrow("capture primer", "reverse transcription")
    sc.labels("capture primer")
    return sc


def targeting_rows(enzyme: str, result: str) -> list[Row]:
    return [
        Row(chunks=[("chromatin ---- [target epitope] ---- chromatin", None, False)]),
        Row(chunks=[("                         ^ primary antibody", "cbc", False)]),
        Row(chunks=[(f"                         | Protein A–{enzyme}", "me", False)]),
        Row(chunks=[(f"                         └─ {result}", None, False)]),
    ]


def duplex_family_rows() -> list[Row]:
    return [
        Row(chunks=[("original duplex tags:       α — insert strand A — β", "umi", False)]),
        Row(chunks=[("complementary partner:      β — insert strand B — α", "umi", False)]),
        Row(chunks=[("PCR descendants → SSCS(A)       PCR descendants → SSCS(B)", None, False)]),
        Row(chunks=[("                    SSCS(A) + SSCS(B) → duplex consensus", None, False)]),
    ]
