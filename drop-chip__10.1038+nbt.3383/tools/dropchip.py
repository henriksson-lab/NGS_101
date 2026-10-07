"""Molecular construct model for Drop-ChIP (Rotem et al., 2015)."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from chemdraw import Construct, Scene, Segment, complement_segments, revcomp

# Representative row 1 from Supplementary Table 2.  Every other barcode adaptor has
# the same core and replaces the eight-base barcode and its reversed copy.
REPRESENTATIVE_ADAPTOR = (
    "TTAAGGGCTTTCGTATCCGGGGGACCTTAATTAAGGTGGGGGGGATACCTTTCGGGTTAA"
)
REPRESENTATIVE_BARCODE = "GGGCTTTC"
ADAPTOR_CORE = "GTATCCGGGGGACCTTAATTAAGGTGGGGGGGATAC"
END = "TTAA"

SC_PCR1 = "TAAGGTGGGGGGGATAC"
SC_PCR2 = "TAAGGTCCCCCGGATAC"
PACI = "TTAATTAA"
BCIVI = "GTATCC"


def _seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


def barcode_adaptor(barcode: str = REPRESENTATIVE_BARCODE) -> Construct:
    """One published 60-nt adaptor strand, with symmetry enforced by construction."""
    if len(barcode) != 8 or any(b not in "ACGT" for b in barcode):
        raise ValueError("Drop-ChIP barcode must be eight DNA bases")
    return Construct([
        _seg("left blunt end", END),
        _seg("barcode", barcode, "cbc"),
        _seg("symmetric adaptor core", ADAPTOR_CORE),
        _seg("reversed barcode", barcode[::-1], "cbc"),
        _seg("right blunt end", END),
    ], name="Drop-ChIP barcode adaptor")


def adaptor_scene(barcode: str = REPRESENTATIVE_BARCODE) -> Scene:
    return Scene.duplex(list(barcode_adaptor(barcode)))


def blunt_ligated_fragment(barcode: str = REPRESENTATIVE_BARCODE) -> Construct:
    """Representative one of four ligation orientations: same top-strand orientation.

    Both adaptor/insert junctions are TTAA-to-genome blunt junctions.  The duplex itself
    is derived by Scene at render time rather than independently entered.
    """
    left = barcode_adaptor(barcode)
    right = barcode_adaptor(barcode)
    return Construct([
        *[_seg("left " + s.name, s.top, s.tag, placeholder=s.placeholder) for s in left],
        _seg("nucleosomal DNA", "XXXXXXXX...XXXXXXXX", placeholder=True),
        *[_seg("right " + s.name, s.top, s.tag, placeholder=s.placeholder) for s in right],
    ], name="blunt-ligated Drop-ChIP fragment")


def paci_product(barcode: str = REPRESENTATIVE_BARCODE) -> Construct:
    """Adaptor-labelled fragment after PacI removes the outer adaptor halves."""
    full = barcode_adaptor(barcode).top()
    site = full.index(PACI)
    # PacI cuts TTAAT^TAA.  The left molecule retains sequence after that top-strand cut;
    # at the right junction it retains the reciprocal prefix.
    left = full[site + 5:]
    right = full[:site + 3]
    if not left.startswith(SC_PCR1) or revcomp(SC_PCR2) not in right:
        raise ValueError("PacI product no longer exposes both published SC-PCR sites")
    if left != SC_PCR1 + barcode[::-1] + END \
            or right != END + barcode + revcomp(SC_PCR2):
        raise ValueError("PacI product no longer preserves the symmetric barcode junctions")
    return Construct([
        _seg("SC-PCR1 site", SC_PCR1),
        _seg("left barcode", barcode[::-1], "cbc"),
        _seg("left PacI half", END),
        _seg("nucleosomal DNA", "XXXXXXXX...XXXXXXXX", placeholder=True),
        _seg("right PacI half", END),
        _seg("right barcode", barcode, "cbc"),
        _seg("SC-PCR2 site'", revcomp(SC_PCR2)),
    ], name="PacI-trimmed Drop-ChIP fragment")


def sc_pcr_scene(barcode: str = REPRESENTATIVE_BARCODE) -> Scene:
    """PacI product with both SC-PCR primers placed on their derived sites."""
    top = list(paci_product(barcode))
    bottom = complement_segments(top)
    sc = Scene()
    sc.strand("top", top, label="template")
    sc.anneal("bottom", bottom, to="top", pair=(top[0].name + "'", top[0].name),
              label="template")
    sc.anneal("SC-PCR1", [_seg("SC-PCR1 primer", SC_PCR1)], to="bottom",
              pair=("SC-PCR1 primer", "SC-PCR1 site'"), label="SC-PCR1")
    sc.anneal("SC-PCR2", [_seg("SC-PCR2 primer", SC_PCR2)], to="top",
              pair=("SC-PCR2 primer", "SC-PCR2 site'"), label="SC-PCR2", above=True)
    sc.arrow("SC-PCR1", "")
    sc.arrow("SC-PCR2", "")
    sc.stack("SC-PCR2", "top", "bottom", "SC-PCR1")
    return sc


def bcivi_product(barcode: str = REPRESENTATIVE_BARCODE) -> tuple[Scene, Construct]:
    """BciVI-cut insert, with a single 3'-A overhang at each end.

    The exact retained sequence is sliced from the symmetric adaptor constructor.  The
    right terminal A is explicit; the reciprocal bottom-strand A is derived as its own
    unpaired segment in the Scene.
    """
    full = barcode_adaptor(barcode).top()
    left_site = full.index(BCIVI)
    right_site = full.index(revcomp(BCIVI))
    left = full[right_site - 5:]       # GGGGGGGATAC · reversed barcode · TTAA
    right = full[:left_site + 12]      # TTAA · barcode · GTATCCGGGGGA
    if not right.endswith("A"):
        raise ValueError("BciVI right product must end in the published single-A overhang")
    paired = [
        _seg("left 11-nt constant", left[:11]),
        _seg("left barcode", left[11:19], "cbc"),
        _seg("left PacI half", left[19:]),
        _seg("nucleosomal DNA", "XXXXXXXX...XXXXXXXX", placeholder=True),
        _seg("right PacI half", right[:4]),
        _seg("right barcode", right[4:12], "cbc"),
        _seg("right 11-nt constant", right[12:-1]),
    ]
    product = Construct([*paired, _seg("right 3' A", right[-1])],
                        name="BciVI-cut Drop-ChIP insert")
    sc = Scene()
    sc.strand("top", list(product), label="BciVI product")
    sc.anneal("bottom", [*complement_segments(paired), _seg("left 3' A", "A")],
              to="top", pair=(paired[0].name + "'", paired[0].name),
              label="BciVI product", unpaired=("left 3' A",))
    sc.mark("top", "right 3' A", "3'-A")
    sc.mark("bottom", "left 3' A", "3'-A")
    return sc, product


def inferred_final_library(barcode: str = REPRESENTATIVE_BARCODE) -> Construct:
    """Only the supported inner architecture; outer Illumina sequences are unpublished."""
    _, inner = bcivi_product(barcode)
    return Construct([
        _seg("left unpublished Illumina arm", "[ILLUMINA ARM]", inferred=True,
             placeholder=True),
        *list(inner),
        _seg("right unpublished Illumina arm", "[ILLUMINA ARM]", inferred=True,
             placeholder=True),
    ], name="Drop-ChIP final library, authoritative boundary")


def _validate_source_transcription() -> None:
    if barcode_adaptor().top() != REPRESENTATIVE_ADAPTOR:
        raise ValueError("representative adaptor no longer matches Supplementary Table 2")
    if barcode_adaptor().top().count(PACI) != 1:
        raise ValueError("barcode adaptor must contain one central PacI site")
    if barcode_adaptor().top().count(BCIVI) != 1 \
            or barcode_adaptor().top().count(revcomp(BCIVI)) != 1:
        raise ValueError("barcode adaptor must expose reciprocal BciVI sites")
    adaptor_scene().rows()
    sc_pcr_scene().rows()
    bcivi_product()[0].rows()


_validate_source_transcription()
