"""NanoSeq molecular states (Abascal et al., Nature 2021, author manuscript Methods).

Source scope (author accepted manuscript, Cambridge repository; published Supplementary
Information Notes 2 and 4):
* HpyCH4V (TG^CA) fragmentation leaves blunt ends beginning 5'-CA (Methods "NanoSeq
  library preparation"; Supplementary Note 2; read-processing Methods "NNNTCA").
* A-tailing uses Klenow fragment (3'->5' exo-) with equimolar dATP/ddBTPs, so a strand
  extended from an internal nick terminates on a ddBTP and becomes unamplifiable
  (Supplementary Note 4).
* xGen Duplex Seq Adapters (IDT 1080799) carry a 3-nt random barcode followed by a 3' T
  overhang, optionally preceded by a spacer base (read-processing Methods). The full
  adapter strands are not printed in the inspected sources.
"""
from chemdraw import Construct, Scene, Segment, feature
from endprep import EndPreparedDuplex, dA_tailed_scene

INSERT_N = 30


def _fragment_body() -> list[Segment]:
    return [Segment("HpyCH4V 5' end", "CA"),
            Segment("genomic insert", "X" * INSERT_N, placeholder=True),
            Segment("HpyCH4V 3' end", "TG")]


def hpych4v_fragment() -> Scene:
    """Blunt HpyCH4V fragment: TG^CA cleavage leaves 5'-CA ... TG-3' on both strands."""
    from chemdraw import complement_segments
    body = _fragment_body()
    bottom = complement_segments(body)
    sc = Scene(); sc.strand("top", body, label="HpyCH4V-generated blunt genomic duplex")
    sc.anneal("bottom", bottom, to="top", pair=(bottom[-1].name, body[0].name),
              label="HpyCH4V-generated blunt genomic duplex")
    sc.note("top", "TG^CA cleavage: blunt ends, each strand starts 5'-CA")
    return sc


def a_tailed_product() -> EndPreparedDuplex:
    """Klenow exo- dA tail on each 3' end; 5'-phosphate from restriction cleavage (derived)."""
    return EndPreparedDuplex(Construct(_fragment_body(), name="HpyCH4V fragment"),
                             five_prime_phosphate=True, three_prime_overhang="A")


def a_tailed_fragment() -> Scene:
    return dA_tailed_scene(a_tailed_product(), "dA-tailed HpyCH4V fragment")


def ddbtp_blocked_nick() -> Scene:
    """Off-product: extension from an internal nick stops at an incorporated ddBTP.

    Drawn as two top-strand pieces: the 3' end at the nick now carries a ddBTP
    (no 3'-OH), so adapter ligation cannot seal it and the strand cannot be copied
    through (Supplementary Note 4). Position and length of the extension are schematic.
    """
    up = [Segment("upstream strand at nick", "X" * 14, placeholder=True),
          Segment("incorporated ddBTP", "X", placeholder=True, inferred=True)]
    sc = Scene(); sc.strand("upstream", up, label="nicked strand")
    sc.anneal("bottom", [Segment("intact partner strand", "x" * 31, placeholder=True)],
              to="upstream", pair=("intact partner strand", "upstream strand at nick"),
              label="partner strand")
    sc.anneal("downstream", [Segment("downstream strand", "X" * 16, placeholder=True)],
              to="bottom", pair=("downstream strand", "intact partner strand"), shift=15,
              above=True, label="nicked strand")
    sc.mark("upstream", "incorporated ddBTP",
            "ddBTP-terminated 3' end: strand unamplifiable (Supplementary Note 4)")
    return sc


def umi(name: str, ident: str) -> Segment:
    return Segment(name, "NNN", "umi", placeholder=True, length_bp=3,
                   feature=feature(ident, "umi", "random"))


def duplex_adapter_end() -> Scene:
    """Source-known adapter end: ...NNN-T-3' overhang; the remaining arm is not printed."""
    top = [Segment("xGen adapter arm (not printed)", "[adapter arm]", placeholder=True,
                   inferred=True),
           umi("3-nt random barcode", "nanoseq_adapter_umi"),
           Segment("adapter 3' T overhang", "T")]
    bottom = [Segment("barcode complement", "nnn", placeholder=True),
              Segment("adapter arm partner (not printed)", "[arm partner]",
                      placeholder=True, inferred=True)]
    sc = Scene(); sc.strand("adapter top", top, label="xGen Duplex Seq Adapter end")
    sc.anneal("adapter bottom", bottom, to="adapter top",
              pair=("barcode complement", "3-nt random barcode"),
              label="xGen Duplex Seq Adapter end",
              unpaired=("adapter arm partner (not printed)",))
    sc.overhang("adapter top", "adapter 3' T overhang", polarity="3-prime", owner="adapter",
                text="3′ T")
    return sc


# qPCR quantification primers, author-manuscript Methods "DNA quantification, dilution and
# PCR amplification" (both HPLC-purified, no modifications printed).
NANOQPCR1 = Construct([Segment("NanoqPCR1", "ACACTCTTTCCCTACACGAC", "r1")], name="NanoqPCR1")
NANOQPCR2 = Construct([Segment("NanoqPCR2", "GTGACTGGAGTTCAGACGTG", "r2")], name="NanoqPCR2")


def oligo_scene(con: Construct) -> Scene:
    sc = Scene(); sc.strand("oligo", list(con), label=con.name)
    return sc
