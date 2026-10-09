"""Quick et al. 2017 PrimalSeq / ARTIC tiled-amplicon workflow."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import illumina as il
import seqprimers as sp
from batch_ngs import seg
from chemdraw import Construct, Row, Scene, feature
from tiled_amplicon import alternating_scheme, scheme_rows

TITLE = "ARTIC / PrimalSeq — two-pool tiled viral amplicon sequencing"
NOTES = "01_artic-tiled-amplicon.html"
SOURCE = ('Defining protocol: <a href="https://doi.org/10.1038/nprot.2017.066">'
          'Quick et al., <i>Nature Protocols</i> (2017)</a>.')
SUMMARY = ("Overlapping amplicons tile a viral genome, but adjacent products alternate between "
           "two multiplex PCR pools so that overlapping primers cannot make preferential short products. "
           "The products then enter distinct MinION or MiSeq library branches.")
CAVEAT = ("The paper names SureSelectXT2/KAPA but does not print the indexing-adapter bases. "
          "The MiSeq adapter geometry is therefore shown as inferred Illumina-compatible sequence; "
          "the two-pool amplicons and both platform branches are source-defined.")

# Compact representative tiling path. Coordinates are diagram units, not a particular virus scheme.
SCHEME = alternating_scheme([(0, 12), (9, 21), (18, 30), (27, 39), (36, 48), (45, 57)])

FINAL_LIBRARY = Construct([
    seg("P5", il.P5, "p5", inferred=True),
    seg("Read 1 arm", il.TRUSEQ_READ1, "r1", inferred=True),
    seg("one tiled viral amplicon", "X" * 42, placeholder=True),
    seg("dA junction", "A", inferred=True),
    seg("Index 1 / Read 2 arm", il.INDEX1_PRIMER, "r2", inferred=True),
    seg("sample index reverse complement", "I" * 8, "cbc", placeholder=True, inferred=True,
        feature=feature("sample_i7", "sample_index", "unknown")),
    seg("P7 reverse complement", il.P7_RC, "p7", inferred=True),
], name="inferred SureSelectXT2-indexed MiSeq branch")
SEQ_PRIMERS = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["R2"])
if problems := sp.verify(FINAL_LIBRARY, SEQ_PRIMERS,
                         required_roles=tuple(p.role for p in SEQ_PRIMERS)):
    raise ValueError("ARTIC MiSeq branch: " + "; ".join(problems))

FINAL_CAPTION = ("INFERRED — MiSeq branch after KAPA end repair/dA-tailing and SureSelectXT2 "
                 "index-adapter ligation; one approximately 400-base tiled amplicon is shown.")
SEQUENCING_INTRO = ("The protocol recommends paired 250-base MiSeq reads for 400-base amplicons. "
                    "These standard primer sites express the inferred Illumina-compatible geometry.")


def illumina_ligation_scene() -> Scene:
    sc = Scene.duplex(list(FINAL_LIBRARY), label="MiSeq branch library")
    sc.junction("top", "Read 1 arm", "one tiled viral amplicon", "adapter ligation")
    sc.junction("top", "dA junction", "Index 1 / Read 2 arm", "adapter ligation")
    sc.labels("top")
    return sc


def nanopore_ligation_scene() -> Scene:
    parts = [
        seg("sequencing adapter and motor", "[sequencing adapter + motor]", "me", placeholder=True),
        seg("native barcode", "[native barcode]", "cbc", placeholder=True,
            feature=feature("sample_native", "sample_index", "whitelist", whitelist="Oxford Nanopore native barcode set")),
        seg("tiled amplicon", "[tiled viral amplicon]", placeholder=True),
    ]
    sc = Scene(); sc.strand("nanopore library", parts, label="MinION branch")
    sc.junction("nanopore library", "sequencing adapter and motor", "native barcode", "adapter ligation")
    sc.junction("nanopore library", "native barcode", "tiled amplicon", "barcode ligation")
    sc.labels("nanopore library")
    return sc


def sections():
    return [
        ("Design one alternating tiling path", scheme_rows(SCHEME),
         "Neighboring overlapping amplicons are assigned to different pools by construction."),
        ("Amplify the two pools separately", [
            Row(chunks=[("cDNA or viral DNA + pool 1 -> odd tiled amplicons", "cbc", False)]),
            Row(chunks=[("cDNA or viral DNA + pool 2 -> even tiled amplicons", "cbc", False)]),
            Row(chunks=[("98 °C denaturation / 65 °C combined anneal-extension; clean each pool", None, False)]),
        ], "Low per-primer concentration and a five-minute anneal/extension support a highly multiplexed reaction."),
        ("MinION branch", nanopore_ligation_scene().rows(),
         "After end repair/dA-tailing, pools may receive separate native barcodes. The paper used 2D SQK-LSK208 and also specifies compatibility with 1D SQK-LSK108."),
        ("INFERRED — MiSeq branch", illumina_ligation_scene().rows(),
         "KAPA Hyper performs end repair/dA-tailing and amplification; SureSelectXT2 indexing adapters replace the KAPA adapters."),
    ]
