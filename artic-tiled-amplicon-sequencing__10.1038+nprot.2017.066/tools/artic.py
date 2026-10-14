"""Quick et al. 2017 PrimalSeq / ARTIC tiled-amplicon workflow."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import illumina as il
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry
from batch_ngs import seg
from chemdraw import Construct, Row, Scene, Segment, feature, oligo
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
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

# ZikaAsian pool 1, first tiled amplicon; Quick et al. Supplementary Table 1.
ZIKA_400_1_LEFT = "GACAGTTCGAGTTTGAAGCGAAAG"
ZIKA_400_1_RIGHT = "AGTATGCACTCCCACGTCTAGT"

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


def oligos():
    yield oligo("ZikaAsian ZIKA_400_1_LEFT (pool 1)",
                [Segment("viral forward primer", ZIKA_400_1_LEFT)])
    yield oligo("ZikaAsian ZIKA_400_1_RIGHT (pool 1)",
                [Segment("viral reverse primer", ZIKA_400_1_RIGHT)])


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
    ont = AdapterLigationGeometry.unavailable("Quick et al. name SQK-LSK208/SQK-LSK108 and native barcoding, but Oxford Nanopore does not disclose the motor/barcode ligation-facing sequences or overhangs for these kits.")
    miseq = AdapterLigationGeometry.unavailable("Quick et al. name SureSelectXT2 indexing adapters but do not print their ligation-facing strands or overhang geometry.")
    return [
        ("Design one alternating tiling path", scheme_rows(SCHEME),
         "Neighboring overlapping amplicons are assigned to different pools by construction."),
        ("Amplify the two pools separately", [
            Row(chunks=[("cDNA or viral DNA + pool 1 -> odd tiled amplicons", "cbc", False)]),
            Row(chunks=[("cDNA or viral DNA + pool 2 -> even tiled amplicons", "cbc", False)]),
            Row(chunks=[("98 °C denaturation / 65 °C combined anneal-extension; clean each pool", None, False)]),
        ], "Low per-primer concentration and a five-minute anneal/extension support a highly multiplexed reaction."),
        ("MinION branch", ont.rows()+nanopore_ligation_scene().rows(),
         "After end repair/dA-tailing, pools may receive separate native barcodes. The paper used 2D SQK-LSK208 and also specifies compatibility with 1D SQK-LSK108."),
        ("INFERRED — MiSeq branch", miseq.rows()+illumina_ligation_scene().rows(),
         "KAPA Hyper performs end repair/dA-tailing and amplification; SureSelectXT2 indexing adapters replace the KAPA adapters."),
    ]


def protocol_graph():
    graph=ProtocolGraph(TITLE,authority=SOURCE)
    viral=graph.input("viral_template","viral DNA or reverse-transcribed cDNA",
                      MaterialKind.MOLECULAR,primary=True)
    pool1_template=graph.product("pool1_template","aliquot for odd amplicons",
                                 MaterialKind.MOLECULAR)
    pool2_template=graph.product("pool2_template","aliquot for even amplicons",
                                 MaterialKind.MOLECULAR)
    graph.react("split_pools","Split template between alternating primer pools",
                StepKind.SPLIT,((viral,InputRole.TEMPLATE),),
                (pool1_template,pool2_template),changes_molecule=False)
    pool1_f=graph.input("pool1_forward","odd-pool forward-primer members including source-printed ZIKA_400_1_LEFT",MaterialKind.OLIGO)
    pool1_r=graph.input("pool1_reverse","odd-pool reverse-primer members including source-printed ZIKA_400_1_RIGHT",MaterialKind.OLIGO)
    pool2_f=graph.input("pool2_forward","even-pool forward-primer members",MaterialKind.OLIGO)
    pool2_r=graph.input("pool2_reverse","even-pool reverse-primer members",MaterialKind.OLIGO)
    polymerase=graph.input("amplicon_polymerase","multiplex PCR polymerase",
                           MaterialKind.CATALYST)
    odd=graph.product("odd_amplicons","odd tiled viral amplicons",MaterialKind.MOLECULAR)
    even=graph.product("even_amplicons","even tiled viral amplicons",MaterialKind.MOLECULAR)
    tiled_template=Construct([seg("representative tiled viral target","X"*42,
                                  placeholder=True)],name="variable tiled viral targets")
    graph.react("amplify_odd","Amplify odd tiled amplicons",StepKind.PCR,
                ((pool1_template,InputRole.TEMPLATE),(pool1_f,InputRole.STRUCTURAL),
                 (pool1_r,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),
                (odd,),geometry=PCRGeometry.unavailable(
                    tiled_template,reason="Each pool contains many source-defined inward primer pairs rather than one PCR geometry."))
    graph.react("amplify_even","Amplify even tiled amplicons",StepKind.PCR,
                ((pool2_template,InputRole.TEMPLATE),(pool2_f,InputRole.STRUCTURAL),
                 (pool2_r,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),
                (even,),geometry=PCRGeometry.unavailable(
                    tiled_template,reason="Each pool contains many source-defined inward primer pairs rather than one PCR geometry."))
    pooled=graph.product("pooled_amplicons","cleaned pooled tiled amplicons",
                         MaterialKind.MOLECULAR)
    graph.react("pool_amplicons","Pool odd and even tiled amplicons",StepKind.POOL,
                ((odd,InputRole.TEMPLATE),(even,InputRole.TEMPLATE)),(pooled,),
                changes_molecule=False)
    ont_input=graph.product("ont_input","amplicons assigned to the MinION branch",
                            MaterialKind.MOLECULAR)
    miseq_input=graph.product("miseq_input","amplicons assigned to the MiSeq branch",
                              MaterialKind.MOLECULAR)
    graph.react("split_platforms","Split pooled amplicons by sequencing platform",
                StepKind.SPLIT,((pooled,InputRole.TEMPLATE),),(ont_input,miseq_input),
                changes_molecule=False)
    ont_endprep=graph.input("ont_endprep","MinION end-repair and dA-tailing enzyme mix",
                            MaterialKind.CATALYST)
    ont_tailed=graph.product("ont_tailed","end-repaired dA-tailed MinION amplicons",
                             MaterialKind.MOLECULAR)
    graph.react("prepare_ont_ends","End-repair and dA-tail the MinION branch",
                StepKind.END_REPAIR,
                ((ont_input,InputRole.TEMPLATE),(ont_endprep,InputRole.CATALYST)),
                (ont_tailed,))
    ont_adapter=graph.input("ont_adapter","native barcode and motor sequencing adapters",
                            MaterialKind.ADAPTER,inferred=True,
                            source=SOURCE+" (kit-facing bases undisclosed)")
    ont_ligase=graph.input("ont_ligase","MinION adapter ligase",MaterialKind.CATALYST)
    ont_library=graph.product("ont_library","motor-adapted MinION tiled-amplicon library",
                              MaterialKind.LIBRARY)
    graph.react("ligate_ont","Ligate native barcode and motor adapters",
                StepKind.LIGATION,
                ((ont_tailed,InputRole.TEMPLATE),(ont_adapter,InputRole.STRUCTURAL),
                 (ont_ligase,InputRole.CATALYST)),(ont_library,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The named MinION kits do not disclose motor/barcode ligation-facing sequences."))
    graph.sequence(ont_library,SequencingPlan(
        limitation="MinION uses a motor-loaded adapter and no sequencing primer; the kit adapter sequence is proprietary."))
    miseq_endprep=graph.input("miseq_endprep","KAPA end-repair and dA-tailing mix",
                              MaterialKind.CATALYST)
    miseq_tailed=graph.product("miseq_tailed","end-repaired dA-tailed MiSeq amplicons",
                               MaterialKind.MOLECULAR)
    graph.react("prepare_miseq_ends","End-repair and dA-tail the MiSeq branch",
                StepKind.END_REPAIR,
                ((miseq_input,InputRole.TEMPLATE),(miseq_endprep,InputRole.CATALYST)),
                (miseq_tailed,))
    miseq_adapter=graph.input("miseq_adapter","SureSelectXT2 indexing adapter",
                              MaterialKind.ADAPTER,inferred=True,
                              source=SOURCE+" (adapter bases not printed)")
    miseq_ligase=graph.input("miseq_ligase","MiSeq-branch adapter ligase",
                             MaterialKind.CATALYST)
    miseq_adapted=graph.product("miseq_adapted","SureSelectXT2-adapted tiled amplicons",
                                MaterialKind.MOLECULAR,inferred=True,
                                source=SOURCE+" (adapter bases not printed)")
    graph.react("ligate_miseq","Ligate SureSelectXT2 indexing adapters",
                StepKind.LIGATION,
                ((miseq_tailed,InputRole.TEMPLATE),(miseq_adapter,InputRole.STRUCTURAL),
                 (miseq_ligase,InputRole.CATALYST)),(miseq_adapted,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The paper names SureSelectXT2 adapters without printing their ligation-facing geometry."))
    fwd=graph.input("miseq_forward","MiSeq library PCR forward primer",MaterialKind.OLIGO,
                    inferred=True,source=SOURCE+" (kit oligo not printed)")
    rev=graph.input("miseq_reverse","MiSeq indexed library PCR reverse primer",
                    MaterialKind.OLIGO,inferred=True,
                    source=SOURCE+" (kit oligo not printed)")
    pcr_polymerase=graph.input("miseq_polymerase","KAPA library PCR polymerase",
                               MaterialKind.CATALYST)
    miseq_library=graph.product("miseq_library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,
                                state=FINAL_LIBRARY,inferred=True,
                                source=SOURCE+" (inferred Illumina-compatible endpoint)")
    graph.react("miseq_pcr","Amplify the indexed MiSeq branch",StepKind.PCR,
                ((miseq_adapted,InputRole.TEMPLATE),(fwd,InputRole.STRUCTURAL),
                 (rev,InputRole.STRUCTURAL),(pcr_polymerase,InputRole.CATALYST)),
                (miseq_library,),geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY,reason="The named kit primers are not printed in the defining protocol."))
    graph.sequence(miseq_library,SequencingPlan(library=FINAL_LIBRARY,
                                                primers=tuple(SEQ_PRIMERS)))
    graph.declare_reviewed_scope((
        "viral_template","split_pools","pool1_forward","pool1_reverse",
        "pool2_forward","pool2_reverse",
        "amplicon_polymerase","amplify_odd","amplify_even","pool_amplicons",
        "split_platforms","ont_endprep","prepare_ont_ends","ont_adapter","ont_ligase",
        "ligate_ont","ont_library","miseq_endprep","prepare_miseq_ends",
        "miseq_adapter","miseq_ligase","ligate_miseq","miseq_forward",
        "miseq_reverse","miseq_polymerase","miseq_pcr","miseq_library"),SOURCE)
    return graph.require_publishable()
