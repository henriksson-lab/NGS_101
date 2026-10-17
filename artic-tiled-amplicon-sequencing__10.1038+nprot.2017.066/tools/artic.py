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
                            SourceClaim, SourceDiscrepancy,
                            StepKind)
from tiled_amplicon import alternating_scheme, scheme_rows

TITLE = "ARTIC / PrimalSeq — two-pool tiled viral amplicon sequencing"
NOTES = "01_artic-tiled-amplicon.html"
SOURCE = ('Defining protocol: <a href="https://doi.org/10.1038/nprot.2017.066">'
          'Quick et al., <i>Nature Protocols</i> (2017)</a>.')
SUMMARY = ("Overlapping amplicons tile a viral genome, but adjacent products alternate between "
           "two multiplex PCR pools so that overlapping primers cannot make preferential short products. "
           "The products then enter distinct MinION or MiSeq library branches.")
CAVEAT = ("The paper names SureSelectXT2 indexes (Agilent G9622A) and the KAPA Hyper Prep Kit but does not print the indexing-adapter bases. "
          "The MiSeq adapter geometry is therefore shown as inferred Illumina-compatible sequence; "
          "the two-pool amplicons and both platform branches are source-defined. The two pools stay separate through "
          "cleanup and are barcoded or indexed individually, as the protocol recommends; MinION users may instead combine pools "
          "under one barcode per sample.")

# Compact representative tiling path. Coordinates are diagram units, not a particular virus scheme.
SCHEME = alternating_scheme([(0, 12), (9, 21), (18, 30), (27, 39), (36, 48), (45, 57)])

# ZikaAsian pool 1, first tiled amplicon; Quick et al. Supplementary Table 1.
ZIKA_400_1_LEFT = "GACAGTTCGAGTTTGAAGCGAAAG"
ZIKA_400_1_RIGHT = "AGTATGCACTCCCACGTCTAGT"

FINAL_LIBRARY = Construct([
    # Single-index adapter: P5 and the Read 1 arm are one oligo, so their shared ACAC
    # occurs once (canonical universal adapter, illumina.TRUSEQ_P5_FULL).
    seg("P5", il.TRUSEQ_P5_FULL[:-len(il.TRUSEQ_READ1)], "p5", inferred=True),
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
SEQUENCING_INTRO = ("The protocol recommends paired 250-base MiSeq reads (MiSeq Reagent Kit v2, 500 cycle) for 400-base amplicons. "
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
        ], "Q5 Hot Start polymerase, 0.015 µM per primer and a five-minute anneal/extension support a highly multiplexed reaction."),
        ("MinION branch", ont.rows()+nanopore_ligation_scene().rows(),
         "After Ultra II end prep, each pool is recommended to receive its own native barcode (Blunt/TA ligase); barcoded pools are combined before adapter ligation. The paper developed the method with 2D EXP-NBD002/SQK-LSK208 and lists 1D EXP-NBD103/SQK-LSK108."),
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
    graph.react("split_pools","Split template between alternating primer pools (2.5 µl cDNA or up to 10 µl DNA each)",
                StepKind.SPLIT,((viral,InputRole.TEMPLATE),),
                (pool1_template,pool2_template),changes_molecule=False)
    polymerase=graph.input("amplicon_polymerase","Q5 Hot Start High-Fidelity DNA polymerase",
                           MaterialKind.CATALYST)
    tiled_template=Construct([seg("representative tiled viral target","X"*42,
                                  placeholder=True)],name="variable tiled viral targets")
    ampure=graph.input("ampure","AMPure XP beads (1.0×/0.8×/0.6× by amplicon length), elute 30 µl EB",
                       MaterialKind.SELECTION)
    scope=["viral_template","split_pools","amplicon_polymerase","ampure"]
    cleaned={}
    for pool,parity,template in ((1,"odd",pool1_template),(2,"even",pool2_template)):
        p=f"pool{pool}"
        printed=(" including source-printed ZIKA_400_1_LEFT" if pool==1 else "",
                 " including source-printed ZIKA_400_1_RIGHT" if pool==1 else "")
        fwd=graph.input(f"{p}_forward",f"{parity}-pool forward-primer members (0.015 µM each){printed[0]}",MaterialKind.OLIGO)
        rev=graph.input(f"{p}_reverse",f"{parity}-pool reverse-primer members (0.015 µM each){printed[1]}",MaterialKind.OLIGO)
        amp=graph.product(f"{parity}_amplicons",f"{parity} tiled viral amplicons",MaterialKind.MOLECULAR)
        graph.react(f"amplify_{parity}",f"Pool {pool} PCR: 98 °C 30 s; cycles 2–40 of 98 °C 15 s, 65 °C 5 min",StepKind.PCR,
                    ((template,InputRole.TEMPLATE),(fwd,InputRole.STRUCTURAL),
                     (rev,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),
                    (amp,),geometry=PCRGeometry.unavailable(
                        tiled_template,reason="Each pool contains many source-defined inward primer pairs rather than one PCR geometry."))
        clean=graph.product(f"{parity}_cleaned",f"bead-cleaned {parity} amplicons (kept separate)",MaterialKind.MOLECULAR)
        waste=graph.product(f"{parity}_bead_waste",f"pool {pool} PCR supernatant and washes",MaterialKind.WASTE)
        graph.react(f"cleanup_{parity}",f"AMPure cleanup of pool {pool} (Step 9)",StepKind.CLEANUP,
                    ((amp,InputRole.TEMPLATE),(ampure,InputRole.SELECTION)),(clean,waste))
        graph.discard(waste,"unbound PCR reaction removed by the bead cleanup")
        ont=graph.product(f"ont_{parity}",f"pool {pool} amplicons for the MinION branch",MaterialKind.MOLECULAR)
        miseq=graph.product(f"miseq_{parity}",f"pool {pool} amplicons for the MiSeq branch (50 ng)",MaterialKind.MOLECULAR)
        graph.react(f"split_{parity}",f"Assign pool {pool} amplicons to a platform branch",StepKind.SPLIT,
                    ((clean,InputRole.TEMPLATE),),(ont,miseq),changes_molecule=False)
        cleaned[pool]=(ont,miseq)
        scope+= [f"{p}_forward",f"{p}_reverse",f"amplify_{parity}",f"cleanup_{parity}",f"split_{parity}"]

    # MinION (Step 12A): per-pool end prep and native barcodes, then pooling.
    ont_endprep=graph.input("ont_endprep","NEBNext Ultra II End Prep (20 °C 5 min, 65 °C 5 min)",
                            MaterialKind.CATALYST)
    barcode=graph.input("native_barcodes","native barcodes NB01–NB12 (one per pool recommended)",
                        MaterialKind.ADAPTER,inferred=True,
                        source=SOURCE+" (kit-facing bases undisclosed)")
    blunt_ta=graph.input("blunt_ta_ligase","NEB Blunt/TA Ligase Master Mix",MaterialKind.CATALYST)
    barcoded=[]
    for pool,(ont,_miseq) in cleaned.items():
        tailed=graph.product(f"ont_tailed_{pool}",f"end-repaired dA-tailed pool {pool} amplicons",
                             MaterialKind.MOLECULAR)
        graph.react(f"prepare_ont_ends_{pool}",f"End-repair and dA-tail pool {pool}, then SPRI cleanup",
                    StepKind.END_REPAIR,
                    ((ont,InputRole.TEMPLATE),(ont_endprep,InputRole.CATALYST)),(tailed,))
        bc=graph.product(f"ont_barcoded_{pool}",f"native-barcoded pool {pool} amplicons",
                         MaterialKind.MOLECULAR,inferred=True,
                         source=SOURCE+" (barcode bases undisclosed)")
        graph.react(f"ligate_barcode_{pool}",f"Ligate a native barcode to pool {pool} (20 °C 10 min, 65 °C 10 min)",
                    StepKind.LIGATION,
                    ((tailed,InputRole.TEMPLATE),(barcode,InputRole.STRUCTURAL),
                     (blunt_ta,InputRole.CATALYST)),(bc,),
                    geometry=AdapterLigationGeometry.unavailable(
                        "The native barcoding kit does not disclose barcode ligation-facing sequences or overhangs."))
        barcoded.append(bc)
        scope+=[f"prepare_ont_ends_{pool}",f"ligate_barcode_{pool}"]
    ont_pooled=graph.product("ont_pooled","pooled native-barcoded amplicons",MaterialKind.MOLECULAR,
                             inferred=True,source=SOURCE+" (barcode bases undisclosed)")
    graph.react("pool_barcoded","Pool barcode ligations and SPRI-clean (elute 30 µl)",StepKind.POOL,
                tuple((b,InputRole.TEMPLATE) for b in barcoded),(ont_pooled,),changes_molecule=False)
    ont_adapter=graph.input("ont_adapter","native-barcoding-kit sequencing adapter (motor)",
                            MaterialKind.ADAPTER,inferred=True,
                            source=SOURCE+" (kit-facing bases undisclosed)")
    ont_library=graph.product("ont_library","motor-adapted MinION tiled-amplicon library",
                              MaterialKind.LIBRARY)
    graph.react("ligate_ont","Ligate the barcoding adapter per kit protocol, with Blunt/TA ligase",
                StepKind.LIGATION,
                ((ont_pooled,InputRole.TEMPLATE),(ont_adapter,InputRole.STRUCTURAL),
                 (blunt_ta,InputRole.CATALYST)),(ont_library,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The named MinION kits do not disclose motor adapter ligation-facing sequences."))
    graph.sequence(ont_library,SequencingPlan(
        limitation="MinION uses a motor-loaded adapter and no sequencing primer; the kit adapter sequence is proprietary."))
    scope+=["ont_endprep","native_barcodes","blunt_ta_ligase","pool_barcoded","ont_adapter",
            "ligate_ont","ont_library"]

    # MiSeq (Step 12B): one KAPA Hyper library per pool, then equimolar pooling.
    miseq_endprep=graph.input("miseq_endprep","KAPA Hyper end-repair and dA-tailing mix",
                              MaterialKind.CATALYST)
    miseq_adapter=graph.input("miseq_adapter","SureSelectXT2 indexing adapter (one index per pool library)",
                              MaterialKind.ADAPTER,inferred=True,
                              source=SOURCE+" (adapter bases not printed)")
    miseq_ligase=graph.input("miseq_ligase","KAPA Hyper ligase",MaterialKind.CATALYST)
    fwd=graph.input("miseq_forward","KAPA library amplification primer",MaterialKind.OLIGO,
                    inferred=True,source=SOURCE+" (kit oligo not printed)")
    rev=graph.input("miseq_reverse","KAPA library amplification primer (partner)",
                    MaterialKind.OLIGO,inferred=True,
                    source=SOURCE+" (kit oligo not printed)")
    pcr_polymerase=graph.input("miseq_polymerase","KAPA library PCR polymerase",
                               MaterialKind.CATALYST)
    spri=graph.input("miseq_spri","0.8× SPRI post-amplification cleanup",MaterialKind.SELECTION)
    pool_libraries=[]
    for pool,(_ont,miseq) in cleaned.items():
        tailed=graph.product(f"miseq_tailed_{pool}",f"end-repaired dA-tailed pool {pool} amplicons",
                             MaterialKind.MOLECULAR)
        graph.react(f"prepare_miseq_ends_{pool}",f"KAPA Hyper end-repair and dA-tailing of pool {pool}",
                    StepKind.END_REPAIR,
                    ((miseq,InputRole.TEMPLATE),(miseq_endprep,InputRole.CATALYST)),(tailed,))
        adapted=graph.product(f"miseq_adapted_{pool}",f"SureSelectXT2-adapted pool {pool} amplicons",
                              MaterialKind.MOLECULAR,inferred=True,
                              source=SOURCE+" (adapter bases not printed)")
        graph.react(f"ligate_miseq_{pool}",f"Ligate SureSelectXT2 indexing adapters to pool {pool} in place of KAPA adapters",
                    StepKind.LIGATION,
                    ((tailed,InputRole.TEMPLATE),(miseq_adapter,InputRole.STRUCTURAL),
                     (miseq_ligase,InputRole.CATALYST)),(adapted,),
                    geometry=AdapterLigationGeometry.unavailable(
                        "The paper names SureSelectXT2 adapters without printing their ligation-facing geometry."))
        amplified=graph.product(f"miseq_amplified_{pool}",f"amplified pool {pool} library",
                                MaterialKind.MOLECULAR,state=FINAL_LIBRARY,inferred=True,
                                source=SOURCE+" (inferred Illumina-compatible endpoint)")
        graph.react(f"miseq_pcr_{pool}",f"Amplify the pool {pool} library per KAPA Hyper instructions",StepKind.PCR,
                    ((adapted,InputRole.TEMPLATE),(fwd,InputRole.STRUCTURAL),
                     (rev,InputRole.STRUCTURAL),(pcr_polymerase,InputRole.CATALYST)),
                    (amplified,),geometry=PCRGeometry.unavailable(
                        FINAL_LIBRARY,reason="The named kit primers are not printed in the defining protocol."))
        lib=graph.product(f"miseq_library_{pool}",f"0.8×-cleaned pool {pool} library",
                          MaterialKind.MOLECULAR,state=FINAL_LIBRARY,inferred=True,
                          source=SOURCE+" (inferred Illumina-compatible endpoint)")
        dimer=graph.product(f"adapter_dimer_{pool}",f"candidate adapter dimers and small fragments (pool {pool})",
                            MaterialKind.WASTE,inferred=True,
                            source=SOURCE+" (0.8× cleanup is stated to remove potential adaptor–dimers)")
        graph.react(f"spri_{pool}",f"0.8× SPRI cleanup of pool {pool} instead of 1×",StepKind.CLEANUP,
                    ((amplified,InputRole.TEMPLATE),(spri,InputRole.SELECTION)),(lib,dimer))
        graph.discard(dimer,"0.8× SPRI supernatant removes potential adaptor–dimers")
        pool_libraries.append(lib)
        scope+=[f"prepare_miseq_ends_{pool}",f"ligate_miseq_{pool}",f"miseq_pcr_{pool}",f"spri_{pool}"]
    miseq_library=graph.product("miseq_library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,
                                state=FINAL_LIBRARY,inferred=True,
                                source=SOURCE+" (inferred Illumina-compatible endpoint)")
    graph.react("pool_miseq","TapeStation QC, KAPA qPCR quantification and equimolar pooling",StepKind.POOL,
                tuple((l,InputRole.TEMPLATE) for l in pool_libraries),(miseq_library,),
                changes_molecule=False)
    graph.sequence(miseq_library,SequencingPlan(library=FINAL_LIBRARY,
                                                primers=tuple(SEQ_PRIMERS)))
    graph.annotate(SourceDiscrepancy(
        "amplify_odd","Published multiplex-PCR denaturation temperature differs",
        (SourceClaim("Experimental design, Primal Scheme primer design (published PDF p.1263)",
                     "2-step protocol: 95 °C denaturation, 65 °C combined annealing and extension"),
         SourceClaim("Procedure Step 8 cycling table (published PDF p.1269)",
                     "98 °C 30 s, then 98 °C 15 s / 65 °C 5 min")),
        "The model uses the Step 8 cycling table (98 °C), the executable procedure.",
        "Both statements are in the published Nature Protocols article.",
        "The article does not say whether 95 °C in the design rationale is a generic value or an alternative."))
    scope+=["miseq_endprep","miseq_adapter","miseq_ligase","miseq_forward","miseq_reverse",
            "miseq_polymerase","miseq_spri","pool_miseq","miseq_library"]
    graph.declare_reviewed_scope(tuple(scope),SOURCE)
    return graph.require_publishable()
