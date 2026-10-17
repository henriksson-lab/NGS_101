"""SureSelect XT HS2 target-capture model."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from adapter_ligation import AdapterLigationGeometry
from batch_ngs import seg
from chemdraw import Construct, MolecularState, Row, Scene, Workflow
from endprep import dA_tailed_scene, repair_and_dA_tail
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from batch_ngs import index_window, seg
import illumina as il
import seqprimers as sp
from chemdraw import complement_segments, feature

def hybrid_capture_scene() -> Scene:
    """An indexed library hybridized to one complementary biotinylated RNA bait."""
    library = [seg("left adapter", "X" * 12, placeholder=True),
               seg("off-target flank", "N" * 8, placeholder=True),
               seg("target interval", "ACGTTGCAACGTTGCAACGT"),
               seg("right flank / adapter", "Y" * 16, placeholder=True)]
    bait = [seg("biotinylated RNA bait", "ACGTTGCAACGTTGCAACGT")]
    sc = Scene()
    sc.strand("library", library, label="indexed library")
    sc.anneal("bait", complement_segments(bait), to="library",
              pair=("biotinylated RNA bait'", "target interval"),
              label="5'-biotin RNA bait")
    sc.mark("bait", "biotinylated RNA bait'", "streptavidin capture")
    return sc


G9957 = "Agilent G9957-90000 Version B0 (May 2026)"
MBC_NOTE = (G9957 + " Step 4 'Process and analyze the reads', MBC-tagged libraries: inline MBCs "
            "are added to both insert ends; one approach trims 5 bases at the start of each "
            "read, takes the first 3 as MBC and discards the remaining 2.")


def _mbc_segments(side: str) -> list:
    """The five read-start bases at one insert end: 3-nt MBC and 2 discarded bases.

    Written as they appear on the top strand, so the right end is mirrored: Read 2 meets
    the MBC first, then the two discarded bases, then the insert.
    """
    mbc = seg(f"{side} MBC", ("U" if side == "left" else "V") * 3, "umi", placeholder=True,
              inferred=True,
              feature=feature(f"molecular_barcode_{side}", "umi", "random",
                              group="molecular_barcode", part=side, note=MBC_NOTE))
    spacer = seg(f"{side} non-MBC read-start bases", "SS", placeholder=True, inferred=True, length_bp=2)
    return [mbc, spacer] if side == "left" else [spacer, mbc]


def sureselect_library():
    """HS2 MBC-tagged form: Read 1 and Read 2 each start with 5 adaptor bases (3 MBC + 2).

    The adaptor arms are not published; canonical TruSeq arms are inferred for orientation.
    The insert's dA/T ligation junction is not drawn at the insert edge because the manual
    assigns the first five read bases to the adaptor tag region before the insert.
    """
    lib = Construct([
        seg("P5", il.P5, "p5", inferred=True),
        index_window("i5", "J", 8, "sample_index_i5", inferred=True),
        seg("Read 1 arm", il.TRUSEQ_READ1, "r1", inferred=True),
        *_mbc_segments("left"),
        seg("captured target insert", "N" * 34, placeholder=True),
        *_mbc_segments("right"),
        # Last base of the canonical TruSeq Read 2 primer site; here an inferred adaptor-arm
        # base, not the insert dA tail (the manual puts five tag bases before the insert).
        seg("Read 2 site terminal base", "A", "r2", inferred=True),
        seg("Index 1 / Read 2 arm", il.INDEX1_PRIMER, "r2", inferred=True),
        index_window("i7 reverse complement", "I", 8, "sample_index_i7", inferred=True),
        seg("P7 reverse complement", il.P7_RC, "p7", inferred=True),
    ], name="SureSelect XT HS2 captured library (MBC-tagged)")
    primers = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["I2"], sp.TRUSEQ["R2"])
    problems = sp.verify(lib, primers, required_roles=tuple(p.role for p in primers))
    if problems:
        raise ValueError("; ".join(problems))
    return lib, primers


TITLE = "SureSelect XT HS2 DNA target enrichment"
NOTES = "01_sureselect-xt-hs2.html"
SOURCE = 'Agilent <a href="https://www.agilent.com/cs/library/usermanuals/public/G9957-90000.pdf">SureSelect XT HS2 DNA Kits, G9957-90000 B0 (2026)</a>.'
SUMMARY = "A dual-indexed, optionally duplex-MBC-tagged Illumina library is pooled, hybridized to biotinylated RNA baits, recovered on streptavidin beads and amplified after capture."
CAVEAT = "The B0 guide (G9957-90000) publishes the MBC read layout (five read-start bases at each end, of which one suggested approach keeps the first three as MBC and discards two) and 384 paired eight-base indexes, but not the HS2 adaptor, indexing-primer or post-capture primer bases. Canonical TruSeq arms in the final construct are inferred for orientation and primer landing, not established kit oligos; the selected index pair is unknown. The MBC-free adaptor branch is not drawn."
FINAL_LIBRARY, SEQ_PRIMERS = sureselect_library()
FINAL_CAPTION = "INFERRED — canonical TruSeq arms place the manual's five read-start bases at each insert end (3-nt MBC + 2 discarded bases) and the eight-base index windows for orientation only; Agilent's exact arm sequences are unresolved."
SEQUENCING_INTRO = "The manual specifies standard Illumina paired-end and dual-index primers with no custom sequencing primers, and 8-bp i7 and i5 index reads. Read 1/Read 2 length is user-chosen: Table 46 lists 2 × 100 or 2 × 150 bp runs and Table 47 gives 151 + 8 + 8 + 151 cycles as an example (one extra cycle per read)."

def workflow():
    insert=Construct([seg("sheared genomic insert","X"*40,placeholder=True)])
    w=Workflow(MolecularState("sheared genomic DNA",tuple(Scene.duplex(list(insert)).rows())))
    tailed=repair_and_dA_tail(insert)
    w.react("End-repair and add 3′ dA",dA_tailed_scene(tailed).rows(),name="dA-tailed genomic DNA",note="The manual explicitly specifies the end-repair/dA-tailing step.")
    geometry=AdapterLigationGeometry.unavailable("Agilent's current SureSelect XT HS2 manual specifies HS2 adaptor ligation after dA-tailing but does not publish the adaptor oligo or its terminal bases; a complementary T is not fabricated here.")
    product=Scene.duplex(list(FINAL_LIBRARY),label="indexed HS2 library",unpaired=tuple(s.name+"'" for s in FINAL_LIBRARY if s.is_role_token()))
    w.ligate_adapter("INFERRED — ligate the HS2 adaptor",geometry,product.rows(),name="HS2 adapter-ligated library",note="MBC-tagged adaptors put five read-start bases (3-nt MBC + 2 discarded) at both insert ends; the exact adapter arms remain undisclosed.")
    w.react("Hybridize the target interval to a biotinylated bait",hybrid_capture_scene().rows(),name="bait-bound target library",note="Libraries are pooled (8 or 16 indexes, 4 µg) and pre-incubated with the XT HS2 Blocker Mix before the probe is added.")
    w.react("Capture, wash and amplify",[
            Row(chunks=[("bait—library duplex + streptavidin bead → retained", None, False)]),
            Row(chunks=[("unbound library → stringent wash → discarded", None, False)]),
            Row(chunks=[("bead-bound target library → post-capture PCR", None, False)])],name="captured HS2 library",note="Affinity selection changes abundance, not adapter orientation.")
    return w


def protocol_graph():
    src = lambda where: f"{G9957}, {where}"
    graph=ProtocolGraph(TITLE,authority=SOURCE)
    dna=graph.input("genomic_dna","10–200 ng genomic DNA (intact or FFPE)",MaterialKind.MOLECULAR,primary=True)
    fragmentation=graph.input("fragmentation_condition",
                              "Covaris shearing (duration by read length) or SureSelect enzymatic fragmentation, 50 µL",
                              MaterialKind.CONDITION,source=src("Chapter 2 Step 2"))
    fragmented=graph.product("fragmented","fragmented genomic DNA",MaterialKind.MOLECULAR)
    graph.react("fragment","Fragment genomic DNA",StepKind.FRAGMENTATION,
                ((dna,InputRole.TEMPLATE),(fragmentation,InputRole.CONDITION)),
                (fragmented,))
    endprep=graph.input("endprep_mix","End Repair-A Tailing Buffer and Enzyme Mix",
                        MaterialKind.CATALYST,source=src("Chapter 3 Step 2, Table 19"))
    tailed=graph.product("tailed","end-repaired 3-prime-dA-tailed DNA",
                         MaterialKind.MOLECULAR)
    graph.react("endprep","Repair fragment ends and add 3-prime dA (20 °C 15 min, 72 °C 15 min)",
                StepKind.END_REPAIR,
                ((fragmented,InputRole.TEMPLATE),(endprep,InputRole.CATALYST)),
                (tailed,))
    adapter=graph.input("hs2_adapter","XT HS2 Adaptor Oligo Mix (MBC-tagged; the MBC-Free Adaptor Oligo Mix branch is not drawn)",
                        MaterialKind.ADAPTER,inferred=True,
                        source=SOURCE+" (complete oligo bases undisclosed)")
    ligase=graph.input("adapter_ligase","T4 DNA Ligase in Ligation Buffer",MaterialKind.CATALYST,
                       source=src("Chapter 3 Steps 1 and 3"))
    ligated=graph.product("ligated","HS2 adapter-ligated molecules",
                          MaterialKind.MOLECULAR,inferred=True,
                          source=SOURCE+" (hidden adapter sequence)")
    geometry=AdapterLigationGeometry.unavailable(
        "Agilent specifies ligation after dA-tailing but does not publish the HS2 adapter terminal bases.")
    graph.react("ligate_adapter","Ligate the HS2 adapter mix (20 °C 30 min)",StepKind.LIGATION,
                ((tailed,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),
                 (ligase,InputRole.CATALYST)),(ligated,),geometry=geometry)
    ampure=graph.input("ampure_xp","AMPure XP beads",MaterialKind.SELECTION,
                       source=src("Tables 22, 26 and 42"))
    ligated_clean=graph.product("ligated_clean","purified adapter-ligated molecules",
                                MaterialKind.MOLECULAR,inferred=True)
    lig_waste=graph.product("ligation_cleanup_waste","AMPure supernatant (unligated adaptor, reaction components)",
                            MaterialKind.WASTE)
    graph.react("ligation_cleanup","AMPure XP clean-up (80 µL beads, elute 35 µL)",StepKind.CLEANUP,
                ((ligated,InputRole.TEMPLATE),(ampure,InputRole.SELECTION)),
                (ligated_clean,lig_waste))
    graph.discard(lig_waste,"bead supernatant removed")
    i5=graph.input("i5_primer","SureSelect XT HS2 Index Primer Pair, P5 member (8-nt i5)",MaterialKind.OLIGO,
                   source=src("Tables 61–68 (index bases only)"))
    i7=graph.input("i7_primer","SureSelect XT HS2 Index Primer Pair, P7 member (8-nt i7)",MaterialKind.OLIGO,
                   source=src("Tables 61–68 (index bases only)"))
    polymerase=graph.input("pcr_polymerase","Herculase II Fusion DNA Polymerase",MaterialKind.CATALYST,
                           source=src("Chapter 3 Step 5 and Chapter 4 Step 5"))
    precapture=graph.product("precapture","indexed pre-capture HS2 library",
                             MaterialKind.MOLECULAR,state=FINAL_LIBRARY)
    graph.react("precapture_pcr","Index the library by PCR (98 °C 2 min; 8–14 cycles by input amount and quality; 72 °C 5 min), then AMPure (50 µL)",StepKind.PCR,
                ((ligated_clean,InputRole.TEMPLATE),(i5,InputRole.STRUCTURAL),
                 (i7,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),
                (precapture,),geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY,reason="Index sequences are published, but the hidden adapter oligos prevent an exact precursor geometry."))
    others=graph.input("other_indexed_libraries","other uniquely dual-indexed HS2 libraries",
                       MaterialKind.MOLECULAR,source=src("Chapter 4 Step 1, Table 30"))
    pool=graph.product("capture_pool","pre-capture pool of 8 or 16 indexed libraries (4 µg total, dried to 4 µL)",
                       MaterialKind.MOLECULAR,state=FINAL_LIBRARY)
    graph.react("pool","Pool indexed libraries before capture",StepKind.POOL,
                ((precapture,InputRole.TEMPLATE),(others,InputRole.TEMPLATE)),(pool,),
                changes_molecule=False)
    blocker=graph.input("blocker_mix","SureSelect XT HS2 Blocker Mix (95 °C 5 min, then held at the hybridization temperature)",
                        MaterialKind.MOLECULAR,source=src("Chapter 4 Step 2"))
    bait=graph.input("rna_baits","target-specific SureSelect probe (biotinylated RNA baits) with RNase Block in hybridization buffer",
                     MaterialKind.MOLECULAR)
    hybridized=graph.product("hybridized","bait-bound target library molecules",
                             MaterialKind.MOLECULAR,
                             state=tuple(hybrid_capture_scene().rows()))
    graph.react("hybridize","Hybridize target intervals to the probe, 65 °C or 67.5 °C for 16–24 h",
                StepKind.HYBRIDIZATION,
                ((pool,InputRole.TEMPLATE),(blocker,InputRole.STRUCTURAL),
                 (bait,InputRole.STRUCTURAL)),
                (hybridized,))
    beads=graph.input("streptavidin","SureSelect Streptavidin Beads; Wash Buffer 1, then 6 × Wash Buffer 2 at 70 °C",
                      MaterialKind.SELECTION,source=src("Chapter 4 Steps 3–4"))
    captured=graph.product("captured","bead-bound target-enriched library",
                           MaterialKind.MOLECULAR)
    unbound=graph.product("unbound","unbound off-target library molecules",
                          MaterialKind.WASTE)
    graph.react("capture","Capture bait-bound molecules (30 min, room temperature) and wash off unbound library",
                StepKind.AFFINITY_CAPTURE,
                ((hybridized,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),
                (captured,unbound))
    graph.discard(unbound,"unbound molecules are removed with the supernatant and stringent washes")
    # The kit supplies one "Post-Capture Primer Mix"; its two end-specific members are
    # implied by exponential amplification of both library ends, not listed by Agilent.
    post_f=graph.input("post_forward","SureSelect Post-Capture Primer Mix: P5-end primer (bases not published)",
                       MaterialKind.OLIGO,source=src("Chapter 4 Step 5"))
    post_r=graph.input("post_reverse","SureSelect Post-Capture Primer Mix: P7-end primer (bases not published)",
                       MaterialKind.OLIGO,source=src("Chapter 4 Step 5"))
    amplified=graph.product("amplified","post-capture PCR product in supernatant",
                            MaterialKind.MOLECULAR,state=FINAL_LIBRARY)
    spent_beads=graph.product("spent_streptavidin_beads","streptavidin beads with bound probe",
                              MaterialKind.WASTE)
    graph.react("postcapture_pcr","Amplify the captured library on the beads (98 °C 2 min; 10–16 cycles by probe size; 72 °C 5 min), then remove the beads on a magnet",StepKind.PCR,
                ((captured,InputRole.TEMPLATE),(post_f,InputRole.STRUCTURAL),
                 (post_r,InputRole.STRUCTURAL),
                 (polymerase,InputRole.CATALYST)),
                (amplified,spent_beads),geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY,reason="Post-capture PCR preserves the final construct, but the manual does not disclose the Post-Capture Primer Mix bases."))
    graph.discard(spent_beads,"beads removed on the magnetic stand; the supernatant is kept")
    library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,
                          state=FINAL_LIBRARY)
    final_waste=graph.product("final_cleanup_waste","AMPure supernatant",MaterialKind.WASTE)
    graph.react("final_cleanup","AMPure XP clean-up of the final pool (50 µL beads)",StepKind.CLEANUP,
                ((amplified,InputRole.TEMPLATE),(ampure,InputRole.SELECTION)),
                (library,final_waste))
    graph.discard(final_waste,"bead supernatant removed")
    graph.sequence(library,SequencingPlan(library=FINAL_LIBRARY,
                                          primers=tuple(SEQ_PRIMERS),
                                          read_lengths={"Index 1 (i7)":8,"Index 2 (i5)":8}))
    graph.declare_reviewed_scope((
        "genomic_dna","fragmentation_condition","fragment","endprep_mix","endprep",
        "hs2_adapter","adapter_ligase","ligate_adapter","ampure_xp","ligation_cleanup",
        "ligation_cleanup_waste","i5_primer","i7_primer","pcr_polymerase","precapture_pcr",
        "other_indexed_libraries","pool","blocker_mix","rna_baits","hybridize","streptavidin",
        "capture","post_forward","post_reverse","postcapture_pcr","spent_streptavidin_beads",
        "final_cleanup","final_cleanup_waste","library"),SOURCE)
    return graph.require_publishable()
