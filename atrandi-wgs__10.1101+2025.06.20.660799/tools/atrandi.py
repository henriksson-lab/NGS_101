"""
Segment and oligo definitions for the Atrandi / PTA single-cell WGS chemistry.

Everything the page draws derives from this file. Sources for each sequence are in the
notes: 01_barcoding_kit.md, 02_library_prep.md, 04_nebnext_illumina.md,
05_barcode_cassette_model.md.

Evidence marking mirrors the notes: segments carrying inferred=True are guesses and render
inside <inf>. Placeholders (barcodes, linkers, insert) are stand-ins, never complemented.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import illumina as il
import pta as pta_chem
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, feature, revcomp
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_graph import (AnnotationSlot, InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, SourceClaim, SourceDiscrepancy, StepKind)

TITLE = "Atrandi semi-permeable capsules + PTA — single-microbe WGS"
SOURCE = 'Defining source: <a href="https://doi.org/10.1101/2025.06.20.660799">Gourlé et al., preprint v2</a>; local PTA implementation with Atrandi DGPM02323198001 and DGPM02323206001, V3 (3 November 2023) source comparison.'

# ------------------------------------------------------- canonical Illumina (by reference)
# The sequences live once, in lib/illumina.py; these are local aliases, never copies.

P5_SEQ = il.P5
P7_SEQ = il.P7
TRUSEQ_R1 = il.TRUSEQ_READ1
TRUSEQ_R2 = il.TRUSEQ_READ2
TRIM_SEEN_IN_R1 = il.TRIM_SEEN_IN_READ1
TRIM_SEEN_IN_R2 = il.TRIM_SEEN_IN_READ2

# ------------------------------------------- Atrandi ligation adapter (DGPM...206001)
# As printed by Atrandi (citation only -- the values below are derived, and the selftest
# pins the derivation to this text):
#   Top:    /5Phos/GATCGGAAGAGCGTCGTGTAGGGAAAGAGTG*T
#   Bottom: /5AmMC6/GCTCTTCCGATCT
# i.e. the top strand is Illumina's Read-2 trim sequence minus the insert's dA, and the
# bottom strand is the last 13 nt of TruSeq Read 1.
LIGADAPT_TOP = TRIM_SEEN_IN_R2[1:]             # 32 nt
LIGADAPT_STEM = LIGADAPT_TOP[:len(il.STEM)]    # 12 bp duplex stem == il.STEM
LIGADAPT_ARM = LIGADAPT_TOP[len(il.STEM):]     # 20 nt single-stranded 3' arm
LIGADAPT_BOT = TRUSEQ_R1[-13:]                 # 13 nt: 12 paired + 1-nt 3'-T overhang

# ------------------------------------------------ Atrandi indexing primers
# Atrandi's i7 primer = P7 + 6-nt index + the first 27 nt of TruSeq Read 2; its P5 primer =
# P5 + TruSeq Read 1 up to the length of the adapter arm (the two share P5's final ACAC).
ATRANDI_I7_INDEX = "AACCTG"                              # 6 nt, i7 only
ATRANDI_I7_ANNEAL = TRUSEQ_R2[:27]                       # 27 nt
ATRANDI_I7 = P7_SEQ + ATRANDI_I7_INDEX + ATRANDI_I7_ANNEAL
ATRANDI_P5_ANNEAL = TRUSEQ_R1[:len(LIGADAPT_ARM)]        # 20 nt == length of LIGADAPT_ARM
_P5_R1_OVERLAP = 4                                       # "ACAC", end of P5 = start of Read 1
ATRANDI_P5 = P5_SEQ[:-_P5_R1_OVERLAP] + ATRANDI_P5_ANNEAL   # 45 nt
ATRANDI_P5_TAIL = ATRANDI_P5[: len(ATRANDI_P5) - len(ATRANDI_P5_ANNEAL)]  # 25 nt 5' overhang

# ------------------------------------- our indexing primers (as ordered, ref/our_index_primers.tsv)
# Used for both our MDA/PTA scWGS and florian-PTA-rnaseq. Atrandi's i7 design with the 6-nt
# AACCTG swapped for a 10-nt IDT UDP i7 index; the oligo carries the revcomp of the
# sample-sheet i7. P5 is Atrandi's, unchanged. Both end in a 3'-terminal phosphorothioate.
# The selftest checks every modelled oligo against the ordered TSV.
OUR_I7_INDEX = {                                         # index as it sits in the oligo
    "UDP0005": "TAATCTCGTC",
    "UDP0006": "GCGCGATGTT",
    "UDP0007": "AGAGCACTAG",
    "UDP0008": "TGCCTTGATC",
}
OUR_I7_INDEX_LEN = 10
OUR_I7 = {k: P7_SEQ + v + ATRANDI_I7_ANNEAL for k, v in OUR_I7_INDEX.items()}
OUR_P5 = ATRANDI_P5

BARCODE_LEN = 8
LINKER_LEN = 4

_D_ARM_TAIL = TRUSEQ_R2[len(ATRANDI_I7_ANNEAL):]          # "CCGATCT", 7 nt


def _seg(name, top, tag=None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


def ligation_adapter() -> Scene:
    """The source-known V3 adapter, shared by its graph input and detail view."""
    sc = Scene()
    sc.strand("top", [_seg("stem", LIGADAPT_STEM, "s5"),
                      _seg("arm", LIGADAPT_ARM, "s5",
                           note="V3 p4: terminal phosphorothioate linkage before the last base")],
              label="", mod5="/5Phos/")
    sc.anneal("bottom", [_seg("stem'", LIGADAPT_BOT[:-1], "s5"),
                          _seg("T", LIGADAPT_BOT[-1:])],
              to="top", pair=("stem'", "stem"), label="", mod5="/5AmMC6/")
    sc.mark("bottom", "T", "3'-T overhang (the ligation end)")
    sc.mark("top", "stem", "12-bp stem")
    sc.mark("top", "arm", "20-nt single-stranded arm = the P5 primer landing site")
    return sc


def read2_arm_segments() -> list[Segment]:
    """The Read-2 / P7 landing site that the round-D cassette must present."""
    return [
        _seg("TruSeq Read 2", ATRANDI_I7_ANNEAL, "t7"),
        _seg("", _D_ARM_TAIL, "t7", inferred=True,
             note="required for the standard Read 2 primer to bind and extend"),
    ]


def barcode_block() -> list[Segment]:
    """D - linker - C - linker - B - linker - A, as read in R2."""
    out: list[Segment] = []
    for i, letter in enumerate("DCBA"):
        out.append(_seg(
            f"BC-{letter}", letter * BARCODE_LEN, "cbc", placeholder=True,
            feature=feature(
                f"cell_barcode_{letter.lower()}", "cell_barcode", "combinatorial",
                group="cell_barcode", part=f"round {letter}",
                note="Atrandi split-pool ligation barcode",
            ),
        ))
        if i < 3:
            out.append(_seg(f"{letter}-{chr(ord(letter) - 1)} linker",
                            "L" * LINKER_LEN, "r2", placeholder=True, inferred=True,
                            length_bp=LINKER_LEN,
                            note="4-nt cohesive overhang in this model; exact bases not established by this declaration"))
    return out


def dA_junction() -> Segment:
    """The FS-fragment's dA tail opposite the adapter's 3'-T overhang (P5 side). A real base
    pair, drawn explicitly so A + LIGADAPT_TOP == Illumina's Read-2 trim sequence and the
    TruSeq Read 1 primer's 3'-T has its partner (it was formerly hidden in the insert)."""
    return _seg("dA", "A", None, note="insert's dA tail / adapter's 3'-T overhang")


def final_library() -> Construct:
    """The finished Illumina library, top strand 5'->3' from the P7 end (Atrandi Fig. 2)."""
    segs: list[Segment] = [
        _seg("Illumina P7", P7_SEQ, "p7"),
        _seg("i7", ATRANDI_I7_INDEX, None,
             feature=feature("sample_i7", "sample_index", "fixed")),
        *read2_arm_segments(),
        *barcode_block(),
        _seg("ligation junction", "T", None, placeholder=True, inferred=True,
             bottom="A", length_bp=1,
             note="dA/dT ligation junction; Bascet's '+1 to account for ligation'"),
        _seg("insert", "XXX...XXX", None, placeholder=True),
        dA_junction(),
        _seg("TruSeq Read 1", LIGADAPT_TOP, "s5"),
        _seg("Illumina P5", revcomp(ATRANDI_P5_TAIL), "p5"),
    ]
    return Construct(segs, name="Atrandi/PTA scWGS final library")


def barcode_read() -> Construct:
    """The barcode-side read: barcode block onward. Offsets must match Bascet."""
    lib = final_library()
    names = [s.name for s in lib]
    start = names.index("BC-D")
    return Construct(lib.segments[start:], name="barcode-side read")


# ----------------------------------------------------------------- PTA amplicon (06_pta.md)
# Research-grade PTA: exo-resistant random hexamer 5'-NpNpNpNpsNpsN-3' (Dean 2002), extended
# by phi29 and terminated by an alpha-thio-dideoxynucleotide. The commercial kit uses a 6-9mer
# and an undisclosed terminator. 5' end is a hydroxyl, so it needs kinasing before ligation;
# 3' end is a dead end -- no 3'-OH, phosphorothioate linkage, exonuclease-resistant.

PTA_PRIMER_LEN = 6          # 6-9 in the commercial kit (BioSkryb TAS-007 Fig. 2B)
PTA_PRIMER_PS_LINKAGES = 2  # two 3'-terminal phosphorothioates (Thermo SO181 / Dean 2002)


def pta_amplicon() -> Construct:
    """One PTA amplicon as drawn on the page. Both termini are chemically distinctive."""
    shared = pta_chem.amplification_product(
        primer_len=PTA_PRIMER_LEN,
        primer_note="5'-OH; last two linkages are phosphorothioate -- needs T4 PNK before ligation")
    return Construct([
        _seg("random primer", shared.segments[0].top, "tso", placeholder=True, inferred=True,
             note="Illustrative research-grade primer; exact commercial ResolveDNA primer chemistry is not established here. " + shared.segments[0].note),
        _seg("amplified genomic DNA", "XXX...XXX", None, placeholder=True),
        _seg("terminator", shared.segments[-1].top, "w1", placeholder=True, inferred=True,
             note="Illustrative research-grade alpha-thio-ddN terminator; exact commercial ResolveDNA terminator identity is not established here"),
    ], name="PTA amplicon")


def pre_pcr() -> Construct:
    """The barcoded, adapter-ligated molecule as it enters indexing PCR (no P7/i7/P5 yet)."""
    segs = [*read2_arm_segments(), *barcode_block(),
            _seg("ligation junction", "T", None, placeholder=True, inferred=True,
                 bottom="A", length_bp=1),
            _seg("insert", "XXX...XXX", None, placeholder=True),
            dA_junction(),
            _seg("TruSeq Read 1", LIGADAPT_TOP, "s5")]
    return Construct(segs, name="pre-PCR molecule")


def adapter_both_ends() -> Construct:
    """A fragment that received the ligation adapter at both ends and no barcode."""
    return Construct([
        _seg("Read 1 site", TRUSEQ_R1, "s5"),
        _seg("insert", "XXX...XXX", None, placeholder=True),
        _seg("Read 1 site ", "A" + LIGADAPT_TOP, "s5"),
    ], name="adapter at both ends")


# ------------------------------------------------------- sequencing primers (lib/seqprimers.py)
# The run acquires both mates and i7, but not i5.  Which standard primers bind, their
# strands and their read directions are inferred centrally from the final sequence.
ACQUIRED_ROLES = ("Read 1", "Index 1 (i7)", "Read 2")
SEQ_PRIMERS = sp.standard_primers(final_library(), ACQUIRED_ROLES)


def protocol_graph():
    graph=ProtocolGraph(TITLE,authority=SOURCE)
    cells=graph.input("microbial_cells","single microbes in semi-permeable capsules",
                      MaterialKind.BIOLOGICAL,primary=True)
    lysis=graph.input("lysis_condition","enzymatic, proteinase-K and chaotropic lysis",
                      MaterialKind.CONDITION)
    genome=graph.product("encapsulated_dna","capsule-retained microbial genomic DNA",
                         MaterialKind.MOLECULAR)
    debris=graph.product("lysis_waste","solubilized cellular debris",MaterialKind.WASTE)
    graph.react("lyse","Lyse microbes while retaining genomic DNA in capsules",
                StepKind.EXTRACTION,
                ((cells,InputRole.TEMPLATE),(lysis,InputRole.CONDITION)),(genome,debris))
    graph.discard(debris,"solubilized cellular material passes out of the capsule")
    random_primer=graph.input("pta_primer","illustrative exo-resistant random PTA primer",MaterialKind.OLIGO, inferred=True)
    phi29=graph.input("phi29","illustrative phi29 strand-displacing polymerase",MaterialKind.CATALYST, inferred=True)
    amplified=graph.product("pta_product","terminated PTA amplicons",MaterialKind.MOLECULAR,
                            state=pta_amplicon(), inferred=True, source=SOURCE+"; research-grade molecular illustration, not verified commercial reagent identities")
    graph.react("pta_amplify","INFERRED molecular illustration: primary template-directed amplification",
                StepKind.EXTENSION,
                ((genome,InputRole.TEMPLATE),(random_primer,InputRole.STRUCTURAL),
                 (phi29,InputRole.CATALYST)),(amplified,))
    # One Atrandi End Prep reaction (DGPM02323198001 V3 steps 25-30: End Prep Enzyme Mix +
    # Reaction Buffer, 20 C 30 min then 65 C 30 min). The guide does not name the enzymatic
    # activities; the repair/dA-tail split below is this model's interpretation.
    endprep=graph.input("endprep_enzymes",
                        "Atrandi End Prep Enzyme Mix and Reaction Buffer (CRP-EPE1/EPB1; 20 °C 30 min, 65 °C 30 min)",
                        MaterialKind.CATALYST,
                        source="Atrandi DGPM02323198001 V3 steps 25–30")
    repaired=graph.product("repaired_pta","blunt phosphorylated PTA products",MaterialKind.MOLECULAR,
                           inferred=True,source=SOURCE+"; activities of the single End Prep mix are not named by the V3 guide")
    graph.react("end_repair","INFERRED activity split: blunt and phosphorylate extendable PTA ends",StepKind.END_REPAIR,
                ((amplified,InputRole.TEMPLATE),(endprep,InputRole.CATALYST)),(repaired,))
    tailed=graph.product("tailed_pta","dA-tailed PTA products",MaterialKind.MOLECULAR,
                         inferred=True,source=SOURCE+"; dA-tailing inferred from the same End Prep reaction and the T/A barcode junction")
    graph.react("da_tail","INFERRED activity split: add 3′ dA in the same End Prep reaction",StepKind.TAILING,
                ((repaired,InputRole.TEMPLATE),(endprep,InputRole.CATALYST)),(tailed,))
    ligase=graph.input("dna_ligase","Atrandi Ligation Enzyme and 10X Ligation Buffer (CRP-LGE1/LGB1; 22 °C 15 min per round, then STOP Buffer)",
                       MaterialKind.CATALYST,source="Atrandi DGPM02323198001 V3 steps 31–43")
    current=tailed
    for letter in "ABCD":
        cassette=graph.input("barcode_"+letter.lower(),"round "+letter+" barcode cassette",
            MaterialKind.ADAPTER,inferred=True,
            source=SOURCE+" (cassette strand/end geometry is not established in this declaration)")
        product=graph.product("barcoded_"+letter.lower(),"barcode rounds A–"+letter+" product",
            MaterialKind.MOLECULAR,inferred=True,
            source=SOURCE+" (cassette strand/end geometry is not established in this declaration)")
        graph.react("ligate_"+letter.lower(),"Ligate split-pool barcode round "+letter,
            StepKind.LIGATION,
            ((current,InputRole.TEMPLATE),(cassette,InputRole.STRUCTURAL),
             (ligase,InputRole.CATALYST)),(product,),
            geometry=AdapterLigationGeometry.unavailable(
                "Sequential ligation is source-known; the exact ordered cassette strands and end geometry have not been established in this model."))
        current=product
    release=graph.input("release_condition","capsule release reagent",MaterialKind.CONDITION)
    released=graph.product("released_dna","released four-round-barcoded DNA",MaterialKind.MOLECULAR)
    graph.react("release","Release barcoded DNA from the capsules",StepKind.EXTRACTION,
                ((current,InputRole.TEMPLATE),(release,InputRole.CONDITION)),(released,))
    # Preprint v2 Methods "Library preparation for sequencing": 0.8X AMPure XP; V3 library
    # guide section I steps 6-17: 80 uL beads to 100 uL, two 80% ethanol washes, 27 uL elution.
    beads=graph.input("ampure_xp","AMPure XP paramagnetic beads",MaterialKind.SELECTION,
                      source="Preprint v2 Methods; Atrandi DGPM02323206001 V3 sections I, IV and VI")
    purified=graph.product("purified_dna","0.8X AMPure XP-purified released DNA",MaterialKind.MOLECULAR)
    release_waste=graph.product("release_supernatant","0.8X bead supernatant and ethanol washes after release",
                                MaterialKind.WASTE)
    graph.react("release_cleanup","Purify released DNA with 0.8X AMPure XP",StepKind.CLEANUP,
                ((released,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),(purified,release_waste),
                source="Preprint v2 Methods; Atrandi DGPM02323206001 V3 section I steps 6–17")
    graph.discard(release_waste,"bead supernatant and washes are removed")
    fs=graph.input("fs_mix","NEBNext Ultra II FS fragmentation mix",MaterialKind.CATALYST)
    fragments=graph.product("fs_fragments","barcoded 300–700-bp DNA fragments",MaterialKind.MOLECULAR)
    fs_waste=graph.product("fs_waste","distal fragmentation products without the barcode end",
                           MaterialKind.WASTE)
    graph.react("fragment","Fragment released barcoded DNA",StepKind.FRAGMENTATION,
                ((purified,InputRole.TEMPLATE),(fs,InputRole.CATALYST)),(fragments,fs_waste))
    graph.discard(fs_waste,"only barcode-proximal fragments can form the intended library")
    fs_tailed=graph.product("fs_tailed","end-repaired and dA-tailed barcode-proximal fragments",
                            MaterialKind.MOLECULAR)
    graph.react("fs_endprep","End-repair and dA-tail the FS fragments",StepKind.END_REPAIR,
                ((fragments,InputRole.TEMPLATE),(fs,InputRole.CATALYST)),(fs_tailed,))
    lig_adapter=graph.input("ligation_adapter","Atrandi T-overhang ligation adapter",
                            MaterialKind.ADAPTER, state=ligation_adapter(),
                            source="Atrandi DGPM02323206001 V3 p4 adapter table")
    ultra_lig=graph.input("ultraii_ligation",
        "NEBNext Ultra II Ligation Master Mix and Ligation Enhancer (20 °C 15 min; 2.5 µL of 1.5 µM adapter)",
        MaterialKind.CATALYST,source="Atrandi DGPM02323206001 V3 section III")
    ligated=graph.product("ligated",pre_pcr().name+" in ligation mix",MaterialKind.MOLECULAR,
                          state=pre_pcr(),inferred=True,source=SOURCE+" (barcode-side arm inferred)")
    graph.react("adapter_ligation","Ligate the non-barcode sequencing adapter",
                StepKind.LIGATION,
                ((fs_tailed,InputRole.TEMPLATE),(lig_adapter,InputRole.STRUCTURAL),
                 (ultra_lig,InputRole.CATALYST)),(ligated,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The ligation adapter is source-known; the opposite barcode-cassette end is not established in this model."),
                source="Atrandi DGPM02323206001 V3 section III")
    pre=graph.product("pre_pcr",pre_pcr().name,MaterialKind.MOLECULAR,state=pre_pcr(),
                      inferred=True,source=SOURCE+" (barcode-side arm inferred)")
    lig_waste=graph.product("ligation_supernatant","0.8X bead supernatant and ethanol washes after ligation",
                            MaterialKind.WASTE)
    graph.react("ligation_cleanup","Purify the ligation with 0.8X AMPure XP",StepKind.CLEANUP,
                ((ligated,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),(pre,lig_waste),
                source="Atrandi DGPM02323206001 V3 section IV")
    graph.discard(lig_waste,"bead supernatant and washes are removed; their composition is not characterized by the source")
    p5=graph.input("p5_primer","Atrandi V3 example P5 library primer",MaterialKind.OLIGO,
        state=Construct([_seg("P5 primer", ATRANDI_P5,
            note="V3 p4: terminal phosphorothioate linkage before the last base")]),
        source="Atrandi DGPM02323206001 V3 p4 example primer table")
    i7=graph.input("i7_primer","Atrandi V3 example indexed P7 library primer",MaterialKind.OLIGO,
        state=Construct([_seg("P7", P7_SEQ, "p7"),
            _seg("example i7", ATRANDI_I7_INDEX,
                 feature=feature("sample_i7", "sample_index", "fixed")),
            _seg("Read 2 arm", ATRANDI_I7_ANNEAL, "t7",
                 note="V3 p4: terminal phosphorothioate linkage before the last base")]),
        source="Atrandi DGPM02323206001 V3 p4 example primer table")
    pcr_pol=graph.input("pcr_polymerase",
        "NEBNext Ultra II Q5 Master Mix (E7805S component; 98 °C 20 s / 54 °C 30 s / 72 °C 20 s, 8–14 cycles)",
        MaterialKind.CATALYST,source="Atrandi DGPM02323206001 V3 section V")
    final=final_library()
    amplified_lib=graph.product("pcr_product",final.name+" in PCR mix",MaterialKind.LIBRARY,
                                state=final,inferred=True,source=SOURCE+" (barcode-side arm inferred)")
    graph.react("library_pcr","Install P5, P7 and the sample index",StepKind.PCR,
                ((pre,InputRole.TEMPLATE),(p5,InputRole.STRUCTURAL),
                 (i7,InputRole.STRUCTURAL),(pcr_pol,InputRole.CATALYST)),(amplified_lib,),
                geometry=PCRGeometry.unavailable(
                    final,reason="This declaration infers part of the barcode-side landing arm; exact cassette strand evidence remains under source review."))
    library=graph.product("library",final.name,MaterialKind.LIBRARY,state=final,inferred=True,
                          source=SOURCE+" (barcode-side arm inferred)")
    large=graph.product("large_fraction","material bound at 0.6X (beads discarded)",MaterialKind.WASTE)
    small=graph.product("final_supernatant","0.8X bead supernatant and ethanol washes",MaterialKind.WASTE)
    graph.react("double_sided_cleanup","Double-sided 0.6X/0.8X AMPure XP size selection",
                StepKind.CLEANUP,((amplified_lib,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),
                (library,large,small),
                source="Preprint v2 Methods (double-sided AMPure); Atrandi DGPM02323206001 V3 section VI")
    graph.discard(large,"0.6X-bound beads are left behind when the supernatant is transferred")
    graph.discard(small,"bead supernatant and washes are removed")
    graph.sequence(library,SequencingPlan(library=final,primers=tuple(SEQ_PRIMERS),
        acquired_roles=ACQUIRED_ROLES,read_lengths={"Index 1 (i7)":6}))
    graph.annotate(SourceDiscrepancy(
        target="library", title="Barcode A extraction offset in the published kit PDF",
        claims=(
            SourceClaim("Atrandi DGPM02323206001 V3, 3 November 2023, pp. 13–14 barcode extraction pseudocode",
                        "Barcode A is read[36:42], a six-base slice."),
            SourceClaim("The same V3 PDF, p. 13 barcode description/Figure 6 and pp. 16–17 barcode tables",
                        "Every barcode round has eight bases; barcode A occupies positions 37–44.")),
        model_choice="The read-layout model uses the eight-base slice read[36:44].",
        evidence="The stated eight-base barcode length and positions 37–44 agree with that slice.",
        uncertainty="The quoted V3 pseudocode conflicts with its eight-base description; a later corrected revision or intended slice has not been established here.",
        slot=AnnotationSlot.INTRO))
    graph.annotate(SourceDiscrepancy(
        target="library", title="V3 Figure 6 example and barcode tables",
        claims=(
            SourceClaim("Atrandi DGPM02323206001 V3 (3 November 2023), printed p. 13, Figure 6 image",
                        "The base-resolved example shows GAACTGGA at barcode D positions 1–8 and ATGCGCCC at barcode C positions 13–20, with intervening linker tracks."),
            SourceClaim("The same V3 PDF, printed pp. 16–17, Tables 2–3",
                        "Neither example eight-mer nor its reverse complement occurs among the 96 listed barcode sequences.")),
        model_choice="Barcode and linker segments remain placeholders; Figure 6 is not used to assign kit-wide linker sequences or ordered cassette strands.",
        evidence="The two readable example barcode spans were checked directly in the original figure and against all original table entries on both strands.",
        uncertainty="The relationship of this example to the listed kit barcode sets is unresolved. The example contains bases, but does not establish that its linker bases apply to the represented kit or experiment.",
        slot=AnnotationSlot.INTRO))
    graph.annotate(SourceDiscrepancy(
        target="barcoded_a", title="Round-A barcode whitelist disagreement",
        claims=(
            SourceClaim("Atrandi DGPM02323206001 V3 (3 November 2023), Tables 2–3",
                        "The round-A table contains 24 eight-base sequences, including AGACTCAC, but not GTCCGATT."),
            SourceClaim("Authors' Bascet demultiplexer, crates/bascet-cli/src/barcode/atrandi_barcodes.tsv",
                        "The round-A whitelist contains GTCCGATT, but not AGACTCAC; only 12 of 24 sequences overlap the PDF table, and well mappings also differ.")),
        model_choice="The schematic uses an eight-base placeholder and does not select either whitelist.",
        evidence="Direct set comparison of both original tables; rounds B, C and D have matching sequence sets.",
        uncertainty="A kit revision, plate lot or source error may explain the round-A difference; the used plate is unknown.",
        slot=AnnotationSlot.INTRO))
    scope=["microbial_cells","lysis_condition","lyse","pta_primer","phi29","pta_amplify",
           "endprep_enzymes","end_repair","da_tail","dna_ligase"]
    for letter in "abcd": scope.extend(("barcode_"+letter,"ligate_"+letter))
    scope.extend(("release_condition","release","ampure_xp","release_cleanup","fs_mix",
                  "fragment","fs_endprep","ligation_adapter","ultraii_ligation",
                  "adapter_ligation","ligation_cleanup","p5_primer","i7_primer",
                  "pcr_polymerase","library_pcr","double_sided_cleanup","library"))
    graph.declare_reviewed_scope(scope,SOURCE)
    return graph.require_publishable()
