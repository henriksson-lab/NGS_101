"""Chen et al. single-cell TAPS.

The cell-barcoding transposome oligos (Additional file 5 Table S3) and the custom
sequencing primers (Table S4) belong to the joint scTAPS/scCAPS+ paper.  They are
declared here once; scCAPS+ imports exactly these declarations because the same
published tables define both protocols' front end and readout.
"""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
import illumina as il
import nextera as nx
import seqprimers as sp
from batch_ngs import seg
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "taps__10.1038+s41587-019-0041-2" / "tools"))
from taps import taps_path
from chemdraw import Construct, Scene, feature, revcomp
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from tagmentation import loaded_tn5_adapter

TITLE = "scTAPS — single-cell 5mC + 5hmC sequencing"
NOTES = "01_sctaps.html"
SOURCE = 'Defining preprint: <a href="https://doi.org/10.1101/2024.02.12.579864">Chen et al. (2024)</a>; authoritative version: <a href="https://doi.org/10.1186/s13059-025-03708-1">Genome Biology (2025)</a>.'
TABLES = SOURCE + " Additional file 5 Tables S3–S4"
SUMMARY = "Each cell is fragmented with one of 96 i5/i7-barcoded Tn5 transposome combinations. Ninety-six gap-filled reactions are pooled, then TET oxidation and pyridine–borane reduction make both 5mC and 5hmC read as T."
CAVEAT = ("Transposome oligos and the four custom sequencing primers are drawn from Additional file 5 Tables S3–S4. "
          "The outer P5/P7 ends and the 8-nt Nextera XT sample indices come from the Nextera XT Index Kit (FC-131-1001), whose primer sequences and chosen indices the paper does not print; "
          "they are shown as INFERRED. The 8-cycle index-read lengths are inferred from the 8-nt Tn5 barcodes; the paper states only 120-bp paired-end reads.")
CONVERTS = ("5mC", "5hmC")
TAPS_PATHS = {base: taps_path(base) for base in ("5mC", "5hmC", "C")}

# --- Additional file 5 Table S3: Tn5MEDS transfer strands (all /5Phos/) -----------
# i5 member: S5 + I5_SPACER + 8-nt barcode + C15 connector + ME (63 nt)
# i7 member: S7 + I7_SPACER + 8-nt barcode + D15 connector + ME (69 nt)
I5_SPACER = "TCCACGC"
C15 = "GCGATCGAGGACGGC"            # "Universal_Connector_A_C15"
I7_SPACER = "CTGTCCCTGTCC"
D15 = "CACCGTCTCCGCCTC"            # "Universal_Connector_B_D15"
I5_BARCODES = ("TATAGCCT", "ATAGAGGC", "CCTATCCT", "GGCTCTGA",
               "AGGCGAAG", "TAATCTTA", "CAGGACGT", "GTACTGAC")
I7_BARCODES = ("CGAGTAAT", "TCTCCGGA", "AATGAGCG", "GGAATCTC", "TTCTGAAT", "ACGAATTC",
               "AGCTTCAG", "GCGCATTA", "CATAGCCG", "TTCGCGGA", "GCGCGAGA", "CTATCGCT")
PMENTS = nx.ME_RC                  # Table S3 "pMENTS", /5Phos/CTGTCTCTTATACACATCT
# The source guard admits only files carrying the catalogue (preprint) DOI, so the
# published Additional file 5 workbook is filtered out, and the bioRxiv text splits
# these oligos with printed line numbers.  tools/selftest.py checks every Table S3
# member and Table S4 primer exactly against the published workbook instead.
SOURCE_OLIGO_EXCEPTIONS = {
    m: "Exact match checked by this protocol's selftest against Genome Biology 2025 Additional file 5 Table S3 (PMC12359873 MOESM5); the shared guard excludes that published file by DOI locator and cannot read the line-numbered preprint copy."
    for m in ("i5_transposon", "i7_transposon")}

CELL_I5 = feature("cell_barcode_i5", "cell_barcode", "combinatorial", group="cell_barcode",
                  part="i5", whitelist="Additional file 5 Table S3 (8 i5 members)")
CELL_I7 = feature("cell_barcode_i7", "cell_barcode", "combinatorial", group="cell_barcode",
                  part="i7", whitelist="Additional file 5 Table S3 (12 i7 members)")
XT_I5 = feature("sample_index_i5", "sample_index", "unknown")
XT_I7 = feature("sample_index_i7", "sample_index", "unknown")


def i5_meds(barcode: str | None = None) -> list:
    """Table S3 P5_i5_n transfer strand; ``barcode`` must be a published member."""
    if barcode is not None and barcode not in I5_BARCODES:
        raise ValueError("i5 barcode is not a Table S3 member")
    bc = (seg("Tn5 i5 cell barcode", barcode, "cbc", feature=CELL_I5) if barcode else
          seg("Tn5 i5 cell barcode", "J" * 8, "cbc", placeholder=True, feature=CELL_I5))
    return [seg("S5", nx.S5, "s5"), seg("i5 spacer", I5_SPACER), bc,
            seg("C15 connector", C15), seg("i5 mosaic end", nx.ME, "me")]


def i7_meds(barcode: str | None = None) -> list:
    """Table S3 P7_i7_n transfer strand; ``barcode`` must be a published member."""
    if barcode is not None and barcode not in I7_BARCODES:
        raise ValueError("i7 barcode is not a Table S3 member")
    bc = (seg("Tn5 i7 cell barcode", barcode, "cbc", feature=CELL_I7) if barcode else
          seg("Tn5 i7 cell barcode", "I" * 8, "cbc", placeholder=True, feature=CELL_I7))
    return [seg("S7", nx.S7, "s7"), seg("i7 spacer", I7_SPACER), bc,
            seg("D15 connector", D15), seg("i7 mosaic end", nx.ME, "me")]


def table_s3_oligos() -> dict[str, str]:
    """The 20 transfer strands plus pMENTS, composed from the declared parts."""
    out = {f"P5_i5_{i}": "".join(s.top for s in i5_meds(b)) for i, b in enumerate(I5_BARCODES, 1)}
    out.update({f"P7_i7_{i}": "".join(s.top for s in i7_meds(b)) for i, b in enumerate(I7_BARCODES, 1)})
    out["pMENTS"] = PMENTS
    return out


def loaded(parts: list, label: str) -> Scene:
    sc = loaded_tn5_adapter(parts, mosaic_end=parts[-1].name, label=label,
                            bottom_name="pMENTS", bottom_label="/5Phos/ pMENTS (Tn5MEDS-Rev)",
                            bottom_sequence=PMENTS, mod5="p")
    sc.mark("Tn5 top", next(s.name for s in parts if s.feature), "well-specific Tn5 barcode")
    return sc


def _rc(parts: list) -> list:
    return [seg(s.name + " reverse complement", revcomp(s.top) if not s.placeholder else s.top,
                s.tag, placeholder=s.placeholder, feature=s.feature) for s in reversed(parts)]


def tagged_fragment(insert: str) -> Construct:
    """Both Tn5 adapters on one fragment; the gap-filled duplex is drawn."""
    return Construct([*i5_meds(), seg(insert, "X" * 40, placeholder=True), *_rc(i7_meds())],
                     name="cell-barcoded Tn5 fragment")


# --- Additional file 5 Table S4: custom NovaSeq sequencing primers -----------------
READ1_PRIMER = C15 + nx.ME                     # GCGATCGAGGACGGCAGATGTGTATAAGAGACAG
READ2_PRIMER = D15 + nx.ME                     # CACCGTCTCCGCCTCAGATGTGTATAAGAGACAG
INDEX1_PRIMER = nx.ME_RC + revcomp(D15)        # CTGTCTCTTATACACATCTGAGGCGGAGACGGTG
INDEX2_PRIMER = nx.ME_RC + revcomp(C15)        # CTGTCTCTTATACACATCTGCCGTCCTCGATCGC
TABLE_S4 = {"Read 1": READ1_PRIMER, "Read 2": READ2_PRIMER,
            "Index 1": INDEX1_PRIMER, "Index 2": INDEX2_PRIMER}
SEQ_PRIMERS = (
    sp.custom("Read 1", "custom Read 1 (Table S4)", READ1_PRIMER, TABLES),
    sp.custom("Index 1 (i7)", "custom Index 1 (Table S4)", INDEX1_PRIMER, TABLES),
    sp.custom("Index 2 (i5)", "custom Index 2 (Table S4)", INDEX2_PRIMER, TABLES),
    sp.custom("Read 2", "custom Read 2 (Table S4)", READ2_PRIMER, TABLES),
)
READ_LENGTHS = {"Read 1": 120, "Index 1 (i7)": 8, "Index 2 (i5)": 8, "Read 2": 120}
READ_USES = {"Read 1": "include", "Read 2": "include",
             "Index 1 (i7)": "extract_and_demultiplex", "Index 2 (i5)": "extract_and_demultiplex"}
SHOW_READ_CONTENT = True


def final_library(insert: str, name: str) -> Construct:
    lib = Construct([
        seg("P5", il.P5, "p5", inferred=True),
        seg("Nextera XT i5 sample index", "J" * 8, "cbc", placeholder=True, inferred=True, feature=XT_I5),
        *i5_meds(), seg(insert, "X" * 400, placeholder=True, length_bp=400), *_rc(i7_meds()),
        seg("Nextera XT i7 sample index reverse complement", "I" * 8, "cbc", placeholder=True,
            inferred=True, feature=XT_I7),
        seg("P7 reverse complement", il.P7_RC, "p7", inferred=True),
    ], name=name)
    if problems := sp.verify(lib, SEQ_PRIMERS):
        raise ValueError(name + ": " + "; ".join(problems))
    return lib


FINAL_LIBRARY = final_library("TAPS-converted genomic insert", "scTAPS library")
FINAL_CAPTION = ("Table S3 transposon arms flank the insert; Nextera XT PCR adds P5/P7 and sample indices (INFERRED, kit sequences not printed). "
                 "In the reads, 5mC and 5hmC appear as T.")
SEQUENCING_INTRO = ("The four Table S4 custom primers land on the C15/D15 connector plus mosaic end. Read 1 and Read 2 (120 bp each, published) start in the insert; "
                    "the custom index primers read the Tn5 i7 and i5 cell barcodes, not the Nextera XT sample indices further out. "
                    "Index-read cycles are not published; 8 is inferred from the barcode length.")

_initial = Scene.duplex([seg("single-cell genomic DNA", "D" * 52, placeholder=True)], label="genomic DNA from one FACS-sorted cell")
INITIAL_ROWS = tuple(_initial.rows()); INITIAL_NAME = "Single-cell genomic DNA"


def tagged_scene(label: str) -> Scene:
    sc = Scene.duplex(list(tagged_fragment("genomic fragment")), label=label)
    sc.labels("top")
    return sc


def base_scene(paths: dict, title: str, stage: int, label: str) -> Scene:
    states = {base: path.after(stage) for base, path in paths.items()}
    sc = tagged_scene(label)
    for base, letter in (("5mC", "M"), ("5hmC", "H"), ("C", "C")):
        if base in states:
            sc.strand(base, [seg(f"{base} path: {states[base]}", letter, "w1" if base != "C" else None,
                                 placeholder=True)], label=f"{base}: {states[base]}")
    sc.footer(title, "C"); return sc


def front_end(g: ProtocolGraph, controls_label: str):
    """Shared published front end: sort, lyse, tagment, strip, gap-fill, pool."""
    cell = g.input("single_cell_dna", "genomic DNA of one FACS-sorted cell or nucleus",
                   MaterialKind.BIOLOGICAL, state=INITIAL_ROWS, primary=True)
    protease = g.input("protease_lysis", "4 µl TAPS-NaOH/MgCl2/NEB buffer 4/Triton X-100 lysis buffer with 1 µg Qiagen protease; 50 °C 3 h, 75 °C 30 min, 80 °C 15 min",
                       MaterialKind.CATALYST)
    lysate = g.product("lysate", "protease-released single-cell DNA", MaterialKind.MOLECULAR, state=INITIAL_ROWS)
    g.react("lyse", "Lyse and heat-inactivate the protease", StepKind.DEPROTEINIZATION,
            ((cell, InputRole.TEMPLATE), (protease, InputRole.CATALYST)), (lysate,), changes_molecule=False)
    i5 = g.input("i5_transposon", "Table S3 i5 Tn5MEDS annealed to pMENTS (i5-1 drawn; 8 members)",
                 MaterialKind.ADAPTER, state=loaded(i5_meds(I5_BARCODES[0]), "/5Phos/ P5_i5_1 transfer strand"), source=TABLES)
    i7 = g.input("i7_transposon", "Table S3 i7 Tn5MEDS annealed to pMENTS (i7-1 drawn; 12 members)",
                 MaterialKind.ADAPTER, state=loaded(i7_meds(I7_BARCODES[0]), "/5Phos/ P7_i7_1 transfer strand"), source=TABLES)
    tn5 = g.input("tn5", "in-house pTBX1-Tn5, 1 µl 1× preassembled per well; 1× TAPS buffer, 6% PEG 3350, 50 °C 15 min",
                  MaterialKind.CATALYST)
    tagged = g.product("tagged", "Tn5-bound cell-barcoded fragments with nine-base gaps (gaps not drawn)",
                       MaterialKind.MOLECULAR, state=tagged_scene("i5/i7 transposon arms joined to genomic 5′ ends"))
    same_arm = g.product("same_arm_fragments", "i5/i5 or i7/i7 fragments and unproductive insertions", MaterialKind.WASTE)
    g.react("tagment", "Tagment with one well-specific i5/i7 transposome", StepKind.TAGMENTATION,
            ((lysate, InputRole.TEMPLATE), (i5, InputRole.STRUCTURAL), (i7, InputRole.STRUCTURAL),
             (tn5, InputRole.CATALYST)), (tagged, same_arm))
    g.discard(same_arm, "carried into PCR; lacks one i5 and one i7 arm, so it is not exponentially amplified by the paired XT primers; no physical separation is modelled")
    sds = g.input("sds_strip", "1.11 µl 1.0% SDS 5 min RT, then 1.24 µl 5.0% NP-40 or Triton X-100 5 min RT", MaterialKind.CONDITION)
    stripped = g.product("stripped", "protein-free gapped cell-barcoded fragments", MaterialKind.MOLECULAR,
                         state=tagged_scene("Tn5 released; nine-base gaps remain (not drawn)"))
    g.react("strip", "Strip Tn5 with SDS and neutralize", StepKind.DEPROTEINIZATION,
            ((tagged, InputRole.TEMPLATE), (sds, InputRole.CONDITION)), (stripped,), changes_molecule=False)
    phusion = g.input("gap_fill_polymerase", "12.35 µl 2× Phusion High-Fidelity PCR Master Mix (F531L), 72 °C 5 min", MaterialKind.CATALYST)
    filled = g.product("gap_filled", "gap-filled cell-barcoded duplex", MaterialKind.MOLECULAR,
                       state=tagged_scene("barcoded, gap-filled single-cell fragment"))
    g.react("gap_fill", "Fill the nine-base transposition gaps", StepKind.EXTENSION,
            ((stripped, InputRole.TEMPLATE), (phusion, InputRole.CATALYST)), (filled,))
    peers = g.input("peer_wells", "gap-filled fragments from the other 95 indexed wells", MaterialKind.MOLECULAR,
                    state=tagged_scene("other well-indexed gap-filled fragments"))
    controls = g.input("conversion_controls", controls_label, MaterialKind.MOLECULAR)
    pooled = g.product("pooled", "pooled 96-cell DNA plus filler and spike-ins", MaterialKind.MOLECULAR,
                       state=tagged_scene("pooled DNA from 96 uniquely indexed wells plus controls"))
    g.react("pool_wells", "Pool 96 wells in a 50-ml tube with filler and spike-ins", StepKind.POOL,
            ((filled, InputRole.TEMPLATE), (peers, InputRole.TEMPLATE), (controls, InputRole.TEMPLATE)),
            (pooled,), changes_molecule=False)
    beads = g.input("pool_purification", "DNA purification of the pool (method not specified)", MaterialKind.SELECTION)
    purified = g.product("purified_pool", "purified pooled DNA", MaterialKind.MOLECULAR,
                         state=tagged_scene("purified pooled cell-barcoded DNA"))
    pool_waste = g.product("pool_waste", "SDS, detergent, enzyme and buffer components", MaterialKind.WASTE)
    g.react("purify_pool", "Purify the pooled DNA", StepKind.CLEANUP,
            ((pooled, InputRole.TEMPLATE), (beads, InputRole.SELECTION)), (purified, pool_waste))
    g.discard(pool_waste, "removed by purification")
    return purified


def library_back_end(g: ProtocolGraph, template, library: Construct, reviewed: list):
    """Shared published back end: XT-indexed Uracil+ PCR, bead selection, custom run."""
    i5p = g.input("xt_i5_primer", "Nextera XT Index Kit (FC-131-1001) i5 primer", MaterialKind.OLIGO,
                  state=[seg("P5", il.P5, "p5", inferred=True),
                         seg("XT i5", "J" * 8, "cbc", placeholder=True, inferred=True), seg("S5", nx.S5, "s5", inferred=True)])
    i7p = g.input("xt_i7_primer", "Nextera XT Index Kit (FC-131-1001) i7 primer", MaterialKind.OLIGO,
                  state=[seg("P7", il.P7, "p7", inferred=True),
                         seg("XT i7", "I" * 8, "cbc", placeholder=True, inferred=True), seg("S7", nx.S7, "s7", inferred=True)])
    pol = g.input("pcr_polymerase", "KAPA HiFi HotStart Uracil+ ReadyMix; 98 °C 45 s, 13 × (98 °C 10 s, 60 °C 15 s, 72 °C 1 min), 72 °C 3 min",
                  MaterialKind.CATALYST)
    amplified = g.product("amplified", library.name + " (unselected PCR product)", MaterialKind.MOLECULAR, state=library)
    g.react("library_pcr", "INFERRED — Nextera XT P5/P7 ends and sample indices: amplify 13 cycles with Uracil+ polymerase; DHU is copied as T", StepKind.PCR,
            ((template, InputRole.TEMPLATE), (i5p, InputRole.STRUCTURAL), (i7p, InputRole.STRUCTURAL),
             (pol, InputRole.CATALYST)), (amplified,),
            geometry=PCRGeometry.unavailable(library, reason="The paper names the Nextera XT Index Kit but prints neither its primer sequences nor the indices used; the kit original was not read in this pass."))
    beads = g.input("size_selection", "0.5×–0.25× AMPure XP double-sided selection", MaterialKind.SELECTION)
    final = g.product("library", library.name, MaterialKind.LIBRARY, state=library)
    rejected = g.product("rejected_sizes", "PCR components and fragments outside 500–700 bp", MaterialKind.WASTE)
    g.react("size_select", "Purify and select 500–700 bp", StepKind.SELECTION,
            ((amplified, InputRole.TEMPLATE), (beads, InputRole.SELECTION)), (final, rejected))
    g.discard(rejected, "removed by bead purification and double-sided size selection")
    g.sequence(final, SequencingPlan(library=library, primers=SEQ_PRIMERS,
                                     acquired_roles=tuple(p.role for p in SEQ_PRIMERS),
                                     read_lengths=READ_LENGTHS, read_uses=READ_USES))
    g.declare_reviewed_scope([
        "single_cell_dna", "protease_lysis", "lyse", "i5_transposon", "i7_transposon", "tn5", "tagment",
        "sds_strip", "strip", "gap_fill_polymerase", "gap_fill", "peer_wells", "conversion_controls",
        "pool_wells", "pool_purification", "purify_pool", *reviewed,
        "xt_i5_primer", "xt_i7_primer", "pcr_polymerase", "library_pcr", "amplified",
        "size_selection", "size_select", "library"], SOURCE)
    return g.require_publishable()


def protocol_graph():
    g = ProtocolGraph(TITLE, authority=SOURCE)
    pooled = front_end(g, "200 ng 2-kb filler DNA plus 0.1% spike-ins (methylated lambda, 2-kb unmodified, 5hmC spike-in)")
    mtet1 = g.input("mtet1", "4 µM mTet1CD, Fe(II)/α-KG/ascorbate/DTT/ATP in HEPES pH 8.0, 37 °C 80 min; proteinase K 50 °C 1 h; 1.8× AMPure; repeated once",
                    MaterialKind.CATALYST)
    oxidized = g.product("oxidized_dna", "5caC-containing indexed DNA", MaterialKind.MOLECULAR,
                         state=base_scene(TAPS_PATHS, "two mTET1 oxidation rounds", 1, "after two TET rounds"))
    g.react("oxidize", "Oxidize 5mC and 5hmC to 5caC (two mTET1 rounds)", StepKind.CHEMICAL_CONVERSION,
            ((pooled, InputRole.TEMPLATE), (mtet1, InputRole.CATALYST)), (oxidized,))
    borane = g.input("pyridine_borane", "pyridine–borane reduction; Zymo column with Oligo Binding Buffer", MaterialKind.CONDITION)
    reduced = g.product("reduced_dna", "DHU-containing indexed DNA", MaterialKind.MOLECULAR,
                        state=base_scene(TAPS_PATHS, "pyridine–borane reduction", 2, "after borane reduction"))
    g.react("reduce", "Reduce 5caC to dihydrouracil", StepKind.CHEMICAL_CONVERSION,
            ((oxidized, InputRole.TEMPLATE), (borane, InputRole.CONDITION)), (reduced,))
    return library_back_end(g, reduced, FINAL_LIBRARY,
                            ["mtet1", "oxidize", "oxidized_dna", "pyridine_borane", "reduce"])
