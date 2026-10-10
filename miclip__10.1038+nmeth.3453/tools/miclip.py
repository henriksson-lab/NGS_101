"""Linder et al. m6A-miCLIP library and crosslink signatures."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
import illumina as il
import seqprimers as sp
from batch_ngs import seg
from chemdraw import Construct, MolecularState, Scene, Workflow, circle_rows, feature
from circular import circularize_ssdna

TITLE = "miCLIP — single-nucleotide m6A/m6Am mapping"
NOTES = "01_miclip.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1038/nmeth.3453">Linder et al. (2015)</a>; detailed open protocol: <a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC5562447/">Methods in Molecular Biology</a>.'
SUMMARY = "UV crosslinking fixes an anti-m6A antibody at a modified adenosine. Reverse transcriptase records that protein–RNA adduct as a characteristic substitution or a cDNA truncation before circularization and PCR."
CAVEAT = "The final construct shows the published inline random/barcode layout; antibody-dependent RT outcomes are biological alternatives, not deterministic base conversions."

CLIP_R1 = "GTTCAGAGTTCTACAGTCCGACGATC"
INLINE_LAYOUT = (3, 4, 2)  # random bases, sample barcode, random bases


def inline_identifiers():
    """Published NNN–BBBB–NN inline layout, with identity attached at creation."""
    return [
    seg("UMI bases 1–3", "N" * INLINE_LAYOUT[0], "umi", placeholder=True,
        feature=feature("miclip_umi_5p", "umi", "random", group="miclip_umi", part="bases 1–3")),
    seg("four-base sample barcode", "B" * INLINE_LAYOUT[1], "cbc", placeholder=True,
        feature=feature("inline_sample_barcode", "sample_index", "whitelist", whitelist="published miCLIP barcode set")),
    seg("UMI bases 4–5", "N" * INLINE_LAYOUT[2], "umi", placeholder=True,
        feature=feature("miclip_umi_3p", "umi", "random", group="miclip_umi", part="bases 4–5")),
    ]


FINAL_LIBRARY = Construct([
    seg("P5", il.P5, "p5"), seg("miCLIP Read 1 site", CLIP_R1, "r1"),
    *inline_identifiers(),
    seg("crosslinked RNA-derived cDNA", "R" * 36, placeholder=True),
    seg("3-prime adapter-derived arm", "X" * 28, "r2", placeholder=True),
    seg("P7 reverse complement", il.P7_RC, "p7"),
], name="single-end miCLIP library")
SEQ_PRIMERS = (sp.custom("Read 1", "miCLIP small-RNA Read 1", CLIP_R1,
                         "Linder et al. miCLIP library design"),)
if problems := sp.verify(FINAL_LIBRARY, SEQ_PRIMERS, required_roles=("Read 1",)):
    raise ValueError("miCLIP final library: " + "; ".join(problems))
READ_LENGTHS = {"Read 1": 51}
FINAL_CAPTION = "Read 1 begins with NNN–BBBB–NN: five random bases form a split UMI and the central four bases identify the library; RNA-derived cDNA follows."
SEQUENCING_INTRO = "The inline identifier cycles are computed from the Read 1 landing site. Crosslink-induced substitutions or stops are interpreted within the following RNA-derived cDNA."

_initial = Scene(); _initial.strand("RNA", [seg("fragmented poly(A) RNA", "R" * 42, placeholder=True)], label="fragmented poly(A) RNA")
_initial.mark("RNA", "fragmented poly(A) RNA", "contains m6A or cap-proximal m6Am")
INITIAL_ROWS = tuple(_initial.rows()); INITIAL_NAME = "Fragmented poly(A) RNA"


def workflow():
    cross = Scene(); cross.strand("RNA", [seg("RNA before site", "R" * 18, placeholder=True), seg("m6A crosslink site", "A", "w1"), seg("RNA after site", "R" * 18, placeholder=True)], label="anti-m6A–RNA UV adduct")
    cross.mark("RNA", "m6A crosslink site", "254-nm UV antibody crosslink")
    lig = Scene(); lig.strand("RNA", [seg("crosslinked RNA", "R" * 36, placeholder=True), seg("pre-adenylated 3-prime adapter", "X" * 18, "r2", placeholder=True)], label="immunopurified adapter-ligated RNA")
    lig.junction("RNA", "crosslinked RNA", "pre-adenylated 3-prime adapter", "RNA ligation")
    outcomes = Scene(); outcomes.strand("mutation", [
        seg("miCLIP Read 1 site", CLIP_R1, "r1"), *inline_identifiers(),
        seg("cDNA before site", "D" * 18, placeholder=True),
        seg("crosslink-induced substitution", "X", "w1", placeholder=True),
        seg("continued cDNA", "D" * 18, placeholder=True),
        seg("3-prime adapter-derived arm", "X" * 28, "r2", placeholder=True),
    ], label="RT read-through: mutation signature")
    outcomes.strand("truncation", [
        seg("miCLIP Read 1 site", CLIP_R1, "r1"), *inline_identifiers(),
        seg("truncated RNA-derived cDNA", "D" * 18, placeholder=True),
        seg("3-prime adapter-derived arm", "X" * 28, "r2", placeholder=True),
    ], label="RT stop: cDNA-end signature")
    outcomes.mark("mutation", "crosslink-induced substitution", "misincorporation at/near m6A")
    outcomes.mark("truncation", "truncated RNA-derived cDNA", "3′ end records RT arrest")
    linear_cdna = Construct([
        seg("miCLIP Read 1 site", CLIP_R1, "r1"), *inline_identifiers(),
        seg("RNA-derived cDNA", "D" * 36, placeholder=True),
        seg("3-prime adapter-derived arm", "X" * 28, "r2", placeholder=True),
    ], name="miCLIP first-strand cDNA")
    closed = circularize_ssdna(linear_cdna, five_prime_phosphate=True)
    circle = circle_rows(closed.linear, "CircLigase closure")
    reopened = Scene(); reopened.strand("cDNA", list(linear_cdna), label="restriction-linearized cDNA")
    reopened.mark("cDNA", "miCLIP Read 1 site", "PCR-accessible linear product")
    final = Scene.duplex(list(FINAL_LIBRARY), label="PCR-amplified miCLIP library")
    w = Workflow(MolecularState(INITIAL_NAME, INITIAL_ROWS))
    w.react("Bind anti-m6A and UV-crosslink", cross.rows(), name="Antibody–RNA adduct", note="The antibody is covalently fixed at m6A/m6Am-containing fragments.")
    w.react("Immunopurify and ligate the 3′ adapter", lig.rows(), name="Adapter-ligated RNA", note="Protein A/G enrichment precedes 3′-adapter ligation and size purification.")
    w.react("Reverse-transcribe across the residual peptide adduct", outcomes.rows(), name="Mutation or truncation cDNA", note="Read-through creates substitutions; arrest creates cDNA truncations. Both locate the crosslink.")
    w.react("Circularize first-strand cDNA", circle, name="Covalently closed cDNA", note="CircLigase joins the required 5′ phosphate and 3′ hydroxyl.")
    w.react("Restriction-linearize the cDNA circle", reopened.rows(), name="Re-linearized cDNA", note="A restriction cut opens the circle at a defined adapter site for PCR.")
    w.react("PCR-amplify", final.rows(), name="Sequencing library", note="PCR installs the flow-cell arms while preserving inline UMI and sample-barcode bases.")
    return w
