"""Spatial-CUT&Tag / hsrChST-seq molecular model.

Defining publication: Deng et al., Science 2022 (doi:10.1126/science.abg7216). Oligos,
barcode pools and reaction conditions come from its Supplementary Materials (author
manuscript EMS145971 DOCX, Materials and Methods and Tables S1-S3, read 2026-10-17).
The authors' later detailed protocol, Farzad et al., Nature Protocols 2024 (MOESM2
Supp_Table 1-3), prints identical oligos and barcodes except Tn5ME-A.

The construct is built locally rather than through ``spatial_epigenome``: that shared
recipe draws linker 1's reverse complement directly against genomic DNA and so omits
the strand that pA-Tn5 transfers.  The transferred Tn5ME-A strand's 5' extension is the
part ligation linker 1 splints to barcode A; the 2022 Tn5ME-A provides exactly that
extension (checked below), whereas the 2024 Supp_Table 1 Tn5ME-A does not.
"""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import illumina as il
import nextera as nx
import seqprimers as sp
from adapter_ligation import reported_adapter_end
from batch_ngs import seg
from chemdraw import Construct, Scene, feature, revcomp
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (Disposition, InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, SourceClaim, SourceDiscrepancy,
                            StepKind, Terminal)


SCIENCE = "https://doi.org/10.1126/science.abg7216"
SCIENCE_SUPP = ("Deng et al. 2022 Science Supplementary Materials (EMS145971), Materials and "
                "Methods \u201cTransposome preparation\u201d and \u201cSpatial-CUT&Tag\u201d; Tables S1\u2013S3")
FARZAD = "https://doi.org/10.1038/s41596-024-01013-y"
TITLE = "Spatial-CUT&Tag (hsrChST-seq)"
NOTES = "01_spatial-cut-tag.html"
SOURCE = (f'Defining source: <a href="{SCIENCE}">Deng et al., <i>Science</i> (2022)</a>; '
          'first named hsrChST-seq in the '
          '<a href="https://doi.org/10.1101/2021.03.11.434985">2021 preprint</a>. '
          'Oligos and procedure: its Supplementary Materials. Later detailed protocol: '
          f'<a href="{FARZAD}">Farzad et al., <i>Nature Protocols</i> (2024)</a>.')
SUMMARY = ("Antibody-tethered pA–Tn5 transfers linker-bearing adapters beside a chosen "
           "histone mark. Two perpendicular in-tissue ligations append 8-nt A and B "
           "barcodes whose ordered pair identifies a tissue pixel.")
CAVEAT = ("Oligos, barcode pools and conditions follow the 2022 Science Supplementary "
          "Materials. The authors' 2024 detailed protocol prints the same oligos and "
          "barcodes except Tn5ME-A, whose 2024 sequence does not pair with ligation "
          "linker 1 (see the discrepancy below). Which Tn5ME-B/Tn5ME-A pairing each "
          "transposome carries is not stated; the drawn fragment is the A–B end class "
          "that both ligations and the PCR require.")

# -- Deng et al. 2022 Supplementary Materials, "Transposome preparation" (DOCX text
# lines 911-918) and Table S1 (lines 1377-1391); Farzad et al. 2024 MOESM2 Supp_Table 1
# prints the same bases for every oligo here except Tn5ME-A.
TN5_ME_REV_WRITTEN = "/5Phos/CTGTCTCTTATACACATCT"
TN5_ME_A_WRITTEN = "/5Phos/CATCGGCGTACGACTAGATGTGTATAAGAGACAG"          # 2022 only
TN5_ME_A_TABLE1_WRITTEN = "/5Phos/TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"    # 2024 Supp_Table 1
TN5_ME_B = "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"
LIGATION_LINKER_1 = "AGTCGTACGCCGATGCTGTCTCTTATACACATCT"
LIGATION_LINKER_2 = "CGAATGCTCTGGCCTCTCAAGCACGTGGAT"
CUSTOM_READ_1 = "GCCTGTCCGCGGAAGCAGTGGTATCAACGCAGAGTAC"
# -- 2022 Table S1 (2024 Supp_Table 2).
NEW_P5_SMART_WRITTEN = ("AATGATACGGCGACCACCGAGATCTACACGCCTGTCCGCGGAAGCAG"
                        "TGGTATCAACGCAGAGT*A*C")
AD2_1 = "CAAGCAGAAGACGGCATACGAGATTCGCCTTAGTCTCGTGGGCTCGGAGATGTG"
# -- 2022 Table S2 (2024 Supp_Table 3), first member of each plate.
BARCODE_A_OLIGO_1 = "/5Phos/AGGCCAGAGCATTCGAACGTGATAGATGTGTATAAGAGACAG"  # Tn5ME-A_BC1_1
BARCODE_B_OLIGO_1 = ("GCCTGTCCGCGGAAGCAGTGGTATCAACGCAGAGTACAACGTGAT"
                     "ATCCACGTGCTTGAG")                                   # BC2_R1_1


def _bases(written: str) -> str:
    return written.removeprefix("/5Phos/").replace("*", "")


NEW_P5_SMART = _bases(NEW_P5_SMART_WRITTEN)
if NEW_P5_SMART != il.P5 + CUSTOM_READ_1:
    raise ValueError("New-P5-SMART is not P5 followed by the custom Read 1 site")
if revcomp(_bases(TN5_ME_REV_WRITTEN)) != nx.ME:
    raise ValueError("Tn5MErev is not the mosaic-end complement")
BARCODE_A_BODY = _bases(BARCODE_A_OLIGO_1)
LINKER_2_A_SIDE = BARCODE_A_BODY[:15]
if BARCODE_A_BODY[23:] != nx.ME or not BARCODE_B_OLIGO_1.startswith(CUSTOM_READ_1):
    raise ValueError("barcode oligos do not have the expected mosaic-end / Read 1 ends")
LINKER_2_B_SIDE = BARCODE_B_OLIGO_1[len(CUSTOM_READ_1) + 8:]
if revcomp(LIGATION_LINKER_2) != LINKER_2_B_SIDE + LINKER_2_A_SIDE:
    raise ValueError("ligation linker 2 does not splint barcode B to barcode A")
# Linker 1 = revcomp(barcode-A mosaic end + X): X is what the transferred adapter must
# present at its 5' end for the splinted ligation.  It is derived, not printed.
_L1_READ = revcomp(LIGATION_LINKER_1)
if not _L1_READ.startswith(nx.ME):
    raise ValueError("ligation linker 1 does not pair the barcode-A mosaic end")
LINKER_1_REQUIRED_EXTENSION = _L1_READ[len(nx.ME):]
TABLE1_EXTENSION = _bases(TN5_ME_A_TABLE1_WRITTEN)[:-len(nx.ME)]
TABLE1_TN5_ME_A_PAIRS_LINKER_1 = TABLE1_EXTENSION == LINKER_1_REQUIRED_EXTENSION
TN5_ME_A = _bases(TN5_ME_A_WRITTEN)
if not TN5_ME_A.endswith(nx.ME):
    raise ValueError("2022 Tn5ME-A does not end in the mosaic end")
TN5_ME_A_EXTENSION = TN5_ME_A[:-len(nx.ME)]
if TN5_ME_A_EXTENSION != LINKER_1_REQUIRED_EXTENSION:
    raise ValueError("2022 Tn5ME-A extension is not splinted by ligation linker 1")

SPATIAL_A = feature("spatial_a", "spatial_barcode", "whitelist", group="spatial_coordinate",
                    part="A", whitelist="Deng et al. 2022 Table S2, 1st barcodes (A1–A50)")
SPATIAL_B = feature("spatial_b", "spatial_barcode", "whitelist", group="spatial_coordinate",
                    part="B", whitelist="Deng et al. 2022 Table S2, 2nd barcodes (B1–B50)")
SAMPLE_I7 = feature("sample_index_i7", "sample_index", "whitelist",
                    whitelist="Deng et al. 2022 Table S1, Ad2.1–Ad2.8")
INSERT = "histone-mark-associated genomic fragment"


def _gdna():
    return seg(INSERT, "X" * 46, placeholder=True)


def _adapter():
    return [seg("Tn5ME-A 5′ extension", TN5_ME_A_EXTENSION, "me"),
            seg("Tn5ME-A mosaic end", nx.ME, "me")]


def _right_arm():
    return [seg("Tn5ME-B mosaic end (rc)", nx.ME_RC, "me"),
            seg("S7 (rc)", nx.S7_RC, "s7")]


def tagged_fragment() -> Construct:
    return Construct([*_adapter(), _gdna(), *_right_arm()],
                     name="pA–Tn5 tagmented fragment")


def a_barcoded_fragment() -> Construct:
    return Construct([seg("linker-2 A-side", LINKER_2_A_SIDE),
                      seg("spatial barcode A", "A" * 8, "cbc", placeholder=True,
                          feature=SPATIAL_A),
                      seg("barcode-A mosaic end", nx.ME, "me"),
                      *_adapter(), _gdna(), *_right_arm()],
                     name="A-barcoded tagmented fragment")


def ab_barcoded_fragment() -> Construct:
    return Construct([seg("custom Read 1 site", CUSTOM_READ_1, "r1"),
                      seg("spatial barcode B", "B" * 8, "cbc", placeholder=True,
                          feature=SPATIAL_B),
                      seg("linker-2 B-side", LINKER_2_B_SIDE),
                      *a_barcoded_fragment()], name="A/B-barcoded tagmented fragment")


def library():
    lib = Construct([seg("P5", il.P5, "p5"), *ab_barcoded_fragment(),
                     seg("i7 (rc)", "I" * 8, "cbc", placeholder=True, feature=SAMPLE_I7),
                     seg("P7 (rc)", il.P7_RC, "p7")],
                    name=INSERT + " spatial-CUT&Tag library")
    primers = (sp.custom("Read 1", "Custom Read 1 primer", CUSTOM_READ_1,
                         "Deng et al. 2022 Table S1"),
               sp.NEXTERA["I1"], sp.NEXTERA["R2"])
    problems = sp.verify(lib, primers, required_roles=tuple(p.role for p in primers))
    if problems:
        raise ValueError("invalid Spatial-CUT&Tag library: " + "; ".join(problems))
    return lib, primers


FINAL_LIBRARY, SEQ_PRIMERS = library()
READ_LENGTHS = {"Read 1": 150, "Index 1 (i7)": 8, "Read 2": 150}
SHOW_READ_CONTENT = True
FINAL_CAPTION = ("A–B end class: barcode B, linker-2 halves and barcode A are joined to "
                 "the transferred Tn5ME-A strand; the opposite end carries Tn5ME-B. The "
                 "custom-primer end reads barcode B, linker 2 and barcode A first. ** "
                 "marks the two in-tissue ligations.")
_GENOMIC_START = 1 + sum(len(x) for x in (
    "B" * 8, LINKER_2_B_SIDE, LINKER_2_A_SIDE, "A" * 8, nx.ME, TN5_ME_A_EXTENSION, nx.ME))
SEQUENCING_INTRO = ("HiSeq 4000 paired-end 150 sequencing with a custom Read 1 primer "
                    "(2022 Supplementary Materials; GEO GSE165217). The custom Read 1 "
                    "primer is printed in 2022 Table S1. Raw Read 1 carries barcode B at "
                    "cycles 1–8 and barcode A at cycles 39–46; the authors filter it on "
                    "the two constant linker sequences, and genomic DNA begins at cycle "
                    f"{_GENOMIC_START}. Read 2 reports the opposite genomic end.")


def oligos():
    return (f"<b>Tn5MErev</b> — 5′-{TN5_ME_REV_WRITTEN}-3′",
            f"<b>Tn5ME-A</b> (2022) — 5′-{TN5_ME_A_WRITTEN}-3′",
            f"<b>Tn5ME-A (2024 Supp_Table 1; does not pair with linker 1)</b> — "
            f"5′-{TN5_ME_A_TABLE1_WRITTEN}-3′",
            f"<b>Tn5ME-B</b> — 5′-{TN5_ME_B}-3′",
            f"<b>Ligation linker 1</b> — 5′-{LIGATION_LINKER_1}-3′",
            f"<b>Ligation linker 2</b> — 5′-{LIGATION_LINKER_2}-3′",
            f"<b>Tn5ME-A_BC1_1</b> (barcode A, 1 of 50) — 5′-{BARCODE_A_OLIGO_1}-3′",
            f"<b>BC2_R1_1</b> (barcode B, 1 of 50) — 5′-{BARCODE_B_OLIGO_1}-3′",
            f"<b>New-P5-SMART PCR hybrid oligo</b> — 5′-{NEW_P5_SMART_WRITTEN}-3′",
            f"<b>Ad2.1_TAAGGCGA</b> (1 of 8) — 5′-{AD2_1}-3′",
            f"<b>Custom Read 1 primer</b> — 5′-{CUSTOM_READ_1}-3′")


def final_rows():
    sc = Scene.duplex(list(FINAL_LIBRARY), label="library")
    sc.junction("top", "linker-2 B-side", "linker-2 A-side", "ligation B")
    sc.junction("top", "barcode-A mosaic end",
                "Tn5ME-A 5′ extension", "ligation A")
    sc.labels("top")
    return sc.rows()


def protocol_graph() -> ProtocolGraph:
    auth = SCIENCE + " (" + SCIENCE_SUPP + ")"
    g = ProtocolGraph("Spatial-CUT&Tag", authority=auth)
    tissue = g.input("tissue", "fixed tissue section", MaterialKind.BIOLOGICAL, primary=True)
    npd = g.input("np40_digitonin", "NP40–digitonin wash buffer, 5 min",
                  MaterialKind.CONDITION)
    perm = g.product("permeabilized", "permeabilized tissue chromatin",
                     MaterialKind.BIOLOGICAL)
    g.react("permeabilize", "Wash and permeabilize the section", StepKind.PERMEABILIZATION,
            ((tissue, InputRole.TEMPLATE), (npd, InputRole.CONDITION)), (perm,))
    primary = g.input("primary_antibody", "histone-mark primary antibody (H3K27me3, "
                      "H3K4me3 or H3K27ac), 1:50", MaterialKind.MOLECULAR)
    bound1 = g.product("primary_bound", "primary-antibody-marked chromatin",
                       MaterialKind.BIOLOGICAL)
    unbound1 = g.product("unbound_primary", "removed primary antibody solution",
                         MaterialKind.WASTE)
    g.react("primary_binding", "Bind primary antibody overnight at 4 °C",
            StepKind.BINDING, ((perm, InputRole.TEMPLATE), (primary, InputRole.ANALYTE)),
            (bound1, unbound1))
    g.discard(unbound1, "primary antibody solution is removed and discarded")
    secondary = g.input("secondary_antibody", "secondary antibody, 1:50",
                        MaterialKind.MOLECULAR)
    bound2 = g.product("secondary_bound", "primary/secondary-antibody-marked chromatin",
                       MaterialKind.BIOLOGICAL)
    unbound2 = g.product("unbound_secondary", "secondary antibody and wash-buffer fraction",
                         MaterialKind.WASTE)
    g.react("secondary_binding", "Bind secondary antibody 30 min RT, 5 min wash-buffer wash",
            StepKind.BINDING, ((bound1, InputRole.TEMPLATE), (secondary, InputRole.ANALYTE)),
            (bound2, unbound2))
    g.discard(unbound2, "secondary antibody solution and washes are removed")
    adapter = g.input("pa_tn5_adapter", "Tn5ME-A and Tn5ME-B, each annealed to Tn5MErev, "
                      "loaded per Diagenode guidelines (transferred Tn5ME-A strand shown)",
                      MaterialKind.ADAPTER, state=tuple(_adapter()))
    patn5 = g.input("pa_tn5", "unloaded pA–Tn5 (Diagenode C01070002); complex 1:100 in "
                    "300-wash buffer, 1 h RT", MaterialKind.CATALYST)
    mg = g.input("tagmentation_buffer", "10 mM MgCl2 in 300-wash buffer, 37 °C 1 h; 40 mM "
                 "EDTA 5 min RT; NEBuffer 3.1 wash", MaterialKind.CONDITION)
    tagged = g.product("tagged_fragment", "pA–Tn5 tagmented fragment", MaterialKind.MOLECULAR,
                       state=tagged_fragment())
    unbound3 = g.product("unbound_patn5", "excess pA–Tn5 (5 min 300-wash)",
                         MaterialKind.WASTE)
    same = g.product("same_adapter_fragments", "candidate same-adapter fragments",
                     MaterialKind.WASTE, inferred=True,
                     source=auth + " (end classes derived from a two-adapter transposome; "
                     "not described in the source)")
    g.react("tagment", "Bind pA–Tn5 1 h RT, wash with 300-wash, activate with Mg2+ "
            "37 °C 1 h, stop with EDTA",
            StepKind.TAGMENTATION, ((bound2, InputRole.TEMPLATE),
            (adapter, InputRole.STRUCTURAL), (patn5, InputRole.CATALYST),
            (mg, InputRole.CONDITION)), (tagged, unbound3, same))
    g.discard(unbound3, "unbound pA–Tn5 is washed away before activation")
    g.mark_terminal(Terminal(same.id, Disposition.UNRESOLVED,
        "Candidate same-adapter fragments remain in the tissue through ligation and lysis; "
        "no physical separation is described. Failing to amplify exponentially is not "
        "a physical removal."))
    linker1 = g.input("barcode_a", "barcode A1–A50 annealed to ligation linker 1 "
                      "(2X annealing buffer); 1 µL of 25 µM per channel",
                      MaterialKind.ADAPTER, state=(BARCODE_A_BODY, LIGATION_LINKER_1))
    ligase = g.input("ligase", "T4 DNA ligase (NEB M0202L) in T4 ligase buffer with "
                     "Triton X-100 and NEBuffer 3.1", MaterialKind.CATALYST)
    a_prod = g.product("a_barcoded", "A-barcoded tagmented fragment", MaterialKind.MOLECULAR,
                       state=a_barcoded_fragment())
    g.react("barcode_a_ligation", "Flow barcode A through the 50 first channels, ligate "
            "37 °C 30 min, wash NEBuffer 3.1",
            StepKind.LIGATION, ((tagged, InputRole.TEMPLATE), (linker1, InputRole.STRUCTURAL),
            (ligase, InputRole.CATALYST)), (a_prod,),
            geometry=reported_adapter_end(
                "barcode-A oligo", kind="splinted", polarity="3-prime",
                terminal_sequence=nx.ME,
                statement="Linker 1's 3' half pairs the barcode-A 3' mosaic end; its 5' "
                          "half pairs the 5'-phosphorylated extension of the transferred "
                          "Tn5ME-A strand (2022 sequence)."))
    linker2 = g.input("barcode_b", "barcode B1–B50 annealed to ligation linker 2; "
                      "1 µL of 25 µM per channel",
                      MaterialKind.ADAPTER, state=(BARCODE_B_OLIGO_1, LIGATION_LINKER_2))
    ab_prod = g.product("ab_barcoded", "A/B-barcoded tagmented fragment",
                        MaterialKind.MOLECULAR, state=ab_barcoded_fragment())
    g.react("barcode_b_ligation", "Attach the perpendicular device, flow barcode B, ligate "
            "37 °C 30 min, wash DPBS",
            StepKind.LIGATION, ((a_prod, InputRole.TEMPLATE), (linker2, InputRole.STRUCTURAL),
            (ligase, InputRole.CATALYST)), (ab_prod,),
            geometry=reported_adapter_end(
                "barcode-B oligo", kind="splinted", polarity="3-prime",
                terminal_sequence=LINKER_2_B_SIDE,
                statement="Ligation linker 2 is the reverse complement of the barcode-B "
                          "3' end followed by the 5'-phosphorylated barcode-A 5' end."))
    pk = g.input("lysis", "TAPS washes; 0.1% SDS in 10 mM TAPS, 60 °C 2 h in a wet box, "
                 "then 65 °C 1 h with rotation", MaterialKind.CONDITION)
    lysate = g.product("lysate", "reverse-crosslinked tissue lysate", MaterialKind.MOLECULAR,
                       state=ab_barcoded_fragment())
    g.react("reverse_crosslinks", "Lyse the region of interest and reverse crosslinks",
            StepKind.DEPROTEINIZATION, ((ab_prod, InputRole.TEMPLATE),
            (pk, InputRole.CONDITION)), (lysate,))
    p5 = g.input("p5_primer", "New-P5-SMART PCR hybrid oligo", MaterialKind.OLIGO,
                 state=NEW_P5_SMART)
    p7 = g.input("indexed_p7_primer", "Ad2.x indexed primer (Ad2.1 shown)",
                 MaterialKind.OLIGO, state=AD2_1)
    pol = g.input("pcr_polymerase", "NEBNext High-Fidelity 2X PCR Master Mix (25 µL) with "
                  "15 µL 0.67% Triton X-100 neutralization per 5 µL lysate; 2 µL of 10 µM each primer",
                  MaterialKind.CATALYST)
    amplified = g.product("amplified", FINAL_LIBRARY.name, MaterialKind.LIBRARY,
                          state=FINAL_LIBRARY)
    g.react("library_pcr", "58 °C 5 min, 72 °C 5 min gap fill, 12 cycles (60 °C anneal)",
            StepKind.PCR, ((lysate, InputRole.TEMPLATE), (p5, InputRole.STRUCTURAL),
            (p7, InputRole.STRUCTURAL), (pol, InputRole.CATALYST)), (amplified,),
            geometry=PCRGeometry.exact(
                ab_barcoded_fragment(), PCRPrimer("New-P5-SMART", NEW_P5_SMART),
                PCRPrimer("Ad2.1", AD2_1),
                statement="New-P5-SMART ends in the custom Read 1 site carried by barcode B; "
                          "Ad2 ends in S7 plus the first Tn5ME-B mosaic-end bases. The 2022 "
                          "methods name a 'new P5 PCR primer' and 'i7 primers'; Table S1 "
                          "prints New-P5-SMART and Ad2.1–Ad2.8."))
    beads = g.input("ampure", "1.3× AMPure XP", MaterialKind.SELECTION)
    lib = g.product("library", FINAL_LIBRARY.name, MaterialKind.LIBRARY,
                    state=FINAL_LIBRARY)
    sup = g.product("bead_supernatant", "bead supernatant and ethanol washes",
                    MaterialKind.WASTE)
    g.react("bead_cleanup", "1.3× AMPure XP cleanup to remove primers, elute in 10 mM "
            "Tris-HCl pH 8", StepKind.CLEANUP,
            ((amplified, InputRole.TEMPLATE), (beads, InputRole.SELECTION)), (lib, sup),
            changes_molecule=False)
    g.discard(sup, "supernatant and washes removed from the magnet-held beads")
    g.sequence(lib, SequencingPlan(
        library=FINAL_LIBRARY, primers=SEQ_PRIMERS,
        acquired_roles=tuple(p.role for p in SEQ_PRIMERS), read_lengths=READ_LENGTHS,
        read_uses={"Read 1": "extract_then_include", "Index 1 (i7)": "demultiplex",
                   "Read 2": "include"}))
    g.annotate(SourceDiscrepancy(
        "pa_tn5_adapter", "Tn5ME-A loading oligo: 2022 Science versus 2024 Nature Protocols",
        (SourceClaim("Deng et al. 2022 Supplementary Materials, 'Transposome preparation'",
                     f"Tn5ME-A is {TN5_ME_A_WRITTEN}."),
         SourceClaim("Farzad et al. 2024 Supp_Table 1, row 'Tn5ME-A'",
                     f"Tn5ME-A is {TN5_ME_A_TABLE1_WRITTEN}."),
         SourceClaim("Both sources, ligation linker 1",
                     f"Ligation linker 1 is {LIGATION_LINKER_1}.")),
        "The construct uses the 2022 Tn5ME-A from the defining publication; the 2024 "
        "sequence is listed but not used.",
        "Linker 1 is the exact reverse complement of the barcode-A mosaic end followed by "
        f"{LINKER_1_REQUIRED_EXTENSION}, which is the 2022 Tn5ME-A extension; the 2024 "
        f"extension {TABLE1_EXTENSION} has no complementary stretch in linker 1, so a "
        "fragment carrying it could not be splint-ligated to barcode A. All other oligos "
        "and barcodes are identical in the two tables.",
        "Whether the 2024 row is a table error or the 2024 procedure used other loading "
        "oligos is not stated."))
    g.declare_reviewed_scope((
        "tissue", "np40_digitonin", "permeabilize", "permeabilized", "primary_antibody",
        "primary_binding", "primary_bound", "unbound_primary", "secondary_antibody",
        "secondary_binding", "secondary_bound", "unbound_secondary", "pa_tn5_adapter",
        "pa_tn5", "unbound_patn5", "tagmentation_buffer", "tagment",
        "tagged_fragment", "same_adapter_fragments", "barcode_a", "ligase",
        "barcode_a_ligation", "a_barcoded", "barcode_b", "barcode_b_ligation",
        "ab_barcoded", "lysis", "reverse_crosslinks", "lysate", "p5_primer",
        "indexed_p7_primer", "pcr_polymerase", "library_pcr", "amplified", "ampure",
        "bead_cleanup", "library", "bead_supernatant"), auth)
    return g.require_publishable()
