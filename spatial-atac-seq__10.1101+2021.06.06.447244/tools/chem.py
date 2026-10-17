"""Spatial-ATAC-seq orthogonal in-tissue barcoding (Deng et al., Nature 2022).

All oligos below are transcribed from the defining publication: the transposome
oligos from Methods "Preparation of the transposome" (PMC9452302 JATS) and the
linkers, barcodes, PCR and custom sequencing primers from Supplementary Table 3
(MOESM5 workbook, sheets "DNA oligos" and "DNA barcode sequences").

The construct is built locally rather than through ``spatial_epigenome``: that
shared recipe draws linker 1's reverse complement directly against genomic DNA,
omitting the Tn5ME-A transferred strand.  Here the transferred Tn5ME-A strand
carries a 15-nt 5' extension that ligation linker 1 splints to the 3' mosaic end
of the barcode-A oligo; the import-time checks below derive that geometry.
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


DENG = "https://doi.org/10.1038/s41586-022-05094-1"
FARZAD = "https://doi.org/10.1038/s41596-024-01013-y"
TITLE = "Spatial-ATAC-seq — spatial chromatin accessibility"
NOTES = "01_spatial-atac-seq.html"
SOURCE = (f'Defining source: <a href="{DENG}">Deng et al., <i>Nature</i> (2022)</a>, '
          'Methods and Supplementary Table 3 (oligos and barcodes). The authors\' later '
          f'<a href="{FARZAD}">Farzad et al., <i>Nature Protocols</i> (2024)</a> '
          'describes the same workflow.')
SUMMARY = ("Tn5 loaded with a linker-bearing Tn5ME-A adapter tags accessible chromatin in "
           "a fixed tissue section. Two perpendicular microfluidic ligations append 8-nt A "
           "and B barcodes whose ordered pair identifies a tissue pixel.")

# -- Deng et al. 2022, Methods "Preparation of the transposome", 5'->3' as printed.
TN5_ME_REV_WRITTEN = "/5Phos/CTGTCTCTTATACACATCT"
TN5_ME_A_WRITTEN = "/5Phos/CATCGGCGTACGACTAGATGTGTATAAGAGACAG"
TN5_ME_B = "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"
# -- Supplementary Table 3, sheet "DNA oligos".
LIGATION_LINKER_1 = "AGTCGTACGCCGATGCTGTCTCTTATACACATCT"
LIGATION_LINKER_2 = "CGAATGCTCTGGCCTCTCAAGCACGTGGAT"
NEW_P5_SMART_WRITTEN = ("AATGATACGGCGACCACCGAGATCTACACGCCTGTCCGCGGAAGCAG"
                        "TGGTATCAACGCAGAGT*A*C")
AD2_1 = "CAAGCAGAAGACGGCATACGAGATTCGCCTTAGTCTCGTGGGCTCGGAGATGTG"
CUSTOM_READ_1 = "GCCTGTCCGCGGAAGCAGTGGTATCAACGCAGAGTAC"
# -- Supplementary Table 3, sheet "DNA barcode sequences" (first member of each set).
BARCODE_A_OLIGO_1 = "/5Phos/AGGCCAGAGCATTCGAACGTGATAGATGTGTATAAGAGACAG"  # Tn5ME-A_BC1_1
BARCODE_B_OLIGO_1 = ("GCCTGTCCGCGGAAGCAGTGGTATCAACGCAGAGTACAACGTGAT"
                     "ATCCACGTGCTTGAG")                                   # BC2_R1_1

# Farzad et al. 2024 Supplementary Table 1 prints a different Tn5ME-A.
FARZAD_TN5_ME_A_WRITTEN = "/5Phos/TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"


def _bases(written: str) -> str:
    return written.removeprefix("/5Phos/").replace("*", "")


NEW_P5_SMART = _bases(NEW_P5_SMART_WRITTEN)
if NEW_P5_SMART != il.P5 + CUSTOM_READ_1:
    raise ValueError("New-P5-SMART is not P5 followed by the custom Read 1 site")
TN5_ME_A = _bases(TN5_ME_A_WRITTEN)
if not TN5_ME_A.endswith(nx.ME) or revcomp(_bases(TN5_ME_REV_WRITTEN)) != nx.ME:
    raise ValueError("Tn5ME-A / Tn5MErev do not form the mosaic-end duplex")
TN5_ME_A_EXTENSION = TN5_ME_A[:-len(nx.ME)]          # 5' extension beyond the ME
BARCODE_A_BODY = _bases(BARCODE_A_OLIGO_1)
LINKER_2_A_SIDE = BARCODE_A_BODY[:15]
BARCODE_A_1 = BARCODE_A_BODY[15:23]
if BARCODE_A_BODY[23:] != nx.ME:
    raise ValueError("barcode-A oligo does not end in the mosaic end")
BARCODE_B_1 = BARCODE_B_OLIGO_1[len(CUSTOM_READ_1):len(CUSTOM_READ_1) + 8]
LINKER_2_B_SIDE = BARCODE_B_OLIGO_1[len(CUSTOM_READ_1) + 8:]
if not BARCODE_B_OLIGO_1.startswith(CUSTOM_READ_1):
    raise ValueError("barcode-B oligo does not start with the custom Read 1 site")
# Splint relations: each linker must be the exact reverse complement of the two ends
# it joins.  Linker 1 joins barcode-A's 3' mosaic end to the Tn5ME-A extension;
# linker 2 joins barcode-B's 3' end to barcode-A's 5' end.
if revcomp(LIGATION_LINKER_1) != nx.ME + TN5_ME_A_EXTENSION:
    raise ValueError("ligation linker 1 does not splint barcode A to Tn5ME-A")
if revcomp(LIGATION_LINKER_2) != LINKER_2_B_SIDE + LINKER_2_A_SIDE:
    raise ValueError("ligation linker 2 does not splint barcode B to barcode A")
# The 2024 Table 1 adapter shares no complementarity with linker 1's 5' half.
FARZAD_EXTENSION = _bases(FARZAD_TN5_ME_A_WRITTEN)[:-len(nx.ME)]
FARZAD_ADAPTER_PAIRS_LINKER_1 = revcomp(LIGATION_LINKER_1).endswith(FARZAD_EXTENSION)

SPATIAL_A = feature("spatial_a", "spatial_barcode", "whitelist", group="spatial_coordinate",
                    part="A", whitelist="Deng et al. 2022 Supplementary Table 3, 1st barcodes")
SPATIAL_B = feature("spatial_b", "spatial_barcode", "whitelist", group="spatial_coordinate",
                    part="B", whitelist="Deng et al. 2022 Supplementary Table 3, 2nd barcodes")
SAMPLE_I7 = feature("sample_index_i7", "sample_index", "whitelist",
                    whitelist="Deng et al. 2022 Supplementary Table 3, Ad2.1-Ad2.8")
INSERT = "accessible genomic fragment"


def _gdna():
    return seg(INSERT, "X" * 46, placeholder=True)


def _right_arm():
    return [seg("Tn5ME-B mosaic end (rc)", nx.ME_RC, "me"),
            seg("S7 (rc)", nx.S7_RC, "s7")]


def _a_side():
    return [seg("barcode-A mosaic end", nx.ME, "me"),
            seg("Tn5ME-A extension", TN5_ME_A_EXTENSION, "me"),
            seg("Tn5ME-A mosaic end", nx.ME, "me"),
            _gdna(), *_right_arm()]


def tagged_fragment() -> Construct:
    return Construct([seg("Tn5ME-A extension", TN5_ME_A_EXTENSION, "me"),
                      seg("Tn5ME-A mosaic end", nx.ME, "me"), _gdna(), *_right_arm()],
                     name="Tn5ME-A/Tn5ME-B tagmented fragment")


def a_barcoded_fragment() -> Construct:
    return Construct([seg("linker-2 A-side", LINKER_2_A_SIDE),
                      seg("spatial barcode A", "A" * 8, "cbc", placeholder=True,
                          feature=SPATIAL_A),
                      *_a_side()], name="A-barcoded tagmented fragment")


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
                    name=INSERT + " spatial-ATAC library")
    primers = (sp.custom("Read 1", "Custom Read 1 primer", CUSTOM_READ_1,
                         "Deng et al. 2022 Supplementary Table 3"),
               sp.NEXTERA["I1"], sp.NEXTERA["R2"])
    problems = sp.verify(lib, primers, required_roles=tuple(p.role for p in primers))
    if problems:
        raise ValueError("invalid Spatial-ATAC library: " + "; ".join(problems))
    return lib, primers


FINAL_LIBRARY, SEQ_PRIMERS = library()
READ_LENGTHS = {"Read 1": 150, "Index 1 (i7)": 8, "Read 2": 150}
SHOW_READ_CONTENT = True
def _read1_start(name: str) -> int:
    """1-based raw Read 1 cycle of a segment, counted from the custom primer's 3' end."""
    segs = list(FINAL_LIBRARY)
    names = [s.name for s in segs]
    first = names.index("custom Read 1 site") + 1
    return sum(len(s.top) for s in segs[first:names.index(name)]) + 1


GENOMIC_READ1_START = _read1_start(INSERT)
FINAL_CAPTION = ("The custom-primer end reads barcode B, linker 2, barcode A, the "
                 "barcode-A mosaic end, the Tn5ME-A extension and its mosaic end before "
                 "accessible genomic DNA. ** marks the two in-tissue ligations.")
SEQUENCING_INTRO = ("HiSeq 4000 or NovaSeq 6000 paired-end 150 sequencing with a custom "
                    "Read 1 primer. Raw Read 1 carries barcode B at cycles 1–8 and barcode "
                    f"A at cycles 39–46; genomic DNA starts at cycle {GENOMIC_READ1_START}. "
                    "Read 2 reports the opposite genomic end.")


def oligos():
    return (f"<b>Tn5MErev</b> — 5′-{TN5_ME_REV_WRITTEN}-3′",
            f"<b>Tn5ME-A</b> — 5′-{TN5_ME_A_WRITTEN}-3′",
            f"<b>Tn5ME-B</b> — 5′-{TN5_ME_B}-3′",
            f"<b>Ligation linker 1</b> — 5′-{LIGATION_LINKER_1}-3′",
            f"<b>Ligation linker 2</b> — 5′-{LIGATION_LINKER_2}-3′",
            f"<b>Tn5ME-A_BC1_1</b> (barcode A, 1 of 50) — 5′-{BARCODE_A_OLIGO_1}-3′",
            f"<b>BC2_R1_1</b> (barcode B, 1 of 50) — 5′-{BARCODE_B_OLIGO_1}-3′",
            f"<b>New-P5-SMART PCR hybrid oligo</b> — 5′-{NEW_P5_SMART_WRITTEN}-3′",
            f"<b>Ad2.1_TAAGGCGA</b> (1 of 8) — 5′-{AD2_1}-3′",
            f"<b>Custom Read 1 primer</b> — 5′-{CUSTOM_READ_1}-3′")


def _scene(con: Construct, label: str, junctions=()):
    sc = Scene.duplex(list(con), label=label)
    for left, right, text in junctions:
        sc.junction("top", left, right, text)
    sc.labels("top")
    return sc


_J_A = ("barcode-A mosaic end", "Tn5ME-A extension", "ligation A")
_J_B = ("linker-2 B-side", "linker-2 A-side", "ligation B")


def final_rows():
    return _scene(FINAL_LIBRARY, "library", (_J_B, _J_A)).rows()


def protocol_graph() -> ProtocolGraph:
    g = ProtocolGraph("Spatial-ATAC-seq", authority=DENG + " (Methods; Supplementary Table 3)")
    tissue = g.input("tissue", "fresh-frozen tissue section", MaterialKind.BIOLOGICAL,
                     primary=True)
    fa = g.input("formaldehyde", "0.2% formaldehyde, 5 min; 1.25 M glycine quench",
                 MaterialKind.CONDITION)
    fixed = g.product("fixed_tissue", "fixed tissue section", MaterialKind.BIOLOGICAL)
    g.react("fix", "Fix with 0.2% formaldehyde and quench with glycine",
            StepKind.CROSS_LINKING, ((tissue, InputRole.TEMPLATE), (fa, InputRole.CONDITION)),
            (fixed,))
    lysis = g.input("lysis_buffer", "lysis buffer (NP-40, Tween-20, digitonin), 15 min",
                    MaterialKind.CONDITION)
    perm = g.product("permeabilized", "permeabilized tissue chromatin", MaterialKind.BIOLOGICAL)
    g.react("permeabilize", "Permeabilize the section and wash", StepKind.PERMEABILIZATION,
            ((fixed, InputRole.TEMPLATE), (lysis, InputRole.CONDITION)), (perm,))
    transposome = g.input("transposome", "Diagenode Tn5 loaded with Tn5ME-A/Tn5MErev and "
                          "Tn5ME-B/Tn5MErev", MaterialKind.ADAPTER,
                          state=(TN5_ME_A, TN5_ME_B, _bases(TN5_ME_REV_WRITTEN)))
    tn5 = g.input("tn5", "Tn5 transposase (Diagenode C01070010)", MaterialKind.CATALYST)
    tagged = g.product("tagged_fragment", "Tn5ME-A/Tn5ME-B tagmented fragment",
                       MaterialKind.MOLECULAR, state=tagged_fragment())
    same = g.product("same_adapter_fragments",
                     "candidate Tn5ME-A/A and Tn5ME-B/B fragments", MaterialKind.WASTE,
                     inferred=True, source=DENG + " (end classes derived from the mixed "
                     "transposome; not described in the source)")
    g.react("tagment", "Bulk in-tissue tagmentation, 37 °C 30 min; stop with 40 mM EDTA",
            StepKind.TAGMENTATION, ((perm, InputRole.TEMPLATE),
            (transposome, InputRole.STRUCTURAL), (tn5, InputRole.CATALYST)), (tagged, same))
    g.mark_terminal(Terminal(same.id, Disposition.UNRESOLVED,
        "Candidate same-adapter fragments stay in the tissue through both ligations and "
        "lysis; no physical separation is described. They lack the P5/P7 primer pair "
        "needed for exponential amplification, which is not a physical removal."))
    linker1 = g.input("barcode_a", "barcode A1–A50 annealed to ligation linker 1",
                      MaterialKind.ADAPTER, state=(BARCODE_A_BODY, LIGATION_LINKER_1))
    ligase = g.input("ligase", "T4 DNA ligase (NEB M0202L)", MaterialKind.CATALYST)
    a_prod = g.product("a_barcoded", "A-barcoded tagmented fragment", MaterialKind.MOLECULAR,
                       state=a_barcoded_fragment())
    g.react("barcode_a_ligation", "Flow barcode A through channels 1–50 and ligate, 37 °C 30 min",
            StepKind.LIGATION, ((tagged, InputRole.TEMPLATE), (linker1, InputRole.STRUCTURAL),
            (ligase, InputRole.CATALYST)), (a_prod,),
            geometry=reported_adapter_end(
                "barcode-A oligo", kind="splinted", polarity="3-prime",
                terminal_sequence=nx.ME,
                statement="Ligation linker 1 is the reverse complement of the barcode-A "
                          "3' mosaic end followed by the Tn5ME-A 15-nt 5' extension, so it "
                          "splints barcode A's 3'-OH to the 5'-phosphorylated transferred "
                          "Tn5ME-A strand."))
    linker2 = g.input("barcode_b", "barcode B1–B50 annealed to ligation linker 2",
                      MaterialKind.ADAPTER, state=(BARCODE_B_OLIGO_1, LIGATION_LINKER_2))
    ab_prod = g.product("ab_barcoded", "A/B-barcoded tagmented fragment",
                        MaterialKind.MOLECULAR, state=ab_barcoded_fragment())
    g.react("barcode_b_ligation", "Rotate the device, flow barcode B and ligate, 37 °C 30 min",
            StepKind.LIGATION, ((a_prod, InputRole.TEMPLATE), (linker2, InputRole.STRUCTURAL),
            (ligase, InputRole.CATALYST)), (ab_prod,),
            geometry=reported_adapter_end(
                "barcode-B oligo", kind="splinted", polarity="3-prime",
                terminal_sequence=LINKER_2_B_SIDE,
                statement="Ligation linker 2 is the reverse complement of the barcode-B "
                          "3' end followed by the 5'-phosphorylated barcode-A 5' end."))
    pk = g.input("proteinase_k", "proteinase K in 1% SDS, 58 °C 2 h then 65 °C overnight",
                 MaterialKind.CATALYST)
    lysate = g.product("lysate", "reverse-crosslinked tissue lysate", MaterialKind.MOLECULAR,
                       state=ab_barcoded_fragment())
    g.react("reverse_crosslinks", "Digest the region of interest and reverse crosslinks",
            StepKind.DEPROTEINIZATION, ((ab_prod, InputRole.TEMPLATE),
            (pk, InputRole.CATALYST)), (lysate,))
    zymo = g.input("zymo_column", "Zymo DNA Clean & Concentrator-5", MaterialKind.SELECTION)
    purified = g.product("purified_dna", "purified spatial fragments", MaterialKind.MOLECULAR,
                         state=ab_barcoded_fragment())
    flow = g.product("column_flowthrough", "column flow-through and washes (protein, SDS)",
                     MaterialKind.WASTE)
    g.react("column_cleanup", "Purify DNA and elute in 20 µl", StepKind.CLEANUP,
            ((lysate, InputRole.TEMPLATE), (zymo, InputRole.SELECTION)), (purified, flow),
            changes_molecule=False)
    g.discard(flow, "column flow-through and washes are discarded")
    p5 = g.input("p5_primer", "New-P5-SMART PCR hybrid oligo", MaterialKind.OLIGO,
                 state=NEW_P5_SMART)
    p7 = g.input("indexed_p7_primer", "Ad2.x indexed primer (Ad2.1 shown)",
                 MaterialKind.OLIGO, state=AD2_1)
    pol = g.input("pcr_polymerase", "NEBNext High-Fidelity 2× PCR Master Mix",
                  MaterialKind.CATALYST)
    amplified = g.product("amplified", FINAL_LIBRARY.name, MaterialKind.LIBRARY,
                          state=FINAL_LIBRARY)
    g.react("library_pcr", "72 °C 5 min gap fill, 5 cycles, then qPCR-determined extra cycles",
            StepKind.PCR, ((purified, InputRole.TEMPLATE), (p5, InputRole.STRUCTURAL),
            (p7, InputRole.STRUCTURAL), (pol, InputRole.CATALYST)), (amplified,),
            geometry=PCRGeometry.exact(
                ab_barcoded_fragment(), PCRPrimer("New-P5-SMART", NEW_P5_SMART),
                PCRPrimer("Ad2.1", AD2_1),
                statement="New-P5-SMART ends in the custom Read 1 site carried by barcode B; "
                          "Ad2 ends in S7 plus the first bases of the Tn5ME-B mosaic end. "
                          "The 72 °C step fills the Tn5 gaps before cycling."))
    beads = g.input("ampure", "1× AMPure XP beads", MaterialKind.SELECTION)
    lib = g.product("library", FINAL_LIBRARY.name, MaterialKind.LIBRARY, state=FINAL_LIBRARY)
    sup = g.product("bead_supernatant", "bead supernatant (PCR primer residues)",
                    MaterialKind.WASTE)
    g.react("bead_cleanup", "1× AMPure XP cleanup, elute in 20 µl", StepKind.CLEANUP,
            ((amplified, InputRole.TEMPLATE), (beads, InputRole.SELECTION)), (lib, sup),
            changes_molecule=False)
    g.discard(sup, "supernatant removed from the magnet-held beads")
    g.sequence(lib, SequencingPlan(
        library=FINAL_LIBRARY, primers=SEQ_PRIMERS,
        acquired_roles=tuple(p.role for p in SEQ_PRIMERS), read_lengths=READ_LENGTHS,
        read_uses={"Read 1": "extract_then_include", "Index 1 (i7)": "demultiplex",
                   "Read 2": "include"}))
    g.annotate(SourceDiscrepancy(
        "transposome", "Tn5ME-A loading oligo",
        (SourceClaim("Deng et al. 2022 Methods, 'Preparation of the transposome'",
                     f"Tn5ME-A is {TN5_ME_A_WRITTEN}."),
         SourceClaim("Farzad et al. 2024 Supplementary Table 1 (used in Procedure "
                     "Steps 16–24, transposome assembly for spatial-ATAC-seq)",
                     f"Tn5ME-A is {FARZAD_TN5_ME_A_WRITTEN}.")),
        "The model uses the defining 2022 Methods oligo.",
        "Ligation linker 1 (printed identically in both sources) is the exact reverse "
        "complement of the barcode-A mosaic end plus CATCGGCGTACGACT, the 5' extension of "
        "the 2022 Tn5ME-A; the 2024 Table 1 extension TCGTCGGCAGCGTC has no complementary "
        "stretch in linker 1.",
        "Whether the 2024 Table 1 entry describes a different experiment or is a table "
        "error is not stated."))
    g.declare_reviewed_scope((
        "tissue", "formaldehyde", "fix", "fixed_tissue", "lysis_buffer", "permeabilize",
        "permeabilized", "transposome", "tn5", "tagment", "tagged_fragment",
        "same_adapter_fragments", "barcode_a", "ligase", "barcode_a_ligation", "a_barcoded",
        "barcode_b", "barcode_b_ligation", "ab_barcoded", "proteinase_k",
        "reverse_crosslinks", "lysate", "zymo_column", "column_cleanup", "purified_dna",
        "column_flowthrough", "p5_primer", "indexed_p7_primer", "pcr_polymerase",
        "library_pcr", "amplified", "ampure", "bead_cleanup", "library",
        "bead_supernatant"), g.authority)
    return g.require_publishable()
