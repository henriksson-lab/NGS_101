"""Light-directed covalent spatial barcoding and Light-Seq library chemistry.

Light-Seq's final library is a deliberately asymmetric Nextera product: a custom
S5XX.GATE primer selects the original GATE end while a normal N7XX primer selects a
tagmented end.  Keeping that assembly here makes the read structure a consequence of
the published oligos rather than a prose assertion in the page builder.
"""
from __future__ import annotations

from types import SimpleNamespace

import illumina as il
import nextera as nx
import seqprimers as sp
from batch_ngs import seg
from chemdraw import Construct, feature, inline_sequence, revcomp
from protocol_model import chain, molecule, role


RT_5N_3G = "TTTACACGATTGAGTTATNNNNNGGG"
GATE = "GGAGTTGGAGTGAGTGGATGAGTGATG"
GATC = "GAGAATGTGAGTGAAGATGTATGGTGA"
GATC_20T = GATC + "T" * 20

# Verbatim base strings from the protocols.io oligo table. Modifications are kept in
# the display strings below because cnvK and the 3-prime dyes are not ordinary bases.
BARCODE_1 = GATE + "D" * 12 + "TATGGATGAGTTATATAACTCA"
BARCODE_2 = GATE + "D" * 12 + "GTTAGGTGAGTTATATAACTCA"
BARCODE_3 = GATE + "D" * 12 + "AGGGTATGAGTTATATAACTCA"
BARCODE_1_AS_WRITTEN = BARCODE_1 + "[cnvK]TCGTGTAAAT[Cy5-3]"
BARCODE_2_AS_WRITTEN = BARCODE_2 + "[cnvK]TCGTGTAAAT[Cy3-3]"
BARCODE_3_AS_WRITTEN = BARCODE_3 + "[cnvK]TCGTGTAAAT[Fl-3]"

READ1_PRIMER = "CGCC" + GATE
INDEX2_PRIMER = revcomp(READ1_PRIMER)
S502_INDEX = "CTCTCTAT"
S502_GATE = il.P5 + S502_INDEX + READ1_PRIMER
N701_INDEX = "TCGCCTTA"
N701 = il.P7 + N701_INDEX + nx.S7


UMI = feature("light_seq.umi", "umi", "random",
              note="12 IUPAC D positions in the published barcode strand")
ROI = feature("light_seq.roi", "spatial_barcode", "fixed",
              note="six-base ROI identity parsed by the authors' analysis code")
I5 = feature("light_seq.i5", "sample_index", "fixed",
             note="representative S502 index")
I7 = feature("light_seq.i7", "sample_index", "fixed",
             note="representative N701 index")


def final_library() -> tuple[Construct, tuple]:
    """Representative S502/N701 library, checked against all four run primers."""
    lib = Construct([
        seg("P5", il.P5, "p5"),
        seg("i5 S502", S502_INDEX, "cbc", feature=I5),
        seg("custom Read 1 landing prefix", "CGCC", "r1"),
        seg("GATE", GATE, "r1"),
        seg("12-nt UMI", "D" * 12, "umi", placeholder=True, feature=UMI),
        seg("ROI barcode 1", "TATGGA", "cbc", feature=ROI),
        seg("barcode constant", "TGAGTTAT"),
        seg("cross-junction P domain", "ATAACTCA"),
        role("captured transcript cDNA", length_bp=40),
        seg("right mosaic end reverse complement", nx.ME_RC, "me"),
        seg("S7 reverse complement", nx.S7_RC, "s7"),
        seg("i7 N701 reverse complement", revcomp(N701_INDEX), "cbc", feature=I7),
        seg("P7 reverse complement", il.P7_RC, "p7"),
    ], name="representative Light-Seq S502/N701 library")
    primers = (
        sp.custom("Read 1", "P5.GATE custom Read 1", READ1_PRIMER,
                  "Light-Seq protocols.io protocol, oligo table"),
        sp.NEXTERA["I1"],
        sp.custom("Index 2 (i5)", "GATE*.P5* custom i5", INDEX2_PRIMER,
                  "Light-Seq protocols.io protocol, oligo table",
                  "The protocol reports inefficient i5 index sequencing on NovaSeq and recommends unique i7 indexes."),
        sp.NEXTERA["R2"],
    )
    problems = sp.verify(lib, primers, required_roles=tuple(p.role for p in primers))
    if problems:
        raise ValueError("Light-Seq final library: " + "; ".join(problems))
    return lib, primers


def workflow():
    spatial = feature("light_seq.roi.workflow", "spatial_barcode", "fixed")
    umi = feature("light_seq.umi.workflow", "umi", "random")
    return chain(molecule("fixed, permeabilized tissue RNA", (
        role("RNA template"),)), (
        ("Prime total RNA in situ with RT.5N.3G and reverse-transcribe",
         molecule("RNA–cDNA hybrid", (
             role("RNA template"), role("cDNA"), role("5-prime docking site"))),
         "Five random bases followed by three G bases prime RNAs without requiring poly(A); the primer's 5-prime constant becomes the barcode docking site."),
        ("Add a 3-prime poly(dA) tail to cDNA",
         molecule("A-tailed RNA–cDNA hybrid", (
             role("RNA template"), role("cDNA"), role("poly(dA)"),
             role("5-prime docking site"))),
         "Terminal transferase creates the handle used by GATC.20T during cross-junction synthesis."),
        ("Hybridize a cnvK barcode and illuminate only the selected ROI",
         molecule("photocrosslinked barcode–cDNA complex", (
             role("GATE PCR handle"),
             role("12-nt UMI", length_bp=12, molecular_feature=umi),
             role("ROI barcode", length_bp=6, molecular_feature=spatial),
             role("cnvK covalent photocrosslink"), role("docking site"),
             role("cDNA"), role("RNA template"))),
         "UV makes an interstrand covalent bond only inside the optical mask. A stringent formamide wash removes hybridized but non-crosslinked barcode strands; orthogonal barcodes can label further ROIs."),
        ("Digest RNA with RNase H and collect crosslinked barcode–cDNA",
         molecule("released crosslinked DNA complex", (
             role("barcode strand"), role("cnvK covalent photocrosslink"),
             role("cDNA"), role("poly(dA)"))),
         "The sample remains on the slide. RNase H releases the selected cDNAs by cleaving their RNA templates."),
        ("Extend GATC.20T across the cross-junction",
         molecule("single stitched readout strand", (
             role("GATC handle"), role("captured transcript cDNA"),
             role("cross-junction P domain"),
             role("ROI barcode", length_bp=6, molecular_feature=spatial),
             role("12-nt UMI", length_bp=12, molecular_feature=umi),
             role("GATE handle"))),
         "Bst extends to the crosslink, then the shared P domain forms a three-way junction and branch-migrates onto the barcode template. The new product contains barcode and cDNA on one uncrosslinked strand."),
        ("Amplify with GATE and GATC primers",
         molecule("Light-Seq amplicon", (
             role("GATE"), role("12-nt UMI", length_bp=12, molecular_feature=umi),
             role("ROI barcode", length_bp=6, molecular_feature=spatial),
             role("captured transcript cDNA"), role("GATC"))),
         "The authors choose the endpoint from a small qPCR aliquot to avoid over-amplification."),
        ("Tagment, then select GATE-to-Nextera fragments with S5XX.GATE and N7XX",
         molecule("indexed Light-Seq library", (
             role("P5 and i5"), role("custom GATE Read 1 site"),
             role("12-nt UMI", length_bp=12, molecular_feature=umi),
             role("ROI barcode", length_bp=6, molecular_feature=spatial),
             role("captured transcript cDNA"), role("Nextera S7, i7 and P7"))),
         "Unlike ordinary Nextera XT, the custom i5-side primer anchors PCR at the original GATE end. The opposite end is a random Nextera tagmentation boundary."),
    ))


def oligos():
    return (
        f"RT.5N.3G — {inline_sequence(RT_5N_3G, 'RT.5N.3G')}",
        f"GATE.D12.B1, as ordered — {BARCODE_1_AS_WRITTEN}",
        f"GATE.D12.B2, as ordered — {BARCODE_2_AS_WRITTEN}",
        f"GATE.D12.B3, as ordered — {BARCODE_3_AS_WRITTEN}",
        f"GATC.20T — {inline_sequence(GATC_20T, 'GATC.20T')}",
        f"GATE PCR primer — {inline_sequence(GATE, 'GATE PCR primer')}",
        f"GATC PCR primer — {inline_sequence(GATC, 'GATC PCR primer')}",
        f"P5.GATE custom Read 1 — {inline_sequence(READ1_PRIMER, 'P5.GATE')}",
        f"GATE*.P5* custom i5 — {inline_sequence(INDEX2_PRIMER, 'GATE*.P5*')}",
        f"S502.GATE example — {inline_sequence(S502_GATE, 'S502.GATE')}",
        f"Next.N701 example — {inline_sequence(N701, 'Next.N701')}",
    )


def get():
    final, primers = final_library()
    return SimpleNamespace(
        TITLE="Light-Seq — light-directed spatial transcriptome barcoding",
        NOTES="01_light_seq.html",
        SOURCE='<a href="https://doi.org/10.1038/s41592-022-01604-1">Kishi et al., <i>Nature Methods</i> (2022)</a>; detailed workflow and oligos in <a href="https://doi.org/10.17504/protocols.io.x54v9jno4g3e/v1">the authors\' protocol</a>.',
        SUMMARY="A random-primed cDNA docking site receives a cnvK-bearing barcode only where UV illuminates the specimen. Cross-junction synthesis converts the covalently linked barcode and cDNA into one strand; a custom asymmetric Nextera XT PCR preserves the GATE end so Read 1 reports UMI, ROI and transcript.",
        CAVEAT="The final molecule is one informative S502/N701 example. S503 and S505–S508 change only the i5 sample index; N702–N706 change only i7. The protocol reports weak custom-i5 sequencing on NovaSeq and recommends unique i7 indexes for demultiplexing.",
        oligos=oligos, workflow=workflow,
        FINAL_LIBRARY=final,
        FINAL_CAPTION="Representative S502/N701 Light-Seq molecule. The custom Read 1 primer ends immediately before the 12-base UMI; the authors' parser then recognizes the six-base ROI code, fixed TGAGTTAT and eight-base P domain before transcript sequence. The right boundary is random Nextera tagmentation.",
        SEQ_PRIMERS=primers,
        READ_LENGTHS=None,
        SEQUENCING_INTRO="Libraries were sequenced on HiSeq, NovaSeq and MiSeq with the custom P5.GATE Read 1 primer and, where supported, GATE*.P5* i5 primer. Read 1 begins UMI–ROI barcode–cross-junction constants–transcript; Read 2 approaches the same insert from the random Nextera end.",
        SEQUENCING_ENDING="", SEQUENCING_UNAVAILABLE="",
    )
