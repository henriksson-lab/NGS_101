"""Zhu et al. 2021 PairedTag DNA/RNA library architecture."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import illumina as il
import nextera as nx
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, MolecularState, Row, Scene, Segment, Workflow, feature, revcomp


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name, top, tag, **kw)


TITLE = "PairedTag — joint histone-mark and transcriptome profiling"
NOTES = "01_pairedtag.html"
SOURCE = ('Defining source: <a href="https://doi.org/10.1038/s41592-021-01060-3">'
          'Zhu et al., <i>Nature Methods</i> (2021)</a>.')
SUMMARY = ("Antibody-tethered, barcoded pA–Tn5 marks chromatin while barcoded reverse "
           "transcription copies nuclear RNA. Two subsequent split-pool ligations give "
           "both modalities the same three-part nuclear identity before separate library "
           "finishing reactions.")
CAVEAT = ("The final panels preserve the published read geometry and barcode lengths. "
          "Barcode letters are role placeholders because each round uses a plate of many "
          "published sequences; canonical Illumina primer sites are shown exactly.")

PMENTS = "CTGTCTCTTATACACATCT"
ADAPTOR_A = nx.ADAPTOR_S5
LINKER_R02 = "CGAATGCTCTGGCCTCTCAAGCACGTGGAT"
LINKER_R03 = "GGTCTGAGTTCGCACCGAAACATCGGCCAC"
P5_FOKI = il.TRUSEQ_READ1
PA_F = "CAGACGTGTGCTCTTCCGATCT"
PA_R = "AAGCAGTGGTATCAACGCAGAGT"
PUBLISHED_READ2_BARCODE_WINDOWS = ((10, 13), (47, 50), (84, 87))


def barcode_cassette(modality: str) -> list[Segment]:
    """Top-strand order; Read 2 traverses this list from right to left."""
    return [
        seg(f"round-1 {modality} barcode", "B" * 6, "cbc", placeholder=True,
            feature=feature(f"cell_bc1_{modality.lower()}", "cell_barcode", "combinatorial",
                            group="cell_id", part="round 1")),
        seg("round-2 linker plus ligation boundary", "L" * 31, placeholder=True),
        seg("round-2 barcode", "B" * 6, "cbc", placeholder=True,
            feature=feature("cell_bc2", "cell_barcode", "combinatorial",
                            group="cell_id", part="round 2")),
        seg("round-3 linker", "L" * 30, placeholder=True),
        seg("round-3 barcode", "B" * 6, "cbc", placeholder=True,
            feature=feature("cell_bc3", "cell_barcode", "combinatorial",
                            group="cell_id", part="round 3")),
        seg("10-nt UMI", "U" * 10, "umi", placeholder=True,
            feature=feature("umi", "umi", "random")),
    ]


def _right_end() -> list[Segment]:
    return [seg("adapter junction", "A"),
            seg("Index 1 / Read 2 arm", il.INDEX1_PRIMER, "r2"),
            seg("i7 reverse complement", "I" * 6, "cbc", placeholder=True,
                feature=feature("sample_i7", "sample_index", "unknown")),
            seg("P7 reverse complement", il.P7_RC, "p7")]


DNA_LIBRARY = Construct([
    seg("P5", il.P5, "p5"), seg("TruSeq Read 1", il.TRUSEQ_READ1, "r1"),
    seg("histone-mark DNA", "X" * 34, placeholder=True),
    *barcode_cassette("DNA"), *_right_end()], name="PairedTag DNA library")
DNA_PRIMERS = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["R2"])

RNA_LIBRARY = Construct([
    seg("P5", il.P5, "p5"),
    seg("N5 carried index", "J" * 8, "cbc", placeholder=True,
        feature=feature("sample_i5", "sample_index", "unknown")),
    seg("S5", nx.S5, "s5"), seg("mosaic end", nx.ME, "me"),
    seg("transcript cDNA", "X" * 34, placeholder=True),
    *barcode_cassette("RNA"), *_right_end()], name="PairedTag RNA library")
RNA_PRIMERS = (sp.NEXTERA["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["R2"])

for library, primers in ((DNA_LIBRARY, DNA_PRIMERS), (RNA_LIBRARY, RNA_PRIMERS)):
    errors = sp.verify(library, primers, required_roles=tuple(p.role for p in primers))
    if errors:
        raise ValueError(f"{library.name}: " + "; ".join(errors))


def read2_barcode_starts(library: Construct, modality: str) -> tuple[int, int, int]:
    """1-based first cycles for BC3, BC2 and BC1, derived from the final construct."""
    hit = sp.locate(library, sp.TRUSEQ["R2"])
    if hit is None or hit.strand != "top":
        raise ValueError(f"{library.name}: Read 2 must extend leftward from the top strand")
    names = ("round-3 barcode", "round-2 barcode", f"round-1 {modality} barcode")
    return tuple(hit.start - library.span(name)[1] + 1 for name in names)


for library, modality in ((DNA_LIBRARY, "DNA"), (RNA_LIBRARY, "RNA")):
    starts = read2_barcode_starts(library, modality)
    if not all(lo <= got <= hi for got, (lo, hi) in
               zip(starts, PUBLISHED_READ2_BARCODE_WINDOWS)):
        raise ValueError(f"{library.name}: barcode geometry {starts} misses published windows")

FINAL_LIBRARIES = (
    ("Histone-mark DNA library", DNA_LIBRARY, DNA_PRIMERS,
     "Read 1 maps the antibody-targeted genomic fragment. Read 2 reports UMI and the "
     "three ligated nuclear barcodes; Index 1 reports the PCR-added sample index.",
     "The DNA branch receives a P5-FokI adapter after SbfI/FokI digestion."),
    ("Transcriptome library", RNA_LIBRARY, RNA_PRIMERS,
     "Read 1 maps transcript cDNA from the N5/ME end. Read 2 reports UMI and the same "
     "three-round nuclear identity; Index 1 reports the sample index.",
     "The RNA branch receives its second sequencing end by N5-loaded Tn5 after NotI digestion."),
)


def sections():
    wf = Workflow(MolecularState("Nuclei with chromatin and RNA", (Row(chunks=[("chromatin + nuclear RNA", None, False)]),)))
    wf.react("Mark chromatin and RNA with matched round-1 identities", [
            Row(chunks=[("antibody — pA–Tn5 — DNA adapter/BC1 ** histone-mark DNA", "me", False)]),
            Row(chunks=[("RNA primer/BC1 ** poly(A) or random-primed nuclear RNA → cDNA", "r1", False)]),
        ], note="Twelve first-round wells receive a matching DNA transposome barcode and RNA RT-primer barcode.")
    for n, linker in ((2, LINKER_R02), (3, LINKER_R03)):
        rc = revcomp(linker)
        sc = Scene(); sc.strand("splint", [seg(f"linker R0{n} new half", linker[15:]), seg(f"linker R0{n} old half", linker[:15])], label=f"Linker-R0{n}")
        sc.anneal("barcode ends", [seg(f"round-{n} 5-prime arm", rc[:15]), seg("previous barcode 3-prime arm", rc[15:])], to="splint", pair=(f"round-{n} 5-prime arm", f"linker R0{n} new half"), label="adjacent barcode oligos", mod5="p")
        sc.ligation_end("barcode ends", f"round-{n} 5-prime arm", kind="splinted", polarity="5-prime", owner="adapter", text="5′-phosphorylated nick partner")
        geometry = AdapterLigationGeometry.known("splinted", sc, statement=f"The supplement prints Linker-R0{n}; its two 15-base halves anneal the adjacent barcode oligos across the ligation nick.")
        wf.ligate_adapter(f"Ligate barcode round {n}", geometry, [Row(chunks=[(f"molecule — BC1 ** " + " ** ".join(f"BC{x}" for x in range(2, n + 1)), "cbc", False)])], note="Nuclei are pooled and redistributed between barcode plates.")
    wf.react("Amplify together, then split by modality", [
            Row(chunks=[("TdT dC tail → anchor extension → common pre-amplification", None, False)]),
            Row(chunks=[("RNA: NotI → N5-Tn5 tagmentation", "r2", False)]),
        ], note="The shared barcode cassette is copied before the material is divided; distinct finishing reactions create the two libraries.")
    wf.ligate_adapter("Ligate the P5 adapter to the SbfI/FokI-cut DNA branch",
        AdapterLigationGeometry.unavailable("The PairedTag supplement identifies the P5-FokI finishing adapter and restriction steps, but the current source transcription does not establish a ligation-ready adapter overhang sequence and polarity."),
        [Row(chunks=[("DNA: SbfI + FokI → P5-adapter ** histone-mark insert", "r1", False)])],
        note="The unavailable state is kept separate from the RNA branch's ligation-free N5-Tn5 finishing reaction.")
    return [(r.action, r.after.rows, r.note) for r in wf.reactions]
