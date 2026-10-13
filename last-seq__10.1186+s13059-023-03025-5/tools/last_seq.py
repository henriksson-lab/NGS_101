"""Molecular model for LAST-seq (Lyu and Chen, 2023)."""
from __future__ import annotations

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import illumina as il
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry
from chemdraw import (Construct, MolecularState, Row, Scene, Segment, Workflow,
                      annotation_rows, feature, oligo)
from rna_amplification import direct_ssrna_ivt, poly_a_capture
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from protocol_model import molecule, role
from topology import HairpinTopology


def seg(name, bases, tag=None, **kw):
    return Segment(name, bases, tag, **kw)


T7_PROMOTER = "TAATACGACTCACTATAGGG"
RANDOM_RT_HANDLE = "TACACGACGCTCTTCCGATCT"
RANDOM_RT = RANDOM_RT_HANDLE + "NNNNNN"
DTRU_LINKER = "CCACCTTTCATTCACCCT"
DTRU_DNA = "TTTTTT"
DTRU_RNA = "U" * 19
DTRU_TAIL = DTRU_LINKER + DTRU_DNA + DTRU_RNA
HAIRPIN_LINKER_ARM = "GGGTGAATGAAAGGTGG"
HAIRPIN_CORE = "CTCTAATACGACTCACTATAGGGAGACGTGTGCTCTTCCGATCT"
HAIRPIN_CLOSING_ARM = "TCGTATTAGAG"
CELL_BARCODE_LEN = 6
UMI_LEN = 8
BC2 = "GTATCC"
BC2_HAIRPIN = HAIRPIN_LINKER_ARM + HAIRPIN_CORE + "N" * UMI_LEN + BC2 + HAIRPIN_CLOSING_ARM
BC2_HAIRPIN_TOPOLOGY = HairpinTopology(
    "LAST-seq BC2 hairpin module", BC2_HAIRPIN,
    first_stem_start=len(HAIRPIN_LINKER_ARM),
    second_stem_start=len(BC2_HAIRPIN) - len(HAIRPIN_CLOSING_ARM),
    stem_length=len(HAIRPIN_CLOSING_ARM), colour="var(--c-cbc)")

CELL_FEATURE = feature("last_cell_barcode", "cell_barcode", "whitelist",
                       whitelist="Supplementary Table S3")
UMI_FEATURE = feature("last_umi", "umi", "random")


def capture_primer():
    return [
        seg("hairpin module", "[hairpin: PCR handle + T7 promoter + UMI + cell barcode]",
            "r2", placeholder=True),
        seg("complementary linker", DTRU_LINKER, "tso"),
        seg("dT7", DTRU_DNA),
        seg("rU19", DTRU_RNA),
    ]


def capture_scene():
    sc = poly_a_capture(capture_primer(), annealing_segment="rU19",
                        primer_label="3′-ddC-blocked LAST-seq primer")
    sc.junction("primer", "hairpin module", "complementary linker", "primer assembly ligation")
    sc.note("primer", "the dT/rU tract is nicked by RNase H; 3′ ddC blocks extension from the primer end")
    return sc


def primer_ligation_geometry():
    sc = Scene(); sc.strand("hairpin", [seg("hairpin linker arm", HAIRPIN_LINKER_ARM),
                                         seg("folded barcode/promoter hairpin", "[folded hairpin]", placeholder=True)],
                                    label="hairpin module")
    sc.anneal("rU-dT module", [seg("phosphorylated paired linker", DTRU_LINKER[:-1]),
                                seg("linker 3-prime dT overhang", DTRU_LINKER[-1])],
              to="hairpin", pair=("phosphorylated paired linker", "hairpin linker arm"),
              label="rU-dT module", mod5="p")
    sc.ligation_end("rU-dT module", "phosphorylated paired linker", kind="splinted",
                    polarity="5-prime", owner="adapter", text="5′-phosphorylated nick end")
    sc.overhang("rU-dT module", "linker 3-prime dT overhang", polarity="3-prime",
                owner="adapter", text="one unpaired dT")
    return AdapterLigationGeometry.known(
        "splinted", sc,
        statement="Supplementary Table S2 gives a 17-base hairpin arm complementary to the first 17 bases of the 18-base rU-dT linker; E. coli DNA ligase seals their annealed nick and one linker dT remains unpaired.")


def promoter_scene():
    promoter = seg("T7 promoter", T7_PROMOTER, "t7")
    template = [seg("cell barcode", "B" * CELL_BARCODE_LEN, "cbc", placeholder=True,
                    feature=CELL_FEATURE),
                seg("UMI", "U" * UMI_LEN, "umi", placeholder=True, feature=UMI_FEATURE),
                seg("mRNA template", "X" * 34, placeholder=True)]
    return direct_ssrna_ivt(template, promoter)


def arna_scene():
    sc = Scene(); sc.strand("aRNA", [
        seg("aRNA transcript", "X" * 34, placeholder=True),
        seg("copied UMI", "U" * UMI_LEN, "umi", placeholder=True, feature=UMI_FEATURE),
        seg("copied cell barcode", "B" * CELL_BARCODE_LEN, "cbc", placeholder=True,
            feature=CELL_FEATURE)], label="antisense RNA copies")
    sc.note("aRNA", "hundreds of linear-amplification products from one original ssRNA")
    return sc


def random_rt_scene():
    arna = [seg("aRNA body", "X" * 34, placeholder=True)]
    primer = [seg("Read 1 handle", RANDOM_RT_HANDLE, "r1"),
              seg("random N6", "N" * 6, placeholder=True)]
    sc = Scene(); sc.strand("aRNA", arna, label="amplified antisense RNA")
    sc.anneal("random RT primer", primer, to="aRNA", pair=("random N6", "aRNA body"),
              shift=len(arna[0]) - 6, label="random RT primer",
              unpaired=("Read 1 handle", "random N6", "aRNA body"))
    sc.arrow("random RT primer", "SuperScript IV")
    sc.labels("random RT primer")
    return sc


FINAL_LIBRARY = Construct([
    seg("P5", il.P5, "p5"), seg("Read 1 arm", il.TRUSEQ_READ1, "r1"),
    seg("random-primer N6", "N" * 6, placeholder=True, length_bp=6),
    seg("aRNA-derived cDNA", "X" * 38, placeholder=True),
    seg("cell barcode", "B" * CELL_BARCODE_LEN, "cbc", placeholder=True,
        feature=CELL_FEATURE),
    seg("UMI", "U" * UMI_LEN, "umi", placeholder=True, feature=UMI_FEATURE),
    seg("Read 2 arm completion base", "A", "r2"),
    seg("Index 1 / Read 2 arm", il.INDEX1_PRIMER, "r2"),
    seg("i7 reverse complement", "I" * 8, "cbc", placeholder=True,
        feature=feature("last_sample_index", "sample_index", "unknown")),
    seg("P7 reverse complement", il.P7_RC, "p7"),
], name="LAST-seq library")

SEQ_PRIMERS = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["R2"])
_errors = sp.verify(FINAL_LIBRARY, SEQ_PRIMERS,
                    required_roles=tuple(p.role for p in SEQ_PRIMERS))
if _errors:
    raise ValueError("LAST-seq final library: " + "; ".join(_errors))


def final_rows():
    sc = Scene.duplex(list(FINAL_LIBRARY), label="LAST-seq")
    sc.labels("top")
    return [*sc.rows(), *annotation_rows(FINAL_LIBRARY)]


def workflow():
    wf = Workflow(MolecularState("Separate hairpin and rU-dT primer modules",
                                 (Row(chunks=[("hairpin module + 5′-phosphorylated rU-dT module", None, False)]),)))
    wf.ligate_adapter("Anneal and nick-ligate the LAST-seq capture primer",
                      primer_ligation_geometry(), capture_scene().rows(),
                      name="assembled hairpin/rU-dT primer",
                      note="The gel-purified product is then annealed to mRNA poly(A).")
    wf.react("RNase H nicking and Klenow exo− extension", promoter_scene().rows(),
             note="A short rA/dT hybrid is nicked and extended so the original RNA gains a duplex T7 promoter without reverse transcription.")
    wf.react("T7 in-vitro transcription directly from ssRNA", arna_scene().rows(),
             note="The original single-stranded RNA is the transcription template; amplification is linear.")
    wf.react("Random-prime the pooled aRNA and reverse-transcribe", random_rt_scene().rows(),
             note="Reverse transcription occurs only after direct RNA amplification.")
    wf.react("P5/P7 PCR and library enrichment", final_rows(),
             note="PCR completes the Illumina library; the inline cell barcode and UMI originate in the capture hairpin.")
    return wf


TITLE = "LAST-seq — direct linear amplification of single-stranded RNA"
NOTES = "01_last-seq.html"
SOURCE = ('Defining source: <a href="https://doi.org/10.1186/s13059-023-03025-5">'
          'Lyu and Chen, <i>Genome Biology</i> (2023)</a>.')
SUMMARY = ("A ligated hairpin/rU-dT primer installs a short duplex T7 promoter at an mRNA’s "
           "3′ end. T7 polymerase then amplifies the original RNA directly, before any reverse transcription.")
FINAL_CAPTION = "Representative PCR-completed LAST-seq molecule. Read 2 reaches the hairpin-derived UMI and cell barcode; Read 1 enters the randomly primed aRNA-derived insert."
SEQUENCING_INTRO = "NextSeq 550: Read 1, 50 cycles; Read 2, 25 cycles. The paper reports 6-nt cellular barcodes and 8-nt UMIs."
READ_LENGTHS = {"Read 1": 50, "Read 2": 25}
SHOW_READ_CONTENT = True


def oligos():
    yield oligo("rU-dT module", [seg("linker", DTRU_LINKER, "tso"),
                                  seg("dT7", DTRU_DNA), seg("rU19", DTRU_RNA)],
                mods="5′ phosphate; 3′ ddC")
    yield oligo("RandomRT primer6", [seg("Read 1 handle", RANDOM_RT_HANDLE, "r1"),
                                      seg("N6", "N" * 6, placeholder=True)])


def protocol_graph():
    graph = ProtocolGraph(TITLE, authority=SOURCE)
    rna = graph.input("rna", "polyadenylated cellular RNA", MaterialKind.BIOLOGICAL,
                      state=molecule("mRNA", (role("transcript body"), role("poly(A) tail"))),
                      primary=True)
    hairpin = graph.input("hairpin", "barcode/UMI/T7-promoter hairpin module",
                          MaterialKind.MOLECULAR)
    rudt = graph.input("rudt", "5-prime-phosphorylated rU-dT module",
                       MaterialKind.ADAPTER)
    ligase = graph.input("primer_ligase", "E. coli DNA ligase", MaterialKind.CATALYST)
    assembled = graph.product("capture_primer", "assembled hairpin/rU-dT capture primer",
                              MaterialKind.MOLECULAR, state=capture_scene())
    graph.react("assemble_primer", "Nick-ligate the LAST-seq capture primer",
                StepKind.LIGATION,
                ((hairpin, InputRole.TEMPLATE), (rudt, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (assembled,),
                geometry=primer_ligation_geometry())
    hybrid = graph.product("captured_rna", "capture primer annealed to mRNA poly(A)",
                           MaterialKind.MOLECULAR, state=capture_scene())
    graph.react("capture", "Anneal the assembled primer to mRNA poly(A)",
                StepKind.HYBRIDIZATION,
                ((rna, InputRole.TEMPLATE), (assembled, InputRole.STRUCTURAL)), (hybrid,))
    rnase = graph.input("rnase_h", "RNase H", MaterialKind.CATALYST)
    klenow = graph.input("klenow", "Klenow exo-minus polymerase", MaterialKind.CATALYST)
    promoter_rna = graph.product("promoter_rna", "ssRNA carrying a duplex T7 promoter",
                                 MaterialKind.MOLECULAR, state=promoter_scene())
    graph.react("install_promoter", "Nick the rA/dT tract and extend the duplex promoter",
                StepKind.EXTENSION,
                ((hybrid, InputRole.TEMPLATE), (rnase, InputRole.CATALYST),
                 (klenow, InputRole.CATALYST)), (promoter_rna,))
    t7 = graph.input("t7_polymerase", "T7 RNA polymerase", MaterialKind.CATALYST)
    arna = graph.product("arna", "linearly amplified antisense RNA copies",
                         MaterialKind.MOLECULAR, state=arna_scene())
    graph.react("transcribe", "Transcribe antisense RNA directly from the original ssRNA",
                StepKind.TRANSCRIPTION,
                ((promoter_rna, InputRole.TEMPLATE), (t7, InputRole.CATALYST)), (arna,))
    random_primer = graph.input("random_rt_primer", "RandomRT primer6", MaterialKind.OLIGO)
    rt = graph.input("reverse_transcriptase", "SuperScript IV", MaterialKind.CATALYST)
    cdna = graph.product("cdna", "random-primed cDNA from pooled aRNA",
                         MaterialKind.MOLECULAR, state=random_rt_scene())
    graph.react("reverse_transcribe", "Random-prime and reverse-transcribe pooled aRNA",
                StepKind.REVERSE_TRANSCRIPTION,
                ((arna, InputRole.TEMPLATE), (random_primer, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (cdna,))
    forward = graph.input("pcr_forward", "P5-side library PCR primer", MaterialKind.OLIGO)
    reverse = graph.input("pcr_reverse", "P7-side library PCR primer", MaterialKind.OLIGO)
    polymerase = graph.input("pcr_polymerase", "library PCR polymerase", MaterialKind.CATALYST)
    library = graph.product("library", FINAL_LIBRARY.name, MaterialKind.LIBRARY,
                            state=FINAL_LIBRARY)
    graph.react("library_pcr", "Complete P5/P7 arms by PCR", StepKind.PCR,
                ((cdna, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (library,), geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY, reason="The pooled random-RT precursor is heterogeneous; the completed construct has checked primer sites."))
    graph.sequence(library, SequencingPlan(
        library=FINAL_LIBRARY, primers=tuple(SEQ_PRIMERS),
        read_lengths=READ_LENGTHS))
    graph.declare_reviewed_scope([
        "rna", "hairpin", "rudt", "primer_ligase", "assemble_primer", "capture",
        "rnase_h", "klenow", "install_promoter", "t7_polymerase", "transcribe",
        "random_rt_primer", "reverse_transcriptase", "reverse_transcribe",
        "pcr_forward", "pcr_reverse", "pcr_polymerase", "library_pcr", "library"], SOURCE)
    return graph.require_publishable()
