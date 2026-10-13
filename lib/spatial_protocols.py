"""Shared molecular models for spatial capture and sample hashing protocols."""
from __future__ import annotations

from dataclasses import dataclass
from types import SimpleNamespace

import illumina as il
import nextera as nx
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry
from batch_ngs import joined_scene, seg, truseq_library
from chemdraw import (Construct, MolecularState, Scene, Segment, Workflow, duplex_rows,
                      feature, oligo, revcomp)
from multimodal_spatial import HybridizationLigationRound
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from protocol_model import chain, molecule, role, role_library


def ident(name, identifier, kind, *, length=None, tag="cbc"):
    return role(name, tag=tag, length_bp=length,
                molecular_feature=feature(identifier, kind, "unknown"))


SPATIAL_A = ident("spatial barcode A", "spatial_a", "spatial_barcode")
SPATIAL_B = ident("spatial barcode B", "spatial_b", "spatial_barcode")
SPATIAL_BEAD = ident("decoded bead barcode", "spatial_bead", "spatial_barcode")
UMI = ident("UMI", "umi", "umi", tag="umi")
SAMPLE = ident("sample index", "sample_index", "sample_index", length=8)
HASH = ident("nuclear hash", "sample_hash", "sample_index")
CELL_A = ident("first-round cell barcode", "cell_round1", "cell_barcode")
CELL_B = ident("second-round cell barcode", "cell_round2", "cell_barcode")
SPACE_HASH = ident("printed spatial hash", "sci_space.grid", "spatial_barcode")
SECTION_HASH = ident("section hash", "sci_space.section", "sample_index")


# Li et al., Nat Protoc 2025 Supplementary Table 1.  The first member of each
# 100-oligo pool is retained as a representative; the molecular feature records
# the full published plate as the identifier source.
SPATIAL_RT_HANDLE = "CATCGGCGTACGACT"
SPATIAL_RT_UMI = "N" * 10
SPATIAL_RT_POLYT = "T" * 15 + "VN"
SPATIAL_LINKER_1 = "AGTCGTACGCCGATGCGAAACATCGGCCAC"
SPATIAL_LINKER_2 = "CGAATGCTCTGGCCTCTCAAGCACGTGGAT"
SPATIAL_A_PREFIX = "AGGCCAGAGCATTCG"
SPATIAL_A1 = "AACGTGAT"
SPATIAL_A_SUFFIX = "GTGGCCGATGTTTCG"
SPATIAL_B_PREFIX = "CAAGCGTTGGCTTCTCGCATCT"
SPATIAL_B1 = "AACGTGAT"
SPATIAL_B_SUFFIX = "ATCCACGTGCTTGAG"
SPATIAL_TN5_A = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
SPATIAL_TN5_B = "CATCGGCGTACGACTAGATGTGTATAAGAGACAG"
SPATIAL_TN5_REV = "CTGTCTCTTATACACATCT"
SPATIAL_N501 = "AATGATACGGCGACCACCGAGATCTACACTAGATCGCTCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
SPATIAL_N701 = "CAAGCAGAAGACGGCATACGAGATTCGCCTTAGTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGCAAGCGTTGGCTTCTCGCATCT"

_SPATIAL_A_FEATURE = feature(
    "spatial_a", "spatial_barcode", "whitelist",
    whitelist="100-member barcode A plate (Li et al. Supplementary Table 1)")
_SPATIAL_B_FEATURE = feature(
    "spatial_b", "spatial_barcode", "whitelist",
    whitelist="100-member barcode B plate (Li et al. Supplementary Table 1)")

SPATIAL_A_ROUND = HybridizationLigationRound(
    "barcode A", SPATIAL_RT_HANDLE, SPATIAL_A_PREFIX, SPATIAL_A_SUFFIX,
    SPATIAL_LINKER_1, _SPATIAL_A_FEATURE, phosphorylated_5=True,
    barcode_length=8)
SPATIAL_B_ROUND = HybridizationLigationRound(
    "barcode B", SPATIAL_A_PREFIX, SPATIAL_B_PREFIX, SPATIAL_B_SUFFIX,
    SPATIAL_LINKER_2, _SPATIAL_B_FEATURE, barcode_length=8)


def _spatial_branch_scene(kind: str, stage: int) -> Scene:
    """Representative source-exact strand after RT/tagmentation and A/B ligation."""
    if kind == "RNA":
        core = [seg("RT handle", SPATIAL_RT_HANDLE, "r1"),
                seg("10-nt UMI", "U" * 10, "umi", placeholder=True,
                    feature=feature("spatial_rna_umi", "umi", "random")),
                seg("poly(dT)VN", SPATIAL_RT_POLYT, placeholder=True),
                seg("transcript cDNA", "X" * 28, placeholder=True),
                seg("TSO-derived end", "CCCATTC ACTCTGCGTTGATACCACTGCTT".replace(" ", ""))]
    else:
        core = [seg("Tn5ME-B RT handle", SPATIAL_RT_HANDLE, "r1"),
                seg("Tn5ME-B mosaic end", "AGATGTGTATAAGAGACAG", "me"),
                seg(f"{kind} genomic DNA", "X" * 30, placeholder=True),
                seg("Tn5ME-A mosaic end reverse complement",
                    "CTGTCTCTTATACACATCT", "me"),
                seg("Tn5ME-A outer arm reverse complement", "TCGTCGGCAGCGTCAG")]
    parts = core
    junctions = []
    if stage >= 1:
        a = [seg("barcode A outgoing arm", SPATIAL_A_PREFIX),
             seg("barcode A", "A" * 8, "cbc", placeholder=True,
                 feature=_SPATIAL_A_FEATURE),
             seg("barcode A ligating arm", SPATIAL_A_SUFFIX)]
        parts = [*a, *parts]
        junctions.append(("barcode A ligating arm", core[0].name, "ligation A"))
    if stage >= 2:
        b = [seg("barcode B PCR arm", SPATIAL_B_PREFIX),
             seg("barcode B", "B" * 8, "cbc", placeholder=True,
                 feature=_SPATIAL_B_FEATURE),
             seg("barcode B ligating arm", SPATIAL_B_SUFFIX)]
        parts = [*b, *parts]
        junctions.insert(0, ("barcode B ligating arm", "barcode A outgoing arm",
                             "ligation B"))
    return joined_scene(parts, tuple(junctions), label=f"{kind} branch", duplex=False)


def _spatial_rows(stage: int, dna_kind: str) -> list:
    return [*_spatial_branch_scene(dna_kind, stage).rows(),
            *_spatial_branch_scene("RNA", stage).rows()]


def _spatial_oligos():
    return (
        oligo("RT primer (5′-phosphate; one dT is internal biotin)", [
            seg("RT handle", SPATIAL_RT_HANDLE),
            seg("10-nt UMI", SPATIAL_RT_UMI, "umi", placeholder=True),
            seg("biotin-dT / poly(dT)VN", SPATIAL_RT_POLYT, placeholder=True)]),
        oligo("Tn5ME-B", [seg("RT-compatible handle", SPATIAL_RT_HANDLE),
                           seg("mosaic end", "AGATGTGTATAAGAGACAG")]),
        oligo("Barcode A-1 (5′-phosphate)", [seg("A outgoing arm", SPATIAL_A_PREFIX),
              seg("A-1", SPATIAL_A1, "cbc", feature=_SPATIAL_A_FEATURE),
              seg("A ligating arm", SPATIAL_A_SUFFIX)]),
        oligo("Ligation linker 1", [seg("splint", SPATIAL_LINKER_1)]),
        oligo("Barcode B-1", [seg("PCR arm", SPATIAL_B_PREFIX),
              seg("B-1", SPATIAL_B1, "cbc", feature=_SPATIAL_B_FEATURE),
              seg("B ligating arm", SPATIAL_B_SUFFIX)]),
        oligo("Ligation linker 2", [seg("splint", SPATIAL_LINKER_2)]),
    )


def _spatial_dna_libraries(dna_kind: str):
    """Exact representative DNA branch before and after N501/N701 PCR.

    The source draws the barcode-bearing end P7-to-insert.  Constructs use the
    canonical P5-to-P7 top-strand orientation, hence the reverse complements of
    the B/A/RT chain on the right side of the insert.
    """
    i5 = feature("spatial_i5", "sample_index", "whitelist",
                 whitelist="N501 index from Li et al. Supplementary Table 1")
    i7 = feature("spatial_i7", "sample_index", "whitelist",
                 whitelist="N701 index from Li et al. Supplementary Table 1")
    insert = [
        seg("Tn5ME-A outer arm", nx.S5, "s5"),
        seg("Tn5ME-A mosaic end", nx.ME, "me"),
        seg(f"{dna_kind} genomic DNA", "X" * 30, placeholder=True),
        seg("Tn5ME-B mosaic end reverse complement", nx.ME_RC, "me"),
        seg("shared RT/Tn5ME-B handle reverse complement", revcomp(SPATIAL_RT_HANDLE)),
        seg("barcode A ligating arm reverse complement", revcomp(SPATIAL_A_SUFFIX)),
        seg("barcode A reverse complement", "A" * 8, "cbc", placeholder=True,
            feature=_SPATIAL_A_FEATURE),
        seg("barcode A outgoing arm reverse complement", revcomp(SPATIAL_A_PREFIX)),
        seg("barcode B ligating arm reverse complement", revcomp(SPATIAL_B_SUFFIX)),
        seg("barcode B reverse complement", "B" * 8, "cbc", placeholder=True,
            feature=_SPATIAL_B_FEATURE),
        seg("barcode B PCR arm reverse complement", revcomp(SPATIAL_B_PREFIX)),
    ]
    template = Construct(insert, name=f"{dna_kind} branch before indexed PCR")
    final = Construct([
        seg("P5", il.P5, "p5"),
        seg("i5 (N501)", "TAGATCGC", "cbc", feature=i5),
        *insert,
        seg("Nextera Read 2 mosaic end reverse complement", nx.ME_RC, "me"),
        seg("S7 reverse complement", nx.S7_RC, "s7"),
        seg("i7 (N701) reverse complement", revcomp("TCGCCTTA"), "cbc", feature=i7),
        seg("P7 reverse complement", il.P7_RC, "p7"),
    ], name=f"representative {dna_kind} sequencing library")
    geometry = PCRGeometry.exact(
        template,
        PCRPrimer("N501", SPATIAL_N501, anneal_nt=len(nx.ADAPTOR_S5)),
        PCRPrimer("N701", SPATIAL_N701, anneal_nt=len(SPATIAL_B_PREFIX)),
        statement="The 2025 protocol oligo table gives the complete N501 and N701 primers; their 3′ arms land exactly on Tn5ME-A and the barcode-B PCR arm.")
    primers = sp.standard_primers(final)
    return template, final, geometry, primers


def _dbit(*, cuttag: bool):
    enzyme = "antibody-tethered pA–Tn5" if cuttag else "Tn5"
    dna_kind = "histone-mark DNA" if cuttag else "accessible DNA"
    start = Scene(); start.strand("targets", [role(dna_kind), role("poly(A) RNA")],
                                 label="fixed permeabilized tissue")
    wf = Workflow(MolecularState("Fixed tissue targets", tuple(start.rows())))
    wf.react(f"Tag {dna_kind} with {enzyme} and reverse-transcribe RNA",
             _spatial_rows(0, dna_kind), name="shared ligation handle on both branches",
             note="Tn5ME-B and the 5′-phosphorylated RT primer carry the same exact CATCGGCGTACGACT handle.")
    g1 = AdapterLigationGeometry.known(
        "splinted", SPATIAL_A_ROUND.adapter_scene(),
        statement="Linker 1 is exactly the reverse complement of the shared RT/Tn5ME-B handle followed by the barcode-A ligating arm.")
    wf.ligate_adapter("Ligate barcode A in the first microfluidic channel set", g1,
                      _spatial_rows(1, dna_kind), name="A-barcoded DNA and RNA branches",
                      note="The A-oligo 3′-OH is sealed to the shared 5′-phosphate; A-1 is shown from the 100-member plate.")
    g2 = AdapterLigationGeometry.known(
        "splinted", SPATIAL_B_ROUND.adapter_scene(),
        statement="Linker 2 is exactly the reverse complement of the barcode-A outgoing arm followed by the barcode-B ligating arm.")
    wf.ligate_adapter("Rotate the device and ligate barcode B", g2,
                      _spatial_rows(2, dna_kind), name="two-dimensional pixel address",
                      note="The B-oligo 3′-OH is sealed to the barcode-A 5′-phosphate; the ordered A–B pair identifies one channel intersection.")
    wf.react("Recover the biotin-selected RNA branch and amplify both libraries",
             _spatial_rows(2, dna_kind), name="separate matched DNA and RNA libraries",
             note="Both physical libraries retain the same exact B–A spatial chain; only the RNA branch carries the RT UMI and poly(dT).")
    _template, final, geometry, _primers = _spatial_dna_libraries(dna_kind)
    wf.pcr("Add the indexed Illumina shell to the DNA branch", geometry,
           duplex_rows(final, label="DNA library", unpaired=(f"{dna_kind} genomic DNA'",)),
           name="indexed DNA-branch library",
           note="N501 supplies P5/i5 at the Tn5ME-A end; N701 supplies P7/i7 at the spatial-barcode end. The independently recovered RNA branch is Nextera-tagmented after pre-amplification and is not collapsed into this fixed cut site.")
    return wf


def _slide_seq():
    return chain(molecule("10-µm bead puck", (SPATIAL_BEAD, UMI, role("poly(dT) capture oligo"))), (
        ("Decode each bead position by SOLiD sequencing-by-ligation",
         molecule("mapped bead puck", (SPATIAL_BEAD, role("known x–y coordinate"))),
         "The decoded bead sequence is a spatial address before tissue is applied."),
        ("Transfer a fresh-frozen tissue section and release RNA",
         molecule("RNA captured on bead", (SPATIAL_BEAD, UMI, role("poly(dT)"),
                                             role("poly(A) RNA"))),
         "RNA inherits position by hybridizing to an oligo on a mapped bead."),
        ("Reverse-transcribe and prepare a 3′ barcoded RNA-seq library",
         molecule("Slide-seq cDNA", (SPATIAL_BEAD, UMI, role("3-prime transcript cDNA"))),
         "The original method follows the Drop-seq-style 3′ capture logic."),
        ("Amplify, tagment and add the sample index",
         molecule("Slide-seq sequencing library", (SPATIAL_BEAD, UMI,
             role("3-prime transcript insert"), SAMPLE)),
         "Sequencing joins each transcript molecule to the bead position decoded on the puck."),
    ))


def _sciplex_atac():
    return chain(molecule("permeabilized nuclei separated by sample", (role("accessible chromatin"),)), (
        ("Let nuclei absorb a sample-specific unmodified ssDNA hash oligo",
         molecule("hashed nucleus", (HASH, role("accessible chromatin"))),
         "The hash records specimen or treatment; it is not a cell barcode."),
        ("Extend the polyadenylated hash and tagment accessible chromatin in indexed wells",
         molecule("first-round indexed nucleus", (HASH, CELL_A, role("tagmented accessible DNA"))),
         "Hash extension makes the label recoverable alongside the ATAC assay."),
        ("Pool, redistribute and add the second combinatorial index",
         molecule("two-level cell identity", (HASH, CELL_A, CELL_B,
                                                role("accessible-DNA fragment"))),
         "The pair of well indices identifies a nucleus; the independent hash identifies its sample."),
        ("Amplify ATAC fragments and hash reads",
         molecule("sciPlex-ATAC readouts", (HASH, CELL_A, CELL_B,
                                             role("ATAC insert"), SAMPLE)),
         "Two-level experiments sequence hash and ATAC products separately; later three-level variants can intermix them."),
    ))


def _sci_space():
    return chain(molecule("fresh-frozen tissue section", (
        role("poly(A) RNA in intact nuclei"), role("anatomical position"))), (
        ("Permeabilize the section and add one polyadenylated section hash",
         molecule("section-labelled nuclei in tissue", (SECTION_HASH,
             role("poly(A) RNA in intact nuclei"))),
         "The section hash identifies the slide after sections are pooled; it is distinct from the printed coordinate hashes."),
        ("Juxtapose tissue with the agarose slide carrying a gridded hash pattern",
         molecule("spatially hashed nuclei", (SPACE_HASH, SECTION_HASH,
             role("poly(A) RNA in intact nuclei"))),
         "Single-stranded poly(A)-tailed hash oligos diffuse from each printed spot into nearby permeabilized nuclei. Fluorescent fiducials register the grid to the tissue image."),
        ("Recover nuclei, fix them and perform barcoded sci-RNA-seq reverse transcription",
         molecule("first-round indexed nuclear molecules", (
             ident("RT-well cell barcode", "sci_space.rt", "cell_barcode"),
             UMI, role("transcript cDNA or copied hash"))),
         "The same oligo(dT) reaction copies endogenous mRNA, the grid hash and the section hash, associating all three with a nucleus-specific RT-well barcode."),
        ("Pool, redistribute, second-strand synthesize, tagment and PCR-index",
         molecule("sci-Space sequencing molecules", (
             ident("RT-well cell barcode", "sci_space.rt", "cell_barcode"),
             UMI, SPACE_HASH, role("transcript or hash payload"), SAMPLE)),
         "The downstream combinatorial sci-RNA-seq indices identify one nucleus. Hash molecules and transcript molecules are separate physical reads joined computationally by that shared nucleus identity."),
    ))


def sci_space_library():
    """Representative spatial-hash molecule in the sci-RNA-seq hybrid shell."""
    def s(name, top, tag=None, **kw):
        return Segment(name, top, tag, **kw)
    grid = feature("sci_space.grid", "spatial_barcode", "whitelist",
                   whitelist="sci-Space printed coordinate-hash panel",
                   note="printed coordinate-hash identity registered by fluorescent fiducials")
    cell = feature("sci_space.rt", "cell_barcode", "combinatorial",
                   group="cell_barcode", part="RT round")
    umi = feature("sci_space.umi", "umi", "random")
    i5 = feature("sci_space.i5", "sample_index", "unknown")
    i7 = feature("sci_space.i7", "sample_index", "unknown")
    lib = Construct([
        s("P5", il.P5, "p5"),
        s("i5", "I" * 10, "cbc", placeholder=True, feature=i5),
        s("TruSeq Read 1", il.TRUSEQ_READ1, "r1"),
        role("UMI", tag="umi", molecular_feature=umi),
        role("RT-well cell barcode", tag="cbc", molecular_feature=cell),
        role("copied printed spatial hash", tag="cbc", molecular_feature=grid),
        s("mosaic end reverse complement", nx.ME_RC, "me"),
        s("S7 reverse complement", nx.S7_RC, "s7"),
        s("i7 reverse complement", "J" * 10, "cbc", placeholder=True, feature=i7),
        s("P7 reverse complement", il.P7_RC, "p7"),
    ], name="representative sci-Space spatial-hash library molecule")
    primers = (sp.TRUSEQ["R1"], sp.NEXTERA["I1"],
               sp.TRUSEQ["I2"], sp.NEXTERA["R2"])
    problems = sp.verify(lib, primers, required_roles=tuple(p.role for p in primers))
    if problems:
        raise ValueError("sci-Space library: " + "; ".join(problems))
    return lib, primers


HDST_PRECURSOR_WRITTEN = (
    "/AmC6/UUUUUGACTCGTAATACGACTCACTATAGGGACACGACGCTCTTCCGATCT-"
    "Spatial_barcode_Pool1-Bridge1"
)


def _hdst():
    return chain(molecule("2-micrometre amine-functionalized silica bead", (
        role("surface amines"),)), (
        ("Couple the uracil-bearing T7 / Illumina precursor with barcode pool 1",
         molecule("first-round bead oligo", (role("bead linkage and U5 release site"),
             role("T7 promoter and Illumina handle"),
             ident("15-nt barcode pool 1", "hdst.bc1", "spatial_barcode", length=15),
             role("bridge 1"))),
         "The published precursor skeleton is exact; the individual Illumina-supplied barcode and decoder sequences are proprietary."),
        ("Split, ligate barcode pools 2 and 3, then append UMI and poly(dT) capture",
         molecule("completed bead capture oligo", (
             ident("barcode pool 1", "hdst.bc1", "spatial_barcode", length=15),
             role("bridge 1"), ident("barcode pool 2", "hdst.bc2", "spatial_barcode", length=15),
             role("bridge 2"), ident("barcode pool 3", "hdst.bc3", "spatial_barcode", length=15),
             ident("six-base UMI", "hdst.umi", "umi", length=6, tag="umi"),
             role("poly(dT) capture"))),
         "The three split-and-pool rounds provide 65 × 211 × 211 combinations. Every oligo on one bead shares its three-part spatial identity but has an independent UMI."),
        ("Load beads randomly into 2-micrometre wells and decode each position",
         molecule("decoded HDST array", (role("known x-y well"),
             role("one three-part spatial barcode"), role("many capture oligos"))),
         "Fourteen cycles of FAM, Cy3 or dark decoder hybridization map each randomly placed bead to an image coordinate before tissue is applied."),
        ("Place, image and permeabilize tissue; reverse-transcribe captured poly(A) RNA",
         molecule("bead-bound spatial cDNA", (
             role("bead"), role("three-part spatial barcode"),
             role("six-base UMI"), role("poly(dT)"), role("antisense transcript cDNA"))),
         "Superscript III extends the bead oligo overnight. Barcode, UMI and transcript information are now covalently linked in one first strand."),
        ("Digest tissue, release at the five-uracil site with USER and prepare RNA-seq library",
         molecule("paired-end HDST library", (
             role("Read 1 arm"), role("three barcode blocks and bridges"),
             role("six-base UMI"), role("transcript insert"), role("Read 2 and indexed arm"))),
         "The released cDNA follows the published Spatial Transcriptomics CEL-seq-like amplification and library workflow. Read 1 reports position and UMI; Read 2 reports transcript."),
    ))


def hdst_library():
    bc1 = feature("hdst.bc1", "spatial_barcode", "whitelist", group="hdst.spatial",
                  part="pool 1", whitelist="proprietary HDST pool 1")
    bc2 = feature("hdst.bc2", "spatial_barcode", "whitelist", group="hdst.spatial",
                  part="pool 2", whitelist="proprietary HDST pool 2")
    bc3 = feature("hdst.bc3", "spatial_barcode", "whitelist", group="hdst.spatial",
                  part="pool 3", whitelist="proprietary HDST pool 3")
    umi = feature("hdst.umi", "umi", "random")
    insert = [
        seg("spatial barcode pool 1", "B" * 15, "cbc", placeholder=True, feature=bc1),
        seg("bridge 1", "X" * 15, placeholder=True, inferred=True, length_bp=15),
        seg("spatial barcode pool 2", "C" * 15, "cbc", placeholder=True, feature=bc2),
        seg("bridge 2", "X" * 15, placeholder=True, inferred=True, length_bp=15),
        seg("spatial barcode pool 3", "D" * 15, "cbc", placeholder=True, feature=bc3),
        seg("terminal spacer", "X", placeholder=True, inferred=True, length_bp=1),
        seg("UMI", "U" * 6, "umi", placeholder=True, feature=umi),
        role("transcript cDNA"),
    ]
    lib, all_primers = truseq_library(insert, "representative HDST sequencing library",
                                      dual_index=False, inferred_adapters=True, i7_nt=8)
    primers = tuple(p for p in all_primers if p.role in ("Read 1", "Index 1 (i7)", "Read 2"))
    return lib, primers


@dataclass(frozen=True)
class Definition:
    title: str; doi: str; citation: str; summary: str; endpoint: str; notes: str


DEFS = {
    "hdst": Definition("HDST — High-Definition Spatial Transcriptomics", "10.1038/s41592-019-0548-y",
        "Vickovic et al., Nature Methods (2019)", "A decoded array of 2-micrometre wells holds randomly placed beads bearing three-part spatial barcodes, UMIs and poly(dT); captured transcripts inherit the bead's imaged coordinate.",
        "Read 1 encodes the three barcode blocks and six-base UMI; Read 2 is transcript sequence. The representative construct marks unavailable proprietary bridge and barcode-panel bases.", "01_hdst.html"),
    "sci-space": Definition("sci-Space", "10.1126/science.abb9536",
        "Srivatsan et al., Science (2021)", "Polyadenylated DNA hashes diffuse from a registered printed grid into permeabilized tissue nuclei; sci-RNA-seq then associates each grid identity with a single-nucleus transcriptome.",
        "Separate transcript, printed-grid-hash and section-hash molecules share the same combinatorial nucleus indices; the representative final construct below is a spatial-hash molecule.", "01_sci_space.html"),
    "spatial-atac-rna": Definition("spatial-ATAC–RNA-seq", "10.1038/s41586-023-05795-1",
        "Zhang et al., Nature (2023)", "Tn5 tags accessible DNA and in-situ poly(dT) RT captures RNA before orthogonal microfluidic barcode ligations assign both branches the same tissue pixel.",
        "The exact representative panel is the ATAC branch: Read 1 enters accessible DNA, while Read 2 enters the B–A spatial chain from the N701 end. The separately recovered RNA branch retains the same spatial address plus a UMI and is Nextera-tagmented after pre-amplification.", "01_spatial_atac_rna_seq.html"),
    "spatial-cuttag-rna": Definition("spatial-CUT&Tag–RNA-seq", "10.1038/s41586-023-05795-1",
        "Zhang et al., Nature (2023)", "Antibody-tethered pA–Tn5 tags histone-marked DNA and in-situ poly(dT) RT captures RNA before both branches receive orthogonal spatial barcodes.",
        "The exact representative panel is the CUT&Tag branch: Read 1 enters histone-mark DNA, while Read 2 enters the B–A spatial chain from the N701 end. The separately recovered RNA branch retains the same spatial address plus a UMI and is Nextera-tagmented after pre-amplification.", "01_spatial_cuttag_rna_seq.html"),
    "slide-seq": Definition("Slide-seq (v1)", "10.1126/science.aaw1219",
        "Rodriques et al., Science (2019)", "A pre-decoded monolayer of 10-µm DNA-barcoded beads captures RNA released from a tissue section for spatially registered 3′ RNA sequencing.",
        "The barcode read identifies a bead's decoded position and its UMI; the paired insert read identifies the captured transcript.", "01_slide_seq.html"),
    "sciplex-atac": Definition("sciPlex-ATAC-seq", "10.1101/2023.03.05.531201",
        "Booth et al., bioRxiv (2023; Genome Research 2023)", "Permeabilized nuclei absorb sample-specific ssDNA hashes, then two-level combinatorial ATAC indexing assigns a separate cell identity.",
        "Hash reads identify sample origin; the ordered combinatorial indices identify the nucleus; paired genomic reads report accessibility.", "01_sciplex_atac_seq.html"),
}


def _finish_pcr(graph, prefix, template, final, primers, endpoint, reviewed, *, label,
                geometry=None):
    forward = graph.input(prefix + "_forward", label + " forward PCR primer",
                          MaterialKind.OLIGO)
    reverse = graph.input(prefix + "_reverse", label + " reverse PCR primer",
                          MaterialKind.OLIGO)
    polymerase = graph.input(prefix + "_polymerase", label + " PCR polymerase",
                             MaterialKind.CATALYST)
    library = graph.product(prefix + "_library", label + " sequencing library",
                            MaterialKind.LIBRARY, state=final)
    graph.react(prefix + "_pcr", "Amplify and index the " + label + " library",
                StepKind.PCR,
                ((template, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (library,), geometry=(geometry or PCRGeometry.unavailable(
                    final, reason="The barcode panel or kit revision does not define one universal amplicon geometry.")))
    plan = (SequencingPlan(library=final, primers=tuple(primers)) if primers
            else SequencingPlan(limitation=endpoint))
    graph.sequence(library, plan)
    reviewed.extend((prefix + "_forward", prefix + "_reverse", prefix + "_polymerase",
                     prefix + "_pcr", prefix + "_library"))
    return library


def _hdst_graph(d, final, primers):
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    bead = graph.input("bead", "amine-functionalized silica bead bearing the T7/Illumina precursor",
                       MaterialKind.BIOLOGICAL, state=_hdst().initial, primary=True)
    pools = graph.input("barcode_pools", "three split-and-pool barcode/bridge oligo pools",
                        MaterialKind.ADAPTER)
    ligase = graph.input("ligase", "barcode-assembly DNA ligase", MaterialKind.CATALYST)
    capture_bead = graph.product("capture_bead", "bead bearing three-part spatial barcode, UMI and poly(dT)",
                                 MaterialKind.MOLECULAR,
                                 state=molecule("completed bead capture oligo", (
                                     role("three barcode blocks and bridges"), UMI,
                                     role("poly(dT) capture"))))
    graph.react("assemble_bead", "Ligate the split-and-pool barcode rounds and capture tail",
                StepKind.LIGATION,
                ((bead, InputRole.TEMPLATE), (pools, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (capture_bead,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The public source gives the precursor skeleton and component lengths but withholds the proprietary barcode and decoder sequences."))
    decoder = graph.input("decoder_oligos", "fluorescent decoder oligo series",
                          MaterialKind.OLIGO)
    decoded = graph.product("decoded_array", "position-decoded HDST bead array",
                            MaterialKind.MOLECULAR,
                            state=molecule("decoded HDST array", (
                                role("known x-y well"), role("three-part spatial barcode"))))
    graph.react("decode", "Decode bead positions by iterative hybridization",
                StepKind.HYBRIDIZATION,
                ((capture_bead, InputRole.TEMPLATE), (decoder, InputRole.STRUCTURAL)), (decoded,))
    tissue_rna = graph.input("tissue_rna", "poly(A) RNA released from the overlaid tissue section",
                             MaterialKind.MOLECULAR)
    captured = graph.product("captured_rna", "spatially barcoded bead-bound RNA",
                             MaterialKind.MOLECULAR,
                             state=molecule("RNA captured on bead", (SPATIAL_BEAD, UMI,
                                 role("poly(dT)"), role("poly(A) RNA"))))
    graph.react("capture_rna", "Hybridize released tissue RNA to bead poly(dT)",
                StepKind.HYBRIDIZATION,
                ((decoded, InputRole.TEMPLATE), (tissue_rna, InputRole.STRUCTURAL)), (captured,))
    rt_primer = graph.input("rt_primer", "bead-bound poly(dT) RT primer", MaterialKind.OLIGO)
    rt = graph.input("reverse_transcriptase", "Superscript III reverse transcriptase",
                     MaterialKind.CATALYST)
    cdna = graph.product("cdna", "bead-bound spatial cDNA", MaterialKind.MOLECULAR,
                         state=molecule("bead-bound spatial cDNA", (SPATIAL_BEAD, UMI,
                             role("antisense transcript cDNA"))))
    graph.react("reverse_transcribe", "Reverse-transcribe captured poly(A) RNA",
                StepKind.REVERSE_TRANSCRIPTION,
                ((captured, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (cdna,))
    user = graph.input("user_enzyme", "USER release enzyme", MaterialKind.CATALYST)
    released = graph.product("released_cdna", "cDNA released at the five-uracil site",
                             MaterialKind.MOLECULAR, state=cdna.state)
    graph.react("release", "Release cDNA at the five-uracil site", StepKind.DIGESTION,
                ((cdna, InputRole.TEMPLATE), (user, InputRole.CATALYST)), (released,))
    reviewed = ["bead", "barcode_pools", "ligase", "assemble_bead", "decoder_oligos",
                "decode", "tissue_rna", "capture_rna", "rt_primer", "reverse_transcriptase",
                "reverse_transcribe", "user_enzyme", "release"]
    _finish_pcr(graph, "library", released, final, primers, d.endpoint, reviewed,
                label="HDST")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _slide_graph(d, final):
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    puck = graph.input("puck", "DNA-barcoded 10-micrometre bead puck",
                       MaterialKind.BIOLOGICAL, state=_slide_seq().initial, primary=True)
    decoder = graph.input("decoder", "SOLiD sequencing-by-ligation decoder probes",
                          MaterialKind.OLIGO)
    mapped = graph.product("mapped_puck", "position-decoded bead puck",
                           MaterialKind.MOLECULAR,
                           state=molecule("mapped bead puck", (SPATIAL_BEAD,
                               role("known x-y coordinate"))))
    graph.react("decode", "Decode each bead position by sequencing-by-ligation",
                StepKind.HYBRIDIZATION,
                ((puck, InputRole.TEMPLATE), (decoder, InputRole.STRUCTURAL)), (mapped,))
    rna = graph.input("tissue_rna", "poly(A) RNA released from the tissue section",
                      MaterialKind.MOLECULAR)
    captured = graph.product("captured_rna", "RNA captured on mapped beads",
                             MaterialKind.MOLECULAR,
                             state=molecule("RNA captured on bead", (SPATIAL_BEAD, UMI,
                                 role("poly(A) RNA"))))
    graph.react("capture", "Capture tissue RNA on mapped bead oligo(dT)",
                StepKind.HYBRIDIZATION,
                ((mapped, InputRole.TEMPLATE), (rna, InputRole.STRUCTURAL)), (captured,))
    rt_primer = graph.input("rt_primer", "bead oligo(dT) RT primer", MaterialKind.OLIGO)
    rt = graph.input("reverse_transcriptase", "reverse transcriptase", MaterialKind.CATALYST)
    cdna = graph.product("cdna", "bead-barcoded 3-prime cDNA", MaterialKind.MOLECULAR,
                         state=molecule("Slide-seq cDNA", (SPATIAL_BEAD, UMI,
                             role("3-prime transcript cDNA"))))
    graph.react("reverse_transcribe", "Reverse-transcribe captured RNA",
                StepKind.REVERSE_TRANSCRIPTION,
                ((captured, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (cdna,))
    transposome = graph.input("transposome", "library-adapter-loaded transposome",
                              MaterialKind.ADAPTER)
    tn5 = graph.input("tn5", "Tn5 transposase", MaterialKind.CATALYST)
    tagged = graph.product("tagged_cdna", "tagmented Slide-seq cDNA", MaterialKind.MOLECULAR,
                           state=cdna.state)
    waste = graph.product("tn5_waste", "non-library tagmentation products", MaterialKind.WASTE)
    graph.react("tagment", "Tagment amplified bead-barcoded cDNA", StepKind.TAGMENTATION,
                ((cdna, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged, waste))
    graph.discard(waste, "does not retain an amplifiable spatial transcript insert")
    reviewed = ["puck", "decoder", "decode", "tissue_rna", "capture", "rt_primer",
                "reverse_transcriptase", "reverse_transcribe", "transposome", "tn5", "tagment"]
    _finish_pcr(graph, "library", tagged, final, (), d.endpoint, reviewed,
                label="Slide-seq")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _sci_space_graph(d, final, primers):
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    tissue = graph.input("tissue", "fresh-frozen tissue section containing intact nuclei",
                         MaterialKind.BIOLOGICAL, state=_sci_space().initial, primary=True)
    section_hash = graph.input("section_hash", "polyadenylated section-hash oligo",
                               MaterialKind.OLIGO)
    sectioned = graph.product("sectioned", "section-labelled tissue nuclei",
                              MaterialKind.MOLECULAR,
                              state=molecule("section-labelled nuclei", (SECTION_HASH,
                                  role("poly(A) RNA"))))
    graph.react("label_section", "Hybridize the section hash in permeabilized tissue",
                StepKind.HYBRIDIZATION,
                ((tissue, InputRole.TEMPLATE), (section_hash, InputRole.STRUCTURAL)),
                (sectioned,))
    grid_hash = graph.input("grid_hash", "registered printed spatial-hash grid",
                            MaterialKind.OLIGO)
    spatial = graph.product("spatial_nuclei", "nuclei carrying section and grid hashes",
                            MaterialKind.MOLECULAR,
                            state=molecule("spatially hashed nuclei", (SPACE_HASH,
                                SECTION_HASH, role("poly(A) RNA"))))
    graph.react("transfer_grid", "Transfer printed spatial hashes into nearby nuclei",
                StepKind.HYBRIDIZATION,
                ((sectioned, InputRole.TEMPLATE), (grid_hash, InputRole.STRUCTURAL)),
                (spatial,))
    rt_primer = graph.input("rt_primer", "barcoded sci-RNA-seq oligo(dT) RT primer",
                            MaterialKind.OLIGO)
    rt = graph.input("reverse_transcriptase", "reverse transcriptase", MaterialKind.CATALYST)
    first = graph.product("first_strand", "RT-well-indexed transcript and hash cDNA",
                          MaterialKind.MOLECULAR,
                          state=molecule("first-round indexed nuclear molecules", (
                              CELL_A, UMI, role("transcript cDNA or copied hash"))))
    graph.react("reverse_transcribe", "Copy endogenous RNA and both hash classes",
                StepKind.REVERSE_TRANSCRIPTION,
                ((spatial, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (first,))
    second_primer = graph.input("second_strand_primer", "second-strand synthesis primer",
                                MaterialKind.OLIGO)
    polymerase = graph.input("second_strand_polymerase", "second-strand polymerase",
                             MaterialKind.CATALYST)
    duplex = graph.product("duplex", "double-stranded indexed transcript/hash products",
                           MaterialKind.MOLECULAR, state=first.state)
    graph.react("second_strand", "Synthesize the second strand after pool redistribution",
                StepKind.EXTENSION,
                ((first, InputRole.TEMPLATE), (second_primer, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (duplex,))
    transposome = graph.input("transposome", "indexed library transposome", MaterialKind.ADAPTER)
    tn5 = graph.input("tn5", "Tn5 transposase", MaterialKind.CATALYST)
    tagged = graph.product("tagged", "tagmented sci-Space products", MaterialKind.MOLECULAR,
                           state=duplex.state)
    waste = graph.product("tn5_waste", "non-library tagmentation products", MaterialKind.WASTE)
    graph.react("tagment", "Tagment transcript and hash products", StepKind.TAGMENTATION,
                ((duplex, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged, waste))
    graph.discard(waste, "does not retain an amplifiable transcript or hash molecule")
    reviewed = ["tissue", "section_hash", "label_section", "grid_hash", "transfer_grid",
                "rt_primer", "reverse_transcriptase", "reverse_transcribe",
                "second_strand_primer", "second_strand_polymerase", "second_strand",
                "transposome", "tn5", "tagment"]
    _finish_pcr(graph, "library", tagged, final, primers, d.endpoint, reviewed,
                label="sci-Space")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _sciplex_graph(d, final):
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    nuclei = graph.input("nuclei", "permeabilized nuclei separated by sample",
                         MaterialKind.BIOLOGICAL, state=_sciplex_atac().initial, primary=True)
    hash_oligo = graph.input("hash_oligo", "sample-specific polyadenylated ssDNA hash",
                             MaterialKind.OLIGO)
    hashed = graph.product("hashed", "nuclei carrying absorbed sample hashes",
                           MaterialKind.MOLECULAR,
                           state=molecule("hashed nucleus", (HASH, role("accessible chromatin"))))
    graph.react("hash", "Let nuclei absorb the sample-specific hash oligo",
                StepKind.HYBRIDIZATION,
                ((nuclei, InputRole.TEMPLATE), (hash_oligo, InputRole.STRUCTURAL)), (hashed,))
    hash_primer = graph.input("hash_primer", "hash extension primer", MaterialKind.OLIGO)
    hash_polymerase = graph.input("hash_polymerase", "hash extension polymerase",
                                  MaterialKind.CATALYST)
    extended = graph.product("extended", "extended recoverable hash in each nucleus",
                             MaterialKind.MOLECULAR, state=hashed.state)
    graph.react("extend_hash", "Extend the absorbed sample hash", StepKind.EXTENSION,
                ((hashed, InputRole.TEMPLATE), (hash_primer, InputRole.STRUCTURAL),
                 (hash_polymerase, InputRole.CATALYST)), (extended,))
    transposome = graph.input("transposome", "first-well-indexed ATAC transposome",
                              MaterialKind.ADAPTER)
    tn5 = graph.input("tn5", "Tn5 transposase", MaterialKind.CATALYST)
    tagged = graph.product("tagged", "sample-hashed first-round-indexed nuclei",
                           MaterialKind.MOLECULAR,
                           state=molecule("first-round indexed nucleus", (HASH, CELL_A,
                               role("tagmented accessible DNA"))))
    waste = graph.product("tn5_waste", "non-library tagmentation products", MaterialKind.WASTE)
    graph.react("tagment", "Tagment chromatin with a first-round well index",
                StepKind.TAGMENTATION,
                ((extended, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged, waste))
    graph.discard(waste, "does not retain a recoverable hash or ATAC insert")
    second = graph.input("second_round", "second-round combinatorial-index wells",
                         MaterialKind.CONDITION)
    indexed = graph.product("indexed", "two-level-indexed hashed nuclei",
                            MaterialKind.MOLECULAR,
                            state=molecule("two-level cell identity", (HASH, CELL_A, CELL_B,
                                role("accessible-DNA fragment"))))
    collisions = graph.product("collisions", "ambiguous combinatorial-index collisions",
                               MaterialKind.WASTE)
    graph.react("partition", "Pool and redistribute into second-index wells",
                StepKind.PARTITION,
                ((tagged, InputRole.TEMPLATE), (second, InputRole.CONDITION)),
                (indexed, collisions), changes_molecule=False)
    graph.discard(collisions, "does not have an unambiguous two-level nucleus identity")
    atac = graph.product("atac_branch", "ATAC fragment branch", MaterialKind.MOLECULAR,
                         state=indexed.state)
    hashes = graph.product("hash_branch", "sample-hash readout branch", MaterialKind.MOLECULAR,
                           state=indexed.state)
    graph.react("split_readouts", "Separate ATAC-fragment and hash readouts",
                StepKind.SPLIT, ((indexed, InputRole.TEMPLATE),), (atac, hashes),
                changes_molecule=False)
    reviewed = ["nuclei", "hash_oligo", "hash", "hash_primer", "hash_polymerase",
                "extend_hash", "transposome", "tn5", "tagment", "second_round",
                "partition", "split_readouts"]
    atac_final = role_library("sciPlex-ATAC fragment roles",
                              (HASH, CELL_A, CELL_B, role("ATAC insert")), index=SAMPLE)
    hash_final = role_library("sciPlex-ATAC hash-read roles",
                              (HASH, CELL_A, CELL_B), index=SAMPLE)
    _finish_pcr(graph, "atac", atac, atac_final, (), "ATAC branch: " + d.endpoint,
                reviewed, label="sciPlex-ATAC")
    _finish_pcr(graph, "hashes", hashes, hash_final, (), "hash branch: " + d.endpoint,
                reviewed, label="sciPlex hash")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _dbit_graph(d, final, primers, *, cuttag):
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    dna_kind = "histone-mark DNA" if cuttag else "accessible DNA"
    tissue = graph.input("tissue", "fixed permeabilized tissue containing DNA and poly(A) RNA",
                         MaterialKind.BIOLOGICAL,
                         state=molecule("fixed tissue targets", (role(dna_kind), role("poly(A) RNA"))),
                         primary=True)
    substrate = tissue
    reviewed = ["tissue"]
    if cuttag:
        antibody = graph.input("antibody", "histone-mark primary antibody",
                               MaterialKind.SELECTION)
        bound = graph.product("antibody_bound", "antibody-bound tissue chromatin",
                              MaterialKind.MOLECULAR, state=tissue.state)
        unbound = graph.product("unbound_chromatin", "unbound chromatin", MaterialKind.WASTE)
        graph.react("bind_antibody", "Bind antibody to the selected histone mark",
                    StepKind.BINDING,
                    ((tissue, InputRole.TEMPLATE), (antibody, InputRole.SELECTION)),
                    (bound, unbound))
        graph.discard(unbound, "does not retain the selected histone-mark complex")
        substrate = bound
        reviewed.extend(("antibody", "bind_antibody"))
    transposome = graph.input("transposome", "RT-handle-bearing loaded Tn5 adapter",
                              MaterialKind.ADAPTER)
    tn5 = graph.input("tn5", "antibody-tethered pA–Tn5" if cuttag else "Tn5 transposase",
                      MaterialKind.CATALYST)
    tagged = graph.product("tagged_tissue", "tissue with tagged DNA and retained poly(A) RNA",
                           MaterialKind.MOLECULAR,
                           state=molecule("shared ligation handle on both branches", (
                               role(dna_kind + " carrying the shared handle"), role("poly(A) RNA"))))
    waste = graph.product("tn5_waste", "non-library tagmentation products", MaterialKind.WASTE)
    graph.react("tagment", "Tag " + dna_kind + " with the shared-handle transposome",
                StepKind.TAGMENTATION,
                ((substrate, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged, waste))
    graph.discard(waste, "does not retain a spatial-barcode-compatible genomic fragment")
    rt_primer = graph.input("rt_primer", "phosphorylated shared-handle oligo(dT) RT primer",
                            MaterialKind.OLIGO)
    rt = graph.input("reverse_transcriptase", "reverse transcriptase", MaterialKind.CATALYST)
    paired = graph.product("paired_branches", "DNA and cDNA branches bearing the shared handle",
                           MaterialKind.MOLECULAR,
                           state=molecule("shared-handle DNA and RNA branches", (
                               role(dna_kind + " fragment"), UMI, role("transcript cDNA"))))
    graph.react("reverse_transcribe", "Reverse-transcribe RNA through the shared RT handle",
                StepKind.REVERSE_TRANSCRIPTION,
                ((tagged, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (paired,))
    a_adapter = graph.input("barcode_a", "barcode-A oligo pool and linker 1",
                            MaterialKind.ADAPTER)
    a_ligase = graph.input("a_ligase", "barcode-A DNA ligase", MaterialKind.CATALYST)
    a_product = graph.product("a_product", "A-barcoded DNA and cDNA branches",
                              MaterialKind.MOLECULAR, state=paired.state)
    g1 = AdapterLigationGeometry.known(
        "splinted", SPATIAL_A_ROUND.adapter_scene(),
        statement="Linker 1 is the reverse complement of the shared RT/Tn5ME-B handle followed by the barcode-A ligating arm.")
    graph.react("ligate_a", "Ligate barcode A in the first channel set", StepKind.LIGATION,
                ((paired, InputRole.TEMPLATE), (a_adapter, InputRole.STRUCTURAL),
                 (a_ligase, InputRole.CATALYST)), (a_product,), geometry=g1)
    b_adapter = graph.input("barcode_b", "barcode-B oligo pool and linker 2",
                            MaterialKind.ADAPTER)
    b_ligase = graph.input("b_ligase", "barcode-B DNA ligase", MaterialKind.CATALYST)
    b_product = graph.product("b_product", "two-dimensionally addressed DNA and cDNA branches",
                              MaterialKind.MOLECULAR, state=a_product.state)
    g2 = AdapterLigationGeometry.known(
        "splinted", SPATIAL_B_ROUND.adapter_scene(),
        statement="Linker 2 is the reverse complement of the barcode-A outgoing arm followed by the barcode-B ligating arm.")
    graph.react("ligate_b", "Rotate the device and ligate barcode B", StepKind.LIGATION,
                ((a_product, InputRole.TEMPLATE), (b_adapter, InputRole.STRUCTURAL),
                 (b_ligase, InputRole.CATALYST)), (b_product,), geometry=g2)
    template, _same_final, geometry, _same_primers = _spatial_dna_libraries(dna_kind)
    dna = graph.product("dna_branch", dna_kind + " spatial branch", MaterialKind.MOLECULAR,
                        state=template)
    rna = graph.product("rna_branch", "matched spatial RNA branch", MaterialKind.MOLECULAR,
                        state=molecule("spatial RNA branch", (SPATIAL_A, SPATIAL_B, UMI,
                            role("transcript cDNA"))))
    graph.react("split_branches", "Recover matched DNA and RNA branches",
                StepKind.SPLIT, ((b_product, InputRole.TEMPLATE),), (dna, rna),
                changes_molecule=False)
    reviewed.extend(("transposome", "tn5", "tagment", "rt_primer", "reverse_transcriptase",
                     "reverse_transcribe", "barcode_a", "a_ligase", "ligate_a",
                     "barcode_b", "b_ligase", "ligate_b", "split_branches"))
    _finish_pcr(graph, "dna", dna, final, primers, "DNA branch: " + d.endpoint,
                reviewed, label=d.title + " DNA", geometry=geometry)
    rna_transposome = graph.input("rna_transposome", "RNA-branch library transposome",
                                  MaterialKind.ADAPTER)
    rna_tn5 = graph.input("rna_tn5", "RNA-branch Tn5 transposase", MaterialKind.CATALYST)
    rna_tagged = graph.product("rna_tagged", "tagmented spatial cDNA branch",
                               MaterialKind.MOLECULAR, state=rna.state)
    rna_waste = graph.product("rna_tn5_waste", "non-library RNA-branch tagmentation products",
                              MaterialKind.WASTE)
    graph.react("tagment_rna", "Tagment the pre-amplified RNA branch",
                StepKind.TAGMENTATION,
                ((rna, InputRole.TEMPLATE), (rna_transposome, InputRole.STRUCTURAL),
                 (rna_tn5, InputRole.CATALYST)), (rna_tagged, rna_waste))
    graph.discard(rna_waste, "does not retain a spatially addressed transcript insert")
    reviewed.extend(("rna_transposome", "rna_tn5", "tagment_rna"))
    rna_final = role_library(d.title + " RNA-library roles",
                             (SPATIAL_A, SPATIAL_B, UMI, role("transcript insert")), index=SAMPLE)
    _finish_pcr(graph, "rna", rna_tagged, rna_final, (), "RNA branch: " + d.endpoint,
                reviewed, label=d.title + " RNA")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def get(key: str):
    d = DEFS[key]
    workflows = {
        "hdst": _hdst,
        "sci-space": _sci_space,
        "spatial-atac-rna": lambda: _dbit(cuttag=False),
        "spatial-cuttag-rna": lambda: _dbit(cuttag=True),
        "slide-seq": _slide_seq,
        "sciplex-atac": _sciplex_atac,
    }
    payloads = {
        "hdst": (SPATIAL_BEAD, UMI, role("transcript insert")),
        "sci-space": (SPACE_HASH, SECTION_HASH, CELL_A, CELL_B, UMI,
                      role("transcript or hash insert")),
        "spatial-atac-rna": (SPATIAL_A, SPATIAL_B, UMI, role("ATAC or RNA insert")),
        "spatial-cuttag-rna": (SPATIAL_A, SPATIAL_B, UMI, role("CUT&Tag or RNA insert")),
        "slide-seq": (SPATIAL_BEAD, UMI, role("3-prime transcript insert")),
        "sciplex-atac": (HASH, CELL_A, CELL_B, role("ATAC insert")),
    }
    if key == "sci-space":
        final, primers = sci_space_library()
    elif key == "hdst":
        final, primers = hdst_library()
    elif key in ("spatial-atac-rna", "spatial-cuttag-rna"):
        dna_kind = ("accessible DNA" if key == "spatial-atac-rna"
                    else "histone-mark DNA")
        _template, final, _geometry, primers = _spatial_dna_libraries(dna_kind)
    else:
        final = role_library(d.title + " final library roles", payloads[key], index=SAMPLE)
        primers = ()
    spec = SimpleNamespace(
        TITLE=d.title, NOTES=d.notes,
        SOURCE=f'Defining source: <a href="https://doi.org/{d.doi}">{d.citation}</a>.',
        SUMMARY=d.summary,
        CAVEAT=("The printed hash, section hash and transcript are separate molecules that are associated through shared nucleus indices; the final construct shows the spatial-hash branch, not a fictitious covalent fusion to a transcript." if key == "sci-space" else
                "<inf>Individual barcode pools and decoder sequences were proprietary. Their published lengths and Read 1 coordinates are preserved, while unavailable bridge bases and the compatible completed Illumina shell are marked inferred.</inf>" if key == "hdst" else ""),
        FINAL_LIBRARY=final, SEQ_PRIMERS=primers,
        FINAL_CAPTION=d.endpoint, SEQUENCING_INTRO=d.endpoint, READ_LENGTHS=None,
        SEQUENCING_ENDING=("" if primers else d.endpoint),
        SEQUENCING_UNAVAILABLE=("" if primers else "The defining source establishes the read roles; kit revisions and barcode panels prevent one universal sequencing-primer sequence from being asserted here."),
        workflow=workflows[key],
        oligos=(_spatial_oligos if key in ("spatial-atac-rna", "spatial-cuttag-rna")
                else (lambda: (f"Bead primer precursor, as written — {HDST_PRECURSOR_WRITTEN}",)
                      if key == "hdst" else ())),
    )
    factories = {
        "hdst": lambda: _hdst_graph(d, final, primers),
        "sci-space": lambda: _sci_space_graph(d, final, primers),
        "spatial-atac-rna": lambda: _dbit_graph(d, final, primers, cuttag=False),
        "spatial-cuttag-rna": lambda: _dbit_graph(d, final, primers, cuttag=True),
        "slide-seq": lambda: _slide_graph(d, final),
        "sciplex-atac": lambda: _sciplex_graph(d, final),
    }
    spec.protocol_graph = factories[key]
    return spec
