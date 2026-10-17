"""Bounded molecular constructors for spatial capture and sample hashing."""
from __future__ import annotations

import illumina as il
import nextera as nx
import seqprimers as sp
from batch_ngs import joined_scene, seg, truseq_library
from chemdraw import Construct, Scene, feature, revcomp
from multimodal_spatial import HybridizationLigationRound
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import InputRole, MaterialKind, SequencingPlan, StepKind
from protocol_model import role


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

# Rodriques et al. 2019, supplementary methods/Table S1.  J denotes the
# split-and-pool bead barcode and N the molecule-specific UMI in the source.
SLIDE_BEAD_HANDLE = "CTACACGACGCTCTTCCGATCT"
SLIDE_BARCODE_1 = "J" * 8
SLIDE_BRIDGE = "TCTTCAGCGTTCCCGAGA"
SLIDE_BARCODE_2 = "J" * 7
SLIDE_UMI = "N" * 8
SLIDE_POLYT = "T" * 30
SLIDE_SMART_PCR = "AAGCAGTGGTATCAACGCAGAGT"
SLIDE_TSO_WRITTEN = "AAGCTGGTATCAACGCAGAGTGAATrG+GrG"
SLIDE_TRUSEQ5 = il.TRUSEQ_P5_FULL


# Zhang et al., Nature 2023 (doi:10.1038/s41586-023-05795-1): Methods
# ("Preparation of transposome") for Tn5ME-A/Tn5ME-B/Tn5MErev, Supplementary
# Table 4 for the RT primer, ligation linkers and N501/N701, and Supplementary
# Tables 5 and 6 for the barcode A/B plates.  The first member of each
# 100-oligo pool is retained as a representative; the molecular feature records
# the full published plate as the identifier source.  The letters are identical
# to the later Li et al. 2025 Nature Protocols Supplementary Table 1.  The RT
# primer's internal /iBiodT/ is the first T after the UMI; SPATIAL_RT_POLYT is
# the remaining T15VN.
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
    whitelist="100-member barcode A plate (Zhang et al. 2023 Supplementary Table 5)")
_SPATIAL_B_FEATURE = feature(
    "spatial_b", "spatial_barcode", "whitelist",
    whitelist="100-member barcode B plate (Zhang et al. 2023 Supplementary Table 6)")

SPATIAL_A_ROUND = HybridizationLigationRound(
    "barcode A", SPATIAL_RT_HANDLE, SPATIAL_A_PREFIX, SPATIAL_A_SUFFIX,
    SPATIAL_LINKER_1, _SPATIAL_A_FEATURE, phosphorylated_5=True,
    barcode_length=8)
SPATIAL_B_ROUND = HybridizationLigationRound(
    "barcode B", SPATIAL_A_PREFIX, SPATIAL_B_PREFIX, SPATIAL_B_SUFFIX,
    SPATIAL_LINKER_2, _SPATIAL_B_FEATURE, barcode_length=8)


def _spatial_branch_scene(kind: str, stage: int, *,
                          template_switched: bool | None = None) -> Scene:
    """Representative source-exact strand after RT/tagmentation and A/B ligation.

    ``stage`` counts completed in-situ barcode ligations (0, 1 or 2).  The RNA
    branch carries a TSO-derived end only after the separate on-bead template
    switch, so an RNA caller must say explicitly whether that step has happened.
    """
    if kind == "RNA":
        if template_switched is None:
            raise ValueError("the RNA branch needs an explicit template_switched stage")
        core = [seg("RT handle", SPATIAL_RT_HANDLE, "r1"),
                seg("10-nt UMI", "U" * 10, "umi", placeholder=True,
                    feature=feature("spatial_rna_umi", "umi", "random")),
                seg("internal biotin-dT (/iBiodT/)", "T"),
                seg("poly(dT)VN", SPATIAL_RT_POLYT, placeholder=True),
                seg("transcript cDNA", "X" * 28, placeholder=True)]
        if template_switched:
            # Zhang et al. 2023 Methods TSO, DNA letters of AAGCAGTGGTATCAACGCAGAGTGAATrGrG+G.
            core.append(seg("TSO-derived end", revcomp("AAGCAGTGGTATCAACGCAGAGTGAATGGG")))
    else:
        if template_switched is not None:
            raise ValueError("template switching applies only to the RNA branch")
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


def _spatial_dna_libraries(dna_kind: str):
    """Exact representative DNA branch before and after N501/N701 PCR.

    The source draws the barcode-bearing end P7-to-insert.  Constructs use the
    canonical P5-to-P7 top-strand orientation, hence the reverse complements of
    the B/A/RT chain on the right side of the insert.
    """
    i5 = feature("spatial_i5", "sample_index", "whitelist",
                 whitelist="N501 index from Zhang et al. 2023 Supplementary Table 4")
    i7 = feature("spatial_i7", "sample_index", "whitelist",
                 whitelist="N701 index (of N701-N707) from Zhang et al. 2023 Supplementary Table 4")
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
        statement="Zhang et al. 2023 Supplementary Table 4 gives the complete N501 and N701 primers; their 3′ arms land exactly on Tn5ME-A and the barcode-B PCR arm.")
    primers = sp.standard_primers(final)
    return template, final, geometry, primers


def slide_seq_library():
    barcode = feature("slide_seq.bead", "spatial_barcode", "whitelist",
                      group="slide_seq.bead", whitelist="decoded Slide-seq bead pool")
    umi = feature("slide_seq.umi", "umi", "random")
    lib = Construct([
        seg("P5", il.P5, "p5"),
        seg("TruSeq Read 1 remainder", il.TRUSEQ_READ1[4:], "r1"),
        seg("bead barcode block 1", "B" * 8, "cbc", placeholder=True,
            feature=barcode),
        seg("bead barcode bridge", SLIDE_BRIDGE),
        seg("bead barcode block 2", "B" * 7, "cbc", placeholder=True,
            feature=barcode),
        seg("UMI", "U" * 8, "umi", placeholder=True, feature=umi),
        seg("poly(dT)", SLIDE_POLYT),
        role("3-prime transcript cDNA"),
        seg("Nextera mosaic end reverse complement", nx.ME_RC, "me"),
        seg("S7 reverse complement", nx.S7_RC, "s7"),
        seg("i7 reverse complement", "I" * 8, "cbc", placeholder=True,
            feature=feature("slide_seq.i7", "sample_index", "whitelist",
                            whitelist="Nextera N700 index set")),
        seg("P7 reverse complement", il.P7_RC, "p7"),
    ], name="representative Slide-seq v1 sequencing library")
    # This intentionally is a mixed TruSeq/Nextera endpoint.  The completed molecule
    # also contains an incidental opposite-strand Nextera-R1 match, so protocol context
    # is required to select the source-declared barcode read.
    primers = (sp.TRUSEQ["R1"], sp.NEXTERA["I1"], sp.NEXTERA["R2"])
    problems = sp.verify(lib, primers,
                         required_roles=("Read 1", "Index 1 (i7)", "Read 2"))
    if problems:
        raise ValueError("invalid Slide-seq endpoint: " + "; ".join(problems))
    return lib, primers


HDST_PRECURSOR_WRITTEN = (
    "/AmC6/UUUUUGACTCGTAATACGACTCACTATAGGGACACGACGCTCTTCCGATCT-"
    "Spatial_barcode_Pool1-Bridge1"
)
HDST_BRIDGE12 = "GACTTGTCTAGAGC"
HDST_BRIDGE34 = "TGATGCCACACTACTC"


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
        seg("ligated bridge 1+2", HDST_BRIDGE12),
        seg("spatial barcode pool 2", "C" * 14, "cbc", placeholder=True, feature=bc2),
        seg("ligated bridge 3+4", HDST_BRIDGE34),
        seg("spatial barcode pool 3", "D" * 15, "cbc", placeholder=True, feature=bc3),
        seg("UMI", "U" * 5, "umi", placeholder=True, feature=umi),
        role("transcript cDNA"),
    ]
    lib, all_primers = truseq_library(insert, "representative HDST sequencing library",
                                      dual_index=False, inferred_adapters=True, i7_nt=8)
    primers = tuple(p for p in all_primers if p.role in ("Read 1", "Index 1 (i7)", "Read 2"))
    return lib, primers


def _finish_pcr(graph, prefix, template, final, primers, endpoint, reviewed, *, label,
                geometry=None, forward_state=None, reverse_state=None):
    """PCR to the sequenced library. Primer states are the caller's source-exact
    reagents; None leaves that primer without an exact state."""
    forward = graph.input(prefix + "_forward", label + " forward PCR primer",
                          MaterialKind.OLIGO, state=forward_state)
    reverse = graph.input(prefix + "_reverse", label + " reverse PCR primer",
                          MaterialKind.OLIGO, state=reverse_state)
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
