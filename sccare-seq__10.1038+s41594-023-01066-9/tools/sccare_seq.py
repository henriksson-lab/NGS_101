"""scCARE-seq molecular model."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import illumina as il
import nextera as nx
import seqprimers as sp
from chemdraw import Construct, Row, Segment, feature, oligo, revcomp
from multimodal_spatial import finish_pcr_branch
from pcr import PCRGeometry, PCRPrimer
from protocol_extensions import proximity_product
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind

RT_HANDLE = "ACGACGCTCTTCCGATCT"
RT_BARCODES = (
    "CTTAGGAC", "GTTGGAGC", "TCTTACGC", "GATGGTCA",
    "GAGAATCC", "CGCAACGA", "CGGACCTA", "GAATTAAC",
)
UMI = "N" * 8
POLYT = "T" * 30 + "VN"

# Representative members of the three complete index-primer sets printed in
# Supplementary Table 1. Their invariant arms come from the shared Illumina model.
N501_INDEX = "TAGATCGC"
N701_INDEX = "TCGCCTTA"
T501_INDEX = "TATAGCCT"
N501 = nx.n5xx_primer(N501_INDEX)
N701 = nx.n7xx_primer(N701_INDEX)
T501 = il.nebnext_i5_primer(T501_INDEX)
COMMON_PRIMER = nx.S7 + nx.ME[:13]
DNA_PRIMER = nx.S5 + nx.ME[:13]
RNA_PRIMER = "CCCTA" + RT_HANDLE
TN5_N5 = nx.ADAPTOR_S5
TN5_N7 = nx.ADAPTOR_S7
TN5_REV = nx.ME_RC
CARRIER = "TCAGGTTTTCCTGAA"

CELL = feature(
    "sccare_cell", "cell_barcode", "whitelist",
    whitelist="eight RT barcodes printed in Supplementary Table 1",
)
UMI_FEATURE = feature("sccare_umi", "umi", "random")
I5 = feature("sccare_i5", "sample_index", "whitelist",
             whitelist="N501–N516 or T501–T516 in Supplementary Table 1")
I7 = feature("sccare_i7", "sample_index", "whitelist",
             whitelist="N701–N724 in Supplementary Table 1")


def _dna_library() -> Construct:
    """N501/N701 representative, assembled only from source-printed parts."""
    return Construct([
        Segment("P5", il.P5, "p5"),
        Segment("i5 (N501)", N501_INDEX, "cbc", feature=I5),
        Segment("S5", nx.S5, "s5"),
        Segment("left mosaic end", nx.ME, "me"),
        Segment("contact-junction DNA", "X" * 46, placeholder=True),
        Segment("right mosaic end reverse complement", nx.ME_RC, "me"),
        Segment("S7 reverse complement", nx.S7_RC, "s7"),
        Segment("i7 (N701) reverse complement", revcomp(N701_INDEX), "cbc", feature=I7),
        Segment("P7 reverse complement", il.P7_RC, "p7"),
    ], name="scCARE-seq DNA-contact library (N501/N701 representative)")


def _rna_library() -> Construct:
    """Hybrid TruSeq-i5/Nextera-i7 RNA branch shown in Extended Data Fig. 1."""
    return Construct([
        Segment("P5", il.P5, "p5"),
        Segment("i5 (T501)", T501_INDEX, "cbc", feature=I5),
        Segment("TruSeq Read 1 arm", il.TRUSEQ_READ1, "r1"),
        Segment("UMI", UMI, "umi", placeholder=True, feature=UMI_FEATURE),
        Segment("cell barcode", RT_BARCODES[0], "cbc", feature=CELL),
        Segment("oligo(dT)VN / cDNA boundary", POLYT, placeholder=True),
        Segment("cDNA insert", "X" * 42, placeholder=True),
        Segment("mosaic end reverse complement", nx.ME_RC, "me"),
        Segment("S7 reverse complement", nx.S7_RC, "s7"),
        Segment("i7 (N701) reverse complement", revcomp(N701_INDEX), "cbc", feature=I7),
        Segment("P7 reverse complement", il.P7_RC, "p7"),
    ], name="scCARE-seq RNA library (T501/N701 representative)")


DNA = _dna_library()
DNA_P = tuple(sp.NEXTERA[k] for k in ("R1", "I1", "I2", "R2"))
RNA = _rna_library()
RNA_P = (sp.TRUSEQ["R1"], sp.NEXTERA["I1"], sp.TRUSEQ["I2"], sp.NEXTERA["R2"])
for _library, _primers in ((DNA, DNA_P), (RNA, RNA_P)):
    _problems = sp.verify(_library, _primers, required_roles=tuple(p.role for p in _primers))
    if _problems:
        raise ValueError(f"{_library.name}: " + "; ".join(_problems))

DNA_GEOMETRY = PCRGeometry.exact(
    DNA, PCRPrimer("N501", N501), PCRPrimer("N701", N701),
    statement="Supplementary Table 1 prints the complete N501 and N701 primers.",
)
RNA_GEOMETRY = PCRGeometry.exact(
    RNA, PCRPrimer("T501", T501), PCRPrimer("N701", N701),
    statement=("Supplementary Table 1 prints the complete T501 and N701 primers; "
               "Extended Data Fig. 1 assigns them to the RNA branch."),
)

TITLE = "scCARE-seq — chromatin architecture and RNA expression"
NOTES = "01_sccare-seq.html"
SOURCE = '<a href="https://doi.org/10.1038/s41594-023-01066-9">Deng et al., <i>Nature Structural &amp; Molecular Biology</i> (2023)</a>.'
SUMMARY = "A biotinylated, well-barcoded oligo(dT) labels RNA while in-situ tagmentation and proximity ligation preserve chromatin contacts. DNA and RNA branches are separated and indexed for paired analysis."
CAVEAT = "Supplementary Table 1 prints all eight RT barcodes and the complete N501–N516, N701–N724 and T501–T516 primer sets. The page uses Oligo-dT-1, N501, N701 and T501 as exact representatives; X and N remain variable biological or randomized sequence."
FINAL_LIBRARIES = (
    ("DNA-contact library", DNA, DNA_P, "N501/N701-amplified, proximity-ligated genomic fragments.", "Paired genomic reads report chromatin contacts; the two sample indices identify the indexed library."),
    ("RNA library", RNA, RNA_P, "T501/N701-amplified well-barcoded UMI cDNA.", "Read 1 encounters the UMI and cell barcode; Read 2 enters cDNA from the Nextera end."),
)


def oligos():
    yield oligo("Tn5-N5", [Segment("S5", nx.S5, "s5"), Segment("ME", nx.ME, "me")])
    yield oligo("Tn5-N7", [Segment("S7", nx.S7, "s7"), Segment("ME", nx.ME, "me")])
    yield oligo("Tn5-Rev", [Segment("ME reverse complement", TN5_REV, "me")], mods="/5Phos/")
    yield oligo("Carrier ssDNA", [Segment("carrier", CARRIER)])
    yield oligo("Common primer (N7-PreA)", [Segment("common primer", COMMON_PRIMER)])
    yield oligo("DNA primer (N5-PreA)", [Segment("DNA primer", DNA_PRIMER)])
    yield oligo("RNA primer (A5-PreA)", [Segment("RNA primer", RNA_PRIMER)])
    yield oligo("Oligo-dT-1", [
        Segment("RNA handle", RT_HANDLE, "r1"),
        Segment("UMI", UMI, "umi", placeholder=True, feature=UMI_FEATURE),
        Segment("cell barcode 1", RT_BARCODES[0], "cbc", feature=CELL),
        Segment("oligo(dT)VN", POLYT, placeholder=True),
    ], mods="one dT is biotinylated in the capture version")
    yield oligo("N501", [Segment("complete N501 primer", N501)])
    yield oligo("N701", [Segment("complete N701 primer", N701)])
    yield oligo("T501", [Segment("complete T501 primer", T501)])


def sections():
    return [
        ("Fix and tagment nuclei", [Row(chunks=[("Tn5-N5 ** chromatin ** Tn5-N7", "me", False)])], "Printed Nextera transposome oligos install DNA handles."),
        ("Reverse-transcribe and label RNA", [Row(chunks=[("RNA handle — UMI — well barcode — oligo(dT) → poly(A) RNA", "umi", False)])], "A biotin-dT version permits physical recovery of the RNA branch."),
        ("Generate contact junctions", proximity_product(left="tagged locus A", right="tagged locus B", bridge="proximity-ligation junction").rows(), "Crosslinked fragments are ligated before single-cell separation."),
        ("Separate and amplify", [Row(chunks=[("N501/N701 DNA-contact PCR ← split → T501/N701 RNA PCR", None, False)])], "The supplement assigns an exact, different i5 primer family to each branch; both use N7xx on the other end."),
    ]


def protocol_graph():
    g = ProtocolGraph(TITLE, authority=SOURCE)
    nuclei = g.input("nuclei", "fixed nuclei containing chromatin and poly(A) RNA", MaterialKind.BIOLOGICAL, state=sections()[0][1], primary=True)
    transposome = g.input("transposome", "printed N5/N7-loaded Nextera transposome", MaterialKind.ADAPTER)
    tn5 = g.input("tn5", "Tn5 transposase", MaterialKind.CATALYST)
    tagged = g.product("tagged", "nuclei retaining Tn5-tagged chromatin and RNA", MaterialKind.MOLECULAR, state=sections()[0][1])
    twaste = g.product("tn5_waste", "non-library chromatin tagmentation products", MaterialKind.WASTE)
    g.react("tagment", "Tagment chromatin in fixed nuclei", StepKind.TAGMENTATION, ((nuclei, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL), (tn5, InputRole.CATALYST)), (tagged, twaste))
    g.discard(twaste, "does not retain a selected N5/N7 chromatin fragment")
    rt_primer = g.input("rt_primer", "biotinylated well-barcoded UMI oligo(dT) primer", MaterialKind.OLIGO)
    rt = g.input("rt", "reverse transcriptase", MaterialKind.CATALYST)
    copied = g.product("copied", "nuclei containing barcoded biotin-cDNA and tagged chromatin", MaterialKind.MOLECULAR, state=sections()[1][1])
    g.react("reverse_transcribe", "Reverse-transcribe and label RNA", StepKind.REVERSE_TRANSCRIPTION, ((tagged, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL), (rt, InputRole.CATALYST)), (copied,))
    ligase = g.input("contact_ligase", "in-situ proximity ligase", MaterialKind.CATALYST)
    contacts = g.product("contacts", "nuclei containing proximity-ligated chromatin and biotin-cDNA", MaterialKind.MOLECULAR, state=sections()[2][1])
    g.react("proximity_ligate", "Generate chromatin-contact junctions", StepKind.PROXIMITY_LIGATION, ((copied, InputRole.TEMPLATE), (ligase, InputRole.CATALYST)), (contacts,))
    streptavidin = g.input("streptavidin", "streptavidin affinity matrix", MaterialKind.SELECTION)
    rna = g.product("rna_branch", "biotin-captured cDNA branch", MaterialKind.MOLECULAR, state=RNA)
    dna = g.product("dna_branch", "DNA-contact flow-through branch", MaterialKind.MOLECULAR, state=DNA)
    g.react("separate", "Separate biotin-cDNA from DNA-contact products", StepKind.AFFINITY_CAPTURE, ((contacts, InputRole.TEMPLATE), (streptavidin, InputRole.SELECTION)), (rna, dna), changes_molecule=False)
    reviewed = ["nuclei", "transposome", "tn5", "tagment", "rt_primer", "rt", "reverse_transcribe", "contact_ligase", "proximity_ligate", "streptavidin", "separate"]
    finish_pcr_branch(g, "dna", dna, DNA, DNA_P, reviewed, label="scCARE-seq DNA-contact", geometry=DNA_GEOMETRY)
    finish_pcr_branch(g, "rna", rna, RNA, RNA_P, reviewed, label="scCARE-seq RNA", geometry=RNA_GEOMETRY)
    g.declare_reviewed_scope(reviewed, SOURCE)
    return g.require_publishable()
