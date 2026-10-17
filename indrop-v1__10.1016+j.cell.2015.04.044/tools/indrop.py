"""Molecular constructors for the published inDrop v1 and v2 protocols."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from chemdraw import Construct, Scene, Segment, feature, revcomp
import illumina as il
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_graph import (Disposition, InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, SourceClaim, SourceDiscrepancy, StepKind,
                            Terminal)

V1 = "inDrop v1"
V2 = "inDrop v2"

LEADER = "CGATGACG"
T7_PROMOTER = "TAATACGACTCACTATAGGG"
MID = "ATACCACCATGG"
PE1 = il.TRUSEQ_READ1[3:]
W1_STAR = "AAGGCGTCACAAGCAATCACTC"
W1 = revcomp(W1_STAR)
BC2_LEN = 8
UMI_LEN = 6
DT_LEN = 19
PE2 = il.PE_READ2_ARM

# Cell 2015 Supplemental Information (mmc1) Table S5, "2. Library preparation as per
# ref. (2)", as printed. Table S5 prints no sequencing primers.
V1_RLO_WRITTEN = "/5Phos/AGATCGGAAGAGCGGTTCAGCAGGAATGCC/3SpC3/"
RLO = "AGATCGGAAGAGCGGTTCAGCAGGAATGCC"
V1_RT2 = "GTCTCGGCATTCCTGCTGAAC"
V1_PCR1 = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGA"
V1_PCR2 = "CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAAC"
# Table S5 numbers the barcode-1 positions 1-10 in both the plate template and the
# assembled primer; the 2017 v2 plate sheet instead lists 8-11 nt.
V1_BC1_LEN = 10
V1_READ1_SITE = il.TRUSEQ_READ1[4:]
if RLO != revcomp(PE2[2:]):
    raise ValueError("v1 RLO is not the complement of the PE2 read site")
if revcomp(V1_RT2)[:16] != RLO[-16:]:
    raise ValueError("v1 2nd RT primer 3' 16 nt must pair with the RLO 3' end")
if V1_PCR2 != il.P7 + "CG" + V1_RT2 or V1_PCR1 != il.P5 + V1_READ1_SITE[:15]:
    raise ValueError("v1 PCR primers are not P7/P5 tails on the RT2 and PE1 arms")
if not PE1.endswith(V1_READ1_SITE):
    raise ValueError("v1 PCR primer 1 does not land on the bead PE1 arm")
# Bases the 2nd RT primer and PCR primer 2 add beyond the RLO, before P7'.
V1_RT2_EXTENSION = revcomp(V1_PCR2)[16:-len(il.P7)]
V2_R1 = PE2[2:]
V2_INDEX = revcomp(PE1)
V2_R2 = PE1

SOURCE_ACRYDITE_BASES = "CGATGACGTAATACGACTCACTATAGGGATACCACCATGGCTCTTTCCCTACACGACGCTCTTC"
SOURCE_PE2_N6 = "TCGGCATTCCTGCTGAACCGCTCTTCCGATCTNNNNNN"
BC1_FEATURE = feature("cell_barcode_1", "cell_barcode", "combinatorial",
                      group="cell_barcode", part="barcode 1")
BC2_FEATURE = feature("cell_barcode_2", "cell_barcode", "combinatorial",
                      group="cell_barcode", part="barcode 2")
UMI_FEATURE = feature("umi", "umi", "random")
I7_FEATURE = feature("sample_index_i7", "sample_index", "unknown")
FEATURES = {"barcode 1": BC1_FEATURE, "barcode 1'": BC1_FEATURE,
            "barcode 2": BC2_FEATURE, "barcode 2'": BC2_FEATURE,
            "UMI": UMI_FEATURE, "i7 read": I7_FEATURE}


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    kw.setdefault("feature", FEATURES.get(name))
    return Segment(name=name, top=top, tag=tag, **kw)


def var(name: str, n: int, letter: str, tag: str | None = None,
        inferred: bool = False) -> Segment:
    return seg(name, letter * n, tag, placeholder=True, inferred=inferred)


def bead_primer(bc1_len: int = 8, *, inferred: bool = False) -> Construct:
    if bc1_len not in (8, 9, 10, 11):
        raise ValueError("inDrop barcode 1 must be 8-11 nt")
    return Construct([
        seg("leader", LEADER, inferred=inferred), seg("T7 promoter", T7_PROMOTER, "tso", inferred=inferred),
        seg("middle", MID, inferred=inferred), seg("PE1", PE1, "r2", inferred=inferred),
        var("barcode 1'", bc1_len, "A", "cbc", inferred), seg("W1", W1, "r3", inferred=inferred),
        var("barcode 2'", BC2_LEN, "B", "cbc", inferred), var("UMI", UMI_LEN, "N", "umi", inferred),
        seg("T19", "T" * DT_LEN, inferred=inferred), var("V", 1, "V", None, inferred),
    ], name="finished inDrop bead primer")


SI_LIBRARY = ("Cell 2015 Supplemental Information (mmc1) 'DNA-Library Preparation', "
              "lines 721-736")
SI_BHM = ("Cell 2015 Supplemental Information (mmc1) 'BHM split-pool combinatorial "
          "barcoding', lines 611-634")
SI_DROPLET = ("Cell 2015 Supplemental Information (mmc1) 'Reverse transcription/lysis mix', "
              "lines 567-571")

_V1_PRIMER_SOURCE = ('Illumina "Illumina Adapter Sequences" #1000000002694; Cell 2015 '
                     'Table S5 prints no sequencing primer, so this canonical primer is '
                     'derived from the library ends its printed oligos build')
V1_PRIMERS = (
    sp.custom("Read 1", "Illumina Read 1 (assumed)", il.TRUSEQ_READ1, _V1_PRIMER_SOURCE),
    sp.custom("Read 2", "Illumina PE Read 2 (assumed)", il.PE_READ2_PRIMER, _V1_PRIMER_SOURCE),
)
V2_PRIMERS = (
    sp.custom("Read 1", "Custom Read 1", V2_R1, "Nature Protocols supplementary oligo table"),
    sp.custom("Index 1 (i7)", "Custom Index Read", V2_INDEX, "Nature Protocols supplementary oligo table"),
    sp.custom("Read 2", "Custom Read 2", V2_R2, "Nature Protocols supplementary oligo table"),
)
# No acquired cycle count is stated: neither the Cell 2015 Experimental Procedures nor
# its Supplemental Information (Table S1 lists instruments only) gives read lengths.
# The GEO GSE65525 records name the instruments but no cycle counts; SRA runinfo for
# SRP053052 gives only combined average spot lengths (89/122/133), not per-read lengths.
V1_READ_LENGTHS: dict[str, int] = {}
V1_READ_USES = {"Read 1": "extract", "Read 2": "include"}


def v1_library(bc1_len: int = V1_BC1_LEN) -> Construct:
    """v1 library computed from the Cell 2015 Table S5 bead, RLO, RT and PCR oligos."""
    con = Construct([
        seg("P5", il.P5, "p5"),
        seg("Read 1 site", V1_READ1_SITE, "r1"),
        var("barcode 1'", bc1_len, "A", "cbc"), seg("W1", W1, "r3"),
        var("barcode 2'", BC2_LEN, "B", "cbc"), var("UMI", UMI_LEN, "N", "umi"),
        seg("T19", "T" * DT_LEN), var("V", 1, "V", None),
        var("antisense mRNA", 19, "X", "r1"), seg("RLO", RLO, "r2"),
        seg("RT2 / PCR primer 2 extension", V1_RT2_EXTENSION, "r2"),
        seg("P7'", il.P7_RC, "p7"),
    ], name="inDrop v1 library")
    _validate(con, V1_PRIMERS)
    return con


def v1_plate_oligos(bc1_len: int = V1_BC1_LEN) -> tuple[Construct, Construct]:
    """Cell 2015 Table S5 barcode templates (W1*-bc1-PE1*, T19V*-UMI-bc2-W1*)."""
    return (Construct([seg("W1*", W1_STAR, "r3"), var("barcode 1", bc1_len, "A", "cbc"),
                       seg("PE1*", revcomp(PE1), "r2")], name="barcode 1 template"),
            Construct([var("B", 1, "B"), seg("A19", "A" * DT_LEN),
                       var("UMI", UMI_LEN, "N", "umi"), var("barcode 2", BC2_LEN, "B", "cbc"),
                       seg("W1*", W1_STAR, "r3")], name="barcode 2/UMI template"))


def rlo_ligation_scene() -> Scene:
    """Fragmented aRNA 3' end beside the printed 5'-phosphorylated, 3'-C3-blocked RLO."""
    sc = Scene()
    sc.strand("aRNA", [Segment("aRNA fragment", "X" * 20, placeholder=True)],
              label="fragmented antisense RNA")
    sc.strand("RLO", [Segment("RLO", RLO, "r2")], label="RLO (/5Phos/ … /3SpC3/)",
              at=26, mod5="p", mod3="C3")
    sc.ligation_end("aRNA", "aRNA fragment", kind="single-strand", polarity="3-prime",
                    owner="substrate", text="aRNA 3′-OH")
    sc.ligation_end("RLO", "RLO", kind="single-strand", polarity="5-prime",
                    owner="adapter", text="5′-phosphate RLO end")
    return sc


def v2_library(bc1_len: int = 8) -> Construct:
    if bc1_len not in (8, 9, 10, 11): raise ValueError("barcode 1 must be 8-11 nt")
    con = Construct([
        seg("P5", il.P5, "p5"), seg("GGTC", "GGTC", "r1"), seg("PE2", PE2, "r1"),
        var("sense mRNA", 19, "X", "r1"), seg("A19", "A" * DT_LEN), var("B", 1, "B"),
        var("UMI", UMI_LEN, "N", "umi"), var("barcode 2", BC2_LEN, "B", "cbc"),
        seg("W1*", W1_STAR, "r3"), var("barcode 1", bc1_len, "A", "cbc"),
        seg("PE1*", revcomp(PE1), "r2"), var("i7 read", 6, "I", "cbc"),
        seg("P7'", il.P7_RC, "p7"),
    ], name="inDrop v2 library")
    _validate(con, V2_PRIMERS)
    return con


def _validate(con: Construct, primers) -> None:
    for p in primers:
        if sp.locate(con, p) is None:
            raise ValueError(f"{con.name}: {p.name} site missing")


def library(protocol: str) -> Construct:
    if protocol == V1: return v1_library()
    if protocol == V2: return v2_library()
    raise ValueError(protocol)


def protocol_graph(protocol: str) -> ProtocolGraph:
    doi = "10.1016/j.cell.2015.04.044" if protocol == V1 else "10.1038/nprot.2016.154"
    source = f"Defining source: https://doi.org/{doi}"
    v1 = protocol == V1
    boundary = source + (" (Supplemental Information Table S5 oligos)" if v1 else "")
    graph = ProtocolGraph(protocol, authority=source)
    bead_seed = graph.input("bead_seed", "/5Acryd//iSpPC/ acrydite photocleavable bead primer",
                            MaterialKind.OLIGO, source=boundary,
                            state=Construct([seg("leader", LEADER), seg("T7 promoter", T7_PROMOTER, "tso"),
                                             seg("middle", MID), seg("PE1 partial", PE1[:24], "r2")],
                                            name="hydrogel-incorporated DNA primer") if v1 else None)
    plate_oligos = graph.input("barcode_plate_oligos", "two split-pool barcode plate oligo sets",
                               MaterialKind.OLIGO, source=boundary,
                               state=v1_plate_oligos() if v1 else None)
    bead_polymerase = graph.input(
        "bead_polymerase",
        "Bst 2.0 DNA polymerase (per well: 85 °C 2 min, 60 °C 20 min hybridization, "
        "60 °C 60 min; EDTA stop; 150 mM NaOH strips the plate strand)"
        if v1 else "bead-primer extension polymerase",
        MaterialKind.CATALYST, source=SI_BHM if v1 else None)
    bead = graph.product("barcoded_bead", "finished split-pool barcode bead primer",
                         MaterialKind.OLIGO,
                         state=bead_primer(V1_BC1_LEN) if v1 else bead_primer(),
                         source=boundary)
    graph.react("bead_synthesis", "Build the bead barcode in two split-pool extensions",
                StepKind.EXTENSION,
                ((bead_seed, InputRole.TEMPLATE), (plate_oligos, InputRole.STRUCTURAL),
                 (bead_polymerase, InputRole.CATALYST)), (bead,))
    cells = graph.input("cells", "single-cell suspension", MaterialKind.BIOLOGICAL,
                        primary=True)
    loading = graph.input("droplet_loading", "one-cell/one-bead droplet loading condition",
                          MaterialKind.CONDITION)
    droplets = graph.product("productive_droplets", "droplets containing one cell and one barcode bead",
                             MaterialKind.MOLECULAR)
    excluded = graph.product("excluded_droplets", "empty or multiplet droplets", MaterialKind.WASTE)
    graph.react("partition", "Co-encapsulate cells and beads and photorelease primers",
                StepKind.PARTITION,
                ((cells, InputRole.TEMPLATE), (bead, InputRole.SPIKE_IN),
                 (loading, InputRole.CONDITION)), (droplets, excluded))
    graph.mark_terminal(Terminal(excluded.id, Disposition.UNRESOLVED,
        "Contents of empty and multiplet droplets are pooled when the emulsion is broken "
        "(Cell 2015 Library preparation; not stated for v2 in the inspected supplements); "
        "their exclusion is computational, not a physical discard."))
    if v1:
        graph.annotate(SourceDiscrepancy(
            "partition", "Published carrier-oil flow rate differs",
            (SourceClaim("Cell 2015 Experimental Procedures, 'Microfluidic Operation'",
                         "90 µL/hr carrier oil"),
             SourceClaim("Supplemental Information, 'Microfluidic device operation'",
                         "80 µL/hr carrier oil")),
            "The graph states no oil flow rate; the other three flow rates agree "
            "(100, 100 and 10-20 µL/hr) and both passages give 4 nL drops.",
            "Both values are printed in the published article and its supplement "
            "(mmc6 lines 1636-1640; mmc1 lines 540-543).",
            "Which rate was used for the reported runs is not stated."))
    rt = graph.input("reverse_transcriptase",
                     "SuperScript III reverse transcriptase (in-drop RT/lysis mix)" if v1
                     else "reverse transcriptase", MaterialKind.CATALYST,
                     source=SI_DROPLET if v1 else None)
    first = graph.product("barcoded_cdna", "bead-barcoded first-strand cDNA",
                          MaterialKind.MOLECULAR)
    rt_label = ("In-drop RT from released T19V primers (50 °C 2 h, 70 °C 15 min), "
                "then demulsify aliquots with perfluorooctanol" if protocol == V1 else
                "Reverse-transcribe mRNA from released T19V primers")
    graph.react("reverse_transcription", rt_label,
                StepKind.REVERSE_TRANSCRIPTION,
                ((droplets, InputRole.TEMPLATE), (bead, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (first,))
    if v1:
        graph.annotate(SourceDiscrepancy(
            "reverse_transcription", "Published demulsification volume differs",
            (SourceClaim("Cell 2015 Experimental Procedures, 'Library Preparation'",
                         "0.2X 20% perfluorooctanol / 80% HFE-7500, after splitting into "
                         "aliquots of 100-3,500 cells"),
             SourceClaim("Supplemental Information, 'Barcoding inside droplets'",
                         "1 volume of PFO solution (20% perfluorooctanol, 80% HFE-7500)")),
            "The RT step names perfluorooctanol demulsification without a volume.",
            "Both values are printed in the published article and its supplement "
            "(mmc6 lines 1652-1655; mmc1 lines 578-584). The two passages also name "
            "different UV lamps (UVP B-100; BlackRay Xenon).",
            "Which volume was used is not stated; the breaking chemistry is the same."))
    v1_extra: tuple[str, ...] = ()
    if v1:
        exo = graph.input("exonuclease_i", "Exonuclease I", MaterialKind.CATALYST,
                          source=SI_LIBRARY)
        hinfi = graph.input("hinfi", "HinFI (10 U), added to digest RT primer dimers",
                            MaterialKind.CATALYST, source=SI_LIBRARY)
        treated = graph.product("exo_treated_cdna",
                                "ExoI/HinFI-treated cDNA with digested primers and primer dimers",
                                MaterialKind.MOLECULAR, source=SI_LIBRARY)
        graph.react("exo_hinfi",
                    "Spin out debris and gels (4 °C, 15 min, 14 krcf), then ExoI with 10 U HinFI",
                    StepKind.DIGESTION,
                    ((first, InputRole.TEMPLATE), (exo, InputRole.CATALYST),
                     (hinfi, InputRole.CATALYST)), (treated,), source=SI_LIBRARY)
        spri = graph.input("spri_beads", "SPRI beads", MaterialKind.SELECTION,
                           source=SI_LIBRARY)
        purified = graph.product("purified_cdna", "SPRI-purified first-strand cDNA",
                                 MaterialKind.MOLECULAR, source=SI_LIBRARY)
        spri_waste = graph.product("spri_waste", "SPRI supernatant with digestion products",
                                   MaterialKind.WASTE)
        graph.react("spri_after_exo", "SPRI purification", StepKind.CLEANUP,
                    ((treated, InputRole.TEMPLATE), (spri, InputRole.SELECTION)),
                    (purified, spri_waste), source=SI_LIBRARY)
        graph.mark_terminal(Terminal(spri_waste.id, Disposition.DISCARD,
            "SPRI supernatant; what each enzyme removes is the authors' stated intent "
            "(HinFI 'to digest primer dimers'), not an observed separation product."))
        first = purified
        v1_extra = ("exonuclease_i", "hinfi", "exo_hinfi", "spri_beads",
                    "spri_after_exo", "spri_waste")
    second_pol = graph.input("second_strand_polymerase", "second-strand polymerase",
                             MaterialKind.CATALYST)
    ds = graph.product("t7_templates", "double-stranded T7-bearing cDNA templates",
                       MaterialKind.MOLECULAR)
    graph.react("second_strand", "Make the bead-borne T7 promoter double-stranded"
                + (" (second-strand synthesis), then SPRI" if v1 else ""),
                StepKind.EXTENSION,
                ((first, InputRole.TEMPLATE), (second_pol, InputRole.CATALYST)), (ds,))
    t7 = graph.input("t7_polymerase", "T7 RNA polymerase", MaterialKind.CATALYST)
    arna = graph.product("arna", "linearly amplified antisense RNA", MaterialKind.MOLECULAR)
    graph.react("ivt", "Linearly amplify pooled cDNA by T7 transcription"
                + ("; no DNase step; BioAnalyzer check, then SPRI" if v1 else ""),
                StepKind.TRANSCRIPTION,
                ((ds, InputRole.TEMPLATE), (t7, InputRole.CATALYST)), (arna,))

    if protocol == V1:
        fragmentation = graph.input("fragmentation_condition", "aRNA fragmentation condition",
                                    MaterialKind.CONDITION)
        # Supplemental Information "DNA-Library Preparation": SPRI, RNA fragmentation,
        # SPRI, SAP treatment (added before ligation), primer ligation, RT, enrichment PCR.
        fragmented = graph.product("fragmented_arna", "fragmented antisense RNA",
                                   MaterialKind.MOLECULAR, source=boundary)
        graph.react("fragment", "Fragment the amplified antisense RNA, then SPRI",
                    StepKind.FRAGMENTATION,
                    ((arna, InputRole.TEMPLATE), (fragmentation, InputRole.CONDITION)),
                    (fragmented,))
        sap = graph.input("sap", "Shrimp Alkaline Phosphatase", MaterialKind.CATALYST,
                          source=SI_LIBRARY)
        dephos = graph.product("dephosphorylated_arna", "SAP-treated fragmented aRNA",
                               MaterialKind.MOLECULAR, source=SI_LIBRARY)
        graph.react("sap_treatment",
                    "Treat with Shrimp Alkaline Phosphatase for 30 min before primer ligation "
                    "(purpose not stated)", StepKind.DEPHOSPHORYLATION,
                    ((fragmented, InputRole.TEMPLATE), (sap, InputRole.CATALYST)), (dephos,),
                    source=SI_LIBRARY)
        fragmented = dephos
        rlo = graph.input("rlo", "RNA ligation oligo " + V1_RLO_WRITTEN, MaterialKind.ADAPTER,
                          source=boundary, state=Construct([seg("RLO", RLO, "r2")],
                                                           name="RNA ligation oligo"))
        ligase = graph.input("rna_ligase", "RNA ligase", MaterialKind.CATALYST)
        ligated = graph.product("rlo_ligated_arna", "RLO-ligated aRNA fragments",
                                MaterialKind.MOLECULAR, source=boundary)
        graph.react("rlo_ligation", "Ligate the blocked RLO to fragmented aRNA",
                    StepKind.LIGATION,
                    ((fragmented, InputRole.TEMPLATE), (rlo, InputRole.STRUCTURAL),
                     (ligase, InputRole.CATALYST)), (ligated,),
                    geometry=AdapterLigationGeometry.known(
                        "single-strand", rlo_ligation_scene(),
                        statement="Cell 2015 Table S5 prints the RNA ligation oligo with a 5′ phosphate and a 3′ C3 spacer; the phosphate is the ligation donor for the aRNA 3′ end, and the C3 spacer blocks further ligation at the RLO 3′ end (blocking role inferred). The ligase and conditions follow the cited CEL-Seq/MARS-seq protocol and are not printed."))
        conversion_primer = graph.input("second_rt_primer", "2nd RT primer", MaterialKind.OLIGO,
                                        source=boundary,
                                        state=Construct([seg("2nd RT primer", V1_RT2, "r2")],
                                                        name="2nd RT primer"))
        converted = graph.product("second_cdna", "second-round cDNA", MaterialKind.MOLECULAR,
                                  source=boundary)
        graph.react("second_rt", "Reverse-transcribe from the ligated RLO",
                    StepKind.REVERSE_TRANSCRIPTION,
                    ((ligated, InputRole.TEMPLATE), (conversion_primer, InputRole.STRUCTURAL),
                     (rt, InputRole.CATALYST)), (converted,))
        pcr_template = converted
        reviewed_tail = ("fragmentation_condition", "fragment", "sap", "sap_treatment",
                         "rlo", "rna_ligase",
                         "rlo_ligation", "second_rt_primer", "second_rt")
        primers = V1_PRIMERS
        lengths = V1_READ_LENGTHS
        uses = V1_READ_USES
    else:
        pe2n6 = graph.input("pe2_n6", "PE2-N6 random primer", MaterialKind.OLIGO)
        converted = graph.product("sense_cdna", "PE2-tagged sense cDNA", MaterialKind.MOLECULAR)
        graph.react("second_rt", "Random-prime aRNA with PE2-N6",
                    StepKind.REVERSE_TRANSCRIPTION,
                    ((arna, InputRole.TEMPLATE), (pe2n6, InputRole.STRUCTURAL),
                     (rt, InputRole.CATALYST)), (converted,))
        pcr_template = converted
        reviewed_tail = ("pe2_n6", "second_rt")
        primers = V2_PRIMERS
        lengths = {}
        uses = {"Read 1": "include", "Index 1 (i7)": "demultiplex",
                "Read 2": "extract"}
    forward = graph.input("pcr_forward", "P5 library PCR primer", MaterialKind.OLIGO,
                          source=boundary,
                          state=Construct([seg("P5", il.P5, "p5"),
                                           seg("PE1 anneal", V1_PCR1[len(il.P5):], "r1")],
                                          name="PCR enrichment primer 1") if v1 else None)
    reverse = graph.input("pcr_reverse",
                          "P7 library PCR primer (no sample index printed)"
                          if v1 else "indexed P7 library PCR primer",
                          MaterialKind.OLIGO, source=boundary,
                          state=Construct([seg("P7", il.P7, "p7"),
                                           seg("RT2 anneal", V1_PCR2[len(il.P7):], "r2")],
                                          name="PCR enrichment primer 2") if v1 else None)
    pcr_pol = graph.input("pcr_polymerase", "library PCR polymerase", MaterialKind.CATALYST)
    final = library(protocol)
    out = graph.product("library", final.name, MaterialKind.LIBRARY, state=final,
                        source=boundary)
    graph.react("library_pcr", "Add P5 and P7 by enrichment PCR (10–13 cycles)" if protocol == V1 else
                "Add P5, P7 and sample index by enrichment PCR",
                StepKind.PCR,
                ((pcr_template, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (pcr_pol, InputRole.CATALYST)),
                (out,), geometry=PCRGeometry.unavailable(
                    final, reason="The pre-PCR IVT-conversion intermediate is not represented base-for-base."))
    graph.sequence(out, SequencingPlan(
        library=final, primers=tuple(primers),
        acquired_roles=tuple(p.role for p in primers), read_lengths=lengths,
        read_uses=uses,
        limitation=("The v1 library is assembled from the Cell 2015 Table S5 oligos; Table S5"
                    " prints no sequencing primer, so canonical Illumina Read 1 and PE Read 2"
                    " primers are assumed from the library ends. No acquired cycle counts are"
                    " stated in the Cell 2015 article or its Supplemental Information."
                    if protocol == V1 else
                    "Acquired cycle counts are not stated in the inspected supplements.")))
    graph.declare_reviewed_scope((
        "bead_seed", "barcode_plate_oligos", "bead_polymerase", "bead_synthesis",
        "cells", "droplet_loading", "partition", "reverse_transcriptase",
        "reverse_transcription", *v1_extra, "second_strand_polymerase", "second_strand",
        "t7_polymerase", "ivt", *reviewed_tail, "pcr_forward", "pcr_reverse",
        "pcr_polymerase", "library_pcr", "library"), source)
    return graph.require_publishable()

