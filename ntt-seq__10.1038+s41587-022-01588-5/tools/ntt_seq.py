"""NTT-seq molecular model from Stuart et al. Extended Data Table 2."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import illumina as il
import nextera as nx
import seqprimers as sp
from batch_ngs import seg
from chemdraw import Construct, MolecularState, Scene, Segment, Workflow, feature, oligo, revcomp
from pcr import PCRGeometry, PCRPrimer
from tagmentation import loaded_tn5_adapter
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from protocol_model import molecule, role

TITLE = "NTT-seq — nanobody-tethered transposition with single-cell indexing"
NOTES = "01_ntt-seq.html"
SOURCE = '<a href="https://doi.org/10.1038/s41587-022-01588-5">Stuart et al., <i>Nature Biotechnology</i> (2023)</a>.'
SUMMARY = ("Species-selective nanobody–Tn5 fusions transfer a target-specific MEDSA "
           "barcode and a MEDSB opposite end before Chromium ATAC adds cell identity.")
CAVEAT = ("The target-specific adapters and Tn5 loading geometry are source-defined. "
          "The linked 10x ATAC v1.1 guide (Rev G) prints gel-bead and SI-PCR oligo "
          "scaffolds, drawn here at their binding sites; the 16-nt cell barcode and 8-nt sample index vary. "
          "This drawing follows the chromatin library. The published optional TotalSeq-A "
          "surface-protein branch adds a blocked bridge oligo and makes a separate antibody-tag library "
          "(Methods lines 110–112 and 128–136); its tag products are not represented by this molecule. "
          "The transferred intermediate is drawn closed. The single-cell Methods (JATS line 130) load "
          "the tagmented cells directly into the 10x workflow, so gap filling is assigned to the 72 °C "
          "5-min opening step of the GEM incubation; that assignment is inferred from the 10x program "
          "and the bulk PCR's own '72°C for 5 min (gap filling)' label (line 126).")

MEDSA_BARCODES = ("GGATTGCT", "GTAATGCA", "GTCAAGGA", "GTGAGCGT", "GTGTGACC", "TAAGGTGG")
CONNECTOR = "GCGATCGAGGACGGC"
MEDSA1 = nx.S5 + MEDSA_BARCODES[0] + CONNECTOR + nx.ME
MEDSB = nx.S7 + nx.ME
ME_REV = nx.ME_RC
CUSTOM_R1 = CONNECTOR + nx.ME
CUSTOM_I5 = nx.ME_RC + revcomp(CONNECTOR)


ATAC_GUIDE = "10x Chromium Next GEM Single Cell ATAC v1.1 guide CG000209 Rev G, Appendix A1"
GEL_BEAD_FORWARD = il.P5 + "N" * 16 + nx.S5
SI_PCR_FORWARD = "AATGATACGGCGACCACCGAGA"
SI_PCR_REVERSE = "CAAGCAGAAGACGGCATACGAGAT" + "N" * 8 + "GTCTCGTGGGCTCGG"


def si_pcr_geometry():
    """CG000209 Rev G Appendix A1 SI-PCR Primer B (partial P5) and P7-N8-partial Read2N."""
    return PCRGeometry.exact(
        PRE_PCR,
        PCRPrimer("SI-PCR forward (partial P5)", SI_PCR_FORWARD, anneal_nt=len(SI_PCR_FORWARD)),
        PCRPrimer("SI-PCR reverse (P7-i7-partial Read2N)", SI_PCR_REVERSE, anneal_nt=len(nx.S7)),
        statement="10x CG000209 Rev G Appendix A1 prints the SI-PCR forward and indexed reverse scaffolds; "
                  "the reverse primer anneals at the MEDSB S7 entry site. Rev G is not shown to be the historical kit revision NTT used.")


def modality_feature():
    return feature("ntt_modality", "feature_barcode", "whitelist",
                   whitelist="Stuart et al. Extended Data Table 2")


def medsa_parts():
    return [seg("S5", nx.S5, "s5"),
            seg("MEDSA_1 modality barcode", MEDSA_BARCODES[0], "cbc",
                feature=modality_feature()),
            seg("connector", CONNECTOR), seg("MEDSA mosaic end", nx.ME, "me")]


def medsb_parts():
    return [seg("S7", nx.S7, "s7"), seg("MEDSB mosaic end", nx.ME, "me")]


def loaded(parts, end, label):
    return loaded_tn5_adapter(parts, mosaic_end=end, label=label,
                              bottom_name="ME reverse", bottom_label="5′-phosphorylated ME reverse",
                              bottom_sequence=ME_REV, mod5="p")


TAGGED = Construct([
    *medsa_parts(), seg("targeted chromatin fragment", "X" * 38, placeholder=True),
    seg("MEDSB mosaic end reverse complement", nx.ME_RC, "me"),
    seg("S7 reverse complement", nx.S7_RC, "s7"),
], name="MEDSA_1/MEDSB target-coded fragment")

PRE_PCR = Construct([
    seg("10x gel-bead P5", il.P5, "p5"),
    seg("10x cell barcode", "B" * 16, "cbc", placeholder=True,
        feature=feature("cell_barcode", "cell_barcode", "whitelist",
                        whitelist="10x Chromium ATAC whitelist")),
    *list(TAGGED),
], name="cell-barcoded NTT-seq molecule before SI-PCR")

FINAL_LIBRARY = Construct([
    seg("P5", il.P5, "p5"),
    seg("10x cell barcode", "B" * 16, "cbc", placeholder=True,
        feature=feature("cell_barcode", "cell_barcode", "whitelist",
                        whitelist="10x Chromium ATAC whitelist")),
    *list(TAGGED),
    seg("i7 sample index reverse complement", "I" * 8, "cbc", placeholder=True,
        feature=feature("sample_index_i7", "sample_index", "unknown")),
    seg("P7 reverse complement", il.P7_RC, "p7"),
], name="NTT-seq MEDSA_1 library")

SEQ_PRIMERS = (
    sp.custom("Read 1", "NTT custom Read 1", CUSTOM_R1,
              "Stuart et al. Extended Data Table 2"),
    sp.NEXTERA["I1"],
    sp.custom("Index 2 (i5)", "NTT custom i5", CUSTOM_I5,
              "Stuart et al. Extended Data Table 2"),
    sp.NEXTERA["R2"],
)
_problems = sp.verify(FINAL_LIBRARY, SEQ_PRIMERS,
                      required_roles=tuple(p.role for p in SEQ_PRIMERS))
if _problems:
    raise ValueError("NTT-seq final library: " + "; ".join(_problems))

FINAL_CAPTION = ("Representative MEDSA_1/MEDSB molecule. The MEDSA barcode identifies "
                 "the nanobody target and the 10x barcode identifies the cell.")
SEQUENCING_INTRO = ("Custom Read 1 enters through the MEDSA connector/ME and custom i5 "
                    "reads back through that adapter toward the 10x cell barcode.")
READ_LENGTHS = {"Read 1": 60, "Index 1 (i7)": 8,
                "Index 2 (i5)": 38, "Read 2": 60}
SHOW_READ_CONTENT = True


def oligos():
    yield oligo("MEDSA_1", medsa_parts())
    yield oligo("MEDSB", medsb_parts())
    yield oligo("ME reverse", [Segment("ME reverse", ME_REV, "me")], mods="/5Phos/")
    yield oligo("Custom R1", [seg("binding region", CUSTOM_R1, "r1")])
    yield oligo("Custom i5", [seg("binding region", CUSTOM_I5, "r2")])

    yield oligo("10x gel-bead oligo (CG000209 A1)", [
        seg("P5", il.P5, "p5"),
        seg("10x cell barcode", "N" * 16, "cbc", placeholder=True,
            feature=feature("cell_barcode", "cell_barcode", "whitelist",
                            whitelist="10x Chromium ATAC whitelist")),
        seg("partial Read1N", nx.S5, "s5")])
    yield oligo("10x SI-PCR forward (CG000209 A1)",
                [seg("partial P5", SI_PCR_FORWARD, "p5")])
    yield oligo("10x indexed PCR reverse (CG000209 A1)", [
        seg("P7", il.P7, "p7"),
        seg("sample index", "N" * 8, "cbc", placeholder=True,
            feature=feature("sample_index_i7", "sample_index", "unknown")),
        seg("partial Read2N", "GTCTCGTGGGCTCGG", "r2")])


def workflow():
    start = Scene.duplex([seg("antibody-targeted chromatin", "X" * 38,
                              placeholder=True)], label="primary-antibody-bound chromatin")
    w = Workflow(MolecularState("Antibody-bound chromatin", tuple(start.rows())))
    w.react("Anneal MEDSA_1 and MEDSB separately to pMENT, mix 1:1 and load nanobody–Tn5 (1 h RT)",
            [*loaded(medsa_parts(), "MEDSA mosaic end", "MEDSA_1 nanobody–Tn5").rows(),
             *loaded(medsb_parts(), "MEDSB mosaic end", "MEDSB nanobody–Tn5").rows()],
            name="MEDSA_1/MEDSB nanobody transposomes")
    w.react("Tagment target-proximal DNA in 10 mM MgCl2, 300 mM NaCl, 37 °C 1 h",
            Scene.duplex(list(TAGGED), label="MEDSA_1/MEDSB target-coded fragment").rows(),
            name="target-coded chromatin fragments")
    w.react("GEM incubation: 72 °C 5 min (inferred gap fill), then 12 linear-amplification cycles",
            Scene.duplex(list(PRE_PCR), label="cell-barcoded NTT-seq copies").rows(),
            name="cell-barcoded linear products",
            note="The 10x ATAC v1.1 guide Appendix A1 prints P5–N16–partial Read1N for the gel-bead oligo.")
    geometry = si_pcr_geometry()
    w.pcr("PCR-complete the Chromium ATAC library", geometry,
          Scene.duplex(list(FINAL_LIBRARY), label="NTT-seq sequencing library").rows(),
          name="NTT-seq sequencing library")
    return w


def protocol_graph():
    graph = ProtocolGraph(TITLE, authority=SOURCE)
    nuclei = graph.input(
        "nuclei", "permeabilized nuclei", MaterialKind.BIOLOGICAL,
        state=molecule("permeabilized nuclei", (role("chromatin epitopes"),)),
        primary=True)
    antibody = graph.input("primary_antibodies", "species-resolved primary antibodies",
                           MaterialKind.SELECTION)
    loaded_adapters = graph.input(
        "loaded_medsa_medsb", "MEDSA/MEDSB-loaded nanobody–Tn5 transposomes",
        MaterialKind.ADAPTER,
        state=(loaded(medsa_parts(), "MEDSA mosaic end", "MEDSA_1 nanobody–Tn5"),
               loaded(medsb_parts(), "MEDSB mosaic end", "MEDSB nanobody–Tn5")))
    nanobody_tn5 = graph.input("nanobody_tn5", "species-selective nanobody–Tn5 fusions",
                               MaterialKind.CATALYST)
    targeted = graph.product("targeted_nuclei", "antibody-bound nanobody–Tn5 nuclei",
                             MaterialKind.MOLECULAR, state=nuclei.state)
    wash = graph.product("unbound_targeting_reagents", "unbound antibodies and transposomes",
                         MaterialKind.WASTE)
    graph.react("target_chromatin", "Primary antibodies overnight 4 °C, wash ×2; nb-Tn5 (2.5 µl per target) in 300 mM NaCl 1 h RT, high-salt wash ×2",
                StepKind.AFFINITY_CAPTURE,
                ((nuclei, InputRole.TEMPLATE), (antibody, InputRole.SELECTION),
                 (loaded_adapters, InputRole.STRUCTURAL),
                 (nanobody_tn5, InputRole.CATALYST)), (targeted, wash))
    graph.discard(wash, "is removed before target-proximal tagmentation")
    tagged = graph.product("target_coded_fragments", TAGGED.name,
                           MaterialKind.MOLECULAR, state=TAGGED)
    tn5_waste = graph.product("tn5_waste", "non-amplifiable or mistargeted Tn5 products",
                              MaterialKind.WASTE)
    mg_buffer = graph.input("tagmentation_buffer",
                            "tagmentation buffer (10 mM MgCl2, 300 mM NaCl), 37 °C 1 h",
                            MaterialKind.CONDITION)
    graph.react("tagment", "Transfer MEDSA and MEDSB at target-proximal DNA (Mg2+, 37 °C 1 h)",
                StepKind.TAGMENTATION,
                ((targeted, InputRole.TEMPLATE), (loaded_adapters, InputRole.STRUCTURAL),
                 (nanobody_tn5, InputRole.CATALYST), (mg_buffer, InputRole.CONDITION)),
                (tagged, tn5_waste))
    graph.discard(tn5_waste, "does not retain the productive MEDSA/MEDSB configuration")
    gel_bead = graph.input(
        "gel_bead_oligo", "10x scATAC gel-bead amplification oligo",
        MaterialKind.OLIGO, state=Construct([seg("P5", il.P5, "p5"),
            seg("10x cell barcode", "N" * 16, "cbc", placeholder=True,
                feature=feature("cell_barcode", "cell_barcode", "whitelist",
                                whitelist="10x Chromium ATAC whitelist")),
            seg("partial Read1N", nx.S5, "s5")]), source=ATAC_GUIDE)
    linear_polymerase = graph.input(
        "linear_polymerase", "Chromium barcoding-mix polymerase (gap fill and linear amplification)",
        MaterialKind.CATALYST)
    cell_barcoded = graph.product(
        "cell_barcoded", PRE_PCR.name, MaterialKind.MOLECULAR, state=PRE_PCR,
        source=SOURCE + "; " + ATAC_GUIDE)
    graph.react("linear_amplification",
                "GEM incubation: 72 °C 5 min (inferred gap fill), 98 °C 30 s, 12 cycles linear amplification with gel-bead oligo",
                StepKind.EXTENSION,
                ((tagged, InputRole.TEMPLATE), (gel_bead, InputRole.STRUCTURAL),
                 (linear_polymerase, InputRole.CATALYST)), (cell_barcoded,))
    forward = graph.input("si_pcr_forward", "10x scATAC SI-PCR forward primer",
                          MaterialKind.OLIGO, source=ATAC_GUIDE)
    reverse = graph.input("si_pcr_reverse", "10x scATAC SI-PCR reverse primer",
                          MaterialKind.OLIGO, source=ATAC_GUIDE)
    polymerase = graph.input("pcr_polymerase", "SI-PCR polymerase", MaterialKind.CATALYST)
    library = graph.product("sequencing_library", FINAL_LIBRARY.name,
                            MaterialKind.LIBRARY, state=FINAL_LIBRARY)
    graph.react("library_pcr", "PCR-complete the Chromium ATAC library", StepKind.PCR,
                ((cell_barcoded, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (library,), geometry=si_pcr_geometry())
    graph.sequence(library, SequencingPlan(
        library=FINAL_LIBRARY, primers=tuple(SEQ_PRIMERS),
        acquired_roles=tuple(p.role for p in SEQ_PRIMERS), read_lengths=READ_LENGTHS))
    graph.declare_reviewed_scope((
        "nuclei", "primary_antibodies", "loaded_medsa_medsb", "nanobody_tn5",
        "target_chromatin", "tagmentation_buffer", "tagment",
        "gel_bead_oligo", "linear_polymerase", "linear_amplification",
        "si_pcr_forward", "si_pcr_reverse", "pcr_polymerase", "library_pcr",
        "sequencing_library"), SOURCE)
    return graph.require_publishable()
