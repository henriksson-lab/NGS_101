"""NTT-seq molecular model from Stuart et al. Extended Data Table 2."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import illumina as il
import nextera as nx
import seqprimers as sp
from batch_ngs import seg
from chemdraw import Construct, MolecularState, Scene, Segment, Workflow, feature, oligo, revcomp
from pcr import PCRGeometry
from tagmentation import loaded_tn5_adapter

TITLE = "NTT-seq — nanobody-tethered transposition with single-cell indexing"
NOTES = "01_ntt-seq.html"
SOURCE = '<a href="https://doi.org/10.1038/s41587-022-01588-5">Stuart et al., <i>Nature Biotechnology</i> (2023)</a>.'
SUMMARY = ("Species-selective nanobody–Tn5 fusions transfer a target-specific MEDSA "
           "barcode and a MEDSB opposite end before Chromium ATAC adds cell identity.")
CAVEAT = ("MEDSA_1/MEDSB and their Tn5 loading geometry are exact. The complete 10x "
          "gel-bead and SI-PCR oligos are vendor-supplied and remain explicit unavailable regions.")

MEDSA_BARCODES = ("GGATTGCT", "GTAATGCA", "GTCAAGGA", "GTGAGCGT", "GTGTGACC", "TAAGGTGG")
CONNECTOR = "GCGATCGAGGACGGC"
MEDSA1 = nx.S5 + MEDSA_BARCODES[0] + CONNECTOR + nx.ME
MEDSB = nx.S7 + nx.ME
ME_REV = nx.ME_RC
CUSTOM_R1 = CONNECTOR + nx.ME
CUSTOM_I5 = nx.ME_RC + revcomp(CONNECTOR)


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
    seg("10x gel-bead amplification handle", "H" * 18, placeholder=True, inferred=True),
    seg("10x cell barcode", "B" * 16, "cbc", placeholder=True,
        feature=feature("cell_barcode", "cell_barcode", "whitelist",
                        whitelist="10x Chromium ATAC whitelist")),
    *list(TAGGED),
], name="cell-barcoded NTT-seq molecule before proprietary PCR")

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


def workflow():
    start = Scene.duplex([seg("antibody-targeted chromatin", "X" * 38,
                              placeholder=True)], label="primary-antibody-bound chromatin")
    w = Workflow(MolecularState("Antibody-bound chromatin", tuple(start.rows())))
    w.react("Anneal MEDSA_1 and MEDSB separately to ME reverse and load nanobody–Tn5",
            [*loaded(medsa_parts(), "MEDSA mosaic end", "MEDSA_1 nanobody–Tn5").rows(),
             *loaded(medsb_parts(), "MEDSB mosaic end", "MEDSB nanobody–Tn5").rows()],
            name="MEDSA_1/MEDSB nanobody transposomes")
    w.react("Tagment target-proximal DNA and fill the staggered gaps",
            Scene.duplex(list(TAGGED), label="MEDSA_1/MEDSB target-coded fragment").rows(),
            name="target-coded chromatin fragments")
    w.react("INFERRED — linearly amplify and barcode fragments in Chromium ATAC droplets",
            Scene.duplex(list(PRE_PCR), label="cell-barcoded NTT-seq copies").rows(),
            name="cell-barcoded linear products",
            note="The paper specifies the 10x scATAC workflow, but the complete gel-bead oligo is not public.")
    geometry = PCRGeometry.unavailable(
        PRE_PCR,
        reason="The source follows the 10x scATAC kit but does not print the complete SI-PCR primer sequences.")
    w.pcr("PCR-complete the Chromium ATAC library", geometry,
          Scene.duplex(list(FINAL_LIBRARY), label="NTT-seq sequencing library").rows(),
          name="NTT-seq sequencing library")
    return w
