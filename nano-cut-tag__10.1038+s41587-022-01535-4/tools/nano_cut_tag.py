"""nano-CUT&Tag molecular model from Bartosovic et al. Supplementary Table 1."""
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

TITLE = "nano-CUT&Tag — multimodal nanobody–Tn5 chromatin profiling"
NOTES = "01_nano-cut-tag.html"
SOURCE = '<a href="https://doi.org/10.1038/s41587-022-01535-4">Bartosovic et al., <i>Nature Biotechnology</i> (2023)</a>.'
SUMMARY = ("Nanobody–Tn5 fusions first transfer a target-specific MeA adapter. Chromium "
           "linear amplification adds cell identity; a second MeB-only Tn5 installs the "
           "opposite end before PCR.")
CAVEAT = ("The representative BcdA adapter and both Tn5 duplexes are exact. The 10x "
          "gel-bead oligo and SI-PCR primers are proprietary, so that binding state is "
          "explicitly unavailable rather than reconstructed.")

OUTER = "TCCACGC"
BARCODE_A = "TATAGCCT"
CONNECTOR = "GCGATCGAGGACGGC"
ME_REV = nx.ME_RC
CUSTOM_R1 = CONNECTOR + nx.ME
CUSTOM_I2 = nx.ME_RC + revcomp(CONNECTOR)


def modality_feature():
    return feature("nanoct_modality", "feature_barcode", "whitelist",
                   whitelist="Bartosovic et al. Supplementary Table 1")


def mea_parts():
    return [seg("S5", nx.S5, "s5"), seg("outer connector", OUTER),
            seg("BcdA modality barcode", BARCODE_A, "cbc", feature=modality_feature()),
            seg("inner connector", CONNECTOR), seg("MeA mosaic end", nx.ME, "me")]


def meb_parts():
    return [seg("S7", nx.S7, "s7"), seg("MeB mosaic end", nx.ME, "me")]


def loaded(parts, end, label):
    return loaded_tn5_adapter(parts, mosaic_end=end, label=label,
                              bottom_name="Tn5_Rev", bottom_label="5′-phosphorylated Tn5_Rev",
                              bottom_sequence=ME_REV, mod5="p")


FIRST_TAGGED = Construct([*mea_parts(), seg("targeted genomic DNA", "X" * 38,
                                            placeholder=True)],
                         name="BcdA-tagged target fragment")

PRE_PCR = Construct([
    seg("10x gel-bead amplification handle", "H" * 18, placeholder=True, inferred=True),
    seg("10x cell barcode", "B" * 16, "cbc", placeholder=True,
        feature=feature("cell_barcode", "cell_barcode", "whitelist",
                        whitelist="10x Chromium ATAC whitelist")),
    *mea_parts(), seg("targeted genomic DNA", "X" * 38, placeholder=True),
    seg("MeB mosaic end reverse complement", nx.ME_RC, "me"),
    seg("S7 reverse complement", nx.S7_RC, "s7"),
], name="twice-tagmented nano-CT molecule before proprietary PCR")

FINAL_LIBRARY = Construct([
    seg("P5", il.P5, "p5"),
    seg("10x cell barcode", "B" * 16, "cbc", placeholder=True,
        feature=feature("cell_barcode", "cell_barcode", "whitelist",
                        whitelist="10x Chromium ATAC whitelist")),
    *mea_parts(), seg("targeted genomic DNA", "X" * 38, placeholder=True),
    seg("MeB mosaic end reverse complement", nx.ME_RC, "me"),
    seg("S7 reverse complement", nx.S7_RC, "s7"),
    seg("i7 sample index reverse complement", "I" * 8, "cbc", placeholder=True,
        feature=feature("sample_index_i7", "sample_index", "unknown")),
    seg("P7 reverse complement", il.P7_RC, "p7"),
], name="nano-CUT&Tag BcdA library")

SEQ_PRIMERS = (
    sp.custom("Read 1", "nano-CT custom Read 1", CUSTOM_R1,
              "Bartosovic et al. Supplementary Table 1"),
    sp.NEXTERA["I1"],
    sp.custom("Index 2 (i5)", "nano-CT custom Index 2", CUSTOM_I2,
              "Bartosovic et al. Supplementary Table 1"),
    sp.NEXTERA["R2"],
)
_problems = sp.verify(FINAL_LIBRARY, SEQ_PRIMERS,
                      required_roles=tuple(p.role for p in SEQ_PRIMERS))
if _problems:
    raise ValueError("nano-CT final library: " + "; ".join(_problems))

FINAL_CAPTION = ("Representative BcdA product. Custom R1 traverses MeA into genomic DNA; "
                 "custom Index 2 traverses the modality adapter toward the 10x cell barcode.")
SEQUENCING_INTRO = ("The source specifies a 36-8-48-36 R1-I1-I2-R2 run with custom R1 and "
                    "Index 2 primers; R2 and Index 1 are standard.")
READ_LENGTHS = {"Read 1": 36, "Index 1 (i7)": 8,
                "Index 2 (i5)": 48, "Read 2": 36}
SHOW_READ_CONTENT = True


def oligos():
    yield oligo("Tn5_P5_MeA_BcdA_0N", mea_parts())
    yield oligo("Tn5_P7_MeB", meb_parts())
    yield oligo("Tn5_Rev", [Segment("ME reverse", ME_REV, "me")], mods="/5Phos/")
    yield oligo("R1_seq", [seg("binding region", CUSTOM_R1, "r1")])
    yield oligo("I2_seq", [seg("binding region", CUSTOM_I2, "r2")])


def workflow():
    start = Scene.duplex([seg("antibody-targeted chromatin", "X" * 38,
                              placeholder=True)], label="primary-antibody-bound chromatin")
    w = Workflow(MolecularState("Antibody-bound chromatin", tuple(start.rows())))
    w.react("Anneal BcdA MeA to Tn5_Rev and load nanobody–Tn5",
            loaded(mea_parts(), "MeA mosaic end", "BcdA nano-Tn5").rows(),
            name="target-coded nano-Tn5")
    w.react("Tagment antibody-proximal DNA and fill the staggered gap",
            Scene.duplex(list(FIRST_TAGGED), label="BcdA-tagged genomic template").rows(),
            name="BcdA-tagged genomic DNA")
    linear = Construct([seg("10x gel-bead amplification handle", "H" * 18,
                            placeholder=True, inferred=True),
                        seg("10x cell barcode", "B" * 16, "cbc", placeholder=True,
                            feature=feature("cell_barcode", "cell_barcode", "whitelist",
                                            whitelist="10x Chromium ATAC whitelist")),
                        *list(FIRST_TAGGED)], name="cell-barcoded linear-amplification product")
    w.react("INFERRED — linearly amplify in Chromium droplets",
            Scene.duplex(list(linear), label="cell-barcoded copies").rows(),
            name="cell-barcoded linear products",
            note="The paper establishes the product order, but the complete gel-bead oligo is not public.")
    w.react("Anneal MeB to Tn5_Rev and load standard Tn5",
            loaded(meb_parts(), "MeB mosaic end", "MeB-only Tn5").rows(),
            name="MeB-only transposome")
    w.react("Tagment the linear products with MeB-only Tn5",
            Scene.duplex(list(PRE_PCR), label="twice-tagmented nano-CT template").rows(),
            name="PCR-ready nano-CT molecule")
    geometry = PCRGeometry.unavailable(
        PRE_PCR,
        reason="The source names 10x SI-PCR primer B and Single Index N Set A, but does not print their complete sequences.")
    w.pcr("Amplify with the 10x scATAC SI-PCR reagents", geometry,
          Scene.duplex(list(FINAL_LIBRARY), label="nano-CT sequencing library").rows(),
          name="nano-CT sequencing library")
    return w
