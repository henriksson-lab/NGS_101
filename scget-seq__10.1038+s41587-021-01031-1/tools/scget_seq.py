"""scGET-seq molecular model from Tedesco et al. Supplementary Table 10."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import illumina as il
import nextera as nx
import seqprimers as sp
from batch_ngs import seg
from chemdraw import Construct, MolecularState, Scene, Segment, Workflow, feature, oligo
from pcr import PCRGeometry
from protocol_graph import (Disposition, Terminal, InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, StepKind)
from tagmentation import loaded_tn5_adapter

TITLE = "scGET-seq — open and heterochromatic DNA in the same cell"
NOTES = "01_scget-seq.html"
SOURCE = '<a href="https://doi.org/10.1038/s41587-021-01031-1">Tedesco et al., <i>Nature Biotechnology</i> (2022)</a>.'
SUMMARY = ("Ordinary Tn5 and HP1α–Tn5 transfer different inline adapter barcodes; "
           "Chromium ATAC subsequently adds a shared cell identity.")
CAVEAT = ("The representative Tn5ME-A.1 and TnHME-A.4 top oligos are exact. Their common "
          "mosaic-end reverse strand follows standard Tn5 loading. Complete 10x gel-bead "
          "and SI-PCR oligos remain unassigned in this model; vendor originals are deferred, so this is not a confirmed nondisclosure claim.")

CONNECTOR = "TCCGATCT"
TN5_BARCODE = "CGTACTAG"
TNH_BARCODE = "TAAGGCGA"
TN5MEA1 = nx.S5 + CONNECTOR + TN5_BARCODE + nx.ME
TNHMEA4 = nx.S5 + CONNECTOR + TNH_BARCODE + nx.ME
ME_REV = nx.ME_RC


def modality_feature():
    return feature("scget_modality", "feature_barcode", "whitelist",
                   whitelist="Tedesco et al. Supplementary Table 10")


def adapter_parts(kind, barcode):
    return [seg("S5", nx.S5, "s5"), seg("adapter connector", CONNECTOR),
            seg(kind + " modality barcode", barcode, "cbc", feature=modality_feature()),
            seg(kind + " mosaic end", nx.ME, "me")]


def loaded(kind, barcode):
    return loaded_tn5_adapter(
        adapter_parts(kind, barcode), mosaic_end=kind + " mosaic end",
        label=kind + "ME-A top oligo", bottom_name="ME reverse",
        bottom_label="5′-phosphorylated ME reverse", bottom_sequence=ME_REV, mod5="p")


def tagged(kind, barcode, dna):
    return Construct([*adapter_parts(kind, barcode),
                      seg(dna, "X" * 38, placeholder=True)],
                     name=kind + "-tagged genomic fragment")


TN5_TAGGED = tagged("Tn5", TN5_BARCODE, "accessible genomic DNA")
TNH_TAGGED = tagged("TnH", TNH_BARCODE, "H3K9me3-associated genomic DNA")

PRE_PCR = Construct([
    seg("10x gel-bead amplification handle", "H" * 18, placeholder=True, inferred=True),
    seg("10x cell barcode", "B" * 16, "cbc", placeholder=True,
        feature=feature("cell_barcode", "cell_barcode", "whitelist",
                        whitelist="10x Chromium ATAC whitelist")),
    *list(TNH_TAGGED),
    seg("opposite mosaic end reverse complement", nx.ME_RC, "me"),
    seg("S7 reverse complement", nx.S7_RC, "s7"),
], name="cell-barcoded scGET-seq molecule before proprietary PCR")

FINAL_LIBRARY = Construct([
    seg("P5", il.P5, "p5"),
    seg("10x cell barcode", "B" * 16, "cbc", placeholder=True,
        feature=feature("cell_barcode", "cell_barcode", "whitelist",
                        whitelist="10x Chromium ATAC whitelist")),
    *list(TNH_TAGGED),
    seg("opposite mosaic end reverse complement", nx.ME_RC, "me"),
    seg("S7 reverse complement", nx.S7_RC, "s7"),
    seg("i7 sample index reverse complement", "I" * 8, "cbc", placeholder=True,
        feature=feature("sample_index_i7", "sample_index", "unknown")),
    seg("P7 reverse complement", il.P7_RC, "p7"),
], name="representative TnH scGET-seq library")

SEQ_PRIMERS = (
    sp.custom("Read 1", "S5-side genomic read", nx.S5,
              "scGET-seq adapter and 10x scATAC structure"),
    sp.NEXTERA["I1"],
    sp.NEXTERA["R2"],
)
_problems = sp.verify(FINAL_LIBRARY, SEQ_PRIMERS,
                      required_roles=tuple(p.role for p in SEQ_PRIMERS))
if _problems:
    raise ValueError("scGET-seq final library: " + "; ".join(_problems))

FINAL_CAPTION = ("Representative TnHME-A.4 product. The inline barcode distinguishes "
                 "HP1α-targeted TnH from ordinary Tn5; the 10x barcode identifies the cell.")
SEQUENCING_INTRO = ("Paired genomic reads identify the insertion locus; the inline adapter "
                    "barcode assigns chromatin modality and the 10x barcode assigns cell.")
SEQUENCING_UNAVAILABLE = ("The complete dedicated 10x cell-barcode sequencing-primer oligo "
                          "is vendor-supplied; genomic Read 1, Index 1 and Read 2 sites are shown.")


def oligos():
    yield oligo("Tn5ME-A.1", adapter_parts("Tn5", TN5_BARCODE))
    yield oligo("TnHME-A.4", adapter_parts("TnH", TNH_BARCODE))
    yield oligo("ME reverse", [Segment("ME reverse", ME_REV, "me")], mods="/5Phos/")


def workflow():
    start = Scene.duplex([seg("nuclear chromatin", "X" * 38, placeholder=True)],
                         label="permeabilized nuclear chromatin")
    w = Workflow(MolecularState("Permeabilized nuclei", tuple(start.rows())))
    w.react("Anneal Tn5ME-A.1 and TnHME-A.4 separately to ME reverse and load Tn5/TnH",
            [*loaded("Tn5", TN5_BARCODE).rows(), *loaded("TnH", TNH_BARCODE).rows()],
            name="two modality-coded transposomes")
    w.react("Tagment accessible DNA with Tn5 and H3K9me3-associated DNA with HP1α–Tn5",
            [*Scene.duplex(list(TN5_TAGGED), label="accessible-DNA product").rows(),
             *Scene.duplex(list(TNH_TAGGED), label="heterochromatin product").rows()],
            name="modality-coded genomic fragments")
    w.react("INFERRED — linearly amplify and barcode both product classes in Chromium ATAC droplets",
            Scene.duplex(list(PRE_PCR), label="representative cell-barcoded TnH product").rows(),
            name="cell-barcoded linear products",
            note="The source establishes the 10x workflow and molecule roles, while the complete gel-bead oligo remains unassigned.")
    geometry = PCRGeometry.unavailable(
        PRE_PCR,
        reason="The source follows the 10x scATAC kit but the complete SI-PCR primer assignment remains unchecked.")
    w.pcr("PCR-complete the Chromium ATAC library", geometry,
          Scene.duplex(list(FINAL_LIBRARY), label="scGET-seq sequencing library").rows(),
          name="scGET-seq sequencing library")
    return w


def protocol_graph():
    source = "Defining source: https://doi.org/10.1038/s41587-021-01031-1"
    proprietary = source + " (10x gel-bead and SI-PCR oligos remain unassigned pending vendor originals)"
    graph = ProtocolGraph("scGET-seq", authority=source)
    nuclei = graph.input("nuclei", "permeabilized nuclear chromatin",
                         MaterialKind.BIOLOGICAL, primary=True)
    adapters = graph.input("modality_adapters", "Tn5ME-A.1 and TnHME-A.4 loaded adapters",
                           MaterialKind.ADAPTER)
    transposases = graph.input("transposases", "Tn5 and HP1α–Tn5 transposases",
                               MaterialKind.CATALYST)
    tagged_products = graph.product("modality_tagged_fragments",
                                    "accessible and H3K9me3 modality-coded fragments",
                                    MaterialKind.MOLECULAR)
    tagmentation_waste = graph.product("tagmentation_waste",
                                       "candidate non-library tagmentation products", MaterialKind.MOLECULAR)
    graph.react("dual_tagmentation",
                "Tagment accessible DNA with Tn5 and H3K9me3-associated DNA with HP1α–Tn5",
                StepKind.TAGMENTATION,
                ((nuclei, InputRole.TEMPLATE), (adapters, InputRole.STRUCTURAL),
                 (transposases, InputRole.CATALYST)), (tagged_products, tagmentation_waste))
    graph.mark_terminal(Terminal(tagmentation_waste.id, Disposition.UNRESOLVED, "Candidate tagmentation product classes and their physical recovery are not resolved by the adapter table."))
    gel_bead = graph.input("gel_bead_oligo", "10x Chromium ATAC gel-bead oligo",
                           MaterialKind.ADAPTER, inferred=True, source=proprietary)
    cell_barcoded = graph.product("cell_barcoded_products",
                                  "cell-barcoded linear scGET-seq products",
                                  MaterialKind.MOLECULAR, state=PRE_PCR,
                                  inferred=True, source=proprietary)
    graph.react("droplet_barcoding",
                "Linearly amplify and barcode both product classes in Chromium ATAC droplets",
                StepKind.ADAPTER_ATTACHMENT,
                ((tagged_products, InputRole.TEMPLATE), (gel_bead, InputRole.STRUCTURAL)),
                (cell_barcoded,), source=proprietary)
    i5 = graph.input("si_p5_primer", "Chromium ATAC SI-PCR P5 primer", MaterialKind.OLIGO,
                     inferred=True, source=proprietary)
    i7 = graph.input("si_p7_primer", "Chromium ATAC indexed SI-PCR P7 primer",
                     MaterialKind.OLIGO, inferred=True, source=proprietary)
    polymerase = graph.input("pcr_polymerase", "SI-PCR polymerase", MaterialKind.CATALYST)
    library = graph.product("library", FINAL_LIBRARY.name, MaterialKind.LIBRARY,
                            state=FINAL_LIBRARY, inferred=True, source=proprietary)
    graph.react("library_pcr", "PCR-complete the Chromium ATAC library", StepKind.PCR,
                ((cell_barcoded, InputRole.TEMPLATE), (i5, InputRole.STRUCTURAL),
                 (i7, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)), (library,),
                geometry=PCRGeometry.unavailable(
                    PRE_PCR, reason="The source follows the 10x scATAC kit but complete SI-PCR primer assignment remains unchecked."),
                source=proprietary)
    graph.sequence(library, SequencingPlan(
        library=FINAL_LIBRARY, primers=SEQ_PRIMERS,
        acquired_roles=tuple(p.role for p in SEQ_PRIMERS),
        read_uses={"Read 1": "include", "Index 1 (i7)": "demultiplex",
                   "Read 2": "include"}, limitation=SEQUENCING_UNAVAILABLE))
    graph.declare_reviewed_scope((
        "nuclei", "modality_adapters", "transposases", "dual_tagmentation",
        "gel_bead_oligo", "droplet_barcoding", "si_p5_primer", "si_p7_primer",
        "pcr_polymerase", "library_pcr", "library"), source)
    return graph.require_publishable()
