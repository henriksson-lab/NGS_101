"""ASAP-seq molecular model."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from antibody_tags import bridge_extension_scene, feature_library
from batch_ngs import nextera_library, seg
from chemdraw import Row, feature
from multimodal_spatial import modality_split

TITLE = "ASAP-seq — ATAC and surface-protein profiling"
NOTES = "01_asap-seq.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1038/s41587-021-00927-2">Mimitou et al., <i>Nature Biotechnology</i> (2021)</a>.'
SUMMARY = "Antibody-derived tags are extended with an ATAC-bead capture handle before a 10x single-cell ATAC workflow yields separate chromatin-accessibility and protein-tag libraries."
CAVEAT = "TotalSeq tag and bridge sequences vary by reagent generation. Their published molecular roles are shown without asserting one product's proprietary bases."
CELL=feature("asap_cell","cell_barcode","whitelist",whitelist="10x Chromium ATAC bead barcode")
ATAC, ATAC_P = nextera_library([seg("cell barcode association","C"*16,"cbc",placeholder=True,feature=CELL),seg("accessible genomic DNA", "X" * 42, placeholder=True)], "ASAP-seq ATAC library")
ADT, ADT_P = feature_library("ASAP-seq antibody-tag library", barcode_id="asap_antibody",cell_id="asap_cell")
FINAL_LIBRARIES = (
    ("ATAC library", ATAC, ATAC_P, "10x scATAC bead barcoding and indexed PCR complete the tagmented chromatin library.", "Paired genomic reads use the standard Nextera primer sites."),
    ("Antibody-tag library", ADT, ADT_P, "The extended TotalSeq tag is amplified as a separate Illumina library.", "Read cycles recover antibody barcode, UMI and sample indexes."),
)
READ_LENGTHS = {"ATAC library": {"Read 1": 50, "Read 2": 50}, "Antibody-tag library": {"Read 1": 28, "Read 2": 50}}

def sections():
    return [
        ("Stain fixed, permeabilized cells with DNA-barcoded antibodies", [Row(chunks=[("antibody — feature barcode — UMI — bridge tail", "cbc", False)])], "Surface or intracellular epitopes retain their antibody-derived DNA tags."),
        ("Extend antibody tags across a blocked bridge", bridge_extension_scene().rows(), "The bridge supplies sequence complementary to the scATAC bead oligo; its 3-prime block prevents bridge self-extension."),
        ("Tagment accessible chromatin and partition nuclei", [Row(chunks=[("Tn5 → accessible chromatin fragments; GEM bead → shared cell identity", "me", False)])], "Chromatin and antibody tags from the same cell receive the same droplet barcode."),
        ("Separate and amplify modalities", modality_split("ATAC library", "antibody-tag library"), "Distinct PCRs recover genomic inserts and antibody-derived tags."),
    ]
