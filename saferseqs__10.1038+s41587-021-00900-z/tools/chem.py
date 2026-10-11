"""SaferSeqS copied-duplex-UID and strand-specific target library model."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from batch_ngs import duplex_fragment, seg
from chemdraw import Construct, MolecularState, Workflow, feature
from adapter_ligation import AdapterLigationGeometry
from molecular_consensus import CopiedDuplexUID

TITLE = "SaferSeqS — copied duplex UID with strand-specific enrichment"
NOTES = "01_saferseqs.html"
SOURCE = ('Defining source: <a href="https://doi.org/10.1038/s41587-021-00900-z">'
          'Cohen et al., <i>Nature Biotechnology</i> (2021)</a>.')
SUMMARY = ("A 3′-adapter UID is copied into both strands of each original duplex before "
           "amplification. Separate Watson- and Crick-specific nested PCRs preserve strand "
           "identity, allowing paired family consensuses from the same molecule.")
CAVEAT = ("Target-specific and sample-index sequences vary by panel. The final page models "
          "their published roles and the fixed 14-nt UID; exact panel oligos remain in the "
          "supplementary source table.")

UID = feature("saferseqs_duplex_uid", "umi", "random", group="duplex_uid", part="shared")
FINAL_LIBRARY = Construct([
    seg("P5 graft and custom read arm", "[P5 + custom Read 1 arm]", placeholder=True),
    seg("i5 sample index", "J" * 8, "cbc", placeholder=True,
        feature=feature("sample_index_i5", "sample_index", "unknown")),
    seg("14-nt copied duplex UID", "N" * 14, "umi", placeholder=True, feature=UID),
    seg("strand-specific target amplicon", "X" * 42, placeholder=True),
    seg("i7 sample index", "I" * 8, "cbc", placeholder=True,
        feature=feature("sample_index_i7", "sample_index", "unknown")),
    seg("custom Read 2 arm and P7 graft", "[custom Read 2 arm + P7]", placeholder=True),
], name="SaferSeqS enriched strand library")
SEQ_PRIMERS = ()
FINAL_CAPTION = ("The first 14 bases of Read 1 identify Watson products; the corresponding "
                 "Crick orientation exposes the same UID in Read 2. Sample identity is in "
                 "the dual index reads.")
SEQUENCING_INTRO = ""
SEQUENCING_ENDING = ("Paired-end sequencing uses custom Read 1, index and Read 2 primers "
                     "mixed with the standard Illumina primers. Watson products expose the "
                     "14-nt UID at the start of Read 1; the reversed Crick products expose it "
                     "at the start of Read 2. Dual eight-base index reads identify samples.")
SEQUENCING_UNAVAILABLE = ("Exact custom sequencing-primer oligos are published in "
                          "Supplementary Table 7, which the automated source fetch records "
                          "as a gated manual download; no bases are inferred here.")


def workflow():
    chemistry = CopiedDuplexUID(14)
    w = Workflow(MolecularState("End-repaired genomic fragment",
                                tuple(duplex_fragment("blunt dephosphorylated DNA").rows())))
    uid_end = AdapterLigationGeometry.unavailable(
        "The defining article establishes ligation to the dephosphorylated 3′ fragment ends, "
        "but the fetched article does not disclose Supplementary Table 7 adapter-end bases.")
    w.ligate_adapter("Ligate the UID-bearing 3′ adapter", uid_end,
            chemistry.three_prime_adapter().rows(), name="3′-adapter-ligated duplex",
            note="The adapter is installed at the two 3′ fragment ends after dephosphorylation and blunting.")
    five_prime_end = AdapterLigationGeometry.unavailable(
        "The defining article establishes the later 5′-adapter ligation, but the fetched "
        "article does not disclose Supplementary Table 7 adapter-end bases.")
    w.ligate_adapter("Extend the UID and ligate the 5′ adapter", five_prime_end,
            chemistry.copied_barcode().rows(), name="duplex-barcoded library",
            note="Extension copies one random UID into both strands before redundant amplification.")
    w.react("Partition into Watson- and Crick-specific hemi-nested PCRs",
            chemistry.strand_families().rows(), name="oriented target families",
            note="A second nested PCR adds sample indices and further target specificity.")
    return w


def oligos():
    return ()
