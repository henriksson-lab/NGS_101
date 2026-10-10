"""AIRR Community bulk UMI 5-prime-RACE BCR library."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from batch_ngs import seg
from chemdraw import Construct, MolecularState, Scene, Workflow, feature

TITLE = "Bulk UMI 5′-RACE AIRR-seq — B-cell receptor repertoires"
NOTES = "01_airr-umi-5-race.html"
SOURCE = 'Authoritative protocol: <a href="https://doi.org/10.1007/978-1-0716-2115-8_19">Gupta et al., AIRR Community (2022)</a>.'
SUMMARY = "Oligo(dT)-primed SMART reverse transcription places a 12-nt UMI and universal handle at the transcript 5′ end. Semi-nested constant-region PCR then captures complete BCR V(D)J regions and adds dual Illumina indices."
CAVEAT = "The kit's gene-specific, SMART and complete sequencing-primer sequences are proprietary. Published molecular roles, lengths and Read 1/Read 2 orientation are shown without inventing bases or primer landing sites."
UMI_LENGTH = 12
SMART_LINKER_LENGTH = 7

FINAL_LIBRARY = Construct([
    seg("P5 arm", "[P5 arm]", "p5", placeholder=True),
    seg("eight-base i5 sample index", "J" * 8, "cbc", placeholder=True,
        feature=feature("sample_index_i5", "sample_index", "whitelist",
                        whitelist="published BCR indexing reverse-primer set")),
    seg("Read 1 entry site", "[Read 1 site]", "r1", placeholder=True),
    seg("short constant-region sequence", "C" * 60, placeholder=True, length_bp=60),
    seg("BCR V(D)J sequence", "V" * 300, placeholder=True, length_bp=300),
    seg("seven-base SMART linker", "L" * SMART_LINKER_LENGTH, placeholder=True),
    seg("twelve-base molecular UMI", "U" * UMI_LENGTH, "umi", placeholder=True,
        feature=feature("bcr_molecule_umi", "umi", "random")),
    seg("Read 2 entry site", "[Read 2 site]", "r2", placeholder=True),
    seg("eight-base i7 sample index", "I" * 8, "cbc", placeholder=True,
        feature=feature("sample_index_i7", "sample_index", "whitelist",
                        whitelist="published BCR PCR2 universal-forward set")),
    seg("P7 arm", "[P7 arm]", "p7", placeholder=True),
], name="dual-index bulk BCR AIRR-seq library")
SEQ_PRIMERS = ()
FINAL_CAPTION = "Read 1 starts in the BCR constant-region side; Read 2 starts at the SMART end, reads the 12-base UMI and then the seven-base linker before entering the upstream V region."
SEQUENCING_ENDING = "Paired reads cover constant-region identity and the rearranged variable region. The first 19 Read 2 bases are the 12-nt UMI followed by the seven-base SMART/template-switch linker. Dual eight-base indexes identify the sample."
SEQUENCING_UNAVAILABLE = "The source explicitly says the vendor primers are unavailable. It gives the read orientation and index/UMI lengths, but not exact sequencing-primer sequences from which binding could be computed."

_initial = Scene(); _initial.strand("mRNA", [seg("BCR transcript", "R" * 54, placeholder=True), seg("poly(A)", "A" * 18)], label="polyadenylated BCR mRNA")
INITIAL_ROWS = tuple(_initial.rows()); INITIAL_NAME = "Bulk BCR mRNA"


def workflow():
    cdna = Scene(); cdna.strand("cDNA", [seg("dT-primer handle", "X" * 14, placeholder=True),
                                         seg("BCR cDNA", "D" * 54, placeholder=True),
                                         seg("twelve-base UMI", "U" * UMI_LENGTH, "umi", placeholder=True,
                                             feature=feature("bcr_molecule_umi", "umi", "random")),
                                         seg("SMART universal handle", "L" * 18, placeholder=True)], label="template-switched first-strand cDNA")
    cdna.mark("cDNA", "twelve-base UMI", "one UMI per parental cDNA")
    pcr1 = Scene.duplex([seg("universal SMART end", "L" * 18, placeholder=True),
                         seg("complete V(D)J", "V" * 54, placeholder=True),
                         seg("constant-region segment", "C" * 32, placeholder=True)], label="PCR1 BCR amplicon")
    pcr1.mark("top", "constant-region segment", "IgG / IgM / IgK / IgL reverse primer")
    final = Scene.duplex(list(FINAL_LIBRARY), label="semi-nested PCR2 library",
                         unpaired=tuple(s.name + "'" for s in FINAL_LIBRARY if s.is_role_token()))
    w = Workflow(MolecularState(INITIAL_NAME, INITIAL_ROWS))
    w.react("Prime reverse transcription from poly(A)", cdna.rows(), name="UMI-tagged first-strand cDNA",
            note="SMARTScribe adds non-templated bases; the SMART UMI oligo template-switches and copies a 12-nt UMI plus universal handle.")
    w.react("PCR1 from universal handle to constant region", pcr1.rows(), name="Chain-specific PCR1 amplicon",
            note="Separate IgG, IgM, kappa and lambda reverse primers retain complete V(D)J plus part of the constant region.")
    w.react("Semi-nested PCR2 and dual indexing", final.rows(), name="Sequencing library",
            note="PCR2 uses chain-specific nested reverse primers and indexed universal forward primers to add P5, P7, i5 and i7.")
    return w
