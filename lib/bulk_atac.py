"""Shared chemistry for the named bulk ATAC-seq variants.

Omni-ATAC and Fast-ATAC differ before and during transposition, but both hand the
same TDE1-tagged molecule to the original Buenrostro single-index amplification.
Keeping that boundary here prevents sample-preparation variants from drifting into
different library geometries.
"""
from __future__ import annotations

from dataclasses import dataclass
from types import SimpleNamespace

import illumina as il
import nextera as nx
from batch_ngs import nextera_library, seg
from protocol_model import chain, molecule, role


# Buenrostro et al., Current Protocols 2015, Basic Protocol 1 reagent table.
# Only seven mosaic-end bases occur at the 3' ends of the amplification primers;
# the complete 19-bp mosaic end was already installed by Tn5.
ATAC_AD1_NOMX = il.P5 + nx.S5 + nx.ME[:7]
ATAC_INDEX_1 = "TAAGGCGA"
ATAC_AD2_1 = il.P7 + "TCGCCTTA" + nx.S7 + nx.ME[:6]


@dataclass(frozen=True)
class AtacVariant:
    key: str
    title: str
    source: str
    summary: str
    caveat: str
    first_step: str
    first_state: str
    first_note: str
    transposition_action: str
    transposition_note: str
    cleanup_note: str
    sequencing_intro: str
    read_lengths: dict[str, int] | None


VARIANTS = {
    "omni-atac": AtacVariant(
        key="omni-atac",
        title="Omni-ATAC — low-background chromatin accessibility sequencing",
        source='<a href="https://doi.org/10.1038/nmeth.4396">Corces et al., Nature Methods (2017)</a> and Supplementary Protocol 1.',
        summary="NP-40, digitonin and Tween permeabilize cells, a Tween wash removes mitochondria, and PBS is retained during Tn5 transposition; the resulting accessible-DNA fragments enter the original single-index ATAC-seq amplification.",
        caveat="This page follows the defining 2017 Omni-ATAC implementation with Illumina TDE1 and the original Ad1/Ad2 single-index primers. Later protocols often substitute home-made Tn5 or dual-index primers; those substitutions are not silently folded into this construct.",
        first_step="Lyse 50,000 cells in multi-detergent ATAC-RSB, then wash nuclei",
        first_state="washed, permeabilized nuclei depleted of free mitochondria",
        first_note="The defining recipe uses 0.1% NP-40, 0.1% Tween-20 and 0.01% digitonin for 3 minutes on ice, followed by a 0.1% Tween-20 wash before pelleting.",
        transposition_action="Transpose accessible chromatin in the PBS-containing Omni mix",
        transposition_note="The 50-microlitre reaction contains 25 microlitres 2x TD buffer, 2.5 microlitres TDE1, 16.5 microlitres PBS, 0.01% digitonin and 0.1% Tween-20; incubate 30 minutes at 37 degrees C with 1,000-rpm mixing.",
        cleanup_note="A silica cleanup stops transposition. Five initial PCR cycles fill the nine-base Tn5 gaps and install flow-cell/index sequence; qPCR on an aliquot determines the minimum additional cycles.",
        sequencing_intro="The defining method uses paired-end Illumina sequencing so the two reads enter accessible genomic DNA from opposite mosaic ends. It does not prescribe one immutable genomic-read length; the single eight-base Index 1 read reports the Ad2 sample index.",
        read_lengths={"Index 1 (i7)": 8},
    ),
    "fast-atac": AtacVariant(
        key="fast-atac",
        title="Fast-ATAC — one-step accessibility sequencing for primary blood cells",
        source='<a href="https://doi.org/10.1038/ng.3646">Corces et al., Nature Genetics (2016)</a> and its Online Methods.',
        summary="Digitonin permeabilizes 5,000 sorted primary blood cells inside the Tn5 reaction itself, eliminating a separate lysis step while preserving the original single-index ATAC-seq library architecture.",
        caveat="Fast-ATAC names the one-step digitonin sample preparation, not a distinct adapter system. The page therefore shares the published Ad1/Ad2 single-index construct with Omni-ATAC while keeping the defining Fast-ATAC reaction conditions separate.",
        first_step="Sort 5,000 primary blood cells and pellet them cold",
        first_state="small intact-cell pellet",
        first_note="The defining workflow sorts purified blood populations, pellets them at 500 r.c.f. for 5 minutes at 4 degrees C, and removes the supernatant without a separate detergent lysis.",
        transposition_action="Permeabilize and transpose together with digitonin",
        transposition_note="The 50-microlitre reaction contains 25 microlitres 2x TD buffer, 2.5 microlitres TDE1, 0.5 microlitre 1% digitonin and 22 microlitres water; incubate 30 minutes at 37 degrees C with 300-rpm mixing.",
        cleanup_note="MinElute cleanup stops the one-tube reaction. Five initial PCR cycles fill the Tn5 gaps and install the complete library arms; an aliquot qPCR chooses the remaining cycles.",
        sequencing_intro="The defining study sequenced Fast-ATAC libraries as 50-base paired-end reads on HiSeq 4000. Standard Nextera Read 1 and Read 2 primers enter accessible DNA from opposite mosaic ends; the eight-base Index 1 read reports the Ad2 sample index.",
        read_lengths={"Read 1": 50, "Index 1 (i7)": 8, "Read 2": 50},
    ),
}


def library(name: str):
    """The historical Ad1_noMX/Ad2 single-index ATAC library."""
    return nextera_library(
        [seg("accessible genomic DNA", "X" * 44, placeholder=True)],
        name + " library", dual_index=False, i7_nt=8)


def workflow(spec: AtacVariant):
    return chain(molecule("fresh cells", (role("intact cells with nuclear and mitochondrial DNA"),)), (
        (spec.first_step,
         molecule(spec.first_state, (role("permeabilized nuclei"), role("accessible chromatin"))),
         spec.first_note),
        (spec.transposition_action,
         molecule("heterologously tagmented accessible-DNA fragment", (
             role("S5 plus mosaic-end transferred strand"),
             role("accessible genomic DNA"),
             role("S7 plus mosaic-end transferred strand")),
             junctions=(("S5 plus mosaic-end transferred strand", "accessible genomic DNA", "Tn5 strand transfer"),
                        ("accessible genomic DNA", "S7 plus mosaic-end transferred strand", "Tn5 strand transfer"))),
         spec.transposition_note + " S5/S7 is the amplifiable end combination; S5/S5 and S7/S7 products lack one flow-cell arm."),
        ("Purify, gap-fill and pre-amplify for five cycles",
         molecule("single-index ATAC library", (
             role("P5 and S5/left mosaic end"), role("accessible genomic insert"),
             role("right mosaic end/S7, eight-base i7 and P7"))),
         spec.cleanup_note),
        ("Use aliquot qPCR to choose only the necessary extra cycles",
         molecule("amplified, complexity-preserving ATAC library", (
             role("P5/S5/ME"), role("accessible genomic insert"), role("ME/S7/i7/P7"))),
         "The fluorescence threshold is used to avoid a fixed over-amplification programme; the completed molecule is ready for paired-end size-aware analysis."),
    ))


def oligos():
    return (
        f"Ad1_noMX — 5′-{ATAC_AD1_NOMX}-3′",
        f"Ad2.1_TAAGGCGA — 5′-{ATAC_AD2_1}-3′",
        f"Ad2.1 sample-sheet Index 1 value — 5′-{ATAC_INDEX_1}-3′",
        f"TDE1 mosaic-end sequence installed at each genomic junction — 5′-{nx.ME}-3′",
    )


def get(key: str):
    spec = VARIANTS[key]
    final, primers = library(spec.title.split(" — ", 1)[0])
    return SimpleNamespace(
        TITLE=spec.title,
        NOTES="01_" + key.replace("-", "_") + ".html",
        SOURCE=spec.source,
        SUMMARY=spec.summary,
        CAVEAT=spec.caveat,
        oligos=oligos,
        workflow=lambda: workflow(spec),
        FINAL_LIBRARY=final,
        FINAL_CAPTION="The original Ad1_noMX and Ad2 architecture is shown: no i5 index is present. Tn5 supplies the complete 19-bp mosaic ends, while PCR appends P5 on the S5 side and P7 plus one eight-base i7 sample index on the S7 side.",
        SEQ_PRIMERS=primers,
        READ_LENGTHS=spec.read_lengths,
        SEQUENCING_INTRO=spec.sequencing_intro,
        SEQUENCING_ENDING="",
        SEQUENCING_UNAVAILABLE="",
    )
