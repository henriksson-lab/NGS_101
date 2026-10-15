"""Bounded chemistry shared by ATAC-seq sample-preparation variants.

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
from pcr import PCRGeometry
from protocol_graph import (Disposition, InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, StepKind, Terminal)
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
