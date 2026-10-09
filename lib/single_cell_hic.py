"""Reusable, deliberately small model for single-cell chromosome-conformation pages.

The workflow object is the policy boundary: it rejects combinations that are not real
protocol classes (for example, streptavidin capture without biotin marking), while each
paper-specific module remains responsible for its enzyme, ordering and amplification.
"""
from __future__ import annotations

from dataclasses import dataclass

import illumina as il
import seqprimers as sp
from chemdraw import Construct, Segment, feature, revcomp
from restriction import Digest, RestrictionEnzyme


@dataclass(frozen=True)
class ContactWorkflow:
    name: str
    enzyme: RestrictionEnzyme
    isolate: str
    marking: str                 # biotin-fill, biotin-blunt, none, or bridge-adaptor
    amplification: str           # PCR, MDA, META, or bisulfite-PCR
    barcoding_rounds: int = 0
    capture: bool = False

    def __post_init__(self) -> None:
        if self.marking not in {"biotin-fill", "biotin-blunt", "none", "bridge-adaptor"}:
            raise ValueError(f"unknown end-marking mode: {self.marking}")
        if self.amplification not in {"PCR", "MDA", "META", "bisulfite-PCR"}:
            raise ValueError(f"unknown amplification mode: {self.amplification}")
        if self.capture and not self.marking.startswith("biotin"):
            raise ValueError("streptavidin capture requires a biotin-marked workflow")
        if self.marking == "bridge-adaptor" and self.barcoding_rounds != 2:
            raise ValueError("indexed bridge-adaptor Hi-C requires two barcode rounds")

    def digest(self) -> Digest:
        return Digest(self.enzyme, "X" * 16, "X" * 16)

    def junction(self) -> Construct:
        d = self.digest()
        if self.marking == "biotin-fill":
            middle = d.fill_in(biotin_base="A").junction
        elif self.enzyme.end == "blunt":
            middle = self.enzyme.site
        else:
            # Cohesive ends re-form one recognition site when compatible ends ligate.
            middle = self.enzyme.site
        return Construct([
            Segment("locus A", "X" * 22, placeholder=True),
            Segment("contact junction", middle, tag="me"),
            Segment("locus B", "X" * 22, placeholder=True),
        ], name=f"{self.name} contact product")


SEQ_PRIMERS = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["R2"])


def inferred_truseq_library(insert: Construct, name: str) -> Construct:
    """Canonical paired-end shell when a paper names, but does not print, the kit arms."""
    lib = Construct([
        Segment("P5", il.P5, "p5", inferred=True),
        Segment("Read 1 arm", il.TRUSEQ_READ1[4:], "r1", inferred=True),
        *list(insert), Segment("dA junction", "A", inferred=True),
        Segment("Index 1 / Read 2 arm", il.INDEX1_PRIMER, "r2", inferred=True),
        Segment("i7 reverse complement", "I" * 8, "cbc", placeholder=True, inferred=True,
                feature=feature("sample_i7", "sample_index", "unknown")),
        Segment("P7 reverse complement", il.P7_RC, "p7", inferred=True),
    ], name=name)
    errors = sp.verify(lib, SEQ_PRIMERS,
                       required_roles=("Read 1", "Index 1 (i7)", "Read 2"))
    if errors:
        raise ValueError("invalid Hi-C library: " + "; ".join(errors))
    return lib


def unresolved_illumina_library(insert: Construct, name: str, *, indexed: bool | None) -> Construct:
    """Library whose source identifies Illumina prep but not the adapter oligos.

    Role-labelled placeholders deliberately prevent a named kit generation, index count,
    or modern primer set from being smuggled into a historical protocol.
    """
    index = ("indexed adapter (sequence not reported)" if indexed is True else
             "adapter/index region (not reported)" if indexed is None else
             "adapter arm (sequence not reported)")
    semantic = (feature("sample_index_region", "sample_index", "unknown",
                        note="index sub-boundary is not reported")
                if indexed is True else None)
    return Construct([
        Segment("Illumina left arm (not reported)", "X" * 20,
                placeholder=True, inferred=True),
        *list(insert),
        Segment(index, "X" * 20, placeholder=True, inferred=True, feature=semantic),
    ], name=name)


# Nagano 2013 and Stevens 2017 print this historical paired-end primer in their
# supplementary methods.  It is not the modern TruSeq Read 2 primer.
HISTORICAL_PE_READ2 = "CGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT"


def historical_inline_pe_library(insert: Construct, name: str, barcode: str = "CAA") -> Construct:
    """Exact custom 3-bp inline-tag architecture shared by Nagano and Stevens."""
    if len(barcode) != 3 or any(b not in "ACGT" for b in barcode):
        raise ValueError("historical Hi-C identification tag must be three DNA bases")
    return Construct([
        Segment("P5", il.P5, "p5"),
        Segment("historical Read 1 arm", il.TRUSEQ_READ1[4:], "r1"),
        Segment("3-bp identification tag", barcode, "cbc",
                feature=feature("cell_inline_a", "cell_barcode", "whitelist",
                                whitelist="published identification-tag set", group="cell_id", part="read 1 end")),
        *list(insert),
        Segment("opposite 3-bp identification tag", revcomp(barcode), "cbc",
                feature=feature("cell_inline_b", "cell_barcode", "whitelist",
                                whitelist="published identification-tag set", group="cell_id", part="read 2 end")),
        Segment("historical Read 2 site", revcomp(HISTORICAL_PE_READ2), "r2"),
        Segment("P7 reverse complement", il.P7_RC, "p7"),
    ], name=name)
