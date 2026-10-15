"""Reusable molecular states for UMI-family and strand-aware consensus methods.

The constructors keep two meanings separate: a single-strand UID identifies one
amplified template strand, whereas a duplex UID is copied into both strands before
amplification and can therefore reconnect Watson and Crick families.
"""
from __future__ import annotations

from dataclasses import dataclass

from batch_ngs import seg
from chemdraw import Construct, Scene, Segment, feature


def uid_segment(name: str, length: int, identity: str, *, group: str | None = None,
                part: str | None = None) -> Segment:
    """Construct a random molecular identifier with validated positive length."""
    if length < 1:
        raise ValueError("a molecular UID must contain at least one base")
    return seg(name, "N" * length, "umi", placeholder=True,
               feature=feature(identity, "umi", "random", group=group, part=part))


@dataclass(frozen=True)
class SingleStrandUID:
    """A primer-borne UID assigned before redundant target amplification."""

    length: int = 14

    def assignment(self) -> Scene:
        target = [seg("target body", "X" * 20, placeholder=True),
                  seg("target-primer site", "T" * 12)]
        primer = [seg("universal PCR handle", "H" * 14, placeholder=True),
                  uid_segment("random UID", self.length, "safe_seqs_uid"),
                  seg("target-specific primer", "A" * 12, placeholder=True)]
        sc = Scene(); sc.strand("template", target, label="one target strand")
        sc.anneal("UID primer", primer, to="template",
                  pair=("target-specific primer", "target-primer site"),
                  label="UID-assignment primer",
                  unpaired=("universal PCR handle", "random UID"))
        sc.arrow("UID primer", "limited-cycle UID assignment")
        sc.mark("UID primer", "random UID", "one UID defines one strand family")
        sc.labels("UID primer")
        return sc

    def family(self) -> Scene:
        parts = [seg("universal PCR handle", "H" * 14, placeholder=True),
                 uid_segment("copied UID", self.length, "safe_seqs_uid"),
                 seg("target amplicon", "X" * 32, placeholder=True),
                 seg("opposite universal handle", "K" * 14, placeholder=True)]
        sc = Scene.duplex(parts, label="redundant copies in one UID family")
        sc.mark("top", "copied UID", "group reads before family consensus")
        sc.note("top", "a variant is accepted only when it dominates the UID family")
        sc.labels("top")
        return sc


@dataclass(frozen=True)
class CopiedDuplexUID:
    """A one-ended UID copied across a duplex before library amplification."""

    length: int = 14

    def three_prime_adapter(self) -> Scene:
        parts = [seg("blunt dephosphorylated insert", "X" * 32, placeholder=True),
                 uid_segment("3-prime adapter UID", self.length,
                             "saferseqs_duplex_uid"),
                 seg("3-prime adapter arm", "A" * 14, placeholder=True)]
        sc = Scene(); sc.strand("ligated end", parts, label="one 3′ end of the original duplex")
        sc.junction("ligated end", "blunt dephosphorylated insert",
                    "3-prime adapter UID", "3′ UID-adapter ligation")
        sc.note("ligated end", "the opposite 3′ fragment end receives the same adapter class")
        sc.labels("ligated end")
        return sc

    def copied_barcode(self) -> Scene:
        parts = [seg("5-prime adapter", "H" * 14, placeholder=True),
                 seg("genomic insert", "X" * 32, placeholder=True),
                 uid_segment("copied duplex UID", self.length,
                             "saferseqs_duplex_uid"),
                 seg("3-prime adapter arm", "K" * 14, placeholder=True)]
        sc = Scene.duplex(parts, label="duplex-barcoded original molecule")
        sc.junction("top", "5-prime adapter", "genomic insert", "5′ adapter ligation")
        sc.junction("top", "genomic insert", "copied duplex UID", "3′ UID-adapter ligation")
        sc.mark("top", "copied duplex UID", "same molecular identity on Watson and Crick")
        sc.note("top", "extension copies the UID; 5′-adapter ligation seals the duplex")
        sc.labels("top")
        return sc

    def strand_families(self) -> Scene:
        sc = Scene()
        sc.strand("Watson family", [uid_segment("Watson UID", self.length,
                                                "saferseqs_duplex_uid"),
                                     seg("Watson target", "X" * 30, placeholder=True)],
                  label="Watson-specific hemi-nested PCR")
        sc.strand("Crick family", [uid_segment("Crick UID", self.length,
                                               "saferseqs_duplex_uid"),
                                    seg("Crick target", "Y" * 30, placeholder=True)],
                  label="Crick-specific hemi-nested PCR")
        sc.note("Watson family", "paired by the shared 14-nt UID after family consensus")
        sc.labels("Watson family"); sc.labels("Crick family")
        return sc
