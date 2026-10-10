"""Bivalent bridges that covalently join an RNA end to a nearby DNA end.

The bridge is deliberately more structured than a generic adapter: one oligo has a
single-stranded RNA ligation arm followed by a DNA duplex arm, while its partner is DNA.
Construction validates the duplex, the pre-adenylated RNA end, and the positions of the
internal molecular barcode and biotinylated base.
"""
from __future__ import annotations

from dataclasses import dataclass

from chemdraw import Scene, Segment, feature, is_real_dna, revcomp


@dataclass(frozen=True)
class BivalentRNADNABridge:
    name: str
    dna_strand: str
    rna_fixed_5: str
    barcode_length: int
    rna_terminal: str
    duplex_arm: str
    biotin_offset: int
    preadenylated_rna_5: bool
    barcode_id: str

    def __post_init__(self) -> None:
        if (not is_real_dna(self.dna_strand) or not is_real_dna(self.duplex_arm)
                or "U" in self.dna_strand.upper() or "U" in self.duplex_arm.upper()):
            raise ValueError("bridge duplex arms must be explicit DNA")
        if self.duplex_arm != revcomp(self.dna_strand):
            raise ValueError("bridge DNA arms must be exact reverse complements")
        if not self.rna_fixed_5 or "U" not in self.rna_fixed_5 or any(
                b not in "ACGUN" for b in self.rna_fixed_5 + self.rna_terminal):
            raise ValueError("bridge RNA ligation arm must be written as RNA with U")
        if self.barcode_length < 1:
            raise ValueError("bridge molecular barcode must have positive length")
        if not self.preadenylated_rna_5:
            raise ValueError("ATP-free RNA ligation requires a pre-adenylated bridge end")
        if not 0 <= self.biotin_offset < len(self.duplex_arm):
            raise ValueError("bridge biotin position lies outside the duplex arm")
        if self.duplex_arm[self.biotin_offset] != "T":
            raise ValueError("the structured internal biotin must modify a thymidine")

    def _hybrid_segments(self) -> list[Segment]:
        before, bio, after = (self.duplex_arm[:self.biotin_offset],
                              self.duplex_arm[self.biotin_offset],
                              self.duplex_arm[self.biotin_offset + 1:])
        out = [Segment("RNA-arm fixed bases", self.rna_fixed_5, "me"),
               Segment("RNA-arm molecular barcode", "N" * self.barcode_length, "umi",
                       placeholder=True,
                       feature=feature(self.barcode_id, "umi", "random")),
               Segment("RNA-arm terminal bases", self.rna_terminal, "me")]
        if before:
            out.append(Segment("duplex before biotin", before, "me"))
        out.append(Segment("internal biotin-dT", bio, "me",
                           note="internal biotinylated thymidine"))
        if after:
            out.append(Segment("duplex after biotin", after, "me"))
        return out

    @property
    def duplex_names(self) -> tuple[str, ...]:
        return tuple(x for x in ("duplex before biotin", "internal biotin-dT",
                                 "duplex after biotin")
                     if any(s.name == x for s in self._hybrid_segments()))

    def scene(self, *, at: int = 0) -> Scene:
        """The free bridge, with the entire mixed strand pairing-validated."""
        sc = Scene()
        sc.strand("bridge DNA", [Segment("bridge DNA strand", self.dna_strand, "me")],
                  label="5′-phosphorylated DNA strand", at=at)
        hybrid = self._hybrid_segments()
        unpaired = ("RNA-arm fixed bases", "RNA-arm molecular barcode",
                    "RNA-arm terminal bases")
        sc.anneal("bridge hybrid", hybrid, to="bridge DNA",
                  pair=(self.duplex_names[-1], "bridge DNA strand"),
                  label="5′-App mixed RNA/DNA strand", unpaired=unpaired)
        sc.mark("bridge hybrid", "RNA-arm molecular barcode", "molecular barcode")
        sc.mark("bridge hybrid", "internal biotin-dT", "streptavidin capture")
        return sc

    def contact_scene(self, *, rna_length: int = 20, dna_length: int = 20) -> Scene:
        """Bridge after both ligations, retaining RNA, both bridge strands and DNA."""
        top = [Segment("bridge DNA strand", self.dna_strand, "me"),
               Segment("AluI genomic end", "D" * dna_length, placeholder=True)]
        hybrid = [Segment("chromatin RNA", "R" * rna_length, placeholder=True),
                  *self._hybrid_segments()]
        sc = Scene()
        sc.strand("DNA side", top, label="bridge DNA strand ligated to genomic DNA")
        sc.anneal("RNA side", hybrid, to="DNA side",
                  pair=(self.duplex_names[-1], "bridge DNA strand"),
                  label="chromatin RNA ligated to mixed RNA/DNA strand",
                  unpaired=("chromatin RNA", "RNA-arm fixed bases",
                            "RNA-arm molecular barcode", "RNA-arm terminal bases"))
        sc.junction("RNA side", "chromatin RNA", "RNA-arm fixed bases", "RNA ligation")
        sc.junction("DNA side", "bridge DNA strand", "AluI genomic end", "DNA ligation")
        sc.mark("RNA side", "internal biotin-dT", "streptavidin capture")
        return sc

    def unligated_contact_scene(self, *, rna_length: int = 20,
                                dna_length: int = 20) -> Scene:
        """Fixed RNA and DNA together with the free, annealed bridge reagent."""
        sc = Scene()
        sc.strand("chromatin RNA", [Segment("chromatin-associated RNA", "R" * rna_length,
                                             placeholder=True)],
                  label="crosslinked chromatin RNA")
        sc.strand("genomic DNA", [Segment("nearby genomic DNA", "D" * dna_length,
                                           placeholder=True)],
                  label="nearby AluI-cut genomic DNA")
        sc.strand("bridge DNA", [Segment("bridge DNA strand", self.dna_strand, "me")],
                  label="5′-phosphorylated bridge DNA", at=max(rna_length, dna_length) + 8)
        hybrid = self._hybrid_segments()
        sc.anneal("bridge hybrid", hybrid, to="bridge DNA",
                  pair=(self.duplex_names[-1], "bridge DNA strand"),
                  label="5′-App mixed RNA/DNA bridge strand",
                  unpaired=("RNA-arm fixed bases", "RNA-arm molecular barcode",
                            "RNA-arm terminal bases"))
        sc.mark("bridge hybrid", "RNA-arm molecular barcode", "molecular barcode")
        sc.mark("bridge hybrid", "internal biotin-dT", "streptavidin capture")
        return sc
