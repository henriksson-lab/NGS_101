"""Reusable molecular states for RNA-specialized sequencing workflows."""
from __future__ import annotations

from chemdraw import Construct, Row, Scene, Segment, circle_rows, feature
from dumbbell import Dumbbell, dumbbell_rows


def role(name: str, width: int, tag: str | None = None, *, inferred: bool = False,
         **kw) -> Segment:
    """A length-preserving unknown region whose evidence state cannot be forgotten."""
    if width < 1:
        raise ValueError(f"{name}: molecular roles must have positive width")
    return Segment(name, "X" * width, tag, placeholder=True, inferred=inferred, **kw)


def rna_fragment(name: str = "RNA fragment", width: int = 32) -> list[Segment]:
    return [role(name, width)]


def adapter_ligation_scene(*, fragment_name: str = "RNA fragment",
                           adapter_name: str = "3-prime RNA adapter",
                           inferred: bool = False) -> Scene:
    """Single-stranded 3-prime RNA adapter ligation with a named exact junction."""
    parts = [role(fragment_name, 32), role(adapter_name, 18, "r2", inferred=inferred)]
    sc = Scene(); sc.strand("ligated RNA", parts, label="adapter-ligated RNA")
    sc.junction("ligated RNA", fragment_name, adapter_name, "ligation")
    sc.labels("ligated RNA")
    return sc


def template_switch_scene(*, body: str = "RNA body", tail: str = "poly(A)",
                          switch: str = "template-switch handle",
                          inferred: bool = False, umi: bool = False) -> Scene:
    """Poly(A)-primed RT with a template-switch handle at the completed cDNA 5-prime end."""
    mrna = [role(body, 30), Segment(tail, "A" * 14)]
    primer = [role("RT/PCR handle", 14, "r2", inferred=inferred), Segment("poly(dT)", "T" * 14)]
    sc = Scene(); sc.strand("RNA", mrna, label="polyadenylated RNA")
    sc.anneal("RT primer", primer, to="RNA", pair=("poly(dT)", tail),
              label="RT primer", unpaired=("RT/PCR handle",))
    sc.arrow("RT primer", "reverse transcription to the RNA 5-prime end")
    switch_parts = [role(switch, 15, "r1", inferred=inferred)]
    if umi:
        switch_parts.append(role("UMI", 8, "umi", inferred=inferred,
                                 feature=feature("umi", "umi", "random",
                                                 note="exact length is proprietary")))
    switch_parts.append(Segment("rGrGrG", "GGG"))
    sc.strand("switch oligo", switch_parts, label="template-switch oligo")
    if umi:
        sc.mark("switch oligo", "UMI", "UMI")
    return sc


def mutational_read_rows() -> list[Row]:
    return [
        Row(chunks=[("unmodified RNA   ———— paired / flexible nucleotides ————", None, False)]),
        Row(chunks=[("modified RNA     ———— ● ——— ●● ———— ● ——————————————", "umi", False)]),
        Row(chunks=[("MaP reverse transcription: chemical adducts → substitutions in cDNA", "r1", False)]),
    ]


def nanopore_entry_rows(molecule: str, *, amplified: bool) -> list[Row]:
    source = "PCR-amplified double-stranded cDNA" if amplified else "PCR-free double-stranded cDNA"
    return [
        Row(chunks=[(source + " ** rapid/sequencing adapter ** motor", "me", False)]),
        Row(chunks=[(f"motor presents one strand of {molecule} to the pore  ————————>", None, False)]),
    ]


def isoseq_smrtbell(insert: Construct) -> tuple[Dumbbell, list[Row]]:
    """Build a closed Iso-Seq SMRTbell; topology is enforced by :class:`Dumbbell`."""
    mol = Dumbbell(insert, "N" * 18, "N" * 18, name="Iso-Seq SMRTbell",
                   insert_overhang="A", adapter_overhang="T")
    return mol, dumbbell_rows(mol)


def circular_cdna_rows(parts: list[Segment]) -> list[Row]:
    return circle_rows(Construct(parts, name="circular cDNA"), "cDNA circularisation")
