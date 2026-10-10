"""Reusable molecular states for newer multi-omic and end-mapping protocols.

These helpers encode the stable chemistry-level distinctions shared by several methods.
Protocol modules still own their exact oligos and reaction order.
"""
from __future__ import annotations

from chemdraw import Row, Scene, Segment, feature


def molecule(name: str, bases: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name, bases, tag, **kw)


def analyte_split(*branches: str, source: str = "one cell") -> list[Row]:
    """A physical split: branches share provenance, not a covalent barcode."""
    if len(branches) < 2:
        raise ValueError("an analyte split needs at least two branches")
    rows = [Row(chunks=[(source, None, False)])]
    for i, branch in enumerate(branches):
        rows.append(Row(chunks=[(("├─ " if i + 1 < len(branches) else "└─ ") + branch,
                                 None, False)]))
    return rows


def coamplified_analytes(*analytes: str, handle: str = "common amplification handle") -> list[Row]:
    """Different analytes acquire one handle before a pool is split downstream."""
    if len(analytes) < 2:
        raise ValueError("co-amplification needs at least two analytes")
    return [Row(chunks=[(f"{name} ** [{handle}]", "me", False)]) for name in analytes]


def indexed_molecule(*, payload: str, barcode_parts: tuple[tuple[str, int], ...],
                     umi: int = 0, target_barcode: int = 0,
                     label: str = "indexed molecule") -> Scene:
    """One strand carrying explicit, structured identifiers.

    Barcode length and multipart identity are constructor inputs, so downstream rendering
    and read-span derivation do not infer them from colour or placeholder letters.
    """
    parts = [molecule(payload, "X" * 30, placeholder=True)]
    for name, length in barcode_parts:
        if length < 1:
            raise ValueError("barcode parts must have positive lengths")
        key = name.lower().replace(" ", "_").replace("-", "_")
        parts.append(molecule(name, "B" * length, "cbc", placeholder=True,
                              feature=feature(f"cell_barcode_{key}", "cell_barcode",
                                              "combinatorial", group="cell_barcode",
                                              part=name)))
    if umi:
        parts.append(molecule("UMI", "U" * umi, "umi", placeholder=True,
                              feature=feature("umi", "umi", "random")))
    if target_barcode:
        parts.append(molecule("target barcode", "T" * target_barcode, "cbc",
                              placeholder=True,
                              feature=feature("target_barcode", "feature_barcode",
                                              "unknown")))
    sc = Scene(); sc.strand("product", parts, label=label)
    previous = payload
    for part in parts[1:]:
        sc.junction("product", previous, part.name, "barcode transfer")
        previous = part.name
    sc.labels("product")
    return sc


def multiplex_target_fragment(*, target_barcode_length: int = 8) -> Scene:
    """CUT&Tag-like product in which the tethered transposome names its target."""
    return indexed_molecule(payload="chromatin fragment", barcode_parts=(),
                            target_barcode=target_barcode_length,
                            label="target-tagged chromatin fragment")


def clicked_nascent_rna(*, barcode_length: int = 8, umi_length: int = 8) -> Scene:
    """Run-on RNA terminated in an alkyne nucleotide and clicked to an azide oligo."""
    rna = [molecule("nascent RNA", "R" * 28, placeholder=True),
           molecule("3-prime O-propargyl NTP", "N", placeholder=True)]
    tag = [molecule("click linker", "C" * 6, placeholder=True),
           molecule("cell barcode", "B" * barcode_length, "cbc", placeholder=True,
                    feature=feature("cell_barcode", "cell_barcode", "unknown")),
           molecule("UMI", "U" * umi_length, "umi", placeholder=True,
                    feature=feature("umi", "umi", "random")),
           molecule("RT handle", "H" * 16, placeholder=True)]
    sc = Scene(); sc.strand("RNA", [*rna, *tag], label="clicked nascent-RNA product")
    sc.junction("RNA", "3-prime O-propargyl NTP", "click linker", "triazole click junction")
    sc.mark("RNA", "cell barcode", "cell identity")
    sc.mark("RNA", "UMI", "molecule identity")
    sc.labels("RNA")
    return sc


def capped_full_length_cdna(*, retain_poly_a: bool = True,
                            cell_barcode: int = 0, umi: int = 0) -> Scene:
    """Full-length cDNA selected simultaneously by its copied cap and RNA 3-prime end."""
    parts = [molecule("5-prime cap-derived end", "C" * 8, placeholder=True),
             molecule("full-length cDNA", "X" * 34, placeholder=True)]
    if retain_poly_a:
        parts.append(molecule("copied poly(A) tail", "T" * 16))
    if cell_barcode:
        parts.append(molecule("cell barcode", "B" * cell_barcode, "cbc", placeholder=True,
                              feature=feature("cell_barcode", "cell_barcode", "unknown")))
    if umi:
        parts.append(molecule("long UMI", "U" * umi, "umi", placeholder=True,
                              feature=feature("umi", "umi", "random")))
    sc = Scene(); sc.strand("cDNA", parts, label="full-length cDNA")
    sc.mark("cDNA", "5-prime cap-derived end", "cap selected")
    if retain_poly_a:
        sc.mark("cDNA", "copied poly(A) tail", "native tail retained")
    sc.labels("cDNA")
    return sc


def native_three_prime_end(*, captured: str = "native DNA 3-prime OH") -> Scene:
    """A native DNA end joined to an end-specific adapter before bulk fragmentation."""
    parts = [molecule("genomic DNA upstream", "X" * 28, placeholder=True),
             molecule(captured, "E", placeholder=True),
             molecule("end-capture adapter", "A" * 18, "r1", placeholder=True)]
    sc = Scene(); sc.strand("captured end", parts, label="native-end capture product")
    sc.junction("captured end", captured, "end-capture adapter", "end-specific ligation")
    sc.labels("captured end")
    return sc


def proximity_product(*, left: str, right: str, bridge: str,
                      complex_barcode: int = 0) -> Scene:
    """A covalent contact product with an explicit bridge and optional complex ID."""
    parts = [molecule(left, "X" * 22, placeholder=True),
             molecule(bridge, "L" * 16, "me", placeholder=True)]
    if complex_barcode:
        parts.append(molecule("complex barcode", "B" * complex_barcode, "cbc",
                              placeholder=True,
                              feature=feature("complex_barcode", "molecular_complex_barcode",
                                              "random")))
    parts.append(molecule(right, "Y" * 22, placeholder=True))
    sc = Scene(); sc.strand("contact", parts, label="contact chimera")
    sc.junction("contact", left, bridge, "proximity ligation")
    tail = "complex barcode" if complex_barcode else bridge
    if complex_barcode:
        sc.junction("contact", bridge, tail, "barcode ligation")
    sc.junction("contact", tail, right, "proximity ligation")
    sc.labels("contact")
    return sc


def overload_then_unpack(*, droplet_barcode: int, well_barcode: int) -> list[Row]:
    """UDA-style overloaded droplets followed by a second index after unpacking."""
    if droplet_barcode < 1 or well_barcode < 1:
        raise ValueError("both overload and unpacking barcodes must be present")
    return [
        Row(chunks=[(f"many cells per droplet → [{droplet_barcode}-nt droplet barcode]", "cbc", False)]),
        Row(chunks=[(f"unpack pool into wells → [{well_barcode}-nt well barcode]", "cbc", False)]),
        Row(chunks=[("cell identity = droplet barcode × well barcode", None, False)]),
    ]
