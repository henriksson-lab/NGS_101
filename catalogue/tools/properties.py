"""Authoritative, validated properties used by the protocol finder.

The TSV is deliberately separate from prose metadata.  Every protocol must classify every
facet, and this loader refuses missing rows, empty cells, unknown values and duplicate tags.
That makes completeness a property of the function used by the site build rather than a
parallel test assertion.
"""
from __future__ import annotations

import csv
from collections import OrderedDict
from pathlib import Path

HERE = Path(__file__).resolve().parent
TSV = HERE.parent / "properties.tsv"


FACETS = OrderedDict([
    ("index_introduction", ("Index introduction", (
        "PCR", "ligation", "transposition", "reverse transcription", "extension / fill-in",
        "combinatorial indexing", "pre-indexed adaptor", "no index"))),
    ("index_architecture", ("Index architecture", (
        "TruSeq", "Nextera", "custom Illumina", "DNBSEQ", "PacBio barcode", "Nanopore barcode",
        "inline barcode", "UDI", "CDI", "single index", "dual index", "no barcode"))),
    ("assay", ("Assay type", (
        "RNA-seq", "RNA detection", "DNA-seq / WGS", "ATAC / accessibility", "DNA methylation",
        "chromatin conformation", "chromatin protein mapping", "CRISPR screening", "protein / feature detection",
        "translation profiling", "nascent transcription", "RNA-protein interaction",
        "promoter / RNA 5-prime mapping", "RNA structure", "targeted DNA sequencing",
        "amplicon profiling", "metagenomics", "spatial / in situ", "telomere", "multiomic"))),
    ("platform", ("Sequencing platform", (
        "Illumina", "DNBSEQ", "PacBio", "Oxford Nanopore", "in situ imaging", "flow cytometry"))),
    ("partitioning", ("Partitioning", (
        "bulk", "plate / well", "droplet", "microwell", "combinatorial indexing",
        "single nucleus", "single cell", "spatial / in situ"))),
    ("starting_material", ("Starting material", (
        "genomic DNA", "total RNA", "poly(A) RNA", "small RNA", "chromatin",
        "fixed cells / nuclei", "damaged DNA", "native RNA", "adapter-ligated library",
        "amplicon", "pre-amplified cDNA", "microbial community DNA",
        "ribosome-protected RNA", "nascent RNA"))),
    ("fragmentation", ("Fragmentation or entry", (
        "mechanical", "enzymatic", "restriction digest", "Tn5 / tagmentation",
        "RNase / MNase", "no fragmentation"))),
    ("adapter_installation", ("Adapter installation", (
        "ligation", "Tn5", "PCR-added", "reverse transcription", "template switching", "tailing + priming",
        "padlock / circularization", "splint ligation", "hairpin / dumbbell", "full-length adaptor",
        "hybridization scaffold"))),
    ("amplification", ("Amplification", (
        "PCR-free", "endpoint PCR", "linear amplification", "rolling-circle amplification",
        "whole-genome amplification", "whole-transcriptome amplification", "two-stage PCR",
        "branched-DNA signal amplification"))),
    ("topology", ("Molecular topology", (
        "linear", "circular", "hairpin", "dumbbell / SMRTbell", "concatemer",
        "proximity-ligation junction", "RNA-DNA hybrid", "DNA nanoball",
        "branched hybridization tree"))),
    ("strand_handling", ("Strand handling", (
        "unstranded", "directional RNA", "strand displacement", "second-strand destruction",
        "single-stranded library", "duplex sequencing", "native RNA sequencing",
        "RNA hybridization"))),
    ("identifiers", ("Molecular identifiers", (
        "UMI", "cell barcode", "sample index", "guide barcode", "spatial barcode",
        "lineage barcode", "no UMI"))),
    ("read_structure", ("Read structure", (
        "paired-end", "single-end", "index reads", "inline barcode", "barcode in Read 1",
        "barcode in Read 2", "custom sequencing primer", "long-read consensus",
        "microscopy", "flow cytometry"))),
    ("selection", ("Selection or enrichment", (
        "poly(A) selection", "rRNA depletion", "size selection", "hybrid capture",
        "restriction-site selection", "spatial capture", "affinity / antibody capture",
        "guide-specific enrichment", "circularization selection", "ribosome-footprint selection",
        "cap selection", "none"))),
    ("conversion", ("Base conversion or marking", (
        "bisulfite", "enzymatic methyl conversion", "GpC methyltransferase",
        "adenine methyltransferase", "dUTP strand marking", "chemical conversion", "none"))),
    ("availability", ("Availability", (
        "published academic protocol", "commercial kit", "discontinued / historical",
        "exact oligos public", "proprietary sequences", "schematic complete",
        "inferred regions present"))),
])

COLUMNS = ["dir", "profile", "overrides"]


def _profile(**values: str) -> dict[str, tuple[str, ...]]:
    missing = set(FACETS) - set(values)
    extra = set(values) - set(FACETS)
    if missing or extra:
        raise ValueError(f"bad property profile: missing={sorted(missing)}, extra={sorted(extra)}")
    return {key: tuple(x.strip() for x in value.split("|") if x.strip())
            for key, value in values.items()}


_ILLUMINA = dict(
    index_introduction="PCR", index_architecture="custom Illumina|dual index",
    assay="DNA-seq / WGS", platform="Illumina", partitioning="bulk",
    starting_material="genomic DNA", fragmentation="mechanical",
    adapter_installation="ligation", amplification="endpoint PCR", topology="linear",
    strand_handling="unstranded", identifiers="sample index|no UMI",
    read_structure="paired-end|index reads", selection="none", conversion="none",
    availability="published academic protocol|exact oligos public|schematic complete")


def _from(base: dict[str, str], **changes: str) -> dict[str, tuple[str, ...]]:
    return _profile(**(base | changes))


PROFILES = {
    "illumina_dna": _profile(**_ILLUMINA),
    "illumina_rna": _from(_ILLUMINA, assay="RNA-seq", starting_material="total RNA",
                           fragmentation="enzymatic"),
    "plate_scrna": _from(
        _ILLUMINA, assay="RNA-seq", partitioning="plate / well|single cell",
        starting_material="poly(A) RNA", fragmentation="Tn5 / tagmentation",
        adapter_installation="template switching|Tn5|PCR-added",
        amplification="whole-transcriptome amplification|endpoint PCR",
        identifiers="sample index|no UMI", read_structure="paired-end|index reads"),
    "droplet_scrna": _from(
        _ILLUMINA, assay="RNA-seq", partitioning="droplet|single cell",
        starting_material="poly(A) RNA", fragmentation="Tn5 / tagmentation",
        adapter_installation="template switching|Tn5|PCR-added",
        amplification="whole-transcriptome amplification|endpoint PCR",
        identifiers="UMI|cell barcode|sample index",
        read_structure="paired-end|index reads|barcode in Read 1"),
    "combi_scrna": _from(
        _ILLUMINA, index_introduction="reverse transcription|ligation|PCR|combinatorial indexing",
        index_architecture="custom Illumina|inline barcode|dual index", assay="RNA-seq",
        partitioning="combinatorial indexing|single cell", starting_material="fixed cells / nuclei",
        fragmentation="no fragmentation", adapter_installation="ligation|PCR-added",
        amplification="endpoint PCR", identifiers="UMI|cell barcode|sample index",
        read_structure="paired-end|index reads|inline barcode"),
    "atac_plate": _from(
        _ILLUMINA, index_introduction="transposition|PCR", index_architecture="Nextera|dual index",
        assay="ATAC / accessibility", partitioning="plate / well|single cell",
        starting_material="chromatin", fragmentation="Tn5 / tagmentation",
        adapter_installation="Tn5|PCR-added"),
    "atac_droplet": _from(
        _ILLUMINA, index_introduction="transposition|extension / fill-in|PCR",
        index_architecture="Nextera|inline barcode|dual index", assay="ATAC / accessibility",
        partitioning="droplet|single nucleus", starting_material="chromatin",
        fragmentation="Tn5 / tagmentation", adapter_installation="Tn5|PCR-added",
        identifiers="cell barcode|sample index|no UMI",
        read_structure="paired-end|index reads|barcode in Read 1"),
    "atac_combi": _from(
        _ILLUMINA, index_introduction="transposition|PCR|combinatorial indexing",
        index_architecture="Nextera|inline barcode|dual index", assay="ATAC / accessibility",
        partitioning="combinatorial indexing|single nucleus", starting_material="chromatin",
        fragmentation="Tn5 / tagmentation", adapter_installation="Tn5|PCR-added",
        identifiers="cell barcode|sample index|no UMI",
        read_structure="paired-end|index reads|inline barcode"),
    "hic_bulk": _from(
        _ILLUMINA, assay="chromatin conformation", starting_material="chromatin",
        fragmentation="restriction digest", adapter_installation="ligation",
        topology="linear|proximity-ligation junction", selection="affinity / antibody capture"),
    "hic_single": _from(
        _ILLUMINA, assay="chromatin conformation", partitioning="single nucleus|single cell",
        starting_material="chromatin", fragmentation="restriction digest",
        adapter_installation="ligation", amplification="whole-genome amplification|endpoint PCR",
        topology="linear|proximity-ligation junction", identifiers="sample index|no UMI"),
    "methyl_dna": _from(
        _ILLUMINA, assay="DNA methylation", fragmentation="no fragmentation",
        conversion="bisulfite"),
    "crispr_screen": _from(
        _ILLUMINA, index_introduction="PCR", index_architecture="custom Illumina|single index",
        assay="CRISPR screening", starting_material="genomic DNA", fragmentation="no fragmentation",
        adapter_installation="PCR-added", identifiers="sample index|guide barcode|no UMI",
        read_structure="single-end|index reads|custom sequencing primer",
        selection="guide-specific enrichment"),
    "pacbio": _from(
        _ILLUMINA, index_introduction="pre-indexed adaptor", index_architecture="PacBio barcode",
        platform="PacBio", fragmentation="mechanical", adapter_installation="hairpin / dumbbell",
        amplification="PCR-free", topology="dumbbell / SMRTbell", identifiers="sample index|no UMI",
        read_structure="long-read consensus"),
    "nanopore": _from(
        _ILLUMINA, index_introduction="ligation|pre-indexed adaptor",
        index_architecture="Nanopore barcode", platform="Oxford Nanopore",
        fragmentation="no fragmentation", adapter_installation="ligation|full-length adaptor",
        amplification="PCR-free", identifiers="sample index|no UMI",
        read_structure="single-end"),
    "in_situ": _from(
        _ILLUMINA, index_introduction="reverse transcription", index_architecture="inline barcode",
        assay="RNA-seq|spatial / in situ", platform="in situ imaging",
        partitioning="spatial / in situ", starting_material="total RNA",
        fragmentation="no fragmentation", adapter_installation="padlock / circularization",
        amplification="rolling-circle amplification", topology="circular|concatemer",
        identifiers="spatial barcode|no UMI", read_structure="single-end|inline barcode"),
    "bdna_imaging": _from(
        _ILLUMINA, index_introduction="no index", index_architecture="no barcode",
        assay="RNA detection|spatial / in situ", platform="in situ imaging",
        partitioning="spatial / in situ", starting_material="total RNA",
        fragmentation="no fragmentation", adapter_installation="hybridization scaffold",
        amplification="branched-DNA signal amplification",
        topology="branched hybridization tree", strand_handling="RNA hybridization",
        identifiers="no UMI", read_structure="microscopy", selection="none",
        availability="commercial kit|proprietary sequences|schematic complete"),
    "bdna_flow": _from(
        _ILLUMINA, index_introduction="no index", index_architecture="no barcode",
        assay="RNA detection|protein / feature detection|multiomic", platform="flow cytometry",
        partitioning="single cell", starting_material="total RNA",
        fragmentation="no fragmentation", adapter_installation="hybridization scaffold",
        amplification="branched-DNA signal amplification",
        topology="branched hybridization tree", strand_handling="RNA hybridization",
        identifiers="no UMI", read_structure="flow cytometry", selection="none",
        availability="commercial kit|proprietary sequences|schematic complete"),
}


def _values(cell: str) -> tuple[str, ...]:
    return tuple(x.strip() for x in cell.split("|") if x.strip())


def _overrides(cell: str, path: Path, line: int) -> dict[str, tuple[str, ...]]:
    out = {}
    if not cell.strip():
        return out
    for assignment in cell.split(";"):
        if "=" not in assignment:
            raise ValueError(f"{path}:{line}: malformed override {assignment!r}")
        key, value = (x.strip() for x in assignment.split("=", 1))
        if key not in FACETS or key in out:
            raise ValueError(f"{path}:{line}: unknown or duplicate facet {key!r}")
        out[key] = _values(value)
    return out


def load(expected_dirs: set[str] | None = None, path: Path = TSV) -> dict[str, dict[str, tuple[str, ...]]]:
    """Load complete protocol properties, rejecting ambiguity before rendering the site."""
    with path.open(encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames != COLUMNS:
            raise ValueError(f"{path}: columns are {reader.fieldnames}, expected {COLUMNS}")
        rows = list(reader)
    out: dict[str, dict[str, tuple[str, ...]]] = {}
    for line, row in enumerate(rows, 2):
        d = row["dir"]
        if not d or d in out:
            raise ValueError(f"{path}:{line}: missing or duplicate directory {d!r}")
        if row["profile"] not in PROFILES:
            raise ValueError(f"{path}:{line}: {d}: unknown profile {row['profile']!r}")
        rec = dict(PROFILES[row["profile"]])
        rec.update(_overrides(row["overrides"], path, line))
        for key, (_, allowed) in FACETS.items():
            values = rec[key]
            if not values:
                raise ValueError(f"{path}:{line}: {d}: empty {key}")
            if len(values) != len(set(values)):
                raise ValueError(f"{path}:{line}: {d}: duplicate {key} value")
            bad = set(values) - set(allowed)
            if bad:
                raise ValueError(f"{path}:{line}: {d}: unknown {key}: {sorted(bad)}")
            rec[key] = values
        out[d] = rec
    if expected_dirs is not None and set(out) != expected_dirs:
        missing, extra = expected_dirs - set(out), set(out) - expected_dirs
        raise ValueError(f"{path}: directory mismatch; missing={sorted(missing)}, extra={sorted(extra)}")
    return out


def option_counts(records: list[dict]) -> dict[str, dict[str, int]]:
    """Counts shown beside filters, derived from exactly the records being rendered."""
    return {key: {value: sum(value in r["properties"][key] for r in records)
                  for value in allowed}
            for key, (_, allowed) in FACETS.items()}
