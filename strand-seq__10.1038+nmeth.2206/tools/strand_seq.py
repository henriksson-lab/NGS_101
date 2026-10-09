"""Falconer et al. Strand-seq parental-template and historical PE library model."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import illumina as il
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, feature, revcomp
from chromatin_epigenetics import strand_selection_rows


def seg(name, top, tag=None, **kw):
    return Segment(name, top, tag, **kw)


TITLE = "Strand-seq — parental template-strand sequencing"
NOTES = "01_strand-seq.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1038/nmeth.2206">Falconer et al., <i>Nature Methods</i> (2012)</a>.'
SUMMARY = "Cells replicate once in BrdU; after single-cell isolation, Hoechst/UV selectively nicks the BrdU-containing nascent strands so only the parental template strands amplify."
CAVEAT = "The paper prints the custom indexed PCR primer and index-read primer and names Illumina PE adapters and PE&nbsp;1.0. The adapter and standard read-primer bases below come from Illumina’s authoritative obsolete paired-end oligo table."

# Falconer Online Methods. The apparent seven-N run is explicitly described as a
# fault-tolerant hexamer. The seven-cycle index read therefore continues into P7'.
INDEX_NT = 6
PE_ADAPTER_BOTTOM = il.TRUSEQ_READ1
INDEX_READ = "GATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCG"
# The printed index-read primer spans the adapter oligo and four adjacent bases.
# The old Read 2 primer is its reverse complement with the dA-junction base included.
PE_ADAPTER_TOP = INDEX_READ[:-4]
PE_READ2 = revcomp("A" + INDEX_READ)


def adapter_scene() -> Scene:
    sc = Scene()
    sc.strand("bottom", [
        seg("Read 1 fork arm", PE_ADAPTER_BOTTOM[:-13], "r1"),
        seg("PE stem complement", PE_ADAPTER_BOTTOM[-13:-1], "r1"),
        seg("3′ T", PE_ADAPTER_BOTTOM[-1], "r1"),
    ], label="PE adapter oligo")
    sc.anneal("top", [
        seg("PE stem", PE_ADAPTER_TOP[:12], "r2"),
        seg("PE fork arm", PE_ADAPTER_TOP[12:], "r2"),
    ], to="bottom", pair=("PE stem", "PE stem complement"),
              label="PE adapter oligo", mod5="p", unpaired=("PE fork arm",))
    sc.mark("bottom", "3′ T", "ligation overhang")
    return sc


def pe1_segments():
    return [seg("P5", il.P5, "p5"),
            seg("Read 1 remainder", il.TRUSEQ_READ1[4:], "r1")]


def indexed_primer_segments():
    return [seg("P7", il.P7, "p7"),
            seg("hexamer index", "N" * INDEX_NT, "cbc", placeholder=True),
            seg("PE Read 2 arm", PE_READ2, "r2")]


SEQ_PRIMERS = (
    sp.custom("Read 1", "obsolete PE Read 1", il.TRUSEQ_READ1,
              'Illumina "Oligonucleotide Sequences for Paired End DNA" (obsolete)'),
    sp.custom("Index 1 (i7)", "Strand-seq custom index primer", INDEX_READ,
              "Falconer et al. 2012 Online Methods",
              "Seven cycles report the six-base index and the adjacent first base of P7'."),
    sp.custom("Read 2", "obsolete PE Read 2", PE_READ2,
              'Illumina "Oligonucleotide Sequences for Paired End DNA" (obsolete)'),
)


def final_library() -> Construct:
    lib = Construct([
        *pe1_segments(),
        seg("retained parental template fragment", "X" * 36, placeholder=True),
        seg("PE Read 2 site", revcomp(PE_READ2), "r2"),
        seg("i7 reverse complement", "N" * INDEX_NT, "cbc", placeholder=True,
            feature=feature("cell_i7", "cell_barcode", "whitelist", whitelist="published hexamer index set")),
        seg("P7 reverse complement", il.P7_RC, "p7"),
    ], name="Strand-seq library")
    problems = sp.verify(lib, SEQ_PRIMERS,
                         required_roles=("Read 1", "Index 1 (i7)", "Read 2"))
    if problems:
        raise ValueError("invalid Strand-seq library: " + "; ".join(problems))
    return lib


FINAL_LIBRARY = final_library()
FINAL_CAPTION = "PCR-completed, single-index historical paired-end library. The six-base i7 identifies the cell; there is no i5 index."
SEQUENCING_INTRO = "Paired 76-nt reads used the obsolete Illumina PE primers. A third seven-cycle read used the printed custom primer to report the six-base barcode plus one adjacent fixed base."


def sections():
    return [
        ("Label nascent DNA during one S phase", strand_selection_rows()[:2],
         "Semiconservative replication leaves each chromatid with one parental template and one BrdU-containing newly synthesized strand."),
        ("Ligate the historical Illumina PE adapter", adapter_scene().rows(),
         "The 5′-phosphorylated PE oligo and Read 1 oligo form the adapter named by the defining paper; ligation precedes strand destruction."),
        ("Destroy the nascent strand before PCR", strand_selection_rows()[2:],
         "Hoechst sensitization and UV nick BrdU DNA. PE 1.0 and the custom indexed primer then amplify only intact parental templates."),
    ]
