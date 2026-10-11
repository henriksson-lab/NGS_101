"""Sequencing readouts of distance-dependent covalent labelling in fixed cells."""
from __future__ import annotations

from types import SimpleNamespace

import seqprimers as sp
from batch_ngs import seg, truseq_library
from protocol_model import chain, molecule, role


def tsa_library():
    """Compatible Kapa Hyper Prep shell; undisclosed kit oligos stay inferred."""
    lib, all_primers = truseq_library([
        seg("TSA-biotin-labelled genomic fragment", "X" * 46, placeholder=True),
    ], "TSA-seq pulldown library", dual_index=False,
       inferred_adapters=True, i7_nt=6)
    primers = tuple(p for p in all_primers if p.role != "Read 2")
    problems = sp.verify(lib, primers, required_roles=tuple(p.role for p in primers))
    if problems:
        raise ValueError("TSA-seq library: " + "; ".join(problems))
    return lib, primers


def tsa_workflow():
    return chain(molecule("formaldehyde-fixed, permeabilized cells", (
        role("target nuclear compartment"), role("chromosomes at varying distances"))), (
        ("Immunostain the compartment marker and recruit horseradish peroxidase",
         molecule("localized HRP source in an intact nucleus", (
             role("antibody-marked nuclear compartment"), role("HRP at the target"),
             role("surrounding chromosomal DNA"))),
         "The defining demonstration used SON for nuclear speckles and lamin antibodies for the nuclear lamina; a no-primary sample controls nonspecific labelling."),
        ("Generate short-lived tyramide-biotin radicals around the HRP source",
         molecule("distance-decaying DNA labelling field", (
             role("HRP source"), role("strongly biotin-labelled nearby DNA"),
             role("weakly labelled distant DNA"))),
         "Radical concentration and therefore direct DNA biotinylation decay with cytological distance. This is a spatial concentration gradient, not a DNA contact or ligation junction."),
        ("Reverse crosslinks, purify DNA and sonicate to 200–800 bp",
         molecule("sheared labelled genomic DNA", (
             role("random sonication end"), role("genomic DNA carrying sparse biotin"),
             role("random sonication end"))),
         "Most fragments were about 500 bp. Ten percent input was retained before affinity purification."),
        ("Capture biotinylated fragments on streptavidin beads and release DNA",
         molecule("TSA-seq pulldown DNA", (
             role("random sonication end"), role("biotin-enriched genomic fragment"),
             role("random sonication end"))),
         "Stringent washes enrich fragments in proportion to their direct tyramide-biotin labelling; pulldown coverage is compared with the matched input."),
        ("Blunt-end, dA-tail, ligate indexed Kapa adapters and PCR for 8–16 cycles",
         molecule("single-index TSA-seq library", (
             role("P5 and Read 1 arm", inferred=True),
             role("biotin-enriched genomic insert"),
             role("Read 2 arm, six-base sample index and P7", inferred=True))),
         "The defining paper names Kapa Hyper Prep and six-base indexed adapters but does not publish the kit oligo sequences."),
    ))


def get():
    final, primers = tsa_library()
    return SimpleNamespace(
        TITLE="TSA-seq — distance mapping from a nuclear compartment",
        NOTES="01_tsa_seq.html",
        SOURCE='<a href="https://doi.org/10.1083/jcb.201807108">Chen et al., <i>Journal of Cell Biology</i> (2018)</a>.',
        SUMMARY="Antibody-localized HRP generates a tyramide-biotin radical gradient in fixed nuclei. Nearby DNA is labelled more strongly than distant DNA, then sonicated, captured with streptavidin and sequenced against matched input.",
        CAVEAT="TSA-seq does not encode distance in a molecule-level barcode or junction: distance is inferred from population enrichment. <inf>The final drawing uses compatible standard Illumina geometry for the Kapa Hyper Prep adapters because their exact six-base indexed oligos were not published.</inf>",
        oligos=lambda: (),
        workflow=tsa_workflow,
        FINAL_LIBRARY=final,
        FINAL_CAPTION="One enriched genomic fragment inside the compatible single-index Kapa/Illumina shell. Biotin enabled pulldown before library construction and is not a sequenced base. <inf>Adapter bases show compatible geometry rather than a source transcription of the proprietary kit oligos.</inf>",
        SEQ_PRIMERS=primers,
        READ_LENGTHS={"Read 1": 100, "Index 1 (i7)": 6},
        SEQUENCING_INTRO="The defining libraries were sequenced for 101 instrument cycles and yielded 100-nt single-end genomic reads plus a six-base sample index. Read density in streptavidin pulldown versus input supplies the distance-dependent signal.",
        SEQUENCING_ENDING="",
        SEQUENCING_UNAVAILABLE="",
    )
