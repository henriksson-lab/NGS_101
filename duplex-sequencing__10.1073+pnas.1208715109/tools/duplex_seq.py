"""Original Schmitt et al. Duplex Sequencing library model."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from batch_ngs import duplex_family_rows, seg, truseq_library
from chemdraw import Row, Scene

TITLE = "Duplex Sequencing — complementary-strand consensus"
NOTES = "01_duplex-sequencing.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1073/pnas.1208715109">Schmitt et al., <i>PNAS</i> (2012)</a>.'
SUMMARY = "Complementary double-stranded random tags identify both strands of each original DNA duplex. PCR descendants first form two single-strand consensuses, which are paired by their swapped end tags to form a duplex consensus."
CAVEAT = "α and β denote the two independently random 12-nt end tags on one original fragment. They are molecule identifiers, not fixed sequences; complementary strands carry the reciprocal αβ/βα relationship by construction."

FINAL_LIBRARY, SEQ_PRIMERS = truseq_library([
    seg("duplex tag α", "N" * 12, "umi", placeholder=True),
    seg("genomic insert", "X" * 38, placeholder=True),
    seg("duplex tag β", "N" * 12, "umi", placeholder=True)],
    "Duplex Sequencing library", dual_index=False)
FINAL_CAPTION = "Each end contributes a 12-nt duplex tag, giving a 24-nt αβ identity for one strand family. The complementary family is recognized by the reciprocal βα tag relationship."
SEQUENCING_INTRO = "Paired-end sequencing observes the two end tags and the intervening insert. Consensus construction is part of the molecular design: strand families are linked only when their tags are complementary and swapped."

def sections():
    adapter = [seg("asymmetric PCR arm", "X" * 12, placeholder=True),
               seg("12-nt random tag", "N" * 12, "umi", placeholder=True),
               seg("3-prime dA", "A")]
    return [
        ("Synthesize a double-stranded random-tag adapter",
         Scene.duplex(adapter, label="Duplex Tag adapter").rows(),
         "A single-stranded N12 region is copied to make a complementary duplex tag. Extended polymerase incubation completes the adapter’s 3′ dA."),
        ("Ligate tagged adapters to both ends of T-tailed DNA",
         [Row(chunks=[("adapter—N12(α)—A:T ** genomic duplex ** T:A—N12(β)—adapter", "umi", False)])],
         "The original implementation T-tailed the sheared insert and A-tailed the adapters. ** marks the sealed adapter–insert boundaries."),
        ("Build paired strand families and require duplex agreement",
         duplex_family_rows(),
         "Reads sharing αβ form one SSCS; reads sharing the reciprocal βα tags form its complementary SSCS. Only agreeing bases survive into the DCS."),
    ]
