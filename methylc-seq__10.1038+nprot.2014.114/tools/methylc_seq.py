"""Conventional ligation-first MethylC-seq / WGBS model."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from chemdraw import Row
from chromatin_epigenetics import illumina_ligation_library, ligated_insert_scene

TITLE = "MethylC-seq — conventional whole-genome bisulfite sequencing"
NOTES = "01_methylc-seq.html"
SOURCE = 'Protocol source: <a href="https://doi.org/10.1038/nprot.2014.114">Urich et al., <i>Nature Protocols</i> (2015)</a>.'
SUMMARY = "Fragment genomic DNA, repair and dA-tail it, ligate fully methylated adapters, then denature and bisulfite-convert before limited uracil-tolerant PCR."
CAVEAT = "The protocol requires fully methylated adapter cytosines but does not make an ordinary canonical adapter sequence itself a scientific result. Canonical TruSeq arms are therefore shown as inferred."
FINAL_LIBRARY, SEQ_PRIMERS = illumina_ligation_library("C→T-converted genomic insert", inferred=True)
FINAL_CAPTION = "INFERRED — PCR-restored duplex library. Unmethylated genomic C is read as T; protected 5mC remains C. Dotted adapter bases denote the canonical inferred TruSeq layout."
SEQUENCING_INTRO = "Paired reads traverse the converted insert. Mapping compares the C/T outcome against the unconverted reference sequence."

def sections():
    return [
        ("Fragment, end-repair and dA-tail", [Row(chunks=[("genomic DNA → ~200-bp fragments → blunt 5′-phosphorylated ends → 3′ dA", None, False)])],
         "Mechanical fragmentation is followed by conventional end preparation."),
        ("Ligate protected adapters before conversion", ligated_insert_scene("C→T-converted genomic insert").rows(),
         "INFERRED — canonical TruSeq placement represents the paper's fully methylated Y-adapter; ** marks both ligations."),
        ("Convert and amplify", [Row(chunks=[("unmethylated C → U → T in PCR     5mC → C", "w1", False)]),
                                  Row(chunks=[("adapter-ligated duplex → denatured non-complementary single strands → PCR duplex", None, False)])],
         "Bisulfite conversion occurs after adapter ligation and can break already tagged molecules."),
    ]
