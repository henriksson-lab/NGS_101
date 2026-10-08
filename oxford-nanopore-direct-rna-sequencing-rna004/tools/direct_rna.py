"""Oxford Nanopore SQK-RNA004 direct-RNA library architecture."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from batch_ngs import seg
from chemdraw import Row, Scene

TITLE = "Oxford Nanopore Direct RNA Sequencing — SQK-RNA004"
NOTES = "01_direct-rna.html"
SOURCE = 'Commercial protocol: Oxford Nanopore <a href="https://nanoporetech.com/document/direct-rna-sequencing-sqk-rna004">Direct RNA Sequencing Kit SQK-RNA004</a>; architecture paper: <a href="https://doi.org/10.1038/nmeth.4577">Garalde et al. (2018)</a>.'
SUMMARY = "A duplex reverse-transcription adapter is ligated to a poly(A) RNA 3′ end, primes a stabilizing cDNA, and receives a motor-loaded sequencing adapter. The nanopore reads the native RNA strand from its 3′ end toward 5′."
CAVEAT = "Current RTA and sequencing-adapter sequences, motor identity and attachment chemistry are proprietary. Bracketed components below are structural roles, not invented bases."
FINAL_LIBRARY = None
SEQ_PRIMERS = ()
SEQUENCING_ENDING = "No sequencing primer extends this library. The motor-loaded sequencing adapter captures the RNA–cDNA hybrid at the pore and controls passage of the native RNA strand, beginning at its poly(A)/3′ end and proceeding toward the RNA 5′ end. The cDNA stabilizes the RNA but is not the strand base-called in direct-RNA mode."

def sections():
    mrna = [seg("RNA body", "X" * 30, placeholder=True), seg("poly(A)", "A" * 18)]
    rta = [seg("RTA duplex / RT handle", "X" * 14, "r2", placeholder=True,
               inferred=True), seg("oligo(dT)", "T" * 10)]
    sc = Scene(); sc.strand("RNA", mrna, label="native RNA", mod5="")
    sc.anneal("RTA", rta, to="RNA", pair=("oligo(dT)", "poly(A)"),
              label="reverse-transcription adapter", unpaired=("RTA duplex / RT handle",))
    sc.arrow("RTA", "reverse transcription stabilizes RNA")
    return [
        ("Ligate the reverse-transcription adapter to the RNA 3′ end", sc.rows(),
         "The RTA’s oligo(dT) region anneals to poly(A), positioning ligation at the RNA 3′ end and providing the primer/handle for optional reverse transcription."),
        ("Synthesize the stabilizing cDNA strand",
         [Row(chunks=[("native RNA     5′—RNA body—poly(A)—RTA—3′", None, False)]),
          Row(chunks=[("                    ||||||||||||||", None, False)]),
          Row(chunks=[("stabilizing cDNA 3′←——————————————5′", "r2", False)])],
         "Reverse transcription reduces structure-related pore blockage. Only the RNA strand is sequenced."),
        ("Attach the motor-loaded sequencing adapter",
         [Row(chunks=[("[motor + membrane tether] ** [sequencing adapter] ** RTA—RNA:cDNA", "me", False)]),
          Row(chunks=[("                         pore entry → RNA 3′ to 5′", None, False)])],
         "** marks vendor-specified ligation/assembly boundaries whose molecular sequences are undisclosed."),
    ]
