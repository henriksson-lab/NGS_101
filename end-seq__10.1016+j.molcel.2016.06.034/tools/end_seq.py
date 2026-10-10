"""Canela et al. END-seq model."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from batch_ngs import joined_scene, seg, truseq_library
from chemdraw import Row
from end_capture import UserHairpin, captured_break_scene, user_opening_rows

TITLE = "END-seq — direct capture of DNA double-strand-break ends"
NOTES = "01_end-seq.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1016/j.molcel.2016.06.034">Canela et al., <i>Molecular Cell</i> (2016)</a>.'
SUMMARY = "Native DNA ends are blunted and trapped in agarose with a biotinylated P5 hairpin; after shearing and enrichment, a P7 hairpin captures the opposite library end and USER opens both loops for indexed PCR."
CAVEAT = "Blunting deliberately erases the original overhang chemistry. Read 1 starts at the first base retained at the processed break end."
A1 = UserHairpin("ENDseq adaptor 1", "GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGUUTUTUUACACTCTTTCCCTACACGACGCTCTTCCGATCT", (33,35))
A2 = UserHairpin("ENDseq adaptor 2", "GATCGGAAGAGCACACGTCUUUUUUUUAGACGTGTGCTCTTCCGATCT")
FINAL_LIBRARY, SEQ_PRIMERS = truseq_library([seg("break-derived genomic insert", "X"*42, placeholder=True)], "END-seq library", inferred_adapters=False)
FINAL_CAPTION = "USER-opened hairpin arms are completed by TruSeq index PCR; the captured DSB boundary is adjacent to the Read 1 arm."
SEQUENCING_INTRO = "Single-end END-seq uses Read 1 to enter the captured genomic end; the complete construct also shows the standard TruSeq index and opposite-end primer sites created by PCR."

def sections():
    captured = [seg("P5 hairpin", "X"*20, "r1", placeholder=True),
                seg("break-derived fragment", "X"*34, placeholder=True),
                seg("P7 hairpin", "X"*20, "r2", placeholder=True)]
    return [
        ("Deproteinize DNA in agarose, blunt breaks and dA-tail", [Row(chunks=[("genomic DNA ---- processed DSB end-A", None, False)])],
         "Embedding before manipulation limits adventitious shearing; polymerases make the captured end blunt before one dA is added."),
        ("Ligate biotinylated P5 hairpin adaptor 1 in the plug", captured_break_scene().rows(),
         "The 3-prime T overhang accepts the break-end dA; two biotin-dU residues later support streptavidin capture and USER opening."),
        ("Shear, enrich captured molecules and ligate P7 hairpin adaptor 2", joined_scene(captured, (("P5 hairpin", "break-derived fragment", "adaptor 1 ligation"), ("break-derived fragment", "P7 hairpin", "adaptor 2 ligation")), label="two-hairpin captured fragment").rows(),
         "Sonication creates the distal end; its repair and dA-tail precede the second hairpin ligation on streptavidin beads."),
        ("Open both hairpins with USER and index-amplify", user_opening_rows(A1,A2),
         "USER cleavage at deoxyuridines exposes the two PCR-addressable Illumina arms."),
    ]
