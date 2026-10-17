"""
Pooled CRISPR screening: the shared sequence landmarks and guide cloning.

Sequences here were extracted from real Addgene maps in lenticrispr-v1-screening__10.1126+science.1247005/ref/plasmids/ using
lib/plasmid.py, not retyped from papers. The cloning simulation below is the check that
matters: digest a real vector, ligate a guide, and confirm the product reconstitutes
U6 -> spacer -> scaffold contiguously.
"""

from __future__ import annotations

from chemdraw import revcomp
from plasmid import BsmBI, Plasmid, digest

# ----------------------------------------------------------------- landmarks
# The human U6 promoter's 3' end, NOT including the +1 base.
U6_3PRIME = "GACGAAACACC"
# The +1 transcription start. It is supplied by the spacer when the spacer already begins
# with G, and appended otherwise -- so U6_3PRIME + G + spacer double-counts a leading G.
# That off-by-one is why the two are separate constants.
U6_PLUS1 = "G"
U6_3PRIME_WITH_G = U6_3PRIME + U6_PLUS1       # the usual "GACGAAACACCG" search motif
# The 23 nt of human U6 immediately 5' of U6_3PRIME (lentiGuide-Puro / lentiCRISPR maps).
# Readout primers anneal across U6_UPSTREAM + U6_3PRIME at various lengths.
U6_UPSTREAM = "GGCTTTATATATCTTGTGGAAAG"

# The original Cong/Ran scaffold (lentiCRISPR v1, v2, lentiGuide-Puro). 76 nt.
SCAFFOLD_V1 = ("GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTT"
               "GAAAAAGTGGCACCGAGTCGGTGC")

# Diagnostic prefix of the v1 scaffold in a map or a read.
SCAFFOLD_V1_HEAD = SCAFFOLD_V1[:23]      # GTTTTAGAGCTAGAAATAGCAAG

# What BsmBI leaves behind in the Zhang-lineage backbones, verified by digestion.
OVERHANG_UPSTREAM = "CACC"               # insert supplies this on its top strand
OVERHANG_DOWNSTREAM = "GTTT"             # insert supplies revcomp(=AAAC) on its bottom
FILLER_LEN = 1885                        # bp between the two BsmBI nicks

SPACER_LEN = 20


# ------------------------------------------------------------ guide oligos
def guide_oligos(spacer: str) -> tuple[str, str]:
    """The annealed oligo pair for cloning a spacer into a BsmBI-cut lentiCRISPR vector.

    Returns (top, bottom), both 5'->3'. U6 needs a G at +1, so a spacer that does not
    already start with G gets one prepended -- the standard Zhang-lab convention.
    """
    spacer = spacer.upper()
    if spacer.startswith("G"):
        return OVERHANG_UPSTREAM + spacer, revcomp(OVERHANG_DOWNSTREAM) + revcomp(spacer)
    return (OVERHANG_UPSTREAM + "G" + spacer,
            revcomp(OVERHANG_DOWNSTREAM) + revcomp(spacer) + "C")


def clone_guide(p: Plasmid, spacer: str, enz=BsmBI) -> Plasmid:
    """Simulate the digest-and-ligate, returning the recombinant sequence.

    The backbone is everything OUTSIDE the two cuts (wrapping the origin); the filler
    between them is replaced by G + spacer.
    """
    cuts = digest(p, enz)
    if len(cuts) != 2:
        raise ValueError(f"{p.name}: expected 2 {enz.name} sites, found {len(cuts)}")
    up, down = cuts
    n = len(p)
    backbone = p.sub(down.top_cut, up.top_cut + n)       # the long way round
    spacer = spacer.upper()
    insert = (OVERHANG_UPSTREAM + ("" if spacer.startswith("G") else U6_PLUS1) + spacer)
    return Plasmid(name=f"{p.name}+{spacer}", seq=backbone + insert,
                   circular=True, features=[])


def guide_context(recombinant, spacer: str) -> bool:
    """Does the clone reconstitute U6 -> spacer -> scaffold without a seam?"""
    recombinant = recombinant.seq if isinstance(recombinant, Plasmid) else recombinant
    spacer = spacer.upper()
    g = "" if spacer.startswith("G") else U6_PLUS1
    expect = U6_3PRIME + g + spacer + SCAFFOLD_V1_HEAD
    doubled = recombinant + recombinant
    return expect in doubled
