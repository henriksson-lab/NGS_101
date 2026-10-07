"""
CRISPR-UMI (Schmierer) -- a Random Sequence Label carried in the guide library itself.

Schmierer B, Botla SK, Zhang J, Turunen M, Kivioja T, Taipale J.
"CRISPR/Cas9 screening using unique molecular identifiers."
Mol Syst Biol 2017;13(10):945. doi:10.15252/msb.20177834  (PMC5658704, CC BY)

Not to be confused with CRISPR-MIP (`../crispr-mip__10.1101+2024.03.28.587082/`), which also puts a UMI in a CRISPR
screen readout. The two are opposite in design and the names collide:

    Schmierer      UMI is cloned INTO the library, one per virus. It labels a
                   TRANSDUCED CELL LINEAGE, and is read out of genomic DNA. Counting
                   distinct RSLs counts clones.
    CRISPR-MIP     UMI is carried on a capture PROBE applied to gDNA afterwards. It
                   labels a CAPTURED MOLECULE. Counting distinct UMIs counts molecules.

Both are "CRISPR + UMI"; only one of them changes the plasmid.
"""

from __future__ import annotations

import sys
from pathlib import Path

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE.parents[1] / "lib"))
# Appended, not inserted: crispr-mip__10.1101+2024.03.28.587082/tools has a build_page.py of its own, and putting it in
# front of sys.path would let `import build_page` pick up the wrong page's module.
sys.path.append(str(_HERE.parents[1] / "crispr-mip__10.1101+2024.03.28.587082" / "tools"))

import crisprmip as cm
import illumina as il
import seqprimers as sp
from chemdraw import Construct, Segment, revcomp
import crispr
from crispr import SCAFFOLD_V1, SPACER_LEN, U6_3PRIME_WITH_G, clone_guide
from plasmid import Plasmid, read_genbank

# ------------------------------------------------------------------- the vector
# pLenti-Puro-AU-flip-3xBsmBI: lentiGuide-Puro (Addgene #52963) with two edits, quoted
# verbatim from Materials and Methods:
#   "replacing the sequence gttttagagctagaaatagcaagttaaaa......TTTTTT with
#    gtttAagagctagaaatagcaagttTaaa......TTTTTTcgtctct to create an AU-flip
#    (Chen et al, 2013) and an additional BsmBI site downstream of the tracrRNA"
PARENT_VECTOR = "lentiGuide-Puro (Addgene #52963)"
VECTOR_NAME = "pLenti-Puro-AU-flip-3xBsmBI"

# The AU-flip is two substitutions in the scaffold's first 29 nt.
SCAFFOLD_AU_FLIP = (SCAFFOLD_V1.replace("GTTTT", "GTTTA", 1)[:25] + "T"
                    + SCAFFOLD_V1.replace("GTTTT", "GTTTA", 1)[26:])
EXTRA_BSMBI = "CGTCTCT"            # added after the terminator

# ------------------------------------------------------------- the library insert
# "The full insert sequence is ..." -- the array fragment plus a 119-bp overlapping oligo
# that carries the RSL and the Illumina i7 primer site.
# The oligo's U6 overlap. Its 3' end is lib's U6_3PRIME_WITH_G (GACGAAACACC + the Pol III
# +1 G); the 23 nt upstream are not in lib.
U6_UPSTREAM = crispr.U6_UPSTREAM
U6_OVERLAP = U6_UPSTREAM + U6_3PRIME_WITH_G              # ends on the Pol III +1 G
TERMINATOR = "TTTTTT"          # 6 T, as in lentiGuide-Puro (verified: gives the published 288 bp)
# Note what this is: the Illumina Read-2 / Index-1 adapter, built INTO the construct. It is
# the canonical 33-nt form -- the same oligo Illumina calls the Index-1 sequencing primer --
# and NOT the 34-nt revcomp of TruSeq Read 2, which carries one extra 5' A. That missing A
# is why the stock Read-2 primer has no site here; see SEQ_PRIMERS.
ILLUMINA_ADAPTER = il.INDEX1_PRIMER
RSL_LEN = 6                                              # the "Random Sequence Label"
DOWNSTREAM = "AAGCTTGGCGTAACTAGATCTTGAGACAAA"            # lentiGuide-Puro, past the HindIII
DOWNSTREAM_PCR = "TGGCAG"      # the next 6 vector nt; PCR2-R / PCR3-R's 3' end sits here
DEMO_SPACER = "ATCG" * 5       # a stand-in guide, used wherever a page or test needs one


def insert(spacer: str = "N" * SPACER_LEN, rsl: str = "N" * RSL_LEN) -> str:
    """The full library insert, 5'->3', as published."""
    return (U6_OVERLAP + spacer + SCAFFOLD_AU_FLIP + TERMINATOR
            + ILLUMINA_ADAPTER + rsl + DOWNSTREAM)


# The vector, rebuilt from the real parent map. One definition, used by the page, the
# selftest and the padlock design script alike.
#
# The map is third-party (Addgene's / the depositor's annotation) and so is not committed
# here; see ref/MANIFEST.md. Two files carry the same plasmid and either will do -- they are
# the same 10,183 bp at different origins, which is immaterial because every measurement
# below is made with find_both()/amplify() on a circle:
#   1. the Addgene GenBank deposit kept with the lentiCRISPR screen notes -- the citable one;
#   2. a SnapGene .dna of the same plasmid, if a local copy happens to be around.
PARENT_MAP_ORIGIN = ("Addgene #52963 (lentiGuide-Puro), full sequence -- see "
                     "ref/MANIFEST.md in this directory")
PARENT_MAP_CANDIDATES = (
    _HERE.parents[1] / "lenticrispr-gecko-screen__10.1126+science.1247005" / "ref"
    / "plasmids" / "addgene-52963_lentiGuide-Puro.gb",
    _HERE.parents[1] / "to_debug" / "lentiGuide-Puro.dna",
)


class MissingMap(FileNotFoundError):
    """The parent plasmid map is not on disk.

    Deliberately an exception and not a silent fallback: without the real map nothing here
    can be computed, and a guessed sequence would be worse than no answer. Catchable, so a
    page can refuse to build with one clear line and the selftest can skip (not fail) the
    checks that need it.
    """


def parent_map() -> Path | None:
    """The first candidate map that is actually on disk, or None on a fresh clone."""
    for p in PARENT_MAP_CANDIDATES:
        if p.exists():
            return p
    return None


# Kept as a module attribute because the page and the selftest both name it. Resolved at
# import time, but nothing is read: import never touches the filesystem beyond exists().
PARENT_MAP = parent_map() or PARENT_MAP_CANDIDATES[0]


def parent_vector() -> Plasmid:
    """lentiGuide-Puro, from the Addgene map. Raises MissingMap if the map is not here."""
    p = parent_map()
    if p is None:
        raise MissingMap(f"missing {PARENT_MAP_CANDIDATES[0].name} -- obtain it from: "
                         f"{PARENT_MAP_ORIGIN}")
    return read_genbank(str(p))


def rebuild_vector(cloned: bool = False, spacer: str = DEMO_SPACER) -> Plasmid:
    """pLenti-Puro-AU-flip-3xBsmBI: the parent map with the published edits spliced in."""
    p = parent_vector()
    i = p.seq.find(SCAFFOLD_V1 + TERMINATOR)
    if i < 0:
        raise ValueError(f"{PARENT_MAP.name} does not carry the original scaffold")
    v = Plasmid(VECTOR_NAME,
                p.seq[:i] + SCAFFOLD_AU_FLIP + TERMINATOR + ILLUMINA_ADAPTER
                + "N" * RSL_LEN + p.seq[i + len(SCAFFOLD_V1) + len(TERMINATOR):], True, [])
    return clone_guide(v, spacer) if cloned else v


LIBRARY = dict(genes=2325, guides=23279, non_targeting=101,
               source="guide sequences from Wang et al. 2014", synthesis="CustomArray")

# What the RSL buys, in the paper's own terms.
ANALYSES = (
    ("TCA", "total count analysis", "conventional read counting, ignoring the RSL"),
    ("LDA", "lineage dropout analysis", "count distinct RSLs per guide: how many clones "
                                        "survived, not how many reads"),
    ("IRA", "internal replicate analysis", "treat RSL sublineages within one guide as "
                                           "replicates generated inside a single screen"),
)


# ------------------------------------------------- a padlock probe for this vector
# Designed with tools/design_padlock.py against the rebuilt vector. The extension arm is
# unchanged (U6 is untouched by the AU-flip); only the ligation arm moves.
#
# The trap: the CRISPR-MIP backbone carries AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC, and this
# vector carries GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -- a substring of it. A capture that
# spans the vector's copy therefore puts TWO Read-2 primer sites in one amplicon unless the
# probe backbone drops its own.
# WARNING: the RSL sits 141 nt from the end of the extension arm, and 115 nt of that is
# fixed construct (scaffold + terminator + the built-in Illumina adapter). The extension arm
# must be upstream of the guide, so 141 nt is the FLOOR for capturing guide and RSL together.
# The published CRISPR-MIP probe fills 112 nt. This design is therefore 1.26x the only gap
# length the chemistry is known to work at, and has not been tested. Treat as a candidate.
PADLOCK_LIG_ARM = DOWNSTREAM[:26]                 # 26 nt, Tm 58.5, 4.4 C above the ext arm
PADLOCK_GAP = 142
PADLOCK_CIRCLE = 271
RSL_DISTANCE = 141                                # irreducible guide-to-RSL span
VALIDATED_GAP = 112                               # what the published probe actually fills

# An alternative arm inside the scaffold's 3' half, clear of both AU-flip positions. It
# captures the guide only (gap 66), but works unchanged on lentiCRISPR v1, lentiCRISPRv2,
# lentiGuide-Puro AND this vector -- the same trick that makes the Joung readout primer
# vector-independent. Needs its own PCR primers: the capture is too short to hold the
# existing pair.
PADLOCK_LIG_ARM_UNIVERSAL = SCAFFOLD_V1[45:68]    # scaffold positions 46-68 (1-based)
PADLOCK_GAP_UNIVERSAL = 66


# --------------------------------------------- the drafted probes (see 02_padlock_protocol.md)
# Backbone for option A drops the Read-2 site, because the vector already carries it.
def probe_A_backbone(index: str = cm.PROBE_INDICES[0]) -> str:
    return cm.UMI + index + cm.READ1_SITE


PROBE_LEN_A = len(cm.EXT_ARM) + len(probe_A_backbone()) + len(PADLOCK_LIG_ARM)
CIRCLE_LEN_A = PADLOCK_GAP + PROBE_LEN_A
PCR_FWD_A = il.P5 + ILLUMINA_ADAPTER[-20:]
PCR_REV_A = cm.P7_TRACR_REV                       # the CRISPR-MIP reverse primer, unchanged
CUSTOM_SEQ_PRIMER_A = ILLUMINA_ADAPTER[-20:]      # = 3' 20 nt of the TruSeq Read-2 primer
# The amplicon is the circle minus the arc outside the primer pair (45 nt, measured in
# selftest), plus the two non-templated P5/P7 tails. Derived rather than typed, because it
# moves with the CRISPR-MIP backbone's own lengths.
OUTSIDE_PRIMERS_A = 45
LIBRARY_LEN_A = CIRCLE_LEN_A - OUTSIDE_PRIMERS_A + len(il.P5) + len(il.P7)


# --------------------------------------------------------- the PUBLISHED readout (Methods) 🟢
# Schmierer et al. 2017, Mol Syst Biol 13:945, "Library preparation and sequencing".
# Three nested PCRs off 200 ug gDNA: PCR1 in 40 parallel reactions (5 ug each) x 14 cycles,
# pooled; PCR2 19 cycles off 5 ul of that pool; PCR3 14 cycles off 2 ul of PCR2.
# KAPA HiFi HotStart throughout. The 288-bp product is gel purified.
# PCR1 anneals outside the cassette, in plain vector sequence.
# TODO(lib): neither of these two is an Illumina or CRISPR part, so neither belongs in lib/.
PCR1_FW = "GGACTATCATATGCTTACCGTAACTTGAAAGTATTTCG"
PCR1_REV = "CTTTAGTTTGTATGTCTGTTGCTATTATGTCTACTATTCTTTCC"

# PCR2 and PCR3 build the Illumina ends out of lib constants. Every tail below is a slice of
# il.*; the only typed sequence is the vector anneal, which is a slice of DOWNSTREAM.
READ1_SITE = il.TRUSEQ_READ1[4:]                  # 29 nt: TruSeq Read 1 less its 5' ACAC,
#                                                 which P5's own 3' ACAC supplies
VECTOR_ANNEAL_REV = revcomp(DOWNSTREAM[-21:] + DOWNSTREAM_PCR)
PCR2_FW = READ1_SITE + U6_OVERLAP[-23:-2]         # ...stops 2 nt short of the +1 G
PCR2_REV = il.P7[5:] + VECTOR_ANNEAL_REV          # P7 less its 5' CAAGC
PCR3_REV = il.P7 + VECTOR_ANNEAL_REV[:10]
LIBRARY_LEN_PUBLISHED = 288                       # the paper's own figure, reproduced in selftest

# Martin's CRISPR_PCR1-F / -R, as transcribed from the workbook in to_debug/. Defined here
# once, and used by both the page and the selftest: -F is Schmierer's PCR1-F with a 3-nt 5'
# extension, -R shares only its 5' 15 nt with Schmierer's PCR1-R and then diverges.
OUR_PCR1_FW = "AAT" + PCR1_FW
OUR_PCR1_REV = PCR1_REV[-15:] + "CCTGCACTGT" + "TGTGGGCGATGTGCGCTCTG"

# PCR3's forward primer carries the i5 sample index between the P5 tail and the Read-1 site.
def pcr3_fw(i5: str = "N" * 6) -> str:
    """P5 + i5 sample index + the Read-1 site. The paper writes the index as a blank.

    Note what it omits: the canonical dual-index P5 arm is P5 + i5 + TruSeq Read 1 *whole*,
    i.e. with the ACAC. Here the index is followed by READ1_SITE, so the library has no ACAC
    in front of the Read-1 site -- which is why the reverse-complement-workflow i5 primer
    (il.INDEX2_PRIMER_RC, 3' end GTGT) has no site. See SEQ_PRIMERS.
    """
    return il.P5 + i5 + READ1_SITE[:-4]           # ...less the ATCT: the vector supplies it


# The custom read primer -- the paper spells the name "CRIPSRSEQ" in its own primer table.
# It ends ...GACGAAACACC + the Pol III +1 G, so its 3' end abuts the spacer exactly: the
# guide is read from cycle 1 with no leader to skip.
CUSTOM_SEQ_PRIMER = READ1_SITE[-6:] + U6_OVERLAP[-23:]

# Sequencing on a HiSeq4000: 20 cycles Read 1 (the guide), then TWO index reads.
INSTRUMENT = "HiSeq 4000"
READ1_CYCLES = 20
I5_CYCLES = 6                                     # i5 = the Illumina sample index
I7_CYCLES = 6                                     # i7 = the RSL  <-- note which is which


# ----------------------------------------------------- the finished library, as sequenced
def _seg(name, top, tag=None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


def final_library(spacer: str = "N" * SPACER_LEN, i5: str = "I" * I5_CYCLES,
                  rsl: str = "N" * RSL_LEN) -> Construct:
    """The gel-purified PCR3 product, top strand -- the thing that goes on the flow cell.

    Every segment is a reference: nothing here is typed twice. selftest checks this
    construct base-for-base against the PCR3 product amplified off the rebuilt vector, so
    it cannot drift from the simulation, and the sequencing-primer table is computed off it.
    """
    return Construct([
        _seg("Illumina P5", il.P5, "p5"),
        _seg("i5 sample index", i5, None, placeholder=True,
             note="carried by PCR3's forward primer"),
        # No ACAC in front of the Read-1 site: PCR3-F joins P5 straight to READ1_SITE.
        _seg("TruSeq Read 1 site", READ1_SITE, "p5"),
        _seg("U6 3' end", U6_OVERLAP[-23:], "r3"),
        _seg("sgRNA spacer", spacer, "cbc", placeholder=True),
        _seg("AU-flip scaffold", SCAFFOLD_AU_FLIP, "r2"),
        _seg("term", TERMINATOR, "me"),
        _seg("Illumina Read 2 / Index 1 site", ILLUMINA_ADAPTER, "t7",
             note="templated by the vector, not grafted on by a primer"),
        _seg("RSL", rsl, "umi", placeholder=True),
        _seg("vector", DOWNSTREAM[:-4], "w1"),
        _seg("PCR3-R site", revcomp(PCR3_REV[len(il.P7):]), "w1"),
        _seg("Illumina P7", il.P7_RC, "p7"),
    ], name="CRISPR-UMI (Schmierer) library")


# Sequencing primers, by reference: the sequences live in lib/illumina.py, except the
# paper's own read primer above. Everything else -- where each one lands, what it reads
# first -- is computed from final_library() by lib/seqprimers.py.
_PAPER = "Schmierer et al. 2017, Mol Syst Biol 13:945, primer table"

SEQ_PRIMERS = [
    sp.custom("Read 1", "CRIPSRSEQ (custom Read 1; the paper's own spelling)",
              CUSTOM_SEQ_PRIMER, _PAPER,
              f"Spiked in and used instead of the stock primer: its 3' base is the Pol III "
              f"+1 G, so cycle 1 is spacer base 1 and all {READ1_CYCLES} cycles are guide."),
    # Not a dud, but useless here: an important distinction, so the reason spells it out.
    sp.custom("Read 1", sp.TRUSEQ["R1"].name + " (for contrast: NOT what is loaded)",
              sp.TRUSEQ["R1"].seq, sp.TRUSEQ["R1"].source,
              "It does prime (its 5' ACAC is a flap: PCR3-F omits it), but read 1 would "
              f"start in U6, {len(U6_OVERLAP[-23:])} nt before the spacer, so all "
              f"{READ1_CYCLES} cycles would be constant vector. That is why the paper "
              "spikes in CRIPSRSEQ."),
    sp.TRUSEQ["I1"],
    sp.custom("Index 2 (i5)", f"flow-cell P5 oligo ({INSTRUMENT}, forward-strand workflow)",
              il.P5, 'Illumina "Indexed Sequencing Overview Guide" #15057455',
              "A forward-strand instrument has no free i5 primer: it extends the grafted "
              "P5 oligo, so the i5 read begins on the index itself."),
    sp.mismatching(sp.TRUSEQ["I2"],
                   "the reverse-complement-workflow primer ends ...GAGTGT, the complement "
                   "of the ACAC that PCR3's forward primer omits -- so its 3' end is "
                   "unpaired and this library cannot be i5-read on a NovaSeq-style "
                   "instrument without re-designing PCR3-F"),
    sp.mismatching(sp.TRUSEQ["R2"],
                   "no Read 2 is run (the RSL is read as the i7), and the stock primer "
                   "could not run anyway: the templated site is the canonical 33-nt "
                   "Index-1 oligo, one A short of the 34-nt revcomp of Read 2, and that "
                   "missing A is exactly where the primer's 3' base would pair"),
]
