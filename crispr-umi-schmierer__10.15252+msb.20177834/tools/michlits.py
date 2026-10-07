"""
CRISPR-UMI (Michlits) -- the second of the two library-encoded UMI methods.

Michlits G, et al. "CRISPR-UMI: single-cell lineage tracing of pooled CRISPR-Cas9 screens."
Nat Methods 2017;14(12):1191-1197. doi:10.1038/nmeth.4466

Distinct from the Schmierer method in `crisprumi.py` (same directory): the barcode is
cloned in a SEPARATE STEP BEFORE the guide, so each ligation pairs a guide with a different
barcode and the UMI is the sgRNA-barcode PAIR. And the readout fights PCR bias by enriching
the cassette 10^3-10^4 fold with a PacI digest and size selection, rather than by replacing
the PCR.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import illumina as il
# Same directory: the two papers share the U6 3' end verbatim, so it is defined once, there.
from crisprumi import U6_OVERLAP

VECTORS = {"pLenti-UMI": dict(addgene=222694, bp=10694, enzyme="BsmBI"),
           "pRetro-UMI": dict(addgene=222686, bp=6991, enzyme="BbsI")}

BARCODE_LEN = 10              # the random stretch cloned in step 1
EXP_INDEX_LEN = 6             # the experimental (sample) index, added by the PCR primer
ILLUMINA_I7_SITE = il.INDEX1_PRIMER      # built into both vectors, as in Schmierer

# P7 as carried in the vector: one C->A so the adaptor contains no BbsI site, which would
# otherwise be cut during the step-2 Golden Gate. The readout primer restores the C.
P7_RC_CANONICAL = il.P7_RC
BBSI_FIX_POS = 17                        # 0-based, in P7_RC: the C that becomes an A
P7_RC_IN_VECTOR = (P7_RC_CANONICAL[:BBSI_FIX_POS] + "A"
                   + P7_RC_CANONICAL[BBSI_FIX_POS + 1:])

PACI = "TTAATTAA"
PACI_FRAGMENT_BP = 589        # the paper's own construct; the 2024 deposit gives 596
ENRICHMENT_FOLD = (1000, 10000)
CASSETTE_FRACTION_PPM = 0.1   # of total gDNA, before enrichment

# Readout. A single short PCR on the enriched template; its product IS the library.
# TODO(lib): the two vector anneals below ("CGAGGGCCT..." and "ACCGTTGATGAGTAG") are
# pLenti-UMI sequence, not Illumina parts, so they stay here.
PCR_FWD = il.P5 + "N" * EXP_INDEX_LEN + "CGAGGGCCTATTTCCCATGATTCCTTC"
PCR_REV = il.P7 + "ACCGTTGATGAGTAG"
# U6 again, plus 9 nt of vector 5' of it; ends on the Pol III +1 G like Schmierer's.
CUSTOM_READ1_PRIMER = "CGATTTCTT" + U6_OVERLAP

# Guides carrying any of these are excluded, because the two-step Golden Gate would cut them.
FORBIDDEN_IN_GUIDE = ("GAAGAC", "GTCTCC", "CTCGAG", "CGTCTC", "GAGACG")
FORBIDDEN_PREFIX, FORBIDDEN_SUFFIX = "AAGAC", "CTCGA"

LIBRARY = dict(genes=6560, guides_per_gene=4, non_targeting=112,
               cloning_events_per_sgrna=(954, 8776), complexity=83_500_000, skew_fold=4)


def guide_allowed(spacer: str) -> bool:
    """Does this spacer survive the paper's exclusion rules?"""
    s = spacer.upper()
    return (not any(m in s for m in FORBIDDEN_IN_GUIDE)
            and not s.startswith(FORBIDDEN_PREFIX) and not s.endswith(FORBIDDEN_SUFFIX))
