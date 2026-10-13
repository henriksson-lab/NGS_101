from dsb_capture import BLESS_DISTAL, BLESS_PROXIMAL, get
from topology import HairpinTopology
globals().update(vars(get("bless")))

STEM_NT = 42
PROXIMAL_TOPOLOGY = HairpinTopology(
    "BLESS L1 proximal hairpin", BLESS_PROXIMAL, 0,
    len(BLESS_PROXIMAL) - STEM_NT, STEM_NT,
    colour="var(--c-r1)")
DISTAL_TOPOLOGY = HairpinTopology(
    "BLESS L3 distal hairpin", BLESS_DISTAL, 0,
    len(BLESS_DISTAL) - STEM_NT, STEM_NT,
    colour="var(--c-r2)")
