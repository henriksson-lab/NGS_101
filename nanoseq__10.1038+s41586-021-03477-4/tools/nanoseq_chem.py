"""NanoSeq-specific restrictive end-repair view; source audit remains separate."""
from chemdraw import Scene, Segment


def _p(name: str, base: str, n: int, tag: str | None = None) -> Segment:
    return Segment(name, base * n, tag, placeholder=True)


def nanoseq_restrictive_end_repair() -> Scene:
    """NanoSeq end state: intact restriction ends, with nick extension suppressed."""
    top = [_p("restriction fragment", "X", 34), Segment("blocked nick", "H", placeholder=True)]
    sc = Scene(); sc.strand("top", top, label="restriction-generated genomic duplex")
    sc.anneal("bottom", [_p("restriction fragment partner", "x", 34)], to="top",
              pair=("restriction fragment partner", "restriction fragment"),
              label="restriction-generated genomic duplex")
    sc.mark("top", "blocked nick", "ddNTP terminates nick extension")
    return sc
