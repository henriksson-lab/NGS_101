"""
Plasmid handling: GenBank parsing, restriction digestion, primer binding and amplicons.

Built so that claims about a protocol ("this primer pair gives a 263 bp product",
"BsmBI leaves an ACCG overhang") can be *checked against a real map* rather than copied
out of a paper. Everything here is circular-aware, because lentiviral vectors are circular
and amplicons across the cloning site routinely wrap the origin.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

from chemdraw import revcomp


# ------------------------------------------------------------------ GenBank
@dataclass
class Feature:
    key: str
    start: int          # 0-based inclusive
    end: int            # 0-based exclusive
    strand: int         # +1 / -1
    qualifiers: dict = field(default_factory=dict)

    @property
    def label(self) -> str:
        for k in ("label", "gene", "product", "note"):
            if k in self.qualifiers:
                return self.qualifiers[k]
        return self.key

    def __len__(self) -> int:
        return self.end - self.start


@dataclass
class Plasmid:
    name: str
    seq: str
    circular: bool = True
    features: list[Feature] = field(default_factory=list)

    def __len__(self) -> int:
        return len(self.seq)

    def sub(self, start: int, end: int) -> str:
        """Sequence from start to end, wrapping the origin if the plasmid is circular."""
        n = len(self.seq)
        if 0 <= start <= end <= n:
            return self.seq[start:end]
        if not self.circular:
            raise ValueError(f"{start}..{end} runs off a linear sequence of {n}")
        return (self.seq * 3)[start + n:end + n]

    def feature(self, pattern: str) -> list[Feature]:
        rx = re.compile(pattern, re.I)
        return [f for f in self.features if rx.search(f.label) or rx.search(f.key)]


_LOC = re.compile(r"(complement\()?(?:join\()?(\d+)\.\.(\d+)")


def read_genbank(path: str) -> Plasmid:
    """Minimal GenBank reader: LOCUS, FEATURES and ORIGIN. Coordinates are 1-based
    inclusive in the file and converted to 0-based half-open here."""
    text = open(path).read()
    m = re.search(r"^LOCUS\s+(\S+)\s+(\d+) bp(.*)$", text, re.M)
    name = m.group(1) if m else path
    circular = bool(m and "circular" in m.group(3).lower())

    seq = ""
    if "ORIGIN" in text:
        seq = re.sub(r"[^acgtnACGTN]", "", text.split("ORIGIN")[-1].split("//")[0]).upper()
    if m and int(m.group(2)) != len(seq):
        raise ValueError(f"{path}: LOCUS says {m.group(2)} bp but ORIGIN holds {len(seq)}")

    features: list[Feature] = []
    if "FEATURES" in text:
        block = text.split("FEATURES", 1)[1].split("ORIGIN")[0]
        cur: Feature | None = None
        qual_key = None
        for line in block.split("\n"):
            if re.match(r"^ {5}\S", line):                      # new feature
                key, loc = line[5:21].strip(), line[21:].strip()
                lm = _LOC.search(loc)
                if lm:
                    cur = Feature(key, int(lm.group(2)) - 1, int(lm.group(3)),
                                  -1 if lm.group(1) else 1)
                    features.append(cur)
                qual_key = None
            elif cur is not None and re.match(r"^ {21}/", line):  # qualifier
                q = line.strip()[1:]
                if "=" in q:
                    qual_key, v = q.split("=", 1)
                    cur.qualifiers[qual_key] = v.strip('"')
                else:
                    qual_key = None
            elif cur is not None and qual_key and re.match(r"^ {21}\S", line):
                cur.qualifiers[qual_key] += line.strip().strip('"')
    return Plasmid(name=name, seq=seq, circular=circular, features=features)


# ------------------------------------------------------------------ searching
def find_all(seq: str, motif: str, circular: bool = True) -> list[int]:
    """0-based start positions of `motif`, wrapping the origin when circular."""
    hay = seq + (seq[:len(motif) - 1] if circular and len(motif) > 1 else "")
    out, i = [], hay.find(motif)
    while i != -1:
        out.append(i % len(seq))
        i = hay.find(motif, i + 1)
    return sorted(set(out))


def find_both(seq: str, motif: str, circular: bool = True) -> list[tuple[int, int]]:
    """(position, strand) for a motif on either strand. Position is always plus-strand."""
    hits = [(p, 1) for p in find_all(seq, motif, circular)]
    rc = revcomp(motif)
    if rc != motif:
        hits += [(p, -1) for p in find_all(seq, rc, circular)]
    return sorted(hits)


# -------------------------------------------------------------- restriction
@dataclass(frozen=True)
class Enzyme:
    """A Type IIS (or other) restriction enzyme.

    `top_cut` / `bottom_cut` are offsets from the END of the recognition site on the
    top strand. BsmBI is CGTCTC(1/5): top strand cut 1 nt downstream, bottom 5 nt,
    leaving a 4-nt 5' overhang.
    """
    name: str
    site: str
    top_cut: int
    bottom_cut: int

    @property
    def overhang_len(self) -> int:
        return abs(self.bottom_cut - self.top_cut)


BsmBI = Enzyme("BsmBI/Esp3I", "CGTCTC", 1, 5)
BbsI = Enzyme("BbsI", "GAAGAC", 2, 6)
BsaI = Enzyme("BsaI", "GGTCTC", 1, 5)


@dataclass
class Cut:
    enzyme: str
    site_start: int        # 0-based start of the recognition site, plus-strand coords
    strand: int
    top_cut: int           # 0-based position of the top-strand nick
    bottom_cut: int
    overhang: str          # the 5' overhang sequence left behind, 5'->3'


def digest(p: Plasmid, enz: Enzyme = BsmBI) -> list[Cut]:
    """Every cut this enzyme makes, with the 5' overhang each one leaves."""
    cuts: list[Cut] = []
    n = len(p.seq)
    for pos, strand in find_both(p.seq, enz.site, p.circular):
        if strand == 1:
            end = pos + len(enz.site)                 # first base after the site
            top, bot = end + enz.top_cut, end + enz.bottom_cut
            overhang = p.sub(top % n, bot % n)
        else:
            start = pos                               # site reads leftward on the minus strand
            top, bot = start - enz.top_cut, start - enz.bottom_cut
            overhang = p.sub(bot % n, top % n)
            top, bot = bot, top
        cuts.append(Cut(enz.name, pos, strand, top % n, bot % n, overhang))
    return sorted(cuts, key=lambda c: c.top_cut)


# -------------------------------------------------------------------- PCR
@dataclass
class Binding:
    primer: str
    strand: int            # +1 binds the plus strand (extends rightward)
    three_prime: int       # 0-based position of the primer's 3'-most templated base
    annealed: int          # how many 3' bases matched
    tail: str              # unmatched 5' tail (adapters, indices)


def bind(p: Plasmid, primer: str, min_anneal: int = 15) -> list[Binding]:
    """Where does this primer anneal? Matching is anchored at the 3' end, so 5' tails
    (P5/P7, indices, stagger) are allowed and reported."""
    primer = primer.upper().replace(" ", "")
    out: list[Binding] = []
    for strand in (1, -1):
        probe_full = primer if strand == 1 else revcomp(primer)
        best = 0
        for k in range(len(primer), min_anneal - 1, -1):
            probe = probe_full[-k:] if strand == 1 else probe_full[:k]
            for pos in find_all(p.seq, probe, p.circular):
                if k < best:
                    continue
                best = k
                three = (pos + k - 1) % len(p) if strand == 1 else pos
                tail = primer[:len(primer) - k]
                out.append(Binding(primer, strand, three, k, tail))
            if best:
                break
    return out


@dataclass
class Amplicon:
    length: int
    seq: str
    fwd: Binding
    rev: Binding
    wraps_origin: bool


def amplify(p: Plasmid, fwd: str, rev: str, min_anneal: int = 15) -> list[Amplicon]:
    """Predicted PCR products, including the primers' 5' tails. Circular-aware."""
    fs = [b for b in bind(p, fwd, min_anneal) if b.strand == 1]
    rs = [b for b in bind(p, rev, min_anneal) if b.strand == -1]
    out: list[Amplicon] = []
    n = len(p)
    for f in fs:
        for r in rs:
            start = f.three_prime - f.annealed + 1
            end = r.three_prime + r.annealed
            if p.circular:
                # only wrap coordinates on a circle -- doing it on a linear template
                # sends an amplicon that ends on the last base to 0 and loses it
                start %= n
                end %= n
                wraps = end <= start
            else:
                if end <= start:
                    continue
                wraps = False
            inner = p.sub(start, end + n if wraps else end)
            seq = f.tail + inner + revcomp(r.tail)
            out.append(Amplicon(len(seq), seq, f, r, wraps))
    return sorted(out, key=lambda a: a.length)
