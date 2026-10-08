#!/usr/bin/env python3
"""
Harvest the protocol catalogue from Teichlab/scg_lib_structs into one table.

scg_lib_structs is the reference collection this repo's house style is modelled on. Its
README lists every method it documents, grouped by assay type, plus a TODO list of methods
it has not drawn yet -- and those TODO entries already carry a paper link each. Each
documented method page cites its own paper(s) as links.

So the catalogue is scraped, never typed:

    README            -> category, display name, page URL (or, for TODO rows, paper URL)
    each method page  -> the paper link(s) it cites, with the citation text
    Crossref / NCBI   -> title, journal, year, PMID for each DOI

Output: catalogue/scg_lib_structs.tsv, one row per (protocol, paper) pair -- a protocol
documented in two papers gets two rows, which is why `protocol` is not a key. The `slug`
column is the proposed directory name, `<name>_pmid<PMID>` (or `_doi<...>` when a paper has
no PMID, e.g. a preprint).

Downloads are cached, so re-running is cheap and offline-safe:

    python3 catalogue/tools/fetch_scg_lib_structs.py              # scrape + resolve
    python3 catalogue/tools/fetch_scg_lib_structs.py --cache DIR   # cache elsewhere
    python3 catalogue/tools/fetch_scg_lib_structs.py --offline     # cache only, no network
    python3 catalogue/tools/fetch_scg_lib_structs.py --refresh     # re-resolve every paper

Papers already resolved in the existing table are reused, not looked up again, so a
re-run only costs requests for pages and papers that are new.

The cache holds third-party HTML. It is only ever read as text and parsed with regexes
here; nothing from it is executed, and it is kept outside the repo.
"""

from __future__ import annotations

import argparse
import csv
import html
import json
import re
import subprocess
import sys
import time
import urllib.parse
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "catalogue" / "scg_lib_structs.tsv"
OURS = ROOT / "catalogue" / "ours.tsv"
README_URL = "https://raw.githubusercontent.com/Teichlab/scg_lib_structs/master/README.md"
PAGE_BASE = "https://teichlab.github.io/scg_lib_structs/methods_html/"
REPO = "https://github.com/Teichlab/scg_lib_structs"
UA = "chem-catalogue/1.0 (+https://github.com/henriksson-lab; mailto:he.johan@gmail.com)"

# A link on a method page is a paper if it points at one of these.
PAPER_HOSTS = re.compile(
    r"(doi\.org|pubmed\.ncbi|ncbi\.nlm\.nih\.gov/pmc|biorxiv\.org|medrxiv\.org|"
    r"researchsquare|nature\.com|science\.org|sciencemag\.org|sciencedirect\.com|"
    r"cell\.com|genome\.cshlp|genesdev\.cshlp|genomebiology|biomedcentral|"
    r"academic\.oup\.com|pnas\.org|embopress|elifesciences|frontiersin|"
    r"journals\.plos|wiley|springer|tandfonline|mdpi)", re.I)
# ...but not these: navigation, data files, tools.
NOT_PAPER = re.compile(r"scg_lib_structs|readthedocs|seqspec|/data/|\.pdf$|\.fa$|\.txt$|"
                       r"support\.illumina|illumina\.com|10xgenomics\.com|teichlab\.github", re.I)

DOI_IN_URL = re.compile(r"(10\.\d{4,9}/[^\s\"'<>&?#]+)")


def fetch(url: str, cache: Path, offline: bool = False, pause: float = 0.34) -> str:
    """GET `url`, memoised on disk by a sanitised filename.

    Uses curl rather than urllib: the system Python here ships without a CA bundle, and
    curl uses the OS trust store, so certificates are still verified.
    """
    key = re.sub(r"[^A-Za-z0-9._-]", "_", url)[-180:]
    hit = cache / key
    if hit.exists():
        return hit.read_text(encoding="utf-8", errors="replace")
    if offline:
        raise RuntimeError(f"--offline and not cached: {url}")
    for wait in (pause, 2, 5, 15):                      # be polite; back off on HTTP 429
        time.sleep(wait)
        done = subprocess.run(["curl", "-sSL", "--fail", "--max-time", "60",
                               "-A", UA, url], capture_output=True, text=True, check=False)
        if "429" not in done.stderr:
            break
    if done.returncode != 0:
        raise RuntimeError(f"curl {done.returncode} for {url}: {done.stderr.strip()[:120]}")
    hit.write_text(done.stdout, encoding="utf-8")
    return done.stdout


# --------------------------------------------------------------------- the README

def parse_readme(md: str) -> list[dict]:
    """-> [{category, name, page_url | paper_url, documented}] in README order."""
    rows, category = [], ""
    for line in md.split("\n"):
        h = re.match(r"-\s+###\s+(.*)", line.strip())
        if h:
            category = h.group(1).strip().rstrip(":")
            continue
        if not re.match(r"\s*-\s+\[", line):
            continue
        links = re.findall(r"\[([^\]]+)\]\(([^)]+)\)", line)
        if not links:
            continue
        # A TODO row may carry several named papers: "[Slide-seq](url) / [Slide-seqV2](url)".
        if links[0][1].startswith(PAGE_BASE):
            name, url = links[0]
            rows.append(dict(category=category, name=name.strip(), page_url=url,
                             paper_urls=[], documented=True))
        else:
            for name, url in links:
                rows.append(dict(category=category, name=name.strip(), page_url="",
                                 paper_urls=[url], documented=False))
    return rows


def family_members(name: str) -> list[str]:
    """"SMART-seq family (including SMART-seq, SMART-seq2/3 ...)" -> the listed members."""
    m = re.search(r"\(including ([^)]+)\)", name)
    if not m:
        return []
    return [p.strip() for p in re.split(r",| and ", m.group(1)) if p.strip()]


def short_name(name: str) -> str:
    return re.sub(r"\s*\(including[^)]*\)", "", name).strip()


# ----------------------------------------------------------------- a method page

# Every method page is: a preamble that introduces the method and cites ITS OWN paper(s),
# then the oligo list and the numbered steps. Papers linked from inside a step are
# incidental technique citations -- the semi-suppressive-PCR paper, the method a step was
# adapted from -- and must not be recorded as the protocol's paper. Position decides it.
# Match only HEADINGS -- the preamble prose says "step-by-step" too, which is how the
# Quartz-seq page first fooled this into calling its own two papers incidental.
BODY_START = re.compile(r"<h[1-4][^>]*>[^<]*?(?:adapter and primer|step-by-step|\(1\))", re.I)
CITATION_TEXT = re.compile(r"\b(?:19|20)\d{2}\b")        # "Science 360, 176-182 (2018)"


def _norm(t: str) -> str:
    return re.sub(r"[^a-z0-9]", "", t.lower())


def page_papers(page_html: str, names: list[str]) -> list[tuple[str, str, str]]:
    """-> [(anchor text, url, role)] for the paper links a method page carries, in order.

    "primary" = the method's own paper, by any of three signals: the link sits in the
    preamble, its anchor text names the method (or a family member), or the anchor reads
    like a citation. "cited" = an incidental technique reference from inside a step, e.g.
    the semi-suppressive-PCR paper, which 20 pages link and none of them is about.
    """
    b = BODY_START.search(page_html)
    boundary = b.start() if b else len(page_html)
    h1 = re.search(r"<h1[^>]*>(.*?)</h1>", page_html, re.S | re.I)
    blob = _norm(" ".join(names) + " " + re.sub(r"<[^>]+>", " ", h1.group(1) if h1 else ""))
    out, seen = [], set()
    for m in re.finditer(r'<a[^>]+href="([^"]+)"[^>]*>(.*?)</a>', page_html, re.S | re.I):
        url = html.unescape(m.group(1)).strip()
        text = html.unescape(re.sub(r"<[^>]+>", " ", m.group(2)))
        text = re.sub(r"\s+", " ", text).strip()
        if NOT_PAPER.search(url) or not PAPER_HOSTS.search(url):
            continue
        if url in seen:
            continue
        seen.add(url)
        n = _norm(text)
        named = len(n) >= 3 and n in blob
        primary = m.start() < boundary or named or bool(CITATION_TEXT.search(text))
        out.append((text, url, "primary" if primary else "cited"))
    return out


# ------------------------------------------------------------------- identifiers

# Publishers whose article URL encodes the DOI by a fixed rule, so no lookup is needed.
URL_DOI_RULES = (
    # Modern Nature ids map straight to a DOI; legacy ones (nbt0400_424, cr201782) do not,
    # which the repair pass below catches when the DOI returns no metadata.
    (re.compile(r"nature\.com/articles/([A-Za-z0-9._-]+?)(?:\.pdf)?/?$"), "10.1038/{0}"),
    (re.compile(r"science(?:mag)?\.org/content/early/[\d/]+/(science\.[a-z0-9]+)"),
     "10.1126/{0}"),
    (re.compile(r"elifesciences\.org/articles/(\d+)"), "10.7554/eLife.{0}"),
    (re.compile(r"(?:www\.)?researchsquare\.com/article/(rs-[\w.-]+)/v(\d+)"),
     "10.21203/rs.3.{0}/v{1}"),
    (re.compile(r"protocolexchange\.researchsquare\.com/article/(pex-[\w.-]+)/v(\d+)"),
     "10.21203/rs.3.{0}/v{1}"),
    (re.compile(r"genomebiology\.biomedcentral\.com/articles/(10\.\d+/[\w.-]+)"), "{0}"),
)


def doi_of(url: str) -> str:
    """The DOI a paper URL carries or implies, without a network lookup."""
    if "biorxiv.org/content/" in url or "medrxiv.org/content/" in url:
        m = re.search(r"content/(10\.\d{4,9}/[^\s\"'<>?#]+?)(v\d+)?(\.full)?/?$", url)
        if m:
            return m.group(1)
    m = DOI_IN_URL.search(urllib.parse.unquote(url))
    if m:
        return re.sub(r"[.,;)]+$", "", m.group(1))
    for rx, tmpl in URL_DOI_RULES:
        m = rx.search(url)
        if m:
            return tmpl.format(*m.groups())
    return ""


# <meta name="citation_doi" content="..."> and friends: what is left after the rules are
# publisher landing pages (Science, ScienceDirect/Cell, CSHL), which all carry one.
META_DOI = re.compile(
    r"""<meta[^>]+(?:name|property)=["'](?:citation_doi|dc\.identifier|DC\.Identifier|"""
    r"""prism\.doi|og:url)["'][^>]+content=["']\s*(?:https?://(?:dx\.)?doi\.org/)?"""
    r"""(10\.\d{4,9}/[^"'\s]+)""", re.I)


# Paywalled publishers answer a landing-page GET with 403, so their DOI is recovered from
# metadata services instead: an Elsevier PII is a Crossref alternative-id, and a
# journal/volume/page citation is a PubMed query.
PII_URL = re.compile(r"(?:sciencedirect\.com/science/article/pii/|cell\.com/(?:[\w-]+/)*"
                     r"(?:abstract|fulltext)/)(S[\dX()\-.]+)", re.I)
VOL_PAGE_URL = re.compile(
    r"(science\.sciencemag\.org|science\.org|www\.pnas\.org|pnas\.org|genome\.cshlp\.org|"
    r"genesdev\.cshlp\.org)/content/(\d+)/[\w.-]+/(\d+)")
JOURNAL_OF_HOST = {"science.sciencemag.org": "Science", "science.org": "Science",
                   "www.pnas.org": "Proc Natl Acad Sci U S A",
                   "pnas.org": "Proc Natl Acad Sci U S A",
                   "genome.cshlp.org": "Genome Res", "genesdev.cshlp.org": "Genes Dev"}


def doi_from_pii(url: str, cache: Path, offline: bool) -> str:
    """Elsevier/Cell PII -> DOI, via Crossref's alternative-id index."""
    m = PII_URL.search(url)
    if not m:
        return ""
    pii = re.sub(r"[^A-Za-z0-9]", "", m.group(1)).upper()
    q = (f"https://api.crossref.org/works?filter=alternative-id:{pii}"
         "&select=DOI&rows=1&mailto=he.johan%40gmail.com")
    try:
        items = json.loads(fetch(q, cache, offline))["message"]["items"]
    except Exception:                                    # noqa: BLE001
        return ""
    return items[0]["DOI"] if items else ""


def pmid_from_citation(url: str, cache: Path, offline: bool) -> str:
    """A journal/volume/page URL -> PMID, via a PubMed field query."""
    m = VOL_PAGE_URL.search(url)
    if not m:
        return ""
    host, vol, page = m.groups()
    journal = JOURNAL_OF_HOST.get(host, "")
    if not journal:
        return ""
    term = urllib.parse.quote(f'"{journal}"[jour] AND {vol}[volume] AND {page}[page]')
    q = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&retmode=json"
         "&tool=chem-catalogue&email=he.johan%40gmail.com&term=" + term)
    try:
        ids = json.loads(fetch(q, cache, offline))["esearchresult"]["idlist"]
    except Exception:                                    # noqa: BLE001
        return ""
    return ids[0] if len(ids) == 1 else ""


def pmid_of_doi(doi: str, cache: Path, offline: bool) -> str:
    """DOI -> PMID by a PubMed field query. The PMC id converter only knows PMC-indexed
    records, so this catches everything else that is in PubMed at all."""
    term = urllib.parse.quote(f"{doi}[doi]")
    q = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&retmode=json"
         "&tool=chem-catalogue&email=he.johan%40gmail.com&term=" + term)
    try:
        ids = json.loads(fetch(q, cache, offline))["esearchresult"]["idlist"]
    except Exception:                                    # noqa: BLE001
        return ""
    return ids[0] if len(ids) == 1 else ""


def doi_of_pmids(pmids: list[str], cache: Path, offline: bool) -> dict[str, str]:
    """PMID -> DOI, from the articleids esummary returns."""
    out: dict[str, str] = {}
    uniq = sorted({str(p) for p in pmids if p})
    for i in range(0, len(uniq), 200):
        q = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed"
             "&retmode=json&tool=chem-catalogue&email=he.johan%40gmail.com&id="
             + ",".join(uniq[i:i + 200]))
        try:
            res = json.loads(fetch(q, cache, offline)).get("result", {})
        except Exception:                                # noqa: BLE001
            continue
        for pmid in uniq[i:i + 200]:
            for aid in (res.get(pmid) or {}).get("articleids", []):
                if aid.get("idtype") == "doi" and aid.get("value"):
                    out[pmid] = str(aid["value"]).strip()
    return out


def doi_from_page(url: str, cache: Path, offline: bool) -> str:
    """Fetch a paper's landing page and read the DOI out of its citation metadata."""
    try:
        page = fetch(url, cache, offline)
    except Exception as e:                               # noqa: BLE001 -- paywalls, 403s
        print(f"    no DOI for {url[:70]}: {e}", file=sys.stderr)
        return ""
    m = META_DOI.search(page)
    return re.sub(r"[.,;)]+$", "", m.group(1)) if m else ""


def pmid_of(url: str) -> str:
    m = re.search(r"pubmed\.ncbi\.nlm\.nih\.gov/(\d+)", url)
    return m.group(1) if m else ""


def resolve_dois(dois: list[str], cache: Path, offline: bool) -> dict[str, dict]:
    """DOI -> {pmid, pmcid} via NCBI's ID converter, 200 at a time."""
    out: dict[str, dict] = {}
    uniq = sorted({str(d) for d in dois if d})
    for i in range(0, len(uniq), 200):
        chunk = uniq[i:i + 200]
        url = ("https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?format=json"
               "&tool=chem-catalogue&email=he.johan%40gmail.com&ids="
               + urllib.parse.quote(",".join(chunk)))
        try:
            data = json.loads(fetch(url, cache, offline))
        except Exception as e:                           # noqa: BLE001
            print(f"  idconv failed ({e}); continuing without PMIDs for this chunk",
                  file=sys.stderr)
            continue
        for rec in data.get("records", []):
            if rec.get("doi"):
                # idconv returns some ids as JSON numbers
                out[str(rec["doi"]).lower()] = {"pmid": str(rec.get("pmid", "") or ""),
                                                "pmcid": str(rec.get("pmcid", "") or "")}
    return out


def crossref(doi: str, cache: Path, offline: bool) -> dict:
    """Title, journal, year for a DOI. Works for preprints too."""
    url = f"https://api.crossref.org/works/{urllib.parse.quote(doi)}"
    try:
        msg = json.loads(fetch(url, cache, offline))["message"]
    except Exception:                                    # noqa: BLE001
        return {}
    title = (msg.get("title") or [""])[0]
    date = (msg.get("issued", {}).get("date-parts") or [[""]])[0]
    return {"title": re.sub(r"\s+", " ", html.unescape(title)).strip(),
            "journal": (msg.get("container-title") or msg.get("institution", [{}])[0]
                        .get("name", "") if msg.get("institution") else
                        (msg.get("container-title") or [""])[0]),
            "year": str(date[0]) if date and date[0] else "",
            "type": msg.get("type", "")}


def esummary(pmids: list[str], cache: Path, offline: bool) -> dict[str, dict]:
    """PMID -> {title, journal, year} from NCBI, for papers Crossref did not answer for."""
    out: dict[str, dict] = {}
    uniq = sorted({str(p) for p in pmids if p})
    for i in range(0, len(uniq), 200):
        chunk = uniq[i:i + 200]
        url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed"
               "&retmode=json&tool=chem-catalogue&email=he.johan%40gmail.com&id="
               + ",".join(chunk))
        try:
            res = json.loads(fetch(url, cache, offline)).get("result", {})
        except Exception as e:                           # noqa: BLE001
            print(f"  esummary failed ({e})", file=sys.stderr)
            continue
        for pmid in chunk:
            rec = res.get(pmid)
            if not rec:
                continue
            out[pmid] = {"title": re.sub(r"\s+", " ", rec.get("title", "")).strip(". "),
                         "journal": rec.get("fulljournalname") or rec.get("source", ""),
                         "year": (rec.get("pubdate", "") or "")[:4]}
    return out


# ------------------------------------------------------- the directory naming scheme
# A protocol directory is named "<protocol>__<doi>", e.g.
#     smart-seq-family__10.1038+nmeth.2639
# Every "/" in the DOI becomes "+", which is the only substitution needed: "+" is
# filesystem- and shell-safe and does not occur in DOIs (the selftest checks that for
# every DOI in the table), so decoding is a plain reverse. Replacing only the first "/"
# would not do -- a Research Square DOI has two (10.21203/rs.3.rs-3210240/v1) -- and
# using "_" would be ambiguous, because DOI suffixes contain underscores themselves
# (10.1038/nbt0400_424). "__" separates the name from the DOI; a slugified name never
# contains a run of underscores.
#
# DOI rather than PMID because every preprint has one from the day it is posted, and
# because far more rows here carry a DOI than a PMID. A protocol with no paper at all
# (a vendor kit) is named by its name alone.

def encode_doi(doi: str) -> str:
    return doi.replace("/", "+")


def decode_doi(encoded: str) -> str:
    return encoded.replace("+", "/")


def name_slug(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "-", short_name(name)).strip("-").lower()


def slug(name: str, doi: str, pmid: str = "") -> str:
    """The proposed directory name for a protocol, given its defining paper."""
    base = name_slug(name)
    if doi:
        return f"{base}__{encode_doi(doi)}"
    if pmid:                              # no DOI anywhere: fall back to the PMID
        return f"{base}__pmid{pmid}"
    return base                           # vendor protocol, no publication


# Vendor kits whose page cites only the published techniques the kit reads out, none of
# them the kit's own paper. Like the other vendor pages they get a paperless defining row
# (named by name alone); the cited papers stay in the table as associated rows.
VENDOR_KITS = {
    "10x Chromium Single Cell 3' FeatureBarcoding",   # cites CRISPR-screen and CITE-seq
    "10x Chromium Single Cell ATAC",                  # cites Amini 2014 transposition
}


def own_names(protocol: str, members: str) -> tuple[str, list[str]]:
    """-> (head name, every name) for a protocol, normalised.

    "SPLiT-seq / microSPLiT" -> head "splitseq"; "SMART-seq family" -> head "smartseq",
    the first family member. The head is the method itself; the rest are later variants.
    """
    parts = [p for p in re.split(r"\s*/\s*|\s+and\s+", protocol) if p]
    fam = [m for m in members.split("; ") if m]
    head = fam[0] if fam else re.sub(r"\s+family$", "", parts[0])
    names = [head, *parts, *fam, *(p for f in fam for p in f.split("/"))]
    return _norm(head), [n for n in dict.fromkeys(map(_norm, names)) if len(n) >= 3]


def defining_paper(rows: list[dict]) -> dict[str, dict]:
    """protocol -> the row naming its directory: the method's OWN paper.

    A method page's preamble cites its own paper alongside the methods it is built from --
    ISSAAC-seq cites ATAC-seq, scNMT-seq cites Smart-seq2 and scBS-seq, SNARE-seq cites
    Drop-seq -- and those are always older, so "earliest primary" names the directory
    after a component. Rank instead, earliest first within a tier:

      0  the anchor text or title names the method itself (its head name)
      1  it names a later variant of it, reads like a citation, or names nothing known
      2  it names a DIFFERENT method in the catalogue: a component, kept only as fallback

    Every other paper stays in the table as a further row against the same protocol.
    """
    names = {r["protocol"]: own_names(r["protocol"], r["family_members"]) for r in rows}
    known = {n for _, ns in names.values() for n in ns}
    best: dict[str, dict] = {}
    for i, r in enumerate(rows):
        if r["role"] != "primary":
            continue
        key = r["protocol"]
        if key in VENDOR_KITS and r["paper_url"]:
            continue
        head, mine = names[key]
        text, title = _norm(r["citation_text"]), _norm(r["title"])
        if head in text or head in title:
            tier = 0
        elif any(n in text for n in mine) or not any(n in text for n in known):
            tier = 1
        else:
            tier = 2
        year = int(r["year"]) if r["year"].isdigit() else 9999
        rank = (tier, year, i)
        if key not in best or rank < best[key]["_rank"]:
            best[key] = {"_rank": rank, "row": r}
    out = {k: v["row"] for k, v in best.items()}
    # The same paper as a preprint first: name the directory for the earlier DOI.
    for k, win in out.items():
        for r in rows:
            if (r["protocol"] == k and r["role"] == "primary" and r is not win
                    and PREPRINT.search(r["paper_url"] + " " + r["journal"])
                    and r["year"] <= win["year"] and same_paper(r["title"], win["title"])):
                out[k] = r
                break
    return out


PREPRINT = re.compile(r"biorxiv|medrxiv|researchsquare|research square|arxiv", re.I)


def same_paper(a: str, b: str) -> bool:
    """Preprint and journal titles of one paper: mostly the same words."""
    wa, wb = set(re.findall(r"[a-z0-9]+", a.lower())), set(re.findall(r"[a-z0-9]+", b.lower()))
    return bool(wa and wb) and len(wa & wb) / min(len(wa), len(wb)) >= 0.6


COLUMNS = ["source", "category", "protocol", "family_members", "documented", "role",
           "scg_page", "paper_url", "citation_text", "doi", "pmid", "pmcid", "title",
           "journal", "year", "slug", "is_defining", "our_dir", "our_status", "note"]


RESOLVED = ("doi", "pmid", "pmcid", "title", "journal", "year")


# ours.tsv: one row per directory of ours. `section` says where the website lists it
# (published protocols, searchable; or wip, our own unpublished work) -- it is not
# copied into the scraped table, but a bad value stops the build here rather than later.
OURS_COLUMNS = ["dir", "protocol", "doi", "status", "section", "modality", "note"]
OURS_VOCAB = {"status": ("documented", "draft", "notes"),
              "section": ("published", "wip"),
              "modality": ("DNA", "RNA", "multi")}


def read_ours(path: Path = OURS) -> list[dict]:
    """catalogue/ours.tsv, with its header and closed vocabularies validated."""
    with path.open(encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        mine = list(reader)
        header = reader.fieldnames or []
    if header != OURS_COLUMNS:
        sys.exit(f"{path}: columns are {header}, expected {OURS_COLUMNS}")
    for m in mine:
        for col, allowed in OURS_VOCAB.items():
            if m[col] not in allowed:
                sys.exit(f"{path}: {m['dir']}: {col}={m[col]!r}, expected one of {allowed}")
    return mine


def known_papers(tsv: Path) -> dict[str, dict]:
    """paper URL -> its identifiers and metadata, from a previous run's table.

    Only fully resolved papers (a title came back) are reused; anything else is retried.
    A paper keeps its identifiers whichever protocol cites it, so the URL is the key.
    """
    if not tsv.exists():
        return {}
    with tsv.open(encoding="utf-8") as fh:
        return {r["paper_url"]: {k: r[k] for k in RESOLVED}
                for r in csv.DictReader(fh, delimiter="\t")
                if r["paper_url"] and r["title"] and r["doi"]}


def build(cache: Path, offline: bool, known: dict[str, dict] | None = None) -> list[dict]:
    readme = fetch(README_URL, cache, offline)
    entries = parse_readme(readme)
    print(f"README: {len(entries)} entries "
          f"({sum(e['documented'] for e in entries)} documented)", file=sys.stderr)

    for e in entries:
        if not e["documented"]:
            continue
        page = fetch(e["page_url"], cache, offline)
        e["papers"] = page_papers(page, [short_name(e["name"]), *family_members(e["name"])])
        n_pri = sum(1 for *_, role in e["papers"] if role == "primary")
        print(f"  {e['name'][:48]:48s} {n_pri} primary + "
              f"{len(e['papers']) - n_pri} cited", file=sys.stderr)

    rows: list[dict] = []
    for e in entries:
        papers = e.get("papers") or [("", u, "primary") for u in e["paper_urls"]]
        if not papers or short_name(e["name"]) in VENDOR_KITS:
            papers = [("", "", "primary"), *papers]
        for text, url, role in papers:
            rows.append({"source": "scg_lib_structs",
                         "category": e["category"], "protocol": short_name(e["name"]),
                         "family_members": "; ".join(family_members(e["name"])),
                         "documented": "yes" if e["documented"] else "todo",
                         "scg_page": e["page_url"], "paper_url": url,
                         "citation_text": text, "role": role,
                         "doi": doi_of(url), "pmid": pmid_of(url)})

    # Papers already resolved by a previous run are taken as they are, and every lookup
    # below sees only the new ones -- so a re-run costs a request per NEW paper.
    known = known or {}
    for r in rows:
        r.update(known.get(r["paper_url"], {}))
    done = sum(1 for r in rows if r["paper_url"] in known)
    print(f"{done} of {len(rows)} rows already resolved; looking up the rest",
          file=sys.stderr)
    old = rows
    rows = [r for r in old if r["paper_url"] not in known]

    # Resolve, cheapest first: PII -> Crossref, citation -> PubMed, else the landing page.
    missing = sorted({r["paper_url"] for r in rows if r["paper_url"] and not r["doi"]})
    if missing:
        print(f"resolving {len(missing)} paper URL(s) with no DOI in the URL...",
              file=sys.stderr)
        found: dict[str, str] = {}
        by_pmid: dict[str, str] = {}
        for u in missing:
            d = doi_from_pii(u, cache, offline)
            if not d:
                pm = pmid_from_citation(u, cache, offline)
                if pm:
                    by_pmid[u] = pm
                    continue
            if not d:
                d = doi_from_page(u, cache, offline)
            found[u] = d
        pm_dois = doi_of_pmids(list(by_pmid.values()), cache, offline)
        for r in rows:
            if r["doi"]:
                continue
            pm = by_pmid.get(r["paper_url"], "")
            if pm:
                r["pmid"] = r["pmid"] or pm
                r["doi"] = pm_dois.get(pm, "")
            else:
                r["doi"] = found.get(r["paper_url"], "")

    # A row that gave us only a PMID (a PubMed link) still needs its DOI.
    only_pmid = [r["pmid"] for r in rows if r["pmid"] and not r["doi"]]
    if only_pmid:
        pm_dois = doi_of_pmids(only_pmid, cache, offline)
        for r in rows:
            if r["pmid"] and not r["doi"]:
                r["doi"] = pm_dois.get(r["pmid"], "")

    ids = resolve_dois([r["doi"] for r in rows], cache, offline)
    for r in rows:
        got = ids.get(r["doi"].lower(), {})
        r["pmid"] = r["pmid"] or got.get("pmid", "")
        r["pmcid"] = got.get("pmcid", "")

    gaps = sorted({r["doi"] for r in rows if r["doi"] and not r["pmid"]})
    if gaps:
        print(f"looking up {len(gaps)} PMID(s) by DOI in PubMed...", file=sys.stderr)
        more = {d: pmid_of_doi(d, cache, offline) for d in gaps}
        for r in rows:
            if not r["pmid"]:
                r["pmid"] = more.get(r["doi"], "")

    meta_by_pmid = esummary([r["pmid"] for r in rows], cache, offline)
    for r in rows:
        meta = {}
        if r["doi"]:
            meta = crossref(r["doi"], cache, offline)
        if not meta.get("title") and r["pmid"]:
            meta = meta_by_pmid.get(r["pmid"], {})
        r["title"] = meta.get("title", "")
        r["journal"] = meta.get("journal", "")
        r["year"] = meta.get("year", "")

    # Repair pass. A DOI that no metadata service recognises is wrong, not just obscure --
    # a URL rule mis-firing on a legacy article id. Re-read it from the landing page.
    broken = sorted({r["paper_url"] for r in rows
                     if r["paper_url"] and r["doi"] and not r["title"]})
    if broken:
        print(f"repairing {len(broken)} DOI(s) that returned no metadata...", file=sys.stderr)
        for u in broken:
            real = doi_from_page(u, cache, offline)
            for r in rows:
                if r["paper_url"] != u:
                    continue
                if real and real.lower() != r["doi"].lower():
                    r["doi"] = real
                meta = crossref(r["doi"], cache, offline)
                r["title"] = meta.get("title", "")
                r["journal"] = meta.get("journal", "")
                r["year"] = meta.get("year", "")
        ids = resolve_dois([r["doi"] for r in rows if not r["pmid"]], cache, offline)
        for r in rows:
            if not r["pmid"]:
                r["pmid"] = ids.get(r["doi"].lower(), {}).get("pmid", "")

    for r in rows:
        if not r["title"] and r["pmid"]:
            m = esummary([r["pmid"]], cache, offline).get(r["pmid"], {})
            r.update({k: m.get(k, r[k]) for k in ("title", "journal", "year")})
    rows = old

    # One directory per protocol, named for its defining paper -- so every row of a
    # protocol carries the same slug, and `is_defining` marks which paper named it.
    define = defining_paper(rows)
    for r in rows:
        d = define.get(r["protocol"])
        src = d or r
        r["slug"] = slug(r["protocol"], src["doi"], src["pmid"])
        r["is_defining"] = "yes" if d is not None and r is d else "no"

    # ---- our own coverage, from catalogue/ours.tsv (the only record of it) -------------
    # Joined by directory name first: a directory IS the protocol's slug, so the join is
    # exact even for a vendor kit with no DOI, and for one paper defining two protocols
    # (HyDrop-RNA / HyDrop-ATAC), which a DOI join would map to a single directory. A
    # An exact protocol-name match comes next. A directory named otherwise falls back to
    # the DOI only for a single-paper upstream protocol: DOI matching every citation would
    # wrongly merge component methods (for example SMART-seq2) into every multi-omics
    # protocol that uses them, and would prevent us splitting an upstream "family" page.
    # A protocol we cover that upstream does not list is appended, so the table is the
    # whole worklist rather than only upstream's part of it.
    mine = read_ours()
    by_dir = {m["dir"]: m for m in mine}
    matched_protocols = {r["protocol"]: by_dir[r["slug"]] for r in rows if r["slug"] in by_dir}
    by_name = {m["protocol"].casefold(): m for m in mine}
    for r in rows:
        if r["protocol"] not in matched_protocols and r["protocol"].casefold() in by_name:
            matched_protocols[r["protocol"]] = by_name[r["protocol"].casefold()]
    # An upstream aggregate page can contain several protocols that we document
    # separately.  Its paper-anchor text names the individual method; when that name and
    # DOI both agree with ours, join that paper row to the individual directory without
    # assigning the whole aggregate to one member.
    matched_rows = {}
    for i, r in enumerate(rows):
        m = by_name.get(r["citation_text"].casefold())
        aggregate = (r["protocol"] + ";" + r["family_members"]).casefold()
        if m is None:
            candidates = [x for x in mine if x["doi"]
                          and x["doi"].casefold() == r["doi"].casefold()
                          and x["protocol"].casefold() in aggregate]
            m = candidates[0] if len(candidates) == 1 else None
        # Some upstream aggregate names omit the explicit version names (notably
        # ``inDrop`` for the v1/v2 papers).  A DOI-exact ours row whose slugged name is
        # the aggregate head plus a version suffix is still a row-level match.  Keep
        # this row-level: assigning the entire aggregate would collapse its versions.
        member_named = (m and (m["protocol"].casefold() in aggregate
                               or name_slug(m["protocol"]).startswith(
                                   name_slug(r["protocol"]) + "-")))
        if (m and r["protocol"] not in matched_protocols and member_named and m["doi"]
                and m["doi"].casefold() == r["doi"].casefold()):
            matched_rows[i] = m
    used = ({m["dir"] for m in matched_protocols.values()}
            | {m["dir"] for m in matched_rows.values()})
    doi_mine = {}
    for m in mine:
        if m["doi"] and m["dir"] not in used:
            doi_mine.setdefault(m["doi"], []).append(m)
    # A DOI shared by two separately documented protocols is not enough to select one.
    by_doi = {doi: ms[0] for doi, ms in doi_mine.items() if len(ms) == 1}
    primary_count = {p: sum(x["role"] == "primary" for x in rows if x["protocol"] == p)
                     for p in {x["protocol"] for x in rows}}
    for r in rows:
        if r["protocol"] in matched_protocols or primary_count[r["protocol"]] != 1 \
                or r["is_defining"] != "yes":
            continue
        head = re.split(r"\s*/\s*", r["protocol"], maxsplit=1)[0].casefold()
        m = by_name.get(head)
        if m and m["doi"].casefold() == r["doi"].casefold():
            matched_protocols[r["protocol"]] = m
    for r in rows:
        if (r["protocol"] not in matched_protocols and r["is_defining"] == "yes"
                and primary_count[r["protocol"]] == 1 and r["doi"] in by_doi):
            matched_protocols[r["protocol"]] = by_doi[r["doi"]]
    for i, r in enumerate(rows):
        m = matched_rows.get(i) or matched_protocols.get(r["protocol"])
        r["our_dir"] = m["dir"] if m else ""
        r["our_status"] = m["status"] if m else ""
        if m and not r.get("note"):
            r["note"] = m["note"]
    hit = ({m["dir"] for m in matched_protocols.values()}
           | {m["dir"] for m in matched_rows.values()})
    for m in mine:
        if m["dir"] in hit:
            continue
        url = f"https://doi.org/{m['doi']}" if m["doi"] else ""
        if url in known and known[url]["doi"] == m["doi"]:
            meta, pm = known[url], known[url]["pmid"]
        else:
            meta = crossref(m["doi"], cache, offline) if m["doi"] else {}
            pm = pmid_of_doi(m["doi"], cache, offline) if m["doi"] else ""
        rows.append({"source": "ours", "category": "Ours, not in scg_lib_structs",
                     "protocol": m["protocol"], "family_members": "", "documented": "ours",
                     "role": "primary", "scg_page": "", "paper_url": url,
                     "citation_text": "", "doi": m["doi"], "pmid": pm, "pmcid": "",
                     "title": meta.get("title", ""), "journal": meta.get("journal", ""),
                     "year": meta.get("year", ""), "slug": m["dir"], "is_defining": "yes",
                     "our_dir": m["dir"], "our_status": m["status"], "note": m["note"]})

    # Never assert an identifier nothing confirms: a DOI no metadata service knows is
    # probably a URL rule mis-firing, so drop it and say so. The URL is still recorded.
    for r in rows:
        if r["doi"] and not r["title"]:
            r["note"] = f"DOI unconfirmed by Crossref/NCBI, dropped: {r['doi']}"
            r["doi"] = ""
        elif r["documented"] == "yes" and not (define.get(r["protocol"], r)["doi"]
                                               or define.get(r["protocol"], r)["pmid"]):
            r["note"] = "vendor/kit protocol: no defining publication"
        elif r["protocol"] == "SMART-seq family" and not r["our_dir"]:
            r["note"] = "upstream aggregate; documented here as five separate protocols"
        else:
            r["note"] = r.get("note", "")        # keep what ours.tsv said
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--cache", type=Path, help="download cache (default: tools/.cache)")
    ap.add_argument("--offline", action="store_true", help="use the cache only")
    ap.add_argument("-o", "--out", type=Path, default=OUT)
    ap.add_argument("--refresh", action="store_true",
                    help="re-resolve every paper, not only those missing from --out")
    args = ap.parse_args()

    cache = args.cache or Path(__file__).resolve().parent / ".cache"
    cache.mkdir(parents=True, exist_ok=True)

    rows = build(cache, args.offline, {} if args.refresh else known_papers(args.out))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", encoding="utf-8") as fh:
        fh.write("\t".join(COLUMNS) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")).replace("\t", " ") for c in COLUMNS) + "\n")

    print(f"  ours: {len({r['our_dir'] for r in rows if r['our_dir']})} protocols "
          f"documented here, over {sum(1 for r in rows if r['our_dir'])} rows",
          file=sys.stderr)
    pri = sum(1 for r in rows if r["role"] == "primary")
    done = sum(1 for r in rows if r["documented"] == "yes")
    print(f"\n{args.out}: {len(rows)} rows, "
          f"{len({r['protocol'] for r in rows})} protocols", file=sys.stderr)
    print(f"  {done} rows from documented pages, {len(rows) - done} from the TODO list",
          file=sys.stderr)
    print(f"  {pri} primary (a method's own paper), {len(rows) - pri} incidental citations",
          file=sys.stderr)
    print(f"  with DOI: {sum(1 for r in rows if r['doi'])}; "
          f"with PMID: {sum(1 for r in rows if r['pmid'])}; "
          f"with title: {sum(1 for r in rows if r['title'])}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
