#!/usr/bin/env python3
"""Fetch the readable source documents for a protocol into $CHEM_DATA. Never into the repo.

For each paper of a protocol (its rows in catalogue/scg_lib_structs.tsv, or a DOI given
directly), fetch whatever is legitimately reachable without a login or a browser:

  upstream    the scg_lib_structs method page, and -- from the Teichlab/scg_lib_structs
              GitHub repository -- the vendor user guides, oligo tables, barcode
              whitelists and Markdown/FASTA sources behind it (tools/fetchers/scg_github.py)
  PMC         JATS full text, PDF and (open-access subset / CC author manuscripts)
              supplements from the PMC Cloud bucket (tools/fetchers/pmc_cloud.py); else
              Europe PMC full text + supplementary zip; else the PMC article page
  bioRxiv     the full-text PDF and every supplementary file of the latest version
              (DOIs 10.1101/YYYY.MM.DD.NNNNNN, 10.64898/YYYY.MM.DD.NNNNNN or
              10.1101/NNNNNN; other 10.1101/ DOIs are
              Cold Spring Harbor journals -- Genome Research, Genes & Dev -- and go the
              journal route)
  linked      a preprint's published version and a journal paper's preprint (Crossref
              relations, the bioRxiv API), fetched the same way
  Springer    supplementary files of Nature-family papers (static-content.springer.com)
  fetchers    special handlers in tools/fetchers/ (eLife API, protocols.io, vendor URLs,
              the scg_lib_structs repo, ...), triggered by DOI/URL pattern or declared
              per protocol in catalogue/extra_sources.tsv

For a protocol directory in catalogue/ours.tsv, its defining DOI is always fetched; a
scg catalogue paper with a different DOI is fetched only when a preprint relation or
agreeing Crossref title/authors links it to that defining DOI (else a (note) row says why).

Source identity: a cache name already registered to another URL is never reused or
overwritten; a byte-identical download from the other URL (a publisher mirror) is recorded
as a `(duplicate)` row of the registered file, a different one as `(manual)`. A second
name for an already saved URL, and an identical copy of a file already here, become
`(duplicate)` rows; a URL-less hand-placed row gets the URL of an identical download;
scg_lib_structs copies of vendor files stay under their scg_ name beside upstream_* twins.

Every download is validated before it is kept (a PDF must be a whole PDF, an xlsx a whole
zip, an HTML page real text and not a captcha/"Preparing to download" gate); what fails
is recorded in MANIFEST.tsv as `(manual) NAME` with its URL and the reason. Each document
gets its plain-text twin FILE.txt (tools/doctext.py) as soon as it is saved -- read, grep
and cite that. Re-running skips what is on disk. Everything lands in
$CHEM_DATA/sources/<slug>/ (default /husky/henriksson/NGS101_papers/sources/<slug>/):
third-party material stays outside the repository, with text twins beside originals.

    python3 tools/get_sources.py "SPLiT-seq / microSPLiT"     # a protocol, by name
    python3 tools/get_sources.py splitseq                     # ...or a fragment of its slug
    python3 tools/get_sources.py --slugs a__10.1+x,b          # exact slugs (comma-separated)
    python3 tools/get_sources.py --doi 10.1101/829960         # one paper
    python3 tools/get_sources.py --todo                       # every untackled protocol
    python3 tools/get_sources.py --all-ours                   # every protocol in ours.tsv
    python3 tools/get_sources.py --revalidate                 # re-check files on disk
    python3 tools/get_sources.py --list                       # what is already fetched
    python3 tools/get_sources.py --fetchers                   # registered special handlers
"""
from __future__ import annotations

import argparse
import html
import re
import shutil
import sys
import urllib.parse
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "catalogue" / "tools"))
sys.path.insert(0, str(ROOT / "tools"))

import catalogue as cat  # noqa: E402
import fetchers  # noqa: E402
from fetchers import Job, handlers_for_doi, protocolsio, run_auto, run_extra  # noqa: E402
from fetchers.base import (Folder, data_dir, gate_reason, get, get_json, is_twin,  # noqa: E402
                           md5, safe_name, say, validate)
from fetchers.pmc_cloud import fetch_article as pmc_cloud_article  # noqa: E402
from fetchers.pmc_cloud import settle_manual_rows as settle_pmc_rows  # noqa: E402

EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest"
SCG_CACHE = ROOT / "catalogue" / "tools" / ".cache"
MAX_LINKED_PROTOCOLS = 6

BIORXIV_DOI = re.compile(r"^10\.(?:1101|64898)/\d{4}\.\d{2}\.\d{2}\.\d{5,6}(?:v\d+)?$"
                         r"|^10\.1101/\d{6}(?:v\d+)?$")


# ---------------------------------------------------------------------- routing

def doi_kind(doi: str) -> str:
    """'biorxiv' | 'protocolsio' | 'researchsquare' | 'journal'. Offline, by pattern:
    10.1101/ is Cold Spring Harbor's prefix, shared by bioRxiv/medRxiv (dated or 6-digit
    suffixes) and its journals (gr., gad., rna., lm., pdb., cshperspect...). Newer
    bioRxiv preprints use 10.64898/ with a dated suffix."""
    d = doi.strip().lower()
    if BIORXIV_DOI.match(d):
        return "biorxiv"
    if d.startswith("10.17504/protocols.io."):
        return "protocolsio"
    if d.startswith("10.21203/rs."):
        return "researchsquare"
    return "journal"


def doi_tag(doi: str) -> str:
    return re.sub(r"[^\w.-]", "_", doi.split("/", 1)[1])[:40]


def norm_doi(doi: str) -> str:
    return re.sub(r"^(?:https?://(?:dx\.)?doi\.org/|doi:)", "", doi.strip(), flags=re.I)


# -------------------------------------------------------------------- per source

def upstream(folder: Folder, scg_page: str) -> None:
    if not scg_page:
        return
    name = "upstream_" + scg_page.rsplit("/", 1)[-1]
    if folder.have(name) and not validate(folder.path / name, name, 500):
        return
    cached = SCG_CACHE / re.sub(r"[^A-Za-z0-9._-]", "_", scg_page)[-180:]
    if cached.exists() and not validate(cached, name, 500):
        folder.keep(name, cached.read_bytes(), scg_page, "scg_lib_structs method page")
    else:
        folder.save(name, scg_page, "scg_lib_structs method page", min_text=500)


def epmc_record(query: str) -> dict:
    q = urllib.parse.quote(query)
    res = get_json(f"{EPMC}/search?query={q}&format=json&resultType=core&pageSize=5", retries=2)
    hits = res.get("resultList", {}).get("result", [])
    return hits[0] if hits else {}


def supplement_names(folder: Folder, tag: str, pmcid: str) -> list[str]:
    """Supplementary file names of a PMC article, from its JATS XML or its article page."""
    names: list[str] = []
    xml = folder.path / f"{tag}_{pmcid}.xml"
    if xml.exists():
        t = xml.read_text(errors="replace")
        for block in re.findall(r"<supplementary-material\b.*?(?:</supplementary-material>|/>)", t, re.S):
            names += re.findall(r'xlink:href="([^"]+)"', block)
    page = folder.path / f"{tag}_{pmcid}.html"
    if page.exists():
        names += [h.rsplit("/", 1)[-1] for h in
                  re.findall(r'href="(?:https://pmc\.ncbi\.nlm\.nih\.gov)?(/articles/instance/\d+/bin/[^"]+)"',
                             page.read_text(errors="replace"))]
    return [n for n in dict.fromkeys(html.unescape(n) for n in names) if "/" not in n and "." in n]


def pmc(folder: Folder, tag: str, pmcid: str, open_access: bool) -> bool:
    """-> True if a full text (JATS XML with a body, or the article page) was saved."""
    cloud = pmc_cloud_article(folder, tag, pmcid)
    got_text = cloud["xml"]
    if not got_text:
        got_text = bool(folder.save(f"{tag}_{pmcid}.xml", f"{EPMC}/{pmcid}/fullTextXML",
                                    "full-text JATS XML (Europe PMC)", need=rb"<body"))
    if not got_text:
        got_text = bool(folder.save(f"{tag}_{pmcid}.html", f"https://pmc.ncbi.nlm.nih.gov/articles/{pmcid}/",
                                    "PMC article page: full text", min_text=5000))
        if got_text:
            folder.rows.pop(f"(manual) {tag}_{pmcid}.xml", None)
    if got_text:      # an older run's failed attempts at the same text are moot now
        folder.rows.pop(f"(manual) {tag}_{pmcid}.xml", None)
        folder.rows.pop(f"(manual) {tag}_{pmcid}.html", None)
    if cloud["supp"]:
        folder.rows.pop(f"(manual) {tag}_{pmcid}_supplementary.zip", None)
    if open_access and not cloud["supp"]:
        folder.save(f"{tag}_{pmcid}_supplementary.zip", f"{EPMC}/{pmcid}/supplementaryFiles",
                    "supplementary files (Europe PMC zip)")
        z = folder.path / f"{tag}_{pmcid}_supplementary.zip"
        if z.exists():
            folder.expand_zip(z.name, f"{tag}_supp_", f"{EPMC}/{pmcid}/supplementaryFiles",
                              "supplementary file")
    # supplements still missing: PMC serves them behind a download gate. Try one; if the
    # gate is up, list the rest for hand download without hammering the server.
    num = pmcid.removeprefix("PMC")
    lic = cloud["meta"].get("license_code") if cloud["meta"] else None
    for n in supplement_names(folder, tag, pmcid):
        name = safe_name(f"{tag}_supp_{n}")
        if (folder.path / name).exists() or f"(duplicate) {name}" in folder.rows:
            continue
        url = f"https://pmc.ncbi.nlm.nih.gov/articles/instance/{num}/bin/{urllib.parse.quote(n)}"
        if "pmc.ncbi.nlm.nih.gov" in folder.gated_hosts:
            folder.manual(name, url, "PMC supplementary file",
                          "PMC download gate" + (f" (not in the PMC Cloud bucket: licence {lic})" if lic else ""))
        else:
            folder.save(name, url, "PMC supplementary file", retries=0)
    return got_text


def biorxiv_details(doi: str) -> tuple[str, list[dict]]:
    for server in ("biorxiv", "medrxiv"):
        info = get_json(f"https://api.biorxiv.org/details/{server}/{doi}", retries=2)
        if info.get("collection"):
            return server, info["collection"]
    return "", []


def biorxiv(folder: Folder, tag: str, doi: str) -> tuple[bool, str]:
    """-> (full text saved?, the published version's DOI or '')."""
    server, vers = biorxiv_details(doi)
    if not vers:
        folder.manual(f"{tag}.full.pdf", f"https://doi.org/{doi}", "preprint",
                      "not found in the bioRxiv/medRxiv API")
        return False, ""
    last = vers[-1]
    v = last["version"]
    published = last.get("published", "")
    published = "" if published in ("", "NA") else published
    base = f"https://www.{server}.org/content/{doi}v{v}"
    ok = bool(folder.save(f"{tag}_v{v}.full.pdf", f"{base}.full.pdf", f"{server} full text, v{v}"))
    if ok:          # an earlier run's "not found in the API" row is moot now
        folder._clear_manual(f"{tag}.full.pdf", f"https://doi.org/{doi}")
    code, body = get(f"{base}.supplementary-material")
    gate = gate_reason(body[:65536])
    if code != 200 or gate:
        folder.manual(f"{tag}_supplementary", f"{base}.supplementary-material",
                      "supplementary page", gate or f"HTTP {code}")
        return ok, published
    links = re.findall(rb'href="((?:https://www\.(?:bio|med)rxiv\.org)?/content/[^"]+/media-\d+[^"]*'
                       rb'|(?:https://www\.(?:bio|med)rxiv\.org)?/highwire/filestream/[^"]+)"', body)
    for url in dict.fromkeys(links):
        url = html.unescape(url.decode())
        if url.startswith("/"):
            url = f"https://www.{server}.org{url}"
        fname = re.sub(r"\?.*", "", url.rsplit("/", 1)[-1])
        folder.save(f"{tag}_{fname}", url, f"{server} supplementary file, v{v}")
    return ok, published


def springer(folder: Folder, tag: str, doi: str) -> list[str]:
    """Nature-family supplements: listed on the article page, served from static-content."""
    code, body = get(f"https://www.nature.com/articles/{doi.split('/', 1)[1]}")
    if code != 200:
        return []
    links = re.findall(rb'https://(?:static-content\.springer\.com|media\.springernature\.com/'
                       rb'original/springer-static)/esm/[^"\s]+?_ESM\.\w+', body)
    saved = []
    for url in sorted(set(links)):
        url = html.unescape(url.decode())
        if folder.save(f"{tag}_" + url.rsplit("/", 1)[-1], url, "Springer supplementary file"):
            saved.append(url)
    return saved


STOP = {"a", "an", "the", "of", "in", "and", "for", "with", "by", "to", "on", "at", "from",
        "single", "cell", "cells", "sequencing", "seq", "using", "via", "high", "throughput"}


def _similar(a: str, b: str) -> bool:
    """Same paper, judged on the distinctive words -- field vocabulary like 'single-cell
    sequencing' is shared by half the catalogue and says nothing."""
    wa = set(re.findall(r"[a-z0-9]+", a.lower())) - STOP
    wb = set(re.findall(r"[a-z0-9]+", b.lower())) - STOP
    return bool(wa) and len(wa & wb) / len(wa | wb) >= 0.6


def _year(item: dict) -> int:
    for k in ("posted", "issued", "created"):
        parts = (item.get(k) or {}).get("date-parts") or [[0]]
        if parts[0] and parts[0][0]:
            return int(parts[0][0])
    return 0


def crossref_relations(doi: str) -> dict:
    return get_json(f"https://api.crossref.org/works/{urllib.parse.quote(doi)}",
                    retries=1).get("message", {}).get("relation", {}) or {}


def related_dois(rel: dict, kind: str) -> list[str]:
    return [r["id"] for r in rel.get(kind, []) if r.get("id-type") == "doi" and r.get("id")]


def preprints_of(doi: str, title: str = "", year: str = "", search: bool = False,
                 first_author: str = "") -> list[str]:
    """Preprint DOIs of a journal paper: Crossref has-preprint, the bioRxiv publication
    index; with `search`, a title search in Europe PMC / Crossref as a last resort."""
    found = related_dois(crossref_relations(doi), "has-preprint")
    pubs = get_json(f"https://api.biorxiv.org/pubs/biorxiv/{doi}")
    found += [c.get("preprint_doi", "") for c in pubs.get("collection", []) if c.get("preprint_doi")]
    if not found and search and title:
        found += [d for d in [preprint_by_title(doi, title, year, first_author)] if d]
    return list(dict.fromkeys(d.lower() for d in found))


def _surname(name: str) -> str:
    """'Li Y' / 'Li, Yun' / 'Yun Li' (catalogue first_author is a surname) -> 'li'."""
    name = name.strip().split(",")[0].strip()
    parts = name.split()
    if len(parts) > 1 and len(parts[-1]) <= 3 and parts[-1].isupper():   # 'Li Y', 'Smith JA'
        parts = parts[:-1]
    return re.sub(r"[^a-z]", "", " ".join(parts).lower())


def _same_first_author(want: str, got: str) -> bool:
    """A title match also needs the same first author when both are known."""
    a, b = _surname(want), _surname(got)
    return not a or not b or a == b or a.endswith(b) or b.endswith(a)


def preprint_by_title(doi: str, title: str, year: str = "", first_author: str = "") -> str:
    rel = get_json("https://api.crossref.org/works?rows=5&filter=relation.type:is-preprint-of,"
                   f"relation.object:{urllib.parse.quote(doi)}")
    for it in rel.get("message", {}).get("items", []):
        if doi_kind(it.get("DOI", "")) == "biorxiv":
            return it["DOI"]
    words = " ".join(re.findall(r"[A-Za-z0-9]+", title)[:10])
    y = int(year) if year.isdigit() else 0
    ok_year = (lambda py: not y or not py or y - 2 <= py <= y)
    rec = epmc_record(f'SRC:PPR AND TITLE:"{words}"')
    if (doi_kind(rec.get("doi", "")) == "biorxiv" and _similar(title, rec.get("title", ""))
            and ok_year(int(rec.get("pubYear") or 0))
            and _same_first_author(first_author, rec.get("authorString", "").split(",")[0])):
        return rec["doi"]
    res = get_json("https://api.crossref.org/works?rows=5&filter=type:posted-content"
                   f"&query.title={urllib.parse.quote(title)}")
    for it in res.get("message", {}).get("items", []):
        lead = ((it.get("author") or [{}])[0]).get("family", "")
        if (doi_kind(it.get("DOI", "")) == "biorxiv"
                and _similar(title, (it.get("title") or [""])[0]) and ok_year(_year(it))
                and _same_first_author(first_author, lead)):
            return it["DOI"]
    return ""


_WORKS: dict[str, dict] = {}


def crossref_work(doi: str) -> dict:
    """Crossref metadata of a DOI ({} if unavailable), cached for the run."""
    key = doi.lower()
    if key not in _WORKS:
        _WORKS[key] = get_json(f"https://api.crossref.org/works/{urllib.parse.quote(doi)}",
                               retries=1).get("message", {}) or {}
    return _WORKS[key]


def works_agree(a: dict, b: dict) -> str:
    """Why two Crossref records describe the same work or a companion of it ('' if they
    do not): a similar title, or a shared first/last author within three years. A
    shared method name alone (CAP-seq 2010 vs CAP-seq 2024) is not agreement."""
    if not a or not b:
        return ""
    ta, tb = (a.get("title") or [""])[0], (b.get("title") or [""])[0]
    if ta and tb and _similar(ta, tb):
        return "similar title"
    names = lambda w: [_surname(x.get("family", "")) for x in w.get("author") or []  # noqa: E731
                       if x.get("family")]
    na, nb = names(a), names(b)
    ends = lambda n: {n[0], n[-1]} if n else set()  # noqa: E731
    ya, yb = _year(a), _year(b)
    shared = (ends(na) & set(nb)) | (ends(nb) & set(na))
    if shared and ya and yb and abs(ya - yb) <= 3:
        return f"shared author {sorted(shared)[0]}, {ya} vs {yb}"
    return ""


def linked_to_defining(defining: str, doi: str) -> str:
    """Why a catalogue paper DOI belongs with a protocol defined by `defining` ('' if no
    evidence): the same DOI, a Crossref/bioRxiv preprint relation, or agreeing metadata."""
    if norm_doi(doi).lower() == norm_doi(defining).lower():
        return "defining DOI"
    d, p = norm_doi(defining).lower(), norm_doi(doi).lower()
    if p in preprints_of(d) or d in published_of(p):
        return "preprint relation"
    return works_agree(crossref_work(d), crossref_work(p))


def published_of(doi: str, hint: str = "") -> list[str]:
    """Journal versions of a preprint: the bioRxiv API's 'published', Crossref is-preprint-of."""
    found = [hint] if hint else []
    found += related_dois(crossref_relations(doi), "is-preprint-of")
    return list(dict.fromkeys(d.lower() for d in found if d))


def paper(job: Job, doi: str, title: str = "", year: str = "", pmcid: str = "",
          link: bool = True, first_author: str = "") -> None:
    """Fetch one paper the generic way, plus its linked preprint/published version."""
    doi = norm_doi(doi)
    if not doi or doi.lower() in job.seen:
        return
    job.seen.add(doi.lower())
    folder = job.folder
    tag = doi_tag(doi)
    old = folder.rows.get(f"(manual) {tag}")          # older code's catch-all row
    if old and old[4] == "no open copy found: publisher page only":
        del folder.rows[f"(manual) {tag}"]
    kind = doi_kind(doi)
    say(f"  {doi}  [{kind}]  {title[:60]}")
    for h in handlers_for_doi(doi):
        h.fn(job, doi, "")
    if kind == "protocolsio":
        return
    if kind == "researchsquare":
        url = f"https://doi.org/{doi}"
        if not settle_researchsquare(folder, doi) and not (
                f"(manual) {tag}" in folder.rows and folder.rows[f"(manual) {tag}"][1] == url):
            folder.manual(tag, url, "Research Square preprint / Protocol Exchange",
                          "Research Square is behind Cloudflare: fetch by hand")
        return
    if kind == "biorxiv":
        got_full_text, pub = biorxiv(folder, tag, doi)
        if got_full_text:
            # A dated DOI once routed through the journal fallback may have left
            # this obsolete manual row in an existing source manifest.
            folder.rows.pop(f"(manual) {tag} full text", None)
        for d in (published_of(doi, pub) if link else []):
            say(f"    published version {d}")
            paper(job, d, link=False)
        return
    rec = epmc_record(f"DOI:{doi}")
    pmcid = rec.get("pmcid") or pmcid
    got_text = False
    if pmcid:
        got_text = pmc(folder, tag, pmcid, rec.get("isOpenAccess") == "Y")
    esm = springer(folder, tag, doi) if doi.startswith("10.1038/") else []
    if esm:   # PMC's gated author-manuscript copies are then usually redundant
        for k, r in folder.rows.items():
            if k.startswith(f"(manual) {tag}_supp_") and "publisher" not in r[4]:
                r[4] += "; the publisher's ESM files were fetched, so probably redundant"
    if link:
        got_supp = any(k.startswith(tag) and "supp" in k.lower() and not k.startswith("(")
                       for k in folder.rows)
        for d in preprints_of(doi, title, year, search=not pmcid or not got_supp,
                              first_author=first_author):
            say(f"    preprint {d}")
            paper(job, d, link=False)
    have_text = got_text or any(k.startswith(tag) and k.endswith((".xml", ".html", ".pdf"))
                                and not k.startswith("(") and "supp" not in k and "ESM" not in k
                                and (folder.path / k).exists() for k in folder.rows)
    url, key = f"https://doi.org/{doi}", f"(manual) {tag} full text"
    copy = "" if have_text else article_printing_doi(folder, doi)
    if have_text or copy:   # an earlier run's failed attempt at the same full text is moot now
        row = folder.rows.get(key)
        folder._clear_manual(f"{tag} full text", url)
        if copy and row and row[1] == url:   # keep a hand annotation of how it was obtained
            folder.note(f"{tag} full text", url, f"{row[4]}; full text on disk as {copy}, "
                        "which prints this DOI on its first pages")
    elif not (key in folder.rows and folder.rows[key][1] == url):   # keep annotations
        folder.manual(f"{tag} full text", url, "journal article",
                      "no open full text found (publisher page only)")
    folder.write()


def settle_researchsquare(folder, doi: str) -> str:
    """A Research Square / Protocol Exchange DOI whose document is on disk (e.g. fetched from
    its asset URL by an extra-source row) needs no '(manual)' row: replace that row by a note
    keeping its annotation. -> the on-disk file name, or ''."""
    tag, url = doi_tag(doi), f"https://doi.org/{doi}"
    copy = article_printing_doi(folder, doi)
    if not copy:
        return ""
    row = folder.rows.get(f"(manual) {tag}")
    row = row if row and row[1] == url else None
    folder._clear_manual(tag, url)
    if row or f"(note) {tag}" not in folder.rows:   # keep an earlier hand annotation
        folder.note(tag, url, (f"{row[4]}; " if row else "")
                    + f"on disk as {copy}, which prints this DOI on its first pages")
    folder.write()
    return copy


ARTICLE_TEXT = 10000    # bytes of extracted text: correction notices are shorter, protocols longer


def article_printing_doi(folder, doi: str) -> str:
    """A validated on-disk original under any name (e.g. an Elsevier mmcN article-plus-SI
    attachment) that prints `doi` near its start -> its name, or ''. Supplements named as
    such are excluded; the DOI must appear in the first pages, as an article prints its own,
    and the text must be article-length (a correction notice also prints the DOI it corrects).
    Another version that cites this DOI is not it: a JATS file must declare it as its own
    article-id, and a document printing a preprint-server DOI is a preprint (which may
    print its published version's DOI) unless that DOI is the requested one."""
    want = doi.lower()
    for name, row in folder.rows.items():
        if (name.startswith("(") or not name.lower().endswith((".pdf", ".xml", ".html"))
                or "supp" in name.lower() or "ESM" in name or not row[3]):
            continue
        path, twin = folder.path / name, folder.path / (name + ".txt")
        if not (path.is_file() and twin.is_file() and twin.stat().st_size >= ARTICLE_TEXT):
            continue
        with twin.open(errors="replace") as fh:
            head = fh.read(6000).lower()
        if want not in head:
            continue
        if name.lower().endswith(".xml"):
            with path.open("rb") as fh:
                own = re.findall(rb'<article-id pub-id-type="doi">\s*([^<\s]+)', fh.read(200_000))
            if want.encode() not in [d.lower() for d in own]:
                continue
        elif doi_kind(want) not in PREPRINT_KINDS and any(
                doi_kind(d) in PREPRINT_KINDS for d in re.findall(r"10\.\d{4,9}/[^\s\"'<>]+", head)):
            continue
        return name
    return ""


PREPRINT_KINDS = {"biorxiv", "researchsquare"}


def linked_protocols(job: Job) -> None:
    """protocols.io protocols the fetched papers cite."""
    refs: list[str] = []
    for t in sorted(job.folder.path.glob("*.txt")):
        if is_twin(t) and not t.name.startswith(("upstream_", "scg_", "protocolsio_")):
            refs += protocolsio.find_refs(t.read_text(errors="replace"))
    refs = [r for r in dict.fromkeys(refs) if f"protocolsio:{r}" not in job.seen]
    for r in refs[:MAX_LINKED_PROTOCOLS]:
        say(f"  cited protocols.io: {r}")
        fetchers.run_spec(job, f"protocolsio:{r}", "protocols.io protocol cited by a fetched paper")


# ------------------------------------------------------------------------ driver

def protocol_rows(query: str) -> tuple[str, list[dict]]:
    catalogue_rows = cat.rows()
    our_dirs = {r["dir"] for r in cat.ours()}
    if query in our_dirs:
        rows = [r for r in catalogue_rows if r["our_dir"] == query]
        if rows:
            return query, rows
    slugs = {r["slug"] for r in catalogue_rows}
    if query in slugs:
        rows = [r for r in catalogue_rows if r["slug"] == query]
        names = list(dict.fromkeys(r["protocol"] for r in rows))
        return query, [r for r in catalogue_rows if r["protocol"] == names[0]]
    q = re.sub(r"[^a-z0-9]", "", query.lower())
    flat = lambda s: re.sub(r"[^a-z0-9]", "", s.lower())  # noqa: E731
    names = [p for p in cat.protocols() if q == flat(p)] or \
            [p for p in cat.protocols()
             if q in flat(p) or q in flat(next(r["slug"] for r in catalogue_rows if r["protocol"] == p))]
    if len(names) != 1:
        raise SystemExit(f"{query!r} matches {len(names)} protocols: {names[:10]}")
    rows = [r for r in catalogue_rows if r["protocol"] == names[0]]
    return rows[0]["slug"], rows


def defining_doi(slug: str) -> str:
    """The defining DOI catalogue/ours.tsv names for a protocol directory ('' if none)."""
    return next((o.get("doi", "") for o in cat.ours() if o["dir"] == slug), "")


def fetch_protocol(slug: str, rows: list[dict], cited: bool = False, papers: bool = True) -> None:
    say(f"{rows[0]['protocol']}  ->  {data_dir() / slug}")
    folder = Folder(data_dir() / slug)
    job = Job(slug, rows, folder)
    # An extra-source doi: row names one declared document. Its own preprint/published
    # versions are not this protocol's papers (a companion method's preprint would be
    # attributed here), so they are not followed; declare them explicitly if needed.
    job.paper = lambda d: paper(job, d, link=False)
    job.pmc = lambda p: pmc(folder, p.upper(), p.upper(),
                            epmc_record(p).get("isOpenAccess") == "Y")
    upstream(folder, rows[0]["scg_page"])
    run_auto(job)
    if papers:
        roles = {"primary", "cited"} if cited else {"primary"}
        defining = defining_doi(slug)
        if defining:      # ours.tsv's defining paper, whatever DOI the scg catalogue row has
            own = next((r for r in rows if r["doi"].lower() == defining.lower()), {})
            paper(job, defining, own.get("title", ""), own.get("year", ""), own.get("pmcid", ""),
                  first_author=own.get("first_author", ""))
        for r in sorted(rows, key=lambda r: (r["is_defining"] != "yes", r["role"] != "primary")):
            if r["role"] not in roles or not r["doi"]:
                continue
            if defining and r["role"] == "primary":
                why = linked_to_defining(defining, r["doi"])
                if not why:   # e.g. another method sharing the name: not this protocol's paper
                    folder.note(f"unlinked catalogue paper {r['doi']}", f"https://doi.org/{r['doi']}",
                                f"not fetched: catalogue row '{r['protocol']}' has no preprint "
                                f"relation or agreeing title/authors with the defining DOI {defining}")
                    say(f"  skip {r['doi']}: no link to the defining DOI {defining}")
                    continue
                folder.rows.pop(f"(note) unlinked catalogue paper {r['doi']}", None)
            paper(job, r["doi"], r["title"], r.get("year", ""), r.get("pmcid", ""),
                  first_author=r.get("first_author", ""))
    run_extra(job)
    if papers:
        linked_protocols(job)
    protocolsio.settle_manual_rows(folder)
    settle_pmc_rows(folder)
    for key, row in list(folder.rows.items()):   # Research Square documents fetched by asset URL
        if key.startswith("(manual) ") and doi_kind(norm_doi(row[1])) == "researchsquare" \
                and key == f"(manual) {doi_tag(norm_doi(row[1]))}":
            settle_researchsquare(folder, norm_doi(row[1]))
    folder.write()
    n_ok = sum(1 for k in folder.rows if not k.startswith("("))
    n_man = sum(1 for k in folder.rows if k.startswith("(manual)"))
    say(f"  done: {n_ok} files, {n_man} to fetch by hand (see MANIFEST.tsv)\n")


def revalidate(dirs: list[Path]) -> None:
    """Re-check every top-level file already on disk (fetched by older code or by hand).
    A file that is not what its name says (gate page, truncated zip/PDF) moves to
    _invalid/ with its .txt twin, and its manifest row becomes (manual) with the reason.
    Valid files with no row are adopted; stale (manual) rows for present files dropped;
    missing or raw-XML .txt twins are (re)made."""
    import doctext
    for d in dirs:
        if not (d / "MANIFEST.tsv").exists() and not any(d.iterdir()):
            continue
        folder = Folder(d)
        bad = adopted = remade = 0
        for f in sorted(p for p in d.iterdir() if p.is_file()):
            if f.name in ("MANIFEST.tsv",) or f.name.endswith(".part") or is_twin(f):
                continue
            row = folder.rows.get(f.name)
            need = rb"<body" if f.suffix.lower() == ".xml" and row and "full" in row[4].lower() else None
            reason = validate(f, f.name, min_text=1500 if "article" in (row[4] if row else "").lower()
                              else 300, need=need)
            if reason:
                inv = d / "_invalid"
                inv.mkdir(exist_ok=True)
                f.replace(inv / f.name)
                twin = f.with_name(f.name + ".txt")
                if twin.exists():
                    twin.replace(inv / twin.name)
                folder.rows.pop(f.name, None)
                if row and row[1].startswith("http"):     # something to fetch again by hand
                    folder.rows[f"(manual) {f.name}"] = [f"(manual) {f.name}", row[1], "", "",
                                                         f"{row[4]}; invalid on disk: {reason}"]
                say(f"  {d.name}: demoted {f.name}: {reason}")
                bad += 1
                continue
            if not row:
                folder.rows[f.name] = [f.name, "", str(f.stat().st_size), md5(f),
                                       "(on disk, not fetched by get_sources: added by hand?)"]
                adopted += 1
            folder.rows.pop(f"(manual) {f.name}", None)
            twin = f.with_name(f.name + ".txt")
            stale_xml = (f.suffix.lower() == ".xml" and twin.exists()
                         and twin.read_bytes()[:200].lstrip().startswith(b"<"))
            if f.suffix.lower() in doctext.CONVERT and (not twin.exists() or stale_xml):
                if stale_xml:
                    twin.unlink()
                folder.to_text(f)
                remade += twin.exists()
        protocolsio.settle_manual_rows(folder)
        settle_pmc_rows(folder)
        # rows of older code: "(manual) <tag>  no open copy found" for papers now handled
        for k in [k for k, r in folder.rows.items()
                  if k.startswith("(manual) ") and r[4] == "no open copy found: publisher page only"]:
            tag, url = k.removeprefix("(manual) "), folder.rows.pop(k)[1]
            if not tag.startswith("protocols.io") and \
                    not any(n.startswith(tag) and n.endswith((".xml", ".html", ".pdf"))
                            and not n.startswith("(") and "supp" not in n.lower()
                            and "ESM" not in n for n in folder.rows) \
                    and f"(manual) {tag} full text" not in folder.rows:
                folder.rows[f"(manual) {tag} full text"] = [
                    f"(manual) {tag} full text", url, "", "",
                    "journal article; no open full text found (publisher page only)"]
        folder.write()
        say(f"{d.name:70s} {bad} demoted, {adopted} adopted, {remade} text twins made")


def list_dirs() -> None:
    for d in sorted(p for p in data_dir().glob("*") if p.is_dir() and not p.name.startswith(".")):
        files = [f for f in d.iterdir() if f.is_file() and f.name != "MANIFEST.tsv" and not is_twin(f)]
        rows = (d / "MANIFEST.tsv").read_text(errors="replace").splitlines()[1:] \
            if (d / "MANIFEST.tsv").exists() else []
        manual = sum(1 for ln in rows if ln.startswith("(manual)"))
        dup = sum(1 for ln in rows if ln.startswith("(duplicate)"))
        invalid = len(list((d / "_invalid").glob("*"))) if (d / "_invalid").is_dir() else 0
        invalid -= sum(1 for f in (d / "_invalid").glob("*.txt")) if invalid else 0
        say(f"{d.name:68s} {len(files):3d} files  {manual:3d} manual  {invalid:2d} invalid"
            f"  {dup:2d} dup")


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 epilog=__doc__.split("\n", 2)[2])
    ap.add_argument("protocol", nargs="*", help="protocol name or slug (fragment)")
    ap.add_argument("--doi", action="append", default=[], help="fetch one paper by DOI")
    ap.add_argument("--slugs", default="", help="comma-separated exact catalogue slugs")
    ap.add_argument("--todo", action="store_true", help="every protocol not yet tackled here")
    ap.add_argument("--all-ours", action="store_true", help="every protocol in catalogue/ours.tsv")
    ap.add_argument("--cited", action="store_true", help="also fetch role=cited papers")
    ap.add_argument("--no-papers", action="store_true",
                    help="only the upstream page, repo files and extra sources (no papers)")
    ap.add_argument("--revalidate", action="store_true",
                    help="re-check files on disk (all dirs, or the protocols named)")
    ap.add_argument("--list", action="store_true", help="list what is already fetched")
    ap.add_argument("--fetchers", action="store_true", help="list the special-case fetchers")
    a = ap.parse_args(argv)
    try:
        sys.stdout.reconfigure(line_buffering=True)
    except AttributeError:
        pass

    if a.list:
        list_dirs()
        return 0
    if a.fetchers:
        for h in fetchers.load().values():
            pats = ", ".join(x for x in (h.doi and f"doi {h.doi.pattern}", h.url and f"url {h.url.pattern}",
                                         h.auto and "auto") if x)
            say(f"{h.name:12s} {h.doc}\n{'':12s} [{pats or 'spec only'}]")
        return 0

    targets: list[tuple[str, list[dict]]] = [protocol_rows(q) for q in a.protocol]
    targets += [protocol_rows(s.strip()) for s in a.slugs.split(",") if s.strip()]
    if a.todo:
        targets += [(r["slug"], [x for x in cat.rows() if x["protocol"] == r["protocol"]])
                    for r in cat.untackled()]
    if a.all_ours:
        targets += [protocol_rows(o["dir"]) for o in cat.ours()]
    seen, uniq = set(), []
    for t in targets:
        if t[0] not in seen:
            seen.add(t[0])
            uniq.append(t)

    if a.revalidate:
        dirs = [data_dir() / s for s, _ in uniq] or \
               sorted(p for p in data_dir().glob("*") if p.is_dir() and not p.name.startswith("."))
        revalidate([d for d in dirs if d.is_dir()])
        return 0
    for doi in a.doi:
        rows = [r for r in cat.rows() if r["doi"].lower() == doi.lower()]
        row = rows[0] if rows else {"doi": doi, "title": "", "pmcid": "", "year": ""}
        folder = Folder(data_dir() / ("doi_" + doi.replace("/", "+")))
        job = Job("doi_" + doi.replace("/", "+"), [row], folder)
        job.paper = lambda d, job=job: paper(job, d)
        paper(job, doi, row["title"], row.get("year", ""), row.get("pmcid", ""),
              first_author=row.get("first_author", ""))
        folder.write()
    for slug, rows in uniq:
        fetch_protocol(slug, rows, a.cited, not a.no_papers)
    if not (a.doi or uniq):
        ap.print_help()
    return 0


if __name__ == "__main__":
    import signal
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)      # `| head` ends quietly
    shutil.which("curl") or sys.exit("get_sources.py needs curl")
    sys.exit(main(sys.argv[1:]))
