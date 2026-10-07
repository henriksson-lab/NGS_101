"""scg_lib_structs on GitHub (Teichlab/scg_lib_structs, branch master): the files behind
an upstream method page.

The rendered page (teichlab.github.io/scg_lib_structs/methods_html/X.html) is fetched by
get_sources itself. This handler adds what the page is built from and links to, which
for vendor kits (10x Genomics, BD Rhapsody, Bio-Rad ddSEQ) *is* the primary source:

  data/<dir>/...               vendor user guides (10x CG000*.pdf, BD GMX_*.pdf), oligo
                               tables (xlsx/xls), barcode whitelists (*.txt.gz, BD_CLS*.txt)
                               -- every ../data/ file the page links, plus the data/<dir>
                               whose name matches the page (e.g. STRT-seq_family)
  methods_html/<page>_*.html   sub-pages the page links (e.g. PIP-seq_v1p.html)
  docs/source/<area>/<m>.md    the seqspec-era Markdown source of the method
  chemistries/adapters/<m>.fa  adapter FASTA

Mapping is by the page's own links first, then by normalised name (see `scg_paths`).
Images are skipped; files over SIZE_CAP are recorded as (manual) with their raw URL.
Licensing: the repository has no LICENSE file; its Zenodo release
(doi:10.5281/zenodo.10042390) is CC-BY-4.0. Vendor documents mirrored under data/ remain
the vendors' copyright. Both facts go into MANIFEST.tsv as a (note) row.

Spec form for catalogue/extra_sources.tsv: `scg:<glob>` over repository paths, e.g.
`scg:data/10X-Genomics/CG000185*` or `scg:data/illumina-adapter-sequences-*.pdf`.
"""
from __future__ import annotations

import fnmatch
import json
import re
import time
import urllib.parse
from pathlib import Path

from . import register
from .base import cache_dir, get, say

REPO = "Teichlab/scg_lib_structs"
BRANCH = "master"
TREE_URL = f"https://api.github.com/repos/{REPO}/git/trees/{BRANCH}?recursive=1"
RAW = f"https://raw.githubusercontent.com/{REPO}/{BRANCH}/"
LICENSE_NOTE = ("scg_lib_structs (github.com/Teichlab/scg_lib_structs): no LICENSE file in the "
                "repository; its Zenodo release doi:10.5281/zenodo.10042390 is CC-BY-4.0. "
                "Vendor documents mirrored under data/ remain the vendors' copyright.")
SIZE_CAP = 40_000_000
IMAGE = re.compile(r"\.(png|jpe?g|gif|svg|tiff?|bmp|ico|eps)$", re.I)
TREE_MAX_AGE = 7 * 86400

# Upstream page stem -> repository files that the name rules below cannot find.
ALIASES = {
    "10xChromium_scATAC": ["docs/source/epi/10x_scATAC-seq.md"],
    "10xChromium_multiome": ["docs/source/multi/10x_multiome.md"],
    "10xChromium3": ["docs/source/ge/10xChromium3v2.md", "docs/source/ge/10xChromium3v3.md",
                     "chemistries/adapters/Chromium.fa"],
    "10xChromium3v1": ["chemistries/adapters/Chromium.fa"],
    "sci-RNA-seq_family": ["docs/source/ge/sci-RNA-seq3.md"],
    "SMART-seq_family": ["docs/source/ge/SMART-seq2.md", "chemistries/adapters/SMART_seq.fa"],
    "SCRB-seq": ["docs/source/ge/mcSCRB-seq.md"],
    "SureCell": ["data/illumina-adapter-sequences-1000000002694-14.pdf"],
    "Drop-seq": ["chemistries/adapters/Drop_seq_Seq_Well.fa"],
    "Seq-Well": ["chemistries/adapters/Drop_seq_Seq_Well.fa"],
    "Quartz-seq_family": ["data/Quartz-seq2_cell_barcodes.txt"],
}


def norm(s: str) -> str:
    return re.sub(r"[^a-z0-9]", "", s.lower())


def page_stem(scg_page: str) -> str:
    return urllib.parse.unquote(scg_page.rsplit("/", 1)[-1]).removesuffix(".html")


def page_links(page_html: str) -> list[str]:
    """Repository paths an upstream methods_html page links: ../data/... files and
    sibling methods_html pages."""
    out = []
    for href in re.findall(r'href="([^"#?]+)"', page_html):
        href = urllib.parse.unquote(href)
        if href.startswith("../data/"):
            out.append(href[3:])
        elif m := re.match(r"(?:https://teichlab\.github\.io/scg_lib_structs/)?"
                           r"(?:\.\./)?methods_html/([^/]+\.html)$", href):
            out.append("methods_html/" + m.group(1))
        elif re.match(r"^[^/:]+\.html$", href):        # relative sibling page
            out.append("methods_html/" + href)
    return list(dict.fromkeys(out))


def scg_paths(tree: list[str], stem: str, links: list[str]) -> list[str]:
    """The repository files that belong to the upstream page `stem` (e.g. '10xChromium3').
    Pure: `tree` is every blob path in the repository, `links` the page's own links."""
    blobs = set(tree)
    pages = {Path(p).stem for p in tree if p.startswith("methods_html/") and p.endswith(".html")}
    n = norm(stem)
    picked: list[str] = []

    def add(p: str) -> None:
        if p in blobs and p not in picked and not IMAGE.search(p):
            picked.append(p)

    for p in links:
        if p.startswith("data/"):
            add(p)
        elif p.startswith("methods_html/"):
            sub = Path(p).stem      # a sub-page (PIP-seq_v1p), not another method (10xChromium3v1)
            if sub.startswith(stem + "_"):
                add(p)
    # data/<dir>/ named like the page (STRT-seq_family, HyDrop for HyDrop_RNA, ...)
    dirs = {p.split("/")[1] for p in tree if p.startswith("data/") and p.count("/") >= 2}
    for d in sorted(dirs):
        nd = norm(d)
        if nd == n or (len(nd) >= 5 and n.startswith(nd)):
            for p in sorted(tree):
                if p.startswith(f"data/{d}/"):
                    add(p)
    # docs/source/<area>/<name>.md: exact name, or the page is a family of that name, or
    # the md is a version of the page (10xChromium3 -> ...v2/v3) that has no page of its own
    for p in sorted(tree):
        if not (p.startswith("docs/source/") and p.endswith(".md") and p.count("/") == 3):
            continue
        m = Path(p).stem
        nm = norm(m)
        if nm == n or (len(nm) >= 5 and n.startswith(nm)) or \
                (nm.startswith(n) and m not in pages and len(n) >= 5):
            add(p)
    for p in sorted(tree):
        if p.startswith("chemistries/adapters/") and p.endswith(".fa"):
            na = norm(Path(p).stem)
            if na == n or (len(na) >= 5 and n.startswith(na)):
                add(p)
    for p in ALIASES.get(stem, []):
        add(p)
    return picked


def local_name(path: str) -> str:
    """Repository path -> file name in the source folder."""
    if path.startswith("methods_html/"):
        return "upstream_" + Path(path).name
    if path.startswith("docs/source/"):
        return "scg_docs_" + path.removeprefix("docs/source/").replace("/", "_")
    if path.startswith("chemistries/"):
        return "scg_" + path.removeprefix("chemistries/").replace("/", "_")
    return "scg_" + Path(path).name


def tree(refresh: bool = False) -> dict[str, int]:
    """{path: size} of every blob, from the GitHub API (cached a week: the unauthenticated
    API allows 60 calls an hour)."""
    cache = cache_dir() / "scg_lib_structs_tree.json"
    fresh = cache.exists() and time.time() - cache.stat().st_mtime < TREE_MAX_AGE
    if not fresh or refresh:
        code, body = get(TREE_URL)
        try:
            data = json.loads(body) if code == 200 else {}
        except ValueError:
            data = {}
        if data.get("tree"):
            cache.write_bytes(body)
        elif not cache.exists():
            say(f"    !! GitHub tree for {REPO} unavailable (HTTP {code})")
            return {}
    data = json.loads(cache.read_bytes())
    return {e["path"]: e.get("size", 0) for e in data.get("tree", []) if e.get("type") == "blob"}


def _fetch_paths(job, paths: list[str], sizes: dict[str, int], what: str) -> None:
    f = job.folder
    f.note("scg_lib_structs licence", f"https://github.com/{REPO}", LICENSE_NOTE)
    for p in paths:
        url = RAW + urllib.parse.quote(p)
        name = local_name(p)
        desc = what or f"scg_lib_structs repo file {p}"
        if sizes.get(p, 0) > SIZE_CAP:
            f.manual(name, url, desc, f"{sizes[p] // 1_000_000} MB, over the {SIZE_CAP // 1_000_000} MB cap")
            continue
        f.save(name, url, desc, min_text=200)


def _page_html(job) -> str:
    for name in job.folder.rows:
        if name.startswith("upstream_") and name.endswith(".html") and (job.folder.path / name).exists():
            return (job.folder.path / name).read_text(errors="replace")
    return ""


@register("scg", url=r"^https://(?:raw\.githubusercontent\.com|github\.com)/Teichlab/scg_lib_structs/",
          auto=lambda job: "teichlab.github.io/scg_lib_structs" in (job.rows[0].get("scg_page") or ""),
          doc="scg_lib_structs GitHub repo: vendor guides, oligo tables, whitelists behind a page")
def fetch(job, arg: str, what: str = "") -> None:
    sizes = tree()
    if not sizes:
        return
    if arg.startswith("http"):                     # a raw/blob URL into the repo
        arg = re.sub(r"^https://[^/]+/Teichlab/scg_lib_structs/(?:blob/|raw/)?[^/]+/", "", arg)
    if arg:
        paths = [p for p in sorted(sizes) if fnmatch.fnmatch(p, arg) and not IMAGE.search(p)]
        if not paths:
            say(f"    !! scg:{arg} matches no file in {REPO}")
        _fetch_paths(job, paths, sizes, what)
        return
    stem = page_stem(job.rows[0]["scg_page"])
    paths = scg_paths(list(sizes), stem, page_links(_page_html(job)))
    say(f"  scg_lib_structs repo: {len(paths)} file(s) for {stem}")
    _fetch_paths(job, paths, sizes, what)
