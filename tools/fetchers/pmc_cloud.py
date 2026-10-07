"""PubMed Central articles from the PMC Cloud Service (AWS Open Data bucket
`pmc-oa-opendata`, https://pmc-oa-opendata.s3.amazonaws.com/), no browser, no gate.

Layout: `PMC<id>.<version>/` holding `PMC<id>.<v>.json` (metadata: license_code,
is_pmc_openaccess, is_manuscript), `.xml` (JATS full text), `.txt`, and -- for the open-
access subset and CC-licensed author manuscripts -- `.pdf` plus every supplementary file
under its original name. Listing is the public S3 ListObjectsV2 call
(`?list-type=2&prefix=PMC<id>.`).

What this does and does not solve (investigated 2026-10):
* Full text: the bucket has JATS XML for *every* PMC article in the OA subset and for NIH
  author manuscripts (license_code "TDM" included), which replaces the PMC article HTML
  page that now often returns a reCAPTCHA "Checking your browser" interstitial.
* Supplements of NIHMS author manuscripts that are not CC-licensed (license_code TDM,
  is_pmc_openaccess false; e.g. PMC5894354, sci-RNA-seq Science 2017): the bucket holds
  only json/xml/txt, no supplements. The other routes are closed too:
  pmc.ncbi.nlm.nih.gov/articles/instance/<id>/bin/<file> returns a JavaScript "Preparing
  to download" page (then reCAPTCHA); europepmc.org/articles/PMC<id>/bin/<file> is
  behind Cloudflare (403 "Just a moment..."); Europe PMC REST /supplementaryFiles answers
  "not open access"; the NCBI OA web service (oa.fcgi) covers the OA subset only and now
  404s. So those stay `(manual)` rows -- a person downloads them in a browser (or the
  same file is often on the publisher's own static server, e.g. Springer ESM, which
  get_sources tries for Nature-family DOIs).
"""
from __future__ import annotations

import json
import re

from . import register
from .base import get, say

BUCKET = "https://pmc-oa-opendata.s3.amazonaws.com"
SKIP = re.compile(r"\.(jpe?g|gif|png|tiff?|txt)$", re.I)


def parse_listing(xml: bytes) -> list[tuple[str, int]]:
    """S3 ListObjectsV2 XML -> [(key, size)]."""
    keys = re.findall(rb"<Key>([^<]+)</Key>.*?<Size>(\d+)</Size>", xml, re.S)
    return [(k.decode(), int(s)) for k, s in keys]


def latest(objects: list[tuple[str, int]], pmcid: str) -> list[tuple[str, int]]:
    """Objects of the newest version folder PMC<id>.<v>/ only."""
    vers = {}
    for k, s in objects:
        m = re.match(rf"{pmcid}\.(\d+)/", k)
        if m:
            vers.setdefault(int(m.group(1)), []).append((k, s))
    return vers[max(vers)] if vers else []


def listing(pmcid: str) -> list[tuple[str, int]]:
    code, body = get(f"{BUCKET}/?list-type=2&prefix={pmcid}.")
    return latest(parse_listing(body), pmcid) if code == 200 else []


def fetch_article(folder, tag: str, pmcid: str) -> dict:
    """Save the JATS XML, PDF and supplements of `pmcid` from the bucket.
    -> {'xml': bool, 'pdf': bool, 'supp': [names], 'meta': {...}} (empty if absent)."""
    objs = listing(pmcid)
    got = {"xml": False, "pdf": False, "supp": [], "meta": {}}
    if not objs:
        return got
    say(f"    PMC Cloud: {len(objs)} object(s) for {pmcid}")
    for key, size in objs:
        base = key.split("/", 1)[1]
        url = f"{BUCKET}/{key}"
        if base.endswith(".json"):
            code, body = get(url)
            try:
                got["meta"] = json.loads(body) if code == 200 else {}
            except ValueError:
                pass
            continue
        if SKIP.search(base):
            continue
        if re.fullmatch(rf"{pmcid}\.\d+\.xml", base):
            got["xml"] = bool(folder.save(f"{tag}_{pmcid}.xml", url,
                                          "full-text JATS XML (PMC Cloud)", need=rb"<body"))
        elif re.fullmatch(rf"{pmcid}\.\d+\.pdf", base):
            got["pdf"] = bool(folder.save(f"{tag}_{pmcid}.pdf", url, "article PDF (PMC Cloud)"))
        else:
            if folder.save(f"{tag}_supp_{base}", url, "supplementary file (PMC Cloud)"):
                got["supp"].append(base)
    m = got["meta"]
    if m:
        folder.note(f"{pmcid} licence", f"{BUCKET}/{objs[0][0].split('/')[0]}/",
                    f"PMC licence {m.get('license_code')}; open access {m.get('is_pmc_openaccess')}; "
                    f"author manuscript {m.get('is_manuscript')}")
    return got


@register("pmccloud", doc="PMC Cloud (AWS pmc-oa-opendata): JATS XML, PDF, supplements of a PMCID")
def fetch(job, arg: str, what: str = "") -> None:
    pmcid = arg.upper() if arg.upper().startswith("PMC") else f"PMC{arg}"
    fetch_article(job.folder, pmcid, pmcid)
