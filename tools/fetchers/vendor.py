"""Generic URL fetcher: vendor documents (10x Genomics CDN, BD, Illumina support PDFs, a
raw GitHub file, ...) declared by URL in catalogue/extra_sources.tsv.

The file name normally uses the URL's last path segment. Accession-based CGI URLs
retain their record identity and scope, and an explicit export format supplies its
suffix. The download is validated
like every other (a vendor's bot wall or 404 page becomes a `(manual)` row, not a .pdf).
Vendor documents stay the vendor's copyright: they are fetched into $CHEM_DATA for
reading, never into the repository.
"""
from __future__ import annotations

import urllib.parse
import re

from . import register


def url_filename(url: str) -> str:
    parsed = urllib.parse.urlparse(url)
    path = parsed.path
    # Addgene has many distinct /<plasmid>/sequences/ pages. Preserve the ID so
    # two linked plasmids do not silently share one cached sequences.html file.
    parts = path.strip("/").split("/")
    if parsed.netloc.lower() in {"addgene.org", "www.addgene.org"} and len(parts) == 2 and parts[0].isdigit() and parts[1] == "sequences":
        return f"addgene-{parts[0]}-sequences.html"
    # GEO identifies the original in its query, not the shared acc.cgi path.
    # Keep the requested record scope too: series-only and family metadata differ.
    if parsed.netloc.lower() in {"ncbi.nlm.nih.gov", "www.ncbi.nlm.nih.gov"} and path == "/geo/query/acc.cgi":
        query = dict(urllib.parse.parse_qsl(parsed.query))
        accession = query.get("acc", "").upper()
        if re.fullmatch(r"G(?:SE|SM|PL|DS)\d+", accession):
            target = re.sub(r"[^a-z0-9_-]", "_", query.get("targ", "self").lower())
            view = re.sub(r"[^a-z0-9_-]", "_", query.get("view", "full").lower())
            form = query.get("form", "html").lower()
            suffix = "txt" if form == "text" else "html"
            return f"geo-{accession}-{target}-{view}.{suffix}"
    name = urllib.parse.unquote(path.rstrip("/").rsplit("/", 1)[-1]) or "index.html"
    # Documentation projects reuse these generic chapter names. Preserve their
    # project identity while keeping established descriptive oligo filenames stable.
    if (parsed.netloc.lower() == "support-docs.illumina.com"
            and len(parts) >= 3 and parts[0] == "SHARE"
            and name.lower() in {"overview.htm", "revisionhistory.htm"}):
        return parts[1] + "-" + name
    if parsed.netloc.lower() in {"ebi.ac.uk", "www.ebi.ac.uk"} and path.startswith("/ena/browser/api/xml/"):
        return name if name.lower().endswith(".xml") else name + ".xml"
    if "." not in name:
        requested = dict(urllib.parse.parse_qsl(parsed.query)).get("format", "").lower()
        suffix = requested if requested in {"csv", "docx", "pdf", "tsv", "xlsx", "xml", "zip"} else "html"
        name += "." + suffix
    return name


def vendor_prefix(url: str) -> str:
    host = urllib.parse.urlparse(url).netloc.lower()
    for key, tag in (("10xgenomics", "10x"), ("bdbiosciences", "bd"), ("bd.com", "bd"),
                     ("illumina", "illumina"), ("bio-rad", "biorad"), ("githubusercontent", "gh")):
        if key in host:
            return tag + "_"
    return ""


@register("url", doc="any URL: download as its last path segment, validated by type")
def fetch(job, arg: str, what: str = "") -> None:
    name = vendor_prefix(arg) + url_filename(arg)
    job.folder.save(name, arg, what or "extra source (catalogue/extra_sources.tsv)")

