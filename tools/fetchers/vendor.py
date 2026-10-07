"""Generic URL fetcher: vendor documents (10x Genomics CDN, BD, Illumina support PDFs, a
raw GitHub file, ...) declared by URL in catalogue/extra_sources.tsv.

The file name is the URL's last path segment (query dropped). The download is validated
like every other (a vendor's bot wall or 404 page becomes a `(manual)` row, not a .pdf).
Vendor documents stay the vendor's copyright: they are fetched into $CHEM_DATA for
reading, never into the repository.
"""
from __future__ import annotations

import urllib.parse

from . import register


def url_filename(url: str) -> str:
    path = urllib.parse.urlparse(url).path
    name = urllib.parse.unquote(path.rstrip("/").rsplit("/", 1)[-1]) or "index.html"
    if "." not in name:
        name += ".html"
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

