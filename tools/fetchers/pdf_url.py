"""Explicit PDF originals whose download URLs have no file extension."""
from hashlib import sha256
from urllib.parse import urlsplit

from . import register


def filename(url: str) -> str:
    if urlsplit(url).scheme not in {"http", "https"}:
        raise ValueError("an original PDF URL must use HTTP or HTTPS")
    # Query parameters may be the complete identity of an opaque export endpoint.
    return "source_pdf_" + sha256(url.encode()).hexdigest()[:20] + ".pdf"


@register("pdf", doc="explicit PDF URL, including opaque exports; validate PDF bytes and preserve full URL identity")
def fetch(job, arg: str, what: str = "") -> None:
    job.folder.save(filename(arg), arg, what or "explicit PDF original")
