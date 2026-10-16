"""protocols.io protocols (DOIs 10.17504/protocols.io.*, or protocols.io/view/... URLs).

Without credentials, protocols.io serves its own PDF export of a public protocol at
`https://www.protocols.io/view/<slug>.pdf`: resolve the DOI (doi.org redirects to the
view URL) and fetch that. Its printed integer ID also identifies an anonymous public
page API response, which can expose attachments omitted from the PDF. Fetch only the
exported revision's public metadata. The v4 JSON API needs a bearer token;
if PROTOCOLS_IO_TOKEN is set in the environment (a free "client access token" from
https://www.protocols.io/developers), the JSON -- steps, materials, attachments -- is
saved as well. Private links (protocols.io/private/...) cannot be fetched by script.

get_sources also scans the text of every fetched paper for protocols.io DOIs/links and
hands them to this fetcher (capped, see get_sources.MAX_LINKED_PROTOCOLS).
"""
from __future__ import annotations

import json
import os
import re
import subprocess
from urllib.parse import urlencode

from . import register
from .base import UA, get_json, safe_name, say

DOI_RE = re.compile(r"10\.17504/protocols\.io\.([a-z0-9]+)", re.I)
VIEW_RE = re.compile(r"https?://(?:www\.)?protocols\.io/view/([a-z0-9-]+?)(?:/v\d+)?(?:\.pdf|\.html)?(?=[\s\"'<>)\],;]|$)", re.I)
ATTACHMENT_RE = re.compile(
    rb"https://content\.protocols\.io/files/[A-Za-z0-9._~-]+\.(?:pdf|xlsx|xlsm|xls|docx|doc|pptx|ppt|zip|csv|tsv|txt|png|jpeg|jpg)(?![A-Za-z0-9])",
    re.I,
)


def find_refs(text: str) -> list[str]:
    """protocols.io DOIs and view URLs mentioned in a text, normalised, in order."""
    out = [f"10.17504/protocols.io.{m.group(1).lower().rstrip('.')}" for m in DOI_RE.finditer(text)]
    out += [f"https://www.protocols.io/view/{m.group(1).lower()}" for m in VIEW_RE.finditer(text)]
    return list(dict.fromkeys(out))


def attachment_urls(pdf: bytes) -> list[str]:
    """Public source attachments linked by a protocols.io PDF export.

    The anonymous PDF is sometimes the only public index of its supplementary oligo
    workbook: the attachment URL lives in a PDF link annotation and is not repeated in
    extracted text.  Preserve order and collapse the several identical annotations that
    a single attachment commonly receives.
    """
    return list(dict.fromkeys(m.group(0).decode("ascii")
                              for m in ATTACHMENT_RE.finditer(pdf)))


def protocol_integer_id(text: str) -> str:
    """The exported revision's ID; do not substitute the latest public revision."""
    match = re.search(r"Protocol\s+Integer\s+ID\s*:\s*(\d+)", text, re.I)
    return match.group(1) if match else ""


def public_api_url(ident: str) -> str:
    fields = ("guidelines", "documents", "description", "materials", "steps")
    return (f"https://www.protocols.io/api/v1/protocols/{ident}?"
            + urlencode([("fields[]", field) for field in fields]))


def matching_public_protocol(data: dict, ident: str) -> dict:
    """Only use attachment metadata for the requested, public protocol revision."""
    if not isinstance(data, dict):
        return {}
    protocol = data.get("protocol", {})
    if (not isinstance(protocol, dict) or str(protocol.get("id")) != ident
            or protocol.get("public") not in (True, 1, "1")):
        return {}
    return protocol


def public_attachments(folder, ident: str) -> list[str]:
    """Recover files omitted from PDF link annotations via the anonymous page API."""
    url = public_api_url(ident)
    path = folder.save(f"protocolsio_{ident}_public.json", url,
                       "protocols.io anonymous public metadata for the exported revision")
    if not path:
        return []
    data = json.loads(path.read_text())
    protocol = matching_public_protocol(data, ident)
    if not protocol:
        folder.manual(f"protocolsio_{ident}_attachments", url, "protocols.io public attachments",
                      "public API did not identify the requested public revision; attachments not followed")
        return []
    return attachment_urls(json.dumps(protocol, ensure_ascii=False).encode())


def resolve(arg: str) -> str:
    """DOI or URL -> the protocol's view URL (no version suffix), or ''."""
    if arg.startswith("http") and "/view/" in arg:
        return re.sub(r"(/v\d+)?(\.pdf|\.html)?$", "", arg.rstrip("/"))
    doi = arg.removeprefix("https://doi.org/").removeprefix("https://dx.doi.org/")
    done = subprocess.run(["curl", "-sSIL", "--max-time", "60", "-A", UA, "-o", "/dev/null",
                           "-w", "%{url_effective}", f"https://doi.org/{doi}"],
                          capture_output=True, text=True, check=False)
    url = done.stdout.strip()
    if "/view/" not in url:
        return ""
    return re.sub(r"(/v\d+)?/?$", "", url.split("?")[0])


@register("protocolsio", doi=r"^10\.17504/protocols\.io\.", url=r"^https?://(?:www\.)?protocols\.io/",
          doc="protocols.io: public PDF, revision-specific public attachments, optional authenticated v4 JSON")
def fetch(job, arg: str, what: str = "") -> None:
    f = job.folder
    api = re.search(r"/api/v1/protocols/(\d+)(?:[/?]|$)", arg)
    if api:
        for url in public_attachments(f, api.group(1)):
            f.save("protocolsio_attachment_" + url.rsplit("/", 1)[-1], url,
                   what or "public attachment identified by protocols.io metadata")
        return
    if "/private/" in arg:
        f.manual("protocolsio_" + safe_name(arg.rsplit("/", 1)[-1]), arg, what or "protocols.io",
                 "private protocols.io link: needs a browser session")
        return
    view = resolve(arg)
    m = DOI_RE.search(arg)
    short = m.group(1).lower() if m else ""
    if not view:
        f.manual(f"protocolsio_{short or safe_name(arg)}.pdf", arg if arg.startswith("http")
                 else f"https://doi.org/{arg}", what or "protocols.io protocol",
                 "DOI did not resolve to a protocols.io page")
        return
    slug = view.rsplit("/", 1)[-1]
    say(f"  protocols.io {slug}")
    desc = what or f"protocols.io protocol {arg} (PDF export)"
    pdf = f.save(f"protocolsio_{slug}.pdf", view + ".pdf", desc)
    if pdf:
        urls = attachment_urls(pdf.read_bytes())
        twin = pdf.with_name(pdf.name + ".txt")
        if twin.exists():
            text = twin.read_text()
        else:
            # --no-text still needs the exported ID to find missing attachment links.
            import doctext
            text = doctext.read(pdf) or ""
        ident = protocol_integer_id(text)
        if ident:
            urls.extend(public_attachments(f, ident))
        for i, url in enumerate(dict.fromkeys(urls), 1):
            ext = url.rsplit(".", 1)[-1].lower()
            ident = url.rsplit("/", 1)[-1].rsplit(".", 1)[0]
            f.save(f"protocolsio_{slug}_attachment_{i}_{ident}.{ext}", url,
                   f"{desc}; public attachment linked from the PDF export")
    token = os.environ.get("PROTOCOLS_IO_TOKEN")
    if token:
        ident = short or slug
        data = get_json(f"https://www.protocols.io/api/v4/protocols/{ident}?content_format=markdown",
                        headers={"Authorization": f"Bearer {token}"})
        if data.get("payload"):
            f.keep(f"protocolsio_{slug}.json", json.dumps(data, indent=1).encode(),
                   f"https://www.protocols.io/api/v4/protocols/{ident}", "protocols.io API JSON")
