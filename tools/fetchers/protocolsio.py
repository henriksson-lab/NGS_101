"""protocols.io protocols (DOIs 10.17504/protocols.io.*, or protocols.io/view/... URLs).

Without credentials, protocols.io serves its own PDF export of a public protocol at
`https://www.protocols.io/view/<slug>.pdf`: resolve the DOI (doi.org redirects to the
view URL) and fetch that. Explicit `/vN` references keep that revision in the request
and cache filename; the PDF must print the requested revision before attachments
are followed. A byte-identical cached PDF is reused only after checking its digest
and printed revision, with the historical URL retained in the manifest. Its printed integer ID also identifies an anonymous public
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
from .base import UA, get_json, md5, safe_name, say

DOI_RE = re.compile(r"10\.17504/protocols\.io\.([a-z0-9]+)(/v\d+)?", re.I)
VIEW_RE = re.compile(r"https?://(?:www\.)?protocols\.io/view/([a-z0-9-]+)(/v\d+)?(?:\.pdf|\.html)?(?=[\s\"'<>)\],;]|$)", re.I)
ATTACHMENT_RE = re.compile(
    rb"https://content\.protocols\.io/files/[A-Za-z0-9._~-]+\.(?:pdf|xlsx|xlsm|xls|docx|doc|pptx|ppt|zip|csv|tsv|txt|png|jpeg|jpg)(?![A-Za-z0-9])",
    re.I,
)


def find_refs(text: str) -> list[str]:
    """protocols.io DOIs and view URLs mentioned in a text, normalised, in order."""
    out = [f"10.17504/protocols.io.{m.group(1).lower().rstrip('.')}{m.group(2) or ''}" for m in DOI_RE.finditer(text)]
    out += [f"https://www.protocols.io/view/{m.group(1).lower()}{m.group(2) or ''}" for m in VIEW_RE.finditer(text)]
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


def requested_version(arg: str) -> str:
    match = re.search(r"/v(\d+)(?=[/.?#]|$)", arg, re.I)
    return match.group(1) if match else ""


def protocol_version(text: str) -> str:
    """Read the printed revision label, not a latest-version link in the PDF."""
    match = re.search(r"^\s*Version\s+(\d+)\s*$", text, re.I | re.M)
    return match.group(1) if match else ""


def _normalized_view(url: str, version: str = "") -> str:
    url = url.split("?")[0].split("#")[0].rstrip("/")
    url = re.sub(r"(?:\.pdf|\.html)$", "", url)
    url = re.sub(r"/v\d+$", "", url).rstrip("/")
    return url + ("/v" + version if version else "")


def resolve(arg: str) -> str:
    """DOI/URL to a view URL; retain an explicitly requested historical revision."""
    version = requested_version(arg)
    if arg.startswith("http") and "/view/" in arg:
        return _normalized_view(arg, version)
    doi = arg.removeprefix("https://doi.org/").removeprefix("https://dx.doi.org/")
    done = subprocess.run(["curl", "-sSIL", "--max-time", "60", "-A", UA, "-o", "/dev/null",
                           "-w", "%{url_effective}", f"https://doi.org/{doi}"],
                          capture_output=True, text=True, check=False)
    url = done.stdout.strip()
    if "/view/" not in url:
        return ""
    return _normalized_view(url, version)


def _duplicate_pdf(folder, name: str):
    """An explicit download may equal a legacy cached file; verify its digest."""
    row = getattr(folder, "rows", {}).get("(duplicate) " + name)
    if not row or not row[3]:
        return None
    existing = getattr(folder, "by_md5", {}).get(row[3])
    if not existing:
        return None
    path = folder.path / existing
    return path if path.is_file() and md5(path) == row[3] else None


def _reject_revision(folder, pdf, url: str, desc: str, requested: str, observed: str) -> None:
    """Retain a misidentified response outside the accepted original cache."""
    invalid = pdf.parent / "_invalid"
    invalid.mkdir(exist_ok=True)
    target = invalid / pdf.name
    count = 1
    while target.exists():
        target = invalid / f"{pdf.stem}_{count}{pdf.suffix}"
        count += 1
    twin = pdf.with_name(pdf.name + ".txt")
    pdf.rename(target)
    if twin.exists():
        twin.rename(target.with_name(target.name + ".txt"))
    if hasattr(folder, "rows"):
        folder.rows.pop(pdf.name, None)
        for digest, name in list(folder.by_md5.items()):
            if name == pdf.name:
                folder.by_md5.pop(digest)
    folder.manual(pdf.name, url, desc,
                  f"requested revision v{requested}, but printed Version is {observed or 'unverified'}; "
                  f"response retained at {target.relative_to(pdf.parent)}; attachments not followed")


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
    version = requested_version(arg)
    slug = re.sub(r"/v(\d+)$", r"_v\1", view.split("/view/", 1)[-1])
    say(f"  protocols.io {slug}")
    desc = what or f"protocols.io protocol {arg} (PDF export)"
    cache_name = f"protocolsio_{slug}.pdf"
    pdf = f.save(cache_name, view + ".pdf", desc)
    duplicate = False
    if pdf is None and version:
        pdf = _duplicate_pdf(f, cache_name)
        duplicate = pdf is not None
    if pdf:
        urls = attachment_urls(pdf.read_bytes())
        twin = pdf.with_name(pdf.name + ".txt")
        if twin.exists():
            text = twin.read_text()
        else:
            # --no-text still needs the exported ID to find missing attachment links.
            import doctext
            text = doctext.read(pdf) or ""
        if version and protocol_version(text) != version:
            if duplicate:
                f.manual(f"protocolsio_{slug}_revision_identity", view + ".pdf", desc,
                         f"requested revision v{version}, but identical cached payload prints "
                         f"Version {protocol_version(text) or 'unverified'}; original cache retained; attachments not followed")
            else:
                _reject_revision(f, pdf, view + ".pdf", desc, version, protocol_version(text))
            return
        ident = protocol_integer_id(text)
        if version and hasattr(f, "note"):
            f._clear_manual(f"protocolsio_{slug}_revision_identity", view + ".pdf")
            f.note(cache_name + " revision", view + ".pdf",
                   f"printed Version {version} verified; exported integer ID {ident or 'not printed'}; local original {pdf.name}")
            f.write()
        if ident:
            urls.extend(public_attachments(f, ident))
        for i, url in enumerate(dict.fromkeys(urls), 1):
            ext = url.rsplit(".", 1)[-1].lower()
            ident = url.rsplit("/", 1)[-1].rsplit(".", 1)[0]
            f.save(f"protocolsio_{slug}_attachment_{i}_{ident}.{ext}", url,
                   f"{desc}; public attachment linked from the PDF export")
    elif version:
        f.manual(f"protocolsio_{slug}_revision_identity", view + ".pdf", desc,
                 "historical revision could not be independently verified from a returned PDF; "
                 "an unversioned or duplicate cache file is not substituted")
        return
    token = os.environ.get("PROTOCOLS_IO_TOKEN")
    if token:
        ident = protocol_integer_id(text) if version else (short or slug)
        if not ident:
            return
        data = get_json(f"https://www.protocols.io/api/v4/protocols/{ident}?content_format=markdown",
                        headers={"Authorization": f"Bearer {token}"})
        if version and (not isinstance(data.get("payload"), dict)
                        or str(data["payload"].get("id")) != ident):
            f.manual(f"protocolsio_{slug}_authenticated_identity",
                     f"https://www.protocols.io/api/v4/protocols/{ident}",
                     "protocols.io authenticated revision metadata",
                     "API payload did not identify the verified PDF's protocol integer ID; not retained")
            return
        if data.get("payload"):
            f.keep(f"protocolsio_{slug}.json", json.dumps(data, indent=1).encode(),
                   f"https://www.protocols.io/api/v4/protocols/{ident}", "protocols.io API JSON")
