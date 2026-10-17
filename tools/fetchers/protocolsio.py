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


PRINTED_DOI_RE = re.compile(
    r"^\s*DOI:\s*https?://(?:dx\.)?doi\.org/(10\.17504/protocols\.io\.[a-z0-9]+(?:/v\d+)?)",
    re.I | re.M)


def printed_doi(text: str) -> str:
    """The export's own DOI (its 'DOI: ...' header line), not a DOI it cites."""
    match = PRINTED_DOI_RE.search(text)
    return match.group(1).lower() if match else ""


def printed_revision(text: str, accession: str = "") -> str:
    """The export's own revision: its printed Version label or, when it prints none, the /vN
    of its own printed DOI -- only when that DOI names `accession`, so a cited or
    latest-version DOI of another protocol cannot supply the revision."""
    label = protocol_version(text)
    doi = printed_doi(text)
    doi_version = requested_version(doi)
    if label and doi_version and label != doi_version:
        return ""  # contradictory identity is not a verified revision
    if label or not accession or re.sub(r"/v\d+$", "", doi) != f"10.17504/protocols.io.{accession.lower()}":
        return label
    return requested_version(doi)


def cited_export(cited: str, printed: str, version: str) -> bool:
    """Does an export printing DOI `printed` (and Version label `version`) identify the
    cited DOI? The accession must match; a cited /vN must match the printed /vN or Version
    label -- a printed different revision is not the cited one. An unversioned citation
    names no revision, so the export of that accession answers it."""
    split = lambda d: (re.sub(r"/v\d+$", "", d.lower()), requested_version(d))  # noqa: E731
    (c_acc, c_ver), (p_acc, p_ver) = split(cited), split(printed)
    if p_ver and version and p_ver != version:
        return False
    return bool(p_acc) and c_acc == p_acc and (not c_ver or c_ver == (p_ver or version))


def cached_export(folder, cited: str):
    """A protocols.io PDF already in the folder whose printed identity is the cited DOI."""
    for name in sorted(getattr(folder, "rows", {})):
        if not (name.startswith("protocolsio_") and name.endswith(".pdf")):
            continue
        pdf = folder.path / name
        twin = pdf.with_name(name + ".txt")
        if not (pdf.is_file() and twin.is_file()):
            continue
        row = folder.rows[name]
        if not row[3] or md5(pdf) != row[3]:
            continue
        text = twin.read_text(errors="replace")
        if cited_export(cited, printed_doi(text), protocol_version(text)):
            return pdf
    return None


def _answered(folder, short: str, version: str, doi_url: str):
    """Record that the folder's cached export answers a cited DOI; clear its manual row."""
    cited = f"10.17504/protocols.io.{short}{version}"
    cached = cached_export(folder, cited)
    if not cached:
        return None
    text = cached.with_name(cached.name + ".txt").read_text(errors="replace")
    folder._clear_manual(f"protocolsio_{short}.pdf", doi_url)
    folder.note(f"protocols.io {short}{version}", doi_url,
                f"cited DOI answered by {cached.name}, which prints DOI {printed_doi(text)}"
                + (f", Version {protocol_version(text)}" if protocol_version(text) else ""))
    folder.write()
    return cached


def settle_manual_rows(folder) -> int:
    """Clear '(manual) protocolsio_<accession>.pdf' rows of earlier runs whose DOI is now
    answered by an export on disk (fetched under its view-URL name). -> rows cleared."""
    cleared = 0
    for key, row in list(getattr(folder, "rows", {}).items()):
        m = DOI_RE.search(row[1]) if key.startswith("(manual) protocolsio_") else None
        if m and key == f"(manual) protocolsio_{m.group(1).lower()}.pdf":
            cleared += bool(_answered(folder, m.group(1).lower(), m.group(2) or "", row[1]))
    return cleared


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


def _reject_revision(folder, pdf, url: str, desc: str, requested: str, observed: str,
                     *, identity_problem: str = "") -> None:
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
    problem = identity_problem or (
        f"requested revision v{requested}, but printed Version is {observed or 'unverified'}")
    folder.manual(pdf.name, url, desc,
                  problem + "; "
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
    m = DOI_RE.search(arg)
    short = m.group(1).lower() if m else ""
    doi_url = arg if arg.startswith("http") else f"https://doi.org/{arg}"
    cached = _answered(f, short, m.group(2) or "", doi_url) if m else None
    if cached:      # the DOI's own export is already here (fetched under its view-URL slug)
        say(f"  protocols.io {short}: already here as {cached.name}")
        return
    view = resolve(arg)
    if not view:
        f.manual(f"protocolsio_{short or safe_name(arg)}.pdf", doi_url, what or "protocols.io protocol",
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
        accession = short or re.sub(r"/v\d+$", "", view.split("/view/", 1)[-1]).rsplit("-", 1)[-1]
        observed = printed_revision(text, accession)
        export_doi = printed_doi(text)
        identity_problem = ""
        label = protocol_version(text)
        doi_version = requested_version(export_doi)
        if label and doi_version and label != doi_version:
            identity_problem = (f"export prints contradictory DOI {export_doi} "
                                f"and Version {label}")
        elif short and export_doi and not cited_export(
                f"10.17504/protocols.io.{short}" + (f"/v{version}" if version else ""),
                export_doi, protocol_version(text)):
            identity_problem = (f"requested {arg}, but export prints DOI {export_doi} "
                                f"and Version {protocol_version(text) or 'unlabelled'}")
        if identity_problem or (version and observed != version):
            if duplicate:
                f.manual(f"protocolsio_{slug}_revision_identity", view + ".pdf", desc,
                         (identity_problem or f"requested revision v{version}, but identical cached payload prints "
                          f"Version {observed or 'unverified'}")
                         + "; original cache retained; attachments not followed")
            else:
                _reject_revision(f, pdf, view + ".pdf", desc, version, observed,
                                 identity_problem=identity_problem)
            return
        ident = protocol_integer_id(text)
        if version and hasattr(f, "note"):
            f._clear_manual(f"protocolsio_{slug}_revision_identity", view + ".pdf")
            how = (f"printed Version {version}" if protocol_version(text)
                   else f"printed DOI {printed_doi(text)} (no Version label)")
            f.note(cache_name + " revision", view + ".pdf",
                   f"{how} verified; exported integer ID {ident or 'not printed'}; local original {pdf.name}")
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
