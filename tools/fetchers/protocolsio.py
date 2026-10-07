"""protocols.io protocols (DOIs 10.17504/protocols.io.*, or protocols.io/view/... URLs).

Without credentials, protocols.io serves its own PDF export of a public protocol at
`https://www.protocols.io/view/<slug>.pdf`: resolve the DOI (doi.org redirects to the
view URL) and fetch that. The JSON API (`/api/v4/protocols/<id>`) needs a bearer token;
if PROTOCOLS_IO_TOKEN is set in the environment (a free "client access token" from
https://www.protocols.io/developers), the JSON -- steps, materials, attachments -- is
saved as well. Private links (protocols.io/private/...) cannot be fetched by script.

get_sources also scans the text of every fetched paper for protocols.io DOIs/links and
hands them to this fetcher (capped, see get_sources.MAX_LINKED_PROTOCOLS).
"""
from __future__ import annotations

import os
import re
import subprocess

from . import register
from .base import UA, get_json, safe_name, say

DOI_RE = re.compile(r"10\.17504/protocols\.io\.([a-z0-9]+)", re.I)
VIEW_RE = re.compile(r"https?://(?:www\.)?protocols\.io/view/([a-z0-9-]+?)(?:/v\d+)?(?:\.pdf|\.html)?(?=[\s\"'<>)\],;]|$)", re.I)


def find_refs(text: str) -> list[str]:
    """protocols.io DOIs and view URLs mentioned in a text, normalised, in order."""
    out = [f"10.17504/protocols.io.{m.group(1).lower().rstrip('.')}" for m in DOI_RE.finditer(text)]
    out += [f"https://www.protocols.io/view/{m.group(1).lower()}" for m in VIEW_RE.finditer(text)]
    return list(dict.fromkeys(out))


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
          doc="protocols.io: the public PDF export (and the JSON API with PROTOCOLS_IO_TOKEN)")
def fetch(job, arg: str, what: str = "") -> None:
    f = job.folder
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
    f.save(f"protocolsio_{slug}.pdf", view + ".pdf", desc)
    token = os.environ.get("PROTOCOLS_IO_TOKEN")
    if token:
        ident = short or slug
        data = get_json(f"https://www.protocols.io/api/v4/protocols/{ident}?content_format=markdown",
                        headers={"Authorization": f"Bearer {token}"})
        if data.get("payload"):
            import json
            f.keep(f"protocolsio_{slug}.json", json.dumps(data, indent=1).encode(),
                   f"https://www.protocols.io/api/v4/protocols/{ident}", "protocols.io API JSON")
