"""eLife (DOIs 10.7554/eLife.N): the article PDF, JATS XML and every supplementary file,
from eLife's public API (api.elifesciences.org/articles/N), served by cdn.elifesciences.org.

PMC's copies of eLife supplements sit behind the PMC download gate; eLife's CDN does not.
"""
from __future__ import annotations

import re

from . import register
from .base import get_json, say


@register("elife", doi=r"^10\.7554/elife\.\d+",
          doc="eLife API: PDF, XML and supplementary files from cdn.elifesciences.org")
def fetch(job, arg: str, what: str = "") -> None:
    m = re.search(r"elife\.(\d+)", arg, re.I)
    if not m:
        return
    n = m.group(1)
    art = get_json(f"https://api.elifesciences.org/articles/{n}", retries=2)
    if not art:
        job.folder.manual(f"eLife.{n}_supplementary", f"https://elifesciences.org/articles/{n}",
                          "eLife article", "eLife API unavailable")
        return
    v = art.get("version", "")
    say(f"  eLife {n} v{v}")
    f = job.folder
    tag = f"eLife.{n}"
    if art.get("pdf"):
        f.save(f"{tag}_v{v}.pdf", art["pdf"], f"eLife article PDF, v{v}")
    if art.get("xml"):
        f.save(f"{tag}_v{v}.xml", art["xml"], f"eLife JATS XML, v{v}", need=rb"<body")
    files = list(art.get("additionalFiles") or [])
    # figure-supplement source data hang off the figures in the body
    def walk(o):
        if isinstance(o, dict):
            for s in o.get("sourceData") or []:
                files.append(s)
            for val in o.values():
                walk(val)
        elif isinstance(o, list):
            for val in o:
                walk(val)
    walk(art.get("body") or [])
    for x in files:
        uri, name = x.get("uri"), x.get("filename")
        if uri and name:
            f.save(f"{tag}_{name}", uri, f"eLife supplementary file {x.get('id', '')}: "
                   f"{(x.get('title') or x.get('label') or '')[:80]}")
