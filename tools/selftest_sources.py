#!/usr/bin/env python3
"""Offline self-test for the source fetcher (tools/get_sources.py, tools/fetchers/,
tools/doctext.py XML handling, catalogue/extra_sources.tsv). No network: every input is
synthetic and built here.

Run:  python3 tools/selftest_sources.py
"""
from __future__ import annotations

import contextlib
import copy
import gzip
import io
import os
import sys
import tempfile
import zipfile
import zlib
from types import SimpleNamespace
from unittest.mock import patch
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / "lib"))
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "catalogue" / "tools"))

from checks import Check  # noqa: E402

import doctext  # noqa: E402
import scrape_primers  # noqa: E402
import audit_protocol_sources as source_reviews  # noqa: E402
import fetchers  # noqa: E402
import get_sources as gs  # noqa: E402
from fetchers import base, pdf_url, pmc_cloud, protocolsio, scg_github, vendor  # noqa: E402
from source_oligo_check import (_exact_hit, check_graph, declared_locators,
                                unrepresented_input_ids,
                                unrepresented_inputs)  # noqa: E402
from selftest_source_oligos import graph_module_paths  # noqa: E402
from protocol_graph import MaterialKind, ProtocolGraph  # noqa: E402
from chemdraw import Scene, Segment  # noqa: E402
from storage import data_root, legacy_reference, sources_root  # noqa: E402

check = Check()
check.section("External source storage")
with patch.dict(os.environ):
    os.environ.pop("CHEM_DATA", None)
    check("fetcher defaults to the external article collection",
          base.data_dir(), Path("/husky/henriksson/NGS101_papers/sources"))
    check("legacy references use external storage without repository fallback",
          legacy_reference("pdf/paper.pdf"), data_root() / "legacy/pdf/paper.pdf")
with tempfile.TemporaryDirectory() as folder, patch.dict(os.environ, {"CHEM_DATA": folder}):
    check("an explicit storage override applies to fetchers and readers together",
          (base.data_dir(), sources_root(), legacy_reference("ref/paper.pdf")),
          (Path(folder) / "sources", Path(folder) / "sources", Path(folder) / "legacy/ref/paper.pdf"))
check.raises("a reference path cannot escape the selected external root",
             lambda: legacy_reference("../paper.pdf"))
check("a graph owned by a page builder enters the source transcription guard",
      [path.name for path in graph_module_paths(
          ROOT / "quartz-seq__10.1186+gb-2013-14-4-r31")],
      ["build_page.py"])

# ------------------------------------------------------------------------------
check.section("DOI routing")
check("bioRxiv dated DOI", gs.doi_kind("10.1101/2019.12.17.879304"), "biorxiv")
check("bioRxiv 6-digit DOI", gs.doi_kind("10.1101/003236"), "biorxiv")
check("bioRxiv DOI with a version suffix", gs.doi_kind("10.1101/2023.05.11.540245v2"), "biorxiv")
check("bioRxiv DOI with newer prefix", gs.doi_kind("10.64898/2026.07.07.736965"), "biorxiv")
check("bioRxiv newer DOI with a version suffix", gs.doi_kind("10.64898/2026.07.07.736965v1"), "biorxiv")
check("Genome Research (10.1101/gr.) is a journal", gs.doi_kind("10.1101/gr.110882.110"), "journal")
check("Genes & Dev (10.1101/gad.) is a journal", gs.doi_kind("10.1101/gad.350434.123"), "journal")
check("CSH Protocols (10.1101/pdb.) is a journal", gs.doi_kind("10.1101/pdb.prot5384"), "journal")
check("protocols.io", gs.doi_kind("10.17504/protocols.io.b4xwqxpe"), "protocolsio")
check("protocols.io PDF attachment discovery",
      protocolsio.attachment_urls(
          b"/URI (https://content.protocols.io/files/a1.xlsx) /URI "
          b"(https://content.protocols.io/files/a1.xlsx) "
          b"https://content.protocols.io/files/figure-2.pdf"),
      ["https://content.protocols.io/files/a1.xlsx",
       "https://content.protocols.io/files/figure-2.pdf"])
check("protocols.io legacy workbook link is not lost or truncated to xls",
      protocolsio.attachment_urls(b"https://content.protocols.io/files/order.xls "
                                 b"https://content.protocols.io/files/order.xlsm"),
      ["https://content.protocols.io/files/order.xls",
       "https://content.protocols.io/files/order.xlsm"])
check("protocols.io exported revision is read from the printed integer ID",
      protocolsio.protocol_integer_id("Protocol Integer ID:\n51484\nLatest version: 117877"),
      "51484")
check("protocols.io attachment metadata cannot silently select another revision",
      protocolsio.matching_public_protocol({"protocol": {"id": 117877, "public": 1}}, "51484"), {})
check("protocols.io private metadata is not followed",
      protocolsio.matching_public_protocol({"protocol": {"id": 51484, "public": 0}}, "51484"), {})
with tempfile.TemporaryDirectory() as folder:
    pdf_path = Path(folder) / "protocol.pdf"
    pdf_path.write_bytes(b"%PDF-1.4\n(no attachment annotations)")
    pdf_path.with_name(pdf_path.name + ".txt").write_text("Protocol Integer ID: 51484")
    calls = []
    class ProtocolFolder:
        def save(self, name, url, what):
            calls.append((name, url))
            return pdf_path if name.endswith(".pdf") else None
    workbook = "https://content.protocols.io/files/order.xls"
    with patch.object(protocolsio, "resolve", return_value="https://www.protocols.io/view/example"), \
         patch.object(protocolsio, "public_attachments", return_value=[workbook]) as public_metadata, \
         patch.dict("os.environ", {}, clear=True):
        protocolsio.fetch(SimpleNamespace(folder=ProtocolFolder()), "https://www.protocols.io/view/example")
    check("PDF without annotations still fetches the public original workbook",
          ([args.args[1] for args in public_metadata.call_args_list], calls[-1][1]),
          (["51484"], workbook))
# Explicit historical protocol revisions must not resolve to a newer cached export.
check("protocols.io explicit view revision survives PDF normalization",
      protocolsio.resolve("https://www.protocols.io/view/example/v3.pdf"),
      "https://www.protocols.io/view/example/v3")
check("protocols.io unversioned view keeps legacy export route",
      protocolsio.resolve("https://www.protocols.io/view/example.pdf"),
      "https://www.protocols.io/view/example")
with patch.object(protocolsio.subprocess, "run", return_value=SimpleNamespace(
        stdout="https://www.protocols.io/view/example/v4?redirect=true")):
    check("protocols.io versioned DOI cannot silently select redirected v4",
          protocolsio.resolve("10.17504/protocols.io.example/v3"),
          "https://www.protocols.io/view/example/v3")
check("protocols.io reference discovery preserves explicit DOI and view revisions",
      protocolsio.find_refs("10.17504/protocols.io.example/v3 and "
                            "https://www.protocols.io/view/example/v3.pdf)"),
      ["10.17504/protocols.io.example/v3", "https://www.protocols.io/view/example/v3"])
check("protocols.io printed revision does not use a latest-version link",
      protocolsio.protocol_version("Version 3\nLatest version: 4\nProtocol Integer ID: 123"), "3")
for observed in ("3", "4", ""):
    with tempfile.TemporaryDirectory() as folder:
        f = base.Folder(Path(folder), convert=False)
        current = Path(folder) / "protocolsio_example.pdf"
        current.write_bytes(b"%PDF-current-v4")
        calls = []
        def version_export(name, url, what):
            calls.append((name, url))
            pdf = Path(folder) / name
            pdf.write_bytes(b"%PDF-historical-response")
            pdf.with_name(pdf.name + ".txt").write_text(
                ("Version " + observed + "\n" if observed else "") +
                "Protocol Integer ID: 123\nLatest version: 999")
            return pdf
        with patch.object(f, "save", side_effect=version_export), \
             patch.object(protocolsio, "public_attachments", return_value=[]) as metadata, \
             patch.dict("os.environ", {}, clear=True):
            protocolsio.fetch(SimpleNamespace(folder=f), "https://www.protocols.io/view/example/v3")
        check("protocols.io v3 cache and requested PDF are distinct from current v4: " + (observed or "unlabelled"),
              (calls[0], current.read_bytes()),
              (("protocolsio_example_v3.pdf", "https://www.protocols.io/view/example/v3.pdf"), b"%PDF-current-v4"))
        if observed == "3":
            check("protocols.io metadata follows the verified historical PDF integer ID",
                  [call.args[1] for call in metadata.call_args_list], ["123"])
        else:
            check("protocols.io wrong or unverified revision cannot supply attachments: " + (observed or "unlabelled"),
                  (metadata.call_count, (Path(folder) / "protocolsio_example_v3.pdf").exists(),
                   bool(list((Path(folder) / "_invalid").glob("*.pdf"))),
                   any("requested revision v3" in row[4] for row in f.rows.values())),
                  (0, False, True, True))
            check("mandatory directory scraper excludes the actual rejected revision and text twin: " + (observed or "unlabelled"),
                  scrape_primers.files([folder]), [current])
            rejected_twin = next((Path(folder) / "_invalid").glob("*.pdf.txt"))
            check("explicit rejected-document inspection remains possible: " + (observed or "unlabelled"),
                  scrape_primers.files([str(rejected_twin)]), [rejected_twin])
with tempfile.TemporaryDirectory() as folder:
    f = base.Folder(Path(folder), convert=False)
    with patch.object(f, "save", return_value=None), \
         patch.object(protocolsio, "public_attachments") as metadata:
        protocolsio.fetch(SimpleNamespace(folder=f), "https://www.protocols.io/view/example/v3")
    check("protocols.io failed or duplicate historical download is not a verified original",
          (metadata.call_count, any("not substituted" in row[4] for row in f.rows.values())), (0, True))
for observed in ("3", "4"):
    with tempfile.TemporaryDirectory() as folder:
        f = base.Folder(Path(folder), convert=False)
        current = Path(folder) / "protocolsio_example.pdf"
        current.write_bytes(b"%PDF-existing-export")
        current.with_name(current.name + ".txt").write_text(
            "Version " + observed + "\nProtocol Integer ID: 123")
        digest = base.md5(current)
        f.by_md5[digest] = current.name
        f.rows["(duplicate) protocolsio_example_v3.pdf"] = [
            "(duplicate) protocolsio_example_v3.pdf", "https://www.protocols.io/view/example/v3.pdf", "", digest, "identical cached original"]
        with patch.object(f, "save", return_value=None), \
             patch.object(protocolsio, "public_attachments", return_value=[]) as metadata, \
             patch.dict("os.environ", {}, clear=True):
            protocolsio.fetch(SimpleNamespace(folder=f), "https://www.protocols.io/view/example/v3")
        check("protocols.io duplicate PDF must verify printed historical identity: v" + observed,
              (metadata.call_count, current.exists(), (Path(folder) / "_invalid").exists()),
              (1 if observed == "3" else 0, True, False))
        if observed == "3":
            check("protocols.io verified duplicate records original filename and revision",
                  any("Version 3 verified" in row[4] and current.name in row[4] for row in f.rows.values()))
for returned_id in ("123", "999"):
    with tempfile.TemporaryDirectory() as folder:
        f = base.Folder(Path(folder), convert=False)
        pdf = Path(folder) / "protocolsio_example_v3.pdf"
        pdf.write_bytes(b"%PDF-historical-export")
        pdf.with_name(pdf.name + ".txt").write_text("Version 3\nProtocol Integer ID: 123")
        with patch.object(f, "save", return_value=pdf), \
             patch.object(f, "keep") as keep, \
             patch.object(protocolsio, "public_attachments", return_value=[]), \
             patch.object(protocolsio, "get_json", return_value={"payload": {"id": returned_id}}) as api, \
             patch.dict("os.environ", {"PROTOCOLS_IO_TOKEN": "offline-test-token"}, clear=True):
            protocolsio.fetch(SimpleNamespace(folder=f), "https://www.protocols.io/view/example/v3")
        check("protocols.io authenticated metadata must match historical printed ID: " + returned_id,
              (api.call_args.args[0], keep.call_count),
              ("https://www.protocols.io/api/v4/protocols/123?content_format=markdown", 1 if returned_id == "123" else 0))
check("Research Square", gs.doi_kind("10.21203/rs.3.rs-3210240/v1"), "researchsquare")
check("Nature journal", gs.doi_kind("10.1038/ncomms14049"), "journal")
check("doi.org prefix stripped", gs.norm_doi("https://doi.org/10.1038/nmeth.1470"), "10.1038/nmeth.1470")
check("eLife DOIs go to the eLife fetcher",
      [h.name for h in fetchers.handlers_for_doi("10.7554/eLife.73971")], ["elife"])
check("protocols.io DOIs go to the protocols.io fetcher",
      [h.name for h in fetchers.handlers_for_doi("10.17504/protocols.io.nrkdd4w")], ["protocolsio"])
check("Crossref relation parsing",
      gs.related_dois({"has-preprint": [{"id-type": "doi", "id": "10.1101/126268"}]}, "has-preprint"),
      ["10.1101/126268"])

# ------------------------------------------------------------------------------
check.section("download validation")
RECAPTCHA = (b'<!doctype html><html lang="en-US" dir="ltr"><head><base href="https://www.google.com/'
             b'recaptcha/challengepage/"><script>window.ppConfig={productName:"RecaptchaChallengePageUi"}'
             b'</script></head><body>' + b" " * 20000 + b"</body></html>")
PREPARING = (b"<html><head><meta name=\"viewport\"><title>Preparing to download ...</title></head>"
             b"<body>Preparing to download</body></html>")
CLOUDFLARE = b"<!DOCTYPE html><html><head><title>Just a moment...</title></head><body></body></html>"
ARTICLE = (b"<html><head><title>A paper</title></head><body>" +
           b"<p>We describe a single-cell method. " * 300 + b"</p></body></html>")
PDF = b"%PDF-1.4\n1 0 obj<<>>endobj\ntrailer<<>>\n%%EOF\n"


def zbytes(entries: dict[str, bytes]) -> bytes:
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        for k, v in entries.items():
            z.writestr(k, v)
    return buf.getvalue()


XLSX = zbytes({"xl/workbook.xml": b"<workbook/>", "[Content_Types].xml": b"<Types/>"})

with tempfile.TemporaryDirectory() as tmp:
    t = Path(tmp)

    def v(name: str, data: bytes, **kw) -> str:
        p = t / name
        p.write_bytes(data)
        return base.validate(p, name, **kw)

    check("a real PDF passes", v("a.pdf", PDF), "")
    check("a reCAPTCHA page saved as .pdf is rejected", "reCAPTCHA" in v("b.pdf", RECAPTCHA), True)
    check("a truncated PDF (no %%EOF) is rejected", "truncated" in v("c.pdf", PDF[:-7]), True)
    check("PMC 'Preparing to download' saved as .xlsx is rejected",
          "Preparing to download" in v("d.xlsx", PREPARING), True)
    check("a reCAPTCHA page saved as .docx is rejected", v("e.docx", RECAPTCHA) != "", True)
    check("a real xlsx (zip) passes", v("f.xlsx", XLSX), "")
    check("a truncated zip (cut before the central directory) is rejected",
          "truncated or corrupt zip" in v("g.zip", zbytes({"a.txt": b"x" * 5000})[:-40]), True)
    check("a zip of zero-byte placeholders is rejected", v("h.zip", zbytes({"x.gif": b""})),
          "zip holds no files")
    check("a reCAPTCHA page saved as .html is rejected", "reCAPTCHA" in v("i.html", RECAPTCHA), True)
    check("a Cloudflare interstitial is rejected", v("j.html", CLOUDFLARE) != "", True)
    check("an HTML page with article text passes", v("k.html", ARTICLE), "")
    check("an HTML page with almost no text is rejected",
          "no article content" in v("l.html", b"<html><body>Hi</body></html>"), True)
    topic_body = b"Actual documented index pooling guidance. " * 5
    topic_start = b'<html data-mc-runtime-file-type="Topic;Default"><body><div id="mc-main-content">'
    short_topic = topic_start + b"<h1>Pooling</h1><div><p>" + topic_body + b"</p></div></div></body></html>"
    check("a short complete documentation topic is accepted",
          v("short-topic.htm", short_topic), "")
    check("topic structure cannot stand in for real main text",
          bool(v("empty-topic.htm", topic_start + b'</div><nav>' + topic_body + b'</nav></body></html>')), True)
    check("main-content marker alone is not a documentation topic",
          bool(v("untyped-topic.htm", short_topic.replace(b'Topic;Default',b'WebHelp2'))), True)
    check("an unclosed short chapter is rejected",
          bool(v("cut-topic.htm", topic_start + topic_body)), True)
    check("short documentation cannot bypass a challenge page",
          bool(v("gated-topic.htm", RECAPTCHA + short_topic)), True)
    check("short documentation still needs requested source markers",
          bool(v("wrong-topic.htm", short_topic, need=b"different-source-marker")), True)
    check("Europe PMC errorBean XML is rejected",
          v("m.xml", b'<?xml version="1.0"?><ns4:errorBean><errCode>404</errCode></ns4:errorBean>') != "", True)
    check("JATS without <body> fails when full text is required",
          "front matter" in v("n.xml", b"<article><front>t</front></article>", need=rb"<body"), True)
    check("JATS with a body passes", v("o.xml", b"<article><body><p>x</p></body></article>",
                                       need=rb"<body"), "")
    check("legacy .xls must be OLE2", v("p.xls", b"<html><head></head><body>x</body></html>") != "", True)
    check("an OLE2 .xls passes", v("q.xls", bytes.fromhex("D0CF11E0A1B11AE1") + b"\0" * 600), "")
    check("a complete .gz passes", v("r.txt.gz", gzip.compress(b"ACGT\n" * 1000)), "")
    check("a truncated .gz is rejected", "gzip" in v("s.txt.gz", gzip.compress(b"ACGT\n" * 1000)[:-12]), True)
    check("an HTML page saved as a barcode .txt is rejected", v("t.txt", ARTICLE) != "", True)
    check("a plain barcode list passes", v("u.txt", b"AAACCCAAGAAACACT\nAAACCCAAGAAACCAT\n"), "")
    _big = t / "cut.mp4"
    with _big.open("wb") as _fh:
        _fh.truncate(50_000_000)
    check("a file of exactly the old 50 MB cap is truncated", "truncated" in base.validate(_big), True)
    _big.unlink()
    check("unknown extension: only the gate check applies", "reCAPTCHA" in v("w.bad", RECAPTCHA), True)

    # Folder: invalid download becomes (manual); identical files are not stored twice
    _quiet = contextlib.ExitStack()
    _quiet.enter_context(contextlib.redirect_stdout(io.StringIO()))
    _quiet.enter_context(contextlib.redirect_stderr(io.StringIO()))
    f = base.Folder(t / "proto", convert=False)
    f.keep("a.pdf", PDF, "https://x/a.pdf", "paper")
    f.keep("b.pdf", PDF, "https://y/b.pdf", "same paper elsewhere")
    f.keep("c.xlsx", PREPARING, "https://x/c.xlsx", "supplement")
    check("keep: a valid file gets a row with its md5", f.rows["a.pdf"][3], base.md5(t / "proto" / "a.pdf"))
    check("keep: an identical second copy is a (duplicate) row, not a file",
          ((t / "proto" / "b.pdf").exists(), "(duplicate) b.pdf" in f.rows), (False, True))
    check("keep: a gate page becomes (manual) with the reason",
          ((t / "proto" / "c.xlsx").exists(), "Preparing to download" in f.rows["(manual) c.xlsx"][4]),
          (False, True))
    check("manifest is written as rows arrive",
          (t / "proto" / "MANIFEST.tsv").read_text().count("\n"), 4)
    check("same source URL reuses its validated cached original",
          f.save("a.pdf", "https://x/a.pdf", "paper"), t / "proto" / "a.pdf")
    check("different source URL cannot reuse a cached original by filename",
          f.save("a.pdf", "https://x/another.pdf", "another paper"), None)
    check("different source bytes cannot overwrite a registered original",
          f.keep("a.pdf", PDF + b"\n", "https://x/third.pdf", "third paper"), None)
    check("source collisions retain original file and provenance",
          ((t / "proto" / "a.pdf").read_bytes(), f.rows["a.pdf"][1]),
          (PDF, "https://x/a.pdf"))
    check("source collision records requested URL and reason",
          any(row[1] == "https://x/another.pdf" and "cache filename a.pdf" in row[4]
              for row in f.rows.values()), True)
    check("duplicate cache row also retains its source identity",
          f.save("b.pdf", "https://x/different.pdf", "different paper"), None)
    f.manual("hub.txt", "https://one/hub.txt", "first hub", "HTTP 404")
    f.manual("hub.txt", "https://two/hub.txt", "second hub", "HTTP 503")
    f.manual("hub.txt", "https://two/hub.txt", "second hub", "HTTP 404")
    failures = [row for row in f.rows.values() if row[0].startswith("(manual) hub")]
    check("failed URLs sharing a basename keep separate retryable manifest entries",
          sorted((row[1],row[4].split("; ")[-1]) for row in failures),
          [("https://one/hub.txt","HTTP 404"),("https://two/hub.txt","HTTP 404")])
    f.keep("hub.txt", b"hub actual_source\n", "https://two/hub.txt", "second hub")
    check("source recovery clears its own failure without erasing another original's failure",
          [(row[1],row[4].split("; ")[-1]) for row in f.rows.values()
           if row[0].startswith("(manual) hub")],
          [("https://one/hub.txt","HTTP 404")])

    # zip expansion pulls documents out into tracked top-level files
    (t / "proto" / "z.zip").write_bytes(zbytes({"S1/Table_S1.xlsx": XLSX, "fig.png": b"\x89PNG",
                                                "__MACOSX/._x": b"junk"}))
    f.expand_zip("z.zip", "z_", "https://x/z.zip", "supp")
    check("expand_zip keeps documents, skips images and resource forks",
          sorted(p.name for p in (t / "proto").glob("z_*")), ["z_Table_S1.xlsx"])

    # revalidate demotes a gate page that an old run saved as a supplement
    (t / "proto" / "old_supp.pdf").write_bytes(RECAPTCHA)
    (t / "proto" / "old_supp.pdf.txt").write_text("x")
    f.rows["old_supp.pdf"] = ["old_supp.pdf", "https://pmc/x.pdf", "1", "0", "PMC supplementary file"]
    f.write()
    gs.revalidate([t / "proto"])
    _quiet.close()
    g = base.Folder(t / "proto", convert=False)
    check("revalidate moves a gate page and its twin to _invalid/",
          sorted(p.name for p in (t / "proto" / "_invalid").iterdir()), ["old_supp.pdf", "old_supp.pdf.txt"])
    check("...and turns its row into (manual) keeping the URL",
          (g.rows["(manual) old_supp.pdf"][1], "old_supp.pdf" in g.rows), ("https://pmc/x.pdf", False))

    # revalidate upgrades the old ambiguous paper row without mistaking a supplement
    # for the article itself.
    legacy = t / "legacy"
    legacy.mkdir()
    tag = "s41587-023-01685-z"
    (legacy / f"{tag}_supp_table.xlsx").write_bytes(XLSX)
    (legacy / "MANIFEST.tsv").write_text(
        base.Folder.HEADER
        + f"{tag}_supp_table.xlsx\thttps://x/table.xlsx\t1\t0\tsupplementary file\n"
        + f"(manual) {tag}\thttps://doi.org/10.1038/{tag}\t\t\t"
          "no open copy found: publisher page only\n")
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        gs.revalidate([legacy])
    migrated = base.Folder(legacy, convert=False)
    check("revalidate renames an old catch-all manual row to a full-text row",
          (f"(manual) {tag}" in migrated.rows, f"(manual) {tag} full text" in migrated.rows),
          (False, True))
    check("a fetched supplement does not count as fetched article full text",
          migrated.rows[f"(manual) {tag} full text"][1],
          f"https://doi.org/10.1038/{tag}")

# ------------------------------------------------------------------------------
check.section("extra sources table")
SAMPLE = ("# comment\nslug\tsource\twhat\twhy\n"
          "a__x\thttps://cdn.example.com/x/CG000331_Guide_RevE.pdf\t10x guide\tvendor kit\n"
          "a__x\tscg:data/BD/*\tBD guides\t\n"
          "b__y\tdoi:10.1038/s41598-017-16546-4\tjournal version\n")
rows = fetchers.parse_extra_sources(SAMPLE)
check("rows parse, comments skipped, short rows padded",
      [(r["slug"], r["source"], r["why"]) for r in rows],
      [("a__x", "https://cdn.example.com/x/CG000331_Guide_RevE.pdf", "vendor kit"),
       ("a__x", "scg:data/BD/*", ""), ("b__y", "doi:10.1038/s41598-017-16546-4", "")])
try:
    fetchers.parse_extra_sources("slug\tsource\twhat\twhy\nonly-a-slug\n")
    _err = ""
except ValueError as e:
    _err = str(e)
check("a row without a source is an error", "needs slug and source" in _err, True)
check("spec split: handler form", fetchers.split_spec("scg:data/BD/*"), ("scg", "data/BD/*"))
check("spec split: a URL is not a handler name", fetchers.split_spec("https://x.org/a.pdf"),
      ("", "https://x.org/a.pdf"))
_route_calls = []
_route_job = fetchers.Job("x", [], None)
_route_job.paper = lambda doi: _route_calls.append(("doi", doi))
_route_job.pmc = lambda pmcid: _route_calls.append(("pmc", pmcid))
fetchers.run_spec(_route_job, "pmc:PMC8238917")
fetchers.run_spec(_route_job, "doi:10.1038/example")
check("PMC and DOI extra-source routes stay distinct", _route_calls,
      [("pmc", "PMC8238917"), ("doi", "10.1038/example")])
check("URL routing: protocols.io URL -> protocolsio",
      fetchers.handler_for_url("https://www.protocols.io/view/hydrop-atac-v1-0-bxsbpnan").name, "protocolsio")
check("URL routing: raw GitHub of scg_lib_structs -> scg",
      fetchers.handler_for_url("https://raw.githubusercontent.com/Teichlab/scg_lib_structs/master/data/BD/BD_CLS1.txt").name,
      "scg")
check("URL routing: anything else -> url", fetchers.handler_for_url("https://cdn.10xgenomics.com/a.pdf").name, "url")
check("vendor file name from a URL with %20 and a query",
      vendor.vendor_prefix("https://cdn.10xgenomics.com/x") +
      vendor.url_filename("https://cdn.10xgenomics.com/raw/upload/v1/in-line%20documents/Dual_Index_Kit_TT_Set_A.csv?x=1"),
      "10x_Dual_Index_Kit_TT_Set_A.csv")
check("vendor file name respects an explicit export format",
      vendor.url_filename("https://docs.google.com/spreadsheets/d/example/export?format=xlsx"),
      "export.xlsx")
check("explicit ZIP export uses archive validation and extraction",
      vendor.url_filename("https://www.ebi.ac.uk/europepmc/webservices/rest/PMC4691826/supplementaryFiles?format=zip"),
      "supplementaryFiles.zip")
check("unknown export formats do not bypass HTML validation",
      vendor.url_filename("https://example.org/export?format=unknown"), "export.html")
check("generic documentation chapter filenames retain their project identity",
      [vendor.url_filename(f"https://support-docs.illumina.com/SHARE/{project}/Content/{name}")
       for project,name in (("IndexAdaptersPooling","Overview.htm"),
                            ("AdapterSequences","Overview.htm"),
                            ("AdapterSequences","RevisionHistory.htm"))],
      ["IndexAdaptersPooling-Overview.htm", "AdapterSequences-Overview.htm",
       "AdapterSequences-RevisionHistory.htm"])
check("descriptive oligo-topic filenames remain stable",
      vendor.url_filename("https://support-docs.illumina.com/SHARE/AdapterSequences/Content/TruSeq-SmallRNA.htm"),
      "TruSeq-SmallRNA.htm")
check("Addgene sequence listings retain distinct plasmid IDs",
      [vendor.url_filename(f"https://www.addgene.org/{plasmid}/sequences/")
       for plasmid in (52963, 196709)],
      ["addgene-52963-sequences.html", "addgene-196709-sequences.html"])
check("original experiment XML export keeps its type",
      vendor.url_filename("https://www.ebi.ac.uk/ena/browser/api/xml/SRX2134246?format=xml"),
      "SRX2134246.xml")
check("ENA XML endpoint keeps its type without an export query",
      vendor.url_filename("https://www.ebi.ac.uk/ena/browser/api/xml/SRX2134246"),
      "SRX2134246.xml")
check("GEO records with one CGI path retain distinct accessions and scope",
      [vendor.url_filename("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?" + query)
       for query in ("acc=GSE60456", "acc=GSE43087",
                     "view=full&form=text&targ=all&acc=GSE60456",
                     "acc=GSE60456&targ=self&form=text&view=brief")],
      ["geo-GSE60456-self-full.html", "geo-GSE43087-self-full.html",
       "geo-GSE60456-all-full.txt", "geo-GSE60456-self-brief.txt"])
import catalogue as cat  # noqa: E402
_slugs = cat.source_slugs()
_real = fetchers.extra_sources()
check("catalogue/extra_sources.tsv parses", len(_real) > 0, True)
check("every extra_sources slug is a catalogue slug", [r["slug"] for r in _real if r["slug"] not in _slugs], [])
check("every extra_sources handler name is registered",
      [r["source"] for r in _real
       if fetchers.split_spec(r["source"])[0] not in ("", "doi", "pmc", *fetchers.load())], [])

# ------------------------------------------------------------------------------
check.section("scg_lib_structs repo mapping (fixture tree)")
TREE = [
    "README.md", "chemistries/adapters/Chromium.fa", "chemistries/adapters/STRT_seq.fa",
    "chemistries/adapters/SureCell.fa", "chemistries/whitelists/3M-february-2018.csv.gz",
    "data/10X-Genomics/CG000108_AssayConfiguration_SC3v2.pdf", "data/10X-Genomics/737K-august-2016.txt.gz",
    "data/10X-Genomics/CG000183_ChromiumSingleCell3__v3_UG_Rev-A.pdf",
    "data/10X-Genomics/CG00026_Chromium_Single_Cell_3__Reagent_Kits_User_Guide_RevB.PDF",
    "data/BD/BD_CLS1.txt", "data/BD/GMX_BD-Rhapsody-WTA-alpha-Protocol_UG_EN.pdf",
    "data/STRT-seq_family/STRT_GenomeRes_2011_SI.pdf", "data/STRT-seq_family/STRT_bc.fa",
    "data/HyDrop/elife-73971-supp1-v4.docx", "data/HyDrop/HyDrop_20200130_plate-1-96_barcode.csv",
    "data/PIP-seq/fb_v3_bc1.tsv", "data/illumina-adapter-sequences-1000000002694-14.pdf",
    "data/itChIP-seq_Table3.xlsx", "docs/source/ge/10xChromium3v1.md", "docs/source/ge/10xChromium3v2.md",
    "docs/source/ge/10xChromium3v3.md", "docs/source/ge/STRT-seq.md", "docs/source/ge/sci-RNA-seq3.md",
    "methods_html/10xChromium3.html", "methods_html/10xChromium3v1.html", "methods_html/PIP-seq.html",
    "methods_html/PIP-seq_v1p.html", "methods_html/STRT-seq_family.html", "methods_html/BD_Rhapsody.html",
    "methods_html/HyDrop_RNA.html", "methods_html/itChIP-seq.html", "methods_html/SureCell.html",
    "methods_html/figs/x.png", "data/BD/bead.png",
]
PAGE_10X = ('<a href="../data/10X-Genomics/CG000108_AssayConfiguration_SC3v2.pdf">v2</a>'
            '<a href="../data/10X-Genomics/737K-august-2016.txt.gz">wl</a>'
            '<a href="../data/10X-Genomics/CG000183_ChromiumSingleCell3__v3_UG_Rev-A.pdf">v3</a>'
            '<a href="https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium3v1.html">v1</a>'
            '<a href="../style_related/page_format.css">css</a>')
check("page links: ../data files and methods_html pages",
      scg_github.page_links(PAGE_10X),
      ["data/10X-Genomics/CG000108_AssayConfiguration_SC3v2.pdf", "data/10X-Genomics/737K-august-2016.txt.gz",
       "data/10X-Genomics/CG000183_ChromiumSingleCell3__v3_UG_Rev-A.pdf", "methods_html/10xChromium3v1.html"])
check("10xChromium3: linked guides + whitelist + v2/v3 md, not the v1 page or v1 md",
      scg_github.scg_paths(TREE, "10xChromium3", scg_github.page_links(PAGE_10X)),
      ["data/10X-Genomics/CG000108_AssayConfiguration_SC3v2.pdf", "data/10X-Genomics/737K-august-2016.txt.gz",
       "data/10X-Genomics/CG000183_ChromiumSingleCell3__v3_UG_Rev-A.pdf", "docs/source/ge/10xChromium3v2.md",
       "docs/source/ge/10xChromium3v3.md", "chemistries/adapters/Chromium.fa"])
check("BD_Rhapsody: linked files only, never images",
      scg_github.scg_paths(TREE, "BD_Rhapsody", ["data/BD/BD_CLS1.txt", "data/BD/bead.png",
                                                 "data/BD/GMX_BD-Rhapsody-WTA-alpha-Protocol_UG_EN.pdf"]),
      ["data/BD/BD_CLS1.txt", "data/BD/GMX_BD-Rhapsody-WTA-alpha-Protocol_UG_EN.pdf"])
check("STRT-seq_family: the unlinked data/STRT-seq_family/ dir, its md and adapters, by name",
      scg_github.scg_paths(TREE, "STRT-seq_family", []),
      ["data/STRT-seq_family/STRT_GenomeRes_2011_SI.pdf", "data/STRT-seq_family/STRT_bc.fa",
       "docs/source/ge/STRT-seq.md", "chemistries/adapters/STRT_seq.fa"])
check("HyDrop_RNA: data/HyDrop/ by name prefix",
      scg_github.scg_paths(TREE, "HyDrop_RNA", []),
      ["data/HyDrop/HyDrop_20200130_plate-1-96_barcode.csv", "data/HyDrop/elife-73971-supp1-v4.docx"])
check("PIP-seq: its own sub-page PIP-seq_v1p is fetched",
      scg_github.scg_paths(TREE, "PIP-seq", ["methods_html/PIP-seq_v1p.html", "methods_html/inDrop.html"]),
      ["methods_html/PIP-seq_v1p.html", "data/PIP-seq/fb_v3_bc1.tsv"])
check("SureCell: adapters by name plus the aliased Illumina adapter PDF",
      scg_github.scg_paths(TREE, "SureCell", []),
      ["chemistries/adapters/SureCell.fa", "data/illumina-adapter-sequences-1000000002694-14.pdf"])
check("local names", [scg_github.local_name(p) for p in
                      ("data/BD/BD_CLS1.txt", "docs/source/ge/STRT-seq.md", "chemistries/adapters/STRT_seq.fa",
                       "methods_html/PIP-seq_v1p.html")],
      ["scg_BD_CLS1.txt", "scg_docs_ge_STRT-seq.md", "scg_adapters_STRT_seq.fa", "upstream_PIP-seq_v1p.html"])

# ------------------------------------------------------------------------------
check.section("PMC Cloud listing, protocols.io references")
LISTING = (b'<?xml version="1.0"?><ListBucketResult><Contents><Key>PMC123.1/PMC123.1.xml</Key>'
           b'<Size>10</Size></Contents><Contents><Key>PMC123.2/PMC123.2.xml</Key><Size>11</Size></Contents>'
           b'<Contents><Key>PMC123.2/mmc1.pdf</Key><Size>99</Size></Contents></ListBucketResult>')
check("S3 listing parsed and reduced to the newest version",
      pmc_cloud.latest(pmc_cloud.parse_listing(LISTING), "PMC123"),
      [("PMC123.2/PMC123.2.xml", 11), ("PMC123.2/mmc1.pdf", 99)])
check("protocols.io DOIs and view links found in text",
      protocolsio.find_refs("see dx.doi.org/10.17504/protocols.io.be8mjhu6. and "
                            "https://www.protocols.io/view/hydrop-atac-v1-0-bxsbpnan)"),
      ["10.17504/protocols.io.be8mjhu6", "https://www.protocols.io/view/hydrop-atac-v1-0-bxsbpnan"])

# ------------------------------------------------------------------------------
check.section("doctext: JATS XML")
JATS = ('<?xml version="1.0"?><!DOCTYPE article><article><front><article-meta><title-group>'
        '<article-title>A <italic>method</italic></article-title></title-group><contrib><name>'
        '<surname>Li</surname><given-names>Yun</given-names></name></contrib></article-meta></front>'
        '<body><sec><title>Methods</title><p>RT primer <bold>AAGCAGTGG</bold>TATCAACGCAGAGT &amp; 37 &#x000b0;C.</p>'
        '<table-wrap><table><tr><td>P5</td><td>AATGATACGG</td></tr></table></table-wrap></sec></body></article>')
_lines = [ln for ln in doctext.jats_text(JATS).splitlines() if ln.strip()]
check("titles and paragraphs on their own lines; inline markup does not split a sequence",
      _lines, ["A method", "Li Yun", "Methods", "RT primer AAGCAGTGGTATCAACGCAGAGT & 37 °C.", "P5\tAATGATACGG"])
check(".xml is converted", ".xml" in doctext.CONVERT, True)

check.section("opaque original PDF URLs")
_pdf_a = "https://example.org/download?artifactId=first"
_pdf_b = "https://example.org/download?artifactId=second"
check("opaque PDF query identity survives cache naming",
      (pdf_url.filename(_pdf_a) != pdf_url.filename(_pdf_b),
       pdf_url.filename(_pdf_a).endswith(".pdf")), (True, True))
_pdf_calls = []
pdf_url.fetch(SimpleNamespace(folder=SimpleNamespace(save=lambda *args: _pdf_calls.append(args))),
              _pdf_a, "original exhibit")
check("explicit PDF route preserves the original URL and validation suffix",
      _pdf_calls, [(pdf_url.filename(_pdf_a), _pdf_a, "original exhibit")])
check.raises("explicit PDF route rejects a non-HTTP source", lambda: pdf_url.filename("file:///tmp/source.pdf"))

check.section("doctext: optional EPS conversion")
with tempfile.TemporaryDirectory() as eps_folder:
    eps = Path(eps_folder) / "source figure.eps"
    eps.write_bytes(b"%!PS-Adobe-3.0 EPSF-3.0\n%%BoundingBox: 0 0 40 40\nshowpage\n")
    with patch.object(doctext, "GHOSTSCRIPT", None), contextlib.redirect_stderr(io.StringIO()):
        check("EPS without optional converter skips without a fake text twin",
              (doctext.convert(eps), doctext.text_path(eps).exists()), (None, False))
    with patch.object(doctext, "GHOSTSCRIPT", "gs"), \
            patch.object(doctext.subprocess, "run", return_value=SimpleNamespace(returncode=1)), \
            contextlib.redirect_stderr(io.StringIO()):
        check("failed EPS conversion does not claim a readable original", doctext.read(eps), None)
    def fake_eps_conversion(args, **kwargs):
        output = next(arg.split("=", 1)[1] for arg in args if arg.startswith("-sOutputFile="))
        Path(output).write_bytes(b"%PDF-1.4\nsynthetic converted figure")
        return SimpleNamespace(returncode=0)
    with patch.object(doctext, "GHOSTSCRIPT", "gs"), \
            patch.object(doctext.subprocess, "run", side_effect=fake_eps_conversion) as eps_run, \
            patch.object(doctext, "_read", return_value="Source figure oligo ACGT\n") as pdf_read:
        check("EPS text passes through the ordinary PDF reader", doctext.read(eps),
              "Source figure oligo ACGT\n")
        check("EPS conversion uses safer mode and absolute source path",
              ("-dSAFER" in eps_run.call_args.args[0], eps_run.call_args.args[0][-1],
               pdf_read.call_args.args[0].suffix), (True, str(eps.resolve()), ".pdf"))
    check("directory conversion discovers EPS originals", doctext.documents([eps_folder]), [eps])

check.section("doctext: optional RTF conversion")
with tempfile.TemporaryDirectory() as rtf_folder:
    rtf = Path(rtf_folder) / "source protocol.rtf"
    rtf.write_bytes(b"{\\rtf1\\ansi Original protocol.}")
    with patch.object(doctext, "SOFFICE", None), contextlib.redirect_stderr(io.StringIO()):
        check("RTF without optional converter skips without a fake text twin",
              (doctext.convert(rtf), doctext.text_path(rtf).exists()), (None, False))
    def fake_rtf_conversion(args, **kwargs):
        output_dir = Path(args[args.index("--outdir") + 1])
        (output_dir / "source protocol.docx").write_bytes(b"converted document fixture")
        return SimpleNamespace(returncode=0)
    with patch.object(doctext, "SOFFICE", "soffice"), \
            patch.object(doctext.subprocess, "run", side_effect=fake_rtf_conversion) as rtf_run, \
            patch.object(doctext, "_read", return_value="Original protocol.\n") as docx_read:
        check("RTF uses ordinary DOCX extraction after conversion", doctext.read(rtf),
              "Original protocol.\n")
        check("RTF conversion passes a distinct profile and absolute source path",
              (any(arg.startswith("-env:UserInstallation=file://") for arg in rtf_run.call_args.args[0]),
               rtf_run.call_args.args[0][-1], docx_read.call_args.args[0].suffix),
              (True, str(rtf.resolve()), ".docx"))

check.section("doctext: filled XFA PDF form")
def xfa_stream(number: int, xml: bytes) -> bytes:
    packed = zlib.compress(xml)
    return (f"{number} 0 obj\n<</Filter/FlateDecode/Length {len(packed)}>>\nstream\n".encode()
            + packed + b"\nendstream\nendobj\n")

blank_xfa = b'<xfa:datasets xmlns:xfa="http://www.xfa.org/schema/xfa-data/1.0/"><xfa:data><form><sampleprep/></form></xfa:data></xfa:datasets>'
filled_xfa = b'<xfa:datasets xmlns:xfa="http://www.xfa.org/schema/xfa-data/1.0/"><xfa:data><form><sampleprep>30 min tagmentation; 40 mM EDTA</sampleprep><index>ACGT</index></form></xfa:data></xfa:datasets>'
sample_xfa_pdf = b"%PDF-1.7\n" + xfa_stream(2, blank_xfa) + xfa_stream(3, filled_xfa) + b"%%EOF\n"
check("filled XFA dataset chosen over blank template", doctext.xfa_text(sample_xfa_pdf),
      "sampleprep: 30 min tagmentation; 40 mM EDTA\nindex: ACGT\n")
check("ordinary PDF without XFA has no form text", doctext.xfa_text(PDF), None)
template_xfa = (b'<template xmlns="http://www.xfa.org/schema/xfa-template/2.8/">'
                b'<draw><value><text>Consumable: an adapter</text></value></draw>'
                b'<draw><value><text>/5Phos/ACGTACGTACGT</text></value></draw>'
                b'<draw><value><exData contentType="text/html"><body><p>Supplier reference</p></body></exData></value></draw>'
                b'<event><script>not_source_text()</script></event>'
                b'<draw><value><image>encoded_image</image></value></draw></template>')
template_pdf = b"%PDF-1.7\n" + xfa_stream(2, template_xfa) + xfa_stream(3, blank_xfa) + b"%%EOF\n"
template_text = ("## XFA template text (form layout and scripts not evaluated)\n"
                 "Consumable: an adapter\n/5Phos/ACGTACGTACGT\nSupplier reference\n")
check("static form template survives an empty XFA dataset; scripts and images excluded",
      doctext.xfa_text(template_pdf), template_text)
check("filled dataset takes precedence over static form labels",
      doctext.xfa_text(template_pdf + xfa_stream(4, filled_xfa)),
      "sampleprep: 30 min tagmentation; 40 mM EDTA\nindex: ACGT\n")
with tempfile.TemporaryDirectory() as folder:
    source = Path(folder) / "dynamic.pdf"
    source.write_bytes(template_pdf)
    with patch.object(doctext, "_read", return_value="Please wait...\n"):
        check("PDF reader replaces placeholder with labelled template content",
              doctext.read(source), template_text)

check.section("source oligo transcription")
check("protocols.io source locator retains the accession before /v1",
      "rm7vzqr92vx1" in declared_locators(
          "hybrid-bag-seq__10.1186+s13059-025-03790-5"), True)
check("source match accepts a printed exact oligo",
      _exact_hit("primer 5′-ACGTACGTACGTACGT-3′", "ACGTACGTACGTACGT"), True)
check("source match rejects a one-base copy error",
      _exact_hit("primer 5′-ACGTACGTACGTACGA-3′", "ACGTACGTACGTACGT"), False)
with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    (folder / "paper.pdf.txt").write_text(
        "first 5′-ACGTACGTAACCGGTT-3′\nsecond 5′-GATTACAGGCTAGCTA-3′\n")
    scene = Scene()
    scene.strand("first", [Segment("left", "ACGTACGT"), Segment("right", "AACCGGTT")])
    scene.strand("second", [Segment("complete", "GATTACAGGCTAGCTA")], rev=True)
    graph = ProtocolGraph("scene pool")
    graph.input("adapters", "two independent ordering oligos", MaterialKind.ADAPTER,
                state=scene)
    check("scene checks separate whole strands in ordering orientation",
          check_graph(graph, folder), ([], [], 2))
    (folder / "paper.pdf.txt").write_text(
        "first 5′-ACGTACGTAACCGGTT-3′\nsecond 5′-GATTACAGGCTAGCTG-3′\n")
    check("a scene strand transcription error is detected",
          (len(check_graph(graph, folder)[0]), check_graph(graph, folder)[2]), (1, 2))
    scene.strand("unassigned", [Segment("unknown", "[unassigned strand]", placeholder=True)])
    check("partly represented scene still reports its incomplete input",
          unrepresented_input_ids(graph), ("adapters",))
    check("unknown scene strand does not hide the known whole strands",
          check_graph(graph, folder)[2], 2)
with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    (folder / "paper.pdf.txt").write_text("primer 5′-ACGTACGTACGTACGT-3′")
    graph = ProtocolGraph("example")
    graph.input("primer", "PCR primer", MaterialKind.OLIGO,
                state="ACGTACGTACGTACGT")
    check("declared graph oligo is checked against local text",
          check_graph(graph, folder), ([], [], 1))
    (folder / "paper.pdf.txt").write_text("primer 5′-ACGTACGTACGTACGA-3′")
    failures, skips, count = check_graph(graph, folder)
    check("copy error fails the source check", (len(failures), skips, count), (1, [], 1))
    check("reasoned material exception opts out", check_graph(graph, folder, {"primer": "source image"}),
          ([], [], 0))
    check("stale exception fails", len(check_graph(graph, folder, {"other": "old"})[0]), 2)
    check("missing text skips", check_graph(graph, folder / "absent"),
          ([], [f"example: no extracted source text in {folder / 'absent'}"], 0))
    graph.input("unnamed_sequence", "an unsequenced oligo", MaterialKind.OLIGO)
    check("unrepresented graph input is reported", unrepresented_inputs(graph), 1)
    check("coverage details name the missing input",
          unrepresented_input_ids(graph), ("unnamed_sequence",))

check.section("independent source-review acceptance gate")
_review_pass = {"reviewer": "reader-one", "date": "2026-10-15",
                "locations": "Methods paragraph 2", "finding": "Primer list checked"}
_second_pass = {**_review_pass, "reviewer": "reader-two",
                "locations": "Methods paragraph 2 and supplement table"}
_review = {
    "dir": "example", "status": "reviewed",
    "source_checklist": [{"id": "paper", "role": "defining", "locator": "DOI example",
                          "availability": "read", "relevance": "relevant",
                          "pass_1": _review_pass, "pass_2": _second_pass}],
    "linked_information": {"result": "none_found", "pass_1": _review_pass,
                           "pass_2": _second_pass},
    "missing_claims": [], "conclusion": "No missing source detail claimed",
}
_known = {"example": {"dir": "example"}}
check("two independent source and link checks pass",
      source_reviews.validate("example", _review, _known), [])
for _pending_field in ("pending_sources", "pending_linked_targets"):
    _pending_review = copy.deepcopy(_review)
    _pending_review[_pending_field] = [{"locator": "named original awaiting inspection"}]
    check("completed review cannot bypass " + _pending_field,
          bool(source_reviews.validate("example", _pending_review, _known)), True)
    _pending_review["status"] = "in_progress"
    check("partial review can retain explicit " + _pending_field,
          source_reviews.validate("example", _pending_review, _known), [])
_vendor_review = copy.deepcopy(_review)
_vendor_review["source_checklist"][0]["role"] = "vendor"
check("vendor guide defines a DOI-free commercial protocol",
      source_reviews.validate("example", _vendor_review, _known), [])
check("vendor guide alone cannot replace a DOI article",
      bool(source_reviews.validate("example", _vendor_review,
                                   {"example": {"doi": "10.1/example"}})), True)
_one_pass = copy.deepcopy(_review)
del _one_pass["source_checklist"][0]["pass_2"]
check("one source check cannot complete review",
      bool(source_reviews.validate("example", _one_pass, _known)), True)
_same_reviewer = copy.deepcopy(_review)
_same_reviewer["linked_information"]["pass_2"]["reviewer"] = "reader-one"
check("same reader cannot supply both link checks",
      bool(source_reviews.validate("example", _same_reviewer, _known)), True)
_unavailable = copy.deepcopy(_review)
_unavailable["source_checklist"][0].update(availability="unavailable", reason="paywall")
check("unavailable relevant source blocks completed review",
      bool(source_reviews.validate("example", _unavailable, _known)), True)
_available = copy.deepcopy(_review)
_available["source_checklist"][0].update(availability="available", pass_1=None, pass_2=None)
check("downloaded but unread source cannot complete review",
      bool(source_reviews.validate("example", _available, _known)), True)
_available["status"] = "in_progress"
check("downloaded source can await reading without an access blocker",
      source_reviews.validate("example", _available, _known), [])
_excluded_once = copy.deepcopy(_review)
_excluded_once["source_checklist"].append({
    "id": "supplement", "role": "supplement", "locator": "figure-only supplement",
    "availability": "read", "relevance": "not_relevant", "reason": "no methods",
    "pass_1": _review_pass})
check("excluded source also needs two independent checks",
      bool(source_reviews.validate("example", _excluded_once, _known)), True)
_unlisted = copy.deepcopy(_review)
_unlisted["linked_information"]["result"] = "enumerated"
check("linked material must appear in source checklist",
      bool(source_reviews.validate("example", _unlisted, _known)), True)
_blocked_claim = copy.deepcopy(_review)
_blocked_claim["missing_claims"] = [{"claim": "index primer", "disposition": "blocked",
                                     "evidence": "linked table inaccessible"}]
check("blocked claim prevents completed review",
      bool(source_reviews.validate("example", _blocked_claim, _known)), True)
_confirmed_claim = copy.deepcopy(_review)
_confirmed_claim["missing_claims"] = [{"claim": "index primer",
                                       "disposition": "confirmed_missing",
                                       "evidence": "two source checks"}]
check("confirmed absence requires accepted_missing status",
      bool(source_reviews.validate("example", _confirmed_claim, _known)), True)
_known_protocols = source_reviews.catalogue()
_review_errors = []
for _path in source_reviews.REVIEWS.glob("*.json"):
    try:
        _record = source_reviews.json.loads(_path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as _exc:
        _review_errors.append(f"{_path.name}: {_exc}")
        continue
    _review_errors.extend(source_reviews.validate(_path.stem, _record, _known_protocols))
check("checked-in source reviews satisfy the acceptance gate", _review_errors, [])
_accepted_path = source_reviews.REVIEWS / "10x-chromium-single-cell-3-ge-v1__10.1038+ncomms14049.json"
_accepted = source_reviews.json.loads(_accepted_path.read_text(encoding="utf-8"))
_unrendered = copy.deepcopy(_accepted)
_unrendered["source_checklist"].append({**copy.deepcopy(_accepted["source_checklist"][0]),
                                         "id": "unrendered_checked_source"})
check("accepted absence page must enumerate every checked source",
      any("public page omits source unrendered_checked_source" in issue for issue in
          source_reviews.validate(_accepted_path.stem, _unrendered, _known_protocols)), True)

check.report()
