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
check("Nature article and figure pages are named by article id with an .html suffix",
      [vendor.url_filename("https://www.nature.com/articles/" + a) for a in (
          "nprot.2013.038", "s41596-019-0202-2/figures/2", "s41586-020-1234-5.pdf")],
      ["nprot.2013.038.html", "s41596-019-0202-2-figure-2.html", "s41586-020-1234-5.pdf"])
check("protocols.io unlabelled export takes its revision only from its own printed DOI",
      [protocolsio.printed_revision(t, "abc123") for t in (
          "DOI: https://dx.doi.org/10.17504/protocols.io.abc123/v1",
          "DOI: https://dx.doi.org/10.17504/protocols.io.other9/v1",
          "Version 2\nDOI: https://dx.doi.org/10.17504/protocols.io.abc123/v1",
          "cites https://dx.doi.org/10.17504/protocols.io.abc123/v1")], ["1", "", "", ""])
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
    check("a publisher 'Client Challenge' interstitial is reported as a gate",
          base.gate_reason(b"<html><head><title>Client Challenge</title></head>"),
          "challenge/gate page (Client Challenge)")
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
    _offline = patch.object(base, "download", return_value=(0, "offline test: no network"))
    with _offline:
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
    with _offline:
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

    # An article saved under another name (Elsevier mmcN) counts as the full text only
    # when it prints its DOI at the start and is article-length; a correction does not.
    other = t / "other_name"
    other.mkdir()
    doi = "10.1016/j.molcel.2019.07.015"
    files = {"1-s2.0-S1-mmc7.pdf": "Article\nhttps://doi.org/" + doi + "\n" + "x" * 25000,
             "correction.pdf": "Correction to https://doi.org/" + doi + "\n",
             "1-s2.0-S1-mmc1.pdf": "Supplemental Information\n" + "y" * 25000 + doi}
    rows = ""
    for name, text in files.items():
        (other / name).write_bytes(PDF)
        (other / (name + ".txt")).write_text(text)
        rows += f"{name}\thttps://x/{name}\t1\tabc\tfile\n"
    (other / "MANIFEST.tsv").write_text(base.Folder.HEADER + rows)
    found = base.Folder(other, convert=False)
    check("an article-length original printing the DOI first is its full text; a correction is not",
          gs.article_printing_doi(found, doi), "1-s2.0-S1-mmc7.pdf")
    (other / "1-s2.0-S1-mmc7.pdf.txt").unlink()
    check("...and without it the DOI-citing correction or late mention does not count",
          gs.article_printing_doi(found, doi), "")
    for name, own, text in (("pre.xml", "10.1101/2025.03.14.642805", "bioRxiv 10.1101/2025.03.14.642805 published as " + doi),
                            ("pre.pdf", "", "doi: 10.1101/2025.03.14.642805 now published in " + doi)):
        (other / name).write_bytes(f'<article><article-id pub-id-type="doi">{own}</article-id>'.encode()
                                   if own else PDF)
        (other / (name + ".txt")).write_text(text + "\n" + "w" * 25000)
        found.rows[name] = [name, "https://x/" + name, "1", "abc", "file"]
    check("a preprint that prints its published version's DOI is not that article's full text",
          gs.article_printing_doi(found, doi), "")
    rs = "10.21203/rs.3.pex-1/v1"
    (other / "asset.pdf").write_bytes(PDF)
    (other / "asset.pdf.txt").write_text("DOI: https://doi.org/" + rs + "\n" + "z" * 25000)
    found.rows["asset.pdf"] = ["asset.pdf", "https://assets/asset.pdf", "1", "abd", "asset"]
    found.rows["(manual) rs.3.pex-1_v1"] = ["(manual) rs.3.pex-1_v1", "https://doi.org/" + rs, "", "",
                                             "hand annotation"]
    with contextlib.redirect_stdout(io.StringIO()):
        gs.settle_researchsquare(found, rs)
    check("a Research Square document on disk replaces its manual row by a note keeping the annotation",
          ("(manual) rs.3.pex-1_v1" in found.rows, found.rows.get("(note) rs.3.pex-1_v1", [""] * 5)[4]),
          (False, "hand annotation; on disk as asset.pdf, which prints this DOI on its first pages"))

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

check.section("stalled downloads")
with tempfile.TemporaryDirectory() as tmp, contextlib.redirect_stdout(io.StringIO()):
    folder = base.Folder(Path(tmp), convert=False)
    stalled = []
    def _stall(args):
        stalled.append(args)
        return SimpleNamespace(returncode=28, stdout=b"200", stderr=b"")
    with patch.object(base, "_curl", side_effect=_stall), patch.object(base.time, "sleep"), \
         patch.dict(os.environ, {"CHEM_DATA": tmp}):
        got = folder.save("manual.pdf", "https://vendor.example/manual.pdf", "vendor manual")
    check("a stalled vendor download is tried once and becomes a (manual) row",
          (got, len(stalled), folder.rows.get("(manual) manual.pdf", [""] * 5)[4].endswith("stalled (download incomplete)"),
           "--speed-time" in stalled[0], list(Path(tmp).glob("*.part"))),
          (None, 1, True, True, []))

check.section("source identity: mirrors, back-fill, repository mirrors, same URL")
def _serve(data: bytes):
    """A fake curl download writing `data`; records each requested URL."""
    calls = []
    def fake(url, dest, max_bytes=base.MAX_BYTES, max_time=900):
        calls.append(url)
        Path(dest).write_bytes(data)
        return 200, ""
    return fake, calls
with tempfile.TemporaryDirectory() as tmp, contextlib.redirect_stdout(io.StringIO()):
    folder = base.Folder(Path(tmp), convert=False)
    static = "https://static-content.springer.com/esm/x/MOESM3_ESM.xlsx"
    media = "https://media.springernature.com/original/springer-static/esm/x/MOESM3_ESM.xlsx"
    folder.keep("t_MOESM3_ESM.xlsx", XLSX, static, "Springer supplementary file")
    fake, calls = _serve(XLSX)
    with patch.object(base, "download", side_effect=fake):
        got = folder.save("t_MOESM3_ESM.xlsx", media, "Springer supplementary file")
        again = folder.save("t_MOESM3_ESM.xlsx", media, "Springer supplementary file")
    check("a byte-identical mirror download returns the registered original",
          (got, again, len(calls)), (Path(tmp) / "t_MOESM3_ESM.xlsx",) * 2 + (1,))
    check("...records the alternate URL, not a source collision",
          ([r[1] for k, r in folder.rows.items() if k.startswith("(duplicate) t_MOESM3")],
           any("source-collision" in k for k in folder.rows), folder.rows["t_MOESM3_ESM.xlsx"][1]),
          ([media], False, static))
    other = zbytes({"xl/workbook.xml": b"<workbook>other</workbook>"})
    fake, calls = _serve(other)
    with patch.object(base, "download", side_effect=fake):
        first = folder.save("t_MOESM3_ESM.xlsx", "https://elsewhere/MOESM3_ESM.xlsx", "other")
        second = folder.save("t_MOESM3_ESM.xlsx", "https://elsewhere/MOESM3_ESM.xlsx", "other")
    check("different content at another URL keeps the identity guard (checked once)",
          (first, second, len(calls), (Path(tmp) / "t_MOESM3_ESM.xlsx").read_bytes() == XLSX,
           any("content differs" in r[4] for r in folder.rows.values())),
          (None, None, 1, True, True))
    # a hand-placed original without a URL gets its URL from a byte-identical download
    folder2 = base.Folder(Path(tmp) / "b", convert=False)
    (folder2.path / "hand_table.xlsx").write_bytes(XLSX)
    folder2.rows["hand_table.xlsx"] = ["hand_table.xlsx", "", "1", base.md5(folder2.path / "hand_table.xlsx"),
                                       "(on disk, not fetched by get_sources: added by hand?)"]
    folder2.by_md5[folder2.rows["hand_table.xlsx"][3]] = "hand_table.xlsx"
    folder2.keep("pmc_supp_Table_5.xlsx", XLSX, "https://pmc/Table_5.xlsx", "PMC supplement")
    check("a URL-less hand row is back-filled from an identical download",
          (folder2.rows["hand_table.xlsx"][1], "(duplicate) pmc_supp_Table_5.xlsx" in folder2.rows,
           (folder2.path / "pmc_supp_Table_5.xlsx").exists()),
          ("https://pmc/Table_5.xlsx", True, False))
    # repository mirrors of vendor originals stay visible under scg_ names
    folder3 = base.Folder(Path(tmp) / "c", convert=False)
    (folder3.path / "upstream_data_Guide.pdf").write_bytes(PDF)
    folder3.rows["upstream_data_Guide.pdf"] = ["upstream_data_Guide.pdf", "", "1", base.md5(folder3.path / "upstream_data_Guide.pdf"),
                                               "(on disk, not fetched by get_sources: added by hand?)"]
    folder3.by_md5[folder3.rows["upstream_data_Guide.pdf"][3]] = "upstream_data_Guide.pdf"
    folder3.rows["(duplicate) scg_Guide.pdf"] = ["(duplicate) scg_Guide.pdf", "https://raw/data/Guide.pdf", "",
                                                 folder3.rows["upstream_data_Guide.pdf"][3], "identical to upstream_data_Guide.pdf"]
    fake, calls = _serve(PDF)
    with patch.object(base, "download", side_effect=fake):
        kept = folder3.save("scg_Guide.pdf", "https://raw/data/Guide.pdf", "scg_lib_structs repo file")
    folder3.keep("scg_docs_ge_Guide.md", b"same bytes\n", "https://raw/docs/Guide.md", "repo markdown")
    folder3.keep("upstream_page.html", ARTICLE, "https://page", "scg_lib_structs method page")
    folder3.keep("scg_page_copy.html", ARTICLE, "https://raw/page", "repo copy of a page")
    check("an scg_ repository mirror of an upstream_* vendor file is kept under its scg_ name",
          (kept, "(duplicate) scg_Guide.pdf" in folder3.rows, (folder3.path / "upstream_data_Guide.pdf").exists()),
          (folder3.path / "scg_Guide.pdf", False, True))
    check("...but never for the drawn page itself",
          ((folder3.path / "scg_page_copy.html").exists(), "(duplicate) scg_page_copy.html" in folder3.rows),
          (False, True))
    check("mirror-name rule excludes the repository's own docs and adapters",
          [base.Folder.mirror_name_kept(n, "upstream_data_x", "")
           for n in ("scg_x.pdf", "scg_docs_ge_x.md", "scg_adapters_x.fa")], [True, False, False])
    # the same URL requested under a second name is not a second original
    folder4 = base.Folder(Path(tmp) / "d", convert=False)
    folder4.keep("p_v1.full.pdf", PDF, "https://biorxiv/p.v1.full.pdf", "biorxiv full text")
    fake, calls = _serve(PDF + b"% stamped differently\n%%EOF\n")
    with patch.object(base, "download", side_effect=fake):
        folder4.save("pv1.full.pdf", "https://biorxiv/p.v1.full.pdf", "extra source")
    check("a second name for an already saved URL becomes a (duplicate) row without a download",
          (calls, (folder4.path / "pv1.full.pdf").exists(), folder4.rows["(duplicate) pv1.full.pdf"][3]),
          ([], False, folder4.rows["p_v1.full.pdf"][3]))
    check("held() resolves a (duplicate) name to the file holding its bytes",
          (folder4.held("pv1.full.pdf"), folder4.held("absent.pdf")), (folder4.path / "p_v1.full.pdf", None))

with tempfile.TemporaryDirectory() as tmp, contextlib.redirect_stdout(io.StringIO()):
    folder = base.Folder(Path(tmp), convert=False)
    folder.manual("PMC4224642 PMC Cloud source", f"{pmc_cloud.BUCKET}/?list-type=2&prefix=PMC4224642.",
                  "defining manuscript", "no usable file returned")
    folder.manual("PMC5024639 PMC Cloud source", f"{pmc_cloud.BUCKET}/?list-type=2&prefix=PMC5024639.",
                  "other manuscript", "no usable file returned")
    folder.keep("nature13802_PMC4224642.pdf", PDF, f"{pmc_cloud.BUCKET}/PMC4224642.1/PMC4224642.1.pdf",
                "article PDF (PMC Cloud)")
    check("a stale PMC Cloud manual row goes once that article's bucket file is on disk",
          (pmc_cloud.settle_manual_rows(folder), "(manual) PMC4224642 PMC Cloud source" in folder.rows,
           "(manual) PMC5024639 PMC Cloud source" in folder.rows), (1, False, True))

check.section("legacy Office content, error bodies, HTTP/2 fallback")
OLE = base.OLE_MAGIC + b"\0" * 600
with tempfile.TemporaryDirectory() as tmp, contextlib.redirect_stdout(io.StringIO()):
    t = Path(tmp)
    (t / "w.docx").write_bytes(OLE)
    (t / "z.zip").write_bytes(OLE)
    check("an OLE2 Word file under a .docx name is a genuine Office file; a .zip is not",
          (base.validate(t / "w.docx"), bool(base.validate(t / "z.zip"))), ("", True))
    check("content name gives OLE2 bytes their legacy suffix",
          (base.content_name("supp.docx", t / "w.docx"), base.content_name("supp.pdf", t / "w.docx")),
          ("supp.doc", "supp.pdf"))
    folder = base.Folder(t / "o", convert=False)
    stored = folder.keep("supp_R2.docx", OLE, "https://pmc/supp_R2.docx", "supplementary file")
    check("a legacy .doc served as .docx is stored as .doc with a (duplicate) row for the served name",
          (stored, folder.rows["(duplicate) supp_R2.docx"][3] == folder.rows["supp_R2.doc"][3],
           folder.held("supp_R2.docx")),
          (folder.path / "supp_R2.doc", True, folder.path / "supp_R2.doc"))
    (t / "expired.html").write_bytes(b"<html><body>PDF link has expired or is invalid" + b" x" * 1000 + b"</body></html>")
    check("a Europe PMC 'PDF link has expired or is invalid' body is no document",
          "expired" in base.validate(t / "expired.html", min_text=10), True)
runs = []
def _curl92(args, **kwargs):
    runs.append(args)
    return SimpleNamespace(returncode=92 if len(runs) == 1 else 0, stdout=b"200", stderr=b"")
with tempfile.TemporaryDirectory() as tmp, patch.object(base.subprocess, "run", side_effect=_curl92), \
        patch.object(base.time, "sleep"), patch.dict(os.environ, {"CHEM_DATA": tmp}):
    base.download("https://vendor.example/manual.pdf", Path(tmp) / "m.pdf")
check("curl HTTP/2 stream error (92) is retried once over HTTP/1.1",
      (len(runs), "--http1.1" in runs[0], "--http1.1" in runs[1]), (2, False, True))

check.section("protocols.io: cited DOI answered by a cached export")
_export = ("Version 2\nDiMeLo-seq V.2\nDOI\nhttps://dx.doi.org/10.17504/protocols.io.b2u8qezw\n"
           "DOI: https://dx.doi.org/10.17504/protocols.io.b2u8qezw\nExternal link: https://doi.org/10.1038/x\n"
           "see also dx.doi.org/10.17504/protocols.io.zcpf2vn\nProtocol Integer ID: 55936\n")
check("printed DOI is the export's own 'DOI:' line, not a DOI it cites",
      protocolsio.printed_doi(_export), "10.17504/protocols.io.b2u8qezw")
check("cited export identity: accession and revision",
      [protocolsio.cited_export(c, p, v) for c, p, v in (
          ("10.17504/protocols.io.b2u8qezw", "10.17504/protocols.io.b2u8qezw", "2"),
          ("10.17504/protocols.io.rm7vzqr92vx1/v1", "10.17504/protocols.io.rm7vzqr92vx1/v1", ""),
          ("10.17504/protocols.io.dm6gp3o68vzp", "10.17504/protocols.io.dm6gp3o68vzp/v1", ""),
          ("10.17504/protocols.io.abc/v3", "10.17504/protocols.io.abc/v4", ""),
          ("10.17504/protocols.io.abc/v3", "10.17504/protocols.io.abc", "4"),
          ("10.17504/protocols.io.zcpf2vn", "10.17504/protocols.io.b2u8qezw", "2"))],
      [True, True, True, False, False, False])
check("a cached export with conflicting DOI/Version labels answers neither citation",
      [protocolsio.cited_export(c, "10.17504/protocols.io.abc/v1", "2")
       for c in ("10.17504/protocols.io.abc/v1", "10.17504/protocols.io.abc")],
      [False, False])
with patch.object(base, "cooling_down", return_value=""), \
     patch.object(base, "throttle"), \
     patch.object(base, "_curl", return_value=SimpleNamespace(
         stdout=b"<html><title>Just a moment</title></html>\n200")), \
     patch.object(base, "mark_blocked") as blocked:
    base.get("https://publisher.example/api")
    check("small API/HTML requests also put a challenged host in cool-down",
          (blocked.call_count, blocked.call_args.args[0]), (1, "https://publisher.example/api"))
with tempfile.TemporaryDirectory() as tmp, contextlib.redirect_stdout(io.StringIO()):
    f = base.Folder(Path(tmp), convert=False)
    pdf = Path(tmp) / "protocolsio_wrong_v1.pdf"
    pdf.write_bytes(PDF)
    pdf.with_name(pdf.name + ".txt").write_text(
        "Version 1\nDOI: https://doi.org/10.17504/protocols.io.wrong/v1\n"
        "Protocol Integer ID: 123\n")
    with patch.object(protocolsio, "resolve", return_value="https://www.protocols.io/view/wrong/v1"), \
         patch.object(f, "save", return_value=pdf), \
         patch.object(protocolsio, "public_attachments") as metadata:
        protocolsio.fetch(SimpleNamespace(folder=f), "10.17504/protocols.io.expected/v1")
    check("matching revision of a different protocol cannot supply attachments",
          (metadata.call_count, pdf.exists(), (Path(tmp) / "_invalid" / pdf.name).exists()),
          (0, False, True))
with tempfile.TemporaryDirectory() as tmp, contextlib.redirect_stdout(io.StringIO()):
    f = base.Folder(Path(tmp), convert=False)
    pdf = Path(tmp) / "protocolsio_dimelo-seq-n2bvjxe4wlk5.pdf"
    pdf.write_bytes(PDF)
    pdf.with_name(pdf.name + ".txt").write_text(_export)
    f.rows[pdf.name] = [pdf.name, "https://www.protocols.io/view/dimelo-seq-n2bvjxe4wlk5.pdf", "1", base.md5(pdf), "protocol"]
    f.manual("protocolsio_b2u8qezw.pdf", "https://doi.org/10.17504/protocols.io.b2u8qezw", "cited", "DOI did not resolve")
    with patch.object(protocolsio, "resolve", side_effect=AssertionError("no network")):
        protocolsio.fetch(SimpleNamespace(folder=f), "10.17504/protocols.io.b2u8qezw")
    check("a cited DOI printed by a cached export clears its stale manual row without a request",
          ("(manual) protocolsio_b2u8qezw.pdf" in f.rows,
           any(pdf.name in r[4] for k, r in f.rows.items() if k == "(note) protocols.io b2u8qezw")),
          (False, True))
    check("an export printing another DOI does not answer a citation",
          protocolsio.cached_export(f, "10.17504/protocols.io.zcpf2vn"), None)
    _original_bytes = pdf.read_bytes()
    pdf.write_bytes(_original_bytes + b"changed")
    check("a cached export with changed bytes does not clear an unavailable citation",
          protocolsio.cached_export(f, "10.17504/protocols.io.b2u8qezw"), None)
    pdf.write_bytes(_original_bytes)
    f.manual("protocolsio_b2u8qezw.pdf", "https://doi.org/10.17504/protocols.io.b2u8qezw", "cited", "DOI did not resolve")
    f.manual("protocolsio_zcpf2vn.pdf", "https://doi.org/10.17504/protocols.io.zcpf2vn", "cited", "DOI did not resolve")
    check("a sweep clears only manual DOI rows that a cached export answers",
          (protocolsio.settle_manual_rows(f), "(manual) protocolsio_b2u8qezw.pdf" in f.rows,
           "(manual) protocolsio_zcpf2vn.pdf" in f.rows), (1, False, True))

check.section("doctext: legacy Word, streamed XLSX, OCR fallback")
with tempfile.TemporaryDirectory() as tmp:
    t = Path(tmp)
    doc = t / "supp.docx"
    doc.write_bytes(OLE)
    with patch.object(doctext, "_legacy", return_value="legacy text\n") as legacy:
        check("an OLE2 Word file under a .docx name is read as legacy .doc",
              (doctext.read(doc), legacy.call_args.args[1]), ("legacy text\n", ".doc"))
    with patch.object(doctext, "SOFFICE", None), patch.object(doctext, "CATDOC", "catdoc"), \
            patch.object(doctext.subprocess, "run",
                         return_value=SimpleNamespace(returncode=0, stdout=b"Primer ACGT\n")) as run:
        check("catdoc reads a .doc when LibreOffice is absent",
              (doctext._legacy(doc, ".doc"), run.call_args.args[0][0]), ("Primer ACGT\n", "catdoc"))
    with patch.object(doctext, "SOFFICE", None), patch.object(doctext, "CATDOC", None), \
            contextlib.redirect_stderr(io.StringIO()):
        check("without converters a .doc is skipped, not faked", doctext._legacy(doc, ".doc"), None)
    sheet = (b'<worksheet><sheetData><row r="1"><c r="A1" t="s"><v>0</v></c><c r="B1" t="inlineStr">'
             b'<is><t>ACGT &amp; T</t></is></c></row>' + b'<row r="2" s="3" customFormat="1"/>' * 3000 +
             b'<row r="3001"><c r="A3001"><v>42</v></c></row></sheetData></worksheet>')
    book = t / "big.xlsx"
    book.write_bytes(zbytes({"xl/workbook.xml": b'<workbook><sheets><sheet name="Oligos" sheetId="1"/></sheets></workbook>',
                             "xl/sharedStrings.xml": b"<sst><si><t>Name</t></si></sst>",
                             "xl/worksheets/sheet1.xml": sheet}))
    check("streamed XLSX skips self-closing empty rows and keeps cell content",
          doctext.xlsx_text(book), "## sheet Oligos\nName\tACGT & T\n42\n")
    with patch.object(doctext, "_read", side_effect=AssertionError("unstreamed")):
        check(".xlsx goes through the streamed reader", doctext.read(book).splitlines()[0], "## sheet Oligos")
    pdf = t / "scan.pdf"
    pdf.write_bytes(PDF)
    with patch.object(doctext, "TESSERACT", None), contextlib.redirect_stderr(io.StringIO()):
        check("without OCR tools an empty text layer is returned unchanged", doctext.with_ocr(pdf, "\f\f"), "\f\f")
    with patch.object(doctext, "TESSERACT", "tesseract"), patch.object(doctext, "PDFTOPPM", "pdftoppm"), \
            patch.object(doctext, "_page_sizes", return_value={}), \
            patch.object(doctext, "ocr_page", side_effect=lambda p, n, longest=792: f"page {n} 5'-ACGT-3'\n"):
        full = doctext.with_ocr(pdf, "\f\f")
        layer = "Real text layer " * 40 + "\f" + "Header only\f"
        with patch.object(doctext, "raster_pages", return_value=[2]) as raster:
            partial = doctext.with_ocr(pdf, layer)
    check("an image-only PDF twin starts with the OCR header, then page text",
          (full.splitlines()[0].startswith("## OCR TEXT"), "page 1 5'-ACGT-3'" in full, "page 2" in full),
          (True, True, True))
    check("an image-only page in a text PDF is appended after the unchanged text layer",
          (partial.startswith(layer.rstrip("\n")), "## OCR TEXT of image-only page 2" in partial,
           raster.call_args.args[1]),
          (True, True, [2]))

check.section("catalogue paper linking")
_cap2010 = {"title": ["Orphan CpG Islands Identify Numerous Conserved Promoters in the Mammalian Genome"],
            "author": [{"family": "Illingworth"}, {"family": "Bird"}], "issued": {"date-parts": [[2010]]}}
_cap2024 = {"title": ["Large-scale single-cell long-read genomics enables high-resolution microbiome profiling"],
            "author": [{"family": "Li"}, {"family": "Zheng"}], "posted": {"date-parts": [[2024]]}}
_science = {"title": ["Droplet-based forward genetic screening of astrocyte-microglia cross-talk"],
            "author": [{"family": "Wheeler"}, {"family": "Quintana"}], "issued": {"date-parts": [[2023]]}}
_protocol = {"title": ["SPEAC-seq: a protocol for screening astrocyte-microglia interactions"],
             "author": [{"family": "Lee"}, {"family": "Wheeler"}, {"family": "Quintana"}],
             "issued": {"date-parts": [[2024]]}}
check("a shared method name with no shared title or authors is no link",
      gs.works_agree(_cap2010, _cap2024), "")
check("a companion paper with shared senior author within three years is linked",
      bool(gs.works_agree(_protocol, _science)), True)
check("unknown metadata never establishes a link", gs.works_agree({}, _science), "")
check("the defining DOI itself needs no lookup",
      gs.linked_to_defining("10.1371/journal.pgen.1001134", "https://doi.org/10.1371/JOURNAL.PGEN.1001134"),
      "defining DOI")
check("first-author agreement for title-matched preprints",
      [gs._same_first_author(a, b) for a, b in (("Li", "Li Y"), ("Illingworth", "Li Y"), ("", "Li"))],
      [True, False, True])

check.section("Polite access: per-host interval and challenge cool-down")
with tempfile.TemporaryDirectory() as folder, patch.dict(os.environ, {"CHEM_DATA": folder}):
    check("publisher hosts get the slow interval, others the default",
          (base.host_interval("link.springer.com"), base.host_interval("www.ncbi.nlm.nih.gov")),
          (base.SLOW_INTERVAL, base.HOST_INTERVAL))
    slept = []
    with patch.object(base.time, "sleep", side_effect=slept.append):
        base.throttle("https://link.springer.com/a.pdf")
        base.throttle("https://link.springer.com/b.pdf")
    check("a second request to the same host waits for the interval", bool(slept) and slept[-1] > 7, True)
    base.mark_blocked("https://www.science.org/doi/pdf/x", "challenge/gate page (Just a moment)")
    with patch.object(base, "_curl") as curl:
        code, problem = base.download("https://www.science.org/doi/suppl/y", Path(folder) / "y.pdf")
    check("a host in cool-down is not requested again", (curl.called, "cool-down" in problem), (False, True))
    check("cool-down is per host", base.cooling_down("https://www.ncbi.nlm.nih.gov/x"), "")
    with patch.object(base.time, "time", return_value=base.time.time() + base.COOLDOWN_SECONDS + 1):
        check("cool-down expires", base.cooling_down("https://www.science.org/doi/z"), "")
    fold = base.Folder(Path(folder) / "sources" / "slug", convert=False, expand=False)
    for reason, url in (("PMC 'Preparing to download' gate", "https://pmc.ncbi.nlm.nih.gov/a.pdf"),
                        ("challenge/gate page (Client Challenge)", "https://link.springer.com/c.pdf")):
        with patch.object(base, "download", return_value=(200, "")) as dl, \
             patch.object(base, "validate", return_value=reason), \
             patch.object(base.time, "sleep"):
            fold._fetch(url, Path(folder) / "part", "c.pdf", 1500, None, base.MAX_BYTES, 2)
        check(f"{reason.split()[0]} gate: fetched once, host cool-down {'pmc' not in url}",
              (dl.call_count, bool(base.cooling_down(url))), (1, "pmc" not in url))

check.report()
