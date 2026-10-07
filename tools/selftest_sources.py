#!/usr/bin/env python3
"""Offline self-test for the source fetcher (tools/get_sources.py, tools/fetchers/,
tools/doctext.py XML handling, catalogue/extra_sources.tsv). No network: every input is
synthetic and built here.

Run:  python3 tools/selftest_sources.py
"""
from __future__ import annotations

import contextlib
import gzip
import io
import sys
import tempfile
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / "lib"))
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "catalogue" / "tools"))

from checks import Check  # noqa: E402

import doctext  # noqa: E402
import fetchers  # noqa: E402
import get_sources as gs  # noqa: E402
from fetchers import base, pmc_cloud, protocolsio, scg_github, vendor  # noqa: E402

check = Check()

# ------------------------------------------------------------------------------
check.section("DOI routing")
check("bioRxiv dated DOI", gs.doi_kind("10.1101/2019.12.17.879304"), "biorxiv")
check("bioRxiv 6-digit DOI", gs.doi_kind("10.1101/003236"), "biorxiv")
check("bioRxiv DOI with a version suffix", gs.doi_kind("10.1101/2023.05.11.540245v2"), "biorxiv")
check("Genome Research (10.1101/gr.) is a journal", gs.doi_kind("10.1101/gr.110882.110"), "journal")
check("Genes & Dev (10.1101/gad.) is a journal", gs.doi_kind("10.1101/gad.350434.123"), "journal")
check("CSH Protocols (10.1101/pdb.) is a journal", gs.doi_kind("10.1101/pdb.prot5384"), "journal")
check("protocols.io", gs.doi_kind("10.17504/protocols.io.b4xwqxpe"), "protocolsio")
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
import catalogue as cat  # noqa: E402
_slugs = {r["slug"] for r in cat.rows()}
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

check.report()
