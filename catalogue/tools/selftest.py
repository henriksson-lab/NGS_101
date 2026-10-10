#!/usr/bin/env python3
"""
Self-test for the protocol catalogue.

The table is scraped, so these are not checks on typing but on the scrape's invariants:
that the many-to-many shape holds, that the directory-naming scheme round-trips on real
DOIs, and that no row asserts an identifier nothing confirmed.

Run:  python3 catalogue/tools/selftest.py
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "lib"))

import catalogue as cat  # noqa: E402
from checks import Check  # noqa: E402
from fetch_scg_lib_structs import (COLUMNS, decode_doi, doi_of, encode_doi,  # noqa: E402
                                   family_members, name_slug, page_papers, short_name,
                                   slug)

check = Check()
rows = cat.rows()

check.section("the table loads and is shaped as declared")
check("the TSV exists", cat.TSV.exists())
check("every declared column is present", list(rows[0].keys()), COLUMNS)
check("every row names its source", sorted({r["source"] for r in rows}),
      ["ours", "scg_lib_structs"])
check("no row is missing a category", [r for r in rows if not r["category"]], [])
check("no row is missing a protocol", [r for r in rows if not r["protocol"]], [])
check("every row is documented=yes, todo or ours",
      sorted({r["documented"] for r in rows}), ["ours", "todo", "yes"])
check("every row has a role", sorted({r["role"] for r in rows}), ["cited", "primary"])
check("there are rows from both the pages and the TODO list",
      bool(cat.documented()) and bool(cat.todo()))

check.section("one row per (protocol, paper); the relation is many-to-many")
check("protocol is NOT a key -- some have several papers", cat.n_rows() > cat.n_protocols())
check("...and that is listed, not hidden", len(cat.multi_paper()) > 0)
check("a paper can define more than one protocol", len(cat.shared_papers()) > 0)
pairs = [(r["protocol"], r["doi"], r["paper_url"]) for r in rows]
check("no duplicate (protocol, paper) pair", len(pairs), len(set(pairs)))
check("every upstream-documented row links its upstream page",
      [r["protocol"] for r in cat.documented()
       if not r["scg_page"].startswith("https://teichlab.github.io/scg_lib_structs/")], [])

check.section("directory naming round-trips on every real DOI")
dois = sorted({r["doi"] for r in rows if r["doi"]})
check("there are DOIs to test", len(dois) > 100)
bad = [d for d in dois if decode_doi(encode_doi(d)) != d]
check("encode -> decode is the identity for all of them", bad, [])
check("a DOI whose suffix contains an underscore still round-trips",
      decode_doi(encode_doi("10.1038/nbt0400_424")), "10.1038/nbt0400_424")
check("...because '_' is left alone and only '/' is substituted",
      encode_doi("10.1038/nbt0400_424"), "10.1038+nbt0400_424")
check("a DOI with TWO slashes round-trips as well (Research Square)",
      encode_doi("10.21203/rs.3.rs-3210240/v1"), "10.21203+rs.3.rs-3210240+v1")
check("no DOI in the table contains the '+' the scheme reserves",
      [d for d in dois if "+" in d], [])
check("slugs are filesystem-safe",
      [s for s in {r["slug"] for r in rows} if not re.fullmatch(r"[A-Za-z0-9._+-]+", s)], [])
check("the name and the DOI are separated by '__'",
      [r["slug"] for r in cat.defining() if r["doi"] and "__" not in r["slug"]], [])
check("a slugified name never contains '__', so the split is unambiguous",
      [r["protocol"] for r in rows if "__" in name_slug(r["protocol"])], [])

check.section("one directory per protocol, named for its defining paper")
multi = [p for p in cat.protocols()
         if len({r["slug"] for r in rows if r["protocol"] == p}) > 1]
check("every row of a protocol carries the same slug", multi, [])
check("distinct slugs == distinct protocols", len({r["slug"] for r in rows}),
      cat.n_protocols())
check("exactly one defining row per protocol that has a paper",
      len(cat.defining()), len({r["protocol"] for r in cat.primary()}))
check("a defining row is always a primary one",
      [r["protocol"] for r in cat.defining() if r["role"] != "primary"], [])
for r in cat.defining():
    if r["doi"]:
        want = f"{name_slug(r['protocol'])}__{encode_doi(r['doi'])}"
        if r["slug"] != want:
            check(f"{r['protocol']}: slug is built from its defining DOI", r["slug"], want)
check("slug built from the defining DOI (checked above for every defining row)", True)
check("each SMART-seq protocol is named for its own defining paper",
      [slug(name, doi) for name, doi in
       [("SMART-seq", "10.1038/nbt.2282"), ("SMART-seq2", "10.1038/nmeth.2639"),
        ("SMART-seq3", "10.1038/s41587-020-0497-0"),
        ("SMART-seq3xpress", "10.1038/s41587-022-01311-4"),
        ("FLASH-seq", "10.1038/s41587-022-01312-3")]],
      ["smart-seq__10.1038+nbt.2282", "smart-seq2__10.1038+nmeth.2639",
       "smart-seq3__10.1038+s41587-020-0497-0",
       "smart-seq3xpress__10.1038+s41587-022-01311-4",
       "flash-seq__10.1038+s41587-022-01312-3"])

# A page's preamble cites the methods a protocol is built from next to its own paper, and
# those are older -- so "earliest" once named ISSAAC-seq after ATAC-seq and three
# multi-omics methods after Smart-seq2. One paper defining two protocols is rare and real.
JOINT = {"10.1038/ncomms14049",          # 10x 3' GE V1 and V2-V4: one Zheng 2017 paper
         "10.7554/eLife.73971",          # HyDrop-RNA and HyDrop-ATAC
         "10.1038/s41587-021-00962-z",   # s3-ATAC and s3-WGS
         "10.1038/s41587-019-0147-6",    # dscATAC-seq and dsciATAC-seq
         "10.1101/2024.02.12.579864",    # scTAPS and scCAPS+
         "10.1038/s41587-021-00927-2"}   # ASAP-seq and DOGMA-seq
by_doi: dict[str, list[str]] = {}
for r in cat.defining():
    if r["doi"]:
        by_doi.setdefault(r["doi"], []).append(r["protocol"])
check("no two protocols share a defining paper, except the known joint papers",
      {d: ps for d, ps in by_doi.items() if len(ps) > 1 and d not in JOINT}, {})
check("...and each joint paper really is shared",
      sorted(d for d in JOINT if len(by_doi.get(d, [])) < 2), [])
_defd = {r["protocol"]: r for r in cat.defining()}
for proto, doi in [("ISSAAC-seq", "10.1038/s41592-022-01601-4"),
                   ("scNMT-seq", "10.1038/s41467-018-03149-4"),
                   ("SNARE-seq", "10.1038/s41587-019-0290-0"),
                   ("scDamID", "10.1016/j.cell.2015.08.040")]:
    check(f"{proto} is named for its own paper, not a component's", _defd[proto]["doi"], doi)
check("the background papers stay in the table as associated rows",
      sorted(r["protocol"] for r in rows if r["doi"] == "10.1038/nmeth.2639"
             and r["is_defining"] == "no"
             and r["protocol"] not in ("SMART-seq family", "SMART-seq2")),
      ["scM&T-seq", "scMT-seq", "scNMT-seq"])
check("a preprint of the method's own paper names it over the journal version",
      _defd["scifi-RNA-seq"]["journal"], "bioRxiv")
check("a vendor kit citing only techniques it reads out has no defining paper",
      [p for p in ("10x Chromium Single Cell ATAC",
                   "10x Chromium Single Cell 3' FeatureBarcoding") if _defd[p]["doi"]], [])

check.section("identifiers are confirmed, never guessed")
check("no row carries a DOI without a resolved title",
      [r["protocol"] for r in rows if r["doi"] and not r["title"]], [])
check("a dropped DOI leaves a note saying so",
      all(r["note"] for r in rows if r["paper_url"] and not r["doi"] and not r["pmid"]))
check("every DOI looks like a DOI",
      [r["doi"] for r in rows if r["doi"] and not re.fullmatch(r"10\.\d{4,9}/\S+", r["doi"])],
      [])
check("every PMID is digits",
      [r["pmid"] for r in rows if r["pmid"] and not r["pmid"].isdigit()], [])
check("a row with no paper at all explains why -- a vendor kit, or unpublished work of ours",
      [r["protocol"] for r in rows if not r["paper_url"] and not r["note"]], [])
check("...and the only unpublished one is ours",
      [r["protocol"] for r in rows
       if not r["paper_url"] and "vendor" not in r["note"] and r["source"] != "ours"], [])

check.section("the parsing rules that classify a page's links")
check("a DOI inside a publisher URL is read straight out of it",
      doi_of("https://www.science.org/doi/10.1126/science.aam8999"),
      "10.1126/science.aam8999")
check("a Nature article id implies its DOI",
      doi_of("https://www.nature.com/articles/s41592-018-0107-y"),
      "10.1038/s41592-018-0107-y")
check("a bioRxiv URL gives the preprint DOI without its version suffix",
      doi_of("https://www.biorxiv.org/content/10.1101/2023.01.12.523500v1"),
      "10.1101/2023.01.12.523500")
check("a family page's members are parsed from its title",
      family_members("SMART-seq family (including SMART-seq, SMART-seq2/3 and FLASH-seq)"),
      ["SMART-seq", "SMART-seq2/3", "FLASH-seq"])
check("the display name drops the '(including ...)' part",
      short_name("Quartz-seq family (including Quartz-seq and Quartz-seq2)"),
      "Quartz-seq family")

_PAGE = """<h1>Demo-seq</h1><p>Demo-seq (<a href="https://doi.org/10.1/own">Nat 2020</a>)
is a method. See <a href="https://doi.org/10.9/fam">Demo-seq2</a>.</p>
<h2>Adapter and primer sequences:</h2>
<p>uses <a href="https://doi.org/10.2/aside">semi-suppressive PCR</a></p>"""
roles = {u: role for _, u, role in page_papers(_PAGE, ["Demo-seq", "Demo-seq2"])}
check("a paper cited in the preamble is the method's own",
      roles["https://doi.org/10.1/own"], "primary")
check("a paper whose anchor names a family member is too",
      roles["https://doi.org/10.9/fam"], "primary")
check("a paper cited from inside a step is incidental",
      roles["https://doi.org/10.2/aside"], "cited")
check("the real table classifies the semi-suppressive-PCR paper as incidental everywhere",
      sorted({r["role"] for r in rows if r["doi"] == "10.1038/nmeth.1470"}), ["cited"])
check("...and it is not the defining paper of anything",
      [r["protocol"] for r in cat.defining() if r["doi"] == "10.1038/nmeth.1470"], [])

check.section("our own coverage (catalogue/ours.tsv)")
check("ours.tsv exists", cat.OURS.exists())
check("every row has a directory, a protocol and a status",
      [r for r in cat.ours() if not (r["dir"] and r["protocol"] and r["status"])], [])
check("every directory named in ours.tsv exists on disk",
      [r["dir"] for r in cat.ours() if not (ROOT / r["dir"]).is_dir()], [])
check("status is one of documented / draft / notes",
      sorted({r["status"] for r in cat.ours()} - {"documented", "draft", "notes"}), [])
check("a 'notes' directory has at least one reference note",
      [r["dir"] for r in cat.ours() if r["status"] == "notes"
       and not list((ROOT / r["dir"]).glob("0*.md"))], [])
check("every documented or draft directory has a build script and a selftest",
      [r["dir"] for r in cat.ours() if r["status"] != "notes"
       and not ((ROOT / r["dir"] / "tools" / "build_page.py").exists()
               and (ROOT / r["dir"] / "tools" / "selftest.py").exists())], [])
# ``dir`` is the explicit stable identifier.  Commercial protocols intentionally use
# concise product/version identifiers (for example ``3-ge-v3`` rather than spelling out
# “Gene Expression”), so reconstructing it from display prose would be a second and
# conflicting source of truth.  DOI-bearing identifiers are still checked reversibly.
check("a directory carrying a DOI encodes it reversibly",
      [r["dir"] for r in cat.ours()
       if r["doi"] and decode_doi(r["dir"].split("__", 1)[1]) != r["doi"]], [])
check("no two of our directories collide",
      len({r["dir"] for r in cat.ours()}), len(cat.ours()))
check("every protocol directory on disk is declared in ours.tsv",
      sorted(p.name for p in ROOT.iterdir()
             if p.is_dir() and (p / "tools" / "build_page.py").exists()
             and p.name not in {r["dir"] for r in cat.ours()}), [])
check.section("website sections and modality (catalogue/ours.tsv)")
from fetch_scg_lib_structs import OURS_COLUMNS, OURS_VOCAB, read_ours  # noqa: E402

check("ours.tsv has exactly the declared columns, in order",
      list(cat.ours()[0].keys()), cat.OURS_COLUMNS)
check("the fetcher and catalogue.py agree on those columns", OURS_COLUMNS, cat.OURS_COLUMNS)
check("...and on each closed vocabulary",
      OURS_VOCAB, {"status": cat.STATUSES, "section": cat.SECTIONS,
                   "modality": cat.MODALITIES})
check("the fetcher reads ours.tsv without complaint", len(read_ours()), cat.n_ours())
check("section is one of published / wip",
      sorted({r["dir"] for r in cat.ours() if r["section"] not in cat.SECTIONS}), [])
# Declared here on purpose: moving a directory between website sections is a decision,
# so it must be made twice -- in ours.tsv and in this list.
WIP = {"crispr-mip__10.1101+2024.03.28.587082"}   # our own method in progress
check("the work-in-progress section is exactly the declared set",
      sorted(r["dir"] for r in cat.wip()), sorted(WIP))
check("every other directory is a published protocol",
      len(cat.published()) + len(cat.wip()), cat.n_ours())
check("modality is one of DNA / RNA / multi",
      sorted({r["dir"] for r in cat.ours() if r["modality"] not in cat.MODALITIES}), [])
_cat_of = {r["our_dir"]: r["category"] for r in cat.covered()}
check("modality agrees with the scg_lib_structs category wherever that implies one",
      [(r["dir"], r["modality"], _cat_of.get(r["dir"])) for r in cat.ours()
       if _cat_of.get(r["dir"]) in cat.CATEGORY_MODALITY
       and cat.CATEGORY_MODALITY[_cat_of[r["dir"]]] != r["modality"]], [])
try:
    cat.ours_in("drafts")
    _refused = False
except ValueError:
    _refused = True
check("ours_in() refuses a section that does not exist", _refused)

check.section("our own coverage, continued")
check("the five SMART-seq protocols are separate entries in ours.tsv",
      {"SMART-seq", "SMART-seq2", "SMART-seq3", "SMART-seq3xpress", "FLASH-seq"}
      <= {r["protocol"] for r in cat.ours()})
check("the table records our coverage", len(cat.covered()) > 0)
check("a directory named for a catalogue slug is joined to that protocol",
      [r["dir"] for r in cat.ours() if r["dir"] in {x["slug"] for x in cat.rows()}
       and r["dir"] not in {x["our_dir"] for x in cat.rows() if x["slug"] == r["dir"]}], [])
check("...for every directory in ours.tsv",
      sorted({r["our_dir"] for r in cat.covered()}),
      sorted(r["dir"] for r in cat.ours()))
check("a covered row's our_status comes from ours.tsv",
      sorted({(r["our_dir"], r["our_status"]) for r in cat.covered()}),
      sorted((r["dir"], r["status"]) for r in cat.ours()))
check("protocols only we cover are appended with source=ours",
      {r["source"] for r in cat.rows() if r["documented"] == "ours"}, {"ours"})
check("every appended row is its own defining row",
      [r["protocol"] for r in cat.rows()
       if r["source"] == "ours" and r["is_defining"] != "yes"], [])
check("an appended row's slug IS its directory",
      [r["protocol"] for r in cat.rows()
       if r["source"] == "ours" and r["slug"] != r["our_dir"]], [])
check("the overlap with upstream is by DOI and non-empty", len(cat.overlap()) > 0)
check("untackled + covered defining rows == all defining rows",
      len(cat.untackled()) + len([r for r in cat.defining() if r["our_dir"]]),
      len(cat.defining()))
check("SPLiT-seq, the archived style exemplar, is in the catalogue",
      any(p.startswith("SPLiT-seq") for p in cat.protocols()))

sys.exit(check.report())
