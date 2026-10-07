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
      [r["slug"] for r in rows if r["doi"] and "__" not in r["slug"]], [])
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
check("SMART-seq family is named after the 2012 Nature Biotech paper, its earliest",
      slug("SMART-seq family", "10.1038/nbt.2282"), "smart-seq-family__10.1038+nbt.2282")

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
check("every one of those directories has a build script and a selftest",
      [r["dir"] for r in cat.ours()
       if not ((ROOT / r["dir"] / "tools" / "build_page.py").exists()
               and (ROOT / r["dir"] / "tools" / "selftest.py").exists())], [])
check("each directory name is rebuilt exactly by the naming scheme",
      [r["dir"] for r in cat.ours() if r["dir"] != slug(r["protocol"], r["doi"])], [])
check("a directory carrying a DOI encodes it reversibly",
      [r["dir"] for r in cat.ours()
       if r["doi"] and decode_doi(r["dir"].split("__", 1)[1]) != r["doi"]], [])
check("no two of our directories collide",
      len({r["dir"] for r in cat.ours()}), len(cat.ours()))
check("every protocol directory on disk is declared in ours.tsv",
      sorted(p.name for p in ROOT.iterdir()
             if p.is_dir() and (p / "tools" / "build_page.py").exists()
             and p.name not in {r["dir"] for r in cat.ours()}), [])
check("SMART-seq family is in both the catalogue and ours.tsv",
      "SMART-seq family" in cat.protocols()
      and "SMART-seq family" in {r["protocol"] for r in cat.ours()})
check("the table records our coverage", len(cat.covered()) > 0)
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
