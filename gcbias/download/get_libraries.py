#!/usr/bin/env python3
"""Fetch the published sgRNA library whitelists into $CHEM_DATA. Never into the repo.

Why a whitelist matters. Guide reads carry basecalling errors, and the error rate rises
steeply with GC (Spearman +0.59 on the Schmierer input). Reads whose guide is miscalled do
not match the library and are silently dropped, so GC-rich guides lose more reads -- which
looks exactly like reduced amplification. Correcting that needs the real library.

Deriving the whitelist from the data instead (keep every 20-mer above a read threshold) is
CIRCULAR, because the threshold itself is depth-dependent and depth is what we are
measuring. Checked against the published Schmierer list: a >=100-read threshold recovered
23,272 of 23,279 guides, but the 7 it missed had median GC 0.70 while the 19 false ones it
kept had median GC 0.35. Small, but biased in precisely the direction under study.

Sources (recorded in MANIFEST.json alongside the files):

  Schmierer 2017, Mol Syst Biol 13:945, Dataset EV1 (MSB-13-945-s002.csv)
      via the EuropePMC REST supplementaryFiles endpoint for PMC5658704.
      Format is already the lib.csv the authors' RSLC pipeline expects:
      GuideID,GuideSequence,TargetGene -- no header row.

  Michlits 2017, Nat Methods 14:1191, Supplementary Tables 2/3/4
      via Springer static-content (the article itself is not open access, but these
      supplementary objects are served without a paywall).
      MOESM3 = the guide library (group, 20nt guide, oligo name, 77nt oligo)
      MOESM4 = PCR and sequencing primers
      MOESM5 = sample <-> 6-bp experimental index map, incl. the original BAM filenames.
               Useful only if you obtain those BAMs: the index read is not in the public
               FASTQ. See gcbias/datasets/crispr-umi-michlits.md.

Usage
    python3 get_libraries.py --list
    python3 get_libraries.py schmierer
    python3 get_libraries.py --all
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import os
import re
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path

SOURCES = {
    "schmierer": dict(
        paper="Schmierer 2017, Mol Syst Biol 13:945",
        what="Dataset EV1 -- the RSL-guide library, 23,279 unique guides",
        url="https://www.ebi.ac.uk/europepmc/webservices/rest/PMC5658704/supplementaryFiles",
        member="MSB-13-945-s002.csv",
        out="schmierer_library.csv",
    ),
    "michlits-library": dict(
        paper="Michlits 2017, Nat Methods 14:1191",
        what="Supplementary Table 2 -- the CRISPR_DigDeep guide library",
        url="https://static-content.springer.com/esm/art%3A10.1038%2Fnmeth.4466/"
            "MediaObjects/41592_2017_BFnmeth4466_MOESM3_ESM.xlsx",
        out="michlits_library.xlsx",
    ),
    "michlits-primers": dict(
        paper="Michlits 2017, Nat Methods 14:1191",
        what="Supplementary Table 3 -- PCR and sequencing primers",
        url="https://static-content.springer.com/esm/art%3A10.1038%2Fnmeth.4466/"
            "MediaObjects/41592_2017_BFnmeth4466_MOESM4_ESM.xlsx",
        out="michlits_primers.xlsx",
    ),
    "michlits-samples": dict(
        paper="Michlits 2017, Nat Methods 14:1191",
        what="Supplementary Table 4 -- sample <-> 6-bp index map and BAM filenames",
        url="https://static-content.springer.com/esm/art%3A10.1038%2Fnmeth.4466/"
            "MediaObjects/41592_2017_BFnmeth4466_MOESM5_ESM.xlsx",
        out="michlits_samples.xlsx",
    ),
}


def data_dir() -> Path:
    d = Path(os.environ.get("CHEM_DATA", Path(__file__).resolve().parents[2] / "_data"))
    (d / "libraries").mkdir(parents=True, exist_ok=True)
    return d / "libraries"


def get(url: str) -> bytes:
    if not shutil.which("curl"):
        raise SystemExit("curl not found on PATH")
    # curl rather than urllib: urllib fails behind a TLS-inspecting proxy
    p = subprocess.run(["curl", "-sfL", "-A", "Mozilla/5.0", url], capture_output=True)
    if p.returncode != 0 or not p.stdout:
        raise SystemExit(f"download failed ({p.returncode}): {url}")
    return p.stdout


def xlsx_rows(blob: bytes) -> list[list[str]]:
    """Minimal xlsx reader: shared strings + the first sheet, in row order."""
    z = zipfile.ZipFile(io.BytesIO(blob))
    sst: list[str] = []
    if "xl/sharedStrings.xml" in z.namelist():
        xml = z.read("xl/sharedStrings.xml").decode("utf8", "replace")
        sst = [re.sub(r"<[^>]+>", "", m) for m in
               re.findall(r"<si>(.*?)</si>", xml, re.S)]
    sheets = sorted(n for n in z.namelist() if re.match(r"xl/worksheets/sheet\d+\.xml$", n))
    rows: list[list[str]] = []
    xml = z.read(sheets[0]).decode("utf8", "replace")
    cell_re = re.compile(r"<c\b([^>]*)>(.*?)</c>|<c\b([^>]*)/>", re.S)
    for rm in re.finditer(r"<row[^>]*>(.*?)</row>", xml, re.S):
        cells = []
        for cm in cell_re.finditer(rm.group(1)):
            attrs = cm.group(1) or cm.group(3) or ""
            body = cm.group(2) or ""
            tm = re.search(r'\bt="([^"]+)"', attrs)
            typ = tm.group(1) if tm else ""
            if typ == "inlineStr":
                t = re.findall(r"<t[^>]*>(.*?)</t>", body, re.S)
                cells.append("".join(t))
                continue
            v = re.search(r"<v>(.*?)</v>", body, re.S)
            if v is None:
                cells.append("")
            elif typ == "s":
                i = v.group(1)
                cells.append(sst[int(i)] if i.isdigit() and int(i) < len(sst) else "")
            else:
                cells.append(v.group(1))
        rows.append(cells)
    return rows


def normalise(key: str, path: Path, out: Path) -> dict | None:
    """Write a uniform guide_id/guide/gene TSV that extract.py --whitelist can read."""
    guides: list[tuple[str, str, str]] = []
    if key == "schmierer":
        for r in csv.reader(path.open()):
            if len(r) >= 3 and re.fullmatch(r"[ACGT]{20}", r[1]):
                guides.append((r[0], r[1], r[2]))
    elif key == "michlits-library":
        rows = xlsx_rows(path.read_bytes())
        hdr = [h.strip().lower() for h in rows[0]] if rows else []
        gi = next((i for i, h in enumerate(hdr) if "20nt" in h or h == "guide"), None)
        ni = next((i for i, h in enumerate(hdr) if "oligo name" in h), None)
        gr = next((i for i, h in enumerate(hdr) if h == "group"), None)
        if gi is None:
            return None
        for r in rows[1:]:
            if len(r) > gi and re.fullmatch(r"[ACGT]{20}", r[gi].strip().upper()):
                guides.append((r[ni].strip() if ni is not None and len(r) > ni else "",
                               r[gi].strip().upper(),
                               r[gr].strip() if gr is not None and len(r) > gr else ""))
    else:
        return None
    if not guides:
        return None
    uniq = {g[1] for g in guides}
    with out.open("w") as fh:
        fh.write("guide_id\tguide\tgene\n")
        for a, b, c in guides:
            fh.write(f"{a}\t{b}\t{c}\n")
    return dict(path=str(out), rows=len(guides), unique_guides=len(uniq))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("which", nargs="*", choices=list(SOURCES))
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--list", action="store_true")
    a = ap.parse_args()

    if a.list or (not a.which and not a.all):
        print("published sgRNA libraries and sample maps\n")
        for k, v in SOURCES.items():
            print(f"  {k:18s} {v['paper']}")
            print(f"  {'':18s} {v['what']}")
        print(f"\ndata dir: {data_dir()}   (set $CHEM_DATA to move it)")
        print("\nUse the normalised *.whitelist.tsv with:")
        print("  extract.py --whitelist <file> --correct-guide")
        return 0

    out = data_dir()
    recs = []
    for key in (list(SOURCES) if a.all else a.which):
        src = SOURCES[key]
        dest = out / src["out"]
        print(f"{key}: {src['what']}\n   {src['url']}")
        blob = get(src["url"])
        if src.get("member"):
            z = zipfile.ZipFile(io.BytesIO(blob))
            if src["member"] not in z.namelist():
                raise SystemExit(f"{src['member']} not in the archive; "
                                 f"got {z.namelist()[:5]}")
            blob = z.read(src["member"])
        dest.write_bytes(blob)
        rec = dict(key=key, paper=src["paper"], what=src["what"], url=src["url"],
                   member=src.get("member"), path=str(dest), bytes=len(blob),
                   md5=hashlib.md5(blob).hexdigest())
        wl = normalise(key, dest, out / f"{key}.whitelist.tsv")
        if wl:
            rec["whitelist"] = wl
            print(f"   -> {dest.name}  {len(blob):,} B")
            print(f"   -> {Path(wl['path']).name}  {wl['rows']:,} rows, "
                  f"{wl['unique_guides']:,} unique guides")
        else:
            print(f"   -> {dest.name}  {len(blob):,} B  (no whitelist extracted)")
        recs.append(rec)

    man = out / "MANIFEST.json"
    old = json.loads(man.read_text()) if man.exists() else []
    man.write_text(json.dumps(old + recs, indent=2))
    print(f"\nmanifest: {man}")
    print("Reminder: nothing under this directory is ever committed (see README.md).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
