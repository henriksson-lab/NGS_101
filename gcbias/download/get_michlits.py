#!/usr/bin/env python3
"""Fetch Michlits 2017 (SRA PRJNA383356) into $CHEM_DATA.

READ THIS FIRST: the public FASTQ for this study does NOT contain the UMI. The SRA runs
report nreads=1 -- one 50-nt read, of which only the first 20 bases are the sgRNA and the
rest is run-off. The 10-nt UMI was index read 1 and the 6-nt experimental index was index
read 2, and neither was deposited. See gcbias/datasets/crispr-umi-michlits.md.

So this script gets you guide-level read counts only. That is still useful -- the plasmid
library run is a no-selection baseline for a GC-vs-abundance check -- but there is no
molecule counter, and the pooled arms inside a run cannot be separated.

The original per-lane BAMs probably do carry the index reads, but they are
`access_type="Use Cloud Data Delivery"` (requester-pays); a plain HTTPS request to the
sra-pub-src bucket returns 403. --show-bams prints them for whoever pursues that route.

Usage
    python3 get_michlits.py --list
    python3 get_michlits.py library --max-reads 25000000
    python3 get_michlits.py --show-bams
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

STUDY = "PRJNA383356"
BASE = "https://ftp.sra.ebi.ac.uk/vol1/fastq"

RUNS = {
    "library":  ("SRR5559298", "CRISPR_DigDeep_KO_Library: Plasmid Library Sequencing",
                 625_912_869),
    "pilot":    ("SRR5559297", "pilot screen, mESC, etoposide hypersensitivity",
                 292_762_401),
    "etoposide": ("SRR5559299", "negative selection, mESC, etoposide (control+treated POOLED)",
                  563_804_278),
    "ipsc":     ("SRR5559300", "positive selection, MEF->iPSC reprogramming roadblocks",
                 73_462_865),
}

BAMS = {
    "SRR5559298": [("Plasmid_NDS_C89KRANXX_6_20160127B_20160128.bam", 7.36),
                   ("Plasmid_NDS_C9PLEANXX_1_20161111B_20161114.bam", 12.58)],
    "SRR5559297": [("Pilot_screen_C9B8KANXX_7_20160406B_20160408.bam", 2.07),
                   ("Pilot_screen_C8LABANXX_2_20160414B_20160416.bam", 7.08)],
    "SRR5559299": [("30k_Etoposide_CAF25ANXX_1_07-02-2017.bam", 6.89),
                   ("30k_Etoposide_CAF25ANXX_2_07-02-2017.bam", 3.93),
                   ("30k_Etoposide_CAR6CANXX_7_20170310B_20170311.bam", 4.84),
                   ("30k_Etoposide_CAR6CANXX_8_20170310B_20170311.bam", 4.23)],
    "SRR5559300": [("30k_iPSC_H5F3TBCXX_1_20160927B_20160928.bam", 2.90)],
}


def data_dir() -> Path:
    d = Path(os.environ.get("CHEM_DATA", Path(__file__).resolve().parents[2] / "_data"))
    (d / "michlits").mkdir(parents=True, exist_ok=True)
    return d / "michlits"


def url_for(run: str) -> str:
    return f"{BASE}/{run[:6]}/{run[-3:].zfill(3)}/{run}/{run}.fastq.gz"


def fetch(key: str, max_reads: int | None, out: Path) -> dict:
    run, title, total = RUNS[key]
    url = url_for(run)
    dest = out / (f"{key}.{max_reads}.fq.gz" if max_reads else f"{key}.fq.gz")
    if dest.exists():
        print(f"  {key}: already at {dest} -- skipping")
        return dict(key=key, run=run, url=url, path=str(dest), skipped=True)
    if not shutil.which("curl"):
        raise SystemExit("curl not found on PATH")
    print(f"  {key}: {title}\n       {url}\n       {total:,} reads total; "
          + (f"first {max_reads:,}" if max_reads else "ALL of it"))
    t0, md5, n = time.time(), hashlib.md5(), 0
    tmp = dest.with_suffix(".part")
    curl = subprocess.Popen(["curl", "-sfL", url], stdout=subprocess.PIPE)
    gun = subprocess.Popen(["gzip", "-dc"], stdin=curl.stdout, stdout=subprocess.PIPE,
                           stderr=subprocess.DEVNULL)
    curl.stdout.close()
    try:
        with gzip.open(tmp, "wb", compresslevel=1) as fh:
            while max_reads is None or n < max_reads:
                rec = [gun.stdout.readline() for _ in range(4)]
                if not rec[0]:
                    break
                for ln in rec:
                    fh.write(ln)
                    md5.update(ln)
                n += 1
                if n % 2_000_000 == 0:
                    print(f"       {n:,} reads ({time.time() - t0:.0f}s)", flush=True)
    finally:
        for p in (gun, curl):
            if p.poll() is None:
                p.terminate()
        gun.stdout.close()
    if n == 0:
        tmp.unlink(missing_ok=True)
        raise SystemExit(f"no reads retrieved from {url}")
    tmp.rename(dest)
    print(f"       -> {dest}  {n:,} reads, {dest.stat().st_size / 1e6:.1f} MB")
    return dict(key=key, run=run, url=url, path=str(dest), reads=n,
                md5_of_stream=md5.hexdigest(), capped=max_reads, skipped=False)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("samples", nargs="*", choices=list(RUNS))
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--max-reads", type=int, default=None)
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--show-bams", action="store_true",
                    help="print the original BAMs that may hold the index reads")
    a = ap.parse_args()

    if a.show_bams:
        print(f"Original per-lane BAMs in SRA for {STUDY}.")
        print("These are the only place the 10-nt UMI may survive. NOT anonymously")
        print("downloadable: access_type='Use Cloud Data Delivery' (requester-pays);")
        print("a direct https GET to sra-pub-src-* returns 403.\n")
        tot = 0.0
        for run, items in BAMS.items():
            for fn, gb in items:
                print(f"  {run}  {gb:>6.2f} GB  {fn}")
                tot += gb
        print(f"\n  total {tot:.1f} GB across {sum(len(v) for v in BAMS.values())} files")
        print("\n  Before paying to move all of it, verify on the smallest file")
        print("  (30k_iPSC_H5F3TBCXX_1, 2.90 GB) that the index reads are actually")
        print("  present as BAM tags. Asking the authors is probably faster.")
        return 0

    if a.list or (not a.samples and not a.all):
        print(f"Michlits 2017, SRA {STUDY}\n")
        print("*** the public FASTQ has NO UMI and NO experimental index ***")
        print("*** guide-level counts only; see --show-bams              ***\n")
        print(f"{'key':11s} {'run':12s} {'reads':>16s}  description")
        for k, (run, title, n) in RUNS.items():
            print(f"{k:11s} {run:12s} {n:>16,}  {title}")
        print(f"\ndata dir: {data_dir()}   (set $CHEM_DATA to move it)")
        print("'library' is the plasmid run -- the no-selection baseline.")
        return 0

    keys = list(RUNS) if a.all else a.samples
    out = data_dir()
    if a.max_reads is None:
        print("WARNING: no --max-reads; full files are tens of GB.\n")
    recs = [fetch(k, a.max_reads, out) for k in keys]
    man = out / "MANIFEST.json"
    old = json.loads(man.read_text()) if man.exists() else []
    man.write_text(json.dumps(old + recs, indent=2))
    print(f"\nmanifest: {man}")
    print("Reminder: nothing under this directory is ever committed (see README.md).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
