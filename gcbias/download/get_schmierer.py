#!/usr/bin/env python3
"""Fetch Schmierer 2017 (ENA PRJEB18436) FASTQ into $CHEM_DATA. Never into the repo.

IMPORTANT: this downloads the *author-submitted* files, not the ENA-generated ones. The
RSL (the lineage UMI) lives in the index field of the FASTQ header, and the generated
`ERR*.fastq.gz` renames every record and throws that away. See
gcbias/datasets/crispr-umi-schmierer.md.

Usage
    python3 get_schmierer.py --list
    python3 get_schmierer.py input --max-reads 2000000
    python3 get_schmierer.py --all --max-reads 500000

`--max-reads N` streams only as much of the gzip member as needed, so a 49 GB file costs a
few hundred MB. Omit it to pull the whole thing (127 GB for --all; you probably do not
want that).
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import sys
import shutil
import subprocess
import time
from pathlib import Path

BASE = "https://ftp.sra.ebi.ac.uk/vol1/run"
STUDY = "PRJEB18436"

# run -> (author filename, sample title, reads, approx gz bytes)
RUNS = {
    "input":  ("ERR2114695", "input.fq.gz",     "RKO_NTU_INPUT",          2_503_174_977, 49_169_994_748),
    "r1_d4":  ("ERR2114696", "ntu_r1_d4.fq.gz", "RKO_NTU_Replicate1_D4",    660_790_803, 14_157_292_526),
    "r1_d28": ("ERR2114697", "ntu_r1_d28.fq.gz", "RKO_NTU_Replicate1_D28", 1_189_493_924, 23_691_745_671),
    "r2_d4":  ("ERR2114698", "ntu_r2_d4.fq.gz", "RKO_NTU_Replicate2_D4",    794_719_413, 16_808_250_341),
    "r2_d28": ("ERR2114699", "ntu_r2_d28.fq.gz", "RKO_NTU_Replicate2_D28", 1_145_941_439, 22_928_184_615),
}


def data_dir() -> Path:
    d = Path(os.environ.get("CHEM_DATA", Path(__file__).resolve().parents[2] / "_data"))
    (d / "schmierer").mkdir(parents=True, exist_ok=True)
    return d / "schmierer"


def url_for(key: str) -> str:
    run, fname, *_ = RUNS[key]
    return f"{BASE}/{run[:6]}/{run}/{fname}"


def fetch_capped(key: str, max_reads: int | None, out: Path) -> dict:
    """Stream the gzip member via curl and stop after max_reads records.

    curl rather than urllib on purpose: urllib honours the system trust store, which fails
    behind a TLS-inspecting proxy ("self-signed certificate in certificate chain"), and
    this script is meant to work without the caller debugging their certificates.
    """
    run, fname, title, total_reads, gz_bytes = RUNS[key]
    url = url_for(key)
    dest = out / (f"{key}.{max_reads}.fq.gz" if max_reads else f"{key}.fq.gz")
    if dest.exists():
        print(f"  {key}: already at {dest} -- skipping")
        return dict(key=key, run=run, url=url, path=str(dest), skipped=True)
    if not shutil.which("curl"):
        raise SystemExit("curl not found on PATH")

    print(f"  {key}: {title}\n       {url}\n       full {gz_bytes / 1e9:.1f} GB / "
          f"{total_reads:,} reads; "
          + (f"first {max_reads:,} reads" if max_reads else "ALL of it"))
    t0 = time.time()
    tmp = dest.with_suffix(".part")
    md5 = hashlib.md5()
    n = 0
    curl = subprocess.Popen(["curl", "-sfL", url], stdout=subprocess.PIPE)
    gunzip = subprocess.Popen(["gzip", "-dc"], stdin=curl.stdout,
                              stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    curl.stdout.close()
    try:
        with gzip.open(tmp, "wb", compresslevel=1) as fh:
            while max_reads is None or n < max_reads:
                rec = [gunzip.stdout.readline() for _ in range(4)]
                if not rec[0]:
                    break
                for ln in rec:
                    fh.write(ln)
                    md5.update(ln)
                n += 1
                if n % 2_000_000 == 0:
                    print(f"       {n:,} reads ({time.time() - t0:.0f}s)", flush=True)
    finally:
        # we deliberately stop early, so SIGPIPE on the producers is expected
        for p in (gunzip, curl):
            if p.poll() is None:
                p.terminate()
        gunzip.stdout.close()
    if n == 0:
        tmp.unlink(missing_ok=True)
        raise SystemExit(f"no reads retrieved from {url} -- check connectivity")
    tmp.rename(dest)
    print(f"       -> {dest}  {n:,} reads, {dest.stat().st_size / 1e6:.1f} MB, "
          f"{time.time() - t0:.0f}s")
    return dict(key=key, run=run, url=url, path=str(dest), reads=n,
                md5_of_stream=md5.hexdigest(), seconds=round(time.time() - t0, 1),
                capped=max_reads, skipped=False)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("samples", nargs="*", choices=list(RUNS) + [], help="which samples")
    ap.add_argument("--all", action="store_true", help="every sample")
    ap.add_argument("--max-reads", type=int, default=None,
                    help="stop after this many reads (strongly recommended)")
    ap.add_argument("--list", action="store_true", help="show what is available and exit")
    a = ap.parse_args()

    if a.list or (not a.samples and not a.all):
        print(f"Schmierer 2017, ENA {STUDY} -- author-submitted FASTQ (RSL in the header)\n")
        print(f"{'key':8s} {'run':12s} {'sample':26s} {'reads':>16s} {'gz':>9s}")
        for k, (run, fn, title, nr, gb) in RUNS.items():
            print(f"{k:8s} {run:12s} {title:26s} {nr:>16,} {gb / 1e9:>7.1f} GB")
        print(f"\ndata dir: {data_dir()}   (set $CHEM_DATA to move it)")
        print("the plasmid INPUT sample is the one to start with -- no selection, so a GC")
        print("trend in it is amplification rather than biology")
        return 0

    keys = list(RUNS) if a.all else a.samples
    out = data_dir()
    if a.max_reads is None:
        print("WARNING: no --max-reads; this will download the full files.\n")
    recs = [fetch_capped(k, a.max_reads, out) for k in keys]

    man = out / "MANIFEST.json"
    old = json.loads(man.read_text()) if man.exists() else []
    man.write_text(json.dumps(old + recs, indent=2))
    print(f"\nmanifest: {man}")
    print("Reminder: nothing under this directory is ever committed (see README.md).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
