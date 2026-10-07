#!/usr/bin/env python3
"""Fetch our own CRISPR-MIP data (ArrayExpress E-MTAB-14379) into $CHEM_DATA.

"CRISPR-MIP replaces PCR and reveals GC and oversampling bias in pooled CRISPR screens".
Deposited after the preprint, which still says "will be deposited".

The interesting set for GC work is the **probe-concentration series**: Brunello library,
padlock protocol, non-digested gDNA, 7 ug and 1e6 cells throughout, three probe
concentrations spanning 100x, three replicates each -- everything else held constant. Plus
two digested samples at the same 7 ug, giving the digest contrast at fixed input.

Read structure (from the authors' own subsamp_mip.py):
    R1[20:40] = the 20-nt sgRNA spacer
    R2[0:13]  = the 13-nt UMI

Usage
    python3 get_crisprmip.py --list
    python3 get_crisprmip.py probe --max-reads 4000000
    python3 get_crisprmip.py --all --max-reads 4000000
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

BASE = "https://ftp.sra.ebi.ac.uk/vol1/fastq"
STUDY = "E-MTAB-14379"

# condition -> (run, probe_uM or None, digested)
RUNS = {
    "digested_rep1":   ("ERR13510408", None,  True),
    "digested_rep2":   ("ERR13510409", None,  True),
    "probe0.002_rep1": ("ERR13510410", 0.002, False),
    "probe0.002_rep2": ("ERR13510411", 0.002, False),
    "probe0.002_rep3": ("ERR13510412", 0.002, False),
    "probe0.02_rep1":  ("ERR13510413", 0.02,  False),
    "probe0.02_rep2":  ("ERR13510414", 0.02,  False),
    "probe0.02_rep3":  ("ERR13510415", 0.02,  False),
    "probe0.2_rep1":   ("ERR13510416", 0.2,   False),
    "probe0.2_rep2":   ("ERR13510417", 0.2,   False),
    "probe0.2_rep3":   ("ERR13510418", 0.2,   False),
}
# Matched PCR-vs-padlock pairs on the SAME material. This is the comparison the
# preprint's GC claim is actually about: our conventional PCR readout versus
# CRISPR-MIP, not plasmid versus genomic DNA.
# PCR samples carry no UMI, and their guide is found by ANCHOR (20 bp before
# GTTTTAGAGC) because the Broad P5 primers carry 0-8 nt staggers.
RUNS.update({
    # Brunello-kinome plasmid library, 0.2 ug -- no biology, purest contrast
    "kin_lib_pcr_rep1": ("ERR13510381", None, False), "kin_lib_pcr_rep2": ("ERR13510382", None, False),
    "kin_lib_pcr_rep3": ("ERR13510383", None, False), "kin_lib_pcr_rep4": ("ERR13510384", None, False),
    "kin_lib_mip_rep1": ("ERR13510419", None, False), "kin_lib_mip_rep2": ("ERR13510420", None, False),
    "kin_lib_mip_rep3": ("ERR13510421", None, False),
    # Brunello-kinome normoxia 3d gDNA, 10.1 ug -- realistic screen sample
    "kin_n3d_pcr_rep1": ("ERR13510391", None, False), "kin_n3d_pcr_rep2": ("ERR13510392", None, False),
    "kin_n3d_pcr_rep3": ("ERR13510393", None, False),
    "kin_n3d_mip_rep1": ("ERR13510400", None, False), "kin_n3d_mip_rep2": ("ERR13510401", None, False),
})
RUNS.update({
    # Brunello-kinome normoxia 14d and hypoxia 3d -- the remaining matched pairs
    "kin_n14_pcr_rep1": ("ERR13510385", None, False), "kin_n14_pcr_rep2": ("ERR13510386", None, False),
    "kin_n14_pcr_rep3": ("ERR13510387", None, False),
    "kin_n14_mip_rep1": ("ERR13510402", None, False), "kin_n14_mip_rep2": ("ERR13510403", None, False),
    "kin_n14_mip_rep3": ("ERR13510404", None, False),
    "kin_hyp_pcr_rep1": ("ERR13510405", None, False), "kin_hyp_pcr_rep2": ("ERR13510406", None, False),
    "kin_hyp_pcr_rep3": ("ERR13510407", None, False),
    "kin_hyp_mip_rep1": ("ERR13510397", None, False), "kin_hyp_mip_rep2": ("ERR13510398", None, False),
    "kin_hyp_mip_rep3": ("ERR13510399", None, False),
})
PCR_KEYS = {k for k in RUNS if "_pcr_" in k}

GROUPS = {
    "probe": [k for k in RUNS if k.startswith("probe")],
    "digest": [k for k in RUNS if k.startswith("digested")] + [k for k in RUNS if "0.2_" in k],
    "kinlib": [k for k in RUNS if k.startswith("kin_lib")],
    "kin3d":  [k for k in RUNS if k.startswith("kin_n3d")],
    "kin14d": [k for k in RUNS if k.startswith("kin_n14")],
    "kinhyp": [k for k in RUNS if k.startswith("kin_hyp")],
}


def data_dir() -> Path:
    d = Path(os.environ.get("CHEM_DATA", Path(__file__).resolve().parents[2] / "_data"))
    (d / "crisprmip").mkdir(parents=True, exist_ok=True)
    return d / "crisprmip"


def url_for(run: str, mate: int) -> str:
    """ENA FASTQ path. The subdirectory is NOT the last 3 characters.

    ENA's rule: with a 6-digit accession there is no subdirectory; with more, the
    subdirectory is the digits beyond the first six, zero-padded to 3. So ERR13510410
    (8 digits) -> '10' -> '010', NOT '410'. Getting this wrong 404s silently.
    Cross-checked against the FASTQ_URI column of the E-MTAB-14379 SDRF.
    """
    num = "".join(c for c in run if c.isdigit())
    sub = "" if len(num) <= 6 else f"/{num[6:].zfill(3)}"
    return f"{BASE}/{run[:6]}{sub}/{run}/{run}_{mate}.fastq.gz"


def fetch_one(key: str, mate: int, max_reads: int | None, out: Path) -> dict:
    run, probe, dig = RUNS[key]
    url = url_for(run, mate)
    dest = out / (f"{key}_R{mate}.{max_reads}.fq.gz" if max_reads
                  else f"{key}_R{mate}.fq.gz")
    if dest.exists():
        return dict(key=key, mate=mate, run=run, url=url, path=str(dest), skipped=True)
    if not shutil.which("curl"):
        raise SystemExit("curl not found on PATH")
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
    finally:
        for p in (gun, curl):
            if p.poll() is None:
                p.terminate()
        gun.stdout.close()
    if n == 0:
        tmp.unlink(missing_ok=True)
        raise SystemExit(f"no reads from {url}")
    tmp.rename(dest)
    print(f"    R{mate}: {n:,} reads, {dest.stat().st_size / 1e6:.0f} MB, "
          f"{time.time() - t0:.0f}s")
    return dict(key=key, mate=mate, run=run, url=url, path=str(dest), reads=n,
                md5_of_stream=md5.hexdigest(), capped=max_reads, skipped=False)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("which", nargs="*")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--max-reads", type=int, default=None)
    ap.add_argument("--list", action="store_true")
    a = ap.parse_args()

    if a.list or (not a.which and not a.all):
        print(f"{STUDY} -- CRISPR-MIP, Brunello library, padlock protocol\n")
        print("All samples: 7 ug gDNA, 1e6 cells, NextSeq 500/550 -- "
              "only probe concentration\nand digest state vary.\n")
        print(f"{'key':18s} {'run':14s} {'probe uM':>9s} {'digested':>9s}")
        for k, (run, probe, dig) in RUNS.items():
            print(f"{k:18s} {run:14s} {('-' if probe is None else probe):>9} "
                  f"{str(dig):>9s}")
        print(f"\ngroups: {', '.join(GROUPS)}")
        print(f"data dir: {data_dir()}   (set $CHEM_DATA to move it)")
        print("\nread structure: R1[20:40] = sgRNA, R2[0:13] = UMI")
        return 0

    keys: list[str] = []
    for w in (list(RUNS) if a.all else a.which):
        keys.extend(GROUPS.get(w, [w]))
    keys = sorted(dict.fromkeys(keys))
    out = data_dir()
    recs = []
    for k in keys:
        if k not in RUNS:
            raise SystemExit(f"unknown sample {k}")
        run, probe, dig = RUNS[k]
        print(f"  {k} ({run}, probe={probe} uM, digested={dig})")
        # PCR samples have no UMI, so only R1 is needed
        mates = (1,) if k in PCR_KEYS else (1, 2)
        for mate in mates:
            recs.append(fetch_one(k, mate, a.max_reads, out))
    man = out / "MANIFEST.json"
    old = json.loads(man.read_text()) if man.exists() else []
    man.write_text(json.dumps(old + recs, indent=2))
    print(f"\nmanifest: {man}")
    print("Reminder: nothing under this directory is ever committed (see README.md).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
