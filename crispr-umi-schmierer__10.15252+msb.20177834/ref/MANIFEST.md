# CRISPR-UMI (Schmierer / Michlits) — reference sequences this directory reads

**None of these files is committed.** They are third-party sequence data — Addgene's and the
depositors' deposits and annotation — and not ours to redistribute, so a fresh clone will not
have them. Nothing here is lost by that: every fact derived from them lives in
`tools/crisprumi.py`, `tools/michlits.py` and `tools/selftest.py`, and each group of checks
that reads one of these files is wrapped in `checks.have(...)`, so the self-test **skips** it,
names the file and says where to get it, rather than failing as though the chemistry were
wrong. `tools/build_page.py` and `tools/design_padlock.py` do not skip: a page or a design made
without the real map would be a guess, so they print one line and exit 1.

- **Format:** FASTA (`.fa`, full plasmid, one record) and GenBank (`.gb`, circular, annotated)
- **Papers:** Schmierer *et al.* 2017, [doi:10.15252/msb.20177834](https://doi.org/10.15252/msb.20177834);
  Michlits *et al.* 2017, [doi:10.1038/nmeth.4466](https://doi.org/10.1038/nmeth.4466)
- **Where to get plasmid sequence:** `https://www.addgene.org/<id>/sequences/` — note that
  Addgene put its sequence downloads behind a login in October 2025; the routes that still work
  are written up in
  [`../../lenticrispr-gecko-screen__10.1126+science.1247005/ref/plasmids/MANIFEST.md`](../../lenticrispr-gecko-screen__10.1126+science.1247005/ref/plasmids/MANIFEST.md)

---

## In this directory (`ref/`)

### `addgene-222694_pLenti-UMI.fa`

The lentiviral CRISPR-UMI vector of Michlits *et al.*, 10,694 bp, full sequence.
**Source:** Addgene **#222694** (pLenti-UMI), full sequence, saved as plain FASTA.

Depended on by `tools/selftest.py`:

- *pLenti-UMI as deposited* — length, the built-in Illumina i7 site, the BbsI-repaired P7
  (one C→A vs canonical), presence of the CRISPR-MIP extension arm and absence of its
  ligation arm.
- *the PacI enrichment fragment* — exactly two `TTAATTAA` sites, and the fragment between them
  within 10 bp of the paper's 589 bp (the 2024 deposit gives 596).

### `addgene-222686_pRetro-UMI.fa`

The retroviral sibling vector, 6,991 bp, full sequence.
**Source:** Addgene **#222686** (pRetro-UMI), full sequence, saved as plain FASTA.

Depended on by `tools/selftest.py`: *pRetro-UMI as deposited* — the same five checks as above.

---

## Read from the sibling screen directory

These are **not copies**: the Addgene maps for the pooled-screening vectors are kept once, with
the lentiCRISPR screen notes, in
`../../lenticrispr-gecko-screen__10.1126+science.1247005/ref/plasmids/`. That directory has its
own manifest describing how each file was obtained.

### `addgene-52963_lentiGuide-Puro.gb` — the parent vector

The plasmid Schmierer *et al.* edited into pLenti-Puro-AU-flip-3xBsmBI. Everything
vector-shaped on this page is rebuilt from it by `crisprumi.rebuild_vector()`, so it is the one
file the directory genuinely cannot do without.
**Source:** Addgene **#52963** (lentiGuide-Puro), full sequence.

Depended on by:

- `tools/build_page.py` → `crisprumi.html`: the cloned-vector and PCR1 sizes in *What changes,
  relative to the parent vector*. Missing → the page refuses to build.
- `tools/design_padlock.py`: every candidate ligation arm is scored against this sequence.
  Missing → the design refuses to run.
- `tools/selftest.py`, six check groups: *the parent vector and its PCR1 sizes* (the +39 bp
  edit, 557 → 596 bp PCR1), *capture of both probes on the rebuilt vector*, *probe A's circle,
  library and cycle positions*, *the three nested PCRs and the 288-bp product*, *the same
  readout run on plain lentiGuide-Puro* (249 bp, 39 bp short of the published 288), and
  *universal arm on lentiGuide-Puro*.

A SnapGene map of the same plasmid, `to_debug/lentiGuide-Puro.dna`, is accepted as a fallback
if a local copy happens to be present — it is the same 10,183 bp at a different origin, which
is immaterial because every measurement is made on a circle. `*.dna` is gitignored (the
annotation is SnapGene's), so do not rely on it.

### `addgene-49535_lentiCRISPRv1.gb`, `addgene-52961_lentiCRISPRv2.gb`

**Sources:** Addgene **#49535** (lentiCRISPR v1) and **#52961** (lentiCRISPRv2), full sequences.

Depended on by `tools/selftest.py`, check groups *universal arm on lentiCRISPR v1* and
*universal arm on lentiCRISPRv2*: the whole point of the alternative ligation arm
(`PADLOCK_LIG_ARM_UNIVERSAL`, scaffold positions 46–68) is that it gives one capture at the
same 66-nt gap on these vectors as on Schmierer's, so it needs their maps to be worth claiming.
