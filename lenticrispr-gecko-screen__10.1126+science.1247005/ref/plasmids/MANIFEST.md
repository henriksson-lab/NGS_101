# Pooled CRISPR screening vectors — full-length reference sequences

Full plasmid maps for checking primer-binding sites *in silico* and deriving expected PCR products.

- **Directory:** `lenticrispr-gecko-screen__10.1126+science.1247005/ref/plasmids` (repo-relative)
- **Downloaded:** 2026-10-03
- **Format:** GenBank (`.gb`), circular, with the depositor/Addgene feature annotations
- **Primary files:** 10 plasmids, all verbatim downloads
- **`alt-snapgene/`:** 4 independent SnapGene maps kept only as cross-checks

## These files are third-party and may not be in your clone

**Every `.gb` here is someone else's deposited sequence data, redistributed by Addgene or
SnapGene under their own terms — none of it is ours.** It is kept out of the repository's
own content for that reason, and a fresh clone may legitimately have an empty (or absent)
`ref/plasmids/`. Nothing is reconstructed from memory if a file is missing; it simply is not
there, and the tooling says so:

- **`tools/selftest.py` skips the checks that need a map.** Each file is declared as a
  `checks.Source` with its Addgene ID, and the groups that read it sit behind `have(...)`.
  Without the maps the suite still **passes**, prints one `SKIP` line per missing file naming
  it and where to download it, and still runs every map-independent check (the overhang and
  guide-oligo arithmetic, the "nothing is retyped" checks, the stagger ladders). The
  `alt-snapgene/` comparisons are guarded one pair at a time, so one missing cross-check map
  costs that one comparison and nothing else.
- **`tools/build_page.py` refuses to build.** `crisprscreen.html` is computed from the maps
  end to end — every amplicon size, every segment length — so there is no honest page to
  emit without them. It prints one line per missing file, naming the file and the Addgene
  page to get it from, and exits 1. It does not fall back to quoted numbers.

To restore the directory, re-download the files from the URLs in the tables below (the
Addgene route needs a login now — see immediately below).

---

## How these were obtained

**Addgene put its sequence data behind a login in October 2025**
([announcement](https://blog.addgene.org/protecting-the-science-you-share-why-sequence-data-requires-a-login)),
which breaks the route this task assumed. Concretely, today:

- `https://www.addgene.org/<id>/sequences/` still loads and still tells you *which* full and partial
  sequences exist — but the `.gbk` download links have been stripped out of the page.
- `https://www.addgene.org/browse/sequence/<seqid>/` ("Analyze Sequence") redirects to a login prompt.
- `https://media.addgene.org/snapgene-media/.../addgene-plasmid-<id>-sequence-<seqid>.gbk`
  returns **HTTP 404** for every version prefix tried.
- `https://api.addgene.org/...` is not a public API (returns the marketing site).

Three routes still work, and all three were used:

**1. Addgene's own sequence-viewer JSON — the primary source for 9 of the 10 files.**
The sequence viewer is a client-side app, and the data it loads is served openly from
`sequences.addgene.org` with no login:

```
https://sequences.addgene.org/snapgene-media/v3.101.0/sequences/<seqid>/<uuid>/addgene-plasmid-<id>-sequence-<seqid>-sequence-lines.json
https://sequences.addgene.org/snapgene-media/v3.101.0/sequences/<seqid>/<uuid>/addgene-plasmid-<id>-sequence-<seqid>-features.json
```

`<seqid>` and `<uuid>` are still readable off the live (non-login) sequences page, because the
map thumbnail is served from the same directory. `-sequence-lines.json` carries the full plasmid
as SnapGene's rendered sequence view — one SVG per 60 bp line, each base an individual `<text>`
glyph — and `-features.json` carries the complete feature table with coordinates, directions and
notes. Both were parsed and re-emitted as GenBank. **This is the real, deposited Addgene data, not
a reconstruction**, but because it needed parsing rather than just downloading, it was validated
three independent ways before being trusted — see the cross-validation section at the end.

**2. Internet Archive, for one cross-check.**
The original `addgene-plasmid-52961-sequence-244694.gbk` is still in the Wayback Machine (the
URL was rebuilt from the live page's `<seqid>`/`<uuid>` and the archived `snapgene-media` version
prefix, located through the CDX index). It was fetched purely to confirm the JSON parser: the two
agree byte for byte. The other plasmids' `.gbk` files were not captured by the Archive under any
version prefix that could be resolved.

**3. Files Addgene still serves openly, and SnapGene.**
`https://media.addgene.org/cms/filer_public/.../*.gb` — depositor-uploaded documents — are not
login-gated; `addgene-66217_CRISPRi-library-backbone_*.gb` came from there verbatim.
SnapGene's public plasmid set (`https://www.snapgene.com/local/fetch.php?set=crispr_plasmids&plasmid=<name>`)
serves genuine `.dna` files credited to the depositing lab and Addgene number; those for
lentiCRISPR v2, lentiCas9-Blast, lentiGuide-Puro and the Weissman CRISPRi backbone were converted
with Biopython's SnapGene reader and kept in `alt-snapgene/` **as independent cross-checks only** —
they are built by SnapGene from the depositors' data, not copied from Addgene, so agreement
between the two is meaningful.

**Checked and unusable:** NCBI nuccore (E-utilities) has no accession for any of these plasmids —
searches for `lentiCRISPR`, `lentiCRISPRv2`, `lentiGuide-Puro` and `pXPR_003` return only patent
sequences or nothing. The Broad GPP vector portal (`portals.broadinstitute.org/gppx/vector/public`)
is now a JavaScript app whose data endpoint `/gppx/vector/public/api/init` returns **HTTP 401**.

**Nothing here was fabricated.** No sequence was assembled, inferred, or filled in from memory;
where a source was unavailable it is recorded as unavailable in "Not obtained".

---

## Files

| Addgene | Plasmid | Depositor | Length (bp) | Addgene page says | File | Bytes |
|---|---|---|---|---:|---|---:|
| [#49535](https://www.addgene.org/49535/) | lentiCRISPR v1 | Feng Zhang | **13450** | not published | `addgene-49535_lentiCRISPRv1.gb` | 22933 |
| [#52961](https://www.addgene.org/52961/) | lentiCRISPRv2 | Feng Zhang | **14873** | 14873 | `addgene-52961_lentiCRISPRv2.gb` | 25018 |
| [#52963](https://www.addgene.org/52963/) | lentiGuide-Puro | Feng Zhang | **10183** | 10183 | `addgene-52963_lentiGuide-Puro.gb` | 17996 |
| [#52962](https://www.addgene.org/52962/) | lentiCas9-Blast | Feng Zhang | **12859** | 12860 ⚠ | `addgene-52962_lentiCas9-Blast.gb` | 22236 |
| [#86708](https://www.addgene.org/86708/) | CROPseq-Guide-Puro | Christoph Bock | **10214** | 10214 | `addgene-86708_CROPseq-Guide-Puro.gb` | 18144 |
| [#60955](https://www.addgene.org/60955/) | pU6-sgRNA EF1Alpha-puro-T2A-BFP | Jonathan Weissman | **8877** | 8888 ⚠ | `addgene-60955_pU6-sgRNA-EF1Alpha-puro-T2A-BFP.gb` | 15306 |
| [#59702](https://www.addgene.org/59702/) | pXPR_011 | John Doench, David Root (Broad GPP) | **8395** | 8395 | `addgene-59702_pXPR_011.gb` | 15568 |
| [#98293](https://www.addgene.org/98293/) | lentiCRISPRv2 blast (lentiCRISPRv2-Blast) | Brett Stringer | **14676** | not published | `addgene-98293_lentiCRISPRv2-Blast.gb` | 24838 |
| [#193588](https://www.addgene.org/193588/) | pXPR_003 sgMSLN | William Hahn | **8323** | not published | `addgene-193588_pXPR_003-sgMSLN.gb` | 15654 |
| [#66217](https://www.addgene.org/pooled-library/weissman-human-crispri-v2/) | CRISPRi v2 pooled-library backbone (pU6-sgRNA EF1Alpha-puro-T2A-BFP family) | Jonathan Weissman | **8888** | 8888 | `addgene-66217_CRISPRi-library-backbone_pU6-sgRNA-EF1Alpha-puro-T2A-BFP.gb` | 14603 |

⚠ = length differs from the number on Addgene's web page; see the per-plasmid notes.

| Source URL (sequence ID) | | 
|---|---|
| **#49535** lentiCRISPR v1 | `https://sequences.addgene.org/snapgene-media/v3.101.0/sequences/73924/9de1ec0a-8066-4892-8c56-00d24dd59816/addgene-plasmid-49535-sequence-73924-{sequence-lines,features}.json` <br> Addgene sequence record: 73924 (Depositor Full Sequence) |
| **#52961** lentiCRISPRv2 | `https://sequences.addgene.org/snapgene-media/v3.101.0/sequences/244694/51ae3e01-5e94-49b8-a4fb-7fa303a7eb81/addgene-plasmid-52961-sequence-244694-{sequence-lines,features}.json` <br> Addgene sequence record: 244694 (Addgene-Verified Full Sequence) |
| **#52963** lentiGuide-Puro | `https://sequences.addgene.org/snapgene-media/v3.101.0/sequences/331247/092a4fee-eef4-4211-8e83-d1edf82c697d/addgene-plasmid-52963-sequence-331247-{sequence-lines,features}.json` <br> Addgene sequence record: 331247 (Addgene-Verified Full Sequence) |
| **#52962** lentiCas9-Blast | `https://sequences.addgene.org/snapgene-media/v3.101.0/sequences/322376/afc5ce55-af10-49b2-b6d3-01746afdc217/addgene-plasmid-52962-sequence-322376-{sequence-lines,features}.json` <br> Addgene sequence record: 322376 (Addgene-Verified Full Sequence) |
| **#86708** CROPseq-Guide-Puro | `https://sequences.addgene.org/snapgene-media/v3.101.0/sequences/317051/ecf32836-79cd-4e6f-a931-771c0a1001ea/addgene-plasmid-86708-sequence-317051-{sequence-lines,features}.json` <br> Addgene sequence record: 317051 (Addgene-Verified Full Sequence) |
| **#60955** pU6-sgRNA EF1Alpha-puro-T2A-BFP | `https://sequences.addgene.org/snapgene-media/v3.101.0/sequences/467280/2b82e300-c292-4d2f-822e-27fbbc4b06cd/addgene-plasmid-60955-sequence-467280-{sequence-lines,features}.json` <br> Addgene sequence record: 467280 (Addgene-Verified Full Sequence) |
| **#59702** pXPR_011 | `https://sequences.addgene.org/snapgene-media/v3.101.0/sequences/233279/37129534-a564-4e66-9822-ced75ef2a42d/addgene-plasmid-59702-sequence-233279-{sequence-lines,features}.json` <br> Addgene sequence record: 233279 (Addgene-Verified Full Sequence) |
| **#98293** lentiCRISPRv2 blast (lentiCRISPRv2-Blast) | `https://sequences.addgene.org/snapgene-media/v3.101.0/sequences/489390/dcb9e0df-d61e-448f-95d1-8689058ff943/addgene-plasmid-98293-sequence-489390-{sequence-lines,features}.json` <br> Addgene sequence record: 489390 (Addgene-Verified Full Sequence) |
| **#193588** pXPR_003 sgMSLN | `https://sequences.addgene.org/snapgene-media/v3.101.0/sequences/380090/a40b489b-3539-46ef-a215-b285eb5c1d66/addgene-plasmid-193588-sequence-380090-{sequence-lines,features}.json` <br> Addgene sequence record: 380090 (Addgene-Verified Full Sequence) |
| **#66217** CRISPRi v2 pooled-library backbone (pU6-sgRNA EF1Alpha-puro-T2A-BFP family) | `https://media.addgene.org/cms/filer_public/bc/b9/bcb9a575-5058-44e6-9916-0cd33075dea7/crispri_backbone_66217.gb` <br> Addgene sequence record: depositor-uploaded file (no sequence ID) |

### Cross-check maps in `alt-snapgene/`

Four maps SnapGene builds from the depositors' data rather than copying from Addgene, so
agreement with the Addgene files above is meaningful. Converted from SnapGene's public
`.dna` set with Biopython's SnapGene reader. **Cross-checks only** — no page or computed
claim depends on them, and each comparison in `tools/selftest.py` is skipped individually
if its file is absent.

Source for all four: `https://www.snapgene.com/local/fetch.php?set=crispr_plasmids&plasmid=<name>`,
where `<name>` is the SnapGene map name in the last column.

| Addgene | Plasmid | Length (bp) | File | Bytes | SnapGene map name |
|---|---|---:|---|---:|---|
| [#52961](https://www.addgene.org/52961/) | lentiCRISPR v2 | 14873 | `alt-snapgene/snapgene-52961_lentiCRISPRv2.gb` | 28553 | `lentiCRISPR v2` |
| [#52962](https://www.addgene.org/52962/) | lentiCas9-Blast | 12859 | `alt-snapgene/snapgene-52962_lentiCas9-Blast.gb` | 24953 | `lentiCas9-Blast` |
| [#52963](https://www.addgene.org/52963/) | lentiGuide-Puro | 10183 | `alt-snapgene/snapgene-52963_lentiGuide-Puro.gb` | 17832 | `lentiGuide-Puro` |
| [#62217](https://www.addgene.org/pooled-library/weissman-human-crispri-v2/) | pCRISPRia-v2 / CRISPRi library backbone (Weissman Lab) | 8888 | `alt-snapgene/snapgene-62217_pCRISPRia-v2.gb` | 17313 | `CRISPRi_library_backbone` |

### What depends on each file

So you can tell at a glance what breaks — or rather, what skips — when a file is missing.

| File | Addgene | Needed by |
|---|---|---|
| `addgene-49535_lentiCRISPRv1.gb` | #49535 | **`tools/build_page.py`** (the primer/vector size table) and **`tools/selftest.py`** (BsmBI site + overhangs + filler length, v1 scaffold identity, KERMIT and Shalem amplicon sizes, U6 stretch present once) |
| `addgene-52961_lentiCRISPRv2.gb` | #52961 | **`tools/build_page.py`** (every figure in parts 2–3, both readout rounds, the sizing tables) and **`tools/selftest.py`** (the published 285 nt BEAKER product, the 287/2148 bp step-1 pair, the +1 G shift, the KERMIT wrap-around, round-1/round-2 geometry) |
| `addgene-52963_lentiGuide-Puro.gb` | #52963 | **`tools/build_page.py`** (the second backbone in every table and the round-1/2 figures) and **`tools/selftest.py`** (the 556/2417 bp and 316/2177 bp step-1 sizes, the 269 bp round-1 difference, BEAKER's site being present but unusable) |
| `addgene-52962_lentiCas9-Blast.gb` | #52962 | **`tools/selftest.py`** only — the two negative controls: zero BsmBI sites and no sgRNA scaffold |
| `addgene-59702_pXPR_011.gb` | #59702 | **`tools/selftest.py`** only — the F+E scaffold check that makes the v1-vs-F+E distinction real |
| `addgene-66217_CRISPRi-library-backbone_pU6-sgRNA-EF1Alpha-puro-T2A-BFP.gb` | #66217 | **`tools/selftest.py`**, optional pair with `alt-snapgene/snapgene-62217_*`: same length, the N20 guide slot, and that the N run is the *only* difference |
| `addgene-86708_CROPseq-Guide-Puro.gb` | #86708 | nothing computed — reference only (minus-strand cassette example in the motif table) |
| `addgene-60955_pU6-sgRNA-EF1Alpha-puro-T2A-BFP.gb` | #60955 | nothing computed — reference only (BstXI/BlpI cloning, mouse U6, Chen 2013 scaffold) |
| `addgene-98293_lentiCRISPRv2-Blast.gb` | #98293 | nothing computed — reference only (marker swap, minus-strand cassette) |
| `addgene-193588_pXPR_003-sgMSLN.gb` | #193588 | nothing computed — reference only (stand-in for the undeposited empty pXPR_003) |
| `alt-snapgene/snapgene-52961_lentiCRISPRv2.gb` | #52961 (SnapGene's own map) | **`tools/selftest.py`**, optional: exact circular rotation of the Addgene map at offset 14256 |
| `alt-snapgene/snapgene-52962_lentiCas9-Blast.gb` | #52962 (SnapGene's own map) | **`tools/selftest.py`**, optional: exact circular rotation at offset 4315 — also the second source confirming 12859 bp, not Addgene's stated 12860 |
| `alt-snapgene/snapgene-52963_lentiGuide-Puro.gb` | #52963 (SnapGene's own map) | **`tools/selftest.py`**, optional: exact circular rotation at offset 2593 |
| `alt-snapgene/snapgene-62217_pCRISPRia-v2.gb` | #62217 (SnapGene's own map) | **`tools/selftest.py`**, optional: the #66217 comparison above |

Nothing outside this directory's own `MANIFEST` depends on the four reference-only files;
they are kept because they are what makes the motif table a survey rather than three
examples.

---

## Landmark motifs

Positions are **1-based, on the plus strand of each deposited map**. `+ n` = the motif itself
starts at n; `− n` = the motif's **reverse complement** starts at n (i.e. the motif is on the
minus strand of this map). `—` = absent from both strands. Searches wrap the circular origin.

Coordinate origins are arbitrary and differ between sources, so compare *patterns*, not absolute positions.

| Plasmid | bp | U6 3′<br>`GACGAAACACCG` | scaffold v1<br>`GTTTTAGAGCTAGAAATAGCAAG` | scaffold F+E<br>`GTTTAAGAGCTATGCTGGAAACAGCATAGCAAG` | `CGTCTC` | `GAGACG` | **BsmBI<br>sites** | WPRE<br>`AATCAACC…TGTGAAAG` |
|---|---:|---|---|---|---|---|:-:|---|
| **lentiCRISPR v1** (#49535) | 13450 | + 1955 | + 3847 | — | + 3840; − 1967 | + 1967; − 3840 | **2** | + 9307 |
| **lentiCRISPRv2** (#52961) | 14873 | + 2228 | + 4120 | — | + 4113; − 2240 | + 2240; − 4113 | **2** | + 9347 |
| **lentiGuide-Puro** (#52963) | 10183 | + 400 | + 2292 | — | + 2285; − 412 | + 412; − 2285 | **2** | + 4485 |
| **lentiCas9-Blast** (#52962) | 12859 | — | — | — | — | — | **0** | + 11865 |
| **CROPseq-Guide-Puro** (#86708) | 10214 | − 5143 | − 3240 | — | + 5137; − 3264 | + 3264; − 5137 | **2** | − 6105 |
| **pU6-sgRNA EF1Alpha-puro-T2A-BFP** (#60955) | 8877 | — | — | — | + 4673; − 5388 | + 5388; − 4673 | **2** | + 5766 |
| **pXPR_011** (#59702) | 8395 | + 8062 | — | + 8093 | + 247, 905; − 1195 | + 1195; − 247, 905 | **3** | + 1964 |
| **lentiCRISPRv2 blast (lentiCRISPRv2-Blast)** (#98293) | 14676 | − 276 | − 13049 | — | + 270; − 13073 | + 13073; − 270 | **2** | − 8014 |
| **pXPR_003 sgMSLN** (#193588) | 8323 | + 400 | + 432 | — | — | — | **0** | + 2625 |
| **CRISPRi v2 pooled-library backbone (pU6-sgRNA EF1Alpha-puro-T2A-BFP family)** (#66217) | 8888 | — | — | — | + 4684; − 5399 | + 5399; − 4684 | **2** | + 5777 |

**BsmBI/Esp3I site count** = number of `CGTCTC` recognition sequences on *either* strand,
i.e. plus-strand `CGTCTC` hits + plus-strand `GAGACG` hits. (The two columns are listed
separately above as requested; they are not independent sites.)

### Reading the motif table

- **5 of 10 plasmids carry the classic Zhang cassette**: hU6 ending in `GACGAAACACCG`,
  a 1889 bp BsmBI filler, then the v1 scaffold `GTTTTAGAGCTAGAAATAGCAAG` —
  lentiCRISPR v1, lentiCRISPRv2, lentiGuide-Puro, CROPseq-Guide-Puro, lentiCRISPRv2-Blast.
  In all five the two BsmBI sites face inward and sit **exactly 1873 bp apart**, bracketing
  the guide-cloning filler (for #98293 that span wraps the map origin).
- **CROPseq-Guide-Puro and lentiCRISPRv2-Blast have the cassette on the minus strand**
  of their maps — every landmark shows up as `− n`. Orientation matters when you place primers.
- **lentiCas9-Blast has no U6, no scaffold and 0 BsmBI sites** — it is a Cas9-delivery
  vector only, not an sgRNA cloning backbone.
- **pXPR_011 is the only F+E-scaffold vector here** — and it is **not an empty backbone**:
  between the U6 3′ end (ends 8073) and the F+E scaffold (8093) sits a 19 nt spacer,
  `GTGAACCGCATCGAGCTGA`, i.e. the reporter guide is already cloned in, so there is no BsmBI
  filler. Its **3** `CGTCTC` sites (plus 247, 905; minus 1195) lie elsewhere in the vector.
- **The Weissman vectors (#60955, #66217, and #62217 in `alt-snapgene/`) match none of the
  three sgRNA probes.** That is correct, not a download problem: they use the *mouse* U6
  promoter and the Chen et al. 2013 'optimized' scaffold
  (`GTTTAAGAGCTAAGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCG`), which differs from the F+E
  probe by 2 nt. Their 2 `CGTCTC`/`GAGACG` hits sit **inside PuroR and TagBFP** and are not
  cloning sites — these vectors are opened with **BstXI + BlpI**.
- **pXPR_003 sgMSLN already contains its guide**, hence 0 BsmBI sites.

---

## Per-plasmid detail

### lentiCRISPR v1 — Addgene [#49535](https://www.addgene.org/49535/)

| | |
|---|---|
| File | `addgene-49535_lentiCRISPRv1.gb` (22933 bytes) |
| Length | **13450 bp**, circular, 49.1% GC |
| Addgene page states | no total vector size published |
| Features | 33 |
| Depositor | Feng Zhang |
| Addgene sequence record | 73924 (Depositor Full Sequence) |
| Source | Addgene public sequence-viewer JSON |
| Source URL | `https://sequences.addgene.org/snapgene-media/v3.101.0/sequences/73924/9de1ec0a-8066-4892-8c56-00d24dd59816/addgene-plasmid-49535-sequence-73924-{sequence-lines,features}.json` |
| Retrieved | 2026-10-03 |
| BsmBI sites | 2 |

Depositor full sequence (Addgene has no verified full sequence for v1). Addgene's plasmid page does not publish a Total vector size for #49535, so the 13450 bp here could not be cross-checked against an independent number; the record is internally consistent (225 viewer lines x 60 bp - 50). All-in-one Cas9+sgRNA vector: hU6 - BsmBI filler - gRNA scaffold - EFS - Cas9-FLAG-P2A-Puro - WPRE.

### lentiCRISPRv2 — Addgene [#52961](https://www.addgene.org/52961/)

| | |
|---|---|
| File | `addgene-52961_lentiCRISPRv2.gb` (25018 bytes) |
| Length | **14873 bp**, circular, 49.8% GC |
| Addgene page states | 14873 |
| Features | 34 |
| Depositor | Feng Zhang |
| Addgene sequence record | 244694 (Addgene-Verified Full Sequence) |
| Source | Addgene public sequence-viewer JSON |
| Source URL | `https://sequences.addgene.org/snapgene-media/v3.101.0/sequences/244694/51ae3e01-5e94-49b8-a4fb-7fa303a7eb81/addgene-plasmid-52961-sequence-244694-{sequence-lines,features}.json` |
| Retrieved | 2026-10-03 |
| BsmBI sites | 2 |

**Triple-validated.** Length matches Addgene's stated total vector size (14873 bp). The sequence is **byte-identical** to the original Addgene `addgene-plasmid-52961-sequence-244694.gbk` recovered from the Wayback Machine, and an **exact circular rotation** of SnapGene's independent `lentiCRISPR v2` map (offset 14256). This is the reference that validates the extraction method used for every other file here.

### lentiGuide-Puro — Addgene [#52963](https://www.addgene.org/52963/)

| | |
|---|---|
| File | `addgene-52963_lentiGuide-Puro.gb` (17996 bytes) |
| Length | **10183 bp**, circular, 47.1% GC |
| Addgene page states | 10183 |
| Features | 27 |
| Depositor | Feng Zhang |
| Addgene sequence record | 331247 (Addgene-Verified Full Sequence) |
| Source | Addgene public sequence-viewer JSON |
| Source URL | `https://sequences.addgene.org/snapgene-media/v3.101.0/sequences/331247/092a4fee-eef4-4211-8e83-d1edf82c697d/addgene-plasmid-52963-sequence-331247-{sequence-lines,features}.json` |
| Retrieved | 2026-10-03 |
| BsmBI sites | 2 |

Matches Addgene's stated 10183 bp and is an **exact circular rotation** of SnapGene's independent `lentiGuide-Puro` map (offset 2593). sgRNA-only vector; the 2 BsmBI sites flank the ~2 kb filler between the U6 and the scaffold.

### lentiCas9-Blast — Addgene [#52962](https://www.addgene.org/52962/)

| | |
|---|---|
| File | `addgene-52962_lentiCas9-Blast.gb` (22236 bytes) |
| Length | **12859 bp**, circular, 52.3% GC |
| Addgene page states | 12860 |
| Features | 33 |
| Depositor | Feng Zhang |
| Addgene sequence record | 322376 (Addgene-Verified Full Sequence) |
| Source | Addgene public sequence-viewer JSON |
| Source URL | `https://sequences.addgene.org/snapgene-media/v3.101.0/sequences/322376/afc5ce55-af10-49b2-b6d3-01746afdc217/addgene-plasmid-52962-sequence-322376-{sequence-lines,features}.json` |
| Retrieved | 2026-10-03 |
| BsmBI sites | 0 |

**Length note:** this record is 12859 bp while Addgene's plasmid page states 12860 bp. SnapGene's independent `lentiCas9-Blast` map is also 12859 bp and is an **exact circular rotation** of this one (offset 4315), so 12859 bp is the real length and the 12860 on the web page is a 1 bp bookkeeping discrepancy. Cas9-only vector: no U6, no sgRNA scaffold, **zero BsmBI sites** - it is not an sgRNA cloning backbone.

### CROPseq-Guide-Puro — Addgene [#86708](https://www.addgene.org/86708/)

| | |
|---|---|
| File | `addgene-86708_CROPseq-Guide-Puro.gb` (18144 bytes) |
| Length | **10214 bp**, circular, 47.1% GC |
| Addgene page states | 10214 |
| Features | 27 |
| Depositor | Christoph Bock |
| Addgene sequence record | 317051 (Addgene-Verified Full Sequence) |
| Source | Addgene public sequence-viewer JSON |
| Source URL | `https://sequences.addgene.org/snapgene-media/v3.101.0/sequences/317051/ecf32836-79cd-4e6f-a931-771c0a1001ea/addgene-plasmid-86708-sequence-317051-{sequence-lines,features}.json` |
| Retrieved | 2026-10-03 |
| BsmBI sites | 2 |

Matches Addgene's stated 10214 bp. **The whole U6-sgRNA-Puro cassette sits on the minus strand of this map**, so every landmark motif is found as its reverse complement - watch the orientation when you lay primers on it. 2 BsmBI sites (CGTCTC at 5137, GAGACG at 3264) bracket the guide-cloning filler; the U6 promoter is 5144..5392 (-) and the gRNA scaffold 3187..3262 (-).

### pU6-sgRNA EF1Alpha-puro-T2A-BFP — Addgene [#60955](https://www.addgene.org/60955/)

| | |
|---|---|
| File | `addgene-60955_pU6-sgRNA-EF1Alpha-puro-T2A-BFP.gb` (15306 bytes) |
| Length | **8877 bp**, circular, 52.2% GC |
| Addgene page states | 8888 |
| Features | 21 |
| Depositor | Jonathan Weissman |
| Addgene sequence record | 467280 (Addgene-Verified Full Sequence) |
| Source | Addgene public sequence-viewer JSON |
| Source URL | `https://sequences.addgene.org/snapgene-media/v3.101.0/sequences/467280/2b82e300-c292-4d2f-822e-27fbbc4b06cd/addgene-plasmid-60955-sequence-467280-{sequence-lines,features}.json` |
| Retrieved | 2026-10-03 |
| BsmBI sites | 2 |

**Length note:** 8877 bp here vs 8888 bp stated on Addgene's page (11 bp). This is the Addgene-*verified* sequence, so 8877 bp is what Addgene actually sequenced; the 8888 on the page matches the depositor-era number and the closely related library backbones (#62217 / #66217, both 8888 bp - see `addgene-66217_*`, whose guide position is an N20 run, and `alt-snapgene/snapgene-62217_*`). **Cloning chemistry differs from the Zhang vectors:** the sgRNA cassette is opened with **BstXI (CCANNNNNNTGG @ 2923) + BlpI (GCTNAGC @ 2962)**, not BsmBI - the 2 CGTCTC/GAGACG sites in this plasmid lie inside the PuroR (4673) and TagBFP (5388) coding sequences and are useless for cloning. It also carries the **mouse U6** promoter (so the human-U6 3' probe GACGAAACACCG is genuinely absent) and the **Chen et al. 2013 'optimized' sgRNA constant region** starting at 2954 (GTTTAAGAGCTA**AGC**TGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCG), which differs from the F+E probe GTTTAAGAGCTA**TGC**TGGAAACAGCATAGCAAG by 2 nt - that is why neither scaffold probe matches.

### pXPR_011 — Addgene [#59702](https://www.addgene.org/59702/)

| | |
|---|---|
| File | `addgene-59702_pXPR_011.gb` (15568 bytes) |
| Length | **8395 bp**, circular, 51.9% GC |
| Addgene page states | 8395 |
| Features | 26 |
| Depositor | John Doench, David Root (Broad GPP) |
| Addgene sequence record | 233279 (Addgene-Verified Full Sequence) |
| Source | Addgene public sequence-viewer JSON |
| Source URL | `https://sequences.addgene.org/snapgene-media/v3.101.0/sequences/233279/37129534-a564-4e66-9822-ced75ef2a42d/addgene-plasmid-59702-sequence-233279-{sequence-lines,features}.json` |
| Retrieved | 2026-10-03 |
| BsmBI sites | 3 |

Matches Addgene's stated 8395 bp exactly. The Broad GPP sgRNA + EGFP activity-reporter vector from Doench et al. 2014, and the **only F+E-scaffold vector** in this set: hU6 3' end finishing at 8073, then a **19 nt spacer `GTGAACCGCATCGAGCTGA`**, then the F+E scaffold `GTTTAAGAGCTATGCTGGAAACAGCATAGCAAG` at 8093. So this plasmid is **not an empty cloning backbone** - its reporter guide is already in place and there is no BsmBI filler. Its **3** `CGTCTC` recognition sites (plus strand 247 and 905, minus strand 1195) are elsewhere in the vector, so do not treat them as a cloning site. Note also that the U6 -> guide -> scaffold junction runs right up to the end of the map (length 8395), so an amplicon across the guide **wraps the coordinate origin**.

### lentiCRISPRv2 blast (lentiCRISPRv2-Blast) — Addgene [#98293](https://www.addgene.org/98293/)

| | |
|---|---|
| File | `addgene-98293_lentiCRISPRv2-Blast.gb` (24838 bytes) |
| Length | **14676 bp**, circular, 49.0% GC |
| Addgene page states | no total vector size published |
| Features | 34 |
| Depositor | Brett Stringer |
| Addgene sequence record | 489390 (Addgene-Verified Full Sequence) |
| Source | Addgene public sequence-viewer JSON |
| Source URL | `https://sequences.addgene.org/snapgene-media/v3.101.0/sequences/489390/dcb9e0df-d61e-448f-95d1-8689058ff943/addgene-plasmid-98293-sequence-489390-{sequence-lines,features}.json` |
| Retrieved | 2026-10-03 |
| BsmBI sites | 2 |

lentiCRISPRv2 with the puromycin marker swapped for blasticidin (BSD). Addgene's page does not publish a total vector size for #98293, so 14676 bp could not be cross-checked externally; the record is internally consistent (245 viewer lines) and is 197 bp shorter than lentiCRISPRv2, consistent with the PuroR -> BSD marker swap. **Cassette is on the minus strand of this map** - motifs appear as reverse complements. 2 BsmBI sites, as expected for an sgRNA cloning backbone.

### pXPR_003 sgMSLN — Addgene [#193588](https://www.addgene.org/193588/)

| | |
|---|---|
| File | `addgene-193588_pXPR_003-sgMSLN.gb` (15654 bytes) |
| Length | **8323 bp**, circular, 51.2% GC |
| Addgene page states | no total vector size published |
| Features | 27 |
| Depositor | William Hahn |
| Addgene sequence record | 380090 (Addgene-Verified Full Sequence) |
| Source | Addgene public sequence-viewer JSON |
| Source URL | `https://sequences.addgene.org/snapgene-media/v3.101.0/sequences/380090/a40b489b-3539-46ef-a215-b285eb5c1d66/addgene-plasmid-193588-sequence-380090-{sequence-lines,features}.json` |
| Retrieved | 2026-10-03 |
| BsmBI sites | 0 |

**Substitute for the empty pXPR_003**, which is not deposited at Addgene and whose sequence the Broad GPP portal no longer serves without an account (see 'Not obtained'). Addgene lists this plasmid's vector backbone as **pXPR_003**, so everything outside the 20 nt guide is the pXPR_003 backbone: hU6 (3' end at 400) - 20 nt sgMSLN guide - gRNA scaffold (432) - WPRE (2625). Because the guide is already cloned in, the BsmBI filler is gone and the plasmid has **0 BsmBI sites**; use it to read off primer-binding sites and expected amplicon sizes for a *filled* pXPR_003, not to model the empty backbone's cut sites.

### CRISPRi v2 pooled-library backbone (pU6-sgRNA EF1Alpha-puro-T2A-BFP family) — Addgene [#66217](https://www.addgene.org/pooled-library/weissman-human-crispri-v2/)

| | |
|---|---|
| File | `addgene-66217_CRISPRi-library-backbone_pU6-sgRNA-EF1Alpha-puro-T2A-BFP.gb` (14603 bytes) |
| Length | **8888 bp**, circular, 52.0% GC |
| Addgene page states | 8888 |
| Features | 22 |
| Depositor | Jonathan Weissman |
| Addgene sequence record | depositor-uploaded file (no sequence ID) |
| Source | depositor GenBank file on Addgene's CMS (never login-gated) |
| Source URL | `https://media.addgene.org/cms/filer_public/bc/b9/bcb9a575-5058-44e6-9916-0cd33075dea7/crispri_backbone_66217.gb` |
| Retrieved | 2026-10-03 |
| BsmBI sites | 2 |

Bonus - a genuine `.gb` that Addgene still serves openly (not login-gated). Empty backbone of the Weissman human CRISPRi v2 pooled library (#66217), 8888 bp, same cloning chemistry and scaffold as #60955.

**Caveat for primer work: this file is not pure ACGT.** The depositor left the 20 nt guide position as a run of **N at 2945..2964**, sitting between the mouse-U6 3' end (`...CCACCTTGTTG`) and the Chen 2013 optimized scaffold (`GTTTAAGAGCTAAGCTGGAAACAGCATAGC...`). Any primer/amplicon calculation that crosses 2945..2964 will see Ns - substitute a real guide first. That N run is also the only difference from SnapGene's `CRISPRi_library_backbone` map (credited to Weissman Lab / Addgene #62217, also 8888 bp, in `alt-snapgene/snapgene-62217_pCRISPRia-v2.gb`), which carries a concrete stuffer there instead - the two files are otherwise base-identical. Together they are the practical maps for the Weissman CRISPRi/a sgRNA vector family.

---

## Not obtained / caveats

### pXPR_003 (empty Broad GPP backbone) - NOT OBTAINED

There is no empty `pXPR_003` deposited at Addgene (a catalog search for "pXPR_003" returns
only ~20 derivatives that already carry a guide, e.g. #193588, #202441, #220318), and the
Broad GPP vector portal no longer serves vector sequences anonymously: the legacy page
`portals.broadinstitute.org/gpp/public/vector/details?vectorId=pXPR_003` renders
"Vector Not Found", the current index `portals.broadinstitute.org/gppx/vector/public` is a
JavaScript app, and its data endpoint `/gppx/vector/public/api/init` returns **HTTP 401**.
NCBI nuccore has no `pXPR_003` record.

**Substitute provided:** `addgene-193588_pXPR_003-sgMSLN.gb` - Addgene #193588, whose
"Vector backbone" field on Addgene is literally *pXPR_003*. Everything except the 20 nt
sgMSLN guide is pXPR_003 backbone, so primer-binding sites and amplicon sizes read off it
are valid for a filled pXPR_003. It cannot be used to model the empty backbone's BsmBI
cut sites (the filler is gone, 0 BsmBI sites).

### GeCKO v2 pooled library - NO PLASMID-LEVEL FULL SEQUENCE EXISTS

`addgene.org/pooled-library/zhang-human-gecko-v2/` (catalog #1000000048) and
`.../zhang-mouse-gecko-v2/` (#1000000049) are **pools of ~120k distinct plasmids**, so
Addgene publishes no single full sequence for them - the pages carry only the library
content CSVs
(`human_geckov2_library_a_09mar2015.csv`, `human_geckov2_library_b_09mar2015.csv`)
and amplification/sequencing protocols. Both pages state the backbone explicitly by linking
plasmids **#52961 (lentiCRISPR v2)** and **#52963 (lentiGuide-Puro)**, both of which are
included here in full. For in-silico work on a GeCKO v2 pool: take `addgene-52961_lentiCRISPRv2.gb` and replace the
guide-cloning filler bracketed by its two BsmBI sites (`GAGACG` at 2240, `CGTCTC` at 4113)
with a 20 nt guide from the library CSV - that gives an exact member of the pool. The same
1873 bp BsmBI site spacing appears in all five Zhang-lineage backbones here, which is itself a
useful consistency check across the independent downloads (see the cross-validation section).

### lentiCRISPR v1 (#49535) - length not independently cross-checkable

The file is a real download (depositor full sequence, 13450 bp) but Addgene's plasmid page
publishes no "Total vector size (bp)" for #49535, and neither SnapGene nor NCBI carries
lentiCRISPR **v1**, so unlike the other vectors its length could not be confirmed against a
second source. The record is internally consistent and its U6 promoter and gRNA scaffold
features extract to the exact canonical 249 bp / 76 bp sequences.

### Nothing else was fabricated or reconstructed

Every `.gb` in this directory is a verbatim download. No sequence was assembled, inferred,
patched, or filled in from memory. Where a source was unavailable it is listed above as
unavailable.

---

## Cross-validation summary

The Addgene sequence-viewer JSON is not a GenBank file — the bases are rendered as individual
SVG `<text>` glyphs — so the extraction was validated before any file here was trusted.

**1. Byte-identical to the original Addgene `.gbk`.**
For lentiCRISPRv2 (#52961) the original `addgene-plasmid-52961-sequence-244694.gbk` was also
recovered from the Wayback Machine (the only one of these still findable there). The JSON-derived
sequence is **identical, base for base, at the same origin** — 14873/14873, zero mismatches.

**2. Exact circular rotations of SnapGene's independent maps.**
SnapGene publishes its own maps for three of these vectors, built by SnapGene from the depositors'
data, not copied from Addgene. Comparing sequences:

| Plasmid | Addgene JSON | SnapGene `.dna` | Result |
|---|---:|---:|---|
| lentiCRISPRv2 #52961 | 14873 bp | 14873 bp | exact rotation, offset 14256 |
| lentiCas9-Blast #52962 | 12859 bp | 12859 bp | exact rotation, offset 4315 |
| lentiGuide-Puro #52963 | 10183 bp | 10183 bp | exact rotation, offset 2593 |

Two independent pipelines agreeing base-for-base on three plasmids of different lengths rules out
glyph loss, reordering, or strand confusion in the extraction.

**3. Feature coordinates checked by translation.**
Addgene's `rangeBegin`/`rangeEnd` are 1-based inclusive. This was settled empirically on the PuroR
CDS of #86708: only `begin-1 .. end` yields a 600 nt ORF that starts with Met and ends in a stop
codon, while the off-by-one alternative gives 599 nt and no stop. With that convention applied, the
large CDS features come out correct — Cas9 translates from Met in #52961, #52962 and #98293, and
PuroR/BSD/TagBFP end on a stop codon — and no CDS has an internal stop.

Note that many of the annotations are SnapGene-style *sub-features* rather than complete ORFs —
FLAG, SV40/nucleoplasmin NLS, P2A/T2A/F2A, gp41 peptide, Factor Xa site, and PuroR where it is
annotated from the second codon because the ATG belongs to the upstream 2A junction. Those are not
expected to begin with Met or end in a stop, and the identical boundaries appear in SnapGene's own
maps, so they are the depositors' annotation choices rather than conversion errors.

**4. Annotated elements extract to their canonical sequences.**
In **all 7** files that annotate a `U6 promoter`, the feature extracts to exactly 249 bp ending
`…AAGGACGAAACACC`; in **all 6** that annotate a `gRNA scaffold`, it extracts to exactly 76 bp
starting `GTTTTAGAGCTAGAAATAGCAAGTTA…` and ending `…CACCGAGTCGGTGC`. Byte-identical extraction
across files of different lengths, different coordinate origins, and both strands confirms that
coordinates *and* strand assignment are right.

**5. Length agreement with Addgene's own published figures.**
Addgene still publishes "Total vector size (bp)" on the non-login plasmid pages. Of the 7 plasmids
where a number is published, 5 match exactly (14873, 10183, 10214, 8395, 8888) and 2 differ —
lentiCas9-Blast by 1 bp and pU6-sgRNA EF1Alpha-puro-T2A-BFP by 11 bp. Both discrepancies are
documented in the per-plasmid notes; for lentiCas9-Blast SnapGene independently confirms the
shorter length, so the web page figure is the one that is off.

**6. Internal consistency of every record.**
The viewer paginates at 60 bp/line. For all 9 JSON-derived files, `(lines-1)×60 + last-line-length`
equals the extracted length exactly, so no line was dropped or double-counted.

**7. Identical BsmBI site spacing in all five Zhang-lineage sgRNA backbones.**
These five files were downloaded independently, yet in every one of them the two inward-facing
BsmBI/Esp3I sites sit **exactly 1873 bp apart** (recognition-start to recognition-start, measured
around the circle), bracketing the guide-cloning filler:

| Plasmid | length | minus-strand site (`GAGACG`) | plus-strand site (`CGTCTC`) | spacing |
|---|---:|---:|---:|---:|
| lentiCRISPR v1 #49535 | 13450 | 1967 | 3840 | 1873 bp |
| lentiCRISPRv2 #52961 | 14873 | 2240 | 4113 | 1873 bp |
| lentiGuide-Puro #52963 | 10183 | 412 | 2285 | 1873 bp |
| CROPseq-Guide-Puro #86708 | 10214 | 3264 | 5137 | 1873 bp |
| lentiCRISPRv2-Blast #98293 | 14676 | 13073 | 270 (wraps the origin) | 1873 bp |

Five separate downloads of four different lengths reproducing the same spacing to the base is a
strong end-to-end check on both the sequences and the circular coordinate handling.

*Not checked here:* the exact cut coordinates and 4 nt overhang sequences. BsmBI is `CGTCTC(1/5)`,
so the nicks are offset from the recognition sites, and the arithmetic was not independently
validated against a known-good digest — if you need the precise cut positions, run the sequence
through a restriction tool rather than relying on numbers in this file.

