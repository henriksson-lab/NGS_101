# CRISPR-MIP — third-party reference material

**These files are deliberately not committed.** They are other people's copyrighted material —
a preprint's supplementary table, a plasmid map — and are not ours to redistribute. A fresh
clone will not have them, and that is expected: `tools/selftest.py` guards every check that
reads one, records a SKIP naming the file and where to get it, and **passes** on the remaining
checks. A missing source file means the clone is incomplete, not that the chemistry is wrong.

Page builders are stricter: `tools/build_page.py` and `tools/show_probe.py` compute everything
from the real vector map, so without it they print one line saying which file is missing and
exit 1 rather than emitting a half-computed page.

- **Protocol:** Selinger M, Yakovenko I, Nazir I, Henriksson J. "CRISPR-MIP replaces PCR and
  reveals GC and oversampling bias in pooled CRISPR screens." *bioRxiv* 2024.03.28.587082,
  doi:[10.1101/2024.03.28.587082](https://doi.org/10.1101/2024.03.28.587082)

---

## Files expected in this directory

### `TableS2_primers_and_probes.xlsx`

| | |
|---|---|
| **What it is** | Supplementary Table S2, "primers & probes": 38 named oligos in three columns (`name`, `sequence (5' -> 3')`, `notes`) — `MIP_probe_1..9`, `P5_tracrRNA_fwd`, `P7_tracrRNA_rev`, the `CRISPR_PCR1/PCR2` primers and the `Brunello-P7-i*` index primers. |
| **Where to get it** | The preprint's supplementary file **`media-2.xlsx`**, from the bioRxiv page for doi:10.1101/2024.03.28.587082 (Supplementary Material → Table S2). Save it here under this name. |
| **Needed by** | `tools/selftest.py`, the Table S2 cross-check — the four checks that compare the probe and primer sequences **as assembled from `lib/` constants** against the table verbatim: all nine probes, `P5_tracrRNA_fwd`, `P7_tracrRNA_rev`, and the test that the per-probe variable stretch is 8 nt rather than 6. That last one is why the table matters: it is what caught the i7 index being 8 nt, which the text does not state. Nothing else reads the file — `crisprmip.py` carries the derived constants, and the pages are built from those. |
| **Also needs** | `openpyxl` (optional, `pip3 install openpyxl`). If the file is present but the module is not, the same checks skip, with that reason. |

### Nothing else

`ref/` holds only the one file above. Two further sources are used from outside this directory
and are listed here because the checks here depend on them:

| Source | Used for | Where to get it |
|---|---|---|
| `lenticrispr-v1-screening__10.1126+science.1247005/ref/plasmids/addgene-52961_lentiCRISPRv2.gb` | The real vector the screens used. Everything measured rather than quoted: the single capture site in the 13 kb vector, the 112-nt gap fill, the closed circle, the 269-bp inverse-PCR product, the restriction sites that spare the probe footprint, and all of `crisprmip.html`'s drawings. Guarded by `HAVE_VEC` in `tools/selftest.py`; fatal for the two page builders. | Addgene **#52961** (lentiCRISPRv2, Zhang lab). That directory's own `ref/plasmids/MANIFEST.md` records the exact route, which is no longer a plain download. |
| `ref/concepts/padlock-circularization.md` (repo root) | Prose only: the general padlock chemistry, linked from `01_crispr-mip.md`. No check reads it. | Written in this repo, not third-party. |

The Brunello kinome library itself (Addgene **#75314**) is named in the page as the library
that was screened, but no file for it is needed: only its backbone, #52961, is used.
