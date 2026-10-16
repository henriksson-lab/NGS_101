# scDamID: nuclear-lamina contacts mapped in single cells by Dam methylation

> **Evidence marking.** 🟢 verbatim from the source · 🟡 derived or inferred · 🔴 not
> published, or published only in a file that could not be fetched. Relationships marked
> 🟡 *(computed)* were worked out with `lib/`; the model is `tools/scdamid.py`.
>
> **The scg_lib_structs drawn page is not a source here.** It is inspiration only: no
> sequence, length or reaction in this note or in the model is taken from it, and
> agreement with it is not evidence.
>
> **Where the bases do come from.** No oligo sequence appears in any document *from these
> authors* that could be fetched: Figure S1E (adaptor and primers), the Extended
> Experimental Procedures and Table S2 (indexing primers) are PMC supplements behind the
> download gate 🔴. What the accessible text does say is that the methylation-detection
> steps are "similar to that of conventional DamID (Vogel et al., 2007) and use largely
> the same reagents" 🟢, and the Figure S1E legend states the primer's "NNNN: random
> nucleotide tail added at the 5′ end" and the resulting **51-nt** read 🟢. The reused
> conventional-DamID oligos — **AdRt, AdRb, AdR_PCR** — are themselves published, and this
> repository holds a fetched original that prints them verbatim with their Vogel-2007
> attribution (see the source table). They are therefore 🟡 for scDamID: the reuse is
> stated, the per-oligo identity in *this* paper is not, so `tools/scdamid.py` creates
> every one of them through `predecessor()`, which forces inferred styling.

**scDamID**: Kind J, Pagie L, de Vries SS, Nahidiazar L, Dey SS, Bienko M, Zhan Y,
Lajoie B, de Graaf CA, Amendola M, Fudenberg G, Imakaev M, Mirny LA, Jalink K, Dekker J,
van Oudenaarden A, van Steensel B. "Genome-wide maps of nuclear lamina interactions in
single human cells." *Cell* 163:134–147 (2015).
doi:[10.1016/j.cell.2015.08.040](https://doi.org/10.1016/j.cell.2015.08.040),
PMID 26365489, PMC4583798 (NIH author manuscript). Data: GEO GSE69423 (single-cell and
conventional DamID).

Builds on **DamID**: van Steensel B, Henikoff S. "Identification of in vivo DNA targets
of chromatin proteins using tethered Dam methyltransferase." *Nat Biotechnol*
18:424–428 (2000). doi:[10.1038/74487](https://doi.org/10.1038/74487). The paper says
its methylation-detection steps resemble, and use largely the same reagents as, the
conventional mammalian DamID protocol of Vogel, Peric-Hupkes and van Steensel, *Nat Protoc* 2:1467–1478 (2007),
doi:[10.1038/nprot.2007.148](https://doi.org/10.1038/nprot.2007.148). That paper is not
in the catalogue and was not fetched.

Sources read (fetched by `tools/get_sources.py`, into
`/husky/henriksson/NGS101_papers/sources/scdamid__10.1016+j.cell.2015.08.040/`, never committed):

| File | What | Used for |
|---|---|---|
| `j.cell.2015.08.040_PMC4583798.html.txt` | PMC full text (author manuscript) | workflow (Results "Single-cell DamID methodology", Experimental Procedures, Fig. 1A and S1 legends), read length, QC numbers |
| `upstream_scDamID.html.txt` | scg_lib_structs page | **inspiration only, not a source** — nothing is cited from or verified against it |
| `/husky/henriksson/NGS101_papers/sources/epidamid__10.1101+2021.10.26.465688/j.molcel.2022.03.009_supp_mmc4.pdf.txt:2206-2224` and `…_PMC9153956.xml.txt:574-583` | EpiDamID (Markodimitraki et al., *Mol Cell* 2022) key-resources oligo table | prints **AdRt** `CTAATACGACTCACTATAGGGCAGCGTGGTCGCGGCCGAGGA`, **AdRb** `TCCTCGGCCGCG` and **AdR_PCR** `GGTCGCGGCCGAGGATC`, each credited to "Vogel et al., 2007" 🟢 — the oligos scDamID says it reuses |
| `j.cell.2015.08.040_supp_NIHMS718969-supplement-7.pdf.recaptcha.html` | a reCAPTCHA page saved under a `.pdf` name (renamed by hand) | nothing |

Could not be fetched (PMC download gate; Europe PMC says the article is "not open
access"; no Cell / ScienceDirect copy was tried successfully):

| Wanted | URL | Why it matters |
|---|---|---|
| Supplement 6: supplementary figures, including **Figure S1E** (adaptor and primer sequences, 51-nt read structure) | https://pmc.ncbi.nlm.nih.gov/articles/instance/4583798/bin/NIHMS718969-supplement-6.pdf | primary source for every scDamID oligo |
| Supplement 7: probably the **Extended Experimental Procedures** (580.7 KB PDF, no caption on PMC) | https://pmc.ncbi.nlm.nih.gov/articles/instance/4583798/bin/NIHMS718969-supplement-7.pdf | lysis buffer, volumes, temperatures, cycle numbers, library prep, read processing |
| Supplement 8: Table S1 (GO) and **Table S2** (indexing primers) | https://pmc.ncbi.nlm.nih.gov/articles/instance/4583798/bin/NIHMS718969-supplement-8.doc | the per-cell index sequences |
| Supplements 1–5: OE-score data files | https://pmc.ncbi.nlm.nih.gov/articles/instance/4583798/bin/NIHMS718969-supplement-1.pdf (…-2 to -5) | not chemistry |
| van Steensel & Henikoff 2000 | https://doi.org/10.1038/74487 | original DamID; not needed for oligos |
| Vogel et al. 2007 *Nat Protoc* | https://doi.org/10.1038/nprot.2007.148 | the conventional DamID protocol this paper says it reuses, and the primary print of AdRt / AdRb / AdR_PCR. Declared in `catalogue/extra_sources.tsv`; until it is fetched, the EpiDamID table above is the original we actually read |

---

## 1. What it is

DamID marks where a protein sits on DNA without antibodies or crosslinking. The protein
of interest, here **Lamin B1**, is fused to *E. coli* **Dam** methyltransferase and
expressed in living cells. Wherever the fusion touches DNA it methylates the adenine of
`GATC` (G-m6A-TC), a mark that human cells otherwise lack. The restriction enzyme
**DpnI** cuts only methylated `GATC`, in the middle and blunt (`GA|TC`). So the fragments
that can be amplified are the ones between two methylated GATCs, and they come from the
regions that contacted the lamina. 🟢 (principle and DpnI specificity: paper, Results)
🟡 (DpnI cut position: standard enzymology, not stated in the text that could be read)

The single-cell version changes the handling, not the chemistry. 🟢 (Experimental
Procedures):

- a clonal line with **inducible Dam-LmnB1** (Shield1-stabilised, per the Fig. S1D
  legend) plus the **Fucci** reporter, induced for 15 h;
- single cells **FACS-sorted at the G1/S transition** (red+green, "yellow" Fucci) into
  lysis buffer in a 96-well plate;
- a different lysis solution 🔴 (composition is in the Extended Experimental Procedures);
- **every enzymatic step in the same well** by adding reagents one after another, with no
  purification or concentration in between;
- **4–6 more PCR cycles** than conventional DamID (26 in total);
- **Illumina sequencing** instead of microarrays.

There is **no cell barcode in the chemistry before pooling**. Each cell is its own well,
and cell identity enters only with the Illumina index added during library prep (Table
S2, not fetched). There is **no UMI**. Duplicates are collapsed by mapping position. 🟢
(Fig. S1G legend: "all duplicate reads are removed")

| | New thing here | Where else it turns up |
|---|---|---|
| 1 | DamID **in one tube per cell**: lysis, DpnI, ligation and PCR without cleanup | the same no-cleanup idea appears in most plate-based single-cell genome methods |
| 2 | A sequencing readout for DamID: Illumina adapters ligated **onto the DamID PCR product** | later scDamID&T-seq (Rooijers et al. 2019) adds a cell barcode to the DamID adapter and a T7 promoter for joint RNA 🟡 (from memory; not in the fetched sources) |

## 2. Oligos

The DamID half 🟡 — the **published conventional-DamID oligos the paper says it reuses**,
read in the EpiDamID key-resources table (each row credited to Vogel et al. 2007), not in
any scDamID document:

```
DamID adaptor (annealed AdRt + AdRb):
AdRt     CTAATACGACTCACTATAGGGCAGCGTGGTCGCGGCCGAGGA                42 nt
AdRb     TCCTCGGCCGCG                                              12 nt
AdR_PCR  GGTCGCGGCCGAGGATC                                         17 nt

scDamID pre-amplification primer = NNNN + AdR_PCR:
         NNNNGGTCGCGGCCGAGGATC                                     21 nt
```

The `NNNN` is the one part of that primer the defining paper states itself 🟢 (Fig. S1E
legend: "NNNN: random nucleotide tail added at the 5′ end of the primer"); the legend also
states the **51-nt** read structure 🟢. Any 5′ phosphate or phosphorothioate on these
oligos is 🔴 — the figure that would show the modifications is the gated one.

The Illumina half 🔴 — the paper names no kit, and Table S2 (the indexing primers) is
gated. The model therefore draws canonical `illumina` constants (`P5`,
`TRUSEQ_READ1[4:]`, `revcomp(TRUSEQ_READ2)`, `P7_RC`) and a 6-nt index **placeholder**,
all inferred, and gives the indexed-adapter ligation an `unavailable` geometry. Earlier
revisions of this note printed a full Y-adaptor, "PCR Primer 1.0" and a "Multiplexing PCR
Primer" with phosphorothioate ends; those reached it only through the drawn page and are
gone.

### How they interlock 🟡 *(computed, and enforced in `_validate()`)*

**DamID adaptor**

- AdRb is the reverse complement of the **last 12 nt of AdRt** (`CGCGGCCGAGGA`). The
  annealed adaptor is therefore **blunt** at the end that ligates to the DpnI fragment,
  with a 30-nt single-stranded 5' tail on AdRt. (The drawn page's AdRb is two bases
  shorter than the published oligo; the model now uses the published 12-mer.)
- AdRt ends in `…GAGGA`. When it is ligated to a DpnI end (which starts `TC…`), the
  junction reads `GAGGA|TC…` and **rebuilds the GATC**. So every genomic fragment in the
  read starts with `GATC`. That is the motif the read filter checks: on average 92 % of
  the mapped locations started with GATC 🟢 (Results).
- `AdR_PCR` = **last 15 nt of AdRt + `TC`**, i.e. it spans the adaptor/GATC junction and
  its 3′ end is the half-site the DpnI cut supplies. It primes only on a fragment that
  carries the adaptor **and** has the regenerated junction. Because the same adaptor sits
  on both ends of a fragment, a **single primer** amplifies it. `tools/scdamid.py` asserts
  both relations (`revcomp(ADRB) == ADRT[-12:]`,
  `ADR_PCR == ADRT[-15:] + JUNCTION`) in its construction code, so a mistyped oligo
  cannot be drawn.
- AdRt[1:20] contains `TAATACGACTCACTATAG`, the **T7 promoter**. Nothing in scDamID uses
  it, and the PCR primer does not include it, so it is not in the final library. It is
  left over from the conventional adaptor. 🟡

**Illumina adaptor and PCR primers** — 🔴 unestablished

The paper says both "ligation of indexed adaptors" (Results) and "indexing primers as
listed in Table S2" (Experimental Procedures) 🟢, which cannot be reconciled without
Table S2. No kit, adaptor or primer is named. The model therefore asserts no Illumina
oligo at all: it ligates an "indexed Illumina adaptor" with
`AdapterLigationGeometry.unavailable(...)`, runs the library PCR with
`PCRGeometry.unavailable(...)`, and draws the finished arms from canonical `illumina`
constants (`P5`, `TRUSEQ_READ1[4:]`, `revcomp(TRUSEQ_READ2)`, `P7_RC`) plus a 6-nt index
placeholder, every one of them inferred 🟡 *(the single-index TruSeq shell is the model's
assumption, not a published layout)*.

## 3. Step by step

Volumes, buffers, temperatures and incubation times are in the Extended Experimental
Procedures, which could not be fetched 🔴. The order is 🟢 from the paper (Results, Fig.
1A legend, Experimental Procedures). The molecular detail is 🟡 *(computed from the
reused DamID oligos)*.

1. **Express and sort.** Induce Dam-LmnB1 for 15 h. FACS-sort single early-S-phase
   (Fucci double-positive) KBM7 cells directly into a small volume of lysis buffer in a
   96-well plate. 🟢 The paper's text names lysis only; whether a protease step is used
   is in the gated Extended Experimental Procedures 🔴.
2. **DpnI digestion** in the same well. This cuts Dam-methylated `G-m6A|TC`, blunt, and
   leaves unmethylated GATC uncut. Fig. 1A has a **heat-inactivation (HI)** step 🟢 (legend
   abbreviation). Exactly where it falls is shown only in the figure image 🔴.
   Conventional DamID adds a **DpnII** digestion after ligation to cut fragments with
   unmethylated internal GATCs. Whether scDamID keeps it is not stated in the text that
   could be read 🔴. The main text lists only DpnI → ligation → PCR.
3. **Ligate the DamID adaptor** (AdRt/AdRb) to the blunt DpnI ends. 🟢 (step), 🟡
   (adaptor identity). Product: `AdRt · GATC-fragment · rc(AdRt)` on one strand.
   Which strand of the adaptor is joined is not established 🔴: AdRb's 5′ phosphorylation
   state is in the gated figure, and without a 5′ phosphate AdRb cannot ligate, so only
   AdRt's 3′ end would join the fragment's DpnI 5′-phosphate end and the opposite-strand
   adaptor site would be made by extension in the first PCR cycle, after the short AdRb
   melts away 🟡. Either way the final product is the same, which is why the model gives
   this ligation an `unavailable` end geometry rather than drawing a junction.
4. **Pre-amplify with the single primer** (`NNNN` + AdR_PCR), 26 cycles in total, in the
   same well. 🟢 (cycle count). About half of the wells give a visible smear 🟢. Product
   (🟡): `NNNN · GGTCGCGGCCGAGGATC · <gDNA> · GATCCTCGGCCGCGACC · NNNN`.
   The T7 part of AdRt is lost because the primer starts inside AdRt.
5. **Illumina library prep from the PCR product**: purify, end-repair, A-tail, ligate an
   indexed Illumina adaptor, PCR. 🟢 (steps, Results + Experimental Procedures) / 🔴 (every
   reagent: no kit, adaptor or primer is named, and the paper's "ligation of indexed
   adaptors" and "indexing primers as listed in Table S2" cannot be reconciled without the
   gated Table S2).
6. **Pool and sequence.** 🟢

## 4. Final library 🟡 *(computed)* — DamID interior reused, Illumina shell assumed

Top strand, 5'→3':

```
illumina.P5 (29) · TRUSEQ_READ1[4:] (29) · NNNN · GGTCGCGGCCGAGGA (AdRt[-15:]) · TC
  · <genomic insert, Dam-methylated GATC to GATC> · GA
  · CCTCGGCCGCGACC (rc) · NNNN · A (A-tail)
  · rc(TRUSEQ_READ2)[1:] (33) · i7 (placeholder) · illumina.P7_RC (24)
```

As a segment list: `P5 · Read1 · N4 · DamID stub · GATC · insert · GATC · rc(DamID stub)
· N4 · A · rc(Read2) · i7 · P7'`. The A-tail base is also the first base of
rc(TRUSEQ_READ2). The insert can be in either orientation, because the same adaptor is on
both ends. On the page every one of those regions is dotted.

## 5. Read layout and sequencing

- 🟢 **Single-end 51-nt reads** (Fig. S1E, S1G legends). Read count: 1.2 × 10⁶ raw reads
  per cell on average. Then 68 % pass the adaptor-structure filter, 40 % map uniquely,
  37 % map and start with GATC, and 13 % remain after duplicate removal. Median
  5 × 10⁵ uniquely mappable reads per cell. The instrument is not named in the text read 🔴.
- 🟡 *(computed)* Read 1 starts after TruSeq Read 1: `NNNN` (4) + AdR_PCR (17, ending in
  the restored `TC`) → read positions 1–4 random, 5–21 adaptor, with the regenerated
  `GATC` at positions 18–21 and genomic sequence from position 22 to 51 (**30 nt**). The
  paper's "filtered on the correct adaptor structure" step is consistent with a fixed
  21-nt prefix 🟡; the actual 51-nt read model in Fig. S1E was not seen 🔴.
- 🔴 Index read 1: an i7 carrying the cell identity, drawn 6 nt wide as a **placeholder**;
  Table S2 would give the real width and bases. No index 2 is modelled (the run is
  single-end, single-index).

## 6. Open questions

- 🔴 The **actual Figure S1E** sequences. Do they match the reused conventional-DamID
  AdRt / AdRb / AdR_PCR exactly, including any 5' phosphate on AdRb and the length of the
  N tail? Fetching Vogel et al. 2007 (declared in `catalogue/extra_sources.tsv`) would
  settle the predecessor oligos; only the gated figure can settle what *this* paper used.
- 🔴 **Table S2**: the index sequences and their width, and whether the index is in a
  ligated adaptor ("indexed adaptors", Results) or in a PCR primer ("indexing primers",
  Experimental Procedures).
- 🔴 Lysis buffer, protease step, DpnI and ligation conditions, the place of the heat
  inactivation, whether a DpnII step is kept, and how the 26 PCR cycles are split.
- 🔴 Whether AdRb is 5'-phosphorylated, i.e. whether it is ligated at all (§3 step 3).

## 7. How this note was made (tool evaluation)

`tools/get_sources.py` got the PMC full text (and the drawn page, which is not used).
Seven of the eight PMC supplements were reported as `(manual)` (download gate). The eighth
(`supplement-7.pdf`) was **saved, but it was a Google reCAPTCHA HTML page**, not the PDF.
`doctext.py` then failed on it ("pdftotext failed"). The file was renamed to
`…supplement-7.pdf.recaptcha.html` so that a later run fetches it again. The Europe PMC
supplementaryFiles API and article-bin URLs were tried by hand: the API says the article
is not open access, and the bin URLs return 403. `tools/scrape_primers.py` found **no
oligo in the paper** — the scDamID sources that count contain none. The reused
conventional-DamID oligos were instead read in the EpiDamID key-resources table already
in `/husky/henriksson/NGS101_papers/sources/` (`--find GGTCGCGGCCGAGGATC` hits it), which prints them with their
Vogel-2007 attribution; `catalogue/extra_sources.tsv` now declares both that table and
the Vogel paper for this slug.
