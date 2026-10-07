# scDamID: nuclear-lamina contacts mapped in single cells by Dam methylation

> **Evidence marking.** 🟢 verbatim from the source · 🟡 derived or inferred · 🔴 not
> published, or published only in a file that could not be fetched. Relationships marked
> 🟡 *(computed)* were worked out with `lib/` while writing this note. This protocol has
> no `tools/` module yet, so no self-test asserts them.
>
> **Note on sources.** No oligo sequence appears in any document from the authors that
> could be fetched. The figure with the adaptor and primers (Figure S1E), the Extended
> Experimental Procedures and the indexing-primer table (Table S2) are all in PMC
> supplements behind the download gate. Every sequence below therefore comes from the
> upstream scg_lib_structs page and is marked 🟡 (secondary source), not 🟢.

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
`_data/sources/scdamid__10.1016+j.cell.2015.08.040/`, never committed):

| File | What | Used for |
|---|---|---|
| `j.cell.2015.08.040_PMC4583798.html.txt` | PMC full text (author manuscript) | workflow (Results "Single-cell DamID methodology", Experimental Procedures, Fig. 1A and S1 legends), read length, QC numbers |
| `upstream_scDamID.html.txt` | scg_lib_structs page | **every oligo sequence**, the step drawings, final library, read layout |
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
| Vogel et al. 2007 *Nat Protoc* | https://doi.org/10.1038/nprot.2007.148 | where the AdRt / AdRb adaptor and the DamID PCR primer come from |

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

🟡 **All from the upstream scg_lib_structs page** (`upstream_scDamID.html.txt`, lines
22–67), written as that page gives them. The paper's own list (Fig. S1E, Table S2) was not
available to check against. `/Phos/` = 5' phosphate; `*` = phosphorothioate bond.

```
DamID adaptor (Vogel et al. 2007 design; annealed AdRt + AdRb):
AdRt   CTAATACGACTCACTATAGGGCAGCGTGGTCGCGGCCGAGGA                  42 nt
AdRb   TCCTCGGCCG                                                  10 nt

scDamID preAmp PCR primer (single primer):
       NNNNGTGGTCGCGGCCGAGGATC                                     23 nt

Illumina Y adaptor (annealed top + bottom):
top    /Phos/ GATCGGAAGAGCACACGTCT                                 20 nt
bottom ACACTCTTTCCCTACACGACGCTCTTCCGATCT                           33 nt

Illumina PCR Primer 1.0:
       AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATC*T  58 nt
Illumina Multiplexing PCR Primer (i7 = 6 nt):
       CAAGCAGAAGACGGCATACGAGAT [i7] GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC*T   64 nt
```

The upstream page also lists the standard sequencing primers: TruSeq Read 1
`TCTTTCCCTACACGACGCTCTTCCGATCT`, TruSeq Read 2
`GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT`, and the index read primer
`GATCGGAAGAGCACACGTCTGAACTCCAGTCAC`. 🟡

What the paper itself says about these oligos 🟢 (Fig. S1E legend): the adaptor and
primer sequences are in Fig. S1E; **NNNN is a random nucleotide tail added at the 5' end
of the primer**; and the PCR product gives a **51-nt read** with a defined structure. That
agrees with the upstream preAmp primer, which has 4 N at its 5' end. Upstream says the Ns
are there for base diversity on the sequencer. 🟡

### How they interlock 🟡 *(computed)*

**DamID adaptor**

- AdRb is the reverse complement of the **last 10 nt of AdRt** (AdRt[32:42],
  `CGGCCGAGGA`). The annealed adaptor is therefore **blunt** at the end that ligates to the
  DpnI fragment, with a 32-nt single-stranded 5' tail on AdRt.
- AdRt ends in `…GAGGA`. When it is ligated to a DpnI end (which starts `TC…`), the
  junction reads `GAGGA|TC…` and **rebuilds the GATC**. So every genomic fragment in the
  read starts with `GATC`. That is the motif the read filter checks: on average 92 % of
  the mapped locations started with GATC 🟢 (Results).
- The preAmp primer (minus NNNN) = **AdRt[25:42] + `TC`**, i.e. the last 17 nt of AdRt
  (`GTGGTCGCGGCCGAGGA`) followed by the two restored bases of the GATC. It primes only on
  a fragment that carries the adaptor **and** starts with the regenerated `ATC`. Because
  the same adaptor sits on both ends of a fragment, a **single primer** amplifies it.
- AdRt[1:20] contains `TAATACGACTCACTATAG`, the **T7 promoter**. Nothing in scDamID uses
  it, and the PCR primer does not include it, so it is not in the final library. It is
  left over from the conventional adaptor. 🟡

**Illumina adaptor and PCR primers** (checked against `lib/illumina.py`)

- Y-adaptor bottom = **`illumina.TRUSEQ_READ1`** exactly (33 nt).
- Y-adaptor top = the first 20 nt of **`illumina.INDEX1_PRIMER`**. Its first 12 nt are
  **`illumina.STEM`** (`GATCGGAAGAGC`), which pairs with **`illumina.STEM_COMPLEMENT`** in
  the bottom oligo (bottom[-13:-1]), leaving the bottom oligo's 3'-**T** overhang for
  T/A ligation. The other 8 nt (`ACACGTCT`) and the 5' 20 nt of the bottom oligo stay
  unpaired: the fork of the Y.
- PCR Primer 1.0 = **`illumina.P5` + `TRUSEQ_READ1[4:]`** (=
  `illumina.NEBNEXT_UNIVERSAL_PRIMER`, 58 nt).
- Multiplexing PCR Primer = **`illumina.P7`** + 6-nt i7 + **`illumina.TRUSEQ_READ2`**.
  The reverse complement of `TRUSEQ_READ2` is `A` + the Y-adaptor top oligo + `GAACTCCAGTCAC`.
  So this primer anneals to the copy of the top strand of the adaptor; Primer 1.0 cannot
  anneal until that first copy has been made (upstream's own remark). It is the standard
  TruSeq single-index (i7 only) arrangement.
- The 3' end of PCR Primer 1.0 and of the Multiplexing PCR primer (`…GATC*T`) carries a
  phosphorothioate bond, which protects it from 3' exonuclease activity. Whether the
  authors' oligos had it 🔴.

## 3. Step by step

Volumes, buffers, temperatures and incubation times are in the Extended Experimental
Procedures, which could not be fetched 🔴. The order is 🟢 from the paper (Results, Fig.
1A legend, Experimental Procedures). The molecular detail is 🟡 from the upstream page
and from computing with the oligos.

1. **Express and sort.** Induce Dam-LmnB1 for 15 h. FACS-sort single early-S-phase
   (Fucci double-positive) KBM7 cells directly into a small volume of lysis buffer in a
   96-well plate. 🟢 Upstream adds "digest protein" (a protease step) 🟡. The paper's text
   names lysis only.
2. **DpnI digestion** in the same well. This cuts Dam-methylated `G-m6A|TC`, blunt, and
   leaves unmethylated GATC uncut. Fig. 1A has a **heat-inactivation (HI)** step 🟢 (legend
   abbreviation). Exactly where it falls is shown only in the figure image 🔴.
   Conventional DamID adds a **DpnII** digestion after ligation to cut fragments with
   unmethylated internal GATCs. Whether scDamID keeps it is not stated in the text that
   could be read 🔴. The main text lists only DpnI → ligation → PCR.
3. **Ligate the DamID adaptor** (AdRt/AdRb) to the blunt DpnI ends. 🟢 (step), 🟡
   (adaptor identity). Product: `AdRt · GATC-fragment · rc(AdRt)` on one strand.
   Upstream's drawing shows **AdRb joined to the fragment** on the other strand. The
   upstream list gives AdRb with **no 5' phosphate**, though, so AdRb cannot be ligated.
   Only AdRt's 3' end can join the fragment's 5'-phosphate end from DpnI. The adaptor
   site on the opposite strand would then be made by extension in the first PCR cycle,
   after the short AdRb melts away. 🟡 This is an inferred **disagreement with the
   upstream drawing**. It does not change the final product.
4. **Pre-amplify with the single preAmp primer**, 26 cycles in total, in the same well.
   🟢 (cycle count). About half of the wells give a visible smear 🟢. Product (🟡):
   `NNNN · GTGGTCGCGGCCGAGGATC · <gDNA> · GATCCTCGGCCGCGACCAC · NNNN`.
   The T7 part of AdRt is lost because the primer starts inside AdRt.
5. **Illumina library prep from the PCR product**: purify, end-repair, A-tail, ligate the
   Y adaptor, PCR with PCR Primer 1.0 + an indexed Multiplexing PCR Primer. 🟡 (upstream).
   The paper says both "ligation of indexed adaptors" (Results) and "indexing primers as
   listed in Table S2" (Experimental Procedures) 🟢. These cannot be reconciled without
   Table S2. Upstream puts the index in the PCR primer.
6. **Pool and sequence.** 🟢

## 4. Final library 🟡 (computed; identical to the upstream final structure)

Top strand, 5'→3':

```
illumina.P5 (29) · TRUSEQ_READ1[4:] (29) · NNNN · GTGGTCGCGGCCGAGGA (AdRt[25:42]) · TC
  · <genomic insert, Dam-methylated GATC to GATC> · GA
  · TCCTCGGCCGCGACCAC (rc AdRt[25:42]) · NNNN · A (A-tail)
  · rc(TRUSEQ_READ2)[1:] (33) · i7 (6, as rc of the primer's [i7]) · illumina.P7_RC (24)
```

The same as a segment list: `P5 · Read1 · N4 · DamID-adaptor-stub · GATC · insert ·
GATC · rc(DamID-adaptor-stub) · N4 · A · rc(Read2) · i7 · P7'`. 77 nt of fixed and random
sequence sit 5' of the first `GATC`, and 83 nt sit 3' of the second `GATC`. The A-tail
base is also the first base of rc(TRUSEQ_READ2). The insert can be in either orientation, because the same adaptor is on both ends.

Built from the oligos with `lib/`, the sequence matches the upstream "Final library
structure" letter for letter outside the `XXX…XXX` insert. 🟡

## 5. Read layout and sequencing

- 🟢 **Single-end 51-nt reads** (Fig. S1E, S1G legends). Read count: 1.2 × 10⁶ raw reads
  per cell on average. Then 68 % pass the adaptor-structure filter, 40 % map uniquely,
  37 % map and start with GATC, and 13 % remain after duplicate removal. Median
  5 × 10⁵ uniquely mappable reads per cell. The instrument is not named in the text read 🔴.
- 🟡 *(computed)* Read 1 starts after TruSeq Read 1:
  `NNNN` (4) + `GTGGTCGCGGCCGAGGA` (17) + `TC` → read positions 1–4 random, 5–21 adaptor,
  **positions 20–23 are `GATC`**, and genomic sequence runs from position 24 to 51
  (**28 nt**). The paper's "filtered on the correct adaptor structure" step is consistent
  with this fixed 23-nt prefix 🟡. Its read model (51-nt structure in Fig. S1E) was not seen.
- 🟡 Index read 1: i7, **6 nt** (upstream), primed by `illumina.INDEX1_PRIMER`. This index
  is the cell identity. No index 2 (single-index TruSeq).
- 🟡 Upstream also draws a Read 2 from the TruSeq Read 2 primer. With the 51-nt
  single-end design reported in the paper, Read 2 was presumably not run.

## 6. Open questions

- 🔴 The **actual Figure S1E** sequences. Do they match the upstream AdRt / AdRb / preAmp
  primer exactly, including any 5' phosphate on AdRb and the length of the N tail?
- 🔴 **Table S2**: the index sequences, their length (upstream says 6 nt), and whether the
  index is in a ligated adaptor ("indexed adaptors", Results) or in a PCR primer
  ("indexing primers", Experimental Procedures; upstream).
- 🔴 Lysis buffer, protease step, DpnI and ligation conditions, the place of the heat
  inactivation, whether a DpnII step is kept, and how the 26 PCR cycles are split.
- 🔴 Whether the paired Read 2 in the upstream drawing was ever sequenced. The paper
  reports only 51-nt reads.
- 🟡 AdRb ligation: upstream draws it ligated, but it is listed without a phosphate (§3,
  step 3).

## 7. How this note was made (tool evaluation)

`tools/get_sources.py` got the upstream page and the PMC full text. Seven of the eight
PMC supplements were reported as `(manual)` (download gate). The eighth
(`supplement-7.pdf`) was **saved, but it was a Google reCAPTCHA HTML page**, not the PDF.
`doctext.py` then failed on it ("pdftotext failed"). The file was renamed to
`…supplement-7.pdf.recaptcha.html` so that a later run fetches it again. The Europe PMC
supplementaryFiles API and article-bin URLs were tried by hand: the API says the article
is not open access, and the bin URLs return 403. `tools/scrape_primers.py` found 19
sequences, **all in the upstream page**. The paper text contains no oligo. Every
`--find` for the oligos above hits only `upstream_scDamID.html.txt`, which is why no
sequence in this note is 🟢.
