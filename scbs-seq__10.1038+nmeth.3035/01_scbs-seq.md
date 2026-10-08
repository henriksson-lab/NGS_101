# scBS-seq — single-cell whole-genome bisulfite sequencing by post-bisulfite random priming

> **Evidence marking.** 🟢 verbatim from the source · 🟡 derived or inferred · 🔴 not
> published. Relationships marked 🟡 *(computed)* were worked out with `lib/` while
> writing this note; they are not yet asserted in a self-test, because this protocol has
> no `tools/` module yet (status `notes` in `catalogue/ours.tsv`).

**scBS-seq** — Smallwood SA, Lee HJ, Angermueller C, Krueger F, Saadeh H, Peat J,
Andrews SR, Stegle O, Reik W, Kelsey G. "Single-cell genome-wide bisulfite sequencing
for assessing epigenetic heterogeneity." *Nature Methods* 11:817–820 (2014).
doi:[10.1038/nmeth.3035](https://doi.org/10.1038/nmeth.3035) · PMC4117646 · PMID 25042786.
Data: GEO GSE56879.

Built on: **PBAT** (post-bisulfite adaptor tagging), Miura F *et al.*, *Nucleic Acids
Res.* 40:e136 (2012), doi:[10.1093/nar/gks454](https://doi.org/10.1093/nar/gks454) —
the paper's ref. 7. Indexing primer: "iPCRTag", Quail MA *et al.*, *Nat. Methods* 9:10–11
(2012) — ref. 20, not fetched.

Later uses of the same chemistry (separate catalogue entries, all citing this paper):
scM&T-seq (`scm-t-seq__10.1038+nmeth.3728`), scNMT-seq
(`scnmt-seq__10.1038+s41467-018-03149-4`), scNOMe-seq / scCOOL-seq
(`scnome-seq-sccool-seq__10.7554+eLife.23203`). A step-by-step *Nature Protocols* version
(Clark SJ *et al.* 2017, doi:10.1038/nprot.2016.187) exists but was not fetched. 🟡

Sources read (fetched by `tools/get_sources.py "scbs-seq__10.1038+nmeth.3035"`, into
`_data/sources/scbs-seq__10.1038+nmeth.3035/`, never committed):

| File | What | Used for |
|---|---|---|
| `nmeth.3035_PMC4117646.html.txt` | author manuscript, full text (PMC) | **all methods and oligos** (Online Methods, line 1287) |
| `upstream_scBS-seq.html.txt` | scg_lib_structs page | second account: oligos, read primers, iPCRtag, library drawing |
| `nmeth.3035_41592_2014_BFnmeth3035_MOESM643_ESM.pdf.txt` | Supplementary Figures + Note (Springer) | QC context (poly-T artefact, Supp. Fig. 2) |
| `nmeth.3035_supp_NIHMS59406-supplement-2.pdf.txt` | Supplementary Note 1 (statistics) | not chemistry |
| `nmeth.3035_supp_NIHMS59406-supplement-1.doc` | Supplementary Figures (Word; no `.txt` twin, read via `catdoc`) | nothing new — figure legends only |
| `*_MOESM644_ESM.xlsx.txt`, `*_MOESM645_ESM.xlsx.txt` (= `Supp_Table_1/2`) | library statistics, region coverage | not chemistry |
| `nmeth.3035_PMC4117646.xml`, `*.zip`, `*.gif/jpg` | duplicates / figure images | not used |

Not fetched / not available:

| What | Why it matters |
|---|---|
| Quail *et al.* 2012 (iPCRTag primer), doi:10.1038/nmeth.1814 | fetched after the first review; Supplementary Table 1 supplies the PE adapter, index-read primer and 24 indexed PCR primers 🟢 |
| Clark *et al.* 2017 *Nat. Protoc.* doi:10.1038/nprot.2016.187 | later, detailed protocol; may give oligo vendors/purifications |
| Miura *et al.* 2012 PBAT (doi:10.1093/nar/gks454) | the bisulfite clean-up "as described previously" |

---

## 1. What it is

Bulk whole-genome bisulfite sequencing ligates adapters first and then bisulfite-treats,
and the bisulfite step destroys most of the adapter-tagged fragments. With one cell's
worth of DNA that loss is fatal. scBS-seq turns the order round 🟢: **bisulfite first**
(which converts unmethylated C to U *and* fragments the DNA), **adapters after**, added
by **random priming** with oligos that carry an Illumina adapter tail and a 3' random
9-mer. The paper calls it "a modification of PBAT" (Miura 2012) 🟢. Its two defining
features, as described in the Online Methods 🟢 (which of them were already in the
original PBAT protocol was not checked, since Miura 2012 was not fetched 🟡 — PBAT is
generally described as also using a biotinylated first-strand primer and bead capture):

| | Feature | Consequence |
|---|---|---|
| 1 | **Five rounds** of oligo1 random priming + Klenow exo- extension, with a 95 °C denaturation between rounds | more of the few template strands get tagged, and each is copied several times; but oligo1 can now also prime on oligo1 copies, so the library is **non-directional** 🟢 (Online Methods, mapping) |
| 2 | Oligo1 is **5'-biotinylated**; tagged strands are captured on streptavidin beads, which are washed with **0.1 N NaOH** (🟡 presumably to strip the non-biotinylated template) | second-adapter priming (oligo2) and PCR are done **on the beads** |

There is **no cell barcode and no UMI in the chemistry** (none is described 🟡): one cell
per tube, identified by the index of the "indexed iPCRTag reverse primer" at the final
PCR, 12–14 cells pooled per sequencing run 🟢; the index length (8 nt) comes only from
the upstream page 🟡. None of the repository concept pages (Tn5, template switching,
ligation, RT, padlock) applies — this is DNA-polymerase random priming on bisulfite-
converted single-stranded DNA. 🟡

## 2. Oligos

🟢 Verbatim from the Online Methods (PMC text line 1287). `[Btn]` = 5' biotin as the
paper writes it; the upstream page writes it `/Bio/`.

```
oligo1             [Btn]CTACACGACGCTCTTCCGATCTNNNNNNNNN                              31 nt
oligo2             TGCTGAACCGCTCTTCCGATCTNNNNNNNNN                                   31 nt
PE1.0 forward      AATGATACGGCGACCACCGAGATCTACACTCTTTC-CCTACACGACGCTCTTCCGATCT      58 nt
iPCRTag reverse    "indexed iPCRTag reverse primer" — sequence not given, cites ref. 20
```

The hyphen in PE1.0 is a line-break artefact of the manuscript (`TCTTTC-CCTACAC`); the
bases are `AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT`, which is also
exactly how the upstream page writes it. 🟡

🟢 From Quail *et al.* 2012 Supplementary Table 1, the source cited by Smallwood *et al.*:

```
iPCRtag reverse primer      CAAGCAGAAGACGGCATACGAGAT [8-bp sample index] GAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT   74 nt
Read 1 sequencing primer    ACACTCTTTCCCTACACGACGCTCTTCCGATCT
Index sequencing primer     AAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTC
Read 2 sequencing primer    CGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT
```

The upstream page itself notes that this iPCRtag design follows an **obsolete** Illumina
adapter layout (the old PE / multiplexing read-2 arm, not the TruSeq one).

### Upstream vs. paper

| Point | Paper | Upstream | Verdict |
|---|---|---|---|
| oligo1, oligo2, PE1.0 bases | as above | identical | agree 🟢 |
| 5' modification of oligo1 | `[Btn]` | `/Bio/` | same thing, different notation |
| iPCRTag sequence | not given (cites Quail 2012) | given | confirmed in cited Quail Supplementary Table 1 🟢 |
| read / index primers | not given | given | confirmed in cited Quail Supplementary Table 1 🟢 |
| rounds of oligo1 priming | 1 + 1 + 3 = 5 | step "repeated five times" | agree (upstream wording loosely) |
| exonuclease I + 0.8× AMPure before capture | yes | not mentioned | upstream omits |
| capture, NaOH washes, oligo2 and PCR **on the beads** | yes | "purify products" between oligo2 and PCR | upstream simplifies |
| sample index called | "indexed" primer, pools of 12–14 | "this is the cell barcode" | agree in effect |

### How the oligos interlock — 🟡 (computed)

- **Oligo1 handle** (its 22 nt before the N9), `CTACACGACGCTCTTCCGATCT`, is
  `illumina.TRUSEQ_READ1[11:]` — the 3' 22 nt of the TruSeq Read 1 primer — and is the
  3' end of PE1.0.
- **PE1.0 = `illumina.P5` + `illumina.TRUSEQ_READ1[4:]`** (P5 ends and Read 1 begins
  with the same `ACAC`, so they overlap by 4 nt); it is identical to
  `illumina.NEBNEXT_UNIVERSAL_PRIMER`. PE1.0 therefore reaches the oligo1 handle with
  22 nt of overlap and extends it out to P5.
- **Oligo2 handle**, `TGCTGAACCGCTCTTCCGATCT`, is the 3' 22 nt of the upstream Read 2
  primer and of the iPCRtag 3' arm: iPCRtag 3' arm (42 nt) = `GAGAT` + Read 2 primer
  (37 nt), and the Read 2 primer ends with the oligo2 handle.
- **iPCRtag 5' arm = `illumina.P7`** exactly.
- The **index sequencing primer** is the dedicated 35-nt Quail primer. It anneals within
  the reverse-complemented iPCRtag arm; it is not the full 42-nt reverse complement.
- The two handles share their 3' 14 nt, `CGCTCTTCCGATCT`, which contains
  `illumina.STEM_COMPLEMENT` (`GCTCTTCCGATC`). So oligo1 and oligo2 differ only in their
  5' 8 nt (`CTACACGA` vs `TGCTGAAC`). Neither handle is the TruSeq Read 2 arm
  (that shares only the 13-nt `GCTCTTCCGATCT`).

## 3. Step by step

Volumes, temperatures and enzymes 🟢 from the Online Methods unless marked.

1. **Isolate and lyse.** MII oocytes by mouth pipette; ESCs by FACS (BD Influx, single
   drop, ToPro-3 / Hoechst gating for live G0/G1). Lysis: oocytes in 2× buffer (10 mM
   Tris-Cl pH 7.4, 2 % SDS) + 0.5 µL proteinase K, 12 µL final; ESCs in 12 µL of
   10 mM Tris-Cl pH 7.4, 0.6 % SDS, 0.5 µL proteinase K. 37 °C 1 h. Negative controls:
   empty tubes or sorted Accudrop beads.
2. **Bisulfite conversion** on the lysate: Imprint DNA Modification Kit (Sigma), all
   volumes halved; after chemical denaturation 65 °C 90 min, 95 °C 3 min, 65 °C 20 min.
   Clean-up as in PBAT (ref. 7); elute in 10 mM Tris-Cl pH 8.5. Unmethylated C → U,
   5mC stays C; DNA is now single-stranded and fragmented.
3. **Oligo1 round 1.** 24 µL: 0.4 mM dNTPs, 0.4 µM oligo1, 1× Blue Buffer; 65 °C 3 min,
   4 °C; add 50 U Klenow exo-; 4 °C 5 min, ramp +1 °C / 15 s to 37 °C, 37 °C 30 min.
   The N9 anneals anywhere on a converted strand and Klenow copies it, so the new strand
   is `5'-biotin · oligo1 handle · N9 · copy of the converted strand-3'`. Because the
   template is U-containing, the copy has A opposite U (read as T after PCR). 🟡
4. **Oligo1 rounds 2–5.** 95 °C 1 min, straight onto ice; add fresh oligo1 (10 pmol),
   Klenow exo- (25 U), dNTPs (1 nmol) in 2.5 µL; same annealing ramp and 37 °C 30 min.
   Four times in all after round 1 → five rounds. Strand displacement by Klenow exo-
   lets primers downstream be displaced (upstream). In later rounds oligo1 can prime on a
   strand made in an earlier round, which is why the library is non-directional 🟢, and
   a strand primed by oligo1 on such a copy can run into the earlier copy's handle and
   end with the **reverse complement of the oligo1 handle** (upstream's "product 2";
   these carry handle1 at both ends and cannot be amplified by PE1.0 + iPCRtag) 🟡.
5. **Remove primer.** 40 U exonuclease I, 37 °C 1 h (digests single-stranded oligo1 so
   free biotin-oligo does not occupy the beads 🟡), then 0.8× AMPure XP; elute in 10 mM
   Tris-Cl pH 8.5.
6. **Capture** on washed M-280 Streptavidin Dynabeads, 20 min, rotating, room
   temperature. Wash twice with **0.1 N NaOH** (strips the non-biotinylated bisulfite
   template off; 🟡 the paper gives only the washes) and twice with 10 mM Tris-Cl.
7. **Oligo2, one round, on the beads.** 48 µL: 0.4 mM dNTPs, 0.4 µM oligo2, 1× Blue
   Buffer; 95 °C 45 s, ice; 100 U Klenow exo-; 4 °C 5 min, +1 °C / 15 s to 37 °C,
   37 °C **90 min**. Oligo2's N9 primes on the captured strand and extends to the 5'
   end, copying the oligo1 handle. Product: a duplex with handle1 at one end and handle2
   at the other. 🟡
8. **Index PCR, on the beads.** Wash beads; 50 µL: 0.4 mM dNTPs, 0.4 µM PE1.0, 0.4 µM
   indexed iPCRTag reverse primer, 1 U KAPA HiFi HotStart, 1× HiFi Fidelity Buffer.
   95 °C 2 min; **12–13 ×** (94 °C 80 s, 65 °C 30 s, 72 °C 30 s); 72 °C 3 min.
   0.8× AMPure XP; Bioanalyser HS chip; KAPA library quantification.
9. Bulk controls (not single cell): one round of oligo1 at 0.8 mM dNTPs, 4 µM oligo1;
   oligo2 at 4 µM; 9–12 PCR cycles.

## 4. Final library — 🟡 (assembled from the oligos above and the upstream drawing)

Top strand, 5'→3':

```
5'- P5 · TruSeq Read 1 site · N9 (oligo1) · <bisulfite-converted insert> · rc N9 (oligo2) · rc oligo2 handle · rest of rc iPCRtag 3' arm · <8-nt index, rc of primer> · P7' -3'
```

Equivalently, written by parts:

```
AATGATACGGCGACCACCGAGATCTACAC      P5 (illumina.P5; its last 4 nt ACAC are also the first 4 of Read 1)
TCTTTCCCTACACGACGCTCTTCCGATCT      rest of the TruSeq Read 1 primer site (TRUSEQ_READ1[4:]); last 22 nt = oligo1 handle
NNNNNNNNN                          oligo1 random 9-mer (not genomic sequence)
<insert>                           bisulfite-converted DNA, either strand orientation
NNNNNNNNN                          oligo2 random 9-mer, as its reverse complement
AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTC   rc of the iPCRtag 3' arm (= index read primer site)
<8-nt index>                       reverse complement of the index written in iPCRtag
ATCTCGTATGCCGTCTTCTGCTTG           P7' (illumina.P7_RC)
```

Fixed length excluding the insert, random 9-mers and index: 58 + 42 + 24 = 124 nt;
with both N9 and the 8-nt index, 150 nt + insert. 🟡 (computed)

## 5. Sequencing

🟢 HiSeq 2500, rapid-run mode, **100 bp paired-end**, pools of 12–14 single-cell
libraries. Both reads are clipped by 9 nt before mapping (Trim Galore
`--clip_r1 9 --clip_r2 9 --paired`) 🟢; that these 9 nt are the oligo1 / oligo2 N9 is
inferred from the construct 🟡. Mapping with Bismark `--non_directional` 🟢, because the
multiple oligo1 rounds make the library non-directional 🟢; hence all four bisulfite
strands (OT, OB, CTOT, CTOB) can appear in Read 1 🟡. Reads were first mapped to human (to
remove contaminants; 1.4 % mapped there), then the rest to mouse in **single-end** mode.

Read primers 🟢: Read 1 and Read 2 are the obsolete Illumina PE primers, and the dedicated
35-nt Quail index-read primer reads the eight-base iPCRTag. This is single-indexed: there
is no i5 index or Index 2 read.

## 6. Open questions

- 🔴 Oligo purification (HPLC / PAGE) and whether the N9 is hand-mixed are not given.
- 🔴 The **poly-T** artefact (reads with ≥50 T, present in negative controls too, main
  cause of low mapping; Supplementary Fig. 2) — its origin is called "a contaminant" and
  "inherent to our current methodology"; no mechanism is given. A guess: oligo1 / Klenow
  extension on a contaminating poly(A) template. 🟡
- 🟡 Read 1 vs. Read 2 methylation information: whether the first ~9 nt *after* the
  clipped N9 still carry primer-induced mismatch bias is not discussed.
- The bulk-cell protocol uses only one oligo1 round; whether bulk libraries are therefore
  directional is not stated 🔴 (mapping parameters are given only for single cells).

## 7. How this note was made (tool evaluation)

`tools/get_sources.py` fetched the PMC author manuscript, the Europe PMC supplementary
zip and the Springer supplements (the two Supp. Table xlsx are duplicated under both
names). `tools/doctext.py` made no `.txt` twin of the Word supplement
`nmeth.3035_supp_NIHMS59406-supplement-1.doc` (legacy binary .doc, "not a format this
reads"); it was read with `catdoc` and contains only figure legends. All oligos are in a
single paragraph of the Online Methods, which `tools/scrape_primers.py` found; it named
the oligo1 hit "Btn" (from the `[Btn]` prefix) without parsing `[Btn]` as a 5'
modification, whereas the upstream `/Bio/` form was parsed. The iPCRTag primer exists
only by citation, so it came from the upstream page.
