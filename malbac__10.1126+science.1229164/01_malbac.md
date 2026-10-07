# MALBAC — quasi-linear whole-genome amplification of a single cell

> **Evidence marking.** 🟢 verbatim from the source · 🟡 derived or inferred · 🔴 not
> published (or not in any source that could be fetched). Relationships marked 🟡
> *(computed)* were worked out with `lib/` while writing this note; they are not yet
> asserted in a self-test, because this protocol has no `tools/` module yet.
>
> **Source caveat.** The primer sequences and all reaction conditions beyond the main
> text live in the Supplementary Online Material (SOM), which could not be fetched (see
> below). Every oligo sequence in this note therefore comes from the upstream
> scg_lib_structs page and is 🟡 (secondary source), not 🟢.

**MALBAC** (Multiple Annealing and Looping Based Amplification Cycles) — Zong C, Lu S,
Chapman AR, Xie XS. "Genome-Wide Detection of Single Nucleotide and Copy Number
Variations of a Single Human Cell." *Science* 2012 Dec 21;338(6114):1622–1626.
doi:[10.1126/science.1229164](https://doi.org/10.1126/science.1229164) · PMID 23258894 ·
PMC3600412 (NIH author manuscript).

Related entries in `catalogue/scg_lib_structs.tsv` that cite the same paper or build on it:

- **DR-seq** (`dr-seq__10.1038+nbt.3129`) — genome + transcriptome from one cell, uses a
  MALBAC-style preamplification; separate note.
- **MALBAC-DT** (`malbac-dt__10.1101+2019.12.31.892190`, bioRxiv 2020) — a transcriptome
  adaptation, on upstream's TODO list; not covered here.

Sources read (fetched by `tools/get_sources.py malbac__10.1126+science.1229164`, into
`_data/sources/malbac__10.1126+science.1229164/`, never committed):

| File | What | Used for |
|---|---|---|
| `science.1229164_PMC3600412.html(.txt)` | PMC author manuscript, full main text + figure legends | principle, temperatures, primer architecture (lengths), results |
| `upstream_MALBAC.html(.txt)` | scg_lib_structs MALBAC page | primer sequences, step diagrams, example downstream libraries |

Could not be fetched 🔴:

| File | URL | Why |
|---|---|---|
| SOM, `NIHMS445444-supplement-supplementary_data.docx` (1.3 MB) | https://pmc.ncbi.nlm.nih.gov/articles/instance/3600412/bin/NIHMS445444-supplement-supplementary_data.docx | PMC served a reCAPTCHA page (saved under the `.docx` name, 21 kB HTML); Europe PMC says the article is not open access |
| Science SOM PDF | https://www.science.org/doi/suppl/10.1126/science.1229164/suppl_file/zong.sm.pdf (guessed path) | publisher bot wall (HTML returned) |

The SOM holds the Materials and Methods: lysis buffer, polymerase(s), reagent
concentrations, exact cycling program, PCR cycle number, library kit and sequencer.

---

## 1. What it is

A **single-cell whole-genome amplification (WGA)** method, not a library chemistry. It
turns the picograms of DNA in one lysed cell into micrograms of amplified DNA, which is
then made into a sequencing library with any standard kit. There is no cell barcode,
no UMI, no index of its own: one cell = one tube = one library. 🟡 (inferred: the main
text and Fig. 1 legend describe no barcode or UMI step)

The idea 🟢 (main text): exponential WGA amplifies whatever got ahead early, so bias
compounds. MALBAC instead runs **five quasi-linear cycles** in which only the original
genome and first-generation copies serve as templates, and only then switches to PCR.
The trick that makes the preamplification linear is **looping**: a copy that carries the
common primer sequence at both ends (a "full amplicon") folds into a hairpin at 58 °C
and is withdrawn from further priming.

| | New thing here | Where else it turns up |
|---|---|---|
| 1 | Random primer = **common 27-nt tag + 8 random nt**, so every copy of a copy has complementary ends | DOP-PCR / PicoPLEX-type primers also carry a 5' tag on a degenerate 3' end; MALBAC adds the looping step |
| 2 | **Looping** of full amplicons removes them from the template pool → quasi-linear preamplification | DR-seq (built on this); MALBAC-DT |
| 3 | Strand-displacing polymerase at elevated temperature, **re-melted every cycle** (not isothermal like MDA) | MDA is the isothermal phi29 alternative the paper compares against |

None of the repo's concept pages (Tn5 tagmentation, template switching, ligation, RT,
padlock) is used *by MALBAC itself*; the downstream library in §5 is a plain Nextera or
TruSeq prep — see [Tn5 tagmentation](../ref/concepts/tn5-tagmentation.md) for the
Nextera one.

## 2. Oligos

🟡 (secondary: upstream scg_lib_structs page; the paper's SOM could not be read). As
written there, no modifications given:

```
MALBAC random primers  5'- GTGAGTGATGGTTGAGGTAGTGTGGAGNNNNNNNN -3'
MALBAC PCR primer      5'- GTGAGTGATGGTTGAGGTAGTGTGGAG -3'
```

What the paper itself says about them 🟢 (main text): the random primers are a pool,
"each having a common 27-nucleotide sequence and 8 variable nucleotides"; the PCR uses
"oligos with the common 27-nucleotide sequence" as primers. No sequence is printed in
the main text.

### How they interlock — 🟡 (computed)

- **Lengths agree with the paper**: the common part is 27 nt and the random primer is
  27 + N8 = 35 nt — matching "27-nucleotide sequence and 8 variable nucleotides". The
  PCR primer is exactly the 27-nt common part.
- N8 → 4⁸ = 65,536 possible 3' ends.
- **Base composition of the 27-nt tag: G 14, T 8, A 5, C 0** — it contains **no C**
  (GC 52 %). A C-free sequence cannot pair with itself in G·C pairs, and no 6-mer of the
  tag has its reverse complement elsewhere in the tag (checked), so the tag forms no
  internal hairpin and primer copies cannot cross-pair through the tag. Whether this
  was the design intent is not stated in the main text. 🔴
- Reverse complement of the tag: `CTCCACACTACCTCAACCATCACTCAC`. A full amplicon is
  `tag · N8-primed copy · … · revcomp(tag)`, i.e. its 3' end is complementary to its
  5' end — which is what lets it loop, and what makes the single PCR primer amplify it
  from both ends.
- **No overlap with any Illumina / Nextera sequence in `lib/`** (`illumina.P5/P7`,
  TruSeq read primers, `nextera.ME/S5/S7` checked: no shared run ≥ 10 nt on either
  strand). The tag is a pure amplification handle; sequencing adapters are added later.

## 3. Step by step

Conditions 🟢 where they are in the main text; the upstream page adds the rest 🟡.

1. **Pick and lyse** one cell. 🟢 Lysis composition, volume and whether DNA is
   fragmented by protease/heat 🔴 (SOM). The paper notes the templates are DNA fragments
   of ~10–100 kb. 🟢
2. **Melt** the genomic DNA at **94 °C**. 🟢
3. **Anneal** the random primers at **0 °C** ("evenly hybridize to the templates at
   0 °C"). 🟢
4. **Extend** with a **strand-displacing polymerase** at an elevated temperature —
   **65 °C** in the paper 🟢; the upstream page says **68 °C** and names "Bst or Phi29"
   🟡. Products are **semi-amplicons** of 0.5–1.5 kb: `5'-tag · N8 · genomic copy -3'`. 🟢
   (lengths); 🟡 (structure). Note: phi29 is a mesophilic enzyme and would not work at
   65–68 °C, so "Phi29" on the upstream page is doubtful; the enzyme actually used is
   in the SOM. 🔴
5. **Melt** at 94 °C to release the semi-amplicons. 🟢
6. **Cycles 2–5**: anneal 0 °C → extend → melt 94 °C → **loop at 58 °C**. 🟢 Genomic DNA
   keeps producing new semi-amplicons; primers landing on a semi-amplicon copy it
   through to its 5' tag, giving a **full amplicon** `tag · insert · revcomp(tag)`
   🟡 *(computed)*. At 58 °C the two complementary 27-nt ends of a full amplicon pair
   intramolecularly and the looped molecule is no longer a template. 🟢 (mechanism,
   Fig. 1 legend) "Five cycles of preamplification" in total 🟢; exact hold times 🔴.
7. **PCR** with the 27-nt common primer to reach micrograms. 🟢 Only full amplicons have
   the primer site at both ends, so only they amplify exponentially 🟢 (Fig. 1 legend).
   Cycle number, polymerase 🔴.
8. **Library prep** of the amplified DNA with a standard kit. The paper's kit is not
   named in the main text 🔴; the upstream page shows Nextera or TruSeq as examples 🟡.

### Upstream page vs paper — agreements and disagreements

| Point | Paper (main text) | Upstream page | Verdict |
|---|---|---|---|
| Primer = 27-nt common + 8 random | 🟢 stated | sequence given, 27 + 8 | agree (lengths computed 🟡) |
| Melt / anneal / loop temps | 94 / 0 / 58 °C | 94 / 0 / 58 °C | agree |
| Extension temperature | 65 °C | 68 °C | **disagree**; SOM would decide |
| Polymerase | "with strand displacement activity" | "Bst or Phi29" | upstream more specific; phi29 implausible at these temps 🟡 |
| Number of linear cycles | five | 5 | agree |
| PCR primer = common 27-mer | 🟢 | given | agree |
| Final product strands | — | step 8 bottom strand written `3'- revcomp(tag) … tag -5'` | **notation slip** 🟡 *(computed)*: that line is not base-complementary to the top strand; the correct bottom strand read 3'→5' begins `CACTCACTACC…` (complement of the tag). The meaning (tag at both 5' ends) is right. |

## 4. Amplified product — 🟡 (assembled)

Double-stranded, each strand:

```
5'- tag(27: GTGAGTGATGGTTGAGGTAGTGTGGAG) · <genomic insert, starting with the N8-primed bases> · revcomp(tag)(27) -3'
```

The insert's first 8 nt on each side were contributed by the random N8 and may carry
priming mismatches. 🟡 (inferred; a general property of random priming, not discussed in
the main text)

## 5. Final library — 🟡 (example structures from the upstream page, checked with `lib/`)

The amplified DNA is fragmented and adapted by the downstream kit. The 27-nt tags
survive only at the ends of each amplicon, so most reads are pure genomic sequence;
reads from fragments that kept an amplicon end begin or end in the tag. 🟡

**Nextera** (upstream example; computed equal to `illumina.P5 · i5 · nextera.S5 ·
nextera.ME` and `nextera.ME_RC · nextera.S7_RC · i7 · illumina.P7_RC`):

```
5'- P5 · i5(8) · s5 · ME · <MALBAC DNA fragment> · ME' · s7' · i7'(8) · P7' -3'
```

**TruSeq** (upstream example, single-indexed; top-strand 5' part computed equal to
`illumina.NEBNEXT_UNIVERSAL_PRIMER`, 3' part to `revcomp(illumina.TRUSEQ_READ2) · i7 ·
illumina.P7_RC`):

```
5'- P5 · TruSeq Read 1 · <MALBAC DNA fragment> · TruSeq Read 2' · i7'(8) · P7' -3'
```

Which kit the paper actually used: 🔴 (SOM).

## 6. Sequencing

🔴 Platform and read length are in the SOM. From the main text 🟢: single cells
sequenced to ~25× mean depth gave ~85 % and up to 93 % genome coverage at ≥1×; CNV
calling used only 0.8× per cell (bin size in the SOM 🔴); allele dropout ~1 % vs ~65 % for MDA.
Standard read primers for whichever kit is used (upstream: "standard sequencing
workflow") 🟡.

## 7. Open questions

- 🔴 Exact polymerase(s) for the linear phase and for PCR, buffer, and the time at each
  temperature — all in the SOM. Upstream's "Bst or Phi29" is unverified.
- 🔴 Extension temperature: 65 °C (paper) vs 68 °C (upstream).
- 🔴 Lysis chemistry and PCR cycle number.
- 🔴 Library kit and sequencing platform / read length used in the paper.
- 🔴 Whether the C-free tag was designed deliberately against self-pairing (computed
  property only).
- 🟡 Whether the published primer sequence on the upstream page is exactly that in the
  SOM — cannot be checked without the SOM. (A commercial MALBAC kit exists from Yikon Genomics 🟡, outside knowledge;
  Yikon is S. Lu's current address in the PMC record 🟢.)

## 8. How this note was made (tool evaluation)

`tools/get_sources.py` fetched the PMC full text and the upstream page; the supplement
link returned a reCAPTCHA page that was saved under the `.docx` name, so `doctext`
skipped it ("File is not a zip file") and the tool did not flag it as `(manual)`.
`tools/scrape_primers.py` found only the upstream page's sequences (the main text has
none). Primer lengths, composition and the library segment identities were computed
with `lib/chemdraw.revcomp`, `lib/illumina` and `lib/nextera`.
