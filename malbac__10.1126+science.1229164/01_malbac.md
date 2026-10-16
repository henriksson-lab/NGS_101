# MALBAC — quasi-linear whole-genome amplification of a single cell

> **Evidence marking.** 🟢 verbatim from the source · 🟡 derived or inferred · 🔴 not
> published (or not in any source that could be fetched). Relationships marked 🟡
> *(computed)* were worked out with `lib/`; the model is `tools/malbac.py`.
>
> **The scg_lib_structs drawn page is not a source here.** It is inspiration only: no
> sequence, temperature, enzyme or library structure in this note or in the model is taken
> from it, and agreement with it is not evidence. An earlier revision of this note carried
> that page's reconstruction — a 27-nt handle `GTGAGTGATGGTTGAGGTAGTGTGGAG`, an extension
> temperature of 68 °C, "Bst or Phi29" as the polymerase, and two worked example
> Nextera/TruSeq libraries — and built a base-composition argument on the handle. The
> inventor patent discussed below independently prints that 27-mer for a related MALBAC
> embodiment, but the original paper's SOM remains unavailable and the other details
> remain unverified. None of those details is used as established model input.
>
> **What is left, and it is enough for the schematic.** The accessible PMC author
> manuscript states the architecture and the whole temperature programme in words, and
> `tools/malbac.py` draws exactly that and nothing more: the handle is a 27-character
> `inferred=True` placeholder (`H` × 27), the random foot an 8-character `N` placeholder,
> the insert an `X` placeholder, and the downstream platform arms the literal
> `[PLATFORM ARM]` token. **The model contains no MALBAC base at all** — there is nothing
> for a secondary source to have contaminated.

**MALBAC** (Multiple Annealing and Looping Based Amplification Cycles) — Zong C, Lu S,
Chapman AR, Xie XS. "Genome-Wide Detection of Single Nucleotide and Copy Number
Variations of a Single Human Cell." *Science* 2012 Dec 21;338(6114):1622–1626.
doi:[10.1126/science.1229164](https://doi.org/10.1126/science.1229164) · PMID 23258894 ·
PMC3600412 (NIH author manuscript).

Related entries in `catalogue/scg_lib_structs.tsv` that cite the same paper or build on it:

- **DR-seq** (`dr-seq__10.1038+nbt.3129`) — genome + transcriptome from one cell, uses a
  MALBAC-style preamplification; separate note.
- **MALBAC-DT** (`malbac-dt__10.1101+2019.12.31.892190`, bioRxiv 2020) — a transcriptome
  adaptation; separate note.

Sources read (fetched by `tools/get_sources.py malbac__10.1126+science.1229164`, into
`/husky/henriksson/NGS101_papers/sources/malbac__10.1126+science.1229164/`, never committed):

| File | What | Used for |
|---|---|---|
| `science.1229164_PMC3600412.html(.txt)`, `…xml(.txt)` | PMC author manuscript, full main text + figure legends | **everything asserted here**: principle, all four temperatures, cycle count, primer architecture (lengths only), results |
| `upstream_MALBAC.html(.txt)` | scg_lib_structs MALBAC page | **inspiration only, not a source** — nothing in this note or the model is cited from, or verified against, it |

Could not be fetched 🔴:

| File | URL | Why |
|---|---|---|
| SOM, `NIHMS445444-supplement-supplementary_data.docx` (1.3 MB) | https://pmc.ncbi.nlm.nih.gov/articles/instance/3600412/bin/NIHMS445444-supplement-supplementary_data.docx | PMC served a reCAPTCHA page; it is quarantined in `_invalid/` and the manifest carries it as `(manual)`. PMC3600412 is a gated NIH author manuscript (licence `TDM`, open access `False`). |
| Science SOM PDF | https://www.science.org/doi/suppl/10.1126/science.1229164 | publisher bot wall (HTML returned) |

The SOM holds the Materials and Methods: **the 27-nt handle sequence**, lysis buffer,
polymerase(s), reagent concentrations, exact cycling program, PCR cycle number, library kit
and sequencer. It is now requested through the PMC Cloud bucket in
`catalogue/extra_sources.tsv` (`pmccloud:PMC3600412`), the route that already served this
article's JATS XML when the gated article page would not.

An alternate primary source narrows, but does not close, the handle gap: the original
inventors Xie, Zong and Lu filed [US20140200146A1, *Methods of amplifying whole genome
of a single cell*](https://patents.google.com/patent/US20140200146A1/en) before the Science
paper. Its Example II lists the common 27-mer `GTGAGTGATGGTTGAGGTAGTGTGGAG` as the PCR
primer (paragraph around lines 1279–1286 in the linked text). Its worked preamplification
primers, however, have **5 random bases followed by GGG or TTT**, whereas the accessible
Science main text describes an **8-variable-base** foot. The patent therefore supports
this sequence for an inventor MALBAC embodiment, but does not establish that the Science
experiment used the same exact handle. The paper's SOM is still needed; the model keeps
the handle inferred and its bases unknown.

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

## 2. Oligos — 🔴 the Science experiment's exact handle remains unverified

Only two oligos exist, and the main text describes both in words without printing either
🟢 (`science.1229164_PMC3600412.xml.txt`):

> "…initiated with a pool of random primers, each having a common 27-nucleotide sequence
> and 8 variable nucleotides that can evenly hybridize to the templates at 0 °C."
>
> "After the five cycles of linear preamplification, only the full amplicons can be
> exponentially amplified in the following PCR using the common 27-nucleotide sequence as
> the primer."

So the architecture is 🟢 and the bases are 🔴. What the model draws:

| Oligo | `tools/malbac.py` | Established | Not established |
|---|---|---|---|
| random primer | `random_primer()` = `handle()` + `random_foot()` = `H` × 27 (`placeholder`, `inferred=True`, note "27-nt handle length is published; exact bases are unavailable") + `N` × 8 | the 27 + 8 architecture and the 35-nt total 🟢 | every base of the handle 🔴 |
| PCR primer | `pcr_primer()` = the same 27-character handle placeholder | that it *is* the common 27-mer 🟢 | its bases 🔴 |

`handle()` is written so that it **cannot** be constructed as established: `inferred=True`
is set inside the helper, not by its callers, so the uncertainty is a property of the
constructor rather than of each use site.

### What can still be said without the bases — 🟡 *(computed)*

- N8 → 4⁸ = 65,536 possible 3' ends.
- A full amplicon is `handle · N8-primed copy · … · revcomp(handle)`, so its 3' end is
  complementary to its 5' end whatever the bases are. That is what lets it loop at 58 °C
  and what lets one primer species amplify it from both ends — and it is the whole
  mechanism of the method. `malbac.loop_topology()` encodes it as a `HairpinTopology`
  with two 27-nt stems and a computed loop length, so the stem geometry is checked by
  construction instead of being drawn by hand.
- Because the handle has no bases, `malbac.protocol_graph()` gives the PCR
  `PCRGeometry.unavailable(..., reason="The paper gives the common-primer length but not
  its bases, so exact landing geometry cannot be asserted.", symmetric_species=True)`.
  The single-species claim is still recorded; the landing coordinates are not invented.
- The earlier revision's base-composition argument (a C-free, hairpin-free tag) rested
  entirely on the drawn page's handle and is withdrawn 🔴; whether the real handle has
  that property is unknown.

## 3. Step by step

Conditions 🟢 where they are in the main text; everything else is 🔴 and left unstated.
The four temperatures and the cycle count below are the model's constants
(`MELT_C`, `ANNEAL_C`, `EXTEND_C`, `LOOP_C`, `PREAMPLIFICATION_CYCLES`) and each is
verbatim in the sentence quoted beside it.

1. **Pick and lyse** one cell. 🟢 Lysis composition, volume and whether DNA is
   fragmented by protease/heat 🔴 (SOM). The paper notes the templates are DNA fragments
   of ~10–100 kb. 🟢
2. **Melt** the genomic DNA at **94 °C**. 🟢
3. **Anneal** the random primers at **0 °C** ("evenly hybridize to the templates at
   0 °C"). 🟢
4. **Extend** at **65 °C** 🟢 ("At an elevated temperature of 65 °C, DNA polymerases with
   strand displacement activity are used to generate semiamplicons with variable lengths
   (0.5–1.5 kb)"). Products are **semi-amplicons** of 0.5–1.5 kb 🟢 (lengths):
   `5'-handle · N8 · genomic copy -3'` 🟡 *(structure, computed from the primer
   architecture)*. **Which** strand-displacing polymerase is 🔴 — the main text says only
   "DNA polymerases with strand displacement activity", and the model's
   `strand_displacing_polymerase` input is a bare `MaterialKind.CATALYST` with no identity.
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
8. **Library prep** of the amplified DNA with a standard kit. **The kit is not named**
   anywhere we hold 🔴 — see §5.

## 4. Amplified product — 🟡 *(assembled from the published architecture)*

Double-stranded, each strand:

```
5'- handle(27, bases 🔴) · <genomic insert, starting with the N8-primed bases> · revcomp(handle)(27) -3'
```

`malbac.full_amplicon()` builds exactly this, with the right-hand copy derived by
`complement_segments([handle(...)])` rather than typed, so the two ends stay
complementary by construction however the handle is later resolved.

The insert's first 8 nt on each side were contributed by the random N8 and may carry
priming mismatches. 🟡 (inferred; a general property of random priming, not discussed in
the main text)

## 5. Final library — 🔴 unestablished

MALBAC is a whole-genome amplification method, not a library chemistry: the paper ends at
micrograms of amplified DNA and sequences it, **without naming the library kit, the
adapters or the platform** anywhere in the accessible text. The model says exactly that and
no more:

- `library_kit` is an `ADAPTER` input with `inferred=True` and the source note
  "(kit identity and adapter chemistry not reported)".
- `unresolved_sequencing_library()` draws the two ends as the literal token
  `[PLATFORM ARM]`, `placeholder=True, inferred=True`, with the note
  "downstream library-preparation kit not reported".
- The step itself is named `"INFERRED — convert amplified DNA with the unreported
  sequencing-library kit"`, so the page caption carries the required `INFERRED — …`
  prefix, and the `SequencingPlan` carries `limitation="The defining paper leaves the
  downstream sequencing-library kit and primer sites unspecified."` with no declared
  primer.

What *can* be said 🟡 *(computed)*: the 27-nt handles survive only at the two ends of each
amplicon, so after fragmentation most reads are pure genomic sequence and only reads
spanning an amplicon end contain handle bases.

An earlier revision of this note printed worked Nextera and TruSeq libraries here as
"example structures". They were the drawn page's choice of kit, not the paper's, and are
withdrawn.

## 6. Sequencing

🔴 Platform and read length are in the SOM. From the main text 🟢: single cells
sequenced to ~25× mean depth gave ~85 % and up to 93 % genome coverage at ≥1×; CNV
calling used only 0.8× per cell (bin size in the SOM 🔴); allele dropout ~1 % vs ~65 % for MDA.
**No sequencing primer is declared** in the model, and the page renders no primer-binding
view: with the kit unknown there is no site to locate 🔴.

## 7. Open questions

- 🔴 **Which 27-nt handle the Science experiment used.** The inventors' earlier patent
  gives a candidate for a related embodiment, but its 5N+GGG/TTT primer feet differ from
  the paper's 8-variable-base description. The paper's SOM must resolve that difference.
  This is the gap that keeps the whole schematic at placeholder resolution, and it
  is the reason `catalogue/extra_sources.tsv` now requests `pmccloud:PMC3600412`.
- 🔴 Exact polymerase(s) for the linear phase and for PCR, buffer, and the time at each
  temperature — all in the SOM. The main text gives only "DNA polymerases with strand
  displacement activity".
- 🔴 Lysis chemistry and PCR cycle number.
- 🔴 Library kit and sequencing platform / read length used in the paper (§5).
- 🟡 A commercial MALBAC kit exists from Yikon Genomics (outside knowledge; Yikon is
  S. Lu's current address in the PMC record 🟢). A vendor manual, if one is published,
  would be a legitimate original for the handle — but we have not located one, and no URL
  is guessed here.

## 8. How this note was made (tool evaluation)

`tools/get_sources.py` fetched the PMC full text and XML plus the drawn page; the
supplement link returned a reCAPTCHA page that was saved under the `.docx` name. Download
validation has since caught it: the file now sits in `_invalid/` and the manifest carries
it as `(manual)` with its URL and reason, which is how it should behave.

`tools/scrape_primers.py` finds **no oligo at all** in the main text, correctly — the paper
prints none. There is therefore nothing for
`tools/selftest_source_oligos.py` to compare either: every oligo input in the graph is
`inferred=True`, so the guard reports `0 checked` rather than a false pass. That is the
honest state of this protocol, and it will change the moment the SOM arrives.

## Bounded primary-source checkpoint, 2026-10-16

The original main article (PMC3600412 text36–42) confirms the27-nt common handle, eight variable bases, 0/65/94/58°C operations and five preamplification cycles before handle PCR. Exact-handle scraper hits in the upstream drawing were not used as evidence. The named original SOM remains uninspected in this pass. The page now distinguishes unassigned handle/library/run choices from published absence and shows an explicit downstream warning without assuming a preparation chemistry. The separate patent comparison and complete molecular topology are deferred.
