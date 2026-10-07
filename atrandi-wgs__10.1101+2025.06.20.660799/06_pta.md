# PTA (Primary Template-directed Amplification) — quick reference

> **Evidence marking.** 🟢 = verbatim from a paper, patent or vendor document · 🟡 = derived or
> inferred (reasoning given) · 🔴 = proprietary or unmeasured.

**Primary paper.** Gonzalez-Pena V, Natarajan S, Xia Y, Klein D, Carter R, Pang Y, Shaner B, Annu K,
Putnam D, Chen W, Connelly J, Pruett-Miller S, Chen X, Easton J, Gawad C. "Accurate genomic variant
detection in single cells with primary template-directed amplification." *PNAS* 2021;**118**(24):
e2024176118. doi:[10.1073/pnas.2024176118](https://doi.org/10.1073/pnas.2024176118). PMCID PMC8214697.
Open access.

**Kit used in our protocol.** BioSkryb ResolveDNA Whole Genome Single-Cell Core Kit, **#100954**;
30 °C / 2.5 h → 65 °C / 5 min (see `03_our_protocol.md` §5).

---

## 1. The three components

### Polymerase — phi29 🟢
> "PTA … takes advantage of the processivity, strong strand displacement activity, and low error
> rate of **phi29 polymerase** used in MDA." — PNAS, Introduction

🔴 The enzyme in the *shipped kit* is undisclosed (SDS lists only "glycerol / proprietary").
🟡 BioSkryb filed WO2021163052A3 *"Phi29 mutants and use thereof"*, so an engineered variant is
plausible but unconfirmed. The reaction also contains a single-stranded binding protein
(WO2023215524A2 ¶[0126]).

### Primer — exonuclease-resistant random hexamer 🟢
> "1 µL **500 µM exonuclease-resistant random primer**" — PNAS Methods

🟢 Length: "random primers **6-9 bases in length**" (WO2021022085A2 ¶[00188]); BioSkryb's own
figure labels it "Primer (6-9mer)" (TAS-007 Fig. 2B).
🟢 Structure, from the MDA literature it is inherited from — Dean et al. 2002, Methods:
> "Thiophosphate-modified random hexamer (**5'-NpNpNpNpsNpsN-3'**)"

i.e. six random bases where **the last two internucleotide linkages are phosphorothioate** and the
first three are ordinary phosphodiester. 🟡 Almost certainly Thermo Exo-Resistant Random Primer
**SO181** — the patent names ThermoFisher, and PNAS's 500 µM stock is exactly SO181's spec.

### Terminators — α-thio-ddNTPs 🟢
The detail that matters most. PNAS Methods:
> "**alpha-thio-ddNTPs (Trilink Bio Technologies) at equal ratios at a concentration of 1,200 µM**"

Resolved against the TriLink catalogue, these are **2′,3′-dideoxyribonucleoside
5′-O-(1-thiotriphosphate)** — N-8009 (ddA), N-8010 (ddC), N-8011 (ddG), N-8012 (ddT).

**Two distinct chemical features, acting at two different atoms. Both are needed:**

| Feature | Atom | Consequence |
|---|---|---|
| **2′,3′-dideoxy** — no 3′-OH | sugar | chain extension irreversibly blocked |
| **α-thio** — S on the α-phosphate | backbone | the new linkage is a **phosphorothioate**, resisting phi29's 3′→5′ proofreading exonuclease |

🟢 Why the second one is essential — PNAS Results:
> "when the amplicons underwent sequencing, we found that the reactions had created poor quality
> products. The mapping rates … were 15.0% ± 2.2 … we hypothesized that amplicons … could reprime
> similar repetitive regions … However, for this model to be accurate, **the polymerase would need
> to remove the terminator from the amplicons prior to extension.** To test this model, **we
> incorporated an alpha-thio group into the terminators, which created an exonuclease-resistant
> phosphorothioate bond when amplification was terminated.** … the percent of reads mapped
> increased from 15.0 ± 2% to **97.9 ± 0.62%**"

So plain ddNTPs give termination *only*; phi29 simply chews them off and repriming produces
chimeras. The α-thio group is what makes termination stick.

> ⚠ **Two terminology traps.**
> 1. The paper calls plain ddNTPs "reversible terminators". Chemically that is wrong — they are
>    irreversible chain terminators. The authors mean *removable by phi29's exonuclease, hence
>    functionally reversible in this reaction*. Gloss it if quoting.
> 2. PTA does **not** use Illumina-style 3′-O-blocked reversible terminators. The patents list
>    those as alternatives; the reduction to practice is α-thio-ddNTPs.

🔴 **The shipped kit's terminator is a declared trade secret.** All three current BioSkryb manuals
say only: *"the polymerase incorporates **proprietary nucleic acid bases** which result in the
termination of the extension of the amplicon."* The published chemistry (α-thio-ddNTPs, 1200 µM in
a Qiagen REPLI-g mix) and the commercial formulation are **different claims** — state them separately.
🟡 The patents' preferred range (10–200 µM) is 6–100× below the paper's 1200 µM, so the research
protocol and the kit are not the same formulation.

## 2. The amplicon, base-by-base

This is what the page has to draw.

```
 5'-N p N p N p N ps N ps N -- [newly synthesised DNA, normal phosphodiester] -- ddN -3'
    +------ random hexamer ------+                                              +- alpha-thio
       2 terminal PS linkages                                                      dideoxy
```

### The 5′ end 🟡
The random primer, incorporated as synthesised. **Standard oligo synthesis gives a 5′-hydroxyl,
not a 5′-phosphate**, and SO181's datasheet lists no 5′ modification.

> **Consequence for our protocol:** PTA amplicon 5′ ends are **not ligatable without kinasing**.
> This is why end-prep cannot be skipped even though fragmentation can — T4 polynucleotide kinase
> is present in the NEBNext Ultra II End Prep mix (and in Atrandi's equivalent). 🟡 No publication
> states this; it is reasoning from the reagent's chemistry.

### The 3′ end 🟢 — a dead end
🟢 US12559794B2 claim 1(b): *"a **3′ terminal alpha-thio dideoxy nucleotide**."*

- **No 3′-OH** → cannot be extended, **cannot be dA-tailed**, cannot donate to a ligation.
- **Phosphorothioate linkage** → 🟢 the patents claim 3′→5′ exonuclease activity is *"reduced by at
  least 85%"* (US12559794B2 cl. 16), elsewhere up to 99.5%.

**You cannot rely on an end-repair exonuclease to polish this end off.**

### ⚠ The unresolved "terminator removal" step
The sources contradict each other, and the page should show this rather than smooth it over:

| Source | Says |
|---|---|
| 🟢 US11643682B2 Fig. 1D | workflow box: *"**removing the terminator**, repairing ends, and performing A-tailing prior to adapter ligation"* |
| 🟢 WO2019148119A1 claim 80 | recites *"**removing at least one terminator nucleotide**"* as a method step |
| 🟢 Same spec, ¶[0101] | *"**irreversible terminators are not capable of substantial removal by an exonuclease**"* |
| 🟢 PNAS | drops the clause entirely: *"the terminated amplification products can undergo **direct ligation of adapters**"* |
| 🟢 ResolveDNA v1 protocol | no such step: amplification → bead cleanup → "ready to use" |

🔴 **No source names an enzyme that removes an α-thio-dideoxy 3′ terminus, or demonstrates that it
happens.** The step is asserted, never enabled.

🟡 **Most likely resolution.** A blunt end whose 3′ strand ends in a PS-linked ddN is unusable: it
cannot be A-tailed and cannot be exonuclease-rescued. So the molecules that enter the library are
almost certainly ligated at the ends that *do* carry a normal 3′-OH — produced by polymerase
dissociation before terminator incorporation, by nicks, and by fill-in from a normal 3′-OH at
recessed/staggered ends. **Each ~1 kb amplicon has two ends and only needs one.** This is reasoning,
not a published claim, and it is the most interesting open experiment in the area.

### Strandedness and branching
🟢 "creating **smaller double-stranded amplification products**" (asserted in PNAS and all three
preprint versions); BioSkryb TAS-097: *"PTA reaction products are double-stranded."*
🟡 The only evidence is indirect — Qubit **dsDNA** HS quantitation and native TapeStation traces.
A dsDNA assay shows dsDNA is present, not that ssDNA is absent. No denaturing gel, no S1/ExoI
sensitivity test.

🔴 **Branching has never been characterised, in either direction.** Zero occurrences of
"branch(ed)", "hyperbranched", "debranch", "T7 endonuclease" or "S1 nuclease" in the PNAS paper, all
three preprint versions, four patent specifications, the STAR Protocols paper, or any of 19 BioSkryb
PDFs. 🟡 It *should* be far less branched than MDA: branch points in MDA are displacement forks, and
PTA terminates extension after a few hundred to ~2000 nt, so each polymerase displaces far less
strand and exposes far shorter single-stranded tails. But nobody has shown it.

## 3. Quasi-linear amplification

🟢 PNAS, Fig. 1A legend:
> "in PTA, the incorporation of exonuclease-resistant terminators in the reaction result in smaller
> double-stranded amplification products that **undergo limited subsequent amplification**,
> resulting in a quasilinear process with more amplification originating from the primary template."

> ⚠ **The preprint said something stronger and the authors retreated from it.** v1 read
> *"cannot undergo further amplification"*; v3 and PNAS read *"undergo limited subsequent
> amplification"*. **Do not draw terminated amplicons as absolute dead ends** — draw them as poor
> templates.

The chain, as the sources state it:

1. 🟢 Termination caps extension length — intended to *"reduce the average amplicon product length
   by at least … 65%"* (US11643682B2 ¶[0072]).
2. 🟢 Short duplex amplicons are poor templates: *"amplicons comprising terminator nucleotides
   **form loops or hairpins which reduce a polymerase's ability to use such amplicons as
   templates**"* (¶[0073]).
3. 🟢 So the primary template stays preferred — BioSkryb Fig. 2B labels the reprimed genomic strand
   *"Primary template reprimed for additional amplification (favored)"*.
4. 🟢 Copies-of-copies become a minority: claim 1 requires *"at least 5% of the polynucleotides are
   direct copies"*, with the spec ranging to 95%.

🟡 **Why this reduces allelic dropout.** In MDA the first priming events at a locus seed an
exponential cascade, so whichever allele is primed first wins by orders of magnitude — that
stochastic head start *is* allelic skew. In PTA no site can run away, so copy number tracks how
often the *original* template was primed, a much lower-variance quantity. Coverage follows genome
content rather than priming luck, equally for both alleles.

🟢 Empirically — PNAS Results: *"**PTA had significantly diminished allelic dropout and skewing**
at those heterozygous sites (Fig. 3D)."* ⚠ PNAS reports **no numeric ADO percentage**; Fig. 3D is
distributional. Any "PTA ADO = X%" cited to that paper is a downstream re-derivation.

🟢 A diagnostic consequence: *"no template control reactions with the irreversible terminators
**did not have detectable amplification products**"* (Fig. 1B).

🟢 And a happy one for microbes: *"by relying more on the primary template for amplification, PTA
enabled **smaller circular templates, such as the mitochondrial genome, to compete** with the larger
nuclear chromosomes."* 🟡 Predicts PTA should favour small replicons — plasmids, phage, small
bacterial chromosomes.

## 4. Amplicon size

⚠ **PNAS gives no numbers.** SI Fig. S1A is a TapeStation trace with no quantification. Every
commonly quoted figure comes from patents or vendor documents. (The "400–700 bp" in PNAS Methods is
a *sequencing library*, not an amplicon.)

| Source | Size |
|---|---|
| 🟢 US11643682B2 ¶[0072] | *"average length of **50-2000 nucleotides** … as compared to an average product length of **>10,000 nucleotides for MDA**"* |
| 🟢 ResolveDNA v2.0 (TAS-097) | *"Average fragment size in this sample is **1275 bp**, which is typical"* |
| 🟢 ResolveDNA v1 manual | *"typically create fragments from a size range **200-4000 bp**"* |
| 🟢 BioSkryb bacteria app note | **centred ~900 bp**, range ~100–1500 bp — *"slightly smaller than observed for eukaryotic cells"* |
| 🟢 McGowan 2026 | *"**PTA amplification products are shorter than traditional MDA**, potentially limiting their utility for long-read sequencing"* |

**Defensible summary:** a broad distribution peaking near **1–1.3 kb**, bulk between **~250 bp and
~2 kb**, tail to ~4 kb. On **bacteria, shifted smaller — centred ~900 bp**. MDA, for contrast,
averages **>10 kb**.

Yields: >1 µg per human cell (v1); **~100–200 ng per single bacterium** (app note).

## 5. MDA vs PTA

| | MDA | PTA |
|---|---|---|
| Polymerase | phi29 | phi29 (same) |
| Primer | PS random hexamer | same class, 6–9mer |
| Extra component | — | **α-thio-ddNTPs** |
| Kinetics | exponential, locally runaway | **quasi-linear** |
| Template of record | copies-of-copies | **primary template** (≥5% direct copies) |
| Product size | **>10 kb** avg | **~250–2000 bp**, avg 1275 bp |
| Architecture | **hyperbranched network**, many ssDNA 5′ tails | discrete duplexes; branching uncharacterised |
| **Needs debranching?** | **Yes** | **No** — see §6 |
| Needs fragmentation? | yes | no in principle (PNAS: "direct ligation"); BioSkryb's current kit still runs one |
| NTC amplification | comparable to a real cell | **undetectable** |
| Reaction | 30 °C, 8–16 h | 30 °C, 2.5 h (v2.0) / 10 h (v1) |

**Where PTA loses** (for balance): 🟢 *"MDA/REPLI-g displayed improved coverage in high GC content
regions compared to PTA/ResolveDNA"* (Hernández-Hernández 2025, Fig. 2C); PTA underestimates ploidy
for CNV work (Lin 2026); and it still produces *"hundreds to thousands of false-positive single base
substitutions and indels in each amplification reaction"* with a C>T and homopolymer-indel signature
— hence the dedicated callers **SCAN2** (Luquette 2022) and **PTATO** (Middelkamp 2023).

## 6. ⚠ No debranching on the PTA path — now independently confirmed

🟢 **Negreira GH, Monsieurs P, Dujardin J-C, Domagalska MA.** *bioRxiv* 2025.09.10.675331 — PTA
substituted into the **Atrandi SPC + split-pool workflow**, i.e. the closest published precedent to
our protocol. Methods, verbatim:

> "For samples SPC-STD1-4, amplified DNA was debranched by combining … **10X Debranching Buffer
> (Atrandi)** and **Debranching Enzyme (Atrandi)** followed by an incubation at 37 °C for 1 hour …
> **The samples submitted to PTA SPC-PTA1-2 did not require debranching and were directly submitted
> to barcoding instead.**"

🟢 They **did** still run end-prep: *"End Prep Buffer (Atrandi) and End Prep Enzyme (Atrandi) …
20 °C for 30 minutes and inactivation at 65 °C for 30 minutes."*

**This independently confirms the omission we inferred from our own preprint** (`03_our_protocol.md`
§5): the MDA arm names "debranching *and* end preparation", the PTA arm names end preparation only.

🟡 The Atrandi kit's debranching enzyme runs **37 °C / 1 h** — identical to the T7 Endonuclease I
conditions in ONT's and NEB's protocols — so it is likely T7 Endo I or a functional equivalent, and
it belongs to that kit's **MDA** chemistry. 🟡 Likewise the Atrandi WGA polymerase runs at **45 °C**,
outside wild-type phi29's range but squarely within **EquiPhi29**'s (30–45 °C, optimum 42 °C).

**For the page: no debranching step in the PTA path. End-prep stays.**

## 7. ⚠ Amplicon leakage from capsules — the central unquantified risk

🟢 **Mullaney DB, et al.** "Capsule-Based Single-Cell Genome Sequencing." *bioRxiv* 2025.03.14.643253,
Supplement:

> "we also tested the retention of DNA during isothermal whole genome amplification using primary
> template amplification (PTA), a reaction performed at 30 °C. Analysis of the supernatant and
> capsule fractions after amplification showed **the presence of DNA ranging from 150bp to 1,500bp
> in the supernatant indicating the leakage of amplified template** even in the absence of
> thermocycling and elevated temperatures (Fig. S11)."

Amplicons were also found entering **empty** capsules (Fig. S12). Set against Atrandi's stated
retention cutoff — SPCs *"retain DNA fragments >500 bp"* (Ling et al. 2025) — **PTA's design goal
(small amplicons) collides directly with capsule-based compartmentalisation.** Worth stating on the
page as a caveat of the chemistry.

## 8. PTA in microbes

🟢 **Bowers RM, Gonzalez-Pena V, Wardhani K, et al.** "scMicrobe PTA: near complete genomes from
single bacterial cells." *ISME Communications* 2024;**4**(1):ycae085.
doi:[10.1093/ismeco/ycae085](https://doi.org/10.1093/ismeco/ycae085)

| | PTA | MDA | WGA-X |
|---|---|---|---|
| *B. subtilis* completeness | **94%** | 60% | failed |
| *E. coli* completeness | **91%** | 62% | — |
| Aquatic, median @1M reads | **83%** | 17% | 11% |
| Medium/high-quality drafts | **78%** | <10% | <10% |
| Median contamination | 1.5% | <0.1% | <0.1% |

Used the **ResolveDNA Bacteria kit** (SL-B lysis reagent, RT 30 min; 12 h at 30 °C), library prep by
Nextera tagmentation, **no barcoding**. 🔴 That kit has no public catalogue number or protocol.

**Prior art for PTA + capsules + split-pool is two preprints, neither bacterial** (Negreira:
*Leishmania*; Mullaney: human/mouse). 🟡 **Our preprint is the only bacterial PTA + SPC + split-pool
work on record.**

## 9. Decisions this settles for the page

1. **Name the chemistry at two levels.** Academic PTA = α-thio-ddNTPs, equal ratios, 1200 µM, TriLink,
   in Qiagen REPLI-g mix. The BioSkryb kit formulation is a **declared trade secret**. Different claims.
2. **Draw the 3′ end as a dead end** (`…-ddN(α-S)-3'`: no 3′-OH, PS linkage) and show the
   terminator-removal tension rather than smoothing it.
3. **No debranching in the PTA path.** End-prep stays.
4. **Flag amplicon leakage** as a caveat of PTA-in-capsules.

## 10. Key references

1. Gonzalez-Pena V, et al. *PNAS* 2021;118(24):e2024176118. doi:[10.1073/pnas.2024176118](https://doi.org/10.1073/pnas.2024176118)
2. Preprint: *bioRxiv* 2020.11.20.391961. doi:[10.1101/2020.11.20.391961](https://doi.org/10.1101/2020.11.20.391961)
3. Patents: WO2019148119A1 · **US11643682B2** (clean text — use for concentrations) · US11905553B2 · US12559794B2 · EP3746564B1
4. BioSkryb: WO2021022085A2 · WO2023215524A2 / US2026/0002203A1 · WO2021163052A3 (phi29 mutants)
5. TriLink N-8009/N-8010/N-8011/N-8012 (α-thio-ddNTPs); Thermo SO181 (exo-resistant random primer)
6. Dean FB, et al. *PNAS* 2002;99(8):5261–5266. doi:[10.1073/pnas.082089499](https://doi.org/10.1073/pnas.082089499)
7. Lasken RS, Stockwell TB. *BMC Biotechnol* 2007;7:19. doi:[10.1186/1472-6750-7-19](https://doi.org/10.1186/1472-6750-7-19)
8. Luquette LJ, et al. *Nat Genet* 2022;54(10):1564–1571. doi:[10.1038/s41588-022-01180-2](https://doi.org/10.1038/s41588-022-01180-2) (SCAN2)
9. Middelkamp S, et al. *Cell Genomics* 2023;3(9):100389. doi:[10.1016/j.xgen.2023.100389](https://doi.org/10.1016/j.xgen.2023.100389) (PTATO)
10. Bowers RM, et al. *ISME Commun* 2024;4(1):ycae085. doi:[10.1093/ismeco/ycae085](https://doi.org/10.1093/ismeco/ycae085)
11. **Negreira GH, et al.** *bioRxiv* 2025.09.10.675331. doi:[10.1101/2025.09.10.675331](https://doi.org/10.1101/2025.09.10.675331) — PTA in Atrandi SPCs
12. **Mullaney DB, et al.** *bioRxiv* 2025.03.14.643253. doi:[10.1101/2025.03.14.643253](https://doi.org/10.1101/2025.03.14.643253) — capsule leakage
13. Ling M, et al. *Front Microbiol* 2025;15:1516656. doi:[10.3389/fmicb.2024.1516656](https://doi.org/10.3389/fmicb.2024.1516656) — Atrandi SPC reference
14. Hernández-Hernández A, et al. *bioRxiv* 2025.10.30.685509. doi:[10.1101/2025.10.30.685509](https://doi.org/10.1101/2025.10.30.685509) — independent PTA vs MDA
