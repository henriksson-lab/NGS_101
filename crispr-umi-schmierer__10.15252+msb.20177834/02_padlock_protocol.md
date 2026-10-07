# Draft protocol — padlock capture of guides from pLenti-Puro-AU-flip-3xBsmBI

> **Status: design, not a validated protocol.** Every sequence and size below is computed
> from the reconstructed vector (`tools/design_padlock.py`, asserted in `tools/selftest.py`).
> Nothing here has been run. The section on what is likely to go wrong is the one to read
> first.

Companion to [`01_crispr-umi-schmierer.md`](01_crispr-umi-schmierer.md). Wet-lab steps are
adapted from CRISPR-MIP ([`../crispr-mip__10.1101+2024.03.28.587082/01_crispr-mip.md`](../crispr-mip__10.1101+2024.03.28.587082/01_crispr-mip.md)).

---

## 0. Read this before ordering anything

Two options, and they are **not** equivalent in risk.

| | **A — guide + RSL** | **B — guide only** |
|---|---|---|
| gap the polymerase must fill | **{{= crisprumi.PADLOCK_GAP }} nt** | **{{= crisprumi.PADLOCK_GAP_UNIVERSAL }} nt** |
| vs the {{= crisprumi.VALIDATED_GAP }} nt CRISPR-MIP is known to fill | **{{= round(crisprumi.PADLOCK_GAP / crisprumi.VALIDATED_GAP, 2) }}× — untested** | {{= round(crisprumi.PADLOCK_GAP_UNIVERSAL / crisprumi.VALIDATED_GAP, 2) }}× — comfortably shorter |
| captures the lineage label | **yes** | no |
| works on other guide vectors | no | **v1, v2, lentiGuide, this one** |
| sequencing | needs a **custom read primer** | standard |
| PCR primers both inside the capture | yes | **no** — one must sit in the backbone |

**Neither is a drop-in.** A risks the gap length; B forfeits the RSL *and* one of the four
selection filters. If the aim is to find out whether padlock capture works on this vector at
all, **run B first** — it is the lower-risk experiment and it isolates the gap-length
question. Then try A.

## 1. Oligos to order

### Probe A — captures guide and RSL ({{= crisprumi.PROBE_LEN_A }} nt)

```
5'-/5Phos/ AAGCTTGGCGTAACTAGATCTTGAGA NNNNNNNNNNNNN ATCACG ACTCTTTCCCTACACGACGCTCTTCCGATCT GTGGAAAGGACGAAACACC -3'
           \___ ligation arm, 26 ___/ \__ UMI 13 __/ \i7 6/ \______ Read-1 site, 31 ______/ \__ ext arm, 19 ___/
```

🔴 **Note what is *missing*: the probe carries no Read-2 site.** The vector already contains
`GATCGGAAGAGCACACGTCTGAACTCCAGTCAC`, which is a substring of the CRISPR-MIP backbone's
Read-2 site; including both would put two copies of one primer site in a single amplicon.
The vector's copy serves instead.

Order PAGE-purified, 5'-phosphorylated. `N` = hand-mixed equimolar.

### Probe B — captures guide only, works on any guide vector ({{= crisprmip.PROBE_LEN }} nt)

```
5'-/5Phos/ TCAACTTGAAAAAGTGGCACCGA NNNNNNNNNNNNN AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC ATCACG ACTCTTTCCCTACACGACGCTCTTCCGATCT GTGGAAAGGACGAAACACC -3'
           \__ ligation arm, 23 __/ \__ UMI 13 __/ \______ Read-2 site, 34 __________/ \i7 6/ \______ Read-1 site, 31 ______/ \__ ext arm, 19 ___/
```

Backbone identical to the published CRISPR-MIP probe; only the ligation arm differs.

### PCR primers (option A)

```
fwd (P5)  AATGATACGGCGACCACCGAGATCTACAC CACGTCTGAACTCCAGTCAC      3' Tm 55.7 C
rev (P7)  CAAGCAGAAGACGGCATACGAGAT      GTTGATAACGGACTAGCC        3' Tm 50.0 C
```

⚠ The reverse primer's 3' Tm is **50.0 °C**, noticeably cooler than the forward's 55.7. It is
the existing `P7_tracrRNA_rev` tail, kept so one fewer thing changes — but if PCR is
inefficient, lengthening it into the scaffold is the first thing to try.

### Custom sequencing primer (option A only)

```
CACGTCTGAACTCCAGTCAC          Tm 55.7 C
```

This is the **3' 20 nt of the standard TruSeq Read-2 / Index-1 primer**, so it is not exotic —
but the run still has to be set up to use it. **Agree this with the facility before you
start.** See §5.

## 2. The reaction

As for CRISPR-MIP, with one addition at the front.

| Step | Conditions |
|---|---|
| **0. Linearise** 🟡 | A supercoiled plasmid barely denatures &mdash; its strands are topologically interlinked and reanneal on cooling, so the arms never get a window to bind. Cut once outside the capture before hybridising. ⚠ The vector has a **third BsmBI site**; pick an enzyme that does not cut between the extension and ligation arms. |
| **1. Denature & hybridise** | 1–10 µg DNA, 0.2–0.002 µM probe, 1.5× Ampligase buffer. 94 °C 5 min → ramp to 60 °C at −0.1 °C/s → **overnight at 60 °C** |
| **2. Extend & ligate** | add Ampligase + 4 U Phusion HF + 0.2 mM dNTPs, pre-heated to 60 °C. **60 °C, 1 h** |
| **3. Exonuclease** | 37 °C, 45 min: 10 U Exo I + 50 U Exo III. Inactivate 80 °C, 20 min |
| **4. PCR** | cycle number by qPCR on a pilot |

🟡 **Step 2 is the one to extend if option A underperforms.** {{= crisprumi.PADLOCK_GAP }} nt is {{= round(crisprumi.PADLOCK_GAP / crisprumi.VALIDATED_GAP, 2) }}× the fill the
chemistry is known to manage; a longer extension (2 h) and more polymerase are the obvious
first variables, and are cheap to test.

## 3. Expected products ✅

| | option A | option B |
|---|---|---|
| probe | {{= crisprumi.PROBE_LEN_A }} nt | {{= crisprmip.PROBE_LEN }} nt |
| gap filled | {{= crisprumi.PADLOCK_GAP }} nt | {{= crisprumi.PADLOCK_GAP_UNIVERSAL }} nt |
| closed circle | **{{= crisprumi.CIRCLE_LEN_A }} nt** | {{= crisprmip.PROBE_LEN + crisprumi.PADLOCK_GAP_UNIVERSAL }} nt |
| final library | **{{= crisprumi.LIBRARY_LEN_A }} bp** | needs its own primer pair |

### Option A final library, computed

```
   0- 28   29 nt  Illumina P5
  29- 48   20 nt  vector adapter, 3' remnant     <- the custom read primes here
  49- 54    6 nt  RSL                            <- lineage label
  55- 80   26 nt  ligation arm
  81- 93   13 nt  UMI                            <- molecule label
  94- 99    6 nt  i7
 100-130   31 nt  Read-1 site
 131-149   19 nt  extension arm
 150        1 nt  Pol III +1 G
 151-170   20 nt  sgRNA spacer                   <- the guide
 171-220   50 nt  scaffold
 221-244   24 nt  Illumina P7
```

## 4. Controls

| control | expectation |
|---|---|
| **no probe** | nothing survives the exonuclease; no PCR product |
| **no ligase** | nothing survives — isolates ligation from hybridisation |
| **no Exo** | product plus a large background of un-circularised probe |
| **uncut plasmid** vs linearised | 🟡 if uncut gives nothing and cut works, step 0 is confirmed as necessary |
| **empty (filler) vector** | ✅ no product — the gap would be 1,972 nt |
| **plain lentiGuide-Puro** with probe B | ✅ should work (gap 66); with probe A, nothing — its ligation arm is absent |

That last row is a genuinely useful specificity control: **probe B works on both vectors,
probe A only on this one.**

## 5. Sequencing

| read | cycles | reads |
|---|---|---|
| Read 1 (standard primer) | **50** | extension arm 19 + G, then **the guide at cycles 21–40** |
| custom read, primer `CACGTCTGAACTCCAGTCAC` | **50** | **RSL at cycles 1–6**, ligation arm 7–32, **UMI at cycles 33–45** |

⚠ **This is the awkward part of option A.** Both reads run in the same direction on the same
strand, from different start points, so they do **not** map onto a standard paired-end
Read 1 / Read 2 pair. In practice that means a custom run recipe. Confirm the facility will
do it before committing — if they will not, option A does not work as drawn and the fallback
is B.

🟡 Shortening the ligation arm would bring the UMI closer and let a shorter custom read reach
both labels, at the cost of Tm margin. 26 nt was chosen for Tm, not for read economy.

## 6. What is most likely to go wrong

1. **The {{= crisprumi.PADLOCK_GAP }} nt gap** — the headline risk. Diagnostic: run probe B alongside. If B gives a
   product and A does not, the gap is the problem, not the probe design.
2. **Plasmid topology** — if step 0 is skipped, expect nothing regardless of probe.
3. **The cool reverse primer** (3' Tm 50.0 °C) — lengthen into the scaffold if PCR is weak.
4. **The custom read** — see §5.
5. **The reconstructed vector** 🔴 — everything here is computed against a vector *rebuilt*
   from lentiGuide-Puro plus the published edits, because we do not have the real map. If the
   facility's plasmid differs anywhere between the arms, these numbers move.
   **Get the real sequence before ordering.**

## 7. Before you order

- [ ] Obtain the actual plasmid map from the facility, and re-run `tools/design_padlock.py`
      against it rather than the reconstruction
- [ ] Confirm the facility will run a custom sequencing primer (option A only)
- [ ] Choose a single-cutter for step 0 that does not cut between the arms
- [ ] Decide A or B — or order both probes and settle the gap-length question in one
      experiment
