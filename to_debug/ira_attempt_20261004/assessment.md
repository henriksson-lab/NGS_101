# Assessment of the 2026-10-04 Brunello attempt

Source material in this folder (`.pdf`/`.docx` converted to `.txt` alongside):

| file | what it is |
|---|---|
| `broadgpp-sequencing-protocol.txt` | Broad GPP "PCR of sgRNAs for Illumina sequencing" — the one-step protocol |
| `new_design_padlock.txt` | three padlock designs: **MS1/MS2** for lentiCRISPRv2, **MS3** for lentiGuide-Puro |
| `email.txt` | the two-step PCR primer pair being considered |

Everything below is computed against the real maps in `../` — see [`../ira1.md`](../ira1.md)
§8 for why the Brunello pool identity is the crux.

## Verdict in one line

🟢 **The MS3 padlock design is correct and should work on Brunello #73178.**
⚠ **The second PCR builds a TruSeq *Small RNA* library, not a standard TruSeq one** — whether
that reads on a given run depends on the instrument's primer mix, which is **unverified**
(see §2). 🔴 **And the primers labelled "Brunello" in the workbook are
the lentiCRISPRv2 ones**, which are wrong for #73178.

---

## 1. ✅ MS3 padlock — this should work

| | MS1/MS2 (lentiCRISPRv2) | **MS3 (lentiGuide-Puro)** |
|---|---|---|
| extension arm | `GTGGAAAGGACGAAACACC` | `GTGGAAAGGACGAAACACC` (same) |
| | 19 nt, Tm 54.1 °C | 19 nt, Tm 54.1 °C |
| ligation arm | `AGCTAGGTCTTGAAAGGAGTGGG` | **`CGTAACTAGATCTTGAGACAAATGGCA`** |
| | 23 nt, Tm 58.5, ΔTm +4.4 | 27 nt, Tm 57.7, **ΔTm +3.7** |
| capture on lentiCRISPRv2 | ✅ gap 111 | 🔴 none (lig arm absent) |
| capture on lentiGuide-Puro | 🔴 none (lig arm absent) | ✅ **gap 110** |
| circle | 236 nt | **239 nt** |
| final library | ✅ **267 bp** — doc says 267 | ✅ **272 bp** — doc says 272 |

🟢 Both final library sizes reproduce the design document **exactly**, which is a strong
end-to-end check: arms, backbone, UMI, index and both PCR primers all have to be right for
that number to come out.

✅ The arm design is sound on its own terms too — the ligation arm is the hotter of the two
by **+3.7 °C**, mirroring the published probe's +4.4 °C, which is the ordering a padlock
wants (the ligation arm must stay annealed while the polymerase fills toward it).

✅ The two designs are **correctly mutually exclusive**: MS3 cannot capture on lentiCRISPRv2
and MS1/MS2 cannot capture on lentiGuide-Puro. So the probe is now a positive identifier of
the backbone, not just a reagent.

🟡 The captured region is **110–111 nt**, against **112 nt** for the published probe — i.e.
essentially the validated gap length, so no new gap-length risk is being taken. (The doc's
"111" counts the Pol III +1 G; the code's "110" counts strictly between the arms. Same
molecule, different convention — do not "fix" either.)

🟡 The captured sequence is written with the **original** scaffold
(`GTTTTAGAGCTAGAAATAGCAAGTT…`), which is correct for Brunello/lentiGuide-Puro and must
**not** be replaced with the AU-flip variant — that belongs to the Schmierer vector
(`../ira1.md` §7).

---

## 2. 🔴 The email's second PCR uses the wrong Illumina chemistry

```
PCR2 fwd  AATGATACGGCGACCACCGAGATCTACAC GTTCAGAGTTCTACAGTCCGACGATC ttgtggaaaggacgaaacaccg
                     P5                 ^^^^^^ TruSeq SMALL RNA RA5 ^^^^^^
PCR2 rev  CAAGCAGAAGACGGCATACGAGAT TGCCGA GTGACTGGAGTTCCTTGGCACCCGAGAATTCCA tctactattc...
                     P7              i7           ^^^^ revcomp of SMALL RNA RA3 ^^^^
```

✅ Checked explicitly:

| | present? |
|---|---|
| standard TruSeq Read 1 (`ACACTCTTTCCCTACACGACGCTCTTCCGATCT`) | 🔴 **no** |
| standard TruSeq Read 2 (`GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT`) | 🔴 **no** |
| TruSeq **Small RNA** RA5 (`GTTCAGAGTTCTACAGTCCGACGATC`) | ✅ yes |
| revcomp of TruSeq **Small RNA** RA3 (`CCTTGGCACCCGAGAATTCCA`) | ✅ yes |

✅ **The construct is a faithful TruSeq Small RNA library, not a chimera.** The P7 end is
exactly `P7 + i7 + GTGACTGGAGTTCCTTGGCACCCGAGAATTCCA`, which is the authentic small-RNA stub
(it shares only a 13-nt 5' stub with the standard Read 2 primer before diverging). So this is
a deliberate choice of chemistry, not a construction error.

✅ What is certain: **the standard TruSeq Read 1 primer cannot prime on it.** Only its last
6 nt (`CGATCT`) appear anywhere in the amplicon; at 8 nt and beyond there is no match, so
there is no 3' end for it to extend from. Read 1 must be primed off **RA5**.

🔴 **What is NOT established — and I previously overstated this:** whether a standard Illumina
run supplies a small-RNA-compatible primer. Illumina reagent cartridges have historically
carried *mixes* covering several library types, so a standard run may well read this library
without any custom spike-in. I could not verify it — Illumina's documentation site has
restructured (the URL the design document itself cites now 404s) and this session's web-search
budget is exhausted.

**Action:** ask the facility which Read 1 / index primers their run supplies, and whether the
cartridge mix covers TruSeq Small RNA. Also check how the i7 index read is primed — the 6-nt
`TGCCGA` sits between P7 and RA3, which is not where a standard TruSeq i7 sits. 🟡 This is
worth settling because it is cheap, **not** because it is established as the fault.

✅ The **first** PCR in the email is fine: `TCTTGTGGAAAGGACGAAACACCG` /
`TCTACTATTCTTTCCCCTGCACTGT` gives **229 bp** on cloned lentiGuide-Puro (and a 12.7 kb
wrap-around on lentiCRISPRv2, so it is also backbone-diagnostic).

---

## 3. ⚠ The P7 primer is vector-specific, and the naming in the workbook is a trap

🟢 The Broad protocol gives **two different P7 primers**, differing only in the 3' vector-
binding sequence:

| Broad P7 | 3' vector sequence | stated product |
|---|---|---|
| "for use with **lentiGuide**" | `TCTACTATTCTTTCCCCTGCACTGT` | 354 nt |
| "for use with **lentiCRISPRv2**" | `CCAATTCCCACTCCTTTCAAGACCT` | 285 nt |

✅ Reproduced against the maps (P5 1-nt stagger + the matching P7):

| | lentiGuide P7 | lentiCRISPRv2 P7 |
|---|---|---|
| on **lentiGuide-Puro** | ✅ **352 bp** | ⚠ 556 bp |
| on **lentiCRISPRv2** | ⚠ 12,796 bp | ✅ **287 bp** |

(352 vs Broad's 354 and 287 vs 285 — a 2-bp offset from where the cloning junction is drawn,
not a disagreement.)

> ⚠ **Both P7 binding sites exist in both vectors** — they are lentiviral backbone elements.
> So the wrong P7 does not simply fail; on lentiGuide-Puro it gives a **556 bp** product.
> A band is not evidence the right primer was used. **Only the size tells you.**

🔴 **And in `../Primers combined from Martin.xlsx` the naming is inverted relative to intent:**

| workbook primer | 3' sequence | actually targets |
|---|---|---|
| `Brunello-P7-iA01…iA11` | `CCAATTCCCACTCCTTTCAAGACCT` | **lentiCRISPRv2** |
| `CRISPR_PCR2-R-A05…B10` | `TCTACTATTCTTTCCCCTGCACTGT` | **lentiGuide-Puro** |

The workbook's own note already flags the asymmetry — *"the P5 primers for CRISPR PCR2 are
absolutely the same as for Brunello library PCR check, HOWEVER, the P7 primers are
DIFFERENT"* — but the primers **named after Brunello encode the lentiCRISPRv2 backbone**.
🔴 If the Brunello stock is **#73178 (lentiGuide-Puro)**, the "Brunello" primers are the wrong
ones and give 556 bp instead of 352 bp.

---

## 4. 🟡 Note for the GC-bias work: the Broad readout is Ex Taq

🟢 From the protocol: **Ex Taq DNA polymerase (Clontech RR001A)**, **28 cycles** in a single
step (95 °C 30 s / 53 °C 30 s / 72 °C 30 s), **up to 10 µg gDNA or 200 ng plasmid** per 100 µl,
≥4 parallel reactions per sample, 8 staggered P5 primers pooled for flowcell diversity.

⚠ Ex Taq is **Taq-based and non-proofreading** — the GC-biased end of the range, unlike the
KAPA HiFi / Q5 used by every dataset in
[`../../gcbias/datasets/polymerases.md`](../../gcbias/datasets/polymerases.md). A 28-cycle
single-step Taq amplification is the most GC-bias-prone readout encountered so far, and it is
the protocol the client's own workflow is based on. Added to that table.

---

## What to do, in order

1. 🔴 **Settle the sequencing chemistry.** Ask the facility, in writing, whether the run uses
   TruSeq **Small RNA** primers. If not, PCR2 must be rebuilt with standard TruSeq tails. This
   is the cheapest check and the most likely total-failure cause.
2. 🔴 **Get the Addgene pool number** (#73178 vs #73179). One number decides which P7 and
   which padlock.
3. ✅ **Size the one-step Broad PCR**: 352 bp → lentiGuide-Puro (#73178), 287 bp →
   lentiCRISPRv2 (#73179), 556 bp → lentiGuide amplified with the wrong (lentiCRISPRv2) P7.
4. ✅ **Then order MS3** if the pool is #73178 — the design checks out end to end.
5. 🟡 If GC bias matters for the readout, consider replacing Ex Taq with KAPA HiFi; see the
   polymerase table.
