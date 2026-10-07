# SMART-seq3, SMART-seq3xpress and FLASH-seq — quick reference

> **Evidence marking.** 🟢 verbatim from the source · 🟡 derived or inferred · 🔴 not published.

**SMART-seq3** — Hagemann-Jensen M, Ziegenhain C, Chen P, Ramsköld D, Hendriks G-J,
Larsson AJM, Faridani OR, Sandberg R. "Single-cell RNA counting at allele and isoform
resolution using Smart-seq3." *Nat Biotechnol* 2020;**38**(6):708–714.
doi:[10.1038/s41587-020-0497-0](https://doi.org/10.1038/s41587-020-0497-0)

**SMART-seq3xpress** — Hagemann-Jensen M, Ziegenhain C, Sandberg R. *Nat Biotechnol*
2022;**40**(10):1452–1457. doi:[10.1038/s41587-022-01311-4](https://doi.org/10.1038/s41587-022-01311-4)

**FLASH-seq** — Hahaut V, Pavlinic D, Carbone W, *et al.* *Nat Biotechnol*
2022;**40**(10):1447–1451. doi:[10.1038/s41587-022-01312-3](https://doi.org/10.1038/s41587-022-01312-3)

Structure transcribed from
<https://teichlab.github.io/scg_lib_structs/methods_html/SMART-seq_family.html>.

---

## 1. What changed, and why it matters

🟢 SMART-seq3 is a redesign of the oligos, not an optimisation. Four changes:

1. **Maxima H− reverse transcriptase** instead of MMLV — RNase H−, more thermostable, more processive.
2. **NaCl instead of KCl** during reverse transcription.
3. **5% PEG** as a crowding reagent, as in mcSCRB-seq.
4. **An 11-bp tag and an 8-bp UMI in the TSO**, so that a 5' read can be distinguished from
   internal reads.

Point 4 is the architectural one. SMART-seq/2 made both ends of the cDNA identical; SMART-seq3
deliberately makes them **different**, which buys three things at once:

- **PCR is no longer single-primer** — a forward and a reverse primer, each matching one end.
- **A read carrying the 11-bp tag is known to come from the transcript's 5' end**, because
  template switching only happens where the polymerase ran off the template.
- **The UMI counts molecules**, since it was attached before any amplification.

So one protocol gives **full-length coverage *and* UMI counting at allele and isoform
resolution** — which is the headline result.

🟢 **SMART-seq3xpress** and **FLASH-seq** do not change the library architecture. They cut
cost and time by removing intermediate purification steps and dropping to nanolitre reaction
volumes. The library generation procedure is the same and the final libraries are *almost* the
same: the only difference is a different TSO variant, which makes the 5' fragment slightly
different.

## 2. Oligos

```
Smartseq3_OligodT30VN  5'- /5Biosg/ ACGAGCATCAGCAGCATACGA TTTTT...TTTTT VN -3'   (30 T)

Smartseq3_N8_TSO       5'- /5Biosg/ AGAGACAGATTGCGCAATG [8-bp UMI]       rGrGrG -3'
Smartseq3xpress_TSO    5'- /5Biosg/ AGAGACAGATTGCGCAATG [8-bp UMI] WW    rGrGrG -3'
FLASH-seq_TSO          5'- /5Biosg/ AGAGACAGATTGCGCAATG [8-bp UMI] CTAAC rGrGrG -3'

Fwd_PCR_primer         5'- TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGATTGCGCAA*T*G -3'
Rev_PCR_primer         5'- ACGAGCATCAGCAGCATAC*G*A -3'
```

`*` is a phosphorothioate linkage (3'-exonuclease protection). `/5Biosg/` is 5' biotin.

🟡 Note the TSO handle: `AGAGACAG` is the **3' end of the Nextera mosaic end**, and
`ATTGCGCAATG` is the 11-bp tag. So the forward PCR primer
(`s5 + ME + tag` = `TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG` + `ATTGCGCAATG`) rebuilds a complete
Nextera s5 arm as it amplifies — the TSO is already most of one.

The three family members differ **only** in the spacer between UMI and G-tail: none,
`WW`, or `CTAAC`.

## 3. Steps

| # | Step | Note |
|---|---|---|
| 1 | Anneal Smartseq3_OligodT30VN; reverse transcribe with **Maxima H− MMLV** | NaCl, 5% PEG |
| 2 | Terminal transferase adds untemplated **CCC** | |
| 3 | TSO anneals to the CCC; template switch copies tag + UMI onto the cDNA | the UMI is now attached, pre-amplification |
| 4 | PCR with **Fwd_PCR_primer and Rev_PCR_primer** | two different primers — not single-primer |
| 5 | Purify amplified cDNA | |
| 6 | **Nextera tagmentation** | 9-bp gap |
| 7 | **72 °C gap fill-in** | first cycle of the Nextera PCR |
| 8 | PCR with N/S5xx and N7xx index primers | |

## 4. Two kinds of fragment

This is the point of the design, and it is why the source page draws the final structure twice.

**5' fragments** retained the TSO, so they carry the 11-bp tag and the UMI:

```
P5 [i5] s5 ME ATTGCGCAATG [8-bp UMI] GGG --- cDNA --- revcomp(ME) revcomp(s7) [i7] revcomp(P7)
```

**Internal fragments** did not, and are **character-for-character identical to a SMART-seq2
fragment** (asserted in `tools/selftest.py`). They carry no UMI and no positional information,
and are used for coverage only.

🟡 So the demultiplexer's job is simply: does read 1 begin with the 11-bp tag? If yes, the
read pair is 5'-derived and the following 8 nt are its UMI.

Single cells are identified by the **i5 + i7 index combination**, as in SMART-seq2 — there is
still no in-molecule cell barcode.

## 5. Back end

Identical to SMART-seq2 and to every other Nextera method — see
[`../ref/concepts/tn5-tagmentation.md`](../ref/concepts/tn5-tagmentation.md). The sequencing
primers are unchanged; only read 1's *content* differs, depending on whether that fragment
kept the TSO.
