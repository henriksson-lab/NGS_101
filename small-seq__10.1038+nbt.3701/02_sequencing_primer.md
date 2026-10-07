# Which sequencing primer reads a Small-seq library?

> **Evidence marking.** 🟢 verbatim from the source · 🟡 derived or inferred · 🔴 not
> published · ✅ computed and asserted in `tools/selftest.py`.

Short answer: **the paper never says.** It does not use the words *sequencing primer*
anywhere, and it does not state whether a custom primer has to be supplied. Everything it
does say about sequencing is quoted below in full, so that nobody has to re-read it.

This note exists because the question is live elsewhere in the repo: Small-seq is built on
**TruSeq *Small RNA*** adapters, whose read-1 landing site is not the TruSeq DNA one that
every other protocol here uses.

---

## 1. Everything the paper says about sequencing — verbatim

🟢 Abstract: *"Ready-to-sequence libraries can be generated in 2–3 d, starting from cell
collection"*, and the method *"relies on standard reagents and instruments."*

🟢 *Overview of the protocol* (p. 2408): *"Finally, the libraries are pooled and purified,
**ready for sequencing on an Illumina sequencer such as a HiSeq or NextSeq platform**."*

🟢 Fig. 1 legend: *"Sequences required for Illumina cluster generation and sample indexing
are added through two rounds of PCR."*

🟢 *Experimental design — RT, PCR amplification, and sequencing* (p. 2411): *"**Sequencing
is performed using Illumina platforms.** We normally sequence 1–2 million reads per cell.
Sequencing all 192 barcoded samples (two 96-well plates) in an Illumina NextSeq lane
produces roughly 2 million reads per cell."*

🟢 *Experimental design — Ligation of adapters* (p. 2410): *"**Sequencing starts from the
UMI**, which is immediately followed by a '-CA-' couple, which separates the UMI from the
small RNA sequences (Fig. 1)."*

🟢 *Equipment* (p. 2413): *"Illumina DNA sequencer (e.g., HiSeq, MiSeq, NovaSeq, or
NextSeq)."*

🟢 Step 37: *"Perform a pilot sequencing experiment with a pool of 48 samples, using a
MiSeq or a NextSeq instrument. Optimize the size-selection window or sequencing depth that
results in the highest number of small RNA molecules or genes of interest in the final data.
**Carry out sequencing according to the manufacturer's protocol**."*

🟢 Step 38: *"Perform sequencing for a single read of 51 bp or longer, depending on the
small RNA molecules of interest. Sequencing of 51 bp is enough to detect miRNAs."*

🟢 *Reagent setup — SRX* (p. 2414): *"SRX DNA index primers are **modified from standard
Illumina TruSeq small RNA index primers** to have barcodes of 8 bp for 192 samples. For
sequence information, see Supplementary Table 1."*

That is the complete set.

## 2. What is *not* stated

🔴 The strings **"sequencing primer"**, **"custom primer"**, **"primer mix"**, **"Read 1"**,
**"flow cell"**, **"spike-in"** and **"PhiX"** do not occur anywhere in the article — zero
occurrences each. Nor is there any instruction to load a primer into a reagent cartridge
position, any note on which platforms do or do not need one, and any mention of the index
read beyond "sample indexing".

🔴 Supplementary Table 1 — the one place an oligo list could have settled it — is not in
the extracted text. **It may well contain a sequencing primer.** The *Materials* list does
not: 🟢 *"Primers and adapters (5.8S rRNA-masking oligonucleotide, RA3, RA5, RTP, SRX, and
RP1 …)"* — five oligos plus the mask, and none of them is a sequencing primer.

🟡 **Reading the omission.** Taken with 🟢 *"relies on standard reagents and instruments"*
and 🟢 *"Carry out sequencing according to the manufacturer's protocol"*, the authors
evidently did not consider the primer a protocol-level decision. That is consistent with
their having used stock TruSeq Small RNA sequencing on a stock run — but it is an argument
from silence, and the paper does not license any stronger claim.

## 3. What the construct forces, regardless

The parts of the library the paper *does* print are enough to pin the primer down
structurally. ✅ All of the following are computed in `tools/selftest.py` from the
published sequences alone.

**(a) The read-1 primer's 3' end is fixed by the construct.** 🟢 "Sequencing starts from
the UMI", and ✅ the UMI begins at offset **55** of the library, immediately after
P5 + RA5. So the primer must extend to the last base of RA5 and no further — its 3'
terminus is RA5's 3' terminus, `…GTTCAGAGTTCTACAGTCCGACGATC`. ✅ Everything upstream of
that point is `P5 + RA5` and nothing else, so **any usable read-1 primer is a suffix of
`P5 + RA5`**.

```
5'- AATGATACGGCGACCACCGAGATCTACAC GTTCAGAGTTCTACAGTCCGACGATC HHHHHHHH CA  small RNA  …
    |------------- P5 -----------|------------ RA5 ----------|-- UMI --|
                                                             ^
                                      read 1 must start here (base 56)
```

**(b) The TruSeq *DNA* read-1 primer cannot do it.** ✅ `ACACTCTTTCCCTACACGACGCTCTTCCGATCT`
does not occur in the library in either orientation — nor does the TruSeq read-2 primer
sequence. There is no TruSeq DNA adapter anywhere in this construct. 🟡 So a run configured
only with the standard TruSeq DNA read-1 primer would produce nothing from these libraries,
and that is true of any TruSeq Small RNA library, not a Small-seq peculiarity.

**(c) RP1 is not the answer either.** ✅ `RP1 == P5 + RA5[:21]`: it stops **5 nt short** of
RA5's 3' end, the missing bases being `CGATC`. 🟡 Priming read 1 with RP1 would emit
`CGATC` before the UMI, contradicting 🟢 "Sequencing starts from the UMI". So the primer
used was *not* simply the PCR primer re-used — it reaches 5 bases further.

**(d) Small-seq's RA5 is a stock adapter with the UMI appended at its 3' end.** 🟡 RA5's
first 26 nt are the standard TruSeq Small RNA 5' adapter; the UMI and `CA` are *added
after* it, between adapter and insert. The direct consequence:

> 🟡 **Any primer that reads a stock TruSeq Small RNA library also reads a Small-seq
> library.** It anneals to the identical 26-nt region. The only difference is that it emits
> **10 extra bases — 8 UMI + `CA` — before the insert**, which is precisely the 10 nt the
> pipeline trims (🟢 `-u 2` in cutadapt, after UMI extraction with `NNNNNNNN`).

**(e) The index read is as unresolved as SRX.** 🟡 Whatever primes the i7 index read must
anneal between RA3 and the index — i.e. inside SRX's own linker, which is 🔴 unpublished.
Since SRX is 🟢 *"modified from standard Illumina TruSeq small RNA index primers"* and the
modification is described as being to the **barcode length only** (6 bp → 8 bp), 🟡 a stock
TruSeq Small RNA index read should work. The paper does not say so.

## 4. Conclusion, stated precisely

| Question | Answer |
|---|---|
| Does the paper name a sequencing primer? | 🔴 **No.** Zero occurrences of the term. |
| Does it say a custom primer is needed? | 🔴 **No.** |
| Does it say the standard instrument primers suffice? | 🔴 **No.** The closest is 🟢 *"ready for sequencing on an Illumina sequencer such as a HiSeq or NextSeq platform"*, 🟢 *"relies on standard reagents and instruments"* and 🟢 *"Carry out sequencing according to the manufacturer's protocol."* |
| Can a TruSeq **DNA** read-1 primer read it? | ✅ **No** — that site is absent from the construct. |
| Can a TruSeq **Small RNA** read-1 primer read it? | 🟡 **Yes**, with 10 extra bases (UMI + CA) before the insert, which the pipeline trims. |
| Which primer did the authors actually use? | 🔴 **Unknown.** Possibly in Supplementary Table 1, which is not in the extracted text. |

🟡 So the honest statement for the TruSeq-Small-RNA-compatibility question is: *Small-seq
does not change the read-1 landing site relative to a stock TruSeq Small RNA library; it
only pushes the insert 10 nt further into the read. Whatever primer arrangement a facility
uses for TruSeq Small RNA libraries applies unchanged.* The paper supports the construct
half of that claim and is silent on the primer half.

## 5. How to close this

1. **Get Supplementary Table 1** (doi 10.1038/s41596-018-0049-y). It is the only place in
   the paper that could name a primer, and it also holds the 🔴 SRX sequences.
2. **Add a TruSeq Small RNA block to `lib/illumina.py`** from the Illumina adapter-sequences
   document — RA5, RA3, RTP, RP1, RPI and the small-RNA sequencing primer — verbatim, as
   that file requires. Then point (d) above stops being 🟡: `RA5[:26]` either equals the
   vendor sequence or it does not, and the selftest can assert it.
3. **Check the GitHub pipeline** (<https://github.com/eyay/smallseq>) for the adapter FASTA
   it trims with (`adapters/cutadapt_3prime.fa`, referenced in 🟢 Step 41). That file pins
   what the authors expected to read through into at the 3' end.
