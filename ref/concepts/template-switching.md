# Template switching and the TSO

How a defined handle gets onto the 5' end of a cDNA without any ligation.
Implemented in `lib/rt.py`.

## The mechanism

MMLV-family reverse transcriptases have a **terminal transferase** side activity. When the
enzyme reaches the 5' end of an mRNA template and runs out of template, it does not simply
stop — it adds a few **untemplated nucleotides** to the 3' end of the new cDNA strand,
predominantly **C**, usually three.

A **template-switching oligo (TSO)** ending in **`rGrGrG`** base-pairs with that `CCC`
overhang. The polymerase treats the annealed TSO as a new template, switches to it, and
copies it. The cDNA now carries a known sequence at its 5' end.

```
                                      MMLV adds CCC here
                                              |
 mRNA  5'- XXXXXXXXXXXXXXXXXXXXXXXXXXXX(A)n   v
 cDNA              3'- CCC XXXXXXXXXXXXXXXXNV(T)30 [oligo-dT handle] -5'
                       |||
 TSO       5'- [handle] GGG -3'    ---> polymerase switches template and copies the handle
```

**Why ribo-G.** The contact is only three base pairs, which is far too weak to be useful as
plain DNA. The terminal G's are **ribonucleotides** (`rG`), or in SMART-seq2 the last one is
a **locked nucleic acid** (`+G`). Both raise the duplex stability of that short G:C contact
enough to make switching efficient — the LNA substitution was one of the two changes that
defined SMART-seq2.

## Why it matters

- **Full-length coverage without ligation.** Both ends of the transcript end up flanked by
  defined sequence, so the whole cDNA can be amplified by PCR.
- **It marks the 5' end specifically.** Template switching only happens where the polymerase
  ran off the end of the template — so anything the TSO carries (a tag, a UMI) ends up
  attached to a genuine transcript 5' end, not an internal position. SMART-seq3 exploits
  this directly.
- **It is a counting opportunity.** A UMI placed in the TSO labels an individual mRNA
  molecule *before* amplification.

## Two design patterns

**Symmetric (SMART-seq, SMART-seq2).** The oligo-dT primer and the TSO carry the *same*
handle, so after amplification both ends of the cDNA are identical and a **single primer**
(ISPCR) amplifies the whole thing. Simple, but this is also *semi-suppressive* PCR, and
there is nothing to distinguish a 5' fragment from an internal one.

**Asymmetric (SMART-seq3).** The oligo-dT and the TSO carry *different* handles, so PCR uses
two different primers. The TSO additionally carries an **11-bp tag** and an **8-bp UMI**.
A read containing the tag is known to come from the transcript's 5' end and carries a
molecular count; reads without it are used for coverage only. That single change is what lets
one protocol do both full-length coverage and UMI counting.

## The shared handle

```
SMART / ISPCR handle    AAGCAGTGGTATCAACGCAGAGT     23 nt
```

Reused by SMART-seq, SMART-seq2, Drop-seq, SPLiT-seq and the Takara SMARTer kits. If a
protocol mentions an "ISPCR primer" or a "SMART handle", it is this sequence.

## References

- Zhu YY *et al.* *Reverse transcriptase template switching: a SMART approach for full-length cDNA library construction*. Biotechniques 2001;30(4):892–7. doi:[10.2144/01304pf02](https://doi.org/10.2144/01304pf02)
- Picelli S *et al.* *Smart-seq2 for sensitive full-length transcriptome profiling in single cells*. Nat Methods 2013;10(11):1096–8. doi:[10.1038/nmeth.2639](https://doi.org/10.1038/nmeth.2639)
- Hagemann-Jensen M *et al.* *Single-cell RNA counting at allele and isoform resolution using Smart-seq3*. Nat Biotechnol 2020;38(6):708–714. doi:[10.1038/s41587-020-0497-0](https://doi.org/10.1038/s41587-020-0497-0)
