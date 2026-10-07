# Tn5 tagmentation (Nextera)

The single most reused adapter chemistry in NGS. Implemented in `lib/nextera.py`.

## What it does

Tn5 transposase is a cut-and-paste transposase. *In vitro*, a **hyperactive Tn5 mutant** is
loaded with short synthetic adaptors instead of a transposon, forming a dimer. Each monomer
carries one adaptor; the dimer binds DNA, cuts both strands, and covalently joins its
adaptor's 3' end to the 5' end of the cut. Fragmentation and adapter addition happen in the
same step — hence *tagmentation*.

This replaces two steps (mechanical/enzymatic shearing, then ligation) with one, which is
why it wins on input amount and hands-on time, and why it is the back end of so many
single-cell methods.

## The sequences

```
Mosaic end (ME), the Tn5 binding site   AGATGTGTATAAGAGACAG      19 bp
what you read through into (revcomp)    CTGTCTCTTATACACATCT
s5  N/S5xx primer entry point           TCGTCGGCAGCGTC           14 nt
s7  N7xx primer entry point             GTCTCGTGGGCTCGG          15 nt

loaded adaptor = entry point + ME
```

Index primers rebuild the flowcell ends by PCR:

```
N/S5xx   AATGATACGGCGACCACCGAGATCTACAC [8-bp i5] TCGTCGGCAGCGTC
N7xx     CAAGCAGAAGACGGCATACGAGAT      [8-bp i7] GTCTCGTGGGCTCGG
```

Sequencing primers fall out of the same two pieces:

| Read | Primer |
|---|---|
| Read 1 | `s5 + ME` |
| Read 2 | `s7 + ME` |
| Index 1 (i7) | `revcomp(ME) + revcomp(s7)` |
| Index 2 (i5) | `revcomp(ME) + revcomp(s5)` |

## Two details that trip people up

**The 9-bp gap.** The two Tn5 monomers cut the two strands **9 bp apart**, so after
tagmentation each end has a 9-nt single-stranded gap. Only the transferred strand is joined;
the other is not. This is why every Nextera PCR programme begins with a hold at **72 °C
before any denaturation** — that step is a gap fill-in, not an extension, and skipping it
loses the library.

**Only one product in three amplifies.** A reaction loaded with both adaptors inserts them at
random, so fragments end up with:

| Ends | Amplifies? | Why |
|---|---|---|
| s5 / s5 | no | P5 at both ends; no P7 site |
| s7 / s7 | no | P7 at both ends; no P5 site |
| **s5 / s7** | **yes** | the only product that can bridge P5 to P7 |

This is suppression by design, not an inefficiency — and it is why the scg_lib_structs pages
draw every product and label which survives. When one end of the molecule is *not* tagmented
(e.g. the 3' end of a barcoded cDNA), the product table grows accordingly; SPLiT-seq draws
five.

## Where it shows up

SMART-seq / SMART-seq2 / SMART-seq3 (their separate `smart-seq*__/` protocol directories), SPLiT-seq, ATAC-seq and scATAC, Nextera
DNA Flex, and most "tagmentation-based" single-cell methods. The front ends differ
completely; the back end is this.

## References

- Reznikoff WS. *Transposon Tn5*. Annu Rev Genet 2008;42:269–286. doi:[10.1146/annurev.genet.42.110807.091656](https://doi.org/10.1146/annurev.genet.42.110807.091656)
- Picelli S *et al.* *Tn5 transposase and tagmentation procedures for massively scaled sequencing projects*. Genome Res 2014;24(12):2033–40. doi:[10.1101/gr.177881.114](https://doi.org/10.1101/gr.177881.114)
- Illumina Adapter Sequences, Document # 1000000002694
