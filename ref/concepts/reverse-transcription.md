# Reverse transcription, and what changes between RNA and DNA protocols

The step that decides most of a protocol's downstream shape. Implemented in `lib/rt.py`.

## Priming

**Oligo-dT.** Anneals to the poly(A) tail. The 3'-most two bases are normally **anchored**
(`VN`, where V = A/C/G and N = anything) so the primer sits at the poly(A)/transcript
junction rather than sliding along the tail — without the anchor, the start position is
arbitrary and 3' bias worsens. Selects polyadenylated RNA, which for eukaryotes is mostly
what you want.

**Random priming (`NNNNNN`).** Needed when there is no poly(A): bacterial mRNA, rRNA-depleted
total RNA, degraded material. Also used *alongside* oligo-dT to recover internal and
non-polyadenylated sequence — SPLiT-seq uses half oligo-dT and half random hexamer wells for
exactly this reason.

**Gene-specific priming.** Targeted panels; not covered here.

## The enzyme matters

| Enzyme | Where it shows up | Why |
|---|---|---|
| MMLV (wild type) | SMART-seq, SMART-seq2 | has the terminal transferase activity template switching depends on |
| **Maxima H−** | SMART-seq3, SPLiT-seq | RNase H−, more thermostable, more processive; better full-length yield |
| SuperScript II/IV | many | similar niche |

RNase H− matters because the RNA strand of the RNA:DNA hybrid is not degraded mid-synthesis,
which helps on long transcripts.

## RNA protocols vs DNA protocols

This is the structural fork that determines nearly everything downstream.

| | **RNA (cDNA)** | **DNA (genomic)** |
|---|---|---|
| First step | reverse transcription | none — the template is already DNA |
| Getting a 5' handle | **template switching** (TSO), or ligation | ligation, or tagmentation |
| Getting a 3' handle | the primer carries it | ligation, or tagmentation |
| Amplification | PCR off the two handles | PCR, MDA, or **PTA** |
| Counting unit | mRNA molecules — a **UMI** is meaningful and usual | genome copies — a UMI is usually pointless, since each locus is present once |
| Coverage goal | per-transcript, 3'-biased or full-length | even across the genome |
| Typical artefact | 3' bias, template-switching chimeras | amplification bias, allelic dropout, chimeras |

Two consequences worth holding onto:

- **UMIs belong to RNA protocols.** They count molecules that existed before amplification.
  In single-cell WGS each locus starts at one or two copies, so there is nothing to count —
  which is why `atrandi-wgs__10.1101+2025.06.20.660799/` has no UMI, and why its read structure has four barcodes and
  no UMI field.
- **The 5' end is the hard end for RNA, and both ends are hard for DNA.** RNA protocols get
  the 3' handle for free (the primer carries it) and solve the 5' end with template
  switching. DNA protocols have no primer to hang a handle on, so they use ligation (with
  its end-repair and dA-tailing prerequisites) or tagmentation (which installs both ends at
  once).

## References

- Picelli S *et al.* Nat Methods 2013;10(11):1096–8. doi:[10.1038/nmeth.2639](https://doi.org/10.1038/nmeth.2639)
- Hagemann-Jensen M *et al.* Nat Biotechnol 2020;38(6):708–714. doi:[10.1038/s41587-020-0497-0](https://doi.org/10.1038/s41587-020-0497-0)
- Rosenberg AB *et al.* *Single-cell profiling of the developing mouse brain and spinal cord with split-pool barcoding*. Science 2018;360(6385):176–182. doi:[10.1126/science.aam8999](https://doi.org/10.1126/science.aam8999)
