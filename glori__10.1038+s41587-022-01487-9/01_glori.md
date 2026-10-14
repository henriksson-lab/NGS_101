# GLORI

**Defining source:** Liu *et al.*, *Nature Biotechnology* (2023),
doi:10.1038/s41587-022-01487-9. **Detailed protocol:** Shen *et al.*, *Nature
Protocols* (2024), doi:10.1038/s41596-023-00937-1. The methods, public data record
and protocol description were read.

## Evidence used

- 🟢 Glyoxal protects G and catalyses nitrite conversion of unmodified A to inosine;
  m6A resists conversion and remains A.
- 🟢 The detailed procedure specifies fragmentation, glyoxal protection, 750 mM nitrite
  conversion at 16 °C for 8 h plus 4 °C overnight, heat deprotection and end repair.
- 🟢 Defining libraries used modified eCLIP with a 10-nt inline UMI and paired-end HiSeq
  X Ten sequencing; NEBNext Small RNA and SMARTer alternatives were also demonstrated.
- 🟢 The defining article's Supplementary Data 8 (`41587_2022_1487_MOESM9_ESM.xlsx`)
  publishes sequences for synthetic m6A spike-in **RNA templates** and other assay
  controls. These are not the adapters or PCR primers that form the sequencing library.
- 🟢 The detailed protocol's [Supplementary Tables 1–3](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41596-023-00937-1/MediaObjects/41596_2023_937_MOESM3_ESM.xlsx) give the **3′ RNA adapter, RT primer and 3′ cDNA adapter with a 10-nt UMI**. Table 3's RT/PCR/sequencing primers are specifically for conversion evaluation of synthetic Probe 5, not the complete RNA-seq library.
- 🟡 Complete modified-eCLIP library amplification oligos and ordered sequencing-primer oligos are still not reproduced in the accessible sources; their final roles remain placeholders.

## Oligo source boundary

🟢 Detailed protocol Supplementary Table 1 prints, as ordered: 3′ RNA adaptor
`rAppAGATCGGAAGAGCGTCGTG`, RT primer `ACACGACGCTCTTCCGATCT`, and 3′ cDNA adaptor
`Phos-NNNNNNNNNNAGATCGGAAGAGCACACGTCTG`. The `N` span is the 10-nt UMI. These are
the 2024 detailed protocol's oligos; the original 2023 experiment's precise lot or
version equivalence is not separately documented.

🔴 The complete modified-eCLIP library PCR oligos remain unavailable. The detailed
protocol names NEBNext Small RNA and SMARTer kit alternatives, but the defining branch is
modified eCLIP; those kit oligos cannot be substituted for it. Raw FASTQ from GSE210563
could support a read-visible adapter consensus and UMI boundary, but cannot establish an
ordered oligo's end modifications or unsequenced tail. Keep only the unverified outer
library arms inferred.

## Construct logic

GLORI is a negative conversion assay: an unmodified genomic A should be read as G, while
the retained-A fraction estimates m6A stoichiometry. The model therefore keeps chemical
state distinct from the downstream library branch and carries UMI identity explicitly.
