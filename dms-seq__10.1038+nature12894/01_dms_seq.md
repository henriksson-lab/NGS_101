# DMS-seq

**Defining source:** Rouskin *et al.*, *Nature* (2014), doi:10.1038/nature12894.
The PMC full text and Methods were read.

## Evidence used

- 🟢 DMS modifies solvent-accessible A and C and blocks reverse transcription.
- 🟢 Methods specify Zn²⁺ fragmentation, 60–80-nt RNA selection, PNK end repair,
  modified linker-1, 25–50-nt stopped-cDNA selection, CircLigase, 8–10 PCR cycles and HiSeq 2000.
- 🟢 In-vivo, refolded in-vitro and denatured controls distinguish cellular structure and bias.
- 🟡 PCR and sequencing-primer sequences remain role-level because the paper does not print them.

## Cited library-primer source

🟢 The DMS-seq Methods explicitly call the library preparation a modified version of Ingolia *et al.*'s ribosome-profiling protocol, [doi:10.1038/nprot.2012.086](https://doi.org/10.1038/nprot.2012.086). That protocol (PMC3535016, “Reagents”) prints the same `/5rApp/CTGTAGGCACCATCAAT/3ddC/` linker and a circularization-compatible RT primer:

`5′-(Phos)-AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGC-(SpC18)-CACTCA-(SpC18)-TTCAGACGTGTGCTCTTCCGATCTATTGATGGTGCCTACAG-3′`

It also prints a forward library PCR primer, `5′-AATGATACGGCGACCACCGAGATCTACAC-3′`, and indexed reverse-primer pattern, `5′-CAAGCAGAAGACGGCATACGAGATNNNNNNGTGACTGGAGTTCAGACGTGTGCTCTTCCG-3′` (the source separates `AC GTG` typographically). The N run is the reverse complement of the Illumina index.

🟡 These are exact **Ingolia** reagent sequences, not direct evidence that Rouskin *et al.* used every base or index unchanged in DMS-seq; the DMS-seq paper says “modified version” and prints only the linker. Until direct DMS-seq oligo evidence is found, the RT and PCR molecules remain inferred candidates and the published DMS-seq schematic should keep their sites unresolved.

## Construct logic

The library preserves a termination coordinate rather than a full RNA fragment. A read end
immediately before A or C reports a DMS adduct; relative counts measure accessibility.

## Bounded curation round 4 (2026-10-16)

Inspected: Published PMC3966492 JATS42–60, DMS modification, Library Generation and sequencing.

Changed T4PNK end treatment to 3-prime dephosphorylation rather than phosphorylation; named SuperScriptIII52°C RT and 8–10 PCR cycles. Replaced unsupported broad primer-absence wording with unassigned model geometry.

Deferred: Original predecessor library protocol and source primer list, explicit 60–80nt RNA/25–50nt-above-primer cDNA gel selections, RNA recovery and exact read acquisition. This pass does not establish a complete source review or publication-wide absence.
