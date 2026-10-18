# DMS-seq

**Defining source:** Rouskin *et al.*, *Nature* (2014), doi:10.1038/nature12894.
The PMC full text and Methods were read.

## Evidence used

- 🟢 DMS modifies solvent-accessible A and C and blocks reverse transcription.
- 🟢 Methods specify Zn²⁺ fragmentation, 60–80-nt RNA selection, PNK end repair,
  modified linker-1, 25–50-nt stopped-cDNA selection, CircLigase, 8–10 PCR cycles and HiSeq 2000.
- 🟢 In-vivo, refolded in-vitro and denatured controls distinguish cellular structure and bias.
- 🟡 Selected study-specific PCR and sequencing-primer geometry remain unassigned; the modified predecessor does not establish those choices.

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

## Local validation migration (2026-10-18)

Re-read PMC3966492 Methods text lines 48–58 and ran the source-bundle oligo
scraper. The graph follows the in-vivo 30 °C yeast branch: DMS treatment now precedes
quenching, RNA extraction and poly(A) selection. It explicitly retains the 60–80-nt
RNA gel selection, post-ligation gel purification and stopped-cDNA band 25–50 nt
above RT-primer size, with their physically excluded fractions.

The linker reagent includes its terminal dideoxyC nucleotide as well as the printed
5′ rApp/3′ ddC end modifications. Its single-strand donor geometry and the RNA–linker
junction are explicit. SuperScript III, CircLigase and eight to ten PCR cycles are
local declarations. The actual selected RT/PCR members, circle junction, indexing
status and read acquisition remain unassigned; a modified predecessor method does
not establish those study-specific choices. PCR uses its incoming circle boundary,
without substituting the expected library as the substrate.

Verification: source-exact modified-linker test, native graph/page, source guard,
pairing and scoped chemistry pass. This migration leaves reaction relations and
the remaining source targets open; it is not a full chemistry or absence audit.
