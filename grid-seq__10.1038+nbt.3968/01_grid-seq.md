# GRID-seq

## Purpose

GRID-seq records a chromatin-associated RNA and a nearby genomic DNA end in one sequenceable molecule by joining both to a directional, biotinylated bridge.

## Sources read

- 🟢 Li et al., *GRID-seq reveals the global RNA–chromatin interactome*, Nature Biotechnology (2017), [doi:10.1038/nbt.3968](https://doi.org/10.1038/nbt.3968), especially Methods and Extended Data Fig. 1.
- 🟢 The open full text (PMCID `PMC5953555`) prints the two linker strands, their modifications, ligases, AluI digestion, affinity purification, MmeI tag release, adapter sequences, PCR primers and sequencing primer.
- 🟢 Zhou et al., *GRID-seq for comprehensive analysis of global RNA–chromatin interactions*, Nature Protocols (2019), [doi:10.1038/s41596-019-0172-4](https://doi.org/10.1038/s41596-019-0172-4), PMCID `PMC7721247`: Table 1 reprints the full linker, Y-adapter and PCR oligos; Methods identify MmeI's two-base 3′ overhang and the adapter's matching `NN` end (lines 106–122), and step 151 specifies single-end 100-base sequencing through the linker.

## Construct-changing path

🟢 The linker DNA strand is `5′-/5Phos/GTTGGAGTTCGGTGTGTGGGAGTGAGCTGTGTC-3′`. The partner is a mixed RNA/DNA strand, `5′-/5Phos/rGrUrUrGrGrArUrUrCrNrNrNrGrACACAGC/iBiodT/CACTCCCACACACCGAACTCCAAC-3′`, pre-adenylated before RNA ligation. Its printed `rGrA` after `rNNN` are the first two bases of the arm that pairs with L2 (the reverse complement of L2 begins `GACACAGC`), so the linker is 9 + 3 + 33 = {{= grid_seq.LINKER_LENGTH }} bp, as the 2019 protocol states ("45-bp linker and 20–21-nt RNA and DNA at each end", PMC7721247 JATS line 90). The `rNNN` bases identify PCR duplicates; the internal biotin supports capture.

🔴 **Published linker discrepancy.** The 2017 defining paper, Methods “Construction of GRID-seq library” (PMC5953555 HTML text line 1272), prints the mixed strand ending `.../iBiodT/CACTCCCACACACCGAACTCCAAC`. The same group's 2019 detailed protocol, Table 1 (PMC7721247 JATS line 1781), prints `.../iBiodT/CACTCCACACACCGAACTCCAAC`, one C shorter. Both print L2 as `GTTGGAGTTCGGTGTGTGGGAGTGAGCTGTGTC` (2017 line 1272; 2019 line 1783). Its reverse complement matches the 2017 mixed strand's DNA duplex arm exactly, so the model uses that pairing-compatible 2017 value. The shortened 2019 strand's synthesis/use is not established by these texts.

🟢 Fixed, SDS-permeabilized nuclei are cut with AluI (37 °C, 2 h), then treated with T4 PNK and ATP (2017 Methods, PMC5953555 JATS text line 132). The 2019 protocol (PMC7721247 JATS line 114) states that PNK phosphorylates 5′ ends and dephosphorylates 3′ ends of both DNA and RNA. ATP-free T4 RNA ligase 2 truncated KQ then joins RNA 3′ ends to the linker's pre-adenylated RNA end; SuperScript III is added directly, and no separate primer is listed (line 134): the 2019 text says the complementary linker strand is extended into the ligated RNA. Nuclei are washed to remove free linker before T4 DNA ligase joins the duplex end to proximal AluI-cut genomic DNA (16 °C overnight).

🟢 Proteinase K/SDS at 65 °C reverses crosslinks; extracted DNA is bound to streptavidin beads, and 150 mM NaOH releases the non-biotinylated strand (line 136; 2019 line 118 calls it the cDNA-containing strand). Random hexamers and Klenow (3′→5′ exo−) make it double-stranded. MmeI with SAM cuts about 20 nt beyond each of two outward-facing linker sites; a 12% native gel selects the 85-bp RNA-tag–linker–DNA-tag band, while a 65-bp band diagnoses a linker joined to only one partner. The 2019 protocol (line 122) states that each MmeI product carries a two-base 3′ overhang, matched by the Y adapter's random 3′-`NN` end.

🟢 Shrimp alkaline phosphatase dephosphorylates the tags; the Y adapter (`/5Phos/AGATCGGAAGAGCACACGTCT` annealed to `ACACTCTTTCCCTACACGACGCTCTTCCGATCTNN`) is ligated; T4 PNK and a second ligase addition seal the remaining nick; a 165–185-bp band is gel-selected; Phusion PCR (16 cycles) with Primer #1 `AATGATACGGCGACCACCGAGATCTACACBBBBBACACTCTTTCCCTACACGACGCTCTTCCGATCT` and Primer #2 `CAAGCAGAAGACGGCATACGAGACGTGTGCTCTTCCGATCT` gives a 194-bp band that is gel-recovered (line 138).

🟡 Primer #2 ends in the exact reverse complement of the phosphorylated Y strand; its remaining 5′ tail is `{{= grid_seq.GRID_P7_TAIL }}` on the library strand, not the canonical TruSeq `ATCTCGTATG…` P7 string. With 20-nt tags the modelled library is {{= len(grid_seq.FINAL_LIBRARY) }} bp; with the 2019 protocol's 20–21-nt tags the published 194-bp band lies within the modelled range.

**Model corrections, 2026-10-17 bounded pass.** The earlier model added two spurious 6-nt "MmeI boundary" blocks inside the insert (MmeI sites lie within the linker ends) and appended the canonical P7 string after the Y arm, giving a 210-bp library; both are removed. The graph previously placed AluI after reverse transcription, added a separate RT primer and omitted PNK, strand release, gel selections, SAP and nick sealing; the MmeI end was drawn as a 5′ overhang. These are repository corrections, not source discrepancies. A later 2026-10-17 correction removed an extra unpaired G that the model had placed between `rNNN` and the L2-paired arm (the printed `rG` is the first base of `GACACAGC`); the model previously drew a 46-bp linker and a 194-bp library with 20-nt tags. The shared bridge helper now takes the caller's declared ribonucleotide positions inside the duplex arm, so the paired `rGrA` are drawn as an RNA segment of the L1 strand (2026-10-17 follow-up); earlier they were drawn in DNA letters. Random hexamers, the Y adapter (with its printed random 3′-`NN`) and PCR Primer #1 (with its `BBBBB` barcode window) are now exact graph input states.

## Read layout

🟢 The paper reports single-end 100-nt HiSeq 2500 sequencing with `ACACTCTTTCCCTACACGACGCTCTTCCGATCT`. Linker sequence and orientation divide each read into an RNA-derived tag and a genomic-DNA tag. Because the Y adapter is ligated at both MmeI ends, either tag can lie next to Read 1; the linker orientation assigns them. PCR Primer #1 contains a five-base multiplex barcode between P5 and the Read 1 site.

🔴 Neither the 2017 paper nor the 2019 detailed protocol states an index-read cycle recipe or a separate primer for that five-base barcode. Its physical position is modeled from PCR Primer #1, but an index-primer binding site is not inferred.

## Open questions (bounded pass 2026-10-17)

- Index-read recipe for the five-base Primer #1 barcode (as above).
- Exact MmeI tag lengths (2019 protocol: 19–23 bp); the model uses 20 nt per tag.
- The model does not draw primer-binding geometry for the Y-ligated PCR template.
- Not inspected in this pass: Extended Data Fig. 1 rendered image, 2019 protocol procedure steps beyond Table 1 and the design sections.
