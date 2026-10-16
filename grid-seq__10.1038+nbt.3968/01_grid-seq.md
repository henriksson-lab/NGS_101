# GRID-seq

## Purpose

GRID-seq records a chromatin-associated RNA and a nearby genomic DNA end in one sequenceable molecule by joining both to a directional, biotinylated bridge.

## Sources read

- 🟢 Li et al., *GRID-seq reveals the global RNA–chromatin interactome*, Nature Biotechnology (2017), [doi:10.1038/nbt.3968](https://doi.org/10.1038/nbt.3968), especially Methods and Extended Data Fig. 1.
- 🟢 The open full text (PMCID `PMC5953555`) prints the two linker strands, their modifications, ligases, AluI digestion, affinity purification, MmeI tag release, adapter sequences, PCR primers and sequencing primer.
- 🟢 Zhou et al., *GRID-seq for comprehensive analysis of global RNA–chromatin interactions*, Nature Protocols (2019), [doi:10.1038/s41596-019-0172-4](https://doi.org/10.1038/s41596-019-0172-4), PMCID `PMC7721247`: Table 1 reprints the full linker, Y-adapter and PCR oligos; Methods identify MmeI's two-base 3′ overhang and the adapter's matching `NN` end (lines 106–122), and step 151 specifies single-end 100-base sequencing through the linker.

## Construct-changing path

🟢 The linker DNA strand is `5′-/5Phos/GTTGGAGTTCGGTGTGTGGGAGTGAGCTGTGTC-3′`. The partner is a mixed RNA/DNA strand, `5′-/5Phos/rGrUrUrGrGrArUrUrCrNrNrNrGrACACAGC/iBiodT/CACTCCCACACACCGAACTCCAAC-3′`, pre-adenylated before RNA ligation. The `rNNN` bases identify PCR duplicates; the internal biotin supports capture.

🔴 **Published linker discrepancy.** The 2017 defining paper, Methods “Construction of GRID-seq library” (PMC5953555 HTML text line 1272), prints the mixed strand ending `.../iBiodT/CACTCCCACACACCGAACTCCAAC`. The same group's 2019 detailed protocol, Table 1 (PMC7721247 JATS line 1781), prints `.../iBiodT/CACTCCACACACCGAACTCCAAC`, one C shorter. Both print L2 as `GTTGGAGTTCGGTGTGTGGGAGTGAGCTGTGTC` (2017 line 1272; 2019 line 1783). Its reverse complement matches the 2017 mixed strand's DNA duplex arm exactly, so the model uses that pairing-compatible 2017 value. The shortened 2019 strand's synthesis/use is not established by these texts.

🟢 Fixed nuclei are cut with AluI. T4 RNA ligase 2 truncated KQ first joins RNA to the linker's pre-adenylated RNA end; SuperScript III copies that RNA from the linker primer. T4 DNA ligase then joins the duplex end to proximal genomic DNA.

🟢 After reverse crosslinking, streptavidin purification and random-primed second-strand synthesis, MmeI cuts about 20 nt beyond each of two oppositely oriented linker sites. The 85-bp RNA-tag–linker–DNA-tag product is selected, ligated to printed Y adapters and PCR amplified.

## Read layout

🟢 The paper reports single-end 100-nt HiSeq 2500 sequencing with `ACACTCTTTCCCTACACGACGCTCTTCCGATCT`. Linker sequence and orientation divide each read into an RNA-derived tag and a genomic-DNA tag. PCR Primer #1 contains a five-base multiplex barcode between P5 and the Read 1 site.

🔴 Neither the 2017 paper nor the 2019 detailed protocol states an index-read cycle recipe or a separate primer for that five-base barcode. Its physical position is modeled from PCR Primer #1, but an index-primer binding site is not inferred.
