# REPLI-g Single Cell

## Defining source

🟢 QIAGEN, *REPLI-g Single Cell Handbook* (October 2012; original PDF mirrored by a reagent supplier), Introduction pp. 6–7 and *Amplification of Genomic DNA from Single Cells* pp. 12–14.

🟡 Not read from an original: a search-index record of the QIAGEN RG35 *Nextera NGS Library Preparation from DNA/cDNA amplified with REPLI-g Kits* August 2014 revision lists the named REPLI-g Single Cell Kit (cat. 150343/150345); its full PDF remains inaccessible. The November 2018 RG35 revision is recorded in the source manifest as listing the **Advanced DNA Single Cell** kit, WTA and WGA/WTA kits; that revision was not retrieved either.

## Molecular workflow

🟢 A single cell is lysed, genomic DNA is denatured and the reaction is neutralized. Single-cell protocol pp. 12–14 (handbook text lines 541–655, steps 1–13): 4 µl cells in PBS plus 3 µl Buffer D2 (33 µl reconstituted Buffer DLB + 3 µl 1 M DTT per 12 reactions), 10 min at 65 °C; 3 µl Stop Solution, ice; 40 µl master mix (9 µl H2O sc, 29 µl REPLI-g sc Reaction Buffer, 2 µl REPLI-g sc DNA Polymerase) to the 10 µl denatured DNA; 8 h at 30 °C; polymerase inactivation 3 min at 65 °C. The handbook does not say which kit component carries the primers, nor that Buffer D2 is alkaline; the graph shows each step without a pH claim. The 2012 handbook pp. 6–7 names exonuclease-resistant primers and phi29 DNA polymerase for multiple-displacement amplification. The 2012 handbook does not specify the primers' bases, length, modifications or exact formulation (a current product-page "random hexamer" description is not an on-disk original and is not used).

🟢 The handbook cites Dean et al. 2002, whose original [PNAS Methods](https://doi.org/10.1073/pnas.082089499) prints the degenerate protected predecessor pattern `5′-NpNpNpNpsNpsN-3′`. That predecessor pattern does not establish the protection positions or formulation of the QIAGEN kit primer. 🟡

**Unavailable (published, not obtainable):** Downstream library: the 2012 handbook (step 13, p. 14; text lines 648–650) says only that amplified DNA can be used for next-generation sequencing. The cited named-kit library route, QIAGEN RG35 *Nextera NGS Library Preparation* (August 2014), and the Illumina *Nextera DNA Sample Preparation Guide* 15027987 Rev B (October 2012) it cites are not available on disk; earlier search-index excerpts (50 ng input, indexed PCR, MinElute cleanup) are not an inspected original and are not used. The November 2018 RG35 identifies the REPLI-g Advanced DNA Single Cell kit and is not a revision substitute. The graph therefore ends MDA with an aliquot/stored-remainder split and a single unresolved (MAGIC) library-preparation boundary; no tagmentation adapters, PCR primers, index scheme or widths, adapter ends or read plan are modelled.

🟢 The same 2012 QIAGEN handbook, single-cell protocol p. 12, explicitly reports high-molecular-weight DNA in no-template controls from random extension of primer dimers. The graph shows this parallel control artifact as discarded material; it is not an input to the sample library.

**Unverified:** Exact REPLI-g primer bases/modifications are not printed in the 2012 handbook and remain unverified; the 2024 handbook and RG35 revisions were not retrieved, so no absence claim is made. No read plan is modelled: the REPLI-g handbook does not prescribe one and the cited library originals are unavailable.

## Review trail

2026-10-17 bounded pass (round breadth-first-16): re-read handbook Introduction/Principle lines 224–315 and single-cell protocol lines 477–640; added Buffer D2/65 °C, Stop Solution neutralization, master-mix components, 8 h/30 °C and 65 °C/3 min inactivation to the graph; removed unsupported "alkaline" and "proprietary random hexamers" wording from the page summary. Current 2024 handbook route still redirects to a QIAGEN search page; full 2014 RG35 and October 2012 Illumina Nextera guide remain unread.

2026-10-19 native-graph migration (bounded, 1 original): re-read the 2012 handbook lines 1–6, 224–315 and 477–655. Removed the queued `batch_ngs.nextera_library` endpoint (dual 8-base indexes, canonical arms, 2 × 150 run plan). The module now owns the MDA graph, an explicit aliquot/stored-remainder split (steps 11–12) and an unresolved library-preparation boundary naming the unavailable August 2014 RG35 and Nextera 15027987 Rev B originals.

2026-10-19 independent implementation review (bounded, same original; lines 1–30, 224–330, 470–660 and reference/footer lines 1178, 1320): step order, volumes, temperatures and times confirmed. Corrections: step 11 storage is conditional ("if not being used directly"), so the stored portion is no longer described as "the rest"; the no-template-control extension now uses the same reaction buffer as the sample reaction; RG35 kit-listing statements re-marked 🟡 because neither RG35 revision was read; protocol line range corrected; self-test now checks that primer-dimer and stored material never reach the library.
