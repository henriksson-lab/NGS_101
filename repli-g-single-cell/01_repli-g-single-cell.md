# REPLI-g Single Cell

## Defining source

🟢 QIAGEN, *REPLI-g Single Cell Handbook* (October 2012; original PDF mirrored by a reagent supplier), Introduction pp. 6–7 and *Amplification of Genomic DNA from Single Cells* pp. 12–14. The QIAGEN RG35 *Nextera NGS Library Preparation from DNA/cDNA amplified with REPLI-g Kits* August 2014 indexed original explicitly lists the named REPLI-g Single Cell Kit (cat. 150343/150345). Its full PDF remains inaccessible; the November 2018 RG35 revision instead lists the **Advanced DNA Single Cell** kit, WTA and WGA/WTA kits.

## Molecular workflow

🟢 A single cell is lysed, genomic DNA is denatured and the reaction is neutralized. Single-cell protocol pp. 12–14 (handbook text lines 559–630): 4 µl cells in PBS plus 3 µl Buffer D2 (33 µl reconstituted Buffer DLB + 3 µl 1 M DTT per 12 reactions), 10 min at 65 °C; 3 µl Stop Solution, ice; 40 µl master mix (9 µl H2O sc, 29 µl REPLI-g sc Reaction Buffer, 2 µl REPLI-g sc DNA Polymerase) to the 10 µl denatured DNA; 8 h at 30 °C; polymerase inactivation 3 min at 65 °C. The handbook does not say which kit component carries the primers, nor that Buffer D2 is alkaline; the graph shows each step without a pH claim. The 2012 handbook pp. 6–7 names exonuclease-resistant primers and phi29 DNA polymerase for multiple-displacement amplification. The current QIAGEN Single Cell product page describes random hexamers; the 2012 handbook does not specify their bases, modifications or exact formulation.

🟢 The handbook cites Dean et al. 2002, whose original [PNAS Methods](https://doi.org/10.1073/pnas.082089499) prints the degenerate protected predecessor pattern `5′-NpNpNpNpsNpsN-3′`. That predecessor pattern does not establish the protection positions or formulation of the QIAGEN kit primer. 🟡

🟡 The indexed August 2014 RG35 original advises 50 ng amplified double-stranded DNA, Nextera DNA Sample Preparation and indexed PCR, with MinElute cleanup. The full original and the cited October 2012 Illumina Nextera guide have not been locally inspected; the Nextera junction geometry is a family-level model.

🟢 The same 2012 QIAGEN handbook, single-cell protocol p. 12, explicitly reports high-molecular-weight DNA in no-template controls from random extension of primer dimers. The graph shows this parallel control artifact as discarded material; it is not an input to the sample library.

🔴 Exact REPLI-g primer bases/modifications remain unverified across the 2024 handbook and linked material. A conventional 2 × 150 run is only an example; the REPLI-g kit itself does not prescribe read length.

## Review trail

2026-10-17 bounded pass (round breadth-first-16): re-read handbook Introduction/Principle lines 224–315 and single-cell protocol lines 477–640; added Buffer D2/65 °C, Stop Solution neutralization, master-mix components, 8 h/30 °C and 65 °C/3 min inactivation to the graph; removed unsupported "alkaline" and "proprietary random hexamers" wording from the page summary. Current 2024 handbook route still redirects to a QIAGEN search page; full 2014 RG35 and October 2012 Illumina Nextera guide remain unread.
