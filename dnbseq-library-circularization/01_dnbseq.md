# DNBSEQ library circularization and DNA-nanoball preparation

## 1. What it is

DNBSEQ converts an adapter-bearing linear library into single-stranded circles and then into compact rolling-circle products containing many tandem copies of each original library molecule.

## 2. Sources and evidence

- 🟢 MGI, *MGIEasy Fast FS Library Prep Set User Manual 2.0*, including the circularization and DNB-preparation chapter.
- 🟢 The vendor's [manual 3.0](https://mgi-tech.com/uploads/files/20250818/f9e41c9729cb0df7ff07f9ddade2324b.pdf), section 4.1, likewise identifies the separate Dual Barcode Circularization Kit and names Dual Barcode Splint Buffer plus DNA Rapid Ligase for ssDNA closure. It gives reagent amounts and conditions, but no splint bases or index oligo ordering sequences; these cannot be recovered from the workflow description alone.
- 🟢 Manual 7.0 (part SOP-013-B02-185, kit version V1.0, Sep. 2024) gives the whole kit procedure numerically: input 1 pmol of dual-barcode PCR product (minimum 0.5 pmol) in 48 µL TE (section 2.1, Table 7, Formula 1); denaturation at 95 °C for 3 min then 2 min on ice (section 3.1.2, Table 9); 11.6 µL Dual Barcode Splint Buffer plus 0.5 µL DNA Rapid Ligase at 37 °C for 30 min (section 3.1.3, Tables 10–11); 1.4 µL Digestion Buffer plus 2.6 µL Digestion Enzyme at 37 °C for 30 min, stopped with 7.5 µL Digestion Stop Buffer (section 3.2, Tables 12–14); 170 µL DNA Clean Beads, two 80% ethanol washes, elution in 22 µL TE with 20 µL recovered (section 3.3); and Qubit ssDNA quantification with a yield of at least 80 fmol (section 3.4, Formula 3, Table 16).
- 🟢 Sections 1.1–1.2 restrict the input to MGI Duplex UMI / UDB adapter PCR products or an "MGI Dual-barcode standard product" and state the kit is **not** suitable for single-barcode libraries; section 1.3 says libraries made with it run on any MGI sequencer.
- 🔴 Manual 7.0 contains no DNB, rolling-circle, polymerase or sequencing-primer text at all: its endpoint is the purified, QC'd ssCir library. The nanoball steps on the page therefore describe the DNBSEQ platform rather than this kit, and are marked as such. Likewise the manual names only "Digestion Enzyme" — its identity and activity (and hence the explicit claim that it is an exonuclease acting on non-circular DNA) are not disclosed.
- 🟢 DNB generation uses rolling-circle replication; the product contains repeated copies of one circular template and loads onto the patterned DNBSEQ array.
- 🟢 MGI's [DNBSEQ Dual Barcode Adapter & Barcode Sequences 3.0](https://en.mgi-tech.com/Uploads/Temp/file/20240530/66583952cae5b.pdf), p. 4, prints UDB/UMI library-prep adapter trimming arms: forward `5′-AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA-3′` and reverse `5′-AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG-3′`. Its following tables print kit-specific index choices. These are not a universal circularization input or a complete adapter-junction description.
- 🟢 The [Dual Barcode Circularization Kit Manual 7.0](https://global-mgitech.com/wp-content/uploads/Documents/SOP-013-B02-185-7.0-MGIEasy-Dual-Barcode-Circularization-Kit-User-Manual-2.pdf), pp. 1–2, names Dual Barcode Splint Buffer and DNA Rapid Ligase but does not print the splint bases.
- 🔴 The manual does not say which of the two denatured strands is circularized, or what becomes of the complementary strand; the model keeps that fate unresolved.
- 🔴 The exact circle junction for a selected upstream library kit and the sequencing-primer bases remain unresolved. The adapter guide refers to a separate DNBSEQ Platform Sequencing Primer Reference Guide, which has not yet been reviewed.

## 3. Molecular path

Circularization changes the topology before sequencing. Rolling-circle replication then changes copy number without making independent PCR molecules: all repeated units in one nanoball derive from one closed library strand.
