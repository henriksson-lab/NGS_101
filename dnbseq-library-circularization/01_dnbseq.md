# DNBSEQ library circularization and DNA-nanoball preparation

## 1. What it is

DNBSEQ converts an adapter-bearing linear library into single-stranded circles and then into compact rolling-circle products containing many tandem copies of each original library molecule.

## 2. Sources and evidence

- 🟢 MGI, *MGIEasy Fast FS Library Prep Set User Manual 2.0*, including the circularization and DNB-preparation chapter.
- 🟢 The vendor's [manual 3.0](https://mgi-tech.com/uploads/files/20250818/f9e41c9729cb0df7ff07f9ddade2324b.pdf), section 4.1, likewise identifies the separate Dual Barcode Circularization Kit and names Dual Barcode Splint Buffer plus DNA Rapid Ligase for ssDNA closure. It gives reagent amounts and conditions, but no splint bases or index oligo ordering sequences; these cannot be recovered from the workflow description alone.
- 🟢 The documented route denatures the PCR product, uses a dual-barcode splint and DNA ligase for ss-circle formation, digests residual linear DNA, and performs DNB preparation.
- 🟢 DNB generation uses rolling-circle replication; the product contains repeated copies of one circular template and loads onto the patterned DNBSEQ array.
- 🟢 MGI's [DNBSEQ Dual Barcode Adapter & Barcode Sequences 3.0](https://en.mgi-tech.com/Uploads/Temp/file/20240530/66583952cae5b.pdf), p. 4, prints UDB/UMI library-prep adapter trimming arms: forward `5′-AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA-3′` and reverse `5′-AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG-3′`. Its following tables print kit-specific index choices. These are not a universal circularization input or a complete adapter-junction description.
- 🟢 The [Dual Barcode Circularization Kit Manual 7.0](https://global-mgitech.com/wp-content/uploads/Documents/SOP-013-B02-185-7.0-MGIEasy-Dual-Barcode-Circularization-Kit-User-Manual-2.pdf), pp. 1–2, names Dual Barcode Splint Buffer and DNA Rapid Ligase but does not print the splint bases.
- 🔴 The exact circle junction for a selected upstream library kit and the sequencing-primer bases remain unresolved. The adapter guide refers to a separate DNBSEQ Platform Sequencing Primer Reference Guide, which has not yet been reviewed.

## 3. Molecular path

Circularization changes the topology before sequencing. Rolling-circle replication then changes copy number without making independent PCR molecules: all repeated units in one nanoball derive from one closed library strand.
