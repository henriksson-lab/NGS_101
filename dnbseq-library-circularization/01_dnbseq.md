# DNBSEQ library circularization and DNA-nanoball preparation

## 1. What it is

DNBSEQ converts an adapter-bearing linear library into single-stranded circles and then into compact rolling-circle products containing many tandem copies of each original library molecule.

## 2. Sources and evidence

- 🟢 MGI, *MGIEasy Fast FS Library Prep Set User Manual 2.0*, including the circularization and DNB-preparation chapter.
- 🟢 The documented route denatures the PCR product, uses a dual-barcode splint and DNA ligase for ss-circle formation, digests residual linear DNA, and performs DNB preparation.
- 🟢 DNB generation uses rolling-circle replication; the product contains repeated copies of one circular template and loads onto the patterned DNBSEQ array.
- 🔴 Current adapter, splint and sequencing-primer sequences are proprietary in the cited manuals.

## 3. Molecular path

Circularization changes the topology before sequencing. Rolling-circle replication then changes copy number without making independent PCR molecules: all repeated units in one nanoball derive from one closed library strand.
