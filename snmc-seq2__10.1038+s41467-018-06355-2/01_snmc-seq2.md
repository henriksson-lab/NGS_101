# snmC-seq2

## 1. What it is

snmC-seq2 is a plate-scalable single-nucleus methylome protocol that improves library complexity and combines inline molecular barcodes with dual unique indexes.

## 2. Authoritative source

🟢 Luo et al., *Robust single-cell DNA methylome profiling with snmC-seq2*, Nature Communications (2018), DOI `10.1038/s41467-018-06355-2`.

## 3. Construct-changing operations

🟢 The workflow performs bisulfite conversion and random-primer-based adapter incorporation on single-nucleus DNA, followed by indexed amplification.

## 4. Barcode geometry

🟢 Supplementary Methods list eight random primers. Each has a 5′ spacer, a TruSeq Read 1 tail, a 6-nt inline cell barcode, and nine 3′ random bases. `P5L_AD002_H` carries `CGATGT` as the representative inline barcode drawn on the page.

🟢 The 768-primer table supplies paired P5/i5 and P7/i7 amplification primers. The proprietary Accel-NGS Adaptase-added bases are not printed, so that short region remains explicitly inferred on the schematic.
