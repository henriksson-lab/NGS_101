# snmC-seq2

## 1. What it is

snmC-seq2 is a plate-scalable single-nucleus methylome protocol that improves library complexity and combines inline molecular barcodes with dual unique indexes.

## 2. Authoritative source

🟢 Luo et al., *Robust single-cell DNA methylome profiling with snmC-seq2*, Nature Communications (2018), DOI `10.1038/s41467-018-06355-2`.

## 3. Construct-changing operations

🟢 The workflow performs bisulfite conversion and random-primer-based adapter incorporation on single-nucleus DNA, followed by indexed amplification.

## 4. Barcode geometry

🟢 Supplementary Methods list eight random primers. Each has a 5′ spacer, a partial
TruSeq Read 1 tail, a 6-nt inline cell barcode, and nine 3′ H bases (each A, C or T;
G is excluded). The paper explicitly contrasts this RP-H pool with the earlier RP-N
pool, so an unrestricted N9 representation would change the chemistry.
`P5L_AD002_H` carries `CGATGT` as the representative inline barcode drawn on the page.

🟢 Supplementary Table 2 supplies 768 paired P5/i5 and P7/i7 amplification-primer
pairs. Pair 1 is shown base by base on the page: its P5 index is `ACGATCAGAC` and its
P7 index, read from the library, is `CATTGCCTAT`. The 25-nt and 27-nt 3′ annealing
regions land on the two final-library arms by computation. The proprietary Accel-NGS
Adaptase-added bases are not printed, so only that short region remains explicitly
inferred on the schematic.

Sources read: the defining Nature Communications article, its Supplementary Methods
(`s41467-018-06355-2_supp_41467_2018_6355_MOESM1_ESM.pdf`), Supplementary Table 2
(`294355_294355-3.xlsx`, identical material linked from the associated preprint), and
the associated bioRxiv manuscript (doi `10.1101/294355`).
