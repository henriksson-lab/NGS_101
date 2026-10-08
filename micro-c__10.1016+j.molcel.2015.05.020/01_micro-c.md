# Micro-C

Hsieh T-HS, Weiner A, Lajoie B, Dekker J, Friedman N, Rando OJ. “Mapping Nucleosome
Resolution Chromosome Folding in Yeast by Micro-C.” *Molecular Cell* 58, 1208–1218
(2015). doi:[10.1016/j.molcel.2015.05.020](https://doi.org/10.1016/j.molcel.2015.05.020).

Evidence: 🟢 stated by the defining paper · 🟡 computed molecular consequence · 🔴 unknown.

## 1. What it is

Micro-C is a chromosome-conformation library protocol in which MNase replaces the
restriction enzyme used by Hi-C. Crosslinked chromatin is digested mostly to
mononucleosomes, the heterogeneous DNA ends are repaired and biotin-labelled, and nearby
nucleosomes are proximity-ligated. Paired reads from the two genomic ends identify a
contact at nucleosome-scale resolution. 🟢

This page models the original budding-yeast Micro-C protocol. Micro-C XL adds an extended
protein–protein crosslinker and is a separate 2016 variant.

## 2. Sources read

- 🟢 The open primary article, PMCID **PMC4509605**, including Methods sections
  “Chromatin digestion and end repair,” “Proximity ligation,” and “Library preparation
  and sequencing.”
- 🟢 The paper’s Figure 1 overview and description of the two-stage ligation-product
  purification.
- 🟡 Canonical TruSeq arms from `lib/illumina.py`; the paper names Illumina adaptors and
  paired-end primers but does not print their sequences.

## 3. Workflow

1. 🟢 Fix budding yeast with **3% formaldehyde for 15 min**, quench with 125 mM glycine,
   spheroplast, and digest with MNase to **>95% mononucleosomes**.
2. 🟢 Stop MNase, concentrate the chromatin and dephosphorylate with Antarctic phosphatase.
3. 🟢 Treat crosslinked chromatin with **T4 DNA polymerase plus ATP** to leave 5′
   single-stranded termini. MNase is sequence-independent, so these are heterogeneous;
   Micro-C has no fixed restriction-junction motif.
4. 🟢 Add **biotin-dCTP, biotin-dATP, dTTP and dGTP** to make biotinylated blunt dsDNA.
   🟡 Which positions acquire biotin is determined by each genomic overhang; the shared
   end-repair function derives those positions rather than assigning a fixed motif.
5. 🟢 Dilute 0.5–1 µg crosslinked chromatin to 10 mL and proximity-ligate with T4 DNA
   ligase. Contact products contain an internal, biotin-labelled blunt ligation boundary.
6. 🟢 Heat-inactivate ligase, reconcentrate, then treat with **100 U exonuclease III for
   5 min**. This eliminates labelled ends of unligated DNA while internal junctions are
   protected.
7. 🟢 Reverse crosslinks overnight with proteinase K, purify DNA, treat with RNase A and
   gel-select **250–350 bp** products.
8. 🟢 End-repair with End-It, dA-tail with exonuclease-minus Klenow and ligate Illumina
   adapters.
9. 🟢 Capture adapter-ligated DNA on streptavidin beads. This second selection separates
   biotin-bearing ligation products from undigested dinucleosomal DNA.
10. 🟢 PCR-amplify on beads for approximately **12–15 cycles** using Illumina paired-end
    primers, then sequence paired-end on HiSeq.

## 4. Molecular distinctions from restriction Hi-C

- Micro-C fragmentation is by MNase around protected nucleosomes, not by a sequence-
  specific restriction enzyme.
- T4 DNA polymerase first resects ends to 5′ overhangs and then fills them to blunt ends.
- Both dATP and dCTP are biotinylated; there is no invariant `GATCGATC`-like junction.
- Specificity comes from two selections: exonuclease III against unligated terminal
  labels, then streptavidin capture of internal labelled junctions after adapter ligation.

## 5. Final library and sequencing

🟡 The paper specifies Illumina indexed adapters without their bases. The page therefore
marks a canonical single-index TruSeq structure as inferred:

```
P5 · Read 1 · nucleosome A DNA · variable repaired/ligated junction · nucleosome B DNA · Read 2' · i7' · P7'
```

Read 1 and Read 2 begin at opposite genomic ends. They are mapped separately to recover
the contacting nucleosomes; the variable proximity junction generally lies inside the
insert rather than defining either read start. 🟢
