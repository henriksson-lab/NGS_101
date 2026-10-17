# Telo-seq

## 1. What it is

Telo-seq selectively places an Oxford Nanopore motor adapter at the native human telomere end, so reads begin on the telomeric C-rich strand and proceed toward the subtelomere.

## 2. Sources read

- 🟢 Schmidt et al. (2024), *Nat. Commun.* 15:5149, [doi:10.1038/s41467-024-48917-7](https://doi.org/10.1038/s41467-024-48917-7): published PMC11189484 JATS (`s41467-024-48917-7_PMC11189484.xml.txt`) Results line 61 and Methods "Telo-seq library preparation" line 169; Supplementary Information (`..._MOESM1_ESM.pdf.txt`) Supplementary Table 7, lines 1780–1796. The bioRxiv preprint 10.1101/2023.11.28.569082 is superseded and not used.
- 🟢 Oxford Nanopore, [TELO-seq requirements](https://nanoporetech.com/document/requirements/TELO-seq), current commercial implementation (additional, not defining).
- 🟡 Makarova et al. (2025), [doi:10.3389/fmolb.2025.1725112](https://doi.org/10.3389/fmolb.2025.1725112), comparative review.

## 3. Oligos and reaction order

🟢 Six 5′-phosphorylated telorettes share `AGCAATACGTAACTGAACGAAGT` and carry the six phases of an 18-nt telomere-complement arm. S1 is `ACTTCGTTCAGTTACGTATTGCTAGCAAT`.

The published Supplementary Table 7 (`MOESM1_ESM.pdf.txt`, lines 1780–1796) labels its third and fourth distinct telorette sequences both `T3`, then continues with `T5` and `T6`; there is no `T4` row. The Methods (JATS `.xml.txt`, line 169) say six phased adapters. The model keeps the fourth printed sequence as the fourth member because it supplies a distinct phase; whether the authors intended the label `T4` is unconfirmed. 🟢 (printed labels and sequences) / 🟡 (ordinal assignment)

🟡 Derived with the model: the first {{= len(teloseq.CORE) }} nt of S1 are the reverse complement of the telorette core (model check `S1 prefix == revcomp(core)`: {{= teloseq.S1_CORE_PAIRING }}), leaving a 3′ extension `{{= teloseq.S1_EXTENSION }}`. The Methods say S1 has complementarity to the AMII adapter; because the AMII sequence is not printed, whether this extension is the AMII-pairing region is not verified.

Reaction order and conditions (🟢 Methods line 169):

1. 🟢 15 µg genomic DNA + 20 µL telorette pre-mix (six oligos in equal parts, 1 µM) with 10,000 U T4 DNA ligase (NEB M0202) in 200 µL 1x rCutSmart + 0.4 mM ATP, 35 °C 16 h, then 65 °C 10 min. Telorettes anneal to the G-overhang and are ligated to the C-strand (Results line 61).
2. 🟢 80 U EcoRV-HF (NEB R3195), 37 °C 30 min; 65 °C 20 min.
3. 🟢 45 U Klenow Fragment (3′→5′ exo−) (NEB E6053) in 250 µL 1x NEBNext dA-tailing buffer, split into five 50 µL aliquots, 37 °C 30 min. The dA-tail on the blunt ends reduces concatemer ligation (Results line 61).
4. 🟢 1x v/v AMPure XP; elution in 190 µL water, 37 °C 15 min.
5. 🟢 Anneal S1 (100 nM) in 50 mM NaCl, 50 °C 1 h.
6. 🟢 0.5x v/v AMPure XP; elution in 30 µL water, 37 °C 15 min.
7. 🟢 From the sequencing-adapter ligation step onward, follow the ONT Native Barcoding Kit (SQK-NBD111) instructions with the AMII adapter. 🔴 The kit's ligation reagents and the AMII terminal bases are not given in the article; the kit document was not inspected.

## 4. Final library and sequencing

There is no sequencing primer. The telomere-side motor controls entry of the C-rich strand, producing a telomere-to-subtelomere read. 🟢 GridION, R9.4.1 flow cells, one flow cell per library; basecalling with a customized Bonito (0.6.2) telomere model (Methods lines 169, 173).

## 5. Open questions

- The fourth telorette's intended label (printed `T3`, no `T4`).
- SQK-NBD111 kit document (AMII adapter assembly, ligation reagents and conditions): not inspected in round 2026-10-17-breadth-first-20.
- Lai et al. (reference 24; ligation conditions cited by the Methods) and the S1 3′-extension pairing with AMII.
- The ONT TELO-seq requirements page is a later commercial implementation, not the defining protocol.
