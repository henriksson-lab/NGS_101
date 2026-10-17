# NEBNext Single Cell / Low Input RNA

## 1. What it is

A commercial full-length cDNA workflow for individual cells or very low RNA inputs, followed by Ultra II FS Illumina library conversion.

## 2. Sources and evidence

- 🟢 NEB E6420 manual version 6.0_01/24 (`manuale6420.pdf.txt`; Section 1 cells lines ~77–1016, Oligo Sequences/Adaptor Trimming/Read Length lines 1840–1870) and [product page](https://www.neb.com/products/e6420-nebnext-single-cell-low-input-rna-library-prep-kit-for-illumina).
- 🟢 The E6420 manual's "Oligo Sequences" table (p. 23, Rev 6.0) prints all three cDNA-stage reagents, including the three terminal riboguanosines on the TSO:
  - RT primer mix: `5′-AAGCAGTGGTATCAACGCAGAGTAC(T)30V-3′` (`V` is A/C/G). The table writes out all 30 T bases.
  - Template switching oligo: `5′-GCTAATCATTGCAAGCAGTGGTATCAACGCAGAGTACATrGrGrG-3′`.
  - cDNA PCR primer: `5′-AAGCAGTGGTATCAACGCAGAGT-3′`.
- 🟢 Amplified cDNA is processed by Ultra II FS fragmentation/end preparation, adaptor ligation and indexed PCR.
- 🟢 The same manual prints the Read 1 and Read 2 adaptor-trimming sequences, `AGATCGGAAGAGCACACGTCTGAACTCCAGTCA` and `AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT`.
- 🔴 The exact kit adaptor body, sample-index sequences and corresponding library PCR primers are not in that table. The trimming sequences are partial library motifs, not the complete adapter oligos.

- 🟢 Required Materials Not Included and step 1.11 make the adaptor, USER Enzyme and index primers part of a separately purchased NEBNext oligo kit; single-index kits (one index primer plus the universal primer) and dual-index plates are both allowed. The drawn dual 8-base index is therefore one option, not the generic E6420 layout.

## 3. Molecular path (Section 1, cells; Section 2 total RNA is parallel)

Template switching supplies the second amplification handle without ligating an adaptor to the RNA. The resulting full-length cDNA is amplified first; only then is it fragmented for an ordinary paired-end sequencing library.

1. 🟢 Lysis in 1X Cell Lysis Buffer with Murine RNase Inhibitor (1.2); RT primer annealing 5 min at 70 °C (1.3).
2. 🟢 RT and template switching with Single Cell RT Buffer, TSO and RT Enzyme Mix, 90 min at 42 °C then 10 min at 70 °C (1.4).
3. 🟢 cDNA PCR with Single Cell cDNA PCR Master Mix and the single cDNA PCR primer, 11–21 cycles, 62 °C annealing, 3 min extension (1.5). 🟡 The 23-nt primer is the 5′ end of both the TSO DNA body and the RT-primer handle, so one species primes inward from both ends; `PCRGeometry.symmetric` locates both sites on the model template. The 12-nt TSO 5′ extension `GCTAATCATTGC` lies outside the primer and is absent from PCR products. The GGG after `ACAT` (copied from rGrGrG) and the base pairing the primer's terminal V are derived.
4. 🟢 Double bead cleanup: 0.6X beads, then re-binding with Bead Reconstitution Buffer (1.6).
5. 🟢 Ultra II FS fragmentation/end prep, 25 min at 37 °C then 30 min at 65 °C (1.8).
6. 🟢 Ligation of the 1:25-diluted NEBNext Adaptor for Illumina with Ultra II Ligation Master Mix and Ligation Enhancer, 15 min at 20 °C; then 3 µl USER Enzyme, 15 min at 37 °C, "only required for use with NEBNext Adaptors" (1.9). The manual does not state the adaptor structure or what USER cleaves, so the page does not draw that mechanism.
7. 🟢 0.8X bead cleanup (1.10); Ultra II Q5 index PCR, typically 8 cycles (3–12 by input), 65 °C annealing/extension (1.11); 0.9X bead cleanup (1.12). An ~80 bp primer or 128–140 bp adaptor-dimer peak prompts a repeat cleanup (1.13.3).

Model corrections (2026-10-17 bounded pass): the amplified-cDNA state previously wrote the RT-primer handle in primer orientation on the sense strand, which left no site for the reverse-direction primer; it now carries the handle's reverse complement after the poly(A) copy, and the cDNA PCR is an exact symmetric single-primer geometry instead of an unavailable one. The TSO 5′ extension is kept only on the template-switched template. USER treatment and the three bead cleanups were added, and the final library is marked inferred because its adaptor and primer sequences come from an unprinted oligo kit.

## 4. Reads

🟢 Typical read length 2×75 (Read Length and Depth); about 500,000 reads per single cell. 🟢 Libraries "resemble TruSeq libraries" and are trimmed with the two motifs above.

## 5. Open questions

- 🔴 Adaptor, USER-site and index-primer sequences of the chosen oligo kit (the specific kit and well are user choices; the oligo-kit manuals were not read in this pass).
- RT and cDNA PCR enzyme identities are not named in the manual.
