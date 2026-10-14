# Nano3P-seq — research notes

Defining source: Begik *et al.*, [10.1038/s41592-022-01714-w](https://doi.org/10.1038/s41592-022-01714-w).

## Evidence used

- 🟢 Nano3P-seq is Nanopore 3′-end-capture cDNA sequencing and accepts total, ribodepleted or poly(A)-selected RNA.
- 🟢 Template switching captures RNA-derived molecules from the native 3′ end without RNA-adapter ligation and without PCR amplification.
- 🟢 Long reads support simultaneous abundance, tail-length and tail-composition analysis for polyadenylated and non-polyadenylated RNAs.
- 🟢 The Methods print the complete blocked `R_RNA`, 5′-phosphorylated
  `D_DNA`, and `CompA_DNA` oligos. TGIRT template-switches from the preannealed
  RNA/DNA starter duplex directly onto the native RNA 3′ end.
- 🟢 After RNA digestion, `CompA_DNA` is annealed to the cDNA handle before
  ligation of Oxford Nanopore Adapter Mix (AMX) from SQK-DCS109.
- 🔴 The paper names the commercial AMX motor adapter but does not disclose its
  full sequence. This is the remaining sequence limitation; it does not apply to the
  three custom Nano3P-seq oligos.

## Protocol boundary

This is distinct from Oxford Nanopore Direct RNA Sequencing: Nano3P-seq sequences cDNA produced from the RNA 3′ end.
