# RNAscope

## 1. What it is

An in situ RNA-hybridization assay that uses adjacent double-Z target probes and branched-DNA signal amplification to visualize individual RNA molecules while retaining tissue morphology.

## 2. Sources and evidence

- 🟢 Wang F, Flanagan J, Su N, et al. *RNAscope: A Novel In Situ RNA Analysis Platform for Formalin-Fixed, Paraffin-Embedded Tissues.* Journal of Molecular Diagnostics 14:22–29 (2012). [doi:10.1016/j.jmoldx.2011.08.002](https://doi.org/10.1016/j.jmoldx.2011.08.002).
- 🟢 [Bio-Techne's current technology description](https://www.bio-techne.com/spatial/reagents/rnascope-ish-technology/get-started) confirms 18–25 target-complementary bases and a 14-base tail per Z probe, a 28-base paired-tail site, and 20 probe pairs in a standard set.
- 🟢 The defining paper states 20 amplifier sites per preamplifier and 20 label-probe sites per amplifier: 400 labels per double-Z pair and a theoretical 8,000 per 20-pair target set.
- 🟢 [US7709198B2](https://patents.google.com/patent/US7709198B2/en), the core double-Z patent assigned to Advanced Cell Diagnostics, gives exact human GAPD probe examples in Tables 1–2 and an example preamplifier as SEQ ID NO:1. In its double-Z pair, each probe contains a target-binding region, a five-thymidine spacer and one 14-nt tail; the two tails together pair to the terminal 28 bases of the preamplifier.
- 🟡 The patent is direct evidence for a working double-Z sequence implementation, but not evidence that current commercial RNAscope products retain those exact universal sequences.
- 🔴 The exact oligo pools and amplification-tree sequences in current commercial kits remain undisclosed. The schematic therefore uses patent bases only in a separately labelled historical example and does not assign them to current products.

## 3. Molecular path

Two target probes bind adjacent RNA regions. Their separate 14-base tails become one stable preamplifier docking site only when juxtaposed. The preamplifier recruits amplifiers, each amplifier recruits labeled probes, and microscopy records the resulting fluorescent or chromogenic dot. This is signal amplification rather than RNA or library amplification.

## 4. Public sequence example

The patent's double-Z GAPD probes use the tails `GAAGTTACCGTTTT` and `CTGAGTCAAAGCAT`. Their concatenation is complementary to the example preamplifier's terminal 28-nt docking site. This makes the cooperative junction reconstructable base by base even though the current commercial reagent sequences are not established.
