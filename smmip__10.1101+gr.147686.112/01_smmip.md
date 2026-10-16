# Single-molecule molecular inversion probes (smMIPs)

## 1. What it is

smMIP combines gap-fill molecular inversion probes with a twelve-base random molecular tag. Two probe arms hybridize on either side of a genomic interval, polymerase copies the gap, ligase closes the probe into a circle, exonuclease removes all linear material, and universal PCR converts the surviving circles into an Illumina library. Reads sharing one tag can be collapsed to a single-molecule consensus.

## 2. Defining paper

🟢 Hiatt JB, Pritchard CC, Salipante SJ, O'Roak BJ, Shendure J. *Single molecule molecular inversion probes for targeted, high-accuracy detection of low-frequency variation.* Genome Research 2013. DOI: [10.1101/gr.147686.112](https://doi.org/10.1101/gr.147686.112).

🟢 O'Roak et al. *Targeted Capture and High-Throughput Sequencing Using Molecular Inversion Probes (MIPs).* Methods in Molecular Biology 2017. DOI: [10.1007/978-1-4939-6442-0_6](https://doi.org/10.1007/978-1-4939-6442-0_6). This author protocol supplies the complete PCR and run-primer sequences. The original paper's publisher supplement could not be fetched automatically and remains a manual source in the manifest.

## 3. Published probe architecture

🟢 The original panel uses an 80-nt probe: 40 nt of targeting arms in total (16–20 nt extension arm and 20–24 nt ligation arm), a constant 28-nt backbone, and a 12-nt degenerate molecular tag. Its gap-fill target is fixed at **112 nt**. The model's 18+22 arm split is one permitted representative; none of the target-specific bases is asserted. Source: Hiatt et al. 2013, Methods, `gr.147686.112_PMC3638140.pdf.txt` lines 949–965. The tag is attached during capture, before universal PCR.

🔴 The exact 28-nt invariant backbone is assigned to Supplemental Table S9 in the original paper. That table is absent from the locally obtained PMC image bundle and remains a manual source in `ref/MANIFEST.md`; the capture circle therefore keeps those bases symbolic. The 2017 author protocol separately supplies the universal PCR and sequencing primer sequences.

The synthesized probe is initially unmodified. 🟢 The pooled probes are 5′-phosphorylated by T4 polynucleotide kinase before target capture (same source, lines 958–967). The capture-ready probe is one 5′-phosphorylated single strand:

`5′-p [ligation arm]—[12-nt molecular tag]—[constant backbone]—[extension arm]-3′`

The extension arm supplies the polymerase 3′-OH. The ligation arm supplies the phosphorylated 5′ acceptor that closes the copied product.

🟡 **Published reaction-volume discrepancy.** Hiatt et al. 2013, Methods “Capture and library construction” (`gr.147686.112_PMC3638140.xml.txt` line 147), says to add water “to 10 mL” while also stating a **20 µL** total gap-fill/ligation reaction. Those volumes cannot both apply to one reaction. The construct model does not use a reaction volume; the intended water amount has not been established.

## 4. Selection logic

1. Both arms must bind the same genomic molecule in the correct orientation.
2. Polymerase copies the intervening target bases.
3. Ligase closes the nick, creating a covalently closed ssDNA circle containing the tag and the 112-nt copied target. The graph derives this circle from the gap-filled strand; its strand-order check is conditional because the targeting arms, backbone and captured bases are not supplied as one exact source-written probe.
4. Exonuclease destroys unreacted probes and genomic DNA while the circle survives.
5. Universal PCR adds the sequencing library arms.

## 5. Read layout and scope

🟢 Universal forward PCR primer: `AATGATACGGCGACCACCGAGATCTACACATACGAGATCCGTAATCGGGAAGCTGAAG`.

🟢 Universal reverse PCR primer: `CAAGCAGAAGACGGCATACGAGATNNNNNNNNACACGCACGATCCGACGGTAGTGT`.

🟢 Custom forward, reverse and index sequencing primers are respectively `CATACGAGATCCGTAATCGGGAAGCTGAAG`, `ACACGCACGATCCGACGGTAGTGT` and `ACACTACCGTCGGATCGTGCGTGT`.

The original assay collects paired reads and one eight-base i7 index read; it has no i5 read. Read 2 encounters the twelve-base molecular tag first. 🟡 Target-specific arms and captured-gap lengths remain variable, so the diagram shows their roles without inventing target sequence.
