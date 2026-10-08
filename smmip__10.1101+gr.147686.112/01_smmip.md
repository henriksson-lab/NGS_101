# Single-molecule molecular inversion probes (smMIPs)

## 1. What it is

smMIP combines gap-fill molecular inversion probes with a twelve-base random molecular tag. Two probe arms hybridize on either side of a genomic interval, polymerase copies the gap, ligase closes the probe into a circle, exonuclease removes all linear material, and universal PCR converts the surviving circles into an Illumina library. Reads sharing one tag can be collapsed to a single-molecule consensus.

## 2. Defining paper

🟢 Hiatt JB, Pritchard CC, Salipante SJ, O'Roak BJ, Shendure J. *Single molecule molecular inversion probes for targeted, high-accuracy detection of low-frequency variation.* Genome Research 2013. DOI: [10.1101/gr.147686.112](https://doi.org/10.1101/gr.147686.112).

🟢 O'Roak et al. *Targeted Capture and High-Throughput Sequencing Using Molecular Inversion Probes (MIPs).* Methods in Molecular Biology 2017. DOI: [10.1007/978-1-4939-6442-0_6](https://doi.org/10.1007/978-1-4939-6442-0_6). This author protocol supplies the complete PCR and run-primer sequences. The original paper's publisher supplement could not be fetched automatically and remains a manual source in the manifest.

## 3. Published probe architecture

🟢 The paper describes two 16–24 nt targeting arms joined by a constant 28 nt backbone and a 12 nt degenerate molecular tag. The tag is attached during capture, before universal PCR.

The probe is one 5′-phosphorylated single strand:

`5′-p [ligation arm]—[12-nt molecular tag]—[constant backbone]—[extension arm]-3′`

The extension arm supplies the polymerase 3′-OH. The ligation arm supplies the phosphorylated 5′ acceptor that closes the copied product.

## 4. Selection logic

1. Both arms must bind the same genomic molecule in the correct orientation.
2. Polymerase copies the intervening target bases.
3. Ligase closes the nick, creating a covalently closed ssDNA circle containing the tag and copied target.
4. Exonuclease destroys unreacted probes and genomic DNA while the circle survives.
5. Universal PCR adds the sequencing library arms.

## 5. Read layout and scope

🟢 Universal forward PCR primer: `AATGATACGGCGACCACCGAGATCTACACATACGAGATCCGTAATCGGGAAGCTGAAG`.

🟢 Universal reverse PCR primer: `CAAGCAGAAGACGGCATACGAGATNNNNNNNNACACGCACGATCCGACGGTAGTGT`.

🟢 Custom forward, reverse and index sequencing primers are respectively `CATACGAGATCCGTAATCGGGAAGCTGAAG`, `ACACGCACGATCCGACGGTAGTGT` and `ACACTACCGTCGGATCGTGCGTGT`.

The original assay collects paired reads and one eight-base i7 index read; it has no i5 read. Read 2 encounters the twelve-base molecular tag first. 🟡 Target-specific arms and captured-gap lengths remain variable, so the diagram shows their roles without inventing target sequence.
