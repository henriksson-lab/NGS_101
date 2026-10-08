# Single-molecule molecular inversion probes (smMIPs)

## 1. What it is

smMIP combines gap-fill molecular inversion probes with a twelve-base random molecular tag. Two probe arms hybridize on either side of a genomic interval, polymerase copies the gap, ligase closes the probe into a circle, exonuclease removes all linear material, and universal PCR converts the surviving circles into an Illumina library. Reads sharing one tag can be collapsed to a single-molecule consensus.

## 2. Defining paper

🟢 Hiatt JB, Pritchard CC, Salipante SJ, O'Roak BJ, Shendure J. *Single molecule molecular inversion probes for targeted, high-accuracy detection of low-frequency variation.* Genome Research 2013. DOI: [10.1101/gr.147686.112](https://doi.org/10.1101/gr.147686.112).

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

🟡 Target-specific arms and gap lengths vary across a probe panel. The diagram therefore shows the invariant architecture rather than inventing one universal target sequence. The paper's original library is represented with its universal paired-end Illumina roles; exact modern index sequences depend on the implementation.
