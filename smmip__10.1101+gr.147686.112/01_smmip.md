# Single-molecule molecular inversion probes (smMIPs)

## 1. What it is

smMIP combines gap-fill molecular inversion probes with a twelve-base random molecular tag. Two probe arms hybridize on either side of a genomic interval, polymerase copies the gap, ligase closes the probe into a circle, exonuclease removes all linear material, and universal PCR converts the surviving circles into an Illumina library. Reads sharing one tag can be collapsed to a single-molecule consensus.

## 2. Defining paper

🟢 Hiatt JB, Pritchard CC, Salipante SJ, O'Roak BJ, Shendure J. *Single molecule molecular inversion probes for targeted, high-accuracy detection of low-frequency variation.* Genome Research 2013. DOI: [10.1101/gr.147686.112](https://doi.org/10.1101/gr.147686.112).

🟢 Cantsilieris S, Stessman HA, Shendure J, Eichler EE. *Targeted Capture and High-Throughput Sequencing Using Molecular Inversion Probes (MIPs).* Methods in Molecular Biology 2017 (PMC5484527 author manuscript). DOI: [10.1007/978-1-4939-6442-0_6](https://doi.org/10.1007/978-1-4939-6442-0_6). (Earlier versions of this note misattributed it to O'Roak et al.; the JATS front matter lists Cantsilieris, Stessman, Shendure and Eichler.) This later protocol prints complete PCR and run-primer sequences (Notes 8 and 18, JATS lines 362 and 402), but for a different design: 70-mer MIPs (line 40), 40–45-nt arms, 5–8-nt tags, 162-nt captures and a 276-bp product (Notes 15–17). 🟡 Their identity with the 2013 panel's primers is therefore inferred, not established. The original paper's publisher supplement (Table S9) is behind a Cloudflare challenge at <https://genome.cshlp.org/content/23/5/843/suppl/DC1> (checked 2026-10-17) and remains a manual source.

## 3. Published probe architecture

🟢 The original panel uses an 80-nt probe: 40 nt of targeting arms in total (16–20 nt extension arm and 20–24 nt ligation arm), a constant 28-nt backbone, and a 12-nt degenerate molecular tag. Its gap-fill target is fixed at **112 nt**. The model's 18+22 arm split is one permitted representative; none of the target-specific bases is asserted. Source: Hiatt et al. 2013, Methods, `gr.147686.112_PMC3638140.pdf.txt` lines 949–965. The tag is attached during capture, before universal PCR.

🔴 The exact 28-nt invariant backbone is assigned to Supplemental Table S9 in the original paper. That table is absent from the locally obtained PMC image bundle and remains a manual source in `ref/MANIFEST.md`; the capture circle therefore keeps those bases symbolic. The 2013 article also assigns its custom sequencing primers and 8-nt sample indexes to Table S9 (Methods “Sequencing and primary analysis”, xml.txt lines 161–163); the primers drawn from the 2017 protocol are inferred stand-ins.

The synthesized probe is initially unmodified. 🟢 The pooled probes are 5′-phosphorylated by T4 polynucleotide kinase before target capture (same source, lines 958–967). The capture-ready probe is one 5′-phosphorylated single strand:

`5′-p [ligation arm]—[12-nt molecular tag]—[constant backbone]—[extension arm]-3′`

The extension arm supplies the polymerase 3′-OH. The ligation arm supplies the phosphorylated 5′ acceptor that closes the copied product.

🟡 **Published reaction-volume discrepancy.** Hiatt et al. 2013, Methods “Capture and library construction” (`gr.147686.112_PMC3638140.xml.txt` line 147), says to add water “to 10 mL” while also stating a **20 µL** total gap-fill/ligation reaction. Those volumes cannot both apply to one reaction. The construct model does not use a reaction volume; the intended water amount has not been established.

🟢 **Reaction conditions** (Hiatt et al. 2013, Methods, xml.txt lines 141–155): pooled probes phosphorylated with 50 U T4 PNK in 1× T4 DNA ligase buffer, 37°C 45 min, 80°C 20 min. Capture: 500 ng gDNA, 330 fmol probe, 1× Ampligase buffer in 10 µL; 98°C 3 min, 85°C 30 min, 60°C 60 min, 56°C 120 min. Gap-fill/ligation: dNTPs, betaine, NAD+, 5 U Ampligase, 3.2 U Phusion; 56°C 60 min, 72°C 20 min. Exonuclease I (20 U) + Exonuclease III (100 U), 37°C 45 min, 80°C 20 min. Two 50 µL PCRs per capture with Phusion Hot-Start II, forward and indexed reverse primers at 500 nM, SYBR green, 26 cycles (98°C 15 s, 65°C 15 s, 72°C 45 s). 1.8× Ampure XP per sample; pooled products run on a 6% PAGE gel and the ~280 bp capture band excised to remove a low-molecular-weight artifact (an Ampure-based size enrichment alternative is also described). The graph now carries these named reagents and the cleanup/gel steps.

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

🟡 These five sequences are verbatim from Cantsilieris et al. 2017 (Notes 8 and 18), not from the 2013 paper; see §2.

🟢 The original assay collected two 101-nt reads plus one 8-nt index read on HiSeq 2000, or two 152-nt reads plus one 8-nt read on MiSeq (xml.txt line 161); it has no i5 read. The first 12 nt of the reverse read are the molecular tag (line 163). The page draws the HiSeq read lengths. 🟡 Target-specific arms and captured-gap lengths remain variable, so the diagram shows their roles without inventing target sequence.
