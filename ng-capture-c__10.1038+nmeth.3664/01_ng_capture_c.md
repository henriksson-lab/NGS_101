# NG Capture-C — research notes

## Citation

Davies JOJ et al. *Multiplexed analysis of chromosome conformation at vastly improved sensitivity.* Nature Methods 13, 74–80 (2016). [doi:10.1038/nmeth.3664](https://doi.org/10.1038/nmeth.3664)

## Sources read

- 🟢 The defining article and Methods in [PMC4724891](https://pmc.ncbi.nlm.nih.gov/articles/PMC4724891/).
- 🟢 Online Methods of the accepted author manuscript, PMC supplementary file `NIHMS65861-supplement-01.pdf` (text lines 1–271: 3C preparation, controls, adaptor addition, oligo capture, multiplexing, sequencing, low-cell adaptations, data analysis). The publisher's Online Methods are on the paywalled article page and were not inspected; the publisher supplementary figure PDF (`MOESM217`) has no methods text.
- 🟢 The later modular Capture-C protocol, [PMC7613269](https://pmc.ncbi.nlm.nih.gov/articles/PMC7613269/), used only to clarify the three modules and the distinction from later NuTi Capture-C.

## 1. What it is

🟢 NG Capture-C is a many-viewpoint-versus-all assay. A 3C proximity-ligation library is converted to an ordinary indexed Illumina library, then selected twice with biotinylated DNA oligos directed to the chosen restriction fragments.

## 2. Construct-changing steps

Online Methods (author manuscript supplement-01) unless stated.

1. 🟢 Fix cells in 2% formaldehyde for 10 min, quench with glycine, lyse to nuclei (lines 13–27).
2. 🟢 Digest in three parallel reactions after SDS then Triton X-100, with three 500 U additions of restriction enzyme over 16–24 h (lines 28–43). The enzyme is DpnII: the digestion-control qPCR is named "DpnII digestion control", Fig. 1a says four-cutters "eg Dpn II", and Supplementary Fig. 1a puts GATC in the probes. The page previously offered "DpnII or NlaIII"; NlaIII is not used in this paper, so the model now names DpnII only.
3. 🟢 Heat-inactivate at 65 °C, then ligate overnight at 16 °C with high-concentration T4 DNA ligase (Thermo Scientific) (lines 46–49).
4. 🟢 Proteinase K de-crosslinking at 65 °C overnight, RNase, phenol-chloroform and ethanol precipitation give the "3C library" (lines 50–59).
5. 🟢 Sonicate 5 µg to 200 bp on a Covaris S220 (lines 80–83).
6. 🟢 NEBNext reagents (E6000/E6040/E7335/E7500; E7370 Ultra for low cell numbers) for end repair, dA-tailing, adaptor ligation and PCR indexing with Illumina TruSeq indexed adapters; AMPure XP 1:1.8 cleanups; 6–8 Herculase II cycles "when adding the Truseq indices"; libraries prepared in duplicate and pooled after indexing (lines 83–94, 204–207). 🔴 Adapter strands, end geometry and PCR primer sequences are not printed. 🟡 The wording leaves open whether the index sits on the ligated adapter or is added by the PCR primer; the graph keeps the PCR primers and adapter geometry unavailable.
7. 🟢 Capture round 1: 1.5–2 µg library (or 1–2 µg each of exactly equimolar pooled indexed libraries) with 5 µg COT DNA and 1,000 pmol each NimbleGen HE universal and index-specific blocking oligo, dried, resuspended in NimbleGen hybridization buffer/component A, denatured 95 °C, hybridized with 13 pmol biotinylated oligo pool at 47 °C for 64–72 h (lines 95–116, 161–172).
8. 🟢 Bind to M270 streptavidin beads (47 °C, 45 min); SeqCap EZ Wash Buffer I, two stringent washes at 47 °C, Wash Buffers I–III (lines 120–140). The uncaptured library is discarded in the supernatants.
9. 🟢 PCR directly from the beads with SeqCap EZ Post-Capture LM-PCR master mix and Post LM-PCR oligos (18 cycles) or the KAPA mix of SeqCap EZ accessory kit v2 (14 cycles), in two reactions; AMPure XP cleanup (lines 141–146). Fig. 1a legend: "PCR from beads using adapter P5&7 sequences". 🔴 Those oligo sequences are vendor reagents not printed here.
10. 🟢 Capture round 2 on 75% of the amplified material (≤2 µg) with 24 h hybridization, "following the protocol outlined above" (lines 147–149); for multiplexed experiments all round-1 PCRs are pooled into one second capture (lines 174–178). Fig. 1a legend says steps 8–10 (hybridization, pull-down, PCR from beads) are repeated, so the model ends with a second post-capture PCR.

🟢 The main text says the second capture raised the total PCR cycles from 20 to 34 (JATS Results, single vs double capture paragraph).

## 3. Final library and reads

🟢 Paired-end reads identify a viewpoint-containing sonication fragment and the reporter across its 3C ligation junction. The defining workflow does not use Hi-C-style biotin fill-in to select all ligation junctions.

🟢 The original Supplementary Data workbook (`nmeth3664_MOESM218_ESM.xlsx`, ‘Oligonucleotide sequences’ sheet, row 94) prints a 120-base Hba-1 viewpoint probe as 5′-`/5Biosg/GAT CCA TGG TAG CAC AGG GCA GCT AAG ATG CAA GTC TGA AGG AGG AGT CTG GCG AGC TGC TCC TGC AGT TCC CTG GAC CCA GAA GGA TGA GCT AGC AGA TTC ACT TGA GCC AAA GGA TTC`-3′. The workbook lists further viewpoint-specific probes and positions. The public page shows this one source-written example; the graph retains the capture-oligo role because each viewpoint uses a different sequence.

🟢 Sequencing (lines 179–188): libraries diluted to 4 nM; most material on MiSeq v2 300-cycle kits as 150-nt paired-end reads; one HiSeq run of 100-nt paired-end reads. Index-read length is not stated, so the page gives no acquired cycle counts.

🟢 Duplicate handling is computational (lines 209–249): Trim Galore, FLASH read reconstruction, in silico DpnII digestion (DpnII2E.pl), bowtie1, then CCanalyser2.pl removes PCR duplicates as reads with identical start and end coordinates of each sub-fragment (plus unique interacting-read coordinates for long-range cis/trans). This is not an experimental separation.

🟡 The native model retains the named TruSeq/NEBNext route as role-level adapter ends. The earlier canonical single-index shell and eight-base index default have been removed; naming TruSeq does not identify the complete historical adapter or run geometry.

🟢 Published inconsistency, shown on the page: the Results give "up to 1,000,000 fold enrichment" for double capture, while the Fig. 1a and Fig. 2a legends give up to 3,000,000-fold (both in the PMC4724891 accepted manuscript). The model encodes no enrichment value.

## Open questions

- Publisher-version Online Methods (paywalled) not compared with the accepted-manuscript supplement.
- TruSeq adapter/index identity, NEBNext kit protocol and NimbleGen SeqCap EZ (HE blockers, Post LM-PCR oligos) vendor originals not read.
- Index read length and acquired cycles for the deposited GEO GSE67959 runs not checked.
- Capture-C predecessor (Hughes et al. 2014) and the CapSequm/captureC tools not followed.

## 4. Boundary with related protocols

🟢 Later Nuclear, Low-Input, Titrated and NuTi Capture-C variants are not silently folded into this page.

## Bounded native-graph migration, 2026-10-18

Three originals were inspected: accepted-manuscript PMC4724891 JATS Results 59 and Figure 1/2 legends 320/330; its Online Methods supplement-01 (1–209, fixation through sequencing and reduced-cell variant boundary); and published Supplementary Data workbook, *Oligonucleotide sequences* sheet, rows 94–96. The mandatory bundle scraper ran and the Hba-1 example was checked directly against the workbook text. The later modular protocol was not used for this native migration.

Independent chromatin DNA regions and probe molecules occupy separate Scene strands. A local DpnII contact constructor preserves GATC; a supplied explicit dA-tailed duplex replaces the generic endpoint helper. The native ligation product precedes index PCR and remains distinct from the expected indexed endpoint. Both post-capture PCR geometries use their preceding denatured, probe-bound captured material; original probe/template fractions remain unresolved through PCR and cleanup rather than disappearing by amplification. The source-written 120-base 5′-biotin Hba-1 probe is now a graph reagent, alongside other independently synthesized pool members; a generic target cannot assign its exact binding coordinate.

The page preserves the reported enrichment discrepancy without encoding a yield. It declares source alternatives of MiSeq 150-base and HiSeq 100-base paired-end reads, with index acquisition unassigned. Native historical TruSeq/NEBNext and SeqCap EZ oligo geometry, exact index placement/widths, publisher-version Methods comparison and the separately cited predecessor/vendor originals remain deferred. This capped migration is not a complete protocol review or an accepted absence audit.
