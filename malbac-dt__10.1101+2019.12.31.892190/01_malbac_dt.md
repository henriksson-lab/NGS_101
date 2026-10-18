# MALBAC-DT — research notes

Authority (published version first): Chapman, Lee, Cai, Ma, Li, Sun, Xie, "Correlated gene
modules uncovered by high-precision single-cell transcriptomics", *PNAS* 119, e2206938119
(2022), [10.1073/pnas.2206938119](https://doi.org/10.1073/pnas.2206938119),
[PMC9907105](https://pmc.ncbi.nlm.nih.gov/articles/PMC9907105/). The directory is named for
the superseded preprint (Chapman, Lee *et al.*, bioRxiv
[10.1101/2019.12.31.892190](https://doi.org/10.1101/2019.12.31.892190), v1 posted 2 Jan 2020);
an earlier version of these notes misattributed it to "Chen *et al.*".

## Sources read (2026-10-17 bounded pass)

- PNAS JATS (PMC Cloud) Materials and Methods "MALBAC-DT Protocol" (`.xml.txt` lines 195–247),
  "Automated MALBAC-DT Protocol" (lines 249–251), "Sequence Processing" (lines 257–259); Results
  "Development of MALBAC-DT" (lines 115–121); Fig. 1D legend (line 109).
- PNAS SI Appendix (`supp_pnas.2206938119.sapp.pdf.txt`): Figs S1–S10 and Tables S1–S3 only
  (performance and module tables, lines 165–276); no methods or oligo table.
- Preprint v1 Methods (lines 697–871) only to confirm supersession. The preprint says Ix-Tru2
  replaces the Nextera **S5XX** primer; the published text says **S7XX**. This is an ordinary
  revision change and the published value is used.
- `tools/scrape_primers.py` over the whole bundle: 0 sequences.

## Evidence

- 🟢 Cells are sorted into 3 µL lysis buffer (SSIV buffer, ICA-630, betaine, SUPERase In,
  0.2 µL 50 µM RT-An primer, dNTP); 72 °C 3 min, then SuperScript IV RT at 55 °C 10 min.
- 🟢 Fig. 1D: the RT primer is poly-T with a cell-specific barcode and UMI. Sequence
  Processing: UMIs follow the (HBDV)<sub>5</sub> RT-An pattern (20 nt); consensus UMIs must
  differ from the (VDBH)<sub>5</sub> RT-Bn pattern at ≥3 bases.
- 🟢 Excess RT primers are degraded with Exo I (37 °C 30 min, 80 °C 20 min); the Exo mix also
  contains 0.2 µL 50 µM RT-Bn. 🔴 RT-Bn's molecular role is not stated.
- 🟢 Amplification: GAT-7N + GAT-COM, Deep Vent (exo-) in ThermoPol; 10 MALBAC looping cycles
  (4→65 °C ramp, 95 °C, 58 °C), then PCR steps 12–16 written "Goto 13 17×"; then Tru2-Gn-RT and
  five more cycles of steps 12–15.
- 🟢 Original plate: 8 row RT3-An primers and 12 column Tru2-Gn-RT primers (inner/outer
  barcodes demultiplexed per cell). Automated plate (K562): 96 well-specific RT3-An and one
  Tru2-G-RT; 15 amplification cycles + 5.
- 🟢 1 µL per well is pooled and purified with 0.8× AMPure; the Nextera kit tagments, and Ix-Tru2
  primers replace the Nextera S7XX primer to select the barcoded 3′ ends. The automated
  plate used Vazyme TruePrep and NovaSeq 2 × 150 bp with a custom Read 2 primer.
- 🔴 Exact native RT-An, RT-Bn, GAT-7N, GAT-COM, Tru2-Gn-RT and Ix-Tru2 strands remain unassigned in this model; the linked automated-protocol identity is unresolved. No complete absence claim is accepted.

## Model (derived)

- 🟡 Source-named native Nextera/S5XX plus Ix-Tru2 regions define the original route; a full canonicalP5/i5/S5/ME shell is not assigned. AutomatedTruePrep/customRead2 remains distinct. Transcript and identifier read roles are interpreted with physical primer sites unassigned.
- 🟡 Displayed internal identifier/handle region order remains inferred. No physical P7 end, fixed index width or exact native primer sequence is assigned from primer names alone.

## Open questions

- Oligo sequences and barcode lengths; the cited automated protocol (reference 33 in the
  Methods resolves to an unrelated p53 paper in the reference list).
- Whether "Goto 13 17×" means 17 or 18 PCR cycles; read lengths of the original-plate run.

## Protocol boundary

MALBAC-DT applies MALBAC-style amplification to digital transcriptomics and is separate from
the MALBAC single-cell DNA protocol.

## Native graph migration (2026-10-18, batch04)

Inspected two published originals: PMC9907105 JATS originalMALBAC-DT Methods195–247,
automated249–251, SequenceProcessing257–259 and earlierFig1D/Results109–121; original
SI Appendix contents/FigsS1–S10/TablesS1–S3 headings/performance tables165–276.
Mandatory bundle scraper returned zerohits, which does not establish publication-wide
absence. Superseded preprint was not re-read or used. Linked automated-protocol
reference33 identity remains unresolved; the original cited protocol stays deferred.

The local graph separates originalU2OS/HEK293T Nextera/S5XX+Ix-Tru2 from automatedK562
TruePrep/customRead2. Originalenzyme/primer identities are not copied into unnamed
automatedRT/Exomix/PCRmastermix formulations. Original10loopingcycles precede a distinct
PCR phase (printedGoto13 seventeen times); actual looping-cycle/sscDNA populations
enter PCR, followed by source5cycle Tru2 completion. Automated15+5cycles remain separate.
RNA templates, cDNA and primer/hydrolysis/nontarget populations remain distinct and
carry through Tn5; finalPCR exposes unresolved carried-species fates. Nonamplification
and computational UMI/barcode filters do not physically discard molecules.

Source-known inner/outer barcode and consensusUMI feature identities persist through
intermediate populations. Their displayed region order is inferred, not established
from primer names. The20nt length derives from the reported consensus(HBDV)5pattern;
that does not establish the physicalRTprimer sequence/orientation. Native kit ends,
index widths/P7 shell and physicalRead1/Read2 sites remain unassigned. Automated2×150
is positive; originalreadmask and both index acquisitions are unassigned. Original
0.8×pooling is not silently assigned to the automatedsummary; its recovery is bounded.
No prior broad missing-oligo statement is promoted to confirmedabsence.
