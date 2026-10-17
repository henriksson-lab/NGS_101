# CEL-Seq2 — source note

The detailed source audit, oligo transcriptions, protocol steps, and comparison with the
original CEL-Seq method are retained in the companion
[CEL-Seq / CEL-Seq2 research note](../cel-seq__10.1016+j.celrep.2012.08.003/01_cel-seq.md).

**Defining source:** Hashimshony T *et al.* “CEL-Seq2: sensitive highly-multiplexed
single-cell RNA-Seq.” *Genome Biology* 17, 77 (2016).
doi:[10.1186/s13059-016-0938-8](https://doi.org/10.1186/s13059-016-0938-8).

The schematic uses the CEL-Seq2 RT-primer and randomhexRT sequences printed in
Supplementary Table S2 (Additional file 5, MOESM5) and the 96-member primer list in
the detailed protocol (Additional file 4, MOESM4). Additional file 6 (MOESM6) is the
separate Fluidigm C1 workflow, not an oligo table; this schematic follows the manual
well workflow. The detailed protocol calls for Illumina
RP1 and indexed RPIX, directing readers to Illumina for their sequences. The
[official Illumina TruSeq Small RNA adapter table](https://support-docs.illumina.com/SHARE/AdapterSequences/Content/TruSeq-SmallRNA.htm)
(document 1000000002694 v22) prints both oligos and the RPI index table. The model uses
RP1 and representative RPI1; the actual indexed RPIX member used in each experiment is
not specified. Their 3′ ends land on 16 and 21 known bases, respectively,
of the selected CEL-Seq2 cDNA template. The transcript insert and cell UMI/barcode are
variable; sequencing-primer placement remains unverified.

## Bounded source pass, 2026-10-17 (round breadth-first-07)

Read in the published PMC4848782 package (`$CHEM_DATA/sources/cel-seq2__10.1186+s13059-016-0938-8/`):
JATS Results (primer shortening, reagent changes, ligation-free library, C1) and
Methods (CEL-Seq2 modifications 1–4, Barcode design, Sequencing); Additional file 4
detailed protocol, docx text lines 1–49 (reagents, primer design) and 166–339 (cell
isolation through library QC and sequencing); Additional file 5 Table S2 (lines 1–21).
The oligo scraper was run on the bundle; the RT-primer family (lines 51–146 of
Additional file 4), Table S2 designs and randomhexRT were checked against those lines.

Changes to the schematic, all 🟢 from Additional file 4:

- The tube/plate workflow is now declared step by step in `tools/chem.py` instead of
  through the shared barcoded-IVT helper: SuperScript II RT (42 °C 1 h) with ERCC
  spike-in, 70 °C 10 min heat inactivation, second strand with E. coli DNA polymerase I,
  E. coli DNA ligase and RNase H (16 °C 2 h), pooling of the cells for one IVT, AMPure XP
  with PEG/NaCl bead buffer, MEGAscript T7 IVT (37 °C 13 h), ExoSAP-IT "to remove
  primers" (37 °C 15 min), Mg²⁺ fragmentation (94 °C 3 min) with EDTA stop, 1.8× RNAClean XP,
  randomhexRT library RT (25 °C 10 min, 42 °C 1 h), Phusion HF PCR with RP1/RPIX
  (11 cycles, up to 15) and two 1× AMPure XP cleanups.
- Unincorporated RT primers are carried as a separate material from the first RT to
  the ExoSAP-IT step; their digestion products are only removed physically by the
  RNAClean XP cleanup. Whether the earlier cDNA bead cleanup already removes some of
  them is not stated. 🟡
- Apparent source inconsistency on the barcode set, shown on the page: the protocol
  says "If not using all 48 primers" (line 48) but lists 96 primers (lines 51–146), and
  the Methods say 168 barcodes were selected. The schematic draws one variable 6-nt
  barcode and does not choose a set size.

Still open: the protocol allows replacing the cDNA bead cleanup by 65 °C for 20 min (not
drawn as a branch); internal aRNA fragments without the RA5/UMI/barcode end are not
drawn as an off-product; the actual RPIX indices; the Fluidigm C1 branch
(Additional file 6); and the CEL-Seq 2012 predecessor paper, which remains
inaccessible.
