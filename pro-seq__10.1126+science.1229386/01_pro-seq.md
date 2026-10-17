# PRO-seq — original 2013 source scope

PRO-seq maps engaged RNA polymerases through biotin-labelled nuclear run-on products. Sequencing begins at the reverse-complement nascent-RNA 3′ boundary, which reports the active-site position.

## Sources and original supplement identity

- 🟢 [Kwak et al., Science (2013)](https://doi.org/10.1126/science.1229386), original author manuscript PMC3974810.
- 🟢 The [author-hosted original supplementary PDF](https://www.cecerelab.com/_files/ugd/594009_926d2d49ebf94caa8e6cf707c2f4846a.pdf) prints the same paper title, Science 339:950, February 22, 2013 and defining DOI. It contains **both PRO-seq and PRO-cap**, not a separate PRO-cap-only protocol. This resolves the original Methods source through an available author copy; the PMC download route's gate remains a route-specific issue.
- 🟢 [Mahat et al., Nature Protocols (2016)](https://doi.org/10.1038/nprot.2016.086), Table 1, is a later detailed source. Its explicit RPI-n sequence is retained as a comparator, not silently assigned to every 2013 ordered oligo.

## Exact original adapter/RT option

The original supplement permits two reverse 5′ adapters with corresponding GX2/GX1 or RP1/RPI choices. The page represents its **TruSeq/RP1 option**. Sequences below use DNA T for computational pairing; the ordered adapters are RNA, printed with `rN` bases in the source.

- Reverse 3′ RNA adapter: `5′-p-{{= ''.join('r'+b for b in pro_seq.VRA3.replace('T','U')) }}-/3InvdT/-3′`.
- TruSeq reverse 5′ RNA option: `5′-{{= ''.join('r'+b for b in pro_seq.VRA5.replace('T','U')) }}-3′`.
- RP1: `5′-{{= pro_seq.RP1 }}-3′`.

🟢 Original SOM Methods print the complete modified reverse 3′ adapter, both reverse 5′ options and GX2/RP1. They name GX1/RPI PCR choices. 🟡 The displayed endpoint uses the later 2016 Table 1 RPI-n ordering family as an explicitly conditional comparison; the historical selected RPI oligo and sample assignment remain unassigned.

## Molecular order

The graph now preserves the original PRO-seq sequence:

1. Biotin-NTP run-on, RNA extraction, base hydrolysis and desalting.
2. First streptavidin enrichment and RNA extraction **before** reverse 3′ adapter ligation.
3. Reverse 3′ adapter ligation, then second streptavidin enrichment.
4. TAP and PNK processing of the RNA 5′ end.
5. Reverse 5′ adapter ligation, then third streptavidin enrichment.
6. RP1-primed RT (25 pmol), trial amplification of an RT portion to set cycles, then Phusion PCR for 12–18 cycles with excess RP1 as one primer and an RPI-index primer (12.5 pmol), and PAGE purification (SOM lines 139–148).

The SOM names T4 RNA ligase I (6 h, 20 °C) for the 3′ ligation (lines 121–124) but no enzyme for the 5′ ligation (lines 135–138). Round 2026-10-17-breadth-first-15 corrected the 5′-ligation geometry statement, which had cited the later 2016 protocol and a T4 RNA ligase for a step the 2013 SOM leaves unnamed; this was a repository attribution error, not a publication error. The SOM analysis clips adapters, applies `fastx_trimmer` to the first 26 bases and drops reads shorter than 16 bases (lines 169–171); acquired cycles remain unassigned.

RP1 anneals through its 21-base suffix within the 26-base reverse 3′ adapter. Extension copies the adapter's remaining five bases before entering the nascent RNA. Those five copied adapter bases are present in the cDNA state. Omitting them was a repository modelling error, not a source discrepancy. RNA and cDNA inserts have unknown-length role tokens, without an invented fragment length.

## PRO-cap boundary and sequencing

PRO-cap omits base-hydrolysis fragmentation, uses standard adapters, and applies Antarctic phosphatase before TAP while omitting PNK. Those steps must not replace the PRO-seq branch. The separate Hsp70 length-measurement experiment also uses specialized GX2short/biotin primers; those are not standard-library primers.

Original SOM sequencing uses GenomeAnalyzer 2 or HiSeq 2000 and describes reverse-complementing reads to recover nascent-RNA sense orientation. Exact loaded run-primer oligos and acquired cycles are not declared here. A PCR-site sequence does not establish a physically loaded sequencing primer. The six-base index belongs to the later explicit RPI comparator; no i5 read is introduced.

## Audit work remaining

The missing-original-supplement claim has been withdrawn. The remaining original SOM figures/tables, separate analysis-code supplement and chemistry-relevant cited precursor methods still need bounded source/link checks. The original supplement's available title and Methods do not by themselves complete every linked-source obligation. No missing-information claim is accepted.
