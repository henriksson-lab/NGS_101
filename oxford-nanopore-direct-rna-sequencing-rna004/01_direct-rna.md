# Oxford Nanopore Direct RNA Sequencing (SQK-RNA004)

## 1. What it is

SQK-RNA004 prepares native RNA for nanopore sequencing. The sequenced molecule is never converted into amplified DNA.

## 2. Sources and evidence

- 🟢 Oxford Nanopore, *Direct RNA sequencing (SQK-RNA004)*, DRS_9195_v4_revL_29June2026: https://nanoporetech.com/document/direct-rna-sequencing-sqk-rna004 (local `.txt` twin). Sections read: overview line 286; workflow table lines 296–316; kit contents lines 440–470; RCS note line 498; materials lines 504–519; library prep lines 585–900 (RTA ligation, RT, RNAClean XP cleanup, RLA ligation, WSB washes, REB elution).
- 🟢 Oxford Nanopore, *Chemistry Technical Document* CHTD_500_v1_revAV_24July2026, "Direct RNA Sequencing Kits" paragraph (lines 448–453 of the copy cached under the direct-cDNA sources).
- Garalde DR et al. *Highly parallel direct RNA sequencing on an array of nanopores.* Nature Methods 15:201–206 (2018). [doi:10.1038/nmeth.4577](https://doi.org/10.1038/nmeth.4577). This paper describes the earlier chemistry. It was not read in this pass and is not used for RNA004 facts.

## 3. Library preparation as specified (🟢 unless marked)

1. Input: 300 ng poly(A)-tailed RNA or 1 µg total RNA in 8 µl (line 389). Total RNA, or an enriched sample (poly(A) or ribo-depleted), is accepted (line 208).
2. RTA ligation (lines 625–706): RNA + NEBNext Quick Ligation Reaction Buffer + diluted RNA CS (RCS) or water + Murine RNase Inhibitor + RT Adapter (RTA) + T4 DNA Ligase (2M U/ml; Quick T4 Ligase not recommended, line 567); 15 µl, 10 min at room temperature.
3. RNA CS (RCS) is an optional control strand: yeast *S. cerevisiae* Enolase II, YHR174W (line 498). It is added to the RTA ligation mix (lines 592–654), so when used it goes through the library as a control molecule. This corrects an earlier note claim that RCS does not enter the sample library.
4. Reverse transcription (lines 711–742): water, dNTPs and 5x Induro RT buffer, then 2 µl Induro Reverse Transcriptase; 60 °C for 30 min, then 70 °C for 10 min. The protocol lists RT as a required step and the only pause point (line 305). The technical document says that bypassing RT shortens the workflow but reduces output (line 453); it does not describe that as an option in this protocol.
5. Cleanup (lines 757–810): 72 µl Agencourt RNAClean XP beads, 70% ethanol washes, elution in 23 µl nuclease-free water.
6. RLA ligation (lines 820–846): RT-RNA + NEBNext Quick Ligation Reaction Buffer + RNA Ligation Adapter (RLA) + T4 DNA Ligase; 40 µl, 10 min at room temperature.
7. Cleanup (lines 855–900): 16 µl RNAClean XP beads, two 150 µl Wash Buffer (WSB) washes with resuspension (free-adapter removal), elution in 13 µl RNA Elution Buffer (REB).
8. 🟢 Only the RNA strand is sequenced. The complementary cDNA strand is not sequenced but improves output (line 286; tech doc line 452).

## 4. What the inspected sources do not give

- 🔴 RTA and RLA bases, ligating-end geometry, motor identity and its attachment are not printed in the inspected sections. The scraper found no sequences in the protocol text.
- 🟡 That RTA targets the poly(A) tail is inferred from the input requirement and from the technical document's description of the previous kit, SQK-RNA002 (line 1055: the adapter "ligates to the 3’ poly(A) tail"). The RNA004 protocol does not print an oligo(dT) sequence.

## 5. Molecular path

The RNA 3′ end is adapted first, so it is the end presented to the pore. The native RNA travels through the pore from 3′ toward 5′, while the cDNA strand serves as structural support.

## 6. Open questions

- RTA/RLA architecture and the junctions they form. The Garalde 2018 paper (earlier chemistry) and its supplement are deferred; they cannot be assumed to describe RNA004.
- How the total-RNA input is selected at the RTA step (poly(A) dependence).
