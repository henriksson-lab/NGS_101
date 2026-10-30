# NSCR

Defining DOI: https://doi.org/10.1101/gr.176099.114.

## Source scope

- `gr.176099.114_PMC4381527.xml.txt`: Published original main39/45/49/55 and Methods151/155–157. RNA-leader release with separate assay/NGS branches and experiment-specific library/acquisition.
- `gr.176099.114_supp_supp_gr.176099.114_Supplemental_Information.docx.txt`: Published Supplemental Methods49–54 including six-stepNSCR and control oligos51. Selected cellular native preparation and control roles read; synthetic RNA-DNA controls are not library adapters.

## Preparation and physical readout

🟢 Published main151/SupplementalMethods54: sourceMEF lysis5ml per~25millioncells,1%SDS/40mMEDTA/500mMNaCl/10mMTrispH8/0.5mg/mlPK,37C18h→phenol/chloroform. TRIzol DNA recovery removes bulkRNA while preserving the targeted RNA-DNA chains. Heat denaturation,16mlsucrosegradient and pooled400–2000ntfraction ethanol recovery precede5′thiophosphate labeling byT4PNK and biotin-maleimide coupling (SP1501/MB9001). No current/vendor end recipe is substituted.

🟢 StreptavidinNEBS1420S1hRT captures labeled RNA-DNA chains and directly tagged brokenDNA backgrounds. F flowthrough is collected. Source BufferA/B,4Murea70C5min,20xSSC twice70C10min,TE/RNB washes occur before a200µlRNB15minRT mock wash W. RNaseI0.02U/µl200µl15minRT acts before physical R-supernatant collection. Biotin remains on upstream RNA-containing fragments; no transfer onto releasedDNA. Exact released5′end/RNAremnant/partition is unassembled. Residual bead material contributes sourceC afterTEwash/0.4NNaOH5minRT/0.4NHCl neutralization and physical supernatant collection; residual beads remain separate; no unchangedRNA or guaranteed bead dissociation is inferred.

🟢 Methods151Klenow/randompriming is only the yield assay and Methods155 explicitly omits it forNGS. An unamplified fractionR control aliquot is saved beforeWGA, with no fabricated fraction size. Source experiment1 usesWGAIISEQX→NEBNextChIPSeq→HiSeq2000PE50; experiment2 usesSEQXE→IntegenXPrepXDNAChIP→PE72 (main55/155–157). Internal WGA/kit/primer polymers remain unassembled, without inferringPhi29 from the nameWGA. The local experiment_graph(1/2) keeps these source branches distinct.

🟢 Control oligos (SI51), verbatim: forward `rArUrCrArGrArUrCrArCrUrUCAGCTCCTCGGCTTAAAATG`; reverse `rGrGrUrGrCrArUrUrUrCrUrGTCCTCGCAGTAACGCAATAG`. These validate an RNA12-DNA20 control substrate; neither is a native sequencing adapter.

Open targets: original WGAII/NEBNextChIPSeq/IntegenXPrepX guide identities, exact capture/tag/release geometry, and cited later detailedNSCR protocol applicability. The later protocol is not inspected and its read length is not adopted. Failed PMC supplement index was not a failed original: the actualDOCX is available and read.

## Review limits

Capped first source/model pass, not a missing-information audit. Mandatory full-bundle scraper was run; hits are source pointers, not facts or evidence of absence. Historical kit/whole-polymer/read geometry remains open. Graph/private preview has explicit unreviewedSourceItem and no finished marker. No accepted missing-information claim.
