# MethylSaferSeqS

Detection of rare mutations, copy number alterations, and methylation in the same template DNA molecules. PNAS (2023), [doi:10.1073/pnas.2220704120](https://doi.org/10.1073/pnas.2220704120).

🟢 Sources: published PMC10104560 main244/250–274/398–410; publishedSI TableS1 text and originalp4. TWO originals counted before opening. Mandatory scraper reviewed18hits/1SIfile; modifications/line wrapping were checked in the original table. This selected graph follows paired LA/LM WGS, not the48-target or hemi-nested validation branches. Main250/400 explicitly adapts SaferSeqS(10.1038/s41587-021-00900-z) with methylated adapters and strand copying/separation; unreviewed base-kit/end chemistry is not blindly imported.

🟢 The same original DNA duplex yields physically separate libraries: copied strands (LA) carry genetic sequence, while native strands (LM) retain original methylation until bisulfite conversion. There is no physical concatenation of the libraries. Single-primer copying is not exponential PCR. Recovery, carryover and conversion need not be perfect; computational UID matching does not physically join or discard DNA.

🟢 Modified adapter/primer examples, all5′→3′, fromTableS1:

```
3' adapter: /5Phos/ATAAAA/iMe-dC/GA/iMe-dC/GG/iMe-dC/NNNNNNNNNNNNNNAGAT/iMe-dC/GGAAGAG/iMe-dC/A/iMe-dC/A/iMe-dC/GT/iMe-dC/TGAA/iMe-dC/T/iMe-dC//iMe-dC/AG*T*/iMe-dC//3SpC3/
5' adapter: A/iMe-dC/A/iMe-dC/T/iMe-dC/TTT/iMe-dC//iMe-dC//iMe-dC/TA/iMe-dC/A/iMe-dC/GA/iMe-dC/G/iMe-dC/T/iMe-dC/TT/iMe-dC//iMe-dC/GAT/iMe-dC/T
copy primer: /52-Bio/iSp18/TTTTTTTTTTUGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
LA1: ACACTCTTTCCCTACACGACGCTCTTCCGATCT
LA2: /52-Bio/GACTGGAGTTCAGACGTGTGCTCTTCCGATCT
LM1: ACACTCTTTCCCTACACGACGC
LM2: /52-Bio/GACTGGAGTTCAGACGTGTGC
```

`/52-Bio/` is dual biotin, `/iSp18/` an18-atomhexaethyleneglycol spacer, U is DNAdeoxyuridine, `*` marks bonds and `/3SpC3/` blocks the3′end. The14N UID is not a sample index and its bases are not chemically declared methylated. AdapterCs outside theUID are source5mC. WGSprimers add8N well/sample indexes and two terminal phosphorothioate bonds; full raw modified fixtures and structured chemical states are in the module.

🟢 After the source-delegated SaferSeqS preparation reaches ligation with the modifiedYadapters, use one5µMdualbiotin/dUprimer for three denaturation/copying rounds with KAPAHiFiHotStartUracilPlus7959052001. Original methylation is not copied; hybridized original/copied strands and two earlier copied strands coexist in the source description. StreptavidinMyOneDynabeadsT1 captures copied strands. Heat beads in22.5µlwater75°C2min, retain original-enriched supernatant, repeat with22.5µlwater and combine45µl. Remaining beads/copies follow the separateLA path. Unbound/carryover classes remain unquantified.

🟢 LM: after denaturation add100µlD5004CTconversionreagent,50°C3h dark/4°C10min. SiliconAbindingplate/600µlMbinding/2800g5min;500µlMwash.200µlMdesulphonationRT15min,2800g5min;wash600µl then400µlMwash;elute2×21µl(~42µl). DNAU is the chemical conversion product; PCR later reportsT. LM100µlPCR uses5µMshortGC-endingprimers and KAPAHiFiHotStartUracilPlus. Short primers favor original methylated-adapter DNA after conversion; contaminating copied species are not asserted perfectly absent.

🟢 LA: USER3µl+EB37µl at37°C60min cleaves primerdU on bead-associated copies. Magnetic separation retains40µlsupernatant, removing bead/primer remnants only at this physical step. LA100µlPCR uses5µMfullprimers and KAPAHiFiHotStartK2601, distinct fromLMuracil-tolerant mix. Both branch PCR programs:98°C45s,8cycles98°C15s/60°C30s/72°C30s.180µlSPRISelect/elute100µlEB.100ng of each branch is allocated to WGS while remainder is unresolved;50µlUltraQ5/primers2µM,98°C30s then7cycles98°C10s/65°C75s. NovaSeq6000PE100 sequences the two separately indexed libraries; source index cycle counts remain unassigned.

Five narrowly documented whole-text exceptions cover one modified adapter and the same two wrapped WGS primer rows in each separate branch. Whole modified fixtures retain every fixed base; the local3adapter test normalizes the complete actual source ordering, and the original page is independently reviewable. Remaining exact reaction/primer landing/strand geometry and base/vendor originals are unresolved; no full molecular validation, final source audit or accepted missing-information claim.
