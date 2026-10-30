# AHARIBO

Minati et al.2021, DOI https://doi.org/10.7554/eLife.59303 . Private capped selected-source curation; no completed schematic or missing-information claim.

🟢 Published eLife article v2 XML mainMethods335–365/405–407/453–465 inspected; DOI/title and version identity checked against publisher original. Its appended peer-review correspondence is not the operative recipe. Protein AHA labeling/click capture is distinct from direct labeling of RNA. RNA-seq and protected-footprint branches have source-specific preparations; the proteomics Mag-DBCO branch is not substituted for the soluble DBCO-biotin RNA capture route. Historical library-kit internals remain unreviewed, with no source-family import.

## Chemistry and physical fractions

🟢 Common RNA-seq source preparation methionine-depletes cells, pulses250µM AHA30min and uses anisomycin-containing sBlock10min before ice/PBS wash and hypotonic lysis. AHA is a noncanonical amino-acid building block incorporated into nascent protein; it does not become an RNA base. Clarification removes nuclei/cellular debris. Soluble50µM DBCO-PEG4-biotin reacts with lysate1h4°C before30min eMagSi-cN capture, magnet removal of supernatant and two washes. The separate proteomic route uses Mag-DBCO beads and different washes; it is not substituted into the RNA library graph.

🟢 RNA-bound ribosome/nascent-peptide complexes remain physical associations, not a new covalent RNA-biotin bond. SDS/ProteinaseK37°C75min produces protein/peptide remnants that remain until acid phenol/chloroform/isoamyl-alcohol separation. Aqueous material is transferred and isopropanol/GlycoBlue precipitated, then recovered with ethanol washing and water resuspension. No perfect capture specificity or exclusive peptide targeting is predicted; the source discusses background. The no-AHA control omits metabolic AHA input but does not assert zero DBCO reactivity.

🟢 RNA-seq uses5ng in SMART-Seq Stranded634443 and HiSeq2500 acquisition, with analytical BioAnalyzer/Qubit fractions separate from production. Exact historical kit-internal adapters, RT/PCR and read geometry remain unreviewed. The source TUG1-BOAT expression/qPCR validation sequences found by the mandatory scraper are not adopted as library primers.

🟢 RPF preparation is a branch: mainMethods453 specifies cycloheximide pretreatment/sourceRiboLace lysis, RNaseI digestion45minRT/sBlock and SuperaseIn quench. Its exact relationship to the common metabolic-pulse/sBlock timing is not unified or silently filled in. Capture uses the same locally printed soluble DBCO-biotin sequence, then ProteinaseK/phenol recovery and25–35nt RNA gel selection. Historical RiboLace kitmodule2RL001_mod2 produces the source library for HiSeq2500single-end50bp. This does not establish complete2018RiboLace equivalence or disclose that inaccessible original's entire recipe. Proteomics-only endpoints are not counted as NGS branches.

## Local model and limits

`NascentComplex` keeps peptide and RNA identities independent. A caller-chosen AHA site can be consumed once by click coupling, leaving a peptide-associated biotin handle; proteolysis leaves peptide/biotin remnants rather than transferring biotin ontoRNA. Selected synthetic labels test only these declared physical relations. The default graph is RNA-seq; RPF and no-AHA conditions are branches, not separate catalogue methods. Actual modifications, complete kit polymers, sample indexes, protein-adduct structures and read geometry remain unassembled.

One published original was chemistry-read in this capped pass; its fullmainMethods provide positive preparation evidence. Available data/transparent-reporting supplements are not represented as detailed chemistry reading. The large proteomics data-table text was skipped by the scraper's size limit; no absent-sequence claim follows. Direct historical vendor originals, detailed linked methods and remaining figure/data scopes are deferred. No full source review, accepted absence or finished public schematic.

Independent selected-source implementation review corrected the RPF branch to fork through its own CYH/rapid-ice-harvest/RL001-1 lysis and clarification boundary. It bypasses common RNA-seq cell washing, hypotonic DNaseI/protease-inhibitor lysis and the duplicated clarification. This preserves the unresolved relationship between common AHA/sBlock timing and the RPF-specific schedule without imposing an extra physical preparation. A focused branch-provenance check covers the correction.
