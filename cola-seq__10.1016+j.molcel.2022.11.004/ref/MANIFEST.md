# Original sources and version limits

Run `OMP_THREAD_LIMIT=1 python3 tools/get_sources.py --slugs cola-seq__10.1016+j.molcel.2022.11.004`. Sources and extracted text live in the external source store; none belongs in the checkout.

- [Accepted author manuscript PMC10448999/NIHMS1854138](https://pmc.ncbi.nlm.nih.gov/articles/PMC10448999/), acquired as `j.molcel.2022.11.004_PMC10448999.xml` via the PMC Cloud `.1` object. The manifest and original header identify an author manuscript with text-mining access, not a final publisher PDF. Inspected introduction72/workflow90–94/STAR Methods216–230 and Key Resources1436–1438.
- [Definitive article DOI](https://doi.org/10.1016/j.molcel.2022.11.004): actual attempted route returned HTTP403. Final/accepted-version equivalence remains unverified.
- Actual Tables S4 and S5 are cited for ncRNA depletion probes and library oligos. The manuscript supplement list has unavailable numbered PDF/CSV/XLSX assets; its Table-to-filename mapping was not verified. Those are genuine cited originals, not invented reagent-disclosure documents. Obtain the actual source index before assigning a numbered asset to a table.
- The automatically retrieved2021preprint and its oligo workbook are superseded and excluded. Their scraper hits do not supply final or accepted-manuscript sequences.

The SVG next to the protocol note is our generated **invented-sequence topology illustration**, not a downloaded source figure. It does not add a synthetic sequence to the real preparation. See the source-review checklist for closed inspection scope and unresolved originals.
