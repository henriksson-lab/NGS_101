# Original SMART-seq: source-era oligo follow-up

31 October 2026. Bounded follow-up of the original 2012 protocol, not SMART-seq2.
The local bundle scraper was run before reading and rerun after the two vendor
originals were retrieved. Source documents remain outside the checkout. This is
one positive-evidence pass, not a two-reader missing-information audit.

## Originals and inspected scope

| Original | Identity and locator | Inspected sections / result |
|---|---|---|
| Ramsköld et al., [doi:10.1038/nbt.2282](https://doi.org/10.1038/nbt.2282) | Existing `nbt.2282_PMC3467340.pdf.txt`, lines400–458, PDF pp7–8. The cached PDF explicitly identifies itself as an author manuscript; the Nature version-of-record route was attempted but redirected to inaccessible authentication. | Generation/amplification and sequencing-library Methods. The article material points to the newly commercialized **SMARTer Ultra Low RNA Kit for Illumina sequencing**, and names the **modified Epicentre Nextera** branch. Do not describe this pass as a fresh published-version Methods comparison. |
| Clontech, [PT5163-1 manual](https://www.takara.co.kr/file/manual/pdf/PT5163-1.pdf) | Retrieved vendor original `PT5163-1.pdf`; cover printed **March31,2011**, Cat634935, version **PR0X3693**. Cached in the original SMART-seq source bundle. | Cover, component page3 and its rendered original; text lines104–150. This is the named Ultra Low RNA kit, not older SMART PCR, SMARTer PCR, RACE or a later HV kit. |
| Epicentre, [Nextera Illumina-compatible manual, Lit307](https://ptacts.uspto.gov/ptacts/public-informations/petitions/1554159/download-documents?artifactId=FpvkPISs46Kk23Mbf9eK52zOgdd-cCdROHz5-4spQn6JMjT1X3-7xO0) | Retrieved `source_pdf_e4d8cd11a3afee8c8274.pdf`, an original nine-page vendor PDF reproduced as **PARSE EX1048**. Header lists GA09115, GA091120, GA0911-50, GA0911-96 and GABC0950; sequence notice has copyright2006–2010. No dated revision line was identified in the inspected pages. | Original rendered page2/Figure1/Sequences; text lines155–172, barcode design246–253, PCR348–408. A separate9/2011 web transcription was a discovery aid only; it does not establish revision equivalence to this PDF. |

All three live-source locators and the two new download registrations are preserved
in `catalogue/extra_sources.tsv` and the external bundle's `MANIFEST.tsv`. The
2011 Clontech cover/date and the Epicentre sequence panel were checked as rendered
originals, not accepted from scraper hits.

## Positive Clontech findings

The component page prints the following base identities, 5′→3′:

- **IS PCR**: `AAGCAGTGGTATCAACGCAGAGT`.
- **SMART CDS Primer II A**: `AAGCAGTGGTATCAACGCAGAGTACT(30)N-1N`, where
  the penultimate class is A/G/C and the terminal class is A/C/G/T. This agrees
  with the article's `T(30)VN` notation.
- **SMARTer II A**: `AAGCAGTGGTATCAACGCAGAGTACXXXXX`, explicitly identifying
  X as a proprietary undisclosed base. Preserve that pattern; do not replace it
  with the older SMART II A GCGGG tail or the later SMART-seq2 ACAT/rGrG+G oligo.

The IS-PCR footnote says “Please do not substitute. The oligo has been modified.”
It does **not** identify a chemical modification, terminal group or modified
position. It may describe redesign relative to another reagent. Record the notice;
do not turn it into a new claim of chemical or terminal modification.

The original article material says its libraries predated commercialization and its
cDNA protocol is reflected in the named kit manual. This supports a **scoped
source-era component correspondence**, not blanket identity of every kit reaction,
lot or proprietary base. Exact IS-PCR bases are now positively sourced; complete
SMARTer-IIA chemistry remains unresolved in this pass.

## Exact Epicentre inventory

These are the original PDF's page2 sequence cells, all 5′→3′:

| Reagent | Sequence |
|---|---|
| Transposon end / ME | `AGATGTGTATAAGAGACAG` |
| Primer1 | `AATGATACGGCGACCACCGA` |
| Adaptor1 | `AATGATACGGCGACCACCGAGATCTACACGCCTCCCTCGCGCCATCAG` |
| Primer2 | `CAAGCAGAAGACGGCATACGA` |
| Adaptor2, no barcode | `CAAGCAGAAGACGGCATACGAGATCGGTCTGCCTTGCCAGCCCGCTCAG` |
| Nextera Read1 | `GCCTCCCTCGCGCCATCAGAGATGTGTATAAGAGACAG` |
| Nextera Read2 | `GCCTTGCCAGCCCGCTCAGAGATGTGTATAAGAGACAG` |
| Nextera Index Read | `CTGTCTCTTATACACATCTCTGAGCGGGCTGGCAAGGCAGACCG` |

The two entry tags are consequently `GCCTCCCTCGCGCCATCAG` and
`GCCTTGCCAGCCCGCTCAG`, followed by ME on their transferred strands. These are
**not** modern Nextera s5/s7 entry tags. The source diagram and original inventory
provide the source-era counterparts to the obsolete Epicentre family subsequently
listed in Illumina's adapter document.

The manual explicitly describes a **four-primer** reaction. The50X cocktail contains
Primer1 and Primer2, each10µM, plus Adaptor1 at0.5µM; separate50X Adaptor2 is0.5µM.
Thus the long adapter primers add complete bridge-PCR-compatible ends, and the
short universal primers amplify those ends. Do not collapse these into a modern
N5xx/N7xx two-primer recipe.

Its optional barcode design is:

`CAAGCAGAAGACGGCATACGAGAT-[BAR CODE]-CGGTCTGCCTTGCCAGCCCGCTCAG`

The manual marks the barcode as reverse complement of the sequencing read. Which
optional index member was used in a particular original SMART-seq sample remains
an unassigned study-specific choice here; no dual index or8-base index is inferred.

## Scoped inheritance and remaining boundaries

The original SMART-seq Methods material names Epicentre,5XHMW buffer,55°C for5min,
a20µl reaction and9 PCR cycles. Those are specific matches to the old vendor
family. Its modifications—approximately5ng cDNA,0.25µl transposase, PB stripping
and SPRI recovery—remain the article's choices, rather than the manual's standard
50ng/1µl/Zymo recipe. The old vendor PCR programme begins with72°C for3min **before**
95°C denaturation, then nine cycles. This supplies a source-era repair-before-PCR
operation rather than permission to draw a pre-filled duplex spontaneously.

No exact article kit catalogue number, optional index member, full proprietary TSO
sequence, or occupied target sequence was established by this pass. Do not call
those publication-wide confirmed absences. The article's alternative PE-adapter
branch and later SMART-seq2/Nextera kits remain separate.
