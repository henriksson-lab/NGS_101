# Linking T-cell receptor sequence to functional phenotype at the single-cell level

Defining publication: https://doi.org/10.1038/nbt.2938.

## Molecular distinction and scope

🟢 A sorted single T cell is assayed for TCRα/TCRβ and cytokine/transcription-factor transcripts. Separate amplicons are associated computationally through a common plate–row–column cell identity; they are not covalently fused. Source is the published study's final accepted manuscript in PMC and the publisher's published supplement. This is a bounded source pass, not a completed disclosure audit.

## Sources inspected

- Defining article, Online Methods: amplification reactions, primer design and sequencing analysis (XML text lines112–133).
- Published supplementary PDF MOESM24: SupplementaryFigure1 and its caption; rendered original page1 independently checked against the text layer. SupplementaryTables1/3 are present, including appended OCR, but their full pools were not transcribed or validated here.
- Article links Corrigendum PMC5800857 (2015,33:210; PMID25658286). Direct PMC retrieval returned a browser challenge. Its applicability remains unresolved; no claim that the correction leaves chemistry unchanged.

## Reaction and identities

🟢 OneStep RT-PCR (Qiagen) first reverse-transcribes at50°C30min and amplifies the multiplex TCR/phenotype panel25cycles. Two separate1µl aliquots enter TCR nestedPCR25cycles or phenotype nestedPCR35cycles (HotStarTaq). Oneµl of each nested reaction enters the corresponding barcodingPCR36cycles. The model keeps transferred and untransferred material distinct; the source does not specify the remainder's disposal.

🟢 The5′ barcoding primer carries plate5 and row5 components; the3′ primer carries column5. All three together identify the original cell. They are neither UMIs nor sequencer sample indices. The model draws these as separate components, never as a fictitious contiguous strand. Tumor/colon combinedα/β TCRPCR is the represented branch; the separateα/β variant is deferred.

🟢 PCR products are combined in equal volumes, run on1.2%agarose, and350–380bp gel material is extracted using QIAquick. MiSeq paired-end sequencing is reported. Run-cycle allocation and complete primer landing remain unassembled.

## Source literals and discrepancies

🟢 Common5′ handle: `CCAGGGTTTTCCCAGTCACGAC` (22bases). Phenotype3′ handle: `AGCGGATAACAATTTCACACAGGA`. Model transcription tests read the published supplement.

The main manuscript calls the former a23-base common sequence, whereas its printed literal and rendered SupplementaryFigure1 contain22bases. This model retains the literal; no additional base is invented. Main Methods also describes Cα as threefold excess while the paired numeric concentrations point in the opposite direction (0.15µM versus0.5µM respectively); that quantitative ratio remains unassigned. These are source-text discrepancies pending the linked corrigendum, not silently repaired publication values. The private preview repeats the notices.

## Deferrals and lineage

Full study-specific primer panels, barcode whitelist transcription, complete final strands/read primers, linked corrigendum and the article's cited three-reaction predecessor remain follow-up targets. Unknowns are unassembled, not confirmed absent. Source supports conceptual nested RT-PCR ancestry, but no exact predecessor component import is authorized in this pass.

## Completion and verification

Private incomplete source-local graph. Required source inventory explicitly retains unassembled full-strand chemistry; a valid local reaction contract does not establish source-review completion. Mandatory primer scraper was run against fetched sources and relevant hits checked in originals. No two-reader missing-information pass has been completed.
