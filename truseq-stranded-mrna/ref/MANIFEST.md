# Source manifest

Retrieve the registered originals with:

```sh
TMPDIR="$PWD/_data/tmp" python3 tools/get_sources.py --slugs truseq-stranded-mrna --no-papers
python3 tools/scrape_primers.py _data/sources/truseq-stranded-mrna
```

`get_sources.py` validates downloads and creates text twins. Check the returned original's title/revision and `MANIFEST.tsv` URL before reading or attributing it to a historical version.

- **1000000040498 v00, October 2017**, *TruSeq Stranded mRNA Reference Guide*: registered official PDF, local `illumina_truseq-stranded-mrna-workflow-reference-1000000040498-00.pdf(.txt)`.
- **1000000002694 v22, September 2025**, *Illumina Adapter Sequences*: registered official PDF, local `illumina_1000000002694_22_ilmn-adapter-sequences.pdf(.txt)`. Its positive CD templates support a conditional reference, not an inferred historical lot identity.
- **1000000041074 v17**, current *Index Adapters Pooling Guide*, [TruSeq Kits](https://support-docs.illumina.com/SHARE/IndexAdaptersPooling/Content/TruSeqPoolingGuidelines.htm): local `illumina_TruSeqPoolingGuidelines.htm(.txt)`. The body was recovered through the registered route after the sandbox DNS failure; there is no remaining content-retrieval barrier for this topic.

The actual older pooling annotation, named supporting manuals/bulletin, software planning documentation and customer-sequence-letter endpoint are recorded as pending linked-source scope in the review. The `150XXXXX` checklist target is a literal source placeholder, not a known absent document. No acquired run, historical adapter revision or selected D-index member is implied by downloading a source.

No third-party documents or extracted text are stored in git. Source tests skip when these cached originals are absent.
