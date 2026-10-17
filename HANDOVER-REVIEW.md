# Review handover: uncommitted edits since 8176a87

This is for whoever reviews the working tree against the last commit, `8176a87` (2026-10-17 12:49).
Delete this file once the edits are reviewed and committed.

## Baseline

When this session took over, `git status` showed only three entries:

- `M meta3c__10.7554+eLife.03318/tools/chem.py`
- `M meta3c__10.7554+eLife.03318/tools/selftest.py`
- `?? AUDIT.md`

The two meta3C files came from the previous owner, and they did not import at that point. They were then fixed in this session. `AUDIT.md` was already the curation worklist, untracked. Every other uncommitted change comes from this session: one coordinator plus parallel sub-workers, none of whom committed.

## How to see everything

```sh
git diff HEAD -M --stat                     # staged + unstaged, renames detected
git ls-files --others --exclude-standard    # 22 untracked files (new ledgers, AUDIT.md, ...)
```

Some changes are staged, so use `git diff HEAD`, not plain `git diff`. The staged changes are `git mv` directory renames and three deleted `lib/` modules; nothing else is staged. At hand-off the diff covers about 1230 files with roughly 57k insertions and 18k deletions. About 820 of those files are in protocol directories and 320 in `catalogue/source_reviews/`.

Sources are not in the repository. They live under `/husky/henriksson/NGS101_papers/sources/<slug>/`, so set `export CHEM_DATA=/husky/henriksson/NGS101_papers`. Several fixes depend on files that were fetched into those bundles, and checks skip when the files are absent.

## What changed, by kind

`AUDIT.md` is the narrative record. Its sections are "Ownership and outstanding checks", one section per round, "Deferred protocols" and "Shared-code and catalogue follow-up queue". The machine record is in `catalogue/curation_rounds/*.tsv`, with one row per protocol pass, and in the `bounded_checks` entries appended to `catalogue/source_reviews/<slug>.json`.

1. **Capped audit passes over all 410 protocols** (round 5 completion, rounds 6–21, deferred follow-ups 01 and 02).
   - 375 protocols are inspected: 368 with corrections and 7 checked with no change.
   - 35 are deferred because their defining original could not be obtained (table in `AUDIT.md`).
   - Corrections went into protocol modules (`*/tools/*.py`), notes (`0*.md`), `ref/MANIFEST.md` files and review JSONs.
   - Review points:
     - Every changed base should trace to a source line cited in the ledger or the note.
     - Inferred material should stay marked inferred.
     - Appended review-JSON fields should not rewrite earlier passes.
2. **Shared-code fixes from the queue** (fix rounds 1–3, file-disjoint packages):
   - `tools/scrape_primers.py`, `tools/source_oligo_check.py`, `tools/selftest_source_oligos.py`: guard and scraper parsing and locator rules.
   - `lib/mdfacts.py`: module lookup checks the note's own directory first.
   - `tools/get_sources.py`, `tools/fetchers/*`, `tools/doctext.py`: fetching and validation.
   - `lib/seqprimers.py`: read-table spans for windows of unstated length.
   - `lib/batch_ngs.py`: the single-index TruSeq P5/Read 1 overlap.
   - `lib/adapter_ligation.py` and `lib/dumbbell.py`: inferred adapter ends.
   - `lib/rna_dna_bridge.py`: optional RNA positions in the bridge duplex.
   - `lib/multimodal_spatial.py` and other helpers: neutral default labels and explicit state hooks.
   - `protocol_families/*` and `lib/`: legacy recipes retired or moved into their owning protocol. `lib/split_pool_atac.py` became `spatac-seq__…/tools/spatac_model.py`, and `lib/antibody_tags.py` and `lib/spatial_epigenome.py` were deleted.
   - `tools/lint_p5_read1.py` (new): a P5/Read 1 junction lint, with exemptions that must give a reason.
   - Review points:
     - Shared code must stay protocol-neutral (`AGENTS.md`).
     - New behaviour is covered in `tools/selftest.py` and `tools/selftest_sources.py`.
3. **Catalogue identity.** These directories were renamed through `catalogue/ours.tsv`:
   - `chip-seq__10.1038+nmeth.1068` → `chip-seq__10.1038+nmeth1068` (the dotted DOI never resolved).
   - `ribomethseq__10.1002+anie.201408362` → `ribomethseq-illumina__10.1093+nar+gkw547`. The page draws Marchand 2016 (Illumina); Birkedal 2015 (Ion Torrent) is kept as the undrawn original.
   - Four preprint-named directories → their published versions: easySHARE-seq (eLife), GoT-ChA (Nature), Phospho-seq (Nat Commun), sciPlex-ATAC-seq (BMC Genomics).

   Each rename updated, in the same step:
   - the review JSON;
   - the ledgers and `AUDIT.md`;
   - the `extra_sources.tsv`, `properties.tsv`, `protocol_oligo_audit.tsv` and `paper_overrides` rows;
   - `citations.tsv`, via `tools/fetch_citations.py`, which refreshed every row's count and date;
   - `scg_lib_structs.tsv`, regenerated;
   - the external source directory.

   Old names may remain in `scg_lib_structs.tsv` in the `slug` and `paper_url` columns. Those columns describe the upstream scg_lib_structs paper, and the `our_dir` column maps it to the new directory.
4. **Catalogue tables.**
   - `extra_sources.tsv`: about 190 rows added and about 25 changed or retired (with `# retired` comments).
   - `chemistry_decoupling.tsv`: caller lists recomputed from static imports, and retired symbols removed.
   - Smaller edits in `properties.tsv` (smol-seq), `oligo_source_audit.tsv` and `protocol_oligo_audit.tsv` (inDrop v1 is now verified against Cell 2015 Table S5).
5. **Privacy.** The maintainer's email was embedded in the fetchers' User-Agent (`tools/fetchers/base.py`) and in Crossref/NCBI query parameters (`catalogue/tools/fetch_scg_lib_structs.py`). It is now removed; `grep -rn "he\.johan"` should find nothing in tracked code.
6. **Polite access.** `tools/fetchers/base.py` now throttles requests per host across processes (8 s for publishers). After a challenge or HTTP 429 it puts the host in a six-hour cool-down; the markers are in `$CHEM_DATA/sources/.cache/ratelimit/`. `AGENTS.md` documents this and states that challenges are never worked around.

## Things worth a skeptical look

- **Bulk worker edits.** About 50 parallel workers edited protocol modules under a 10-minute, three-source cap. They were told to verify against sources and to keep inferred values inferred. Spot-check transcriptions against the cited source lines, especially where a value changed: the GRID-seq extra G, ISSAAC dT32 → T30VN, ProBac probe length, MARS-seq N anchor, scifi representative barcode wells (scifi uses one representative well's A/A bases; other wells differ, as its note says), and Microwell TSO.
- **Second readers.** Review JSONs record worker passes. Where a JSON shows two passes, the second may also be a worker. The two-independent-reader gate for confirmed-missing claims is not met by any capped pass.
- **Old pass text kept on purpose.** Some review JSONs store a later pass under an extra key such as `pass_3_second_reader_2026-10-17`, because both pass slots were already full. The audit validator does not read those keys.
- **Read-table wording.** The read-table rewording (`lib/seqprimers.py`) changed the generated text on seven pages: disco, olink-explore, sci-space, scone-seq, sprite, taps and umi-4c.
- **Self-corrected worker mistakes.**
  - One worker staged with `git add -N .` and then removed exactly those intent-to-add entries.
  - One worker ran `git mv` and then reset it, so `lib/split_pool_atac.py` shows as deleted and `spatac_model.py` as untracked.
  - One worker ran a broad `pkill -f tesseract`.

## Verification recorded at the previous hand-off

```sh
export CHEM_DATA=/husky/henriksson/NGS101_papers
for t in */tools/selftest.py; do python3 "$t" >/dev/null || echo "FAIL $t"; done   # 402 suites
for d in $(grep -l PAGE_FILE */tools/*.py | cut -d/ -f1 | sort -u); do python3 tools/protocol.py test $d >/dev/null || echo "FAIL $d"; done  # wrapper-free
python3 tools/selftest.py            # 1262 checks
python3 tools/selftest_sources.py    # 222 checks
python3 catalogue/tools/selftest.py  # 83 checks
python3 tools/check_chemistry.py     # contradictory=0
python3 tools/audit_protocol_sources.py --check $(ls catalogue/source_reviews/*.json | xargs -n1 basename | sed 's/\.json$//')
python3 tools/selftest_source_oligos.py   # slow; run in chunks of directories: 1000 checked, 5 failed
python3 tools/lint_pairing.py        # 0 problems (after build)
python3 tools/lint_p5_read1.py       # 0 findings, 1 exempted (PEM-seq, source-printed)
python3 build_docs.py --check        # all up to date
python3 build_site.py                # 407/407 diagram pages, links resolve
```

The 5 source-guard failures recorded at that hand-off:

- **Atrandi WGS (3):** the vendor guides are behind a reCAPTCHA.
- **ChAR-seq bridge (2):** the bridge is printed only in a figure and the bioRxiv supplement. Whether that counts as a source is an open maintainer decision.

## Still open

These are listed in `AUDIT.md`:

- the follow-up queue: maintainer decisions, the legacy-recipe migration, blocked retrievals and round-3 residuals;
- the 35 deferred protocols, which need originals downloaded by hand;
- "How to resume".

## Critical diff review — 2026-10-17

The baseline/handoff measurements above are historical. This follow-up inspected the
high-impact shared changes: source guard provenance and variable strands;
protocols.io identity/revision handling; fetcher cache identity and access handling;
read-cycle calculations; graph branch rendering; sequencing-plan rendering;
canonical adapter overlap; note module lookup; and catalogue DOI/directory changes.
It did not independently reread all 410 protocols or complete the legacy migration.

Established defects corrected:

- The protocols.io DOI-answering cache path now checks its recorded digest; export
  acceptance checks consistent DOI/Version labels and the cited accession. Newly retrieved wrong-accession or
  contradictory exports cannot supply attachments. Existing cache files are not
  overwritten to satisfy a different identity.
- The oligo guard no longer automatically accepts a linked bioRxiv/medRxiv version
  when the publication is available. Published byte-identical supplement aliases
  remain eligible, and Genome Research `10.1101/gr.*` is not treated as bioRxiv.
- Read-content and feature-cycle arithmetic now uses source-declared molecular
  lengths instead of schematic widths. Features beyond an unknown-length insert
  are not assigned exact cycles.
- A contributing branch's preparation is rendered even when another template
  shares the output's drawn state.
- Mixed exact/unresolved sequencing plans receive separate terminal sections;
  an exact branch cannot hide a neighbouring unresolved branch or supply its sites.
  Pages whose plans are all unresolved retain their existing source-declared view.
- Small API/HTML requests trigger the same challenge cool-down as downloads.
- Index check summaries require a successful process exit; a printed pass banner
  followed by a failure cannot count as passing.

Evidence decision: ChAR-seq does not need a preprint exception. I visually read the
published `elife-27024-fig1-figsupp1.jpg` inside the existing PMC5962340 supplement
ZIP and verified the model's two ordering strands and App/Sp3/TEG labels. The local
transcription checks already pin both strings. A bridge-only image exemption is
now explicit in the source guard and in its output; neither strand is counted as a
text match. The preprint checklist row is historical/not relevant; prior passes
remain intact. Exact TEG linkage and the historical NEBNext kit remain unresolved.

Scoped source checks covered ChAR-seq, Atrandi, scifi-RNA, scifi-ATAC, scRRBS,
scTrio-seq, Spatial-CUT&Tag, GRID-seq, ISSAAC-seq, ProBac-seq and Microwell-seq:
63 strands checked, 3 failed (the unchanged Atrandi locator exclusions), 3 strands
exempted (ChAR-seq x2 and the existing scifi exception), 6 inputs unrepresented,
0 source-bundle skips. Published GRID-seq Methods line 129 and ISSAAC Supplementary
Table 4 rows TruseqR1/R2_oligo_dT were also spot-read: the corrected GRID mixed
RNA/DNA linker and ISSAAC T30VN remain consistent with those originals. These
spot checks establish their stated transcription scope, not complete chemistry.

Atrandi's hand-migrated V3 guides are present and their text contains the three
strings, but their exact original retrieval URLs are unlogged and the source
locator filter excludes them. No blanket trust rule or source exception was added;
those three guard failures remain visible. The 35 deferred defining-original passes,
all-unresolved multi-terminal presentation and the remaining vocabulary/migration
queue remain open. No missing-information claim was accepted.

Final verification on the reviewed working tree:

- 401 legacy protocol suites and 6 wrapper-free suites passed (process exit checked).
- Shared suite: 1,272 checks; source-fetcher suite: 226; catalogue suite: 83.
- All 410 review JSONs pass the schema gate; historical review statuses unchanged.
- The scoped source guards above completed without timeouts. Atrandi's three known
  locator failures remain failures; exemptions and unrepresented inputs are separate.
- Default chemistry pilot checks: 39 established, 39 conditional, 157 unresolved,
  0 contradictory claims. This is the tool's pilot scope, not all 410 recipes.
- Full final site rebuild: 407/407 diagrams, 402/402 documented NGS graph contracts,
  all relative links resolve. Mixed-terminal output spot-checked on ARTIC,
  Multi-CUT&Tag and spatial-CUT&Tag–RNA-seq.
- Pairing lint: 0 problems. P5/Read 1 lint: 402 modules scanned, 0 findings,
  1 existing source-qualified PEM-seq exemption.
- Whitespace checks are clean for the files changed in this review. Existing
  whitespace findings elsewhere in the larger diff were left alone.

No changes were staged or committed during this review. Generated HTML/site output
is ignored; all original source documents remain outside the checkout.

## After the critical review: deferred follow-up 03 (2026-10-17, evening)

The coordinator session's work overlapped with the review above. Everything below changed
after the review's measurements. Check it the same way.

- **Hand-downloaded originals.** The maintainer downloaded 24 files, and the coordinator moved them into
  `$CHEM_DATA/sources/<slug>/`, after checking each one's printed DOI, title or authors. Each has
  a `MANIFEST.tsv` row with its original URL and the note "hand-downloaded by the maintainer
  … identity checked". Four exact duplicates were discarded.
- **Ledger.** `catalogue/curation_rounds/2026-10-17-deferred-followup-03.tsv` re-audits the 31
  protocols whose defining originals are now on disk (29 `checked_with_corrections`). The
  per-protocol paragraphs are in the `AUDIT.md` section "Deferred follow-up 03". Workers edited protocol
  modules, notes and review JSONs (appended `bounded_checks`). No file was reported misplaced.
- **Still deferred (6):** CapSeq, HiRES, Microbe-seq, ONT direct cDNA, TUC-seq (main article
  missing) and RiboMethSeq (the Birkedal 2015 article is missing; its defining Marchand 2016 source is read).
- **Coordinator catalogue edits:**
  - The `ours.tsv` CH-ATAC note.
  - The ChIP-seq `protocol_oligo_audit.tsv` row.
  - A sci-CAR `pmccloud:PMC5866720` row in `extra_sources.tsv`.
  - `chemistry_decoupling.tsv` caller lists: cfMeDIP and scCARE-seq no longer call legacy recipes.
  - The inDrop v1 note's stale "v2 main text not fetched" statements now point to the v2 note.
- **Verification after both sets of edits:**
  - All 402 suites and the 6 wrapper-free suites pass. Shared checks 1274, source fetcher 226, catalogue 83.
  - The full source guard checked 1027 strands: 3 failed (Atrandi), 31 exempted with stated
    reasons, 1001 inputs without exact states.
  - Chemistry 0 contradictory; all review schemas pass.
  - Site 407/407, links resolve; both lints clean.

  Review the exemption list (`selftest_source_oligos.py` prints each one with its reason). Each
  one is a protocol-declared `SOURCE_OLIGO_EXCEPTIONS` entry, not a text match.


## Critical re-review after deferred follow-up 03 — 2026-10-17

Screened the latest 31-protocol batch, concentrating on source promotion, adapter
modifications, acquired reads and product fates. Ten protocols needed corrections;
the scope and findings are in
`catalogue/curation_rounds/2026-10-17-critical-rereview-followup-03.tsv`
and appended per-protocol bounded checks. This was not a complete independent
reading of all 31 source bundles or an absence audit.

- TAIL-seq: the 3′ adapter now includes the two biotin-dT nucleotides and terminal
  ddC as TTC, with the blocked end shown. Native tail lengths are variable. RNA
  intermediates are single-stranded; RT has a distinct limited cDNA state and PCR
  geometry refers to that actual input. The final amplicon does not acquire the
  modified suffix beyond the RT binding arm. The shared scraper now retains a
  nucleotide-bearing terminal modification run; a regression covers this failure.
- SMART-seq2: removed acquired R1=43 inferred from demultiplexed/deposited reads.
  SureSelect XT HS2: the two source-known non-MBC read-start bases have explicit
  length 2 despite their unknown identities.
- SMART-seq2, sci-CAR, CoBATCH, Slide-seq and scCARE: candidate non-selected end
  classes have unresolved amplification/recovery, rather than invented physical
  discard. The legacy SMART-seq2 end-class validator still uses MaterialKind.WASTE
  for its non-selected output; its terminal explicitly states recovery unresolved.
- inDrop v2: ExoI hydrolysis products are distinguished from physical discard.
  Single-stranded ancient DNA: released original templates stay unresolved after
  entering PCR. cfMeDIP: both denatured strand polarities stay in the same tube and
  enter the control/IP split; the complementary polarity is not discarded.

Verification: all 30 applicable local suites pass (Simple-seq is notes-only);
all 31 scoped chemistry checks have zero contradictions; source guards check
66 exact strands with zero failures, 8 declared exemptions and 60 inputs without
exact states. No source-guard timeouts. Shared 1275 and source-fetcher 226 checks
pass. Review schemas and both pairing/P5 lints pass. The cfMeDIP correction initially
failed the two-output strand-separation contract; both outputs are now retained
and rejoined as a population-accounting step, and its page/guard/chemistry pass.
Full-site rebuild result is recorded below after the final run.

Earlier source deferrals, including the six unobtained defining originals and
Atrandi's three historical locator failures, were not resolved by this review.
No files were staged or committed.

Final rebuild: **407/407 diagram pages built**, 402/402 documented NGS graph contracts pass, and every relative site link resolves. The index reports 5525 checks total. All 777 Markdown outputs were current before this final verification entry; rebuilt after adding it. No protocol page is omitted.
