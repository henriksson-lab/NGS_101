# Protocol source reviews

One `<protocol-dir>.json` per audited protocol. This is a source-reading record,
not a store for third-party material. Run `python3 tools/audit_protocol_sources.py`
for the repository census or `--check <protocol-dir>` for one review. Missing files
remain pending; historical rows in the older oligo ledgers do not pass this gate.

The checklist includes the defining source, every supplement and detailed method
relevant to a missing claim, vendor documents and predecessor protocols named by
those sources, and chemistry-relevant linked material. A source can be `available`
(a validated original is retrieved but reading is pending), `read`, or `unavailable`,
and `relevant` or `not_relevant`. Record the local original's locator for an
available source; downloading it supplies no reading pass. State a reason for either kind of
exclusion. `linked_information.result` says `enumerated`, `none_found`, or `blocked`;
both reviewers describe how they checked links. A linked item is another source
row with `role: "linked_material"` and `discovered_from` naming its referring row.

Name an identified original document, repository or file in each source row.
A desired reagent specification or an unassigned experimental index is an open
question, not evidence that a disclosure document exists but is inaccessible.
Keep those questions separate from actual unavailable originals. A cited source
whose identity or version remains unresolved still blocks its missing-information
claim. Reclassifying a question never establishes absence or supplies a reading pass.

The metadata-only [identity queue](../source_review_identity_queue.tsv) flags rows
needing this distinction. Its questions are not source-content findings, accepted
exclusions or completed reviews.

Two independent people or agents check the **original** source text separately,
including sources judged irrelevant to the missing claim so that exclusions do
not depend on one reader.
Each pass records a reviewer, date, exact locations inspected, and findings. A
scraper run can point a reviewer to a location but is never a pass itself. The
second reviewer checks the source afresh, including the linked-material inventory.
If any relevant source or link is inaccessible or unexamined, its
missing-information conclusion remains `blocked`. For the overall review status,
use `in_progress` while the next action is ordinary downloading, reading or an
independent check. Use `blocked` when an attempted next action establishes a
specific access barrier or an unresolved source identity, and record the attempted
route and reason. A failed downloader request alone does not establish that the
document is inaccessible: check the original publisher, repository and available
author-copy routes. Distinguish a network/environment failure from a publisher
restriction. These work statuses never relax the two-pass acceptance gate.
Only `accepted_missing` claims confirmed absence, and it also requires the public
page to display the checked source scope and linked-information disposition.
`reviewed` means two passes found no accepted missing claim. `in_progress` allows
partial entries; it cannot establish absence.
Named source or link obligations may be tracked in `pending_sources` and
`pending_linked_targets` during an audit. Either nonempty queue prevents a
`reviewed` or `accepted_missing` status, even if every current checklist row has
two passes. Resolve each item against its original and retain the located checks
in the source checklist before clearing the queue. Removing a pending item without
that evidence does not establish source completeness.

Example shape (replace every value with a located, protocol-specific review):

```json
{
  "dir": "example__10.0000+example",
  "status": "in_progress",
  "source_checklist": [
    {
      "id": "paper",
      "role": "defining",
      "locator": "DOI and section or file name",
      "availability": "read",
      "relevance": "relevant",
      "pass_1": {
        "reviewer": "reviewer-1",
        "date": "YYYY-MM-DD",
        "locations": "Methods paragraphs 3-7; supplement references",
        "finding": "What was and was not printed"
      },
      "pass_2": {
        "reviewer": "reviewer-2",
        "date": "YYYY-MM-DD",
        "locations": "Independently checked Methods paragraphs 3-7",
        "finding": "Independent finding and any disagreement"
      }
    }
  ],
  "linked_information": {
    "result": "none_found",
    "pass_1": {
      "reviewer": "reviewer-1",
      "date": "YYYY-MM-DD",
      "locations": "Methods citations and supplement links",
      "finding": "No other chemistry-relevant linked material"
    },
    "pass_2": {
      "reviewer": "reviewer-2",
      "date": "YYYY-MM-DD",
      "locations": "Independently followed Methods citations and links",
      "finding": "No linked protocol or vendor guide found"
    }
  },
  "missing_claims": [],
  "conclusion": "Review in progress"
}
```

For `accepted_missing`, add each missing claim with `claim`, `disposition` set
to `confirmed_missing`, and `evidence`; add `public_page` (relative HTML path)
and `public_checklist_marker` (exact visible text in that page). The validator
rejects acceptance if any relevant source is unavailable, the two reviewers
are the same, linked information is blocked, or the public checklist is absent.
