# Yeast Toolkit: parts, cassettes and multigene plasmids

Lee ME, DeLoache WC, Cervantes B, Dueber JE (2015).
[A Highly Characterized Yeast Toolkit for Modular, Multipart Assembly](https://doi.org/10.1021/sb500366v).
ACS Synthetic Biology 4:975–986.

## Source scope

Bounded readings, 2026-10-21; the direct primary-type cassette route is now modelled.
The published article PDF was retrieved from the authors' laboratory website:
`Pub_2015_Highly-Characterized.pdf`. Its own printed DOI is 10.1021/sb500366v.
This author-hosted copy has unfinished journal pagination (`XXXX`); its exact
publication-version identity is not yet established. Exact assembly
boundaries below are now checked against the publisher's separately identified SI.
Read Figure 1 caption and Results “Definition of an Assembly Standard for Yeast”,
local PDF text lines 120–190. Read the article's Supporting Information inventory
at lines 930–942 to identify the next source obligations.

The initial oligo-scraper run returned zero candidate sequences before the
supplements were retrieved. No sequence-absence inference follows. The expanded
published supplement bundle was subsequently scraped and read as described below.

## Source-supported architecture

🟢 The authors adapt modular cloning to yeast using numbered part types. Source
DNA enters a universal part vector by BsmBI assembly. Typed parts are combined
using BsaI into a cassette containing a transcription unit. Cassettes can then be
combined using BsmBI into a multigene plasmid.

🟢 Assembly Connectors, part types 1 and 5, determine cassette order at the higher
assembly level. They are molecular parts with distinct boundary requirements, not
display-only labels or implicit separators supplied by a shared family recipe.

🟢 The article distinguishes selectable vector backbones at successive assembly
levels. The graph must distinguish transferred cassettes from vector backbones,
selection markers and removed dropout material. The whole input plasmid does not
become the insert sequence of its product.

The original MoClo protocol uses a different enzyme progression. Reuse of the
assembly idea does not justify importing MoClo's BpiI step, fusion-site vocabulary
or vector geometry into this protocol.

## Published supplementary follow-up

Retrieved the publisher's two supplement records through its public Figshare API:

- Supplement DOI `10.1021/sb500366v.s001`, article 2054070 version 1, file
  `sb500366v_si_001.pdf`. The retrieved original matches publisher MD5
  `b77fcc4261f4f2954586ccb9f8c9f72c`.
- Supplement DOI `10.1021/sb500366v.s002`, article 2054064 version 1, file
  `sb500366v_si_002.zip`. The retrieved archive matches publisher MD5
  `260662b165a4d6b402015c95b202ba9f`. It includes the original `YTK_Parts.xls`
  sequence workbook and the GenBank plasmid records.

Ran the oligo scraper again on the expanded bundle. Read the SI's part-type
overview (original PDF page 12, Figure S6), hierarchical-assembly description
(pages 20–22), multi-gene example (page 22, Figure S21) and “Differences from
MoClo” (page 23). Inspected the sequence workbook's headings, representative
part identities and sequence fields. Individual GenBank maps are retrieved but
not independently read in this pass.

🟢 Figure S6 gives the primary BsaI part boundaries on the reference strand:

| Part type | Upstream fusion | Downstream fusion |
|---|---|---|
| 1 | CCCT | AACG |
| 2 | AACG | TATG |
| 3 | TATG | ATCC |
| 4 | ATCC | GCTG |
| 5 | GCTG | TACA |
| 6 | TACA | GAGT |
| 7 | GAGT | CCGA |
| 8 | CCGA | CCCT |

🟡 Independently computed BsaI cuts in the published workbook's representative
pYTK002, 009, 032, 051, 067, 074, 081 and 083 sequences agree with those boundaries.
This checks the named source sequences; it does not establish cloning yield or
compatibility of every possible part combination.

🟢 The SI explicitly states that its fusion vocabulary differs from original
MoClo. Subtypes and composite types are separate choices and must not inherit
the primary boundaries by name alone.

🟢 Figure S21 uses ConLS/ConR1, ConL1/ConR2 and ConL2/ConRE to order its example
transcription units. The receiving vector uses ConLS′/ConRE′. These reversed
connectors are distinct molecular inputs. The final connector junctions retain
barcode scars used for assembly verification, rather than sample demultiplexing.
Their exact retained sequence and structured identifier meaning remain to be assigned.

🟢 The preassembled cassette-vector branch retains BsaI sites in the Type 234
GFP dropout. It is not interchangeable with the direct part-to-cassette route.
The SI describes different reaction handling and incorrect concatenated-plasmid
products in that branch. Both the retained sites and those fate distinctions must
remain visible if that branch is modelled.

## Output and implementation queue

The current public schematic models the direct eight-primary-part cassette route,
using the published workbook identities pYTK002, 009, 032, 051, 067, 074, 081 and
083. These illustrate source-defined architecture; this combination is not asserted
to be an experimental clone. Type 8 supplies the AmpR–ColE1 backbone, not an extra
ninth receiving vector. BsaI cleavage, cohesive-end pairing/sealing, and candidate
colony selection are explicit graph steps. Donor backbones and the mRFP1 dropout
are independent molecular classes. Selection does not prove a clone's genotype.

The shared operation checks only the declared cohesive ends and circular strand
order. Unknown part cores keep the assembly conditional. Release from the complete
source plasmids is not established by this relation, and the model makes no claim
about yield, gene expression, off-product frequencies or pooled representation.
No NGS library preparation or Illumina read recipe is assigned.

🟢 Follow-up reading of the original SI Figures S7 and S15, including the rendered
page 13, establishes the different ConLS and ConR1 internal architectures. The
declared 143-bp barcode concatenations and ConLS 21-bp cassette barcode have
structured inline-barcode features; they are not UMIs or demultiplexing indexes.
Their known molecular lengths and recognition/fusion motifs are independently
checked against the published workbook sequences. Barcode bases are not expanded
in the diagram. This does not assign the later multigene scar sequence or length.

🟢 Published discrepancy, checked in rendered SI page 18: Figure S16B prints Type
6 with TGGC at the left and CGAC on the lower right strand (corresponding to GCTG
on the reference strand). Figure S16A, its text, and overview Figure S6 instead
give TACA/GAGT. The model uses TACA/GAGT, also supported by independent BsaI cuts
in published pYTK074. The public page identifies these locations and values. The
authors' intended correction to panel B has not been confirmed.

🟢 The author-hosted article's Methods “Golden Gate Assembly Protocol”, text
lines 723–739, specify T7 DNA ligase and T4 DNA ligase buffer as distinct reagents.
Because that copy's exact publication-version identity is unverified, the current
published-SI-scoped graph leaves ligase identity unassigned. This is a source
identity question, not a claim that the publication omits the reagent.

- The published SI and sequence archive are available and identified. Read the
  remaining detailed connector/subtype figures and source-owned vector records
  needed for each selected branch before defining its retained sequences.
- Resolve the final article's identity before assigning its detailed reagent choice.
- Add part-entry, cassette-vector and cassette-to-multigene branches locally with
  their own retained sequences, geometry and molecular fates. They are not covered
  by publication of the direct cassette route.
- Keep construct sequencing/QC separate from assembly and expression claims.
- Do not infer yields, expression levels or representation of combinatorial products
  from compatible part boundaries.
- A second independent reader and complete linked-source inventory remain pending.
  This note makes no confirmed-missing-information claim.
