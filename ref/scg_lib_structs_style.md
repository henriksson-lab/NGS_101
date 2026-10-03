# scg_lib_structs page style (reference)

Target format for the Atrandi single-cell WGS chemistry page, reverse-engineered from
<https://teichlab.github.io/scg_lib_structs/methods_html/SPLiT-seq.html>
(source + stylesheet archived alongside this file as `SPLiT-seq.html`, `page_format.css`).

## Mechanism

Hand-written HTML, no JavaScript, no classes on sequence text. Unknown element names are
used as inline color spans, defined in `../style_related/page_format.css`.

| Tag | Color | Meaning |
|---|---|---|
| `<p5>` / `<p7>` | `#08519c` / `#a50f15` | Illumina P5 / P7 flowcell adapters |
| `<s5>` / `<s7>` | `#6baed6` / `#fc9272` | Nextera s5 / s7 primer entry points |
| `<me>` | `#969696` | Tn5 mosaic end |
| `<t7>` | `blue` | TruSeq read-2 / small-RNA scaffold |
| `<cbc>` | `#f768a1` | cell barcode |
| `<umi>` | `#807dba` | UMI |
| `<r1>` / `<r2>` / `<r3>` | `#4a1486` / `#6a51a3` / `#807dba` | round-1/2/3 linker regions |
| `<tso>` | `#2ca25f` | template-switch oligo |
| `<w1>` | `#f03b20` | Drop-seq style W1 adapter |
| `<bulge>` `<nexus>` `<cs1>` `<cs2>` `<hairp>` `<lstem>` `<ustem>` `<pe1>` `<pe2>` `<pe3>` `<doab>` | various | method-specific extras |

Typography: body Computer Modern Sans; all sequence Fira Mono.

- `<seq>` — 0.9em, `text-indent: 2.7em`. Used for the oligo list.
- `<align class="small">` — 0.9em, `line-height: 65%`.
- `<align class="long">` — 0.8em, **`line-height: 70%`**. The crushed line-height is what
  makes stacked duplex strands read as one ladder.

## Document skeleton

```
<h1>  method name (with #anchors if several variants share a page)
<p><info>  prose preamble: what the method is, which publication the walkthrough
           follows, link to the supplementary oligo table, protocol-version caveats
<h2>  Adapter and primer sequences:
      <seq> + one <p> per oligo, named exactly as the publication/manufacturer names it
<h2>  Step-by-step library generation
      <h3>(1) ... (2) ...</h3>  each followed by <pre><align class="long"> ... </align></pre>
<h3>  (N) Final library structure:
<h2>  Library sequencing:
      <h3>(1) Read 1 ... (2) Index 1 ... (3) Read 2 ...  with template strand + cycle counts
```

Numbered `<h3>` headings are full imperative sentences naming enzyme and intent, e.g.
*"Anneal Oligo-dTVN and Oligo-randN primers to mRNA and reverse transcription using
Maxima H Minus Reverse Transcriptase `<i>in situ</i>` to add Round 1 barcodes"*.

Every diagram is `<pre>` wrapping `<align>`, so alignment is literal whitespace —
columns are kept by hand, counting characters.

## Drawing conventions

- Top strand 5'->3' left to right as `5'- ... -3'`; bottom strand beneath as
  `3'- ... -5'`, i.e. written right-to-left so complementary bases share a column.
- Intermediate products label the real 5' modification where it sits:
  `/5Phos/-5'`, `/5Biosg/-5'`.
- `---->` / `<----` mark polymerase extension and its direction; arrow length is padded
  to reach the template.
- `--|` at a strand end = tethered to a streptavidin bead.
- Insert is `XXX...XXX`; poly-A as `(pA)` or `(A)<sub>n</sub>`, complement `(dT)`;
  degenerate bases verbatim (`B`, `V`, `N`, `NV(T)15`).
- Placeholders in the oligo list are bracketed prose —
  `<cbc>[8-bp Round2 barcode]</cbc>` — collapsing to literal `NNNNNNNN` in the final
  library diagram.
- `<i>...</i>` captions label sub-products inside a `<pre>` block, e.g.
  *"Product 5 (s5 at one end, 3' of cDNA at the other end, the only amplifiable fragment)"*.
  Tagmentation is drawn as **all five** Tn5 products, each annotated with whether it can
  amplify — a big part of why these pages convince.
- The final-structure block carries two or three annotation lines under the duplex,
  aligning region names and lengths beneath their columns:

```
  <p5>Illumina P5</p5>   ...  <cbc>8 bp</cbc>   <r2>Round2 linker</r2>   <cbc>8 bp</cbc>  ...  <p7>Illumina P7</p7>
                              <cbc>Round1</cbc>                <cbc>Round2</cbc>
                              <cbc>Barcode</cbc>               <cbc>Barcode</cbc>
```

- The sequencing section re-draws the same final library once per read, changing only the
  primer and the arrow, and states template strand plus cycle count in the heading:
  *"(bottom strand as template, cDNA reads, 66 cycles)"*.

## Asset layout

- Page: `methods_html/<Method>.html`
- Stylesheet: `../style_related/page_format.css` (fonts `cmunss.ttf`, `cmunso.ttf`,
  `cmunsx.ttf`, `cmunsi.ttf`, `FiraMono-Regular.ttf`, `FiraMono-Bold.ttf` live beside it)
- Oligo tables: `../data/<Method>/<supp_table>.xlsx`
- Shared figures reused across pages, e.g. `../data/tn5_dimer.svg`

---

## Marking inferred content (project-specific extension)

The upstream scg_lib_structs pages document published chemistries, so they have no convention for
uncertainty. Ours has three genuinely inferred regions (the barcode cassette design, its 4 nt
cohesive overhangs, and the length of the round-D Read-2 arm), so we add one.

**Three layers, all of which must be present:**

1. **A standing caveat in the `<info>` preamble**, listing exactly what is inferred and what is
   documented. A reader who reads only the top of the page must not come away thinking the
   cassette is published chemistry.

2. **A caption on every affected step**, using the `<i>…</i>` convention the upstream pages already
   use for sub-products:

   ```html
   <h3>(5) Pool and split to the Barcode B plate to add the second barcode:</h3>
   <pre>
   <align class="long">
   <i>INFERRED — Atrandi do not publish the barcode cassette design. The 4-nt cohesive
   overhang drawn here is deduced from the read layout (8+4+8+4+8+4+8+1) and from the
   dA-tailed input; the real chemistry may differ.</i>
   ...
   ```

3. **Visual marking of the inferred bases themselves**, so a reader skimming only the diagrams still
   sees it. Add one element to the stylesheet:

   ```css
   inf { border-bottom: 1px dotted #999; }
   ```

   and write inferred placeholder bases in **lowercase** inside it:

   ```html
   <t7>GTGACTGGAGTTCAGACGTGTGCTCTT</t7><inf><t7>ccgatct</t7></inf><cbc>DDDDDDDD</cbc>
   ```

   Lowercase-for-uncertain is a familiar genomics convention (soft-masking) and survives
   copy-paste; the dotted underline carries it visually.

**Rule of thumb:** colour tags (`<p5>`, `<t7>`, `<cbc>`…) say *what a region is*; `<inf>` says
*how well we know it*. The two are orthogonal and nest freely.

**Evidence marking in the notes files** (`01`–`05`): 🟢 verbatim from source · 🟡 derived or
inferred · 🔴 not published. Every 🟡 and 🔴 in the notes must map to one of the three layers above
when it reaches the page.
