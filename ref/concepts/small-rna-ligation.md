# Sequential adapter ligation, and masking oligos

How a 22-nt RNA with no poly(A) tail, no cap and no handle gets two defined ends — with no
polymerase involved. The front end of every small-RNA method (Small-seq, NEBNext Small RNA,
QIAseq miRNA, TruSeq Small RNA, the CLIP family), and the only alternative in this repo to
template switching and tagmentation. Implemented in `small-seq__10.1038+nbt.3701/tools/smallseq.py`.

## The problem

Tn5 needs ~100 bp of duplex to tagment; a miRNA is 22 nt of single-stranded RNA. Oligo-dT
needs a poly(A) tail; small RNAs have none. Template switching needs a reverse transcriptase
to run off the 5' end of a template, which gives you *one* handle and only after RT has
already been primed from somewhere. So: ligate both ends, before anything else happens.

Two ligations, and they cannot be done at once, because an adapter in the reaction is both
a donor and an acceptor. Left to itself the mix produces **adapter dimers** — 5' adapter
ligated straight to 3' adapter, no insert — which are shorter than the real product, amplify
better, and will dominate the run.

## The four tricks

**1. A pre-adenylated 3' adapter, and a ligase that cannot activate anything.**
The 3' adapter is supplied already carrying its own activation: `5'-rApp…`, chemically
adenylated. It is ligated by **T4 RNA ligase 2, truncated** in a reaction with **no ATP**.
With no ATP the ligase cannot adenylate anything else, so the only possible donor in the
tube is the adapter — small RNAs cannot be joined to each other, and the adapter cannot be
joined to itself. The adapter's 3' end is blocked (`ddC`) for the same reason.

> Order chemically adenylated adapter. Enzymatic adenylation after synthesis is incomplete,
> and the unadenylated fraction is dead weight that competes for the acceptor.

**2. Destroy the leftover 3' adapter before adding the 5' adapter.**
This is the step that separates a working small-RNA protocol from a tube of dimers:

- a **5' deadenylase** strips `rApp` from unligated adapter, leaving a plain 5'-phosphate;
- the **RT primer** — which is the exact reverse complement of the 3' adapter — anneals to
  it, making a blunt duplex;
- **lambda exonuclease** digests the 5'-phosphorylated strand of a duplex, so it eats the
  adapter. The RT primer survives because its own 5' end is blocked (biotin).

Adapter already ligated to a small RNA is protected: its 5'-phosphate is now an internal
phosphodiester. So the reaction removes exactly the failure mode and nothing else. The RT
primer is left in the tube and reused twice more — to prime RT, and as the reverse primer
of the first PCR.

**3. A 5' adapter that can only ligate onto a 5'-phosphate.**
The 5' adapter is RNA, blocked at its 5' end (an amino-linker), with a free 3'-OH, ligated
by **T4 RNA ligase 1 with ATP**. Blocking its 5' end stops 5'-adapter concatemers. Requiring
a 5'-phosphate on the acceptor is a selectivity feature, not an accident: **capped mRNA
cannot be ligated**, so the method excludes mRNA by chemistry rather than by size.

**4. The UMI goes in the 5' adapter.**
The 5' adapter is attached before any amplification, so a random stretch inside it labels an
individual molecule. Because it sits at the very start of the read, every read carries a
count. Compare the two other places a UMI can live:

| Where | Method | Labels |
|---|---|---|
| In the ligated 5' adapter | Small-seq | every molecule that got an adapter; read first |
| In a template-switching oligo | SMART-seq3 | only molecules that template-switched, i.e. genuine 5' ends |
| On a capture probe | CRISPR-MIP | every molecule the probe captured |

## Why the dimer never fully goes away

Run the arithmetic: an adapter dimer is the full library construct with a zero-length insert.
For Small-seq that is **130 bp** against **148 bp** for the shortest insert anyone analyses
(18 nt). An 18-bp difference at 130 bp is beyond a 3% agarose cassette, so size selection
*reduces* dimers and never removes them. This is a structural property of ligation-based
library prep, and it is why every such protocol lists adapter dimers as a limitation.

## Masking oligos — depletion by doing nothing

Ligation is indiscriminate: the most abundant 5'-phosphorylated, 3'-OH RNA in the cell wins,
and in a mammalian cell that is a ribosomal RNA fragment. The usual answer is to pull the
rRNA out (biotinylated capture probes, RNase H digestion). A **masking oligonucleotide**
does something cheaper: a DNA oligo complementary to the rRNA's **3' end**, added to the
lysis buffer, which anneals during the lysis incubation and leaves the rRNA's 3' end inside
a duplex. The 3' adapter cannot be ligated there; with no 3' adapter there is no RT primer
site; so the rRNA never enters the library. Nothing is removed from the tube.

Design points, all of which Small-seq's 5.8S mask illustrates:

- **Target the end the adapter needs**, not an internal stretch.
- **Make it long and GC-rich enough to stay annealed through the lysis incubation.** The
  Small-seq mask is 76 nt with a Tm of ~83 °C, against a 72 °C lysis step.
- **Block its own 3' end** (3'-biotin here) so it cannot be extended by the reverse
  transcriptase and become a product in its own right.
- **It is species-specific.** The Small-seq mask is designed for human and mouse; another
  organism needs its own.
- **Leaving it out is a useful positive control.** The masked species reappears as a sharp,
  predictable band — 270–290 bp for 5.8S in Small-seq — which proves a cell was there.

The idea generalises to any ligation protocol with one dominant unwanted species, and it is
cheap: one more oligo in the lysis buffer, no extra handling step, no material lost.
