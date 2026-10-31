# Substrate-driven reaction engine invariants

31 October2026. Mechanistic review for the selected animation examples. No simulator
or protocol code was changed. These are scoped construction requirements, not an
all-enzyme specificity catalogue or a claim to predict yields, rates or stopping lengths.

## The common contract

An operation receives **actual current molecules**, their end groups, base/sugar states,
strand connectivity, explicit binding relations, an enzyme profile and supplied
reagents. It derives new bonds, nucleotide additions, removed products and retained
molecules. It must not receive an intended product sequence as its source of bases.
An optional expected-product assertion can compare the result afterwards; changing
that assertion must never change the generated reaction.

1. Give each physical strand/base interval stable identity. Copying allocates new
   bases and a new daughter identity, while retaining the primer's original bases.
2. Keep base identity, sugar, modifications and annotations separate. U does not
   imply RNA; methylation is not copied as a methyl group into a new DNA daughter.
3. Keep binding separate from covalent connectivity and screen positions. Moving,
   reframing or folding cannot create a bond or unpair a duplex by itself.
4. A supplied enzyme name alone grants no generic capability. Declare relevant
   substrate polymers, direction, end requirements, activation/cofactor policy,
   displacement ability and selected reaction conditions explicitly.
5. Unknown substrate/end/occupancy values remain unresolved. A selected successful
   example can be constructed under declared conditions; it does not establish
   population conversion or physical removal of other species.

## Copying, gap filling and terminal addition

A DNA-polymerase operation starts at an actual extendable3′OH, with an explicit
primer–template binding map covering its3′ end. Traverse the template3′→5′ and add
complementary DNA5′→3′ from supplied dNTPs. Preserve any unpaired5′ primer tail.
The next base comes from the accessible current template site, never from a planned
output string. An unresolved class may generate an unresolved complementary class,
but cannot generate a secretly chosen A/C/G/T.

Primer mismatch acceptance, proofreading, RNA-template use, LNA-template use and
strand displacement are enzyme/condition-specific policies, not universal defaults.
For a strict exact-copy example, fail unsupported/mismatched inputs instead of
inventing a repaired match. An existing downstream strand must be stopped at,
displaced or processed by an explicitly supplied operation/profile.

A gap is missing nucleotides. A nick is a missing bond between adjacent nucleotides.
Polymerase can fill a declared templated gap from its upstream3′ end but does not
seal a remaining nick. A nontemplated addition is a different operation with a
supplied nucleotide class/policy; it cannot borrow its base from an intended endpoint.
The scRRBS extra3′A therefore remains separate from copying its recessed MspI end.

## Ligation and circularization

A ligase joins actual3′OH and5′phosphate ends, or a supplied preadenylated5′ donor
under the appropriate ligase policy. The ends' original base and sugar identities
are retained; no nucleotide is synthesized. Consumed end groups become an internal
backbone junction, not lingering free-end labels. If activation/leaving-group
chemistry is represented, adenylation/AMP release must be a separate tracked state,
not an adenine quietly inserted into the sequenced backbone.

- **Nick sealing:** require adjacent appropriately oriented termini on the declared
  duplex/splint. A nucleotide gap cannot be sealed merely because its ends are drawn
  close together. Compatible sequence/end geometry does not imply every ligase
  accepts every polymer combination or cofactor condition.
- **MIP:** the existing probe3′ extension arm supplies the primer. Copy only the
  interval between its bound end and the original5′ ligation arm. The captured DNA
  sequence is derived from that target. Only after filling can the new3′OH join the
  original probe5′P. A folded loop before sealing remains a linear molecule.
- **BOLORAMIS:** SplintR joins the DNA probe ends on an **RNA splint**. It does not
  reverse-transcribe the RNA or incorporate that RNA into the circle. The conditional
  RNA match in the current example remains inferred. Joining the last nick creates
  closed DNA topology; neither RNA binding nor a circular drawing establishes it.
- **Small-seq3′ ligation:** RA3 is5′App DNA with a blocked3′ddC. The acceptor is the
  small RNA's3′OH. The truncated Rnl2 route uses preactivation and no added ATP.
  A generic ATP-dependent activation rule must not be substituted.
- **Small-seq5′ ligation:** RA5 is RNA with a blocked5′ end and a free3′OH. Rnl1+ATP
  joins that3′OH to the small RNA5′P. Its RNA UMI is an existing adapter interval,
  not newly synthesized by ligase or polymerase.
- **Gibson:** annealing overlaps creates a paired joint, not a sealed backbone.
  Derive any actual gaps from retained substrate intervals. The repository toy has
  zero gaps, so polymerase must add zero bases there. Seal its remaining compatible
  nicks without inventing repair sequence. The source's3′-resection alternatives
  must not be conflated with the selected5′-resection method.

## Exonuclease resection and adapter degradation

Exonuclease removal begins at an accessible end of a permitted current substrate and
progresses in the declared direction. New removed nucleotide products are distinct
from a physically excised intact strand. Remaining molecule coordinates, end groups
and pairing state must follow the actual removed interval. A closed strand has no
free end until an appropriate nick/cut operation creates one.

- Gibson5′ resection exposes3′ single-strand overlaps. Its exact toy stopping points
  are authored illustrative parameters, not enzyme-predicted lengths. Do not remove
  the surviving3′ overlap merely to match a desired joined product.
- Small-seq first deadenylates **unligated RA3**, producing5′P; RTP binds it to form
  the lambda-exonuclease substrate. Lambda preferentially degrades the phosphorylated
  adapter strand; the source's blocked5′ RTP end protects that separate primer.
  The ligated RNA–RA3 species is a different substrate population, not the same free
  adapter selected for digestion. The existing animation's n−1 endpoint is illustrative,
  not a source-defined one-base residual product.
- Do not take `lib/padlock.py:survives_exonuclease()` as a general enzymology engine.
  Its circular-type predicate does not model enzyme identity, ss/ds specificity,
  end accessibility, nicked circles, modifications or rates. Each digestion requires
  its own source-qualified enzyme/substrate policy.

## Phi29 rolling-circle amplification

Require an actual covalently closed DNA template and an existing correctly annealed
primer with free3′OH. Derive the template coordinate at the primer3′ end; repeatedly
step backward modulo the circle length, copying complementary DNA. Preserve every
primer base once and create each later nucleotide by that traversal. Laps are a
consequence of traversal, not concatenation of a supplied intended repeat string.
A nucleotide/time budget can stop the illustrative simulation; three drawn repeats
are not an experimental yield or native termination mechanism.

As synthesis revisits the circle, previously synthesized DNA is displaced and remains
as an ordered linear daughter. The template circle persists. A drawn loop without a
sealed return bond cannot authorize circular wrapping; a nick must not be hidden by
modulo arithmetic. BOLORAMIS aminoallyl-dUTP incorporation is a supplied modified
nucleotide pool; the current model has not assigned native incorporation positions
or occupancy. Do not manufacture a complete modification pattern or DNA–protein
crosslink map from the known reagent alone.

## T7 transcription

Transcription is de novo RNA synthesis, not extension from a DNA primer and not
replacement of DNA with RNA. Require a source-specified **paired DNA promoter**, its
orientation and transcription-initiation coordinate, T7 polymerase, and supplied
rNTPs. Read the DNA template3′→5′ and produce RNA5′→3′ with U rather than T. Specify
the newly initiated5′ chemical end through the reaction; do not inherit a default
DNA-primer end group. Exact promoter recognition/tolerance and initiation heterogeneity
are not determined by a generic sequence-label match.

Retain both DNA template strands, with a locally changing transcription bubble if
shown. After release, allow another initiation on the same promoter while earlier RNA
molecules remain. Repeated copies do not require DNA replication.

The current CEL-Seq source-owned `native.py:rna_parts()` includes the promoter's
terminal **GGG** in the transcript before downstream adapter/barcode/insert sequence.
Starting after the entire displayed20-base promoter would silently omit these bases.
Use the explicit initiation coordinate, not an assumed boundary after a label called
“promoter”. The initialGs are a sequence-derived/modelled interpretation; the cached
bench protocol independently specifies the T7/rNTP incubation, not a new mapped TSS.

## Bisulfite conversion and later copying

This is chemical treatment of existing DNA, not an enzyme-mediated copy. For an
explicit selected converted site, track unmethylatedC through the sulfonated/deaminated
state and desulfonation to U while preserving **DNA sugar and molecular identity**.
Retain declared protected5mC in the selected example. Unknown methylation, conversion
probability, incomplete conversion, damage and recovery are not filled from a desired
C/T readout. Former genomic partners can cease to be complementary after both are
converted; do not silently redraw them as complements.

Subsequent uracil-compatible DNA synthesis creates A opposite original U. A later
complementary DNA copy can contain T. The U-bearing original remains a distinct
molecule. Likewise, the methyl group on original5mC is not physically copied into a
new conventional DNA daughter. `lib/base_conversion.py:bisulfite_path()` combines
chemical/reporting states; its C→U→T or protected5mC→C path is not permission to
rewrite the original molecule directly to the final sequenced base.

The scRRBS source names PfuTurboCx for the first PCR and Phusion for the second.
Do not silently give every polymerase uracil tolerance or erase original U-bearing
material when a PCR daughter is created.

## Inspected evidence and current limits

The required bundle primer scraper was run for BOLORAMIS, Small-seq, Gibson, CEL-Seq
and scRRBS. Upstream schematic hits were excluded as evidence. Cached originals were
read in these exact locations; no downloads, source edits or new absence claims:

| Original | External text twin and inspected scope |
|---|---|
| [BOLORAMIS, published article](https://doi.org/10.1093/nar/gkab120) | `boloramis__10.1093+nar+gkab120/nar_gkab120_PMC8191787.xml.txt`167–181; direct RNA-dependent ligation, source probe/primer cohybridization, SplintR, Phi29/dNTP/aminoallyl-dUTP |
| [Small-seq detailed protocol2018](https://doi.org/10.1038/s41596-018-0049-y) | `small-seq__10.1038+nbt.3701/hagemann_s41596-018-0049-y.pdf.txt`327–364;3′ ligation, adapter deadenylation/duplex digestion,5′ ligation/UMI. This is explicitly the detailed2018 version. |
| [Gibson2009 published supplementary text](https://doi.org/10.1038/nmeth.1318) | `gibson-assembly__10.1038+nmeth.1318/nmeth.1318_41592_2009_BFnmeth1318_MOESM319_ESM.pdf.txt`815–925; source3′-resection alternatives, gap-fill/nick-seal distinctions, circularity/RCA and comparison withT5. Exact selected toy bases/stops are repository inventions. |
| [CEL-Seq2012 detailed protocol](https://ars.els-cdn.com/content/image/1-s2.0-S2211124712002288-mmc4.pdf) | `cel-seq__10.1016+j.celrep.2012.08.003/j.celrep.2012.08.003_mmc4.pdf.txt`172–214; pooled-cDNA recovery and A/G/C/U+T7 buffer/enzyme37°C13h IVT |
| [scRRBS published article](https://doi.org/10.1101/gr.161679.113) | `scrrbs__10.1101+gr.161679.113/gr.161679.113_PMC3847781.xml.txt`105; MspI, fill/A-tail, premethylated adapters, bisulfite/desulfonation, first PfuTurboCx and later Phusion |

Construction review also inspected `lib/circular.py`, `lib/padlock.py`,
`lib/base_conversion.py`, `lib/assembly_geometry.py` and the selected protocol modules.
Those code contracts identify current modelling choices; they do not establish new
source chemistry or general enzyme specificity. Detailed catalyst structures, exact
rates, mismatch spectra, occupancy, native stopping distributions and whole-population
fates remain outside this bounded engine proposal.
