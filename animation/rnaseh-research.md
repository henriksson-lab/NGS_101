# RNase H and CEL-Seq second-strand synthesis

31 October2026. Bounded reading of three originals: the CEL-Seq2012 detailed protocol,
the contemporaneous MessageAmpII manual, and a primary E.coli RNaseHI study. No core
or animation code was changed. An unresolved enzyme-vial identity is not evidence
that the information is unpublished.

## Evidence inspected

| Original | Identity and exact inspected location |
|---|---|
| [CEL-Seq2012 detailed protocol](https://ars.els-cdn.com/content/image/1-s2.0-S2211124712002288-mmc4.pdf) | Cached external `cel-seq__10.1016+j.celrep.2012.08.003/j.celrep.2012.08.003_mmc4.pdf.txt`90–162;16July2012 protocol, actual primer/RT/second-strand recipe. Mandatory bundle primer scraper rerun. |
| [MessageAmpII aRNA Amplification Kit,1751M RevE](https://documents.thermofisher.com/TFS-Assets/LSG/manuals/1751ME.pdf) | Vendor original, printed02/2011, AM1751; printedp3 Procedure Overview, p4 supplied reagents, pp14–15 sectionII.D first-round second-strand mix. Optional second-round sectionIV was also inspected to distinguish its separate RNaseH treatment. Read through web tool; no local source copy is claimed. |
| [Lee et al., RNase H is an exo- and endoribonuclease…](https://doi.org/10.1093/nar/gkab1064) | Primary paper, Nucleic Acids Research50(2022):1801–1814, online12November2021; Introduction, Results “Binding affinities…”/Figure2 and Discussion. Read through publisher's accessible article text. |

Primary and vendor findings below are paraphrased within their inspected scope. The
cached CEL-Seq original remains under the external source root; no source text/PDF
or rendered source page was added to the repository.

## What is established for original CEL-Seq

🟢 The source adds its barcoded oligo(dT) primer and ERCC RNA before ArrayScript RT.
The substrate entering second-strand synthesis is therefore first-strand DNA made
from RNA, with retained RNA lineage and the source primer's5′ handle. It is not an
already fabricated dsDNA molecule.

🟢 Its second-strand table adds water6.3µl, buffer1µl, dNTP0.4µl, **DNA Pol0.2µl**
and **RNaseH0.1µl**, followed by16°C2h. It does not separately list DNA ligase, and
it does not identify its DNA Pol vial as E.coli polymeraseI. This is a positive
transcription of the inspected table, not a conclusion that other documents lack
additional component information.

🟢 The2011 MessageAmpII manual independently describes coupled RNA degradation and
second-strand synthesis. Its100µl first-round mix has the corresponding10-fold reagent
amounts. Its component list names DNA Polymerase and RNaseH as separate supplied
vials. Neither enzyme identity nor a separately added ligase can be inferred merely
from similarity to another second-strand recipe.

The manual's optional **second amplification round** instead adds RNaseH after reverse
transcription before adding fresh T7 oligo(dT) primer. That distinct procedure must
not be substituted for the source's first-round coupled reaction.

## RNaseHI/H1 substrate rule

🟢 Lee et al. study E.coli RNaseHI, also termed H1-class RNaseH. Cleavage recognition
requires four consecutive RNA nucleotides, two on either side of the cleaved bond,
in an RNA:DNA hybrid. Their measured product binding supports that requirement.
The paper distinguishes HI/H1 from HII/H2 and finds that substrate-junction geometry
can change processivity; a3′DNA-overhang substrate does not behave like every blunt
hybrid. This supports an explicit selected endonucleolytic mode, not a universal claim
that RNaseHI cleavage is always uniformly random or has no directionality.

For a strict H1-class demonstration, the engine should therefore derive eligible
bonds from the actual bound substrate:

- The scissile strand must carry RNA sugars across a contiguous four-nucleotide
  neighborhood, with two nucleotides on each side of the candidate bond.
- Those nucleotides must map consecutively and antiparallel to the same DNA partner,
  with an intact local hybrid footprint. A collection of unrelated pairing marks is
  not that substrate. Merely placing
  an RNA label next to DNA does not create the substrate.
- Free RNA, RNA:RNA and DNA:DNA are not this substrate. An isolated ribonucleotide
  embedded in DNA must not pass the H1 four-RNA rule.
- Break the RNA phosphodiester bond. Retain the DNA strand and its bases unchanged.
  A cleavage event creates separate RNA fragments with persistent base ownership;
  dissociation, subsequent cleavage and physical removal are later events.
- Recompute eligibility after every cut/dissociation. A four-RNA recognition patch
  is not an instruction to divide the entire RNA into fixed four-base chunks.
- Choose a site only among eligible bonds. A deterministic seed/selection policy can
  make the illustrative run reproducible, but does not establish native cleavage
  coordinates, equal probabilities, rate or quantitative completion.

The usual hydrolytic fragment-end model is new3′OH/new5′phosphate. This is a mechanistic
end-state convention, not evidence for a particular source fragment length or recovered
primer. Cleaved RNA remains RNA; it does not become DNA or vanish as a side effect of
cutting. Short fragments may remain transiently bound until an explicit dissociation
policy releases them.

HII/H2 needs a **separate enzyme profile**. The present three-source reading does not
parameterize its single-ribonucleotide-excision specificity or kinetics. It would be
incorrect either to assign H2 the H1 four-RNA minimum, or to allow the H1 profile to
remove an isolated embedded ribonucleotide because a differently named RNaseH can.
The MessageAmp vial's exact enzyme identity also remains a distinct source question;
an E.coliHI demonstration is not proof that every commercial RNaseH vial is identical.

## RNA fragments as primers: what the animation can and cannot infer

🟡 A retained, bound RNA fragment with an accessible3′OH is a potential primer for a
polymerase profile that explicitly accepts RNA primers on DNA. Extension adds DNA to
that **existing RNA3′ end**; the product is initially an RNA–DNA chimera. It is not a
new all-DNA strand born with the intended full product sequence. The primer's original
RNA sugars/base identities remain until a later permitted removal/replacement event.

A mechanistic DNA-polymeraseI/nick-translation example would additionally require
explicit RNA-primer extension and downstream5′-exonuclease/processing capabilities.
A continuing nick would require an actual compatible ligation operation if sealing is
shown. Those are meaningful generic operations, but the inspected CEL-Seq table and
MessageAmp manual do not establish that exact commercial-vial implementation. Do not
supply an invented CEL-Seq ligase, fixed RNA-priming site or complete nick-translation
route just to arrive at a desired fully sealed DNA product.

A source-qualified sequence of views can show:

1. The retained RNA–first-strand-DNA hybrid from actual RT.
2. Eligible RNA cleavage by a clearly named H1-class mechanism illustration.
3. Persistent RNA fragments and unchanged DNA; explicitly dissociate only selected
   fragments, preserving any still-bound candidate primers.
4. If the polymerase capability is declared/established, new DNA growing from a
   retained RNA3′OH; show sugar boundaries and any unresolved nicks/remnant ends.
5. The source's observed paired-promoter/dsDNA endpoint as a separate scoped product
   boundary if detailed kit priming/repair has not been mechanistically reconstructed.

That final boundary must not masquerade as a fully substrate-generated base-exact
reaction trajectory. The source establishes successful second-strand preparation;
it does not by itself specify an individual RNA cleavage/priming history. A requested
chemistry test should expose this limitation rather than pass a desired product into
the enzyme operation and animate toward it.

## Concrete implementation checks

- A camera move or oligo displacement cannot change chemical bond state.
- H1 cuts zero bonds on an unhybridized RNA, DNA:DNA duplex or three-RNA patch.
- A paired four-RNA patch supplies a candidate central bond, not a prescribed cut time.
- Total RNA/DNA base ownership is conserved across cleavage; only RNA connectivity
  changes. The DNA sequence cannot change during RNaseH treatment.
- A released fragment no longer supplies a DNA-hybrid recognition patch or a bound
  primer. A retained RNA primer cannot be silently recast as an all-DNA primer.
- RNA-primer extension is rejected unless the selected polymerase profile supports it.
- Sealing is rejected when no ligase is supplied or a nucleotide gap remains.
- Source/native and generic illustrative enzyme identities stay distinct in metadata.

These are scoped mechanistic invariants, not new source facts about CEL-Seq product
stoichiometry, native end repair, mismatch tolerance or all RNaseH family members.
