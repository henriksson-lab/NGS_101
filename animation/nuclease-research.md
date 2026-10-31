# Nuclease substrates and product conservation

Focused source check, 31 October 2026. No simulator code changes. The CRISPR-MIP
bundle scraper was run; the v1 Methods treatment was read directly in the cached
original. Manufacturer pages/guides below were inspected for the named products;
they are not proof of complete digestion of every molecular geometry.

## Source CRISPR-MIP treatment

[Selinger et al., doi:10.1101/2024.03.28.587082v1](https://doi.org/10.1101/2024.03.28.587082)
Methods “CRISPR-MIP sgRNA abundance quantification”, Steps1–4, cached
`2024.03.28.587082_v1.full.pdf.txt` lines471–493: following extension/ligation, add
10U ExonucleaseI **EN0581** and50U ExonucleaseIII **M0206L**,37°C45min, followed by
80°C20min inactivation. The stated purposes are removing non-circularized MIPs and
genomic DNA. That is not a demonstrated complete-recovery/removal claim for all
hybridized, nicked, blocked, damaged or partially extended inputs.

## Exact product profiles

| Product and original source | Inspected scope | Relevant behavior |
|---|---|---|
| [Thermo ExonucleaseI, EN0581, MAN0012007 RevC00](https://documents.thermofisher.com/TFS-Assets/LSG/manuals/MAN0012007_Exonuclease_I_Exo_I_UG.pdf), printed28May2020 | Page1 Description, Applications, Note | Stepwise3′→5′ digestion of ssDNA with free3′OH, releasing5′dNMPs. The5′terminal dinucleotide remains. Phosphoryl- or acetyl-blocked3′OH prevents this activity. The guide says it is unsuitable for removing dsDNA3′overhangs. |
| [NEB ExonucleaseIII, M0206](https://www.neb.com/en-gb/products/m0206-exonuclease-iii-e-coli) | Product Information, substrate summary, Product Notes | Preferred entry is blunt/recessed duplex3′OH; duplex nicks also initiate progressive3′→5′ deletion. Four-base-or-longer3′protrusions are essentially resistant. Phosphorothioate bonds resist cleavage. Activity depends on sequence, helix, conditions and enzyme/substrate ratio. |
| [NEB RNaseH, M0297](https://www.neb.com/en-us/products/m0297-rnase-h) and [cut chemistry FAQ](https://www.neb.com/en-gb/faqs/which-side-of-the-rna-base-does-rnase-h-cut) | Description/Product Information; FAQ answer | Endoribonuclease cleavage of RNA in RNA–DNA hybrids; DNA is retained. Products carry upstream3′OH and downstream5′phosphate. This is not a requirement to enter from a free RNA terminus. |

M0206 is also reported to possess AP-endonuclease,3′phosphatase and RNaseH
activities. Do not silently activate these from the exonuclease name: they require
explicitly represented substrates and activity scope. In particular,3′phosphate
is not a universal protection claim for every activity of ExoIII. The paper's
actual80°C inactivation remains distinct from the supplier's standard recommendation.
RNaseH isoforms are also distinct; do not extend the M0297 profile to an HII/H2
single-ribonucleotide cleavage rule or the thermostable enzyme's temperature range.

## Simulator and animation invariants

1. **Travel along the attacked backbone.** For3′→5′ exonuclease, select an eligible
   actual3′terminus and repeatedly hydrolyze the bond connecting the last nucleotide
   to its predecessor. The enzyme moves toward the original5′end. A top strand
   drawn5′left→3′right is therefore chewed right-to-left; its antiparallel partner
   is chewed left-to-right. No complementary synthesis or reverse polymerization
   occurs, and no dNTP supply is consumed by this hydrolysis.
2. **Retain products and chemical identities.** Released nucleotides keep their
   original IDs/base/sugar and become explicit5′monophosphate products. The
   remaining strand exposes a new3′OH. EN0581's terminal dinucleotide must stay
   linked in its single-activity model. Endonuclease cleavage instead creates two
   oligonucleotide components unless repeated cleavage actually produces monomers.
3. **Check access again after each reaction.** ExoI's substrate profile is ssDNA;
   a hybridized3′arm is not interchangeable with free linear probe. ExoIII's
   preferred duplex3′end can exist on a partially annealed probe or at a nick.
   A strand may expose a different substrate as its partner is digested. Hairpins,
   blocked ends and long3′overhangs need explicit handling, not unconditional
   “linear=true” selection. Low/resistant activity outside the modeled profile is
   not experimentally proven zero activity.
4. **A ring label does not confer immunity.** An intact covalently closed DNA
   backbone has no exonuclease entry end in this selected exo-only model. A nicked
   duplex circle supplies an ExoIII entry. An endonuclease or appropriately modeled
   damage-processing activity can introduce an entry into a previously closed
   molecule. The successfully sealed padlock's survival must follow its actual
   connectivity and permitted activities, not its name.
5. **Keep endonuclease cleavage separate from hybrid melting.** RNaseH attacks
   eligible internal RNA bonds in a hybrid and leaves the DNA backbone intact.
   Cleaved RNA fragments can initially remain base-paired; do not erase every
   duplex bond or destroy the RNA IDs merely because a covalent bond was cut.
   Subsequent fragment dissociation is an explicit physical step. Short-hybrid
   footprint or exact site choices need the selected enzyme's scope; this pass
   does not establish cleavage of a lone ribonucleotide in DNA.

For the cleanup movie, show actual sealed padlock, unligated probe and a represented
linear target fraction together. Digest only eligible ends, keeping both monomers
and residual oligonucleotides visible before they leave the viewport as a chosen
presentation fraction. Unassembled native target ends or alternative access states
remain unresolved; the source purpose does not authorize drawing all material gone.
