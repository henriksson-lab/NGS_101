# Padlock probes and circularisation

How to capture a specific sequence, label it, and throw away everything else — without a
gel, a bead cleanup or a size selection. Implemented in `lib/padlock.py`.

## The mechanism

A **padlock probe** is a single linear ssDNA whose two *ends* are complementary to two
places on a target, a short distance apart. Both ends anneal, leaving the probe arched over
a gap. A polymerase fills the gap from the extension arm's 3'-OH; a ligase joins the new
3' end to the probe's 5' phosphate. The probe is now a **covalently closed circle** that
threads the target.

```
            5'-p [ligation arm]---[ backbone: UMI, primers, index ]---[extension arm] 3'-OH
                      |||                                                   |||
 target  ----------[ anneal ]--------------- gap ----------------------[ anneal ]----------
                                     polymerase fills ->
                                                         <- ligase closes here
```

Used for quantification it is called a **molecular inversion probe (MIP)** — "inversion"
because the captured target ends up flanked by probe backbone, so the probe's UMI and index
travel with it.

## Why it is worth the trouble

- **The circle is the selection.** Exonucleases I and III degrade anything with a free end.
  Unreacted probe, genomic DNA, and probes that annealed but failed to ligate are all
  linear and all destroyed; only correctly closed circles survive. The failure modes are
  removed by enzymology rather than by a cleanup step.
- **Two arms, not one primer.** A PCR primer needs one match; a padlock needs two, in the
  right order, on the same strand, a fixed distance apart. That is a far stricter condition
  — in practice a single site in a whole plasmid.
- **The UMI is attached before amplification.** So the readout counts *molecules*, not
  reads, which is what makes it immune to PCR efficiency differences and sequencing depth.

## Design rules

**Make the ligation arm hotter than the extension arm.** The polymerase travels toward the
ligation arm; if that arm is the weaker duplex, the polymerase displaces it and the circle
never closes. A few °C of Tm margin is the usual fix.

**Keep the gap short.** It must be crossed by the polymerase within the reaction, and
everything in it is synthesised rather than probed — so sequence you need to read goes in
the gap, and sequence you need only for binding goes in the arms.

**Put the amplification primer sites inside the captured region.** Then a probe that was
annealed but not extended cannot amplify, adding a fourth filter after the three above.

## Where it shows up

`crispr-mip__10.1101+2024.03.28.587082/` (sgRNA quantification), targeted resequencing panels, *in situ* sequencing and
rolling-circle readouts, and SNP genotyping — the original MIP application.

## References

- Nilsson M *et al.* *Padlock probes: circularizing oligonucleotides for localized DNA detection*. Science 1994;265(5181):2085–8. doi:[10.1126/science.7522346](https://doi.org/10.1126/science.7522346)
- Hardenbol P *et al.* *Multiplexed genotyping with sequence-tagged molecular inversion probes*. Nat Biotechnol 2003;21(6):673–8. doi:[10.1038/nbt821](https://doi.org/10.1038/nbt821)
- Selinger M *et al.* *CRISPR-MIP…* bioRxiv 2024.03.28.587082. doi:[10.1101/2024.03.28.587082](https://doi.org/10.1101/2024.03.28.587082)
