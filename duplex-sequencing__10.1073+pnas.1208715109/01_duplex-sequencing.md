# Duplex Sequencing

## 1. What it is

Duplex Sequencing attaches complementary random tags to both strands of each original DNA duplex so PCR and sequencing errors can be rejected by agreement between independently copied strand families.

## 2. Sources and evidence

- 🟢 Schmitt MW et al. *Detection of ultra-rare mutations by next-generation sequencing.* PNAS 109:14508–14513 (2012). [doi:10.1073/pnas.1208715109](https://doi.org/10.1073/pnas.1208715109).
- 🟢 A randomized adapter strand is copied to produce a complementary double-stranded tag, then A-tailed adapters are ligated to T-tailed sheared DNA.
- 🟢 Each end contributes 12 random nucleotides. A strand family is identified by the 24-nt αβ pair; its complementary partner carries the reciprocal βα relationship.
- 🟢 PCR descendants yield single-strand consensus sequences before complementary SSCS pairs yield a duplex consensus sequence.
- 🟡 The authors' later detailed protocol, [Kennedy et al., *Nature Protocols* 9:2586–2606 (2014)](https://doi.org/10.1038/nprot.2014.170), prints adapter oligos MWS51/MWS55 and PCR primers MWS13/MWS20/MWS21 in Table 3. It explicitly states that these **redesigned** adapters terminate in a 3′ dT and pair with dA-tailed fragments, whereas the defining 2012 paper used 3′-dA adapters with dT-tailed DNA. Their published Table 3 bases therefore cannot be assigned to the 2012 adapter in this page.
- 🔴 The 2012 paper's SI Materials and Methods PDF (`1208715109_pnas.201208715SI.pdf`) is listed by PMC but its download is gated; the PNAS supplement URL also returned 403 in this audit. The exact 2012 adapter-strand and PCR-primer oligos remain unverified, despite the disclosed 12+12-nt tag geometry.

## 3. Molecular path

The essential object is not an ordinary UMI family but a pair of related UMI families. The tag relationship retains which two sequenced strands originated as the two halves of one physical duplex.
