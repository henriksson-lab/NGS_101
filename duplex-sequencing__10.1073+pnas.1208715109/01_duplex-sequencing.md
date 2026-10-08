# Duplex Sequencing

## 1. What it is

Duplex Sequencing attaches complementary random tags to both strands of each original DNA duplex so PCR and sequencing errors can be rejected by agreement between independently copied strand families.

## 2. Sources and evidence

- 🟢 Schmitt MW et al. *Detection of ultra-rare mutations by next-generation sequencing.* PNAS 109:14508–14513 (2012). [doi:10.1073/pnas.1208715109](https://doi.org/10.1073/pnas.1208715109).
- 🟢 A randomized adapter strand is copied to produce a complementary double-stranded tag, then A-tailed adapters are ligated to T-tailed sheared DNA.
- 🟢 Each end contributes 12 random nucleotides. A strand family is identified by the 24-nt αβ pair; its complementary partner carries the reciprocal βα relationship.
- 🟢 PCR descendants yield single-strand consensus sequences before complementary SSCS pairs yield a duplex consensus sequence.

## 3. Molecular path

The essential object is not an ordinary UMI family but a pair of related UMI families. The tag relationship retains which two sequenced strands originated as the two halves of one physical duplex.
