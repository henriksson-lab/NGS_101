# SLAM-seq

## 1. What it is

SLAM-seq labels newly synthesized RNA with 4-thiouridine and converts the label into a sequencing-detectable base substitution through thiol alkylation.

## 2. Authoritative source

🟢 Herzog et al., *Thiol-linked alkylation of RNA to assess expression dynamics*, Nature Methods (2017), DOI `10.1038/nmeth.4435`.

## 3. Construct-changing operations

🟢 Iodoacetamide alkylates 4-thiouridine. Reverse transcription across the modified base generates characteristic T-to-C mismatches; the downstream RNA library method is otherwise conventional.

## 4. Chemistry and library endpoint

🟢 The paper reports standard alkylation as 10 mM iodoacetamide in 50% DMSO and 50 mM sodium phosphate, pH 8, for 15 min at 50 °C, followed by DTT quenching. IAA attaches a carboxyamidomethyl group to s4U.

🟢 Reverse transcription of alkylated s4U produces the diagnostic T-to-C readout. The defining mRNA experiment used the Lexogen QuantSeq 3′ mRNA library kit and Illumina single-read sequencing.

🟡 The paper names the commercial kit but does not print its adapter oligos. The schematic therefore marks the Illumina-compatible single-index adapter skeleton as inferred rather than presenting kit bases as published facts.
