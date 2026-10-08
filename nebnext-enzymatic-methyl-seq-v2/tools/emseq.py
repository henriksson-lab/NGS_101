"""NEBNext Enzymatic Methyl-seq v2, E8015 manual v2.0 6/25."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from chemdraw import Construct, Row, Segment
from endprep import dA_tailed_scene, repair_and_dA_tail
from nebnext import STANDARD_PRIMERS, indexed_library
TITLE="NEBNext Enzymatic Methyl-seq v2"
NOTES="01_nebnext-em-seq-v2.html"
SOURCE='Commercial protocol: <a href="https://www.neb.com/en-us/products/e8015-nebnext-enzymatic-methyl-seq-v2-kit">NEB E8015 manual v2.0 (June 2025)</a>.'
SUMMARY="Adapters are ligated before enzymatic conversion. TET2 and T4-BGT protect 5mC/5hmC; APOBEC converts unmodified C to U; Q5U then amplifies the converted library."
OUTCOME={"C":"U","5mC":"C","5hmC":"C"}
DNA=Construct([Segment("fragmented DNA", "X"*42, placeholder=True)], name="fragmented DNA")
END_PREP=repair_and_dA_tail(DNA)
PROTECTION=[Row(chunks=[("C → C     5mC → 5caC     5hmC → 5gmC",None,False)])]
DEAMINATION=[Row(chunks=[("unmodified C → U     protected 5caC / 5gmC → C in sequence",None,False)])]
STEPS=(("End-repair, dA-tail and ligate EM-seq adaptors", dA_tailed_scene(END_PREP), "The conversion-resistant EM-seq adaptor is attached before base conversion."),
       ("Protect modified cytosines", PROTECTION, "TET2 oxidizes 5mC through 5caC; T4-BGT glucosylates 5hmC. Both products resist APOBEC."),
       ("Deaminate unmodified cytosine", DEAMINATION, "APOBEC changes unprotected C to U while protected derivatives remain readable as C."),
       ("Amplify converted strands with Q5U", Construct([Segment("converted insert: T where C was unmethylated", "X"*46, placeholder=True)], name="converted library"), "Q5U copies uracil-containing templates; dual-index primers complete the flow-cell arms."))
SEQ_PRIMERS=STANDARD_PRIMERS; REQUIRED_ROLES=tuple(p.role for p in SEQ_PRIMERS)
READOUT="Paired reads sequence the converted insert; C calls represent protected 5mC/5hmC and T calls at genomic C positions represent unmodified cytosine. i5 and i7 are separate index reads."
FINAL_CAPTION="PCR-completed, dual-index EM-seq library. The manual prints TruSeq-compatible trimming sequences, so primer binding is computed from the canonical arms."
def final_library(): return indexed_library(Construct([Segment("enzymatically converted insert", "X"*42, placeholder=True)]), "EM-seq v2 library")
