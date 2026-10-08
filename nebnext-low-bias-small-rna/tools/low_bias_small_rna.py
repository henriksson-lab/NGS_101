"""NEBNext Low-bias Small RNA library preparation, E3420 manual v2.0 9/25."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from chemdraw import Construct, Row, Segment
from nebnext import STANDARD_PRIMERS, indexed_library
TITLE="NEBNext Low-bias Small RNA"
NOTES="01_nebnext-low-bias-small-rna.html"
SOURCE='Commercial protocol: <a href="https://www.neb.com/en-us/products/e3420-nebnext-low-bias-small-rna-library-prep-kit">NEB E3420 manual v2.0 (September 2025)</a>.'
SUMMARY="Small RNAs with 5′ phosphate and 3′ hydroxyl are captured through randomized splint ligation, followed by selective adaptor removal, reverse transcription and UDI PCR."
CAVEAT="The current manual does not disclose adaptor bases. The final panel uses the documented TruSeq-compatible primer-site boundary while rendering every undisclosed base as a placeholder."
INPUT=Construct([Segment("5′-phosphorylated small RNA", "R"*28, placeholder=True)], name="small RNA")
THREE=Construct([Segment("small RNA", "R"*28, placeholder=True), Segment("3′ adaptor", "X"*14, "r2", placeholder=True, inferred=True)], name="3′-ligated RNA")
SPLINT=[Row(chunks=[("5′ adaptor ** small RNA ** 3′ adaptor",None,False)]), Row(chunks=[("       randomized splint pairs across the ligation substrate", "me",True)])]
CDNA=Construct([Segment("5′ adaptor remnant", "X"*12, "r1", placeholder=True, inferred=True), Segment("small-RNA cDNA", "X"*28, placeholder=True), Segment("3′ adaptor remnant", "X"*12, "r2", placeholder=True, inferred=True)], name="small-RNA cDNA")
STEPS=(("Ligate the 3′ adaptor", THREE, "The acceptor is an RNA 3′ hydroxyl. Excess free adaptor is enzymatically removed before the next ligation."),
       ("Randomized-splint 5′ adaptor ligation", SPLINT, "The splint reduces sequence-dependent ligation bias. Exact adaptor and splint bases are not disclosed in the manual."),
       ("Modify the 3′ adaptor and reverse-transcribe", CDNA, "During 5′ ligation the existing 3′ adaptor is enzymatically converted into an RT-priming substrate; RT copies the captured RNA."),
       ("Size-select and UDI-amplify", Construct([Segment("selected small-RNA cDNA", "X"*40, placeholder=True)], name="selected cDNA"), "Bead selection can emphasize miRNAs or retain a broader small-RNA range before dual-index PCR."))
SEQ_PRIMERS=STANDARD_PRIMERS; REQUIRED_ROLES=tuple(p.role for p in SEQ_PRIMERS)
READOUT="A 56-base Read 1 is normally sufficient for the small-RNA insert; i5 and i7 are read separately from the NEBNext LV unique-dual-index library."
FINAL_CAPTION="INFERRED — TruSeq-compatible UDI boundary. Proprietary ligated-adaptor remnants remain dotted; the standard sequencing-primer sites are computed from the completed PCR product."
def final_library():
    ins=Construct([Segment("5′ kit remnant", "X"*8, placeholder=True, inferred=True), Segment("small-RNA insert", "X"*28, placeholder=True), Segment("3′ kit remnant", "X"*8, placeholder=True, inferred=True), Segment("3′ adaptor boundary A", "A", inferred=True, note="required by the documented TruSeq-compatible Read 2 site; exact proprietary adaptor sequence is unavailable")])
    return indexed_library(ins, "NEBNext low-bias small-RNA library",
                           inferred_adaptors=True, dA_junction=False)
