"""NEBNext Ultra II Directional RNA library preparation, E7760."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from chemdraw import Construct, Row, Segment
from nebnext import STANDARD_PRIMERS, indexed_library
TITLE="NEBNext Ultra II Directional RNA"
NOTES="01_nebnext-ultra-ii-directional-rna.html"
SOURCE='Commercial protocol: <a href="https://www.neb.com/en-us/products/next-generation-sequencing-library-preparation/library-preparation-for-illumina/rna-library-prep-for-illumina">NEB E7760 product manual and official RNA workflow</a>.'
SUMMARY="Fragmented RNA is copied into double-stranded cDNA with dUTP in strand two. USER selectively removes that marked strand before PCR, preserving transcript orientation."
RNA=Construct([Segment("fragmented RNA", "R"*42, placeholder=True)], name="RNA fragment")
FIRST=Construct([Segment("first-strand cDNA", "X"*42, placeholder=True)], name="RNA:cDNA product")
SECOND=[Row(chunks=[("first-strand cDNA  3′-xxxxxxxxxxxxxxxx-5′",None,False)]), Row(chunks=[("second strand      5′-xxxUxxxxUxxxxxxx-3′", "w1",False)])]
USER=[Row(chunks=[("USER cuts the dUTP-labelled second strand → it cannot seed library PCR",None,False)]), Row(chunks=[("the unmarked first strand defines the sequenced orientation",None,False)])]
STEPS=(("Fragment RNA and synthesize the first cDNA strand", FIRST, "Random priming copies enriched or rRNA-depleted RNA fragments."),
       ("Make a dUTP-marked second strand", SECOND, "dUTP substitutes for dT during second-strand synthesis; the mark belongs only to strand two."),
       ("End-prep, ligate adaptor and remove strand two", USER, "USER excises the dU-containing strand. Strand-specificity reagent reinforces the same directional selection."),
       ("PCR-complete the surviving strand", Construct([Segment("orientation-preserving cDNA insert", "X"*42, placeholder=True)], name="directional insert"), "NEBNext dual-index primers amplify molecules derived from the surviving first strand."))
SEQ_PRIMERS=STANDARD_PRIMERS; REQUIRED_ROLES=tuple(p.role for p in SEQ_PRIMERS)
READOUT="Read 1 and Read 2 traverse opposite ends of the cDNA insert; their orientation is interpretable because the dUTP-labelled second strand was removed before amplification."
FINAL_CAPTION="Dual-index directional RNA library. The molecule has ordinary TruSeq-compatible sequencing sites; directionality was imposed earlier by selective second-strand destruction."
def final_library(): return indexed_library(Construct([Segment("directional cDNA insert", "X"*42, placeholder=True)]), "Ultra II Directional RNA library")
