"""NEBNext Ultra II FS DNA PCR-free, E7430, manual v3.0 9/25."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
import illumina as il
from chemdraw import Construct, Segment
from endprep import dA_tailed_scene, repair_and_dA_tail
from nebnext import STANDARD_PRIMERS, indexed_library

TITLE="NEBNext Ultra II FS DNA PCR-free"
NOTES="01_nebnext-ultra-ii-fs-pcr-free.html"
SOURCE='Commercial protocol: <a href="https://www.neb.com/en-us/products/e7435-nebnext-ultra-ii-fs-dna-pcr-free-library-prep-with-sample-purification-beads">NEB E7430/E7435 manual v3.0 (September 2025)</a>.'
SUMMARY="Intact DNA is enzymatically fragmented, end-repaired and dA-tailed in one reaction, then receives full-length dual-index UMI adaptors. There is no PCR."
UMI_NT=12
I7_READ_NT=20
FS_PROGRAMS={350:(10,37,30,65),450:(8,37,30,65)}
DNA=Construct([Segment("genomic DNA", "X"*50, placeholder=True)], name="intact DNA")
END_PREP=repair_and_dA_tail(Construct([Segment("FS fragment", "X"*42, placeholder=True)], name="FS fragment"))
ADAPTOR=Construct([Segment("full-length UDI-UMI adaptor", "X"*18, "cbc", placeholder=True, inferred=True,
    note="current sequence proprietary; 8-base i5, 8-base i7 and 12-base UMI are documented")], name="pre-annealed adaptor")
STEPS=(("Fragment, repair ends and add dA", dA_tailed_scene(END_PREP), "The single FS reaction couples time-dependent fragmentation to end repair, 5′ phosphorylation and one-base 3′ dA tails."),
       ("Ligate full-length UDI-UMI adaptors", ADAPTOR, "The adaptor is already flow-cell competent: it carries both unique indices and a 12-base UMI. Dotted bases remain proprietary."),
       ("Size-select without amplification", Construct([Segment("adapter-ligated molecules only", "X"*50, placeholder=True)], name="PCR-free pool"), "Bead selection removes free adaptor and unwanted sizes; no polymerase amplification completes or copies this library."))
SEQ_PRIMERS=STANDARD_PRIMERS; REQUIRED_ROLES=tuple(p.role for p in SEQ_PRIMERS)
READOUT="Read 1 and Read 2 sequence the insert; Index 2 reports i5, while the 20-cycle Index 1 read contains the 8-base i7 followed by the 12-base UMI."
FINAL_CAPTION="Full-length PCR-free molecule. The adapter sequence was withdrawn as proprietary by NEB; documented functional blocks and computed primer sites are retained without inventing bases."
def final_library(): return indexed_library(Construct([Segment("DNA insert", "X"*42, placeholder=True)]), "Ultra II FS PCR-free library", umi_nt=UMI_NT, inferred_adaptors=True)
