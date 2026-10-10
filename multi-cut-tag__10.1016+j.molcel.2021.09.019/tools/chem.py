"""Multi-CUT&Tag target-coded transposition."""
from chemdraw import Row
from multibranch import illumina_branch
from protocol_extensions import multiplex_target_fragment
TITLE="Multi-CUT&Tag — multiplex chromatin targets in the same cell"
NOTES="01_multi-cut-tag.html"
SOURCE='Defining source: <a href="https://doi.org/10.1016/j.molcel.2021.09.019">Gopalan et al., <i>Molecular Cell</i> (2021)</a>.'
SUMMARY="Each antibody-targeted pA–Tn5 complex carries a target-specific adapter barcode; subsequent cell indexing preserves both target and cellular identity on every chromatin fragment."
FINAL_LIBRARY,SEQ_PRIMERS=illumina_branch("Multi-CUT&Tag","targeted chromatin",architecture="nextera",cell=12,target=8)
FINAL_CAPTION="The target barcode names the tethered antibody/transposome and the cell barcode names the cell."
SEQUENCING_INTRO="Declared Nextera primers bind the final multiplex library."
def sections():
    return [("Bind antibodies to several chromatin targets",[Row(chunks=[("chromatin target ← antibody ← pA–Tn5[target code]", "cbc", False)])],"Distinct loaded transposomes are assigned to distinct antibody targets."),("Tagment and transfer the target code",multiplex_target_fragment(target_barcode_length=8).rows(),"Tn5 cleavage transfers a target-identifying adapter to the adjacent DNA."),("Add cell and sample indexes",[Row(chunks=[("target-coded fragment ** cell index ** sequencing adapters", "cbc", False)])],"Later indexing makes target × cell combinations readable in one pool."),]
