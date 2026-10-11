"""Altemose et al. DiMeLo-seq directed-methylation model."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from adapter_ligation import AdapterLigationGeometry
from chemdraw import MolecularState,Row,Workflow
from chromatin_epigenetics import dimelo_mark_rows

TITLE = "DiMeLo-seq — antibody-directed methylation on long DNA"
NOTES = "01_dimelo-seq.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1038/s41592-022-01475-6">Altemose et al., <i>Nature Methods</i> (2022)</a>.'
SUMMARY = "An antibody recruits Protein A–Hia5 to a chromatin target; SAM activation deposits nearby m6A marks that are read directly on unamplified long DNA molecules."
CAVEAT = "The scientific construct is the marked native DNA molecule. Oxford Nanopore library-adapter bases are kit-specific and are intentionally represented by role rather than invented sequence."
FINAL_LIBRARY = None
SEQ_PRIMERS = ()
SEQUENCING_ENDING = "Extract high-molecular-weight native DNA, attach the chosen Oxford Nanopore ligation-sequencing motor adapter, and feed one strand through the pore. Modified-base calling reports exogenous m6A together with endogenous CpG methylation on the same molecule; no sequencing primer binds."

def sections():
    ont=AdapterLigationGeometry.unavailable("Altemose et al. use Oxford Nanopore ligation sequencing kits, whose motor-adapter bases and ligation-facing overhang geometry are kit-specific and proprietary.")
    return [
        ("Bind antibody and pA–Hia5", dimelo_mark_rows()[:4],
         "The antibody supplies target specificity while the Protein A fusion positions Hia5 near the chromatin feature."),
        ("Activate directed methylation", dimelo_mark_rows()[4:],
         "SAM activates Hia5; deposited m6A is a covalent mark on nearby adenines, not an adapter or cleavage site."),
        ("Extract native long DNA", [Row(chunks=[("intact marked chromatin DNA → deproteinize gently → long native molecules", None, False)])],
         "The DNA is not PCR-amplified, preserving co-occurring marks and haplotype information along each molecule."),
        ("Attach the Oxford Nanopore motor adapter",ont.rows()+[Row(chunks=[("[motor-loaded ONT adapter] ** native marked DNA → pore", "r2", False)])],
         "The kit-defined adapter enables nanopore translocation; its ligation-facing end is not publicly disclosed."),
    ]

def workflow():
 s=sections();w=Workflow(MolecularState(s[0][0],tuple(s[0][1])));w.react(s[1][0],s[1][1],note=s[1][2]);w.react(s[2][0],s[2][1],note=s[2][2])
 g=AdapterLigationGeometry.unavailable("Altemose et al. use Oxford Nanopore ligation sequencing kits, whose motor-adapter bases and ligation-facing overhang geometry are kit-specific and proprietary.")
 w.ligate_adapter(s[3][0],g,s[3][1],note=s[3][2]);return w
