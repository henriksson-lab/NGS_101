"""CapTrap-seq dual full-length selection and long-read endpoints."""
from chemdraw import Construct, Row, Segment
from protocol_extensions import capped_full_length_cdna
TITLE="CapTrap-seq — cap- and poly(A)-selected full-length cDNA"
NOTES="01_captrap-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/s41467-024-49523-3">Glinos et al., <i>Nature Communications</i> (2024)</a>.'
SUMMARY="Anchored oligo(dT) selects the RNA 3-prime end, cap trapping selects cDNA that reached the native 5-prime cap, and end-dependent linker ligations retain molecules complete at both ends for long-read sequencing."
CAVEAT="CapTrap-seq cDNA was sequenced on both ONT and PacBio. The endpoint below shows the ONT adapter geometry; PacBio SMRTbell conversion is an alternative platform step, not a different CapTrap construct."
FINAL_LIBRARY=Construct([Segment("motor-loaded ONT adapter","[motor-loaded adapter]","r2",placeholder=True,inferred=True),Segment("5-prime linker","L"*14,placeholder=True),Segment("cap-selected full-length cDNA","X"*42,placeholder=True),Segment("copied poly(A) tail","T"*16),Segment("3-prime linker","R"*14,placeholder=True)],name="ONT-adapted CapTrap-seq cDNA")
FINAL_CAPTION="INFERRED — full-length-selected cDNA after platform adapter attachment; ONT adapter bases are proprietary."
SEQ_PRIMERS=()
SEQUENCING_ENDING="There is no synthesis sequencing primer on the ONT branch. A motor-loaded ligation adapter meters the full-length cDNA strand through the pore. For PacBio, hairpin ligation instead creates a SMRTbell and the platform sequencing primer anneals in the hairpin adapter."
def sections():
    return [("Prime at the poly(A) tail and synthesize first-strand cDNA",[Row(chunks=[("5′ cap—RNA—poly(A) || anchored oligo(dT) → cDNA", "r1", False)])],"Anchored oligo(dT) is the first 3-prime-end selection."),("Biotinylate the cap and capture RNA–cDNA hybrids",[Row(chunks=[("biotin—cap—RNA/cDNA ● streptavidin", "w1", False)])],"Only cDNA still connected to a capped RNA 5-prime end is retained."),("Release full-length cDNA and ligate both end-dependent linkers",capped_full_length_cdna().rows(),"Sequential 5-prime and 3-prime linker ligations select molecules carrying both native-end signatures."),("INFERRED — attach a long-read platform adapter",[Row(chunks=[("5′ linker—full-length cDNA—poly(A/T)—3′ linker → ONT or PacBio", None, False)])],"LA-PCR enriches doubly linked molecules before platform-specific library preparation; proprietary adapter bases are not disclosed."),]
