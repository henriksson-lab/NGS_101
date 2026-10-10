"""CELLO-seq long-UMI, full-length single-cell cDNA."""
from chemdraw import Construct, Row, Segment, feature
from protocol_extensions import capped_full_length_cdna
TITLE="CELLO-seq — long-UMI single-cell long-read RNA-seq"
NOTES="01_cello-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/s41587-021-01093-1">Berrens et al., <i>Nature Biotechnology</i> (2022)</a>; detailed protocol <a href="https://doi.org/10.1038/s41596-025-01203-2">Nat. Protocols (2025)</a>.'
SUMMARY="High-temperature template-switch RT creates full-length cDNA, then a splint ligation adds a cell barcode and 22-nt RYN-patterned UMI; deliberate duplicate generation enables long-read consensus correction."
CAVEAT="The ONT ligation adapter's molecular role is established, but its proprietary bases are not disclosed; that terminal segment is inferred."
FINAL_LIBRARY=Construct([Segment("motor-loaded ONT adapter","[motor adapter]","r2",placeholder=True,inferred=True),Segment("template-switch end","X"*18,"tso",placeholder=True),Segment("full-length cDNA","X"*42,placeholder=True),Segment("cell barcode","B"*8,"cbc",placeholder=True,feature=feature("cell_barcode","cell_barcode","unknown")),Segment("22-nt RYN UMI","U"*22,"umi",placeholder=True,feature=feature("umi","umi","random",note="constrained RYN repeat design")),Segment("oligo(dT)-derived end","T"*18)],name="ONT-adapted CELLO-seq cDNA")
FINAL_CAPTION="INFERRED — full-length cDNA with splint-ligated cell barcode and 22-nt patterned UMI after proprietary nanopore-adapter attachment."
SEQ_PRIMERS=()
SEQUENCING_ENDING="There is no synthesis sequencing primer. The ONT ligation adapter's motor meters the adapted cDNA through the pore; reads sharing the 22-nt patterned UMI are grouped for consensus correction."
def sections():
    return [("Reverse-transcribe at high temperature and template-switch",capped_full_length_cdna(retain_poly_a=True).rows(),"A 3-prime-amino-blocked TSO suppresses TSO-primed artefact inserts while retaining full-length cDNA."),("Digest unused oligos and splint-ligate identifiers",[Row(chunks=[("TSO—full-length cDNA—oligo(dT) ** cell barcode—22-nt RYN UMI", "umi", False)])],"HiFi Taq ligase joins the annealed splint carrying cell and molecule identity."),("Amplify deliberately to obtain UMI families",[Row(chunks=[("one barcoded cDNA → many PCR copies with the same 22-nt UMI", "umi", False)])],"High duplicate depth is intentional: multiple noisy long reads form one corrected consensus."),("INFERRED — attach ONT ligation adapters",[Row(chunks=[("motor adapter ** full-length barcoded cDNA → pore", "r2", False)])],"The published experiments use ONT ligation kits; the proprietary adapter sequence is not disclosed."),]
