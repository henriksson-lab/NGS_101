"""PAIso-seq end extension, USER release and PacBio SMRTbell conversion."""
from batch_ngs import joined_scene, seg
from chemdraw import Construct, Row, feature
from longread_footprinting import conventional_smrtbell, pacbio_entry, smrtbell_rows
from protocol_extensions import capped_full_length_cdna
TITLE="PAIso-seq — poly(A)-inclusive full-length isoform sequencing"
NOTES="01_paiso-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/s41467-019-13228-9">Liu et al., <i>Nature Communications</i> (2019)</a>; detailed protocol <a href="https://doi.org/10.1038/s41596-022-00704-8">Nat. Protocols (2022)</a>.'
SUMMARY="A dU-interrupted oligo templates end extension across the complete RNA poly(A) tail, USER releases the copied product, template switching captures the 5-prime end, and PacBio CCS reads the entire isoform plus native tail."
INSERT=Construct([seg("5-prime cDNA handle","AAGCAGTGGTATCAACGCAGAGT","tso"),seg("full-length cDNA","X"*38,placeholder=True),seg("copied native poly(A) tail","T"*18),seg("sample barcode","B"*16,"cbc",placeholder=True,feature=feature("sample_barcode","sample_index","whitelist",whitelist="published PAIso-seq sample-barcode set"))],name="PAIso-seq full-length cDNA")
SMRTBELL=conventional_smrtbell("PAIso-seq full-length cDNA")
FINAL_LIBRARY=None; SEQ_PRIMERS=(); SEQUENCING_ENDING=pacbio_entry(SMRTBELL)
def oligos():
    return ["<b>TSO</b> — 5′-AAGCAGTGGTATCAACGCAGAGTACATrGrG+G-3′","<b>RT primer</b> — 5′-AAGCAGTGGTATCAACGCAGAGTAC-3′","<b>PCR primer</b> — 5′-AAGCAGTGGTATCAACGCAGAGT-3′"]
def sections():
    return [("Template an end extension across the native poly(A) tail",[Row(chunks=[("RNA body—poly(A) ← dU-interrupted barcoded oligo", "r1", False)])],"The designed oligo lets polymerase copy the complete tail rather than prime within it."),("Open the dU sites with USER and reverse-transcribe",capped_full_length_cdna(retain_poly_a=True).rows(),"USER removes the templating oligo; RT and template switching place a known handle at the transcript 5-prime end."),("Amplify full-length cDNA",[Row(chunks=[("TSO handle—full-length cDNA—copied native tail—barcode", "tso", False)])],"The tail remains inside the amplicon rather than being replaced by an oligo(dT) tract."),("Ligate PacBio hairpins",smrtbell_rows(SMRTBELL),"Damage repair, end repair and hairpin ligation produce a closed CCS template."),]
