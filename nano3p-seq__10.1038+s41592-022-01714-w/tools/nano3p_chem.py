"""Source-exact custom oligos for Nano3P-seq."""
from chemdraw import Construct, Segment, oligo

R_RNA = "GAAGAUAGAGCGACAGGCAAGUGAUCGGAAG"
D_DNA = "CTTCCGATCACTTGCCTGTCGCTCTATCTTCN"
COMPA_DNA = "GAAGATAGAGCGACAGGCAAGTGATCGGAAGA"


def final_library() -> Construct:
    return Construct([
        Segment("SQK-DCS109 AMX motor adapter", "[SQK-DCS109 AMX motor adapter]",
                placeholder=True, inferred=True),
        Segment("D_DNA-derived cDNA handle", D_DNA),
        Segment("cDNA copied from the native RNA 3-prime end", "X" * 42,
                placeholder=True),
    ], name="Nano3P-seq PCR-free cDNA molecule")


def oligos():
    return (
        oligo("R_RNA", [Segment("RNA starter strand", R_RNA.replace("U", "T"))],
              five="5′-", three="-/3SpC3/-3′"),
        oligo("D_DNA", [Segment("TGIRT starter DNA", D_DNA)],
              five="5′-/5Phos/-", three="-3′"),
        oligo("CompA_DNA", [Segment("adapter-ligation complement", COMPA_DNA)]),
    )
