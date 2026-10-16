"""ZipSeq photocaged spatial barcode model."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from chemdraw import Construct, Row, Segment, feature, oligo, revcomp
import illumina as il
import seqprimers as sp
from pcr import PCRGeometry, PCRPrimer
from multimodal_spatial import finish_pcr_branch
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import molecule, role
ZIP1="GTAGCAACCACAGATCGCACCCGAGAATTCCATGATGC"+"A"*28
# This page represents the tumor experiment's v3 branch. Hu Methods delegates
# GEX to 10x and zipcode amplification to the original CITE-seq protocol.
RNA_PRIMERS = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["R2"])
CELL = feature("cell_barcode", "cell_barcode", "whitelist", whitelist="10x Chromium 3-prime v3")
UMI = feature("umi", "umi", "random")
INDEX = feature("sample_index_i7", "sample_index", "unknown")
REGION = feature("zipseq_region", "spatial_barcode", "whitelist", whitelist="Hu named zipcodes")
O1_PRIME = revcomp("CGATCTGTGGTTGCTAC")
ZIPCODE1 = ZIP1[-36:-28]
ZIP_HANDLE = ZIP1[len(O1_PRIME):-36]
# CITE-seq original supplemental protocol, 2017-07-03, p.6, RPI1.
CITE_RPI1 = "CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
CITE_P5 = il.NEBNEXT_UNIVERSAL_PRIMER
SMALL_RNA_R2_SITE = CITE_RPI1[len(il.P7)+6:]
RPI_FAMILY = il.P7 + "N"*6 + SMALL_RNA_R2_SITE
ZIP_READ_LIMITATION = ("INFERRED — The zipcode endpoint and I1/R2 binding view are reconstructed from the cited CITE-seq RPIx architecture. "
                      "The loaded sequencing-primer sequences and exact RPI member for Hu's run have not been independently verified. "
                      "The PCR-derived reverse-primer footprint overlaps the initial CA of the named eight-base zipcode; its actual sequencing entry is not assigned from that footprint alone.")
ZIP_PRIMERS = ()  # Loaded zipcode sequencing oligos remain under source review.


def _cell_end(*, completed):
    left = ([Segment("P5", il.P5, "p5"),
             Segment("TruSeq Read 1 arm", il.TRUSEQ_READ1[4:], "r1")]
            if completed else [Segment("gel-bead Read 1 handle", il.TRUSEQ_READ1[11:], "r1")])
    return [*left,
        Segment("cell barcode", "B"*16, "cbc", placeholder=True, feature=CELL),
        Segment("UMI", "U"*12, "umi", placeholder=True, feature=UMI)]


def require_zipcode_primer_sites(library, *, native_template=False):
    """Check the endpoint against the separately printed predecessor PCR pair."""
    return PCRGeometry.exact(library, PCRPrimer("CITE-seq 10x P5", CITE_P5, anneal_nt=22 if native_template else len(CITE_P5)),
        PCRPrimer("CITE-seq RPIx family", RPI_FAMILY, anneal_nt=len(ZIP_HANDLE)+2 if native_template else len(SMALL_RNA_R2_SITE)),
        statement="Hu Methods inherits these primer classes from CITE-seq; exact RPI index member remains unspecified.")


def _gene_expression_library():
    library = Construct([*_cell_end(completed=True),
        Segment("poly(dT)30", "T"*30), Segment("cDNA insert", "X"*40, placeholder=True),
        Segment("dA junction", "A"), Segment("Index 1 / Read 2 arm", il.INDEX1_PRIMER, "r2"),
        Segment("i7 reverse complement", "I"*8, "cbc", placeholder=True, feature=INDEX),
        Segment("P7 reverse complement", il.P7_RC, "p7")], "ZipSeq v3 gene-expression library")
    problems = sp.verify(library, RNA_PRIMERS, required_roles=tuple(p.role for p in RNA_PRIMERS))
    if problems: raise ValueError("ZipSeq GEX: "+"; ".join(problems))
    return library


def zipcode_captured_cdna():
    return Construct([*_cell_end(completed=False), Segment("poly(dT)30", "T"*30),
        Segment("Zipcode 1 reverse complement", revcomp(ZIPCODE1), "cbc", feature=REGION),
        Segment("copied partial Small RNA handle", revcomp(ZIP_HANDLE), "r2"),
        Segment("copied O1-prime overhang", revcomp(O1_PRIME))], "captured Zipcode 1 cDNA")


def _zipcode_library():
    native_site = ZIP_HANDLE + ZIPCODE1[:2]
    if not SMALL_RNA_R2_SITE.endswith(native_site):
        raise ValueError("Hu native zipcode no longer supplies the CITE-seq RPI 3-prime binding region")
    library = Construct([*_cell_end(completed=True), Segment("poly(dT)30", "T"*30),
        Segment("Zipcode 1 reverse complement", revcomp(ZIPCODE1), "cbc", feature=REGION),
        Segment("copied partial Small RNA handle", revcomp(ZIP_HANDLE), "r2"),
        Segment("RPI-derived extension of Small RNA arm", revcomp(SMALL_RNA_R2_SITE[:-len(native_site)]), "r2", inferred=True),
        Segment("RPI i7 reverse complement", "I"*6, "cbc", placeholder=True, feature=INDEX, inferred=True),
        Segment("P7 reverse complement", il.P7_RC, "p7", inferred=True)], "INFERRED Zipcode 1 PCR endpoint")
    require_zipcode_primer_sites(library)
    return library


RNA_LIB = _gene_expression_library()
ZIP_LIB = _zipcode_library()
READ_LENGTHS = {"Gene-expression library (v3 example)": {"Read 1":28, "Index 1 (i7)":8, "Read 2":91},
                "Zipcode 1 library (conditional primer-derived view)": {"Read 1":28, "Index 1 (i7)":8, "Read 2":91}}
READ_SPECS = {"Zipcode 1 library (conditional primer-derived view)": (
    sp.read_spec("Read 1",28,"16-base cell barcode;12-base UMI (represented v3 branch)","extract"),
    sp.read_spec("Index 1 (i7)",8,"sample index;exact RPI member/index extent not yet verified","demultiplex"),
    sp.read_spec("Index 2 (i5)","not acquired","not acquired in the represented single-index branch","not_applicable"),
    sp.read_spec("Read 2",91,"zipcode identity;exact entry relative to its eight-base identifier remains unverified","extract"),
)}
FINAL_LIBRARIES = (
    ("Gene-expression library (v3 example)", RNA_LIB, RNA_PRIMERS,
     "10x 3-prime v3 transcript endpoint; the paper also used v2 for other experiments.",
     "The represented v3 run follows CG000183: 28/8/91 cycles. Read 1 reports cell barcode and UMI; Read 2 enters cDNA."),
    ("Zipcode 1 library (conditional primer-derived view)", ZIP_LIB, ZIP_PRIMERS,
     "INFERRED — Zipcode 1 endpoint reconstructed using Hu's native tag and the cited CITE-seq PCR-primer architecture.",
     ZIP_READ_LIMITATION))
TITLE="ZipSeq — sequential photo-uncaging of spatial zipcodes"
NOTES="01_zipseq.html"
SOURCE='<a href="https://doi.org/10.1038/s41592-020-0880-2">Hu et al., <i>Nature Methods</i> (2020)</a>.'
SUMMARY="A photocaged DNA anchor coats cells. Patterned illumination exposes an overhang, and a complementary polyadenylated zipcode is added after each illumination round. 10x RNA-seq reads the accumulated spatial tags with the transcriptome."
CAVEAT="The two drawn libraries represent the published tumor experiment's 10x 3-prime v3 endpoint (16-base cell barcode, 12-base UMI). The paper used v2 for its in-tube, wound-healing and lymph-node experiments (16-base barcode, 10-base UMI). Hu et al., Methods (in-tube, wound, lymph-node, tumor and Library Construction); original 10x CG00052 Rev D Sequencing and CG000183 Rev A Step 4/oligo appendix. The page does not assign one kit version to every ZipSeq run. "+ZIP_READ_LIMITATION
def oligos():
 yield oligo("antibody anchor",[Segment("anchor sequence","ATCGTTTTTTTTTTTTGATGTTGACGG")],mods="Amine-C6 linker")
 yield oligo("O1 caged strand",[Segment("published O1 caged-strand bases","CGATCTGTGGTTGCTACCCGTCAACATCAAAAAAAAAAAA")],mods="NPOM-modified T at positions 4, 8, 11 and 15 in the rendered bioRxiv v1 Methods, p. 17; publisher JATS loses colour")
 yield oligo("O1 blocking strand",[Segment("caged-overhang blocker","CGATCTGTGGTTGCTAC")])
 yield oligo("zipcode 1",[Segment("O1′ complement + read handle","GTAGCAACCACAGATCGCACCCGAGAATTC"),Segment("8-base region zipcode","CATGATGC","cbc",feature=feature("zipseq_region","spatial_barcode","whitelist",whitelist="paper zipcode list")),Segment("poly(A)28","A"*28)])
def sections(): return [
 ("Coat cells with photocaged anchors",[Row(chunks=[("cell membrane — anchor || caged strand [overhang hidden]",None,False)])],"The anchor may be antibody- or lipid-conjugated."),
 ("Illuminate a selected region",[Row(chunks=[("365-nm light → cage removed → overhang exposed", "w1",False)])],"Only cells in the projected region become competent for that round's zipcode."),
 ("Hybridize a poly(A)-tailed zipcode",[Row(chunks=[("exposed overhang || zipcode — poly(A)28", "cbc",False)])],"Repeated uncaging and hybridization give cells a history of illuminated regions."),
 ("Dissociate and run 10x 3′ RNA-seq",[Row(chunks=[("same cell barcode ├─ endogenous transcript\n                  └─ zipcode feature",None,False)])],"Zipcodes and mRNA are captured by poly(dT) in the same droplet."),]

def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 cells=g.input("cells","cells containing poly(A) RNA",MaterialKind.BIOLOGICAL,state=molecule("ZipSeq cells",(role("cell membrane"),role("poly(A) RNA"))),primary=True)
 anchors=g.input("anchors","membrane-conjugated photocaged DNA anchors",MaterialKind.MOLECULAR)
 coated=g.product("coated","cells coated with photocaged anchors",MaterialKind.MOLECULAR,state=molecule("anchor-coated cells",(role("poly(A) RNA"),role("caged anchor overhang"))))
 excess=g.product("excess","unbound photocaged anchors",MaterialKind.WASTE)
 g.react("coat","Coat cells with photocaged DNA anchors",StepKind.BINDING,((cells,InputRole.TEMPLATE),(anchors,InputRole.ANALYTE)),(coated,excess)); g.discard(excess,"is not attached to a retained cell")
 light=g.input("light","patterned 365-nm illumination condition",MaterialKind.CONDITION)
 uncaged=g.product("uncaged","region-selected cells with exposed anchor overhangs",MaterialKind.MOLECULAR,state=molecule("uncaged cells",(role("poly(A) RNA"),role("exposed anchor overhang"))))
 g.react("uncage","Illuminate a selected region to expose anchor overhangs",StepKind.CHEMICAL_CONVERSION,((coated,InputRole.TEMPLATE),(light,InputRole.CONDITION)),(uncaged,))
 zipcode=g.input("zipcode","complementary poly(A)-tailed zipcode oligo",MaterialKind.OLIGO,state=ZIP1)
 tagged=g.product("tagged","cells carrying the hybridized regional zipcode",MaterialKind.MOLECULAR,state=molecule("zipcode-tagged cells",(role("poly(A) RNA"),role("region zipcode"),role("zipcode poly(A)"))))
 g.react("hybridize","Hybridize the round-specific poly(A)-tailed zipcode",StepKind.HYBRIDIZATION,((uncaged,InputRole.TEMPLATE),(zipcode,InputRole.STRUCTURAL)),(tagged,))
 droplets=g.input("droplets","10x 3-prime v3 GEM partitioning condition for represented tumor branch",MaterialKind.CONDITION)
 occupied=g.product("occupied","cell-barcoded transcript and zipcode products",MaterialKind.MOLECULAR,state=tagged.state)
 failed=g.product("failed","empty or unresolved GEM products",MaterialKind.WASTE)
 g.react("partition","Dissociate and partition tagged cells with barcoded beads",StepKind.PARTITION,((tagged,InputRole.TEMPLATE),(droplets,InputRole.CONDITION)),(occupied,failed),changes_molecule=False); g.discard(failed,"does not retain a resolvable cell barcode")
 rna=g.product("rna_branch","gene-expression branch",MaterialKind.MOLECULAR,state=occupied.state); zips=g.product("zip_branch","zipcode-feature branch",MaterialKind.MOLECULAR,state=occupied.state)
 g.react("split","Separate transcript and zipcode products",StepKind.SPLIT,((occupied,InputRole.TEMPLATE),),(rna,zips),changes_molecule=False)
 bead=g.input("bead_primer","cell-barcode/UMI oligo(dT) bead primer",MaterialKind.OLIGO)
 rt=g.input("rt","reverse transcriptase",MaterialKind.CATALYST)
 cdna=g.product("cdna","cell-barcoded transcript cDNA",MaterialKind.MOLECULAR,state=molecule("ZipSeq cDNA",(role("cell barcode"),role("UMI"),role("cDNA"))))
 g.react("reverse_transcribe","Capture and reverse-transcribe endogenous mRNA",StepKind.REVERSE_TRANSCRIPTION,((rna,InputRole.TEMPLATE),(bead,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(cdna,))
 zpol=g.input("zipcode_extension","zipcode barcode-transfer polymerase",MaterialKind.CATALYST)
 copied=g.product("copied","cell-barcoded zipcode products",MaterialKind.MOLECULAR,state=zipcode_captured_cdna())
 g.react("copy_zipcodes","Extend the same bead primer across poly(A)-tailed zipcodes",StepKind.EXTENSION,((zips,InputRole.TEMPLATE),(bead,InputRole.STRUCTURAL),(zpol,InputRole.CATALYST)),(copied,))
 reviewed=["cells","anchors","coat","light","uncage","zipcode","hybridize","droplets","partition","split","bead_primer","rt","reverse_transcribe","zipcode_extension","copy_zipcodes"]
 finish_pcr_branch(g,"rna",cdna,RNA_LIB,RNA_PRIMERS,reviewed,label="ZipSeq gene expression",read_lengths=READ_LENGTHS["Gene-expression library (v3 example)"])
 zip_forward=g.input("zip_forward", "CITE-seq 10x P5 PCR primer", MaterialKind.OLIGO,
     state=Construct([Segment("10x P5 primer", CITE_P5)]), source="CITE-seq original supplemental protocol p.6")
 zip_reverse=g.input("zip_reverse", "CITE-seq RPIx indexed PCR-primer family", MaterialKind.OLIGO,
     state=Construct([Segment("RPIx architecture", RPI_FAMILY)]), inferred=True,
     source="Hu Library Construction inherits CITE-seq primers; exact RPI member not assigned")
 zip_pol=g.input("zip_polymerase", "zipcode PCR polymerase", MaterialKind.CATALYST)
 zip_library=g.product("zip_library", "INFERRED Zipcode 1 PCR endpoint", MaterialKind.LIBRARY,
     state=ZIP_LIB, inferred=True, source="Hu native zipcode and cited CITE-seq PCR architecture")
 g.react("zip_pcr", "INFERRED — amplify zipcode branch with cited CITE-seq primer classes", StepKind.PCR,
     ((copied,InputRole.TEMPLATE),(zip_forward,InputRole.STRUCTURAL),(zip_reverse,InputRole.STRUCTURAL),
      (zip_pol,InputRole.CATALYST)),(zip_library,),geometry=require_zipcode_primer_sites(copied.state,native_template=True))
 g.sequence(zip_library, SequencingPlan(library=ZIP_LIB,primers=ZIP_PRIMERS,
     acquired_roles=("Read 1","Index 1 (i7)","Read 2"),
     read_lengths=READ_LENGTHS["Zipcode 1 library (conditional primer-derived view)"],
     read_uses={"Read 2":"extract"},limitation=ZIP_READ_LIMITATION))
 reviewed.extend((zip_forward.id,zip_reverse.id,zip_pol.id,"zip_pcr",zip_library.id))
 g.declare_reviewed_scope(reviewed,SOURCE)
 return g.require_publishable()
