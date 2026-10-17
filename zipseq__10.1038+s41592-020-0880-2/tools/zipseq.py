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
CAVEAT="The two drawn libraries represent the published tumor experiment's 10x 3-prime v3 endpoint (16-base cell barcode, 12-base UMI). The paper used v2 for its in-tube, wound-healing and lymph-node experiments (16-base barcode, 10-base UMI). Hu et al., Methods (in-tube, wound, lymph-node, tumor and Library Construction); original 10x CG00052 Rev D Sequencing and CG000183 Rev A Step 4/oligo appendix. The page does not assign one kit version to every ZipSeq run. The represented tumor branch uses the internal-Cy5 antibody anchor (oligo list line 126) and a live Cy5+ sort (Methods line 272). Library Construction (line 276) spikes a CITE-seq additive primer (partial Small RNA Read 2) into cDNA amplification and recovers zipcode products from the SPRI supernatant before two 3X SPRI cleanups; the additive-primer sequence is not printed by Hu et al. and was not found on the 2017-07-03 CITE-seq protocol pages inspected (pp. 3–6; the cited protocol version is unconfirmed), so it remains an unsequenced input. "+ZIP_READ_LIMITATION
# Published Methods oligonucleotide list (PMC7891292 JATS lines 118-166).
O1_CAGED = "CGATCTGTGGTTGCTACCCGTCAACATCAAAAAAAAAAAA"   # line 130, NPOM-caged
O1_BLOCK = "CGATCTGTGGTTGCTAC"                          # line 134


def oligos():
 yield oligo("antibody anchor, internal Cy5 (Fig. 3, 4; represented tumor branch)",[Segment("anchor 5′ bases","ATCCAG"),Segment("anchor after internal Cy5","TTTTTTTTTTTTGATGTTGACGG")],mods="Amine-C6 linker; internal Cy5 (\\iCy5\\) between ATCCAG and the T run")
 yield oligo("antibody anchor without Cy5 (Fig. 1, S1)",[Segment("anchor sequence","ATCGTTTTTTTTTTTTGATGTTGACGG")],mods="Amine-C6 linker")
 yield oligo("O1 caged strand",[Segment("published O1 caged-strand bases",O1_CAGED)],mods="NPOM-caged T; published list says caged T are highlighted in red but the accessible JATS drops the colour; positions 4, 8, 11, 15 come from the rendered bioRxiv v1 Methods p. 17 (preprint-supplied detail, same bases; published-PDF marking not inspected)")
 yield oligo("O1 blocking strand",[Segment("caged-overhang blocker",O1_BLOCK)])
 yield oligo("zipcode 1",[Segment("O1′ complement + read handle","GTAGCAACCACAGATCGCACCCGAGAATTC"),Segment("8-base region zipcode","CATGATGC","cbc",feature=feature("zipseq_region","spatial_barcode","whitelist",whitelist="paper zipcode list")),Segment("poly(A)28","A"*28)])
def sections(): return [
 ("Coat cells with photocaged anchors",[Row(chunks=[("cell membrane — anchor || caged strand [overhang hidden]",None,False)])],"The anchor may be antibody- or lipid-conjugated."),
 ("Illuminate a selected region",[Row(chunks=[("365-nm light → cage removed → overhang exposed", "w1",False)])],"Only cells in the projected region become competent for that round's zipcode."),
 ("Hybridize a poly(A)-tailed zipcode",[Row(chunks=[("exposed overhang || zipcode — poly(A)28", "cbc",False)])],"Repeated uncaging and hybridization give cells a history of illuminated regions."),
 ("Dissociate and run 10x 3′ RNA-seq",[Row(chunks=[("same cell barcode ├─ endogenous transcript\n                  └─ zipcode feature",None,False)])],"Zipcodes and mRNA are captured by poly(dT) in the same droplet."),]

def protocol_graph():
 """Represented tumor branch (Methods lines 238, 272, 274-276)."""
 g=ProtocolGraph(TITLE,authority=SOURCE)
 cells=g.input("cells","live tumor sections (150 µm) containing poly(A) RNA",MaterialKind.BIOLOGICAL,state=molecule("ZipSeq cells",(role("cell membrane"),role("poly(A) RNA"))),primary=True)
 conjugate=g.input("anchors","anti-CD45 (30F-11) antibody conjugated to the internal-Cy5 amine anchor strand (Thunder-Link)",MaterialKind.MOLECULAR)
 caged=g.input("caged_strand","O1 NPOM-caged strand",MaterialKind.OLIGO,state=O1_CAGED)
 caged_anchor=g.product("caged_anchor","antibody–anchor conjugate prehybridized to the caged O1 strand",MaterialKind.MOLECULAR,state=molecule("caged anchor",(role("antibody anchor"),role("caged O1 overhang"))))
 g.react("prehybridize","Prehybridize caged strand to antibody–anchor conjugate, 1:1, 37 °C 15 min, cool to RT",StepKind.HYBRIDIZATION,((conjugate,InputRole.TEMPLATE),(caged,InputRole.STRUCTURAL)),(caged_anchor,))
 coated=g.product("coated","sections stained with photocaged anchors",MaterialKind.MOLECULAR,state=molecule("anchor-coated cells",(role("poly(A) RNA"),role("caged anchor overhang"))))
 excess=g.product("excess","unbound photocaged anchors",MaterialKind.WASTE)
 g.react("coat","Block with ssDNA/BSA, then stain sections with caged anchor 1 h 4 °C and wash",StepKind.BINDING,((cells,InputRole.TEMPLATE),(caged_anchor,InputRole.ANALYTE)),(coated,excess)); g.discard(excess,"removed by washes")
 light=g.input("light","patterned 365-nm illumination, 1000 ms pulses through the DMD",MaterialKind.CONDITION)
 uncaged=g.product("uncaged","region-selected cells with exposed anchor overhangs",MaterialKind.MOLECULAR,state=molecule("uncaged cells",(role("poly(A) RNA"),role("exposed anchor overhang"))))
 g.react("uncage","Illuminate a selected region to expose anchor overhangs",StepKind.CHEMICAL_CONVERSION,((coated,InputRole.TEMPLATE),(light,InputRole.CONDITION)),(uncaged,))
 zipcode=g.input("zipcode","Zipcode 1 (poly(A)-tailed, O1′ overhang)",MaterialKind.OLIGO,state=ZIP1)
 hybrid=g.product("hybridized","cells carrying the hybridized regional zipcode",MaterialKind.MOLECULAR,state=molecule("zipcode-tagged cells",(role("poly(A) RNA"),role("region zipcode"),role("zipcode poly(A)"))))
 g.react("hybridize","Hybridize Zipcode 1 and wash (1 µM, 10 min in the lymph-node Methods the tumor study follows)",StepKind.HYBRIDIZATION,((uncaged,InputRole.TEMPLATE),(zipcode,InputRole.STRUCTURAL)),(hybrid,))
 block=g.input("blocking_strand","O1 blocking strand (0.25 µM in the lymph-node Methods)",MaterialKind.OLIGO,state=O1_BLOCK)
 tagged=g.product("tagged","zipcoded cells with carried-over zipcode neutralized",MaterialKind.MOLECULAR,state=hybrid.state)
 neutralized=g.product("neutralized","carried-over Zipcode 1 bound by blocking strand",MaterialKind.WASTE)
 g.react("block","Neutralize carried-over zipcode with O1 blocking strand, 5 min; repeat illumination/hybridization for Zipcode 2",StepKind.HYBRIDIZATION,((hybrid,InputRole.TEMPLATE),(block,InputRole.STRUCTURAL)),(tagged,neutralized)); g.discard(neutralized,"removed by washes")
 sorter=g.input("facs","collagenase I/IV dissociation then FACSAria II live Cy5+ sort",MaterialKind.SELECTION)
 sorted_cells=g.product("sorted","live Cy5+ anchor-labelled cells",MaterialKind.MOLECULAR,state=tagged.state)
 unsorted=g.product("unsorted","dead or Cy5-negative cells and debris",MaterialKind.WASTE)
 g.react("sort","Dissociate sections and sort live Cy5+ cells",StepKind.SELECTION,((tagged,InputRole.TEMPLATE),(sorter,InputRole.SELECTION)),(sorted_cells,unsorted),changes_molecule=False); g.discard(unsorted,"excluded by the sort gate")
 droplets=g.input("droplets","10x 3-prime v3 GEM partitioning, 8k target cells",MaterialKind.CONDITION)
 occupied=g.product("occupied","partitioned cells with barcoded beads",MaterialKind.MOLECULAR,state=sorted_cells.state)
 failed=g.product("failed","empty or unresolved GEM products",MaterialKind.WASTE)
 g.react("partition","Partition sorted cells with barcoded beads",StepKind.PARTITION,((sorted_cells,InputRole.TEMPLATE),(droplets,InputRole.CONDITION)),(occupied,failed),changes_molecule=False); g.discard(failed,"does not retain a resolvable cell barcode")
 bead=g.input("bead_primer","10x v3 cell-barcode/UMI oligo(dT) bead primer",MaterialKind.OLIGO)
 rt=g.input("rt","10x reverse transcription reagents",MaterialKind.CATALYST)
 cdna=g.product("cdna","cell-barcoded transcript cDNA",MaterialKind.MOLECULAR,state=molecule("ZipSeq cDNA",(role("cell barcode"),role("UMI"),role("cDNA"))))
 copied=g.product("copied","cell-barcoded zipcode copies",MaterialKind.MOLECULAR,state=zipcode_captured_cdna())
 g.react("reverse_transcribe","In-GEM bead primer captures poly(A) mRNA and poly(A) zipcodes; RT copies both",StepKind.REVERSE_TRANSCRIPTION,((occupied,InputRole.TEMPLATE),(bead,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(cdna,copied))
 cdna_primers=g.input("cdna_primers","10x cDNA amplification primers",MaterialKind.OLIGO)
 additive=g.input("additive_primer","CITE-seq additive primer (partial Small RNA Read 2)",MaterialKind.OLIGO)
 amp_pol=g.input("cdna_polymerase","10x cDNA amplification mix",MaterialKind.CATALYST)
 amp_cdna=g.product("amp_cdna","amplified full-length cDNA",MaterialKind.MOLECULAR,state=cdna.state)
 amp_zip=g.product("amp_zip","amplified short zipcode products (<200 bp)",MaterialKind.MOLECULAR,state=copied.state)
 g.react("cdna_amp","Amplify cDNA with the additive primer spiked in",StepKind.PCR,
     ((cdna,InputRole.TEMPLATE),(copied,InputRole.TEMPLATE),(cdna_primers,InputRole.STRUCTURAL),(additive,InputRole.STRUCTURAL),(amp_pol,InputRole.CATALYST)),
     (amp_cdna,amp_zip),geometry=PCRGeometry.unavailable(copied.state,reason="Hu et al. name the additive primer only as partial Read 2 Small RNA and follow CITE-seq; neither its sequence nor the 10x cDNA-primer sequences are printed in the article."))
 spri=g.input("spri_split","post-amplification SPRI selection",MaterialKind.SELECTION)
 rna=g.product("rna_branch","bead-bound cDNA fraction",MaterialKind.MOLECULAR,state=amp_cdna.state)
 zips=g.product("zip_branch","supernatant zipcode fraction",MaterialKind.MOLECULAR,state=amp_zip.state)
 g.react("split","SPRI: keep bead-bound cDNA and decant the zipcode-containing supernatant",StepKind.SELECTION,((amp_cdna,InputRole.TEMPLATE),(amp_zip,InputRole.TEMPLATE),(spri,InputRole.SELECTION)),(rna,zips),changes_molecule=False)
 spri3=g.input("spri_3x","two successive 3X SPRI cleanups",MaterialKind.SELECTION)
 zip_clean=g.product("zip_clean","purified zipcode fraction",MaterialKind.MOLECULAR,state=zips.state)
 zip_waste=g.product("zip_spri_waste","3X SPRI supernatants",MaterialKind.WASTE)
 g.react("zip_cleanup","Purify the zipcode fraction with two successive 3X SPRI cleanups",StepKind.CLEANUP,((zips,InputRole.TEMPLATE),(spri3,InputRole.SELECTION)),(zip_clean,zip_waste)); g.discard(zip_waste,"removed with the bead supernatants")
 lib_kit=g.input("tenx_library_reagents","10x v3 fragmentation, end-repair/A-tailing and adapter-ligation reagents",MaterialKind.CATALYST)
 ligated=g.product("rna_ligated","fragmented, adapter-ligated cDNA",MaterialKind.MOLECULAR,state=molecule("ligated cDNA",(role("cell barcode"),role("UMI"),role("cDNA fragment"),role("Read 2 adapter"))))
 g.react("rna_construct","10x library construction of the bead-bound cDNA",StepKind.MAGIC,((rna,InputRole.TEMPLATE),(lib_kit,InputRole.CATALYST)),(ligated,),
     reason="Hu et al. process the cDNA fraction 'as recommended in the 10X user guide'; the CG000183 fragmentation, ligation and adapter geometry are not modelled in this pass.")
 reviewed=["cells","anchors","caged_strand","prehybridize","coat","light","uncage","zipcode","hybridize","blocking_strand","block","facs","sort","droplets","partition","bead_primer","rt","reverse_transcribe","cdna_primers","additive_primer","cdna_polymerase","cdna_amp","spri_split","split","spri_3x","zip_cleanup","tenx_library_reagents","rna_construct"]
 finish_pcr_branch(g,"rna",ligated,RNA_LIB,RNA_PRIMERS,reviewed,label="ZipSeq gene expression",read_lengths=READ_LENGTHS["Gene-expression library (v3 example)"])
 zip_forward=g.input("zip_forward", "CITE-seq 10x P5 PCR primer", MaterialKind.OLIGO,
     state=Construct([Segment("10x P5 primer", CITE_P5)]), source="CITE-seq original supplemental protocol p.6")
 zip_reverse=g.input("zip_reverse", "CITE-seq RPIx indexed PCR-primer family", MaterialKind.OLIGO,
     state=Construct([Segment("RPIx architecture", RPI_FAMILY)]), inferred=True,
     source="Hu Library Construction inherits CITE-seq primers; exact RPI member not assigned")
 zip_pol=g.input("zip_polymerase", "2X KAPA HiFi HotStart Master Mix", MaterialKind.CATALYST)
 zip_library=g.product("zip_library", "INFERRED Zipcode 1 PCR endpoint", MaterialKind.LIBRARY,
     state=ZIP_LIB, inferred=True, source="Hu native zipcode and cited CITE-seq PCR architecture")
 g.react("zip_pcr", "INFERRED — amplify zipcode branch with cited CITE-seq primer classes", StepKind.PCR,
     ((zip_clean,InputRole.TEMPLATE),(zip_forward,InputRole.STRUCTURAL),(zip_reverse,InputRole.STRUCTURAL),
      (zip_pol,InputRole.CATALYST)),(zip_library,),geometry=require_zipcode_primer_sites(zip_clean.state,native_template=True))
 g.sequence(zip_library, SequencingPlan(library=ZIP_LIB,primers=ZIP_PRIMERS,
     acquired_roles=("Read 1","Index 1 (i7)","Read 2"),
     read_lengths=READ_LENGTHS["Zipcode 1 library (conditional primer-derived view)"],
     read_uses={"Read 2":"extract"},limitation=ZIP_READ_LIMITATION))
 reviewed.extend((zip_forward.id,zip_reverse.id,zip_pol.id,"zip_pcr",zip_library.id))
 g.declare_reviewed_scope(reviewed,SOURCE)
 return g.require_publishable()
