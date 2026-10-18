"""10x Chromium Next GEM Single Cell Multiome ATAC + Gene Expression."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from chemdraw import Construct, Row, Scene, Segment, feature, revcomp
import illumina as il
import nextera as nx
import seqprimers as sp
from multimodal_spatial import modality_split
from pcr import PCRGeometry
from adapter_ligation import AdapterLigationGeometry
from protocol_graph import (Disposition, InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, SourceClaim, SourceDiscrepancy, StepKind, Terminal)
from protocol_model import molecule, role

seg=Segment

TITLE = "10x Chromium Single Cell Multiome ATAC + Gene Expression"
NOTES = "01_multiome.html"
SOURCE = 'Commercial protocol: 10x Genomics <a href="https://www.10xgenomics.com/support/instruments/chromium-x-series/chromium-next-gem-single-cell-multiome-atac-plus-gene-expression-reagent-kits-user-guide">CG000338 Rev G</a>.'
SUMMARY = "Nuclei are tagmented before GEM formation. Inside each GEM, one gel-bead barcode is transferred to both accessible-DNA fragments and poly(A)-derived cDNA, yielding linked ATAC and gene-expression libraries from the same nucleus."
CAVEAT = ("CG000338 RevG prints both gel-bead structures and complete final-library layouts. "
          "Selected barcode/UMI/index values stay variable. Full physical transposon, TSO, "
          "preamp/cDNA/index primer and adaptor structures and the ATAC barcode-joining "
          "mechanism remain unassigned in this model. Their inferred regions and source-limited "
          "reactions remain explicit. SI-PCR PrimerB is catalog-linked to PN-2000128 in CG000168 RevA.")

CELL_BARCODE = feature("cell_barcode", "cell_barcode", "whitelist",
                       whitelist="10x-multiome-gem-barcode")
UMI_FEATURE = feature("umi", "umi", "random")
I7_FEATURE = feature("sample_index_i7", "sample_index", "unknown")

GEX = Construct([
    seg("P5",il.P5,"p5"),
    seg("i5 sample index","J"*10,"cbc",placeholder=True,
        feature=feature("sample_index_i5","sample_index","unknown")),
    seg("TruSeq Read 1",il.TRUSEQ_READ1,"r1"),
    seg("cell barcode","B"*16,"cbc",placeholder=True,feature=CELL_BARCODE),
    seg("UMI","U"*12,"umi",placeholder=True,feature=UMI_FEATURE),
    seg("poly(dT)30","T"*30),seg("VN anchor","VN",placeholder=True,length_bp=2),
    seg("cDNA insert","X"*36,placeholder=True),
    seg("Read 2 arm reverse complement",revcomp(il.TRUSEQ_READ2),"r2"),
    seg("i7 sample index reverse complement","I"*10,"cbc",placeholder=True,
        feature=I7_FEATURE),
    seg("P7 reverse complement",il.P7_RC,"p7")],
    name="Multiome gene-expression library")
GEX_PRIMERS=tuple(sp.TRUSEQ[k] for k in ("R1","I1","I2","R2"))
_gex_problems=sp.verify(GEX,GEX_PRIMERS,required_roles=tuple(p.role for p in GEX_PRIMERS))
if _gex_problems:
    raise ValueError("Multiome gene-expression library: "+"; ".join(_gex_problems))
ATAC = Construct([
    seg("P5", il.P5, "p5"),
    seg("10x cell barcode", "B" * 16, "cbc", placeholder=True,
        feature=CELL_BARCODE),
    seg("10x spacer","CGCGTCTG"), seg("S5", nx.S5, "s5"),
    seg("left mosaic end", nx.ME, "me"),
    seg("accessible genomic DNA", "X" * 34, placeholder=True),
    seg("right mosaic end reverse complement", nx.ME_RC, "me"),
    seg("S7 reverse complement", nx.S7_RC, "s7"),
    seg("i7 sample index", "I" * 8, "cbc", placeholder=True, feature=I7_FEATURE),
    seg("P7 reverse complement", il.P7_RC, "p7")], name="Multiome ATAC library")
ATAC_PRIMERS = tuple(sp.NEXTERA[k] for k in ("R1", "I1", "I2", "R2"))
_atac_problems = sp.verify(ATAC, ATAC_PRIMERS,
                           required_roles=tuple(p.role for p in ATAC_PRIMERS))
if _atac_problems:
    raise ValueError("Multiome ATAC library: " + "; ".join(_atac_problems))

FINAL_LIBRARIES = (
    ("Final gene-expression library", GEX, GEX_PRIMERS,
     "Source-exact CG000338 Rev G layout. Read 1 reports the GEM cell barcode and UMI.",
     "The declared Illumina primers are located on the final gene-expression duplex."),
    ("Final ATAC library", ATAC, ATAC_PRIMERS,
     "Source-exact CG000338 Rev G layout. Paired genomic reads flank the accessible fragment; the 24-cycle Index 2 read covers the 8-base spacer and the 16-base GEM barcode, and Index 1 reads the sample index.",
     "The declared Nextera sequencing primers are located on the final ATAC duplex."),
)

def oligos():
    return [
        "GEX gel-bead oligo: 5′-CTACACGACGCTCTTCCGATCT-[16-bp cell barcode]-[12-bp UMI]-T30VN-3′",
        "ATAC gel-bead oligo: 5′-AATGATACGGCGACCACCGAGATCTACAC-[16-bp cell barcode]-CGCGTCTG-3′",
        "Final GEX shell: P5-[10-bp i5]-TruSeq Read 1-[cell barcode]-[UMI]-T30VN-cDNA-Index 1/Read 2 arm-[10-bp i7]-P7′",
        "Final ATAC shell: P5-[cell barcode]-CGCGTCTG-s5-ME-insert-ME′-s7′-[8-bp i7]-P7′",
    ]

def sections():
    return [
        ("Tagment accessible chromatin before partitioning",
         [Row(chunks=[("nucleus chromatin → Tn5 cleavage + adapter transfer → tagged open-DNA fragments", "me", False)])],
         "Transposition (ATAC Buffer B + ATAC Enzyme B, 37 °C 60 min) occurs in bulk nuclei; intact nuclei retain both tagged chromatin and RNA."),
        ("Transfer one GEM barcode to RNA and ATAC molecules",
         [Row(chunks=[("poly(A) RNA + barcoded oligo-dT + Template Switch Oligo → cell barcode + UMI + cDNA", "cbc", False)]),
          Row(chunks=[("tagged genomic DNA + P5/barcode/spacer gel-bead oligo → cell barcode + ATAC fragment", "cbc", False)])],
         "Both reactions occur in the same GEM incubation (37 °C 45 min, 25 °C 30 min, then Quenching Agent) and therefore use the same published 16-base barcode; CG000338 prints each gel-bead structure but not the joining mechanism for the ATAC oligo."),
        ("Pre-amplify both modalities together, then split", modality_split("gene-expression cDNA library", "ATAC DNA library"),
         "After Dynabeads and SPRIselect cleanup, one 7-cycle Pre-Amplification PCR fills gaps and amplifies both products; aliquots then go to ATAC Sample Index PCR and to cDNA amplification."),
    ]

GEX_BEAD = Construct([
    seg("partial TruSeq Read 1", il.TRUSEQ_READ1[-22:], "r1"),
    seg("cell barcode", "B" * 16, "cbc", placeholder=True, feature=CELL_BARCODE),
    seg("UMI", "U" * 12, "umi", placeholder=True, feature=UMI_FEATURE),
    seg("poly(dT)30", "T" * 30),seg("VN anchor","VN",placeholder=True,length_bp=2)], name="Single Cell Multiome Gel Beads A GEX oligo")
ATAC_BEAD = Construct([
    seg("P5", il.P5, "p5"),
    seg("10x cell barcode", "B" * 16, "cbc", placeholder=True, feature=CELL_BARCODE),
    seg("10x spacer", "CGCGTCTG")], name="Single Cell Multiome Gel Beads A ATAC oligo")
# CG000338 Rev G Appendix p. 80 prints these two gel-bead oligos verbatim.
assert GEX_BEAD.get("partial TruSeq Read 1").top == "CTACACGACGCTCTTCCGATCT"
# Rev G Sequencing tables: ATAC 50/8/24/49 (i5 24 = 8-nt spacer + 16-nt barcode;
# custom 8-dark-cycle recipe on NextSeq 500/550); GEX 28/10/10/90.
ATAC_READ_LENGTHS = {"Read 1": 50, "Index 1 (i7)": 8, "Index 2 (i5)": 24, "Read 2": 49}
GEX_READ_LENGTHS = {"Read 1": 28, "Index 1 (i7)": 10, "Index 2 (i5)": 10, "Read 2": 90}
# SI-PCR Primer B PN-2000128 is listed in Multiome Amp Kit A (Rev G p. 6) without bases;
# the same part number is printed as AATGATACGGCGACCACCGAGA in CG000168 Rev A Appendix
# (10x scATAC v1). Catalog-linked assignment, not a Rev G reprint.
SI_PCR_PRIMER_B = il.P5[:22]

def mixed_population(label, *populations):
    """Independent DNA/RNA populations share a container, never a covalent strand."""
    sc=Scene()
    for index,population in enumerate(populations):
        sc.strand(f"population {index + 1}",[role(population)],label=f"{label}: {population}")
    return sc


def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    nucleus=g.input("nucleus","nucleus containing accessible chromatin and poly(A) RNA",MaterialKind.BIOLOGICAL,state=mixed_population("intact nucleus","accessible chromatin","poly(A) RNA"),primary=True)
    transposome=g.input("transposome","ATAC Enzyme B PN-2000265 (loaded transposome)",MaterialKind.ADAPTER,inferred=True,source=SOURCE+" (transposon oligos not printed; Read 1N/Read 2N arms appear only in the final ATAC layout)")
    tn5=g.input("tn5","Tn5 transposase",MaterialKind.CATALYST)
    tagged=g.product("tagged_nucleus","nucleus retaining tagmented chromatin and RNA",MaterialKind.MOLECULAR,state=mixed_population("tagged nucleus","adapter-tagged accessible DNA","poly(A) RNA"))
    waste=g.product("tn5_waste","candidate same-adaptor tagmentation products",MaterialKind.MOLECULAR,inferred=True,source=SOURCE+" (candidate end classes; not described in CG000338)")
    g.react("tagment","Tagment accessible chromatin before partitioning",StepKind.TAGMENTATION,((nucleus,InputRole.TEMPLATE),(transposome,InputRole.STRUCTURAL),(tn5,InputRole.CATALYST)),(tagged,waste))
    g.mark_terminal(Terminal(waste.id,Disposition.UNRESOLVED,"Candidate insertion products stay in the nuclei through partitioning; no physical separation of them is described and their recovery is not assigned."))
    droplets=g.input("droplets","Chip J GEM partitioning condition",MaterialKind.CONDITION)
    occupied=g.product("occupied","occupied GEM products",MaterialKind.MOLECULAR,state=tagged.state)
    failed=g.product("failed","empty or unresolved GEM products",MaterialKind.MOLECULAR)
    g.react("partition","Partition tagged nuclei with barcoded gel beads",StepKind.PARTITION,((tagged,InputRole.TEMPLATE),(droplets,InputRole.CONDITION)),(occupied,failed),changes_molecule=False)
    g.mark_terminal(Terminal(failed.id,Disposition.UNRESOLVED,"Empty or multi-nucleus GEMs are pooled and recovered with the rest; exclusion is computational, not a physical discard."))
    rtprimer=g.input("rt_primer","GEX gel-bead oligo (Gel Beads A PN-2000261)",MaterialKind.OLIGO,state=GEX_BEAD,source=SOURCE+" (Rev G Appendix p. 80)")
    rt=g.input("rt","Barcoding Enzyme Mix PN-2000266 (RT activity)",MaterialKind.CATALYST)
    first=g.product("first_strand","GEM mixture with first-strand cell-barcoded cDNA",MaterialKind.MOLECULAR,state=mixed_population("GEM products","cell barcode + UMI + first-strand cDNA","adapter-tagged accessible DNA"))
    g.react("reverse_transcribe","Capture poly(A) RNA and reverse-transcribe in the GEM",StepKind.REVERSE_TRANSCRIPTION,((occupied,InputRole.TEMPLATE),(rtprimer,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(first,))
    tso=g.input("tso","Template Switch Oligo PN-3000228",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (in the GEM Master Mix; bases not printed in Rev G)")
    switched=g.product("switched","GEM mixture with template-switched cDNA",MaterialKind.MOLECULAR,state=mixed_population("GEM products","cell barcode + UMI + cDNA + TSO complement","adapter-tagged accessible DNA"),inferred=True,source=SOURCE+" (template-switch end not printed)")
    g.react("template_switch","INFERRED — copy the Template Switch Oligo at the cDNA 3′ end",StepKind.TEMPLATE_SWITCHING,((first,InputRole.TEMPLATE),(tso,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(switched,))
    atac_bead=g.input("atac_bead","ATAC gel-bead oligo (Gel Beads A PN-2000261)",MaterialKind.OLIGO,state=ATAC_BEAD,source=SOURCE+" (Rev G Appendix p. 80)")
    barcoded=g.product("barcoded_gem","GEM mixture with barcoded cDNA and barcoded ATAC fragments",MaterialKind.MOLECULAR,state=mixed_population("GEM products","cell-barcoded cDNA","P5 + cell barcode + spacer joined to s5 fragments"))
    g.react("barcode_atac","INFERRED — join the P5/barcode/spacer oligo to the s5 end of tagged fragments",StepKind.MAGIC,((switched,InputRole.TEMPLATE),(atac_bead,InputRole.STRUCTURAL)),(barcoded,),reason="CG000338 Rev G states only that GEM incubation produces 10x Barcoded DNA from transposed DNA; the final library places the CGCGTCTG spacer directly 5′ of s5, but the joining reagent (including any splint) and its geometry have not been reconstructed from checked sources; adjacency in the final layout is not a bond-formation mechanism.")
    quench=g.input("quench","Quenching Agent PN-2000269; Dynabeads MyOne SILANE and SPRIselect cleanup",MaterialKind.CONDITION)
    recovered=g.product("recovered","recovered barcoded GEM products",MaterialKind.MOLECULAR,state=barcoded.state)
    gem_waste=g.product("gem_waste","partitioning oil, leftover reagents and primers removed by bead cleanup",MaterialKind.WASTE)
    g.react("recover","Quench, break GEMs and purify (Dynabeads, SPRIselect)",StepKind.CLEANUP,((barcoded,InputRole.TEMPLATE),(quench,InputRole.CONDITION)),(recovered,gem_waste),changes_molecule=False)
    g.discard(gem_waste,"Rev G Step 3: silane beads purify products from leftover biochemical reagents and primers; supernatants are removed.")
    pre_f=g.input("preamp_primers","Pre-Amp Primers PN-2000271",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (bases not printed)")
    pre_r=g.input("preamp_primers_2","Pre-Amp Primers PN-2000271 (second end)",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (bases not printed; one mix, primer count unknown)")
    pre_pol=g.input("preamp_mix","Pre-Amp Mix PN-2000270",MaterialKind.CATALYST)
    preamp=g.product("preamp","pre-amplified ATAC and cDNA products",MaterialKind.MOLECULAR,state=recovered.state)
    g.react("preamp_pcr","Pre-Amplification PCR: 72 °C 5 min gap fill, then 7 cycles",StepKind.PCR,((recovered,InputRole.TEMPLATE),(pre_f,InputRole.STRUCTURAL),(pre_r,InputRole.STRUCTURAL),(pre_pol,InputRole.CATALYST)),(preamp,),geometry=PCRGeometry.unavailable(recovered.state,reason="Rev G names Pre-Amp Primers PN-2000271 but does not print them; the guide states that this PCR fills gaps and amplifies both modalities."))
    atac=g.product("atac_branch","40 µl pre-amplified aliquot for ATAC",MaterialKind.MOLECULAR,state=preamp.state)
    rna=g.product("rna_branch","35 µl pre-amplified aliquot for cDNA amplification",MaterialKind.MOLECULAR,state=preamp.state)
    g.react("split","Split the pre-amplified product into ATAC and cDNA aliquots",StepKind.SPLIT,((preamp,InputRole.TEMPLATE),),(atac,rna),changes_molecule=False)
    cf=g.input("cdna_forward","cDNA Primers PN-2000089 (forward)",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (bases not printed in Rev G)")
    cr=g.input("cdna_reverse","cDNA Primers PN-2000089 (reverse)",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (bases not printed in Rev G)")
    camp=g.input("cdna_amp","Amp Mix PN-2000047",MaterialKind.CATALYST)
    cdna=g.product("cdna","amplified cell-barcoded cDNA",MaterialKind.MOLECULAR,state=molecule("Multiome cDNA",(role("cell barcode"),role("UMI"),role("cDNA"))))
    g.react("cdna_amplification","Amplify barcoded cDNA",StepKind.PCR,((rna,InputRole.TEMPLATE),(cf,InputRole.STRUCTURAL),(cr,InputRole.STRUCTURAL),(camp,InputRole.CATALYST)),(cdna,),geometry=PCRGeometry.unavailable(rna.state,reason="cDNA Primers PN-2000089 are named but not printed in CG000338 Rev G."))
    frag=g.input("fragmentase","Fragmentation Enzyme/Buffer PN-2000090/2000091",MaterialKind.CATALYST)
    tailed=g.product("tailed","fragmented, end-repaired, A-tailed cDNA",MaterialKind.MOLECULAR)
    g.react("frag_er_at","Fragmentation, end repair and A-tailing (one reaction)",StepKind.FRAGMENTATION,((cdna,InputRole.TEMPLATE),(frag,InputRole.CATALYST)),(tailed,))
    adaptor=g.input("adaptor","Adaptor Oligos PN-2000094",MaterialKind.ADAPTER)
    ligase=g.input("ligase","DNA Ligase PN-220110",MaterialKind.CATALYST)
    adapted=g.product("adapted","adaptor-ligated cDNA fragments",MaterialKind.MOLECULAR,
        state=molecule("Multiome GEX pre-index input",(role("adaptor-ligated cDNA; exact physical ends unassigned"),)))
    g.react("ligate","Adaptor ligation",StepKind.LIGATION,((tailed,InputRole.TEMPLATE),(adaptor,InputRole.STRUCTURAL),(ligase,InputRole.CATALYST)),(adapted,),geometry=AdapterLigationGeometry.unavailable("CG000338 Rev G names Adaptor Oligos PN-2000094 but prints neither its strands nor its ligating ends; the final library alone does not establish reagent geometry."))
    reviewed=["nucleus","transposome","tn5","tagment","droplets","partition","rt_primer","rt","reverse_transcribe","tso","template_switch","atac_bead","barcode_atac","quench","recover","gem_waste","preamp_primers","preamp_primers_2","preamp_mix","preamp_pcr","split","cdna_forward","cdna_reverse","cdna_amp","cdna_amplification","fragmentase","frag_er_at","adaptor","ligase","ligate"]
    # Guide-specific index PCR declarations use the actual precursor states.
    # A complete final-library arm cannot supply an unassigned pre-PCR primer site.
    gex_f=g.input("gex_forward","Dual Index TT Set A i5 primer (physical sequence unassigned)",MaterialKind.OLIGO)
    gex_r=g.input("gex_reverse","Dual Index TT Set A i7 primer (physical sequence unassigned)",MaterialKind.OLIGO)
    gex_pol=g.input("gex_polymerase","Amp Mix PN-2000047/2000131",MaterialKind.CATALYST)
    gex_library=g.product("gex_library","gene-expression sequencing library",MaterialKind.LIBRARY,state=GEX)
    g.react("gex_pcr","Sample Index PCR with one Dual Index TT Set A pair (§7.5)",StepKind.PCR,
        ((adapted,InputRole.TEMPLATE),(gex_f,InputRole.STRUCTURAL),(gex_r,InputRole.STRUCTURAL),
         (gex_pol,InputRole.CATALYST)),(gex_library,),geometry=PCRGeometry.unavailable(adapted.state,
            reason="RevG §7.5 names Dual Index TT SetA; complete physical primers and the ligated precursor ends remain unassigned. Final-library TruSeq sites do not establish these input coordinates."))
    g.sequence(gex_library,SequencingPlan(library=GEX,primers=GEX_PRIMERS,
        acquired_roles=tuple(p.role for p in GEX_PRIMERS),read_lengths=GEX_READ_LENGTHS))
    atac_f=g.input("atac_forward","SI-PCR Primer B PN-2000128 (catalog-linked)",MaterialKind.OLIGO,
        state=Construct([seg("SI-PCR PrimerB",SI_PCR_PRIMER_B)]),
        source="CG000338 RevG PN-2000128; same-PN oligo printed in CG000168 RevA Appendix")
    atac_r=g.input("atac_reverse","Single Index Kit N SetA i7 primer mix (physical member unassigned)",MaterialKind.OLIGO)
    atac_pol=g.input("atac_polymerase","Amp Mix PN-2000047/2000103",MaterialKind.CATALYST)
    atac_library=g.product("atac_library","ATAC sequencing library",MaterialKind.LIBRARY,state=ATAC)
    g.react("atac_pcr","Sample Index PCR with SI-PCR PrimerB and Single Index N SetA (§5.1)",StepKind.PCR,
        ((atac,InputRole.TEMPLATE),(atac_f,InputRole.STRUCTURAL),(atac_r,InputRole.STRUCTURAL),
         (atac_pol,InputRole.CATALYST)),(atac_library,),geometry=PCRGeometry.unavailable(atac.state,
            reason="PrimerB is catalog-linked but the complete preamplified ATAC precursor and physical Single Index N primer mixture are unassigned. Final ATAC primer sites cannot stand in for that input."))
    g.sequence(atac_library,SequencingPlan(library=ATAC,primers=ATAC_PRIMERS,
        acquired_roles=tuple(p.role for p in ATAC_PRIMERS),read_lengths=ATAC_READ_LENGTHS,
        limitation="ATAC24cycle i5 recipe includes spacer8+barcode16; on NextSeq500/550 the guide requires8dark+16sequenced cycles, rather than24reported bases."))
    reviewed.extend(("gex_forward","gex_reverse","gex_polymerase","gex_pcr","gex_library",
                     "atac_forward","atac_reverse","atac_polymerase","atac_pcr","atac_library"))
    g.annotate(SourceDiscrepancy(
        "gex_library", "CG000338 Rev B omits four P5 bases in its gene-expression layout",
        (SourceClaim("CG000338 Rev B, Appendix p. 80, Gene Expression Library",
                     "The 5′ P5 arm is printed as AATGATACGGCGACCACCGAGATCT before the N10 i5 index."),
         SourceClaim("CG000338 Rev B, Appendix p. 80, ATAC gel bead and library; "
                     "Rev F/G, Appendix pp. 80–81, Gene Expression Library",
                     "The P5 arm is printed as AATGATACGGCGACCACCGAGATCTACAC before the index or barcode.")),
        "The gene-expression construct uses the full P5 arm ending ACAC, as printed in Rev F/G.",
        "Rev B's own ATAC layout and later guide revisions agree on the full P5 arm; "
        "the canonical Illumina P5 sequence agrees as well.",
        "Whether Rev B's shortened gene-expression line was a diagram error or represented "
        "a reagent variant is not stated by 10x."))
    g.declare_reviewed_scope(reviewed,SOURCE)
    return g.require_publishable()
