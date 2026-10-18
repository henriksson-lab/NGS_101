"""Macaulay et al. G&T-seq separation and branch model."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from chemdraw import Construct,Scene,Segment,feature,revcomp
from pcr import PCRGeometry
from protocol_graph import InputRole,MaterialKind,ProtocolGraph,SequencingPlan,StepKind,Highlight
from protocol_model import molecule,role
import seqprimers as sp

TITLE="G&T-seq — separate genome and transcriptome from one cell"
NOTES="01_g-t-seq.html"
SOURCE=('Defining source: <a href="https://doi.org/10.1038/nmeth.3370">Macaulay et al., <i>Nature Methods</i> (2015)</a>; '
        '<a href="https://zenodo.org/records/895966">accepted author manuscript with Online Methods (Zenodo 895966)</a>; '
        '<a href="https://www.pure.ed.ac.uk/ws/files/181359881/Separation_and_parallel_sequencing_of_the_genomes_and_transcriptomes_of_single_cells_using_G_T_seq.pdf">authors\' detailed protocol, <i>Nature Protocols</i> (2016)</a>.')
SUMMARY="A biotinylated oligo(dT) captures a cell's poly(A) RNA on streptavidin beads, physically separating it from genomic DNA; the two fractions then undergo Smart-seq2-like WTA and MDA/PicoPLEX WGA, respectively."
CAVEAT="G&T-seq separates one cell into genome and transcriptome branches. The original 2015 Online Methods print the capture, template-switch and PCR oligos. Their TSO is AAGCAGTGGTATCAACGCAGAGTACrGrG+G (accepted-manuscript Online Methods, extracted lines 336–337); the authors' 2016 detailed protocol prints AAGCAGTGGTATCAACGCAGAGTACATrGrG+G (oligo list lines 652–653), adding AT immediately before rGrG+G. The model uses the 2015 value for this defining version; the 2016 recipe is a separately scoped later version, not a correction or replacement of the defining2015 oligo. Downstream selected Nextera XT index members and WGA kit-specific oligos still need source checks. INFERRED — the template-switch panel draws untemplated C additions pairing with the TSO rGrG+G end; that mechanism and its length are not printed. The drawn capture oligo and TSO residues are the printed 2015 oligos."
CAPTURE_HANDLE="AAGCAGTGGTATCAACGCAGAGTAC"
CAPTURE_POLYT="T"*30
TSO_DNA="AAGCAGTGGTATCAACGCAGAGTAC"
TSO_WRITTEN=TSO_DNA+"rGrG+G"
TSO_2016="AAGCAGTGGTATCAACGCAGAGTACATrGrG+G"
ISPCR="AAGCAGTGGTATCAACGCAGAGT"
RNA=Construct([role("NexteraXT transcript-library physical adapter/index architecture unassigned",inferred=True,
    molecular_feature=feature("gt_rna_sample_index","sample_index","unknown",note="Libraries individually barcoded per2015Methods"))],name="G&T2015 transcriptome endpoint (bounded)")
DNA=Construct([role("NexteraXT WGA-library physical adapter/index architecture unassigned",inferred=True,
    molecular_feature=feature("gt_dna_sample_index","sample_index","unknown",note="Libraries individually barcoded per2015Methods"))],name="G&T2015 genome endpoint (bounded)")
RNA_P=DNA_P=()
FINAL_LIBRARIES=(("Transcriptome library",RNA,RNA_P,"2015 bead-captured RT/WTA and NexteraXT library; no later kit/index oligos substituted.","Historical NexteraXT physical read/index geometry remains unassigned."),
 ("Genome library",DNA,DNA_P,"2015 MDA/PicoPLEX WGA followed by NexteraXT library construction.","Historical NexteraXT physical read/index geometry remains unassigned."))
READ_SPECS={name:(sp.read_spec("Sequencing output",None,"amplified "+label+" insert; native priming and exact run layout unassigned","undetermined"),) for name,label in (("Transcriptome library","transcript"),("Genome library","genomic"))}
SEQUENCING_UNAVAILABLE="2015 Methods use NexteraXT/HiSeq2500 fast mode or MiSeq; selected physical index members/lengths, acquired reads and cycle counts are unassigned. Deep-WGA and targeted SureSelect runs are separate source branches."


def capture_scene(*,with_dna=False):
    rna=[Segment("mRNA body","X"*28,placeholder=True),Segment("poly(A)","A"*30)]
    oligo=[Segment("biotin capture handle",CAPTURE_HANDLE,"w1"),
           Segment("poly(dT)",CAPTURE_POLYT),Segment("3′ VN anchor","VN",placeholder=True)]
    sc=Scene(); sc.strand("mRNA",rna,label="cellular mRNA")
    sc.anneal("capture oligo",oligo,to="mRNA",pair=("poly(dT)","poly(A)"),
              label="biotin-oligo(dT)",unpaired=("3′ VN anchor",))
    if with_dna:sc.strand("genomic population",[role("genomic DNA in supernatant")],label="independent genomic DNA")
    sc.mark("capture oligo","biotin capture handle","streptavidin bead capture")
    return sc

def template_switch_local_scene():
    """Bead-captured RT and template switch drawn from the printed G&T oligos.

    The first strand is primed by the capture oligo itself. The untemplated C
    additions that pair with the TSO rGrG+G end are a mechanistic inference; their
    number is not printed by the source and the segment is marked inferred.
    """
    rna=[Segment("mRNA body","X"*30,placeholder=True),Segment("poly(A)","A"*30)]
    first=[Segment("capture handle",CAPTURE_HANDLE,"w1"),Segment("poly(dT)",CAPTURE_POLYT),
           Segment("3′ VN anchor","VN",placeholder=True),
           Segment("reverse-transcribed cDNA","X"*28,placeholder=True),
           Segment("untemplated C additions","CCC",inferred=True)]
    sc=Scene(); sc.strand("mRNA",rna,label="bead-captured mRNA")
    sc.anneal("first strand",first,to="mRNA",pair=("poly(dT)","poly(A)"),
              label="capture-oligo-primed cDNA",
              unpaired=("3′ VN anchor","reverse-transcribed cDNA","untemplated C additions"))
    sc.mark("first strand","capture handle","5′ biotin-TEG capture/RT oligo")
    sc.anneal("TSO",[Segment("TSO DNA handle",TSO_DNA,"w1"),Segment("3′ rGrG+G","GGG")],
              to="first strand",pair=("3′ rGrG+G","untemplated C additions"),label="2015 TSO")
    sc.mark("TSO","3′ rGrG+G","rGrG+G: two RNA G then one LNA G (printed residues)")
    return sc

def oligos():
    return [
        "Original 2015 Online Methods (accepted manuscript), extracted lines 319–342:",
        f"Biotinylated oligo-dT30VN: 5′-Biotin-TEG-{CAPTURE_HANDLE}{CAPTURE_POLYT}VN-3′",
        f"Template switching oligo: 5′-{TSO_WRITTEN}-3′",
        f"PCR primer: 5′-{ISPCR}-3′",
        f"2016 author protocol TSO: 5′-{TSO_2016}-3′ (extra AT before rGrG+G; model uses 2015 TSO; separately scoped later version)",
    ]

def first_parts(*,switched=False):
    parts=[Segment("capture handle",CAPTURE_HANDLE),Segment("poly(dT)30",CAPTURE_POLYT),
           Segment("VN anchor","VN",placeholder=True),role("reverse-transcribed transcript; variable length")]
    if switched:parts.extend((Segment("untemplated C additions","CCC",inferred=True),Segment("copied2015TSO",revcomp(TSO_DNA))))
    return parts


def lysate_scene():
    sc=Scene();sc.strand("genomic population",[role("genomic DNA")],label="independent genomic-DNA population")
    sc.strand("RNA population",[role("poly(A) RNA")],label="independent poly(A)-RNA population")
    return sc


def protocol_graph():
    source="Macaulay etal. defining2015 accepted manuscript OnlineMethods319–368; later2016protocol not substituted"
    g=ProtocolGraph(TITLE,authority=SOURCE)
    lysate=g.input("lysate","single-cell RLTPlus lysate with ERCC spike-in",MaterialKind.BIOLOGICAL,state=lysate_scene(),primary=True)
    capture=g.input("capture_oligo","2015 5′Biotin-TEG oligo(dT)30VN, bead-conjugated",MaterialKind.OLIGO,
        state=molecule("capture/RT oligo",(Segment("handle",CAPTURE_HANDLE),Segment("dT30",CAPTURE_POLYT),Segment("VN","VN",placeholder=True))),source=source)
    bead=g.input("streptavidin","Dynabeads carrying capture oligo",MaterialKind.SELECTION)
    hybrid=g.product("hybrid","bead-captured mRNA alongside genomic-DNA supernatant",MaterialKind.MOLECULAR,state=capture_scene(with_dna=True))
    g.react("hybridize","Capture poly(A)RNA20min with mixing",StepKind.HYBRIDIZATION,
        ((lysate,InputRole.TEMPLATE),(capture,InputRole.STRUCTURAL),(bead,InputRole.SELECTION)),(hybrid,),source=source)
    wash=g.input("wash_condition","fourTris/KCl/MgCl2/DTT/Tween/RNaseinhibitor washes; tipwash10µl",MaterialKind.CONDITION)
    rna=g.product("rna_fraction","bead-retained mRNA/capture-oligo hybrid",MaterialKind.MOLECULAR,state=capture_scene())
    dna=g.product("dna_fraction","pooled genomic-DNA supernatant plus bead/tipwashes",MaterialKind.MOLECULAR,state=molecule("DNA supernatant/washes",(role("genomic DNA; recovery/purity not quantitative"),)))
    g.react("separate","Transfer genomic-DNA supernatant and pool all fourwashes",StepKind.AFFINITY_CAPTURE,
        ((hybrid,InputRole.TEMPLATE),(bead,InputRole.SELECTION),(wash,InputRole.CONDITION)),(rna,dna),changes_molecule=False,source=source)
    rt=g.input("rt","SuperScriptII;42°C60min,50°C30min,60°C10min",MaterialKind.CATALYST,source=source)
    first=g.product("first_strand","single-strand oligo(dT)-primed cDNA",MaterialKind.MOLECULAR,state=molecule("first-strand cDNA",tuple(first_parts())))
    g.react("reverse_transcribe","Reverse-transcribe capturedmRNA using its bound oligo",StepKind.REVERSE_TRANSCRIPTION,
        ((rna,InputRole.TEMPLATE),(capture,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(first,),source=source)
    tso=g.input("tso","defining2015 TSO rGrG+G (twoRNA/LNA terminalG)",MaterialKind.OLIGO,
        state=molecule("2015TSO",(Segment("TSODNA",TSO_DNA),Segment("rGrG+G","GGG",note="twoRNA G then oneLNA G"))),source=source)
    switched=g.product("switched","INFERRED switched single-strand cDNA",MaterialKind.MOLECULAR,state=molecule("switchedssDNA",tuple(first_parts(switched=True))))
    g.react("template_switch","INFERRED — template-switch across untemplated C additions",StepKind.TEMPLATE_SWITCHING,
        ((first,InputRole.TEMPLATE),(tso,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(switched,),source=source)
    ispcr=g.input("rna_wta_primer","one physical2015 universalPCR primer species",MaterialKind.OLIGO,state=Construct([Segment("ISPCR",ISPCR)]),source=source)
    kapa=g.input("rna_wta_polymerase","KAPAHiFiHotStart;18×(98°C15s,67°C20s,72°C6min)",MaterialKind.CATALYST,source=source)
    parts=first_parts(switched=True)
    rwta=g.product("rna_wta","duplex amplified transcript cDNA",MaterialKind.MOLECULAR,
        state=Scene.duplex(parts,unpaired=tuple(s.name+"'" for s in parts if s.is_role_token() or s.placeholder)))
    g.react("amplify_rna","INFERRED copied ends — amplify with one universalprimer species",StepKind.PCR,
        ((switched,InputRole.TEMPLATE),(ispcr,InputRole.STRUCTURAL),(kapa,InputRole.CATALYST)),(rwta,),
        geometry=PCRGeometry.unavailable(switched.state,symmetric_species=True,reason="Actual incoming switchedcDNA is single stranded. Printed universalprimer is one physicalspecies; native first-cycle pairing remains a modelboundary, not a final-library template."),source=source)
    rbeads=g.input("rna_cleanup_beads","1:1AMPure,25µlEB eluate",MaterialKind.SELECTION)
    rclean=g.product("rna_clean","purified WTAcDNA",MaterialKind.MOLECULAR,state=rwta.state)
    rwaste=g.product("rna_cleanup_waste","WTAAMPure supernatant/washes",MaterialKind.WASTE)
    g.react("rna_cleanup","Clean WTA product with1:1AMPure",StepKind.CLEANUP,((rwta,InputRole.TEMPLATE),(rbeads,InputRole.SELECTION)),(rclean,rwaste),changes_molecule=False,source=source)
    g.discard(rwaste,"sourceAMPure supernatant/wash removal; exact component recovery unassigned")
    dbeads=g.input("dna_precipitation_beads","0.6×AMPure precipitation from pooledDNAfractions",MaterialKind.SELECTION)
    dclean=g.product("dna_precipitated","genomicDNA eluted into chosenWGA reaction",MaterialKind.MOLECULAR,state=dna.state)
    dwaste=g.product("dna_precipitation_waste","0.6×AMPure removed supernatant/washes",MaterialKind.WASTE)
    g.react("dna_precipitation","Precipitate retained genomicDNA with0.6×AMPure",StepKind.CLEANUP,((dna,InputRole.TEMPLATE),(dbeads,InputRole.SELECTION)),(dclean,dwaste),changes_molecule=False,source=source)
    g.discard(dwaste,"sourcebead precipitation removes liquidfractions; not rejection by sequence amplification")
    wga_kit=g.input("wga_kit","GenomiPhiV2 MDA or PicoPLEX alternative kit; component-specific primer/enzyme recipes unassigned",MaterialKind.CONDITION)
    wga=g.product("wga","whole-genome-amplified DNA from chosen kit",MaterialKind.MOLECULAR,state=molecule("WGA population",(role("genome-derived amplifiedDNA; kit-dependent branching/ends unassigned"),)))
    g.react("amplify_dna","Amplify genomicDNA using chosen source WGAkit",StepKind.MAGIC,((dclean,InputRole.TEMPLATE),(wga_kit,InputRole.CONDITION)),(wga,),
        reason="Source used both MDA and PicoPLEX alternatives. Native primer/enzyme and amplification architecture are not reconstructed as one interchangeable reaction.",source=source)
    wbeads=g.input("wga_cleanup_beads","1:1AMPure,25µlEB",MaterialKind.SELECTION)
    wclean=g.product("wga_clean","purified WGAproduct",MaterialKind.MOLECULAR,state=wga.state)
    wwaste=g.product("wga_cleanup_waste","WGAAMPure supernatant/washes",MaterialKind.WASTE)
    g.react("wga_cleanup","Purify WGA product with1:1AMPure",StepKind.CLEANUP,((wga,InputRole.TEMPLATE),(wbeads,InputRole.SELECTION)),(wclean,wwaste),changes_molecule=False,source=source)
    g.discard(wwaste,"sourcephysicalAMPure purification")
    for prefix,template,final in (("rna",rclean,RNA),("dna",wclean,DNA)):
        kit=g.input(prefix+"_library_kit","1–5ng amplified input; historicalNexteraXT per2015Methods",MaterialKind.CONDITION)
        library=g.product(prefix+"_library","individually indexed "+prefix+" library",MaterialKind.LIBRARY,state=final)
        g.react(prefix+"_library_prep","Prepare "+prefix+" library using source historicalNexteraXT",StepKind.MAGIC,
            ((template,InputRole.TEMPLATE),(kit,InputRole.CONDITION)),(library,),source=source,
            reason="The historical kit's physical tagmentation/gap-fill/indexPCR/cleanup recipe and selectedoligos are not reconstructed; generic endpoint or terminalPCR is not substituted.")
        g.sequence(library,SequencingPlan(library=final,acquisition="unspecified",limitation=SEQUENCING_UNAVAILABLE))
    g.annotate(Highlight("reverse_transcribe","RNA-template fate is unassigned","The output tracks first-strandcDNA. Subsequent physicalRNA recovery/degradation is not assigned or called a discard."))
    g.declare_reviewed_scope(("lysate","capture_oligo","streptavidin","hybridize","wash_condition","separate","rt","reverse_transcribe","tso","template_switch","rna_wta_primer","rna_wta_polymerase","amplify_rna","rna_cleanup_beads","rna_cleanup","dna_precipitation_beads","dna_precipitation","wga_kit","amplify_dna","wga_cleanup_beads","wga_cleanup","rna_library_kit","rna_library_prep","rna_library","dna_library_kit","dna_library_prep","dna_library"),source)
    return g.require_publishable()
