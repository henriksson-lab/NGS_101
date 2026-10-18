"""Guo2013 native scRRBS populations; no borrowed canonical adapter shell."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from chemdraw import Scene, Segment, feature
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_graph import (InputRole as I, MaterialKind as M, ProtocolGraph,
                            SequencingPlan, StepKind as K)
import seqprimers as sp

TITLE="scRRBS — single-cell reduced-representation bisulfite sequencing"
NOTES="01_scrrbs.html"
AUTHORITY="Guo et al., Genome Research2013, doi:10.1101/gr.161679.113; PMC3847781 Methods103–105 and data processing111–113"
SOURCE='<a href="https://doi.org/10.1101/gr.161679.113">Guo et al., Genome Research (2013)</a>, Methods: Construction of single-cell RRBS sequencing libraries.'
SUMMARY="Single-cell and unmethylated lambda DNA are MspI-digested, filled and A-tailed, ligated to premethylated indexed Illumina adaptors, bisulfite-converted and amplified with PfuTurbo Cx then Phusion. A 200–500-bp native-PAGE fraction enters paired-end HiSeq2000 sequencing."
CAVEAT="Native adaptor strands, protection positions, index member/width, library primers, first-cycle geometry and transfer between PCR rounds remain unassigned in this bounded scope. Original U-bearing DNA and methylated adaptor molecules are distinct from ordinary PCR daughters. No canonical TruSeq fork or six-base index is substituted; component-specific column and gel recovery remains unresolved."
FINAL_LIBRARY=None
SEQ_PRIMERS=()
SEQUENCING_ENDING="Source paired-end HiSeq2000 sequencing; native primer sites, read lengths and index acquisition remain unassigned."
READ_SPECS=(sp.read_spec("Read 1",None,"Bisulfite-derived mouse or lambda genomic DNA; methylation readout C/T versus G/A depends on strand","include"),sp.read_spec("Read 2",None,"Opposite end of same source paired-end library","include"),sp.read_spec("Sequencing output",None,"Native ligated index; width, read role and acquisition unassigned","undetermined"))
INDEX=feature("library_sample_index","sample_index","unknown",group="single_cell_library",note="Source indexed adaptor; native coordinate, member, width and acquired read role unassigned")


def role(name,*,indexed=False):
    return Segment(name,"["+name+"]",placeholder=True,feature=INDEX if indexed else None)


def population(stage,*,carrier=False):
    """Each row is an independent polymer class, not a covalent species mixture."""
    sc=Scene()
    for species in ("mouse genomic","unmethylated lambda spike-in"):
        if stage in {"released","digested","prepared"}:
            for strand in ("a","b"):
                parts=[role(species+" "+strand+" strand; "+("MspI C^CGG-bounded fragment / native5′CG overhang" if stage!="released" else "released genomic DNA"))]
                if stage=="prepared":
                    parts=[role(species+" "+strand+" fragment; native methylation and unmodified-dCTP fill-in distinguished"),Segment("extra3′A","A")]
                sc.strand(species+strand,parts,label=("source3′A-tailed strand; exact genome/end lengths unassigned" if stage=="prepared" else "independent strand population; partner relationship not base-resolved"))
        elif stage in {"ligated","converted"}:
            for strand in ("a","b"):
                text=(species+" "+strand+" native premethylated-indexed-adaptor-ligated DNA; strands/junctions/index coordinate unassigned")
                if stage=="converted":
                    text=(species+" "+strand+" converted original: unmethylatedC→U when converted; "+("original genomic5mC and " if species=="mouse genomic" else "")+"protected-adaptorC retained; native junctions unassigned")
                sc.strand(species+strand,[role(text,indexed=True)],label="independent converted original; former genomic partners need not complement" if stage=="converted" else "native adaptor-ligated original strand class")
        elif stage in {"first_pcr","second_template","second_pcr"}:
            sc.strand(species+" daughters",[role(species+" ordinary PCR DNA daughters; converted-site readoutT/A and retainedC/G; native ends/index coordinate unassigned",indexed=True)],label="ordinary DNA copies; methyl protection is not copied")
            sc.strand(species+" originals",[role(species+" possible carried U-bearing originals / protected methylated adaptors; recovery and second-round transfer unassigned",indexed=True)],label="independent original polymers; not chemically erased by PCR")
        else:
            raise ValueError(stage)
    if carrier:
        sc.strand("trna",[role("10ng tRNA carrier; retained amount and later recovery unassigned")],label="independent RNA polymer population")
    if stage in {"first_pcr","second_template","second_pcr"}:
        sc.note("mouse genomic daughters","Polymerase copying does not reproduce 5mC protection or U chemistry in ordinary DNA daughters. Original modified molecules can remain; their survival and native primer selection are unresolved.")
    if stage=="prepared":
        sc.note("mouse genomica","Source unmodified dCTP supplies fill-in C, distinct from native genomic methylation. Exact base-resolved duplex assembly is not supplied.")
    return sc


def protocol_graph():
    g=ProtocolGraph("scRRBS",authority=AUTHORITY)
    cell=g.input("cell","one picked mouse cell",M.BIOLOGICAL,primary=True)
    lysis=g.input("lysis","5µl lysis buffer:20mM Tris-EDTA pH8,20mM KCl,0.3%TritonX100;1mg/ml Qiagen protease;50°C3h,75°C30min",M.CONDITION)
    lam=g.input("lambda","60fg unmethylated lambda DNA",M.MOLECULAR)
    dna=g.product("released","released mouse DNA plus independent lambda spike-in",M.MOLECULAR,state=population("released"))
    g.react("lyse","Lyse cell with lambda spike-in;heat-inactivate protease",K.EXTRACTION,((cell,I.TEMPLATE),(lam,I.SPIKE_IN),(lysis,I.CONDITION)),(dna,))
    mspi=g.input("mspi","9U MspI;18µl reaction;37°C3h",M.CATALYST)
    digested=g.product("mspi_fragments","MspI-digested mouse and lambda fragments",M.MOLECULAR,state=population("digested"))
    g.react("mspi_digest","Digest at source C^CGG sites",K.DIGESTION,((dna,I.TEMPLATE),(mspi,I.CATALYST)),(digested,))
    klenow=g.input("klenow","5U Klenow exo-minus;20µl reaction",M.CATALYST)
    nucleotides=g.input("nucleotides","1mM dATP +0.1mM dGTP +0.1mM unmodified dCTP",M.MOLECULAR)
    prepared=g.product("prepared","filled fragments with extra3′A;ordinary fill-inC",M.MOLECULAR,state=population("prepared"))
    g.react("end_prepare","Fill recessed MspI ends and add source extra3′A",K.END_REPAIR,((digested,I.TEMPLATE),(klenow,I.CATALYST),(nucleotides,I.STRUCTURAL)),(prepared,))
    adapter=g.input("methylated_adapters","Illumina standard premethylated indexed adaptors;native physical strands/index member unassigned",M.ADAPTER)
    ligase=g.input("ligase","30U high-concentration T4 DNA ligase;25µl reaction",M.CATALYST)
    ligated=g.product("ligated","native premethylated-indexed-adaptor ligation population",M.MOLECULAR,state=population("ligated"))
    g.react("adapter_ligation","Ligate native source premethylated indexed adaptors",K.LIGATION,((prepared,I.TEMPLATE),(adapter,I.STRUCTURAL),(ligase,I.CATALYST)),(ligated,),geometry=AdapterLigationGeometry.unavailable("Methods105 identifies standard premethylated indexed Illumina adaptors, but native physical strands, exact junctions and ligation-facing end geometry remain unassigned in this scope. Canonical TruSeq is not substituted."))
    bis=g.input("bisulfite","MethyCode/Invitrogen kit,150µl;98°C10min denaturation then64°C2.5h conversion",M.CONDITION)
    converted=g.product("converted","separate U-bearing converted original strands;protected adaptorC retained",M.MOLECULAR,state=population("converted"))
    g.react("bisulfite_conversion","Denature and bisulfite-convert native ligated populations",K.CHEMICAL_CONVERSION,((ligated,I.TEMPLATE),(bis,I.CONDITION)),(converted,))
    trna=g.input("trna","10ng tRNA carrier",M.MOLECULAR)
    column=g.input("column","Zymo-Spin columns;on-column desulfonation;30µl eluate",M.SELECTION)
    purified=g.product("purified","30µl recovered converted DNA plus possible carried tRNA",M.MOLECULAR,state=population("converted",carrier=True))
    waste=g.product("column_waste","physically removed column flowthrough/washes;class composition unassigned",M.WASTE)
    g.react("column_recovery","Desulfonate and column-purify converted DNA with tRNA carrier",K.CLEANUP,((converted,I.TEMPLATE),(trna,I.STRUCTURAL),(column,I.SELECTION)),(purified,waste))
    g.discard(waste,"Physical column removal; no class-specific RNA/modified-DNA depletion is assumed.")
    f=g.input("forward","native first-round library forward primer;sequence/modifications unassigned",M.OLIGO)
    r=g.input("reverse","native first-round library reverse primer;sequence/modifications unassigned",M.OLIGO)
    cx=g.input("cx","1U PfuTurbo Cx;50µl firstPCR",M.CATALYST)
    first=g.product("first_pcr","firstPCR ordinary daughters plus possible originalDNA/tRNA",M.MOLECULAR,state=population("first_pcr",carrier=True))
    g.react("library_pcr","FirstPCR:95°C2min;25cycles95°C20s/60°C30s/72°C1min",K.PCR,((purified,I.TEMPLATE),(f,I.STRUCTURAL),(r,I.STRUCTURAL),(cx,I.CATALYST)),(first,),geometry=PCRGeometry.unavailable(purified.state,reason="Actual recovered U-bearing native ligation population plus possible tRNA. Exact primer species, first-cycle sites and native adaptor ends unassigned; completed canonical product is not an input."))
    template=g.product("second_template","source secondPCR input;transfer fraction and preparation unassigned",M.MOLECULAR,state=population("second_template",carrier=True))
    g.react("between_rounds","Transfer first-round material to secondPCR;scope unresolved",K.MAGIC,((first,I.TEMPLATE),),(template,),reason="Methods105 specifies two50µlPCR rounds but transfer fraction, dilution, cleanup and differential recovery of originalU-bearing DNA/tRNA versusordinarycopies remain unassigned. No specific physical aliquot or cleanup is invented.")
    sf=g.input("second_forward","native second-round forward primer;identity/relationship to first-round primer unassigned",M.OLIGO)
    sr=g.input("second_reverse","native second-round reverse primer;identity/relationship to first-round primer unassigned",M.OLIGO)
    phusion=g.input("phusion","1U Phusion high-fidelity;50µl secondPCR",M.CATALYST)
    second=g.product("amplified","secondPCR ordinary daughters;possible carried originalDNA/tRNA",M.MOLECULAR,state=population("second_pcr",carrier=True))
    g.react("second_pcr","SecondPCR:98°C2min;22cycles98°C10s/60°C30s/72°C1min",K.PCR,((template,I.TEMPLATE),(sf,I.STRUCTURAL),(sr,I.STRUCTURAL),(phusion,I.CATALYST)),(second,),geometry=PCRGeometry.unavailable(template.state,reason="Actual source transfer population. Native primer identity and first-cycle geometry unassigned; originalU-bearing molecules are not assumed converted into ordinaryDNA merely by transfer."))
    gel=g.input("gel","12%native polyacrylamide TBE gel;recover200–500bp",M.SELECTION)
    library=g.product("library","recovered200–500bp native indexed scRRBS library",M.LIBRARY,state=population("second_pcr",carrier=True))
    off=g.product("gel_waste","physical material outside recovered gel interval",M.WASTE)
    g.react("size_selection","Recover source200–500bp native-PAGE fraction",K.SELECTION,((second,I.TEMPLATE),(gel,I.SELECTION)),(library,off),changes_molecule=False)
    g.discard(off,"Physically excluded gel fraction; co-migrating originalDNA/tRNA/adaptor classes within the slice remain unresolved. Read filtering is computational, not physical removal.")
    g.sequence(library,SequencingPlan(acquired_roles=("Read 1","Read 2"),limitation=SEQUENCING_ENDING+" Lambda reads estimate conversion; sample-index member/width/readrole are unassigned. FragmentAnalyzer and qPCR analytical sampling are summarized without invented volumes."))
    g.declare_reviewed_scope(("cell","lysis","lambda","lyse","mspi","mspi_digest","klenow","nucleotides","end_prepare","methylated_adapters","ligase","adapter_ligation","bisulfite","bisulfite_conversion","trna","column","column_recovery","forward","reverse","cx","library_pcr","between_rounds","second_forward","second_reverse","phusion","second_pcr","gel","size_selection","library"),AUTHORITY)
    return g.require_publishable()
