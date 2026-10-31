"""Published POINT nascent-RNA isolation and source-specific library branch prefixes."""
from dataclasses import dataclass, replace
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from protocol_graph import ProtocolGraph, MaterialKind, InputRole, StepKind, Terminal, Disposition
SOURCE='Sousa-Luís et al.2021 doi:10.1016/j.molcel.2021.02.034, VOR STAR Methods JATS487–503'
BRANCHES={'seq':('POINT-seq','NEBNext Ultra II Directional RNA E7760S'),
          '5':('POINT-5','SMARTer Stranded RNA-seq634836'),
          'nano':('POINT-nano','ONT direct cDNA SQK-DCS109')}

@dataclass(frozen=True)
class NascentRNA:
    carrier: str='RNA'
    artificial_poly_a: bool=False
    native_end_population: str='unassigned native RNA ends'
    def tail(self):
        if self.artificial_poly_a: raise ValueError('do not silently add the same artificial tail twice')
        return replace(self,artificial_poly_a=True)


def prefix_graph_for_branch(branch):
    if branch not in BRANCHES: raise ValueError('POINT branches are seq,5,nano')
    name,kit=BRANCHES[branch]
    g=ProtocolGraph(name+' — Pol-II-associated intact RNA',authority=SOURCE,representation='published_instance')
    sample=g.input('chromatin','nuclear/chromatin fraction from source cells;'+('4×10^7' if branch=='nano' else '10^7')+' cells',MaterialKind.BIOLOGICAL,primary=True)
    dnase=g.input('dnase','Turbo DNase in modified NUN/Empigen preparation',MaterialKind.CATALYST)
    digest=g.product('digested_chromatin','solubilized chromatin including RNA–PolII complexes and DNA digestion products',MaterialKind.MOLECULAR)
    g.react('digest_dna','Source modified NUN/Empigen chromatin preparation and DNase digestion',StepKind.DIGESTION,
            ((sample,InputRole.TEMPLATE),(dnase,InputRole.CATALYST)),(digest,))
    soluble=g.product('soluble_fraction','soluble digested chromatin supernatant',MaterialKind.MOLECULAR)
    pellet=g.product('centrifuge_pellet','material rejected by source13,000rpm clarification',MaterialKind.WASTE)
    g.react('clarify','Centrifuge digested chromatin and retain supernatant',StepKind.CLEANUP,
            ((digest,InputRole.TEMPLATE),),(soluble,pellet));g.discard(pellet,'Actual clarification pellet')
    antibody=g.input('polii_beads','anti-PolII antibody-conjugated beads',MaterialKind.SELECTION)
    captured=g.product('captured_complex','PolII-associated intact RNA complexes',MaterialKind.MOLECULAR)
    washes=g.product('ip_waste','immunoprecipitation flow-through and washes',MaterialKind.WASTE)
    g.react('capture','Anti-PolII IP4°C1h; six NET-2E washes',StepKind.AFFINITY_CAPTURE,
            ((soluble,InputRole.TEMPLATE),(antibody,InputRole.SELECTION)),(captured,washes));g.discard(washes,'Physical IP reject')
    state=NascentRNA()
    rna=g.product('purified_rna','Direct-zol-purified nascent RNA after further DNase treatment',MaterialKind.MOLECULAR,chemical_state=state)
    extract_waste=g.product('extraction_waste','extraction/column cleanup reject',MaterialKind.WASTE)
    g.react('extract','Source two Direct-zol purifications with intervening DNase treatment',StepKind.CLEANUP,
            ((captured,InputRole.TEMPLATE),),(rna,extract_waste));g.discard(extract_waste,'Physical extraction/column reject')
    endpoint=rna
    if branch=='seq':
        condition=g.input('fragmentation','historical E7760S fragmentation mixture;90°C10min',MaterialKind.CONDITION)
        endpoint=g.product('fragmented_rna','RNA-fragment mixture targeting150–200nt; all molecular fragments retained here',MaterialKind.MOLECULAR,
                           chemical_state=replace(state,native_end_population='mixture includes newly created fragment ends'))
        g.react('fragment','Fragment RNA before directional RNA library preparation',StepKind.FRAGMENTATION,
                ((rna,InputRole.TEMPLATE),(condition,InputRole.CONDITION)),(endpoint,))
    elif branch=='nano':
        pap=g.input('pap','E.coli poly(A) polymerase',MaterialKind.CATALYST)
        atp=g.input('atp','ATP required for poly(A) synthesis; source buffer recipe unreviewed',MaterialKind.MOLECULAR,inferred=True)
        tailed=g.product('poly_a_rna','RNA bearing in-vitro-added poly(A); tail length unassigned',MaterialKind.MOLECULAR,chemical_state=state.tail())
        g.react('tail','In-vitro polyadenylation before ONT direct-cDNA library prep',StepKind.TAILING,
                ((rna,InputRole.TEMPLATE),(pap,InputRole.CATALYST),(atp,InputRole.STRUCTURAL)),(tailed,))
        beads=g.input('size_beads','SPRISelect0.6×',MaterialKind.SELECTION)
        endpoint=g.product('selected_rna','source-selected polyadenylated RNA; removes RNA below500nt',MaterialKind.MOLECULAR,chemical_state=state.tail())
        small=g.product('small_rna_reject','RNA rejected by physical size selection',MaterialKind.WASTE)
        g.react('size_select','Source SPRISelect selection before SQK-DCS109',StepKind.SELECTION,
                ((tailed,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),(endpoint,small));g.discard(small,'Physical size-selection reject')
    g.mark_terminal(Terminal(endpoint.id,Disposition.UNRESOLVED,
        kit+' downstream molecular inputs/oligos/revision and exact read geometry remain unreviewed. '
        +('POINT-5 uses random-N6/template switching without RNA fragmentation.' if branch=='5' else '')
        +' No other branch recipe is substituted.'))
    g.declare_reviewed_scope(('chromatin','dnase','digest_dna','clarify','polii_beads','capture','extract'),SOURCE)
    return g

# Only the defining POINT-seq branch is completed here; other branches retain
# their independently source-qualified RNA prefixes and unreviewed kit boundaries.
PAGE_FILE='point.html'

def graph_for_branch(branch):
    return seq_graph() if branch=='seq' else prefix_graph_for_branch(branch)

def protocol_graph():return seq_graph()

def seq_graph():
    from protocol_graph import (MaterialKind as K,InputRole as R,StepKind as S,
                                SequencingPlan)
    g=ProtocolGraph('POINT-seq — source Pol-II nascent-RNA library',authority=SOURCE,
                    representation='published_instance')
    pellet=g.input('chromatin','Source-prepared nuclear/chromatin pellet from10^7 HeLa/HCT116 cells; initial mNET nuclear isolation is outside this selected boundary',K.BIOLOGICAL,primary=True)
    def material(k,label,kind=K.MOLECULAR):return g.input(k,label,kind)
    def op(k,label,x,extra=(),kind=S.MAGIC):
        out=g.product(k+'_product',label+' material; exact occupied molecules and recovery unassembled',K.MOLECULAR)
        g.react(k,label,kind,((x,R.TEMPLATE),)+tuple(extra),(out,),changes_molecule=kind is not S.GENERIC,reason='Source-positive local preparation; molecular input-to-output geometry unassembled.' if kind is S.MAGIC else '')
        return out
    def clean(k,label,x,supply):
        kept=g.product(k+'_kept','Actual retained '+label+' fraction; purity/recovery unassigned',K.MOLECULAR)
        waste=g.product(k+'_waste','Physical '+label+' rejects',K.WASTE)
        g.react(k,label,S.CLEANUP,((x,R.TEMPLATE),(supply,R.CONDITION if supply.kind is K.CONDITION else R.STRUCTURAL)),(kept,waste));g.discard(waste,'Source physical workup rejects; not a species-specific complete removal claim')
        return kept
    nun=material('nun_supplies','Actual NUN1 and modified NUN2 supplies: source Tris/HEPES/NaCl/EDTA/glycerol/MgCl2/NP40/urea/Empigen and inhibitor contributions; exact input volumes unassigned')
    x=op('nun_treat','Resuspend in NUN1; modified NUN2/Empigen treatment on ice10min with gentle inversion',pellet,((nun,R.STRUCTURAL),))
    pellet=g.product('washed_chromatin','Source400g30s retained chromatin pellet after one PBS wash',K.MOLECULAR)
    rejects=g.product('chromatin_wash_reject','Actual nonretained supernatant/PBS-wash fractions',K.WASTE)
    pbs=material('pbs','Actual source PBS wash contribution')
    g.react('chromatin_recover','400g30s chromatin collection and PBS wash',S.CLEANUP,((x,R.TEMPLATE),(pbs,R.STRUCTURAL)),(pellet,rejects));g.discard(rejects,'Actual source pellet/wash separation')
    dnase=material('dnase','Actual0.2U/ul Turbo DNase,37C15min',K.CATALYST)
    db=material('dnase_supplies','Actual source10mMTris7.5/400mMNaCl/100mMMnCl2/2UperulRiboLock contributions; no normalized recipe substitution')
    x=op('digest_dna','Source chromatin DNA digestion37C15min; DNA products remain in mixture',pellet,((dnase,R.CATALYST),(db,R.STRUCTURAL)),S.DIGESTION)
    x=clean('clarify','13000rpm10min soluble digested-chromatin supernatant',x,material('clarification_supply','Actual source clarification operation; no additional reagent dose',K.CONDITION))
    net=material('net_dilution','Actual NET2E tenfold dilution:50mMTris7.4/150mMNaCl/0.05percentNP40/3percentEmpigen stock')
    x=op('dilute','Dilute soluble fraction tenfold in ice-cold NET2E',x,((net,R.STRUCTURAL),),S.GENERIC)
    antibody=material('polii_beads','Actual10ug antiPolII antibody on200ul Dynabeads anti-mouseIgG for POINT-seq; ready-conjugated source reagent class',K.SELECTION)
    washes=material('net_washes','Actual six1ml ice-cold NET2E wash contributions')
    kept=g.product('captured','Source bead-associated PolII/nascent-RNA complexes',K.MOLECULAR);waste=g.product('ip_reject','Actual IP flow-through and six wash fractions',K.WASTE)
    g.react('capture','Source antiPolII capture4C1h and sixNET2Ewashes',S.AFFINITY_CAPTURE,((x,R.TEMPLATE),(antibody,R.SELECTION),(washes,R.STRUCTURAL)),(kept,waste));g.discard(waste,'Actual nonretained IP/wash fractions')
    x=clean('extract_first','First source Direct-zol/Trizol RNA purification',kept,material('directzol_first','Actual first Direct-zol/Trizol workup contributions; internal formulation unassembled'))
    fresh=material('fresh_dnase','Actual fresh intervening Turbo DNase treatment37C10min; amount unassigned',K.CATALYST)
    x=op('rna_dnase','Fresh source DNase treatment37C10min; products retained until workup',x,((fresh,R.CATALYST),),S.DIGESTION)
    x=clean('extract_second','Second source Direct-zol RNA purification',x,material('directzol_second','Actual second Direct-zol workup contributions'))
    qa=g.product('rna_qa','Actual source RNA size-analysis aliquot; amount unassigned',K.QC);stock=g.product('rna_stock','Isolated nascent-RNA stock after analytical allocation; amount unassigned',K.MOLECULAR)
    g.react('rna_qa_allocate','Separate actual TapeStation RNA analytical contribution',S.SPLIT,((x,R.TEMPLATE),),(qa,stock))
    qa_supply=material('tapestation','Actual4150TapeStation high-sensitivity RNA assay supplies')
    qa_result=g.product('rna_qa_result','Source RNA size-analysis readout; not library production',K.QC)
    g.react('rna_qa_assay','Source4150TapeStation RNA size analysis',S.MAGIC,((qa,R.TEMPLATE),(qa_supply,R.STRUCTURAL)),(qa_result,),reason='Actual analytical assay; detailed assay chemistry unassembled.');g.mark_terminal(Terminal(qa_result.id,Disposition.QC,'Source analytical fraction'))
    kit_input=g.product('kit_input','Actual selected source POINT-seq RNA input; amount unassigned',K.MOLECULAR);rest=g.product('rna_remainder','Unallocated isolated-RNA stock class; amount/presence/fate unassigned',K.MOLECULAR)
    g.react('kit_allocate','Select actual source library input without assuming full-stock transfer',S.SPLIT,((stock,R.TEMPLATE),),(kit_input,rest));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Unallocated RNA stock; no nonzero remainder or disposal assumed'))
    fragment=material('fragmentation','Actual historical E7760S RNA fragmentation contributions; internal formulation unassembled')
    fc=material('fragment_condition','Source90C10min physical fragmentation conditions',K.CONDITION)
    x=op('fragment','90C10min RNA fragmentation targeting150–200nt; physical fragments retained',kit_input,((fragment,R.STRUCTURAL),(fc,R.CONDITION)),S.FRAGMENTATION)
    kit=material('kit_reagents','Actual historical NEBNext UltraII Directional RNA Library Prep E7760S reagent collection; edition, components and stock packaging unassembled')
    library=g.product('dna_library','Actual source POINT-seq PCR-produced DNA-library class; whole adapters, RT/PCR orders and read landing unassembled',K.LIBRARY,chemical_state={'polymer':'DNA class of source-reported RNA-seq PCR library','unknown':'Historical kit inner chemistry/ends/ordered polymers/index/read geometry unassembled; not direct RNA sequencing'})
    kit_reject=g.product('kit_reject','Physical source library-preparation rejects, including source-reported loss of shortRNA; recovery distribution unassigned',K.WASTE)
    g.react('kit_prepare','Source-positive historical directional RNA-to-PCR-DNA library preparation',S.MAGIC,((x,R.TEMPLATE),(kit,R.STRUCTURAL)),(library,kit_reject),reason='Defining published Methods493 explicitly prepares a PCR library from fragmentedRNA using E7760S. Kit inner ligation/RT/PCR stages and complete polymers are not reconstructed from a modern recipe. Source loss of smallRNA occurs during library preparation, not fragmentation itself.');g.discard(kit_reject,'Actual source library-preparation rejection boundary; exact fraction/chemical composition unassigned')
    load=g.product('loaded_library','Actual selected NovaSeq6000 DNA-library loading contribution; amount unassigned',K.LIBRARY);remainder=g.product('library_stock','Untransferred prepared-library stock class; amount/presence/fate unassigned',K.LIBRARY)
    g.react('load_allocate','Select source NovaSeq6000 sequencing contribution',S.SPLIT,((library,R.TEMPLATE),),(load,remainder));g.mark_terminal(Terminal(remainder.id,Disposition.UNRESOLVED,'Library stock outside selected loading contribution'))
    g.sequence(load,SequencingPlan(acquisition='unspecified',limitation='Source POINT-seq NovaSeq6000. Exact acquired roles/cycles, index members and physical primer landing remain unassembled. Named historical kit produces source-positive DNA; no modern kit geometry is inherited.'))
    g.declare_reviewed_scope(SEQ_SCOPE,SOURCE)
    return g.require_publishable()

SEQ_SCOPE=tuple("""chromatin nun_supplies nun_treat chromatin_recover pbs dnase dnase_supplies digest_dna clarification_supply clarify net_dilution dilute polii_beads net_washes capture directzol_first extract_first fresh_dnase rna_dnase directzol_second extract_second rna_qa_allocate tapestation rna_qa_assay kit_allocate fragmentation fragment_condition fragment kit_reagents kit_prepare dna_library load_allocate loaded_library""".split())

def render_page():
    from page import head,info
    from protocol_views import assemble_page,WorkflowSection
    g=protocol_graph()
    return assemble_page(g,(head('POINT-seq'),'<main class="wrap" data-protocol-complete="true"><h1>POINT-seq</h1><p>Pol-II affinity capture enriches nascent RNA from source-prepared nuclear/chromatin material. The selected POINT-seq branch fragments that RNA and produces a directional PCR DNA library for NovaSeq6000.</p><p><a href="https://doi.org/10.1016/j.molcel.2021.02.034">Sousa-Luís et al. (2021)</a> · <a href="01_point.html">Research note</a></p>',info('Published POINT-seq preparation from the nuclear/chromatin pellet boundary. Initial crude-nuclear isolation follows the cited mNET-seq protocol and is not reconstructed here. Historical E7760S edition, internal ligation/RT/PCR chemistry, full library polymers and exact read/index geometry remain unassembled. A named kit handoff does not establish those molecular details.'),'<p>Fresh intervening DNase is distinct from the initial chromatin-digestion enzyme. Physical purification rejects, TapeStation analytical material and unassigned RNA/library stocks remain separate. Small-RNA loss is reported during library preparation, rather than declared to occur by RNA fragmentation alone.</p>','<h2>Other POINT branches</h2><p>The defining article also reports unfragmented template-switch POINT-5 and polyadenylated direct-cDNA POINT-nano. Their source-local RNA prefixes remain separate, unresolved models; this schematic completes only POINT-seq. It does not assign direct RNA sequencing or transfer either alternative kit recipe into this branch.</p>','<details><summary>Source preparation, analytical allocation and library handoff</summary>',WorkflowSection('POINT-seq source-selected preparation'),'</details></main>'))

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph()
    c('POINT-nano remains RNA before direct-cDNA kit',graph_for_branch('nano').materials['selected_rna'].chemical_state.carrier,'RNA')
    c('POINT-5 does not inherit fragmentation','fragment' in graph_for_branch('5').steps,False)
    c('Fresh second DNase is a distinct supply',g.steps['rna_dnase'].inputs[1].material,'fresh_dnase')
    c('RNA QA is not a production library',g.terminals['rna_qa_result'].disposition,Disposition.QC)
    c('Kit produces DNA library before sequencing',g.materials['dna_library'].chemical_state['polymer'],'DNA class of source-reported RNA-seq PCR library')
    c('Actual source gate',g.publication_problems(),[])
    original=Source(sources_root()/'point__10.1016+j.molcel.2021.02.034'/'j.molcel.2021.02.034_PMC8122139.xml.txt','Published defining Methods489–493')
    if have(c,original,label='POINT source branch and analytical scope'):
        text=' '.join(original.path.read_text().split())
        for label,needle in [('POINT-seq kit DNA production','NEBNext Ultra II Directional RNA library prep kit to prepare PCR library'),('Source RNA QA','4150 TapeStation'),('Source second DNase','one Turbo DNase treatment at 37°C for 10 min'),('Source fragmentation','90°C for 10 min')]:c(label,needle in text,True)
    c.report()

if __name__=='__main__': selftest()
