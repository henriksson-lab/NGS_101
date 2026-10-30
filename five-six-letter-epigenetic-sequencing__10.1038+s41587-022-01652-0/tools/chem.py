"""Source-owned private partial five/six-letter natural epigenetic sequencing.

One protocol, two explicit preparation branches. No MASD constructor inheritance.
"""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,revcomp
from illumina import TRUSEQ_READ1,NEBNEXT_ARM_READ2
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry

TITLE='Five- and six-letter epigenetic sequencing'
AUTHORITY='Published PMC10567558 main121/123/185–195/203; published SI Figure7 original p8/control and Tables1–3 identities; two originals. Selected80ng genomicDNA route; cfDNA/comparator preparations separate.'
ADAPTER1='ATGACGATGCGTTCGAGCATCGUCAUT'
Y_TOP=TRUSEQ_READ1
Y_BOTTOM=NEBNEXT_ARM_READ2+'A'
def unknown_scene():
    sc=Scene();sc.strand('unassembled',[Segment('unassembled','[source whole molecule/end geometry unassembled]',placeholder=True)]);return sc
def methylation(seq):return tuple((i,'5mC')for i,b in enumerate(seq)if b=='C')
def hairpin_fold():
    n=10;left=ADAPTER1[:n];right=ADAPTER1[16:26]
    if left!=revcomp(right.replace('U','T')):raise ValueError('Source adapter1 stem cannot pair')
    return {'source_ordering':ADAPTER1,'polymer':'DNA','modified_bases':methylation(ADAPTER1)+tuple((i,'deoxyU')for i,b in enumerate(ADAPTER1)if b=='U'),'stem_length':n,'loop':ADAPTER1[10:16],'three_prime_overhang':ADAPTER1[26:],'U_positions':tuple(i for i,b in enumerate(ADAPTER1)if b=='U'),'five_prime':'unassigned','inferred_fold':True}
def hairpin_scene(*,inferred=True):
    f=hairpin_fold();sc=Scene();sc.strand('left',[Segment('stem5',ADAPTER1[:10],inferred=inferred)])
    sc.anneal('right',[Segment('stem3',ADAPTER1[16:26],inferred=inferred),Segment('three prime T',ADAPTER1[26:],inferred=inferred)],to='left',pair=('stem3','stem5'))
    sc.overhang('right','three prime T',polarity='3-prime',owner='adapter',text='Computed source adapter1 3prime T interface')
    sc.note('left','INFERRED candidate fold of one covalent oligo: GTTCGA loop joins source arms. U sites lie in the paired stem, not the loop; post-USER junction/end paths unassembled. Full raw ordering/modifications retained in chemical state.');return sc
def y_scene():
    n=12
    if Y_TOP[-n-1:-1]!=revcomp(Y_BOTTOM[:n]):raise ValueError('Source Y terminal stem cannot pair')
    sc=Scene();sc.strand('top',[Segment('Y top fork',Y_TOP[:-n-1]),Segment('Y top stem',Y_TOP[-n-1:-1]),Segment('Y three prime T',Y_TOP[-1:])])
    sc.anneal('bottom',[Segment('Y bottom stem',Y_BOTTOM[:n]),Segment('Y bottom fork',Y_BOTTOM[n:])],to='top',pair=('Y bottom stem','Y top stem'),unpaired=('Y bottom fork',))
    sc.overhang('top','Y three prime T',polarity='3-prime',owner='adapter',text='Computed source adapter2 3prime T interface');return sc
def local_geometry(adapter):
    sc=hairpin_scene()if adapter==1 else y_scene()
    return AdapterLigationGeometry.known('cohesive',sc,unknown_scene(),statement='Only source adapter local stem/3primeT is computed; source endprep/polymerase substrate and full sealed junction paths unassembled, not confirmed missing. Adapter1 candidate folding is inferred; no proved substrate compatibility or ligation yield.')

def graph_for_letters(letters=5):
    if letters not in(5,6):raise ValueError('Select source five- or six-letter branch')
    g=ProtocolGraph(f'{TITLE}: {letters}-letter80nggDNA branch',authority=AUTHORITY,representation='published_instance')
    x=g.input('input','Source80ng genomicDNA; no imposed M.SssI methylation of study sample',K.MOLECULAR,primary=True)
    def stock(key,label,kind=K.CATALYST,role=R.CATALYST,chem=None,state=None):return(g.input(key,label,kind,chemical_state=chem,state=state),role)
    def op(key,label,kind,inputs=(),product='Source reaction mixture, residual substrates/free reagents/off-products retained until workup',geometry=None):
        nonlocal x
        kw={}
        if kind==S.GENERIC:kw['changes_molecule']=False
        if kind==S.MAGIC:kw['reason']='Source-specific physical transition/bookkeeping remains unassembled in this selected model.'
        if geometry is not None:kw['geometry']=geometry
        y=g.product(key+'_product',product,K.MOLECULAR);g.react(key,label,kind,((x,R.TEMPLATE),)+tuple(inputs),(y,),**kw);x=y
    def cond(key,label):return stock(key,label,K.CONDITION,R.CONDITION)
    def clean(key,label):
        nonlocal x
        c=g.input(key+'_condition',label,K.CONDITION);a=g.product(key+'_kept','Source physically recovered DNA; recovery/carryover unassigned',K.MOLECULAR);w=g.product(key+'_waste','Source separated SPRI/workup fractions',K.WASTE)
        g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(c,R.CONDITION)),(a,w));g.discard(w,'Physical source workup; no complete purity/recovery inferred.');x=a
    def split(key,label,rest,disp=Disposition.UNRESOLVED):
        nonlocal x
        a=g.product(key+'_used',label,K.MOLECULAR);b=g.product(key+'_rest',rest,K.QC if disp==Disposition.QC else K.MOLECULAR);g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(a,b));g.mark_terminal(Terminal(b.id,disp,rest));x=a
    op('sonicate','Source sonicate extractedgDNA to~250bp; cfDNA unfragmented route separate',S.FRAGMENTATION,(cond('sonicate_condition','Source mechanical sonication'),))
    spike=stock('spike_stock','Source1ul0.5ngul ground-truth spike stock for>=80ng input; selected composition/proportions of controls unassigned',K.MOLECULAR,R.SPIKE_IN,{'scope':'Source185–193 describes synthetic5hmC/pUC19/lambda controls; no all-controls-in-every-sample assertion'})
    op('spike','Source add ground-truth spike before endprep',S.GENERIC,(spike,))
    op('endprep','Source3.5ulEndPrepbuffer+1.5ulenzymemix E7647;20C30min65C30min',S.END_REPAIR,(stock('endprep_enz','Source NEB E7647 EndPrep enzyme mix'),cond('endprep_condition','Source stated buffer/time; detailed historical chemistry/ends unassembled')))
    hp=stock('adapter1','Source3.75ul adapter1 hairpin stock; allC5mC and deoxyU',K.ADAPTER,R.STRUCTURAL,hairpin_fold(),hairpin_scene())
    op('ligate1','Source3.75ul adapter1+0.5ulenhancer+15ulligationmastermix E7647,20C15min',S.LIGATION,(hp,stock('ligase1','Source first fresh E7647 ligation master mix'),cond('ligation1_condition','Source enhancer/ligation temperature')),geometry=local_geometry(1))
    clean('ligate1_clean','Source SPRI after first ligation; elute23.75ul')
    op('user','Source23.75ulDNA+3ul10xrCutSmart+3.25ulUSER37C30min; source deoxyU in paired stem, not generic loop-opening recipe',S.DIGESTION,(stock('user_enzyme','Source USER'),cond('user_condition','Source rCutSmart buffer')),product='Source USER-processed hairpin-linked DNA classes and cleavage fragments; exact strand separation/new end/genomic junction paths unassembled, not assumed ideal full opening')
    clean('user_clean','Source SPRI after USER; elute12ul')
    copying=[stock('klenow','Source2ul Klenow exo-minus'),stock('pnk','Source2ul T4PNK'),stock('copy_nts','Source2ul10mM each dNTP',K.MOLECULAR,R.STRUCTURAL),cond('copy_condition','Source12ulDNA+2ul10xNEBbuffer4,37C30min; source3prime-end copying relation unassembled, no external primer addition')]
    if letters==6:copying.extend([stock('early_bgt','Source extra1ulT4BGT added in same PNK/Klenow reaction'),stock('early_udp','Source extra2.5ulUDPglucose in same reaction',K.MOLECULAR,R.STRUCTURAL)])
    op('copy','Source simultaneous PNK/Klenow'+('/earlyBGT'if letters==6 else '')+' treatment; original plus newly synthesized unmodified copy linked by hairpin',S.EXTENSION,copying,product='Source original/template-copy linked DNA classes; no perfect copy completion/sequence/ligation topology imposed')
    op('reanneal','Source95C2min, cooling-0.1C/sec reanneal',S.GENERIC,(cond('cool_condition','Source controlled cooling/annealing'),))
    ya=stock('adapter2','Source2.5ul TOTAL ready Y-adapter stock; not2.5ul of each free strand',K.ADAPTER,R.STRUCTURAL,{'source_top':Y_TOP,'source_bottom':Y_BOTTOM,'top_modified_bases':methylation(Y_TOP),'bottom_modified_bases':methylation(Y_BOTTOM),'top_phosphorothioate_after':(len(Y_TOP)-2,),'five_prime':'unassigned'},y_scene())
    op('ligate2','Source2.5ul ready adapter2+0.5ulenhancer+15ulligationmix+12ulwater20C15min',S.LIGATION,(ya,stock('ligase2','Source fresh second ligation master mix'),cond('ligation2_condition','Source stated ligation programme')),geometry=local_geometry(2))
    clean('ligate2_clean','Source SPRI after second ligation; five-letter elute30ul; six-letter source uses15ulDNA forDNMT5')
    if letters==6:
        split('dnmt_aliquot','Source15ulDNA into DNMT5 reaction','Untransferred second-ligation eluate; exact amount/fate unassigned')
        op('dnmt5','Source15ulDNA+3.95water+0.75ul100mMATP+1.5ul100mMMgCl2+0.4ulSAM+6ul5xbuffer+2.4ulDNMT5,37C3h',S.CHEMICAL_CONVERSION,(stock('dnmt5_enzyme','Source DNMT5'),stock('dnmt5_sam','Source SAM methyl donor',K.MOLECULAR,R.STRUCTURAL),cond('dnmt5_condition','Source5xbuffer250mMTris8/10mMDTT and ATP/Mg')),product='Source DNMT5 methyl-copying mixture after early5hmC glucosylation; unconverted species/free reagents retained')
        op('dnmt_bridge','Source post-DNMT5 preparation to30ul protection input; transfer/cleanup bookkeeping unassembled, no invented SPRI or perfect enzyme removal',S.MAGIC)
    op('protect','Source30ulDNA+10ulTET2 5xsupplement+1ulUDPglucose+1ulT4BGT+1ul100mMDTT+2ulTET2+5ulfresh1:1250dilution500mMFeSO4,37C60min',S.CHEMICAL_CONVERSION,(stock('tet2','Source TET2'),stock('protect_bgt','Source fresh protection-stage T4BGT'),stock('protect_udp','Source protection-stage UDPglucose',K.MOLECULAR,R.STRUCTURAL),cond('protection_condition','Source5xsupplement10mMalphaKG/0.25MTris8/10mMATP; fresh iron dilution retained literally')),product='Source oxidized/glucosylated protected cytosine derivatives and residual substrates; no complete protection yield assumed')
    clean('protect_clean','Source SPRI after protection; elute31ul')
    op('deaminate','Source31ulDNA+12.75water+17.5ul4xbuffer+1.75ul100mMATP+3.5ul100mMMgCl2+2ulA3A+2.5ulUvrD,37C1.5h',S.CHEMICAL_CONVERSION,(stock('apobec','Source APOBEC-A3A'),stock('uvrd','Source UvrD helicase'),cond('deam_condition','Source4xbuffer200mMBisTris6.1/0.4percentTween stock')),product='Source C-to-U and protected original/copy strand classes with residual unprocessed DNA; chemicalU is distinct from PCR-reportedT')
    clean('deam_clean','Source SPRI after deamination; elute20ul')
    pair=g.input('udi_pair','Source5ul Abclonal UDI primer-pair stock RK21624SetA; whole orderings/selected indexes unassigned',K.OLIGO,chemical_state={'identifier_role':'unique dual sample indexes','sequence_status':'unassembled; no invented canonical full primers'})
    forward=g.product('udi_forward','Forward primer component of same5ul supplied UDI pair',K.OLIGO);reverse=g.product('udi_reverse','Reverse primer component of same5ul supplied UDI pair',K.OLIGO)
    g.react('udi_components','Represent two co-supplied primer components of ONE stock, no physical split or extra additions',S.GENERIC,((pair,R.TEMPLATE),),(forward,reverse),changes_molecule=False)
    op('pcr','Source20ulDNA+5ulUDIprimerpair+25ul2xKAPAHiFiU+KK2802;98C30s,6x98C10s/62C30s/65C60s,65C5min for80nggDNA',S.PCR,((forward,R.STRUCTURAL),(reverse,R.STRUCTURAL),stock('pcr_kapa','Source KAPA HiFi UracilPlus master mix')),product='Source amplified original/copy-encoded products and off-products; source six-cycle80ng route, cfDNA7/8cycles separate',geometry=PCRGeometry.unavailable(unknown_scene(),reason='Source UDI primers/hairpin-copy whole geometry unassembled; no source absence claim.'))
    clean('pcr_clean','Source SPRI after PCR; elute15ul10mMTris8')
    split('tapestation','Production library retained','Source TapeStationD5000 QC aliquot, allocation unassigned',Disposition.QC)
    op('dilute','Source dilute library1.2–1.8nM forNovaSeqS4; pool allocation unassigned',S.GENERIC,(cond('loading_condition','Source NovaSeqS4 loading preparation'),))
    op('phix','Source add8percentPhiX to all converted sequencing libraries',S.GENERIC,(stock('phix_stock','Source PhiX stock',K.MOLECULAR,R.SPIKE_IN),))
    lib=g.product('library','Source converted NovaSeqS4 sequencing preparation',K.LIBRARY);g.react('load','Source sequencing loading',S.GENERIC,((x,R.TEMPLATE),),(lib,),changes_molecule=False)
    g.sequence(lib,SequencingPlan(limitation='Source111/8/8/111 or151/8/8/151 run alternatives; selected run unassigned. R1 P7-primed original; R2 P5-primed newly synthesized copy. Whole primer/cycle geometry unassembled. ComputationalNmask/filtering is not physical DNA loss; chemicalU vs reportedT distinct.'))
    g.add_source_item(SourceItem('pending','Full postUSER/copy/junction/read geometry, exact historical AbclonalUDI orderings and six-letter postDNMT5 transfer bookkeeping unassembled.','Source-positive selected80ng route only; no accepted missing-information or full linked-source audit.'))
    return g
def protocol_graph():return graph_for_letters(5)

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();a=graph_for_letters(5);b=graph_for_letters(6);f=hairpin_fold()
    c('USER sites are in paired source stem',all(16<=p<26 for p in f['U_positions']),True)
    c('Candidate hairpin overhang3primeT',f['three_prime_overhang'],'T')
    c('No imposed study-sample M.SssI',any('msss'in key for key in a.materials),False)
    c('USER consumes physically cleaned ligation',a.steps['user'].inputs[0].material,'ligate1_clean_kept')
    c('Copy consumes physically cleaned USER DNA',a.steps['copy'].inputs[0].material,'user_clean_kept')
    c('No external copy primer or SELF_PRIMER',[p.material for p in a.steps['copy'].inputs if p.role in{R.STRUCTURAL,R.SELF_PRIMER}],['copy_nts'])
    c('Five-letter lacks early BGT','early_bgt'in a.materials,False)
    c('Six-letter earlyBGT in samecopyreaction',[p.material for p in b.steps['copy'].inputs if p.role==R.CATALYST],['klenow','pnk','early_bgt'])
    c('Six-letter DNMT after second ligation cleanup',b.steps['dnmt_aliquot'].inputs[0].material,'ligate2_clean_kept')
    c('Unassembled DNMT transfer precedes protection',b.steps['protect'].inputs[0].material,'dnmt_bridge_product')
    c('Ready Y stock consumed once',sum(p.material=='adapter2'for st in b.steps.values()for p in st.inputs),1)
    c('Ready Y retains both full strands',(b.materials['adapter2'].chemical_state['source_top'],b.materials['adapter2'].chemical_state['source_bottom']),(Y_TOP,Y_BOTTOM))
    c('Top terminal PS retained',b.materials['adapter2'].chemical_state['top_phosphorothioate_after'],(len(Y_TOP)-2,))
    c('Protection BGT is separate fresh addition',[(st.id,p.material)for st in b.steps.values()for p in st.inputs if p.material in{'early_bgt','protect_bgt'}],[('copy','early_bgt'),('protect','protect_bgt')])
    c('Physical clean before PCR',b.steps['pcr'].inputs[0].material,'deam_clean_kept')
    c('PCR consumes same-stock produced UDI components',[p.material for p in b.steps['pcr'].inputs if p.role==R.STRUCTURAL],['udi_forward','udi_reverse'])
    c('PhiX added after source loading dilution',a.steps['phix'].inputs[0].material,'dilute_product')
    path=sources_root()/Path(__file__).resolve().parents[1].name/'s41587-022-01652-0_PMC10567558.xml.txt'
    if have(c,Source(path,'https://doi.org/10.1038/s41587-022-01652-0')):
        compact=''.join(path.read_text().split())
        c('Whole source adapter1 ordering',ADAPTER1 in compact,True)
        c('Whole source Y top with terminalPS',Y_TOP[:-1]+'*'+Y_TOP[-1]in compact,True)
        c('Whole source Y bottom ordering',Y_BOTTOM in compact,True)
    c('Both branch graphs stay private/incomplete',all(g.publication_problems()for g in(a,b)),True)
    c.report()
if __name__=='__main__':selftest()
