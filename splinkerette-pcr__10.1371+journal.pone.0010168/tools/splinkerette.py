"""Potter–Luo splinkerette chemistry on an explicitly synthetic short locus.

Source oligos remain literal; invented template intervals are labelled synthetic.
The page and animation share the same nucleotide-resolved reaction trajectory.
"""
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT),str(ROOT/'lib')]
from lib.nucleotide_simulator import Simulator, Polymerase, Ligase, complement
from nucleotide_view import NucleotideState
from chemdraw import Construct, Segment, Scene, revcomp
from topology import HairpinTopology, FoldLayout, FoldStrand, FoldPairing
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import ProtocolGraph, Disposition, Terminal, MaterialKind as K, InputRole as R, StepKind as S
from protocol_views import assemble_page, WorkflowSection, OligoEntry, OligoSection
from page import head, info, table

PAGE_FILE='splinkerette.html'
TITLE='Splinkerette PCR — recover DNA beside an insertion'
DOI='10.1371/journal.pone.0010168'
SOURCE='Potter and Luo (2010), doi:'+DOI+', Introduction splinkerette design, Figure1B–C, Table1, and published Protocol S1 Steps2–7/Tables2–14.'
EXAMPLE='INFERRED — synthetic short locus carrying the published 3′ P-element primer sites in nested order. Inter-site distances, flanking DNA and restriction-site placement are invented; this is not a recovered Drosophila insertion.'
TOP='GATCCCACTAGTGTCGACACCAGTCTCTAATTTTTTTTTTCAAAAAAA'
BOT='CGAAGAGTAACCGTTGCTAGGAGAGACCGTGGCTGAATGAGACTGGTGTCGACACTAGTGG'
S1='CGAAGAGTAACCGTTGCTAGGAGAGACC'
S2='GTGGCTGAATGAGACTGGTGTCGAC'
T1='CACTCAGACTCAATACGACAC'
T2='GGATGTCTCTTGCCGAC'
T3='CGGGACCACCTTATG'
OLIGOS=(('SPLNK-GATC-TOP',TOP),('SPLNK-BOT',BOT),('SPLNK#1 (S1)',S1),
        ('SPLNK#2 (S2)',S2),('3′SPLNK#1 (T1)',T1),('3′SPLNK#2 (T2)',T2),('3′SPLNK-SEQ (T3)',T3))
# Figure1C specifies the fold; these spans are a complementary 7-bp stem
# with the remaining TTTC loop, not a thermodynamic folding prediction.
FOLD=HairpinTopology('SPLNK-GATC-TOP',TOP,30,41,7)
CORE=TOP[4:27]
LOCUS=T1+'ACGT'+T2+'AGCT'+T3+'TCGA'+'ACGTCATGCTACGACTGCTA'+'A'
PRECURSOR=LOCUS+'GATCC'+'TCGACTAG'
CUT_SITE=len(LOCUS)-1


def adapter_fold():
    layout=FOLD.folded_layout();t=layout.strands[FOLD.name]
    # The 23-base core pairs with the bottom oligo's 3′ end. Its 5′ arm
    # deliberately remains single-stranded, like Figure1C.
    x,y=t.points[4]
    bottom=FoldStrand('SPLNK-BOT',BOT,tuple((x+(len(BOT)-1-i)*12,y+24) for i in range(len(BOT))),
                      'var(--c-r2)',label='SPLNK-BOT · 5′ OH')
    return FoldLayout([t,bottom],[*layout.pairings,
        FoldPairing(t.name,4,bottom.name,len(BOT)-len(CORE),len(CORE),'adapter duplex core')])


def simulate(*,top_five='OH',bottom_five='OH',cofactors=('ATP',)):
    s=Simulator();r={'simulator':s,'stages':{}}
    top=s.add('genome_top',PRECURSOR,features=[(0,len(PRECURSOR),'synthetic')],source=EXAMPLE)
    bottom=s.add('genome_bottom',revcomp(PRECURSOR),features=[(0,len(PRECURSOR),'synthetic')],source=EXAMPLE)
    s.duplex('genome_top','genome_bottom',source=EXAMPLE)
    r.update(parent_top=top,parent_bottom=bottom,genome_initial=s.snapshot())
    r['digest']=s.restrict('genome_top',CUT_SITE,'AGATCC',1,5,source=SOURCE)
    gt=s.chain(top[0]);gb=s.chain(bottom[-1]);offtop=s.chain(top[len(LOCUS)]);offbottom=s.chain(bottom[0])
    r.update(genome_top=gt,genome_bottom=gb,off_top=offtop,off_bottom=offbottom,digested=s.snapshot())
    # Separation is noncovalent, not degradation of discarded flanks.
    bridge=[i for i in offtop if s.paired.get(i) in gb]
    r['bridge']=s.unpair(bridge,source=SOURCE+'; selected restriction fragments separate')
    r['fragment_ready']=s.snapshot()
    at=s.add('adapter_top',TOP,five=top_five,features=[(0,len(TOP),'adapter')],source=SOURCE)
    ab=s.add('adapter_bottom',BOT,five=bottom_five,features=[(0,len(BOT),'adapter')],source=SOURCE)
    r.update(adapter_top=at,adapter_bottom=ab,adapter_unfolded=s.snapshot())
    r['fold']=s.fold('adapter_top',30,41,7,source=SOURCE+' Figure1C; explicit selected fold')
    r['core']=s.bind_spans(at[4:27],tuple(reversed(ab[-23:])),source=SOURCE)
    r['adapter_ready']=s.snapshot()
    r['dock']=s.bind_spans(at[:4],tuple(reversed(gb[:4])),source=SOURCE)
    r['ligation_before']=s.snapshot()
    ligase=Ligase('T4 DNA ligase',('DNA',),('DNA',),('P',),('ATP',),('DNA',))
    r['seal']=s.ligate(ab[-1],gb[0],ligase,mode='nick',cofactors=cofactors,source=SOURCE)
    r['linked']=s.chain('adapter_bottom');r['ligated']=s.snapshot()
    # Keep the unsealed opposite junction as two real termini. An adapter
    # with a 5′ phosphate would be a different experimental reagent.
    if top_five!='OH' or bottom_five!='OH':
        raise ValueError('published splinkerette oligos must retain their unmodified 5′ hydroxyl ends')
    s.supply_nt('DNA','ACGT',source=SOURCE)
    pol=Polymerase('Phusion DNA polymerase',('DNA',),'DNA')
    def primer(name,seq):
        ids=s.add(name,seq,features=[(0,len(seq),'primer')],source=SOURCE)
        r[name+'_primer']=ids
        return ids
    primer('first_t1',T1)
    primer('first_s1',S1)
    r['first_before']=s.snapshot()
    r['melt1']=s.unpair(tuple(s.atoms),source=SOURCE+' Round1 thermal denaturation')
    r['first_open']=s.snapshot()
    # There must be no complementary S1 site before a transposon primer
    # has copied the ligated bottom oligo. Check this in construction too.
    try:s.anneal('first_s1','adapter_bottom',binding=len(S1),source=SOURCE)
    except ValueError as error:
        if '0 physical binding sites' not in str(error):raise
    else:raise ValueError('adapter primer site already exists before selective copying')
    r['t1_bind']=s.anneal('first_t1','adapter_bottom',binding=len(T1),source=SOURCE)
    r['first_bound']=s.snapshot()
    ops=s.extend_available(pol,['first_t1','adapter_bottom'],source=SOURCE)
    if len(ops)!=1:raise ValueError('selected first strand requires exactly one eligible extension')
    r['first_extend']=ops[0];r['first_top']=s.chain('first_t1');r['first_copied']=s.snapshot()
    r['melt2']=s.unpair(r['first_top'],source=SOURCE+' next Round1 denaturation')
    r['s1_before']=s.snapshot()
    r['s1_bind']=s.anneal('first_s1','first_t1',binding=len(S1),source=SOURCE)
    r['s1_bound']=s.snapshot()
    ops=s.extend_available(pol,['first_s1','first_t1'],source=SOURCE)
    if len(ops)!=1:raise ValueError('selected opposite-strand copy requires exactly one eligible extension')
    r['s1_extend']=ops[0];r['round1_bottom']=s.chain('first_s1');r['round1']=s.snapshot()
    # Nested PCR first copies the *physical* templates to their ends. The
    # first two products still carry an outer-primer boundary on one side.
    for cycle,top_ref,bottom_ref in [(1,'first_t1','first_s1'),(2,'nested_t1','nested_s1')]:
        f='nested_t'+str(cycle);q='nested_s'+str(cycle)
        primer(f,T2);primer(q,S2)
        r[f'nested{cycle}_before']=s.snapshot()
        r[f'nested{cycle}_melt']=s.unpair(s.chain(top_ref)+s.chain(bottom_ref),source=SOURCE+' nested PCR denaturation')
        r[f'nested{cycle}_open']=s.snapshot()
        r[f+'_bind']=s.anneal(f,bottom_ref,binding=len(T2),source=SOURCE)
        r[q+'_bind']=s.anneal(q,top_ref,binding=len(S2),source=SOURCE)
        r[f'nested{cycle}_bound']=s.snapshot()
        r[f'nested{cycle}_operations']=s.extend_available(pol,[f,q,top_ref,bottom_ref],source=SOURCE)
        if len(r[f'nested{cycle}_operations'])!=2:
            raise ValueError('both nested-PCR complexes must extend from their actual 3′ ends')
        r[f+'_ids']=s.chain(f);r[q+'_ids']=s.chain(q)
        r[f'nested{cycle}_after']=s.snapshot()
    r['final_top']=s.chain('nested_t2');r['final_bottom']=s.chain('nested_s2')
    r['final_before']=s.snapshot()
    r['final_melt']=s.unpair(r['final_top']+r['final_bottom'],source=SOURCE+' representative PCR denaturation')
    r['final_pair']=s.duplex('nested_t2','nested_s2',source=SOURCE+' endpoint strand reassociation')
    r['final']=s.snapshot()
    r['reagent_operations']={op.details['name']:op for op in s.history if op.kind=='reagent'}
    a=''.join(s.atoms[i].base for i in r['final_top']);b=''.join(s.atoms[i].base for i in r['final_bottom'])
    if a!=revcomp(b) or not a.startswith(T2) or not a.endswith(revcomp(S2)):
        raise ValueError('reaction-generated nested product has incorrect physical boundaries')
    if TOP in a or TOP in b:raise ValueError('hairpin oligo was incorrectly inherited into the amplicon')
    return r


def sequence(r,key):return ''.join(r['simulator'].atoms[i].base for i in r[key])


def state(r,snapshot,*tracks,note=''):
    return NucleotideState(r[snapshot],tuple((label,r[key]) for label,key in tracks),note)


def reference(r,key):
    value=sequence(r,key);parts=[]
    # Synthetic template remains inferred even after it has been copied.
    for j,i in enumerate(r[key]):
        atom=r['simulator'].atoms[i]
        parts.append(Segment('base'+str(j),atom.base,'r1' if 'adapter' in atom.features else None,
                             inferred='synthetic' in atom.features))
    return Construct(parts,'INFERRED synthetic splinkerette reference')


def protocol_graph():
    r=simulate();g=ProtocolGraph(TITLE,authority=SOURCE,required_output=Disposition.CONSTRUCT,
        representation='illustrative_example',example_description=EXAMPLE)
    dna=g.input('dna','INFERRED synthetic insertion-containing duplex',K.MOLECULAR,primary=True,
        state=state(r,'genome_initial',('synthetic top','parent_top'),('synthetic bottom','parent_bottom')),source=EXAMPLE,inferred=True)
    nuc=g.input('bstyi','BstYI; selected AGATCC site',K.CATALYST)
    fragment=g.product('fragment','INFERRED selected restriction fragment',K.MOLECULAR,
        state=state(r,'digested',('genomic top','genome_top'),('genomic bottom; 5′ GATC','genome_bottom')))
    other=g.product('other_fragment','Other cut fragment; retained in reaction mixture',K.MOLECULAR,
        state=state(r,'digested',('other top','off_top'),('other bottom','off_bottom')))
    g.react('digest','INFERRED locus — BstYI cleavage; 5′ GATC ends',S.DIGESTION,[(dna,R.TEMPLATE),(nuc,R.CATALYST)],[fragment,other])
    g.mark_terminal(Terminal(other.id,Disposition.UNRESOLVED,'Spectator fragment remains in the simulation; no purification or degradation is implied.'))
    a=g.input('top_oligo','SPLNK-GATC-TOP; 5′ OH',K.OLIGO,state=Construct([Segment('TOP',TOP)],'TOP'))
    b=g.input('bottom_oligo','SPLNK-BOT; 5′ OH',K.OLIGO,state=Construct([Segment('BOT',BOT)],'BOT'))
    cooling=g.input('cooling','95°C 3min, cool to room temperature',K.CONDITION)
    folded=g.product('folded','Terminal hairpin in TOP',K.MOLECULAR,
        state=NucleotideState(r['fold'].after,(('TOP; folded topology above, linear sequence map here',r['adapter_top']),)),
        chemical_state=FOLD)
    from protocol_graph import MolecularRelation as M
    g.react('fold','Fold TOP; retain its free 3′ end',S.FOLDING,[(a,R.TEMPLATE),(cooling,R.CONDITION)],[folded],relations=(M.HAIRPIN_STEM,))
    adapter=g.product('adapter','Asymmetric splinkerette; unfolded sequence map',K.ADAPTER,
        state=state(r,'adapter_ready',('TOP; terminal hairpin shown above','adapter_top'),('BOT; single-stranded 5′ arm','adapter_bottom'),note='Linear sequence map; the TOP hairpin is folded in the topology view.'))
    g.react('anneal_adapter','Anneal the 23-bp core; cool from 95°C',S.HYBRIDIZATION,[(folded,R.TEMPLATE),(b,R.STRUCTURAL)],[adapter])
    ligase=g.input('ligase','T4 DNA ligase + ATP-containing ligase buffer',K.CATALYST)
    geo=Scene();geo.strand('adapter',[Segment('GATC',TOP[:4]),Segment('core',CORE),Segment('hairpin tail',TOP[27:])],mod5='OH',mod3='OH')
    geo.anneal('bottom',[Segment('BOT arm',BOT[:-23]),Segment('core complement',BOT[-23:])],
        to='adapter',pair=('core complement','core'),unpaired=('BOT arm',),
        label='BOT; recessed 3′ OH',mod5='OH',mod3='OH')
    geo.overhang('adapter','GATC',polarity='5-prime',owner='adapter')
    geo.note('bottom','Unfolded sequence map. Only recessed BOT 3′ OH seals to genomic 5′ P; TOP 5′ OH stays unsealed.')
    geometry=AdapterLigationGeometry.known('cohesive',geo,statement='Source 5′ GATC cohesive end. Only BOT 3′ OH → restriction-derived genomic 5′ P is sealed; genomic 3′ OH → TOP 5′ OH remains unsealed. Both oligos are unmodified.')
    linked=g.product('linked','INFERRED locus — one sealed strand; opposite nick remains',K.MOLECULAR,
        state=state(r,'ligated',('BOT covalently joined to genomic bottom','linked'),('genomic top; free 3′ OH','genome_top'),('TOP; free 5′ OH','adapter_top'),note='One bond made; TOP and genomic top remain distinct physical chains. Hairpin shown separately above.'))
    g.react('ligate','INFERRED junction — ligate only BOT to genomic 5′ phosphate',S.LIGATION,
        [(fragment,R.TEMPLATE),(adapter,R.STRUCTURAL),(ligase,R.CATALYST)],[linked],geometry=geometry)
    pol=g.input('polymerase','Phusion DNA polymerase',K.CATALYST)
    nts=g.input('dntps','dATP, dCTP, dGTP, dTTP',K.MOLECULAR)
    primers={name:g.input(name,name.upper(),K.OLIGO,state=Construct([Segment(name,seq)],name))
             for name,seq in [('t1',T1),('s1',S1),('t2',T2),('s2',S2)]}
    heat=g.input('denaturation','98°C PCR denaturation',K.CONDITION)
    opened=g.product('opened','INFERRED separated ligated template; original TOP retained',K.MOLECULAR,
        state=state(r,'first_open',('ligated BOT/genomic template','linked'),
                    note='Thermal separation breaks base pairs, not the ligated backbone. TOP and the original genomic top remain separate strands in the inventory.'))
    original=g.product('original_strands','INFERRED original genomic top and unligated TOP oligo',K.MOLECULAR,
        state=state(r,'first_open',('original genomic top','genome_top'),('unligated TOP','adapter_top')))
    g.react('open_first','INFERRED locus — denature the duplex and hairpin',S.STRAND_SEPARATION,
        [(linked,R.TEMPLATE),(heat,R.CONDITION)],[opened,original])
    g.mark_terminal(Terminal(original.id,Disposition.UNRESOLVED,'Original strands remain in the mixture; no destruction, clearance or complete background-product model is claimed.'))
    bound=g.product('t1_bound','INFERRED T1 primer–template complex',K.MOLECULAR,
        state=state(r,'first_bound',('ligated template','linked'),('T1 primer','first_t1_primer')))
    g.react('bind_t1','INFERRED locus — T1 binds the known insertion end',S.HYBRIDIZATION,
        [(opened,R.TEMPLATE),(primers['t1'],R.STRUCTURAL)],[bound])
    first=g.product('first','INFERRED first strand creates the S1 complement',K.MOLECULAR,
        state=state(r,'first_copied',('T1 extension','first_top'),('ligated template','linked')))
    g.react('selective_copy','INFERRED locus — denature, bind T1, copy into BOT',S.EXTENSION,
        [(bound,R.TEMPLATE),(pol,R.CATALYST),(nts,R.STRUCTURAL)],[first])
    second_open=g.product('second_open','INFERRED first copy released from its template',K.MOLECULAR,
        state=state(r,'s1_before',('newly copied first strand','first_top'),
                    note='Original ligated template remains in the molecular inventory.'))
    parent=g.product('parent_template','INFERRED original ligated BOT/genomic template retained',K.MOLECULAR,
        state=state(r,'s1_before',('original ligated template','linked')))
    g.react('open_second','INFERRED locus — next thermal cycle separates the copy',S.STRAND_SEPARATION,
        [(first,R.TEMPLATE),(heat,R.CONDITION)],[second_open,parent])
    g.mark_terminal(Terminal(parent.id,Disposition.UNRESOLVED,'Retained parental template and later PCR descendants remain in the simulation; this graph expands one representative successful lineage.'))
    s1_bound=g.product('s1_bound','INFERRED S1 binds only the newly created site',K.MOLECULAR,
        state=state(r,'s1_bound',('first copy','first_top'),('S1 primer','first_s1_primer')))
    g.react('bind_s1','INFERRED locus — bind S1 to the copied BOT complement',S.HYBRIDIZATION,
        [(second_open,R.TEMPLATE),(primers['s1'],R.STRUCTURAL)],[s1_bound])
    double=g.product('round1','INFERRED round-1 duplex; T1/S1 boundaries',K.MOLECULAR,
        state=state(r,'round1',('T1-derived top','first_top'),('S1-derived bottom','round1_bottom')))
    g.react('adapter_copy','INFERRED locus — next cycle binds S1 to the newly made site',S.EXTENSION,
        [(s1_bound,R.TEMPLATE),(pol,R.CATALYST),(nts,R.STRUCTURAL)],[double])
    final=g.product('amplicon','INFERRED nested T2/S2 duplex',K.MOLECULAR,
        state=state(r,'final',('nested top','final_top'),('nested bottom','final_bottom')),chemical_state=reference(r,'final_top'))
    pcr=PCRGeometry.exact(reference(r,'first_top'),PCRPrimer('T2',T2,len(T2)),PCRPrimer('S2',S2,len(S2)),
        statement='Exact sites on the synthetic first-round reference. Both first nested extensions run to their actual template ends; subsequent template copying establishes T2/S2 boundaries.')
    g.react('nested_pcr','INFERRED locus — nested PCR with T2 and S2',S.PCR,
        [(double,R.TEMPLATE),(primers['t2'],R.STRUCTURAL),(primers['s2'],R.STRUCTURAL),(pol,R.CATALYST),(nts,R.STRUCTURAL)],
        [final],geometry=pcr)
    g.constructed(final,'Flanking-sequence amplicon. Source proceeds to cleanup and nested-primer Sanger sequencing; no NGS library or Illumina adapters are assigned.')
    g.declare_reviewed_scope(tuple(g.steps)+tuple(g.materials),SOURCE+' Scope ends at the amplified DNA; Step6 cleanup and Step7 Sanger are described on the page, not reconstructed as an NGS readout.')
    return g.require_publishable()


def render_page():
    r=simulate()
    oligos=OligoSection('Adapter and nested primers',tuple(OligoEntry(name,(Segment(name,seq),)) for name,seq in OLIGOS),
        note_html='All sequences 5′→3′, verbatim from Protocol S1 Table13. TOP and BOT are unmodified and unphosphorylated; do not order 5′ phosphates. T1/T2/T3 target the P-element 3′ end. For other constructs use the source’s matching primer set.',multiline=True)
    return assemble_page(protocol_graph(),(head(TITLE),'<div class="wrap"><h1>'+TITLE+'</h1>',
        '<p>Map an insertion by amplifying from its known transposon end into unknown genomic DNA. A restriction digest supplies a ligatable end; the asymmetric splinkerette supplies the second PCR handle. No circularization or cut inside the transposon is required.</p>',
        info('<a href="https://doi.org/'+DOI+'">Potter and Luo, PLoS ONE5:e10168 (2010)</a>; Figure1, Table1 and <a href="https://journals.plos.org/plosone/article/file?id=10.1371/journal.pone.0010168.s002&amp;type=supplementary">published Protocol S1</a>.'),
        info(EXAMPLE+' Published oligos are exact; highlighted inferred template bases and spacing are a teaching example, not a validated experimental insertion.'),oligos,
        '<h2>Why the adapter does not prime immediately</h2><p>S1 and S2 have the same polarity as the unpaired BOT arm. They cannot prime that arm. T1 first binds the known transposon sequence and copies the covalently joined genomic/BOT template. Only that newly synthesized strand contains the complementary adapter-primer sites. S1 can then copy back in the next cycle. TOP’s terminal hairpin reduces background end-repair priming; it is not a universal chemical block or a claim of zero background.</p>',
        adapter_fold().svg(930,300,'Splinkerette: 23-bp core, terminal hairpin and single-stranded BOT arm'),
        '<p>The folded view specifies a complementary 7-bp stem and TTTC loop within TOP. This is a geometric rendering of Figure1C, not a folding-energy prediction. The workflow below uses unfolded sequence maps; the simulator retains the actual self-pairs until thermal denaturation.</p>',
        '<h2>Design and experimental order</h2><ol><li>Choose a transposon end and its three source primers, ordered T1 → T2 → T3 toward the unknown flank. Confirm their sites on the actual construct. T1/S1 perform round1, T2/S2 are nested inside that product, and T3 reads from the transposon into the recovered flank.</li><li>Digest with a compatible 5′ GATC enzyme. The source recommends BstYI first; BglII, BamHI and BfuCI are alternatives. Use a site beyond the insertion junction; the enzyme need not cut inside the transposon.</li><li>Anneal TOP/BOT at95°C for3min and cool to room temperature. Ligate to the restriction digest using T4 DNA ligase for at least2h at room temperature. Only BOT3′ OH joins the genomic5′ phosphate. The opposite junction stays nicked because TOP5′ is OH.</li><li>Round1 uses S1 plus T1 and Phusion HF chemistry. Round2 uses no more than1µl of round1 with S2 plus T2. The first nested extensions reach the physical outer-template ends; later cycles establish both nested boundaries.</li><li>Gel-purify the desired product, or treat the remaining round2 reaction with Antarctic phosphatase and ExonucleaseI at37°C for2h, then80°C for15min. Sequence with T3 using Sanger chemistry and map the genomic sequence beyond the transposon boundary.</li></ol>',
        table(('Source recipe','Published setting'),(
            ('Genomic digest','35µl; ≥2h. BstYI60°C, then80°C20min; BglII37°C; BamHI37°C and purify; BfuCI37°C, then80°C20min.'),
            ('Ligation','35µl digest +2.5µl water +5µl10×ligase buffer +6µl annealed adapter +1.5µl T4 ligase(400U/µl);50µl total.'),
            ('Round1, P-element3′ end','25µl:10µl ligation,5µl5×HF buffer,0.5µl10mM dNTP,0.5µl each10µM S1/T1,0.25µl Phusion,8.25µl water. Source annealing61°C.'),
            ('Round1 cycling','98°C75s;2×(98°C20s,64°C15s);30×(98°C20s,61°C15s,72°C2min);72°C7min.'),
            ('Round2, P-element3′ end','25µl:≤1µl round1,5µl5×HF,0.5µl10mM dNTP,0.5µl each10µM S2/T2,0.25µl Phusion; water to25µl.98°C75s;30×(98°C20s,58°C15s,72°C90s);72°C7min.'))),
        '<p>Follow Protocol S1 for enzyme-specific digest handling, alternative transposon ends and controls. The model does not predict background, insertion multiplicity, fragment recovery, PCR efficiency or Sanger quality. Negative-control DNA and checking both insertion ends are recommended in the source.</p>',
        WorkflowSection('Step-by-step selected chemistry'),
        '<h2>Final molecule and read direction</h2><p>The product runs from the T2 boundary through the remaining known transposon segment, the unknown genomic flank, the restriction-derived junction and the copied S2 adapter boundary. TOP and its hairpin are not incorporated. T3 binds within the transposon and reads toward the unknown flank. The read may eventually reach adapter sequence; trim this before genomic mapping. No sample index, UMI or Illumina read layout is present.</p>',
        '<p><a class="research-notes" href="01_splinkerette.html">Source notes and simulation limits</a></p></div>'))


def selftest():
    from checks import Check, Source, have
    from storage import sources_root
    check=Check();r=simulate();s=r['simulator']
    check('opposite TOP/genomic junction remains chemically unsealed',r['ligated'].next.get(r['genome_top'][-1]),None)
    check('BOT is covalently joined to restriction-derived genomic 5′ phosphate',r['ligated'].next[r['adapter_bottom'][-1]],r['genome_bottom'][0])
    check('adapter primer complement is generated by actual first-strand copying',revcomp(S1) in sequence(r,'first_top'),True)
    check('nested products are complementary after physical-end copying',sequence(r,'final_top'),revcomp(sequence(r,'final_bottom')))
    check('first nested top retains an outer adapter boundary',len(r['nested_t1_ids'])>len(r['final_top']),True)
    check('first nested bottom retains an outer transposon boundary',len(r['nested_s1_ids'])>len(r['final_bottom']),True)
    for label,fn in [('unphosphorylated nick cannot be sealed',lambda:s.ligate(r['genome_top'][-1],r['adapter_top'][0],Ligase('T4',('DNA',),('DNA',),('P',),('ATP',),('DNA',)),mode='nick',cofactors=('ATP',),source=SOURCE)),
                     ('ligation requires ATP',lambda:simulate(cofactors=()))]:
        try:fn()
        except ValueError:failed=True
        else:failed=False
        check(label,failed,True)
    # Read the original published supplemental table, not scraper output.
    source=Source(sources_root()/'splinkerette-pcr__10.1371+journal.pone.0010168'/'journal.pone.0010168_supp_pone.0010168.s002.doc.txt',
                  'https://doi.org/'+DOI)
    if have(check,source,label='published Protocol S1'):
        text=source.path.read_text()
        for name,value in OLIGOS:check('source transcription '+name,value in text,True)
    from chemistry_validation import check_graph,Status
    g=protocol_graph()
    check('no contradictory chemistry declarations',[c for c in check_graph(g) if c.status is Status.CONTRADICTORY],[])
    check('explicit constructed-DNA endpoint',g.required_output,Disposition.CONSTRUCT)
    check('publication gate',g.publication_problems(),[])
    check.report()
