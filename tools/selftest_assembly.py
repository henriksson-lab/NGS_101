"""Negative controls for explicit assembly relations and illustrative provenance."""
from dataclasses import replace
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'lib'))
from checks import Check
from chemdraw import Construct,Scene,Segment,revcomp,feature
from cohesive_assembly import CohesiveFragment,CircularCohesiveAssembly
from assembly_geometry import StrandRetention,BridgedAssembly,OverlapExtension,NucleotideBoundary,pcr_product,DuplexEndState,ActivatedEndTransfer,LinearFivePrimeAssembly,CircularOverlapReference,DuplexPartition
from pcr import PCRGeometry,PCRPrimer
from protocol_graph import ProtocolGraph,Disposition,MaterialKind as K,InputRole as R,StepKind as S,MolecularRelation as M
from chemistry_validation import Status,check_graph


def run(check):
    check.section('explicit assembly geometry and negative controls')
    core=(Segment('core','GCTA'),)
    fragment=CohesiveFragment('fragment','AC',core,'TG','5-prime','5-prime')
    retention=StrandRetention('parent','ready','ACGCTATG',revcomp('ACGCTATG'),(0,6),(0,6),fragment,'fixture intervals')
    check.raises('wrong retained strand boundary is rejected',
                 lambda:replace(retention,top_span=(1,6)))
    check.raises('opposite end polarity cannot reuse the same retained intervals',
                 lambda:replace(retention,fragment=replace(fragment,left_polarity='3-prime')))
    check('removed prefix and suffix remain separate reference intervals',
          replace(retention,top='GGACGCTATGCC',bottom='GGCATAGCGTCC',top_span=(2,8),bottom_span=(2,8)).removed(),
          ('GG','TGCC','GG','GTCC'))
    g=ProtocolGraph('retention fixture',required_output=Disposition.CONSTRUCT)
    before=Scene.duplex((Segment('parent','ACGCTATG'),))
    a=g.input('parent','parent',K.MOLECULAR,primary=True,state=before)
    enzyme=g.input('enzyme','caller enzyme',K.CATALYST)
    b=g.product('ready','retained',K.MOLECULAR,state=fragment,chemical_state=fragment)
    g.react('cut','explicit intervals',S.DIGESTION,[(a,R.TEMPLATE),(enzyme,R.CATALYST)],[b],geometry=retention,relations=(M.STRAND_RETENTION,))
    check('independent interval check accepts the fixture',any(c.status is Status.ESTABLISHED for c in check_graph(g)),True)
    g.materials['ready']=replace(b,chemical_state=replace(fragment,core=(Segment('mutated','GCTT'),)))
    check('independent interval check rejects a substituted product',any(c.status is Status.CONTRADICTORY for c in check_graph(g)),True)
    a=Construct((Segment('a','ACGTTG'),));b=Construct((Segment('b','TACGGA'),))
    geometry=BridgedAssembly((('a',a),('b',b)),(('ab',revcomp(a.top()[-3:]+b.top()[:3])),('ba',revcomp(b.top()[-3:]+a.top()[:3]))),3,True,True,'explicit fixture end chemistry')
    labelled=replace(geometry,fragments=(('a',Construct((Segment('a',a.top(),feature=feature('fixture.index','sample_index','fixed')),))),('b',b)))
    check('joining fragments preserves structured identifier meaning',
          tuple(labelled.product('labelled'))[0].feature.id,'fixture.index')
    check.raises('incompatible bridge is rejected',
                 lambda:replace(geometry,bridges=(('ab','AAAAAA'),geometry.bridges[1])))
    check.raises('missing phosphate is not repaired by a geometry default',
                 lambda:replace(geometry,five_prime_phosphorylated=False))
    check.raises('blocked hydroxyl prevents a ligatable declaration',
                 lambda:replace(geometry,three_prime_hydroxyl=False))
    g=ProtocolGraph('bridge fixture',required_output=Disposition.CONSTRUCT)
    inputs=[]
    for identifier,con in geometry.fragments:inputs.append(g.input(identifier,identifier,K.MOLECULAR,primary=True,state=con))
    bridges=[]
    for identifier,sequence in geometry.bridges:bridges.append(g.input(identifier,identifier,K.OLIGO,state=Construct((Segment(identifier,sequence),))))
    ligase=g.input('ligase','fixture ligase',K.CATALYST)
    out=g.product('out','joined',K.MOLECULAR,state=geometry.product('out'))
    g.react('ligate','bridge-directed ligation',S.LIGATION,[(inputs[0],R.TEMPLATE),(inputs[1],R.STRUCTURAL)]+[(m,R.STRUCTURAL) for m in bridges]+[(ligase,R.CATALYST)],[out],geometry=geometry,relations=(M.BRIDGED_ASSEMBLY,))
    check('independent bridge relation accepts the fixture',any(c.status is Status.ESTABLISHED for c in check_graph(g)),True)
    g.materials['out']=replace(out,state=Construct((Segment('incorporated bridge',a.top()+geometry.bridges[0][1]+b.top()),)))
    check('bridge incorporation is independently rejected',any(c.status is Status.CONTRADICTORY for c in check_graph(g)),True)
    check.raises('illustrative graph requires an invented-component description',
                 lambda:ProtocolGraph('unqualified example',representation='illustrative_example'))
    check('legacy graph provenance stays unspecified',ProtocolGraph('legacy').representation,'unspecified')
    h=ProtocolGraph('recombination shape')
    parent=h.input('parent','template',K.MOLECULAR)
    cat=h.input('cat','catalyst',K.CATALYST)
    product=h.product('product','product',K.MOLECULAR)
    check.raises('site-specific recombination cannot omit its second template',
                 lambda:h.react('exchange','exchange',S.SITE_SPECIFIC_RECOMBINATION,[(parent,R.TEMPLATE),(cat,R.CATALYST)],[product]))

    check.section('reference PCR copying and overlap extension negative controls')
    left=Construct((Segment('left','ACGTCAGTACGATCGTGACTACGA'),))
    right=Construct((Segment('right','GATCTACGAGTCTCGAGTACCTGACATGTCAGATCA'),))
    extended_left=Construct((Segment('left','ACGTCAGTACGATCGTGACTACGAGATCTACGAGTC'),))
    overlap=OverlapExtension('left','right',extended_left,right,12,'explicit fixture overlap')
    check.raises('unequal terminal overlap is rejected',lambda:replace(overlap,overlap_nt=11))
    g=ProtocolGraph('overlap fixture',required_output=Disposition.CONSTRUCT)
    a=g.input('left','left fragment',K.MOLECULAR,primary=True,state=extended_left)
    b=g.input('right','right fragment',K.MOLECULAR,primary=True,state=right)
    polymerase=g.input('polymerase','polymerase',K.CATALYST)
    product=g.product('product','product',K.MOLECULAR,state=overlap.product('product'))
    g.react('extension','extension',S.EXTENSION,[(a,R.TEMPLATE),(b,R.TEMPLATE),(polymerase,R.CATALYST)],[product],geometry=overlap,relations=(M.OVERLAP_EXTENSION,))
    check('overlap relation accepts one retained copy',any(c.status is Status.ESTABLISHED for c in check_graph(g)),True)
    g.materials['product']=replace(product,state=Construct((Segment('duplicated',extended_left.top()+right.top()),)))
    check('independent overlap relation rejects duplicated junction DNA',any(c.status is Status.CONTRADICTORY for c in check_graph(g)),True)
    f=Construct((Segment('forward','GGAT'+left.top()[:12],inferred=True),))
    r=Construct((Segment('reverse','TTCA'+revcomp(left.top()[-12:]),inferred=True),))
    geometry=PCRGeometry.exact(left,PCRPrimer('f',f.top(),12,inferred=True),PCRPrimer('r',r.top(),12,inferred=True),statement='explicit fixture primer tails')
    copied=pcr_product(geometry,f,r,name='reference')
    check('reference amplicon retains independently supplied primer tails',copied.top(),'GGAT'+left.top()+'TGAA')
    check.raises('amplicon derivation rejects a substituted primer',lambda:pcr_product(geometry,Construct((Segment('wrong','AAAA'),)),r,name='bad'))
    g=ProtocolGraph('PCR fixture',required_output=Disposition.CONSTRUCT)
    template=g.input('template','template',K.MOLECULAR,primary=True,state=left)
    fm=g.input('f','forward',K.OLIGO,state=f);rm=g.input('r','reverse',K.OLIGO,state=r)
    polymerase=g.input('polymerase','polymerase',K.CATALYST)
    product=g.product('product','product',K.MOLECULAR,state=copied)
    carry=g.product('carry','unchanged parent class',K.MOLECULAR,state=left)
    geometry=replace(geometry,amplicon_id='product')
    g.react('pcr','PCR',S.PCR,[(template,R.TEMPLATE),(fm,R.STRUCTURAL),(rm,R.STRUCTURAL),(polymerase,R.CATALYST)],[product,carry],geometry=geometry,relations=(M.PCR_PRODUCT,))
    check('independent product relation accepts the tailed PCR reference',not any(c.status is Status.CONTRADICTORY for c in check_graph(g)),True)
    g.materials['product']=replace(product,state=left)
    check('independent product relation rejects lost primer tails',any(c.status is Status.CONTRADICTORY for c in check_graph(g)),True)
    g.steps['pcr']=replace(g.steps['pcr'],geometry=replace(geometry,amplicon_id='not_an_output'))
    check('an amplicon target must belong to its actual PCR step',any(c.status is Status.CONTRADICTORY for c in check_graph(g)),True)
    core=(Segment('core','CTAG'),)
    fragment=CohesiveFragment('boundary','AGTA',core,'ACTA','5-prime','5-prime')
    boundary=NucleotideBoundary('p','q','AGTACTAGACTA',revcomp('AGTACTAGACTA'),(0,8),(0,8),fragment,'fixture 3-prime boundary','G')
    check.raises('a different boundary nucleotide cannot silently reuse the declared endpoints',lambda:replace(boundary,boundary_base='C'))
    check.raises('a forbidden base in the removed suffix is rejected',lambda:replace(boundary,top='AGTACTAGAGTA'))

    check.section('activated end-transfer negative controls')
    vf=CohesiveFragment('recipient','A',(Segment('v','TCAGGATC'),),'T','3-prime','3-prime')
    inf=CohesiveFragment('insert','T',(Segment('i','CGATGACC'),),'A','3-prime','3-prime')
    vector=DuplexEndState(vf,'unassigned','phosphotyrosyl','unassigned','phosphotyrosyl')
    insert=DuplexEndState(inf,'OH','OH','OH','OH')
    transfer=ActivatedEndTransfer('vector','insert',vector,insert,'fixture activated ends')
    check.raises('unactivated donors cannot use the transfer relation',lambda:replace(transfer,vector=replace(vector,top3='OH')))
    check.raises('phosphorylated insert is not a hydroxyl acceptor',lambda:replace(transfer,insert=replace(insert,top5='phosphate')))
    g=ProtocolGraph('transfer fixture',required_output=Disposition.CONSTRUCT)
    vm=g.input('vector','vector',K.MOLECULAR,primary=True,state=vector)
    im=g.input('insert','insert',K.MOLECULAR,primary=True,state=insert)
    out=g.product('out','reference',K.MOLECULAR,state=transfer.reference('out'))
    g.react('transfer','transfer',S.ADAPTER_ATTACHMENT,[(vm,R.TEMPLATE),(im,R.STRUCTURAL)],[out],geometry=transfer,relations=(M.ACTIVATED_END_TRANSFER,))
    check('independent activated-end relation accepts the fixture',any(c.status is Status.ESTABLISHED for c in check_graph(g)),True)
    g.materials['insert']=replace(im,state=replace(insert,bottom5='phosphate'))
    check('changed actual acceptor end chemistry is independently rejected',any(c.status is Status.CONTRADICTORY for c in check_graph(g)),True)

    check.section('linear linker retention and conditional host reference controls')
    pf=CohesiveFragment('prefix','AGTC',(Segment('prefix core','ACGTA'),),'TGAC','5-prime','5-prime')
    body=CohesiveFragment('body','TGAC',(Segment('body core','CTAGA'),),'CAGT','5-prime','5-prime')
    sf=CohesiveFragment('suffix','CAGT',(Segment('suffix core','TCGAT'),),'GACT','5-prime','5-prime')
    geometry=LinearFivePrimeAssembly((('prefix',pf),('body',body),('suffix',sf)),True,True,'fixture supplied end chemistry')
    check.raises('linear linker ligation rejects incompatible internal ports',lambda:replace(geometry,fragments=(('prefix',pf),('body',replace(body,left='AAAA')),('suffix',sf))))
    g=ProtocolGraph('linear fixture',required_output=Disposition.CONSTRUCT)
    ms=[g.input(i,i,K.MOLECULAR,primary=True,state=f) for i,f in geometry.fragments]
    enzyme=g.input('ligase','ligase',K.CATALYST)
    product=g.product('product','product',K.MOLECULAR,state=geometry.product('product'))
    g.react('join','linear joining',S.LIGATION,[(ms[0],R.TEMPLATE)]+[(m,R.STRUCTURAL) for m in ms[1:]]+[(enzyme,R.CATALYST)],[product],geometry=geometry,relations=(M.LINEAR_COHESIVE_ASSEMBLY,))
    check('linear product relation accepts retained external ends',not any(c.status is Status.CONTRADICTORY for c in check_graph(g)),True)
    g.materials['product']=replace(product,state=replace(product.state,left='AAAA'))
    check('independent linear relation rejects a changed external tail',any(c.status is Status.CONTRADICTORY for c in check_graph(g)),True)
    a=Construct((Segment('a','AGTCGATTCAGT'),));b=Construct((Segment('b','CAGTACGAAGTC'),))
    ref=CircularOverlapReference((('a',a,4,4),('b',b,4,4)),'fixture intended reference only')
    check('conditional host references default to inferred product segments',all(s.inferred for s in ref.reference('reference')),True)
    check.raises('circular intended-reference constructor rejects unequal overlaps',lambda:replace(ref,parts=(('a',a,4,4),('b',Construct((Segment('bad','AAAATCGTAGTC'),)),4,4))))
    g=ProtocolGraph('host reference fixture',required_output=Disposition.CONSTRUCT)
    ma=g.input('a','a',K.MOLECULAR,primary=True,state=a);mb=g.input('b','b',K.MOLECULAR,primary=True,state=b)
    out=g.product('out','intended reference',K.MOLECULAR,state=ref.reference('out'))
    g.react('host','unmodelled host process',S.MAGIC,[(ma,R.TEMPLATE),(mb,R.TEMPLATE)],[out],geometry=ref,reason='fixture unresolved biology')
    claims=check_graph(g)
    check('matching reference never upgrades the unmodelled host transformation',any(c.status is Status.UNRESOLVED and c.target=='host' for c in claims),True)
    check('reference arithmetic remains explicitly conditional',any(c.status is Status.CONDITIONAL and c.question=='conditional circular reference' for c in claims),True)
    g.materials['out']=replace(out,state=Construct((Segment('duplicated',a.top()+b.top()),)))
    check('independent intended-reference check rejects duplicated overlaps',any(c.status is Status.CONTRADICTORY for c in check_graph(g)),True)
    g.materials['b']=replace(mb,state=None)
    check('unknown contributor remains unresolved rather than confirmed contradictory',any(c.status is Status.UNRESOLVED and c.question=='conditional circular reference' for c in check_graph(g)),True)

    check.section('reference partition accounting independent of nuclease chemistry')
    parent=Construct((Segment('parent','ACGTCAGTGCTA'),))
    geom=DuplexPartition('parent',parent,(('left',0,3),('target',3,8),('right',8,12)),'chosen reference boundaries only')
    check.raises('partition rejects an omitted parent interval',lambda:replace(geom,intervals=(('left',0,3),('right',4,12))))
    check.raises('partition rejects duplicated parent interval',lambda:replace(geom,intervals=(('left',0,8),('right',7,12))))
    check.raises('partition rejects lost suffix',lambda:replace(geom,intervals=(('left',0,3),('right',3,11))))
    g=ProtocolGraph('reference partition fixture',required_output=Disposition.CONSTRUCT)
    source=g.input('parent','parent',K.MOLECULAR,primary=True,state=Scene.duplex(tuple(parent)),chemical_state=parent)
    enzyme=g.input('enzyme','unassigned nuclease',K.CATALYST)
    products=[g.product(i,i,K.MOLECULAR,state=Scene.duplex(tuple(con)),chemical_state=con) for i,con in geom.products()]
    g.react('partition','chosen reference partition',S.DIGESTION,[(source,R.TEMPLATE),(enzyme,R.CATALYST)],products,geometry=geom,relations=(M.DUPLEX_PARTITION,))
    check('reference partition remains conditional even on explicit bases',any(c.status is Status.CONDITIONAL and c.question==M.DUPLEX_PARTITION.value for c in check_graph(g)),True)
    target=g.materials['target']
    g.materials['target']=replace(target,chemical_state=Construct((Segment('changed','CAGTGA'),)))
    check('independent partition relation rejects changed output',any(c.status is Status.CONTRADICTORY for c in check_graph(g)),True)
    g.materials['target']=target
    g.materials['parent']=replace(source,chemical_state=Construct((Segment('changed','TCGTCAGTGCTA'),)))
    check('independent partition relation rejects changed parent',any(c.status is Status.CONTRADICTORY for c in check_graph(g)),True)

    check.section('explicit duplex reference views across recursive assembly')
    af=CohesiveFragment('a','AGTC',(Segment('a body','ACGTA'),),'TGAC','5-prime','5-prime')
    bf=CohesiveFragment('b','TGAC',(Segment('b body','TCGAT'),),'AGTC','5-prime','5-prime')
    geom=CircularCohesiveAssembly((('a',af),('b',bf)),'fixture circular order')
    g=ProtocolGraph('recursive view fixture',required_output=Disposition.CONSTRUCT)
    a=g.input('a','a',K.MOLECULAR,primary=True,state=af);b=g.input('b','b',K.MOLECULAR,primary=True,state=bf)
    ligase=g.input('ligase','ligase',K.CATALYST)
    scene=Scene.duplex(tuple(geom.product('reference')))
    product=g.product('product','product',K.MOLECULAR,state=scene)
    g.react('seal','seal',S.LIGATION,[(a,R.TEMPLATE),(b,R.STRUCTURAL),(ligase,R.CATALYST)],[product],geometry=geom,relations=(M.COHESIVE_ASSEMBLY,))
    check('cohesive relation reads an explicitly supplied full duplex reference',any(c.status is Status.ESTABLISHED and c.question==M.COHESIVE_ASSEMBLY.value for c in check_graph(g)),True)
    del scene.strands['bottom']
    check('one surviving strand cannot silently become a full duplex reference',any(c.status is Status.UNRESOLVED and c.question==M.COHESIVE_ASSEMBLY.value for c in check_graph(g)),True)


if __name__=='__main__':
    check=Check();run(check);check.report()
