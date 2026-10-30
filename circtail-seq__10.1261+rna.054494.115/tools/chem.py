"""Protocol-owned source-qualified partial chemistry. No public schematic endpoint."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'lib'))
from chemdraw import Scene, Segment
from protocol_graph import ProtocolGraph, MaterialKind as K, InputRole as R, StepKind as S, SourceItem, SequencingPlan, Terminal, Disposition
from pcr import PCRGeometry
from circular import CircularizationGeometry
from adapter_ligation import AdapterLigationGeometry

def scene(name, text):
    sc = Scene(); sc.strand(name, [Segment(name, text, placeholder=True)]); return sc

def op(g, key, x, label, kind=S.MAGIC, *, reagent='Source mixture; internal chemistry unassembled', product=None, structural=(), catalyst_stock=None):
    catalyst = kind in {S.DIGESTION, S.END_REPAIR, S.EXTENSION, S.TAILING, S.REVERSE_TRANSCRIPTION, S.PCR, S.LIGATION, S.DEPHOSPHORYLATION, S.PHOSPHORYLATION, S.CIRCULARIZATION}
    stock = catalyst_stock or g.input(key + '_reagent', reagent, K.CATALYST if catalyst else K.CONDITION)
    inputs = [(x, R.TEMPLATE), (stock, R.CATALYST if catalyst else R.CONDITION)]
    inputs.extend((a, R.STRUCTURAL) for a in structural)
    kwargs = {}
    if kind == S.MAGIC:
        kwargs['reason'] = 'Source operation established; exact molecular transition/kit internals unassembled. No accepted absence claim.'
    if kind == S.PCR:
        kwargs['geometry'] = PCRGeometry.unavailable(Scene(), reason='Exact template/primer endpoint geometry unassembled; selected oligos do not establish whole amplicon.')
        if len(structural) < 2:
            for arm in ('f', 'r'):
                inputs.append((g.input(key + '_' + arm, 'Source PCR primer; exact identity/landing unassembled', K.OLIGO, state=scene(key + '_' + arm, '[source primer]')), R.STRUCTURAL))
    if kind == S.REVERSE_TRANSCRIPTION and not structural:
        inputs.append((g.input(key + '_primer', 'Source RT primer; geometry unassembled', K.OLIGO, state=scene(key + '_primer', '[source RT primer]')), R.STRUCTURAL))
    if kind == S.CIRCULARIZATION:
        kwargs['geometry'] = CircularizationGeometry.unavailable('Source RNA substrates and closure end geometry unassembled; no source absence claim.')
    if kind == S.LIGATION:
        kwargs['geometry'] = AdapterLigationGeometry('unavailable', 'unavailable', (scene(key + '_ends', '[source ligation; complete substrate/end geometry unassembled]'),), 'This incomplete model does not claim publication absence.')
        if not structural:
            inputs.append((g.input(key + '_adapter', 'Source ligating adapter; full geometry unassembled', K.ADAPTER, state=scene(key + '_adapter', '[source adapter]')), R.STRUCTURAL))
    y = g.product(key + '_product', product or label + '; physical companions retained until separation', K.MOLECULAR)
    g.react(key, label, kind, inputs, (y,), **kwargs); return y

def clean(g, key, x, label):
    condition = g.input(key + '_condition', label, K.CONDITION)
    kept = g.product(key + '_kept', 'Source physically recovered fraction; carryover/yield unquantified', K.MOLECULAR)
    waste = g.product(key + '_waste', 'Actual separated cleanup fractions; no perfect classification assumed', K.WASTE)
    g.react(key, label, S.CLEANUP, ((x, R.TEMPLATE), (condition, R.CONDITION)), (kept, waste))
    g.discard(waste, 'Actual physical separation rather than computational removal.'); return kept

def qc(g, key, x, label):
    kept = g.product(key + '_kept', 'Source retained preparation', K.MOLECULAR)
    assay = g.product(key + '_aliquot', 'Separate analytical aliquot; no return to library assumed', K.QC)
    g.react(key, label, S.SPLIT, ((x, R.TEMPLATE),), (kept, assay))
    g.mark_terminal(Terminal(assay.id, Disposition.QC, 'Source analytical aliquot.')); return kept

from illumina import P5,P7,TRUSEQ_READ1,TRUSEQ_READ2
from chemdraw import feature
SOURCE='PublishedPMC4748824 Results60–78/Methods150–160; publishedSupplementaryMethodsteps1–12/TableS1p2rendered. LinkedRiley2014recoveryoriginalunreviewed.'
RT={'co1':'GTAATGAGTACGTTGTAAAACTG','co3p':'AAACCCTCCTCCTCTGGCGGT'}
REV={'co1':'ACTTAAGCACACAAGACATAG','co3p':'GCCTTTTTATTCAATTTTAAACAATCC'}
FWD={'co1':'TTTTATATTTTTTTTTGACAAG','co3p':'AGGGTTTAGTTGGAATGAA'}
INDEX={('co1','a'):'CGATGT',('co1','b'):'TGACCA',('co3p','a'):'ACAGTG',('co3p','b'):'CAGATC'}

def primer_scenes(gene,replicate):
    if gene not in RT or replicate not in {'a','b'}:raise ValueError('Source-definedgene/primerpair required')
    rev=Scene();rev.strand('rev',[Segment('P5',P5),Segment('SourceRead1arm',TRUSEQ_READ1[4:]),Segment('Source5Ndiversityoffset','N'*5,placeholder=True),Segment('Sourcegene-specificreverse',REV[gene])])
    fwd=Scene();fwd.strand('fwd',[Segment('P7',P7),Segment('Source6ntsampleindex',INDEX[gene,replicate],feature=feature('circtail.sample','sample_index','fixed',note='TableS1boldsamplebarcode;notUMI')),Segment('SourceRead2arm',TRUSEQ_READ2),Segment('SourceNNdiversityoffset','N'*(4 if replicate=='a' else 6),placeholder=True),Segment('Sourcegene-specificforward',FWD[gene])])
    return rev,fwd

def aliquot(g,key,x,used,rest):
    a=g.product(key+'_used',used,K.MOLECULAR);b=g.product(key+'_rest',rest,K.MOLECULAR)
    g.react(key,'Physicalsourcealiquot;remainderallocationunassigned',S.SPLIT,((x,R.TEMPLATE),),(a,b));g.mark_terminal(Terminal(b.id,Disposition.UNRESOLVED,'Unusedsourcefraction,notchemicalwaste'));return a

def protocol_graph_for(gene='co1',replicate='a'):
    rev,fwd=primer_scenes(gene,replicate)
    g=ProtocolGraph('circTAIL-seq '+gene+' primer'+replicate+' — source-selected partial path',authority=SOURCE,representation='published_instance')
    x=g.input('rna','Source5ugDNase-treatedtotalTbruceiRNA;oneindependentbiologicalsample',K.MOLECULAR,primary=True)
    x=op(g,'circle',x,'T4RNAligase200ulovernightRT;25uMATP',S.CIRCULARIZATION,reagent='T4RNAligase/RNaseinhibitor/Tris-acetate/ammonium-acetate/Mg/DTT/ATP',product='RNAcirclespluslinearRNA/intermolecularchimeras/reagents; no perfectintramolecularpuritypredicted')
    x=clean(g,'rna_clean',x,'pH5.2phenol-chloroform/overnightprecipitation/wash/20ulwater')
    x=aliquot(g,'rt_input',x,'Source1ugrecoveredRNAinto20ulRT','RemainingrecoveredRNA;mass/yieldnotinferredfrom5uginput')
    primers=[]
    for target,bases in RT.items():
        sc=Scene();sc.strand(target,[Segment('SourceRTprimer',bases)]);primers.append(g.input('rt_'+target,'Published'+target+'cDNAprimer;allgene-specificprimerssourcepoolmembershipbeyondTableS1unassigned',K.OLIGO,state=sc))
    x=op(g,'rt',x,'EpiscriptRTwithgene-specificprimerpool;sourceconditionsbeyondnamedRTunreviewed',S.REVERSE_TRANSCRIPTION,reagent='Episcript/sourceRTbuffer;2pmoleachprimer',structural=tuple(primers),product='cDNAcopiesfromRNAcirclejunctionsandotherphysicaltemplatesplusRNA/primers/reagents')
    x=aliquot(g,'production_input',x,'SourcecDNAproductioninput;0.25ulperstatedreaction/scalingto100ulunassigned','Remaining20ulRTpreparationafteroptimization/productionallocations;exactvolumenotassigned')
    fp=g.input('fwd','Published'+gene+'3primefwd'+replicate+'indexedprimer',K.OLIGO,state=fwd);rp=g.input('rev','Published'+gene+'5primereverseprimer',K.OLIGO,state=rev)
    aliquots=tuple(g.product('reaction_'+str(i),'One20ulaliquotofsource100ulproductionPCR',K.MOLECULAR) for i in range(5))
    c=g.input('pcr_mix','SourceKAPA2GRobust/buffer/weighted8mMdGTP,dCTP/12mMdTTP,dATP;whetherprintedstockorfinalunassigned',K.CATALYST)
    mix=g.product('production_mix','Source100ulPCRmastermixwithtemplateandpairedprimers',K.MOLECULAR)
    g.react('assemble_pcr','SourceproductionPCRassembly',S.GENERIC,((x,R.TEMPLATE),(fp,R.STRUCTURAL),(rp,R.STRUCTURAL),(c,R.CATALYST)),(mix,))
    g.react('split_reactions','SupplementMethod12divides100ulinto20ulaliquots',S.SPLIT,((mix,R.TEMPLATE),),aliquots)
    products=[]
    for i,a in enumerate(aliquots):
        products.append(op(g,'pcr_'+str(i),a,'Source'+('33'if gene=='co1'else'32')+'cycles95/62/72C15seach;enzymeandprimersalreadycarriedinthisaliquot',S.PCR,reagent='KAPA2GRobust/knownsourceprimersalreadyinmastermix',structural=(fp,rp),catalyst_stock=c))
    x=g.product('recombined','Recombinedfive20ulPCRaliquots',K.MOLECULAR);g.react('recombine','Sourcepost-PCRrecombination',S.POOL,tuple((a,R.TEMPLATE)for a in products),(x,))
    x=clean(g,'page',x,'6percentnativePAGE/selectednominal0–150nt-tailrange;Riley2014recoverywith400ulelution/twocolumns/ethanol/15ulwater')
    x=qc(g,'bioanalyzer',x,'SeparateBioanalyzerquality/quantityassay')
    x=qc(g,'kapaqc',x,'SeparateKapaQCtestofamplifiability')
    other=g.input('other_libraries','Otherfivelibrariesfortheactualsix-libraryrun;identities/amountsunassigned',K.LIBRARY,primary=True)
    pool=g.product('pooled','Six-libraryapproximatelyequalBioanalyzerquantitypool',K.MOLECULAR);g.react('pool','Physicalsix-librarypool',S.POOL,((x,R.TEMPLATE),(other,R.TEMPLATE)),(pool,))
    phix=g.input('phix','IndependentdiversePhiXlibrary;sourcefractionunassigned',K.LIBRARY,primary=True);mixed=g.product('with_phix','PhysicalcircTAIL/PhiXloadingmixture',K.LIBRARY)
    g.react('spike','SourcePhiXadditiontounderclusteredrun',S.POOL,((pool,R.TEMPLATE),(phix,R.TEMPLATE)),(mixed,))
    g.sequence(mixed,SequencingPlan(limitation='PublishedMiSeqV2paired-end150bp;actualindexed/offsetprimercoordinatesnotassembled. SourceNNNN/NNNNNNdiversityoffsetsaretrimmedcomputationally,notmolecularUMIs.'))
    g.add_source_item(SourceItem('pending','FullgeneRNA/circle/PCRjunctiongeometry,sourceoptimizationreactionallocations,remainingrunmembers andRiley2014original','Privatepartialgraph;source-tableoligosonly,nofullsource/absenceaudit.'));return g

def protocol_graph():return protocol_graph_for()
def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph();primer_scenes('co3p','b')
    c('ProductionaliquotsrecombinedonlyafterPCR',len(g.steps['recombine'].inputs),5)
    c('PCRtemplateisonephysicalreaction',g.steps['pcr_0'].inputs[0].material,'reaction_0')
    c('ProductionpoolfollowsindependentQCfractions',g.steps['pool'].inputs[0].material,'kapaqc_kept')
    c('PhiXisindependentphysicalspike',g.steps['spike'].inputs[1].material,'phix')
    c('RNAremainderretained',g.terminals['rt_input_rest'].disposition,Disposition.UNRESOLVED)
    c('Noexactampliconclaim',bool(g.publication_problems()),True)
    src=Source(sources_root()/'circtail-seq__10.1261+rna.054494.115'/'rna.054494.115_supp_supp_054494.115_SuppMaterial.pdf.txt','get_sources --slugs circtail-seq__10.1261+rna.054494.115')
    if have(c,src,label='publishedTableS1'):
        t=''.join(src.path.read_text().split()).upper()
        for gene in RT:c('WholeRTprimer'+gene,RT[gene]in t,True)
        for gene in RT:
            for rep in ('a','b'):
                whole=P7+INDEX[gene,rep]+TRUSEQ_READ2+'N'*(4 if rep=='a' else 6)+FWD[gene]
                c('Wholeindexedprimer'+gene+rep,whole in t,True)
    c.report()
if __name__=='__main__':selftest()
