"""Protocol-owned source-qualified partial chemistry. No public schematic endpoint."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'lib'))
from chemdraw import Scene, Segment, revcomp
from protocol_graph import ProtocolGraph, MaterialKind as K, InputRole as R, StepKind as S, SourceItem, SequencingPlan, Terminal, Disposition
from pcr import PCRGeometry
from adapter_ligation import AdapterLigationGeometry

def scene(name, text):
    sc = Scene(); sc.strand(name, [Segment(name, text, placeholder=True)]); return sc

def oligo(g, key, name, bases, **mods):
    sc = Scene()
    sc.strand(key, [Segment(key, bases, placeholder=bool(set(bases)-set('ACGTU')))], **mods)
    return g.input(key, name, K.OLIGO, state=sc)

def op(g, key, x, label, kind=S.MAGIC, *, reagent='Source mixture; internal chemistry unassembled', product=None, structural=()):
    catalyst = kind in {S.DIGESTION, S.END_REPAIR, S.EXTENSION, S.TAILING, S.REVERSE_TRANSCRIPTION, S.PCR, S.LIGATION, S.DEPHOSPHORYLATION, S.PHOSPHORYLATION}
    stock = g.input(key + '_reagent', reagent, K.CATALYST if catalyst else K.CONDITION)
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

SOURCE = 'Published PMC4509939 Results70/Methods130–134; published SuppMethods steps1–7/reagent list. Selected detailed SI route; main differences explicit in note.'
END_EXTEND = 'CAGACGTGTGCTCTTCCGATCTTTTTTTTTTTTTTTTTT'
SPLINT_A_DNA = 'CCCTACACGACGCTCTTCCG'
SPLINT_A_RNA_AS_WRITTEN = 'ATCT'
SPLINT_B = 'NNNNAGATCGGAAGAGCGTCGTGTAGGG'
UNIVERSAL = 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCG'

def splint_a():
    sc=Scene();sc.strand('splint_a',[Segment('DNA arm',SPLINT_A_DNA),Segment('Source ribonucleotide block',SPLINT_A_RNA_AS_WRITTEN,note='Source (rA)(rT)(rC)(rT) retained literally; rT not silently changed to rU')]);return sc

def splint_scene():
    # This establishes only complementary base letters, not the uncertain rT identity.
    if SPLINT_B != 'NNNN'+revcomp(SPLINT_A_DNA+SPLINT_A_RNA_AS_WRITTEN):
        raise ValueError('Published full splint B must complement the full declared A arm')
    sc=splint_a()
    sc.anneal('splint_b',[Segment('Source degenerate overhang','NNNN',placeholder=True),Segment('Pairs to source RNA block',revcomp(SPLINT_A_RNA_AS_WRITTEN)),Segment('Pairs to DNA arm',revcomp(SPLINT_A_DNA))],to='splint_a',pair=('Pairs to DNA arm','DNA arm'),unpaired=('Source degenerate overhang',))
    return sc

def protocol_graph():
    splint_scene()  # Construction checks the source-known duplex before exposing its oligos.
    g=ProtocolGraph('PAT-seq — selected published supplementary workflow',authority=SOURCE,representation='published_instance')
    x=g.input('total_rna','Source1ug total yeastRNA; target and other RNAs retained',K.MOLECULAR,primary=True)
    ee=g.input('end_extend','Published5primeBio end-extend DNA; template then bead-bound RT primer, not a covalent RNA-biotin tag',K.OLIGO,state=(lambda sc:(sc.strand('ee',[Segment('Published Bio-EE DNA',END_EXTEND)],mod5="biotin"),sc)[1])(Scene()))
    x=op(g,'extend',x,'Klenow exo-minus3prime RNA extension on annealed DNA template37C30min;80C5min',S.EXTENSION,reagent='Klenow exo-minus/dNTP/SSIIIbuffer/DTT/RNaseOUT',structural=(ee,),product='RNA-DNA chimeric3prime-extended classes annealed to separate biotinylated DNA templates; unextendedRNA/freeoligos retained')
    x=op(g,'digest',x,'SI RNaseT1 dilution1:1000 digestion on ice2min; main says1min',S.DIGESTION,reagent='Source RNaseT1/digestbuffer',product='PhysicalRNAfragments plus retainedDNA/polyAtails; no guaranteed perfectpolyAselection')
    x=clean(g,'phenol',x,'Immediateacidphenol/chloroform quench and chloroform extraction')
    beads=g.input('beads','Streptavidinmagneticbeads',K.SELECTION)
    bound=g.product('bound','Biotin-DNA associatedRNA/otherboundmaterial',K.MOLECULAR);w=g.product('capture_unbound','Actualunboundphysicalfraction',K.WASTE)
    g.react('capture','DNA:RNAduplexcapture15min37C then magneticseparation',S.AFFINITY_CAPTURE,((x,R.TEMPLATE),(beads,R.SELECTION)),(bound,w));g.discard(w,'Physicalunboundfraction;not computationalfilter')
    x=clean(g,'capture_wash',bound,'Source55C dilutedSSIIIbuffer wash and37C ligasebufferwash')
    x=op(g,'pnk',x,'T4PNK5primephosphorylation37C30min;70C10min',S.PHOSPHORYLATION,reagent='T4PNK/ATP-containingT4DNAligasebuffer')
    x=clean(g,'pnk_wash',x,'Sourcebeadwash afterPNK')
    a=g.input('splint_a','PublishedmixedDNA/RNA splintA; source rT retained',K.ADAPTER,state=splint_a());b=oligo(g,'splint_b','PublishedsplintB in5prime-to3prime orientation;NNNNannealingwindow meaningunassigned',SPLINT_B)
    x=op(g,'ligate',x,'PreannealedsplintA/B;65C10min→16C5min thenT4RNA ligase2 16C16h',S.LIGATION,reagent='T4RNAligase2/sourcePEGbuffers',structural=(a,b),product='SplintA-adaptedRNAclasses onbiotin-DNAmatrix plusunligatedRNA/adapters/splintB/reagents')
    x=clean(g,'ligation_wash',x,'MagneticwashremovesexcesssplintbeforeRT')
    x=op(g,'rt',x,'Bead-boundBio-EE primesSSIII55C60min;70C5min',S.REVERSE_TRANSCRIPTION,reagent='SSIII/dNTP/DTT/RNaseOUT',structural=(ee,),product='BiotinylatedfirststrandDNAcopies plusRNA/matrix/primer/reagentremnants')
    x=clean(g,'rt_wash',x,'MagneticwashafterRT;optionalstorageunmodelled')
    released=g.product('released','Formamide80C5minreleasedDNA-containingliquid',K.MOLECULAR);matrix=g.product('matrix','Sedimentedmagneticmatrixandassociatedmaterial',K.WASTE)
    c=g.input('formamide','2xformamideloadingbuffer80C5min/ice/pulsespin',K.CONDITION)
    g.react('release','Denatureandpulsespinbeforeloadinggel',S.STRAND_SEPARATION,((x,R.TEMPLATE),(c,R.CONDITION)),(released,matrix));g.discard(matrix,'Actualsedimentedmatrix;biotin-bearingDNAcanreleasefromstreptavidinunderharshconditions,no stoichiometricprediction')
    x=clean(g,'page',released,'6percentureaPAGE120–300baseyeastwindow/crushsoak/SpinX;sourceSIisopropanolprecipitation')
    production=g.product('pcr_input','10ulof30ulselectedcDNA',K.MOLECULAR);remaining=g.product('remaining_cdna','Remaining20ulcDNA;fateunassigned',K.MOLECULAR)
    g.react('aliquot','Physicalone-thirdPCRinput',S.SPLIT,((x,R.TEMPLATE),),(production,remaining));g.mark_terminal(Terminal(remaining.id,Disposition.UNRESOLVED,'SourcePCRuses10of30ul;remainingfatenotassigned'))
    f=oligo(g,'universal','Publisheduniversalprimer',UNIVERSAL)
    x=op(g,'pcr',production,'AmpliTaqGold36016cyclePCR/indexintroduction',S.PCR,reagent='AmpliTaqGold360/sourceScriptSeqindexreverseprimer',structural=(f,g.input('index_primer','ScriptSeqIndexPCRreverseprimer;exactchosenindex/landingunassembled',K.OLIGO,state=scene('idx','[source indexed primer]'))))
    x=qc(g,'agarose',x,'Separate10percentPCRproducthighresolution2percentagaroseassay;intensityinforms pool')
    other=g.input('other_libraries','OtherindependentlypreparedPATlibraries;sourcepoolmembersunassigned',K.LIBRARY,primary=True)
    pool=g.product('pooled','PhysicalapproximatelyequimolarPATlibrarypool',K.MOLECULAR);g.react('pool','Sourceintensity-guidedpoolbeforefacilitycleanup',S.POOL,((x,R.TEMPLATE),(other,R.TEMPLATE)),(pool,))
    x=clean(g,'facility_cleanup',pool,'Sequencingfacilitybeadcleanupofprimers')
    x=qc(g,'quantify',x,'Facilityquantitationaliquot;volumeunassigned')
    lib=g.product('library','PAT-seqphysicalDNApool;fullmixed-sugar/junction/readprimersunassembled',K.LIBRARY)
    g.react('library_identity','Physicalsequencinginterface',S.GENERIC,((x,R.TEMPLATE),),(lib,),changes_molecule=False)
    g.sequence(lib,SequencingPlan(limitation='MainactualHiSeq1500rapid100basesingleend;SIworkflowoverview150bpnotassignedtomainexperiment. Exactprimer/indexgeometryunassembled; no alignmentfilterphysicallyremovesmolecules.'))
    g.add_source_item(SourceItem('pending','HistoricalScriptSeqindex/completejunctionandmixedrTidentity','Privatepartialmodel;notpublishedabsence.'))
    return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph()
    c('Bound-DNAprimer retained as same material throughRT',g.steps['rt'].inputs[2].material,'end_extend')
    c('Capture separates actual unbound fraction',g.terminals['capture_unbound'].disposition,Disposition.DISCARD)
    c('Unusedtwo-thirds retained explicitly',g.terminals['remaining_cdna'].disposition,Disposition.UNRESOLVED)
    c('PCR receives selected physicalaliquot',g.steps['pcr'].inputs[0].material,'pcr_input')
    c('LibrarypoolusesproductionaftergelQCsplit',g.steps['pool'].inputs[0].material,'agarose_kept')
    c('Facilitycleanup followsphysicalpool',g.steps['facility_cleanup'].inputs[0].material,'pooled')
    c('Partialmodelcannotpublish',bool(g.publication_problems()),True)
    src=Source(sources_root()/'pat-seq__10.1261+rna.048355.114'/'rna.048355.114_supp_supp_048355.114_SuppMethods.docx.txt','get_sources --slugs pat-seq__10.1261+rna.048355.114')
    if have(c,src,label='publishedworkflow'):
        t=src.path.read_text();c('WholeBioEEtranscription',END_EXTEND in t,True);c('Whole5primeBtranscription',SPLINT_B in t,True);c('Wholeuniversaltranscription',UNIVERSAL in t,True);c('SourceexplicitRNAthymidine', '(rA)(rT)(rC)(rT)' in t,True)
    c.report()
if __name__=='__main__':selftest()
