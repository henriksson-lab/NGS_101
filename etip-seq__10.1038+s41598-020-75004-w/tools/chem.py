"""Private eTIP source-positive partial graph; no complete molecular/absence claim."""
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[1]/'lib'))
from chemdraw import Scene,Segment,revcomp
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,Terminal,Disposition,SequencingPlan,SourceItem
AUTHORITY='PublishedPMC7596052 Methods164–176/Introduction/Results/Fig3actualPDFp5 and modifiedordersPDFp11; publishedMOESM6notesmethodologyp1/FigS1p2. SelectedCGNday2eTIPb default; complementaryeTIPa separate actual preparation, one method.'
SOURCE_ORDERS={'top':'GTTGGATCCGATA[Bio-dT]CGC','bottom':'GGCCGCGATATCGGATCCAACT'}
TOP=SOURCE_ORDERS['top'].replace('[Bio-dT]','T')
BOTTOM=SOURCE_ORDERS['bottom']
ADAPTOR_CHEMICAL_STATE={'source_whole_orders':SOURCE_ORDERS.copy(),'backbone':'DNA','top_internal_biotin_dT_zero_based':13,'ordered_terminal_phosphates':'unassigned','local_duplex':'17bp core/bottom3T/bottom5GGCC','physical_identity':'same annealed adaptor stock; aliquoting does not alter ordering'}
SOURCE_OLIGO_EXCEPTIONS={'adaptor_top':'Published main174/actualPDFp11 prints GTTGGATCCGATA[Bio-dT]CGC: internal Bio-dT is ONE physical thymidine at zero-based13, not a missing linker/base. Generic source parser truncates at bracketed Bio-dT. Whole modified original fixture and exact physical-polymer restoration checked locally; narrow transcription-format exception only, not source absence or shared parser migration.'}

def strand_scene(which):
    sc=Scene()
    if which=='top':
        sc.strand('top',[Segment('before_biotin',TOP[:13]),Segment('physical_biotin_dT',TOP[13:14],note='Physical thymidine nucleotide with biotin, not an extra/linker-only base'),Segment('after_biotin',TOP[14:])])
    elif which=='bottom':sc.strand('bottom',[Segment('joining_5GGCC',BOTTOM[:4]),Segment('core_rc',BOTTOM[4:-1]),Segment('insert_3T',BOTTOM[-1:])])
    else:raise ValueError('Select source top or bottom')
    return sc

def adaptor_scene(end='insert'):
    if end not in ('insert','joining'):raise ValueError('Select insert or joining end')
    if BOTTOM[4:-1]!=revcomp(TOP) or BOTTOM[:4]!=revcomp(BOTTOM[:4]):raise ValueError('Source whole orders must pair with17bp core/self-complementaryGGCC')
    sc=Scene();sc.strand('top',[Segment('core',TOP,note='Source Bio-dT at zero-based13; exact phosphate state of ordered oligo unassigned')])
    sc.anneal('bottom',[Segment('joining_5GGCC',BOTTOM[:4]),Segment('core_rc',BOTTOM[4:-1]),Segment('insert_3T',BOTTOM[-1:])],to='top',pair=('core_rc','core'))
    sc.ligation_end('bottom','insert_3T' if end=='insert' else 'joining_5GGCC',kind='cohesive',polarity='3-prime' if end=='insert' else '5-prime',owner='adapter',text='Known source local3T/sourceDNA3A class; fullsealedmoleculesunassembled' if end=='insert' else 'Known source5GGCC self-compatible joining overhang; sourcePNK intended5P, occupancy/sealedjunctionunassembled')
    return sc

def unknown_scene():
    sc=Scene();sc.strand('unassembled',[Segment('unassembled','[unassembled]',placeholder=True)]);return sc

class Builder:
    def __init__(self,g,x):self.g=g;self.x=x
    def inp(self,key,label,kind=K.CONDITION,role=R.CONDITION,**kw):return(self.g.input(key,label,kind,**kw),role)
    def cat(self,key,label):return self.inp(key,label,K.CATALYST,R.CATALYST)
    def op(self,key,label,ports=(),kind=S.MAGIC,product=None,**kw):
        y=self.g.product(key+'_product',product or 'Actual source preparation; whole heterogeneous molecules/end distributions/conversion/byproduct fates unassembled',K.MOLECULAR)
        if kind==S.MAGIC:kw.setdefault('reason','Source positive operation; exact occupied substrates/full molecular products, conversion and off-products remain unassembled, not an absence claim.')
        if kind==S.GENERIC:kw['changes_molecule']=False
        self.g.react(key,label,kind,((self.x,R.TEMPLATE),)+tuple(ports),(y,),**kw);self.x=y;return y
    def cleanup(self,key,label,ports=()):
        y=self.g.product(key+'_kept','Actual retained workup fraction; recovery/carryover unassigned',K.MOLECULAR);w=self.g.product(key+'_excluded','Actual excluded workup fraction',K.WASTE)
        self.g.react(key,label,S.CLEANUP,((self.x,R.TEMPLATE),)+tuple(ports)+(self.inp(key+'_conditions',label),),(y,w));self.g.discard(w,'Source workup fraction excluded, not quantitative molecular purity.');self.x=y
    def allocate(self,key,label,selected,rest,qc=False):
        a=self.g.product(key+'_selected',selected,K.MOLECULAR);r=self.g.product(key+'_remainder',rest,K.MOLECULAR)
        self.g.react(key,label,S.SPLIT,((self.x,R.TEMPLATE),),(a,r))
        self.g.mark_terminal(Terminal(a.id if qc else r.id,Disposition.QC if qc else Disposition.UNRESOLVED,selected if qc else rest));self.x=r if qc else a;return self.x
    def components(self,key,label,items):
        stock=self.g.input(key,label,K.MOLECULAR);out=[];ports=[]
        for suffix,kind,role,text in items:
            y=self.g.product(key+'_'+suffix,text,kind);out.append(y);ports.append((y,role))
        self.g.react(key+'_components','Inventory decomposition of source supplies; not physical separation or extra fresh additions',S.GENERIC,((stock,R.TEMPLATE),),tuple(out),changes_molecule=False)
        return tuple(ports)

def protocol_graph():
    return preparation_graph('b')

def preparation_graph(branch='b'):
    if branch not in ('a','b'):raise ValueError('Select eTIPa or eTIPb, not a synthetic combined-IP split')
    g=ProtocolGraph('eTIP'+branch+'-seq selectedCGNday2 preparation',authority=AUTHORITY,representation='published_instance')
    x=g.input('cgn','Actual sourceprepared day2CGN culture100mm dish; upstreamanimal/culture outside graph',K.BIOLOGICAL,primary=True);b=Builder(g,x)
    g.add_source_item(SourceItem('fraction_identity','PublishedmainResults/Fig1 andSI6p1–2: P1mostlycovalentGsegments; P2containsnoncovalentTANDresealedG; incomplete preservation/ligationclasses.','No perfectspeciespurity/singlecomplex-onlyproduct orfullcrosslinkatommap. eTIPa/b are separateactualpreparations, notinventedaliquots ofoneIP.'))
    g.add_source_item(SourceItem('modification','Main174/actualPDFp11 topGTTGGATCCGATA[Bio-dT]CGC andbottomGGCCGCGATATCGGATCCAACT.','PhysicalBio-dT retainedasT; computed17bpcore/3T/5GGCC. Wholemodifiedfixturesverify narrowparserexception; orderedterminalphosphates unassigned, sourcePNKintendedclass explicit.'))
    g.add_source_item(SourceItem('filtering','Main188–190unique-mapping/adaptor-startreadfilter andResultsinter-complexchimerainference.','Computationalreadselectiondoesnotphysicallypurifypreparedmolecules; genomicdistanceisnotphysicalchimeraDNAlength orUMI.'))
    b.op('etoposide','Source0.5mMetoposide15minserumfreemedium',ports=(b.inp('etoposide_reagent','Actual0.5mMetoposide',K.MOLECULAR,R.STRUCTURAL),b.inp('etoposide_conditions','Actual15minsourceprogramme')),product='Sourceetoposide-trappedtopoIIbeta/DNAmixedcomplexclasses; occupancy/unreactedmaterialunassembled')
    b.op('lysis','Source1percentsarkosyl/23Gneedle/CsClfinal.5M',ports=(b.inp('lysis_conditions','Actualsarkosyl/CsCl sourceconditions'),))
    b.op('sonicate','Sourcemicroprobesonication/dilute3volTEST100',ports=(b.inp('sonication_conditions','ActualTEST100Tris7.5/EDTA100mMNaCl/Triton/PIC; no inventedsonicationdose'),))
    sup=g.product('ip_supernatant','SourceactualsupernatantforIP',K.MOLECULAR);other=g.product('nonselected_lysate','NonselectedmaterialbeforeIP; sourceyield/exactseparationunassigned',K.MOLECULAR)
    g.react('select_supernatant','Source selects actual supernatant; separating operation parameters unassigned',S.SPLIT,((b.x,R.TEMPLATE),),(sup,other));g.mark_terminal(Terminal(other.id,Disposition.UNRESOLVED,'No sourcequantityorfateassigned'));b.x=sup
    beads=g.input('ip_beads','Actualanti-topoIIbeta3B6coatedProteinGmagneticbeads; preparedselectioncomplex/factoryunassembled',K.SELECTION)
    bound=g.product('ip_bound','ActualselectedtopoIIbeta-DNAbeadassociatedpreparation',K.MOLECULAR);waste=g.product('ip_unbound','ActualunboundIPfraction',K.WASTE)
    g.react('ip','Sourceovernight4Canti-topoIIbetaProteinGIP',S.AFFINITY_CAPTURE,((b.x,R.TEMPLATE),(beads,R.SELECTION)),(bound,waste));g.discard(waste,'Excludedsourceunboundfraction, not complete recovery');b.x=bound;b.cleanup('ip_wash','Source3xTEST200wash')
    def recover(prefix):
        b.op(prefix+'_digest','SourceRNaseA/proteinaseK recovery; exactdosesunassigned',ports=(b.cat(prefix+'_rnase_a','ActualRNaseA'),b.cat(prefix+'_pk','ActualproteinaseK')),kind=S.DEPROTEINIZATION)
        b.cleanup(prefix+'_organic','Sourcephenol/chloroformextraction')
        b.cleanup(prefix+'_ethanol','Sourceethanolprecipitationwithglycogencarrier',ports=(b.inp(prefix+'_glycogen','Actualsourceglycogen molecularcarrier',K.MOLECULAR,R.STRUCTURAL),))
    def finish(prefix,cycles,is_a=False):
        if is_a:
            inv=b.components(prefix+'_kit','ActualTruSeqDNALTsupplies aggregateinventory; stockcount/packaging/formulationsunassigned',[
                ('adapters',K.OLIGO,R.STRUCTURAL,'Actualkitadapterclass, wholeoligos/modificationsunassembled'),('forward',K.OLIGO,R.STRUCTURAL,'ActualkitPCRforwardprimerclass; wholebodyunassembled'),('reverse',K.OLIGO,R.STRUCTURAL,'ActualkitPCRreverseprimerclass; wholebody/indexunassigned'),('enzymes',K.CATALYST,R.CATALYST,'ActualPCRpolymeraseclass; otherkitcatalysts/exactcompositionunassembled'),('nts',K.MOLECULAR,R.STRUCTURAL,'Actualkitnucleotidecomponentclasses, exactstocks/dosesunassigned'),('conditions',K.CONDITION,R.CONDITION,'Actualkitconditions unassembled')])
            b.op(prefix+'_library_pcr','SourceTruSeqDNALTpreparation/'+str(cycles)+'PCRcycles; no inheritedvendorrecipe',ports=inv)
        else:
            inv=b.components(prefix+'_endprep_supplies','Actualon-beadendrepair/dA-tail/paired-endadaptor supplies aggregateinventory; packaging/dosesunassigned',[
                ('adapters',K.OLIGO,R.STRUCTURAL,'Actualpaired-endsequencingadapterclass; wholebodies/ends/modsunassembled'),('catalysts',K.CATALYST,R.CATALYST,'Actualendrepair/tailing/ligationcatalystclasses; exactcomponentsunassigned'),('nts',K.MOLECULAR,R.STRUCTURAL,'Actualendprep/tailing nucleotideinventoryclass; exactstocks/dosesunassigned'),('conditions',K.CONDITION,R.CONDITION,'Actualsourceon-beadendprepconditionsunassembled')])
            b.op(prefix+'_library_adapter','Sourceon-beadendrepair/dA-tail/paired-endadapterligation',ports=inv)
            master=b.components(prefix+'_pcr_master','ONEactualIlluminaPCRmastermix; sourcevolumeunassigned',[
                ('enzyme',K.CATALYST,R.CATALYST,'PCRenzymecomponentofsameactualmaster'),('nts',K.MOLECULAR,R.STRUCTURAL,'PCRnucleotidecomponentofsameactualmaster'),('buffer',K.CONDITION,R.CONDITION,'PCRbuffercomponentofsameactualmaster')])
            primers=b.components(prefix+'_primer_cocktail','ONEactualIlluminaPCRprimercocktail; volumeunassigned',[
                ('forward',K.OLIGO,R.STRUCTURAL,'Actualcocktailforwardprimerclass, exactwholebodyunassembled'),('reverse',K.OLIGO,R.STRUCTURAL,'Actualcocktailreverseprimerclass, exactwholebody/indexunassembled')])
            b.op(prefix+'_pcr','Source'+str(cycles)+'PCRcycles',ports=master+primers,kind=S.PCR,geometry=PCRGeometry.unavailable(unknown_scene(),reason='SourceactualPCRprimercocktail, wholebodies/occupiedprimersites/fullPCRmoleculesunassembled'))
            b.cleanup(prefix+'_ampure','SourceAMPureXPpurification; ratio/elutionvolumeunassigned')
        b.allocate(prefix+'_bioanalyzer','ActualBioAnalyzerQCaliquot; amountunassigned','ActualanalyticalQCcontribution','SourcepostQCproductionstock',qc=True)
        b.allocate(prefix+'_loading','Actuallibraryallocationforsequencing; amountunassigned','Actualsequencingallocation','Untransferredlibrary stock/fateunassigned')
        lib=g.product(prefix+'_library','ActualsourceHiSeq2000PE100libraryclass; wholeduplex/index/readlandingunassembled',K.LIBRARY)
        g.react(prefix+'_load','SourceHiSeq2000pairedend100bpeach',S.GENERIC,((b.x,R.TEMPLATE),),(lib,),changes_molecule=False)
        g.sequence(lib,SequencingPlan(read_lengths={'Read1':100,'Read2':100},limitation='SourceHiSeq2000PE100; historicaladapter/index/wholeprimers/readlandingunassembled; UMI meaning not assigned in selected model. Mappingfilters computationalonly.'))
    if branch=='a':
        p1=g.product('p1','ActualTEST500-retainedbeadP1; mostlycovalentG, notpure',K.MOLECULAR);p2=g.product('p2','ActualTEST500releasedP2 includingTandresealedG, notpure',K.MOLECULAR)
        g.react('salt_fractionate','SourceTEST5004CphysicalP1/P2fractionation',S.CLEANUP,((b.x,R.TEMPLATE),b.inp('test500','ActualTEST500.5MNaCl')),(p1,p2))
        for prefix,fraction in [('p1',p1),('p2',p2)]:
            b.x=fraction;recover(prefix)
            b.op(prefix+'_covaris','SourceCovarisS2~200bp',ports=(b.inp(prefix+'_shear','ActualsourceCovarisconditions; exactdoseunassigned'),))
            b.allocate(prefix+'_library_input','SI6FigS1p2 reports485ng DNAusedforsequencing; actualkitinputdilution/allocationunassembled','Actualsource485ngsequencingpreparationclass; not a new reagent dose','UntransferredDNAstock/fateunassigned')
            finish(prefix,10,True)
    else:
        b.op('blunt','SourceT4DNApolymerasebluntingofshearedends',ports=(b.cat('blunt_t4pol','ActualT4DNApolymerase'),b.inp('blunt_conditions','Sourcebluntingconditions/dosesunassembled')),kind=S.END_REPAIR)
        tail=b.components('tail_module','ActualNEBNextdA-tailingmodule; quantity/formulationunassigned',[
            ('klenow',K.CATALYST,R.CATALYST,'SourceKlenowfragment inactualmodule'),('nucleotide',K.MOLECULAR,R.STRUCTURAL,'dA-tailingnucleotideclass ofactualsourceoperation; stock/doseunassigned'),('buffer',K.CONDITION,R.CONDITION,'SourceNEBreactionbuffer/PICactualco-components')])
        b.op('tail','Source3dA-tailingonbeads',ports=tail,kind=S.TAILING)
        at=g.input('adaptor_top','Sourcewholemodifiedtopoligo, terminalorderedphosphateunassigned',K.OLIGO,state=strand_scene('top'),chemical_state={'source_ordering':SOURCE_ORDERS['top'],'backbone':'DNA','internal_biotin_dT_zero_based':13,'five_prime_ordered_modification':'unassigned'})
        ab=g.input('adaptor_bottom','Sourcewholebottomoligo; notabarcodeorUMI',K.OLIGO,state=strand_scene('bottom'),chemical_state={'source_ordering':SOURCE_ORDERS['bottom'],'backbone':'DNA','five_prime_ordered_modification':'unassigned'})
        ad=g.product('adaptor_duplex','Actualequimolarboil5min/slowcoolDNAadaptorstock; nucleotidelocalgeometryknown',K.OLIGO,state=adaptor_scene(),chemical_state=ADAPTOR_CHEMICAL_STATE.copy())
        g.react('adaptor_anneal','Sourceequimolaroligosboilingwater5min/slowcoolRT',S.BINDING,((at,R.TEMPLATE),(ab,R.STRUCTURAL),b.inp('anneal_conditions','Actualsourceannealprogramme')),(ad,))
        aq=g.product('adaptor_qc','ActualPAGEannealingQCaliquot',K.MOLECULAR);used=g.product('adaptor_ligation_stock','Source1.5uMadaptorligationallocation; amountunassigned',K.OLIGO,state=adaptor_scene(),chemical_state=ADAPTOR_CHEMICAL_STATE.copy());rest=g.product('adaptor_remainder','Actualuntransferredadaptorstock; amountunassigned',K.OLIGO,state=adaptor_scene(),chemical_state=ADAPTOR_CHEMICAL_STATE.copy())
        g.react('adaptor_allocations','ActualPAGEQC/adaptorligationallocation/remainder, amountsunassigned',S.SPLIT,((ad,R.TEMPLATE),),(aq,used,rest));g.mark_terminal(Terminal(aq.id,Disposition.QC,'ActualsourcePAGEannealingcheck'));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Noinventedadaptorstockvolume'))
        b.op('first_ligation','Source1.5uMadaptor/12000UT4ligase/5percentPEG/22C6h',ports=((used,R.STRUCTURAL),b.cat('first_t4ligase','Actual12000UT4DNAligasefirstreaction'),b.inp('first_ligation_conditions','Source5percentPEG8000/NEBbuffer/PIC')),kind=S.LIGATION,geometry=AdapterLigationGeometry.known('cohesive',adaptor_scene(),statement='SourceDNA3A/sourceadapterbottom3T class known; actualgenomicends/Poccupancy/sealedwholejunctionunassembled'))
        pnk_product=g.product('pnk_product','Sourceintended5P-adaptorendclasseswithinbeadcomplex; fulloccupancy/unreactedendsunassembled',K.MOLECULAR)
        joining_ends=g.product('copresent_joining_ends','SAMEsourceadaptorendinventorycopresentinPNKproduct; nofreshadditionorphysicalseparation',K.OLIGO,chemical_state={'copresent_with':'pnk_product','source_class':'intended5P/GGCCadaptorends','occupancy':'unassembled'})
        g.react('pnk','SourceT4PNKphosphorylatesadaptor5endsAFTERfirstligation',S.PHOSPHORYLATION,((b.x,R.TEMPLATE),b.cat('pnk_enzyme','ActualT4PNK'),b.inp('pnk_conditions','ActualsourcePNKconditions/doseunassigned')),(pnk_product,joining_ends));b.x=pnk_product
        b.op('second_ligation','SourcefreshT4DNAligase16Covernightrotating adaptorjoining',ports=((joining_ends,R.STRUCTURAL),b.cat('second_t4ligase','ActualsecondT4DNAligasereaction; doseunassigned'),b.inp('second_ligation_conditions','Actual16Covernightcontinuousrotation')),kind=S.LIGATION,geometry=AdapterLigationGeometry.known('cohesive',adaptor_scene('joining'),statement='KnownselfcompatibleGGCC5overhang/sourcePNKintended5P; wholeoccupiedsame/intercomplexsealedproductsunassembled'))
        b.op('stop','Sourceice10mMEDTAstop',ports=(b.inp('stop_conditions','Actual10mMEDTAonice'),),kind=S.GENERIC)
        recover('b_recover')
        b.op('terminal_biotin_removal','SourceT4DNApolymerasewithdATP/dTTPremovesterminal-adaptorbiotinylatednucleotides',ports=(b.cat('terminal_t4pol','ActualfreshT4DNApolymeraseoperation'),b.inp('terminal_dntp','ActualdATPanddTTP',K.MOLECULAR,R.STRUCTURAL),b.inp('terminal_conditions','Sourceincubationdose/timeunassembled')),product='Sourceintendedterminal-adaptorbiotinremoval/internalbiotinretainedclasses; exactbackbone/byproducts/offproductsunassembled')
        b.op('b_covaris','SourceCovarisshearing~200bp',ports=(b.inp('b_shear_conditions','ActualCovarisprogrammeunassigned'),))
        sb=g.input('streptavidin_beads','ActualDynabeadsMyOneStreptavidinC1',K.SELECTION);yes=g.product('biotin_captured','Actualstreptavidinbead-associatedDNA; notpurejunction-only',K.MOLECULAR);no=g.product('biotin_unbound','ActualexcludedunboundDNAfraction',K.WASTE)
        g.react('biotin_capture','SourcebiotinylatedfragmentcaptureAFTERterminalprocessing/Covaris',S.AFFINITY_CAPTURE,((b.x,R.TEMPLATE),(sb,R.SELECTION)),(yes,no));g.discard(no,'Actualsourceunboundfractionexcluded; recoveryunassigned');b.x=yes;b.cleanup('biotin_wash','Source3xTweenwash')
        finish('b',15)
    return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph();a=preparation_graph('a')
    c('ModifiedTisONEphysicalbase',len(TOP),17)
    c('Bio-dTsitephysicalT',TOP[13],'T')
    c('Wholemodifiedsourcepolymer',SOURCE_ORDERS['top'].replace('[Bio-dT]','T'),TOP)
    c('Wholebottomsource',len(BOTTOM),22)
    c('Literalcorepairs',revcomp(TOP),BOTTOM[4:-1]);adaptor_scene();adaptor_scene('joining')
    c('SourceGGCCjoiningoverhang',BOTTOM[:4],revcomp(BOTTOM[:4]))
    c('Sourcebottom3T',BOTTOM[-1],'T')
    c.raises('Rejectunknownroute',lambda:preparation_graph('synthetic'))
    c('Firstligusesactualadaptorallocation',g.steps['first_ligation'].inputs[1].material,'adaptor_ligation_stock')
    c('PNKAFTERfirstligation',g.steps['pnk'].inputs[0].material,'first_ligation_product')
    c('Secondligfreshenzyme',g.steps['second_ligation'].inputs[2].material,'second_t4ligase')
    c('Secondligusescopresentadaptorendsnotfreshstock',g.steps['second_ligation'].inputs[1].material,'copresent_joining_ends')
    c('JoiningendinventoryproducedinactualPNK','copresent_joining_ends' in g.steps['pnk'].outputs,True)
    c('NoTEST500inselectedB','salt_fractionate'in g.steps,False)
    c('TerminalprocessingbeforeCovaris',g.steps['b_covaris'].inputs[0].material,'terminal_biotin_removal_product')
    c('CaptureactualCovarisproduct',g.steps['biotin_capture'].inputs[0].material,'b_covaris_product')
    c('Aphysicalfractionoutputs',tuple(a.steps['salt_fractionate'].outputs),('p1','p2'))
    c('P2recoveryusesactualeluate',a.steps['p2_digest'].inputs[0].material,'p2')
    c('P1andP2notinventedpool',any(s.kind==S.POOL for s in a.steps.values()),False)
    c('BsourcePCRtyped',g.steps['b_pcr'].kind,S.PCR)
    c('ActualPCRmastercomponents',g.steps['b_pcr_master_components'].inputs[0].material,'b_pcr_master')
    c('Actualprimercocktailcomponents',g.steps['b_primer_cocktail_components'].inputs[0].material,'b_primer_cocktail')
    c('PAGEQCnotproduction',g.terminals['adaptor_qc'].disposition,Disposition.QC)
    p=sources_root()/HERE.parent.name/'s41598-020-75004-w_PMC7596052.xml.txt'
    if have(c,Source(p,'https://doi.org/10.1038/s41598-020-75004-w')):
        text=p.read_text()
        for key,raw in SOURCE_ORDERS.items():c('Wholeoriginalmodifiedorder '+key,raw in text,True)
    c('Privatepartial',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
