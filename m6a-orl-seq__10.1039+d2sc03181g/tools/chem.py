"""Source-qualified, private partial m6A-ORL-Seq preparation.

Residue markers are qualitative source states, not an atom map or yield model.
"""
from dataclasses import dataclass
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / 'lib'))
from chemdraw import Scene, Segment, revcomp, feature
from illumina import INDEX2_PRIMER_RC, INDEX1_PRIMER, TRUSEQ_READ2
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_graph import (ProtocolGraph, MaterialKind as K, InputRole as R,
    StepKind as S, Terminal, Disposition, SequencingPlan, SourceItem)

AUTHORITY = 'Published PMC9600483 mainFig1/Methods73–121; publishedESI10S20–S22/TableS1. ThirdanalyticalXLSXmetadataonly; closed3originalpass.'
RNA_LINKER = INDEX2_PRIMER_RC[:19]
RNA_LINKER_LENGTH = len(RNA_LINKER)
RT_PRIMER = 'A' + revcomp(RNA_LINKER)
DNA_ADAPTER_BODY = 'A' + INDEX1_PRIMER[:21]
DNA_ADAPTER_LENGTH = 10 + len(DNA_ADAPTER_BODY)
QPCR_FORWARD = TRUSEQ_READ2[-22:]

@dataclass(frozen=True)
class ResidueMarker:
    name: str = 'm6A'
    def __post_init__(self):
        if self.name not in ('m6A', 'Nm6A', 'Am6A', 'Bm6A', 'Tm6A'):
            raise ValueError('Only source-named qualitative residue states are declared')
    @property
    def biotin_on_rna(self): return self.name == 'Bm6A'
    @property
    def retained_thiolated_adduct(self): return self.name == 'Tm6A'
    def undergo(self, operation):
        transition = {('m6A','nitrosate'): 'Nm6A', ('Nm6A','reduce'): 'Am6A',
                      ('Am6A','label_BA'): 'Bm6A', ('Bm6A','cleave_disulfide'): 'Tm6A'}
        if (self.name, operation) not in transition:
            raise ValueError('Operation does not match the selected source marker pathway')
        return ResidueMarker(transition[self.name, operation])

def marker_pathway():
    states = [ResidueMarker()]
    for operation in ('nitrosate', 'reduce', 'label_BA', 'cleave_disulfide'):
        states.append(states[-1].undergo(operation))
    return states

def released_rna_descriptor():
    marker=marker_pathway()[-1]
    return {'selected_marker':marker.name,
            'retained_thiolated_rna_adduct':marker.retained_thiolated_adduct,
            'biotin_on_selected_released_rna':marker.biotin_on_rna,
            'scope':'Selected source released class only; whole RNA/atom map and other classes unassembled'}

def rna_linker_scene():
    sc=Scene(); sc.strand('linker', [Segment('RNA3 linker DNA body', RNA_LINKER)],
                         mod5='rApp activation', mod3='SpC3 non-nucleotide blocker'); return sc

def dna_adapter_scene():
    sc=Scene(); sc.strand('adapter', [Segment('source random N10', 'NNNNNNNNNN', placeholder=True,
        feature=feature('m6aorl.fragment_barcode', 'inline_barcode', 'random', note='Source fragment barcode; no UMI/deduplication/cell/sample meaning assigned')),
        Segment('DNA5 adaptor body', DNA_ADAPTER_BODY)], mod5='P', mod3='SpC3 non-nucleotide blocker'); return sc

def rt_primer_scene():
    sc=Scene(); sc.strand('RT', [Segment('source unpaired5A', 'A'), Segment('RT complement', revcomp(RNA_LINKER))]); return sc

def rt_binding_scene():
    """Postligation partial template: activation AMP has already left the donor."""
    sc=Scene()
    sc.strand('template', [Segment('unassembled upstream RNA', '[unassembled RNA]', placeholder=True, inferred=True),
                           Segment('RNA3 linker DNA body', RNA_LINKER)],
              mod3='SpC3 non-nucleotide blocker')
    sc.anneal('RT', [Segment('source unpaired5A', 'A'), Segment('RT complement', revcomp(RNA_LINKER))],
        to='template', pair=('RT complement', 'RNA3 linker DNA body'), unpaired=('source unpaired5A',))
    return sc

def rna3_mixture_descriptor():
    return {'copresent_byproducts':'Released activation AMP plus unreacted linker/RNA/off-products; amounts and fates unassigned',
            'scope':'One physical reaction mixture until source workup; AMP is separate from linker backbone, not physically separated'}

def qforward_scene():
    sc=Scene(); sc.strand('FP', [Segment('source qPCR FP', QPCR_FORWARD)]); return sc

def inferred_acceptor_scene(kind):
    if kind not in ('RNA','cDNA'): raise ValueError('Only the selected RNA/cDNA acceptor classes are represented')
    sc=Scene(); sc.strand('acceptor', [Segment('unassembled '+kind, '[unassembled '+kind+']', placeholder=True, inferred=True)],
                          mod3='INFERRED intended3primeOH')
    sc.ligation_end('acceptor', 'unassembled '+kind, kind='single-strand', polarity='3-prime', owner='substrate', text='INFERRED intended3OH acceptor class')
    return sc

def ligation_geometry(which):
    if which == 'RNA3':
        sc=rna_linker_scene(); strand='linker'; first='RNA3 linker DNA body'; kind='RNA'
        statement='Source19nt preadenylated donor5rApp/3SpC3. IntendedRNA3OH population inferred/unassembled; activation AMP is a leaving group, not backbone nucleotide.'
    elif which == 'DNA5':
        sc=dna_adapter_scene(); strand='adapter'; first='source random N10'; kind='cDNA'
        statement='Source32nt donor5P/3SpC3 with randomN10 fragmentbarcode; intendedcDNA3OH inferred/unassembled. SpC3 is a non-nucleotide blocker.'
    else: raise ValueError('Unknown source ligation branch')
    sc.ligation_end(strand, first, kind='single-strand', polarity='5-prime', owner='adapter')
    return AdapterLigationGeometry.known('single-strand', sc, inferred_acceptor_scene(kind), statement=statement)

def unavailable_scene():
    sc=Scene(); sc.strand('unknown', [Segment('whole products', '[unassembled products]', placeholder=True)]); return sc

def sequence_of(sc):
    return ''.join(s.top for strand in sc.strands.values() for s in strand.segs)

def protocol_graph():
    g=ProtocolGraph('m6A-ORL-Seq selected5ug293TpolyA branch', authority=AUTHORITY, representation='published_instance')
    x=g.input('polya_rna', 'Source5ug293TpolyA RNA; actual occupied sequences/m6A/m6Am and end populations unassigned', K.MOLECULAR, primary=True)
    def reagent(key,label,kind=K.CONDITION,role=R.CONDITION,**kw): return g.input(key,label,kind,**kw),role
    def op(key,label,kind,inputs=(),product='Source mixture including unreacted/background/off-products; whole molecules/occupancy/yield unassembled', chemical_state=None, **kw):
        nonlocal x
        y=g.product(key+'_product',product,K.MOLECULAR,chemical_state=chemical_state)
        if kind==S.GENERIC:kw['changes_molecule']=False
        g.react(key,label,kind,((x,R.TEMPLATE),)+tuple(inputs),(y,),**kw);x=y
    def cleanup(key,label):
        nonlocal x
        y=g.product(key+'_kept','Actual retained RNA/DNA fraction; recovery/carryover not class-pure',K.MOLECULAR)
        w=g.product(key+'_waste','Actual discarded workup fractions; imperfect molecular separation',K.WASTE)
        g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),reagent(key+'_conditions',label)),(y,w));g.discard(w,'Physical workup, not computational filtering');x=y
    def transfer(key,label):
        nonlocal x
        y=g.product(key+'_used',label,K.MOLECULAR);r=g.product(key+'_remainder','Untransferred stock; actual subsequent fate unassigned',K.MOLECULAR)
        g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(y,r));g.mark_terminal(Terminal(r.id,Disposition.UNRESOLVED,'Remainder not automatically discarded'));x=y
    def qc(key,label,extra=()):
        nonlocal x
        y=g.product(key+'_stock','Production stock after actual QC aliquot allocation; amount unassigned',K.MOLECULAR);a=g.product(key+'_aliquot',label,K.MOLECULAR)
        g.react(key+'_split',label,S.SPLIT,((x,R.TEMPLATE),),(y,a))
        if extra:
            z=g.product(key+'_assay','Source QC assay reaction; not production library',K.QC)
            g.react(key+'_measure',label,S.MAGIC,((a,R.TEMPLATE),)+tuple(extra),(z,),reason='QC only; full assay molecular states/landing unassembled.')
            g.mark_terminal(Terminal(z.id,Disposition.QC,label))
        else:g.mark_terminal(Terminal(a.id,Disposition.QC,label+'; assay internals unassembled'))
        x=y

    op('fragment','Source70C15minRNAfragmentation60–200nt',S.MAGIC,(reagent('fragmentation_reagent','Source ThermoRNAfragmentation prepared reagent',K.MOLECULAR,R.STRUCTURAL),),reason='Sourcefragmentation positive; exact cleavages/ends/populations unassembled.')
    op('nitrosate','SourceNOexposure/3Macetate4.5/2mMEDTA;ice6min then4C20h',S.CHEMICAL_CONVERSION,
       (reagent('no_reactant','SourceNO/O2 chemistry reactants, prepared prior to graph',K.MOLECULAR,R.STRUCTURAL),reagent('nitrosation_conditions','Source200ulRNA nitrosation conditions')),
       product='SourceNm6A-bearing RNA classes plus unreactedRNA and AtoI/GtoX/CtoU side-change classes; not nucleoside88percent-yield assignment')
    cleanup('nitrosation_rcc','SourceRCC after nitrosation')
    op('tdo_reduce','Source100mMTDO/300mMcarbonate10/2mMEDTA/30ul;50C10min1000rpm',S.CHEMICAL_CONVERSION,
       (reagent('tdo','SourcepreparedTDO consumed reactant, not catalytic default',K.MOLECULAR,R.STRUCTURAL),reagent('tdo_conditions','SourceRNA-specific30ul reduction conditions')),
       product='SourceAm6A marker class plus unreacted/otherRNA populations; wholeproducts/yields unassigned')
    cleanup('tdo_rcc','SourceRCC after TDO reduction')
    op('ba_label','Source2mMBA/100mMMES5.5/10percentDMF50ul;37C1h1000rpm',S.CHEMICAL_CONVERSION,
       (reagent('ba_probe','PreparedBA physical biotinylated aldehyde reactant withdisulfide; manufactureoutsidegraph',K.MOLECULAR,R.STRUCTURAL),reagent('ba_conditions','Source50ul labelingconditions')),
       product='SourceBm6A biotin-bearing RNA classes plus unreactedRNA/m6Am/background/offproducts; no exclusivem6A specificity claim')
    cleanup('ba_rcc','SourceRCC after BA labeling')
    op('carrier','Sourceadd50ugYeasttRNA AFTERchemicaltreatments',S.POOL,(reagent('yeast_trna','Actual50ugYeasttRNA carrier, not sequenced-target identity',K.MOLECULAR,R.TEMPLATE),))
    beads=g.input('c1_beads','Actual25ulC1streptavidinbeads:two200ulNaOH100mM/NaCl50mMprewashes/twowaterwashes/resuspend100ulbindingmixture',K.SELECTION)
    op('capture','SourceC1binding1h thenmagnet; bindingmixture10ulhighsaltbuffer+90ulwater',S.MAGIC,((beads,R.SELECTION),reagent('binding_conditions','SourceTris/EDTA/NaCl/Tweenbindingconditions')),
       reason='Sourcebiotin-affinitybinding positive; whole occupiedbead/RNA classes unassembled.')
    cleanup('bead_wash','Source three1ml highsaltwashes, retainbeads; unbound/washfractions physically removed')
    # Consecutive releases consume the SAME residual bead fraction; supernatants pool later.
    eluate1=g.product('dtt_eluate1','Actualfirst50ulDTTsupernatant:Tm6Athiolated-adductRNA plus possiblecarryover',K.MOLECULAR,chemical_state=released_rna_descriptor())
    residual1=g.product('dtt_residual1','Residualbeads including biotin-side linker and incompletelyreleasedRNA classes',K.MOLECULAR)
    g.react('dtt_release1','Source50ul100mMDTT37C15min1000rpm; collectsupernatant',S.CHEMICAL_CONVERSION,
        ((x,R.TEMPLATE),reagent('dtt1','Actualfirst50ul100mMDTT reducingreactant',K.MOLECULAR,R.STRUCTURAL),reagent('dtt1_conditions','Source37C15min1000rpmfirstrelease')),(eluate1,residual1))
    eluate2=g.product('dtt_eluate2','Actualsecond50ulDTTsupernatant:thiolated-adductRNA andpossiblecarryover',K.MOLECULAR,chemical_state=released_rna_descriptor())
    residual2=g.product('dtt_residual2','Actualremainingbeads/biotin-side material and incomplete-releaseclasses; fateunassigned',K.MOLECULAR)
    g.react('dtt_release2','Source second50ul100mMDTT50C5min onresidualbeads; collectsupernatant',S.CHEMICAL_CONVERSION,
        ((residual1,R.TEMPLATE),reagent('dtt2','Actualsecond50ul100mMDTT reducingreactant',K.MOLECULAR,R.STRUCTURAL),reagent('dtt2_conditions','Source50C5minsecondrelease')),(eluate2,residual2))
    g.mark_terminal(Terminal(residual2.id,Disposition.UNRESOLVED,'Beads excluded from collectedRNA; notbiotin-streptavidinbondcleavage or restorationofunmodifiedm6A'))
    pooled=g.product('pooled_eluates','Actual combinedfirstandsecondsupernatants; Tm6A-bearing classes retainthiolatedadduct',K.MOLECULAR,chemical_state=released_rna_descriptor())
    g.react('eluate_pool','Physicallycombine TWOactualDTTeluates',S.POOL,((eluate1,R.TEMPLATE),(eluate2,R.TEMPLATE)),(pooled,));x=pooled
    cleanup('enriched_rcc','SourceRCC pooledeluates')
    qc('rna_qc','Source Qubit4 analyticalaliquot')
    transfer('rna_library_transfer','Source50ngRNA aliquot forlibrary, prepared5ulRNA forPNK; untransferredstockunassigned')
    op('pnk','Source10ul:5RNA/1buffer/1PNK/1inhibitor/2water;37C1h65C20min',S.MAGIC,
       (reagent('pnk_enzyme','Source1ulPNK catalyst',K.CATALYST,R.CATALYST),reagent('pnk_conditions','SourcePNKbuffer/inhibitor/water; NOATPaddition listed')),
       reason='SourcePNKoperation positive; do not assert creationofRNA5P withoutATP. Intended3OHacceptorpopulation/wholeends unassembled.')
    linker=g.input('rna3_linker','Source1ul10uM preadenylated19ntDNA3linker/3SpC3',K.ADAPTER,state=rna_linker_scene(),chemical_state={'five_prime':'rApp','three_prime':'SpC3','backbone_length':RNA_LINKER_LENGTH})
    lig=g.product('rna3_ligated','OnephysicalRNA3ligationmixture: intendedRNA-linkerproducts pluscopresentreleasedactivationAMP/unreactedlinker/RNA/offproducts; wholejunctionunassembled',K.MOLECULAR,chemical_state=rna3_mixture_descriptor())
    g.react('rna3_ligate','SourceRNA10+linker1/preheat65C5minice+mix14.5→25.5ul;25C2h16Covernight',S.LIGATION,
       ((x,R.TEMPLATE),(linker,R.STRUCTURAL),reagent('rnl2','Source1ulT4RNAligase2truncatedKQ',K.CATALYST,R.CATALYST),reagent('rna3_conditions','Source2.5buffer/2.5DTT/7.5PEG/1inhibitor; noadditionalATPlisted')),
       (lig,),geometry=ligation_geometry('RNA3'));x=lig
    op('deadenylase','Source1ul5primeDeadenylase30C1h beforeRecJf',S.MAGIC,
       (reagent('deadenylase_enzyme','Source5primeDeadenylase',K.CATALYST,R.CATALYST),),reason='Sourceexcessadapterdeactivation positive; actualcompleteoccupancy/offproductstatesunassembled.',chemical_state=rna3_mixture_descriptor())
    op('recjf','Source1ulRecJf37C1h AFTERdeadenylase',S.MAGIC,
       (reagent('recjf_enzyme','SourceRecJf',K.CATALYST,R.CATALYST),),reason='Sourceexcesslinkerdigestionpositive; allsubstrateclasses/specificity/yieldunassembled.',chemical_state=rna3_mixture_descriptor())
    cleanup('rna3_occ','SourceOCC cleanup/elute12ul; copresentAMP/linker/byproductpartitionandrecoveryunassigned')
    transfer('rt_transfer','Source10ul from12ulRNAeluate')
    rt=g.input('rt_primer','Source1ul2uMRTprimer20nt',K.OLIGO,state=rt_primer_scene())
    op('rt','SourceRNA10+primer1/preheat65C5minice+mix9→20ul;25C3min42C7min52C30min4Chold',S.MAGIC,
       ((rt,R.STRUCTURAL),reagent('ssiii','Source1ulSSIII',K.CATALYST,R.CATALYST),reagent('rt_nts','Actual2ul10mM-eachdNTP',K.MOLECULAR,R.STRUCTURAL),reagent('rt_conditions','Source4ulbuffer/1DTT/1inhibitor; immediatefollowing≤37Chandlingwarning scopeambiguous')),
       reason='SourceRTstops/copyingpositive; 19bpprimer-linkerlocalduplex known, wholeRNA/cDNA/stoppositionsunassembled.')
    op('exoi','Source1ulExoI37C30min afterRT',S.MAGIC,(reagent('exoi_enzyme','SourceExonucleaseI',K.CATALYST,R.CATALYST),),reason='SourceRTprimerdigestion positive; actualfullsubstrate/offproductpopulationsunassembled.')
    op('rna_hydrolysis','Source+10ul1MNaOH+10ul.5MEDTA+10ulwater65C15min',S.MAGIC,
       (reagent('hydroxide','SourceNaOHhydrolysisreactant',K.MOLECULAR,R.STRUCTURAL),reagent('hydrolysis_conditions','SourceNaOH/EDTA/waterheatingconditions')),
       reason='SourceRNAhydrolysispositive; physicalRNAadducts distinctfromDNAcopies/stops, wholeproductsunassembled.')
    cleanup('cdna_occ','SourceOCC/elute7ul')
    transfer('dna5_transfer','Source6ulcDNA from7uleluate')
    dna5=g.input('dna5_adapter','Source.3ul80uM5P/N10/22ntbody/3SpC3DNAadaptor',K.ADAPTER,state=dna_adapter_scene(),chemical_state={'five_prime':'P','three_prime':'SpC3','backbone_length':DNA_ADAPTER_LENGTH})
    op('dna5_ligate','Source6cDNA+.3adaptor+1DMSO/preheat65C5minice+mix12.7→20ul;25Covernight',S.LIGATION,
       ((dna5,R.STRUCTURAL),reagent('rnl1','Source1.5ulRNAligase1',K.CATALYST,R.CATALYST),reagent('dna5_atp','Actual.2ul.1MATP',K.MOLECULAR,R.STRUCTURAL),reagent('dna5_conditions','Source2ulbuffer/9ul50percentPEG8000/1ulDMSO')),
       geometry=ligation_geometry('DNA5'),product='IntendedcDNA3OH-to-DNAadaptor5Pproducts withN10fragmentbarcode, freeadaptor/nativecDNA/offproducts; wholejunctionunassembled')
    cleanup('dna5_occ','SourceOCC/elute30ul')
    f=g.input('qpcr_forward','Source1ul10uMFP-Library-qPCR',K.OLIGO,state=qforward_scene())
    r=g.input('qpcr_reverse','Source1ul10uMRP-Library-qPCR equalsRTprimer',K.OLIGO,state=rt_primer_scene())
    qm=g.input('qpcr_master','Actualprepared5ulHieffPCRSYBRGreenMasterMixstock',K.MOLECULAR)
    qe=g.product('qpcr_enzyme','Enzymecomponentofsame5ulSYBRmasterstock; identity/doseunassigned',K.CATALYST)
    qn=g.product('qpcr_nts','Nucleotidecomponentofsame5ulSYBRmasterstock; doseunassigned',K.MOLECULAR)
    qd=g.product('qpcr_dye','SYBRGreenfluorescentcomponentofsame5ulmasterstock; doseunassigned',K.MOLECULAR)
    qb=g.product('qpcr_buffer','Buffercomponentofsame5ulSYBRmasterstock; compositionunassigned',K.CONDITION)
    g.react('qpcr_components','Representcopresentcomponentsofoneactual5ulpreparedSYBRmasterstock, noindependentadditions',S.GENERIC,((qm,R.TEMPLATE),),(qe,qn,qd,qb),changes_molecule=False)
    qc('cycle_qc','Actual3ulof30ulcDNA allocated to10ulqPCR withFP/RP and5ulSYBRmasterstock',
       ((f,R.STRUCTURAL),(r,R.STRUCTURAL),(qe,R.CATALYST),(qn,R.STRUCTURAL),(qd,R.STRUCTURAL),(qb,R.CONDITION)))
    transfer('pcr_transfer','Actual23ulfromremaining27ul; nominal4ul untransferred')
    master=g.input('q5_master','Actual25ulNEBNextUltraIIQ5MasterMixM0544L preparedstock',K.MOLECULAR)
    enz=g.product('q5_enzyme','Catalystcomponentofsame25ulmasterstock',K.CATALYST)
    nts=g.product('q5_nts','Nucleotidecomponentofsame25ulmasterstock; doseunassigned',K.MOLECULAR)
    buf=g.product('q5_buffer','Buffercomponentofsame25ulmasterstock',K.CONDITION)
    g.react('q5_components','Componentrepresentationofactualmasterstock, notindependentfreshadditions',S.GENERIC,((master,R.TEMPLATE),),(enz,nts,buf),changes_molecule=False)
    pf=g.input('neb_universal','Source1ulNEBuniversalprimer; wholebodyunassembled',K.OLIGO)
    pr=g.input('neb_index','Source1ulNEBindexprimer; wholebody/indexchoiceunassembled',K.OLIGO)
    op('pcr','Source50ul:23cDNA/25master/1universal/1index;98C30s,98C10s65C75scyclescuedbyqPCR,65C5min4Chold',S.PCR,
       ((enz,R.CATALYST),(nts,R.STRUCTURAL),(buf,R.CONDITION),(pf,R.STRUCTURAL),(pr,R.STRUCTURAL)),
       geometry=PCRGeometry.unavailable(unavailable_scene(),reason='WholeNEBprimerbodies/indexchoices/occupiedtemplates/fullcopiesunassembled. Native3SpC3blockerandRNAadductare notassignedtonewPCRcopies.'))
    cleanup('pcr_clean','Source45ulVAHTSbeads on50ulPCR(.9x),two80percentEtOHwashes/elute23ul.1xTE')
    transfer('final_transfer','Source21ulfrom23uleluate')
    transfer('loading','Actualloadingaliquot; untransferredlibrarystockfateunassigned')
    lib=g.product('library','Selectedm6AORLHiSeqXTenlibrary; wholemolecules/index/readlandingunassembled',K.LIBRARY)
    g.react('load','SourceHiSeqXTenPE150',S.GENERIC,((x,R.TEMPLATE),),(lib,),changes_molecule=False)
    g.sequence(lib,SequencingPlan(limitation='SourcePE150; exactNEBprimer/index/runmask/landingunassembled. N10fragmentbarcode notassumedUMI/cell/sample. N10sourceversushexamerprocessingdiscrepancy unresolved; analysisfiltersnotphysicalcleanup.'))
    g.add_source_item(SourceItem('barcode_discrepancy','MainDataPreprocessing119 says randomdecamer(N10), then says trimrandomhexamer; ESI10dprintsN10.','Physicalmodel retainsN10randominlinefragmentbarcode only; no undocumented deduplication/indexinterpretation or N6replacement.'))
    g.add_source_item(SourceItem('temperature_scope','ESI10cprogram25/42/52/4C is immediatelyfollowed by donotraiseabove37C handlingwarning.','Programretained; warninghandling-scopeambiguous, notrewrittenasglobaltemperaturecap.'))
    g.add_source_item(SourceItem('chemistry_limits','MainFig1Bm6A→Tm6A DTTdisulfidecleavage retainsRNAthiolatedadduct; sidechanges/unreacted/m6Amclassesremain.','Noatommap/yieldprediction/exclusivem6Aspecificity; comparison/synthesis/predecessor/vendorbranchesnotmerged.'))
    return g

def selftest():
    from checks import Check, Source, have
    from storage import sources_root
    c=Check();g=protocol_graph();root=sources_root()/HERE.parent.name
    p=root/'d2sc03181g_supp_SC-013-D2SC03181G-s002.pdf.txt'
    if have(c,Source(p,'https://doi.org/10.1039/d2sc03181g')):
        t=''.join(p.read_text().split()).upper()
        for name,ordering in [('RNA3modified','5’rAPP-'+RNA_LINKER+'-3SpC3'),('DNA5modified','5Phos-NNNNNNNNNN'+DNA_ADAPTER_BODY+'-3SpC3'),('RT',RT_PRIMER),('qPCRFP',QPCR_FORWARD)]:
            c('WholeESIsourceordering '+name,ordering.upper() in t,True)
    c('RT localpairing19bpwith5Aoverhang',len(rt_binding_scene().strands['RT'].segs[-1].top),19)
    c('BoundRTtemplatehasreleasedactivationAMP',rt_binding_scene().strands['template'].mod5,'')
    c('RNA3C3blockeraddsnoC',len(sequence_of(rna_linker_scene())),19)
    c('DNA5C3blockeraddsnoC',len(sequence_of(dna_adapter_scene())),32)
    c('Sourcebarcodeonlyinline',dna_adapter_scene().strands['adapter'].segs[0].feature.role,'inline_barcode')
    c('Acceptorinferreddefault',inferred_acceptor_scene('cDNA').strands['acceptor'].segs[0].inferred,True)
    path=marker_pathway()
    c('DTTretainsRNAthiolatedadduct',path[-1].retained_thiolated_adduct,True)
    c('DTTdoesnotrestorem6A',path[-1].name,'Tm6A')
    c('BiotinnotonreleasedTmRNA',path[-1].biotin_on_rna,False)
    try:ResidueMarker('m6A').undergo('cleave_disulfide');bad=False
    except ValueError:bad=True
    c('RejectreleasebeforeBAlabel',bad,True)
    c('Carrierafterchemistry',g.steps['carrier'].inputs[0].material,'ba_rcc_kept')
    c('SecondDTTusesactualresidualbeads',g.steps['dtt_release2'].inputs[0].material,'dtt_residual1')
    c('Poolactualtwoeluates',tuple(p.material for p in g.steps['eluate_pool'].inputs),('dtt_eluate1','dtt_eluate2'))
    c('DeadenylasebeforeRecJf',g.steps['recjf'].inputs[0].material,'deadenylase_product')
    c('RTafteractual10ultransfer',g.steps['rt'].inputs[0].material,'rt_transfer_used')
    c('PCRfrompostQC23ultransfer',g.steps['pcr'].inputs[0].material,'pcr_transfer_used')
    c('qPCRnotproduction',g.terminals['cycle_qc_assay'].disposition,Disposition.QC)
    c('QCmasteronephysicalstock',g.steps['qpcr_components'].inputs[0].material,'qpcr_master')
    c('NoNEBbodyborrowedbyname',g.materials['neb_universal'].state,None)
    c('Untransferredstocknotwaste',g.terminals['pcr_transfer_remainder'].disposition,Disposition.UNRESOLVED)
    c('Privateincompletegraph',bool(g.publication_problems()),True)
    c.report()

if __name__=='__main__':selftest()
