"""Private source-qualified partial preparation; no complete molecular/absence claim."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,feature,revcomp
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
class B:
    def __init__(self,title,authority,label,kind=K.MOLECULAR):
        self.g=ProtocolGraph(title,authority=authority,representation='published_instance');self.x=self.g.input('source',label,kind,primary=True)
    def stock(self,key,label,kind=K.MOLECULAR,role=R.STRUCTURAL,**kw):return(self.g.input(key,label,kind,**kw),role)
    def cat(self,key,label):return self.stock(key,label,K.CATALYST,R.CATALYST)
    def cond(self,key,label):return self.stock(key,label,K.CONDITION,R.CONDITION)
    def op(self,key,label,kind=S.MAGIC,inputs=(),product='Actual source preparation; whole occupied substrates/products/side products and conversion unassembled',**kw):
        y=self.g.product(key+'_product',product,K.MOLECULAR)
        if kind==S.MAGIC:kw.setdefault('reason','Positive source operation; whole occupied molecules, end distributions, detailed molecular products and conversion unassembled, not a source-absence claim.')
        if kind==S.GENERIC:kw['changes_molecule']=False
        self.g.react(key,label,kind,((self.x,R.TEMPLATE),)+tuple(inputs),(y,),**kw);self.x=y;return y
    def clean(self,key,label,inputs=()):
        y=self.g.product(key+'_kept','Actual retained physical source fraction, carryover/recovery unassigned',K.MOLECULAR);w=self.g.product(key+'_waste','Actual excluded workup fraction',K.WASTE)
        self.g.react(key,label,S.CLEANUP,((self.x,R.TEMPLATE),self.cond(key+'_condition',label))+tuple(inputs),(y,w));self.g.discard(w,'Actual source workup fraction excluded; no quantitative purification claim.');self.x=y;return y
    def aliquot(self,key,label,kept='Selected source aliquot',rest='Untransferred stock, amount/fate unassigned',qc=False):
        y=self.g.product(key+'_selected',kept,K.MOLECULAR);r=self.g.product(key+'_remainder',rest,K.MOLECULAR)
        self.g.react(key,label,S.SPLIT,((self.x,R.TEMPLATE),),(y,r))
        if qc:self.g.mark_terminal(Terminal(y.id,Disposition.QC,kept));self.x=r;return r
        self.g.mark_terminal(Terminal(r.id,Disposition.UNRESOLVED,rest));self.x=y;return y




from adapter_ligation import AdapterLigationGeometry
from illumina import TRUSEQ_READ1,TRUSEQ_READ2
import re
AUTHORITY='PublishedPMC7026646 main76–126/150/170–176 andsameoriginalFig1; actualTableS1s-oligos/qPCR andSI1XF2blueprint. SelectedQKIHis→Strep/PK/SSIII/RNaseH XF2 only; noXF1/3/4/FLAG/transientbranchmerge.'
LEFT=revcomp(TRUSEQ_READ1[-18:]);RIGHT=TRUSEQ_READ2[-19:]
# Exact whole modified ordering fixtures from original TableS1, not inferred hairpin sequences.
SOURCE_ORDERINGS={'r2A-P': '/5Phos/rNrNrCrArCrUrUrGrNrYrYrNrNAGATCGGAAGAGCGTCGT/iSp18/ACGTGTGCTCTTCCGATCT/3Phos/', 'r2B-P': '/5Phos/rNrNrCrArCrUrUrGrNrRrRrNrNAGATCGGAAGAGCGTCGT/iSp18/ACGTGTGCTCTTCCGATCT/3Phos/', 'r4A-P': '/5Phos/rNrNrGrGrCrUrArArNrYrYrNrNAGATCGGAAGAGCGTCGT/iSp18/ACGTGTGCTCTTCCGATCT/3Phos/', '4B-P': '/5Phos/rNrNrGrGrCrUrArArNrRrRrNrNAGATCGGAAGAGCGTCGT/iSp18/ACGTGTGCTCTTCCGATCT/3Phos/', '6A-P': '/5Phos/rNrNrUrUrCrCrArCrNrYrYrNrNAGATCGGAAGAGCGTCGT/iSp18/ACGTGTGCTCTTCCGATCT/3Phos/', '6B-P': '/5Phos/rNrNrUrUrCrCrArCrNrRrRrNrNAGATCGGAAGAGCGTCGT/iSp18/ACGTGTGCTCTTCCGATCT/3Phos/', '11A-P': '/5Phos/rNrNrCrUrCrArCrArNrYrYrNrNAGATCGGAAGAGCGTCGT/iSp18/ACGTGTGCTCTTCCGATCT/3Phos/', '11B-P': '/5Phos/rNrNrCrUrCrArCrArNrRrRrNrNAGATCGGAAGAGCGTCGT/iSp18/ACGTGTGCTCTTCCGATCT/3Phos/', '13A-P': '/5Phos/rNrNrGrUrUrCrCrArNrYrYrNrNAGATCGGAAGAGCGTCGT/iSp18/ACGTGTGCTCTTCCGATCT/3Phos/', '13B-P': '/5Phos/rNrNrGrUrUrCrCrArNrRrRrNrNAGATCGGAAGAGCGTCGT/iSp18/ACGTGTGCTCTTCCGATCT/3Phos/', '15A-P': '/5Phos/rNrNrGrUrArArCrCrNrYrYrNrNAGATCGGAAGAGCGTCGT/iSp18/ACGTGTGCTCTTCCGATCT/3Phos/', '15B-P': '/5Phos/rNrNrGrUrArArCrCrNrRrRrNrNAGATCGGAAGAGCGTCGT/iSp18/ACGTGTGCTCTTCCGATCT/3Phos/', '17A-P': '/5Phos/rNrNrGrUrGrArUrGrNrYrYrNrNAGATCGGAAGAGCGTCGT/iSp18/ACGTGTGCTCTTCCGATCT/3Phos/', '17B-P': '/5Phos/rNrNrGrUrGrArUrGrNrRrRrNrNAGATCGGAAGAGCGTCGT/iSp18/ACGTGTGCTCTTCCGATCT/3Phos/', '410a': 'TACACGACGCTCTTCCGATCT', '410b': 'GACGTGTGCTCTTCCGATCT'}

def polymer(name):
    raw=re.sub(r'/[^/]+/','',SOURCE_ORDERINGS[name]);return ''.join(re.findall(r'r?([ACGUTNRY])',raw))
def adapter_parts(name='r2A-P'):
    seq=polymer(name);tag=seq[:13]
    return [Segment('random_first',tag[:2],feature=feature('flash.molecule','umi','random',group='flash.molecule',note='SelectedXF1–8 main150 fiveNdeduplicationpositions, twohere')),Segment('pulldown_condition',tag[2:8],feature=feature('flash.pulldown','sample_index','whitelist',whitelist='PublishedTableS1sevenRNAfixedsixconditioncodes')),Segment('random_middle',tag[8:9],feature=feature('flash.molecule','umi','random',group='flash.molecule',note='Oneoffiverandomdeduppositions')),Segment('replicate_code',tag[9:11],placeholder=True,feature=feature('flash.replicate','sample_index','degenerate',note='SourceRY-spacebiologicalreplicatetags; notselectedpipelineUMI')),Segment('random_last',tag[11:13],feature=feature('flash.molecule','umi','random',group='flash.molecule',note='Twooffiverandomdeduppositions')),Segment('left_stem',LEFT[:13]),Segment('left_loop_flank',LEFT[13:]),Segment('right_loop_flank',RIGHT[:6],note='Non-nucleotideiSp18connectspreviousDNAflankandthisflank; itisnot18DNA/RNAbases'),Segment('right_stem',RIGHT[6:])]
def oligo_scene(name='r2A-P'):
    s=Scene()
    if name in ('410a','410b'):
        seq=SOURCE_ORDERINGS[name];canonical=('T'+TRUSEQ_READ1[-20:])if name=='410a'else('G'+RIGHT)
        if canonical!=seq:raise ValueError('CanonicalqPCRbodymustpreservewholesourceordering')
        s.strand('oligo',[Segment('sourceqPCRprimer',canonical)]);return s
    parts=adapter_parts(name)
    if LEFT[:13]!=revcomp(RIGHT[6:]):raise ValueError('Sources-oligo13bpstemdoesnotpair')
    if ''.join(p.top for p in parts)!=polymer(name):raise ValueError('Sourcepolymerorderingchanged')
    s.strand('oligo',parts[:7],mod5='P',mod3='iSp18 continuation');s.strand('oligo_after_spacer',parts[7:],mod5='iSp18 continuation',mod3='P');s.note('oligo','ONE covalent source oligo: RNA13 plus leftDNA18 → nonnucleotide iSp18 → rightDNA19. Separate linker-connected nucleotide tracts, not free oligos or direct nucleotide adjacency across spacer.');s.ligation_end('oligo','random_first',kind='single-strand',polarity='5-prime',owner='adapter',text='KnownsourceRNA5P/s-oligolocalend; fulloccupiedRNA/hairpinjunctionunassembled')
    return s

def protocol_graph():
    b=B('FLASH selectedXF2 dual-affinity/RNAhairpin/selfprimedcDNAcircle preparation',AUTHORITY,'SourcepreparedinducedQKI3xFlag-HBHcells; constructengineering/upstreamcultureunassembled',K.BIOLOGICAL);g=b.g
    g.add_source_item(SourceItem('tag_meaning','Main172 genericNNXXXXXXNNNNN describes7randomribonucleotides. Main150 prints processing pattern NNB1B2NT1T2T3T4T5T6NN with fiveNdedup/twoRYbioreplicate/sixcondition roles; TableS1 raw RNA order is NN+[fixed6]+N+[YY/RR]+NN.','Both literal source orderings retained, not silently normalized into one coordinate map. TableS1literalYY/RR retainedasbioreplicatefeature, notuniformrandom7ntUMI. Roleassignment source-qualified; actualprocessingcoordinates/readorientation/wholeprimerlandingunassembled; no inferredreadmask.'))
    g.add_source_item(SourceItem('spacer','Main86calls8-atomhexaethyleneglycol; main172/176C18; TableS1orders/iSp18/.','PreserveliteralnonnucleotideiSp18and13bppairedDNAstemwith5/6ntDNAflanks; not18basesandnotperfect18/19bpstem. Preciseatomnomenclatureunassembled.'))
    g.add_source_item(SourceItem('source_programmes','Main96literal0.15mJ/cm2UVC; main126PCRlistedcomponents20+1+1+9+8ul.','LiteralUVvalueandlisted39ultotalretainedwithoutrepair; actualUVcalibration/volumeconventionunassigned. 410a/bareqPCRprequantificationprimers, notproductionP5/P3wholebodies.'))
    b.op('uv','SourcePBSwash/0.15mJcm2UVCicewatercrosslink',inputs=(b.cond('uv_conditions','ActualsourceUVprogramme'),),product='SourceUVcrosslinkedcell/RNPclasses; crosslinksites/yieldunassembled')
    b.clean('cell_pellet','Actual500g/PBSwash/re-pellet/snapfreeze/thaw')
    b.op('lysis','Source550ulNLB/sonicate5min30sONOFFLOW5cycles4C',inputs=(b.cond('nlb','ActualNLBsourceformulation'),));b.clean('clarify','Source20000g10min4Cinsolublematerialexcluded')
    def capture(key,label):
        bound=g.product(key+'_bound','ActualselectedaffinitybeadassociatedRNPpreparation',K.MOLECULAR);sup=g.product(key+'_supernatant','Actualremovedunboundfraction',K.WASTE)
        g.react(key,label,S.AFFINITY_CAPTURE,((b.x,R.TEMPLATE),b.stock(key+'_beads','Actual25ulsourceaffinitybeads washed/resuspended500ulNLB',K.SELECTION,R.SELECTION)),(bound,sup));g.discard(sup,'Actualsourceunboundfractionremoved; noquantitativecaptureclaim');b.x=bound
    capture('his','Source25ulHisbeads5minmax10min');b.clean('his_wash','Source800ulNLBwash')
    eluate=g.product('his_eluate','Actual250mMimidazoleeluate',K.MOLECULAR);res=g.product('his_residual','ActualresidualHisbeads/retainedmaterial',K.MOLECULAR)
    g.react('imidazole_release','Source250mMimidazoleNLB10minice',S.MAGIC,((b.x,R.TEMPLATE),b.cond('imidazole','Actual250mMimidazoleNLB')),(eluate,res),reason='Sourceimidazoleaffinityrelease; fullmolecularclassesandresidualyieldunassembled.');g.mark_terminal(Terminal(res.id,Disposition.UNRESOLVED,'SourceHisresidualfateunassigned'));b.x=eluate
    capture('strep','Source25ulC1streptavidin4C1h');b.clean('strep_wash','SourceLDS/PLB/HSB/NDBserialwash')
    b.op('trim','Source1mlNDB/2ulTurboDNase/10ulRNaseI1to2000/37C3minice2min',S.DIGESTION,(b.cat('turbo_dnase','Actual2ulTurboDNase'),b.cat('rnase_i','Actual10uldilutedRNaseI1to2000'),b.cond('trim_conditions','Actualsourcebuffer/time')));b.clean('trim_wash','Actualsupernatantremoval/HSB/NDBwash')
    def pnk(key,label):
        b.op(key,label,S.DEPHOSPHORYLATION,(b.cat(key+'_enzyme','Actualfresh1ulT4PNK10Uul'),b.cond(key+'_buffer','Actual20ulPNKMES6/RNasin/betaME/noATP/37C20min')));b.clean(key+'_wash','SourceHSB/NDBwash')
    pnk('rna_heal','SourceRNA3cyclicphosphatehealingbeforeadapterligation')
    oligo=g.input('s_oligo','SourceTableS1r2A-Pselectedlocalreference, actualstudyconditionassignmentunassigned',K.OLIGO,state=oligo_scene(),chemical_state={'source_ordering':SOURCE_ORDERINGS['r2A-P'],'sugar_blocks':[{'start':0,'end':13,'sugar':'RNA'},{'start':13,'end':50,'sugar':'DNA'}],'five_prime':'P','three_prime':'P','non_nucleotide_linker':{'after_polymer_index':30,'ordering':'/iSp18/','atom_count_nomenclature':'unassembled'},'fold_candidate':{'inferred':True,'paired_stem':13,'loop_flanks':[5,6],'non_nucleotide_spacer':'iSp18'}})
    b.op('adapter_ligate','Source20ulRnl1/1ul10uMs-oligo/2ul1mMATP/PEG/25C1h',S.LIGATION,((oligo,R.STRUCTURAL),b.cat('rnl1','Actual1ulT4RNAligase1'),b.stock('lig_atp','Actual2ul1mMATP'),b.cond('lig_conditions','Actual2ul10xbuffer/4ulPEG/.5ulRNasin/9.5ulwater')),geometry=AdapterLigationGeometry.known('single-strand',oligo_scene(),statement='KnownsourceRNA5P localadapterend; RNA3OHconditionalafterPNK, fulloccupiedhairpin/sealedRNAjunctionunassembled.'),product='Sourcebeadassociateds-oligoligatedRNP/freeadapter/unreactedclasses; fullproductsunassembled');b.clean('lig_wash','SourceHSB/NDBwash')
    b.aliquot('tagged_pool_input','Actualselecteduniquelytaggedsamplepoolcontribution; amountunassigned')
    others=g.input('other_tagged_samples','Otheractualuniquelytaggedsamplebeadcontributions; count/allocations/ordersunassigned',K.MOLECULAR);pool=g.product('tagged_pool','Actualuniquelytaggedbeadpool',K.MOLECULAR);g.react('tagged_pooling','SourcepoolrelevantsamplesAFTERuniqueadaptertagging',S.POOL,((b.x,R.TEMPLATE),(others,R.TEMPLATE)),(pool,));b.x=pool
    pnk('adapter_deprotect','SourcefreshPNKremovess-oligo3Paftertaggedpooling')
    b.op('pk','SelectedXF2source100ulPKmix/10ul20mgmlPK37Cdurationunassigned',S.DEPROTEINIZATION,(b.cat('proteinase_k','Actual10ul20mgmlPKin100ulmix'),b.cond('pk_conditions','SourceTris/NaCl/Tween/EDTA/SDS37C')),product='SourceRNAreleasedfromprotein; wholecovalentRNA/s-oligoandresidualclassesunassembled')
    b.clean('rna_column','SourceZymoOligo200ulbinding/400ulEtOH/elute9.5ul')
    b.aliquot('rt_input','Actual9of9.5ulRNAeluate','Actual9ulRNAforRT','Actual.5uluntransferredRNAeluate')
    rt=g.product('rt_product','Sourceintramoleculars-oligoselfprimedcDNA/RNAhybridclasses',K.MOLECULAR,chemical_state={'primer_provenance':'Covalentlyretaineds-oligoDNA3OHaftersecondPNK; nofreshRTprimer','copied_tags':'Source13ntRNAadaptersegmentcopiedintoDNA; fullwholecDNAsequenceunassembled','self_priming_geometry':'Sourcepositivehairpincontext; completeoccupiedstructureunassembled'})
    g.react('rt','Source20ulSSIII42/50/55C10mineach',S.MAGIC,((b.x,R.TEMPLATE),b.cat('ssiii','Actual1ulSuperScriptIII'),b.stock('rt_dntp','Actual1ul10mMdNTP'),b.cond('rt_conditions','Actual2ul10xbuffer/4ul25mMMg/2ul.1MDTT/1ulRNaseOUT')),(rt,),reason='Actualintramolecularhairpinpriming; fullRNA/cDNAhybrid/branchendgeometriesunassembled. Noexternaloligoprimeraddition.');b.x=rt
    b.op('rnase_h','SelectedXF2direct1ulRNaseH37C20minthen54ulwater',S.DIGESTION,(b.cat('rnase_h_enzyme','Actual1ulRNaseH'),b.cond('rnase_h_conditions','Actual37C20min/add54ulwater→75ul')),product='Sourceintended5PcDNAafterRNAhybridandoriginalRNAadaptersegmentdigestion; incompleteconversion/fullproductsunassembled')
    b.clean('cdna_column','SourceZymo150ulbinding/300ulEtOH/elute6.5ul')
    b.aliquot('circle_input','Actual6of6.5ulcDNAeluate','Actual6ulforcircle','Actual.5uluntransferredcDNAeluate')
    b.op('circle','Source10ulCircLigaseII60C~16h',S.MAGIC,(b.cat('circligase','Actual.5ulCircLigaseII'),b.cond('circle_conditions','Actual1ulbuffer/.5ul50mMMn/2ul5Mbetaine')),product='SourcecircularizedcDNA/unreacted/offproductclasses; occupiedendjunctionandconversionunassembled')
    circ=b.x;q=g.product('qpcr_aliquot','Actual1ulcircleforcycle-estimatingqPCR',K.MOLECULAR);prod=g.product('production_aliquot','Actual9ulcircleforproductionPCR',K.MOLECULAR);g.react('circle_allocation','SelectedsourceqPCRguidedroute1ulqPCR/9ulproduction',S.SPLIT,((circ,R.TEMPLATE),),(q,prod))
    qpf=g.input('qpcr_410a','Source410aqPCRonly',K.OLIGO,state=oligo_scene('410a'),chemical_state={'source_ordering':SOURCE_ORDERINGS['410a']});qpr=g.input('qpcr_410b','Source410bqPCRonly',K.OLIGO,state=oligo_scene('410b'),chemical_state={'source_ordering':SOURCE_ORDERINGS['410b']});qdone=g.product('qpcr_readout','SourceqPCRusedforestimatingPCRcyclenumber',K.MOLECULAR)
    g.react('qpcr','Source410aqPCRprequantification; Ctusedasproductioncyclenumber',S.MAGIC,((q,R.TEMPLATE),(qpf,R.STRUCTURAL),(qpr,R.STRUCTURAL),b.stock('qpcr_supplies','ActualqPCRsupplies, reactioncompositionunassembled')),(qdone,),reason='ActualqPCRcontrol; detailedprogramme/formulation/Ctreadoutunassigned. Notproductionprimerbodyevidence.');g.mark_terminal(Terminal(qdone.id,Disposition.QC,'Actualcycle-estimatingqPCRaliquot'))
    b.x=prod
    master=g.input('q5_master','Actual20ul2xNEBNextQ5MasterMix',K.MOLECULAR);enz=g.product('q5_enzyme','SameactualmasterQ5enzymecomponent',K.CATALYST);nt=g.product('q5_nt','SameactualmasterdNTPcomponent',K.MOLECULAR);buffer=g.product('q5_buffer','Sameactualmasterbuffercomponent',K.CONDITION);g.react('q5_components','SamephysicalPCRmastercomponents, nonewseparateadditions',S.GENERIC,((master,R.TEMPLATE),),(enz,nt,buffer),changes_molecule=False)
    p5=g.input('production_p5','Actual1ul10uMuniversalP5productionprimer; wholebodyunassembled',K.OLIGO);p3=g.input('production_p3','Actual1ul10uMbarcodedP3productionprimer; wholebody/indexunassembled',K.OLIGO,chemical_state={'sample_index':feature('flash.production_sample','sample_index','unknown',note='SourcebarcodedP3; sequence/coordinateunassembled')})
    b.op('production_pcr','Source20ulmaster/1ulP5/1ulP3/9ulcircle/8ulwater; listed39ultotal',inputs=((p5,R.STRUCTURAL),(p3,R.STRUCTURAL),(enz,R.CATALYST),(nt,R.STRUCTURAL),(buffer,R.CONDITION),b.cond('pcr_conditions','SourceCtguidedcyclenumber, numericvalueunassigned')),product='SourceproductionPCRlibraryclasses; fullduplex/readlandingunassembled')
    b.clean('ampure_first','Actual1.5xAMPurecleanup');b.clean('ampure_second','Actualsecond1.5xAMPurecleanup/elute10mMTris8/.1percentTween')
    b.aliquot('library_qc','ActualQubitDNAHS/BioanalyzerHSQCaliquot; amountunassigned',qc=True)
    b.aliquot('loading','Actuallibrarysubmittedforsequencing; loadingallocationunassigned')
    lib=g.product('library','SourceFLASHXF2sequencinglibraryclass; wholeduplex/readgeometryunassembled',K.LIBRARY);g.react('load','Sourcepairedendhigh-throughputsequencing',S.GENERIC,((b.x,R.TEMPLATE),),(lib,),changes_molecule=False);g.sequence(lib,SequencingPlan(limitation='SourceXF1–8pairedendprocessing main150; exactinstrument/cycles/productionprimerlandingunassembled. FiveNdedup/twoRYreplicatetags/fixedsixcondition distinct; originalRNAadaptertags copiedbeforeRNaseH.'))
    return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph()
    c('Exact13bpDNAstem',LEFT[:13],revcomp(RIGHT[6:]))
    c('Sourceadapterpolymer50not68nt',len(polymer('r2A-P')),50)
    c('Adapter5Pphysicalmark',oligo_scene().strands['oligo'].mod5,'P')
    c('Spacer separates nucleotide tracts of one source oligo',tuple(''.join(x.top for x in strand.segs)for strand in oligo_scene().strands.values()),(polymer('r2A-P')[:31],polymer('r2A-P')[31:]))
    c('Healingbeforeligation',list(g.steps).index('rna_heal')<list(g.steps).index('adapter_ligate'),True)
    c('TaggedpoolbeforefreshPNK',g.steps['adapter_deprotect'].inputs[0].material,'tagged_pool')
    c('RTnineulcomesfromeluate',g.steps['rt'].inputs[0].material,'rt_input_selected')
    c('NoaddedRTprimer',all(x.material not in('production_p5','production_p3','s_oligo')for x in g.steps['rt'].inputs),True)
    c('RNaseHnotNaOHselectedXF2',g.steps['rnase_h'].kind,S.DIGESTION)
    c('Circle6ulstockremainder',g.materials['circle_input_remainder'].id,'circle_input_remainder')
    c('qPCRandproductiondistinctactualaliquots',g.steps['production_pcr'].inputs[0].material,'production_aliquot')
    c('ProductionprimersnotqPCR410',tuple(x.material for x in g.steps['production_pcr'].inputs[1:3]),('production_p5','production_p3'))
    c('CopiedRNAdegradationqualified','originalRNAadaptersegment'in g.materials['rnase_h_product'].label,True)
    p=sources_root()/Path(__file__).resolve().parents[1].name/'nar_gkz1141_supp_gkz1141_supplemental_files_Supplementary_Table_1.xlsx.txt'
    if have(c,Source(p,'https://doi.org/10.1093/nar/gkz1141')):
        rows={}
        for line in p.read_text().split('## sheet s-oligos')[1].splitlines():
            cols=line.split('\t')
            if len(cols)>1:rows[cols[0]]=cols[1].strip()
        for name,raw in SOURCE_ORDERINGS.items():c('Wholeoriginalmodifiedordering '+name,rows.get(name),raw);oligo_scene(name)
    c('Privatepartial',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
