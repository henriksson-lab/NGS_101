"""Private main-article-qualified ExSeq preparation, not an imported FISSEQ recipe."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
AUTHORITY='PrintedpublishedScience371eaax2656/29Jan2021, authorhostedPDF3a64da_132e82c78d884323a15c487cb07fbf01.pdf; MethodsSummary2136–2148/2329–2342,Results224–240/565–572,actualFig1PDFp3. DetailedSM/tablesremainblocked; no wholeFISSEQ/defaultprimers or final-AMequivalence.'
class B:
    def __init__(self):
        self.g=ProtocolGraph('ExSeq selecteduntargeted expandedRNA→in-situ/ex-situreadout',authority=AUTHORITY,representation='published_instance');self.x=self.g.input('source','Sourceintactbiologicaltissue withnativeRNA; selecteduntargetedworkflow, sample-specificconditionsunassembled',K.BIOLOGICAL,primary=True)
    def stock(self,k,label,kind=K.MOLECULAR,role=R.STRUCTURAL,**kw):return(self.g.input(k,label,kind,**kw),role)
    def cond(self,k,label):return self.stock(k,label,K.CONDITION,R.CONDITION)
    def op(self,k,label,kind=S.MAGIC,ins=(),product='Actualsourcepreparation; wholeoccupiedmolecules/conditions/yieldunassembled',**kw):
        y=self.g.product(k+'_product',product,K.MOLECULAR)
        if kind==S.MAGIC:kw.setdefault('reason','Positivepublishedmainboundary; detailedSMrecipe/wholemolecules/chemicalproductsunassembled, notsourceabsence or borrowedFISSEQchemistry.')
        if kind==S.GENERIC:kw['changes_molecule']=False
        self.g.react(k,label,kind,((self.x,R.TEMPLATE),)+tuple(ins),(y,),**kw);self.x=y

def protocol_graph():
    b=B();g=b.g
    b.op('fix','Sourcefixedtissue; actualfixationrecipeunassembled',ins=(b.stock('fixation_supplies','Actualsourcefixationmolecularsupplies; identity/doseunassembled, notassumedPFAdefault'),))
    b.op('labelx','SourceRNAanchoringLabelX treatment',ins=(b.stock('labelx_reagent','ActualsourceLabelXRNAanchoringreagent; exactstructure/preparation/dose/occupiedRNAadductunassembled'),),product='SourceLabelX-treatedRNA/tissuepreparation; anchoringclassknown, wholeadductsunassembled')
    b.op('polymerize','SourceRNAanchoredintoswellablegel; chargedcarboxylic-acidgroupspresent',ins=(b.stock('swellable_gel','Actualswellablegelmolecularsupplies, chargedcarboxylicacidclasssource-positive; wholemonomers/crosslinker/formulationunassembled'),),product='RNA-bearingchargedhydrogel/tissuepreparation; exactpolymer/adductarchitectureunassembled')
    b.op('digest','Sourcemechanicalsoftening/digestionbeforeexpansion',ins=(b.stock('digestion_supplies','Actualsourcesoftening/digestionsupplies; enzymeidentity/compositionunassembled, no borrowedProteinaseK'),))
    b.op('expand','Sourceexpandswithwater; RNAremainsgelanchored',S.GENERIC,(b.cond('expansion_water','Actualwaterexpansionconditions; sampleexpansionfactor/volumeunassigned'),),product='ActualexpandedRNA-bearinggel, notnewtranscripts/sequenceconversion')
    b.op('reembed','SourcestabilizesexpandedmaterialbyreembeddinginUNCHARGEDgel',ins=(b.stock('uncharged_gel','Actualsourceunchargedreembeddingpolymerinventory; preciseformulationunassembled'),))
    b.op('passivate','Sourcechemicalpassivationenablesenzymaticreactionsinneutralchargeenvironment',ins=(b.stock('passivation_supplies','Actualchemicalpassivationcontributors; exactidentity/dose/atommapunassembled'),),product='ActualpassivatedreembeddedexpandedRNApreparation; notassumedcompletechargeconversion/yield')
    # Aggregate functional inventory is NOT evidence that these supplies were one premix or separately fresh doses.
    inv=g.input('library_inventory','Aggregateactualuntargetedin-situlibrarysupplies; physicalpackaging/stockcounts/internalallocationsunassigned',K.MOLECULAR)
    specs=[('rt_primer',K.OLIGO,'SourcerandomRTprimerclass; wholebody/lengthandchemicalmods unassembled',R.STRUCTURAL),('rt_enzyme',K.CATALYST,'ActualRNA-dependentDNApolymerasefunctionalclass; identity/doseunassembled',R.CATALYST),('rt_nt',K.MOLECULAR,'ActualRTnucleotidecontributors; composition/doseunassembled',R.STRUCTURAL),('circle_supplies',K.MOLECULAR,'ActualcDNAcircularizationfunctionalreactioninventory; identities/wholeendsunassembled',R.STRUCTURAL),('rca_supplies',K.MOLECULAR,'Actualrolling-circle-amplificationfunctionalreactioninventory; nucleotide/polymeraseallocationsandidentitiesunassembled',R.STRUCTURAL),('rca_priming',K.OLIGO,'ActualRCAprimingcomponentclass; wholebody/origin/relationtopriorRTprimerunassembled, notassumedfreshsecondprimerstock',R.STRUCTURAL)]
    out=[g.product(k,label,kind)for k,kind,label,role in specs];g.react('library_components','Functionalinventorybookkeeping; notphysicalseparation/premixorclaimoffreshduplicateadditions',S.GENERIC,((inv,R.TEMPLATE),),tuple(out),changes_molecule=False)
    g.react('reverse_transcribe','Sourcerandom-primedreverse-transcriptioninexpandedpassivatedsample',S.REVERSE_TRANSCRIPTION,((b.x,R.TEMPLATE),(out[0],R.STRUCTURAL),(out[1],R.CATALYST),(out[2],R.STRUCTURAL)),(g.product('cdna','Actualrandom-primedcDNApreparation; completebody/modifications/primersunassembled',K.MOLECULAR),));b.x=g.materials['cdna']
    b.op('size_restrict','PublishedResultsrestrictcDNAfragmentstoapproximately100basestoimprovecircularization',ins=(b.stock('size_supplies','ActualcDNAsize-restrictionprocessinginventory; cleavage/selectionmechanism/orderwithotherSMstepsunassembled'),),product='Sourceapproximately100basecDNAfragmentpreparation; observedmean76isnotimposedselectionlength')
    b.op('circularize','SourcecDNAcircularization; wholemolecularends/ligationgeometryandproductsunassembled',ins=((out[3],R.STRUCTURAL),),product='ActualcircularizedcDNApreparationincludingunreactedclasses; no invented5P/body/circlejunction')
    rcapol=g.product('rca_polymerase','DNApolymerasefunctionalcomponentofSAMEactualRCAinventory; identity/doseunassigned',K.CATALYST);rcant=g.product('rca_nucleotides','NucleotidefunctionalcomponentofSAMEactualRCAinventory; composition/doseunassigned',K.MOLECULAR)
    g.react('rca_components','Copresentfunctionalcomponents ofaggregateRCAinventory, notnewstocks/premixorphysicalseparation',S.GENERIC,((out[4],R.TEMPLATE),),(rcapol,rcant),changes_molecule=False)
    b.op('rca','SourceinsituRCAofcDNAcircles',ins=((rcapol,R.CATALYST),(rcant,R.STRUCTURAL),(out[5],R.STRUCTURAL)),product='Actualin-situampliconconcatemerpreparation; wholecircle/primer/repeatsequenceandchemicalmods unassembled')
    b.op('in_situ_sequencing','Sourceiterativefluorescentsequencingchemistryandopticalimagingrounds',ins=(b.stock('in_situ_oligos','Actualfluorescentsequencingoligoinventory; wholeorderings/mods/roundscheduleunassembled',K.OLIGO),b.stock('in_situ_other','Otheractualin-situsequencingchemicalsupplies; identity/packagingunassembled')),product='SAMEexpandedfixedsampleafterin-situopticalsequencing; imagingdoesnotremoveaphysicalDNAaliquot')
    # Optical sequence reads/data matching do not physically alter DNA/RNA molecules.
    released=g.product('ex_situ_amplicons','ActualcDNAampliconselutedfromexpandedsequencedsampleforex-situresequencing',K.MOLECULAR);residual=g.product('residual_sample','Residualgel/tissue/ampliconpreparationafterelution; completeness/fateunassigned',K.MOLECULAR)
    g.react('extract_amplicons','Sourcephysicalpostin-situcDNAampliconelution; no computationalmatchingstepmisclassifiedasextraction',S.MAGIC,((b.x,R.TEMPLATE),b.stock('elution_supplies','Actualampliconelutioncontributors; exactreaction/recoveryunassembled')),(released,residual),reason='PublishedFig1C/MethodsSummaryestablishactualampliconrecovery; detailedSMmechanismandmolecularendsunassembled.');g.mark_terminal(Terminal(residual.id,Disposition.UNRESOLVED,'Residualphysicalsample; no inventedcompletedigestion/discard.'));b.x=released
    b.op('ex_situ_library','Sourceex-situIlluminasequencinglibraryfromrecoveredamplicons',ins=(b.stock('ex_situ_oligos','ActualNGSlibraryoligocontributors; wholebodies/indexstatusandprimerlandingunassembled',K.OLIGO),b.stock('ex_situ_supplies','Otheractualex-situlibrarymolecularsupplies; completekit/enzymes/NTs/physicalstocksunassembled')),product='Actualex-situampliconsequencinglibrary; completeconstruction/fragmentation/indexingchemistryunassembled')
    selected=g.product('selected_loading','Actualselectedex-situsequencingloadingallocation; amountunassigned',K.MOLECULAR);rem=g.product('loading_remainder','Untransferredex-situlibrarystock; fateunassigned',K.MOLECULAR);g.react('loading_allocate','Actualsequencingallocationseparatefromwholestock',S.SPLIT,((b.x,R.TEMPLATE),),(selected,rem));g.mark_terminal(Terminal(rem.id,Disposition.UNRESOLVED,'Untransferredsequencingstockretained.'))
    lib=g.product('library','Sourceex-situIlluminaresequencinglibrary; wholemolecularconstructunassembled',K.LIBRARY);g.react('load','Sourceex-situNGSreadout; exactinstrument/runlengthunassembled',S.GENERIC,((selected,R.TEMPLATE),),(lib,),changes_molecule=False);g.sequence(lib,SequencingPlan(limitation='PublishedFig1in-situ17roundexample,variablein-situ5–30basestypically; ex-situIlluminareadscontainedcDNArepeats. No syntheticUMI/barcode inferred fromendogenoussequencefingerprints or exactrunlength/kitbody borrowed.'))
    g.add_source_item(SourceItem('source_scope',AUTHORITY,'Mainarticle-supporteduntargetedpreparationonly. ActualLabelX/unchargedreembed/passivationandampliconrecoveryboundariessupportidentitybeyondphysicalresolution/software. DetailedSM/tablecontentsunreviewed; targetedpadlockvariantnotmerged.'))
    g.add_source_item(SourceItem('computational_matching','SourceFig1C/MethodsSummaryusesrecoveredampliconsequenceasdictionaryanduniquematchesin-situreads.','Matching/segmentation/barcode-likeendogenoussequenceidentificationarecomputational, notnewwetligation or syntheticUMIinstallation.'))
    return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph()
    c('LabelXbeforegel',g.steps['polymerize'].inputs[0].material,'labelx_product')
    c('Reembedafterexpansion',g.steps['reembed'].inputs[0].material,'expand_product')
    c('RTafterunchargedreembedandpassivation',g.steps['reverse_transcribe'].inputs[0].material,'passivate_product')
    c('RandomprimerbodyunassemblednotborrowedFISSEQ',g.materials['rt_primer'].state,None)
    c('Sizeprocessingbeforecircle',g.steps['circularize'].inputs[0].material,'size_restrict_product')
    c('RCAcatalystfromsameaggregateinventory',g.producer_of('rca_polymerase').id,'rca_components')
    c('RCAnucleotidescontributewithoutfreshsecondstock',any(i.material=='rca_nucleotides' and i.role==R.STRUCTURAL for i in g.steps['rca'].inputs),True)
    c('Ampliconrecoveryafteropticalseq',g.steps['extract_amplicons'].inputs[0].material,'in_situ_sequencing_product')
    c('Residualgelnotautomaticallydiscarded',g.terminals['residual_sample'].disposition,Disposition.UNRESOLVED)
    c('Exsitupreparationusesactualrecoveredmolecules',g.steps['ex_situ_library'].inputs[0].material,'ex_situ_amplicons')
    c('Finalloadingallocationretainsstock',g.terminals['loading_remainder'].disposition,Disposition.UNRESOLVED)
    c('NoexperimentalNGSkitinvented',g.materials['library'].state,None)
    c('OpticalreadoutnotphysicalQCaliquot',any(t.disposition==Disposition.QC for t in g.terminals.values()),False)
    p=sources_root()/Path(__file__).resolve().parents[1].name/'3a64da_132e82c78d884323a15c487cb07fbf01.pdf.txt'
    if have(c,Source(p,'https://doi.org/10.1126/science.aax2656')):
        t=''.join(p.read_text().split());c('Printedpublishedidentity','Science371,eaax2656(2021)'in t,True);c('Actualunchargedreembedding','re-embeddingtheminunchargedgels'in t,True);c('ActualLabelXsource','RNAanchoringreagentLabelX'in t,True)
    c('Endpointretainsprimarybiologicalprovenance',any('lacks primary-template provenance' in x for x in g.publication_problems()),False)
    c('Privateincompletegraph',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
