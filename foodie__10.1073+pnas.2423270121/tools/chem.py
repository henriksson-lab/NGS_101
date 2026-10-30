"""FOODIE bulk library preparation and separate source single-cell adapter component.

The local single-cell duplexes are not the identities of the bulk Vazyme stock.
Protein production and whole occupied chromatin geometry remain outside the graph.
"""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,revcomp,feature
from nextera import ME
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from pcr import PCRGeometry
P5_TAIL='GGAGTGGTAGGTGTGGG'
P7_TAIL='GGGAGTGAGGTGAGGTG'
TN5_BC1='TATATTGTGAG'
AUTHORITY='PublishedPMC11670102main184–188/282–300; publishedSI TableS2 actualPDFp6(printed38). Two originals, preprint excluded.'

def single_cell_adapter(which):
    """Source TableS2 local duplex, not a full integration state or bulk recipe."""
    if which not in ('p5','p7_1'):raise ValueError('Only the two source-selected TableS2 rows are assembled here')
    tail=P5_TAIL if which=='p5'else P7_TAIL
    segs=[Segment('cytosine-depleted handle',tail)]
    if which=='p7_1':segs.append(Segment('source Tn5 barcode1',TN5_BC1,feature=feature('foodie.tn5.bc1','cell_barcode','fixed',group='foodie.cell',part='tn5_round')))
    segs.append(Segment('ME',ME))
    s=Scene();s.strand('transferred',segs,mod3='OH')
    s.anneal('nontransferred',[Segment('sourceMEcomplement',revcomp(ME))],to='transferred',pair=('sourceMEcomplement','ME'),mod5='P')
    return s

def single_cell_orderings():
    return {'p5':P5_TAIL+ME,'p7_1':P7_TAIL+TN5_BC1+ME,'nontransferred':'/5Phos/'+revcomp(ME)}

def unknown_scene():
    s=Scene();s.strand('unknown',[Segment('whole products','[unassembled source products]',placeholder=True)]);return s

def protocol_graph():
    g=ProtocolGraph('FOODIE selected fresh-cell bulk library',authority=AUTHORITY,representation='published_instance')
    x=g.input('cells','Source30000freshK562orGM12878cells; actualchosenlineunassigned, noheatshock/fixationbranchmerged',K.BIOLOGICAL,primary=True)
    def inp(k,label,kind=K.CONDITION,role=R.CONDITION,**kw):return(g.input(k,label,kind,**kw),role)
    def op(k,label,kind,inputs=(),product='Source mixture with incomplete/unreacted/off-target and carryover classes; quantities/occupancy unassigned',**kw):
        nonlocal x
        y=g.product(k+'_product',product,K.MOLECULAR)
        if kind==S.GENERIC:kw['changes_molecule']=False
        g.react(k,label,kind,((x,R.TEMPLATE),)+tuple(inputs),(y,),**kw);x=y
    def clean(k,label):
        nonlocal x
        y=g.product(k+'_kept','Actualsourcephysicallyretainedfraction; imperfectrecovery/carryover unassigned',K.MOLECULAR);w=g.product(k+'_removed','Actualsourceaspirated/removedworkupfractions',K.WASTE)
        g.react(k,label,S.CLEANUP,((x,R.TEMPLATE),inp(k+'_conditions',label)),(y,w));g.discard(w,'Physicalworkup, notperfectmolecularspecificity');x=y
    clean('cell_wash','Source450g5minpellet/equalvolumecoldPBS/500g5min4C')
    op('permeabilize','Sourcecoldpermeabilization10mMTris7.4/10mMNaCl/3mMMg/.1percentIGEPAL/.1percentTween/.1percentdigitonin',S.MAGIC,(inp('perm_conditions','Sourcepermeabilizationbuffer'),),reason='Sourcepermeabilizationpositive; fullcellular/nuclearcomplexesunassembled.')
    clean('nuclei','Source550g5min4C; carefullyaspiratesupernatant, retainnucleipellet')
    bulk=g.input('transposome_stock','Sourceactual1ulVazymeS111transposomestock; exactloadedoligoidentitiesunassembled',K.MOLECULAR)
    tn=g.product('tn5_enzyme','CatalystcomponentoftheactualS111stock',K.CATALYST);ad=g.product('tn5_adapters','DNAcomponentofsameS111stock, notTableS2singlecellcustomoligos',K.ADAPTER);tb=g.product('tn5_stock_buffer','ConditioncomponentofsameS111stock',K.CONDITION)
    g.react('transposome_components','Representcomponentsofsamesourcephysicalstock, notfreshindependentadditions',S.GENERIC,((bulk,R.TEMPLATE),),(tn,ad,tb),changes_molecule=False)
    op('tagment','Source25ul2xTD/1ulS111/24ulPBS+.1percentdigitonin;37C30min800rpm',S.MAGIC,((tn,R.CATALYST),(ad,R.STRUCTURAL),(tb,R.CONDITION),inp('tag_conditions','Source2xTD20mMTAPS8.5/10mMMg/20percentDMF beforemixing')),
       reason='SourceTn5positive; fullfragmentends/loadedoligos/occupiedchromatinunassembled.',product='SourceTn5tagged/unreactedchromatin and othermixtureclasses; fulladapters/productdistributionunassembled')
    clean('tag_stop','Source100ulicecoldRSB/550g5min4C; carefullyaspiratesupernatant')
    op('deaminate','SourcepreparedDddB inDRB37C10min; selectedfreshbulkroute',S.CHEMICAL_CONVERSION,(inp('dddb','SourcepreparedDddB catalyst; proteinproductionoutsidegraph',K.CATALYST,R.CATALYST),inp('deamination_condition','Source5ul2xDRB+5ulpreparedDddBstock; wholechromatinoccupancy/productdistributionunassembled')),
       product='Sourcecytosine-to-uracil DNAclasses plusprotected/incompletelydeaminatedtemplates; chemicalU isNOT DNA-T')
    clean('dna_clean','SourceZymoDNACleanConcentrator5purification')
    # Actual source prepared PCR mix, including two catalysts and structural nucleotides/primers.
    mix=g.input('pcr_mix','SourcepreparedQ5U/Bst3/Nexteraprimerreactionstock',K.MOLECULAR)
    e=g.product('q5u','Source1xQ5Umastermixpolymerasecomponent',K.CATALYST);bs=g.product('bst3','SourceBst3.0component .08Uul',K.CATALYST);nts=g.product('pcr_nts','NucleotidecomponentofsamepreparedPCRstock; doseunassigned',K.MOLECULAR);f=g.product('pcr_forward','Source1.25uMNexteraforwardprimer, exactbulkwholebodyunassembled',K.OLIGO);rev=g.product('pcr_reverse','Source1.25uMNexterareverseprimer, exactbulkwholebodyunassembled',K.OLIGO);buf=g.product('pcr_buffer','ConditioncomponentofsamePCRstock',K.CONDITION)
    g.react('pcr_components','Componentrepresentationofsameactualpreparedstock, noextramolecularadditions',S.GENERIC,((mix,R.TEMPLATE),),(e,bs,nts,f,rev,buf),changes_molecule=False)
    op('preamp','Source65C5min/80C12min/98C2min then5cycles98C15s/60C30s/72C1min',S.PCR,((e,R.CATALYST),(bs,R.CATALYST),(nts,R.STRUCTURAL),(f,R.STRUCTURAL),(rev,R.STRUCTURAL),(buf,R.CONDITION)),geometry=PCRGeometry.unavailable(unknown_scene(),reason='Wholehistoricalbulkprimers/landing/indexedproductsunassembled; sourceDNAcopyreadoutTdoesnotmeanchemicaldeaminationproducedT.'),product='SourcepreamplifiedDNA/copyclasses andcopresentprimers/nucleotides/catalysts; originalU-containingtemplates mayremain')
    production=g.product('remaining_pcr','Actualremainingpreamplifiedreactionafteraliquotsampling; volumeunassigned',K.MOLECULAR);qc=g.product('qpcr_aliquot','ActualcollectedqPCRaliquot; volumeunassigned',K.MOLECULAR)
    g.react('qpcr_split','Sourcesamplesaliquotafter5cycles; remainingreactionretained',S.SPLIT,((x,R.TEMPLATE),),(production,qc));g.mark_terminal(Terminal(qc.id,Disposition.QC,'Source20cycleqPCRmonitoring determinesadditionalcycles; wholeqPCRpreparationunassembled, notproductionlibrary.'));x=production
    template=g.product('remaining_template','Templatecomponentofsamereactionremainder',K.MOLECULAR);ce=g.product('remaining_enzymes','Co-retainedQ5U/Bst3catalysts fromsamereaction, activityunassigned',K.CATALYST);cn=g.product('remaining_nts','Co-retainednucleotidesfromsamereaction, amountsunassigned',K.MOLECULAR);cf=g.product('remaining_primers','Co-retainedprimercomponentsofsamereaction',K.OLIGO);cb=g.product('remaining_buffer','Co-retainedconditioncomponentofsamereaction',K.CONDITION)
    g.react('remaining_components','Representcopresentcomponentsofremainder, notnewstockadditionsorphysicalpurification',S.GENERIC,((production,R.TEMPLATE),),(template,ce,cn,cf,cb),changes_molecule=False);x=template
    op('complete_pcr','SourceadditionalcyclesguidedbyqPCR; generally10–12TOTALnot10–12extra',S.PCR,((ce,R.CATALYST),(cn,R.STRUCTURAL),(cf,R.STRUCTURAL),(cb,R.CONDITION)),geometry=PCRGeometry.unavailable(unknown_scene(),reason='Wholebulkindexedproducts/landingandactualselectedcyclecountunassembled; originalprimers/catalysts reusedfromretainedmix.'))
    clean('size_select','SourceZymoSelectaSizeDNACleanConcentrator retainabove200bp; recovery/carryoverunassigned')
    selected=g.product('loading_aliquot','Actualselectedloadingfraction; volumeunassigned',K.MOLECULAR);remaining=g.product('library_remainder','Untransferredpurifiedlibrarystock; actualfateunassigned',K.MOLECULAR)
    g.react('loading_split','Sourcephysicalsequencerloadingallocation',S.SPLIT,((x,R.TEMPLATE),),(selected,remaining));g.mark_terminal(Terminal(remaining.id,Disposition.UNRESOLVED,'Sourceuntransferredstocknotautomaticallydiscardedorpooled.'));x=selected
    lib=g.product('library','SelectedfreshbulkFOODIEIlluminalibrary',K.LIBRARY);g.react('load','SourcebulkFOODIEordinary2x150PEreadout',S.GENERIC,((x,R.TEMPLATE),),(lib,),changes_molecule=False)
    g.sequence(lib,SequencingPlan(limitation='Sourceordinarybulk2x150PE, separatecooperativitylibrary2x250PEnotmerged. Fullbulkindexidentities/primerlanding/runmaskunassembled. ProtectedC/unconvertedstatesandcopiedTreadoutsarenotdirectTFoccupancyproof; nakedDNAcontrols/mapping/deduplicationareseparatecontrols/computationaloperations, notphysicalcleanup.'))
    g.add_source_item(SourceItem('branches','CustomTableS2singlecellTn5oligosareknownlocalcomponents butnotidentitiesofbulkS111stock.','High-throughputbranch96tagwells/pool/deaminate/10nucleipersecondwell/indexPCRwithcustomprimers remainsoutside selectedbulkgraph. Mainliteral10mMprimer/dosesnottransferredtobulk.'));return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph()
    single_cell_adapter('p5');single_cell_adapter('p7_1')
    c('Onlycustomtail/indexdepletednotME','C'not in P5_TAIL+P7_TAIL+TN5_BC1 and'C'in ME,True)
    c('SourceTn5barcode1length',len(TN5_BC1),11)
    c('Deaminationfollowssourcepelletworkup',g.steps['deaminate'].inputs[0].material,'tag_stop_kept')
    c('ActualqPCRsampletakenafter5cyclepreamp',g.steps['qpcr_split'].inputs[0].material,'preamp_product')
    c('Productioncontinuesfromactualremainder',g.steps['remaining_components'].inputs[0].material,'remaining_pcr')
    c('NoadditionalpreparedstockincontinuedPCR',g.steps['complete_pcr'].inputs[1].material,'remaining_enzymes')
    c('qPCRaliquotisQC',g.terminals['qpcr_aliquot'].disposition,Disposition.QC)
    c('Bulkadapterisnotborrowedcustomcomponent',g.materials['tn5_adapters'].state,None)
    p=sources_root()/Path(__file__).resolve().parents[1].name/'pnas.2423270121_supp_pnas.2423270121.sapp.pdf.txt'
    if have(c,Source(p,'https://doi.org/10.1073/pnas.2423270121')):
        t=''.join(p.read_text().split())
        for k,v in single_cell_orderings().items():c('Publishedwholelocaladapter '+k,v.replace('/5Phos/','') in t,True)
    c('Privatepartial',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
