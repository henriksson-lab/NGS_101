from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,revcomp
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition

def unknown_scene(label='Unassembled whole molecular product'):
    sc=Scene();sc.strand('unassembled',[Segment('unassembled','[unassembled]',placeholder=True)]);return sc

def add(g,key,label,kind=K.MOLECULAR,role=R.STRUCTURAL,**kw):
    return (g.input(key,label,kind,**kw),role)
def condition(g,key,label):return add(g,key,label,K.CONDITION,R.CONDITION)
def opaque(g,key,label,x,ins=(),kind=S.MAGIC):
    y=g.product(key+'_product',label+'; whole molecular identities/junctions unassembled',K.MOLECULAR)
    kw={'reason':'Source-positive operation; whole kit composition, substrate context or product geometry unassembled, not publication absence.'} if kind==S.MAGIC else {}
    g.react(key,label,kind,((x,R.TEMPLATE),)+tuple(ins),(y,),**kw);return y
def cleanup(g,key,label,x):
    y=g.product(key+'_kept','Source retained preparation; recovery/carryover unassigned',K.MOLECULAR);w=g.product(key+'_removed','Source workup fraction',K.WASTE)
    g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),condition(g,key+'_condition',label)),(y,w));g.discard(w,'Physical workup, not perfect molecular selectivity.');return y
def aliquot(g,key,label,x,remainder_label='Untransferred source stock; later fate unassigned'):
    y=g.product(key+'_chosen',label,K.MOLECULAR);r=g.product(key+'_remainder',remainder_label,K.MOLECULAR)
    g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(y,r));g.mark_terminal(Terminal(r.id,Disposition.UNRESOLVED,remainder_label));return y
def qc(g,key,label,x):
    a=g.product(key+'_analytical',label,K.MOLECULAR);k=g.product(key+'_kept','Retained production stock',K.MOLECULAR)
    g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(a,k));g.mark_terminal(Terminal(a.id,Disposition.QC,'Analytical aliquot not returned; amount unassigned.'));return k
from illumina import P5,TRUSEQ_READ2
TITLE='RoCK and ROI: published selected eGFP capture/ROI preparation'
AUTHORITY='PublishedPMC12695968 main288–352 and SI1Table11p41 original modified oligos. Two originals. Protocols.io explicitly citedv2 and exact BD vendor guides ordinary unreviewed; no historical whole-bead or modern BD constructor inherited.'
SPLINT_SUFFIX='CATACCTACTACGCATA'
EGFP_CAPTURE_TEMPLATE='CATGGTCCTGCTGGAGTTCGTGACC'
TSO_END=revcomp(SPLINT_SUFFIX)
T_PRIMER='ACAGGAAACTCATGGTGCGT'
INDEX_PRIMER=P5+'TCTTTCCCTACAC'+T_PRIMER
SEQ_PRIMER='ACACTCTTTCCCTACAC'+T_PRIMER
ROI=TRUSEQ_READ2[11:]+'GGCAACTACAAG'
SOURCE_ORDERINGS={'splint':'/5Phos/'+EGFP_CAPTURE_TEMPLATE+SPLINT_SUFFIX,'polyA':'/5Phos/'+'A'*18,'tso_protection':'/5Phos/'+SPLINT_SUFFIX,'Tprimer':T_PRIMER,'index':INDEX_PRIMER,'customRead1':SEQ_PRIMER,'eGFP_ROI':ROI}
def strand(seq,name='oligo',mod5=''):
    sc=Scene();sc.strand(name,[Segment(name,seq)],mod5=mod5);return sc
def splint_scene():
    sc=Scene();sc.strand('native_bead_end',[Segment('source_tso_end',TSO_END)])
    sc.anneal('splint',[Segment('capture_template',EGFP_CAPTURE_TEMPLATE),Segment('tso_complement',SPLINT_SUFFIX)],to='native_bead_end',pair=('tso_complement','source_tso_end'),unpaired=('capture_template',),mod5='P');return sc

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=AUTHORITY,representation='published_instance')
    raw=g.input('raw_beads','Source full vial EnhancedCellCaptureBeadsV2; complete bead-linked upstream barcodes/UMIs/dT/TSO constructs unassembled',K.BIOLOGICAL,primary=True)
    tubes=[g.product('bead_tube_'+str(i),'Source500ul bead vial '+str(i+1),K.BIOLOGICAL)for i in range(4)]
    g.react('four_vials','Source full bead vial split into four500ul vials processed in parallel',S.SPLIT,((raw,R.TEMPLATE),),tuple(tubes))
    modified=[]
    for i,x in enumerate(tubes):
        p='v'+str(i)+'_'
        x=cleanup(g,p+'wash','Source WaterBuffer thenTE/TW wash',x)
        x=opaque(g,p+'anneal','Source5min37C300rpm preincubation WITHOUT T4polymerase',x,(add(g,p+'splint','Source5P eGFP splint, representative chosen target; full strand in originalTable11',K.OLIGO,R.STRUCTURAL,state=strand(EGFP_CAPTURE_TEMPLATE+SPLINT_SUFFIX,'splint','P')),add(g,p+'protect','Source5P A18 polyA protective oligo',K.OLIGO,R.STRUCTURAL,state=strand('A'*18,'protect','P')),add(g,p+'nt','Source130ul10mMdNTPs per production vial; final1mM wording retained without per-base inference'),condition(g,p+'buffer','Source260ul5xbuffer/857ulwater/6ul1:1splint:protectmix; stated50uMeach convention unassigned')),S.HYBRIDIZATION)
        # Whole bead-bound input/annealed primer is retained; no second splint/primer stock addition.
        x=opaque(g,p+'extend','Source6ul T4DNApolymerase5U/ul, RT10min; preincubation nucleotide/oligo stocks retained',x,(add(g,p+'t4','Source T4DNApolymerase',K.CATALYST,R.CATALYST),),S.EXTENSION)
        x=opaque(g,p+'inactivate','Source75C10min T4polymerase inactivation',x,(condition(g,p+'thermal_t4','Source75C10min'),),kind=S.THERMAL_TREATMENT)
        x=cleanup(g,p+'post_t4_wash','Source WaterBuffer thenTE/TW wash',x)
        x=opaque(g,p+'lambda','Source21ul lambdaexonuclease5000U/ml in stated95ulbuffer/832ulwater;37C30min',x,(add(g,p+'lambda_stock','Source lambdaexonuclease',K.CATALYST,R.CATALYST),condition(g,p+'lambda_buffer','Source actual lambda buffer/water volumes; buffer-concentration wording retained')),S.DIGESTION)
        x=opaque(g,p+'lambda_stop','Source75C10min lambda inactivation',x,(condition(g,p+'thermal_lambda','Source75C10min'),),kind=S.THERMAL_TREATMENT);x=cleanup(g,p+'storage_wash','Source washes then4C TE/TWstorage',x);modified.append(x)
    # No pooling of the four factory tubes is asserted: source moves from bead factory to experiments without assigning that transfer.
    chosen=g.product('selected_modified_beads','Source RoCKseq beads selected for one eGFP experiment, amount/tube allocation unassigned',K.MOLECULAR)
    remaining=g.product('factory_remaining_stock','Unallocated modified-bead stock across factory tubes',K.MOLECULAR)
    g.react('factory_allocation','Source transfer from prepared-bead stock to assay; exact factory tube allocation unassembled',S.MAGIC,tuple((x,R.TEMPLATE)for x in modified),(chosen,remaining),reason='Source stock/experimental allocation boundary unassembled; not a declared physical pooling.');g.mark_terminal(Terminal(remaining.id,Disposition.UNRESOLVED,'Actual stock not silently all consumed.'))
    cells=g.input('cells','Source single-cell preparation; selected sample not fixed as a particular mixing cohort',K.BIOLOGICAL,primary=True)
    x=opaque(g,'capture_rt','Source selected preparedbeads in680ulSampleBuffer; RNAcapture/cDNAsynthesis delegated to exactBDDoc210966',cells,((chosen,R.STRUCTURAL),add(g,'capture_kit','Source BDcDNA reagents/RT and molecular substrates; whole historical formulation unassembled')))
    roi_stock=g.input('roi_stock','Source1ul100uM eGFP ROIprimer diluted to10ul',K.OLIGO,state=strand(ROI,'roi'))
    a=g.product('roi_four_ul','Source4ul dilutedROI mix for randompriming',K.OLIGO,state=strand(ROI,'roi'));r=g.product('roi_six_ul','Source6ul unused dilutedROI mix',K.OLIGO,state=strand(ROI,'roi'))
    g.react('roi_aliquot','Source4ul from10ul dilutedprimerstock',S.SPLIT,((roi_stock,R.TEMPLATE),),(a,r));g.mark_terminal(Terminal(r.id,Disposition.UNRESOLVED,'Actual diluted primer remainder.'))
    x=opaque(g,'roi_extend','Source randompriming/extension with selectedROIprimer plus standard randomers; no original RNA capture replaced by ROIpriming',x,((a,R.STRUCTURAL),add(g,'rpe_stock','Source BD RandomPrimingMix/randomers/extensionreagents; enzyme/nucleotide components unassembled')))
    x=opaque(g,'rpe_pcr','Source RPEPCR with1ul100uM Tprimer added to standardRPEPCRmix/purifiedRPEproduct',x,(add(g,'t_primer','Source Tprimer',K.OLIGO,R.STRUCTURAL,state=strand(T_PRIMER,'Tprimer')),add(g,'rpe_pcr_stock','Source RPEPCR preparation and standard universal primer; full components unassembled')))
    # Two actual aliquots from common primarylibrary, not two complete copies of its mass.
    wta=g.product('wta_primary','Selected primarylibrary aliquot for WTAindexPCR; amount unassigned',K.MOLECULAR);tso=g.product('tso_primary','Selected primarylibrary aliquot for TSOindexPCR; amount unassigned',K.MOLECULAR);rem=g.product('primary_remaining','Untransferred primarylibrary stock',K.MOLECULAR)
    g.react('primary_split','Source separateWTA andTSO indexingPCRs from same primarylibrary; allocation amounts unassembled',S.SPLIT,((x,R.TEMPLATE),),(wta,tso,rem));g.mark_terminal(Terminal(rem.id,Disposition.UNRESOLVED,'Actual allocation unassigned.'))
    selected=[]
    for label,part in [('wta',wta),('tso',tso)]:
        ins=[add(g,label+'_index_mix','Source standardBDLibraryReverseprimer/indexPCRstock; full formulation unassembled')]
        if label=='tso':ins.append(add(g,'custom_index_primer','Source5ul100uM customTindexprimer replaces5ulLibraryForwardPrimer',K.OLIGO,R.STRUCTURAL,state=strand(INDEX_PRIMER,'index_primer')))
        else:ins.append(add(g,'wta_forward','Source BDLibraryForwardPrimer, whole ordering unassembled',K.OLIGO))
        part=opaque(g,label+'_index','Source separate '+label+'indexPCR; same8bp sample index used for both libraries, value unassigned',part,ins)
        selected.append(aliquot(g,label+'_pool_aliquot','Source selected '+label+' library amount for1:1concentrationpool, volume unassigned',part))
    pool=g.product('paired_library_pool','Source1:1WTA:TSO selectedlibrarypool',K.MOLECULAR);g.react('pool','Source actual two-library pool, not whole-stock consumption',S.POOL,tuple((x,R.TEMPLATE)for x in selected),(pool,))
    lib=g.product('sequencing_library','Source pool plus actualPhiXDNAcontrol',K.LIBRARY)
    g.react('loading','Source3percentPhiX and customRead1primer addition',S.GENERIC,((pool,R.TEMPLATE),add(g,'phix','Source PhiXDNAsequencinglibrary',K.MOLECULAR,R.SPIKE_IN),add(g,'custom_read1','Source customRead1sequencingprimer',K.OLIGO,R.STRUCTURAL,state=strand(SEQ_PRIMER,'custom_read1'))),(lib,),changes_molecule=False)
    g.sequence(lib,SequencingPlan(limitation='Selected cell-line/Pdgfra runNovaSeq6000R1:60/R2:62; leukemia/HeLaR2:180 are separate run branches. Full bead/index/read landing and vendor library structures unassembled; custom primer literal retained without invented landing.'))
    g.add_source_item(SourceItem('unassembled','Protocolsv2 and exactBDDoc210966/23-21711-00 ordinary unreviewed; whole beadfactory/polymerase products and read structures unassembled. FluorescentQC20ul/scaledreactions sourcepositive but not imported as production losses.','Not absence/full chemistry audit.'));return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph();c('Four source factory vials',len(g.steps['four_vials'].outputs),4)
    for i in range(4):
        c('No fresh splint at extension '+str(i),len(g.steps['v'+str(i)+'_extend'].inputs),2)
        c('Lambda follows washed T4 branch '+str(i),list(g.steps).index('v'+str(i)+'_post_t4_wash')<list(g.steps).index('v'+str(i)+'_lambda'),True)
    c('ROI actual unused6ul retained','roi_six_ul'in g.terminals,True)
    c('Both selected pool stocks retain remainders',all(k+'_pool_aliquot_remainder'in g.terminals for k in ('wta','tso')),True)
    c('PhiX contributes actual molecular stock',any(p.material=='phix'and p.role==R.SPIKE_IN for p in g.steps['loading'].inputs),True)
    c('Local nativeTSO:splint pairing computes',bool(splint_scene().strands),True)
    p=sources_root()/Path(__file__).resolve().parents[1].name/'s41467-025-66248-z_supp_41467_2025_66248_MOESM1_ESM.pdf.txt'
    if have(c,Source(p,'https://doi.org/10.1038/s41467-025-66248-z')):
        t=''.join(p.read_text().split())
        for name,seq in SOURCE_ORDERINGS.items():
            raw=seq.replace('/5Phos/','');c('OriginalTable11whole '+name,raw in t,True)
    c('Private partial graph remains incomplete',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
