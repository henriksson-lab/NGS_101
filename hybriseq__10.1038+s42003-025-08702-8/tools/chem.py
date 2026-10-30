"""Private HybriSeq reporter preparation; published source-owned probe/barcode inventory."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,feature,revcomp
from illumina import TRUSEQ_READ1,TRUSEQ_READ2
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
TITLE='HybriSeq: selected standard96well DNA-probe reporter preparation'
AUTHORITY='PublishedPMC12365068main125–190 and MOESM3 Data1/2/6–9 selected originalwholeoligos. Twooriginals; sourceSI1/predecessors/fullproductionPCRorderings unreviewed.'
LEFT=TRUSEQ_READ2[14:]
RIGHT='GCTTTAAGGCCGGTCCTAGC'
L='CAGAAATGCTGCTGCAAGAGTTAGAGCGCC';RGT='TATGTAGCTGTGGCTGGCTCATGAGTTGCT'
LINK1='CGAATGCTCTGGCCTGCTAGGACCGGCCTT';LINK2='CTCAAGCACGTGGATAGTCGTACGCCGATG'
Q1=revcomp(LINK1[15:]);Q2=revcomp(LINK2[15:])
BC1='AGGCCAGAGCATTCG'+'CAGAAGA'+'CATCGGCGTACGACT'
BC2='ATCCACGTGCTTGAG'+'CAGAAGA'+'GTGGCCGATGTTTCG'
BC3=TRUSEQ_READ1[9:]+'GACCAAA'+'CGAAACATCGGCCAC'
SOURCE_ORDERINGS={'left_common':LEFT,'right_common':'N'*8+RIGHT,'linker1':LINK1+'/3Phos/','linker2':LINK2+'/3Phos/','quencher1':Q1,'quencher2':Q2,'bc1_A1':BC1,'bc2_A1':BC2,'bc3_A1':BC3}
def probe_scene():
    s=Scene();s.strand('left',[Segment('left_handle',LEFT),Segment('left_target',L)])
    s.anneal('rna',[Segment('rna_right',revcomp(RGT).replace('T','U'),inferred=True),Segment('rna_left',revcomp(L).replace('T','U'),inferred=True)],to='left',pair=('rna_left','left_target'),unpaired=('rna_right',))
    s.anneal('right',[Segment('right_target',RGT),Segment('probe_umi','N'*8,placeholder=True,feature=feature('hybri_probe_umi','umi','random')),Segment('right_handle',RIGHT)],to='rna',pair=('right_target','rna_right'),unpaired=('probe_umi','right_handle'),mod5='P')
    s.ligation_end('left','left_target',kind='splinted',polarity='3-prime',owner='substrate',text='SourceadjacentDNAprobe3OH; RNAcomplement is inferredfromprovidedprobe arms')
    s.ligation_end('right','right_target',kind='splinted',polarity='5-prime',owner='adapter',text='SourcePNKrightprobe5P; wholebiologicaltarget/junctionunassembled');return s

def barcode_scene(n):
    pre,tail,link=('AGGCCAGAGCATTCG','CATCGGCGTACGACT',LINK1)if n==1 else('ATCCACGTGCTTGAG','GTGGCCGATGTTTCG',LINK2)
    s=Scene();s.strand('barcode',[Segment('bc_prefix',pre),Segment('cbc','N'*7,placeholder=True,feature=feature('hybri_bc'+str(n),'cell_barcode','unknown')),Segment('bc_tail',tail)],mod5='P')
    s.anneal('linker',[Segment('link_stem',link[:15]),Segment('link_overhang',link[15:])],to='barcode',pair=('link_stem','bc_prefix'),unpaired=('link_overhang',),mod3='P')
    s.ligation_end('barcode','bc_prefix',kind='splinted',polarity='5-prime',owner='adapter',text='SourcePNK5P donor; linker3P blocks its own3end, partialacceptorcontextunassembled')
    return s

def unknown():
    s=Scene();s.strand('unassembled',[Segment('unassembled','[unassembled]',placeholder=True)]);return s

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=AUTHORITY,representation='published_instance');x=g.input('cells','Sourcefilteredcellpreparation for standard96wellHybriseq; exactsampleandprobe-panelallocationunassigned',K.BIOLOGICAL,primary=True)
    def inp(k,label,kind=K.CATALYST,role=R.CATALYST,**kw):return(g.input(k,label,kind,**kw),role)
    def cond(k,label):return inp(k,label,K.CONDITION,R.CONDITION)
    def op(k,label,kind,ins=(),geometry=None):
        nonlocal x
        y=g.product(k+'_product','Sourcecopresentcells/RNA/DNAreporters/reagents, residualclasses retained untilworkup',K.MOLECULAR);kw={}
        if geometry is not None:kw['geometry']=geometry
        if kind==S.MAGIC:kw['reason']='Publishedoperationpositive; completeplateallocation/sourcepoolmembers/junctions or stockcomposition unassembled, notabsence.'
        g.react(k,label,kind,((x,R.TEMPLATE),)+tuple(ins),(y,),**kw);x=y
    def clean(k,label):
        nonlocal x
        y=g.product(k+'_kept','Sourcephysicallyretainedpreparation;recovery/carryoverunassigned',K.MOLECULAR);w=g.product(k+'_waste','Sourcephysicalworkupfractions',K.WASTE);g.react(k,label,S.CLEANUP,((x,R.TEMPLATE),cond(k+'_condition',label)),(y,w));g.discard(w,'Physicalwashing,notperfectmolecularspecificity.');x=y
    op('fix','Source4percentPFA30minRT at.5millioncells/ml',S.CROSS_LINKING,(cond('pfa','Sourceparaformaldehyde'),));clean('fix_wash','SourcePBSwashes/filter40um/count')
    probes=inp('prepared_probes','Source10nMperprobe preparedleft/rightDNApool; rightPNK5P/column preparationknown, fullpanelmembershipunassembled;KIF23pair islocalinventoryexample only',K.ADAPTER,R.STRUCTURAL,chemical_state={'source_representative_left':LEFT+L,'source_representative_right_5P':RGT+'N'*8+RIGHT,'umi':'rightprobe8nt','source_pool_scope':'otherpanelstrands/selectedcellRNAcontextunassembled'})
    op('hybridize','Sourcehybridizationbufferpreincubation37C10min then10nMprobes37C18–24h',S.HYBRIDIZATION,(probes,cond('hyb_buffer','Source30percentformamide/BSA/Tween/2xSSC/RNasin')))
    clean('hyb_wash','Source2x20percentformamidewash37C5min thenligationbufferwash')
    # Explicit same-tube retained reporter component, not fresh probes after washing.
    t=g.product('cell_component','Same-tubecellularRNAcomponent',K.MOLECULAR);p=g.product('retained_probes','RetainedhybridizedDNAprobecomponents;occupancyunassigned',K.ADAPTER);g.react('probe_components','Same-tubecomponentsnotnewaddition/physicalfractionation',S.GENERIC,((x,R.TEMPLATE),),(t,p),changes_molecule=False);x=t
    op('probe_ligate','Source2uMSplintR37C1h',S.LIGATION,((p,R.STRUCTURAL),inp('splintr','SourceSplintRligase'),cond('probe_lig_buffer','SourceT4ligasebuffer/.4mMATP/RNasin')),AdapterLigationGeometry.known('splinted',probe_scene(),statement='SourceDNAprobe3OH/PNK5P/adjacenttargetdesign;RNAcomplementcandidateinferred, fullreporterproductunassembled.'))
    for n in (1,2):
        adapters=inp('round'+str(n)+'_adapters','Source96wellprepared5Pbarcode/linkerplate,2.31ul5.4uMperwell; drawnfamilycomponent,fullpooldistributionunassembled',K.ADAPTER,R.STRUCTURAL,state=barcode_scene(n),chemical_state={'representative_A1':BC1 if n==1 else BC2,'linker_3P':LINK1 if n==1 else LINK2,'cell_barcode_nt':7,'source_linker_geometry':'barcode_scene('+str(n)+') localfamilyinterface,notwholeplatemixture'})
        prep=inp('round'+str(n)+'_preparation','SourceactualT4ligation/quenching/washstockpreparations,fullphysicalallocationsunassembled',K.MOLECULAR,R.STRUCTURAL)
        op('round'+str(n),'Source96wellphysicalsplit22.69ulcells/500perwell+2.31uladapter;25C2h/2ul62.5uMquencher30min/25ulEDTAwash10min;physicalpool/pellet/washes',S.MAGIC,(adapters,prep,))
        # These are real pool-wide operations; individual96well/remainingstock fractions are unassembled, not falsely one-input POOL claims.
    op('round3','SourcefinalcellpoolThermoPolnormalized23k/ml;460cells/well,round3well-specificreverse andforwardprimersKAPA20cycles',S.MAGIC,(inp('round3_stock','SourceKAPAReadyMix+forward/round3reverseprimerpreparation; samebulkenzyme/nucleotidecomponentsunassembled',K.MOLECULAR,R.STRUCTURAL,chemical_state={'representative_A1_reverse':BC3,'source_pcr_cycles':20}),))
    supern=g.product('round3_supernatant','SourcecollectedPCRreactionDNA solutionsfromallround3wells; fullpoolallocationunassembled',K.MOLECULAR);pellet=g.product('round3_cells','Source500g1mincellpelletfraction; laterfateunassigned',K.MOLECULAR);g.react('round3_collect','SourcePCRsolutioncollectionaftercellpelleting',S.SPLIT,((x,R.TEMPLATE),),(supern,pellet));g.mark_terminal(Terminal(pellet.id,Disposition.UNRESOLVED,'Separatedcellsnotsequencingreads.'));x=supern
    clean('round3_clean','SourceallcollectedPCRsolutionspooled/columnpurified; detailedmultiwellpoolpartitionunassembled')
    q=g.product('tapestation','SourceTapeStationQCaliquot,amountunassigned',K.MOLECULAR);k=g.product('qc_kept','RetainedpooledPCRstock',K.MOLECULAR);g.react('qc','SourceanalyticalTapeStationQCseparatefromproduction',S.SPLIT,((x,R.TEMPLATE),),(q,k));g.mark_terminal(Terminal(q.id,Disposition.QC,'AnalyticalQCnotreturned.'));x=k
    op('conditional_gel','Sourceconditional3percentagarose200V20min/correctsizegelrecovery whenrequired;252±2or232±2bpdependsonleftprimerbody',S.MAGIC,(cond('gel_choice','Sourcebandpurity/leftprobearchitecturecondition; selectedbranchunassigned'),))
    a=g.product('adapter_pcr_input','Source.5ngpooledround3productforadapterPCR',K.MOLECULAR);r=g.product('pcr_stock_remainder','Untransferredround3DNAstock,fateunassigned',K.MOLECULAR);g.react('adapter_pcr_aliquot','Source.5ngaliquot',S.SPLIT,((x,R.TEMPLATE),),(a,r));g.mark_terminal(Terminal(r.id,Disposition.UNRESOLVED,'Actualstocknotallconsumed.'));x=a
    op('adapter_pcr','SourceP7forward/P5reverse.4uMeach/KAPAReadyMix10cycles',S.MAGIC,(inp('adapter_pcr_stock','SourceKAPAReadyMix/P5/P7primerpreparation; wholeprimerorderingsandcomponentamountsunassembled',K.MOLECULAR,R.STRUCTURAL),))
    clean('final_clean','Source.8xAMPureXP')
    lib=g.product('library','SourceIllumina-compatibleHybriSeqDNAreporterpreparation',K.LIBRARY);g.react('load','Source15pMsequencerloading',S.GENERIC,((x,R.TEMPLATE),cond('loading','Source15pM')),(lib,),changes_molecule=False);g.sequence(lib,SequencingPlan(limitation='SourceMiSeq75ntR1cellbarcode/75ntR2probeandUMI. Fullproductionprimers/entireamplicon/readlandingunassembled; DNAreporters,notRNAsequencevariants,aremeasured.'))
    g.add_source_item(SourceItem('pending','Whole96wellfractions/sealedreporters/allpanelstrands/productionP5P7 andreadlandingunassembled.SI1ordinarydeferred; sourcepositiveprivatepartialscope.','Noabsenceorcompletechemistryclaim.'));return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph();c('ReporterbranchdoesnotinventRT',all(s.kind!=S.REVERSE_TRANSCRIPTION for s in g.steps.values()),True)
    c('Ligaseusesretainedprobes',any(p.material=='retained_probes'for p in g.steps['probe_ligate'].inputs),True)
    c('RNApairingcandidateisinferred',all(x.inferred for x in probe_scene().strands['rna'].segs),True)
    c('AdapterPCRretainsstock',g.terminals['pcr_stock_remainder'].disposition,Disposition.UNRESOLVED)
    c('Cellpelletnotsequencingreadfilter',g.terminals['round3_cells'].disposition,Disposition.UNRESOLVED)
    p=sources_root()/Path(__file__).resolve().parents[1].name/'s42003-025-08702-8_supp_42003_2025_8702_MOESM3_ESM.xlsx.txt'
    if have(c,Source(p,'https://doi.org/10.1038/s42003-025-08702-8')):
        text=p.read_text()
        for name,seq in SOURCE_ORDERINGS.items():c('Originalwholeinventory '+name,seq in text,True)
        c('OriginalKIF23wholepair',L+' ,'+RGT in text,True)
    c('Privateincompletegraph',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
