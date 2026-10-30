"""Private Nm-Mut-seq selected Extended-protocol polyA-RNA downstream preparation."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,revcomp
from illumina import TRUSEQ_READ2,TRUSEQ_SMALL_RNA_RA5
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
TITLE='Nm-Mut-seq: restrictive/permissive Extended-protocol pair'
AUTHORITY='Published SI2 Extended protocols originalPDFp35–38/printed24–27, libraryMethods4197–4254. One original source for recipe/oligos. Main/reagents not automatically inherited from modern NEB kits.'
A3_FIXED='TTAGGC'+revcomp(TRUSEQ_READ2)[:21]
A5_FIXED=TRUSEQ_SMALL_RNA_RA5.replace('T','U')
SOURCE_ORDERINGS={'adapter3':'5rApp-'+'N'*5+A3_FIXED+'-3biotin','adapter5':A5_FIXED+'N'*5}
def adapter(which):
    s=Scene();parts=([Segment('random5','N'*5,placeholder=True),Segment('fixed',A3_FIXED)]if which==3 else[Segment('fixed',A5_FIXED),Segment('random5','N'*5,placeholder=True)])
    s.strand('adapter',parts,mod5='rApp'if which==3 else'',mod3='biotin'if which==3 else'');return s

def ends(which):
    s=adapter(which);s.strand('substrate',[Segment('substrate','[RNA fragment; whole geometry unassembled]',placeholder=True)],mod3='OH'if which==3 else'',mod5='P'if which==5 else'')
    s.ligation_end('adapter','random5'if which==3 else'fixed',through='fixed'if which==3 else'random5',kind='single-strand',polarity='5-prime'if which==3 else'3-prime',owner='adapter',text='Source preadenylated5end/biotin3block'if which==3 else'Source unmodified RNA3end')
    s.ligation_end('substrate','substrate',kind='single-strand',polarity='3-prime'if which==3 else'5-prime',owner='substrate',text='Source PNK repaired RNA; whole sealed junction unassembled');return s

def unknown():
    s=Scene();s.strand('unknown',[Segment('unassembled','[unassembled]',placeholder=True)]);return s

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=AUTHORITY,representation='published_instance');x=g.input('polya_rna','Source1ug polyA RNA at ExtendedDay2 boundary; upstream purification not assembled',K.MOLECULAR,primary=True)
    def inp(k,label,kind=K.CATALYST,role=R.CATALYST,**kw):return(g.input(k,label,kind,**kw),role)
    def cond(k,label):return inp(k,label,K.CONDITION,R.CONDITION)
    def nt(k,label):return inp(k,label,K.MOLECULAR,R.STRUCTURAL)
    def op(k,label,kind,ins=(),geometry=None,product='Sourcecopresent molecular/reagent mixture, residual classes retained until workup'):
        nonlocal x
        y=g.product(k+'_product',product,K.MOLECULAR);kw={}
        if geometry is not None:kw['geometry']=geometry
        if kind==S.MAGIC:kw['reason']='Source operation positive; whole stock identity/product geometry unassembled, not confirmed absent.'
        if kind==S.GENERIC:kw['changes_molecule']=False
        g.react(k,label,kind,((x,R.TEMPLATE),)+tuple(ins),(y,),**kw);x=y
    def clean(k,label):
        nonlocal x
        y=g.product(k+'_kept','Source retained fraction; recovery/carryover unassigned',K.MOLECULAR);w=g.product(k+'_waste','Source separated workup fractions',K.WASTE);g.react(k,label,S.CLEANUP,((x,R.TEMPLATE),cond(k+'_condition',label)),(y,w));g.discard(w,'Physical workup, not perfect molecular selection.');x=y
    op('fragment','Source45ulRNA+5ul1McarbonatepH9.2,95C8min/ice1min',S.CHEMICAL_CONVERSION,(cond('fragment_buffer','SourceNa2CO3/NaHCO3buffer/heat'),))
    clean('fragment_clean','SourceOCCelute40ul')
    op('pnk1','Source40ulRNA+5ulbuffer+1ulPNK37C30min',S.DEPHOSPHORYLATION,(inp('pnk1_enzyme','SourcePNK'),))
    op('pnk2','Sourceadd5ul10mMATP+1uladditionalPNK37C30min/65C10min',S.PHOSPHORYLATION,(inp('pnk2_enzyme','SourceadditionalPNK'),nt('pnk_atp','SourceATPphosphate donor')))
    clean('pnk_clean','SourceOCC;eluteappropriatevolume>100ng/6ul')
    a=g.product('rna3ul','Source3ulRNA selected for ligation',K.MOLECULAR);r=g.product('rna_stock_remainder','Untransferred PNKRNA stock, volume/fate unassigned',K.MOLECULAR);g.react('rna_aliquot','Source3ulselected aliquot',S.SPLIT,((x,R.TEMPLATE),),(a,r));g.mark_terminal(Terminal(r.id,Disposition.UNRESOLVED,'Untransferred stock.'));x=a
    ad=g.input('adapter3','Source0.5ul12uM '+SOURCE_ORDERINGS['adapter3'],K.ADAPTER,state=adapter(3),chemical_state={'source_ordering':SOURCE_ORDERINGS['adapter3'],'five_prime':'rApp','three_prime':'biotin'})
    op('anneal3','Source70C2min/ice with3adapter',S.HYBRIDIZATION,((ad,R.STRUCTURAL),))
    t=g.product('rna_component','Same-tube RNA component',K.MOLECULAR);ad=g.product('retained_adapter3','Same-tube retained3adapter',K.ADAPTER,state=adapter(3));g.react('components3','Represent same-tube components,notphysicalfractions',S.GENERIC,((x,R.TEMPLATE),),(t,ad),changes_molecule=False);x=t
    op('ligate3','Source5ul3ligationbuffer+1.5ulligationenzyme25C2h',S.LIGATION,((ad,R.STRUCTURAL),inp('ligase3','Source historical kit3ligationenzyme; exact identity unassembled')),AdapterLigationGeometry.known('single-strand',ends(3),statement='Known adapter5rApp/RNA3OH local interface; whole sealed product unassembled.'))
    primer=g.input('sr_rt_primer','Source0.5ulSRRTprimer dilutedto2.75ul; wholeordering unassembled',K.OLIGO)
    op('anneal_rt','Source75C5min37C15min25C15min,addSRRTprimer',S.HYBRIDIZATION,((primer,R.STRUCTURAL),))
    ad=g.input('adapter5','Source1ul12uM RNA5adapter separately70C2min/ice '+SOURCE_ORDERINGS['adapter5'],K.ADAPTER,state=adapter(5),chemical_state={'source_ordering':SOURCE_ORDERINGS['adapter5'],'sugar':'RNA'})
    op('ligate5','Sourceadd1ul5adapter/.5ulbuffer/1.25ulligationenzyme25C2horovernight',S.LIGATION,((ad,R.STRUCTURAL),inp('ligase5','Sourcehistoricalkit5ligationenzyme; exactidentityunassembled')),AdapterLigationGeometry.known('single-strand',ends(5),statement='Known PNK-RNA5P and RNA-adapter3OH local interface; whole sealed product unassembled.'))
    clean('ligated_clean','SourceOCCelute21ul; RTprimer occupancy/recovery unassigned')
    arms={a:g.product(a+'_rna','Source10ulRNA preparation for'+a+'RT; actual table10ul versus equal-halves wording retained',K.MOLECULAR)for a in('permissive','restrictive')};rem=g.product('rt_rna_remainder','Nominal1ul from21ul after two10ul transfers; equal-halves wording reconciliation unassigned',K.MOLECULAR);g.react('rt_split','Source table two10ul aliquots; text says equal halves of21ul, not silently corrected',S.SPLIT,((x,R.TEMPLATE),),(*arms.values(),rem));g.mark_terminal(Terminal(rem.id,Disposition.UNRESOLVED,'Source allocation wording unresolved.'))
    for arm,start in arms.items():
        x=start
        t=g.product(arm+'_template','Same-tube ligatedRNA component',K.MOLECULAR);p=g.product(arm+'_retained_primer','Previously added SRRTprimer component; occupancy/recovery unassigned, no new stock',K.OLIGO);g.react(arm+'_components','Same-tube components from anneal_rt provenance; no physicalfractionation',S.GENERIC,((x,R.TEMPLATE),),(t,p),changes_molecule=False);x=t
        ins=[(p,R.STRUCTURAL),inp(arm+'_rt_enzyme','Source2ul10uM evolvedRT41B4'),nt(arm+'_dntp','Source2ul10mMdNTP'if arm=='permissive'else'Source2ul20uMdNTP; source mixture dilution convention retained'),cond(arm+'_rt_buffer','Source4ul5xRTbuffer/.5ulRNaseinhibitor; no new primer addition')]
        if arm=='restrictive':ins.append(nt(arm+'_datp','Sourceadditional2ul10mMdATP'))
        op(arm+'_rt','Source37C1h80C10min,restrictive/permissive dNTP compositions',S.REVERSE_TRANSCRIPTION,ins)
        qc=g.product(arm+'_qc','Source1ulcDNAforanalyticalqPCR',K.MOLECULAR);pc=g.product(arm+'_pcr_input','Source10ulcDNAforproductionPCR',K.MOLECULAR);rem=g.product(arm+'_cdna_remainder','UntransferredRTstock remainder; actualvolumeunassigned',K.MOLECULAR);g.react(arm+'_cdna_split','Source1ulQC/10ulPCRallocations',S.SPLIT,((x,R.TEMPLATE),),(qc,pc,rem));g.mark_terminal(Terminal(qc.id,Disposition.QC,'AnalyticalqPCRnotreturnedtolibrary.'));g.mark_terminal(Terminal(rem.id,Disposition.UNRESOLVED,'Untransferredstock.'));x=pc
        mix=g.input(arm+'_longamp','Source12.5ul2xLongAmpTaqbulkstock',K.MOLECULAR);e=g.product(arm+'_pcr_enzyme','SamebulkLongAmpTaqcatalystcomponent',K.CATALYST);n=g.product(arm+'_pcr_nt','Samebulk molecular nucleotide substrate component',K.MOLECULAR);b=g.product(arm+'_pcr_buffer','Samebulkbuffercomponent',K.CONDITION);g.react(arm+'_pcr_components','Representoneactualbulkaddition,notthreeindependentstocks',S.GENERIC,((mix,R.TEMPLATE),),(e,n,b),changes_molecule=False)
        op(arm+'_pcr','Source.625uleachSR/indexprimer+14.38ulPCRmix; prose15cyclespermissive/21restrictive versus adjacentchart12–15 retained',S.PCR,((e,R.CATALYST),(n,R.STRUCTURAL),(b,R.CONDITION),inp(arm+'_sr','SourceSRprimer,wholeorderingunassembled',K.OLIGO,R.STRUCTURAL),inp(arm+'_index','Sourceindexedprimer,wholeordering/indexunassigned',K.OLIGO,R.STRUCTURAL)),PCRGeometry.unavailable(unknown(),reason='Fullhistoricalprimerorderings/ampliconlandingunassembled; notabsenceclaim.'))
        clean(arm+'_gel','Source2.5percentlowmeltagarose90V90min/select160–200bp/MinElute13ul')
        q=g.product(arm+'_qubit','Source1ulfinalQubit aliquot',K.MOLECULAR);lib=g.product(arm+'_library','Source remainingNmMutlibrary afterQubit',K.LIBRARY);g.react(arm+'_final_split','SourcefinalQCaliquotandremaininglibrary',S.SPLIT,((x,R.TEMPLATE),),(q,lib));g.mark_terminal(Terminal(q.id,Disposition.QC,'SourceQubitnotreturned.'));g.sequence(lib,SequencingPlan(limitation='SourceNextSeq550SE75.Wholelibrary/readlandingandtwoN5randomwindowmeaningunassembled; calibration/reportingnotchemicalyieldvalidation.'))
    g.add_source_item(SourceItem('pending','SelectedExtendedDay2/3only. SIshortMethodsstart50ug→100ugintermediate differs fromExtended50ug; not imported upstream. Actual21ul/equal-halves/table10ul and PCR15/21vs12–15wordings retained.','No confirmedabsence or complete source review.'));return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph();c('RNA5adapterliteralUretained','U'in ''.join(x.top for x in adapter(5).strands['adapter'].segs),True)
    c('RTprimerisretainednotnewstock',all(g.producer_of(p.material)is not None for p in g.steps['restrictive_rt'].inputs if p.role==R.STRUCTURAL and p.material.endswith('retained_primer')),True)
    c('RestrictivedATPismolecular',g.materials['restrictive_datp'].kind,K.MOLECULAR)
    c('PNKsecondadditionexplicit',any(p.material=='pnk2_enzyme'for p in g.steps['pnk2'].inputs),True)
    c('SourceRNAremainderexplicit',g.terminals['rt_rna_remainder'].disposition,Disposition.UNRESOLVED)
    c('QCnotsequencing',g.terminals['restrictive_qc'].disposition,Disposition.QC)
    path=sources_root()/Path(__file__).resolve().parents[1].name/'s41422-023-00836-w_41422_2023_836_MOESM2_ESM.pdf.txt'
    if have(c,Source(path,'https://doi.org/10.1038/s41422-023-00836-w')):
        text=path.read_text();c('Originalwhole3adapter','5’rApp-NN NNN TTA GGC AG ATC GGA AGA GCA CAC GTC T-3’biotin'in text,True);c('Originalwhole5adapter','5’-GU UCA GAG UUC UAC AGU CCG ACG AUC NNN NN-3’'in text,True)
        c('Modeledwhole3ordering','N'*5+A3_FIXED,'NNNNNTTAGGCAGATCGGAAGAGCACACGTCT');c('Modeledwhole5ordering',A5_FIXED+'N'*5,'GUUCAGAGUUCUACAGUCCGACGAUCNNNNN')
    c('Privateincompletegraph',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
