"""Private source-qualified HiCuT contact generation and directed tagmentation."""
from pathlib import Path
import sys,re
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[1]/'lib'))
from chemdraw import Scene,Segment,feature
from illumina import P5,P7
from nextera import S5,S7,ME
from pcr import PCRGeometry
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,Terminal,Disposition,SequencingPlan,SourceItem
AUTHORITY='PublishedPMC8979432main143/147/151/153/157; publishedS4Tableactualp1. Closed2originals.'
I5='TAGATCGC'
I7=('TCGCCTTA','CTAGTACG','TTCTGCCT','GCTCAGGA','AGGAGTCC','CATGCCTA','GTAGAGAG','CCTCTCTG','AGCGTAGC','CAGCCTCG','TGCCTCTT','TCCTCTAC')
SOURCE_OLIGO_EXCEPTIONS={'universal_i5':'ActualpublishedS4p1 universal i5 Sequencecell wraps after55nt; GTGTAT remains same Sequencecell. Whole61nt sourcefixture checked and actualoriginalrender inspected. Generic scraper keeps55nt only; no sourcebody/chemistry waiver, no sharedparser change.'}
def primer_scene(which='i5',index=None):
    if which not in ('i5','i7'):raise ValueError('Select sourcei5/i7')
    if index is not None and(index not in I7 or which!='i7'):raise ValueError('Select a source-listed i7 index')
    sc=Scene();sc.strand('primer',[Segment('flowcellarm',P5 if which=='i5'else P7),Segment('sampleindex',I5 if which=='i5'else(index or'NNNNNNNN'),placeholder=which=='i7'and index is None,feature=feature('hicut.'+which,'sample_index','fixed'if which=='i5'else'unknown')),Segment('primerentry',S5 if which=='i5'else S7),Segment('MEpartial',ME[:10]if which=='i5'else ME[:7])]);return sc
def seq(scene):return ''.join(s.top for st in scene.strands.values()for s in st.segs)
def unknown_scene():
    sc=Scene();sc.strand('unassembled',[Segment('Wholelibraryunassembled','[unassembled]',placeholder=True)]);return sc

def protocol_graph():
    g=ProtocolGraph('HiCuT selected CTCF contact/tagmentation preparation',authority=AUTHORITY,representation='published_instance')
    x=g.input('cells','Source100000cells; selectedCTCFpreparation, exactcell/genomeunassigned',K.BIOLOGICAL,primary=True)
    def inp(k,label,kind=K.CONDITION,role=R.CONDITION,**kw):return g.input(k,label,kind,**kw),role
    def op(k,label,ports=(),kind=S.MAGIC,product='Actualsourceheterogeneousmixture; wholeproducts/ends/occupancy/yieldsunassembled',**kw):
        nonlocal x
        y=g.product(k+'_product',product,K.MOLECULAR)
        if kind==S.MAGIC:kw.setdefault('reason','Sourcepositiveoperation; occupiedmolecules/endstates/fullproductsunassembled, noabsenceclaim.')
        if kind==S.GENERIC:kw['changes_molecule']=False
        g.react(k,label,kind,((x,R.TEMPLATE),)+tuple(ports),(y,),**kw);x=y
    def clean(k,label,ports=()):
        nonlocal x
        y=g.product(k+'_kept','Actualretainedproductionfraction; recovery/carryoverunassigned',K.MOLECULAR);w=g.product(k+'_removed','Actualremovedworkupfraction',K.WASTE)
        g.react(k,label,S.CLEANUP,((x,R.TEMPLATE),inp(k+'_conditions',label))+tuple(ports),(y,w));g.discard(w,'Actualremovedphysicalfraction; noperfectspecificityclaim');x=y
    def components(k,label,spec):
        st=g.input(k,label,K.MOLECULAR);outs=[];ports=[]
        for suffix,kind,role,text in spec:
            m=g.product(k+'_'+suffix,text,kind);outs.append(m);ports.append((m,role))
        g.react(k+'_components','CopresentcomponentsofONEactualsourcepreparedstock; nofreshduplicateadditions',S.GENERIC,((st,R.TEMPLATE),),tuple(outs),changes_molecule=False);return tuple(ports)
    def allocate(k,label,qc=False):
        nonlocal x
        u=g.product(k+'_used',label,K.MOLECULAR);r=g.product(k+'_rest','Actualremainingstock; amount/fateunassigned',K.MOLECULAR)
        g.react(k,label,S.SPLIT,((x,R.TEMPLATE),),(u,r));g.mark_terminal(Terminal(u.id if qc else r.id,Disposition.QC if qc else Disposition.UNRESOLVED,label if qc else'Untransferredstocknotautomaticallydiscarded'));x=r if qc else u
    op('formaldehyde','Source.5ml1percentFA10minRT/.2Mglycine5minRT15minice',(inp('fa_supplies','Actualformaldehyde/glycinechemicalsupplies; packagingunassigned',K.MOLECULAR,R.STRUCTURAL),))
    clean('first_wash','Source.05percentBSA/PBSwash/2000g5min',(inp('first_bsa','ActualBSAwashmolecularcarrier',K.MOLECULAR,R.STRUCTURAL),))
    op('dsg','Source3mMDSG500ulPBS40minRTrotation/.4Mglycine5min',(inp('dsg_supplies','ActualDSG/glycinesupplies; packagingunassigned',K.MOLECULAR,R.STRUCTURAL),))
    clean('second_wash','Source.05percentBSA/PBSwash/2000g5min',(inp('second_bsa','ActualfreshBSAwashmolecularcarrier',K.MOLECULAR,R.STRUCTURAL),))
    clean('nuclei','SourceHiClysis30minice/2500g5min/NEBuffer3.1wash')
    op('permeabilize','Source161ulNEBuffer3.1+19ul1percentSDS65C10minice/21.5ul10percentTriton37C15min900rpm',(inp('permeabilization_supplies','ActualSDS/Tritonchemicalsupplies; packagingunassigned',K.MOLECULAR,R.STRUCTURAL),))
    op('digest','Source20ul10UulDdeI/4ul50UulDpnII/2ulNEBuffer;37C3horovernight/65C20min',(inp('ddei','Actual200UDdeI',K.CATALYST,R.CATALYST),inp('dpnii','Actual200UDpnII',K.CATALYST,R.CATALYST),inp('digest_conditions','Actual30son/4minoff900rpm programme')))
    fill=components('fill_master','Actual35ulendfillingmaster; printedlistedamountsnotnormalized',[
        ('polymerase',K.CATALYST,R.CATALYST,'Actual5ulDNApolymeraseI'),('nts',K.MOLECULAR,R.STRUCTURAL,'Actual18.75ul.4mMPLAINdATP andprinted.75ul dTTP/dGTP/dCTP10mMeach; wording/balanceunassigned; nobiotin'),('other',K.CONDITION,R.CONDITION,'Actualmastermixremainingcomposition/volumebalanceunassigned')])
    op('fill','Sourceplainnucleotidefill37C2–3hrotation; notbiotinlabeling',fill)
    lig=components('ligation_master','ONEactual332.5ulligationmaster:60buffer/50Triton/6BSA/2.5ligase/214water',[
        ('ligase',K.CATALYST,R.CATALYST,'Actual2.5ul400UulT4DNAligase'),('bsa',K.MOLECULAR,R.STRUCTURAL,'Actual6ul10mgmlBSAcopresentcarrier'),('conditions',K.CONDITION,R.CONDITION,'Actual60ul10xbuffer/50ul10percentTriton/214ulwater')])
    op('contact_ligate','SourceRT2–3hendoverendrotation',lig)
    clean('exchange','Sourcepellet/wash200ulexchangebuffer/resuspend100ul')
    u=g.product('nuclei_25','Actual25ulof100ulexchangenuclei',K.MOLECULAR);r=g.product('nuclei_75','Actualremaining75ulexchangenuclei',K.MOLECULAR)
    g.react('nuclei_allocation','Actual100ulstock→25ulproduction/75ulremainder',S.SPLIT,((x,R.TEMPLATE),),(u,r));g.mark_terminal(Terminal(r.id,Disposition.UNRESOLVED,'Remaining75ulnotautomaticallyused/discarded'));x=u
    op('cona_bind','Source10ulConAbeadswashed2x100activationbuffer/resuspend100ul/addto25ulnuclei/RT10min',(inp('cona','ActualactivatedConAbeads; nucleiattachmentnotantibodyIP',K.SELECTION,R.SELECTION),),kind=S.BINDING)
    clean('cona_supernatant','Sourcemagnet/removeConAloadingsupernatant; retainattachednuclei')
    op('primary','Source50ulantibodybuffer/1ugCTCFprimary2hRTorovernight4C',(inp('primary_ab','Actual1ugCTCFantibody',K.MOLECULAR,R.STRUCTURAL),),kind=S.BINDING)
    clean('primary_supernatant','Sourcemagnet/removeprimary-incubationsupernatant')
    op('secondary','Source50ullowsaltdigitonin/.25ugsecondaryRT30min',(inp('secondary_ab','Actual.25ugappropriateanti-primarysecondary; selectedidentityunassigned',K.MOLECULAR,R.STRUCTURAL),),kind=S.BINDING)
    clean('secondary_wash','Source2x200ullowsaltdigitoninwash')
    loaded=components('pag_loaded','Actual2.5ul20xCUTANApAG-Tn5LOADEDcomplex; exactoligosunassembled',[
        ('enzyme',K.CATALYST,R.CATALYST,'ActualpAG-Tn5catalystclassofloadedstock'),('oligos',K.OLIGO,R.STRUCTURAL,'Actualloadedoligoclasses; wholebodies/ends/modsunassembled, no vendorconstructorinheritance')])
    op('target_pag','Source50ulhighsaltdigitonin+2.5ulloadedstockRT1h',loaded,kind=S.BINDING)
    clean('pag_wash','Source2x200ulhighsaltdigitoninwash')
    op('tagment','Source50ultagmentationbuffer/10mMMgCl2/37C1h',(inp('tagment_conditions','ActualHEPES/NaCl/spermidine/proteaseinhibitor/Mgprogramme'),))
    clean('tagment_supernatant','SourceREMOVEactualposttagmentationsupernatant; RETAINmaterialforrelease')
    op('release','Source50ulReleaseBuffer/1ul20mgmlPK58C1h68C2h',(inp('release_pk','Actual1ul20mgmlproteinaseK',K.CATALYST,R.CATALYST),inp('release_buffer','Source10mMTAPS8.5/.5percentSDS/BOTH.2mMEDTAand22.5mMEDTA; intendedconcentrationunassigned')))
    released=g.product('release_supernatant','ActualNEWrelease-supernatantCOLLECTEDforDNApurification',K.MOLECULAR);pellet=g.product('release_pellet','Actualresidualmaterialafterrelease; fateunassigned',K.MOLECULAR)
    g.react('collect_release','PhysicallycollectNEWrelease-supernatant; distinctfromremovedtagmentationsupernatant',S.SPLIT,((x,R.TEMPLATE),),(released,pellet));g.mark_terminal(Terminal(pellet.id,Disposition.UNRESOLVED,'ResidualreleasefractionnotPCRinput'));x=released
    clean('dna_clean','SourceZymoD5210/elute21ulPCRreadyDNA')
    f=g.input('universal_i5','Actual2ulsourceuniversal i5',K.OLIGO,state=primer_scene('i5'));i=g.input('indexed_i7','Actualsource10uMi7primer, VOLUMEunassigned/notassumed2ul; studyindexunassigned',K.OLIGO,state=primer_scene('i7'))
    pm=components('pcr_master','ONEactual25ulCUTANAHighFidelity2xPCRMasterMix',[
        ('polymerase',K.CATALYST,R.CATALYST,'ActualPCRpolymeraseclasscopresentinsame25ulstock'),('nts',K.MOLECULAR,R.STRUCTURAL,'ActualPCRnucleotideclasscopresentinsame25ulstock'),('conditions',K.CONDITION,R.CONDITION,'ActualPCRbufferco-componentsinsame25ulstock')])
    op('pcr','Source2uli5/10uMi7volumeunassigned/25ulmaster;58C5min72C5min98C45s/18cycles98C15s60C10s72C1min',((f,R.STRUCTURAL),(i,R.STRUCTURAL))+pm,kind=S.PCR,geometry=PCRGeometry.unavailable(unknown_scene(),reason='Knownprimerbodiesdo notestablishoccupiedTn5template/fullamplicons/readlanding; loadedoligosunassembled.'))
    clean('size_select','SourceAMPureXPsize-selectionratioUNASSIGNED/elute15ulQiagenbuffer')
    allocate('qubit_qc','ActualQubitfluorometeraliquot; volumeunassigned',qc=True)
    allocate('qpcr_qc','ActualqPCRagainstIlluminaprimersaliquot; volumeunassigned',qc=True)
    allocate('loading','ActualNextSeqloadingallocation; volumeunassigned')
    lib=g.product('library','ActualHiCuTNextSeqDNAsequencinglibrary; wholemoleculargeometryunassembled',K.LIBRARY)
    g.react('load','SourceNextSeq75bppairedend',S.GENERIC,((x,R.TEMPLATE),),(lib,),changes_molecule=False)
    g.sequence(lib,SequencingPlan(limitation='SourceNextSeqPE75. WholeloadedTn5strands/fullmaturelibrary/readlanding/indexcyclesunassembled; no biotinpulldown or physicalantibodyIP invented.'))
    g.add_source_item(SourceItem('release_edta','Main153ReleaseBufferlistsboth.2mMEDTAand22.5mMEDTA.','Bothretained; intendedconcentrationunassignedwithoutnormalization/publicationerrorclaim.'))
    g.add_source_item(SourceItem('universal_wrap','OriginalS4p1universali5Sequencecell55nt+GTGTATwrap=61nt; i7_11/12rowsrepeatidentically.','Samecellwrapnotbarcodecolumn/error.12uniquei7; canonicaldecompositiondoesnotestablishloadedvendoroligos.'))
    return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph();p=sources_root()/HERE.parent.name/'journal.pgen.1010121_supp_pgen.1010121.s010.pdf.txt'
    if have(c,Source(p,'https://doi.org/10.1371/journal.pgen.1010121')):
        t=p.read_text();whole=seq(primer_scene());c('Wholewrappeduniversalfixture',bool(re.search(re.escape(whole[:55])+r'\s+i5\s+'+re.escape(whole[55:]),t)),True)
        for n,index in enumerate(I7,1):c('WholeS4i7ordering '+str(n),seq(primer_scene('i7',index))in t,True)
        c('Printedduplicatei7_11retained',t.count(seq(primer_scene('i7',I7[10]))),2)
        c('Printedduplicatei7_12retained',t.count(seq(primer_scene('i7',I7[11]))),2)
    m=p.parent/'journal.pgen.1010121_PMC8979432.xml.txt'
    if have(c,Source(m,'https://doi.org/10.1371/journal.pgen.1010121')):
        t=m.read_text();c('SourceplaindATP','18.75 μL of 0.4 mM dATP'in t,True);c('BothEDTAdescriptors','0.2 mM EDTA 0.5% SDS, 22.5mM EDTA'in t,True)
    try:primer_scene('i7','AAAAAAAA');bad=False
    except ValueError:bad=True
    c('Rejectunlistedindex',bad,True)
    c('Unassignedindexplaceholder',primer_scene('i7').strands['primer'].segs[1].placeholder,True)
    c('Indexsamplemeaning',primer_scene('i7').strands['primer'].segs[1].feature.role,'sample_index')
    c('Nucleiremainderunresolved',g.terminals['nuclei_75'].disposition,Disposition.UNRESOLVED)
    c('ConAuses25ulallocation',g.steps['cona_bind'].inputs[0].material,'nuclei_25')
    c('Releaseusesretainedposttagmentmaterial',g.steps['release'].inputs[0].material,'tagment_supernatant_kept')
    c('PCRDNAfromnewreleasesupernatant',g.steps['dna_clean'].inputs[0].material,'release_supernatant')
    c('Removedtagmentfractiondiscarded',g.terminals['tagment_supernatant_removed'].disposition,Disposition.DISCARD)
    c('PCRNTsfromsame25ulmaster',g.producer_of('pcr_master_nts').id,'pcr_master_components')
    c('FillNTsfromsameactualmaster',g.producer_of('fill_master_nts').id,'fill_master_components')
    c('LigaseBSAcopresent',g.producer_of('ligation_master_bsa').id,'ligation_master_components')
    c('Qubitnotsequenced',g.terminals['qubit_qc_used'].disposition,Disposition.QC)
    c('qPCRnotsequenced',g.terminals['qpcr_qc_used'].disposition,Disposition.QC)
    c('Wholelibraryunassembled',g.materials['library'].state,None)
    c('Privatepartial',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
