"""Protocol-owned, source-qualified private partial MASD preparation."""
TITLE='MASD-seq selected Appendix cfDNA route'
AUTHORITY='Published PMC13658776 main90/96–100/298–310; Appendix PartI126–502(printed14–20); TableEV1 wholeoriginalXLSX includingred/blue formatting;3originals closed. DetailedAppendix selected where main differs, conflicts explicit in note.'
INITIAL='Source extracted plasma cfDNA, selected20ng-scale Appendix route; paired tissue MspI-RRS/WGS preparations remain separate'
CELL_INPUT=False
READOUT='Source NovaSeqXPlus150bp paired-end; source8bp i5/i7 UDI identifiers plus hairpin randomD4/biotindT/D4 positions. Main310 records4baseUMItag whenavailable; correspondence of parsedtag to full8randompositions unassembled. Chemical U afterdeamination is distinct from DNA T reported afterPCR; no perfect conversion, repair or duplex-error-rate prediction.'
DEFERRED='Source main/Appendix differences: ligase2000UvsAppendix2ul2000Uul(4000U),ligation1hvs30min,formamideelution10ulvs20ul,Mg5uMvsAppendix3.5ul100mMin70ul(5mM); selectedAppendix values explicit. MainFig1caption ordersmethylationafterligation butMethods/Appendix concurrentlydigest/methylatebeforeligation. HistoricalEM-seq kit identities and full end/ligation/copy/read geometry unassembled, no absenceclaim.'
OLIGOS={}
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from chemdraw import Scene,Segment,feature,revcomp
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry

def scene(label):
    item=OLIGOS.get(label)
    if item is None:return None
    sc=Scene();sc.strand('oligo',item['segments'](),mod5=item.get('mod5',''),mod3=item.get('mod3',''));return sc

def chemistry(label):return OLIGOS.get(label,{}).get('chemistry')
def unknown_scene():
    sc=Scene();sc.strand('unassembled',[Segment('unassembled geometry','[source geometry unassembled]',placeholder=True)]);return sc

from illumina import TRUSEQ_READ1,NEBNEXT_ARM_READ2,P5,P7
HAIRPIN_LEFT='CGGCCGCAGTC';HAIRPIN_RIGHT='GACTGCGGC'
Y_TOP=TRUSEQ_READ1+'G'
Y_BOTTOM='CGCA'+NEBNEXT_ARM_READ2
HP_RAW="5'-p-CGGCCGCAGTCDDDD/Int Biotin dT/DDDDGACTGCGGC-3'"
OLIGOS={
 'Source Appendix biotin hairpin':{'segments':lambda:[Segment('hairpin left',HAIRPIN_LEFT),Segment('left random UMI','D'*4,placeholder=True,feature=feature('masd_umi','umi','random'),note='D=A/G/T, not an unrestricted N window'),Segment('internal biotin-dT','T',note='Source internal biotinylated dT'),Segment('right random UMI','D'*4,placeholder=True,feature=feature('masd_umi','umi','random'),note='D=A/G/T; full parsedtag assignment unassembled'),Segment('hairpin right',HAIRPIN_RIGHT)],'mod5':'P','chemistry':{'polymer':'DNA','five_prime':'P','internal_biotin_dT':len(HAIRPIN_LEFT)+4,'random_base_classes':'D=A/G/T;8positions'}},
 'Source methylated Y top':{'segments':lambda:[Segment('source Ytop',Y_TOP)],'chemistry':{'polymer':'DNA','phosphorothioate_after':(len(Y_TOP)-2,),'modified_bases':tuple((i,'5mC')for i,b in enumerate(Y_TOP)if b=='C'),'scope':'Appendix allC instruction; TableEV1blueNotI marking doesnot establishunmodifiedC'}},
 'Source methylated Y bottom':{'segments':lambda:[Segment('source Ybottom',Y_BOTTOM)],'mod5':'P','chemistry':{'polymer':'DNA','five_prime':'P','modified_bases':tuple((i,'5mC')for i,b in enumerate(Y_BOTTOM)if b=='C'),'scope':'Appendix allC instruction; firstCGCblue inTableEV1, not an assertionunmodified'}},
 'Source Appendix UDI forward':{'segments':lambda:[Segment('P5',P5),Segment('i5','N'*8,placeholder=True,feature=feature('masd_i5','sample_index','unknown')),Segment('source forwardarm',TRUSEQ_READ1)],'chemistry':{'polymer':'DNA'}},
 'Source Appendix UDI reverse':{'segments':lambda:[Segment('P7',P7),Segment('i7','N'*8,placeholder=True,feature=feature('masd_i7','sample_index','unknown')),Segment('source reversearm','GT'+revcomp(NEBNEXT_ARM_READ2))],'chemistry':{'polymer':'DNA'}}}

def hairpin_fold():
    # Source-declared hairpin, with its geometry computed from exact fixed arms.
    seq=''.join(seg.top for seg in scene('Source Appendix biotin hairpin').strands['oligo'].segs)
    n=len(HAIRPIN_RIGHT);start=len(HAIRPIN_LEFT)-n;second=len(seq)-n
    if seq[start:start+n]!=revcomp(seq[second:]):raise ValueError('Source hairpin stem cannot pair')
    overhang=seq[:start];loop=seq[start+n:second]
    if overhang!='CG' or loop!='DDDDTDDDD':raise ValueError('Source hairpin interface or loop changed')
    return {'sequence':seq,'first_stem_start':start,'second_stem_start':second,'stem_length':n,'overhang_5prime':overhang,'loop':loop,'inferred':True}

def hairpin_end_scene(*,inferred=True):
    f=hairpin_fold();a=f['first_stem_start'];n=f['stem_length'];b=f['second_stem_start'];seq=f['sequence']
    sc=Scene();sc.strand('hairpin_top',[Segment('hairpin5CG',f['overhang_5prime'],inferred=inferred),Segment('hairpin stem5',seq[a:a+n],inferred=inferred)],mod5='P')
    sc.anneal('hairpin_bottom',[Segment('hairpin stem3',seq[b:b+n],inferred=inferred)],to='hairpin_top',pair=('hairpin stem3','hairpin stem5'))
    sc.overhang('hairpin_top','hairpin5CG',polarity='5-prime',owner='adapter',text='Computed 5prime CG cohesive interface')
    sc.note('hairpin_top','INFERRED local fold of one source oligo: these arms are covalently joined by D4-biotin-dT-D4 loop; full ordering is retained in inventory. No thermodynamic occupancy prediction.')
    return sc

def y_end_scene():
    # The two source Y strands have a terminal stem next to bottom5primeCG.
    overhang=Y_BOTTOM[:2]
    candidates=[n for n in range(1,min(len(Y_TOP),len(Y_BOTTOM)-2)+1) if Y_TOP[-n:]==revcomp(Y_BOTTOM[2:2+n])]
    if not candidates:raise ValueError('Source Y strands have no paired terminal stem')
    n=max(candidates)
    if overhang!=hairpin_fold()['overhang_5prime']:raise ValueError('Source adapters use different cohesive interfaces')
    sc=Scene();sc.strand('y_top',[Segment('Yfork top',Y_TOP[:-n]),Segment('Ystem top',Y_TOP[-n:])])
    sc.anneal('y_bottom',[Segment('Y5CG',overhang),Segment('Ystem bottom',Y_BOTTOM[2:2+n]),Segment('Yfork bottom',Y_BOTTOM[2+n:])],to='y_top',pair=('Ystem bottom','Ystem top'),unpaired=('Yfork bottom',),mod5='P')
    sc.overhang('y_bottom','Y5CG',polarity='5-prime',owner='adapter',text='Computed 5prime CG cohesive interface')
    sc.note('y_top','Local Y terminal stem derived from both whole source orderings; arms remain explicitly unpaired. Source methylation/PS modifications remain in oligo inventory.')
    return sc

def local_ligation_geometry():
    hp=hairpin_end_scene();yy=y_end_scene();f=hairpin_fold()
    if revcomp(f['overhang_5prime'])!='CG':raise ValueError('Hairpin end incompatible with source MspI5primeCG')
    sc=Scene();sc.strand('mspi_top',[Segment('MspI5CG','CG'),Segment('MspI flank','[source cut-fragment flank]',placeholder=True)])
    sc.note('mspi_top','Partial source MspI end only; full flank/complement and completed cfDNA junction paths unassembled.')
    sc.overhang('mspi_top','MspI5CG',polarity='5-prime',owner='substrate',text='Source C|CGG MspI interface; actual cfDNA fragment sequence unassigned')
    return AdapterLigationGeometry.known('cohesive',hp,yy,sc,statement='Computed source hairpin9bpstem/5primeCG and Y14bpterminalstem/5primeCG interfaces can pair with the source MspI5primeCG end. Complete sealed cfDNA products, ligation choices and yields remain unassembled.')

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=AUTHORITY,representation='published_instance')
    x=g.input('input',INITIAL,K.BIOLOGICAL if CELL_INPUT else K.MOLECULAR,primary=True)
    def stock(key,label):return g.input(key,label,K.OLIGO,state=scene(label),chemical_state=chemistry(label))
    def op(key,label,kind,reagent,product,oligos=(),symmetric=False):
        nonlocal x
        enzyme=kind in {S.DIGESTION,S.END_REPAIR,S.EXTENSION,S.TAILING,S.LIGATION,S.PCR,S.DEPROTEINIZATION,S.TAGMENTATION,S.TRANSCRIPTION,S.CHEMICAL_CONVERSION}
        c=g.input(key+'_stock',reagent,K.CATALYST if enzyme else K.CONDITION)
        ins=[(x,R.TEMPLATE),(c,R.CATALYST if enzyme else R.CONDITION)]
        for i,o in enumerate(oligos):ins.append((stock(key+'_oligo_'+str(i),o),R.STRUCTURAL))
        kw={}
        if kind==S.MAGIC:kw['reason']='Source-specific compound transition or molecular geometry remains unassembled in this selected model.'
        if kind==S.GENERIC:kw['changes_molecule']=False
        if kind==S.LIGATION:kw['geometry']=AdapterLigationGeometry('unavailable','unavailable',(unknown_scene(),),'Actual source ligation; full source-specific substrate/junction geometry unassembled, not a source absence claim.')
        if kind==S.PCR:kw['geometry']=PCRGeometry.unavailable(unknown_scene(),reason='Source amplification known; whole source-specific primer landing and product geometry unassembled.',symmetric_species=symmetric)
        y=g.product(key+'_product',product,K.MOLECULAR);g.react(key,label,kind,ins,(y,),**kw);x=y
    def clean(key,label):
        nonlocal x
        c=g.input(key+'_condition',label,K.CONDITION);y=g.product(key+'_kept','Physically recovered source preparation; exact recovery/carryover unassigned',K.MOLECULAR);w=g.product(key+'_waste','Physically separated source workup fractions',K.WASTE)
        g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(c,R.CONDITION)),(y,w));g.discard(w,'Source physical workup, not perfect molecular exclusion.');x=y
    def split(key,label,rest_label,disposition=Disposition.UNRESOLVED):
        nonlocal x
        y=g.product(key+'_used',label,K.MOLECULAR);r=g.product(key+'_rest',rest_label,K.QC if disposition==Disposition.QC else K.MOLECULAR);g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(y,r));g.mark_terminal(Terminal(r.id,disposition,rest_label));x=y
    def qc(key,label):split(key,'Production preparation retained after separate analytical aliquot',label,Disposition.QC)
    m=g.input('msss','Source1ulM.SssI20Uul',K.CATALYST);d=g.input('mspi','Source1ulMspI100Uul',K.CATALYST);sam=g.input('sam','Source1ul8mMSAM methyl donor',K.MOLECULAR);buf=g.input('digest_buffer','Source42ulcfDNA/5ul10xrCutSmart/50ultotal37C1h',K.CONDITION)
    cut=g.product('methyl_digest_product','Source MspI-cut CpG-methylated cfDNA classes with incompletely processed molecules; no quantitative protection/completion assumed',K.MOLECULAR)
    g.react('methyl_digest','Concurrent source M.SssI CpGmethylation and MspI restriction; no invented sequential intermediate',S.CHEMICAL_CONVERSION,((x,R.TEMPLATE),(m,R.CATALYST),(d,R.CATALYST),(sam,R.STRUCTURAL),(buf,R.CONDITION)),(cut,));x=cut
    clean('methyl_clean','Source150ulAMPure per50ulreaction;two200ul70percentEtOHwash;elute41ul')
    split('ligation_aliquot','Source40.5ulof41uleluate forligation','Untransferred0.5ul eluate; fate unassigned')
    # Ready adapter stocks are inventories; full folded/annealed junction geometry is unassembled.
    h=stock('hairpin','Source Appendix biotin hairpin');ya=stock('y_top','Source methylated Y top');yb=stock('y_bottom','Source methylated Y bottom')
    hc=g.input('hairpin_condition','Source annealed hairpin stock8ngul,1ulselected; fullstock folding geometryunassembled',K.CONDITION);yc=g.input('y_condition','Source annealed Yadapter4ngul,1ulselected; exact duplex/stoichiometryunassembled',K.CONDITION)
    hairpin=g.product('hairpin_ready','Source hairpin stock; computed local9bpstem/5primeCG interface, same stock rather than newaddition',K.ADAPTER,state=hairpin_end_scene(),chemical_state=hairpin_fold());yready=g.product('y_ready','Source preannealed Yadapter with computed terminalstem/5primeCG; same stock rather than newfree-strand addition',K.ADAPTER,state=y_end_scene())
    g.react('hairpin_prepare','Source intended hairpin; INFERRED local fold from exact9bpstem, no thermodynamic occupancy claim',S.FOLDING,((h,R.TEMPLATE),(hc,R.CONDITION)),(hairpin,))
    g.react('y_prepare','Source Yadapter stock annealing; exact historical annealing programme unassigned',S.HYBRIDIZATION,((ya,R.TEMPLATE),(yb,R.STRUCTURAL),(yc,R.CONDITION)),(yready,))
    e=g.input('ligase','Source Appendix2ulT4ligase2000Uul (4000U); main298says2000U',K.CATALYST);cond=g.input('ligation_condition','Source50ulwith5ul10xrCutSmart/0.5ul100mMATP,20C30min65C10min; main1hdiffers',K.CONDITION)
    lig=g.product('ligate_product','Source hairpin/Yadapter-ligated original-Watson/Crick cfDNA classes, adapterdimers and unligated molecules; no idealproduct yield',K.MOLECULAR)
    g.react('ligate','Source Appendixligation30min with both prepared adapters; main/Appendix differences retained',S.LIGATION,((x,R.TEMPLATE),(hairpin,R.STRUCTURAL),(yready,R.STRUCTURAL),(e,R.CATALYST),(cond,R.CONDITION)),(lig,),geometry=local_ligation_geometry());x=lig
    op('noti','Source1ulNotI20U added to50ulligation;37C15min',S.DIGESTION,'Source NotIR3189','Source cut/uncut adapter-dimer and intendedproduct classes; digestion is not physical dimer disappearance')
    clean('noti_clean','Source150ulAMPure perpostNotImix;two200ulEtOHwashes;elute40ul. Literal150ulretained ratherthan recalculated3xof51ul')
    beads=g.input('streptavidin','Source~10ulMyOneStreptavidinC1,prewashed3times,resuspended40ul2xwashbuffer',K.SELECTION)
    bound=g.product('captured','Source bead-associated biotin-bearing DNA/reagents with possible nonspecific/unligated companions',K.MOLECULAR);unbound=g.product('capture_unbound','Physical noncaptured supernatant',K.WASTE)
    cc=g.input('capture_condition','Source40ulDNA+40ul2xbufferRT15min;magnet2–3min',K.CONDITION)
    g.react('capture','Source biotin-DNA capture',S.AFFINITY_CAPTURE,((x,R.TEMPLATE),(beads,R.SELECTION),(cc,R.CONDITION)),(bound,unbound));g.discard(unbound,'Appendix explicitly removes/discards capture supernatant; specificity not idealized.');x=bound
    clean('capture_wash','Source three200ul1xwashbuffer magnetic washes; retained beadassociatedpreparation')
    elcond=g.input('formamide','Source Appendix20ul95percentformamide/10mMEDTApH8.2,65C5min;main10uldiffers',K.CONDITION)
    released=g.product('formamide_eluate','Source transferred DNA-containing supernatant after noncovalent biotin/streptavidin elution; no covalent biotin-cleavage claim',K.MOLECULAR);spent=g.product('spent_beads','Source co-separated beads/residualDNA; laterfateunassigned',K.MOLECULAR)
    g.react('release','Source magnetic separation and supernatanttransfer afterformamide incubation',S.CLEANUP,((x,R.TEMPLATE),(elcond,R.CONDITION)),(released,spent));g.mark_terminal(Terminal(spent.id,Disposition.UNRESOLVED,'No perfect DNArelease asserted.'));x=released
    clean('release_clean','Source60ulAMPure/elute28ulwater')
    split('oxidation_aliquot','Source27ulof28ulDNA into45uloxidationmix','Untransferred1ul eluate; fate unassigned')
    op('protect','Source oxidation/protection45ulplus5ulfresh1:1000Feactive reagent;37C1h',S.CHEMICAL_CONVERSION,'Source EM-seqE7125 TET2/βGT protection chemistry; buffer5ul/S1ul/reagent1ul/enzyme1 10ul/enzyme2 1ul; exact kit mappingunreviewed','Source protectedoxidized5mCpG/5hmCpG derivatives and residualunconvertedDNA; no preserved endogenous methylation measurementclaim')
    op('protect_stop','Source1ulStopReagent37C30min',S.QUENCHING,'Source EM-seq StopReagent','Source stoppedprotectionmixture')
    clean('protect_clean','Source150ulAMPure/elute54ul')
    split('deam_aliquot','Source53.55ulof54uleluate to70uldeaminationmix','Untransferred0.45ul eluate; fate unassigned')
    ap=g.input('apobec','Source1ulDEAenzymemixAPOBEC3A',K.CATALYST);uv=g.input('uvrd','Source2.5ulUvrD20ugml (50ng)',K.CATALYST);dc=g.input('deam_condition','Source Appendix70ul:7ulDEAReactionBuffer(labeled1x),0.7ulBSA,1.75ul100mMATP,3.5ul100mMMgCl2;37C90min; no separate initialdenaturation',K.CONDITION)
    deam=g.product('deam_product','Source helicase-separated C-to-chemicalU strands still covalently hairpinlinked where ligated; protected/residualCs and incompleteproducts remain',K.MOLECULAR)
    g.react('deam','Concurrent UvrD strandseparation/APOBECdeamination;AppendixMg5mM vsmain5uM explicit',S.CHEMICAL_CONVERSION,((x,R.TEMPLATE),(ap,R.CATALYST),(uv,R.CATALYST),(dc,R.CONDITION)),(deam,));x=deam
    clean('deam_clean','Source210ulAMPure/elute23.5ul')
    split('pcr_aliquot','Source22.5ulof23.5uleluate into50ulPCR','Untransferred1ul eluate; fate unassigned')
    op('pcr','Source25ulKAPAUracilPlus/2.5ulUDIprimermix;98C1min14cycles98C20s65C30s72C30s;72C1min4Chold',S.PCR,'Source KAPAHiFiHotStartUracilPlusKK2802','Source amplified original-duplex derivatives; U readoutT distinguished from chemistry; no predictedperfect duplexconsensus',('Source Appendix UDI forward','Source Appendix UDI reverse'))
    clean('pcr_clean','Source60ulAMPure1.2x/elute25ul')
    qc('final_qc','Source Qubit4dsDNAHS and fragmentanalyzer aliquots;~300–700bpnominalprofile')
    library=g.product('library','Source sequencing preparation; full source-specific library geometry unassembled',K.LIBRARY)
    g.react('loading','Physical source loading preparation',S.GENERIC,((x,R.TEMPLATE),),(library,),changes_molecule=False)
    g.sequence(library,SequencingPlan(limitation=READOUT))
    g.add_source_item(SourceItem('pending',DEFERRED,'Selected source scope only; incomplete linked-source/whole molecular review, no accepted absence claim.'))
    return g

# Whole source ordering is independently tested without collapsing D(A/G/T) into N.
SOURCE_OLIGO_EXCEPTIONS={'hairpin':'Appendixp15line208 and originalTableEV1B6 give identical complete phosphorylated D4/IntBiotindT/D4 hairpin; model retains D baseclasses and all fixed bases/modifications. Textguard collapses placeholderD toN and cannot match its internal-biotindT shorthand. Local whole rawordering normalization and position fixtures remain mandatory.'}
def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph()
    c('Fold candidate geometry inferred by default',all(seg.inferred for strand in hairpin_end_scene().strands.values()for seg in strand.segs),True)
    c('Known local cohesive interface retained',g.steps['ligate'].geometry.kind,'cohesive')
    c('Whole hairpin D classes preserved in computed fold',hairpin_fold()['sequence'],''.join(x.top for x in scene('Source Appendix biotin hairpin').strands['oligo'].segs))
    c('Concurrent digest/methylation single sourceaddition',[(i.material,i.role)for i in g.steps['methyl_digest'].inputs if i.role==R.CATALYST],[('msss',R.CATALYST),('mspi',R.CATALYST)])
    c('No endrepair imposed on MASD branch',any(st.kind==S.END_REPAIR for st in g.steps.values()),False)
    c('Ligation uses produced adapterstocks',[i.material for i in g.steps['ligate'].inputs if i.role==R.STRUCTURAL],['hairpin_ready','y_ready'])
    c('NotI before physical cleanup',g.steps['noti_clean'].inputs[0].material,'noti_product')
    c('Capture nonbound fraction explicitly discarded',g.terminals['capture_unbound'].disposition,Disposition.DISCARD)
    c('Formamide-separated beads not silently consumed',g.terminals['spent_beads'].disposition,Disposition.UNRESOLVED)
    c('Protection before UvrD/APOBEC',g.steps['deam_aliquot'].inputs[0].material,'protect_clean_kept')
    c('No standalone thermaldenaturation step',any(st.kind==S.THERMAL_TREATMENT for st in g.steps.values()),False)
    c('Deamination is simultaneous enzymeaddition',[i.material for i in g.steps['deam'].inputs if i.role==R.CATALYST],['apobec','uvrd'])
    c('Actual PCR aliquot remainder retained',g.terminals['pcr_aliquot_rest'].disposition,Disposition.UNRESOLVED)
    c('Source full UDIforward pattern',''.join(x.top for x in scene('Source Appendix UDI forward').strands['oligo'].segs),'AATGATACGGCGACCACCGAGATCTACAC'+'N'*8+'ACACTCTTTCCCTACACGACGCTCTTCCGATCT')
    c('Source full UDIreverse pattern',''.join(x.top for x in scene('Source Appendix UDI reverse').strands['oligo'].segs),'CAAGCAGAAGACGGCATACGAGAT'+'N'*8+'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC')
    c('D windows retain source baseclasses',[s.top for s in scene('Source Appendix biotin hairpin').strands['oligo'].segs if s.placeholder],['DDDD','DDDD'])
    c('Source hairpin5P retained',chemistry('Source Appendix biotin hairpin')['five_prime'],'P')
    src=Source(sources_root()/Path(__file__).resolve().parents[1].name/'s44321-026-00509-2_supp_44321_2026_509_MOESM1_ESM.pdf.txt','PublishedAppendixp15line208/doi10.1038/s44321-026-00509-2')
    if have(c,src,label='original wholemodifiedhairpin'):
        c('Raw complete modifiedordering in original',HP_RAW in src.path.read_text(),True)
        normalized=HP_RAW.removeprefix("5'-p-").removesuffix("-3'").replace('/Int Biotin dT/','T')
        c('Whole modifiedhairpin source transcription',''.join(s.top for s in scene('Source Appendix biotin hairpin').strands['oligo'].segs),normalized)
        c('Source biotinposition after D4',chemistry('Source Appendix biotin hairpin')['internal_biotin_dT'],normalized.index('DDDD')+4)
    c('Private partial graph remains incomplete',bool(g.publication_problems()),True)
    c.report()
if __name__=='__main__':selftest()
