"""Private CLAPS-seq: independent source-owned adapters and genomic preparation."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,revcomp
from illumina import NEBNEXT_ARM_READ2,TRUSEQ_READ1
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from adapter_ligation import AdapterLigationGeometry
TITLE='CLAPS-seq selected genomic OG preparation'
AUTHORITY='PublishedPMC8643665main77–129 and published SI TableS1 originalPDFp2/lines33–88; two originals. Independent local chemistry, no PADD/source-family constructor inheritance.'
AD1T=NEBNEXT_ARM_READ2+'A'
AD1B='NNNNNGACTGGTTCCAATTGAAAGTGCTCTTCCGATCT'
AD2T=revcomp(TRUSEQ_READ1)
AD2B=TRUSEQ_READ1+'N'*5
EXTENSION='T'+revcomp(NEBNEXT_ARM_READ2)+'T'
SOURCE_ORDERINGS={'ad1t':'Phos-'+AD1T+'SpC3','ad1b':AD1B[:-1]+'*'+AD1B[-1],'ad2t':'Phos-'+AD2T+'SpC3','ad2b':AD2B+'SpC3','extension':EXTENSION}
def ad1():
    s=Scene();s.strand('top',[Segment('ad1_stem',AD1T[:14]),Segment('ad1_tail',AD1T[14:])],mod5='P',mod3='SpC3')
    s.anneal('bottom',[Segment('ad1_random5','N'*5,placeholder=True),Segment('ad1_btail',AD1B[5:23]),Segment('ad1_bstem',AD1B[23:-1]),Segment('ad1_terminal_T','T')],to='top',pair=('ad1_bstem','ad1_stem'),unpaired=('ad1_random5','ad1_btail','ad1_terminal_T'),mod3='OH; terminal PS bond');return s
def ad1_end():
    s=ad1();s.anneal('substrate',[Segment('source_insert','[source fragment]',placeholder=True),Segment('source_dA','A')],to='bottom',pair=('source_dA','ad1_terminal_T'),unpaired=('source_insert',),mod3='OH')
    s.ligation_end('bottom','ad1_terminal_T',kind='cohesive',polarity='3-prime',owner='adapter',text='Literal source3T cohesive overhang; adjacent top donor is source5P')
    s.ligation_end('substrate','source_dA',kind='cohesive',polarity='3-prime',owner='substrate',text='Source dA-tailed fragment3OH; whole dsDNA/junction unassembled');return s
def ad2():
    s=Scene();s.strand('top',[Segment('ad2_stem',AD2T)],mod5='P',mod3='SpC3')
    s.anneal('bottom',[Segment('ad2_bstem',TRUSEQ_READ1),Segment('ad2_random5','N'*5,placeholder=True)],to='top',pair=('ad2_bstem','ad2_stem'),unpaired=('ad2_random5',),mod3='SpC3');return s
def ad2_end():
    s=ad2();s.anneal('substrate',[Segment('extension_body','[extension product]',placeholder=True),Segment('extension_end5','N'*5,placeholder=True,inferred=True,note='Conditional pairing to a complementary member of sourceN5pool; actual target/matched member unassigned')],to='bottom',pair=('extension_end5','ad2_random5'),unpaired=('extension_body',),mod3='OH')
    s.ligation_end('top','ad2_stem',kind='splinted',polarity='5-prime',owner='adapter',text='SourceAD2top5P/duplex stem; N5pool splint member unassigned')
    s.ligation_end('substrate','extension_end5',kind='splinted',polarity='3-prime',owner='substrate',text='Nascent extension3OH; selected matched terminal bases conditional, whole product unassembled');return s

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=AUTHORITY,representation='published_instance');x=g.input('dna','Source8–10ug extractedHeLa genomicDNA; chosen endogenous or treated cohort unassigned; extraction/ELISA outside boundary',K.MOLECULAR,primary=True)
    def inp(k,label,kind=K.CATALYST,role=R.CATALYST,**kw):return(g.input(k,label,kind,**kw),role)
    def cond(k,label):return inp(k,label,K.CONDITION,R.CONDITION)
    def op(k,label,kind,ins=(),geometry=None,product='Source mixture with residual/unconverted substrates and copresent reagents until workup'):
        nonlocal x
        y=g.product(k+'_product',product,K.MOLECULAR);kw={}
        if geometry is not None:kw['geometry']=geometry
        if kind==S.GENERIC:kw['changes_molecule']=False
        if kind==S.MAGIC:kw['reason']='Source preparation positive; complete historical stock composition/whole product geometry unassembled, not publication absence.'
        outputs=(y,)
        if kind==S.AFFINITY_CAPTURE:
            r=g.product(k+'_unretained','Unretained physical affinity-capture fraction; recovery/fate unassigned',K.MOLECULAR);outputs=(y,r)
        g.react(k,label,kind,((x,R.TEMPLATE),)+tuple(ins),outputs,**kw);x=y
        if kind==S.AFFINITY_CAPTURE:g.mark_terminal(Terminal(r.id,Disposition.UNRESOLVED,'Only bead-bound fraction follows source washes.'))
    def clean(k,label):
        nonlocal x
        y=g.product(k+'_kept','Source physically retained fraction; recovery/carryover unassigned',K.MOLECULAR);w=g.product(k+'_waste','Source separated workup fractions',K.WASTE);g.react(k,label,S.CLEANUP,((x,R.TEMPLATE),cond(k+'_condition',label)),(y,w));g.discard(w,'Physicalworkup,notperfectspecificity.');x=y
    # Source prepares full50ul annealed stocks; selected reaction doses/remainders explicit.
    prepared={}
    for n,scene,dose in [(1,ad1(),'100pmol perug selectedDNA'),(2,ad2(),'40pmol')]:
        top_scene=Scene();top_scene.strand('oligo',[Segment('source_whole',AD1T if n==1 else AD2T)],mod5='P',mod3='SpC3')
        bottom_scene=Scene();bottom_scene.strand('oligo',([Segment('random5','N'*5,placeholder=True),Segment('source_fixed',AD1B[5:])]if n==1 else[Segment('source_fixed',TRUSEQ_READ1),Segment('random5','N'*5,placeholder=True)]),mod3='OH;terminalPS'if n==1 else'SpC3')
        top=g.input('ad'+str(n)+'t','Source100uMtopoligo '+SOURCE_ORDERINGS['ad'+str(n)+'t'],K.OLIGO,state=top_scene,chemical_state={'source_ordering':SOURCE_ORDERINGS['ad'+str(n)+'t']})
        bottom=g.input('ad'+str(n)+'b','Source100uMbottomoligo '+SOURCE_ORDERINGS['ad'+str(n)+'b'],K.OLIGO,state=bottom_scene,chemical_state={'source_ordering':SOURCE_ORDERINGS['ad'+str(n)+'b'],'terminal_PS':n==1})
        stock=g.product('ad'+str(n)+'_stock','Source50ul preparedadapterstock,100uMeacholigo',K.ADAPTER,state=scene);g.react('anneal_ad'+str(n),'Sourceboil2min/slowcool25C in.5mMTris8/100mMNaCl/.05mMEDTA',S.HYBRIDIZATION,((top,R.TEMPLATE),(bottom,R.STRUCTURAL),cond('ad'+str(n)+'_buffer','Sourcehybridizationbuffer/thermalprogramme')),(stock,))
        a=g.product('ad'+str(n)+'_aliquot','Source'+dose+'selectedpreparedadapter',K.ADAPTER,state=scene);r=g.product('ad'+str(n)+'_remainder','Untransferredpreparedadapterstock; actualvolume/fateunassigned',K.MOLECULAR);g.react('ad'+str(n)+'_split','Sourceadapterstock selecteddose',S.SPLIT,((stock,R.TEMPLATE),),(a,r));g.mark_terminal(Terminal(r.id,Disposition.UNRESOLVED,'Untransferred adapter stock.'));prepared[n]=a
    op('label','Source20mMaminePEG11biotin+5mMK2IrBr6in100mMphosphatepH8,ice5min',S.CHEMICAL_CONVERSION,(inp('biotin_probe','SourceAminePEG11biotin molecular tag',K.MOLECULAR,R.STRUCTURAL),cond('label_buffer','SourceK2IrBr6oxidation/phosphatebuffer/ice')),product='SourcebiotinylatedOG classes,unreactedOG/otherDNA/background; no perfectselectivity or complete conversion')
    clean('g50','SourceMicrospinG50 prewashedicecoldwater; columnseparationquenchesreaction')
    clean('ethanol_label','Sourceethanolprecipitation/asabove77;TErecovery')
    op('shear','SourceQ800sonicatoraveraging600bp AFTERchemicalbiotinylation',S.FRAGMENTATION,(cond('sonicator','SourceQ800shearing'),))
    used=g.product('dna_used','Source up to5ug purified-fragment input selected for endprep; intervening purification workup/fractions unassembled',K.MOLECULAR);rest=g.product('dna_remainder','UntransferredlabelledDNAstock; mass/fateunassigned',K.MOLECULAR);g.react('dna_split','Sourceup to5ugaliquotafter8–10uglabeling',S.SPLIT,((x,R.TEMPLATE),),(used,rest));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Stockremainder.'));x=used
    op('endprep','SourceNEBE6040endrepair/dA-tailing; completeinternalreagentstocksunassembled',S.END_REPAIR,(inp('endprep_mix','SourcehistoricalNEBNextendprep enzyme preparation'),),product='SourcedA-tailedfragments plusincompleteclasses; wholeendprep/productgeometryunassembled')
    op('ligate1','SourceAd1100pmolperugDNA/NEBNextT4QuickLigase4Covernight',S.LIGATION,((prepared[1],R.STRUCTURAL),inp('quickligase','SourceNEBNextT4QuickLigase')),AdapterLigationGeometry.known('cohesive',ad1_end(),statement='SourceAd1terminalT/dA-fragment local pair andtop5P available; whole sealed duplexproduct unassembled.'))
    clean('ad1_clean','Source0.7xFragSelectbeads/elute15ul.1xTE')
    op('capture','Source5ulC1beads+10ugsalmonDNAinbindingbuffer4C1h',S.AFFINITY_CAPTURE,(inp('c1','SourceMyOneStreptavidinC1',K.SELECTION,R.SELECTION),inp('salmon','Source10ugsalmonspermDNAphysicalblockingstock',K.MOLECULAR,R.STRUCTURAL),cond('capture_buffer','Source.5MNaCl/5mMMgCl2/Tris8/EDTA/Tween/CA630')))
    clean('capture_wash','Source2xB&W/2xbufferII/2xbufferIII/TE/.1xTE,resuspend6ul.1xTE;retainbeads')
    primer_scene=Scene();primer_scene.strand('primer',[Segment('source_whole',EXTENSION)])
    p=g.input('extension_primer','Source30pmol '+EXTENSION,K.OLIGO,state=primer_scene,chemical_state={'source_ordering':EXTENSION})
    op('primer_anneal','SourceattachExtensionPrimertobead-boundbiotinylatedDNA',S.HYBRIDIZATION,((p,R.STRUCTURAL),))
    t=g.product('extension_template','Same-tubebeadboundDNA component',K.MOLECULAR);p=g.product('retained_extension_primer','Same-tube previously added ExtensionPrimer',K.OLIGO);g.react('extension_components','Same-tube components,no newstock orphysicalfractionation',S.GENERIC,((x,R.TEMPLATE),),(t,p),changes_molecule=False);x=t
    extension_stock=g.input('extension_reaction_stock','Source prepared KAPAHiFi extension reaction reagents; nucleotide formulation and original supply mode unassembled, no separate fresh dose inferred',K.MOLECULAR);extension_cat=g.product('extension_catalyst','Catalyst component of same source extension preparation',K.CATALYST);extension_nt=g.product('extension_nucleotides','Molecular nucleotide substrate component of same source extension preparation; amounts/convention unassigned',K.MOLECULAR);g.react('extension_reagent_components','Represent source extension preparation components, not separate invented additions',S.GENERIC,((extension_stock,R.TEMPLATE),),(extension_cat,extension_nt),changes_molecule=False)
    op('extend','SourceKAPAHiFi98C21s68C2min; complete nucleotide/buffer formulation unassembled',S.EXTENSION,((p,R.STRUCTURAL),(extension_cat,R.CATALYST),(extension_nt,R.STRUCTURAL)),product='Sourcebeadboundtemplate/extensioncomplements,stalled/readthrough/unextendedclasses; no idealstop/yieldassumption')
    op('exoi','Source1.5ulExoI37C10min; modelassay15min not imported',S.DIGESTION,(inp('exoi_enzyme','SourceExoI'),))
    op('sds_stop','Sourceequalvolume1percentSDS',S.QUENCHING,(cond('sds','SourceSDSstop'),))
    clean('extension_wash','Source2xB&W/2xbufferII/2xbufferIII;retainbeads')
    comp=g.product('complement_release','Sourcecollectedreleasedcomplementaryextensionstrands; carryoverunassigned',K.MOLECULAR);bound=g.product('bound_biotin_dna','Biotinylatedtemplate/beadresidualfraction; laterfateunassigned',K.MOLECULAR);g.react('release_split','Sourcefresh.1MNaOH/.1percentCA630RT2min,collect complementary strands ratherthanbiotinylatedOGtemplate',S.STRAND_SEPARATION,((x,R.TEMPLATE),cond('naoh_buffer','SourceNaOHelutionbuffer')),(comp,bound));g.mark_terminal(Terminal(bound.id,Disposition.UNRESOLVED,'Unreleased template/resin not sequenced by assertion.'));x=comp
    op('neutralize','Source.5Maceticacidneutralization',S.QUENCHING,(cond('acetate','Sourceaceticacid'),))
    clean('release_ethanol','Sourceethanolprecipitationandrecovery')
    op('ligate2','Source40pmolAd2/InstantStickendLigaseM0370L4Covernight',S.LIGATION,((prepared[2],R.STRUCTURAL),inp('sticky_ligase','SourceInstantStickendLigaseMasterMix')),AdapterLigationGeometry.known('splinted',ad2_end(),statement='SourceAD2fullfixedduplex/N5splintfamily; selected complementary member/targetbases conditional andwhole sealedproductunassembled.'))
    clean('ad2_clean','Source0.8xFragSelectbeadpurification')
    qc=g.product('ligation_qc','Source2percentligatedDNAforlibraryqualitycheck',K.MOLECULAR);prod=g.product('production','Source98percentligatedDNAforproductionPCR',K.MOLECULAR);g.react('qc_split','Source2percentQC/restproduction',S.SPLIT,((x,R.TEMPLATE),),(qc,prod));g.mark_terminal(Terminal(qc.id,Disposition.QC,'SourceQCnotreturned; noanalyticalcyclecountimported.'));x=prod
    op('pcr','SourceNEBMultiplexE7335LorE7600S/UltraIIQ5M0544; completePCRrecipe/cyclesunassembled',S.MAGIC,(inp('pcr_preparation','SourcehistoricalQ5/MultiplexprimerPCRpreparation; exact stockcomponents/oligosunassembled',K.MOLECULAR,R.STRUCTURAL),),product='Sourceamplifiedproductionlibrary; wholeindexedproducts/readmasksunassembled')
    load_aliquot=g.product('loading_aliquot','Selected qualified library aliquot for pooling/loading; amount unassigned',K.MOLECULAR);load_rem=g.product('loading_remainder','Untransferred qualified library stock; later fate unassigned',K.MOLECULAR);g.react('loading_split','Source qualified libraries pooled; selected aliquot and remaining stock unquantified',S.SPLIT,((x,R.TEMPLATE),),(load_aliquot,load_rem));g.mark_terminal(Terminal(load_rem.id,Disposition.UNRESOLVED,'No whole-stock consumption inferred from pooling.'));x=load_aliquot
    others=g.input('other_qualified_aliquots','Source other qualified sequencing-library aliquots at pool boundary; their stock preparations, identities and volumes unassembled',K.MOLECULAR)
    lib=g.product('library','Source actual pool of qualified sequencing-library aliquots',K.LIBRARY);g.react('pool','Source actual qualified library pooling; no equimolar or volume claim',S.POOL,((x,R.TEMPLATE),(others,R.TEMPLATE)),(lib,));g.sequence(lib,SequencingPlan(limitation='SourceHiSeqXorNovaSeq pairedreads/main129; runlengths/wholePCRprimers/readlandingunassembled. N5purpose unassigned; computationaldedup/guaninecoordinatefiltersnotphysicalpurification.'))
    g.add_source_item(SourceItem('pending','Completeextensionreactionformulation/wholeproducts/readlanding unassembled. OnlygenomicCLAPSbranch; controlOGT/ELISA/inputkitroute notmerged.','Twooriginalcappedsourcepass,noconfirmedabsence.'));return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph();c('Chemicaltagprecedes600bpshear',g.steps['shear'].inputs[0].material,'ethanol_label_kept')
    c('Sourceblockedtopdoesnotaddnucleotide',sum(len(s.top)for s in ad1().strands['top'].segs),len(AD1T))
    c('ExtensioninitialTretained',EXTENSION,'TGACTGGAGTTCAGACGTGTGCTCTTCCGATCT')
    c('Extensionconsumesretainedprimer',any(p.material=='retained_extension_primer'and p.role==R.STRUCTURAL for p in g.steps['extend'].inputs),True)
    c('Extensionnucleotidesfromsamepreparedstock',g.producer_of('extension_nucleotides').id=='extension_reagent_components'and any(p.material=='extension_nucleotides'and p.role==R.STRUCTURAL for p in g.steps['extend'].inputs),True)
    c('Sequencedrouteisreleasedcomplement',g.steps['neutralize'].inputs[0].material,'complement_release')
    c('Boundtemplateunresolved',g.terminals['bound_biotin_dna'].disposition,Disposition.UNRESOLVED)
    c('QCnotproductionPCR',g.steps['pcr'].inputs[0].material,'production')
    c('Poolloadingretainsstockremainder',g.terminals['loading_remainder'].disposition,Disposition.UNRESOLVED)
    c('Actualpoolcombinesselectedandotheraliquots',g.steps['pool'].kind==S.POOL and {p.material for p in g.steps['pool'].inputs}=={'loading_aliquot','other_qualified_aliquots'},True)
    c('Adapterstockremaindersretained',sum(k in g.terminals for k in ['ad1_remainder','ad2_remainder']),2)
    path=sources_root()/Path(__file__).resolve().parents[1].name/'nar_gkab1022_supp_gkab1022_supplemental_file.pdf.txt'
    if have(c,Source(path,'https://doi.org/10.1093/nar/gkab1022')):
        text=path.read_text()
        for name,seq in SOURCE_ORDERINGS.items():c('Originalwholemodifiedordering '+name,seq in text.splitlines(),True)
    c('Privateincompletegraph',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
