from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,revcomp,feature
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
class B:
    def __init__(self,title,authority,start,kind=K.BIOLOGICAL):
        self.g=ProtocolGraph(title,authority=authority,representation='published_instance');self.x=self.g.input('source',start,kind,primary=True)
    def stock(self,key,label,kind=K.CATALYST,role=R.CATALYST,**kw):return(self.g.input(key,label,kind,**kw),role)
    def cond(self,key,label):return self.stock(key,label,K.CONDITION,R.CONDITION)
    def op(self,key,label,kind,inputs=(),geometry=None,product='Source mixture including incomplete/unreacted/off-target classes; occupancy/recovery unassigned',reason=None):
        y=self.g.product(key+'_product',product,K.MOLECULAR);kw={}
        if kind==S.MAGIC:kw['reason']=reason or 'Positive source operation; complete occupied molecular inputs/products and historical formulation unassembled, not a missing-information conclusion.'
        if kind==S.GENERIC:kw['changes_molecule']=False
        if geometry is not None:kw['geometry']=geometry
        self.g.react(key,label,kind,((self.x,R.TEMPLATE),)+tuple(inputs),(y,),**kw);self.x=y;return y
    def clean(self,key,label,extra=()):
        y=self.g.product(key+'_kept','Actual source retained physical fraction; yield/carryover unassigned',K.MOLECULAR);w=self.g.product(key+'_waste','Actual separated source workup fraction',K.WASTE)
        self.g.react(key,label,S.CLEANUP,((self.x,R.TEMPLATE),self.cond(key+'_condition',label))+tuple(extra),(y,w));self.g.discard(w,'Actual workup fraction excluded; no perfect specificity or molecular yield claim.');self.x=y
    def split(self,key,label,selected,remainder,disposition=Disposition.UNRESOLVED):
        y=self.g.product(key+'_selected',selected,K.MOLECULAR);r=self.g.product(key+'_remainder',remainder,K.MOLECULAR)
        self.g.react(key,label,S.SPLIT,((self.x,R.TEMPLATE),),(y,r));self.g.mark_terminal(Terminal(r.id,disposition,remainder));self.x=y;return y

def unknown_scene():
    s=Scene();s.strand('unassembled',[Segment('source products unassembled','[unassembled]',placeholder=True)]);return s

AUTHORITY='Published final PerezGoronzy2024 PDF mainMethods3235–3285; publishedSI1 Methods380–587/originalp26; publishedSupplementaryData1(MOESM4) originalZIP/XMLA1–A97. Three originals; AM/preprint excluded.'
ANTIBODY_ORDERING='/5phos/TGACTTGNNNNNNNNTATTATGGTAGATCGGAAGAGCGTCGTGTACACAGAGTC/3Bio/'
WORKBOOK_FIRST='/5phos/TGACTTGNNNNNNNNAACTGTTCCAGATCGGAAGAGCGTCGTGTACACAGAGTC/3Bio/'
def antibody_scene():
    sc=Scene();sc.strand('antibody_id',[
        Segment('Odd tag-compatible sticky end','TGACTTG'),
        Segment('Antibody oligo deduplication UMI','NNNNNNNN',placeholder=True,feature=feature('chipdip.antibody_umi','umi','random')),
        Segment('Antibody-ID identity barcode','TATTATGGT',feature=feature('chipdip.antibody_identity','feature_barcode','fixed',note='Antibody identity; not cell or sample identity')),
        Segment('Source i5 primer-binding region and spacer; exact functional boundary unassembled','AGATCGGAAGAGCGTCGTGTACACAGAGTC')],mod5='P',mod3='biotin');return sc

def protocol_graph():
    b=B('ChIP-DIP selected equal-bead pool and joint DNA/antibody barcode readout',AUTHORITY,'Actual prepared FA-crosslinked sonicated chromatin lysate; upstream fixation/nuclear factory outside selected boundary');g=b.g;chromatin=b.x
    oligo=g.input('antibody_id_oligo','Actual SI1 example modified antibody-ID oligo, study antibody assignment unassembled',K.OLIGO,state=antibody_scene(),chemical_state={'source_ordering':ANTIBODY_ORDERING,'umi_meaning':'Antibody oligo deduplication only, not genomic DNA dedup','barcode_meaning':'9nt antibody identity, no sample-target assignment invented'})
    b.x=oligo;b.op('oligo_streptavidin','Source20ul10uMoligo+75ulPBS+5ul1mgmlstreptavidin,30minRT',S.BINDING,(b.stock('streptavidin','Actual5ulstreptavidin',K.MOLECULAR,R.STRUCTURAL),));oligocomplex=b.x
    b.op('oligo_dilute','909nMstreptavidinmolecule-basedstock diluted1:4 to227nM; fresh1:40 gives5.675nM',S.GENERIC,(b.cond('oligo_dilution','SourcePBS; coupling uses14ul5.675nM'),));b.split('coupling_oligo_aliquot','Actual14ul from freshly diluted5.675nM stock','Actual14ul coupling input','Actual untransferred oligo-streptavidin stock; not fresh duplicate');oligocomplex=b.x
    bead=g.input('protein_g_beads','Actual ProteinG Dynabeads molecular stock1ml, storagebuffer retained',K.MOLECULAR)
    beads=g.product('separated_beads','Actual original ProteinG beads after storagebuffer separation',K.MOLECULAR);saved=g.product('saved_storage_buffer','Actual original ProteinG storagebuffer saved for later resuspension; composition unassembled',K.MOLECULAR)
    g.react('save_storage_buffer','Physically retain original ProteinG storagebuffer before PBSt washes',S.SPLIT,((bead,R.TEMPLATE),),(beads,saved));b.x=beads
    b.clean('bead_wash','OncePBSt; original storagebuffer retained for subsequent resuspension')
    b.op('biotinylate','Source1mlbeads+20ul5mMSulfoNHSbiotin30minRT',S.CHEMICAL_CONVERSION,(b.stock('nhs_biotin','ActualNHSbiotin',K.MOLECULAR,R.STRUCTURAL),b.cond('nhs_condition','SourceRT30min')))
    b.split('quench_volume','Source magnetic remove500ul buffer before Tris replacement','Actual bead-containing retained fraction','Actual removed500ul buffer; fate unassigned')
    b.op('nhs_quench','Add500ul1MTris7.4,30minRT',S.QUENCHING,(b.cond('tris','Actual500ul1MTris7.4'),));b.clean('biotin_wash','2x1mlPBSt; wash fractions distinct from saved original storagebuffer')
    b.op('storage_resuspend','Resuspend biotinylated beads using actual saved original storagebuffer',S.GENERIC,((saved,R.STRUCTURAL),))
    b.split('selected_beads','Select10ulbeads for one antibody; other factory stock retained','Actual10ulbeads in200ulbindingbuffer','Other preparedbiotinylated beads outside selected factory')
    b.op('couple_oligo','14ul5.675nM oligo-streptavidin30minRT1200rpm',S.BINDING,((oligocomplex,R.STRUCTURAL),b.cond('coupling_buffer','.5xPBSt/5mMTris8/.5mMEDTA/1MNaCl')))
    b.clean('oligo_bead_wash','2xM2/2xPBSt; resuspend200ulPBSt')
    b.split('oligo_qc','Representative well20percent removed for terminal-tag ligation; halfPCR10cycles/Tapestation QC outside production','Actual remaining production oligo-beads; representative-well assignment unassembled','Actual source analytical QC allocation; no return to production',disposition=Disposition.QC)
    b.op('antibody_couple','2.5ugantibody overnight4C,30sshakingevery15min',S.BINDING,(b.stock('antibody','Actual2.5ug selected antibody; identity-to-exampleoligo assignment unassembled',K.MOLECULAR,R.STRUCTURAL),))
    b.clean('antibody_wash','2xPBSt')
    b.op('block_sites','Source200ulPBSt+4mMbiotin+2.5ugHumanIgGFc15minRT',S.BINDING,(b.stock('block_biotin','Actualfreebiotin',K.MOLECULAR,R.STRUCTURAL),b.stock('block_fc','ActualHumanIgGFc',K.MOLECULAR,R.STRUCTURAL)))
    b.clean('block_wash','2x200ulPBSt+2mMbiotin; remove finalsupernatant')
    selected=b.x;others=g.input('other_antibody_beads','Actual distinct oligo/antibody-conjugated bead factories; whole per-antibody preparations outside selected boundary',K.MOLECULAR)
    pool=g.product('bead_pool','Actual equal-amount antibody/ID conjugate pool;10ul per antibody before representativeQC allocations, exact postQC quotas unassembled',K.MOLECULAR)
    g.react('pool_beads','Source selected equal-bead pooling; titrated alternative not merged',S.POOL,((selected,R.TEMPLATE),(others,R.TEMPLATE)),(pool,))
    b.x=chromatin;b.op('pooled_ip','Source labeled beadpool added lysate1hRT',S.BINDING,((pool,R.STRUCTURAL),b.stock('ip_biotin_fc','ActualPBSt+10mMbiotin+1xPIC+2.5ugHumanIgGFc per10ulbeads',K.MOLECULAR,R.STRUCTURAL)))
    b.clean('ip_washes','2xIPwashI150mMNaCl/2xII500mMNaCl/2xM2; retain physical antibody/bead/oligo/chromatin complexes')
    b.op('end_repair','OnbeadNEBNextE6050L T4polymerase/T4PNK20C15min',S.MAGIC,(b.stock('repair_mix','Actual enzyme cocktail/reactionbuffer; complete NT/occupied DNA states unassembled',K.MOLECULAR,R.STRUCTURAL),b.stock('repair_blockers','Actual4mMbiotin/1ugIgGFc per10ulbeads',K.MOLECULAR,R.STRUCTURAL)))
    b.clean('repair_stop','3xvolPBSt+100mMEDTA quench/2xPBSt wash')
    b.op('da_tail','SourceNEBNextE6053L dA-tail37C15min',S.MAGIC,(b.stock('tail_mix','Actualsource dA-tail module/reactionbuffer; source parenthetical3prime→5primeexo versus explicitKlenowexo- retained unresolved',K.MOLECULAR,R.STRUCTURAL),b.stock('tail_blockers','Actual4mMbiotin/1ugIgGFc per10ulbeads',K.MOLECULAR,R.STRUCTURAL)))
    b.clean('tail_stop','3xvolPBSt+100mMEDTA quench/2xPBSt wash')
    b.op('dpm_tag','Source DPM split-and-pool ligation attaches common sticky end to DNA',S.MAGIC,(b.stock('dpm_stock','ActualDPM tag/ligase preparation; full modified strands/geometries unassembled',K.MOLECULAR,R.STRUCTURAL),b.stock('dpm_blockers','Actual2mMbiotin/5.4uMProteinG',K.MOLECULAR,R.STRUCTURAL)),reason='Positive bead-associated DPM ligation; no covalent antibody-ID–genomic DNA fusion inferred. Whole tags/geometry/catalyst stock unassembled.')
    b.op('split_pool','Source >=6 Odd/Even/Terminal split-pool ligation rounds,5minRT each',S.MAGIC,(b.stock('tag_plates','Actual round-specific Odd/Even/Terminal tags and ligase stocks; whole orderings unassembled',K.MOLECULAR,R.STRUCTURAL),b.stock('tag_blockers','Actual2mMbiotin/5.4uMProteinG perround',K.MOLECULAR,R.STRUCTURAL)),reason='Source shared combinatorial bead tags join DNA and antibody-ID reads computationally, not one covalent molecule. Most experiments6rounds24–36tags; per-run plate assignments/wholeproducts unassembled.')
    b.split('library_aliquot','Resuspend1mlPKbuffer; source selected.05–4percent bead aliquot','Actual chosen bead aliquot; unique library-index pair later','Actual remaining barcoded bead stock/other aliquots, no duplicated fullstock')
    b.op('deproteinize','Source8ulPK55C2h1200rpm;65Covernight1600rpm reverse-crosslink',S.MAGIC,(b.stock('pk','Actual8ulProteinaseK'),))
    b.clean('dna_isolate','ZymoIC/6xDNAbindingbuffer/elute21ulwater')
    b.op('pcr','21ulDNA+2uli5(12.5uM)+2uli7(12.5uM)+25ul2xQ5,9–12cycles',S.MAGIC,(b.stock('pcr_mix','ActualQ5HS2xMasterMix25ul; whole NT/catalyst composition unassembled',K.MOLECULAR,R.STRUCTURAL),b.stock('pcr_primers','Actuali5/i7 primerpair, full orderings unassembled; aliquot-specific sample indexes',K.MOLECULAR,R.STRUCTURAL)))
    b.clean('library_clean','1.2xSPRI/elute20ulwater/2percentagarose removes unused primers')
    lib=g.product('joint_library','Actual mixed antibody-ID and genomic DNA read classes carrying shared split-pool tags, NOT covalent pairwise DNA-antibody fusion',K.LIBRARY);g.react('load','Source asymmetrical paired-end acquisition',S.GENERIC,((b.x,R.TEMPLATE),),(lib,),changes_molecule=False)
    g.sequence(lib,SequencingPlan(read_lengths={'R1':100,'R2':200},limitation='Selected SI300cycle run100R1/200R2; alternative50/150 and publishedmainAVITI not merged. Exact tags/primers/landing unassembled. DNA dedup alignment position percluster; antibody-ID dedup UMI.'))
    g.add_source_item(SourceItem('tail_reagent_discrepancy','SI1 originalPDFp29 chromatin dA-tail paragraph gives containingKlenow(3prime→5primeexo) but next sentenceKlenowexo-. Model retains source dA-tail operation without selecting complete exonuclease form.','SI1 renderedp26 color assigns7ntsticky/8ntUMI/9ntantibody identity; antibody assignment to example oligo unassembled.'))
    return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    import zipfile,xml.etree.ElementTree as ET
    c=Check();g=protocol_graph();c('Antibody oligo whole source ordering', ''.join(s.top for s in antibody_scene().strands['antibody_id'].segs),ANTIBODY_ORDERING.removeprefix('/5phos/').removesuffix('/3Bio/'));c('Reuses produced original storagebuffer',g.steps['storage_resuspend'].inputs[1].material,'saved_storage_buffer');c('Initial wash consumes separated beads',g.steps['bead_wash'].inputs[0].material,'separated_beads');c('Actual pooledIP precedes DNA repair',g.steps['end_repair'].inputs[0].material,'ip_washes_kept');c('DPM precedes shared bead tagging',g.steps['split_pool'].inputs[0].material,'dpm_tag_product');c('PK uses actual selected aliquot',g.steps['deproteinize'].inputs[0].material,'library_aliquot_selected');c('No covalent DNA antibody fusion constructor',not any('fusion' in x for x in g.steps),True)
    root=sources_root()/Path(__file__).resolve().parents[1].name;p=root/'s41588-024-02000-5_41588_2024_2000_MOESM4_ESM.xlsx'
    if have(c,Source(p,'https://doi.org/10.1038/s41588-024-02000-5')):
        z=zipfile.ZipFile(p);ns={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'};strings=[''.join(x.itertext()) for x in ET.fromstring(z.read('xl/sharedStrings.xml')).findall('s:si',ns)];c('Original workbook first whole modified oligo',strings[0],WORKBOOK_FIRST)
    p=root/'s41588-024-02000-5_41588_2024_2000_MOESM1_ESM.pdf.txt'
    if have(c,Source(p,'https://doi.org/10.1038/s41588-024-02000-5')):
        t=p.read_text().replace('\n','');c('Original SI example whole modified oligo',ANTIBODY_ORDERING in t.replace(' ',''),True);c('Original genomic dedup positional','alignment position for DNA reads or UMI' in p.read_text(),True)
    c('Privatepartial unready',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
