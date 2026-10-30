"""Private aa-tRNA-seq preparation; qualitative link chemistry, not an atom map."""
from dataclasses import dataclass
from pathlib import Path
import re,sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[1]/'lib'))
from chemdraw import Scene,Segment,feature
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,Terminal,Disposition,SequencingPlan,SourceItem

AUTHORITY='Published PMC12368100 Methods220–246/activation182 and actualPDFp11; publishedSI1 TableS2 originalp2. Closed two-original pass; Benchling/vendor originals unread.'
UNIVERSAL='CCTAAGAGCAAGAAGAAGCrCrUrGrGrN'
CHARGED='/5Phos/rGrGrCrUrUrCrUrUrCrUrUrGrCrUrCrUrUrCrCrArArCrCrUrUrGrCrCrUrUAAAAAAAAAA'
UNCHARGED='/5Phos/rGrGrCrUrUrCrUrUrCrUrUrGrCrUrCrUrUrArUrGrGrArArGrGrUrArGrGrCAAAAAAAAAA'
RAW={'universal':UNIVERSAL,'charged':CHARGED,'uncharged':UNCHARGED}
SOURCE_OLIGO_EXCEPTIONS={k+'_stock':'Published SI TableS2 biological adapter is line-wrapped with intervening table labels and a split ribonucleotide r-prefix. Original p2 rendered and inspected; local wholemodified ordering fixture reconstructs ONLY sequence lines within the named biological adapter row and checks the entire ordering. Narrow text extraction exception, no chemistry or source-absence exemption.' for k in ('charged','uncharged')}

@dataclass(frozen=True)
class Bridge:
    aminoacyl_ester: bool=True
    adapter_five: str='phosphate'
    amino_n_to_adapter_p: bool=False
    def __post_init__(self):
        if self.adapter_five not in ('phosphate','imidazolide','phosphoramidate'):
            raise ValueError('Unrepresented source adapter activation state')
        if self.amino_n_to_adapter_p != (self.adapter_five=='phosphoramidate'):
            raise ValueError('Phosphoramidate descriptor and N–P bond must agree')
        if self.amino_n_to_adapter_p and not self.aminoacyl_ester:
            raise ValueError('Selected source bridge requires a retained aminoacyl ester')
    def activate(self):
        if self.adapter_five!='phosphate':raise ValueError('Activation requires source5P precursor')
        return Bridge(self.aminoacyl_ester,'imidazolide')
    def couple(self):
        if not self.aminoacyl_ester or self.adapter_five!='imidazolide':
            raise ValueError('Selected chemical relation needs charged tRNA and activated adapter')
        return Bridge(True,'phosphoramidate',True)

def charged_product_descriptor():
    b=Bridge().activate().couple()
    return {'native_aminoacyl_ester_retained':b.aminoacyl_ester,
            'amino_n_to_adapter_p':b.amino_n_to_adapter_p,
            'amino_acid_is_nucleotide':False,
            'scope':'Selected qualitative charged capture class only; aminoacid identity,2prime/3primeesterregister, whole RNA/atom map and reaction populations unassembled'}

def rna_bases(raw):return ''.join(re.findall(r'r([ACGUN])',raw)).replace('U','T')
def adapter_scene(which):
    raw=RAW[which];s=Scene()
    if which=='universal':
        dna=raw.split('r')[0]
        parts=[Segment('source DNA universal arm',dna,note='Source DNA sugar'),
               Segment('source RNA terminal arm',rna_bases(raw),placeholder=True,note='RNA rC/rU/rG/rG/rN; T is displayed DNA alphabet for source U, terminal rN is not a UMI')]
        mod=''
    else:
        parts=[Segment('source '+which+' RNA arm',rna_bases(raw),
            note='RNA sugar; T displays source U. Whole adapter-arm marker, exact differential-code subwindow not assigned',
            feature=feature('aatrna.capture_'+which,'inline_barcode','fixed',note='Known adapter class distinguishes charged/uncharged capture; no cell/sample/UMI meaning')),
            Segment('source DNA A10 tail','A'*10,note='Deoxyadenosine tail, not RNA polyadenylation')]
        mod='P'
    s.strand('adapter',parts,mod5=mod);return s

def adapter_descriptor(which):
    if which=='universal':blocks=[{'sugar':'DNA','bases':UNIVERSAL.split('r')[0]},{'sugar':'RNA','bases':''.join(re.findall(r'r([ACGUN])',UNIVERSAL))}]
    else:blocks=[{'sugar':'RNA','bases':''.join(re.findall(r'r([ACGUN])',RAW[which]))},{'sugar':'DNA','bases':'A'*10}]
    return {'source_modified_ordering':RAW[which],'sugar_blocks':blocks,'five_prime':'unspecified' if which=='universal' else 'P','scope':'Whole source oligo only, not occupied hybrid/product geometry'}

def protocol_graph():
    g=ProtocolGraph('aa-tRNA-seq selected biological yeast RNA branch',authority=AUTHORITY,representation='published_instance')
    x=g.input('source','Actual frozen yeast pellet; charged/deacylated RNA populations, aminoacid identities and occupied sequences unassembled',K.BIOLOGICAL,primary=True)
    def reagent(key,label,kind=K.CONDITION,role=R.CONDITION,**kw):return g.input(key,label,kind,**kw),role
    def op(key,label,kind=S.MAGIC,inputs=(),product='Source mixture includes unreacted/off-product classes; whole occupied molecules and recovery unassigned',chemical_state=None,**kw):
        nonlocal x
        y=g.product(key+'_product',product,K.MOLECULAR,chemical_state=chemical_state)
        if kind==S.MAGIC:kw.setdefault('reason','Source operation positive; full substrates/products and detailed molecular outcomes unassembled.')
        if kind==S.GENERIC:kw['changes_molecule']=False
        g.react(key,label,kind,((x,R.TEMPLATE),)+tuple(inputs),(y,),**kw);x=y;return y
    def clean(key,label,extra=()):
        nonlocal x
        y=g.product(key+'_kept','Actual retained source fraction; imperfect recovery/carryover',K.MOLECULAR);w=g.product(key+'_waste','Actual excluded workup fractions',K.WASTE)
        g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),reagent(key+'_conditions',label))+tuple(extra),(y,w));g.discard(w,'Physical source workup; not ideal molecular specificity');x=y
    def aliquot(key,label,qc=False):
        nonlocal x
        y=g.product(key+'_selected',label,K.MOLECULAR);r=g.product(key+'_remainder','Untransferred stock; amount/fate unassigned',K.MOLECULAR)
        g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(y,r))
        if qc:g.mark_terminal(Terminal(y.id,Disposition.QC,label+'; assay internals unassembled'));x=r
        else:g.mark_terminal(Terminal(r.id,Disposition.UNRESOLVED,'Untransferred stock, not implicit waste'));x=y
    # Parallel manufacture starts with source synthesized oligos; no synthetic-tRNA factory.
    u=g.input('universal_stock','Source synthesized biological universal5prime splint stock',K.OLIGO,state=adapter_scene('universal'),chemical_state=adapter_descriptor('universal'))
    c=g.input('charged_stock','Source synthesized biological charged3prime adapter5P stock',K.OLIGO,state=adapter_scene('charged'),chemical_state=adapter_descriptor('charged'))
    n=g.input('uncharged_stock','Source synthesized biological uncharged3prime adapter5P stock',K.OLIGO,state=adapter_scene('uncharged'),chemical_state=adapter_descriptor('uncharged'))
    prepared={}
    for key,m in [('universal',u),('charged',c),('uncharged',n)]:
        p=g.product(key+'_purified','Source gel-purified100uM '+key+' adapter stock; whole purification recovery unassigned',K.OLIGO)
        w=g.product(key+'_purification_waste','Actual rejected adapter gel/workup fractions',K.WASTE)
        g.react(key+'_purify','Source2mMstock→6percentTBU/10nmolperlane/UVshadow→crushsoak4Covernight/EtOH→100uM',S.CLEANUP,
            ((m,R.TEMPLATE),reagent(key+'_purification_conditions','Source gel/crushsoak/ethanol conditions')),(p,w));g.discard(w,'Physical adapter purification');prepared[key]=p
    act=g.product('activated_adapter_mixture','Charged adapter5prime imidazolide and copresent EDC/imidazole-derived/unreacted species; full chemistry unassembled',K.MOLECULAR,
        chemical_state={'selected_adapter_five':'imidazolide','backbone_unchanged':True,'scope':'Qualitative source activation class, no complete atom balance'})
    g.react('activate_adapter','Sourcechargedadapter+500foldEDC/imidazole25C2h',S.CHEMICAL_CONVERSION,
        ((prepared['charged'],R.TEMPLATE),reagent('edc','Source EDC activating reactant',K.MOLECULAR,R.STRUCTURAL),reagent('imidazole','Source imidazole activating reactant',K.MOLECULAR,R.STRUCTURAL),reagent('activation_conditions','Source25C2h')),(act,))
    saved=x;x=act
    clean('activation_workup','Source30ulperchlorate-saturatedcoldacetone+1mlacetone/dryice20min/4Cspin/twoacetone-etherwashes/dry→200uM10mMimidazole7')
    aliquot('activated_transfer','Actualfivefold-excess activated3prime adapter contribution');activated=x;x=saved
    # Universal splint contributes separately to chemical and later enzymatic capture.
    split=[]
    for suffix,label in [('chemical','Actualfivefold-excess universal5prime splint forchemicalcapture'),('enzymatic','Actualfresh9pmol universal5prime splint forlaterRNL2'),('remainder','Untransferred prepared universal stock, fateunassigned')]:split.append(g.product('universal_'+suffix,label,K.OLIGO))
    g.react('universal_allocations','Allocate actual two source splint contributions from preparedstock',S.SPLIT,((prepared['universal'],R.TEMPLATE),),tuple(split));g.mark_terminal(Terminal(split[2].id,Disposition.UNRESOLVED,'Saved universalstock'))
    op('extract','Source400ulcoldAES+400ulacidPCI/vortex15s/ice20min/vortex5minintervals',inputs=(reagent('acid_pci','Actual400ulacidphenol/chloroform/isoamylalcohol',K.MOLECULAR,R.STRUCTURAL),reagent('aes','Actual400ul10mMNaOAc4.5/1mMEDTA/.5percentSDS')))
    clean('aqueous','Source18,000xg10min literal40C; retain aqueousphase, intendedtemperatureunresolved')
    aliquot('aqueous_transfer','Actual100ul aqueousphase selected')
    clean('ampure_supernatant','Source.4xAMPure/RT2min/magnet RETAIN SUPERNATANT',extra=(reagent('ampure04','Actual.4xAMPure beads; excluded bead fraction',K.SELECTION,R.SELECTION),))
    aliquot('nanodrop_rna','ActualNanodropRNA analyticalaliquot; volumeunassigned',qc=True)
    aliquot('r1018_input','Source100ugRNA selected forR1018 smallRNAkit')
    clean('r1018','SourceR1018 with10mMNaOAc4.5 dilution NOTwater; elute30ul10mMNaOAc4.5')
    aliquot('smallrna_qc','SourceNanodrop smallRNA analyticalaliquot',qc=True)
    aliquot('chemical_input','Source15–50pmolsmallRNA selected')
    op('chemical_capture','Sourceactivatedcharged3adapter+universal5splint/100mMMES5.5/2.5mMMg/50mMHEI6.5/25C30min',S.CHEMICAL_CONVERSION,
        ((activated,R.STRUCTURAL),(split[0],R.STRUCTURAL),reagent('hei','SourceHEI catalyst',K.CATALYST,R.CATALYST),reagent('chemical_conditions','Source MES/Mg/HEIbuffer; actualRNAmixture')),
        product='Selected aa-tRNA retains aminoacylester and gains amino-N–adapter-P phosphoramidate bridge; uncharged/native/freeadapter/otherclasses copresent',chemical_state=charged_product_descriptor())
    clean('tbu_selection','Source10percentTBU70–150nt/crushsoak4Covernight/EtOH/GlycoBlue/NaOAc4.5; mixed RNA classes retained',extra=(reagent('glycoblue','Source GlycoBlue physical coprecipitant/carrier',K.MOLECULAR,R.STRUCTURAL),))
    aliquot('gel_qc','SourceNanodrop selectedRNA analyticalaliquot',qc=True)
    aliquot('rnl2_input','Actual20pmol fromfirstligation stock')
    nselected=g.product('uncharged9','Actual9pmol uncharged3prime adapter contribution',K.OLIGO);nrest=g.product('uncharged_rest','Untransferred unchargedadapterstock',K.OLIGO)
    g.react('uncharged_allocations','Select actual9pmol frompreparedstock',S.SPLIT,((prepared['uncharged'],R.TEMPLATE),),(nselected,nrest));g.mark_terminal(Terminal(nrest.id,Disposition.UNRESOLVED,'Untransferred adapterstock'))
    op('rnl2_capture','Source20ul10percentPEG/9pmoluncharged3adapter/FRESH9pmol5splint/2ulRNL2.74mgml/25C30min',
        inputs=((nselected,R.STRUCTURAL),(split[1],R.STRUCTURAL),reagent('rnl2','Source homemadeT4RNAligase2 catalyst',K.CATALYST,R.CATALYST),reagent('rnl2_conditions','Source1xRNL2buffer/1ulRNaseinhibitor/10percentPEG8000')),
        product='Source conventional adapter capture ofdeacylatedRNA and5primeadapterjoining; chemicallycaptured aa-linkedRNA/free/unreacted/offproducts retained; full junctions unassembled')
    clean('rnl2_beads','Source1.8xtRNAbeads/ice4min/discardsupernatant/180ul80percentEtOH/dry/elute13ul')
    aliquot('rnl2_qc','Source quantified adapterligatedeluate analyticalaliquot',qc=True)
    aliquot('rta_input','Actual12.5ul from13ulRNAeluate')
    rta=g.input('rta','Actual1.5ulONT RNA004 RTA; wholeadapter/ends/geometryunassembled',K.ADAPTER)
    op('rta_ligate','Source12.5sample/1.5RTA/.5inhibitor/4buffer/1.5ul600UulT4DNAligase25C30min',inputs=((rta,R.STRUCTURAL),reagent('rta_ligase','ActualT4DNAligase catalyst',K.CATALYST,R.CATALYST),reagent('rta_conditions','Actual4ulT4buffer/.5ulRNaseinhibitor')))
    clean('rta_beads','Source1.35xtRNAbeads RT cleanup/elute26ul')
    aliquot('qubit','SourceQuant-iTQubitdsDNAHS assay aliquot, retained asprinted',qc=True)
    aliquot('tapestation','SourceTapeStationHS-DNA1000 assay aliquot, retained asprinted',qc=True)
    aliquot('rla_input','Actual50–400fmol in23ul fromRTAproduct stock')
    rla=g.input('rla','Actual6ulONT RNA004 RLA; wholeadapter/geometryunassembled',K.ADAPTER)
    op('rla_ligate','Source23sample/6RLA/8buffer/3ul600UulT4DNAligase25C30min,samesequencingday',inputs=((rla,R.STRUCTURAL),reagent('rla_ligase','Actual3ulT4DNAligase catalyst',K.CATALYST,R.CATALYST),reagent('rla_conditions','Actual8ulT4buffer')))
    clean('final_cleanup','Source1.8xAMPure/ONTWSB washes; completeRNA004vendoroperationsunassembled')
    aliquot('loading','Actualsequencerloadingallocation fromfinishedstock')
    lib=g.product('library','Source aa-tRNA directRNA library; wholeoccupiedRNA/adapter/readgeometryunassembled',K.LIBRARY)
    g.react('load','SourcePromethIONP2Solo/P2RNAflowcell',S.GENERIC,((x,R.TEMPLATE),),(lib,),changes_molecule=False)
    g.sequence(lib,SequencingPlan(limitation='DirectRNAcurrentreadout; aminoacidbridge notencodedasabase. WholeRTA/RLA/motor/readlanding/biologicaltRNA sequences and2prime/3primeesterregisterunassembled. Basecaller/chargingclassification computational, no physicaltranscriptconversion asserted.'))
    g.add_source_item(SourceItem('temperature','Publishedmainoriginalp11 literallystates18,000xg10min40C duringacidPCIextraction.','Preserve printed40C; intendedtemperatureunresolved, not silently normalizedto4C or a provenauthorerror.'))
    g.add_source_item(SourceItem('bridge_scope','Mainactivation182/Methodschemicalcapture andSIchemicalbridge: aminoacylesterretained plusamino-N–adapter-Pconnection.','Qualitative selectedclassonly; atommap/occupiedRNA/aminoacididentities/allproductsunassembled. UniversalNisstructuraldegeneracy, notUMI; internaladapterclassesnotcell/sampleidentifiers.'))
    return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph();p=sources_root()/HERE.parent.name/'s41467-025-62545-9_supp_41467_2025_62545_MOESM1_ESM.pdf.txt'
    if have(c,Source(p,'https://doi.org/10.1038/s41467-025-62545-9')):
        t=p.read_text()
        for name,raw in RAW.items():
            if name=='universal':source=''.join(t.split())
            else:
                start="Charged 3' adapter" if name=='charged' else "Uncharged 3' adapter"
                end="Uncharged 3' adapter" if name=='charged' else "Gly-GCC-3'-exon probe"
                block=t.split(start,1)[1].split(end,1)[0]
                source=''.join(line.strip() for line in block.splitlines() if '/5Phos/' in line or re.fullmatch('[ACGUNr]+',line.strip()))
            c('Wholemodified biologicaladapter '+name,raw in source,True)
    c('Nativeesterretained',charged_product_descriptor()['native_aminoacyl_ester_retained'],True)
    c('NewbridgeisNtoP',charged_product_descriptor()['amino_n_to_adapter_p'],True)
    c('Aminoacidnotnucleotide',charged_product_descriptor()['amino_acid_is_nucleotide'],False)
    for name,b in [('unactivated',Bridge()),('deacylated',Bridge(False).activate())]:
        try:b.couple();rejected=False
        except ValueError:rejected=True
        c('Reject selectedchemicalcoupling '+name,rejected,True)
    c('RNAactivationdoesnotaddbase',Bridge().activate().adapter_five,'imidazolide')
    c('UniversalterminalriboseN',adapter_scene('universal').strands['adapter'].segs[-1].top,'CTGGN')
    for name in ['charged','uncharged']:c(name+' hasDNA-A-tail',adapter_scene(name).strands['adapter'].segs[-1].top,'A'*10)
    c('Distinctcapturecodes',rna_bases(CHARGED)==rna_bases(UNCHARGED),False)
    c('Chemicalcaptureusestheactivatedallocation',g.steps['chemical_capture'].inputs[1].material,'activated_transfer_selected')
    c('Freshenzymaticsplintallocation',g.steps['rnl2_capture'].inputs[2].material,'universal_enzymatic')
    c('Enzymaticcaptureafter20pmolselection',g.steps['rnl2_capture'].inputs[0].material,'rnl2_input_selected')
    c('RTAafter12point5ultransfer',g.steps['rta_ligate'].inputs[0].material,'rta_input_selected')
    c('RLAafteractualQCstock',g.steps['rla_ligate'].inputs[0].material,'rla_input_selected')
    c('DNAAssaysareQCnotlibrary',g.terminals['qubit_selected'].disposition,Disposition.QC)
    c('No RTAwholebodyborrowed',g.materials['rta'].state,None)
    c('Privatepartial',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
