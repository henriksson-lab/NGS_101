"""Private source-qualified partial molecular preparation; no full chemistry/absence audit."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from protocol_graph import SourceCoverage
from chemdraw import Scene,Segment
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
class B:
    def __init__(self,title,authority,label):
        self.g=ProtocolGraph(title,authority=authority,representation='published_instance');self.x=self.g.input('source',label,K.BIOLOGICAL,primary=True)
    def stock(self,key,label,kind=K.MOLECULAR,role=R.STRUCTURAL,**kw):return(self.g.input(key,label,kind,**kw),role)
    def cond(self,key,label):return self.stock(key,label,K.CONDITION,R.CONDITION)
    def op(self,key,label,kind=S.MAGIC,inputs=(),product='Source preparation; whole occupied molecules, conversion and side products unassembled',**kw):
        y=self.g.product(key+'_product',product,K.MOLECULAR)
        if kind==S.MAGIC:kw.setdefault('reason','Positive source boundary; detailed occupied molecules, molecular products and conversion unassembled. Not source-absence evidence.')
        if kind==S.GENERIC:kw['changes_molecule']=False
        self.g.react(key,label,kind,((self.x,R.TEMPLATE),)+tuple(inputs),(y,),**kw);self.x=y;return y
    def clean(self,key,label):
        y=self.g.product(key+'_kept','Actual retained source workup fraction; carryover and yield unassigned',K.MOLECULAR);w=self.g.product(key+'_waste','Excluded source workup fraction',K.WASTE)
        self.g.react(key,label,S.CLEANUP,((self.x,R.TEMPLATE),self.cond(key+'_condition',label)),(y,w));self.g.discard(w,'Source workup excludes this fraction; no perfect purification claim.');self.x=y;return y
    def aliquot(self,key,label,qc=False):
        y=self.g.product(key+'_selected','Actual selected allocation: '+label,K.MOLECULAR);r=self.g.product(key+'_remainder','Untransferred stock; amount and later fate unassigned',K.MOLECULAR)
        self.g.react(key,label,S.SPLIT,((self.x,R.TEMPLATE),),(y,r))
        if qc:self.g.mark_terminal(Terminal(y.id,Disposition.QC,label));self.x=r
        else:self.g.mark_terminal(Terminal(r.id,Disposition.UNRESOLVED,'Untransferred stock retained as unresolved, no invented discard.'));self.x=y
        return self.x
    def finish(self,label,limitation):
        y=self.g.product('library',label,K.LIBRARY);self.g.react('load_library','Selected loading allocation, whole product unassembled',S.GENERIC,((self.x,R.TEMPLATE),),(y,),changes_molecule=False);self.g.sequence(y,SequencingPlan(limitation=limitation));return y

def oligo(seq,name):
    sc=Scene();sc.strand('DNA',[Segment(name,seq)]);return sc

PROBES={'18S-1': 'CTGGTTGATCCTGCCAGTAGTCATATGCTTGTC', '18S-2': 'CAAAGATTAAGCCATGCATGTCTAAGTATAAGC', '18S-3': 'CAGTGAAACTGCGAATGGCTCATTAAATCAG', '18S-4': 'CATGGTATAACTGTGGTAATTCTAGAGCTAATACATGC', '18S-5': 'TGGAAGAGATGTATTTATTAGAT', '18S-6': 'TGATGATTCATAATAACTTTTCGAATCGCATGGC', '18S-7': 'TGGCGATGGTTCATTCAAATTTCTGCCCTATCAACTT', '18S-8': 'TCCTAATTCAGGGAGGTAGTGACAATAAATAACGA', '18S-9': 'TACCTTAACGAGGAACAATTGGAGGGCAAGTC', '18S-10': 'GGGCAAGTCTGGTGCCAGCAGCCGCGGTAATTCCAGCTCC', '18S-11': 'TTCCAACGGGGCCTTT', '18S-12': 'TTGAAAAAATTAGAGTGTTCAAAGCAGGCGT', '18S-13': 'TTGCTCGAATATATTAGCATGGAATAATAGAAT', '18S-14': 'ATGATTAATAGGGACGGTCGGGGGCAT', '18S-15': 'TGCCAAGGACGTTTTCATTAATCAAGAACGA', '18S-16': 'ACTATGCCGACTAGGGATCGGGTGGTGTT', '18S-17': 'CCCACTCGGCACCTTACGAGAAATCAAAGTCTT', '18S-18': 'AGGTCCAGACACAATAAGGATTGACAGAT', '18S-19': 'AGTTGGTGGAGTGATTTGT', '18S-20': 'CGCTACACTGACGGAGCCAGCGAGT', '18S-21': 'GATTGAATGGCTTAGTGAGGCCTCAGGATC', '18S-22': 'AGCGGAGAATTTGGACAAACTTGGTCATT', '18S-23': 'ACCTGCGGAAGGATCATTA', '18S-24': 'ACCTGCGGAAGGATCATTAAAAAAA', '25S-1': 'GACCTCAAATCAGGTAGGAGTACCCGCTGAACT', '25S-2': 'GGAAAAGAAACCAACCGGGATTGCCTTAG', '25S-3': 'CTTGGAACAGGACGTCATAGAGGGTGAGAATC', '25S-4': 'AGAGTCGAGTTGTTTGGGAATGCAGCTCTAAGT', '25S-5': 'CAGTGATGGAAAGATGAAAAGAACTTTGAAAAG', '25S-6': 'GCCAGCATCAGTTTTGGTGGCAGGAT', '25S-7': 'AATACTGCCAGCTGGGACTGAGGACTG', '25S-8': 'TCTAACGTCTATGCGAGTGTTTGGGTGTAAAAC', '25S-9': 'CAAGAGGTGCACAATCGACCGATCCTGA', '25S-10': 'ATCGAACCATCTAGTAGCTGGTTCCTGCCGAAGT', '25S-11': 'ACGTAGAGTTAAGGTGCCGGAATACACGCTC', '25S-12': 'GGACGTGGGTTAGTCGATC', '25S-13': 'GGATTCTTCACGGTAACGTAACTGAATG', '25S-14': 'CAGCTTATCACCCCGGAATTGGTTTA', '25S-15': 'CCGGAGATGGGGTCTT', '25S-16': 'GCTTGTGACGGCCCGTGAAAATCCACAGGAAGG', '25S-17': 'TCCACAGGAAGGAATAGTTTTCATGCCAGGTCGTAC', '25S-18': 'TCAAAGTGAAGAAATTCAACCAAGCGCGGGTAA', '25S-19': 'GAGGGTGTAGAATAAGTGGGAGCTTCGGCGC', '25S-20': 'CACGTTCTAGCATTCAAGGTCCC', '25S-21': 'ATACCGAAGCAGAATTCGGTAAGCGTTGGATTGT', '25S-22': 'GATGGATACGAATAAGGCGTC', '25S-23': 'TGCTGGCGAATTGCAATGTC', '25S-24': 'AATCATTTGTATACGACTTAGATGTACAACGGG', '5.8S-1': 'GGGCATGCCTGTTTGAGCGTCATT', '5S-1': 'ACGCGAAACTCAGGTGCTGCAATCT'}
AUTHORITY='Published2016SupplementaryDiscussion2a286–340/Fig/abstract and originalTable3MOESM90 DNA probes. Defining Nature main cached as subscriptionpreview, fullMethods and2017detailedprotocol ordinary unreviewed. Narrow applied production route, not reconstructed full2016recipe.'
def protocol_graph():
    b=B('TCP-seq source-qualified applied footprint frontend and PDD inventory',AUTHORITY,'Source live yeast with translating ribosomes and initiating SSU; biological growth/genotype and exact cell input unassembled');g=b.g
    g.add_source_item(SourceItem('selected_scope','Discussion2a says we elected in-vivo formaldehyde fixation; introduced first sedimentation requiring at least one translating ribosome; used cDNA-based PDD without stringent footprint size selection or RNA-level rRNA removal.','Exact fixation dose, footprint nuclease/digestion/fractionation, adapter/RT/PCR/readout chemistry remain unassembled. Prior report1–3percent is not imported as an executed TCP dose. Observed19–75nt is not a physical selection window.'))
    b.op('fix','Source applied in-vivo formaldehyde crosslinking immobilizes initiating SSU on mRNA',S.CROSS_LINKING,(b.stock('formaldehyde','Actual source formaldehyde molecular contributor; dose/time/source executed conditions unassigned'),b.cond('fix_conditions','Source in-vivo formaldehyde fixation; exact executed concentration/time unassembled')),product='Actual fixed translation-complex preparation; crosslink sites and molecular stoichiometry unassembled')
    sediment=g.product('translation_sediment','Actual first-ultracentrifugation sediment: mRNAs carrying at least one translating ribosome, with associated complexes',K.MOLECULAR)
    other=g.product('initial_nonselected','Nonselected first-sedimentation fraction including completely untranslated mRNA/freeSSU classes; quantities and fate unassigned',K.MOLECULAR)
    g.react('first_sediment','Source introduced first ultracentrifugation to enrich translating mRNPs',S.CLEANUP,((b.x,R.TEMPLATE),b.cond('sediment_conditions','Actual centrifugation setup/fraction cut/conditions unassembled')),(sediment,other));g.mark_terminal(Terminal(other.id,Disposition.UNRESOLVED,'Excluded from selected preparation; physical fate not assigned.'));b.x=sediment
    supplies=g.input('footprint_preparation','Aggregate actual source footprint-generation/fractionation supplies; precise molecules/order/packaging unassembled in selected read scope',K.MOLECULAR)
    ssu=g.product('ssu_footprints','Actual selected SSU-associated RNA footprint population; exact digestion/end classes/fraction cut unassembled',K.MOLECULAR)
    rib=g.product('ribosome_footprints','Other actual ribosome-associated footprint preparation; not pooled with selectedSSU library',K.MOLECULAR)
    g.react('footprint_fractions','Opaque applied footprint-generation/fractionation interface; not a borrowed standard Ribo-seq recipe',S.MAGIC,((b.x,R.TEMPLATE),(supplies,R.STRUCTURAL)),(ssu,rib),reason='Published applied SSU/ribosome populations source-positive; full Methods unavailable in selected scope, so nuclease, processing order, dose and products unassembled.');g.mark_terminal(Terminal(rib.id,Disposition.UNRESOLVED,'Other ribosome-fraction library not expanded or implicitly pooled.'));b.x=ssu
    b.op('cdna_boundary','Source footprint-library cDNA preparation before PDD; whole recipe unassembled',inputs=(b.stock('cdna_supplies','Actual source cDNA/library-preparation molecular supplies; whole oligos, enzymes and internal order unassembled'),),product='Actual source SSU-footprint cDNA library entering PDD; occupied molecules/adapter ends unassembled')
    probes=[]
    for j,(name,seq) in enumerate(PROBES.items()):probes.append(b.stock('pdd_probe_'+str(j),name+' sourceTable3 DNA PDD probe',K.OLIGO,R.STRUCTURAL,state=oligo(seq,name)))
    # Table3 records actual probe sequences but not original physical stock packaging or selected-library per-probe doses.
    b.op('pdd','Source used cDNA-level PDD rRNA depletion; whole reaction and fragment fates unassembled',inputs=tuple(probes)+(b.stock('pdd_supplies','Other actual PDD reaction supplies; source cited predecessor recipe not inherited without reading'),),product='Actual PDD-treated source footprint cDNA preparation; degraded fragments/retained strand states and removal yield unassembled')
    b.op('finish','Source applied TCP-seq library finish, PCR presence from Discussion2a; order/primers/cycles unassembled',inputs=(b.stock('finish_supplies','Actual source library-finishing supplies; whole primer/index inventory and packaging unassembled',K.OLIGO),))
    b.aliquot('loading','Actual selected TCP-seq sequencing allocation; amount and instrument/read plan unassembled')
    b.finish('Actual source selectedSSU TCP-seq library; complete molecular construct unassembled','Inspected2016publishedappliedworkflow/abstract confirmssequencedfootprints; platform, run lengths, wholeprimerlanding and full preparation remain unassembled. No standard Ribo-seq kit,2017recipe or readmask imported.');return _publication_coverage(g)

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    from zipfile import ZipFile
    import xml.etree.ElementTree as E
    c=Check();g=protocol_graph()
    c('First physical sediment follows actual fixation',g.steps['first_sediment'].inputs[0].material,'fix_product')
    c('Nonselected first fraction not invented discard',g.terminals['initial_nonselected'].disposition,Disposition.UNRESOLVED)
    c('SSU cDNA branch takes actual selected footprint fraction',g.steps['cdna_boundary'].inputs[0].material,'ssu_footprints')
    c('Ribosome alternative not pooled',g.terminals['ribosome_footprints'].disposition,Disposition.UNRESOLVED)
    c('All Table3 probes actual contributors',sum(i.material.startswith('pdd_probe_') for i in g.steps['pdd'].inputs),len(PROBES))
    c('PDD targets cDNA preparation not rawRNA',g.steps['pdd'].inputs[0].material,'cdna_boundary_product')
    c('No invented19–75 size selection',any(s.kind==S.SAMPLING for s in g.steps.values()),False)
    c('Loading remainder retained','loading_remainder' in g.terminals,True)
    p=sources_root()/Path(__file__).resolve().parents[1].name/'nature18647_41586_2016_BFnature18647_MOESM90_ESM.xlsx'
    if have(c,Source(p,'https://doi.org/10.1038/nature18647')):
        ns={'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
        with ZipFile(p) as z:
            ss=[''.join(x.itertext()) for x in E.fromstring(z.read('xl/sharedStrings.xml')).findall('m:si',ns)]
            original={}
            for r in E.fromstring(z.read('xl/worksheets/sheet1.xml')).findall('m:sheetData/m:row',ns):
                vals={}
                for cell in r.findall('m:c',ns):
                    v=cell.find('m:v',ns)
                    if v is not None:vals[cell.attrib['r']]=''.join(ss[int(v.text)].split()) if cell.attrib.get('t')=='s' else v.text
                names=[v for v in vals.values() if v.startswith(('18S-','25S-','5.8S-','5S-'))]
                if names:original[names[0]]=next(v for k,v in vals.items() if k.startswith('D'))
            for name,seq in PROBES.items():c('Whole originalTable3 ordering '+name,seq,original[name])
    c('Source-qualified publication gate',g.publication_problems(),[]);c.report()


def _cover_reviewed_targets(g,groups):
 """Assign only the literal curator-reviewed targets declared below."""
 for source,targets in groups:
  for target in targets:
   key='reviewed_'+target
   label=g.materials[target].label if target in g.materials else g.steps[target].label
   g.add_source_item(SourceItem(key,label,source));g.cover_source_item(SourceCoverage(key,'modelled',target))

PAGE_FILE='tcp-seq.html'
SOURCE_TARGETS=(
 ('Published2016SupplementaryDiscussion2a286–340: appliedformaldehydefixation/firstsediment, actualSSU/ribosomefractions andcDNAPDD/PCRpresence. Supplies are explicitunassembledoperationinterfaces, not inherited2017recipes.',
  ('source','formaldehyde','fix_conditions','fix','sediment_conditions','first_sediment','footprint_preparation','footprint_fractions','cdna_supplies','cdna_boundary','pdd_supplies','pdd','finish_supplies','finish','loading','load_library','library')),
 ('Published2016Table3MOESM90 originalA4:D53: all50wholeDNA PDDprobeorders, sourceIDsretained; actualdose/stockpackaging/occupiedsitesunassigned.',
  ('pdd_probe_0', 'pdd_probe_1', 'pdd_probe_2', 'pdd_probe_3', 'pdd_probe_4', 'pdd_probe_5', 'pdd_probe_6', 'pdd_probe_7', 'pdd_probe_8', 'pdd_probe_9', 'pdd_probe_10', 'pdd_probe_11', 'pdd_probe_12', 'pdd_probe_13', 'pdd_probe_14', 'pdd_probe_15', 'pdd_probe_16', 'pdd_probe_17', 'pdd_probe_18', 'pdd_probe_19', 'pdd_probe_20', 'pdd_probe_21', 'pdd_probe_22', 'pdd_probe_23', 'pdd_probe_24', 'pdd_probe_25', 'pdd_probe_26', 'pdd_probe_27', 'pdd_probe_28', 'pdd_probe_29', 'pdd_probe_30', 'pdd_probe_31', 'pdd_probe_32', 'pdd_probe_33', 'pdd_probe_34', 'pdd_probe_35', 'pdd_probe_36', 'pdd_probe_37', 'pdd_probe_38', 'pdd_probe_39', 'pdd_probe_40', 'pdd_probe_41', 'pdd_probe_42', 'pdd_probe_43', 'pdd_probe_44', 'pdd_probe_45', 'pdd_probe_46', 'pdd_probe_47', 'pdd_probe_48', 'pdd_probe_49',)),
)
def _publication_coverage(g):
 g.cover_source_item(SourceCoverage('selected_scope','modelled','pdd'))
 _cover_reviewed_targets(g,SOURCE_TARGETS)
 return g.require_publishable()

def render_page():
 from page import head,caveat,table
 from protocol_views import assemble_page,WorkflowSection
 g=protocol_graph()
 probe_rows=tuple((name,seq)for name,seq in PROBES.items())
 return assemble_page(g,[head('TCP-seq — translation-complex footprint preparation'),
 '<main class="wrap" data-protocol-complete="true"><h1>TCP-seq</h1><p>Formaldehyde fixation preserves translation complexes on mRNA. An initial sedimentation enriches mRNAs carrying a translating ribosome; footprint fractions then report small-subunit and ribosome positions.</p><p><a href="https://doi.org/10.1038/nature18647">Archer et al., Nature (2016)</a> · <a href="01_tcp-seq.html">Research note and evidence scope</a></p>',
 caveat('<strong>Applied preparation with incomplete finishing detail.</strong> This page follows the published 2016 Supplementary Discussion and original PDD probe table. Full 2016 main Methods remain unretrieved. Fixation dose, footprint digestion and fractionation details, complete cDNA/library chemistry and run setup remain unassembled. The later 2017 protocol is not substituted.'),
 '<h2>Preparation</h2><ol><li>Fix live yeast with formaldehyde to preserve initiating small-subunit complexes on mRNA; the executed concentration and time remain unassigned.</li><li>Collect the first ultracentrifugation sediment. This selects mRNAs with at least one translating ribosome; completely untranslated mRNA and free small subunits are outside the selected preparation.</li><li>Prepare footprint fractions and carry the selected small-subunit RNA into cDNA preparation. Digestion, end classes and fraction cuts remain unresolved here.</li><li>Apply probe-directed degradation (PDD) to the cDNA library with the source-table DNA probes, then finish and sequence the library. PCR is source-positive; whole primers, cycles and read geometry remain unassembled.</li></ol>',
 '<h2>Known molecular parts</h2>',table(('Part','Role and limit'),(('Formaldehyde','Source-applied fixation. The cited prior report’s 1–3% range is not assigned as TCP-seq’s executed dose.'),('First sediment','Physical enrichment for translating mRNPs; excluded material remains a separate unresolved fraction.'),('50 DNA PDD probes','Whole source orders are retained. PDD acts on cDNA, not by RNA-level rRNA removal; reaction details and individual probe doses remain unresolved.'))),
 '<details><summary>All 50 source-table PDD probe orders, 5′→3′</summary>',table(('Source ID','Whole DNA order'),probe_rows),'</details>',
 '<h2>Read setup</h2><p>The inspected sources establish sequenced footprints but do not assign an instrument, read lengths, index acquisition or complete primer landing sites in this declaration. These values remain unassigned in the model.</p>',
 '<h2>Source boundaries</h2><p>The applied preparation did not use stringent footprint size selection. The discussion’s approximately 19–75 nt range motivates a possible future selection; it is not an executed size window here. Table3 names the final probe “5S-1” while its target column says “5.8S rRNA”; the original ID and whole order are preserved without claiming a resolved target-site assignment.</p>',
 '<details><summary>Full workflow, physical fractions and unresolved chemistry</summary>',WorkflowSection('Source-qualified preparation'),'</details></main>'])

if __name__=='__main__':selftest()
