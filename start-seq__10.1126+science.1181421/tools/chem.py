"""Nechaev2010 Start-seq owns two independent end-reading library routes."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct
from pcr import PCRGeometry
from protocol_graph import Disposition, Highlight, InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind, Terminal
from protocol_model import role

TITLE="Start-seq — start-site-associated RNA sequencing"
NOTES="01_start_seq.html"
SOURCE=('Defining source: <a href="https://doi.org/10.1126/science.1181421">Nechaev et al., Science (2010)</a>; '
        '<a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM463298">original 5′-read deposit</a> and '
        '<a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM463301">original 3′-read deposit</a>.')
SUMMARY=("Short nuclear RNA is size-selected and enriched for capped species by sequential polyphosphatase and exonuclease treatment. "
         "The two source deposits use different library adapters to read the RNA 5′ or 3′ end in separate preparations.")
CAVEAT=("The displayed branches correspond to the untreated 5′-read deposit and LacZ-depleted 3′-read deposit; they are independent source samples. "
        "The original Science supporting methods remain route-blocked. Ordered adapter/RT/PCR strands, historical kit assignment and acquisition masks "
        "remain unassigned model data. The dated August 2010 manufacturer guide is a comparator, not the demonstrated experimental revision. "
        "No publication absence, exact complete depletion or shared index/RT-PCR recipe is asserted.")
RNA_RANGE=(25,100)
SEQ_PRIMERS=()
FINAL_LIBRARY=None
READ_LENGTHS={}
READ_SPECS=(sp.read_spec("Read 1",None,"Separate 5′-read library: RNA 5′ end; separate 3′-read library: RNA 3′ end. Actual acquisition masks unassigned.","include"),)

SEQUENCING_ENDING=("The deposits report Illumina Genome Analyzer II and separate 5′- and 3′-end libraries. "
                   "Trimming reads to 26 nt is a processing rule, not an acquisition-cycle declaration. "
                   "The ordered library ends, sequencing-primer sites and index acquisition are unassembled in this model.")


def state(name,*parts):
    return Construct([role(p) for p in parts],name=name)


FIVE_LIBRARY=state("5′-read Start-seq library — ordered ends unassembled","standard-library outer primer/adapter region",
                   "variable cDNA encoding the RNA 5-prime end","opposite standard-library primer/adapter region")
THREE_LIBRARY=state("3′-read Start-seq library — ordered ends unassembled","custom-library outer primer/adapter region",
                    "variable cDNA encoding the RNA 3-prime end","opposite custom-library primer/adapter region")
EXTRA_MOLECULAR_STATES=(("5′-read library",FIVE_LIBRARY,"Source-known roles from the standard Illumina small-RNA linker branch; exact bases and run assignment unassembled."),
                        ("3′-read library",THREE_LIBRARY,"Source-known roles from the custom-linker/custom-RT branch; no standard kit sequence is substituted."))


def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    polyphosphatase=g.input("polyphosphatase","5-prime polyphosphatase",MaterialKind.CATALYST)
    reviewed=[polyphosphatase.id]

    def branch(prefix,source_sample,condition,kind,final):
        def mid(name):return prefix+'_'+name
        def reagent(name,label,material_kind):
            m=g.input(mid(name),label,material_kind)
            reviewed.append(m.id)
            return m
        def product(name,label,material_kind=MaterialKind.MOLECULAR,parts=()):
            return g.product(mid(name),label,material_kind,state=state(label,*parts) if parts else None)
        def react(name,label,step_kind,inputs,outputs,**kw):
            reviewed.append(mid(name))
            g.react(mid(name),label,step_kind,inputs,outputs,**kw)
        rna=g.input(mid('nuclear_rna'),f"nuclear RNA from {source_sample}: {condition}",MaterialKind.BIOLOGICAL,
                    primary=True,source=f"GEO {source_sample}, original extraction protocol",
                    state=state('nuclear RNA','short capped and other nuclear RNAs'))
        reviewed.append(rna.id)
        gel=reagent('size_rule',f"{RNA_RANGE[0]}–{RNA_RANGE[1]}-nt nuclear RNA gel selection",MaterialKind.SELECTION)
        short=product('short_rna','size-selected nuclear RNA',parts=('variable 25–100-nt nuclear RNA',))
        excluded=product('excluded_sizes','gel-excluded RNA sizes',MaterialKind.WASTE)
        react('select_size','Gel-select 25–100-nt nuclear RNA',StepKind.SELECTION,
              ((rna,InputRole.TEMPLATE),(gel,InputRole.SELECTION)),(short,excluded))
        g.discard(excluded,'The original deposit selects the excised 25–100-nt RNA interval; other gel material is excluded.')
        poly=polyphosphatase
        treated=product('poly_treated','polyphosphatase-treated short RNA')
        react('treat_polyphosphate','Treat short RNA with 5′ polyphosphatase',StepKind.DEPHOSPHORYLATION,
              ((short,InputRole.TEMPLATE),(poly,InputRole.CATALYST)),(treated,))
        terminator=reagent('terminator','Terminator 5-prime exonuclease',MaterialKind.CATALYST)
        enriched=product('cap_enriched','remaining cap-enriched short RNA',parts=('cap-enriched short RNA; residual background unassigned',))
        digested=product('digested_uncapped','exonuclease-digested RNA products; subsequent recovery unassigned')
        react('digest_uncapped','Digest susceptible RNA with Terminator exonuclease',StepKind.DIGESTION,
              ((treated,InputRole.TEMPLATE),(terminator,InputRole.CATALYST)),(enriched,digested))
        g.mark_terminal(Terminal(digested.id,Disposition.UNRESOLVED,'These are a logical class of enzymatic degradation products. Their downstream recovery/removal and degradation efficiency are unassigned; no immediate physical discard is asserted.'))
        adapter3=reagent('adapter3',kind+' 3-prime small-RNA linker; actual ordered strand unassigned',MaterialKind.ADAPTER)
        ligase3=reagent('ligase3','source kit 3-prime-linker ligase; actual formulation unassigned',MaterialKind.CATALYST)
        linked3=product('linked3','3-prime-linker-bearing RNA')
        react('ligate3','Ligate the branch-specific 3′ linker',StepKind.LIGATION,
              ((enriched,InputRole.TEMPLATE),(adapter3,InputRole.STRUCTURAL),(ligase3,InputRole.CATALYST)),(linked3,),
              geometry=AdapterLigationGeometry.unavailable('The current declaration has not assigned this branch\'s ordered 3-prime linker or modified ligation-end geometry. The original supplement remains access-blocked; no absence is accepted.'))
        cleanup=reagent('post_ligation_gel','gel purification after 3-prime ligation',MaterialKind.SELECTION)
        clean=product('clean_linked3','gel-purified 3-prime-linked RNA')
        gel_waste=product('post_ligation_excluded','gel-excluded material; composition unassigned',MaterialKind.WASTE)
        react('purify_linked3','Gel-purify the 3′-linker product',StepKind.CLEANUP,
              ((linked3,InputRole.TEMPLATE),(cleanup,InputRole.SELECTION)),(clean,gel_waste))
        g.discard(gel_waste,'The source explicitly gel-purifies after 3-prime ligation; exact excluded species and complete removal are unassigned.')
        ap=reagent('alkaline_phosphatase','alkaline phosphatase',MaterialKind.CATALYST)
        dephos=product('dephosphorylated','phosphatase-treated linked RNA')
        react('dephosphorylate','Treat with alkaline phosphatase',StepKind.DEPHOSPHORYLATION,
              ((clean,InputRole.TEMPLATE),(ap,InputRole.CATALYST)),(dephos,))
        tap=reagent('tap','tobacco acid pyrophosphatase',MaterialKind.CATALYST)
        decapped=product('decapped','TAP-treated linked RNA ready for 5-prime-linker ligation')
        react('decap','Remove the cap with tobacco acid pyrophosphatase',StepKind.CHEMICAL_CONVERSION,
              ((dephos,InputRole.TEMPLATE),(tap,InputRole.CATALYST)),(decapped,))
        adapter5=reagent('adapter5',kind+' 5-prime small-RNA linker; actual ordered strand unassigned',MaterialKind.ADAPTER)
        ligase5=reagent('ligase5','source kit 5-prime-linker ligase; actual formulation unassigned',MaterialKind.CATALYST)
        linked5=product('linked5','two-linker RNA',parts=('branch-specific 5-prime linker','variable source RNA','branch-specific 3-prime linker'))
        react('ligate5','Ligate the branch-specific 5′ linker',StepKind.LIGATION,
              ((decapped,InputRole.TEMPLATE),(adapter5,InputRole.STRUCTURAL),(ligase5,InputRole.CATALYST)),(linked5,),
              geometry=AdapterLigationGeometry.unavailable('The ordered 5-prime linker and end geometry remain unassembled; no neighboring protocol supplies them.'))
        rt_label='custom 3-prime-read RT primer' if prefix=='three' else 'standard Illumina small-RNA kit RT primer'
        primer=reagent('rt_primer',rt_label+'; ordered strand unassigned',MaterialKind.OLIGO)
        rt=reagent('rt_enzyme','source kit reverse transcriptase; actual experimental identity unassigned',MaterialKind.CATALYST)
        cdna=product('cdna','branch-specific first-strand cDNA',parts=('linker-derived tract','variable RNA-derived cDNA','opposite linker-derived tract'))
        react('reverse_transcribe','Reverse-transcribe with the branch-specific primer',StepKind.REVERSE_TRANSCRIPTION,
              ((linked5,InputRole.TEMPLATE),(primer,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(cdna,))
        forward=reagent('pcr_forward','source kit PCR forward strand; actual ordered oligo unassigned',MaterialKind.OLIGO)
        reverse=reagent('pcr_reverse','source kit PCR reverse strand; actual ordered oligo unassigned',MaterialKind.OLIGO)
        polymerase=reagent('pcr_enzyme','source kit PCR reagents; actual formulation unassigned',MaterialKind.CATALYST)
        library=g.product(mid('library'),final.name,MaterialKind.LIBRARY,state=final)
        react('library_pcr','Amplify the branch-specific cDNA library',StepKind.PCR,
              ((cdna,InputRole.TEMPLATE),(forward,InputRole.STRUCTURAL),(reverse,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(library,),
              geometry=PCRGeometry.unavailable(cdna.state,reason='The model has not assigned the ordered experimental PCR strands and binding sites. This is not a claim that no source prints a primer.'))
        reviewed.append(library.id)
        g.sequence(library,SequencingPlan(library=final,primers=(),limitation='Separate end-reading library; exact acquired mask, physical primer and molecular entry geometry remain unassigned.'))
        return library

    branch('five','GSM463298','untreated 5-prime-read sample','standard Illumina small-RNA',FIVE_LIBRARY)
    branch('three','GSM463301','LacZ-depleted 3-prime-read sample','custom',THREE_LIBRARY)
    g.annotate(Highlight('five_library_pcr','Separate source preparations, not a paired-end run',
        'The displayed 5′ and 3′ libraries come from different deposited sample conditions. No literal split of one physical RNA sample is asserted. '
        'Both records report Genome Analyzer II; their 26-nt trimming rule is computational, not proof of acquisition length or a sample-index read.'))
    g.annotate(Highlight('three_ligate3','Historical kit and custom branch remain distinct',
        'The August 2010 manufacturer Rev. D guide is later than the original experiment and uses standard 5′-then-3′ ligation. '
        'The original deposits explicitly use 3′-then-5′ ligation with intervening gel, alkaline phosphatase and TAP. '
        'The custom 3′-read branch retains its own linkers and RT primer; no standard or TruSeq Small RNA bases, enzyme choices or index are substituted.'))
    g.declare_reviewed_scope(tuple(reviewed),SOURCE)
    return g.require_publishable()
