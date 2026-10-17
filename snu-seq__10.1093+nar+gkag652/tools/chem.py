"""SNU-seq owns its pulse-label selection and two dated kit references."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
import illumina as il
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, feature, revcomp
from pcr import PCRGeometry
from protocol_graph import (Disposition, Highlight, InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, SourceClaim, SourceDiscrepancy, StepKind, Terminal)
from protocol_model import role

TITLE="SNU-Seq — single-nucleotide-resolution 4sU sequencing"
NOTES="01_snu_seq.html"
SOURCE=('Defining source: <a href="https://doi.org/10.1093/nar/gkag652">Gerlevik et al., NAR (2026)</a>, original Supplementary Methods; '
        '<a href="https://www.lexogen.com/wp-content/uploads/2018/06/012UG036V0130_QuantSeq-IonTorrent.pdf">2018 Ion guide</a> and '
        '<a href="https://www.lexogen.com/wp-content/uploads/2021/08/015UG009V0260_QuantSeq_Illumina_2021-08-16.pdf">2021 FWD guide</a>, conditional references.')
SUMMARY=("RNA containing pulse-incorporated 4sU is enriched using reversible HPDP-biotin capture, recovered with DTT and polyadenylated with bacterial PAP. "
         "Two platform-specific QuantSeq routes use oligo(dT) first-strand synthesis, RNA removal and random-primed second-strand synthesis before PCR.")
CAVEAT=("The graph begins with recovered pulse-labelled total RNA and selects the 2026 +bPAP SNU route without intentional sonication. "
        "Labelled molecules may include pre-existing RNA tracts. Optional rRNA-depletion and no-bPAP/no-4sU controls are distinct variants. "
        "INFERRED — displayed Ion Torrent and Illumina FWD endpoints apply dated manufacturer architectures as conditional references; "
        "actual study guide revisions, ordered kit primers and barcode/index members remain unassigned. Exact RNA and copied poly(A) lengths are unknown. "
        "Reads proceed towards the poly(A) end; reaching the native endpoint depends on insert, read length and processing.")
PCR_CYCLES=13
READ_CYCLES=100
SEQ_PRIMERS=()
FINAL_LIBRARY=None
READ_SPECS=(sp.read_spec("Read 1",READ_CYCLES,"Ion reference: TCAG–inline barcode–CGAT then RNA; FWD reference: RNA towards poly(A). Native endpoint reach is conditional.","include"),
            sp.read_spec("Index 1 (i7)",None,"Ion has no separate index read in the dated guide; FWD study index acquisition/member unassigned","undetermined"))
SEQUENCING_ENDING=("The 2026 supplement reports single-end 100-nt sequencing (SR100): Ion Chef/Ion Proton or Lexogen Illumina NextSeq 2000. "
                   "Poly(A) trimming and aligned-read 3′ coverage are computational operations; they do not establish that every acquired read reaches a native RNA end.")
ION_A="CCATCTCATCCCTGCGTGTCTCCGACTCAG"
ION_P1="ATCACCGACTGCCCATAGAGAGGAAAGCGGAGGCGTAGTGG"
ION_SPACER="CGAT"


def state(name,*parts):
    return Construct([role(p) for p in parts],name=name)


def ion_library(*,inferred=True):
    return Construct([
        Segment("Ion A arm",ION_A,"p5",inferred=inferred),
        role("inline barcode, source family 9–11 nt",tag="cbc",molecular_feature=feature("ion_inline_sample","sample_index","unknown"),inferred=inferred),
        Segment("post-barcode spacer",ION_SPACER,inferred=inferred),
        role("variable copied RNA near its native 3-prime end"),role("copied poly(A) and primer-derived junction; length unassigned"),
        Segment("Ion P1 arm",ION_P1,"p7",inferred=inferred)],name="INFERRED — 2018 Ion Torrent library reference")


def fwd_library(*,inferred=True):
    return Construct([
        Segment("P5-side cap",il.P5[:-4],"p5",inferred=inferred),
        Segment("FWD Read 1 site",il.TRUSEQ_READ1,"r1",inferred=inferred),
        role("variable copied RNA near its native 3-prime end"),role("copied poly(A) and primer-derived junction; length unassigned"),
        Segment("FWD opposite adapter arm",revcomp(il.TRUSEQ_READ2),"r2",inferred=inferred),
        Segment("reference i7 family",'N'*6,"cbc",placeholder=True,inferred=inferred,
                feature=feature("fwd_i7_sample","sample_index","unknown",note="dated guide's six-base family, not selected study member")),
        Segment("P7-side cap",il.P7_RC,"p7",inferred=inferred)],name="INFERRED — 2021 Illumina FWD library reference")


ION_LIBRARY=ion_library()
FWD_LIBRARY=fwd_library()
FWD_R1=sp.custom("Read 1","FWD guide printed Read 1 primer",il.TRUSEQ_READ1,"Lexogen 015UG009V0260 Appendix G; conditional guide assignment")
EXTRA_MOLECULAR_STATES=(("Ion Torrent reference endpoint",ION_LIBRARY,"INFERRED — source-known 2018 architecture, conditional study assignment; no barcode member or physical sequencing primer selected."),
                        ("Illumina FWD reference endpoint",FWD_LIBRARY,"INFERRED — source-known 2021 single-i7 architecture, conditional study assignment; optional UMI/UDI modules not selected."))
TOPOLOGY_PANELS=(("INFERRED — FWD guide Read 1 entry reference",
                  sp.section(FWD_LIBRARY,(FWD_R1,),required_roles=("Read 1",),read_lengths={"Read 1":READ_CYCLES},
                             intro="Conditional 2021 FWD guide assignment. Read 1 starts from the random second-strand side towards poly(A); variable inserts prevent a fabricated fixed-coordinate endpoint.")),)


def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    reviewed=[]
    def inp(id,label,kind,**kw):
        m=g.input(id,label,kind,**kw);reviewed.append(id);return m
    def prod(id,label,kind=MaterialKind.MOLECULAR,**kw):return g.product(id,label,kind,**kw)
    def react(id,label,kind,inputs,outputs,**kw):
        reviewed.append(id);g.react(id,label,kind,inputs,outputs,**kw)
    pool=Scene()
    pool.strand('pulse_labelled',[role('RNA with pulse-incorporated 4sU, potentially including pre-existing tracts')],label='label-containing RNA population')
    pool.strand('without_label',[role('RNA without pulse-installed label')],label='RNA without pulse label')
    rna=inp('total_rna','recovered total RNA after 4sU pulse and source spike-in addition',MaterialKind.BIOLOGICAL,primary=True,
            state=pool)
    hpdp=inp('hpdp','EZ-Link HPDP-biotin',MaterialKind.MOLECULAR)
    hpdp_condition=inp('hpdp_condition','source HPDP reaction incubation',MaterialKind.CONDITION)
    labelled=prod('biotinylated','HPDP-biotinylated pulse-labelled RNA within the RNA pool')
    react('biotinylate','Biotinylate incorporated 4sU thiols',StepKind.CHEMICAL_CONVERSION,
          ((rna,InputRole.TEMPLATE),(hpdp,InputRole.STRUCTURAL),(hpdp_condition,InputRole.CONDITION)),(labelled,))
    preclean=inp('pre_capture_cleanup','phase-lock chloroform extraction and isopropanol recovery',MaterialKind.CONDITION)
    cleaned=prod('cleaned','RNA recovered before affinity capture')
    preclean_waste=prod('preclean_waste','non-retained extraction/precipitation fractions',kind=MaterialKind.WASTE)
    react('remove_free_biotin','Clean RNA before streptavidin capture',StepKind.CLEANUP,
          ((labelled,InputRole.TEMPLATE),(preclean,InputRole.CONDITION)),(cleaned,preclean_waste))
    g.discard(preclean_waste,'The source proceeds with recovered RNA; exact composition of discarded cleanup fractions remains unassigned.')
    beads=inp('streptavidin','µMACS streptavidin beads/columns and source wash',MaterialKind.SELECTION)
    bound=prod('bound','column-associated biotinylated RNA')
    flow=prod('flowthrough','non-retained RNA and wash fractions',MaterialKind.WASTE)
    react('capture','Capture biotinylated RNA and wash the column',StepKind.AFFINITY_CAPTURE,
          ((cleaned,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),(bound,flow))
    g.discard(flow,'The selected library route proceeds with RNA released from the column; the discarded flow-through is not asserted to be all and only old RNA.')
    dtt=inp('dtt','DTT for labelled-RNA release',MaterialKind.CONDITION)
    eluate=prod('eluate','DTT-eluted RNA containing pulse-incorporated label')
    column=prod('column_after_elution','column/beads after elution; residual contents unassigned')
    react('release','Release labelled RNA with DTT',StepKind.CHEMICAL_CONVERSION,
          ((bound,InputRole.TEMPLATE),(dtt,InputRole.CONDITION)),(eluate,column))
    g.mark_terminal(Terminal(column.id,Disposition.UNRESOLVED,'Remaining column contents and downstream disposal are unassigned; no quantitative elution efficiency is asserted.'))
    recovery=inp('rna_recovery','miRNeasy recovery with on-column DNase I',MaterialKind.CONDITION)
    recovered=prod('recovered','recovered enriched RNA',state=state('enriched RNA','intact source RNA containing pulse-incorporated 4sU'))
    rna_waste=prod('rna_recovery_waste','non-retained RNA-column cleanup fractions',kind=MaterialKind.WASTE)
    react('recover_rna','Recover RNA and remove DNA contamination',StepKind.CLEANUP,
          ((eluate,InputRole.TEMPLATE),(recovery,InputRole.CONDITION)),(recovered,rna_waste))
    g.discard(rna_waste,'Non-retained column fractions are excluded from the recovered RNA; molecular composition is unassigned.')
    pap=inp('bpap','NEB E. coli poly(A) polymerase',MaterialKind.CATALYST)
    atp=inp('atp','ATP for bacterial polyadenylation',MaterialKind.MOLECULAR)
    tailed=prod('tailed','bPAP-tailed source RNA',state=state('PAP-tailed RNA','variable source RNA with retained native 3-prime boundary','artificial poly(A) tail, unknown length'))
    react('polyadenylate','Add artificial poly(A) tails to enriched RNA',StepKind.TAILING,
          ((recovered,InputRole.TEMPLATE),(pap,InputRole.CATALYST),(atp,InputRole.STRUCTURAL)),(tailed,))
    precipitation=inp('post_pap_cleanup','post-PAP isopropanol precipitation with 1/10 volume sodium acetate and GlycoBlue; pellet resuspended in water',MaterialKind.CONDITION)
    pap_rna=prod('pap_rna','recovered PAP-tailed RNA',state=tailed.state)
    pap_waste=prod('post_pap_waste','non-retained post-PAP precipitation supernatant',kind=MaterialKind.WASTE)
    react('recover_tailed','Recover polyadenylated RNA before kit preparation',StepKind.CLEANUP,
          ((tailed,InputRole.TEMPLATE),(precipitation,InputRole.CONDITION)),(pap_rna,pap_waste))
    g.discard(pap_waste,'The source resuspends the RNA pellet; the non-retained fraction is not assigned an exact species composition.')
    ion_rna=prod('ion_rna','RNA assigned to the Ion library branch',state=tailed.state)
    fwd_rna=prod('fwd_rna','RNA assigned to the Illumina FWD library branch',state=tailed.state)
    react('platform_variants','Represent separate source-reported platform variants',StepKind.MAGIC,
          ((pap_rna,InputRole.TEMPLATE),),(ion_rna,fwd_rna),changes_molecule=False,
          reason='This comparison boundary represents alternative source-reported platform implementations. Actual sample assignment and aliquot relationships are unassigned; it is not a verified physical split or biochemical transformation.')

    def library_branch(prefix,substrate,final,first_label,second_label,removal_label):
        # These dated guide operations are conditional assignments to the study.
        # Keep that scope on every intermediate, rather than only the endpoint.
        def add(name,label,kind,**kw):
            kw.setdefault('inferred',True)
            return inp(prefix+'_'+name,label,kind,**kw)
        def product(name,label,**kw):
            kw.setdefault('inferred',True)
            return prod(prefix+'_'+name,label,**kw)
        def step(name,label,kind,inputs,outputs,**kw):
            if not label.startswith('INFERRED'):
                label='INFERRED dated guide assignment — '+label
            react(prefix+'_'+name,label,kind,inputs,outputs,**kw)
        primer=add('dt_primer',first_label+'; actual ordered reagent unassigned',MaterialKind.OLIGO)
        rt=add('rt_enzyme','dated guide RT enzyme mix E1',MaterialKind.CATALYST)
        first=product('first_strand','oligo(dT)-primed first-strand cDNA/RNA mixture',state=state('first strand','kit first-strand handle','oligo(dT)-primed variable source cDNA'))
        step('rt','Reverse-transcribe from the artificial poly(A) tail',StepKind.REVERSE_TRANSCRIPTION,
             ((substrate,InputRole.TEMPLATE),(primer,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(first,))
        removal=add('rna_removal',removal_label,MaterialKind.CONDITION)
        cdna=product('cdna','cDNA after RNA-template degradation',state=first.state)
        rna_products=product('rna_products','degraded RNA-template products; subsequent fate unassigned')
        step('remove_rna','Remove the RNA template before second-strand synthesis',StepKind.CHEMICAL_CONVERSION,
             ((first,InputRole.TEMPLATE),(removal,InputRole.CONDITION)),(cdna,rna_products))
        g.mark_terminal(Terminal(rna_products.id,Disposition.UNRESOLVED,'RNA degradation is source-described; recovery of individual degradation products through later cleanup is not assigned.'))
        random=add('random_primer',second_label+'; full ordered strand unassigned',MaterialKind.OLIGO)
        enzyme=add('second_enzyme','dated guide second-strand DNA polymerase mix E2',MaterialKind.CATALYST)
        ds=product('double_stranded','random-primed double-stranded DNA',state=state('second strand','kit second-strand handle','variable RNA-derived DNA interval','first-strand-side handle'))
        step('second_strand','Random-prime the second DNA strand',StepKind.EXTENSION,
             ((cdna,InputRole.TEMPLATE),(random,InputRole.STRUCTURAL),(enzyme,InputRole.CATALYST)),(ds,))
        cleanup=add('cleanup','dated guide magnetic bead purification',MaterialKind.SELECTION)
        clean=product('clean_ds','purified double-stranded DNA',state=ds.state)
        waste=product('cleanup_waste','non-retained cleanup fraction',kind=MaterialKind.WASTE)
        step('purify_ds','Purify DNA before PCR',StepKind.CLEANUP,
             ((ds,InputRole.TEMPLATE),(cleanup,InputRole.SELECTION)),(clean,waste))
        g.discard(waste,'Non-retained fraction is removed during the source-described bead purification; exact molecular composition is unassigned.')
        p1=add('pcr_forward','kit PCR forward reagent; ordered strand unassigned',MaterialKind.OLIGO)
        p2=add('pcr_reverse','kit PCR reverse/index reagent; selected ordered member unassigned',MaterialKind.OLIGO)
        polymerase=add('pcr_mix',f'kit PCR mix; 2026 source selects {PCR_CYCLES} cycles',MaterialKind.CATALYST)
        amplified=product('amplified',final.name, state=final,inferred=True,source=SOURCE+'; dated manufacturer reference assignment remains conditional')
        step('pcr',f'INFERRED — amplify using the dated platform reference; {PCR_CYCLES} source PCR cycles',StepKind.PCR,
             ((clean,InputRole.TEMPLATE),(p1,InputRole.STRUCTURAL),(p2,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(amplified,),
             geometry=PCRGeometry.unavailable(clean.state,reason='The final vendor architecture is source-positive, but actual study ordering primers, selected members and exact annealing geometry remain unassigned; no source absence is accepted.'))
        post_cleanup=add('post_pcr_cleanup','dated guide post-PCR magnetic purification',MaterialKind.SELECTION)
        library=product('library',final.name,kind=MaterialKind.LIBRARY,state=final,inferred=True,source=SOURCE+'; dated guide assignment remains conditional')
        post_waste=product('post_pcr_waste','non-retained PCR cleanup fraction',kind=MaterialKind.WASTE)
        step('purify_library','INFERRED — purify the dated reference library after PCR',StepKind.CLEANUP,
             ((amplified,InputRole.TEMPLATE),(post_cleanup,InputRole.SELECTION)),(library,post_waste))
        g.discard(post_waste,'The dated guides physically remove non-retained post-PCR cleanup fractions; exact species composition is unassigned.')
        reviewed.append(library.id)
        g.sequence(library,SequencingPlan(library=final,primers=(),acquired_roles=('Read 1',),read_lengths={'Read 1':READ_CYCLES},
                                        limitation='SR100 is source-reported; exact study primer/member assignment and native-end reach remain conditional.'))

    library_branch('ion',ion_rna,ION_LIBRARY,'Ion P1-side oligo(dT) primer','Ion A-side random second-strand primer with inline barcode',
                   'Ion guide RNA Removal Solution 1, heat, then Solution 2')
    library_branch('fwd',fwd_rna,FWD_LIBRARY,'Illumina Read-2-side oligo(dT) primer','FWD Read-1-side random second-strand primer',
                   'FWD guide RNA Removal Solution with heat')
    g.annotate(Highlight('polyadenylate','Published route and control choices',
        'The 2026 article Methods and Supplementary Methods report 150 ng thio-labelled RNA polyadenylated with E. coli PAP for 45 min at 37°C, '
        'Lexogen QuantSeq Ion Torrent or FWD-for-Illumina preparation with 13 PCR cycles, and SR100 sequencing. '
        'Superseded bioRxiv versions are not used as authorities. '
        'No-bPAP, no-4sU, optional rRNA-depleted and size-fractionated workflows remain distinct source variants.'))
    g.annotate(Highlight('platform_variants','Platform variants are declarations of separate implementations',
        'The branch point compares source-reported platform routes; it is not evidence that a particular experimental sample was split between platforms. '
        'Ion inline barcodes are introduced in second-strand synthesis; the FWD guide introduces i7 by PCR. No selected barcode/index member is assigned.'))
    g.annotate(Highlight('fwd_library','Readout is conditional on reaching the RNA end',
        'Both dated guides read towards poly(A) from the random-primer side. FWD is distinct from the REV/CSP workflow. '
        'Poly(A) clipping, genomic-A masking, quality/mapping filters and aligned-read 3′ coverage are computational. '
        'Neither 100 acquired nucleotides nor genomecov -3 proves that every read terminates at the native RNA 3′ boundary.'))
    g.annotate(SourceDiscrepancy('capture','Wash salt differs between method versions',
        (SourceClaim('2026 original Supplementary Methods p18, text line338','The wash buffer contains 1 mM NaCl.'),
         SourceClaim('Dölken 2008 original Methods p1970, local text966–968','The wash buffer contains 1 M NaCl.')),
        'The declaration preserves capture and washing without choosing a salt concentration.',
        'Both printed original pages were visually inspected, in addition to their text twins.',
        'The 2026 method explicitly includes modifications. This is a source-version difference; no typographical error or corrected concentration is assumed.'))
    g.declare_reviewed_scope(tuple(reviewed),SOURCE)
    return g.require_publishable()
