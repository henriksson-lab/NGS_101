"""Causal checks for source-qualified CEL-Seq, Drop-seq and scRRBS movies."""
import unittest
from dataclasses import replace
from .synthesis_examples import celseq,celseq_simulation,dropseq,dropseq_simulation,scrrbs,scrrbs_simulation,rc
from .simulator import Polymerase
from .molecular_frame import molecular_frame
from .evaluate import presentation,evaluate


class SynthesisTests(unittest.TestCase):
    def test_t7_is_de_novo_retains_dna_and_copies_source_barcode(self):
        r=celseq_simulation();a,b=r['op1'].before,r['op1'].after
        for id in [*r['top'],*r['bottom']]:self.assertEqual(a.atoms[id],b.atoms[id])
        self.assertEqual(a.next,{k:v for k,v in b.next.items()if k in a.atoms})
        start=r['start'];self.assertEqual(''.join(b.atoms[i].base for i in r['rna1']),''.join(a.atoms[i].base for i in r['top'][start:]).replace('T','U'))
        barcode=''.join(b.atoms[i].base for i in r['rna1']if 'cell barcode'in b.atoms[i].features)
        self.assertEqual(barcode,'CAUCACGC')
        self.assertFalse(any(i in b.paired for i in r['rna1']))
        self.assertEqual(b.five[r['rna1'][0]],'ppp')
        self.assertTrue(set(r['rna1']).isdisjoint(r['rna2']))
        self.assertTrue(all(i in r['after'].atoms for i in r['rna1']))

    def test_t7_requires_paired_promoter_and_transcribes_actual_insert(self):
        r=celseq_simulation();s=r['simulator'];prom=r['top'][slice(*r['promoter_span'])]
        s.unpair(prom,source='synthetic causal control')
        with self.assertRaisesRegex(ValueError,'paired promoter'):
            s.transcribe('template',len(r['bottom'])-1-r['start'],Polymerase('T7',('DNA',),'RNA',de_novo=True),promoter=prom,source='synthetic causal control')
        alternative=celseq_simulation(body='GATACCA')
        self.assertNotEqual(len(alternative['rna1']),len(r['rna1']))
        self.assertEqual(''.join(alternative['after'].atoms[i].base for i in alternative['rna1'])[-7:],'GAUACCA')

    def test_dropseq_copies_rna_from_existing_bead_primer(self):
        r=dropseq_simulation();a,b=r['extend'].before,r['extend'].after
        self.assertEqual(r['product'][:len(r['primer'])],tuple(r['primer']))
        self.assertEqual(len(r['extend'].added),len(r['body']))
        self.assertEqual(r['extend'].details['stop'],'physical_template_end')
        self.assertEqual(b.five[r['primer'][0]],'bead')
        for id in r['primer']:self.assertEqual(a.atoms[id],b.atoms[id])
        for id in r['rna']:self.assertEqual(a.atoms[id],b.atoms[id])
        self.assertTrue(all(b.atoms[i].sugar=='DNA'for i in r['product']))
        alternate=dropseq_simulation(barcode='GCTATGACTGAC')
        self.assertNotEqual(''.join(b.atoms[i].base for i in r['primer']),''.join(alternate['after'].atoms[i].base for i in alternate['primer']))
        self.assertEqual(''.join(b.atoms[i].base for i in r['extend'].added),''.join(alternate['after'].atoms[i].base for i in alternate['extend'].added))

    def test_dropseq_requires_rna_capable_polymerase_and_real_poly_a_contacts(self):
        r=dropseq_simulation();s=r['simulator']
        # On the bound initial state, a DNA-only polymerase cannot copy this RNA.
        from .simulator import Simulator
        t=Simulator();t.supply_nt('DNA','ACGT',source='control')
        t.add('rna','CGAUAC'+'A'*30,sugar='RNA',source='control');t.add('primer','T'*30,source='control')
        t.anneal('primer','rna',binding=30,source='control')
        with self.assertRaisesRegex(ValueError,'template sugar'):
            t.extend('primer',Polymerase('DNA only',('DNA',),'DNA'),source='control')
        t.unpair(t.chain('rna'),source='control')
        with self.assertRaisesRegex(ValueError,'annealed'):
            t.extend('primer',Polymerase('RT',('RNA',),'DNA'),source='control')

    def test_bisulfite_changes_original_bases_not_identity_sugar_or_backbone(self):
        r=scrrbs_simulation();a,b=r['conversion'].before,r['conversion'].after
        self.assertEqual(a.next,b.next);self.assertEqual(set(a.atoms),set(b.atoms))
        for id in r['original']:
            self.assertEqual(a.atoms[id].sugar,b.atoms[id].sugar)
            self.assertEqual(b.atoms[id].base,'U'if id in r['sites']else a.atoms[id].base)
        protected=r['original'][r['protected']]
        self.assertEqual(b.atoms[protected].base,'C');self.assertEqual(b.atoms[protected].modification,'5mC')
        with self.assertRaisesRegex(ValueError,'protected'):
            r['simulator'].convert([protected],from_base='C',to_base='U',protected=('5mC',),source='synthetic causal control')

    def test_bisulfite_U_A_T_are_three_different_molecules(self):
        r=scrrbs_simulation();s=r['after'];n=len(r['original'])
        for j,id in enumerate(r['original']):
            if id not in r['sites']:continue
            opposite=r['dna1'][n-1-j];reported=r['dna2'][j]
            self.assertEqual([s.atoms[i].base for i in [id,opposite,reported]],['U','A','T'])
            self.assertEqual(len({id,opposite,reported}),3)
        self.assertTrue(all(id in s.atoms for id in r['original']))
        self.assertFalse(any(id in s.paired for id in r['original']))
        changed=scrrbs_simulation(protected=3)
        self.assertNotEqual([s.atoms[i].base for i in r['dna2']],[changed['after'].atoms[i].base for i in changed['dna2']])

    def test_all_incoming_transitions_preserve_identical_molecular_frames(self):
        for factory in [celseq,dropseq,scrrbs]:
            m=factory()
            for before,after in zip(m.shots,m.shots[1:]):
                self.assertGreater(after.transition,0)
                self.assertEqual(molecular_frame(before,before.timing()[2]),molecular_frame(after,0),(m.title,after.id))
                state,actors,camera=presentation(m,after,0)
                self.assertEqual(camera,before.camera)
                for st in state.values():self.assertEqual(st.molecular,molecular_frame(after,0))

    def test_every_active_nucleotide_has_an_explicit_retained_track(self):
        for factory in [celseq,dropseq,scrrbs]:
            for shot in factory().shots:
                drawn={i for st in shot.strands for i in st.nucleotide_ids}
                for seconds in [0,shot.establish+shot.action*.5,shot.timing()[2]]:
                    frame=molecular_frame(shot,seconds)
                    self.assertTrue((set(frame.atoms)-frame.released)<=drawn,(factory.__name__,shot.id))

    def test_stock_primers_arrive_intact_then_bind_without_synthesis(self):
        m=scrrbs();shot=next(s for s in m.shots if s.id=='primer');event=shot.events[0]
        start=shot.establish+shot.action*event.fraction_start
        before=molecular_frame(shot,start-1e-6);supplied=molecular_frame(shot,start)
        ids={a['id']for a in event.molecular_operation['added']}
        self.assertFalse(ids&set(before.atoms));self.assertTrue(ids<=set(supplied.atoms))
        self.assertEqual(len(supplied.atoms)-len(before.atoms),5)
        self.assertEqual(len(supplied.backbone)-len(before.backbone),4)


if __name__=='__main__':unittest.main(verbosity=2)
