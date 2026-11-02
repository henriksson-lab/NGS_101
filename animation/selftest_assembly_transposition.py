"""Causal mechanism and continuity checks for the three assembly examples."""
import math
import unittest
from .assembly_transposition import gibson,gibson_simulation,atac,atac_simulation,rnaseh,rnaseh_simulation
from .evaluate import evaluate,presentation
from .molecular_frame import molecular_frame
from .simulator import Simulator,Polymerase


class AssemblyTranspositionTests(unittest.TestCase):
    def test_exact_overlap_stalls_without_manufacturing_gap_synthesis(self):
        m,r=gibson_simulation()
        self.assertFalse(r['gap'].added)
        self.assertEqual(r['gap'].details['stop'],'occupied_template')
        self.assertEqual(r['gap'].before.atoms,r['gap'].after.atoms)
        for op in r['seals']:
            self.assertEqual(op.before.atoms,op.after.atoms)
            self.assertEqual(len(set(op.after.next.items())-set(op.before.next.items())),1)
        removed={id for op in r['resect_ops']for id in op.details['removed']}
        self.assertTrue(removed<=r['simulator'].released)
        self.assertTrue(removed.isdisjoint(r['joined_top']+r['joined_bottom']))

    def test_tn5_preserves_target_and_creates_two_transfer_bonds_per_complex(self):
        r=atac_simulation()
        for op in [r['transfer_left'],r['transfer_right']]:
            self.assertEqual(op.before.atoms,op.after.atoms)
            self.assertEqual(len(set(op.before.next.items())-set(op.after.next.items())),2)
            self.assertEqual(len(set(op.after.next.items())-set(op.before.next.items())),2)
        self.assertTrue(set(r['target_top_ids']+r['target_bottom_ids'])<=set(r['repair_after'].atoms))
        self.assertFalse(set(r['target_top_ids']+r['target_bottom_ids'])&r['repair_after'].released)

    def test_all_loaded_me_duplexes_stay_compact_through_docking(self):
        shot=atac().shots[0]
        duplexes=[p for p in shot.pairings if p.id.endswith('_me')]
        self.assertEqual(len(duplexes),4)
        # No viewport filtering: incoming/outgoing offscreen oligos are still
        # physical paired molecules and must use the same geometry.
        for phase in [0,.05,.1,.2,.34,.4,.6,.7,.85,.95,1]:
            states,_=evaluate(shot,shot.establish+shot.action*phase)
            for p in duplexes:
                a,b=states[p.a],states[p.b]
                for j in range(p.length):
                    ia,ib=p.a_start+j,p.b_start-j
                    pa=a.segment_point(ia+.5,ia);pb=b.segment_point(ib+.5,ib)
                    self.assertAlmostEqual(math.dist(pa,pb),24,places=7,msg=(phase,p.id,j))

    def test_repair_stalls_on_short_me_then_copies_available_adapter(self):
        r=atac_simulation()
        for side in ['top','bottom']:
            gap=r[side+'_gap'];copy=r[side+'_adapter_copy']
            self.assertEqual(gap.details['stop'],'occupied_template')
            self.assertEqual(len(gap.added),9)
            self.assertTrue(copy.added)
            self.assertEqual(copy.details['stop'],'physical_template_end')
            ids=r['repaired_'+side+'_ids'];added=set(gap.added+copy.added)
            self.assertEqual(set(ids)-set(r['tag_'+side+'_ids']),added)
        before=r['tag_before'].atoms
        altered=list(r['target']);altered[45]='A'if altered[45]!='A'else'C'
        changed=atac_simulation(target=''.join(altered))
        self.assertNotEqual(r['repaired_top'],changed['repaired_top'])

    def test_invalid_stagger_or_external_cut_fails_instead_of_drawing_product(self):
        with self.assertRaises((ValueError,IndexError)):
            atac_simulation(cuts=(0,88))

    def test_rnase_cuts_rna_only_and_dissociation_is_separate(self):
        r=rnaseh_simulation()
        for op in r['cuts']:
            self.assertEqual(op.before.atoms,op.after.atoms)
            self.assertEqual(op.before.paired,op.after.paired)
            self.assertEqual(len(set(op.before.next.items())-set(op.after.next.items())),1)
            self.assertFalse(op.added)
        for id in r['dna']:
            self.assertEqual(r['before'].atoms[id],r['after'].atoms[id])
        self.assertEqual(r['cleaved'].next,r['after'].next)
        self.assertEqual(r['cleaved'].atoms,r['after'].atoms)
        self.assertFalse(r['after'].released)
        self.assertEqual(len(r['cleaved'].paired)-len(r['after'].paired),32)
        self.assertTrue(all(r['after'].atoms[id].sugar=='RNA'for id in r['rna']))

    def test_reframing_keeps_the_same_physical_state_and_nucleotide_coordinates(self):
        for factory in [gibson,atac,rnaseh]:
            m=factory()
            for previous,shot in zip(m.shots,m.shots[1:]):
                self.assertEqual(molecular_frame(previous,previous.timing()[2]),molecular_frame(shot,0),(factory.__name__,shot.id))
                old,_=evaluate(previous,previous.timing()[2]);actual,_,camera=presentation(m,shot,0)
                self.assertEqual(camera,previous.camera)
                locations={id:(st,j)for st in old.values()for j,id in enumerate(st.strand.nucleotide_ids)if st.fraction(j)}
                for st in actual.values():
                    for j,id in enumerate(st.strand.nucleotide_ids):
                        if st.fraction(j):
                            other,k=locations[id]
                            self.assertLess(math.dist(st.point(j+.5),other.segment_point(k+.5,k)),1e-6)
                presentation(m,shot,shot.transition/2)

    def test_all_active_atoms_have_tracks_including_nonselected_fragments(self):
        for factory in [gibson,atac,rnaseh]:
            for shot in factory().shots:
                represented={id for track in shot.strands for id in track.nucleotide_ids}
                for t in [0,shot.establish+shot.action*.5,shot.timing()[2]]:
                    frame=molecular_frame(shot,t)
                    self.assertFalse(set(frame.atoms)-frame.released-represented,(factory.__name__,shot.id))

    def test_movies_roundtrip_and_share_source_trace(self):
        from .model import Movie
        import json
        for factory in [gibson,atac,rnaseh]:
            m=factory();d=json.loads(json.dumps(m.to_dict()))
            self.assertEqual(Movie.from_dict(d).to_dict(),Movie.from_dict(Movie.from_dict(d).to_dict()).to_dict())
            self.assertTrue(m.provenance['molecular_simulation'])


if __name__=='__main__':unittest.main(verbosity=2)
