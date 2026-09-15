"""Numerical consistency/regression checks, not experimental validation."""
import unittest
import numpy as np
from xfem import Model, quadrature, enrichment, strike, WIDTH, HEIGHT


class NumericalChecks(unittest.TestCase):
    def test_subcell_area_is_conserved(self):
        model = Model('coarse')
        cracks = [{'y': .073, 'length': .112}, {'y': .142, 'length': .087}]
        area = sum(quadrature(v, cracks, model.h)[1].sum() for v, _, _ in model.geometry)
        self.assertAlmostEqual(area, WIDTH*HEIGHT, places=12)

    def test_tip_gradient_matches_finite_difference(self):
        c = {'y': .1, 'length': .12}
        p = np.array([[.16, .112], [.204, .083]])
        _, gradient = enrichment(p, c, 'tip', .025)
        for axis in range(2):
            delta = np.zeros(2); delta[axis] = 1e-7
            numerical = (enrichment(p+delta, c, 'tip', .025)[0]-enrichment(p-delta, c, 'tip', .025)[0])/(2e-7)
            np.testing.assert_allclose(gradient[:, axis], numerical, rtol=1e-7)

    def test_linear_elastic_scaling_and_equilibrium(self):
        model = Model('coarse')
        cracks = [{'y': .103, 'length': .13}]
        a, b = model.solve(cracks, 20e9, .103), model.solve(cracks, 40e9, .103)
        self.assertGreater(a['enriched_dofs'], 0)
        self.assertLess(a['residual'], 1e-7)
        self.assertAlmostEqual(a['energy'], 2*b['energy'], delta=a['energy']*1e-7)
        np.testing.assert_allclose(a['values'], b['values'], rtol=1e-6)

    def test_intact_reflection_symmetry(self):
        model = Model('medium')
        a, b = model.solve([], 20e9, .06), model.solve([], 20e9, .14)
        # Diagonal orientation introduces a small discretization asymmetry.
        self.assertLess(abs(a['energy']/b['energy']-1), .02)

    def test_growth_multiple_cracks_and_separation(self):
        for mesh in ['coarse', 'medium', 'fine']:
            with self.subTest(mesh=mesh):
                p = dict(mesh=mesh, young=70, poisson=.22, toughness=10, mass=.25, height=.01, aim=.3, cracks=[])
                self.assertFalse(strike(p)['accepted'])
                p.update(mass=1, height=.3)
                r = strike(p); self.assertTrue(r['accepted']); p['cracks'] = r['cracks']
                p['aim'] = .7; r = strike(p); self.assertEqual(len(r['cracks']), 2)
                p['cracks'], p['aim'] = r['cracks'], .3
                lengths = [c['length'] for c in r['cracks']]
                r = strike(p)
                self.assertGreater(r['cracks'][0]['length'], lengths[0])
                self.assertEqual(r['cracks'][1]['length'], lengths[1])
                for _ in range(45):
                    p['cracks'] = r['cracks']
                    r = strike(p)
                    if r['fractured']:
                        break
                self.assertTrue(r['fractured'])
                self.assertEqual(r['cracks'][0]['length'], WIDTH)

    def test_invalid_inputs_fail_cleanly(self):
        p = dict(mesh='medium', young=70, poisson=.22, toughness=10, mass=.5, height=.12, aim=.5, cracks=[])
        for key, value in [('young', 0), ('poisson', .6), ('toughness', float('nan')), ('aim', 1)]:
            with self.subTest(key=key), self.assertRaises(ValueError):
                strike({**p, key: value})


if __name__ == '__main__':
    unittest.main()
