import unittest
import numpy as np
from .alpha_transport_environment import relocate_macro,macro_semantics,compare_photons,photon_records


class EnvironmentControls(unittest.TestCase):
    def test_macro_changes_only_paths(self):
        body='/gagg/run/seed 123\n/gagg/run/outputFile old.root\n/gagg/sipm/configFile old.json\n/gun/energy 8.0 MeV\n/run/beamOn 1\n'
        changed=relocate_macro(body,'new.root','new.json')
        self.assertEqual(macro_semantics(body),macro_semantics(changed))
        self.assertIn('/gagg/run/seed 123\n',changed)
        with self.assertRaises(ValueError):relocate_macro('/run/beamOn 1','x','y')

    def test_optical_identity_distinguishes_order_and_physics(self):
        a=dict(time_ps=np.array([1,2,2]),x_mm=np.array([0.,1.,2.],np.float32),
               y_mm=np.zeros(3,np.float32),energy_eV=np.full(3,2.5,np.float32))
        b={k:v[::-1].copy() for k,v in a.items()}
        result=compare_photons(a,b)
        self.assertFalse(result['ordered_equal']);self.assertTrue(result['sorted_exact_equal'])
        b['energy_eV'][0]=np.nextafter(b['energy_eV'][0],np.float32(np.inf))
        self.assertFalse(compare_photons(a,b)['sorted_exact_equal'])
        b['x_mm']=b['x_mm'][:2]
        with self.assertRaises(ValueError):photon_records(b)


if __name__=='__main__':unittest.main()
