import unittest
import numpy as np
from genomedashboard.ds import ds

class Testds(unittest.TestCase):
    
    def test_HP_intra(self):
        hp_intra = ds.HP_intra(-0.32, -0.58, -0.34, -3.80, 4.48, -0.73)
        self.assertEqual(hp_intra.she,-0.32)
        self.assertEqual(hp_intra.str,-0.58)
        self.assertEqual(hp_intra.sta,-0.34)
        self.assertEqual(hp_intra.buc,-3.80)
        self.assertEqual(hp_intra.pro,4.48)
        self.assertEqual(hp_intra.ope,-0.73)
        
    def test_HP_inter(self):
        hp_inter = ds.HP_inter(-0.87, -0.61, 3.07, 2.35, 0.78, 32.52)
        self.assertEqual(hp_inter.shi,-0.87)
        self.assertEqual(hp_inter.sli,-0.61)
        self.assertEqual(hp_inter.ris,3.07)
        self.assertEqual(hp_inter.til,2.35)
        self.assertEqual(hp_inter.rol,0.78)
        self.assertEqual(hp_inter.twi,32.52)
        
    def test_HP(self):
        hp_intra = ds.HP_intra(-0.32, -0.58, -0.34, -3.80, 4.48, -0.73)
        hp_inter = ds.HP_inter(-0.87, -0.61, 3.07, 2.35, 0.78, 32.52)
        hp = ds.HP(hp_intra, hp_inter,hptype='3DNA')
        self.assertIs(hp.HP_intra,hp_intra)
        self.assertIs(hp.HP_inter,hp_inter)
        self.assertEqual(hp.hptype,'3DNA')
        
    def test_RD(self):
        rd = ds.RD(np.zeros(3),np.eye(3))
        self.assertEqual(rd.r.tolist(), np.zeros(3).tolist())
        self.assertEqual(rd.d.tolist(),np.eye(3).tolist())
        
        
if __name__ == '__main__':
    unittest.main()
