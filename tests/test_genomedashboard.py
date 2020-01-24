#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `genomedashboard` package."""


import unittest
import numpy as np
import genomedashboard as gd
from click.testing import CliRunner


class TestGenomeDashboard(unittest.TestCase):
    """Tests for `genomedashboard` package."""

    def setUp(self):
        """Set up test fixtures, if any."""
        self.hp_intra = gd.HP_intra(-0.32, -0.58, -0.34, -3.80, 4.48, -0.73)
        self.hp_inter = gd.HP_inter(-0.87, -0.61, 3.07, 2.35, 0.78, 32.52)

    def tearDown(self):
        """Tear down test fixtures, if any."""
        pass

    def test_HP_intra(self):
        """Test something."""
        self.assertEqual(self.hp_intra.she, -0.32)
        self.assertEqual(self.hp_intra.str, -0.58)
        self.assertEqual(self.hp_intra.sta, -0.34)
        self.assertEqual(self.hp_intra.buc, -3.80)
        self.assertEqual(self.hp_intra.pro, 4.48)
        self.assertEqual(self.hp_intra.ope, -0.73)

    def test_HP_inter(self):
        self.assertEqual(self.hp_inter.shi, -0.87)
        self.assertEqual(self.hp_inter.sli, -0.61)
        self.assertEqual(self.hp_inter.ris, 3.07)
        self.assertEqual(self.hp_inter.til, 2.35)
        self.assertEqual(self.hp_inter.rol, 0.78)
        self.assertEqual(self.hp_inter.twi, 32.52)

    def test_HP(self):
        hp = gd.HP(self.hp_intra, self.hp_inter, hptype='3DNA')
        self.assertIs(hp.HP_intra, self.hp_intra)
        self.assertIs(hp.HP_inter, self.hp_inter)
        self.assertEqual(hp.hptype, '3DNA')

    def test_RD(self):
        rd = gd.RD(np.zeros(3), np.eye(3))
        self.assertEqual(rd.r.tolist(), np.zeros(3).tolist())
        self.assertEqual(rd.d.tolist(), np.eye(3).tolist())

    def test_Mask_3D(self):
        values = np.array([1.2, 2.7, 8.9])
        self.assertEqual(gd.Mask_3D(values).values.tolist(), values.tolist())
        self.assertEqual(gd.Mask_3D(values).dimension, 3)

    def test_SC(self):
        self.assertIs(gd.SC(HP=gd.HP(self.hp_intra, self.hp_inter)
                            ).HP.HP_intra, self.hp_intra)
        self.assertIs(gd.SC(somevalue=2).somevalue, 2)

    def test_SEQ(self):
        self.assertEqual(gd.SEQ('ACGT').tolist(), ['A-T', 'C-G', 'G-C', 'T-A'])
        self.assertEqual(gd.SEQ('ACGT').tostring(), 'ACGT')
        self.assertEqual(gd.SEQ('ACGT').tostep(), ['A-C', 'C-G', 'G-T'])

    def test_command_line_interface(self):
        """Test the CLI."""
        runner = CliRunner()
        result = runner.invoke(cli.main)
        assert result.exit_code == 0
        # assert 'genomedashboard.cli.main' in result.output.
        help_result = runner.invoke(cli.main, ['--help'])
        assert help_result.exit_code == 0
        assert '--help  Show this message and exit.' in help_result.output


if __name__ == '__main__':
    unittest.main()
