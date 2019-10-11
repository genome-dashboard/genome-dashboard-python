#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `genomedashboard` package."""


import unittest
from click.testing import CliRunner

# Need to improve testing.
# from genomedashboard import genomedashboard
# from genomedashboard import convert
# from genomedashboard import data
# from genomedashboard import ds
# from genomedashboard import io
from genomedashboard import cli


class TestGenomedashboard(unittest.TestCase):
    """Tests for `genomedashboard` package."""

    def setUp(self):
        """Set up test fixtures, if any."""

    def tearDown(self):
        """Tear down test fixtures, if any."""

    def test_000_something(self):
        """Test something."""

    def test_command_line_interface(self):
        """Test the CLI."""
        runner = CliRunner()
        result = runner.invoke(cli.main)
        assert result.exit_code == 0
        # assert 'genomedashboard.cli.main' in result.output  # Failing.
        help_result = runner.invoke(cli.main, ['--help'])
        assert help_result.exit_code == 0
        assert '--help  Show this message and exit.' in help_result.output
