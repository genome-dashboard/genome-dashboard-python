#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `genome_dashboard` package."""


import unittest
from click.testing import CliRunner

from genome_dashboard import genome_dashboard
from genome_dashboard import cli


class TestGenome_dashboard(unittest.TestCase):
    """Tests for `genome_dashboard` package."""

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
        assert 'genome_dashboard.cli.main' in result.output
        help_result = runner.invoke(cli.main, ['--help'])
        assert help_result.exit_code == 0
        assert '--help  Show this message and exit.' in help_result.output
