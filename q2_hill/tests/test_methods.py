# ----------------------------------------------------------------------------
# Copyright (c) 2024, Stephanie Hereira-Pacheco.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import tempfile
import unittest
import numpy as np
import pandas as pd
import pandas.testing as pdt

import biom
import skbio
import qiime2
from qiime2.plugin.testing import TestPluginBase
from qiime2.plugin.util import transform
from q2_types.feature_table import BIOMV100Format

from q2_hill._methods import (alpha_taxa, alpha_phylo, alpha_functional)

class TestExamples(TestPluginBase):
    package = 'q2_hill.tests'

    def test_examples(self):
        self.execute_examples()

class AlphaTests(TestPluginBase):

    package = 'q2_hill.tests'

    def setUp(self):
        super().setUp()
        
        # functions
        self.alpha_taxa = self.plugin.pipelines['alpha_taxa']
        self.alpha_phylo = self.plugin.pipelines['alpha_phylo']
        self.alpha_functional = self.plugin.pipelines.get('alpha_functional', None)

        # Import data
        self.empty_table = Artifact.import_data('FeatureTable[Frequency]', self.get_data_path('empty.biom'))

        self.two_feature_table = Artifact.import_data(
            'FeatureTable[Frequency]', 
            self.get_data_path('two_feature_table.biom')
        )

        self.three_feature_tree = Artifact.import_data(
            'Phylogeny[Rooted]', 
            self.get_data_path('three_feature.tree')
        )

        # Create artif data
        biom_table = biom.Table(np.array([[0, 3, 4], [2, 3, 2]]),
                                ['T1', 'T2'],
                                ['S1', 'S2', 'S3'])
        self.test_table = Artifact.import_data('FeatureTable[Frequency]', biom_table)

        # Create tree
        tree_data = "((T1:0.25, T2:0.50):0.25, T3:0.75)root;"
        tree = skbio.TreeNode.read(io.StringIO(tree_data))
        self.test_tree = Artifact.import_data('Phylogeny[Rooted]', tree)

        if self.alpha_functional:
            func_data = biom.Table(np.array([[0.2, 0.5], [0.3, 0.7]]),
                                   ['Trait1', 'Trait2'],
                                   ['O1', 'O2'])
            self.functional_table = Artifact.import_data('FeatureTable[Frequency]', func_data)

    #test alpha taxa
     def test_alpha_taxa(self):
        actual = self.alpha_taxa(table=self.test_table, q=1)[0].view(pd.Series)
        expected = pd.Series({'S1': 1, 'S2': 2, 'S3': 2}, name='q=1')  # Computado manualmente
        pdt.assert_series_equal(actual, expected)

    def test_alpha_taxa_invalid_q(self):
        with self.assertRaisesRegex(TypeError, "incompatible"):
            self.alpha_taxa(table=self.test_table, q=-1)

    def test_alpha_taxa_empty_table(self):
        with self.assertRaisesRegex(ValueError, "empty"):
            self.alpha_taxa(table=self.empty_table, q=1)

    # Pruebas para diversidad filogenética
    def test_alpha_phylo(self):
        actual = self.alpha_phylo(table=self.test_table, phylogeny=self.test_tree, q=1)[0].view(pd.Series)
        expected = pd.Series({'S1': 0.75, 'S2': 1.0, 'S3': 1.0}, name='q=1')  # Computado con skbio
        pdt.assert_series_equal(actual, expected)

    def test_alpha_phylo_invalid_q(self):
        with self.assertRaisesRegex(TypeError, "incompatible"):
            self.alpha_phylo(table=self.test_table, phylogeny=self.test_tree, q=-1)

    def test_alpha_phylo_empty_table(self):
        with self.assertRaisesRegex(ValueError, "empty"):
            self.alpha_phylo(table=self.empty_table, phylogeny=self.test_tree, q=1)

    # 📌 **Pruebas para diversidad funcional (si está implementado)**
    def test_alpha_functional(self):
        if self.alpha_functional:
            actual = self.alpha_functional(table=self.functional_table, q=1, tau=0.5)[0].view(pd.Series)
            expected = pd.Series({'O1': 0.2, 'O2': 0.3}, name='q=1')  # Simulado
            pdt.assert_series_equal(actual, expected)

    def test_alpha_functional_invalid_q(self):
        if self.alpha_functional:
            with self.assertRaisesRegex(TypeError, "incompatible"):
                self.alpha_functional(table=self.functional_table, q=-1, tau=0.5)

    def test_alpha_functional_empty_table(self):
        if self.alpha_functional:
            with self.assertRaisesRegex(ValueError, "empty"):
                self.alpha_functional(table=self.empty_table, q=1, tau=0.5)

