import io
import unittest
import numpy as np
import pandas as pd
import pandas.testing as pdt

import biom
import skbio
import qiime2
from qiime2.plugin.testing import TestPluginBase
from qiime2 import Artifact

from q2_hill._methods import alpha_taxa, alpha_phylo, alpha_functional

class TestExamples(TestPluginBase):
    package = 'q2_hill.tests'

    def test_examples(self):
        self.execute_examples()

class AlphaTests(TestPluginBase):

    package = 'q2_hill.tests'

    def setUp(self):
        super().setUp()
        
        # Asignar métodos del plugin
        self.alpha_taxa = self.plugin.methods['alpha_taxa']
        self.alpha_phylo = self.plugin.methods['alpha_phylo']
        self.alpha_functional = self.plugin.methods.get('alpha_functional', None)

        # Importar datos de prueba
        self.empty_table = Artifact.import_data('FeatureTable[Frequency]', self.get_data_path('empty.biom'))

        self.two_feature_table = Artifact.import_data(
            'FeatureTable[Frequency]', 
            self.get_data_path('two_feature_table.biom')
        )

        self.three_feature_tree = Artifact.import_data(
            'Phylogeny[Rooted]', 
            self.get_data_path('three_feature.tree')
        )

        # Crear tabla BIOM de prueba
        biom_table = biom.Table(np.array([[0, 3, 4], [2, 3, 2]]),
                                ['C1', 'C2'],
                                ['S1', 'S2', 'S3'])
        self.test_table = Artifact.import_data('FeatureTable[Frequency]', biom_table)

        # Crear árbol filogenético de prueba
        tree_data = "((S1:0.25, S2:0.50):0.25, S3:0.75)root;"
        tree = skbio.TreeNode.read(io.StringIO(tree_data))
        self.test_tree = Artifact.import_data('Phylogeny[Rooted]', tree)

        # Crear datos funcionales en formato Metadata
        self.functional_data = qiime2.Metadata(pd.DataFrame(
                 [[0,3,4], [2,3,2], [3,2,4]], columns=['T1', 'T2', 'T3'],
                 index=pd.Index(['S1', 'S2', 'S3'], name='id')).T)

    # 📌 **Pruebas para diversidad taxonómica**
    def test_alpha_taxa(self):
        actual = self.alpha_taxa(table=self.test_table, q=1)[0].view(pd.Series)
        expected = pd.Series({'C1': 1.979627, 'C2': 2.941714}, name='q=1')  
        pdt.assert_series_equal(actual, expected)

    def test_alpha_taxa_invalid_q(self):
        with self.assertRaisesRegex(ValueError, "q must be a non-negative number"):
            self.alpha_taxa(table=self.test_table, q=-1)  # Validación de q

    def test_alpha_taxa_empty_table(self):
        with self.assertRaisesRegex(ValueError, "empty"):
            self.alpha_taxa(table=self.empty_table, q=1)

    # 📌 **Pruebas para diversidad filogenética**
    def test_alpha_phylo(self):
        actual = self.alpha_phylo(table=self.test_table, phylogeny=self.test_tree, q=1)[0].view(pd.Series)
        expected = pd.Series({'C1': 1.484720, 'C2': 1.641855}, name='q=1')  
        pdt.assert_series_equal(actual, expected)

    def test_alpha_phylo_invalid_q(self):
        with self.assertRaisesRegex(ValueError, "q must be a non-negative number"):
            self.alpha_phylo(table=self.test_table, phylogeny=self.test_tree, q=-1)  # Validación de q

    def test_alpha_phylo_empty_table(self):
        with self.assertRaisesRegex(ValueError, "empty"):
            self.alpha_phylo(table=self.empty_table, phylogeny=self.test_tree, q=1)

    # 📌 **Pruebas para diversidad funcional (si está implementado)**
    def test_alpha_functional(self):
        if self.alpha_functional:
            actual = self.alpha_functional(table=self.test_table, traits=self.functional_data, q=1, tau=0.8)[0].view(pd.Series)
            expected = pd.Series({'C1': 1.1256, 'C2': 1.157701}, name='q=1')  
            pdt.assert_series_equal(actual, expected)

    def test_alpha_functional_invalid_q(self):
        if self.alpha_functional:
            with self.assertRaisesRegex(ValueError, "q must be a non-negative number"):
                self.alpha_functional(table=self.test_table, traits=self.functional_data, q=-1, tau=0.8)  # Validación de q

    def test_alpha_functional_empty_table(self):
        if self.alpha_functional:
            with self.assertRaisesRegex(ValueError, "empty"):
                self.alpha_functional(table=self.empty_table, traits=self.functional_data, q=1, tau=0.8)

