import sys
from pathlib import Path
import numpy as np
import pytest
import pandas as pd
import polars as pl
from typing import Dict, List, Tuple, Optional


from ponderosa.classifiers import RelationshipClassifier, SegmentCountClassifier, HaplotypeScoreClassifier, DegreeClassifier, train_classifiers, run_inference
from ponderosa.data_loading import Pairs
from ponderosa.ibd_tools import Features
from ponderosa.pedigree import PedigreeHierarchy, PedigreeRegistry

FIRST_SECOND_REL_DICT = {
    'FS': {1: [[1, -1]], 2: [[1, -1]], 'sex': False, 'degree': '1st'},
    'PO': {1: [[1]], 'sex': False, 'degree': '1st'},
    'HS': {1: [[1, -1]], 'sex': False, 'degree': '2nd'},
    'AV': {1: [[1, 1, -1], [1, 1, -1]], 'sex': False, 'degree': '2nd', "parent": "GPAV"},
    'GP': {1: [[1, 1]], 'sex': False, 'degree': '2nd', "parent": "GPAV"},
    'MGP': {2: [[1, 1]], 'sex': True, 'degree': '2nd', 'parent': 'GP'},
    'PGP': {1: [[1, 1]], 'sex': True, 'degree': '2nd', 'parent': 'GP'},
    'MHS': {2: [[1, -1]], 'sex': True, 'degree': '2nd', 'parent': 'HS'},
    'PHS': {1: [[1, -1]], 'sex': True, 'degree': '2nd', 'parent': 'HS'},
    "GPAV": {"sex": False, "degree": "2nd"},
    "CO": {1: [[1, 1, -1, -1], [1, 1, -1, -1]], 'sex': False, 'degree': '3rd'},
    "HCO": {1: [[1, 1, -1, -1]], 'sex': False, 'degree': '4th'}
}

def generate_synthetic_data(labels, n_samples, n_features):
    """
    Generate synthetic data with different normal distributions for each label.
    
    Parameters:
    -----------
    labels : list
        List of labels/classes
    n_samples : int
        Number of samples to generate per label
    n_features : int
        Number of features per sample
    
    Returns:
    --------
    X : numpy.ndarray
        Feature matrix of shape (len(labels) * n_samples, n_features)
    y : numpy.ndarray
        Labels array of shape (len(labels) * n_samples,)
    """
    
    total_samples = len(labels) * n_samples
    X = np.zeros((total_samples, n_features))
    y = np.zeros(total_samples, dtype=object)
    
    for i, label in enumerate(labels):
        # Generate random means and variances for this label
        means = np.random.uniform(-2, 2, n_features)  # Random means between -2 and 2
        variances = np.random.uniform(0.5, 3, n_features)  # Random variances between 0.5 and 3
        
        # Calculate start and end indices for this label's samples
        start_idx = i * n_samples
        end_idx = (i + 1) * n_samples
        
        # Generate samples for each feature
        for j in range(n_features):
            X[start_idx:end_idx, j] = np.random.normal(means[j], np.sqrt(variances[j]), n_samples)
        
        # Fill in the labels
        y[start_idx:end_idx] = label
    
    return X, y

class TestDegreeClassifier:

    DEGREES = ["FS", "PO", "2nd", "3rd", "4th"]

    def test_basic(self):

        X, y = generate_synthetic_data(self.DEGREES, 10, 2)

        classif = DegreeClassifier()
        classif.fit(X, y)

        n_test = 5

        X_test, _ = generate_synthetic_data(self.DEGREES, n_test, 2)

        proba = classif.predict_proba(X_test)

        assert proba.shape[0] == len(self.DEGREES)*n_test
        assert proba.shape[1] == len(self.DEGREES)
        assert np.all(proba.sum(axis=1))

        predicted = classif.predict(X_test)

        assert predicted.shape[0] == len(self.DEGREES)*n_test
        assert len(set(predicted) - set(self.DEGREES)) == 0

    def test_loo(self):

        X, y = generate_synthetic_data(self.DEGREES, 10, 2)

        classif = DegreeClassifier()
        proba = classif.fit_predict_proba(X, y)

        assert proba.shape[0] == len(self.DEGREES)*10
        assert proba.shape[1] == len(self.DEGREES)
        assert np.all(proba.sum(axis=1))


class TestSegmentCountClassifier:

    SECOND = ["PHS", "MHS", "PGP", "MGP", "AV"]

    def test_training(self):

        classif = SegmentCountClassifier()

        X = np.array([
            [0, 0, 0],
            [1, 1, 1],
            [2, 2, 2]
        ])

        X_test = np.array([0, 2, 4])

        X_transformed, _ = classif._manipulate_training_data(X, None)

        assert X.shape[1] == 3
        assert X_transformed.shape[1] == 2
        assert np.array_equal(X_transformed[:, 1], X_test)

    def test_basic(self):

        X, y = generate_synthetic_data(self.SECOND, 10, 3)

        classif = SegmentCountClassifier()
        classif.fit(X, y)

        n_test = 5
        X_test, y = generate_synthetic_data(self.SECOND, 5, 3)

        proba = classif.predict_proba(X_test)

        assert proba.shape[0] == len(self.SECOND)*n_test
        assert proba.shape[1] == len(self.SECOND)
        assert np.all(proba.sum(axis=1)) 

        predicted = classif.predict(X_test)

        assert predicted.shape[0] == len(self.SECOND)*n_test
        assert len(set(predicted) - set(self.SECOND)) == 0

    def test_loo(self):

        X, y = generate_synthetic_data(self.SECOND, 10, 3)

        classif = SegmentCountClassifier()
        proba = classif.fit_predict_proba(X, y)

        assert proba.shape[0] == len(self.SECOND)*10
        assert proba.shape[1] == len(self.SECOND)
        assert np.all(proba.sum(axis=1))

class TestHaplotypeScoreClassifier:

    CAT = ["HS", "GPAV"]

    def test_training(self):

        X = np.array([
            [0, 1, 2, 3],
            [4, 5, 6, 7]
        ])

        y = np.array(["GPAV", "HS"])

        classif = HaplotypeScoreClassifier()

        X_train, y_train = classif._manipulate_training_data(X, y)

        assert y_train.shape[0] == y.shape[0] * 2
        assert X_train.shape[0] == y.shape[0] * 2
        assert np.array_equal(y_train, np.array(["GPAV", "HS", "Phase error", "Phase error"]))

        X_expected = np.array([
            [0, 1],
            [4, 5],
            [2, 3],
            [6, 7]
        ])

        assert np.array_equal(X_expected, X_train)

    def test_basic(self):

        X, y = generate_synthetic_data(self.CAT, 10, 4)

        classif = HaplotypeScoreClassifier()
        classif.fit(X, y)

        n_test = 5
        X_test, y = generate_synthetic_data(self.CAT, 5, 4)

        proba = classif.predict_proba(X_test)

        assert proba.shape[0] == len(self.CAT)*n_test
        assert proba.shape[1] == len(self.CAT) + 1
        assert np.all(proba.sum(axis=1)) 

        predicted = classif.predict(X_test)

        assert predicted.shape[0] == len(self.CAT)*n_test
        assert set(predicted) - set(self.CAT) == {"Phase error"}

    def test_loo(self):

        X, y = generate_synthetic_data(self.CAT, 10, 4)

        classif = HaplotypeScoreClassifier()
        proba = classif.fit_predict_proba(X, y)

        assert proba.shape[0] == len(self.CAT)*10
        assert proba.shape[1] == len(self.CAT) + 1
        assert np.all(proba.sum(axis=1))

class TestTraining:

    def test_basic(self):

        n_features = Features.NO_FEATURES

        nodes = ["PHS", "MHS", "PGP", "MGP", "AV", "FS", "PO", "CO", "HCO"]

        X, y = generate_synthetic_data(nodes, n_samples=10, n_features=n_features)

        segment_df = pd.DataFrame(np.abs(X))
        segment_df["id1"] = np.arange(0, segment_df.shape[0]).astype(str)
        segment_df["id2"] = np.arange(segment_df.shape[0], segment_df.shape[0]*2).astype(str)

        pairs = Pairs(pl.from_pandas(segment_df))

        hier = PedigreeHierarchy(FIRST_SECOND_REL_DICT)

        registry = PedigreeRegistry(hier)

        for rel, row in zip(y, segment_df.itertuples()):

            registry.add_pair(row.id1, row.id2, rel)

        
        classifiers = train_classifiers(registry, pairs)

        assert len(classifiers.keys()) == 3


class TestTesting:

    def _generate_data(self, hierarchy: PedigreeHierarchy) -> Tuple[PedigreeRegistry, Pairs]:

        n_features = Features.NO_FEATURES

        nodes = list(hierarchy.nodes)

        X, y = generate_synthetic_data(nodes, n_samples=15, n_features=n_features)

        segment_df = pd.DataFrame(np.abs(X))
        segment_df["id1"] = np.arange(0, segment_df.shape[0]).astype(str)
        segment_df["id2"] = np.arange(segment_df.shape[0], segment_df.shape[0]*2).astype(str)

        pairs = Pairs(pl.from_pandas(segment_df))

        registry = PedigreeRegistry(hierarchy)

        for rel, row in zip(y, segment_df.itertuples()):

            registry.add_pair(row.id1, row.id2, rel)

        return registry, pairs

    def test_basic(self):

        hierarchy = PedigreeHierarchy.from_yaml()

        train_registry, train_pairs = self._generate_data(hierarchy)

        classifiers = train_classifiers(train_registry, train_pairs)

        _, test_pairs = self._generate_data(hierarchy)

        mhier = run_inference(test_pairs, classifiers, hierarchy)

        # TODO: actually test that this works!

        # import pytest; pytest.set_trace()

