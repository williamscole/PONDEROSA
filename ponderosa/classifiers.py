from abc import ABC, abstractmethod
from typing import Dict, List, Tuple, Optional
import numpy as np
from sklearn.model_selection import LeaveOneOut
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
import networkx as nx
import pickle as pkl

from .pedigree import PedigreeRegistry, PedigreeHierarchy
from .prediction import MatrixHierarchy, process_phase_error
from .data_loading import Pairs


MIN_TRAINING_PAIRS = 5

def prepare_classifier_data(X: np.ndarray, *args):

    valid_mask = ~np.isnan(X).any(axis=1)

    out_arrs = [X[valid_mask]]

    for arg in args:
        assert X.shape[0] == arg.shape[0]

        out_arrs.append(arg[valid_mask])

    return tuple(out_arrs)


class RelationshipClassifier(ABC):

    def __init__(self, name: str):
        self.name = name
        self.model = None
        self.nodes = None
        self.features = None
        self.test_features = None

    def _validate_basic_input(self, X: np.ndarray, y: np.ndarray = None):
        if np.isnan(X).any():
            raise ValueError(
                f"{self.name} classifier expects clean data with no NaN values. "
                f"Use mask_nan() to clean your data first."
            )

        if not y is None:
            assert X.shape[0] == y.shape[0]

    def _validate_training_data(self, X: np.ndarray, y: np.ndarray, min_training_pairs=MIN_TRAINING_PAIRS) -> Tuple[np.ndarray, np.ndarray]:
        # TODO how to gracefully solve the issue of few training pairs

        self._validate_basic_input(X, y)

        assert np.unique(y).shape[0] > 1

        unique, counts = np.unique(y, return_counts=True)

        assert np.all(counts > min_training_pairs)

        return X, y

    def get_classes(self) -> List[str]:
        if self.model:
            return list(self.model.classes_)

    def data_from_pairs(self, pairs: Pairs, training: bool, registry: PedigreeRegistry = None):
        if training:
            pair_arr = registry.get_relatives_from(self.nodes)
            y = pair_arr[:,2]
            pair_list = pair_arr[:,:2]

            features = self.features

        else:
            y = np.array([])
            pair_list = np.array([])

            features = self.test_features

        X = pairs.get_pair_data_from(pair_list, *features, output_style="flatten")

        return prepare_classifier_data(X, y) if training else prepare_classifier_data(X)[0]

    def fit(self, X, y, loo: bool = False):

        if not loo: # Leave-one-out: X and y have already been properly processed
            self._validate_training_data(X, y)

            X, y = self._manipulate_training_data(X, y)

        self._train_model(X, y)

    def fit_predict_proba(self, X: np.ndarray, y: np.ndarray):

        self._validate_training_data(X, y, MIN_TRAINING_PAIRS + 1)

        X, y = self._manipulate_training_data(X, y)

        out_probas = []

        for train_idx, test_idx in self._get_loo_data(X, y):

            self.fit(X[train_idx], y[train_idx], loo=True)

            probas = self.predict_proba(X[test_idx], loo=True)

            out_probas.append(probas[0])

        return np.array(out_probas)

    def predict_proba(self, X, loo: bool = False):

        self._validate_basic_input(X)

        if not loo: # If loo, it has been processed
            X = self._manipulate_testing_data(X)

        probas = self._predict_proba(X)

        assert np.allclose(np.sum(probas, axis=1), 1.0)

        return probas

    def predict(self, X):

        self._validate_basic_input(X)

        X = self._manipulate_testing_data(X)

        return self._predict(X)

    def _get_loo_data(self, X: np.ndarray, y: np.ndarray):
        loo = LeaveOneOut()

        return loo.split(X)

    def _predict(self, X: np.ndarray):
        assert self.model is not None

        return self.model.predict(X)

    def _predict_proba(self, X: np.ndarray):
        assert self.model is not None

        return self.model.predict_proba(X)

    @abstractmethod
    def _manipulate_training_data(self, X, y):
        pass

    @abstractmethod
    def _manipulate_testing_data(self, X):
        pass

    @abstractmethod
    def _train_model(self, X, y):
        pass

class SegmentCountClassifier(RelationshipClassifier):
    """Classifier for number of IBD segments to distinguish 2nd degree relationships"""
    
    def __init__(self):
        super().__init__("n_segments")

        self.model = None
        self.nodes = ["AV", "PGP", "MGP", "PHS", "MHS"]
        self.features = ["N", "IBD1", "IBD2"]
        self.test_features = ["N", "IBD1", "IBD2"]

    def _manipulate_training_data(self, X: np.ndarray, y: np.ndarray):
        N = 0
        IBD1 = 1
        IBD2 = 2

        # The second feature is the IBD coverage def as IBD1 + IBD2
        X_train = np.column_stack([
            X[:, N],
            X[:, IBD1] + X[:, IBD2]
        ])

        return X_train, y

    def _manipulate_testing_data(self, X: np.ndarray):
        
        X, _ = self._manipulate_training_data(X, None)

        return X

    def _train_model(self, X: np.ndarray, y: np.ndarray):
        """
        Columns in X correspond to no. of segments, IBD1 prop., IBD2 prop.
        """
        
        self.model = LinearDiscriminantAnalysis()
        
        self.model.fit(X, y)


class HaplotypeScoreClassifier(RelationshipClassifier):

    def __init__(self):
        super().__init__("haplotype_score")

        self.model = None
        self.nodes = ["GPAV", "HS"]
        self.features = ["H1", "H2", "H1_ERR", "H2_ERR"]
        self.test_features = ["H1", "H2"]

    def _get_loo_data(self, X: np.ndarray, y: np.ndarray):
        test_indices = np.where(y != "Phase error")[0]
        all_indices = np.arange(y.shape[0])
        
        for test_idx in test_indices:
            train_idx = all_indices[all_indices != test_idx]
            yield train_idx, np.array([test_idx])

    def _manipulate_training_data(self, X, y):
        """
        Columns in X correspond to h1, h2, h1_error, h2_error
        """
        H1 = 0
        H2 = 1
        H1_ERR = 2
        H2_ERR = 3

        assert X.shape[1] == 4

        n_samples = X.shape[0]

        X_train = np.vstack([X[:,[H1,H2]], X[:,[H1_ERR,H2_ERR]]])

        phase_error_y = np.array(["Phase error"]*n_samples)

        y_train = np.concatenate([y, phase_error_y])

        return X_train, y_train

    def _manipulate_testing_data(self, X):

        return X[:,:2]

    def _train_model(self, X: np.ndarray, y: np.ndarray):

        self.model = LinearDiscriminantAnalysis()
        
        self.model.fit(X, y)


class DegreeClassifier(RelationshipClassifier):
    def __init__(self):
        super().__init__("degree")

        self.model = None
        self.nodes = ["PO", "FS", "2nd", "3rd", "4th"]
        self.features = ["IBD1", "IBD2"]
        self.test_features = ["IBD1", "IBD2"]

    def _manipulate_training_data(self, X, y):
        IBD1 = 0
        IBD2 = 1
        return X[:,[IBD1,IBD2]], y

    def _manipulate_testing_data(self, X):
        return X

    def _train_model(self, X, y):

        self.model = LinearDiscriminantAnalysis()

        self.model.fit(X, y)


CLASSIF_DICT = {
        "degree": DegreeClassifier,
        "no_segments": SegmentCountClassifier,
        "hap_score": HaplotypeScoreClassifier
    }

def train_classifiers(registry: PedigreeRegistry, pairs: Pairs) -> Dict[str, RelationshipClassifier]:

    trained_classifiers = {}

    for classif_name, classif in CLASSIF_DICT.items():

        classif = classif()

        X, y = classif.data_from_pairs(pairs=pairs, training=True, registry=registry)

        classif.fit(X, y)

        trained_classifiers[classif_name] = classif

    return trained_classifiers


def train_load_classifiers(registry: PedigreeRegistry, pairs: Pairs, training_file: str = None, prefix: str = None) -> Dict[str, RelationshipClassifier]:
    
    # Training file has been provided
    if training_file:
        with open(training_file, "rb") as pklf:

            trained_classifiers = pkl.load(pklf)
    else:
        trained_classifiers = train_classifiers(registry, pairs, prefix)

    return trained_classifiers


def run_inference(pairs: Pairs, trained_classifiers: Dict[str, RelationshipClassifier], hierarchy: PedigreeHierarchy) -> MatrixHierarchy:

    matrix_hier = MatrixHierarchy.from_hierarchy(hierarchy, pairs.get_pair_dict(index_to_pair=True), list(trained_classifiers.keys()))

    for classif_name, classif in trained_classifiers.items():

        X = classif.data_from_pairs(pairs=pairs, training=False)

        
        proba = classif.predict_proba(X)

        classes = classif.get_classes()

        if classif_name == "hap_score":
            proba = process_phase_error(proba)
            assert classes[-1] == "Phase error"
            classes.remove("Phase error")

        matrix_hier.add_probs(proba, classes, classif_name)

    matrix_hier.compute_probs()

    return matrix_hier


