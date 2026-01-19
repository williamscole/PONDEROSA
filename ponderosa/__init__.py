"""
PONDEROSA: Relationship inference from IBD segments.

PONDEROSA is a tool for inferring genetic relationships between individuals
using Identity By Descent (IBD) segment data. It uses machine learning
classifiers trained on IBD sharing patterns to distinguish between different
degrees of biological relationships.

Quick Start:
    >>> from ponderosa import PonderosaConfig, run_ponderosa
    >>> config = PonderosaConfig.from_yaml("config.yaml")
    >>> results = run_ponderosa(config)
    >>> print(results.summary())

Command-line usage:
    python -m ponderosa --config config.yaml
    python -m ponderosa --ibd segments.txt --fam samples.fam --ibd-caller phasedibd

For more information, see the documentation at:
    https://github.com/williamscole/PONDEROSA
"""

__version__ = "0.2.0"
__author__ = "Cole Williams"

# Core API
from .config import PonderosaConfig, FilesConfig, AlgorithmConfig, OutputConfig
from .core import run_ponderosa, PonderosaResults

# Data loading
from .data_loading import (
    load_individuals,
    load_pairs,
    Individuals,
    Pairs,
    IBDLoader,
    PhasedIBDLoader,
    HapIBDLoader,
)

# Pedigree analysis
from .pedigree import (
    PedigreeHierarchy,
    PedigreeRegistry,
    build_pedigree,
)

# Classification
from .classifiers import (
    RelationshipClassifier,
    DegreeClassifier,
    SegmentCountClassifier,
    HaplotypeScoreClassifier,
    train_classifiers,
    run_inference,
)

# Prediction
from .prediction import MatrixHierarchy

# Output
from .output import write_files

__all__ = [
    # Version info
    "__version__",
    "__author__",
    # Core API
    "run_ponderosa",
    "PonderosaResults",
    "PonderosaConfig",
    "FilesConfig",
    "AlgorithmConfig",
    "OutputConfig",
    # Data loading
    "load_individuals",
    "load_pairs",
    "Individuals",
    "Pairs",
    "IBDLoader",
    "PhasedIBDLoader",
    "HapIBDLoader",
    # Pedigree
    "PedigreeHierarchy",
    "PedigreeRegistry",
    "build_pedigree",
    # Classification
    "RelationshipClassifier",
    "DegreeClassifier",
    "SegmentCountClassifier",
    "HaplotypeScoreClassifier",
    "train_classifiers",
    "run_inference",
    # Prediction
    "MatrixHierarchy",
    # Output
    "write_files",
]
