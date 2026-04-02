"""
Core workflow orchestration for PONDEROSA relationship inference.

This module provides the main `run_ponderosa` function that coordinates
all analysis steps from data loading through relationship classification
and output generation.

Example:
    >>> from ponderosa.config import PonderosaConfig
    >>> from ponderosa.core import run_ponderosa
    >>> config = PonderosaConfig.from_yaml("config.yaml")
    >>> results = run_ponderosa(config)
"""

import logging
from dataclasses import dataclass
from typing import List, Optional, Dict, Any
from pathlib import Path

from .config import PonderosaConfig
from .data_loading import load_individuals, load_pairs, Individuals, Pairs
from .pedigree import build_pedigree, PedigreeHierarchy, PedigreeRegistry
from .classifiers import train_load_classifiers, run_inference, RelationshipClassifier
from .prediction import MatrixHierarchy
from .output import write_files, write_evaluation_report
from .evaluation import EvaluationResults

logger = logging.getLogger(__name__)


@dataclass
class PonderosaResults:
    """Container for PONDEROSA analysis results."""
    
    individuals: Individuals
    pairs: Pairs
    registry: PedigreeRegistry
    hierarchy: PedigreeHierarchy
    classifiers: Dict[str, RelationshipClassifier]
    matrix_hierarchy: MatrixHierarchy
    files_written: List[str]
    evaluation: EvaluationResults
    
    def summary(self) -> str:
        """Return a summary of the analysis results."""
        n_pairs = self.pairs.n_pairs()
        n_relationships = len([p for node_pairs in self.registry.nodes.values() for p in node_pairs])
        
        return (
            f"PONDEROSA Analysis Summary\n"
            f"==========================\n"
            f"Pairs analyzed: {n_pairs}\n"
            f"Known relationships found: {n_relationships}\n"
            f"Classifiers trained: {list(self.classifiers.keys())}\n"
            f"Files written: {len(self.files_written)}\n"
        )


def run_ponderosa(config: PonderosaConfig) -> PonderosaResults:
    """
    Execute the complete PONDEROSA analysis workflow.
    
    This function orchestrates the entire relationship inference pipeline:
    1. Load the relationship hierarchy definition
    2. Load individual-level data (fam file, ages)
    3. Load pairwise IBD data and compute features
    4. Build pedigree and identify known relationships
    5. Train or load classifiers
    6. Run inference to compute posterior probabilities
    7. Write output files
    
    Args:
        config: PonderosaConfig object containing all configuration settings
        
    Returns:
        PonderosaResults: Container with all analysis outputs
        
    Raises:
        FileNotFoundError: If required input files are missing
        ValueError: If configuration is invalid
        
    Example:
        >>> config = PonderosaConfig.from_yaml("analysis.yaml")
        >>> results = run_ponderosa(config)
        >>> print(results.summary())
    """
    
    # Setup logging based on verbosity
    _setup_logging(config.output.verbose)
    
    logger.info("Starting PONDEROSA analysis")
    logger.info(f"\n{config}")
    
    # Step 1: Load the relationship hierarchy
    logger.info("Loading relationship hierarchy...")
    if config.files.rel_tree:
        hierarchy = PedigreeHierarchy.from_yaml(config.files.rel_tree)
        logger.debug(f"Loaded custom hierarchy from {config.files.rel_tree}")
    else:
        hierarchy = PedigreeHierarchy.from_yaml()
        logger.debug("Using default relationship hierarchy")
    
    # Step 2: Load individual-level data
    logger.info("Loading individual data...")
    individuals = load_individuals(config.files)
    logger.info(f"Loaded {len(individuals.individuals)} individuals")
    
    # Step 3: Load pairwise IBD data
    logger.info("Loading IBD segments and computing pairwise features...")
    pairs = load_pairs(config.files, config.algorithm)
    logger.info(f"Computed features for {pairs.n_pairs()} pairs")
    
    # Step 4: Build pedigree and find known relationships
    logger.info("Building pedigree and identifying relationships...")
    registry = build_pedigree(individuals, pairs, hierarchy)
    _log_registry_summary(registry)
    
    # Step 5: Train or load classifiers
    if config.files.training:
        logger.info(f"Loading pre-trained classifiers from {config.files.training}")
    else:
        logger.info("Training classifiers on known relationships...")
    
    classifiers = train_load_classifiers(
        registry=registry,
        pairs=pairs,
        training_file=config.files.training,
        prefix=config.output.output if config.output.write_training else None
    )
    logger.info(f"Classifiers ready: {list(classifiers.keys())}")
    
    # Step 6: Run inference
    logger.info("Running relationship inference...")
    matrix_hierarchy = run_inference(pairs, classifiers, hierarchy)
    logger.info("Inference complete")

    if config.files.truth:
        logger.info(f"Evaluating predictions against ground truth: {config.files.truth}")
        try:
            evaluation = EvaluationResults.from_ponderosa(
            truth_file=config.files.truth,
            pred_mhier=matrix_hierarchy,
            hierarchy=hierarchy,
            level_names=["Degree", "Specific Relationship"]  # Optional
        )
            # TODO: fix this
            # logger.info(f"Evaluation complete - Accuracy: {evaluation.accuracy:.2%}")
            
            # Write evaluation report
            eval_file = write_evaluation_report(evaluation, config.output.output)
            logger.info(f"Evaluation report written: {eval_file}")
            
        except Exception as e:
            logger.warning(f"Failed to evaluate predictions: {e}")
            logger.debug("Continuing without evaluation...")

    # Step 7: Write output files
    logger.info(f"Writing output files with prefix: {config.output.output}")
    files_written = write_files(
        pairs=pairs,
        registry=registry,
        mhier=matrix_hierarchy,
        classifiers=classifiers,
        output_prefix=config.output.output
    )
    
    for f in files_written:
        logger.info(f"  Written: {f}")
    
    logger.info("PONDEROSA analysis complete!")
    
    # Return results container
    results = PonderosaResults(
        individuals=individuals,
        pairs=pairs,
        registry=registry,
        hierarchy=hierarchy,
        classifiers=classifiers,
        matrix_hierarchy=matrix_hierarchy,
        files_written=files_written,
        evaluation=evaluation
    )
    
    return results


def _setup_logging(verbosity: int) -> None:
    """Configure logging based on verbosity level."""
    if verbosity == 0:
        level = logging.WARNING
    elif verbosity == 1:
        level = logging.INFO
    else:
        level = logging.DEBUG
    
    logging.basicConfig(
        level=level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    )


def _log_registry_summary(registry: PedigreeRegistry) -> None:
    """Log a summary of relationships found in the registry."""
    total = 0
    for rel, pairs in registry.nodes.items():
        if pairs:
            count = len(pairs)
            total += count
            logger.debug(f"  {rel}: {count} pairs")
    
    logger.info(f"Found {total} known relationship pairs")
