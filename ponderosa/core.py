from .config import PonderosaConfig
from .data_loading import load_individuals, load_pairs
from .pedigree import build_pedigree, PedigreeHierarchy
from .classifiers import train_load_classifiers, run_inference
from .output import write_files

def run_ponderosa(config: PonderosaConfig):
    # Load the hierarchy
    hierarchy = PedigreeHierarchy.from_yaml()

    # Loads individual-level data
    individuals = load_individuals(config.files)

    # Loads pairwise data, including IBD information
    pairs = load_pairs(config.files, config.algorithm)

    # Finds all relationships in the dataset
    registry = build_pedigree(individuals, pairs, hierarchy)

    # Train/write or load classifiers
    classifiers = train_load_classifiers(registry,
                                         pairs,
                                         config.files.training)

    # Use the classifiers and compute posterior probability
    matrix_hierarchy = run_inference(pairs, classifiers, hierarchy)

    # Write out results
    files_written = write_files(pairs, registry, matrix_hierarchy, classifiers, config.output.output)



