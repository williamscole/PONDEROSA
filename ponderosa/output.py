import pickle as pkl
import pandas as pd

from .prediction import MatrixHierarchy
from .data_loading import Pairs
from .pedigree import PedigreeRegistry
from .evaluation import EvaluationResults


def write_out_readable(mhier: MatrixHierarchy, pairs: Pairs, registry: PedigreeRegistry, output_prefix: str):

    output_df = mhier.to_dataframe()

    # Temp remove the prob columns so we can append to end later
    prob_columns = output_df.attrs["prob_columns"]

    prob_df = output_df[prob_columns]
    output_df = output_df.drop(columns=prob_columns)

    # Add the known relationship
    output_df["known_rel"] = registry.get_relationships(output_df[["id1", "id2"]].values)

    # Add the pairwise data
    ibd_columns = ["IBD1", "IBD2", "N", "H1", "H2"]

    pair_data = pairs.get_pair_data_from(output_df[["id1", "id2"]].apply(tuple, axis=1).values,
                                         *ibd_columns, output_style="flatten")

    pair_data_df = pd.DataFrame(pair_data, columns=[i.lower() for i in ibd_columns])

    # Merge all the data together
    output_df = output_df.merge(pair_data_df, right_index=True, left_index=True)
    output_df = output_df.merge(prob_df, right_index=True, left_index=True)

    # Write out the data
    output_df.to_csv(f"{output_prefix}_pairs.txt", sep="\t", index=False)

    return f"{output_prefix}_pairs.txt"


def write_pickle(obj, output_file: str):
    with open(output_file, "wb") as f:
        pkl.dump(obj, f)

    return output_file

def write_out_matrix_hierarchy(mhier: MatrixHierarchy, output_prefix: str):

    return write_pickle(mhier, f"{output_prefix}.mhier.pkl")

def write_out_classifier(classifiers: dict, output_prefix: str):

    return write_pickle(classifiers, f"{output_prefix}.classif.pkl")

def write_evaluation_report(eval_results: EvaluationResults, output_prefix: str) -> str:
    """Write evaluation report to text file."""
    output_file = f"{output_prefix}_evaluation.txt"
    
    with open(output_file, 'w') as f:
        f.write(str(eval_results))
    
    return output_file

def write_files(pairs: Pairs, registry: PedigreeRegistry, mhier: MatrixHierarchy, classifiers: dict, output_prefix: str):

    files_written = []

    files_written.append(write_out_readable(mhier, pairs, registry, output_prefix))
    files_written.append(write_out_matrix_hierarchy(mhier, output_prefix))
    files_written.append(write_out_classifier(classifiers, output_prefix))

    return files_written



