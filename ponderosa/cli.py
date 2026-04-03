"""
Command-line interface for PONDEROSA relationship inference.

This module provides the CLI entry point for running PONDEROSA analyses.
It supports both configuration files (YAML) and direct command-line arguments,
with CLI arguments taking precedence over config file values.

Usage:
    # Using a config file
    python -m ponderosa --config config.yaml
    
    # Using command-line arguments
    python -m ponderosa --ibd segments.txt --fam samples.fam --ibd-caller phasedibd
    
    # Mixed: config file with CLI overrides
    python -m ponderosa --config config.yaml --output my_results --verbose
"""

import argparse
import sys
import logging
from pathlib import Path

from .config import PonderosaConfig

logger = logging.getLogger(__name__)


def create_parser() -> argparse.ArgumentParser:
    """
    Create the argument parser for PONDEROSA CLI.
    
    Returns:
        argparse.ArgumentParser: Configured argument parser
    """
    parser = argparse.ArgumentParser(
        prog="ponderosa",
        description="PONDEROSA: Relationship inference from IBD segments",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run with a configuration file
  python -m ponderosa --config analysis.yaml
  
  # Run with command-line arguments
  python -m ponderosa --ibd segments.ibd --fam samples.fam --ibd-caller phasedibd
  
  # Override config file settings
  python -m ponderosa --config config.yaml --output new_results -vv

For more information, visit: https://github.com/williamscole/PONDEROSA or email cole.williams@ucsf.edu
        """
    )
    
    # Config file option
    parser.add_argument(
        '--config', '-c',
        type=Path,
        help='YAML configuration file (CLI arguments override config values)'
    )
    
    # === File arguments ===
    files_group = parser.add_argument_group('Input Files')
    
    files_group.add_argument(
        '--ibd',
        type=Path,
        help='IBD segments file (required unless using --config)'
    )
    files_group.add_argument(
        '--fam',
        type=Path,
        help='PLINK FAM file with individual information (required unless using --config)'
    )
    files_group.add_argument(
        '--ibd-caller',
        choices=['phasedibd', 'hap-ibd', 'ibdseq'],
        default='phasedibd',
        dest='ibd_caller',
        help='IBD calling software used (default: phasedibd)'
    )
    files_group.add_argument(
        '--truth',
        type=Path,
        default=None,
        help='Ground truth relationships file for evaluation'
    )
    files_group.add_argument(
        '--map',
        type=Path,
        dest='mapf',
        help='Genetic map file (required for hap-ibd if segments are in bp)'
    )
    files_group.add_argument(
        '--ages',
        type=Path,
        help='File containing age information for individuals'
    )
    files_group.add_argument(
        '--priors',
        type=Path,
        help='File specifying relationship priors (e.g., age-based priors)'
    )
    files_group.add_argument(
        '--populations',
        type=Path,
        help='Population assignment file for population-specific analysis'
    )
    files_group.add_argument(
        '--training',
        type=Path,
        help='Pre-trained classifiers file (pickle) to skip training'
    )
    files_group.add_argument(
        '--rel-tree',
        type=Path,
        dest='rel_tree',
        help='Custom relationship hierarchy YAML file'
    )
    
    # === Algorithm arguments ===
    algo_group = parser.add_argument_group('Algorithm Parameters')
    
    algo_group.add_argument(
        '--min-segment-length',
        type=float,
        default=None,
        dest='min_segment_length',
        help='Minimum IBD segment length in cM (default: 3.0)'
    )
    algo_group.add_argument(
        '--min-total-ibd',
        type=float,
        default=None,
        dest='min_total_ibd',
        help='Minimum total IBD sharing in cM for a pair (default: 50.0)'
    )
    algo_group.add_argument(
        '--max-gap',
        type=float,
        default=None,
        dest='max_gap',
        help='Maximum gap in cM for stitching adjacent segments (default: 1.0)'
    )
    algo_group.add_argument(
        '--population',
        type=str,
        default=None,
        help='Population identifier for analysis (default: pop1)'
    )
    algo_group.add_argument(
        '--genome-length',
        type=float,
        default=None,
        dest='genome_length',
        help='Total genome length in cM (default: 3545.0)'
    )
    algo_group.add_argument(
        '--train-only',
        action='store_true',
        dest='train_only',
        help='Stop after training classifiers (skip inference)'
    )
    
    # === Output arguments ===
    output_group = parser.add_argument_group('Output Settings')
    
    output_group.add_argument(
        '--output', '-o',
        type=str,
        default=None,
        help='Output file prefix (default: ponderosa_results)'
    )
    output_group.add_argument(
        '--min-probability',
        type=float,
        default=None,
        dest='min_probability',
        help='Minimum probability threshold for reporting (default: 0.5)'
    )
    output_group.add_argument(
        '--write-training',
        action='store_true',
        dest='write_training',
        help='Save trained classifiers to pickle file'
    )
    output_group.add_argument(
        '--create-plots',
        action='store_true',
        dest='create_plots',
        help='Generate visualization plots'
    )
    output_group.add_argument(
        '--verbose', '-v',
        action='count',
        default=0,
        help='Increase verbosity (-v for INFO, -vv for DEBUG)'
    )
    output_group.add_argument(
        '--debug',
        action='store_true',
        help='Show full error tracebacks on failure'
    )
    
    return parser


def validate_args(args: argparse.Namespace) -> None:
    """
    Validate that required arguments are present.
    
    Args:
        args: Parsed command-line arguments
        
    Raises:
        SystemExit: If required arguments are missing
    """
    # If no config file, we need at minimum --ibd and --fam
    if not args.config:
        missing = []
        if not args.ibd:
            missing.append('--ibd')
        if not args.fam:
            missing.append('--fam')
        
        if missing:
            print(f"Error: The following arguments are required: {', '.join(missing)}", file=sys.stderr)
            print("  (Or provide a --config file with these values)", file=sys.stderr)
            sys.exit(1)


def main() -> None:
    """
    Main CLI entry point for PONDEROSA.
    
    Parses arguments, creates configuration, validates inputs,
    and runs the analysis pipeline.
    """
    parser = create_parser()
    args = parser.parse_args()
    
    # Validate arguments
    validate_args(args)
    
    try:
        # Convert argparse Namespace to dict, filtering out None values for non-flags
        args_dict = {}
        for key, value in vars(args).items():
            # Include the value if it's not None, or if it's a boolean flag that was set
            if value is not None:
                args_dict[key] = value
        
        # Create config from CLI args (and optional YAML file)
        config = PonderosaConfig.from_cli_and_yaml(args_dict)
        
        # Validate configuration
        config.validate()
        
        # Run the analysis
        from .core import run_ponderosa
        results = run_ponderosa(config)
        
        # Print summary if verbose
        if args.verbose > 0:
            print("\n" + results.summary())
        
        print(f"\nAnalysis complete. Results written to: {config.output.output}*")
        
    except FileNotFoundError as e:
        print(f"Error: File not found - {e}", file=sys.stderr)
        if args.debug:
            import traceback
            traceback.print_exc()
        sys.exit(1)
        
    except ValueError as e:
        print(f"Error: Invalid configuration - {e}", file=sys.stderr)
        if args.debug:
            import traceback
            traceback.print_exc()
        sys.exit(1)
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        if args.debug:
            import traceback
            print("\nDebug traceback:", file=sys.stderr)
            traceback.print_exc()
        else:
            print("(Run with --debug for full traceback)", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
