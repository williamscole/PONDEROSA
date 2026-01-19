"""
PONDEROSA command-line entry point.

Allows running PONDEROSA as a module:
    python -m ponderosa --config config.yaml
    python -m ponderosa --ibd segments.txt --fam samples.fam --ibd-caller phasedibd
"""

from .cli import main

if __name__ == "__main__":
    main()
