import sys
from pathlib import Path
import pytest


from ponderosa.config import PonderosaConfig

# config1 = PonderosaConfig.from_dict({"files": {"fam": "data/test1.fam",
#                                         "ages": "data/age1.txt",
#                                         "ibd": "data/segments1.txt",
#                                         "ibd_caller": "phasedibd"},
#                                         })


class TestYaml:

    def test_test1(self):

        TEST_DIR = Path(__file__).parent
        yaml_file_path = TEST_DIR / "data" / "test1"

        config1 = PonderosaConfig.from_yaml(yaml_file_path / "args.yaml")
        config1.validate()

        config2 = PonderosaConfig.from_yaml(yaml_file_path / "args_split_chrom.yaml")
        config2.validate()

        # pytest.set_trace()

