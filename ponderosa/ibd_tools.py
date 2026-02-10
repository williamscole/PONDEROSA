import numpy as np
import polars as pl
import itertools as it
from typing import List, Tuple, Any
from scipy.stats import multivariate_normal

####### IBD segment methods #######
class IBD(pl.DataFrame):
    """Polars DataFrame subclass for IBD segments with validation."""
    
    _schema = {
        "id1": pl.Utf8,
        "id2": pl.Utf8,
        "id1_haplotype": pl.Int8,
        "id2_haplotype": pl.Int8,
        "chromosome": pl.Int8,
        "start_cm": pl.Float64,
        "end_cm": pl.Float64,
    }
    
    def __init__(self, data=None, validate: bool = True, **kwargs):
        
        if data is None:
            # Create empty DataFrame with correct schema
            data = pl.DataFrame(schema=self._schema)
        
        # Handle different input types
        if isinstance(data, pl.LazyFrame):
            data = data.collect()
        elif isinstance(data, dict):
            data = pl.DataFrame(data)
        elif isinstance(data, IBD):
            validate = False  # Already validated
        
        if validate and not isinstance(data, IBD):
            data = self._validate_and_clean(data)
        
        # Call parent constructor
        super().__init__(data, **kwargs)
    
    def _validate_and_clean(self, data: pl.DataFrame) -> pl.DataFrame:
        """Validate schema and clean data."""
        
        # Check for missing columns
        expected_cols = list(self._schema.keys())
        actual_cols = data.columns
        missing_cols = [col for col in expected_cols if col not in actual_cols]
        
        if missing_cols:
            raise ValueError(f"Missing required columns: {missing_cols}")
        
        # Remove 'l' column if it exists (we ignore it)
        if 'l' in actual_cols:
            data = data.drop('l')
        
        # Select and reorder to expected columns only
        data = data.select(expected_cols)
        
        # Cast to correct types
        data = data.cast(self._schema)
        
        # Validate haplotypes
        self._validate_haplotypes(data)
        
        return data
    
    def _validate_haplotypes(self, data: pl.DataFrame):
        """Validate haplotype values are 0 or 1."""
        for col in ["id1_haplotype", "id2_haplotype"]:
            invalid_mask = ~data[col].is_in([0, 1])
            if invalid_mask.any():
                invalid_values = data.filter(invalid_mask)[col].unique().to_list()
                raise ValueError(
                    f"Column '{col}' contains invalid haplotype values: {invalid_values}. "
                    f"Values must be 0 or 1."
                )
    
    @property
    def _constructor(self):
        """Return constructor for operations that return new DataFrames."""
        def _c(*args, **kwargs):
            if 'validate' not in kwargs:
                kwargs['validate'] = False
            return IBD(*args, **kwargs)
        return _c

class Features:
    """Centralized feature indexing system"""
    IBD1 = 0
    IBD2 = 1
    H1 = 2
    H2 = 3
    HTOT = 4
    N = 5
    N2 = 6
    H1_ERR = 7
    H2_ERR = 8
    HTOT_ERR = 9
    NO_FEATURES = 10

    SWAPPABLE = {
        "H1": "H2",
        "H2": "H1",
        "H1_ERR": "H2_ERR",
        "H2_ERR": "H1_ERR"
    }
    
    # Alternative: use a class method to get feature names
    @classmethod
    def get_feature_names(cls):
        return ["IBD1", "IBD2", "H1", "H2", "HTOT", "N", "N2", "H1_ERR", "H2_ERR", "HTOT_ERR"]

    @classmethod
    def get_index(cls, feature_name: str) -> int:
        """Get index for a feature by name"""
        return getattr(cls, feature_name.upper())

    @classmethod
    def get_swappable(cls, feature_name: str) -> str:
        return cls.SWAPPABLE.get(feature_name, feature_name)



# Processes IBD segments for a pair of individuals
class ProcessSegments:

    CHROMOSOME = 0
    START = 1
    END = 2
    ID1HAP = 3
    ID2HAP = 4

    def __init__(self, segments_arr: np.array, max_gap: int = 1.0):
        """
        Initialize with an nx5 numpy array where columns are:
        [chromosome, start_cm, end_cm, id1_haplotype, id2_haplotype]
        """

        self.segs = segments_arr.copy()

        self.max_gap = max_gap

        self.segs[:,self.ID1HAP] += 1 # Unique haplotype index for id1
        self.segs[:,self.ID2HAP] += 3 # And for id2

    def split_regions(self, region_dict: dict, new_region: Tuple[float, float, Any]) -> dict:
        """
        Takes as input:
        - region_dict: {(start, stop): [obj1, obj2, obj3, ...]}
        - new_region: [start, stop, obj]
        Returns updated region dict that splits overlaps into new segments
        """

        def overlap(region1: Tuple[float, float], region2: Tuple[float, float]) -> float:
            start1, end1 = region1
            start2, end2 = region2
            return min(end1, end2) - max(start1, start2)

        out_region = dict()
        overlapped = {tuple(new_region[:2]): [new_region[2]]}

        for region in sorted(region_dict):
            if overlap(region, new_region[:2]) > 0:
                if tuple(region) == tuple(new_region[:2]):
                    region_dict[region] += [new_region[2]]
                    return region_dict
                overlapped[region] = region_dict[region]
            else:
                out_region[region] = region_dict[region]
        
        sites = sorted(set(it.chain(*overlapped)))
        for start, stop in zip(sites, sites[1:]):
            info = [j for i, j in overlapped.items() if overlap((start, stop), i) > 0]
            out_region[(start, stop)] = sorted(it.chain(*info))
        
        return out_region

    def compute_ibd1_ibd2(self, regions: dict) -> np.array:

        ibd1, ibd2, n2 = 0, 0, 0

        for (start, end), hap in regions.items():
            l = end - start
            if sum(set(hap)) == 10:  # Covered on all 4 haplotypes (1+2+3+4=10)
                ibd2 += l
                n2 += 1
            else:
                ibd1 += l
        
        return np.array([ibd1, ibd2, n2])

    def compute_h(self, regions: dict, inter_phase: bool = False) -> np.array:
        tot = 0
        temp = {1: 0, 2: 0, 3: 0, 4: 0}

        for (start, end), hapl in regions.items():
            if len(hapl) > 2:  # Present on 1+ haplotype for at least one in pair
                continue
                
            l = end - start
            tot += l

            for h in hapl:
                temp[h] += l

        t1 = temp[1] if inter_phase else max(temp[1], temp[2])
        t2 = temp[3] if inter_phase else max(temp[3], temp[4])

        return np.array([t1, t2, tot])

    def _phase_error_h_score(self):
            mean = [0.65, 0.65]
            cov = [[0.004, 0], [0, 0.004]]
            h = multivariate_normal(mean, cov).rvs(1)
            return np.concatenate((h, np.array([1])))

    def _get_ibd_stats(self, chrom_index: int) -> np.array:
        ibd1, ibd2 = 0, 0
        region1 = dict()

        region2 = set()

        for seg in self.segs[chrom_index, :]:

            ####### for computing ibd1, ibd2 and haplotype score ######
            region1 = self.split_regions(region1, seg[[self.START,self.END,self.ID1HAP]])
            region1 = self.split_regions(region1, seg[[self.START,self.END,self.ID2HAP]])

            ###### For computing n of segments ########
            start, stop = seg[[self.START,self.END]]
            overlapped = {start, stop}
            updated_regions = set()

            for r1, r2 in region2:
                if min(stop, r2) - max(start, r1) > -self.max_gap:
                    overlapped |= {r1, r2}
                else:
                    updated_regions |= {(r1, r2)}

            updated_regions |= {(min(overlapped), max(overlapped))}
            region2 = updated_regions.copy()

        ibd_stats = np.zeros(Features.NO_FEATURES)
        ibd_stats[Features.N] = len(region2) # No. of IBD segments
        ibd_stats[[Features.IBD1, Features.IBD2, Features.N2]] += self.compute_ibd1_ibd2(region1)
        ibd_stats[[Features.H1, Features.H2, Features.HTOT]] += self.compute_h(region1)
        ibd_stats[[Features.H1_ERR, Features.H2_ERR, Features.HTOT_ERR]] += self._phase_error_h_score()

        return ibd_stats

    @classmethod
    def ibd_sharing_stats(cls, pair_df: pl.DataFrame, max_gap: float = 1) -> np.array:

        segments = pair_df.select(["chromosome", "start_cm", "end_cm", "id1_haplotype", "id2_haplotype"]).to_numpy()

        instance = cls(segments, max_gap)

        ibd_stats = np.zeros(Features.NO_FEATURES)

        for chrom in np.unique(segments[:,0]):

            chrom_index = np.where(segments[:,0]==chrom)[0]

            ibd_stats += instance._get_ibd_stats(chrom_index)

        ibd_stats[[Features.H1, Features.H2]] /= ibd_stats[Features.HTOT]
        ibd_stats[[Features.H1_ERR, Features.H2_ERR]] /= ibd_stats[Features.HTOT_ERR]

        return ibd_stats