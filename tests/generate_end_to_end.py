"""
Synthetic data generation for end-to-end testing of PONDEROSA.

Generates IBD segments and fam files with known ground truth relationships
where feature distributions are clearly separable for 100% classification accuracy.
"""

import numpy as np
import itertools as it
from typing import Union
import pandas as pd


# Kinship coefficient bins for each relationship type
# These are non-overlapping (except FS/PO which are distinguished by IBD2)
KINSHIP_BINS = {
    "PO": (0.500, 0.500),
    "FS": (0.400, 0.600),
    "2nd": (0.180, 0.320),
    "3rd": (0.100, 0.150),
    "4th": (0.050, 0.070),
    "5th": (0.01, 0.03)
}

# Segment count bins for 2nd degree subtypes
# These distinguish relationships by number of IBD segments
SEGMENT_COUNT_BINS = {
    "PGP": (25, 30),
    "MGP": (35, 40),
    "PHS": (25, 50),
    "MHS": (50, 75),
    "AV": (50, 100),
    "PO": (22, 22),
    "FS": (50, 100),
    "3rd": (15, 25),
    "4th": (10, 20),
    "5th": (1, 8)
}


def get_kinship_coefficient(degree: Union[str, int]) -> float:
    """
    Generate a random kinship coefficient for an n-th degree relative.
    
    Parameters
    ----------
    degree : str or int
        Either "PO", "FS" for 1st degree relatives, or an integer 2-4 
        for 2nd through 4th degree relatives.
    
    Returns
    -------
    float
        A kinship coefficient sampled uniformly from the appropriate bin.
    
    Examples
    --------
    >>> get_kinship_coefficient("PO")  # Returns value in [0.475, 0.5]
    0.487
    >>> get_kinship_coefficient("FS")  # Returns value in [0.4, 0.6]
    0.523
    >>> get_kinship_coefficient(2)     # Returns value in [0.18, 0.32]
    0.251
    """
    
    if degree not in KINSHIP_BINS:
        raise ValueError(f"Unknown degree '{degree}'. Must be 'PO', 'FS', or integer 2-4.")
    
    low, high = KINSHIP_BINS[degree]
    return np.random.uniform(low, high)


def get_segment_count(relationship: str) -> int:
    """
    Generate a random segment count for a 2nd degree relationship subtype.
    
    Parameters
    ----------
    relationship : str
        One of "PGP", "MGP", "PHS", "MHS", or "AV".
    
    Returns
    -------
    int
        A segment count sampled uniformly from the appropriate bin.
    
    Examples
    --------
    >>> get_segment_count("PGP")  # Returns value in [15, 22]
    18
    >>> get_segment_count("MHS")  # Returns value in [50, 75]
    63
    """
    
    if relationship not in SEGMENT_COUNT_BINS:
        raise ValueError(f"Unknown relationship '{relationship}'. Must be one of {list(SEGMENT_COUNT_BINS.keys())}.")
    
    low, high = SEGMENT_COUNT_BINS[relationship]
    return np.random.randint(low, high + 1)


def generate_segment_coordinates(total_ibd: float, n_segments: int = None, genome_len: float = 3200.0, n_chromosomes: int = 22, min_seg_len: float = 5.0, min_gap: float = 5.0) -> "pd.DataFrame":
    """
    Generate random, non-overlapping segment coordinates along a genome.
    
    Segments are distributed across chromosomes and do not cross chromosome
    boundaries. The total length of all segments equals total_ibd.
    
    Parameters
    ----------
    total_ibd : float
        Total IBD to generate in centimorgans.
    n_segments : int, optional
        Number of segments to generate. If None, calculated automatically.
    genome_len : float, optional
        Total genome length in centimorgans. Default is 3200 cM.
    n_chromosomes : int, optional
        Number of chromosomes (equal length). Default is 22.
    min_seg_len : float, optional
        Minimum segment length in cM. Default is 5.0 cM.
    min_gap : float, optional
        Minimum gap between segments on the same chromosome. Default is 5.0 cM.
    
    Returns
    -------
    pd.DataFrame
        DataFrame with columns ['chromosome', 'start_cm', 'end_cm'], 
        sorted by chromosome and start position. Segments are non-overlapping
        and do not cross chromosome boundaries.
    """
    import pandas as pd
    
    if total_ibd <= 0:
        raise ValueError(f"total_ibd must be positive, got {total_ibd}")
    if total_ibd > genome_len:
        raise ValueError(f"total_ibd ({total_ibd}) cannot exceed genome_len ({genome_len})")
    
    chrom_len = genome_len / n_chromosomes
    
    # Special case: total_ibd == genome_len (e.g., PO relationship)
    # Return one segment per chromosome spanning the full length
    if abs(total_ibd - genome_len) < 1.0:
        segments = []
        for chrom in range(1, n_chromosomes + 1):
            segments.append({
                "chromosome": chrom,
                "start_cm": 0.0,
                "end_cm": chrom_len
            })
        return pd.DataFrame(segments)
    
    # Step 1: Determine number of segments
    if n_segments is None:
        max_n_segments = int(total_ibd / min_seg_len)
        n_segments = max(1, max_n_segments // 2)
    
    # Validate n_segments
    if n_segments < 1:
        raise ValueError(f"n_segments must be >= 1, got {n_segments}")
    if n_segments * min_seg_len > total_ibd:
        raise ValueError(f"Cannot fit {n_segments} segments of min {min_seg_len} cM in {total_ibd} cM total IBD")
    
    # Step 2: Generate segment lengths
    # Initialize all segments at min_seg_len (as integers to avoid floating point issues)
    len_arr = np.full(n_segments, int(min_seg_len), dtype=np.float64)
    
    # Calculate remaining IBD to distribute
    remaining = int(total_ibd - (min_seg_len * n_segments))
    
    # Distribute remaining IBD one cM at a time
    max_seg_len = chrom_len - min_gap  # Leave room for gaps
    
    for _ in range(remaining):
        # Try to find a segment that can accept more
        for attempt in range(n_segments * 2):
            idx = np.random.randint(0, n_segments)
            if len_arr[idx] < max_seg_len:
                len_arr[idx] += 1.0
                break
    
    # Sort segments largest to smallest for placement (ensures large ones get placed first)
    len_arr = np.sort(len_arr)[::-1]
    
    # Step 3: Place segments round-robin across chromosomes
    # Track used regions on each chromosome as list of (start, end) tuples
    chrom_regions = [[] for _ in range(n_chromosomes)]
    
    all_segments = []
    cur_chrom = 0  # 0-indexed
    
    for seg_len in len_arr:
        placed = False
        start_chrom = cur_chrom
        
        while True:
            # Find available gaps on this chromosome
            regions = sorted(chrom_regions[cur_chrom])
            
            # Build list of available gaps
            gaps = []
            prev_end = 0.0
            for (reg_start, reg_end) in regions:
                gap_size = reg_start - prev_end - min_gap
                if gap_size >= seg_len:
                    gaps.append((prev_end, reg_start - min_gap))
                prev_end = reg_end + min_gap
            
            # Check gap at the end of chromosome
            gap_size = chrom_len - prev_end
            if gap_size >= seg_len:
                gaps.append((prev_end, chrom_len))
            
            if gaps:
                # Randomly choose a gap and a position within it
                gap_start, gap_end = gaps[np.random.randint(len(gaps))]
                max_start = gap_end - seg_len
                start_pos = np.random.uniform(gap_start, max_start)
                
                # Round start to avoid floating point issues
                start_pos = round(start_pos, 6)
                
                # Add the segment - store the integer length directly
                all_segments.append({
                    "chromosome": cur_chrom + 1,  # 1-indexed for output
                    "start_cm": start_pos,
                    "end_cm": start_pos + seg_len,
                    "length_cm": seg_len  # Store exact integer length
                })
                chrom_regions[cur_chrom].append((start_pos, start_pos + seg_len))
                placed = True
            
            # Move to next chromosome (wrap around)
            cur_chrom = (cur_chrom + 1) % n_chromosomes
            
            # If placed or we've tried all chromosomes, break
            if placed or cur_chrom == start_chrom:
                break
        
        if not placed:
            raise RuntimeError(f"Could not place segment of length {seg_len:.2f} cM - genome too full")
    
    df = pd.DataFrame(all_segments)
    return df.sort_values(["chromosome", "start_cm"]).reset_index(drop=True)

def test_segment_generator():

    genome_len = 3200

    degrees = ["FS", "PO", "2nd", "3rd", "4th"]
    rels = ["PHS", "PGP", "MHS", "MGP", "AV"]

    for _ in range(1000):

        degree = np.random.choice(degrees)
        rel = np.random.choice(rels) if degree == "2nd" else degree

        total_ibd = 2 * get_kinship_coefficient(degree) * genome_len
        n_segs = get_segment_count(rel)

        if rel == "FS":
            total_ibd = total_ibd / 2

        while True:
            try:
                segments_df = generate_segment_coordinates(total_ibd, n_segs, genome_len, 22)
                break
            except:
                print("\tFailed. Trying again.")
                continue

        segments_df["l"] = segments_df.end_cm - segments_df.start_cm

        print(degree, rel, genome_len, total_ibd, n_segs, segments_df.l.sum(), segments_df.shape[0])

        assert segments_df[segments_df.l<5].shape[0] == 0
        assert (total_ibd - 1) <= segments_df.l.sum() <= (total_ibd + 1)
        assert segments_df.chromosome.nunique() <= 22
        assert segments_df.chromosome.max() <= 22

        assert segments_df.shape[0] == n_segs

        for _, chrom_df in segments_df.groupby("chromosome"):
            chrom_index = list(chrom_df.index)

            for idx1, idx2 in it.combinations(chrom_index, r=2):
                x1, x2 = chrom_df.loc[idx1][["start_cm","end_cm"]]
                y1, y2 = chrom_df.loc[idx2][["start_cm","end_cm"]]
                assert not (y1 < x1 < y2)
                assert not (x1 < y1 < x2)

    print("Tests passed!")

def generate_test_segments(degree: str, rel: str = None, genome_len: int = 3200, fixed1: int = 0, fixed2: int = 0):

    def get_segments(total_ibd, n_segs, genome_len, n_chromosomes):
        while True:
            try:
                df = generate_segment_coordinates(total_ibd, n_segs, genome_len, n_chromosomes)
                break
            except:
                continue
        return df

    kinship = get_kinship_coefficient(degree)
    total_ibd = 2 * kinship * genome_len

    n_segs = get_segment_count(degree) if rel is None else get_segment_count(rel)

    print(degree, rel, kinship, total_ibd, n_segs)

    if degree == "FS":
        total_ibd = total_ibd / 2

    if rel is None:
        if degree == "FS":
            df1 = get_segments(total_ibd, n_segs, genome_len, 22)
            df2 = get_segments(total_ibd, n_segs, genome_len, 22)
            df1["id1_haplotype"] = 0; df1["id2_haplotype"] = 0
            df2["id1_haplotype"] = 1; df2["id2_haplotype"] = 1
            segments_df = pd.concat([df1, df2])
        elif degree == "PO":
            segments_df = get_segments(total_ibd, n_segs, genome_len, 22)
            segments_df["id1_haplotype"] = 0
            segments_df["id2_haplotype"] = np.random.choice([0,1], segments_df.shape[0])
        else:
            segments_df = get_segments(total_ibd, n_segs, genome_len, 22)
            segments_df["id1_haplotype"] = fixed1; segments_df["id2_haplotype"] = fixed2
    else:
        if "HS" in rel:
            segments_df = get_segments(total_ibd, n_segs, genome_len, 22)
            segments_df["id1_haplotype"] = int(rel=="PHS"); segments_df["id2_haplotype"] = int(rel=="PHS")
        else:
            segments_df = get_segments(total_ibd, n_segs, genome_len, 22)
            segments_df["id1_haplotype"] = 0
            segments_df["id2_haplotype"] = np.random.choice([0,1], segments_df.shape[0])

    return segments_df
    

            


        



