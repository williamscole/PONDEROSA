
from generate_files import GeneratePairs




if __name__ == "__main__":

    pairs = GeneratePairs(n_pairs=100, simple_segments=False)

    pairs.generate_some(["po", "fs", "hs", "gp", "av", "co"])

    pairs.write_out("del_now/test", n_chrom=22)

    import ipdb; ipdb.set_trace()