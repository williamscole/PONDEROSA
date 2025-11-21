


def split_regions(region_dict, new_region):
    # returns the overlap of 2 regions (<= 0 if no overlap)
    def overlap(region1, region2):
        start1, end1 = region1
        start2, end2 = region2
        return min(end1,end2) - max(start1,start2)
    # out region will be returned; is a dict of regions mapping to members of region    
    out_region = dict()
    # overlapped keeps track of all the regions that overlap with the new region
    overlapped = {tuple(new_region[:2]):[new_region[2]]}
    # iterate through the existing regions
    for region in sorted(region_dict):
        # if overlap
        if overlap(region, new_region[:2]) > 0:
            # the regions completely overlap, just add the member and return region dict
            if tuple(region) == tuple(new_region[:2]):
                region_dict[region] += [new_region[2]]
                return region_dict
            # bc the region overlaps, add it to overlapped
            overlapped[region] = region_dict[region]
        # no overlap, but add the region to the out_region dict
        else:
            out_region[region] = region_dict[region]
    # all the segments in overlapped overlap, so each consecutive pairs of coordinates in sites should/could have different members
    sites = sorted(set(it.chain(*overlapped)))
    # iterate thru consecutive sites
    for start, stop in zip(sites, sites[1:]):
        # get the members of the regions that overlap the consecutive sites
        info = [j for i, j in overlapped.items() if overlap((start, stop), i) > 0]
        # unpack the membership
        out_region[(start,stop)] = sorted(it.chain(*info))
    return out_region



# perform various computations on ibd segments for a pair of individuals
# takes as input a phasedibd segment data frame
class ProcessSegments:
    def __init__(self, pair_df):
        self.segs = pair_df

    '''Function takes as input a region_dict, which has the following format:
        {(start, stop): [obj1, obj2, obj3, ...]}
        Also takes as input a new_region which has format: [start, stop, obj]
        This function sees if there is overlap between the new region and existing regions
        If there is overlap, it splits the current region into new regions and concats the objs'''
    def split_regions(self, region_dict, new_region):
        # returns the overlap of 2 regions (<= 0 if no overlap)
        def overlap(region1, region2):
            start1, end1 = region1
            start2, end2 = region2
            return min(end1,end2) - max(start1,start2)

        # out region will be returned; is a dict of regions mapping to members of region    
        out_region = dict()
        # overlapped keeps track of all the regions that overlap with the new region
        overlapped = {tuple(new_region[:2]):[new_region[2]]}

        # iterate through the existing regions
        for region in sorted(region_dict):
            # if overlap
            if overlap(region, new_region[:2]) > 0:
                # the regions completely overlap, just add the member and return region dict
                if tuple(region) == tuple(new_region[:2]):
                    region_dict[region] += [new_region[2]]
                    return region_dict
                # bc the region overlaps, add it to overlapped
                overlapped[region] = region_dict[region]
            # no overlap, but add the region to the out_region dict
            else:
                out_region[region] = region_dict[region]
        
        # all the segments in overlapped overlap, so each consecutive pairs of coordinates in sites should/could have different members
        sites = sorted(set(it.chain(*overlapped)))
        # iterate thru consecutive sites
        for start, stop in zip(sites, sites[1:]):
            # get the members of the regions that overlap the consecutive sites
            info = [j for i, j in overlapped.items() if overlap((start, stop), i) > 0]
            # unpack the membership
            out_region[(start,stop)] = sorted(it.chain(*info))
        
        return out_region

    # stitches together segments that are at most max_gap apart
    def segment_stitcher(self, segment_list, max_gap = 1):
        regions = {}

        # iterate through the segments
        for start, stop in segment_list:

            '''Alg works by adding start/stop positions to overlapped
               We init overlapped with the start/stop of the segment
               Next we iterate through the current regions
               If there is overlap, we add the start/stop to overlapped
               At the end, we create a new region taking the min of overlapped and the max of overlapped'''
            overlapped = {start, stop}

            # we rewrite the regions
            updated_regions = set()

            # iterate through the current regions
            for r1, r2 in regions:

                # if there is a overlap with the ibd segment
                if min(stop, r2) - max(start, r1) > -max_gap:
                    # add the start/stop of the region
                    overlapped |= {r1, r2}
                # no overlap, so add the region to the updated regions
                else:
                    updated_regions |= {(r1, r2)}

            # add the new segment/new region to updated regions
            updated_regions |= {(min(overlapped), max(overlapped))}

            # for the next iteration
            regions = updated_regions.copy()

        # return the regions
        return regions

    # returns ibd1, ibd2 values for the pair
    def get_ibd1_ibd2(self):

        # init with ibd1 of 0 cM and ibd2 of 0 cM
        ibd1, ibd2 = 0, 0

        # iterate through the chromosomes
        for chrom, chrom_df in self.segs.groupby("chromosome"):

            # rdata frame
            r = {}

            # iterate through the segments
            for _, row in chrom_df.iterrows():
                # add the segments from the perspective of id1; name hap index as 0 --> 1 and 1 --> 2
                r = split_regions(r, [row["start_cm"], row["end_cm"], row["id1_haplotype"]+1])
                # add the segments from the perspective of id2; rename the haplotype index as 0 --> 3 and 1 --> 4
                r = split_regions(r, [row["start_cm"], row["end_cm"], row["id2_haplotype"]+3])

            # iterate through the regions
            for (start, end), hap in r.items():
                # get the length of the region
                l = end - start
                # r is covered on all 4 haplotypes --> IBD2
                if sum(set(hap)) == 10:
                    ibd2 += l
                # not ibd2
                else:
                    ibd1 += l
        
        return ibd1, ibd2

    # returns the number of IBD segments
    def get_n_segments(self):
        n = 0

        # run segment stitcher for each chromosome
        for _, chrom_df in self.segs.groupby("chromosome"):
            n += len(self.segment_stitcher(chrom_df[["start_cm", "end_cm"]].values))

        return n

    # returns the haplotype score of the pair
    def get_h_score(self):
        hap, tot = {0:0, 1:0}, 0

        # iterate through the chromosome
        for _, chrom_df in self.segs.groupby("chromosome"):
            r= {}
            # iterate through segments
            for _, row in chrom_df.iterrows():
                # on id1
                r = self.split_regions(r, [row["start_cm"], row["end_cm"], row["id1_haplotype"]+1])
                # on id2
                r = self.split_regions(r, [row["start_cm"], row["end_cm"], row["id2_haplotype"]+3])

            # holds the hap information for the chromosome
            temp = {1:0, 2:0, 3:0, 4:0}
            for (start, end), hapl in r.items():
                # present on 1+ haplotype for at least one in the pair
                if len(hapl) > 2:
                    continue
                # get length of region and add to total
                l = end - start
                tot += l

                # add the hap information
                for h in hapl:
                    temp[h] += l
            # for id1
            hap[0] += max(temp[1], temp[2])
            # for id2
            hap[1] += max(temp[3], temp[4])

        # return hap score
        h1, h2 = tot and hap[0]/tot or 0, tot and hap[1]/tot or 0
        return h1, h2

    def ponderosa_data(self, genome_len, empty=False):
        # creates an empty class
        class ponderosa: pass
        ponderosa = ponderosa()

        if empty:
            return ponderosa

        # add ibd1, ibd2 data
        ibd1, ibd2 = self.get_ibd1_ibd2()
        ponderosa.ibd1 = ibd1 / genome_len
        ponderosa.ibd2 = ibd2 / genome_len

        # get the number of ibd segments
        ponderosa.n = self.get_n_segments()

        # get haplotype scores
        h1, h2 = self.get_h_score()
        ponderosa.h1 = h1
        ponderosa.h2 = h2

        # return the data
        return ponderosa
    


# pair_df is a dataframe of a pair of relatives
# mean_d is the mean distance between switch errors
def introduce_phase_error(pair_df, mean_d):
    
    # given a mean distance between switch error returns a list of randomly drawn sites
    def generate_switches(mean_d, index):
        #start the switch at 0
        switches, last_switch = [], 0
        
        # longest chrom is 287 cM 
        while last_switch < 300:
            
            # add the new site to the previous site
            switches.append(np.random.exponential(mean_d) + last_switch)
            
            # previous site is now the new site
            last_switch = switches[-1]
            
        # return
        return [(i, index) for i in switches]
    
    # store the newly create segments
    new_segments = []
    
    for chrom, chrom_df in pair_df.groupby("chromosome"):
    
        # generate the switch locations
        s1 = generate_switches(mean_d, 0)
        s2 = generate_switches(mean_d, 1)
        switches = np.array(sorted(s1 + s2))
        switch_index, switches = switches[:,1], switches[:,0]

        # old segments
        segments = chrom_df[["start_cm", "end_cm", "id1_haplotype", "id2_haplotype", "id1", "id2", "chromosome"]].values

        # iterate through the segments
        for start, stop, hap1, hap2, id1, id2, chrom in segments:

            # get number of switches before the segment
            n_dict = {0: len(np.where(np.logical_and(switches<start, switch_index==0))[0]),
                    1: len(np.where(np.logical_and(switches<start, switch_index==1))[0])}


            # get the index of switches within the segment
            b = np.where(np.logical_and(switches>=start, switches<=stop))[0]

            # iterate through the switches and the switch index in the segment
            for s, index in zip(switches[b], switch_index[b]):
                # add the broken segment as the current start --> s and add n, which is the number of preceding switches
                new_segments.append([chrom, id1, id2, hap1, hap2, start, s, n_dict[0], n_dict[1]])

                # new start
                start = s

                # increase the number of switches by 1 but only on the relevant switch
                n_dict[index] += 1

            # add the final segment
            new_segments.append([chrom, id1, id2, hap1, hap2, start, stop, n_dict[0], n_dict[1]])

    pair_df = pd.DataFrame(new_segments, columns = ["chromosome", "id1", "id2", "id1_haplotype", "id2_haplotype", "start_cm", "end_cm", "n1", "n2"])
    pair_df["l"] = pair_df["end_cm"] - pair_df["start_cm"]
    pair_df = pair_df[pair_df.l >= 2.5]
    pair_df["id1_haplotype"] = pair_df[["id1_haplotype", "n1"]].apply(lambda x: x[0] if x[1]%2 == 0 else (x[0]+1)%2, axis = 1)
    pair_df["id2_haplotype"] = pair_df[["id2_haplotype", "n2"]].apply(lambda x: x[0] if x[1]%2 == 0 else (x[0]+1)%2, axis = 1)

    return pair_df