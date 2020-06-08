import pandas as pd
import sys

def get_hap_score(relative_list,par_file,hap_file):
    pars = [par_file[par] for par in ["match_file","map_file","ped_file","out","num_chr","cm_gap","disc_homoz"]]
    match_file,map_file,ped_file,out,num_chr,cm_gap,disc_homoz = pars
    num_chr,cm_gap,disc_homoz = int(num_chr),float(cm_gap),int(disc_homoz)
    class GenotypeData:
        def __init__(self,map_file,ped_file,chrm):
            self.gts = dict()
            if ped_file != "None":
                for ind_gts in open(ped_file.replace("chr1","chr%s" % chrm)).readlines():
                    ind_gts = ind_gts.split()
                    iid,gts = ind_gts[0],ind_gts[6:]
                    gt_out = list()
                    for alleles in range(0,len(gts),2):
                        pair=(gts[alleles],gts[alleles+1])
                        gt_out.append(pair)
                    self.gts[iid] = gt_out

            self.snp_pos = list()
            self.mb_cm = dict()
            for snps in open(map_file.replace("chr1","chr%s" % chrm)).readlines():
                snps = snps.split()
                cm, mb = float(snps[2]),int(snps[3])
                self.snp_pos.append(cm)
                self.mb_cm[mb] = cm

        def mb_to_cm(self,mb):
            return self.mb_cm[mb]

        def gap_discordance(self,iid1,iid2,gap_start,gap_end,threshold):
            status = False
            if self.gts != {}:
                start=self.snp_pos.index(gap_start)+1
                end=self.snp_pos.index(gap_end)
                num_disc_homoz = 0
                for index in range(start,end):
                    gt1,gt2=self.gts[iid1][index],self.gts[iid2][index]
                    if (gt1[0]==gt1[1]) and (gt2[0]==gt2[1]) and (gt1!=gt2):
                        num_disc_homoz+=1
                    if num_disc_homoz > threshold:
                        status=True
                        break
            return status


    class PairData:
        def __init__(self,relative_list,out):
            self.relative_list,self.out = relative_list,out
            self.pair_data = dict()
            for pairs in relative_list:
                iid1,iid2 = pairs.split("_")
                #pair maps to {hap tots}, total ibd (for hap score)
                self.pair_data[pairs] = [{iid1:0,iid2:0},0,0]

        def finish_chrm(self,pair,iid1,iid2,seg_list,num_segs):
            hap_data,tot_ibd = {iid1:{"0":0,"1":0},iid2:{"0":0,"1":0}},0
            for segs in seg_list:
                seg_len = segs[4]-segs[3]
                tot_ibd += seg_len
                hap_data[segs[1]][segs[5]] += seg_len
                hap_data[segs[2]][segs[6]] += seg_len
            self.pair_data[pair][1] += tot_ibd
            self.pair_data[pair][2] += num_segs
            self.pair_data[pair][0][iid1] += max(hap_data[iid1]["0"],hap_data[iid1]["1"])
            self.pair_data[pair][0][iid2] += max(hap_data[iid2]["0"],hap_data[iid2]["1"])

        def get_scores(self,pair):
            info = self.pair_data[pair]
            iid1,iid2 = pair.split("_")
            h1,h2,total,num = info[0][iid1],info[0][iid2],info[1],info[2]
            if total != 0:
                h1,h2 = h1/total,h2/total
                ratio = min(h1/h2,h2/h1)
            else:
                h1,h2,ratio = 0,0,0
            return [pair,iid1,h1,iid2,h2,ratio,num]

        def get_relatives(self):
            return self.relative_list

        def write_out(self,relative_list,out):
            out_df = list()
            for pairs in self.relative_list:
                out_df.append(self.get_scores(pairs))
            out_df = pd.DataFrame(out_df,columns=["PAIR_ID","IID1","H1","IID2","H2","HSR","N"])
            with open("%s.txt" % self.out,"w") as outfile:
                outfile.write(out_df.to_string(index=False))
            outfile.close()
            return out_df

    if hap_file != "None":
        sys.stdout.write("\rSkipping hap score computation.\n")
        hap_df = open(hap_file).readlines()
        hap_df = [i.split() for i in hap_df]
        hap_df = pd.DataFrame(hap_df[1:],columns=hap_df[0])
        for col in ["H1","H2","HSR","N"]:
            hap_df[col] = hap_df[col].astype(float)
        return hap_df


    hap_data = PairData(relative_list,out)

    for chrm in range(1,num_chr+1):
        sys.stdout.write("\rCalculating hap score...chr %s" % chrm)
        sys.stdout.flush()
        genotype_data = GenotypeData(map_file,ped_file,chrm)

        def create_match_df(relative_list):
            match = pd.read_csv(match_file.replace("chr1","chr%s" % chrm),delim_whitespace=True,header=None)
            match["IID1"] = match[1].apply(lambda x: x[:-2])
            match["IID2"] = match[3].apply(lambda x: x[:-2])
            match["PAIR_ID"] = match[["IID1","IID2"]].min(axis=1) + "_" + match[["IID1","IID2"]].max(axis=1)
            match = match[match["PAIR_ID"].isin(relative_list)]
            match["HAP1"] = match[1].apply(lambda x: x[-1])
            match["HAP2"] = match[3].apply(lambda x: x[-1])
            match["CM_START"] = match[5].apply(genotype_data.mb_to_cm)
            match["CM_END"] = match[6].apply(genotype_data.mb_to_cm)
            match = match[["PAIR_ID","IID1","IID2","CM_START","CM_END","HAP1","HAP2"]].values.tolist()
            match.sort()
            return match

        match = create_match_df(hap_data.get_relatives())

        while match != []:
            pair_segs = list()
            pair,iid1,iid2 = match[0][:3]
            while pair in match[0]:
                pair_segs.append(match[0])
                match = match[1:]
                if match == []:
                    break
            start1,end1,index_list,num_segs = pair_segs[0][3],pair_segs[0][4],[0],0
            for i in range(1,len(pair_segs)):
                start2,end2 = pair_segs[i][3:5]
                if start1 < start2 < end1 and start1 < end2 < end1: #Completely overlap
                    continue
                index_list.append(i)
                if start1 == start2: #Segs have same start but seg2 is as long or longer
                    del index_list[-2]
                    continue
                if start2 > end1 + cm_gap or (start2 <= end1 + cm_gap and genotype_data.gap_discordance(iid1,iid2,end1,start2,disc_homoz)): #New seg
                    num_segs += 1
                    start1,end1 = start2,end2
                    continue
                end1 = end2
            num_segs += 1
            hap_data.finish_chrm(pair,iid1,iid2,[pair_segs[i] for i in index_list],num_segs)            

    sys.stdout.write("\rCalculating hap score...done   \n")
    return hap_data.write_out(relative_list,out)
