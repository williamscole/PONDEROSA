import argparse
import pandas as pd
import numpy as np
import math
import os
import sys
import subprocess

from pedigree_tools import RemoveRelateds


def parse_populations(args):

    # load king file
    king_df = pd.read_csv(args.king, delim_whitespace=True, dtype = {"ID1": str, "ID2": str})

    # will hold subsets of the king file
    king_pops = []
    
    # keeps track of the n of individuals in each pop
    n_indv = []

    # keeps track of pop names
    pop_names = []

    # no pops
    if args.pops == "None":
        king_pops.append(king_df)
        n_indv.append(1)
        pop_names.append("pop1")

    else:
        pops_df = pd.read_csv(args.pops, delim_whitespace=True, header=None)

        # iterate through the populations
        for pop, pop_df in pops_df.groupby(1):

            # sample in proportion to the number of close reltives
            close_rels = king_df[king_df.InfType.isin(["4th", "3rd", "2nd", "PO", "FS"])]
            n_indv.append(len(set(close_rels["ID1"]) | set(close_rels["ID2"])))

            # add pop name
            pop_names.append(pop)

            ids = pop_df[0].values

            subset = king_df.apply(lambda x: x.ID1 in ids and x.ID2 in ids, axis=1)

            king_pops.append(king_df[subset])

            king_df = king_df[~subset]


    # n_populations long, each item a list of unrelated individuals from the population
    unrelateds_pops = []

    # iterate through the different king files
    for df in king_pops:
            
        # get unrelateds
        rr = RemoveRelateds()
        
        G = rr.king_graph(df, lambda PropIBD, a, b, c: PropIBD > args.max_pihat)

        unrelateds_pops.append(rr.get_unrelateds(G).unrelateds)

    # returns a list of unrelateds in each pop and the relative size of each population
    return unrelateds_pops, np.array(n_indv) / sum(n_indv), pop_names

def prep_simulations(args, unrelateds_pops, rel_pop_size, pop_names):

    for unrelateds, rel_size, pop in zip(unrelateds_pops, rel_pop_size, pop_names):

        # get the number of final samples we need
        n_samples = int(args.n_pairs * rel_size)

        # number of sims we're able to do per
        sims_per_gp = math.floor(len(unrelateds) / 4)
        sims_per_other = math.floor(len(unrelateds) / 5)

        # number of runs we have to do
        gp_sims_reqd = math.ceil(n_samples / sims_per_gp)
        other_sims_reqd = math.ceil(n_samples / sims_per_other)

        # write the def files
        def_files = [f"def {pop}_pgrandparent {sims_per_gp} 4\n1 1\n2 0 1 1sM\n3 1\n4 1",
                        f"def {pop}_mgrandparent {sims_per_gp} 4\n1 1\n2 0 1 1sF\n3 1\n4 1",
                        f"def {pop}_avuncular {sims_per_other} 4\n2 1 2 \n3 1 2\n4 1 1",
                        f"def {pop}_phs {sims_per_other} 3 M\n2 1 2 1:1 2:1\n3 1 2",
                        f"def {pop}_mhs {sims_per_other} 3 F\n2 1 2 1:1 2:1\n3 1 2"]
        
        # write the def files in the appropriate population dir
        for f, rel in zip(def_files, ["pgp", "mgp", "av", "phs", "mhs"]):
            out = open(f"{args.directory}/{pop}/{rel}.def", "w")
            out.write(f)
            out.close()
        
        # get a list of all the ped-sim runs to make
        runs = []

        # add all the runs
        runs += [f"{args.pedsim} {args.simmap} {args.intf} {pop} pgp {i} {args.ibdcaller}" for i in range(math.ceil(gp_sims_reqd/2))]
        runs += [f"{args.pedsim} {args.simmap} {args.intf} {pop} mgp {i} {args.ibdcaller}" for i in range(math.ceil(gp_sims_reqd/2))]
        runs += [f"{args.pedsim} {args.simmap} {args.intf} {pop} av {i} {args.ibdcaller}" for i in range(2*math.ceil(other_sims_reqd/2))]
        runs += [f"{args.pedsim} {args.simmap} {args.intf} {pop} phs {i} {args.ibdcaller}" for i in range(math.ceil(other_sims_reqd/2))]
        runs += [f"{args.pedsim} {args.simmap} {args.intf} {pop} mhs {i} {args.ibdcaller}" for i in range(math.ceil(other_sims_reqd/2))]

        ### write out the runs
        out = open(f"{args.directory}/{pop}/runs.txt", "w")
        out.write("\n".join(runs) + "\n")
        out.close()

        ### write out the unrelated samples
        out = open(f"{args.directory}/{pop}/unrelateds.txt", "w")
        out.write("\n".join(unrelateds) + "\n")
        out.close()

    # now write a bash script
    out = open(f"{args.directory}/simulate_pops.sh", "w")

    # start the loop and cd into the pop directory
    script = "for pop in " + " ".join(pop_names) + f"\ndo\n\ncd ${{pop}}\n\n"

    # subset the vcf
    script += f"bcftools view -S unrelateds.txt {args.vcf} -o ${{pop}}_unrelateds.vcf\n\n"

    # run a simple ped-sim; create the .def file
    script += 'echo "def temp 1 2" > temp.def\necho "1 1" >> temp.df\necho "2 1" >> temp.def\n\n'

    # create the command
    script += f"{args.pedsim} -i ${{pop}}_unrelateds.vcf -d temp.def -m {args.simmap} --intf {args.intf} -o temp\n"

    # get the .map file
    script += "grep -v '#' temp.vcf | awk '{print $1, $2, $3}' > variants.map\n\ncd ..\n"

    # create the new variants file
    script += f"python3 pedigree_tools.py {args.map} ${{pop}}/variants.map interpolate\n\nmv sim_chr*map ${{pop}}/\n\ncd ${{pop}}\n\n"

    # how run the simulations
    script += f"cat runs.txt | while read cmd\ndo\n\n"

    # run the sims
    cmd = "bash ../single_sim.sh $cmd"
    script += f"{args.slurm.replace('CMD', cmd)}\n\ndone\n\n"

    # run the IBD
    cmd = f"bash {args.ibdcaller}.sh"

    # go back to the prev directory
    script += "cd ..\n\ndone\n"

    out.write(script)
    out.close()
    
if __name__ == "__main__":

    parser = argparse.ArgumentParser("This code helps run ped-sim simulations!")

    parser.add_argument("--pops", type=str, default="None",
                        help="A file with no header where column 1 if the iid and column2 is the population. If not supplied, assumes only one population.")

    parser.add_argument("--king", type=str, help="KING .seg file", required=True)

    parser.add_argument("--max_pihat", type=float, default=0.10,
                        help="Maximum KING pihat value for both members of a pair to be included as a pedigree founder.")

    parser.add_argument("--n_pairs", default=100, type=int,
                        help="Minimum number of 2nd degree pairs for training.")

    parser.add_argument("--simmap", help="Full path and file name of the simmap.", required=True)

    parser.add_argument("--intf", help = "Full path and file name of the interference file", required=True)

    parser.add_argument("--pedsim", help = "Full path to ped-sim, e.g., /path/to/ped-sim/ped-sim")

    parser.add_argument("--slurm", help = "Slurm-like command to submit each simulation as its own job. Must have CMD in it where a command can be run. E.g., 'sbatch -t 1:00:00 --mem 32gb --wrap 'CMD''", default="CMD")

    parser.add_argument("--vcf", help = "Full path and file name of input vcf.", required=True)

    parser.add_argument("--map", help = "Full path to .map file in PLINK format for chromosome 1 or a file with all chroms. Must contain EXACTLY the coordinates in the .vcf", required=True)
    
    parser.add_argument("--directory", help = "Name of directory to store simulations in.", default = "ponderosa_simulations")
    
    parser.add_argument("--ibdcaller", help = "Prefix for .sh bash script to run IBD on vcfs", default = "phasedibd")

    args = parser.parse_args()

    unrelateds_pops, rel_pop_size, pop_names = parse_populations(args)

    directory = args.directory

    # create all necessary directories for running this
    if not os.path.exists(f"{directory}/"):
        # create the directories
        mkdir = subprocess.Popen(f"mkdir {directory}", shell=True, stdout=subprocess.PIPE)
        mkdir.wait()

        # create the directories for each population
        for pop in pop_names:
            mkdir = subprocess.Popen(f"mkdir {directory}/{pop}", shell=True, stdout=subprocess.PIPE)
            mkdir.wait()

        # move the ibdcaller.sh into the new dir
        if os.path.exists(f"{args.ibdcaller}.sh"):
            cp = subprocess.Popen(f"cp {args.ibdcaller}.sh {directory}/", shell=True, stdout=subprocess.PIPE)
            cp.wait()
        else:
            print("Cannot find [ibdcaller].sh")
            sys.exit()

        # make a dir for the segments
        mkdir = subprocess.Popen(f"mkdir {directory}/segments", shell=True, stdout=subprocess.PIPE)
        mkdir.wait()

        # move the python script into the new folder
        cp = subprocess.Popen(f"cp pedigree_tools.py {directory}/", shell=True, stdout=subprocess.PIPE)
        cp.wait()

        cp = subprocess.Popen(f"cp single_sim.sh {directory}/", shell=True, stdout=subprocess.PIPE)
        cp.wait()

    else:
        print("Directory already exists! Exiting...")
        sys.exit()

    prep_simulations(args, unrelateds_pops, rel_pop_size, pop_names)


