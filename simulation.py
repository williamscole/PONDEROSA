import argparse 
import pickle as pkl
from pedigree_tools import Simulations

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-vcf", help="Path and vcf file name.", type=str)
    parser.add_argument("-pop", default="pop1", help="Population to run simulations on. Defaults to running on a single population called pop.")
    parser.add_argument("-pop_file", help="Path and file name that contains ids (col1) and their associated population (col2). Not required; if not supplied, will run as a single population", default=None)
    parser.add_argument("-king", help="Path and file name of the KING file.")
    parser.add_argument("-k", help="Max PropIBD value of a mating pair or in-law", default=0.1)
    parser.add_argument("-pedsim", help="ped-sim executable", default="/users/cwilli50/ped-sim/ped-sim")
    parser.add_argument("-simmap", help="Sex-specific recombination map supplied by ped-sim.", default="/users/cwilli50/ped-sim/refined_mf.simmap")
    parser.add_argument("-intf", help="Path to interference file.", default="/users/cwilli50/ped-sim/interfere/nu_p_campbell.tsv")
    parser.add_argument("-def_dir", help="Path to the directory containing the def files, e.g., path/to/def/")
    parser.add_argument("-map")
    ### TODO: add path argument

    # different run types
    parser.add_argument("-init", help="Initializes the simulation.", action="store_true")
    parser.add_argument("-simulate", help="Creates the commands to simulate.", action="store_true")
    parser.add_argument("-phasedibd", action="store_true")
    parser.add_argument("-ibd", action="store_true")
    args = parser.parse_args()
    return args

bash_script = "while read cmd; do\n$cmd\ndone < iterITER.txt\n"



if __name__ == "__main__":

    args = parse_args()

    # initialize everything
    if args.init:
        sims = Simulations(args)

    else:
        i = open(f"{args.pop}_sims/sim_obj.pkl", "rb")
        sims = pkl.load(i)


    if args.simulate:

        sims.simulation_iter(args)

        i = open(f"{args.pop}_sims/run_pedsim.sh", "w")
        i.write(bash_script.replace("ITER", str(sims.iter)))
        i.close()

    elif args.phasedibd:

        sims.run_phasedibd(args)

    elif args.ibd:

        sims.analyze_ibd(args)



    with open(f"{args.pop}_sims/sim_obj.pkl", "wb") as pkl_file:
        pkl.dump(sims, pkl_file)
    



