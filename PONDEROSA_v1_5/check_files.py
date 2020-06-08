import os.path
import sys
import pandas as pd
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
import time
import datetime

def start_up(parameter_dict,file_dict,run_type):
    def logo():
        tree = ["\n           P O N D E R O S A",
        "                v.1.2          ",
        "                  &",                   
        "                 ##&",                   
        "                #####",                  
        "              &%#######",              
        "            &%&%#####%###",              
        "          &   &%######&   &",            
        "             &%########&",               
        "            %%%%#########",              
        "          %%%%%###########&",           
        "        %%&  &%#########  &&&",          
        "      *     %%###########.",             
        "          &%%##############",            
        "        &&&%%%############&",         
        "     #  &%%%%##/(###########&/",         
        "      %%%%%%##%((%%############%&",      
        "   &%%%%&&  &%%%&&&%#######&&&&&&&&&",   
        "&&.&%     &%%&   &((&&  ###&&&&&",        
        "                 %((&       ",      
        "                 %*(&",                  
        "                 #((&",                 
        "                .(((%",                  
        "                &%(((%",                 
        "               %%#(#((((",               
        "              &&&&(  &(&)",
        "                          ",
        "               (c) 2019",
        "C.M. Williams, C.R. Gignoux, B.M. Henn\n"]

        def print_color(color,text):
            if color == "bold":
                color = "\033[1m"
            if color == "brown":
                color = "\033[0;33m"
            elif color == "green":
                color = "\033[0;32m"
            colored_text = f"\033[{color}{text}\033[00m"
            return colored_text

        for index,i in enumerate(tree):
                if index in [k for k in range(19,26)]:
                    print(print_color("brown",i))
                elif index in [k for k in range(2,19)]:
                    print(print_color("green",i))
                elif index == 0:
                    print(print_color("bold",i))
                elif index in [1,27,28]:
                    print(i)

    def write_args():
        sys.stdout.write("Run type: %s\n\n" % run_type)
        sys.stdout.write("Arguments used:\n")
        for files in list(dict.fromkeys(file_dict)):
            file_name = file_dict[files]
            if file_name == "None":
                file_name = "None provided"
            sys.stdout.write("    %s: %s\n" % (files,file_name))
        for parameters in list(dict.fromkeys(parameter_dict)):
            parameter_val = parameter_dict[parameters]
            sys.stdout.write("    %s: %s\n" % (parameters,parameter_val))

    def locate_files():
        #Check for files
        file_type = list(dict.fromkeys(file_dict))
        not_found = [file_dict[file] for file in file_type if (not os.path.isfile(file_dict[file]) and file_dict[file] != "None")]
        if not_found != []:
            sys.stdout.write("\nWARNING!\nThe following files were not found:\n")
            for msg in not_found:
                sys.stdout.write("\t%s\n" % msg)
            sys.stdout.write("\nExiting PONDEROSA...\n\n")
            sys.exit()


    logo()
    write_args()
    sys.stdout.write("\nChecking files...done\n")
    locate_files()
    return {**parameter_dict,**file_dict}

class LogFile:
    def __init__(self,parameters,run_type):
        self.start_time = time.time()
        self.logfile = open("%s.log" % parameters["out"],"w")
        self.logfile.write("PONDEROSA v1.5\n\n%s\n" % datetime.datetime.now())
        self.logfile.write("\nRun type: %s\n\nParameters/files provided:\n" % run_type)        
        for pars in list(dict.fromkeys(parameters)):
            self.logfile.write("\t%s: %s\n" % (pars,parameters[pars]))

    def write(self,line):
        self.logfile.write(line)

    def validation(self,second_df):
        pass

    def mz_twins(self,twins):
        if twins != {}:
            self.write("\nMZ twins present. PONDEROSA collapses MZ twins into one individual.\n")
            self.write("The ID on the left has been replaced by the ID on the right.\n")
            for iid in twins:
                self.write("%s %s\n" % (iid,twins[iid]))

    def relative_pairs(self,rel_df):
        pass


    def write_log(self):
        time_elapsed = round(time.time() - time.time(),1)
        self.write("\nTime elapsed: %s\n\n" % time_elapsed)
        self.logfile.close()

    def write_errors(self,error_dict):
        if error_dict[2] != []:
            self.write("\nNon-critical errors detected. Please double check.\n")
            for errors in error_dict[2]:
                self.write(errors[0])
                for pairs in errors[1]:
                    self.write("\t"+" ".join(pairs) + "\n")
        for errors in error_dict[1]:
            self.write("Critical errors detected.\n")
            self.write(errors[0])
            self.write_log()
            sys.stdout.write("\nCritical error detected. See log file.\nExiting PONDEROSA...\n")
            sys.exit()






