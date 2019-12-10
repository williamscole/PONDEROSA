import os.path
import sys

def start_up(parameter_dict):

    tree = ["\n       P  O  N  D  E  R  O  S  A",
    "                v.1.1          ",
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

    sys.stdout.write("Arguments used:\n")

    for files in ["fam file","match file","king file","map file","ped file","age file","haps file"]:
        file_name = parameter_dict[files]
        if file_name == "":
            file_name = "None provided"
        sys.stdout.write("    %s: %s\n" % (files,file_name))

    for parameters in ["out","mhs gap","gp gap","po gap","ilash","cm gap","max discordant homoz.","likelihood"]:
        parameter_val = parameter_dict[parameters]
        sys.stdout.write("    %s: %s\n" % (parameters,parameter_val))

    stop = False
    #Check for files
    not_found = list()
    for files in ["fam file","king file","age file","haps file"]:
        files = parameter_dict[files]
        if files == "":
            continue
        if not os.path.isfile(files):
            not_found.append(files)
    for i in range(1,23):
        if not os.path.isfile(parameter_dict["map file"] % str(i)):
            not_found.append(parameter_dict["map file"] % str(i))
        if not os.path.isfile(parameter_dict["match file"] % str(i)):
            not_found.append(parameter_dict["match file"] % str(i))
        if parameter_dict["ped file"] != "" and not os.path.isfile(parameter_dict["ped file"] % str(i)):
            not_found.append(parameter_dict["ped file"] % str(i))

    if not_found != []:
        sys.stdout.write("\nWARNING!\nThe following files were not found:\n")
        for msg in not_found:
            sys.stdout.write("\t%s\n" % msg)
        sys.stdout.write("\nExiting PONDEROSA...\n\n")
        sys.exit()

    #Check for errors in files
    errors = list()
    if not_found == []:
        if "_" in open(parameter_dict["match file"] % "1").readlines()[0].split()[1] and not parameter_dict["ilash"]:
            errors.append("iLASH file detected but --ilash flag not used.\n")
        for i in range(1,23):
            mapf = open(parameter_dict["map file"] % str(i)).readlines()
            if float(mapf[-1].split()[2]) == 0:
                errors.append(".map file must have distance in cM in 3rd column.\n")
                break
        if parameter_dict["age file"] != "":
            for ppl in open(parameter_dict["age file"]).readlines():
                age = ppl.split()[1]
                try:
                    age = int(age)
                except:
                    errors.append("Ages must be in second column and numerical.\n")
                    break
    if errors != []:
        sys.stdout.write("\nWARNING!\n")
        for msg in errors:
            sys.stdout.write(msg)
        sys.stdout.write("\nExiting PONDEROSA...\n\n")
        sys.exit()

    sys.stdout.write("\nChecking files...done\n")



