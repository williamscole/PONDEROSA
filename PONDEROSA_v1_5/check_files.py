import os.path
import sys

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
