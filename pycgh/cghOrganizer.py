import os
from shutil import copyfile

import pycgh.readers

MAX_FILES = 200

def organize(CGH_DIR, OUTPUT_DIR):
    '''
    This function copies all cgh (agilent) files found in the input folder CGH_DIR
    in different subfolders of a given folder OUTPUT_DIR, based on their resolution (the number
    of probes)
    '''
    
    ### create output directory if it does not exist
    if not os.path.exists(OUTPUT_DIR):
        os.mkdir(OUTPUT_DIR)
            
    i = 0
    
    file_list = sorted(os.listdir(CGH_DIR))
    
    for f in file_list:
        
        i += 1
        
        if not f.endswith('txt'):
            print "Skipping file %s" % f
            continue
        
        input_file = os.path.join(CGH_DIR, f)
        
        ### load the cgh file
        acgh = pycgh.readers.agilent(input_file)
        
        N_probes = len(acgh.M['test_signal'])
        dest_dir = os.path.join(OUTPUT_DIR, str(N_probes)) ## retrieve destination folder
        
        ### if there is no such destination folder, create it
        if not os.path.exists(dest_dir):
            os.mkdir(dest_dir)
            
        
        ### copy the original file to the new folder
        dest_file = os.path.join(dest_dir, f)
        print "copying file %s to %s [%d of %d]" % (input_file, dest_file, i, len(file_list))
        copyfile(input_file, dest_file)
        
        ### dereference the variable acgh to free up memory
        del acgh
        
        if i >= MAX_FILES:
            break

if __name__ == '__main__':
    CGH_DIR = '/home/matteo/projects/poland/aCGH_Poland_probandi_anonimi'
    OUTPUT_DIR = '/home/matteo/projects/poland/aCGH_Poland_probandi_anonimi_split'
    organize(CGH_DIR, OUTPUT_DIR)
