###This is the script using bedtools to filter the ChromoPainter BED archaic blocks to a high confidence set that:
#1. Overlap HMM archaic blocks (> 0.1% overlap)
#2. Overlap S* archaic windows (> 0.1% overlap)
#3. Are not covered by ChromoPainter archaic blocks inferred to be the other archaic (<99.9% overlap)

#It requires four BED files per chromosome copy of per individual - CP archaic target, HMM archaic target, S* and CP archaic outgroup
#Each BED file indicates the blocks inferred to be archaic introgression by the method in question, e.g.

"""
chr1	2903159	2903783
chr1	2907644	2919963
chr1	2920431	2927624
chr1	2931250	2932735
...
"""


###GSJ 06/01/2019

import subprocess
import tempfile
import os

def filter_highconf_binary_hmm_sstar_cparchaicalt999(ind_names, cp_deni_bedfiles, cp_nean_bedfiles, hmm_deni_bedfiles, hmm_nean_bedfiles, sstar_bedfiles, outfile_base = os.getcwd() + '/hconf_%s_%s.bed', tmp_folder = os.getcwd(), bedtools_dir = './bedtools2/bin/'):
    """
    ind_names is the list of individuals to filter
    cp_deni_bedfiles is the list of ChromoPainter archaic Denisovan introgression BED files, same order as ind_names
    cp_nean_bedfiles is the list of ChromoPainter archaic Neanderthal introgression BED files, same order as ind_names
    hmm_deni_bedfiles is the list of HMM archaic Denisovan introgression BED files, same order as ind_names
    hmm_nean_bedfiles is the list of HMM archaic Neanderthal introgression BED files, same order as ind_names
    sstar_bedfiles is the list of SStar archaic introgression BED files, same order as ind_names

    bed files must be sorted

    bedtools_dir is the path to the /bedtools2/bin/ directory
    tmp_folder is a path to a folder where some temporary files can be created (then deleted) during the filtering
    """
    assert len(ind_names) == len(cp_deni_bedfiles)
    assert len(ind_names) == len(cp_nean_bedfiles)
    assert len(ind_names) == len(hmm_deni_bedfiles)
    assert len(ind_names) == len(hmm_nean_bedfiles)
    assert len(ind_names) == len(sstar_bedfiles)
    TMP_FOLDER = tmp_folder
    BEDTOOLS_DIR = bedtools_dir
    assess_files = [cp_deni_bedfiles, cp_nean_bedfiles, hmm_deni_bedfiles, hmm_nean_bedfiles, sstar_bedfiles]
    #Deni first
    for ind in range(len(ind_names)):
        print "deni trimming", ind_names[ind]
        #These work by removing all regions that are > 0.001 overlap.
        #And then subtracting the remaining from the original file.
        #HMM first
        command_intersect_hmm = BEDTOOLS_DIR + 'bedtools subtract -a %s -b %s -N -f %f -sorted -g %s' %(assess_files[ind][0], assess_files[ind][2], 0.001, GENOME_FILE) + PIPE + BEDTOOLS_DIR + 'bedtools subtract -a %s -b stdin -sorted -g %s' %(assess_files[ind][0], GENOME_FILE)
        f_intersect_hmm = tempfile.NamedTemporaryFile(mode = 'w+b', dir = TMP_FOLDER)
        bedtools_process_intersect_hmm = subprocess.Popen(args = command_intersect_hmm, stdout = f_intersect_hmm, stderr = f_intersect_hmm, shell = True)
        bedtools_process_intersect_hmm.communicate(input=None)
        f_intersect_hmm.flush()
        f_intersect_hmm.seek(0)
            
        command_intersect_sstar = BEDTOOLS_DIR + 'bedtools subtract -a %s -b %s -N -f %f -sorted -g %s' %(f_intersect_hmm.name, assess_files[ind][4], 0.001, GENOME_FILE) + PIPE + BEDTOOLS_DIR + 'bedtools subtract -a %s -b stdin -sorted -g %s' %(f_intersect_hmm.name, GENOME_FILE)
        f_intersect_sstar = tempfile.NamedTemporaryFile(mode = 'w+b', dir = TMP_FOLDER)
        bedtools_process_intersect_sstar = subprocess.Popen(args = command_intersect_sstar, stdout = f_intersect_sstar, stderr = f_intersect_sstar, shell = True)
        bedtools_process_intersect_sstar.communicate(input=None)
        f_intersect_sstar.flush()
        f_intersect_sstar.seek(0)

        command_subtract_altarchaic = BEDTOOLS_DIR + 'bedtools subtract -a %s -b %s -N -f %f -sorted -g %s' %(f_intersect_sstar.name, assess_files[ind][1], 0.999, GENOME_FILE)
        f_subtract_altarchaic = tempfile.NamedTemporaryFile(mode = 'w+b', dir = TMP_FOLDER)
        bedtools_process_subtract_altarchaic = subprocess.Popen(args = command_subtract_altarchaic, stdout = f_subtract_altarchaic, stderr = f_subtract_altarchaic, shell = True)
        bedtools_process_subtract_altarchaic.communicate(input=None)
        f_subtract_altarchaic.flush()
        f_subtract_altarchaic.seek(0)
        
        #And write out
        
        with open(outfile_base %('deni', ind_name), 'wb') as f_out_deni:
            for line in f_subtract_altarchaic:
                f_out_deni.write(line)
        f_intersect_hmm.close()
        f_intersect_sstar.close()
        f_subtract_altarchaic.close()
    
    #And Nean
    for ind in range(len(ind_names)):
        print "nean trimming", ind_names[ind]
        #HMM first
        command_intersect_hmm = BEDTOOLS_DIR + 'bedtools subtract -a %s -b %s -N -f %f -sorted -g %s' %(assess_files[ind][1], assess_files[ind][3], 0.001, GENOME_FILE) + PIPE + BEDTOOLS_DIR + 'bedtools subtract -a %s -b stdin -sorted -g %s' %(assess_files[ind][1], GENOME_FILE)
        f_intersect_hmm = tempfile.NamedTemporaryFile(mode = 'w+b', dir = TMP_FOLDER)
        bedtools_process_intersect_hmm = subprocess.Popen(args = command_intersect_hmm, stdout = f_intersect_hmm, stderr = f_intersect_hmm, shell = True)
        bedtools_process_intersect_hmm.communicate(input=None)
        f_intersect_hmm.flush()
        f_intersect_hmm.seek(0)
            
        command_intersect_sstar = BEDTOOLS_DIR + 'bedtools subtract -a %s -b %s -N -f %f -sorted -g %s' %(f_intersect_hmm.name, assess_files[ind][4], 0.001, GENOME_FILE) + PIPE + BEDTOOLS_DIR + 'bedtools subtract -a %s -b stdin -sorted -g %s' %(f_intersect_hmm.name, GENOME_FILE)
        f_intersect_sstar = tempfile.NamedTemporaryFile(mode = 'w+b', dir = TMP_FOLDER)
        bedtools_process_intersect_sstar = subprocess.Popen(args = command_intersect_sstar, stdout = f_intersect_sstar, stderr = f_intersect_sstar, shell = True)
        bedtools_process_intersect_sstar.communicate(input=None)
        f_intersect_sstar.flush()
        f_intersect_sstar.seek(0)

        command_subtract_altarchaic = BEDTOOLS_DIR + 'bedtools subtract -a %s -b %s -N -f %f -sorted -g %s' %(f_intersect_sstar.name, assess_files[ind][0], 0.999, GENOME_FILE)
        f_subtract_altarchaic = tempfile.NamedTemporaryFile(mode = 'w+b', dir = TMP_FOLDER)
        bedtools_process_subtract_altarchaic = subprocess.Popen(args = command_subtract_altarchaic, stdout = f_subtract_altarchaic, stderr = f_subtract_altarchaic, shell = True)
        bedtools_process_subtract_altarchaic.communicate(input=None)
        f_subtract_altarchaic.flush()
        f_subtract_altarchaic.seek(0)

        #And write out
        with open(outfile_base %('nean', ind_name), 'wb') as f_out_nean:
            for line in f_subtract_altarchaic:
                f_out_nean.write(line)
        f_intersect_hmm.close()
        f_intersect_sstar.close()
        f_subtract_altarchaic.close()
    return None
