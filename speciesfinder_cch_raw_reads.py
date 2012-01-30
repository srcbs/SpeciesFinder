#this script will parse the chh directory and will look for input files
''' for each project id the complete genome will be tested and the information about the taxonomy will stored in the directory with the output files '''
import os
import time
import sys
import subprocess

sys.path.append("/home/panfs/cbs/projects/cge/people/salvocos/scripts/lib/python/xmsub_manager/")
import xmsub_manager

####### GLOBAL VARIABLES ###########

chhRoot = "/panfs1/cge/data/collab/chh/data_pw/" #root to the data from chh
outDir = "/panvol1/simon/projects/cge/test/speciesfinder/"
chhDirNames = ["Bacillus_subtilis", "pediococcus", "propionibacterium", "streptococcu_thermophilus", "lactococcus_lactis", "oenococcus_oenii"] 
chhMultiLevDirNames = ["Lactobacillus", "bifidobacterium"] #these subdirectories need to be parsed in a different way

#scriptPath = "/panfs1/cge/people/salvocos/CGE/16S_prediction/scripts/ssu/speciesfinder.py"
scriptPath = "/panfs1/cge/people/salvocos/CGE/16S_prediction/scripts/ssu/speciesfinder4.py"
#idsListFile = "/panfs1/cge/people/salvocos/CGE/16S_prediction/tables/annotated_genomes_test_list.txt" # _tmp.txt
pythonPath = "/tools/bin/python"
forbiddenExt = ["gz", "ini", "tar"]
validExt = ["txt", "fq"]
fileExt = ".fsa"
Mbases = 50
#refTableFile = "/panfs1/cge/people/salvocos/CGE/16S_prediction/tables/final_main_table.txt"

####### FUNCTIONS #########

#this function will get the parameters from the command line
def getParams():
    global Mbases
    if len(sys.argv) < 2:
        print "ERROR: too few arguments"
        sys.exit()
    elif len(sys.argv) > 2:
        print "ERROR: too many arguments"
        sys.exit()
    for index, arg in enumerate(sys.argv):
        print "ARG_" + str(index) + " ::\t" + arg
    Mbases = sys.argv[1]
    try:
        Mbases = float(Mbases)
    except e, ValueError:
        print e
    if Mbases < 1:
        print "ERROR: the reads number must be not lower than 1Mb"
        sys.exit(1)



###### MAIN #####

#XMSUB example
#xmsub -N ssu_test_chh_genomes -m ae -M salvocos@cbs.dtu.dk  -q cge --host=cge-s2 -l walltime=48:00:00,mem=2000m,nodes=1:ppn=1 -d /panfs1/cge/people/salvocos/CGE/16S_prediction/logs/xmsub_logs/ -ro /panfs1/cge/people/salvocos/CGE/16S_prediction/logs/std_out/ssu_test_chh_genomes.out -re /panfs1/cge/people/salvocos/CGE/16S_prediction/logs/std_out/ssu_test_chh_genomess.err -de /tools/bin/python  /panfs1/cge/people/salvocos/CGE/16S_prediction/scripts/test_chh_raw_reads.py


getParams()

#check if the files and paths exist
if not os.path.isdir(chhRoot):
    print "ERROR: " + chhRoot + " is not a valid directory"
    sys.exit(1)
if not os.path.isdir(outDir):
    print "ERROR: " + outDir + " is not a valid directory"
    sys.exit(1)
for el in chhDirNames:
    subDir = chhRoot + el + "/"
    if not os.path.isdir(subDir):
        print "ERROR: " + subDir + " is not a valid directory"
        sys.exit(1)
for el in chhMultiLevDirNames:
    subDir = chhRoot + el + "/"
    if not os.path.isdir(subDir):
        print "ERROR: " + subDir + " is not a valid directory"
        sys.exit(1)


#command example
''' python speciesfinder.py <raw_read.fq> --Mbases 50 --sample <out_dir> --tax '''

#traverse the directories and look for possible input files
for el in chhDirNames:
    subDir = chhRoot + el + "/"
    dirFiles = os.listdir(subDir)
    #fPath = inputFastaDir + el.strip() + fileExt
    #outPath = outDir + el.strip()  + "/"

    print el
    for f in dirFiles:
        fPath = subDir + f
        if os.path.isfile(fPath):
            #print "\t" + fPath
            splitted = fPath.split(".")
            ext = splitted[-1]
            if ext in validExt:
                print "INPUT FILE: " + fPath
                #now I will create the output directories
                outDirRoot = outDir + str(Mbases) + "/"
                print "OUTPUT ROOT :: " + outDirRoot
                outName = outDirRoot  + el + "_" + f
                if not os.path.isdir(outDirRoot):
                    os.makedirs(outDirRoot)
                print outName
                cmd = pythonPath + " " + scriptPath + " " + fPath + " --Mbases " + str(Mbases)  + " --sample " + outName + " --tax"
                subprocess.Popen(cmd, shell=True).wait() # execute the system call


#now I will predict the files in the subdirectories

for el in chhMultiLevDirNames:
    subDir = chhRoot + el + "/"
    dirFiles = os.listdir(subDir)
    for dir in dirFiles:
        speciesDir = subDir + dir + "/"
        files = os.listdir(speciesDir)
        for f in files:
            splitted = f.split(".")
            ext = splitted[-1]
            fPath = speciesDir + f
            if ext in validExt:
                #print ext
                print "INPUT FILE: " + fPath
                            #now I will create the output directories
                outDirRoot = outDir + str(Mbases) + "/"
                print "OUTPUT ROOT :: " + outDirRoot
                outName = outDirRoot  + el + "_" + dir + "_" + f
                if not os.path.isdir(outDirRoot):
                    os.makedirs(outDirRoot)
                print outName
                cmd = pythonPath + " " + scriptPath + " " + fPath + " --Mbases " + str(Mbases)  + " --sample " + outName + " --tax"
                subprocess.Popen(cmd, shell=True).wait() # execute the system call
