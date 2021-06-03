import numpy as np
# import pandas as pd
import os
import json
import sys

# import glob

if sys.argv[1] == "EMT":
    folder = "/work/ritchieteam/EMT_irfinder/"
    shortnames = [["T5moins_rep1", "T5moins_rep2", "T5moins_rep3"], ["T1plus_rep1", "T1plus_rep2", "T1plus_rep3"],
                  ["T5plus_rep1", "T5plus_rep2", "T5plus_rep3"]]
    longnames = ["LONG/irfinderrest/T7_moins", "LONG/irfinderrest/T1_plus", "LONG/irfinderrest/T7_plus"]
if sys.argv[1] == "EMT5m":
    folder = "/work/ritchieteam/EMT_irfinder/"
    shortnames = [["T5moins_rep1", "T5moins_rep2", "T5moins_rep3"]]
    longnames = ["LONG/irfinderrest/T7_moins"]
if sys.argv[1] == "EMT5p":
    folder = "/work/ritchieteam/EMT_irfinder/"
    shortnames = [["T5plus_rep1", "T5plus_rep2", "T5plus_rep3"]]
    longnames = ["LONG/irfinderrest/T7_plus"]
if sys.argv[1] == "EMT1p":
    folder = "/work/ritchieteam/EMT_irfinder/"
    shortnames = [["T1plus_rep1", "T1plus_rep2", "T1plus_rep3"]]
    longnames = ["LONG/irfinderrest/T1_plus"]

else:
    if sys.argv[1] == "GM":
        folder = "/work/ritchieteam/GM12878/"
        shortnames = [["short/irfinderrest/SRR521453", "short/irfinderrest/SRR521447", "short/irfinderrest/SRR521448",
                       "short/irfinderrest/SRR521454", "short/irfinderrest/SRR521455"]]
        #shortnames = [["short/irfinderrest/SRR521447"]]
        longnames = ["long/irfinderrest/mergedDirectRNA.fastq"]

nondirmode=True


resume = open(os.path.join(folder, "resumeclasses.tsv"), "a")

# longloc = np.array([])
# # longarray=np.array([])
# longclass = np.array([])
tnoIRfalsepos = 0
tnoIRtruepos = 0
tIRfalsepos = 0
tIRtruepos = 0
for lit, lname in enumerate(longnames):
    longloc = np.array([])
    longclass = np.array([])

    # if os.path.isfile(os.path.join(folder,lname,"IRFinder-IR-dir-AI.txt")):
    #     _array=open(os.path.join(folder,lname,"IRFinder-IR-dir-AI.txt"), "rt")
    #     print(os.path.join(folder,lname,"IRFinder-IR-dir-AI.txt"))
    # else :
    #     _array=open(os.path.join(folder,lname,"IRFinder-IR-nondir-AI.txt"), "rt")
    #     print(os.path.join(folder, lname, "IRFinder-IR-nondir-AI.txt"))
    # while True:
    #     line=_array.readline()
    #     if not line:
    #         break
    #
    #     a=line.strip().split("\t")
    #     loc=a[0]
    #     # locaa=a[0].split(":")
    #     # locab=locaa[1].split("-")
    #     # loc="{}:{}-{}".format(locaa[0],int(locab[0])-15,int(locab[1])+15)
    #
    #     regiondata =np.array( json.loads(a[1]))
    #
    #     if len(regiondata) > 51& np.min(np.concatenate((regiondata[0:15 ,0],regiondata[-15: ,0])))>10:
    #
    #         if np.min(regiondata[25:-25, 0]) >= 1: #1
    #             longclass= np.append(longclass,"hIR")
    #             longloc= np.append(longloc,loc)
    #         else:
    #             if np.sum(regiondata[25:-25, 0] == 0)>= 1: #reduce windows
    #                 longclass = np.append(longclass, "noIR")
    #                 longloc = np.append(longloc, loc)
    #
    # _array.close()
    if os.path.isfile(os.path.join(folder, lname, "IRFinder-IR-dir.txt")):
        _array = open(os.path.join(folder, lname, "IRFinder-IR-dir.txt"), "rt")
        print(os.path.join(folder, lname, "IRFinder-IR-dir.txt"))
    else:
        _array = open(os.path.join(folder, lname, "IRFinder-IR-nondir.txt"), "rt")
        print(os.path.join(folder, lname, "IRFinder-IR-nondir.txt"))
    # read header
    line = _array.readline()
    clir = 0
    clnoir = 0
    ctclaudiofilter = 0
    while True:
        line = _array.readline()
        if not line:
            break

        a = line.strip().split("\t")
        loc = "{}:{}-{}".format(a[0], int(a[1]) - 15, int(a[2]) + 15)
        if loc=="9:124353694-124356778":
            titi=1
        if "9:" in loc :
            tot=1
        intronDepth = int(a[8])
        spliceExact = int(a[18])
        warning = a[20]
        irratio = float(a[19])
        if ((spliceExact + intronDepth) >= 25) and warning == "-":  # (int(a[2])-int(a[1]) > 51) and

            if (intronDepth >= 1) and (spliceExact >= 3):
                if irratio >= 0.1:# 0.05
                    ctclaudiofilter += 1
                    longclass = np.append(longclass, "hIR")
                    longloc = np.append(longloc, loc)
                    clir += 1
            else:
                if ((spliceExact + intronDepth) >= 50) and irratio == 0:
                    longclass = np.append(longclass, "noIR")
                    longloc = np.append(longloc, loc)
                    clnoir += 1

    _array.close()
    print("IR={} noIR={}".format(clir, clnoir))
    #coroborate nondir no ir
    if nondirmode==True:
        _arraynondir = open(os.path.join(folder, lname, "IRFinder-IR-nondir.txt"), "rt")
        print(os.path.join(folder, lname, "IRFinder-IR-nondir.txt"))
        line = _arraynondir.readline()
        loc="0"
        mask = np.ones(len(longloc), dtype=bool)
        noirrm=0
        hirrm=0
        while True:
            line = _arraynondir.readline()
            if not line:
                break

            a = line.strip().split("\t")
            loc = "{}:{}-{}".format(a[0], int(a[1]) - 15, int(a[2]) + 15)
            if loc == "9:124353694-124356778":
                titi = 1

            irratio = float(a[19])
            dirlocindex = np.where(longloc == loc)[0]
            if dirlocindex.size == 1:
                if irratio!=0:
                    if longclass[dirlocindex]=="noIR":
                        noirrm+=1
                        mask[dirlocindex]=False
                if irratio < 0.05 :
                    if  longclass[dirlocindex] == "hIR":
                        hirrm+=1
                        mask[dirlocindex] = False

        #ctnoIRrm=sum(mask==False)
        print("discordant noIR removed {}".format(noirrm))
        print("discordant IR removed {}".format(hirrm))
        clnoir=clnoir-noirrm
        clir=clir-hirrm
        longloc=longloc[mask]
        longclass=longclass[mask]

    print("IR={} noIR={}".format(clir, clnoir))
    print("claudioIR={}".format(ctclaudiofilter), flush=True)
    ## match with shortread
    # shortresfit = np.array([[]])
    # shortlocfit = np.array([])

    for sit, sname in enumerate(shortnames[lit]):
        shortresfit = np.array([[]])
        shortlocfit = np.array([])
        noIRfalsepos = 0
        noIRtruepos = 0
        IRfalsepos = 0
        IRtruepos = 0
        print(sname)
        resume.write("{}\n".format(sname))
        if False: #os.path.isfile(os.path.join(folder, sname, "IRFinder-IR-dir.txt")):
            sres = open(os.path.join(folder, sname, "IRFinder-IR-dir.txt"), "rt")
            print(os.path.join(folder, sname, "IRFinder-IR-dir.txt"))
        else:
            sres = open(os.path.join(folder, sname, "IRFinder-IR-nondir.txt"), "rt")
            print(os.path.join(folder, sname, "IRFinder-IR-nondir.txt"))
        # sres = open(os.path.join(folder, sname, "IRFinder-IR-dir.txt"), "rt")
        sres.readline()
        ct = 0
        while True:
            line = sres.readline()
            if not line:
                break

            a = line.strip().split("\t")

            locs = "{}:{}-{}".format(a[0], int(a[1]) - 15, int(a[2]) + 15)
            if loc == "9:124353694-124356778":
                1
            if a[20] == '-':
                if np.where(longloc == locs)[0].size == 1:  # (locs in longloc):
                    shortlocfit = np.append(shortlocfit, locs)
                    if shortresfit.size == 0:
                        shortresfit = np.array(a)
                    else:
                        shortresfit = np.vstack((shortresfit, a))
                    ct = ct + 1
                else:
                    tmpclass = "NA"

        sres.close()
        print("commons regions long approve vs short res ", ct, flush=True)
        resume.write("commons regions long approve vs short res {}\n".format(ct))
        # outresarray=np.array([shortlocfit.shape])
        stail = os.path.split(sname)[1]
        ltail = os.path.split(lname)[1]
        outres = open(os.path.join(folder, "classes", ltail + "__" + stail + ".tsv"), "w")
        # sarray = open(os.path.join(folder, sname, "IRFinder-IR-dir-AI.txt"), "rt")
        if False :#os.path.isfile(os.path.join(folder, sname, "IRFinder-IR-dir.txt")):
            sarray = open(os.path.join(folder, sname, "IRFinder-IR-dir-AI.txt"), "rt")
            print(os.path.join(folder, sname, "IRFinder-IR-dir-AI.txt"))
        else:
            sarray = open(os.path.join(folder, sname, "IRFinder-IR-nondir-AI.txt"), "rt")
            print(os.path.join(folder, sname, "IRFinder-IR-nondir-AI.txt"))
        # sarray.readline()

        # shortAIloc = np.array([])
        cthir = 0
        ctnoir = 0
        noirirratio = np.array([])
        while True:
            line = sarray.readline()
            a = line.strip().split("\t")[0].split(":")
            if len(a)<2:
                print(a)
            else:
                a = a[0] + ":" + a[1]
            ishort = np.where(shortlocfit == a)[0]
            if ishort.size == 1:

                tmpclass = longclass[longloc == a][0]
                tabshort = shortresfit[ishort, :][0]
                if tmpclass == "noIR":
                    noirirratio = np.append(noirirratio, float(tabshort[19]))
                    if (tabshort[20] == '-') & (float(tabshort[19]) >=0.01):#0.045 for 0.1long#& (float(tabshort[19]) >= 0.08) # 0.1 for 0.05 irratio long read(10 expre), 0.08 for 0.05 irratio long read(15 expre)
                        tmpclass = "noIR"  # confirm
                        noIRfalsepos += 1
                    else:
                        tmpclass = "NA"
                        noIRtruepos += 1
                # outresarray[ishort] = "{}\t0\t{}\n".format(a[0], tmpclass)

                if tmpclass == "hIR":
                    if (tabshort[20] == '-') & (float(tabshort[19]) >= 0.05):
                        IRtruepos += 1
                    else:
                        IRfalsepos += 1
                        tmpclass = "NA"
                    cthir = cthir + 1
                else:
                    if tmpclass == "noIR":
                        ctnoir = ctnoir + 1
                outres.write("{}\t0\t{}\n".format(a, tmpclass))

            else:
                if not line:
                    break
                # outresarray[ishort]="{}\t0\t{}\n".format(a[0], "NA")
                outres.write("{}\t0\t{}\n".format(a, "NA"))

            if not line:
                break
        print("hIR :", cthir)
        print("noIR :", ctnoir)

        # plt.figure()
        # plt.plot(range(len(noirirratio)),np.sort(noirirratio))
        #
        # plt.savefig( sname+"noIRirratio.png")
        # plt.close()
        tnoIRfalsepos += noIRfalsepos
        tnoIRtruepos += noIRtruepos
        tIRfalsepos += IRfalsepos
        tIRtruepos += IRtruepos

        outres.close()
        sarray.close()
        print("{}----\ttruepos\tfalsepos".format(lname))
        print("noIR\t{}\t{}".format(noIRtruepos, noIRfalsepos))
        print("IR\t{}\t{}".format(IRtruepos, IRfalsepos), flush=True)

        resume.write("{}\n-----\ttruepos\tfalsepos\n".format(lname))
        resume.write("noIR\t{}\t{}\n".format(noIRtruepos, noIRfalsepos))
        resume.write("IR\t{}\t{}\n".format(IRtruepos, IRfalsepos))

resume.write("===================\nTOTAL\n-----\ttruepos\tfalsepos\n".format(lname))
resume.write("noIR\t{}\t{}\n".format(tnoIRtruepos, tnoIRfalsepos))
resume.write("IR  \t{}\t{}\n".format(tIRtruepos, tIRfalsepos))
resume.close()
print("TOTAL\n-----\ttruepos\tfalsepos\n".format(lname))
print("noIR\t{}\t{}\n".format(tnoIRtruepos, tnoIRfalsepos))
print("IR  \t{}\t{}\n".format(tIRtruepos, tIRfalsepos))
# open(fname, "rt")
