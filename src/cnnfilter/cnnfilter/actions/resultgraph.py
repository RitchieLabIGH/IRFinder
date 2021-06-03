import sys,os
import pandas as pd
import numpy as np
#EMT5p
if sys.argv[1] == "EMT5m":

    shortnames = [["T5moins_rep1", "T5moins_rep2", "T5moins_rep3"]]

if sys.argv[1] == "EMT5p":

    shortnames = [["T5plus_rep1", "T5plus_rep2", "T5plus_rep3"]]

if sys.argv[1] == "EMT1p":

    shortnames = [["T1plus_rep1", "T1plus_rep2", "T1plus_rep3"]]


folder = "/work/sylvain/IntronScanner/test/training/EMT_training/depth15_ir0.05_noir0.cov50long_.023allshort_rmALLnondirnoncongruant2021-05-07_16_53/"
folder = sys.argv[1]

pred = pd.read_csv("prediction_for_EMT_test.tsv",sep="\t|-|:",skiprows=1,header=None)
pred[2]=pred[2]+15
pred[3]=pred[3]-15
pred["id"]=pred[1].apply(str)+":"+pred[2].apply(str)+"-"+pred[3].apply(str)
pred["truelab"]="no"
for sit, sname in enumerate(shortnames[0]):
    data = pd.read_csv(sname+".tsv",delimiter="\t")
    data["id"]=data["Chr"].apply(str)+":"+data["Start"].apply(str)+"-"+data["End"].apply(str)
    pred.loc[pred[0]==sname].loc[pred['id'].isin(data["id"])]
    data.loc[data['id'].isin(pred.loc[pred[0] == sname]["id"])]
    pred.loc[pred[0] == sname]["truelab"]=data.loc[data['id'].isin(pred.loc[pred[0] == sname]["id"])]["Warnings"]
    (pred[1]==data["Chr"]) & (pred[2]==data["Start"]) & (pred[3]==data["End"])  & (pred[0]==sname)
    pred[pred[0]==sname]
# for sit, sname in enumerate(shortnames[0]):
#     txt = numpy.loadtxt(os.path.join(folder, sname, "IRFinder-IR-nondir.txt"), delimiter="\t)
#     data = np.genfromtxt(fname=os.path.join(folder, sname, "IRFinder-IR-nondir.txt"), delimiter="\t", skip_header=1)
#     sarray = open(os.path.join(folder, sname, "IRFinder-IR-nondir-AI.txt"), "rt")
#     print(os.path.join(folder, sname, "IRFinder-IR-nondir-AI.txt"))
#
#
#     line = _array.readline()
#
#     while True:
#         line = _array.readline()
#         if not line:
#             break
#
#             irratio = float(a[19])
#             classpred= a[20]
#             irratio=float(a[3])