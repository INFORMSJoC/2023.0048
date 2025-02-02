import math
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import matplotlib.ticker

def Generate_a_p_t_function(T):
    Trad=T*math.pi/180
    omega=random.uniform(2*math.pi/Trad, 4*math.pi/Trad)
    phi = 0
    a = random.uniform(0, 1)
    b = random.uniform(a, 1)
    d = random.uniform(0, 0.1)
    c = random.uniform(0, (1 - d) / Trad)
    normalizefactor = random.uniform(1, 100)
    A ,Au, Al = [], [], []
    array=np.linspace(0,Trad,T)
    for i in array:
        x = a * math.cos(omega * i + phi) + b
        A.append(x)
        Au.append(x * (1 + c * i + d))
        Al.append(x * (1 - c * i - d))
    A = np.array(A) / (normalizefactor * max(Au))
    Al = np.array(Al) / (normalizefactor * max(Au))
    Au = np.array(Au) / (normalizefactor * max(Au))
    return Al,A,Au




# how d1 and d2 changes with k in app1 part 2

matplotlib.rc('text', usetex=True)
matplotlib.rcParams.update({'font.size': 42})
x = np.arange(0,19,1)
MINMAX_D=[0.02291, 0.02122, 0.01309, 0.01291, 0.0104, 0.00923, 0.0092, 0.00886, 0.00797, 0.00791, 0.0073, 0.0073, 0.0073, 0.0073, 0.0073]
TOTALAVG_D=[0.00597, 0.00406, 0.0025, 0.00179, 0.00167, 0.00165, 0.00134, 0.00144, 0.00112, 0.00098, 0.00097, 0.00094, 0.00087, 0.00085, 0.00089]
# percentimprovementMinmax_D=[100*(MINMAX_D[0]-MINMAX_D[m])/MINMAX_D[0] for m in range(len(MINMAX_D))]
percentimprovementTotalavg_D=[100*(TOTALAVG_D[0]-TOTALAVG_D[m])/TOTALAVG_D[0] for m in range(len(TOTALAVG_D))]

MinMaxDistance_from16=[0.02565, 0.03752, 0.03752, 0.04007, 0.03957, 0.03957, 0.04007, 0.03957, 0.04007, 0.04147, 0.04147, 0.04147, 0.04147, 0.04147, 0.04147]
Total_Average_distance_from16= [0.01594, 0.01785, 0.01941, 0.02012, 0.02024, 0.02026, 0.02057, 0.02047, 0.02079, 0.02093, 0.02094, 0.02097, 0.02104, 0.02106, 0.02102]
percentimprovementMinmax_D_from16=[-100*(MinMaxDistance_from16[0]-MinMaxDistance_from16[m])/MinMaxDistance_from16[0] for m in range(len(MinMaxDistance_from16))]
percentimprovementTotalavg_D_from16=[-100*(Total_Average_distance_from16[0]-Total_Average_distance_from16[m])/Total_Average_distance_from16[0] for m in range(len(Total_Average_distance_from16))]
# costtomultiplyapp2=16000000000*(13.6/53)/(26000000)
costtomultiplyapp1=(13.6/20) #app1
# plt.plot(x,percentimprovementMinmax_D,label='Maximum regret',color='k',ls='-')
# plt.plot(x,percentimprovementTotalavg_D,label='Average deviation',color='k',ls='-.')
# plt.plot(x,percentimprovementMinmax_D_from16,label='Maximum regret',color='k',ls='-')
# plt.plot(x,percentimprovementTotalavg_D_from16,label='Average deviation',color='k',ls='-.')

import pickle
# below is app1
with open('totalcostapp1forplotting.data', 'rb') as f:
    dataread=pickle.load(f)
f.close()
# # below is for robust data app 1
with open('totalcostapp1forplottingRobust.data', 'rb') as f:
    datareadrobust=pickle.load(f)
f.close()

# # # below is app2
# with open('totalcostapp2forplottingpart1.data', 'rb') as f:
#     datareadpart1=pickle.load(f)
# with open('totalcostapp2forplotting.data', 'rb') as f:
#     datareadpart2=pickle.load(f)
# f.close()
# dataread=datareadpart1+datareadpart2
#
# # below is for robust data app 2
# import pickle
# with open('totalcostapp2forplottingRobust.data', 'rb') as f:
#     datareadrobust=pickle.load(f)
# f.close()
plt.figure(figsize=(12,8))
plt.plot(x[0:7],np.array(dataread)[0:7]*costtomultiplyapp1,label='Regret',color='k',ls='-.')
plt.scatter(x[0:7],np.array(dataread)[0:7]*costtomultiplyapp1, color='r', s=60)
plt.plot(x[0:7],np.array(datareadrobust)[0:7]*costtomultiplyapp1,label='Minmax',color='b',ls='-.')
plt.scatter(x[0:7],np.array(datareadrobust)[0:7]*costtomultiplyapp1, color='r', s=60)
# plt.plot(x[0:13],np.array(dataread)[0:13]*costtomultiplyapp2,label='Regret',color='k',ls='-.')
# plt.scatter(x[0:13],np.array(dataread)[0:13]*costtomultiplyapp2, color='r', s=60)
# plt.plot(x[0:13],np.array(datareadrobust)[0:13]*costtomultiplyapp2,label='Minmax',color='b',ls='-.')
# plt.scatter(x[0:13],np.array(datareadrobust)[0:13]*costtomultiplyapp2, color='r', s=60)
plt.xlabel(r'$\rm{Allowable\; number\; of\; changes}\; (K)$', usetex=True)
plt.ylabel(r'$\rm{Cost}$ ' + r'$(\rm{in\;} \$1\rm{M})$', usetex=True)
plt.legend()








# # for the infection change of part 2 app 2 in the case study for k=3 or K=2 (same thing)
# k=3
# path_matrix_of_ks=np.array([ 0, 23, 34, 52])
# matplotlib.rc('text', usetex=True)
# matplotlib.rcParams.update({'font.size': 29})
# x = np.arange(0,53,1)
# Y=[[0]*5]*53
# budget_path_3d_matrix_of_ks=[[10.2,  9.5, 21.1,  4.2,  0. ],
#        [ 9.8,  8.9, 20.2,  5.2,  1. ],
#        [10.4,  9.9, 21.5,  2.7,  0.5]]
# m=0
# for i in range(len(Y)):
#     if i <= path_matrix_of_ks[m+1]:
#         Y[i]=budget_path_3d_matrix_of_ks[m]
#     else:
#         m=m+1
#         Y[i] = budget_path_3d_matrix_of_ks[m]
# Y=np.array(Y)
#
# plt.plot(x,Y[:,0],label='$x_{1}^{*}$',color='k')
# plt.plot(x,Y[:,1],label='$x_{2}^{*}$',ls='-.',color='k')
# plt.plot(x,Y[:,2],label='$x_{3}^{*}$',ls='--',color='k')
# plt.plot(x,Y[:,3],label='$x_{4}^{*}$',ls=':',color='k')
# plt.plot(x,Y[:,4],label='$x_{5}^{*}$',marker='*',color='k')
# plt.legend(loc='upper right', fontsize=24, prop={"size":20},ncol=3,bbox_to_anchor=(1, 0.92))
# plt.text(6,11,'HIV',)
# plt.text(6,8,'HBV')
# plt.text(6,19,'HCV')
# plt.text(6,5,'Babesiosis')
# plt.text(6,0.8,'WNV')
# plt.xticks(np.arange(53), (' ', 'Jan', ' ', ' ', 'Feb', ' ', ' ', ' ', 'Mar', ' ', ' ', ' ', 'Apr', ' ', ' ', ' ', ' ','May', ' ', ' ', ' ',' Jun', ' ', ' ', ' ', 'Jul', ' ', ' ', ' ',' ' ,'Aug', ' ', ' ', ' ', 'Sep', ' ', ' ', ' ', 'Oct', ' ', ' ', ' ', ' ' ,'Nov', ' ', ' ', ' ', 'Dec', ' ', ' ', ' ', ' ', ' '))
# plt.grid(False)
# plt.show()







# for the bargraph of part 2 app 1 WNV of the case study

# matplotlib.rc('text', usetex=True)
# matplotlib.rcParams.update({'font.size': 30})
#
# plt.subplot(1,2,1)
# x = np.arange(0,53,1)
# y = 0.001*np.array([6, 50, 2, 1, 3, 3, 0, 0, 7, 2, 11, 4, 4, 6, 8, 12, 10, 18, 29, 51, 47, 107, 134, 177, 294, 429, 841, 1275, 1925, 2734, 4417, 4821, 5769, 5737, 5479, 4672, 3886, 2596, 1831, 1198, 702, 475, 329, 217, 164, 118, 64, 59, 59, 36, 18, 35, 1])
#
# plt.bar(x,y, facecolor='k')
#
# plt.xticks(np.arange(53), ('Jan ', ' ', ' ', ' ', 'Feb', ' ', ' ', ' ', 'Mar', ' ', ' ', ' ', 'Apr', ' ', ' ', ' ', ' ','May', ' ', ' ', ' ',' Jun', ' ', ' ', ' ', 'Jul', ' ', ' ', ' ',' ' ,'Aug', ' ', ' ', ' ', 'Sep', ' ', ' ', ' ', 'Oct', ' ', ' ', ' ', ' ' ,'Nov', ' ', ' ', ' ', 'Dec', ' ', ' ', ' ', ' ', ' '))
# plt.xlabel('Time (weeks)')
# plt.ylabel(r'Reported WNV cases $\times 10^{3}$',labelpad=-60)
# plt.grid(False)
# # plt.xticks(rotation = 45)
# # ax[0].show()
#
#
# plt.subplot(1,2,2)
# # for the graph of part 2 app 1 of the case study
# x = np.arange(0, 20, 1)
# y=np.arange(26,46,1)
# WNVPrevalence_Data_List = 100*np.array([0.000023472,0.000039521,0.0000644458,0.00010147,0.000153794,0.00023497,0.000315464,0.000402724,0.000470892,0.000495,0.000491618,0.000448873,0.000377572,0.000294762,0.000215901,0.000143624,0.0000954769,0.0000613813,0.0000391124,0.0000268997])
# WNVPrevalence_Data_List_Upper=100*np.array([3.74784e-05, 0.00010418250000000001, 0.0001753244, 0.000252516, 0.00033735800000000006, 0.0004758535, 0.0006140085, 0.0007555465, 0.000887538, 0.0009975000000000001, 0.0009586505000000001, 0.00090012, 0.0008273110000000001, 0.0007487475000000001, 0.0005999664999999999, 0.000454477, 0.00032105245, 0.00019465364999999999, 0.0001547977, 0.00011996985000000001])
# WNVPrevalence_Data_List_Lower=100*np.array([1.249111e-05, 2.2236880000000002e-05, 3.642055e-05, 5.66539e-05, 8.453720000000001e-05, 0.00012799715, 0.0001711161, 0.00021761804999999998, 0.00025457405, 0.0002695, 0.000266719, 0.00024425655, 0.00020751605, 0.0001650211, 0.00012238295, 8.303685e-05, 5.575565e-05, 3.5500235e-05, 2.3523285e-05, 1.6574435e-05])
# #
#
#
# plt.plot(y,WNVPrevalence_Data_List_Lower,ls=' ')
# # plt.plot(y,WNVPrevalence_Data_List,color='k',label='WNV',linewidth=0.8)
# plt.plot(y,WNVPrevalence_Data_List_Upper,ls=' ',color='k')
#
# plt.scatter(y,WNVPrevalence_Data_List_Lower,color='k',marker='_',s=500)
# plt.scatter(y,WNVPrevalence_Data_List,color='k')
# plt.scatter(y,WNVPrevalence_Data_List_Upper,color='k',marker='_',s=500)
#
# dummy=[]
# for i in x:
#     dummy.append((WNVPrevalence_Data_List_Lower[i],WNVPrevalence_Data_List_Upper[i]))
#
# plt.plot((y,y),([i for (i,j) in dummy], [j for (i,j) in dummy]),c='black',ls=':')
# plt.xticks(y, (' Jun', ' ', ' ', ' ', 'Jul', ' ', ' ', ' ','Aug' ,' ', ' ', ' ', 'Sep', ' ', ' ', ' ', 'Oct', ' ', ' ', ' '))
# # plt.xticks(rotation=45)
# plt.fill_between(y,WNVPrevalence_Data_List_Lower,WNVPrevalence_Data_List_Upper,alpha=0.1)
# # plt.legend(loc='best',fontsize=40)
# plt.xlabel('Time (weeks)')
# plt.ylabel('Prevalence rate $(\%)$',labelpad=-105)
# # ax[1].show()












# for the graph of part 1 of the case study INFORMS PRESENTATION
x = np.arange(0, 100, 1)
y = [0.00016,0.000479798387335,0.000782427175445,0.00106823676432,0.00133757755393,0.00159079994427,0.00182825433533,0.00205029112708,0.00225726071951,0.00244951351261,0.00262739990635,0.00279127030073,0.00294147509572,0.00307836469132,0.0032022894875,0.00331359988425,0.00341264628155,0.00349977907939,0.00357534867775,0.00363970547661,0.00369319987597,0.0037361822758,0.00376900307608,0.00379201267681,0.00380556147796,0.00380999987953,0.00380567828148,0.00379294708382,0.00377215668651,0.00374365748956,0.00370779989293,0.00366493429662,0.0036154111006,0.00355958070487,0.0034977935094,0.00343039991419,0.00335775031921,0.00328019512444,0.00319808472989,0.00311176953551,0.00302159994131,0.00292792634727,0.00283109915336,0.00273146875957,0.0026293855659,0.00252519997231,0.0024192623788,0.00231192318535,0.00220353279194,0.00209444159857,0.0019850000052,0.00187555841183,0.00176646721844,0.00165807682501,0.00155073763153,0.00144480003798,0.00134061444435,0.00123853125062,0.00113890085677,0.00104207366279,0.00094840006867,0.000858230474381,0.000771915279913,0.00068980488525,0.000612249690376,0.000539600095274,0.00047220649993,0.000410419304327,0.000354588908449,0.00030506571228,0.000262200115805,0.000226342519007,0.000197843321871,0.00017705292438,0.000164321726519,0.000160000128272,0.000164438529623,0.000177987330555,0.000200996931054,0.000233817731102,0.000276800130685,0.000330294529786,0.00039465132839,0.00047022092648,0.00055735372404,0.000656400121055,0.000767710517508,0.000891635313384,0.00102852490867,0.00117872970334,0.00134260009739,0.0015204864908,0.00171273928355,0.00191970887563,0.00214174566702,0.0023792000577,0.00263242244767,0.00290176323689,0.00318757282537,0.00349020161308]
x=x[::2]
y=y[::2]

y_err1 = [(1+0.2)*j for j in y]
y_err2 = [(1-0.2)*j for j in y]

matplotlib.rc('text', usetex=True)
matplotlib.rcParams.update({'font.size': 38})
plt.scatter(x,y,color='k',label='prevalence rate $(p_{t})$', marker='x', s=100)
plt.scatter(x,y_err1,color='k',marker='_',s=300)
plt.scatter(x,y_err2,color='k',marker='_',s=300)
plt.axvline(x=0, ls='--', c='r')
plt.text(30,-0.0005,'$t_{1}$',c='r')
plt.axvline(x=30, ls='--', c='r')
plt.text(65,-0.0005,'$t_{2}$',c='r')
plt.axvline(x=65, ls='--', c='r')
plt.text(90,-0.0005,'$t_{3}$',c='r')
plt.axvline(x=90, ls='--', c='r')
plt.axvline(x=100, ls='--', c='r')
plt.fill_between(x,y_err1,y_err2,alpha=0.2)
plt.legend(loc='upper right')
plt.xlabel('Time')
# plt.ylabel('$p(t)$',fontsize=25)
plt.arrow(0, 0.0006, 30, 0,width=0.00002,length_includes_head=True,head_width=0.0001,head_length=1.2,color='b')
plt.arrow(30, 0.0006, -30, 0,width=0.00002,length_includes_head=True,head_width=0.0001,head_length=1.2,color='b')
plt.text(11,0.0008,'$x_{0,t_{1}}$')
plt.arrow(30, 0.0006, 35, 0,width=0.00002,length_includes_head=True,head_width=0.0001,head_length=1.2,color='b')
plt.arrow(65, 0.0006, -35, 0,width=0.00002,length_includes_head=True,head_width=0.0001,head_length=1.2,color='b')
plt.text(44,0.0008,'$x_{t_{1},t_{2}}$')
plt.arrow(65, 0.0006, 25, 0,width=0.00002,length_includes_head=True,head_width=0.0001,head_length=1.2,color='b')
plt.arrow(90, 0.0006, -25, 0,width=0.00002,length_includes_head=True,head_width=0.0001,head_length=1.2,color='b')
plt.text(73,0.0008,'$x_{t_{2},t_{3}}$')
plt.arrow(90, 0.0006, 10, 0,width=0.00002,length_includes_head=True,head_width=0.0001,head_length=1.2,color='b')
plt.arrow(100, 0.0006, -10, 0,width=0.00002,length_includes_head=True,head_width=0.0001,head_length=1.2,color='b')
plt.text(91,0.0008,'$x_{t_{3},T-1}$')
plt.text(48,0.004,'$K=3$')
dummy=[]
for i,K in enumerate(x):
    dummy.append((y_err1[i],y_err2[i]))
plt.plot((x, x), ([i for (i, j) in dummy], [j for (i, j) in dummy]), c='black', ls=':')
plt.show()









# for the graph of part 2 app 2 of the case study
# x = np.arange(0, 53, 1)
# HIVL = 100*np.array([0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005])
# HIV = 100*np.array([0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007])
# HIVU = 100*np.array([0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01])
# HBVL = 100*np.array([0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025])
# HBV = 100*np.array([0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345])
# HBVU = 100*np.array([0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044])
# HCVL = 100*np.array([0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013])
# HCV = 100*np.array([0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016])
# HCVU = 100*np.array([0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019])
# BabesiosisL = 100*np.array([0.000433926,0.000390934,0.000347943,0.000304951,0.00033009,0.000355229,0.000380369,0.000405508,0.000423816,0.000442124,0.000460432,0.00047874,0.000518963,0.000559186,0.000599408,0.000639631,0.000679854,0.000833696,0.000987538,0.001141379,0.001295221,0.001608916,0.00192261,0.002236305,0.00255,0.002445289,0.002340579,0.002235868,0.002131157,0.002026447,0.001784071,0.001541695,0.00129932,0.001056944,0.00100202,0.000947096,0.000892172,0.000837248,0.000808174,0.0007791,0.000750026,0.000720952,0.000691877,0.000644605,0.000597332,0.000550059,0.000502786,0.000455513,0.00040824,0.000360968,0.000313695,0.000266422,0.000219149])
# Babesiosis = 100*np.array([0.000655144,0.000590234,0.000525325,0.000460416,0.000498371,0.000536327,0.000574282,0.000612237,0.000639879,0.00066752,0.000695162,0.000722803,0.000783532,0.000844261,0.000904989,0.000965718,0.001026447,0.001258717,0.001490988,0.001723259,0.001955529,0.002429147,0.002902765,0.003376382,0.00385,0.003691907,0.003533815,0.003375722,0.00321763,0.003059537,0.002693597,0.002327658,0.001961718,0.001595778,0.001512854,0.001429929,0.001347005,0.001264081,0.001220184,0.001176288,0.001132392,0.001088495,0.001044599,0.000973227,0.000901854,0.000830481,0.000759108,0.000687736,0.000616363,0.00054499,0.000473618,0.000402245,0.000330872])
# BabesiosisU = 100*np.array([0.000833819,0.000751207,0.000668596,0.000585984,0.000634291,0.000682598,0.000730904,0.000779211,0.000814391,0.000849571,0.000884751,0.000919931,0.000997222,0.001074514,0.001151805,0.001229096,0.001306387,0.001602004,0.001897621,0.002193238,0.002488856,0.003091642,0.003694428,0.004297214,0.0049,0.004698791,0.004497583,0.004296374,0.004095165,0.003893956,0.003428215,0.002962473,0.002496732,0.00203099,0.00192545,0.00181991,0.00171437,0.00160883,0.001552962,0.001497094,0.001441226,0.001385358,0.00132949,0.001238652,0.001147814,0.001056976,0.000966138,0.0008753,0.000784462,0.000693624,0.000602786,0.000511948,0.00042111])
# WNVL = 100*np.array([0.0000000296953,0.0000000240391,0.0000000183828,0.0000000127266,0.0000000222715,0.0000000318164,0.0000000413614,0.0000000509063,0.0000000593907,0.000000067875,0.0000000763594,0.0000000848438,0.000000129387,0.00000017393,0.000000218473,0.000000263016,0.000000307559,0.000000608224,0.000000908889,0.00000120955,0.00000151022,0.00000495276,0.0000083953,0.0000118378,0.0000152804,0.0000210243,0.0000267682,0.0000325121,0.0000382561,0.000044,0.00004182,0.0000396401,0.0000374601,0.0000352802,0.0000288649,0.0000224497,0.0000160344,0.00000961917,0.00000793417,0.00000624917,0.00000456417,0.00000287917,0.00000119418,0.00000104358,0.000000892981,0.000000742383,0.000000591786,0.000000441188,0.00000029059,0.000000139992,0.00000014,0.00000014,0.00000014])
# WNV = 100*np.array([0.000000334073,0.00000027044,0.000000206807,0.000000143174,0.000000250554,0.000000357935,0.000000465315,0.000000572696,0.000000668145,0.000000763594,0.000000859044,0.000000954493,0.0000014556,0.00000195671,0.00000245782,0.00000295893,0.00000346004,0.00000684252,0.000010225,0.0000136075,0.00001699,0.0000557185,0.0000944471,0.000133176,0.000171904,0.000236523,0.000301142,0.000365762,0.000430381,0.000495,0.000470475,0.000445951,0.000421426,0.000396902,0.00032473,0.000252559,0.000180387,0.000108216,0.0000892594,0.0000703032,0.0000513469,0.0000323907,0.0000134345,0.0000117403,0.000010046,0.00000835181,0.00000665759,0.00000496336,0.00000326914,0.00000157491,0.00000157,0.00000157,0.00000157])
# WNVU = 100*np.array([0.00000101234,0.000000819514,0.000000626687,0.00000043386,0.000000759256,0.00000108465,0.00000141005,0.00000173544,0.00000202468,0.00000231392,0.00000260316,0.0000028924,0.00000441091,0.00000592943,0.00000744794,0.00000896645,0.000010485,0.0000207349,0.0000309849,0.0000412348,0.0000514848,0.000168844,0.000286203,0.000403562,0.000520922,0.000716737,0.000912553,0.001108369,0.001304184,0.0015,0.001425683,0.001351367,0.00127705,0.001202733,0.000984032,0.00076533,0.000546628,0.000327926,0.000270483,0.00021304,0.000155597,0.0000981537,0.0000407106,0.0000355766,0.0000304425,0.0000253085,0.0000201745,0.0000150405,0.00000990648,0.00000477246,0.00000477,0.00000477,0.00000477])
#
# #
# matplotlib.rc('text', usetex=True)
# matplotlib.rcParams.update({'font.size': 50})
#
# # plt.plot(x,HIVL,ls=' ')
# # plt.plot(x,HIV,color='k',label='$p_{1}(t)$:HIV',marker='^',linewidth=0.5)
# # plt.plot(x,HIVU,ls=' ',color='k')
#
# # plt.plot(x,HBVL,ls=' ')
# # plt.plot(x,HBV,color='k',label='$p_{2}(t)$:HBV',marker='+',linewidth=0.5)
# # plt.plot(x,HBVU,ls=' ',color='k')
# #
# # plt.plot(x,HCVL,ls=' ')
# # plt.plot(x,HCV,color='k',label='$p_{3}(t)$:HCV',marker='o',linewidth=0.5)
# # plt.plot(x,HCVU,ls=' ',color='k')
#
# dummy=[]
# for i in x:
#     dummy.append((BabesiosisL[i],BabesiosisU[i]))
#
# plt.plot((x,x),([i for (i,j) in dummy], [j for (i,j) in dummy]),c='black',ls=':')
# # plt.plot(x,BabesiosisL,ls=' ')
# plt.scatter(x,BabesiosisL,marker='_',s=200,c='k')
# # plt.plot(x,Babesiosis,color='k',label='Babesiosis',ls='--',dashes=[10,12],linewidth=0.8)
# plt.scatter(x,Babesiosis,color='k',label='Babesiosis',marker='x',s=100)
# # plt.plot(x,BabesiosisU,ls=' ',color='k')
# plt.scatter(x,BabesiosisU,marker='_',s=200,c='k')
#
#
# dummy=[]
# for i in x:
#     dummy.append((WNVL[i],WNVU[i]))
#
# plt.plot((x,x),([i for (i,j) in dummy], [j for (i,j) in dummy]),c='black',ls=':')
#
# # plt.plot(x,WNVL,ls=' ')
# plt.scatter(x,WNVL,marker='_',s=200,c='k')
# # plt.plot(x,WNV,color='k',label='WNV',linewidth=0.8)
# plt.scatter(x,WNV,color='k',label='WNV',s=30)
# # plt.plot(x,WNVU,ls=' ',color='k')
# plt.scatter(x,WNVU,marker='_',s=200,c='k')
#
# # plt.fill_between(x,HIVL,HIVU,alpha=0.15)
# # plt.fill_between(x,HBVL,HBVU,alpha=0.15)
# # plt.fill_between(x,HCVL,HCVU,alpha=0.15)
# plt.fill_between(x,BabesiosisL,BabesiosisU,alpha=0.3)
# plt.fill_between(x,WNVL,WNVU,alpha=0.3)
#
# plt.legend(loc='best',fontsize=45)
# plt.xlabel('Time (weeks)',fontsize=45)
# plt.xticks(np.arange(53), ('Jan ', ' ', ' ', ' ', 'Feb', ' ', ' ', ' ', 'Mar', ' ', ' ', ' ', 'Apr', ' ', ' ', ' ', ' ','May', ' ', ' ', ' ',' Jun', ' ', ' ', ' ', 'Jul', ' ', ' ', ' ',' ' ,'Aug', ' ', ' ', ' ', 'Sep', ' ', ' ', ' ', 'Oct', ' ', ' ', ' ', ' ' ,'Nov', ' ', ' ', ' ', 'Dec', ' ', ' ', ' ', ' ', ' '))
# plt.ylabel('Prevalence rate $(\%)$',fontsize=45)
# plt.show()










# # creation of matrices of prevalence rates used for app 2
# Matrix_of_prevalence_rates_of_infections=[[0]*53]*5
# Matrix_of_prevalence_rates_of_infections[0]=HIV
# Matrix_of_prevalence_rates_of_infections[1]=HBV
# Matrix_of_prevalence_rates_of_infections[2]=HCV
# Matrix_of_prevalence_rates_of_infections[3]=Babesiosis
# Matrix_of_prevalence_rates_of_infections[4]=WNV
#
# Prevalence_Data_List=Matrix_of_prevalence_rates_of_infections
#
# Matrix_of_prevalence_rates_of_infectionsL=[[0]*53]*5
# Matrix_of_prevalence_rates_of_infectionsL[0]=HIVL
# Matrix_of_prevalence_rates_of_infectionsL[1]=HBVL
# Matrix_of_prevalence_rates_of_infectionsL[2]=HCVL
# Matrix_of_prevalence_rates_of_infectionsL[3]=BabesiosisL
# Matrix_of_prevalence_rates_of_infectionsL[4]=WNVL
#
# Prevalence_Data_List_Lower=Matrix_of_prevalence_rates_of_infectionsL
#
# Matrix_of_prevalence_rates_of_infectionsU=[[0]*53]*5
# Matrix_of_prevalence_rates_of_infectionsU[0]=HIVU
# Matrix_of_prevalence_rates_of_infectionsU[1]=HBVU
# Matrix_of_prevalence_rates_of_infectionsU[2]=HCVU
# Matrix_of_prevalence_rates_of_infectionsU[3]=BabesiosisU
# Matrix_of_prevalence_rates_of_infectionsU[4]=WNVU
#
# Prevalence_Data_List_Upper=Matrix_of_prevalence_rates_of_infectionsU



T=50
x=np.arange(T)
Al,A,Au=Generate_a_p_t_function(T)


# x=x[::2]
# Al,A,Au=Al[::2],A[::2],Au[::2]
import matplotlib
matplotlib.rc('text', usetex=True)
matplotlib.rcParams.update({'font.size': 45})


dummy=[]
for i in x:
    dummy.append((Al[i],Au[i]))
plt.plot((x,x),([i for (i,j) in dummy], [j for (i,j) in dummy]),c='black',ls=':')

plt.scatter(x,Al,marker='_',s=200,c='k')
plt.scatter(x,A,color='k',marker='x',s=100)
plt.scatter(x,Au,marker='_',s=200,c='k')
plt.fill_between(x,Al,Au,alpha=0.3)
plt.yticks([])
plt.axvline(30,ls='-',linewidth=1.5,c='k',dashes=(3,4))
plt.axvline(40,ls='-',linewidth=1.5,c='k',dashes=(3,5))
plt.arrow(30-1, Al[31], 0.0, Au[40]-Al[31] -0.0005, fc="k", ec="k",head_width=0.5, head_length=0.0003)
plt.arrow(30-1, Au[39], 0.0, -Au[40]+Al[31] +0.00075, fc="k", ec="k",head_width=0.5, head_length=0.0003)
plt.text(24,0.017,'$\mathcal{P}_{30,40}$')
X,Y=[30,40],[Au[40],Au[40]]
plt.plot(X,Y,ls=":",linewidth=2,c='k')
X,Y=[30,40],[Al[30],Al[30]]
plt.plot(X,Y,ls=":",linewidth=2,c='k')
plt.xlim(-1,T+1)
plt.xlabel('time')
plt.show()