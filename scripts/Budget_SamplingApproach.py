import cvxpy as cp
import numpy as np
import matplotlib
from matplotlib.pylab import plt
import math
import random

# def ERM_Problem(Possible_Prevalence_Vector_Cases, k, BT, number_of_viruses):
#     Risk_values_of_Possible_Cases = []
#     for i in range(len(Possible_Prevalence_Vector_Cases)):
#         mu = Possible_Prevalence_Vector_Cases[i]
#         B_d = cp.Variable(number_of_viruses)
#         ERM_obj_func = sum(mu[j] * cp.exp(-k[j] * B_d[j]) for j in range(number_of_viruses))
#         constraints_ERM = [sum(B_d) <= BT, B_d >= 0]
#         ERM_prob = cp.Problem(cp.Minimize(ERM_obj_func), constraints_ERM)
#         ERM_prob.solve()
#         #        print("\nThe optimal risk is", ERM_prob.value)
#         #        print("The optimal B_d is")
#         #        print(B_d.value)
#         #        print(sum(B_d.value))
#         Risk_values_of_Possible_Cases.append(ERM_prob.value)
#     return Risk_values_of_Possible_Cases

def Threshold_needed_for_ERM(k,pvector,infectioni,IstarE):
    dummyarray=[(k[j]*pvector[j])**(1/k[j]) for j in IstarE]
    S=sum((1/k[j]) for j in IstarE)
    return math.log(np.product(dummyarray)) - S*math.log(k[infectioni]*pvector[infectioni])


def Budget_ERM_optimal(k,infectioni,IstarE,BT,pvector):
    S = sum((1 / k[j]) for j in IstarE)
    dummyarray=[(k[j]*pvector[j])**(1/(k[j]*S)) for j in IstarE]
    a=BT/(k[infectioni]*S)
    b=(1/k[infectioni])*math.log(k[infectioni]*pvector[infectioni])
    c=-(1/k[infectioni])*math.log(np.product(dummyarray))
    return a+b+c


def closed_form_formula_for_ERM(k,pvector,BT):
    B=[]
    dummy=[-k[m]*pvector[m]*math.exp(-k[m]*0) for m in range(number_of_viruses)]
    A=[]
    for i in range(number_of_viruses):
        A.append(dummy.index(min(dummy)))
        dummy[dummy.index(min(dummy))]=100000000000
    IstarE=[A[0]]
    j=1
    m=A[j]
    TH=Threshold_needed_for_ERM(k, pvector, m, IstarE)
    while TH<=BT:
        IstarE.append(m)
        if j==number_of_viruses-1:
            break
        j=j+1
        m=A[j]
        TH=Threshold_needed_for_ERM(k,pvector,m,IstarE)
    for m in range(number_of_viruses):
        if m in IstarE:
            B.append(Budget_ERM_optimal(k,m,IstarE,BT,pvector))
        else:
            B.append(0)
    objfunc=sum(pvector[j] * math.exp(-k[j] * B[j]) for j in range(number_of_viruses))
    return objfunc




HIVL = [0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005]
HIVU = [0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01]
HBVL = [0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025]
HBVU = [0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044]
HCVL = [0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013]
HCVU = [0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019]
BabesiosisL = [0.000433926,0.000390934,0.000347943,0.000304951,0.00033009,0.000355229,0.000380369,0.000405508,0.000423816,0.000442124,0.000460432,0.00047874,0.000518963,0.000559186,0.000599408,0.000639631,0.000679854,0.000833696,0.000987538,0.001141379,0.001295221,0.001608916,0.00192261,0.002236305,0.00255,0.002445289,0.002340579,0.002235868,0.002131157,0.002026447,0.001784071,0.001541695,0.00129932,0.001056944,0.00100202,0.000947096,0.000892172,0.000837248,0.000808174,0.0007791,0.000750026,0.000720952,0.000691877,0.000644605,0.000597332,0.000550059,0.000502786,0.000455513,0.00040824,0.000360968,0.000313695,0.000266422,0.000219149]
BabesiosisU = [0.000833819,0.000751207,0.000668596,0.000585984,0.000634291,0.000682598,0.000730904,0.000779211,0.000814391,0.000849571,0.000884751,0.000919931,0.000997222,0.001074514,0.001151805,0.001229096,0.001306387,0.001602004,0.001897621,0.002193238,0.002488856,0.003091642,0.003694428,0.004297214,0.0049,0.004698791,0.004497583,0.004296374,0.004095165,0.003893956,0.003428215,0.002962473,0.002496732,0.00203099,0.00192545,0.00181991,0.00171437,0.00160883,0.001552962,0.001497094,0.001441226,0.001385358,0.00132949,0.001238652,0.001147814,0.001056976,0.000966138,0.0008753,0.000784462,0.000693624,0.000602786,0.000511948,0.00042111]
WNVL = [0.0000000296953,0.0000000240391,0.0000000183828,0.0000000127266,0.0000000222715,0.0000000318164,0.0000000413614,0.0000000509063,0.0000000593907,0.000000067875,0.0000000763594,0.0000000848438,0.000000129387,0.00000017393,0.000000218473,0.000000263016,0.000000307559,0.000000608224,0.000000908889,0.00000120955,0.00000151022,0.00000495276,0.0000083953,0.0000118378,0.0000152804,0.0000210243,0.0000267682,0.0000325121,0.0000382561,0.000044,0.00004182,0.0000396401,0.0000374601,0.0000352802,0.0000288649,0.0000224497,0.0000160344,0.00000961917,0.00000793417,0.00000624917,0.00000456417,0.00000287917,0.00000119418,0.00000104358,0.000000892981,0.000000742383,0.000000591786,0.000000441188,0.00000029059,0.000000139992,0.00000014,0.00000014,0.00000014]
WNVU = [0.00000101234,0.000000819514,0.000000626687,0.00000043386,0.000000759256,0.00000108465,0.00000141005,0.00000173544,0.00000202468,0.00000231392,0.00000260316,0.0000028924,0.00000441091,0.00000592943,0.00000744794,0.00000896645,0.000010485,0.0000207349,0.0000309849,0.0000412348,0.0000514848,0.000168844,0.000286203,0.000403562,0.000520922,0.000716737,0.000912553,0.001108369,0.001304184,0.0015,0.001425683,0.001351367,0.00127705,0.001202733,0.000984032,0.00076533,0.000546628,0.000327926,0.000270483,0.00021304,0.000155597,0.0000981537,0.0000407106,0.0000355766,0.0000304425,0.0000253085,0.0000201745,0.0000150405,0.00000990648,0.00000477246,0.00000477,0.00000477,0.00000477]

dictio={'HIV': [HIVL,HIVU],
        'HBV': [HBVL,HBVU],
        'HCV': [HCVL,HCVU],
        'Babesiosis': [BabesiosisL,BabesiosisU],
        'WNV': [WNVL,WNVU]}


############## Below Prevalence list is just to address the comment of reviewer 1: Major Revision 2############
simulations=5000
Prevalence_Data_List_Simulated=[]
for reps in range(simulations):
    dummyalldiseases=[]
    for M in dictio:
        dummyonedisease=[]
        for i in range(len(HIVU)):
            dummyonedisease.append(np.random.uniform(dictio[M][0][i], dictio[M][1][i],1)[0])
            # dummyonedisease.append(np.random.uniform(random.uniform(0, 0.01) * dictio[M][0][i], random.uniform(0, 0.001) + dictio[M][1][i],1)[0])
        dummyalldiseases.append(dummyonedisease)
    Prevalence_Data_List_Simulated.append(dummyalldiseases)
##save the prev data to be read from the regret problem for a fair comparison
import pickle
with open('Prevalence_Data_ListSimulatedAPP2.data', 'wb') as f:
    pickle.dump(Prevalence_Data_List_Simulated, f)
f.close()




BT = 45
k = [0.2800,0.1600,0.1400,0.3800,0.1850]
r=0
number_of_viruses = 5
matrix_column = 2*number_of_viruses
N=len(HIVL)
T = np.arange(0,len(HIVL),1)


#Continuous
F_continuous_Simulations=[]
sum_Contin_Simulation, max_Contin_Simulation=[],[]
for reps in range(simulations):
    Prevalence_Data_List=Prevalence_Data_List_Simulated[reps]
    Transpose_of_Prevalence_Data_List = np.transpose(Prevalence_Data_List)
    # function_continuous_case = ERM_Problem(Transpose_of_Prevalence_Data_List, k, BT, number_of_viruses)
    function_continuous_case = [closed_form_formula_for_ERM(k,Transpose_of_Prevalence_Data_List[x,:],BT) for x in range(len(Transpose_of_Prevalence_Data_List))]
    sum_Contin_Simulation.append(sum(function_continuous_case))
    max_Contin_Simulation.append(max(function_continuous_case))
    F_continuous_Simulations.append(function_continuous_case)
    print('we are in reps', reps)


# K_allowable_changes_array=np.arange(1,5,1)
K_allowable_changes_array=np.arange(1,13,1)

#regret solution is below
Regret_F_constrained_matrix_of_ks_Simulations=[]
Regret_Total_Average_cost_to_save_Simulations, Regret_Max_cost_to_save_Simulations=[],[]
with open('solutionovertime_RegretAPP2.data', 'rb') as f:
    xstar_regret_solutionKs = pickle.load(f)
f.close()
for K_allowable_changes in K_allowable_changes_array:
    xstar_regret_solution=xstar_regret_solutionKs[K_allowable_changes-1]
    SimulatedFvalues=[]
    Total_Average_cost_to_save, Max_cost_to_save=[],[]
    for reps in range(simulations):
        Prevalence_Data_List = Prevalence_Data_List_Simulated[reps]
        Transpose_of_Prevalence_Data_List = np.transpose(Prevalence_Data_List)
        F_constrained_solution = []
        for i in range(len(HIVL)):
            nu = xstar_regret_solution[i]
            mu = Transpose_of_Prevalence_Data_List[i]
            F_constrained_solution.append(sum(mu[j] * math.exp(-k[j] * nu[j]) for j in range(number_of_viruses)))
        SimulatedFvalues.append(F_constrained_solution)
        Total_Average_cost_to_save.append(sum(F_constrained_solution))
        Max_cost_to_save.append(max(F_constrained_solution))
    # Regret_F_constrained_matrix_of_ks_Simulations.append(SimulatedFvalues)
    Regret_Total_Average_cost_to_save_Simulations.append(Total_Average_cost_to_save)
    Regret_Max_cost_to_save_Simulations.append(Max_cost_to_save)





#robust solution is below
Robust_F_constrained_matrix_of_ks_Simulations=[]
Robust_Total_Average_cost_to_save_Simulations, Robust_Max_cost_to_save_Simulations=[],[]
with open('solutionovertime_RobustAPP2.data', 'rb') as f:
    xstar_robust_solutionKs = pickle.load(f)
f.close()
for K_allowable_changes in K_allowable_changes_array:
    xstar_robust_solution=xstar_robust_solutionKs[K_allowable_changes-1]
    SimulatedFvalues=[]
    Total_Average_cost_to_save,Max_cost_to_save=[],[]
    for reps in range(simulations):
        Prevalence_Data_List = Prevalence_Data_List_Simulated[reps]
        Transpose_of_Prevalence_Data_List = np.transpose(Prevalence_Data_List)
        F_constrained_solution = []
        for i in range(len(HIVL)):
            nu = xstar_robust_solution[i]
            mu = Transpose_of_Prevalence_Data_List[i]
            F_constrained_solution.append(sum(mu[j] * math.exp(-k[j] * nu[j]) for j in range(number_of_viruses)))
        SimulatedFvalues.append(F_constrained_solution)
        Total_Average_cost_to_save.append(sum(F_constrained_solution))
        Max_cost_to_save.append(max(F_constrained_solution))
    # Robust_F_constrained_matrix_of_ks_Simulations.append(SimulatedFvalues)
    Robust_Total_Average_cost_to_save_Simulations.append(Total_Average_cost_to_save)
    Robust_Max_cost_to_save_Simulations.append(Max_cost_to_save)






# #stochastic solution is below
# Stochastic_F_constrained_matrix_of_ks_Simulations=[]
# Stochastic_Total_Average_cost_to_save_Simulations, Stochastic_Max_cost_to_save_Simulations=[],[]
# with open('solutionovertime_StochasticAPP2.data', 'rb') as f:
#     xstar_Stochastic_solutionKs = pickle.load(f)
# f.close()
# for K_allowable_changes in K_allowable_changes_array:
#     xstar_Stochastic_solution=xstar_Stochastic_solutionKs[K_allowable_changes-1]
#     SimulatedFvalues=[]
#     Total_Average_cost_to_save,Max_cost_to_save=[],[]
#     for reps in range(simulations):
#         Prevalence_Data_List = Prevalence_Data_List_Simulated[reps]
#         Transpose_of_Prevalence_Data_List = np.transpose(Prevalence_Data_List)
#         F_constrained_solution = []
#         for i in range(len(HIVL)):
#             nu = xstar_Stochastic_solution[i]
#             mu = Transpose_of_Prevalence_Data_List[i]
#             F_constrained_solution.append(sum(mu[j] * math.exp(-k[j] * nu[j]) for j in range(number_of_viruses)))
#         SimulatedFvalues.append(F_constrained_solution)
#         Total_Average_cost_to_save.append(sum(F_constrained_solution))
#         Max_cost_to_save.append(max(F_constrained_solution))
#     # Stochastic_F_constrained_matrix_of_ks_Simulations.append(SimulatedFvalues)
#     Stochastic_Total_Average_cost_to_save_Simulations.append(Total_Average_cost_to_save)
#     Stochastic_Max_cost_to_save_Simulations.append(Max_cost_to_save)








K_tobeused=3 #########Spcify the K value that needs to be used in plotting

price_of_reg_relative_to_ContSUM=100*abs(np.mean(Regret_Total_Average_cost_to_save_Simulations[K_tobeused]) - np.mean(sum_Contin_Simulation))/np.mean(sum_Contin_Simulation)
price_of_rob_relative_to_ContSUM=100*abs(np.mean(Robust_Total_Average_cost_to_save_Simulations[K_tobeused]) - np.mean(sum_Contin_Simulation))/np.mean(sum_Contin_Simulation)
# price_of_stoch_relative_to_ContSUM=100*abs(np.mean(Stochastic_Total_Average_cost_to_save_Simulations[K_tobeused]) - np.mean(sum_Contin_Simulation))/np.mean(sum_Contin_Simulation)

price_of_reg_relative_to_ContMAX=100*abs(np.mean(Regret_Max_cost_to_save_Simulations[K_tobeused]) - np.mean(max_Contin_Simulation))/np.mean(max_Contin_Simulation)
price_of_rob_relative_to_ContMAX=100*abs(np.mean(Robust_Max_cost_to_save_Simulations[K_tobeused]) - np.mean(max_Contin_Simulation))/np.mean(max_Contin_Simulation)
# price_of_stoch_relative_to_ContMAX=100*abs(np.mean(Stochastic_Max_cost_to_save_Simulations[K_tobeused]) - np.mean(max_Contin_Simulation))/np.mean(max_Contin_Simulation)


######################## PLOTS OF DISTRIBUTION FOR SUM AND MAX FOR THE ALL STRATEGIES ########################## K=[1]
costtomultiplyapp2=16000*(13.6/53)/(26)
matplotlib.rc('text', usetex=True)
matplotlib.rcParams.update({'font.size': 28})
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]

sum_Contin_Simulation=np.array(sum_Contin_Simulation)*costtomultiplyapp2
max_Contin_Simulation=np.array(max_Contin_Simulation)*costtomultiplyapp2
Regret_Total_Average_cost_to_save_Simulations=np.array(Regret_Total_Average_cost_to_save_Simulations)*costtomultiplyapp2
Regret_Max_cost_to_save_Simulations=np.array(Regret_Max_cost_to_save_Simulations)*costtomultiplyapp2
Robust_Total_Average_cost_to_save_Simulations=np.array(Robust_Total_Average_cost_to_save_Simulations)*costtomultiplyapp2
Robust_Max_cost_to_save_Simulations=np.array(Robust_Max_cost_to_save_Simulations)*costtomultiplyapp2
# Stochastic_Total_Average_cost_to_save_Simulations=np.array(Stochastic_Total_Average_cost_to_save_Simulations)*costtomultiplyapp2
# Stochastic_Max_cost_to_save_Simulations=np.array(Stochastic_Max_cost_to_save_Simulations)*costtomultiplyapp2


fig, ax = plt.subplots(1,2, figsize=(14,5))
ax[0].hist((sum_Contin_Simulation), bins='auto', label='Nonrestricted', alpha=0.4)
ax[0].axvline(x=np.mean(sum_Contin_Simulation), ls='--', color='b', linewidth=2)
ax[0].hist((Regret_Total_Average_cost_to_save_Simulations[K_tobeused]), bins='auto', label='Regret', alpha=0.4)
ax[0].axvline(x=np.mean(Regret_Total_Average_cost_to_save_Simulations[K_tobeused]), ls='--', color='orange', linewidth=2)
ax[0].hist((Robust_Total_Average_cost_to_save_Simulations[K_tobeused]), bins='auto', label='Minimax', alpha=0.4)
ax[0].axvline(x=np.mean(Robust_Total_Average_cost_to_save_Simulations[K_tobeused]), ls='--', color='g', linewidth=2)
# ax[0].hist((Stochastic_Total_Average_cost_to_save_Simulations[K_tobeused]), bins='auto', label='Stochastic', alpha=0.4)
# ax[0].axvline(x=np.mean(Stochastic_Total_Average_cost_to_save_Simulations[K_tobeused]), ls='--', color='r', linewidth=2)
ax[0].legend(ncol=1, fontsize=20)
ax[0].set_xlabel(r'$\rm{Total \; cost \;(in \; \$1M)}$')
ax[0].set_ylabel(r'\rm{Distribution}')

ax[1].hist((max_Contin_Simulation), bins='auto', label='Nonrestricted', alpha=0.4)
ax[1].axvline(x=np.mean(max_Contin_Simulation), ls='--', color='b', linewidth=2)
ax[1].hist((Regret_Max_cost_to_save_Simulations[K_tobeused]), bins='auto', label='Regret', alpha=0.4)
ax[1].axvline(x=np.mean(Regret_Max_cost_to_save_Simulations[K_tobeused]), ls='--', color='orange', linewidth=2)
ax[1].hist((Robust_Max_cost_to_save_Simulations[K_tobeused]), bins='auto', label='Minimax', alpha=0.4)
ax[1].axvline(x=np.mean(Robust_Max_cost_to_save_Simulations[K_tobeused]), ls='--', color='g', linewidth=2)
# ax[1].hist((Stochastic_Max_cost_to_save_Simulations[K_tobeused]), bins='auto', label='Stochastic', alpha=0.4)
# ax[1].axvline(x=np.mean(Stochastic_Max_cost_to_save_Simulations[K_tobeused]), ls='--', color='r', linewidth=2)
ax[1].legend(ncol=1, fontsize=20)
ax[1].set_xlabel(r'$\rm{Max \; cost \;(in \; \$1M)}$')
ax[1].set_ylabel(r'\rm{Distribution}')
######################### PLOTS OF DISTRIBUTION FOR SUM AND MAX FOR THE ALL STRATEGIES ########################## K=[1]













########################## PLOTS OF VARIATION WRT K ########################## K=[0,1,2,3,4,5,6,7,8,9,10,11,12]
costtomultiplyapp2=16000*(13.6/53)/(26)
matplotlib.rc('text', usetex=True)
matplotlib.rcParams.update({'font.size': 30})
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
sum_Contin_Simulation=np.array(sum_Contin_Simulation)*costtomultiplyapp2 #########IF NOT MULTIPLIED ALREADY IN THE FIRST FIGURE, OTHERWISE WE ARE DOING THAT TWICE
Regret_Total_Average_cost_to_save_Simulations=np.array(Regret_Total_Average_cost_to_save_Simulations)*costtomultiplyapp2
Robust_Total_Average_cost_to_save_Simulations=np.array(Robust_Total_Average_cost_to_save_Simulations)*costtomultiplyapp2
# Stochastic_Total_Average_cost_to_save_Simulations=np.array(Stochastic_Total_Average_cost_to_save_Simulations)*costtomultiplyapp2
meansReg=[np.mean(Regret_Total_Average_cost_to_save_Simulations[x]) for x in range(len(K_allowable_changes_array))]
maxsReg=[np.max(Regret_Total_Average_cost_to_save_Simulations[x]) for x in range(len(K_allowable_changes_array))]
minsReg=[np.min(Regret_Total_Average_cost_to_save_Simulations[x]) for x in range(len(K_allowable_changes_array))]
meansRob=[np.mean(Robust_Total_Average_cost_to_save_Simulations[x]) for x in range(len(K_allowable_changes_array))]
maxsRob=[np.max(Robust_Total_Average_cost_to_save_Simulations[x]) for x in range(len(K_allowable_changes_array))]
minsRob=[np.min(Robust_Total_Average_cost_to_save_Simulations[x]) for x in range(len(K_allowable_changes_array))]
# meansStoch=[np.mean(Stochastic_Total_Average_cost_to_save_Simulations[x]) for x in range(len(K_allowable_changes_array))]
# maxsStoch=[np.max(Stochastic_Total_Average_cost_to_save_Simulations[x]) for x in range(len(K_allowable_changes_array))]
# minsStoch=[np.min(Stochastic_Total_Average_cost_to_save_Simulations[x]) for x in range(len(K_allowable_changes_array))]


fig, ax = plt.subplots(1,1, figsize=(10,6))
ax.plot([x for x in range(len(K_allowable_changes_array))],meansReg, label='Regret', color='orange', ls=':')
ax.scatter([x for x in range(len(K_allowable_changes_array))],meansReg, marker='D',s=50, color='orange')
ax.plot([x for x in range(len(K_allowable_changes_array))],maxsReg, color='orange')
ax.scatter([x for x in range(len(K_allowable_changes_array))],maxsReg, marker='*',s=75, color='orange')
ax.plot([x for x in range(len(K_allowable_changes_array))],minsReg, color='orange')
ax.scatter([x for x in range(len(K_allowable_changes_array))],minsReg, marker='*',s=75, color='orange')
ax.fill_between([x for x in range(len(K_allowable_changes_array))],minsReg,maxsReg,alpha=0.2, color='orange')

ax.plot([x for x in range(len(K_allowable_changes_array))],meansRob, label='Minimax', color='green', ls=':')
ax.scatter([x for x in range(len(K_allowable_changes_array))],meansRob, marker='D',s=50, color='green')
ax.plot([x for x in range(len(K_allowable_changes_array))],maxsRob, color='green')
ax.scatter([x for x in range(len(K_allowable_changes_array))],maxsRob, marker='*',s=75, color='green')
ax.plot([x for x in range(len(K_allowable_changes_array))],minsRob, color='green')
ax.scatter([x for x in range(len(K_allowable_changes_array))],minsRob, marker='*',s=75, color='green')
ax.fill_between([x for x in range(len(K_allowable_changes_array))],minsRob,maxsRob,alpha=0.2, color='green')

# ax.plot([x for x in range(len(K_allowable_changes_array))],meansStoch, label='Stochastic', color='r', ls=':')
# ax.scatter([x for x in range(len(K_allowable_changes_array))],meansStoch, marker='D',s=50, color='r')
# ax.plot([x for x in range(len(K_allowable_changes_array))],maxsStoch, color='r')
# ax.scatter([x for x in range(len(K_allowable_changes_array))],maxsStoch, marker='*',s=75, color='r')
# ax.plot([x for x in range(len(K_allowable_changes_array))],minsStoch, color='r')
# ax.scatter([x for x in range(len(K_allowable_changes_array))],minsStoch, marker='*',s=75, color='r')
# ax.fill_between([x for x in range(len(K_allowable_changes_array))],minsStoch,maxsStoch,alpha=0.2, color='r')
ax.legend(fontsize=28)
ax.set_ylabel(r'$\rm{Total\; cost \; (in \; \$1M) }$')
ax.set_xlabel(r'$\rm{Allowable \; number\; of \; changes\; (K)}$')
########################## PLOTS OF VARIATION WRT K ########################## K=[0,1,2,3,4,5,6,7,8,9,10,11,12]




















# ax[0,0].hist((sum_Contin_Simulation), bins='auto', label='Nonrestricted')
# ax[0,0].axvline(x=np.mean(sum_Contin_Simulation), ls='--', color='b', linewidth=2)
# ax[0,0].hist((Regret_Total_Average_cost_to_save_Simulations[0]), bins='auto', label='Regret')
# ax[0,0].axvline(x=np.mean(Regret_Total_Average_cost_to_save_Simulations[0]), ls='--', color='orange', linewidth=2)
# ax[0,0].hist((Robust_Total_Average_cost_to_save_Simulations[0]), bins='auto', label='Minmax')
# ax[0,0].axvline(x=np.mean(Robust_Total_Average_cost_to_save_Simulations[0]), ls='--', color='g', linewidth=2)
# ax[0,0].legend(ncol=1, fontsize=20)
# ax[0,0].set_xlabel(r'$\rm{Cost \;(in \; \$1M)}$')
# ax[0,0].set_ylabel(r'\rm{Distribution}')
#
#
# ax[0,1].hist((sum_Contin_Simulation), bins='auto', label='Nonrestricted')
# ax[0,1].axvline(x=np.mean(sum_Contin_Simulation), ls='--', color='b', linewidth=2)
# ax[0,1].hist((Regret_Total_Average_cost_to_save_Simulations[1]), bins='auto', label='Regret')
# ax[0,1].axvline(x=np.mean(Regret_Total_Average_cost_to_save_Simulations[1]), ls='--', color='orange', linewidth=2)
# ax[0,1].hist((Robust_Total_Average_cost_to_save_Simulations[1]), bins='auto', label='Minmax')
# ax[0,1].axvline(x=np.mean(Robust_Total_Average_cost_to_save_Simulations[1]), ls='--', color='g', linewidth=2)
# ax[0,1].legend(ncol=1, fontsize=20)
# ax[0,1].set_xlabel(r'$\rm{Cost \;(in \; \$1M)}$')
# ax[0,1].set_ylabel(r'\rm{Distribution}')
#
# ax[1,0].hist((sum_Contin_Simulation), bins='auto', label='Nonrestricted')
# ax[1,0].axvline(x=np.mean(sum_Contin_Simulation), ls='--', color='b', linewidth=2)
# ax[1,0].hist((Regret_Total_Average_cost_to_save_Simulations[2]), bins='auto', label='Regret')
# ax[1,0].axvline(x=np.mean(Regret_Total_Average_cost_to_save_Simulations[2]), ls='--', color='orange', linewidth=2)
# ax[1,0].hist((Robust_Total_Average_cost_to_save_Simulations[2]), bins='auto', label='Minmax')
# ax[1,0].axvline(x=np.mean(Robust_Total_Average_cost_to_save_Simulations[2]), ls='--', color='g', linewidth=2)
# ax[1,0].legend(ncol=1, fontsize=20)
# ax[1,0].set_xlabel(r'$\rm{Cost \;(in \; \$1M)}$')
# ax[1,0].set_ylabel(r'\rm{Distribution}')
#
# ax[1,1].hist((sum_Contin_Simulation), bins='auto', label='Nonrestricted')
# ax[1,1].axvline(x=np.mean(sum_Contin_Simulation), ls='--', color='b', linewidth=2)
# ax[1,1].hist((Regret_Total_Average_cost_to_save_Simulations[3]), bins='auto', label='Regret')
# ax[1,1].axvline(x=np.mean(Regret_Total_Average_cost_to_save_Simulations[3]), ls='--', color='orange', linewidth=2)
# ax[1,1].hist((Robust_Total_Average_cost_to_save_Simulations[3]), bins='auto', label='Minmax')
# ax[1,1].axvline(x=np.mean(Robust_Total_Average_cost_to_save_Simulations[3]), ls='--', color='g', linewidth=2)
# ax[1,1].legend(ncol=1, fontsize=20)
# ax[1,1].set_xlabel(r'$\rm{Cost \;(in \; \$1M)}$')
# ax[1,1].set_ylabel(r'\rm{Distribution}')