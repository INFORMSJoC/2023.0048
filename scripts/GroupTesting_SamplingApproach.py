from scipy.special import lambertw
import math
from statistics import mean
import numpy as np
import matplotlib
from matplotlib.pylab import plt

absilun=0.0000000001
r8_big = 1.0E+30
# lambda_2 = 5000
lambda_2 = 60 # this 60 is nothing but the 615 we had for all TTIS but for WNV only


def n_0(p,Se,Sp):
    y=(2/(math.log(1-p)))*(lambertw(-0.5*(math.sqrt((math.log(1-p))/((-Se-Sp+1)*(1+lambda_1-lambda_1*Sp)))),0))
    return y

def f_n_0(p,n,Se,Sp):
    return (1/n)+Se-(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*(1-p)**n + (lambda_2*p*(1-Se*Se)) + (lambda_1*(1-Sp)*Se*(1-p))




Prevalence_Data_List_Upper=[3.74784e-05, 0.00010418250000000001, 0.0001753244, 0.000252516, 0.00033735800000000006, 0.0004758535, 0.0006140085, 0.0007555465, 0.000887538, 0.0009975000000000001, 0.0009586505000000001, 0.00090012, 0.0008273110000000001, 0.0007487475000000001, 0.0005999664999999999, 0.000454477, 0.00032105245, 0.00019465364999999999, 0.0001547977, 0.00011996985000000001]
Prevalence_Data_List_Lower=[1.249111e-05, 2.2236880000000002e-05, 3.642055e-05, 5.66539e-05, 8.453720000000001e-05, 0.00012799715, 0.0001711161, 0.00021761804999999998, 0.00025457405, 0.0002695, 0.000266719, 0.00024425655, 0.00020751605, 0.0001650211, 0.00012238295, 8.303685e-05, 5.575565e-05, 3.5500235e-05, 2.3523285e-05, 1.6574435e-05]

simulations=5000
############## Below Prevalence list is just to address the comment of reviewer 1: Major Revision 2############
Sampled_Prevalence_Data_List = [[np.random.uniform(Prevalence_Data_List_Lower[m], Prevalence_Data_List_Upper[m]) for m in range(len(Prevalence_Data_List_Lower))] for reps in range(simulations)]
##save the prev data to be read from the regret problem for a fair comparison
import pickle
with open('Prevalence_Data_ListSimulatedAPP1.data', 'wb') as f:
    pickle.dump(Sampled_Prevalence_Data_List, f)
f.close()


T = np.arange(0,len(Sampled_Prevalence_Data_List[0]),1)
Se=0.95
Sp=0.95
lambda_1=30
p_underline = 1 - math.exp(-(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*0.367879)
p_hat = 1 - math.exp(-(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*0.541341)

N=len(Sampled_Prevalence_Data_List[0])
# K_allowable_changes_array=np.arange(1,5,1)
K_allowable_changes_array=np.arange(1,7,1) # needed for the plot for varying over K, total cost along with variations



#16 and continous solution is below
n_star_continuous_Simulations = [[n_0(i, Se, Sp) for i in Sampled_Prevalence_Data_List[m]] for m in range(simulations)]
MinMaxDistance=[]
Total_Average_distance=[]
Total_Average_cost_to_save=[]
F_of_n_equal_16_Simulations=[]
F_of_n_star_continuous_Simulations=[]
sum_Contin_Simulation, max_Contin_Simulation=[],[]
sum_16_Simulation, max_16_Simulation=[],[]
for reps in range(simulations):
    F_of_n_equal_16=[]
    F_of_n_star_continuous=[]
    for i in range(len(Sampled_Prevalence_Data_List[reps])):
        ceil = math.ceil(n_star_continuous_Simulations[reps][i])
        floor = math.floor(n_star_continuous_Simulations[reps][i])
        min_between_ceil_and_floor = min(f_n_0(Sampled_Prevalence_Data_List[reps][i], ceil, Se, Sp),f_n_0(Sampled_Prevalence_Data_List[reps][i], floor, Se, Sp))
        F_of_n_star_continuous.append(min_between_ceil_and_floor)
        F_of_n_equal_16.append(f_n_0(Sampled_Prevalence_Data_List[reps][i],16,Se,Sp))
    F_of_n_equal_16_Simulations.append(F_of_n_equal_16)
    sum_Contin_Simulation.append(sum(F_of_n_star_continuous))
    max_Contin_Simulation.append(max(F_of_n_star_continuous))
    sum_16_Simulation.append(sum(F_of_n_equal_16))
    max_16_Simulation.append(max(F_of_n_equal_16))
    F_of_n_star_continuous_Simulations.append(F_of_n_star_continuous)
F_of_n_star_continuous_Simulations=np.array(F_of_n_star_continuous_Simulations)
# F_of_n_star_continuous_Simulations_OverTime_Mean=[np.mean(F_of_n_star_continuous_Simulations[:, m]) for m in T]
# F_of_n_star_continuous_Simulations_OverTime_Max=[np.max(F_of_n_star_continuous_Simulations[:, m]) for m in T]
# F_of_n_star_continuous_Simulations_OverTime_Min=[np.min(F_of_n_star_continuous_Simulations[:, m]) for m in T]

F_of_n_equal_16_Simulations=np.array(F_of_n_equal_16_Simulations)
# F_of_n_equal_16_Simulations_Simulations_OverTime_Mean=[np.mean(F_of_n_equal_16_Simulations[:, m]) for m in T]
# F_of_n_equal_16_Simulations_Simulations_OverTime_Max=[np.max(F_of_n_equal_16_Simulations[:, m]) for m in T]
# F_of_n_equal_16_Simulations_Simulations_OverTime_Min=[np.min(F_of_n_equal_16_Simulations[:, m]) for m in T]



#regret solution is below
Regret_F_constrained_matrix_of_ks_Simulations=[]
Regret_Total_Average_cost_to_save_Simulations, Regret_Max_cost_to_save_Simulations=[],[]
with open('solutionovertime_RegretAPP1.data', 'rb') as f:
    n_star_constrained_solutionKs = pickle.load(f)
f.close()
for K_allowable_changes in K_allowable_changes_array:
    n_star_constrained_solution=n_star_constrained_solutionKs[K_allowable_changes-1]
    SimulatedFvalues=[]
    Total_Average_cost_to_save, Max_cost_to_save=[],[]
    for reps in range(simulations):
        F_of_n_star_constrained_solution = []
        for i in range(len(Sampled_Prevalence_Data_List[0])):
            F_of_n_star_constrained_solution.append(f_n_0(Sampled_Prevalence_Data_List[reps][i],n_star_constrained_solution[i],Se,Sp))
        SimulatedFvalues.append(F_of_n_star_constrained_solution)
        Total_Average_cost_to_save.append(sum(F_of_n_star_constrained_solution))
        Max_cost_to_save.append(max(F_of_n_star_constrained_solution))
    Regret_F_constrained_matrix_of_ks_Simulations.append(SimulatedFvalues)
    Regret_Total_Average_cost_to_save_Simulations.append(Total_Average_cost_to_save)
    Regret_Max_cost_to_save_Simulations.append(Max_cost_to_save)
Regret_F_constrained_matrix_of_ks_Simulations=np.array(Regret_F_constrained_matrix_of_ks_Simulations)
# Regret_F_constrained_matrix_of_ks_Simulations_OverTime_Mean=[[np.mean(Regret_F_constrained_matrix_of_ks_Simulations[z-1,:, m]) for m in T] for z in K_allowable_changes_array]
# Regret_F_constrained_matrix_of_ks_Simulations_OverTime_Max=[[np.max(Regret_F_constrained_matrix_of_ks_Simulations[z-1,:, m]) for m in T] for z in K_allowable_changes_array]
# Regret_F_constrained_matrix_of_ks_Simulations_OverTime_Min=[[np.min(Regret_F_constrained_matrix_of_ks_Simulations[z-1,:, m]) for m in T] for z in K_allowable_changes_array]


#robust solution is below
Robust_F_constrained_matrix_of_ks_Simulations=[]
Robust_Total_Average_cost_to_save_Simulations, Robust_Max_cost_to_save_Simulations=[],[]
with open('solutionovertime_RobustAPP1.data', 'rb') as f:
    n_star_constrained_solutionKs = pickle.load(f)
f.close()
for K_allowable_changes in K_allowable_changes_array:
    n_star_constrained_solution=n_star_constrained_solutionKs[K_allowable_changes-1]
    SimulatedFvalues=[]
    Total_Average_cost_to_save,Max_cost_to_save=[],[]
    for reps in range(simulations):
        F_of_n_star_constrained_solution = []
        for i in range(len(Sampled_Prevalence_Data_List[0])):
            F_of_n_star_constrained_solution.append(f_n_0(Sampled_Prevalence_Data_List[reps][i],n_star_constrained_solution[i],Se,Sp))
        SimulatedFvalues.append(F_of_n_star_constrained_solution)
        Total_Average_cost_to_save.append(sum(F_of_n_star_constrained_solution))
        Max_cost_to_save.append(max(F_of_n_star_constrained_solution))
    Robust_F_constrained_matrix_of_ks_Simulations.append(SimulatedFvalues)
    Robust_Total_Average_cost_to_save_Simulations.append(Total_Average_cost_to_save)
    Robust_Max_cost_to_save_Simulations.append(Max_cost_to_save)
Robust_F_constrained_matrix_of_ks_Simulations=np.array(Robust_F_constrained_matrix_of_ks_Simulations)
# Robust_F_constrained_matrix_of_ks_Simulations_OverTime_Mean=[[np.mean(Robust_F_constrained_matrix_of_ks_Simulations[z-1,:, m]) for m in T] for z in K_allowable_changes_array]
# Robust_F_constrained_matrix_of_ks_Simulations_OverTime_Max=[[np.max(Robust_F_constrained_matrix_of_ks_Simulations[z-1,:, m]) for m in T] for z in K_allowable_changes_array]
# Robust_F_constrained_matrix_of_ks_Simulations_OverTime_Min=[[np.min(Robust_F_constrained_matrix_of_ks_Simulations[z-1,:, m]) for m in T] for z in K_allowable_changes_array]




# #stochatic solution is below
# Stochastic_F_constrained_matrix_of_ks_Simulations=[]
# Stochastic_Total_Average_cost_to_save_Simulations, Stochastic_Max_cost_to_save_Simulations=[],[]
# with open('solutionovertime_StochasticAPP1.data', 'rb') as f:
#     n_star_constrained_solutionKs = pickle.load(f)
# f.close()
# for K_allowable_changes in K_allowable_changes_array:
#     n_star_constrained_solution=n_star_constrained_solutionKs[K_allowable_changes-1]
#     SimulatedFvalues=[]
#     Total_Average_cost_to_save,Max_cost_to_save=[],[]
#     for reps in range(simulations):
#         F_of_n_star_constrained_solution = []
#         for i in range(len(Sampled_Prevalence_Data_List[0])):
#             F_of_n_star_constrained_solution.append(f_n_0(Sampled_Prevalence_Data_List[reps][i],n_star_constrained_solution[i],Se,Sp))
#         SimulatedFvalues.append(F_of_n_star_constrained_solution)
#         Total_Average_cost_to_save.append(sum(F_of_n_star_constrained_solution))
#         Max_cost_to_save.append(max(F_of_n_star_constrained_solution))
#     Stochastic_F_constrained_matrix_of_ks_Simulations.append(SimulatedFvalues)
#     Stochastic_Total_Average_cost_to_save_Simulations.append(Total_Average_cost_to_save)
#     Stochastic_Max_cost_to_save_Simulations.append(Max_cost_to_save)
# Stochatic_F_constrained_matrix_of_ks_Simulations=np.array(Stochastic_F_constrained_matrix_of_ks_Simulations)
# Robust_F_constrained_matrix_of_ks_Simulations_OverTime_Mean=[[np.mean(Robust_F_constrained_matrix_of_ks_Simulations[z-1,:, m]) for m in T] for z in K_allowable_changes_array]
# Robust_F_constrained_matrix_of_ks_Simulations_OverTime_Max=[[np.max(Robust_F_constrained_matrix_of_ks_Simulations[z-1,:, m]) for m in T] for z in K_allowable_changes_array]
# Robust_F_constrained_matrix_of_ks_Simulations_OverTime_Min=[[np.min(Robust_F_constrained_matrix_of_ks_Simulations[z-1,:, m]) for m in T] for z in K_allowable_changes_array]








K_tobeused=3 #########Spcify the K value that needs to be used in plotting

price_of_reg_relative_to_ContSUM=100*abs(np.mean(Regret_Total_Average_cost_to_save_Simulations[K_tobeused]) - np.mean(sum_Contin_Simulation))/np.mean(sum_Contin_Simulation)
price_of_rob_relative_to_ContSUM=100*abs(np.mean(Robust_Total_Average_cost_to_save_Simulations[K_tobeused]) - np.mean(sum_Contin_Simulation))/np.mean(sum_Contin_Simulation)
price_of_CP_relative_to_ContSUM=100*abs(np.mean(sum_16_Simulation) - np.mean(sum_Contin_Simulation))/np.mean(sum_Contin_Simulation)
# price_of_stoch_relative_to_ContSUM=100*abs(np.mean(Stochastic_Total_Average_cost_to_save_Simulations[K_tobeused]) - np.mean(sum_Contin_Simulation))/np.mean(sum_Contin_Simulation)

price_of_reg_relative_to_ContMAX=100*abs(np.mean(Regret_Max_cost_to_save_Simulations[K_tobeused]) - np.mean(max_Contin_Simulation))/np.mean(max_Contin_Simulation)
price_of_rob_relative_to_ContMAX=100*abs(np.mean(Robust_Max_cost_to_save_Simulations[K_tobeused]) - np.mean(max_Contin_Simulation))/np.mean(max_Contin_Simulation)
price_of_CP_relative_to_ContMAX=100*abs(np.mean(max_16_Simulation) - np.mean(max_Contin_Simulation))/np.mean(max_Contin_Simulation)
# price_of_stoch_relative_to_ContMAX=100*abs(np.mean(Stochastic_Max_cost_to_save_Simulations[K_tobeused]) - np.mean(max_Contin_Simulation))/np.mean(max_Contin_Simulation)


########################## PLOTS OF DISTRIBUTION FOR SUM AND MAX FOR THE ALL STRATEGIES ########################## K=[3]
costtomultiply=13.6/len(T) #(nb of donors)
#########PLOTTING GOES HERE########## #Histograms for simulated values of the sum#
matplotlib.rc('text', usetex=True)
matplotlib.rcParams.update({'font.size': 26})
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
fig, ax = plt.subplots(1,2, figsize=(14,5))
sum_Contin_Simulation=np.array(sum_Contin_Simulation)*costtomultiply
max_Contin_Simulation=np.array(max_Contin_Simulation)*costtomultiply
Regret_Total_Average_cost_to_save_Simulations=np.array(Regret_Total_Average_cost_to_save_Simulations)*costtomultiply
Regret_Max_cost_to_save_Simulations=np.array(Regret_Max_cost_to_save_Simulations)*costtomultiply
Robust_Total_Average_cost_to_save_Simulations=np.array(Robust_Total_Average_cost_to_save_Simulations)*costtomultiply
Robust_Max_cost_to_save_Simulations=np.array(Robust_Max_cost_to_save_Simulations)*costtomultiply
# Stochastic_Total_Average_cost_to_save_Simulations=np.array(Stochastic_Total_Average_cost_to_save_Simulations)*costtomultiply
# Stochastic_Max_cost_to_save_Simulations=np.array(Stochastic_Max_cost_to_save_Simulations)*costtomultiply
sum_16_Simulation=np.array(sum_16_Simulation)*costtomultiply
max_16_Simulation=np.array(max_16_Simulation)*costtomultiply



ax[0].hist((sum_Contin_Simulation), bins='auto', label='Nonrestricted', alpha=0.4)
ax[0].axvline(x=np.mean(sum_Contin_Simulation), ls='--', color='b', linewidth=2)
ax[0].hist((Regret_Total_Average_cost_to_save_Simulations[K_tobeused]), bins='auto', label='Regret', alpha=0.4)
ax[0].axvline(x=np.mean(Regret_Total_Average_cost_to_save_Simulations[K_tobeused]), ls='--', color='orange', linewidth=2)
ax[0].hist((Robust_Total_Average_cost_to_save_Simulations[K_tobeused]), bins='auto', label='Minimax', alpha=0.4)
ax[0].axvline(x=np.mean(Robust_Total_Average_cost_to_save_Simulations[K_tobeused]), ls='--', color='g', linewidth=2)
# ax[0].hist((Stochastic_Total_Average_cost_to_save_Simulations[K_tobeused]), bins='auto', label='Stochastic', alpha=0.4)
# ax[0].axvline(x=np.mean(Stochastic_Total_Average_cost_to_save_Simulations[K_tobeused]), ls='--', color='r', linewidth=2)
ax[0].hist((sum_16_Simulation), bins='auto', label='CP', alpha=0.4)
ax[0].axvline(x=np.mean(sum_16_Simulation), ls='--', color='r', linewidth=2)
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
ax[1].hist((max_16_Simulation), bins='auto', label='CP', alpha=0.4)
ax[1].axvline(x=np.mean(max_16_Simulation), ls='--', color='r', linewidth=2)
ax[1].legend(ncol=1, fontsize=20)
ax[1].set_xlabel(r'$\rm{Max \; cost \;(in \; \$1M)}$')
ax[1].set_ylabel(r'\rm{Distribution}')
########################## PLOTS OF DISTRIBUTION FOR SUM AND MAX FOR THE ALL STRATEGIES ########################## K=[0,1,2,3]














########################## PLOTS OF VARIATION WRT K ########################## K=[0,1,2,3,4,5,6]
costtomultiply=13.6/len(T) #(nb of donors)
#########PLOTTING GOES HERE########## #Histograms for simulated values of the sum#
matplotlib.rc('text', usetex=True)
matplotlib.rcParams.update({'font.size': 25})
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
Regret_Total_Average_cost_to_save_Simulations=np.array(Regret_Total_Average_cost_to_save_Simulations)*costtomultiply
Robust_Total_Average_cost_to_save_Simulations=np.array(Robust_Total_Average_cost_to_save_Simulations)*costtomultiply
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
ax.legend(ncol=2,fontsize=28)
ax.set_ylabel(r'$\rm{Total\; cost \; (in \; \$1M) }$')
ax.set_xlabel(r'$\rm{Allowable \; number\; of \; changes\; (K)}$')

########################## PLOTS OF VARIATION WRT K ########################## K=[0,1,2,3,4,5,6]










































# sns.distplot(sum_16_Simulation, kde=True, ax=ax[0], label='CP',color='r')
# ax[0].axvline(x=np.mean(sum_16_Simulation), ls='--', color='r', linewidth=2)
# sns.distplot(Robust_Total_Average_cost_to_save_Simulations[2], kde=True, ax=ax[0], label='Minimax', color='g')
# ax[0].axvline(x=np.mean(Robust_Total_Average_cost_to_save_Simulations[2]), ls='--', color='g', linewidth=2)
# sns.distplot(Regret_Total_Average_cost_to_save_Simulations[2], kde=True, ax=ax[0],label='Regret', color='orange')
# ax[0].axvline(x=np.mean(Regret_Total_Average_cost_to_save_Simulations[2]), ls='--', color='orange', linewidth=2)
# sns.distplot(sum_Contin_Simulation, kde=True, ax=ax[0], label='Nonrestricted', color='b')
# ax[0].axvline(x=np.mean(sum_Contin_Simulation), ls='--', color='b', linewidth=2)
# ax[0].legend(ncol=1, fontsize=20)
# ax[0].set_xlabel(r'$\rm{Total \; Cost \;(in \; \$1M)}$')
# ax[0].set_ylabel(r'\rm{Distribution}')
#
# sns.distplot(max_16_Simulation, kde=True, ax=ax[1], label='CP',color='r')
# ax[1].axvline(x=np.mean(max_16_Simulation), ls='--', color='r', linewidth=2)
# sns.distplot(Robust_Max_cost_to_save_Simulations[2], kde=True, ax=ax[1], label='Minimax', color='g')
# ax[1].axvline(x=np.mean(Robust_Max_cost_to_save_Simulations[2]), ls='--', color='g', linewidth=2)
# sns.distplot(Regret_Max_cost_to_save_Simulations[2], kde=True, ax=ax[1],label='Regret', color='orange')
# ax[1].axvline(x=np.mean(Regret_Max_cost_to_save_Simulations[2]), ls='--', color='orange', linewidth=2)
# sns.distplot(max_Contin_Simulation, kde=True, ax=ax[1], label='Nonrestricted', color='b')
# ax[1].axvline(x=np.mean(max_Contin_Simulation), ls='--', color='b', linewidth=2)
# ax[1].legend(ncol=1, fontsize=20)
# ax[1].set_xlabel(r'$\rm{Total \; Cost \;(in \; \$1M)}$')
# ax[1].set_ylabel(r'\rm{Distribution}')









#
# ax[0,0].hist((sum_Contin_Simulation), bins='auto', label='Nonrestricted')
# ax[0,0].axvline(x=np.mean(sum_Contin_Simulation), ls='--', color='b', linewidth=2)
# ax[0,0].hist((Regret_Total_Average_cost_to_save_Simulations[0]), bins='auto', label='Regret')
# ax[0,0].axvline(x=np.mean(Regret_Total_Average_cost_to_save_Simulations[0]), ls='--', color='orange', linewidth=2)
# ax[0,0].hist((Robust_Total_Average_cost_to_save_Simulations[0]), bins='auto', label='Minmax')
# ax[0,0].axvline(x=np.mean(Robust_Total_Average_cost_to_save_Simulations[0]), ls='--', color='g', linewidth=2)
# ax[0,0].hist((sum_16_Simulation), bins='auto', label='CP')
# ax[0,0].axvline(x=np.mean(sum_16_Simulation), ls='--', color='r', linewidth=2)
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
# ax[0,1].hist((sum_16_Simulation), bins='auto', label='CP')
# ax[0,1].axvline(x=np.mean(sum_16_Simulation), ls='--', color='r', linewidth=2)
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
# ax[1,0].hist((sum_16_Simulation), bins='auto', label='CP')
# ax[1,0].axvline(x=np.mean(sum_16_Simulation), ls='--', color='r', linewidth=2)
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
# ax[1,1].hist((sum_16_Simulation), bins='auto', label='CP')
# ax[1,1].axvline(x=np.mean(sum_16_Simulation), ls='--', color='r', linewidth=2)
# ax[1,1].legend(ncol=1, fontsize=20)
# ax[1,1].set_xlabel(r'$\rm{Cost \;(in \; \$1M)}$')
# ax[1,1].set_ylabel(r'\rm{Distribution}')










# ##########PLOTTING GOES HERE########## #Histograms for simulated values of the sum#
# matplotlib.rc('text', usetex=True)
# matplotlib.rcParams.update({'font.size': 30})
# matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
# fig, ax = plt.subplots(2,2)
# ax[0,0].scatter(T, F_of_n_star_continuous_Simulations_OverTime_Mean, marker='*', s=10, color='k')
# ax[0,0].plot(T, F_of_n_star_continuous_Simulations_OverTime_Max, color='k')
# ax[0,0].plot(T, F_of_n_star_continuous_Simulations_OverTime_Min, color='k')
# ax[0,0].fill_between(T,F_of_n_star_continuous_Simulations_OverTime_Max,F_of_n_star_continuous_Simulations_OverTime_Min, alpha=0.1, color='k')
#
# ax[0,0].scatter(T, F_of_n_equal_16_Simulations_Simulations_OverTime_Mean, marker='*', s=10, color='g')
# ax[0,0].plot(T, F_of_n_equal_16_Simulations_Simulations_OverTime_Max, color='g', ls='-', dashes=(3,3))
# ax[0,0].plot(T, F_of_n_equal_16_Simulations_Simulations_OverTime_Min, color='g', ls='-', dashes=(3,3))
# ax[0,0].fill_between(T,F_of_n_equal_16_Simulations_Simulations_OverTime_Max,F_of_n_equal_16_Simulations_Simulations_OverTime_Min, alpha=0.1, color='g')
#
# ax[0,0].scatter(T, Regret_F_constrained_matrix_of_ks_Simulations_OverTime_Mean[0], marker='*', s=10, color='b')
# ax[0,0].plot(T, Regret_F_constrained_matrix_of_ks_Simulations_OverTime_Max[0], color='b', ls='-.')
# ax[0,0].plot(T, Regret_F_constrained_matrix_of_ks_Simulations_OverTime_Min[0], color='b', ls='-.')
# ax[0,0].fill_between(T,Regret_F_constrained_matrix_of_ks_Simulations_OverTime_Max[0],Regret_F_constrained_matrix_of_ks_Simulations_OverTime_Min[0], alpha=0.1, color='b')
#
# ax[0,0].scatter(T, Robust_F_constrained_matrix_of_ks_Simulations_OverTime_Mean[0], marker='*', s=10, color='r')
# ax[0,0].plot(T, Robust_F_constrained_matrix_of_ks_Simulations_OverTime_Max[0], color='r', ls=':')
# ax[0,0].plot(T, Robust_F_constrained_matrix_of_ks_Simulations_OverTime_Min[0], color='r', ls=':')
# ax[0,0].fill_between(T,Robust_F_constrained_matrix_of_ks_Simulations_OverTime_Max[0],Robust_F_constrained_matrix_of_ks_Simulations_OverTime_Min[0], alpha=0.1, color='r')
#
#
#
#
# ax[0,1].scatter(T, F_of_n_star_continuous_Simulations_OverTime_Mean, marker='*', s=10, color='k')
# ax[0,1].plot(T, F_of_n_star_continuous_Simulations_OverTime_Max, color='k')
# ax[0,1].plot(T, F_of_n_star_continuous_Simulations_OverTime_Min, color='k')
# ax[0,1].fill_between(T,F_of_n_star_continuous_Simulations_OverTime_Max,F_of_n_star_continuous_Simulations_OverTime_Min, alpha=0.1, color='k')
#
# ax[0,1].scatter(T, F_of_n_equal_16_Simulations_Simulations_OverTime_Mean, marker='*', s=10, color='g')
# ax[0,1].plot(T, F_of_n_equal_16_Simulations_Simulations_OverTime_Max, color='g', ls='-', dashes=(3,3))
# ax[0,1].plot(T, F_of_n_equal_16_Simulations_Simulations_OverTime_Min, color='g', ls='-', dashes=(3,3))
# ax[0,1].fill_between(T,F_of_n_equal_16_Simulations_Simulations_OverTime_Max,F_of_n_equal_16_Simulations_Simulations_OverTime_Min, alpha=0.1, color='g')
#
# ax[0,1].scatter(T, Regret_F_constrained_matrix_of_ks_Simulations_OverTime_Mean[1], marker='*', s=10, color='b')
# ax[0,1].plot(T, Regret_F_constrained_matrix_of_ks_Simulations_OverTime_Max[1], color='b', ls='-.')
# ax[0,1].plot(T, Regret_F_constrained_matrix_of_ks_Simulations_OverTime_Min[1], color='b', ls='-.')
# ax[0,1].fill_between(T,Regret_F_constrained_matrix_of_ks_Simulations_OverTime_Max[1],Regret_F_constrained_matrix_of_ks_Simulations_OverTime_Min[1], alpha=0.1, color='b')
#
# ax[0,1].scatter(T, Robust_F_constrained_matrix_of_ks_Simulations_OverTime_Mean[1], marker='*', s=10, color='r')
# ax[0,1].plot(T, Robust_F_constrained_matrix_of_ks_Simulations_OverTime_Max[1], color='r', ls=':')
# ax[0,1].plot(T, Robust_F_constrained_matrix_of_ks_Simulations_OverTime_Min[1], color='r', ls=':')
# ax[0,1].fill_between(T,Robust_F_constrained_matrix_of_ks_Simulations_OverTime_Max[1],Robust_F_constrained_matrix_of_ks_Simulations_OverTime_Min[1], alpha=0.1, color='r')
#
#
#
#
# ax[1,0].scatter(T, F_of_n_star_continuous_Simulations_OverTime_Mean, marker='*', s=10, color='k')
# ax[1,0].plot(T, F_of_n_star_continuous_Simulations_OverTime_Max, color='k')
# ax[1,0].plot(T, F_of_n_star_continuous_Simulations_OverTime_Min, color='k')
# ax[1,0].fill_between(T,F_of_n_star_continuous_Simulations_OverTime_Max,F_of_n_star_continuous_Simulations_OverTime_Min, alpha=0.1, color='k')
#
# ax[1,0].scatter(T, F_of_n_equal_16_Simulations_Simulations_OverTime_Mean, marker='*', s=10, color='g')
# ax[1,0].plot(T, F_of_n_equal_16_Simulations_Simulations_OverTime_Max, color='g', ls='-', dashes=(3,3))
# ax[1,0].plot(T, F_of_n_equal_16_Simulations_Simulations_OverTime_Min, color='g', ls='-', dashes=(3,3))
# ax[1,0].fill_between(T,F_of_n_equal_16_Simulations_Simulations_OverTime_Max,F_of_n_equal_16_Simulations_Simulations_OverTime_Min, alpha=0.1, color='g')
#
# ax[1,0].scatter(T, Regret_F_constrained_matrix_of_ks_Simulations_OverTime_Mean[2], marker='*', s=10, color='b')
# ax[1,0].plot(T, Regret_F_constrained_matrix_of_ks_Simulations_OverTime_Max[2], color='b', ls='-.')
# ax[1,0].plot(T, Regret_F_constrained_matrix_of_ks_Simulations_OverTime_Min[2], color='b', ls='-.')
# ax[1,0].fill_between(T,Regret_F_constrained_matrix_of_ks_Simulations_OverTime_Max[2],Regret_F_constrained_matrix_of_ks_Simulations_OverTime_Min[2], alpha=0.1, color='b')
#
# ax[1,0].scatter(T, Robust_F_constrained_matrix_of_ks_Simulations_OverTime_Mean[2], marker='*', s=10, color='r')
# ax[1,0].plot(T, Robust_F_constrained_matrix_of_ks_Simulations_OverTime_Max[2], color='r', ls=':')
# ax[1,0].plot(T, Robust_F_constrained_matrix_of_ks_Simulations_OverTime_Min[2], color='r', ls=':')
# ax[1,0].fill_between(T,Robust_F_constrained_matrix_of_ks_Simulations_OverTime_Max[2],Robust_F_constrained_matrix_of_ks_Simulations_OverTime_Min[2], alpha=0.1, color='r')
#
#
#
# ax[1,1].scatter(T, F_of_n_star_continuous_Simulations_OverTime_Mean, marker='*', s=10, color='k')
# ax[1,1].plot(T, F_of_n_star_continuous_Simulations_OverTime_Max, color='k')
# ax[1,1].plot(T, F_of_n_star_continuous_Simulations_OverTime_Min, color='k')
# ax[1,1].fill_between(T,F_of_n_star_continuous_Simulations_OverTime_Max,F_of_n_star_continuous_Simulations_OverTime_Min, alpha=0.1, color='k')
#
# ax[1,1].scatter(T, F_of_n_equal_16_Simulations_Simulations_OverTime_Mean, marker='*', s=10, color='g')
# ax[1,1].plot(T, F_of_n_equal_16_Simulations_Simulations_OverTime_Max, color='g', ls='-', dashes=(3,3))
# ax[1,1].plot(T, F_of_n_equal_16_Simulations_Simulations_OverTime_Min, color='g', ls='-', dashes=(3,3))
# ax[1,1].fill_between(T,F_of_n_equal_16_Simulations_Simulations_OverTime_Max,F_of_n_equal_16_Simulations_Simulations_OverTime_Min, alpha=0.1, color='g')
#
# ax[1,1].scatter(T, Regret_F_constrained_matrix_of_ks_Simulations_OverTime_Mean[3], marker='*', s=10, color='b')
# ax[1,1].plot(T, Regret_F_constrained_matrix_of_ks_Simulations_OverTime_Max[3], color='b', ls='-.')
# ax[1,1].plot(T, Regret_F_constrained_matrix_of_ks_Simulations_OverTime_Min[3], color='b', ls='-.')
# ax[1,1].fill_between(T,Regret_F_constrained_matrix_of_ks_Simulations_OverTime_Max[3],Regret_F_constrained_matrix_of_ks_Simulations_OverTime_Min[3], alpha=0.1, color='b')
#
# ax[1,1].scatter(T, Robust_F_constrained_matrix_of_ks_Simulations_OverTime_Mean[3], marker='*', s=10, color='r')
# ax[1,1].plot(T, Robust_F_constrained_matrix_of_ks_Simulations_OverTime_Max[3], color='r', ls=':')
# ax[1,1].plot(T, Robust_F_constrained_matrix_of_ks_Simulations_OverTime_Min[3], color='r', ls=':')
# ax[1,1].fill_between(T,Robust_F_constrained_matrix_of_ks_Simulations_OverTime_Max[3],Robust_F_constrained_matrix_of_ks_Simulations_OverTime_Min[3], alpha=0.1, color='r')