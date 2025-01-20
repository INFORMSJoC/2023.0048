import pandas as pd
import numpy as np
import matplotlib.pylab as plt
import matplotlib
import math

matplotlib.rc('text', usetex=True)
matplotlib.rcParams.update({'font.size': 34})
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]


def vl_max_formula_OPTIMIZED_as_a_lowerbd_on_theworstcase(N):
    K=[math.floor((12.5+2*N)/9),math.ceil((12.5+2*N)/9)]
    d=[K[0]-1,K[1]-1]
    a=[N-7+5.5*d[0]+2*d[0]*N-2*d[0]*K[0]-2.5*d[0]**2 ,N-7+5.5*d[1]+2*d[1]*N-2*d[1]*K[1]-2.5*d[1]**2]
    return int(max(a))



N_array=[25,50,100]
K_allowable_changes_array=np.arange(2,11,1)
gap_array=np.linspace(0,1,100)

percentage_vl_upper_bound_array=np.array([(vl_max_formula_OPTIMIZED_as_a_lowerbd_on_theworstcase(j))/((j*(j-1))/2) for j in N_array])



# Read data random and minmax rij
Data_minmax_rij=[]
Data_random_rij=[]
for j in range(len(N_array)):
    df1=pd.read_excel(r'C:\Users\Marwan\OneDrive - American University of Beirut\Desktop\Local_Disk_D Marwan\GEORGE MASON UNIVERSITY\Research 1\python trials\untitled1\untitled\data for new runs complexity\caseStudyNumericalExperimentsminmaxregret.xlsx','N={}'.format(N_array[j]),index_col=0)
    df2=pd.read_csv(r'C:\Users\Marwan\OneDrive - American University of Beirut\Desktop\Local_Disk_D Marwan\GEORGE MASON UNIVERSITY\Research 1\python trials\untitled1\untitled\data for new runs complexity\caseStudyNumericalExperimentsrandomized when N={}'.format(N_array[j])+'.csv', index_col=0)
    df1['average_over_K']=np.average(df1.iloc[:,:],axis=1)
    df2['average_over_K']=np.average(df2.iloc[:,:],axis=1)
    Data_minmax_rij.append(df1)
    Data_random_rij.append(df2)





# Plot data
colors=['g','r','b']
for i in range(len(N_array)):
    plt.plot(gap_array,Data_random_rij[i].loc[:, 'average_over_K'],ls='-',dashes=(5,3),label='randomized, $T={}$'.format(N_array[i]),linewidth=1.5,c=colors[i])
for i in range(len(N_array)):
    plt.plot(gap_array,Data_minmax_rij[i].loc[:, 'average_over_K'],label='with structure, $T={}$'.format(N_array[i]),linewidth=1.5,c=colors[i])
# for i in range(len(N_array)):
#     plt.axhline(percentage_vl_upper_bound_array[i],ls=':',c=colors[i],dashes=(1,3),label=r'$\frac{\underline{v}(T)}{Q}$, ' +'$T={}$'.format(N_array[i]),linewidth=2.5)

plt.ylabel(r'Query ratio')
plt.xlabel('Tolerance '+r'$\varepsilon$')
plt.xlim(0,1)
# plt.legend(bbox_to_anchor=(0.6, 0.25, 0.4, 0.5),fontsize=27,ncol=3)
plt.legend(fontsize=27,ncol=2)
plt.show()

















# NEEDE BY DR HADI TO CHECK IF AT ANY TIME WE EXCEEDED THE LOWER BOUND
# epsarray=[0,0.05,0.1,0.2,0.3]
# Narray=[25,50,100]
# checkarray=[]
# for i in epsarray:
#     m=0
#     for j in N_array:
#         df=pd.read_csv('data of N={}'.format(j) +' epsilon={}'.format(i)+'.csv',index_col=0)
#         df=df/((j*(j-1))/2)
#         df=df.to_numpy()
#         check=df[df>=percentage_vl_upper_bound_array[m]]
#         print(check)
#         if check.size!=0:
#             checkarray.append(1)
#         m=m+1