import numpy as np
import matplotlib.pylab as plt
import matplotlib
matplotlib.rc('text', usetex=True)
matplotlib.rcParams.update({'font.size': 36})
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]


avgcostofTTIs=16000000000*(13.6/53)/26000000


Path_avg_B=np.array([ 0, 21, 30, 35, 52])
Path_constant_B=np.array([ 0, 21, 30, 35, 52])
avg_bd_solution_path=np.round([44.09595089869998, 47.812843888336644, 46.50101313801885, 44.23135391139709],1)
constant_bd_solution_path=[45]*len(avg_bd_solution_path)

X_path_avg_Bnotrounded=np.round([np.array([1.01899295e+01, 9.45562291e+00, 2.10123323e+01, 3.57242965e+00,
         2.44842643e-08]),
  np.array([10.13575997,  9.45554102, 20.9322416 ,  5.94199016,  1.29841588]),
  np.array([10.09891364,  9.43006892, 20.87122912,  4.79183091,  1.16837416]),
  np.array([10.3541545 ,  9.86565357, 21.35672104,  2.53862256,  0.11620221])],1)


X_path_avg_B=[]
for j in range(len(X_path_avg_Bnotrounded)):
    D=[]
    for m in range(len(X_path_avg_Bnotrounded[j])):
        if m!=2:
            D.append(X_path_avg_Bnotrounded[j][m])
        else:
            D.append(round(avg_bd_solution_path[j]-sum(X_path_avg_Bnotrounded[j])+X_path_avg_Bnotrounded[j][m],2))
    X_path_avg_B.append(D)


X_path_constant_Bnotrounded=np.round([np.array([1.03304465e+01, 9.70140078e+00, 2.12933313e+01, 3.67482139e+00,
         2.88095405e-08]),
  np.array([ 9.67977191,  8.69351649, 20.03069925,  5.60189954,  0.99411279]),
  np.array([ 9.87126171,  9.05065173, 20.42157142,  4.63317612,  1.023339  ]),
  np.array([10.48881777, 10.09509927, 21.62398364,  2.62845362,  0.16364569])],1)

BT=45
X_path_constant_B=[]
for j in range(len(X_path_constant_Bnotrounded)):
    D=[]
    for m in range(len(X_path_constant_Bnotrounded[j])):
        if m!=2:
            D.append(X_path_constant_Bnotrounded[j][m])
        else:
            D.append(round(BT-sum(X_path_constant_Bnotrounded[j])+X_path_constant_Bnotrounded[j][m],2))
    X_path_constant_B.append(D)




f_avg_B=avgcostofTTIs*np.array([0.00217685, 0.00216009, 0.00214332, 0.00212656, 0.00213643,
        0.0021463 , 0.00215618, 0.00216605, 0.00217326, 0.00218046,
        0.00218767, 0.00219488, 0.00221101, 0.00222713, 0.00224326,
        0.00225939, 0.00227551, 0.00233866, 0.0024018 , 0.00246495,
        0.00252809, 0.00232147, 0.00240145, 0.00248144, 0.00256142,
        0.00259571, 0.00263   , 0.00266429, 0.00269858, 0.00273287,
        0.00285339, 0.0027744 , 0.0026954 , 0.0026164 , 0.00254484,
        0.00269395, 0.00259171, 0.00248947, 0.00245419, 0.00241891,
        0.00238363, 0.00234834, 0.00231306, 0.0022842 , 0.00225534,
        0.00222648, 0.00219763, 0.00216877, 0.00213991, 0.00211105,
        0.00208384, 0.00205664, 0.00202944])
f_constant_B=avgcostofTTIs*np.array([0.00209297, 0.00207684, 0.00206071, 0.00204459, 0.00205409,
        0.00206359, 0.00207309, 0.00208259, 0.00208952, 0.00209646,
        0.0021034 , 0.00211033, 0.00212586, 0.00214139, 0.00215692,
        0.00217245, 0.00218798, 0.00224885, 0.00230971, 0.00237058,
        0.00243144, 0.00262827, 0.00271685, 0.00280543, 0.00289401,
        0.00292896, 0.00296391, 0.00299886, 0.00303382, 0.00306877,
        0.00302176, 0.00293855, 0.00285533, 0.00277212, 0.00269814,
        0.00260403, 0.00250347, 0.00240291, 0.00236835, 0.00233379,
        0.00229924, 0.00226468, 0.00223012, 0.00220219, 0.00217426,
        0.00214632, 0.00211839, 0.00209046, 0.00206253, 0.0020346 ,
        0.00200831, 0.00198202, 0.00195573])
Ideal_function_avg_B=avgcostofTTIs*np.array([0.002204937953920871, 0.0022048742976612163,0.0022048106321745906, 0.002204746947667144,0.0022048543611030164, 0.002204961766584327, 0.002205069165708776, 0.0022051765622838823,0.00220527202097988, 0.002205367478577314, 0.0022054629362569665, 0.0022055583921534666, 0.0022060595119790852, 0.0022065606322792183, 0.0022070617506179868, 0.0022075628674535623, 0.0022080639831408904, 0.0022114464777542187, 0.0022148289660832783, 0.0022182114712676355, 0.002221593974716134, 0.002260322472553183, 0.0022990510673416035, 0.0023377799597551776, 0.0023765079498632717, 0.0024411269250007336, 0.0025057458901667972, 0.0025703658374599866, 0.0026349847476699593, 0.0026996035583934123, 0.00267507864833629,0.002650554709982507, 0.00262602975349816, 0.002601505783383971, 0.002529333853155271, 0.002457162893741406, 0.0023849909195716736, 0.0023128199367136453, 0.0022938633390369016,0.0022749071408198515, 0.0022559508420612493, 0.002236994642751869, 0.0022180384428756796, 0.002216344236635478, 0.002214649928853772, 0.0022129557289603755, 0.0022112614961902035,0.0022095672495344707, 0.0022078730073284815, 0.002206178746122881, 0.0022061737864175835, 0.0022061736828560548, 0.002206173172352297])

Ideal_function_of_constant_budget=avgcostofTTIs*np.array([0.0020532862129138403, 0.002024658126165226, 0.0019931687223876623, 0.001958117692226675, 0.001979164987190433, 0.0019988749331910248, 0.0020174181244323485, 0.0020349350667255737, 0.0020471307598892807, 0.002058881682523727, 0.0020702224197093, 0.0020811821132106427, 0.0021043430980254773, 0.0021260314913253506, 0.0021464389030190407, 0.0021657219889518963, 0.002184009004365657, 0.0022479517911555215, 0.0023028853490681764, 0.002351279387452701, 0.002394695467918006, 0.0025036936846319687, 0.0026016854741469325, 0.002691824039898285, 0.002776056752589893, 0.0028260532618642727, 0.002875497380823374, 0.0029243428673799503, 0.002972534421256997, 0.003020011631105938, 0.0029526586080121, 0.002879936285685088, 0.0028001406558785936, 0.0027105704574900674, 0.0026218772475510525, 0.0025323813705127496, 0.0024419919937632975, 0.0023506054656570884, 0.002321030812902273, 0.002291120811253816, 0.0022608509369190016, 0.0022301946443910527, 0.002199122776081044, 0.0021767541520974717, 0.0021530281892814087, 0.002127737560699177, 0.0021006210100739473, 0.0020713442762667323, 0.0020394673155499934, 0.0020043969932429283, 0.0019669921747422496, 0.0019243494269558406, 0.0018745681252022458])
T=np.arange(0,53)

for i in range(len(Path_avg_B)-1):
     # I have done one line as both methods resulted in same points in time
     plt.axvline(Path_avg_B[i],c='r',ls='-',dashes=(5,2))



# plt.text(Path_avg_B[0]+0.2,3.08,'$b=\${}$'.format(avg_bd_solution_path[0]),c='b')
# plt.text(Path_avg_B[1]+0.2,2.22,'$b=\${}$'.format(avg_bd_solution_path[1]),c='b')
# plt.text(Path_avg_B[2]+0.2,3.08,'$b=\${}$'.format(avg_bd_solution_path[2]),c='b')
# plt.text(Path_avg_B[3]+0.2,3.08,'$b=\${}$'.format(avg_bd_solution_path[3]),c='b')
#
# plt.text(Path_constant_B[0]+0.2,2.62,r'$b=\${}$'.format(constant_bd_solution_path[0]),c='k')
# plt.text(Path_constant_B[1]+5,2.4,r'$b=\${}$'.format(constant_bd_solution_path[1]),c='k')
# plt.text(Path_constant_B[2]+0.2,2.22,r'$b=\${}$'.format(constant_bd_solution_path[2]),c='k')
# plt.text(Path_constant_B[3]+0.2,2.22,r'$b=\${}$'.format(constant_bd_solution_path[3]),c='k')
#
#
# plt.text(Path_avg_B[0]+0.2, 2.69,r'$\boldsymbol{x}^{*}=$' +'\n'+'['+ '$\$$'+ '{}'.format(X_path_avg_B[0][0]) + '\n' + '$\$$'+ '{}'.format(X_path_avg_B[0][1]) + '\n'+'$\$$'+ '{}'.format(X_path_avg_B[0][2]) + '\n'+'$\$$'+ '{}'.format(X_path_avg_B[0][3])+'\n'+'$\$$'+ '{}'.format(X_path_avg_B[0][4]) +']',c='b',fontsize=28)
# plt.text(Path_avg_B[1]+0.2, 1.83,r'$\boldsymbol{x}^{*}=$' +'\n'+'['+ '$\$$'+ '{}'.format(X_path_avg_B[1][0]) + '\n' + '$\$$'+ '{}'.format(X_path_avg_B[1][1]) + '\n'+'$\$$'+ '{}'.format(X_path_avg_B[1][2]) + '\n'+'$\$$'+ '{}'.format(X_path_avg_B[1][3])+ '\n'+'$\$$'+ '{}'.format(X_path_avg_B[1][4])+']',c='b',fontsize=28)
# plt.text(Path_avg_B[2]+0.2, 2.28, r'$\boldsymbol{x}^{*}=$' +'\n'+'['+ '$\$$'+ '{}'.format(X_path_avg_B[2][0]) + '\n' + '$\$$'+ '{}'.format(X_path_avg_B[2][1]) + '\n'+'$\$$'+ '{}'.format(X_path_avg_B[2][2]) + '\n'+'$\$$'+ '{}'.format(X_path_avg_B[2][3])+ '\n'+'$\$$'+ '{}'.format(X_path_avg_B[2][4])+']',c='b',fontsize=28)
# plt.text(Path_avg_B[3]+0.2, 2.69,r'$\boldsymbol{x}^{*}=$' +'\n'+'['+ '$\$$'+ '{}'.format(X_path_avg_B[3][0]) + '\n' + '$\$$'+ '{}'.format(X_path_avg_B[3][1]) + '\n'+'$\$$'+ '{}'.format(X_path_avg_B[3][2]) + '\n'+'$\$$'+ '{}'.format(X_path_avg_B[3][3])+ '\n'+'$\$$'+ '{}'.format(X_path_avg_B[3][4])+']',c='b',fontsize=28)
#
# plt.text(Path_constant_B[0]+0.2, 2.23,r'$\boldsymbol{x}^{*}=$' +'\n'+'['+ '$\$$'+ '{}'.format(X_path_constant_B[0][0]) + '\n' + '$\$$'+ '{}'.format(X_path_constant_B[0][1]) + '\n'+'$\$$'+ '{}'.format(X_path_constant_B[0][2]) + '\n'+'$\$$'+ '{}'.format(X_path_constant_B[0][3])+ '\n'+'$\$$'+ '{}'.format(X_path_constant_B[0][4])+']',fontsize=28)
# plt.text(Path_constant_B[1]+5.75, 1.99, r'$\boldsymbol{x}^{*}=$' +'\n'+'['+ '$\$$'+ '{}'.format(X_path_constant_B[1][0]) + '\n' + '$\$$'+ '{}'.format(X_path_constant_B[1][1]) + '\n'+'$\$$'+ '{}'.format(X_path_constant_B[1][2]) + '\n'+'$\$$'+ '{}'.format(X_path_constant_B[1][3]) + '\n'+'$\$$'+ '{}'.format(X_path_constant_B[1][4])+']',fontsize=28)
# plt.text(Path_constant_B[2]+0.2, 1.83,r'$\boldsymbol{x}^{*}=$' +'\n'+'['+ '$\$$'+ '{}'.format(X_path_constant_B[2][0]) + '\n' + '$\$$'+ '{}'.format(X_path_constant_B[2][1]) + '\n'+'$\$$'+ '{}'.format(X_path_constant_B[2][2]) + '\n'+'$\$$'+ '{}'.format(X_path_constant_B[2][3])+ '\n'+'$\$$'+ '{}'.format(X_path_constant_B[2][4])+']',fontsize=28)
# plt.text(Path_constant_B[3]+0.2, 1.83,r'$\boldsymbol{x}^{*}=$' +'\n'+'['+ '$\$$'+ '{}'.format(X_path_constant_B[3][0]) + '\n' + '$\$$'+ '{}'.format(X_path_constant_B[3][1]) + '\n'+'$\$$'+ '{}'.format(X_path_constant_B[3][2]) + '\n'+'$\$$'+ '{}'.format(X_path_constant_B[3][3])+ '\n'+'$\$$'+ '{}'.format(X_path_constant_B[3][4])+']',fontsize=28)


table0 = r'\begin{tabular}{ l | r | r } & flx & unf \\ $\$b$ & '+'{:.1f}'.format(avg_bd_solution_path[0]) + '& '+'{:.1f}'.format(constant_bd_solution_path[0]) + r'\\ $x^{*}_{1}$'+'& {:.1f}'.format(X_path_avg_B[0][0])+ '& {:.1f}'.format(X_path_constant_B[0][0]) + r'\\ $x^{*}_{2}$'+'& {:.1f}'.format(X_path_avg_B[0][1])+ '& {:.1f}'.format(X_path_constant_B[0][1]) + r'\\ $x^{*}_{3}$'+'& {:.1f}'.format(X_path_avg_B[0][2])+ '& {:.1f}'.format(X_path_constant_B[0][2]) + r'\\ $x^{*}_{4}$'+'& {:.1f}'.format(X_path_avg_B[0][3])+ '& {:.1f}'.format(X_path_constant_B[0][3]) + r'\\ $x^{*}_{5}$'+'& {:.1f}'.format(X_path_avg_B[0][4])+ '& {:.1f}'.format(X_path_constant_B[0][4]) +' \end{tabular}'
plt.text(Path_avg_B[0]+0.2, 0.0027*avgcostofTTIs,table0,size=21)

table1 = r'\begin{tabular}{ l | r | r } & flx & unf \\ $\$ b$ & '+'{:.1f}'.format(avg_bd_solution_path[1]) + '& '+'{:.1f}'.format(constant_bd_solution_path[1]) + r'\\ $x^{*}_{1}$'+'& {:.1f}'.format(X_path_avg_B[1][0])+ '& {:.1f}'.format(X_path_constant_B[1][0]) + r'\\ $x^{*}_{2}$'+'& {:.1f}'.format(X_path_avg_B[1][1])+ '& {:.1f}'.format(X_path_constant_B[1][1]) + r'\\ $x^{*}_{3}$'+'& {:.1f}'.format(X_path_avg_B[1][2])+ '& {:.1f}'.format(X_path_constant_B[1][2]) + r'\\ $x^{*}_{4}$'+'& {:.1f}'.format(X_path_avg_B[1][3])+ '& {:.1f}'.format(X_path_constant_B[1][3]) + r'\\ $x^{*}_{5}$'+'& {:.1f}'.format(X_path_avg_B[1][4])+ '& {:.1f}'.format(X_path_constant_B[1][4]) +' \end{tabular}'
plt.text(Path_avg_B[1]+0.2,0.00212*avgcostofTTIs,table1,size=21)

table2 = r'\begin{tabular}{ l | r | r } & flx & unf \\ $\$b$ & '+'{:.1f}'.format(avg_bd_solution_path[2]) + '& '+'{:.1f}'.format(constant_bd_solution_path[2]) + r'\\ $x^{*}_{1}$'+'& {:.1f}'.format(X_path_avg_B[2][0])+ '& {:.1f}'.format(X_path_constant_B[2][0]) + r'\\ $x^{*}_{2}$'+'& {:.1f}'.format(X_path_avg_B[2][1])+ '& {:.1f}'.format(X_path_constant_B[2][1]) + r'\\ $x^{*}_{3}$'+'& {:.1f}'.format(X_path_avg_B[2][2])+ '& {:.1f}'.format(X_path_constant_B[2][2]) + r'\\ $x^{*}_{4}$'+'& {:.1f}'.format(X_path_avg_B[2][3])+ '& {:.1f}'.format(X_path_constant_B[2][3]) + r'\\ $x^{*}_{5}$'+'& {:.1f}'.format(X_path_avg_B[2][4])+ '& {:.1f}'.format(X_path_constant_B[2][4]) +' \end{tabular}'
plt.text(Path_avg_B[2]+0.05,0.00212*avgcostofTTIs,table2,size=21)

# table21 = r'\begin{tabular}{ l | l } & flx  \\ $\$b$ & '+'{:.1f}'.format(avg_bd_solution_path[2]) + r'\\ $x^{*}_{1}$'+'& {:.1f}'.format(X_path_avg_B[2][0]) + r'\\ $x^{*}_{2}$'+'& {:.1f}'.format(X_path_avg_B[2][1]) + r'\\ $x^{*}_{3}$'+'& {:.1f}'.format(X_path_avg_B[2][2]) + r'\\ $x^{*}_{4}$'+'& {:.1f}'.format(X_path_avg_B[2][3])+  r'\\ $x^{*}_{5}$'+'& {:.1f}'.format(X_path_avg_B[2][4]) +' \end{tabular}'
# plt.text(Path_avg_B[2]+0.01,2.45,table21,size=22)
# table22 = r'\begin{tabular}{ l | l} & unf  \\ $\$b$' + '& '+'{:.1f}'.format(constant_bd_solution_path[2]) + r'\\ $x^{*}_{1}$'+ '& {:.1f}'.format(X_path_constant_B[2][0]) + r'\\ $x^{*}_{2}$'+ '& {:.1f}'.format(X_path_constant_B[2][1]) + r'\\ $x^{*}_{3}$'+ '& {:.1f}'.format(X_path_constant_B[2][2]) + r'\\ $x^{*}_{4}$'+ '& {:.1f}'.format(X_path_constant_B[2][3]) + r'\\ $x^{*}_{5}$'+ '& {:.1f}'.format(X_path_constant_B[2][4]) +' \end{tabular}'
# plt.text(Path_avg_B[2]+1.5,2.08,table22,size=22)

table3 = r'\begin{tabular}{ l | r | r } & flx & unf \\ $\$b$ & '+'{:.1f}'.format(avg_bd_solution_path[3]) + '& '+'{:.1f}'.format(constant_bd_solution_path[3]) + r'\\ $x^{*}_{1}$'+'& {:.1f}'.format(X_path_avg_B[3][0])+ '& {:.1f}'.format(X_path_constant_B[3][0]) + r'\\ $x^{*}_{2}$'+'& {:.1f}'.format(X_path_avg_B[3][1])+ '& {:.1f}'.format(X_path_constant_B[3][1]) + r'\\ $x^{*}_{3}$'+'& {:.1f}'.format(X_path_avg_B[3][2])+ '& {:.1f}'.format(X_path_constant_B[3][2]) + r'\\ $x^{*}_{4}$'+'& {:.1f}'.format(X_path_avg_B[3][3])+ '& {:.1f}'.format(X_path_constant_B[3][3]) + r'\\ $x^{*}_{5}$'+'& {:.1f}'.format(X_path_avg_B[3][4])+ '& {:.1f}'.format(X_path_constant_B[3][4]) +' \end{tabular}'
plt.text(Path_avg_B[3]+0.5,0.00285*avgcostofTTIs,table3,size=21)


# plt.plot(T,1000*Ideal_function_of_constant_budget,c='k',linewidth=0.8,label='$f^{*}_{t}$')
# plt.plot(T,1000*Ideal_function_avg_B,c='b',linewidth=0.8,label='$g^{*}_{t}$')
plt.plot(T,f_constant_B,c='k',linewidth=0.8,ls='-.',label='uniform')
plt.plot(T,f_avg_B,linewidth=0.8,ls='-.',c='b',label='flexible')
plt.text(10,3,'$K=3$')
# plt.title('$K=3$')
plt.legend(loc='best')
plt.ylabel(r'$\rm{Cost}$ ' + r'$(\rm{in\;} 1\rm{M}\$)$')
plt.xticks(np.arange(53), ('Jan ', ' ', ' ', ' ', 'Feb', ' ', ' ', ' ', 'Mar', ' ', ' ', ' ', 'Apr', ' ', ' ', ' ', ' ','May', ' ', ' ', ' ',' Jun', ' ', ' ', ' ', 'Jul', ' ', ' ', ' ',' ' ,'Aug', ' ', ' ', ' ', 'Sep', ' ', ' ', ' ', 'Oct', ' ', ' ', ' ', ' ' ,'Nov', ' ', ' ', ' ', 'Dec', ' ', ' ', ' ', ' ', ' '))
plt.xlim(-1,53)