import matplotlib.pyplot as plt
import numpy as np
import matplotlib

matplotlib.rc('text', usetex=True)
matplotlib.rcParams.update({'font.size': 35})


def create_n_array(time_changes,group_size_changes,t):
    n_array=[]
    j=0
    for i in range(len(t)):
        if i<time_changes[j]:
            n_array.append(group_size_changes[j])
        elif j==len(time_changes)-1:
            while i<=len(t)-1:
                n_array.append(group_size_changes[-1])
                i=i+1
            return n_array
        else:
            j=j+1
            n_array.append(group_size_changes[j])
    return n_array


T=np.arange(0,20,1)
K_allowable_changes_array =[1, 2, 3, 4]
nstarsolution_continuous=np.round([102.68766640163948,79.16507602633837, 62.02062525408915, 49.45182390893853, 40.19096291768017, 32.538774314969885, 28.098844733945356, 24.882926702827582, 23.020521896535183, 22.455900157596272, 22.5325963967961, 23.57547521121782, 25.6944098304785, 29.06466642972339, 33.940097064745196, 41.58535597857481, 50.97648541597246, 63.54696481847239, 79.57687501448439, 95.93038892239838],0)

path_matrix_of_ks =[np.array([ 0, 19]), np.array([ 0,  4, 19]), np.array([ 0,  4, 15, 19]), np.array([ 0,  4,  14, 17, 19])]
group_size_path_of_ks =[np.array([28]), np.array([61, 38]), np.array([61, 33, 53]), np.array([61, 33, 44, 74])]

solK0=[ group_size_path_of_ks[0][0] ]*len(T)
solK1=[ group_size_path_of_ks[1][0] ]*(path_matrix_of_ks[1][1]) + [ group_size_path_of_ks[1][1] ]*(path_matrix_of_ks[1][-1] - path_matrix_of_ks[1][1] +1)
solK2=[ group_size_path_of_ks[2][0] ]*(path_matrix_of_ks[2][1]) + [ group_size_path_of_ks[2][1] ]*(path_matrix_of_ks[2][2] - path_matrix_of_ks[2][1]) + [ group_size_path_of_ks[2][2] ]*(path_matrix_of_ks[2][-1] - path_matrix_of_ks[2][2] +1)
solK3=[ group_size_path_of_ks[3][0] ]*(path_matrix_of_ks[3][1]) + [ group_size_path_of_ks[3][1] ]*(path_matrix_of_ks[3][2] - path_matrix_of_ks[3][1]) + [ group_size_path_of_ks[3][2] ]*(path_matrix_of_ks[3][3] - path_matrix_of_ks[3][2] ) + [ group_size_path_of_ks[3][-1] ]*(path_matrix_of_ks[3][-1] - path_matrix_of_ks[3][3] +1)

plt.figure(figsize=(12,8))
plt.scatter(T,nstarsolution_continuous,marker='d',label='unrestricted',s=75)
plt.scatter(T,solK0,marker='*',label='$K=0$',s=75)
plt.scatter(T,solK1,marker='P',label='$K=1$',s=75)
plt.scatter(T,solK2,marker='X',label='$K=2$',s=75)
plt.scatter(T,solK3,marker='8',label='$K=3$',s=75)

plt.plot(T,nstarsolution_continuous,linewidth=1.5)
plt.plot(T,solK0,linewidth=1.5)
plt.plot(T,solK1,linewidth=1.5)
plt.plot(T,solK2,linewidth=1.5)
plt.plot(T,solK3,linewidth=1.5)
# plt.xlim(0,len(T)-1)
plt.xticks(np.arange(20), ('Jun', ' ', ' ', ' ', 'Jul', ' ', ' ', ' ', 'Aug', ' ', ' ', ' ', 'Sep', ' ', ' ', ' ', 'Oct', ' ', ' ', ' '))
plt.ylabel('Optimal group size')
plt.legend(loc='best')