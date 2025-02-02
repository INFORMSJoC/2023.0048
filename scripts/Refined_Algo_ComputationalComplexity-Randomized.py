import numpy as np
import math
import random
import os


def Modified_BellmanFord_narrowest(s,t,w,mx,source,dest):
    mx = mx+1
    n=max(t)+1
    ones=np.ones((mx,n))
    distance = 1000*ones
    predecessor =np.zeros((mx,n),dtype=int)-1
    # for m in range(0,mx):
    #     distance[m,source]=0
    distance[:,source]=0
    for i in range(1,mx):
        for j in range(0,len(s)):
            a=s[j]
            b=t[j]
            if max(distance[i-1][a], w[j]) <= distance[i][b]:
                distance[i][b] = max(distance[i-1][a], w[j])
                predecessor[i][b] = a
    d = distance[-1][dest]
    ind = predecessor[-1][dest]
    p = [dest,ind]
    i = mx-2
    while ind != source:
        ind = predecessor[i][ind]
        i = i - 1
        p.append(ind)
    p = np.flip(p)
    return d, p


def Modified_BellmanFord_narrowest_no_numpy(s,t,w,mx,source,dest):
    mx = mx+1
    n=max(t)+1
    ones=np.ones((mx,n))
    distance = 1000*ones
    predecessor =np.zeros((mx,n),dtype=int)-1
    distance[:,source]=0
    distance=distance.tolist()
    predecessor=predecessor.tolist()
    for i in range(1,mx):
        for j in range(0,len(s)):
            a=s[j]
            b=t[j]
            if max(distance[i-1][a], w[j]) <= distance[i][b]:
                distance[i][b] = max(distance[i-1][a], w[j])
                predecessor[i][b] = a
    d = distance[-1][dest]
    ind = predecessor[-1][dest]
    p = [dest,ind]
    i = mx-2
    while ind != source:
        ind = predecessor[i][ind]
        i = i - 1
        p.append(ind)
    p = np.flip(p)
    return d, p


def populate_consecutive_arcs(Mvalues,upperconsec,N):
    Mvalues=np.array(Mvalues)
    for i in range(N-1):
        Mvalues[i][i+1]=random.uniform(0.001,upperconsec)
    return Mvalues


def Take_right_above_diag_from_Vactual(Mvalues,Vactual,N):
    Mvalues=np.array(Mvalues)
    for i in range(N-1):
        Mvalues[i][i+1]=Vactual[i][i+1]
    return Mvalues


def populate_other_arcs_actual(Mvalues,upperremaining,N):
    Mvalues=np.array(Mvalues)
    i=N-2
    while i>=0:
        j=i+2
        while j<N:
            dummyrow,dummycol=Mvalues[i,0: j],Mvalues[i:N, j]
            max1,max2=max(dummyrow),max(dummycol)
            maxneeded=max(max1,max2)
            Mvalues[i][j]=maxneeded + random.uniform(0.001,upperremaining)
            j=j+1
        i=i-1
    return Mvalues


def Create_Vactual(N,upperremaining,upperconsec):
    A = -1 * np.ones((N, N))
    Vactualflags=populate_consecutive_arcs(A,upperconsec,N)
    return populate_other_arcs_actual(Vactualflags,upperremaining,N)


def populate_other_arcs(Mvalues,N):
    Mvalues=np.array(Mvalues)
    i=N-2
    while i>=0:
        j=i+2
        while j<N:
            dummyrow,dummycol=Mvalues[i,0: j],Mvalues[i:N, j]
            max1,max2=max(dummyrow),max(dummycol)
            maxneeded=max(max1,max2)
            Mvalues[i][j]=maxneeded
            j=j+1
        i=i-1
    return Mvalues


def populate_other_arcs_efficient(Mvalues,N):
    Mvalues=np.array(Mvalues)
    i=N-2
    while i>=0:
        Mvalues[0:i+1,i+1:N]=np.maximum(Mvalues[0:i+1,i+1:N], Mvalues[i][i+1])
        i=i-1
    return Mvalues


def extract_weights_needed_for_bellman(Mvalue):
    condition=Mvalue!=-1
    dummy=np.extract(condition,Mvalue)
    return dummy


def Determine_weights_given_pathactual(path,Vactual,M):
    P=[]
    for i in range(len(path)-1):
        a,b=path[i],path[i+1]
        w=Vactual[a][b]
        P.append(w)
        M[a][b]=w
    return P,M


def Update_V_approx_after_exact_calc(path,Mvalue,Mvalueactual,N):
    Mvalue=np.array(Mvalue)
    Mvalueactual=np.array(Mvalueactual)
    for i in range(len(path)-1):
        a,b=path[i],path[i+1]
        Mvalue[a][b]=Mvalueactual[a][b]
    Mvalue=populate_other_arcs(Mvalue,N)
    return Mvalue


def Update_V_approx_after_exact_calc_efficient(path,Mvalue,Mvalueactual,N):
    Mvalue=np.array(Mvalue)
    Mvalueactual=np.array(Mvalueactual)
    for i in range(len(path)-1):
        a,b=path[i],path[i+1]
        Mvalue[a][b]=Mvalueactual[a][b]
        Mvalue[0:a+1, b:N]=np.maximum(Mvalue[0:a+1, b:N],Mvalueactual[a][b])
    return Mvalue


def Create_e1_e2(N):
    e1, e2 = [], []
    for i in range(0, N):
        for j in range(i + 1, N):
            e1.append(i)
            e2.append(j)
    return e1,e2


def vl_max_formula_OPTIMIZED_as_a_lowerbd_on_theworstcase(N):
    K=[math.floor((12.5+2*N)/9),math.ceil((12.5+2*N)/9)]
    d=[K[0]-1,K[1]-1]
    a=[N-7+5.5*d[0]+2*d[0]*N-2*d[0]*K[0]-2.5*d[0]**2 ,N-7+5.5*d[1]+2*d[1]*N-2*d[1]*K[1]-2.5*d[1]**2]
    return int(max(a))


def Dr_Hadi_algo_return_vl_width(epsilon,N,K,Vactual,e1,e2):
    A = -1 * np.ones((N, N))
    Vactualflags = Take_right_above_diag_from_Vactual(A,Vactual,N)
    Vapprox = populate_other_arcs_efficient(Vactualflags, N)
    weights = extract_weights_needed_for_bellman(Vapprox)
    widthapprox, pathapprox = Modified_BellmanFord_narrowest(e1, e2, weights, K, 0, N - 1)
    exactvaluesonpath, Vactualflags = Determine_weights_given_pathactual(pathapprox, Vactual,Vactualflags)
    width = max(exactvaluesonpath)
    matrixtobesavedprogress=[]
    while epsilon>0:
        Vapprox = Update_V_approx_after_exact_calc_efficient(pathapprox, Vapprox, Vactual, N)
        weights = extract_weights_needed_for_bellman(Vapprox)
        widthapprox, pathapprox = Modified_BellmanFord_narrowest(e1, e2, weights, K, 0, N - 1)
        exactvaluesonpath, Vactualflags = Determine_weights_given_pathactual(pathapprox, Vactual,Vactualflags)
        width = max(exactvaluesonpath)
        absolute_gap = abs(width - widthapprox)
        if absolute_gap<=epsilon:
            vl_till_now=len(extract_weights_needed_for_bellman(Vactualflags))
            matrixtobesavedprogress.append([width, widthapprox, (vl_till_now*2)/(N*(N-1))])
            epsilon=absolute_gap
    # vl = len(extract_weights_needed_for_bellman(Vactualflags))
    return width,matrixtobesavedprogress


epsilon = 10000000000
upper_remaining, upper_consecutive = 2, 3
N_array=[25,50,100]
K_allowable_changes_array=np.arange(2,11,1)
nbofreplications=1000

# percentage_vl_upper_bound_array=np.array([(vl_max_formula_OPTIMIZED_as_a_lowerbd_on_theworstcase(j)*100)/((j*(j-1))/2) for j in N_array])
# percentage_vl_upper_bound_array.tofile('lower bound on upper bound for all N.csv',sep=',')

for N in N_array:
    Q = N * (N - 1) / 2
    e1, e2 = Create_e1_e2(N)
    crazy_matrix_of_matrices_K_reps=[[1000000 for j in range(nbofreplications)] for i in range(len(K_allowable_changes_array))]
    for z in range(len(K_allowable_changes_array)):
        K_allowable_changes=K_allowable_changes_array[z]
        for reps in range(nbofreplications):
            V_actual=Create_Vactual(N,upper_remaining,upper_consecutive)
            width,M=Dr_Hadi_algo_return_vl_width(epsilon,N,K_allowable_changes,V_actual,e1,e2)
            # print('the width is',width)
            print('we are in the iteration',N,K_allowable_changes,reps)
            crazy_matrix_of_matrices_K_reps[z][reps]=M
    with open(os.path.join('data for new runs complexity','data when N={}'.format(N)+'.npy'), 'wb') as f:
        np.save(f, crazy_matrix_of_matrices_K_reps,allow_pickle=True, fix_imports=True)




# Read files

# Data_needed=[]
# for i in range(len(N_array)):
#     with open(os.path.join('data for new runs complexity','data when N={}'.format(N_array[i])+'.npy'), 'rb') as f:
#         Data_needed.append(np.load(f,allow_pickle=True, fix_imports=True))