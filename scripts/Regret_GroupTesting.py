from scipy.special import lambertw
import math
from statistics import mean
import numpy as np
import operator

absilun=0.0000000001
r8_big = 1.0E+30
# lambda_2 = 5000
lambda_2 = 60 # this 60 is nothing but the 615 we had for all TTIS but for WNV only

# def BellmanFord_narrowest(s,t,w,mx,source,dest):
#     mx = mx+1
#     n=max(t)+1
#     ones=np.ones((mx,n))
#     distance = 1000*ones
#     predecessor =np.zeros((mx,n),dtype=int)-1
#     for m in range(0,mx):
#         distance[m,source]=0
#     #print(distance)
#     #print(predecessor)
#
#     for i in range(1,mx):
#         for j in range(0,len(s)):
#             a=s[j]
#             b=t[j]
#             if max(distance[i-1][a], w[j]) < distance[i][b]:
#                 distance[i][b] = max(distance[i-1][a], w[j])
#                 predecessor[i][b] = a
#
#     d = distance[-1][dest]
#     ind = predecessor[-1][dest]
#     p = [dest,ind]
#     i = mx-2
#     while ind != source:
#         ind = predecessor[i][ind]
#         i = i - 1
#         p.append(ind)
#
#     p = np.flip(p)
#     print(distance)
#     print(predecessor)
#     print(d)
#     print(p)
#     return d, p

def Modified_BellmanFord_narrowest(s,t,w,ns,mx,source,dest,f_n_0,prevalence_rate_array):
    mx = mx+1
    n=max(t)+1
    ones=np.ones((mx,n))
    distance = 1000*ones
    predecessor =np.zeros((mx,n),dtype=int)-1
    for m in range(0,mx):
        distance[m,source]=0
    #print(distance)
    #print(predecessor)
    for i in range(1,mx):
        for j in range(0,len(s)):
            a=s[j]
            b=t[j]
            if max(distance[i-1][a], w[j]) < distance[i][b]:
                distance[i][b] = max(distance[i-1][a], w[j])
                predecessor[i][b] = a
            if max((distance[i-1][a], w[j])) == distance[i][b] and (predecessor[i][b] != a):
                # old predec from the previous lap (on the same row)
                T_ind1 = predecessor[i][b]
                p_ind1 = [b, T_ind1]
                m = i-1
                # print('m is ', m)
                while T_ind1 != source:
                    T_ind1 = predecessor[m][T_ind1]
                    m = m - 1
                    p_ind1.append(T_ind1)
                p_ind1 = np.flip(p_ind1)
                group_size_path1 = []
                for m in range(0, len(p_ind1) - 1):
                    A = p_ind1[m]
                    k = m + 1
                    B = p_ind1[k]
                    for q in range(0, len(s)):
                        if A == s[q]:
                            if B == t[q]:
                                group_size_path1.append(ns[q])
                n_star_constrained_solution1 = []
                u = 0
                v = 1
                m = 0
                while (u < len(group_size_path1) and v < len(p_ind1)):
                    while (p_ind1[v] > m):
                        n_star_constrained_solution1.append(group_size_path1[u])
                        m = m + 1
                    v = v + 1
                    u = u + 1
                while m < b:
                    n_star_constrained_solution1.append(group_size_path1[-1])
                    m = m + 1
                F_of_n_star_constrained_solution1 = []
                for m in range(len(n_star_constrained_solution1)):
                    F_of_n_star_constrained_solution1.append(f_n_0(prevalence_rate_array[m], n_star_constrained_solution1[m], Se, Sp))
                f_cont = [f_n_0(prevalence_rate_array[u], n_0(prevalence_rate_array[u], Se, Sp), Se, Sp) for u in range(b)]
                avg_diff1 = mean(np.abs(list(map(operator.sub, f_cont, F_of_n_star_constrained_solution1))))

                #check for the other predecessor
                T_ind2 = a
                p_ind2 = [b, T_ind2]
                m = i-1
                # print('m is ', m)
                while T_ind2 != source:
                    T_ind2 = predecessor[m][T_ind2]
                    m = m - 1
                    p_ind2.append(T_ind2)
                p_ind2 = np.flip(p_ind2)
                group_size_path2 = []
                for m in range(0, len(p_ind2) - 1):
                    A = p_ind2[m]
                    k = m + 1
                    B = p_ind2[k]
                    for q in range(0, len(s)):
                        if A == s[q]:
                            if B == t[q]:
                                group_size_path2.append(ns[q])
                n_star_constrained_solution2 = []
                u = 0
                v = 1
                m = 0
                while (u < len(group_size_path2) and v < len(p_ind2)):
                    while (p_ind2[v] > m):
                        n_star_constrained_solution2.append(group_size_path2[u])
                        m = m + 1
                    v = v + 1
                    u = u + 1
                while m < b:
                    n_star_constrained_solution2.append(group_size_path2[-1])
                    m = m + 1
                F_of_n_star_constrained_solution2 = []
                for m in range(len(n_star_constrained_solution2)):
                    F_of_n_star_constrained_solution2.append(
                        f_n_0(prevalence_rate_array[m], n_star_constrained_solution2[m], Se, Sp))
                f_cont = [f_n_0(prevalence_rate_array[u], n_0(prevalence_rate_array[u], Se, Sp), Se, Sp) for u in range(b)]
                avg_diff2=mean(np.abs(list(map(operator.sub, f_cont, F_of_n_star_constrained_solution2))))
                if avg_diff2 < avg_diff1:
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
    print(distance)
    print(predecessor)
    print(d)
    print(p)
    return d, p


def h(p,n):
    y=-n-(1)/((2*math.log(1-p))*(1+lambertw(-0.5*(math.sqrt((math.log(1-p))/((-Se-Sp+1)*(1+lambda_1-lambda_1*Sp)))),0)))
    return y


def partial_der_Regret_wrt_P(p,n):
    y=n*(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*(math.pow(1-p,n-1)) + 1/(2*(1-p)*(lambertw(-0.5*(math.sqrt((math.log(1-p))/((-Se-Sp+1)*(1+lambda_1-lambda_1*Sp)))),0)))
    return y


def n_0(p,Se,Sp):
    y=(2/(math.log(1-p)))*(lambertw(-0.5*(math.sqrt((math.log(1-p))/((-Se-Sp+1)*(1+lambda_1-lambda_1*Sp)))),0))
    return y


def n_1(p,Se,Sp):
    y=(2/(math.log(1-p)))*(lambertw(-0.5*(math.sqrt((math.log(1-p))/((-Se-Sp+1)*(1+lambda_1-lambda_1*Sp)))),-1))
    return y


def f_n_0(p,n,Se,Sp):
    return (1/n)+Se-(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*(1-p)**n + (lambda_2*p*(1-Se*Se)) + (lambda_1*(1-Sp)*Se*(1-p))


def  bisectionH (s1,s2,tolerance,n):
    xl=s1
    xu=s2
    m = (xl + xu) / 2
    while(abs(xl-xu)>tolerance):
        m = (xl+xu)/2
        if(h(xl,n)*h(m,n)<tolerance):
            xu=m
        else:
            if(h(xl,n)*h(m,n)>tolerance):
                xl=m
    return m

def  bisectionP (s1,s2,tolerance,n):
    xl=s1
    xu=s2
    m = (xl + xu) / 2
    while(abs(xl-xu)>tolerance):
        m = (xl+xu)/2
        if(partial_der_Regret_wrt_P(xl,n)*partial_der_Regret_wrt_P(m,n)<tolerance):
            xu=m
        else:
            if(partial_der_Regret_wrt_P(xl,n)*partial_der_Regret_wrt_P(m,n)>tolerance):
                xl=m
    return m


def Regret_Formula(fncn0,p,n,Se,Sp):
    return (1/n)-(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*((1-p)**n) - (1/(fncn0(p,Se,Sp)))+(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*((1-p)**(fncn0(p,Se,Sp)))


def MAX_REGRET_case_1(fncn0,p1,p2,n,Se,Sp):
    if (h(p1,n) * h(p2,n)) > 0:
        if partial_der_Regret_wrt_P(p1,n) * partial_der_Regret_wrt_P(p2,n) > 0:

            y1=(1/n)-(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*(1-p1)**n - (1/(fncn0(p1,Se,Sp)))+(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*(1-p1)**(fncn0(p1,Se,Sp))
            y2=(1/n)-(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*(1-p2)**n - (1/(fncn0(p2,Se,Sp)))+(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*(1-p2)**(fncn0(p2,Se,Sp))
            return max(y1,y2)
        else:
            p_root = bisectionP(p1, p2, absilun,n)
            y1=(1/n)-(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*(1-p1)**n - (1/(fncn0(p1,Se,Sp)))+(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*(1- p1)**(fncn0(p1,Se,Sp))
            y2=(1/n)-(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*(1-p2)**n - (1/(fncn0(p2,Se,Sp)))+(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*(1- p2)**(fncn0(p2,Se,Sp))
            y3=(1/n)-(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*(1-p_root)**n - (1/(fncn0(p_root,Se,Sp)))+(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*(1-p_root)**(fncn0(p_root,Se,Sp))
            return max(y1,y2,y3)
    else:
        p_tilde = bisectionH(p1, p2, absilun,n)
        if partial_der_Regret_wrt_P(p_tilde,n) < 0:

            y1=(1/n)-(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*(1-p1)**n - (1/(fncn0(p1,Se,Sp)))+(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*(1-p1)**(fncn0(p1,Se,Sp))
            y2=(1/n)-(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*(1-p2)**n - (1/(fncn0(p2,Se,Sp)))+(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*(1-p2)**(fncn0(p2,Se,Sp))
            return max(y1,y2)
        if partial_der_Regret_wrt_P(p_tilde,n) == 0:
            if partial_der_Regret_wrt_P(p1,n) * partial_der_Regret_wrt_P((p_tilde - absilun),n) < 0:

                x1 = p_tilde
                x2 = bisectionP(p1, p_tilde - absilun, absilun,n)

                y1=(1/n)-(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*(1-p1)**n - (1/(fncn0(p1,Se,Sp)))+(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*(1-p1)**(fncn0(p1,Se,Sp))
                y2=(1/n)-(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*(1-p2)**n - (1/(fncn0(p2,Se,Sp)))+(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*(1-p2)**(fncn0(p2,Se,Sp))
                y3=(1/n)-(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*(1-x1)**n - (1/(fncn0(x1,Se,Sp)))+(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*(1-x1)**(fncn0(x1,Se,Sp))
                y4=(1/n)-(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*(1-x2)**n - (1/(fncn0(x2,Se,Sp)))+(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*(1-x2)**(fncn0(x2,Se,Sp))
                return max(y1,y2,y3,y4)
            else:
                if partial_der_Regret_wrt_P(p2,n) * partial_der_Regret_wrt_P((p_tilde + absilun),n) < 0:

                    x1 = p_tilde
                    x2 = bisectionP(p_tilde + absilun, p2, absilun,n)

                    y1=(1/n)-(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*(1-p1)**n - (1/(fncn0(p1, Se, Sp)))+(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*(1-p1)**(fncn0(p1,Se,Sp))
                    y2=(1/n)-(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*(1-p2)**n - (1/(fncn0(p2, Se, Sp)))+(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*(1-p2)**(fncn0(p2,Se,Sp))
                    y3=(1/n)-(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*(1-x1)**n - (1/(fncn0(x1,Se,Sp)))+(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*(1-x1)**(fncn0(x1,Se,Sp))
                    y4=(1/n)-(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*(1-x2)**n - (1/(fncn0(x2,Se,Sp)))+(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*(1-x2)**(fncn0(x2,Se,Sp))
                    return max(y1,y2,y3,y4)
                else:

                    x1 = p_tilde
                    y1=(1/n)-(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*(1-p1)**n - (1 /(fncn0(p1,Se,Sp)))+(Se+Sp-1)*(1+lambda_1 -lambda_1*Sp)*(1-p1)**(fncn0(p1,Se,Sp))
                    y2=(1/n)-(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*(1-p2)**n - (1/(fncn0(p2,Se,Sp)))+(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*(1-p2)**(fncn0(p2,Se,Sp))
                    y3=(1/n)-(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*(1-x1)**n - (1/(fncn0(x1,Se,Sp)))+(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*(1-x1)**(fncn0(x1,Se,Sp))
                    return max(y1,y2,y3)
        if partial_der_Regret_wrt_P(p1,n) > 0 and partial_der_Regret_wrt_P(p2,n) > 0:

            y1=(1/n)-(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*(1-p1)**n - (1/(fncn0(p1,Se,Sp)))+(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*(1-p1)**(fncn0(p1,Se,Sp))
            y2=(1/n)-(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*(1-p2)**n - (1/(fncn0(p2,Se,Sp)))+(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*(1-p2)**(fncn0(p2,Se,Sp))
            return max(y1,y2)
        else:
            if partial_der_Regret_wrt_P(p1,n) < 0:
                y1_star = bisectionP(p1, p_tilde, absilun,n)
            else:
                y1_star = bisectionP(p_tilde, p2, absilun,n)
            if partial_der_Regret_wrt_P(p1,n) * partial_der_Regret_wrt_P(y1_star - absilun,n) < 0:
                y1_star_prime = bisectionP(p1, y1_star - absilun, absilun,n)

                x1=y1_star
                x2=y1_star_prime

                y1=(1/n)-(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*(1- p1)**n - (1/(fncn0(p1,Se,Sp)))+(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*(1-p1)**(fncn0(p1,Se,Sp))
                y2=(1/n)-(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*(1- p2)**n - (1/(fncn0(p2,Se,Sp)))+(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*(1-p2)**(fncn0(p2,Se,Sp))
                y3=(1/n)-(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*(1-x1)**n - (1/(fncn0(x1,Se,Sp)))+(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*(1-x1)**(fncn0(x1,Se,Sp))
                y4=(1/n)-(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*(1-x2)**n - (1/(fncn0(x2,Se,Sp)))+(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*(1-x2)**(fncn0(x2,Se,Sp))
                return max(y1,y2,y3,y4)
            else:
                if partial_der_Regret_wrt_P(y1_star + absilun,n) * partial_der_Regret_wrt_P(p2,n) < 0:
                    y1_star_prime = bisectionP(y1_star + absilun, p2, absilun,n)

                    x1 = y1_star
                    x2 = y1_star_prime

                    y1=(1/n)-(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*(1- p1)**n - (1/(fncn0(p1,Se,Sp)))+(Se + Sp - 1)*(1+lambda_1-lambda_1*Sp)*(1- p1)**(fncn0(p1,Se,Sp))
                    y2=(1/n)-(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*(1- p2)**n - (1/(fncn0(p2,Se,Sp)))+(Se + Sp - 1)*(1+lambda_1-lambda_1*Sp)*(1- p2)**(fncn0(p2,Se,Sp))
                    y3=(1/n)-(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*(1-x1)**n - (1/(fncn0(x1,Se,Sp)))+(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*(1-x1)**(fncn0(x1,Se,Sp))
                    y4=(1/n)-(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*(1-x2)**n - (1/(fncn0(x2,Se,Sp)))+(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*(1-x2)**(fncn0(x2,Se,Sp))
                    return max(y1,y2,y3,y4)
                else:
                    y1=(1/n)-(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*(1 - p1)**n - (1/(fncn0(p1,Se,Sp)))+(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*(1-p1)**(fncn0(p1,Se,Sp))
                    y2=(1/n)-(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*(1 - p2)**n - (1/(fncn0(p2,Se,Sp)))+(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*(1-p2)**(fncn0(p2,Se,Sp))
                    y3=(1/n)-(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*(1-y1_star)**n - (1/(fncn0(y1_star,Se,Sp)))+(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*(1-y1_star)**(fncn0(y1_star,Se,Sp))
                    return max(y1,y2,y3)


r=0
number_for_the_subplot=1
Prevalence_Data_List = [0.000023472,0.000039521,0.0000644458,0.00010147,0.000153794,0.00023497,0.000315464,0.000402724,0.000470892,0.000495,0.000491618,0.000448873,0.000377572,0.000294762,0.000215901,0.000143624,0.0000954769,0.0000613813,0.0000391124,0.0000268997]
# Prevalence_Data_List_Upper=[0.0000514848,0.000168844,0.000286203,0.000403562,0.000520922,0.000716737,0.000912553,0.001108369,0.001304184,0.0015,0.001425683,0.001351367,0.00127705,0.001202733,0.000984032,0.00076533,0.000546628,0.000327926,0.000270483,0.00021304]
# Prevalence_Data_List_Lower=[0.00000151022,0.00000495276,0.0000083953,0.0000118378,0.0000152804,0.0000210243,0.0000267682,0.0000325121,0.0000382561,0.000044,0.00004182,0.0000396401,0.0000374601,0.0000352802,0.0000288649,0.0000224497,0.0000160344,0.00000961917,0.00000793417,0.00000624917]
Prevalence_Data_List_Upper=[3.74784e-05, 0.00010418250000000001, 0.0001753244, 0.000252516, 0.00033735800000000006, 0.0004758535, 0.0006140085, 0.0007555465, 0.000887538, 0.0009975000000000001, 0.0009586505000000001, 0.00090012, 0.0008273110000000001, 0.0007487475000000001, 0.0005999664999999999, 0.000454477, 0.00032105245, 0.00019465364999999999, 0.0001547977, 0.00011996985000000001]
Prevalence_Data_List_Lower=[1.249111e-05, 2.2236880000000002e-05, 3.642055e-05, 5.66539e-05, 8.453720000000001e-05, 0.00012799715, 0.0001711161, 0.00021761804999999998, 0.00025457405, 0.0002695, 0.000266719, 0.00024425655, 0.00020751605, 0.0001650211, 0.00012238295, 8.303685e-05, 5.575565e-05, 3.5500235e-05, 2.3523285e-05, 1.6574435e-05]



T = np.arange(0,len(Prevalence_Data_List),1)
# Prevalence_Data_List_Upper=[(1+r)*i for i in Prevalence_Data_List]
# Prevalence_Data_List_Lower=[(1-r)*i for i in Prevalence_Data_List]
Se=0.95
Sp=0.95
lambda_1=30
p_underline = 1 - math.exp(-(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*0.367879)
p_hat = 1 - math.exp(-(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*0.541341)
print('p hat is ',p_hat)
print('p underline is ',p_underline)
N=len(Prevalence_Data_List)
# K_allowable_changes_array=np.arange(1,5,1)
K_allowable_changes_array=np.arange(1,7,1) # needed for the plot for varying over K



n_star_continuous = [n_0(i, Se, Sp) for i in Prevalence_Data_List]
MinMaxDistance=[]
Total_Average_distance=[]
Total_Average_cost_to_save=[]
MinMaxDistance_from16=[]
Total_Average_distance_from16=[]
F_of_n_equal_16=[]
F_of_n_star_continuous = []
for i in range(len(Prevalence_Data_List)):
    ceil = math.ceil(n_star_continuous[i])
    floor = math.floor(n_star_continuous[i])
    min_between_ceil_and_floor = min(f_n_0(Prevalence_Data_List[i], ceil, Se, Sp),
                                     f_n_0(Prevalence_Data_List[i], floor, Se, Sp))
    F_of_n_star_continuous.append(min_between_ceil_and_floor)
    F_of_n_equal_16.append(f_n_0(Prevalence_Data_List[i],16,Se,Sp))
print('the n star continuous solution is ', n_star_continuous)


F_constrained_matrix_of_ks=[[0]*len(Prevalence_Data_List)]*len(K_allowable_changes_array)
path_matrix_of_ks=[]
group_size_path_of_ks=[]
solutionovertime=[]
for K_allowable_changes in K_allowable_changes_array:
    weights=[]
    e1=[]
    e2=[]
    N_star_List = []
    for i in range(0,N):
        sublist_upper = []
        sublist_lower = []
        sublist_upper.append(Prevalence_Data_List_Upper[i])
        sublist_lower.append(Prevalence_Data_List_Lower[i])
        for j in range(i+1,N):
            sublist_upper.append(Prevalence_Data_List_Upper[j])
            sublist_lower.append(Prevalence_Data_List_Lower[j])
            a=min(sublist_lower)
            b=max(sublist_upper)
            p_underline = 1 - math.exp(-(Se + Sp - 1) * (1 + lambda_1 - lambda_1 * Sp) * 0.367879)
            p_hat = 1 - math.exp(-(Se + Sp - 1) * (1 + lambda_1 - lambda_1 * Sp) * 0.541341)
            if b < p_underline:
                V = min(p_hat, b)
                Lower_Bound = math.ceil(n_0(V, Se, Sp))
                Upper_Bound = min(400,math.ceil(n_1(a, Se, Sp)))
                n_possible_values = list(range(Lower_Bound, Upper_Bound + 1, 1))
                W = []
                for m in n_possible_values:
                    W.append(MAX_REGRET_case_1(n_0, a, b, m, Se, Sp))

                Limit=(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*(1-a)**(n_0(a,Se,Sp)) - 1/(n_0(a,Se,Sp))
                min_REGRET = min(min(W), Limit)

                if min_REGRET == Limit:
                    MIN = "infinity"
                else:
                    Index = W.index(min_REGRET)
                    MIN = n_possible_values[Index]
                N_star_List.append(MIN)
                Z=complex(min_REGRET,0)
                weights.append(Z.real)
                e1.append(i)
                e2.append(j)
            if p_underline < b:
                V = min(p_hat, p_underline)
                Lower_Bound = math.ceil(n_0(V, Se, Sp))
                Upper_Bound = min(400,math.ceil(n_1(a, Se, Sp)))
                n_possible_values = list(range(Lower_Bound, Upper_Bound + 1, 1))
                W = []
                for m in n_possible_values:
                    W.append(max(MAX_REGRET_case_1(n_0, a, p_underline, m, Se, Sp),(1/m)-(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*(1-p_underline)** m))

                Limit=(Se+Sp-1)*(1+lambda_1-lambda_1*Sp)*(1-a)**(n_0(a,Se,Sp)) - 1/(n_0(a,Se,Sp))
                min_REGRET = min(min(W), Limit)

                if min_REGRET == Limit:
                    MIN = "infinity"
                else:
                    Index = W.index(min_REGRET)
                    MIN = n_possible_values[Index]
                N_star_List.append(MIN)
                Z = complex(min_REGRET, 0)
                weights.append(Z.real)
                e1.append(i)
                e2.append(j)
    source=0
    dest=max(e2)
    print(e1)
    print(e2)

    print(weights)
    print(N_star_List)
    d, path = Modified_BellmanFord_narrowest(e1,e2,weights,N_star_List,K_allowable_changes,source,dest,f_n_0,Prevalence_Data_List)

    group_size_path=[]

    for j in range(0,len(path)-1):
        A=path[j]
        k=j+1
        B=path[k]
        for m in range(0,len(e1)):
            if A==e1[m]:
                if B==e2[m]:
                    group_size_path.append(N_star_List[m])
    print(group_size_path)


    n_star_constrained_solution=[]
    u=0
    v=1
    m=0
    while(u<len(group_size_path) and v < len(path)):
        while (path[v] > m):
            n_star_constrained_solution.append(group_size_path[u])
            m=m+1
        v=v+1
        u=u+1
    while(m<len(Prevalence_Data_List)):
        n_star_constrained_solution.append(group_size_path[-1])
        m=m+1
    print('the n star constrained solution is ', n_star_constrained_solution)
    n_star_constrained_solution = n_star_constrained_solution
    F_of_n_star_constrained_solution = []
    for i in range(len(Prevalence_Data_List)):
        F_of_n_star_constrained_solution.append(f_n_0(Prevalence_Data_List[i],n_star_constrained_solution[i],Se,Sp))
    #down are for plotting
    MinMaxDistance.append(round(d,5))
    Total_Average_distance.append(round(mean(np.abs(list(map(operator.sub, F_of_n_star_continuous, F_of_n_star_constrained_solution)))),5))
    Total_Average_cost_to_save.append(sum(F_of_n_star_constrained_solution))
    MinMaxDistance_from16.append(round(max(list(map(operator.sub, F_of_n_equal_16,F_of_n_star_constrained_solution))),5))
    Total_Average_distance_from16.append(round(mean(np.abs(list(map(operator.sub, F_of_n_equal_16, F_of_n_star_constrained_solution)))),5))
    F_constrained_matrix_of_ks[K_allowable_changes-1]=F_of_n_star_constrained_solution
    path=np.array(path)
    path_matrix_of_ks.append(path)
    group_size_path=np.array(group_size_path)
    group_size_path_of_ks.append(group_size_path)
    solutionovertime.append(n_star_constrained_solution)

import pickle
with open('F_constrained_matrix_of_ks_RegretAPP1.data', 'wb') as f:
    pickle.dump(F_constrained_matrix_of_ks, f)
f.close()
with open('solutionovertime_RegretAPP1.data', 'wb') as f:
    pickle.dump(solutionovertime, f)
f.close()
with open('group_size_path_of_ksRegretAPP1.data', 'wb') as f:
    pickle.dump(group_size_path_of_ks, f)
f.close()
with open('path_matrix_of_ksRegretAPP1.data', 'wb') as f:
    pickle.dump(path_matrix_of_ks, f)
f.close()


# # distance of CP to ideal "minmax and average"
# MinMaxDistance_current_practice_to_ideal = round(max(list(map(operator.sub, F_of_n_equal_16,F_of_n_star_continuous))),5)
# Total_Average_distance_current_practice_to_ideal = round(mean(np.abs(list(map(operator.sub, F_of_n_equal_16,F_of_n_star_continuous)))),5)


# save the dictio in a json file
# import pickle
# with open('totalcostapp1forplotting.data', 'wb') as f:
#     pickle.dump(Total_Average_cost_to_save, f)
# f.close()