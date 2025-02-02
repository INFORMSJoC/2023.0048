import math
import time
import cvxpy as cp
import numpy as np
import itertools
import warnings
from statistics import mean
import operator
warnings.filterwarnings("ignore")


def CreateFeasibleValuesofthevectoru(levelofdescofB,number_of_viruses,B):
    possiblevaluesofb = np.array([np.linspace(0, B, levelofdescofB) for m in range(number_of_viruses)])
    alluvalues = np.array(list(itertools.product(*possiblevaluesofb)))
    sums = np.array([sum(j) for j in alluvalues]).reshape(len(alluvalues), 1)
    alluvalueswithsum = np.hstack((alluvalues, sums))
    alluvalueswithsum = alluvalueswithsum[alluvalueswithsum[:, -1] <= B]
    return alluvalueswithsum[:, :-1]


def CreatePossibleValuesofthevectorp(levelofdescofP,number_of_viruses, pl, pu):
    possiblevaluesofp = np.array([np.linspace(pl[m],pu[m],levelofdescofP) for m in range(number_of_viruses)])
    return np.array(list(itertools.product(*possiblevaluesofp)))


def NLPFormulation(k, number_of_viruses,B, levelofdescofB, pl,pu, levelofdescofP, timehorizon, K):
    Feasibleuvalues=CreateFeasibleValuesofthevectoru(levelofdescofB,number_of_viruses,B)
    z=cp.Variable(1)
    objfunc=z
    x = [cp.Variable(number_of_viruses) for t in range(timehorizon+1)]
    y = cp.Variable(timehorizon+1, integer=False)
    constraints=[]
    for t in range(timehorizon+1):
        Allpvalues = CreatePossibleValuesofthevectorp(levelofdescofP, number_of_viruses, pl[:, t], pu[:, t])
        print('t=',t)
        constraints.append(sum(x[t]) <= B)
        constraints.append(x[t] >= 0)
        constraints.append(y[t] >= 0)
        constraints.append(y[t] <= 1)
        if t < timehorizon:
            for m in range(number_of_viruses):
                constraints.append(x[t+1][m]>=x[t][m]-B*y[t])
                constraints.append(x[t+1][m]<=x[t][m]+B*y[t])
        for p in range(len(Allpvalues)):
            mu=Allpvalues[p]
            for u in range(len(Feasibleuvalues)):
                U=Feasibleuvalues[u]
                f_u_p=sum(mu[j] * math.exp(-k[j] * U[j]) for j in range(number_of_viruses))
                f_x_p=sum(mu[j] * cp.exp(-k[j] * x[t][j]) for j in range(number_of_viruses))
                constraints.append(f_x_p - f_u_p <= z)
    constraints.append(sum(y[t] for t in range(timehorizon+1)) <= K)
    NLP = cp.Problem(cp.Minimize(objfunc), constraints)
    NLP.solve()
    return z.value


def Modified_BellmanFord_narrowest(s,t,w,bs_matrix,mx,source,dest,ERM_Problem,prevalence_rate_matrix,k,BT):
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

            # if max((distance[i-1][a], w[j])) == distance[i][b] and (predecessor[i][b] != a):
            #     T_ind1 = predecessor[i][b]
            #     p_ind1 = [b, T_ind1]
            #     m = i-1
            #     while T_ind1 != source:
            #         T_ind1 = predecessor[m][T_ind1]
            #         m = m - 1
            #         p_ind1.append(T_ind1)
            #     p_ind1 = np.flip(p_ind1)
            #     Budget_allocation_path1 = []
            #     for m in range(0, len(p_ind1) - 1):
            #         A = p_ind1[m]
            #         f = m + 1
            #         B = p_ind1[f]
            #         for q in range(0, len(s)):
            #             if A == s[q]:
            #                 if B == t[q]:
            #                     Budget_allocation_path1.append(bs_matrix[q])
            #     # VIPPPP: I THINK IN WHAT BELOW, the "len(prevalence_rate_matrix[0]" should be up to the node we are testing and maybe " b" (since here we are not getting mulitple optimal sl so does not matter)
            #     # BTW: it will lead to the same answer by chance, since in both cases we are dividing over "b" the number of nodes up to the testing ones, where the first values here,
            #     # happened to be non-zero, with the remaining being not 0, but the -(f_cont) values, but since in both cases we are dividing by the same deno, and we are substracting
            #     # the same values after the "b" in both cases, then it wouldn't really matter by chance
            #     Budget_allocation_star_constrained_solution1 = [[0] * number_of_viruses] * len(prevalence_rate_matrix[0])
            #     u = 0
            #     v = 1
            #     m = 0
            #     while (u < len(Budget_allocation_path1) and v < len(p_ind1)):
            #         while (p_ind1[v] > m):
            #             Budget_allocation_star_constrained_solution1[m] = Budget_allocation_path1[u]
            #             m = m + 1
            #         v = v + 1
            #         u = u + 1
            #     while m < b:
            #         Budget_allocation_star_constrained_solution1[m] = Budget_allocation_path1[-1]
            #         m = m + 1
            #     function_constrained_case1 = []
            #     Transposeofprevmatrix = np.transpose(prevalence_rate_matrix)
            #     for m in range(len(Budget_allocation_star_constrained_solution1)):
            #         nu = Budget_allocation_star_constrained_solution1[m]
            #         mu = Transposeofprevmatrix[m]
            #         function_constrained_case1.append(sum(mu[e] * math.exp(-k[e] * nu[e]) for e in range(number_of_viruses)))
            #     # VIPPPP: In what below I dont think we have to use all of the "Transposeofprevmatrix" however only the rows/cols up to the node we are testing "b"
            #     f_cont = ERM_Problem(Transposeofprevmatrix, k, BT, number_of_viruses)
            #     avg_diff1=mean(np.abs(list(map(operator.sub,f_cont,function_constrained_case1))))
            #
            #     # check for the other predecessor
            #     T_ind2 = a
            #     p_ind2 = [b, T_ind2]
            #     m = i-1
            #     while T_ind2 != source:
            #         T_ind2 = predecessor[m][T_ind2]
            #         m = m - 1
            #         p_ind2.append(T_ind2)
            #     p_ind2 = np.flip(p_ind2)
            #     Budget_allocation_path2 = []
            #     for m in range(0, len(p_ind2) - 1):
            #         A = p_ind2[m]
            #         f = m + 1
            #         B = p_ind2[f]
            #         for q in range(0, len(s)):
            #             if A == s[q]:
            #                 if B == t[q]:
            #                     Budget_allocation_path2.append(bs_matrix[q])
            #     Budget_allocation_star_constrained_solution2 = [[0] * number_of_viruses] * len(prevalence_rate_matrix[0])
            #     u = 0
            #     v = 1
            #     m = 0
            #     while (u < len(Budget_allocation_path2) and v < len(p_ind2)):
            #         while (p_ind2[v] > m):
            #             Budget_allocation_star_constrained_solution2[m] = Budget_allocation_path2[u]
            #             m = m + 1
            #         v = v + 1
            #         u = u + 1
            #     while m < b:
            #         Budget_allocation_star_constrained_solution2[m] = Budget_allocation_path2[-1]
            #         m = m + 1
            #     function_constrained_case2 = []
            #     Transposeofprevmatrix = np.transpose(prevalence_rate_matrix)
            #     for m in range(len(Budget_allocation_star_constrained_solution2)):
            #         nu = Budget_allocation_star_constrained_solution2[m]
            #         mu = Transposeofprevmatrix[m]
            #         function_constrained_case2.append(
            #             sum(mu[e] * math.exp(-k[e] * nu[e]) for e in range(number_of_viruses)))
            #     f_cont = ERM_Problem(Transposeofprevmatrix, k, BT, number_of_viruses)
            #     avg_diff2= mean(np.abs(list(map(operator.sub, f_cont, function_constrained_case2))))
            #     if avg_diff2 == min(avg_diff1,avg_diff2):
            #         predecessor[i][b] = a

    d = distance[-1][dest]
    ind = predecessor[-1][dest]
    p = [dest,ind]
    i = mx-2
    while ind != source:
        ind = predecessor[i][ind]
        i = i - 1
        p.append(ind)

    p = np.flip(p)
    # print(distance)
    # print(predecessor)
    # print(d)
    # print(p)
    return d, p


def Generate_Cases(V_number_of_viruses_2):
    Possible_Cases = list(itertools.product(*V_number_of_viruses_2))
    return Possible_Cases

def ERM_Problem(Possible_Prevalence_Vector_Cases, k, BT, number_of_viruses):
    Risk_values_of_Possible_Cases = []
    for i in range(len(Possible_Prevalence_Vector_Cases)):
        mu = Possible_Prevalence_Vector_Cases[i]
        B_d = cp.Variable(number_of_viruses)
        ERM_obj_func = sum(mu[j] * cp.exp(-k[j] * B_d[j]) for j in range(number_of_viruses))
        constraints_ERM = [sum(B_d) <= BT, B_d >= 0]
        ERM_prob = cp.Problem(cp.Minimize(ERM_obj_func), constraints_ERM)
        ERM_prob.solve()
        #        print("\nThe optimal risk is", ERM_prob.value)
        #        print("The optimal B_d is")
        #        print(B_d.value)
        #        print(sum(B_d.value))
        Risk_values_of_Possible_Cases.append(ERM_prob.value)
    return Risk_values_of_Possible_Cases

def RRM_Problem(Possible_Prevalence_Vector_Cases, k, BT, number_of_viruses,Risk_values_of_Possible_Cases):
    B_r = cp.Variable(number_of_viruses)
    t = cp.Variable(1)
    RRM_obj_func = t
    constraints_RRM = []
    for m in range(len(Risk_values_of_Possible_Cases)):
        mu = Possible_Prevalence_Vector_Cases[m]
        constraints_RRM.append(
            sum(mu[j] * cp.exp(-k[j] * B_r[j]) for j in range(number_of_viruses)) - Risk_values_of_Possible_Cases[m] <= t)
    #    for m in range(number_of_viruses):
    #        constraints_RRM.append(B_r[m] >= 0)
    constraints_RRM.append(B_r >= 0)
    constraints_RRM.append(sum(B_r) <= BT)
    RRM_prob = cp.Problem(cp.Minimize(RRM_obj_func), constraints_RRM)
    RRM_prob.solve()
    # print("\nThe optimal regret is", RRM_prob.value)
    print("The optimal B_r is")
    # print(B_r.value)
    # print(sum(B_r.value))
    return B_r.value, t.value

def minMaxRegret_multiple_infections(PrevalenceRates_number_of_viruses_2, BT, k):
    Possible_Prevalence_Vector_Cases = Generate_Cases(PrevalenceRates_number_of_viruses_2)
    number_of_viruses = len(PrevalenceRates_number_of_viruses_2)
    Risk_values_of_Possible_Cases = ERM_Problem(Possible_Prevalence_Vector_Cases, k, BT, number_of_viruses)
    t = time.time()
    r =  RRM_Problem(Possible_Prevalence_Vector_Cases, k, BT, number_of_viruses,Risk_values_of_Possible_Cases)
    return r



HIVL = [0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005]
HIV = [0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007,0.007]
HIVU = [0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01]
HBVL = [0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025]
HBV = [0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345,0.00345]
HBVU = [0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044,0.0044]
HCVL = [0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013]
HCV = [0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016,0.016]
HCVU = [0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019,0.019]
BabesiosisL = [0.000433926,0.000390934,0.000347943,0.000304951,0.00033009,0.000355229,0.000380369,0.000405508,0.000423816,0.000442124,0.000460432,0.00047874,0.000518963,0.000559186,0.000599408,0.000639631,0.000679854,0.000833696,0.000987538,0.001141379,0.001295221,0.001608916,0.00192261,0.002236305,0.00255,0.002445289,0.002340579,0.002235868,0.002131157,0.002026447,0.001784071,0.001541695,0.00129932,0.001056944,0.00100202,0.000947096,0.000892172,0.000837248,0.000808174,0.0007791,0.000750026,0.000720952,0.000691877,0.000644605,0.000597332,0.000550059,0.000502786,0.000455513,0.00040824,0.000360968,0.000313695,0.000266422,0.000219149]
Babesiosis = [0.000655144,0.000590234,0.000525325,0.000460416,0.000498371,0.000536327,0.000574282,0.000612237,0.000639879,0.00066752,0.000695162,0.000722803,0.000783532,0.000844261,0.000904989,0.000965718,0.001026447,0.001258717,0.001490988,0.001723259,0.001955529,0.002429147,0.002902765,0.003376382,0.00385,0.003691907,0.003533815,0.003375722,0.00321763,0.003059537,0.002693597,0.002327658,0.001961718,0.001595778,0.001512854,0.001429929,0.001347005,0.001264081,0.001220184,0.001176288,0.001132392,0.001088495,0.001044599,0.000973227,0.000901854,0.000830481,0.000759108,0.000687736,0.000616363,0.00054499,0.000473618,0.000402245,0.000330872]
BabesiosisU = [0.000833819,0.000751207,0.000668596,0.000585984,0.000634291,0.000682598,0.000730904,0.000779211,0.000814391,0.000849571,0.000884751,0.000919931,0.000997222,0.001074514,0.001151805,0.001229096,0.001306387,0.001602004,0.001897621,0.002193238,0.002488856,0.003091642,0.003694428,0.004297214,0.0049,0.004698791,0.004497583,0.004296374,0.004095165,0.003893956,0.003428215,0.002962473,0.002496732,0.00203099,0.00192545,0.00181991,0.00171437,0.00160883,0.001552962,0.001497094,0.001441226,0.001385358,0.00132949,0.001238652,0.001147814,0.001056976,0.000966138,0.0008753,0.000784462,0.000693624,0.000602786,0.000511948,0.00042111]
WNVL = [0.0000000296953,0.0000000240391,0.0000000183828,0.0000000127266,0.0000000222715,0.0000000318164,0.0000000413614,0.0000000509063,0.0000000593907,0.000000067875,0.0000000763594,0.0000000848438,0.000000129387,0.00000017393,0.000000218473,0.000000263016,0.000000307559,0.000000608224,0.000000908889,0.00000120955,0.00000151022,0.00000495276,0.0000083953,0.0000118378,0.0000152804,0.0000210243,0.0000267682,0.0000325121,0.0000382561,0.000044,0.00004182,0.0000396401,0.0000374601,0.0000352802,0.0000288649,0.0000224497,0.0000160344,0.00000961917,0.00000793417,0.00000624917,0.00000456417,0.00000287917,0.00000119418,0.00000104358,0.000000892981,0.000000742383,0.000000591786,0.000000441188,0.00000029059,0.000000139992,0.00000014,0.00000014,0.00000014]
WNV = [0.000000334073,0.00000027044,0.000000206807,0.000000143174,0.000000250554,0.000000357935,0.000000465315,0.000000572696,0.000000668145,0.000000763594,0.000000859044,0.000000954493,0.0000014556,0.00000195671,0.00000245782,0.00000295893,0.00000346004,0.00000684252,0.000010225,0.0000136075,0.00001699,0.0000557185,0.0000944471,0.000133176,0.000171904,0.000236523,0.000301142,0.000365762,0.000430381,0.000495,0.000470475,0.000445951,0.000421426,0.000396902,0.00032473,0.000252559,0.000180387,0.000108216,0.0000892594,0.0000703032,0.0000513469,0.0000323907,0.0000134345,0.0000117403,0.000010046,0.00000835181,0.00000665759,0.00000496336,0.00000326914,0.00000157491,0.00000157,0.00000157,0.00000157]
WNVU = [0.00000101234,0.000000819514,0.000000626687,0.00000043386,0.000000759256,0.00000108465,0.00000141005,0.00000173544,0.00000202468,0.00000231392,0.00000260316,0.0000028924,0.00000441091,0.00000592943,0.00000744794,0.00000896645,0.000010485,0.0000207349,0.0000309849,0.0000412348,0.0000514848,0.000168844,0.000286203,0.000403562,0.000520922,0.000716737,0.000912553,0.001108369,0.001304184,0.0015,0.001425683,0.001351367,0.00127705,0.001202733,0.000984032,0.00076533,0.000546628,0.000327926,0.000270483,0.00021304,0.000155597,0.0000981537,0.0000407106,0.0000355766,0.0000304425,0.0000253085,0.0000201745,0.0000150405,0.00000990648,0.00000477246,0.00000477,0.00000477,0.00000477]

time_horizon=5
BT = 45
pl=np.array([HIVL,HBVL,HCVL,BabesiosisL,WNVL])
pu=np.array([HIVU,HBVL,HCVU,BabesiosisU,WNVU])
Prevalence_Data_List=np.array([Babesiosis])
k = [0.2800,0.1600,0.1400,0.3800,0.1850]
# k=[0.2800,0.1600,0.1400,0.3800]
number_of_viruses=len(k)
levelofdescofB=45
levelofdescofP=100
allowablenbchanges=4



t1=time.time()
optobjfunc_z1=NLPFormulation(k, number_of_viruses, BT, levelofdescofB, pl,pu,levelofdescofP,time_horizon,allowablenbchanges)
print('time needed MINLP is', time.time()-t1)












t1=time.time()
N=time_horizon
number_of_arcs = int(N * (N - 1) / 2)
weights = []
e1 = []
e2 = []
source = 0
dest = N - 1
matrix_column = 2*number_of_viruses
A = np.zeros((number_of_arcs, matrix_column))

p = 0
q = 1
for m in range(number_of_viruses):
    v = 0
    for i in range(0, N):
        sublist_upper = []
        sublist_lower = []
        sublist_upper.append(pu[m][i])
        sublist_lower.append(pl[m][i])
        for j in range(i + 1, N):
            sublist_upper.append(pu[m][j])
            sublist_lower.append(pl[m][j])
            a = min(sublist_lower)
            b = max(sublist_upper)
            A[v][p] = a
            A[v][q] = b
            v = v + 1
    p = p + 2
    q = q + 2
# print(A)
PREVS = np.zeros((number_of_viruses, 2))
B_Regret_star_app_2 = np.zeros((number_of_arcs, number_of_viruses))
counter=0
for i in range(number_of_arcs):
    for m in range(number_of_viruses):
        Q = A[i]
        p = 2 * m
        q = 2 * m + 1
        PREVS[m][0] = Q[p]
        PREVS[m][1] = Q[q]
    # print(PREVS)
    B, maxregret = minMaxRegret_multiple_infections(PREVS, BT, k)
    weights.append(maxregret)
    B_Regret_star_app_2[i] = B
    counter=counter+1
    print('we are in counter ', counter)
# print("the B values to each virus is")
# print(B_Regret_star_app_2)

for i in range(0, N):
    for j in range(i + 1, N):
        e1.append(i)
        e2.append(j)
optobjfunc_z2,_= Modified_BellmanFord_narrowest(e1,e2,weights,B_Regret_star_app_2,allowablenbchanges-1,source,dest,ERM_Problem,Prevalence_Data_List,k,BT)
print('time needed OPD is', time.time()-t1)