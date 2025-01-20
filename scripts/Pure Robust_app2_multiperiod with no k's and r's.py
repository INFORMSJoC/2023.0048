import math
import time
import cvxpy as cp
import numpy as np
import itertools
from statistics import mean
import operator


def Modified_BellmanFord_narrowest(s,t,w,bs_matrix,mx,source,dest,prevalence_rate_matrix,k,BT):
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
                T_ind1 = predecessor[i][b]
                p_ind1 = [b, T_ind1]
                m = i-1
                while T_ind1 != source:
                    T_ind1 = predecessor[m][T_ind1]
                    m = m - 1
                    p_ind1.append(T_ind1)
                p_ind1 = np.flip(p_ind1)
                Budget_allocation_path1 = []
                for m in range(0, len(p_ind1) - 1):
                    A = p_ind1[m]
                    f = m + 1
                    B = p_ind1[f]
                    for q in range(0, len(s)):
                        if A == s[q]:
                            if B == t[q]:
                                Budget_allocation_path1.append(bs_matrix[q])
                # VIPPPP: I THINK IN WHAT BELOW, the "len(prevalence_rate_matrix[0]" should be up to the node we are testing and maybe " b" (since here we are not getting mulitple optimal sl so does not matter)
                # BTW: it will lead to the same answer by chance, since in both cases we are dividing over "b" the number of nodes up to the testing ones, where the first values here,
                # happened to be non-zero, with the remaining being not 0, but the -(f_cont) values, but since in both cases we are dividing by the same deno, and we are substracting
                # the same values after the "b" in both cases, then it wouldn't really matter by chance
                Budget_allocation_star_constrained_solution1 = [[0] * number_of_viruses] * len(prevalence_rate_matrix[0])
                u = 0
                v = 1
                m = 0
                while (u < len(Budget_allocation_path1) and v < len(p_ind1)):
                    while (p_ind1[v] > m):
                        Budget_allocation_star_constrained_solution1[m] = Budget_allocation_path1[u]
                        m = m + 1
                    v = v + 1
                    u = u + 1
                while m < b:
                    Budget_allocation_star_constrained_solution1[m] = Budget_allocation_path1[-1]
                    m = m + 1
                function_constrained_case1 = []
                Transposeofprevmatrix = np.transpose(prevalence_rate_matrix)
                for m in range(len(Budget_allocation_star_constrained_solution1)):
                    nu = Budget_allocation_star_constrained_solution1[m]
                    mu = Transposeofprevmatrix[m]
                    function_constrained_case1.append(sum(mu[e] * math.exp(-k[e] * nu[e]) for e in range(number_of_viruses)))
                # VIPPPP: In what below I dont think we have to use all of the "Transposeofprevmatrix" however only the rows/cols up to the node we are testing "b"
                # f_cont = ERM_Problem(Transposeofprevmatrix, k, BT, number_of_viruses)
                f_cont = [closed_form_formula_for_ERM(k,Transposeofprevmatrix[x],BT)[0] for x in range(len(Transposeofprevmatrix))]
                avg_diff1=mean(np.abs(list(map(operator.sub,f_cont,function_constrained_case1))))

                # check for the other predecessor
                T_ind2 = a
                p_ind2 = [b, T_ind2]
                m = i-1
                while T_ind2 != source:
                    T_ind2 = predecessor[m][T_ind2]
                    m = m - 1
                    p_ind2.append(T_ind2)
                p_ind2 = np.flip(p_ind2)
                Budget_allocation_path2 = []
                for m in range(0, len(p_ind2) - 1):
                    A = p_ind2[m]
                    f = m + 1
                    B = p_ind2[f]
                    for q in range(0, len(s)):
                        if A == s[q]:
                            if B == t[q]:
                                Budget_allocation_path2.append(bs_matrix[q])
                Budget_allocation_star_constrained_solution2 = [[0] * number_of_viruses] * len(prevalence_rate_matrix[0])
                u = 0
                v = 1
                m = 0
                while (u < len(Budget_allocation_path2) and v < len(p_ind2)):
                    while (p_ind2[v] > m):
                        Budget_allocation_star_constrained_solution2[m] = Budget_allocation_path2[u]
                        m = m + 1
                    v = v + 1
                    u = u + 1
                while m < b:
                    Budget_allocation_star_constrained_solution2[m] = Budget_allocation_path2[-1]
                    m = m + 1
                function_constrained_case2 = []
                Transposeofprevmatrix = np.transpose(prevalence_rate_matrix)
                for m in range(len(Budget_allocation_star_constrained_solution2)):
                    nu = Budget_allocation_star_constrained_solution2[m]
                    mu = Transposeofprevmatrix[m]
                    function_constrained_case2.append(sum(mu[e] * math.exp(-k[e] * nu[e]) for e in range(number_of_viruses)))
                # f_cont = ERM_Problem(Transposeofprevmatrix, k, BT, number_of_viruses)
                f_cont = [closed_form_formula_for_ERM(k, Transposeofprevmatrix[x], BT)[0] for x in range(len(Transposeofprevmatrix))]
                avg_diff2= mean(np.abs(list(map(operator.sub, f_cont, function_constrained_case2))))
                if avg_diff2 == min(avg_diff1,avg_diff2):
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
    # print(predecessor)
    # print(d)
    # print(p)
    return d, p



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
    return objfunc, B



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




Matrix_of_prevalence_rates_of_infections=[[0]*53]*5
Matrix_of_prevalence_rates_of_infections[0]=HIV
Matrix_of_prevalence_rates_of_infections[1]=HBV
Matrix_of_prevalence_rates_of_infections[2]=HCV
Matrix_of_prevalence_rates_of_infections[3]=Babesiosis
Matrix_of_prevalence_rates_of_infections[4]=WNV
Prevalence_Data_List= Matrix_of_prevalence_rates_of_infections


Matrix_of_prevalence_rates_of_infectionsL=[[0]*53]*5
Matrix_of_prevalence_rates_of_infectionsL[0]=HIVL
Matrix_of_prevalence_rates_of_infectionsL[1]=HBVL
Matrix_of_prevalence_rates_of_infectionsL[2]=HCVL
Matrix_of_prevalence_rates_of_infectionsL[3]=BabesiosisL
Matrix_of_prevalence_rates_of_infectionsL[4]=WNVL

# Prevalence_Data_List_Lower= Matrix_of_prevalence_rates_of_infectionsL
Prevalence_Data_List_Lower= np.array(Matrix_of_prevalence_rates_of_infectionsL)/0.9

Matrix_of_prevalence_rates_of_infectionsU=[[0]*53]*5
Matrix_of_prevalence_rates_of_infectionsU[0]=HIVU
Matrix_of_prevalence_rates_of_infectionsU[1]=HBVU
Matrix_of_prevalence_rates_of_infectionsU[2]=HCVU
Matrix_of_prevalence_rates_of_infectionsU[3]=BabesiosisU
Matrix_of_prevalence_rates_of_infectionsU[4]=WNVU

# Prevalence_Data_List_Upper= Matrix_of_prevalence_rates_of_infectionsU
Prevalence_Data_List_Upper= np.array(Matrix_of_prevalence_rates_of_infectionsU)*0.9




BT = 45
k = [0.2800,0.1600,0.1400,0.3800,0.1850]
r=0
number_of_viruses = 5
matrix_column = 2*number_of_viruses
# Prevalence_Data_List=[[0.007, 0.007, 0.007, 0.007, 0.007],[0.00345,0.00345,0.00345,0.00345,0.00345],[0.016,0.016,0.016,0.016,0.016],[0.003059537, 0.001595778, 0.001264081, 0.0010445990,0.000759108],[0.000495, 0.000396902, 0.000108216, 1.34345E-05, 4.96336E-06]]
# Prevalence_Data_List_Upper = [[0.01, 0.01, 0.01, 0.01, 0.01],[0.0044, 0.0044, 0.0044, 0.0044, 0.0044],[0.019, 0.019, 0.019, 0.019, 0.019],[0.003893956, 0.00203099, 0.00160883, 0.00132949, 0.000966138],[0.0015, 0.001202733, 0.000327926, 4.07106E-05, 1.50405E-05]]
# Prevalence_Data_List_Lower = [[0.005, 0.005, 0.005, 0.005, 0.005],[0.0025, 0.0025, 0.0025, 0.0025, 0.0025],[0.013, 0.013, 0.013, 0.013, 0.013],[0.002026447, 0.001056944, 0.000837248, 0.000691877, 0.000502786],[0.000044, 3.52802E-05, 9.61917E-06, 1.19418E-06, 4.41188E-07]]
path_matrix_of_ks=[]
budget_path_3d_matrix_of_ks=[]
F_constrained_matrix_of_ks=[]
MinMaxDistance = []
Total_Average_distance = []
Total_Average_cost_to_save=[]

N=len(Prevalence_Data_List_Upper[0])
T = np.arange(0,len(Prevalence_Data_List[0]),1)
Transpose_of_Prevalence_Data_List = np.transpose(Prevalence_Data_List)
function_continuous_case = [closed_form_formula_for_ERM(k,Transpose_of_Prevalence_Data_List[x],BT)[0] for x in range(len(Transpose_of_Prevalence_Data_List))]
# ERM_Problem(Transpose_of_Prevalence_Data_List, k, BT, number_of_viruses)
# K_allowable_changes_array=[1,2,3,4]
# below is needed for the new plot for app2
# K_allowable_changes_array=np.arange(1,5,1)
K_allowable_changes_array=np.arange(1,13,1) # needed for the plot for varying over K

number_of_arcs = int(N * (N - 1) / 2)
weights = []
B_Regret_star_app_2 = np.zeros((number_of_arcs, number_of_viruses))
e1 = []
e2 = []
source = 0
dest = len(Prevalence_Data_List_Upper[0]) - 1

PREVS = np.zeros((number_of_viruses, 2))
counter=0
for i in range(0, N):
    for j in range(i + 1, N):
        for m in range(number_of_viruses):
            PREVS[m][0], PREVS[m][1] = min(Prevalence_Data_List_Lower[m][i:j + 1]), max(Prevalence_Data_List_Upper[m][i:j + 1])
        w, x = closed_form_formula_for_ERM(k, PREVS[:, 1], BT)
        weights.append(w)
        B_Regret_star_app_2[counter] = x
        e1.append(i)
        e2.append(j)
        print("we are in ", counter)
        counter=counter+1

solutionovertime=[]
for kcounter, K_allowable_changes in enumerate(K_allowable_changes_array):

    d, path = Modified_BellmanFord_narrowest(e1,e2,weights,B_Regret_star_app_2,K_allowable_changes,source,dest,Prevalence_Data_List,k,BT)

    MinMaxDistance.append(d)
    Budget_allocation_path=[]

    for j in range(0,len(path)-1):
        A=path[j]
        s=j+1
        B=path[s]
        for m in range(0,len(e1)):
            if A==e1[m]:
                if B==e2[m]:
                    Budget_allocation_path.append(B_Regret_star_app_2[m])
    print("the allocated B's for the multi-period are")
    print(Budget_allocation_path)

    Budget_allocation_star_constrained_solution=[[0]*number_of_viruses]*len(Prevalence_Data_List[0])
    u=0
    v=1
    m=0
    while(u<len(Budget_allocation_path) and v < len(path)):
        while (path[v] > m):
            Budget_allocation_star_constrained_solution[m]=Budget_allocation_path[u]
            m=m+1
        v=v+1
        u=u+1
    while(m<len(Prevalence_Data_List[0])):
        Budget_allocation_star_constrained_solution[m]=Budget_allocation_path[-1]
        m=m+1
    print('the budget allocation star constrained solution is ', Budget_allocation_star_constrained_solution)


    function_constrained_case = []
    for i in range(len(Budget_allocation_star_constrained_solution)):
        nu = Budget_allocation_star_constrained_solution[i]
        mu = Transpose_of_Prevalence_Data_List[i]
        function_constrained_case.append(sum(mu[j] * math.exp(-k[j] * nu[j]) for j in range(number_of_viruses)))
    function_constrained_case=np.array(function_constrained_case)

    # needed for plotting
    avgdist=mean(np.abs(list(map(operator.sub, function_continuous_case, function_constrained_case))))
    Total_Average_distance.append(avgdist)
    Total_Average_cost_to_save.append(sum(function_constrained_case))

    F_constrained_matrix_of_ks.append(function_constrained_case)
    path=np.array(path)
    path_matrix_of_ks.append(path)
    Budget_allocation_path=np.array(Budget_allocation_path)
    Budget_allocation_path_rounded=np.around(Budget_allocation_path, decimals=1)
    budget_path_3d_matrix_of_ks.append(Budget_allocation_path_rounded)
    solutionovertime.append(Budget_allocation_star_constrained_solution)

    print('======================we are in =====================', kcounter, K_allowable_changes)


import pickle
with open('F_constrained_matrix_of_ks_RobustAPP2.data', 'wb') as f:
    pickle.dump(F_constrained_matrix_of_ks, f)
f.close()
with open('solutionovertime_RobustAPP2.data', 'wb') as f:
    pickle.dump(solutionovertime, f)
f.close()
with open('group_size_path_of_ksRobustAPP2.data', 'wb') as f:
    pickle.dump(budget_path_3d_matrix_of_ks, f)
f.close()
with open('path_matrix_of_ksRobustAPP2.data', 'wb') as f:
    pickle.dump(path_matrix_of_ks, f)
f.close()


# import pickle
# with open('totalcostapp2forplottingRobust.data', 'wb') as f:
#     pickle.dump(Total_Average_cost_to_save, f)
# f.close()