import numpy as np
import itertools
import cvxpy as cp
import math


def Modified_BellmanFord_narrowest(s,t,w,mx,source,dest):
    mx = mx+1
    n=max(t)+1
    ones=np.ones((mx,n))
    distance = 1000*ones
    predecessor =np.zeros((mx,n),dtype=int)-1
    for m in range(0,mx):
        distance[m,source]=0
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


def RRM_Problem(Possible_Prevalence_Vector_Cases, k, BT, number_of_viruses,Risk_values_of_Possible_Cases):
    B_r = cp.Variable(number_of_viruses)
    t = cp.Variable(1)
    RRM_obj_func = t
    constraints_RRM = []
    for m in range(len(Risk_values_of_Possible_Cases)):
        mu = Possible_Prevalence_Vector_Cases[m]
        constraints_RRM.append(sum(mu[j] * cp.exp(-k[j] * B_r[j]) for j in range(number_of_viruses)) - Risk_values_of_Possible_Cases[m] <= t)
    constraints_RRM.append(B_r >= 0)
    constraints_RRM.append(sum(B_r) <= BT)
    RRM_prob = cp.Problem(cp.Minimize(RRM_obj_func), constraints_RRM)
    RRM_prob.solve()
    return B_r.value, t.value


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


def Generate_Cases(V_number_of_viruses_2):
    Possible_Cases = list(itertools.product(*V_number_of_viruses_2))
    return Possible_Cases


def minMaxRegret_multiple_infections(PrevalenceRates_number_of_viruses_2, BT, k):
    Possible_Prevalence_Vector_Cases = Generate_Cases(PrevalenceRates_number_of_viruses_2)
    Risk_values_of_Possible_Cases=[closed_form_formula_for_ERM(k,Possible_Prevalence_Vector_Cases[x],BT) for x in range(len(Possible_Prevalence_Vector_Cases))]
    r=RRM_Problem(Possible_Prevalence_Vector_Cases,k,BT,number_of_viruses,Risk_values_of_Possible_Cases)
    return r


def populate_consecutive_arcs(Mvalues,prevalencerateL,prevalencerateU,N,BT,k,numberofviruses):
    PREVS = np.zeros((numberofviruses, 2))
    Mvalues=np.array(Mvalues)
    for i in range(N-1):
        j=i+1
        for m in range(numberofviruses):
            PREVS[m][0], PREVS[m][1] = min(prevalencerateL[m][i:j + 1]), max(prevalencerateU[m][i:j + 1])
        x, w = minMaxRegret_multiple_infections(PREVS, BT, k)
        Mvalues[i][j]=w[0]
        print(w)
    return Mvalues


# def populate_other_arcs_actual(Mvalues,upperremaining,N):
#     Mvalues=np.array(Mvalues)
#     i=N-2
#     while i>=0:
#         j=i+2
#         while j<N:
#             dummyrow,dummycol=Mvalues[i,0: j],Mvalues[i:N, j]
#             max1,max2=max(dummyrow),max(dummycol)
#             maxneeded=max(max1,max2)
#             Mvalues[i][j]=maxneeded + round(random.uniform(0,upperremaining),4)
#             j=j+1
#         i=i-1
#     return Mvalues


def populate_other_arcs_approx(Mvalues,N):
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


def populate_other_arcs_approx_efficient(Mvalues,N):
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


# def Create_random_initial_path_of_length_Kplus1(K,N):
#     path = [0]
#     for i in range(1, K):
#         path.append(random.randint(path[i - 1] + 1, N - 1 - K + i))
#     path.append(N - 1)
#     return path


# def Determine_weights_given_pathapprox(path,Mvalue):
#     P=[]
#     for i in range(len(path)-1):
#         a,b=path[i],path[i+1]
#         P.append(Mvalue[a][b])
#     return P


def Determine_weights_given_pathactual(path,prevalencerateL,prevalencerateU,BT,k,numberofviruses,M):
    PREVS = np.zeros((numberofviruses, 2))
    P=[]
    for i in range(len(path)-1):
        a,b=path[i],path[i+1]
        for m in range(number_of_viruses):
            PREVS[m][0], PREVS[m][1] = min(prevalencerateL[m][a:b + 1]), max(prevalencerateU[m][a:b + 1])
        x, w = minMaxRegret_multiple_infections(PREVS, BT, k)
        print(w)
        P.append(w[0])
        M[a][b]=w
    return P,M


def Update_V_approx_after_exact_calc(path,Mvalue,pathactual,N):
    Mvalue=np.array(Mvalue)
    for i in range(len(path)-1):
        a,b=path[i],path[i+1]
        Mvalue[a][b]=pathactual[i]
    Mvalue=populate_other_arcs_approx(Mvalue,N)
    return Mvalue


def Update_V_approx_after_exact_calc_efficient(path,Mvalue,pathactual,N):
    Mvalue=np.array(Mvalue)
    for i in range(len(path)-1):
        a,b=path[i],path[i+1]
        Mvalue[a][b]=pathactual[i]
        Mvalue[0:a+1, b:N]=np.maximum(Mvalue[0:a+1, b:N],pathactual[i])
    return Mvalue


def Create_e1_e2(N):
    e1, e2 = [], []
    for i in range(0, N):
        for j in range(i + 1, N):
            e1.append(i)
            e2.append(j)
    return e1,e2


BT = 45
k = [0.2800,0.1600,0.1400,0.3800,0.1850]
number_of_viruses = len(k)
matrix_column = 2*number_of_viruses

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


Prevalence_Data_List=[HIV,HBV,HCV,Babesiosis,WNV]
Prevalence_Data_List_Lower=[HIVL,HBVL,HCVL,BabesiosisL,WNVL]
Prevalence_Data_List_Upper=[HIVU,HBVU,HCVU,BabesiosisU,WNVU]
N=len(Prevalence_Data_List[0])
Q=N*(N-1)/2
T = np.arange(0,len(Prevalence_Data_List[0]),1)
Transpose_of_Prevalence_Data_List = np.transpose(Prevalence_Data_List)

K_allowable_changes=6
number_of_arcs=int(N*(N-1)/2)




epslion = 0.001
e1, e2 = Create_e1_e2(N)
A = -1 * np.ones((N, N))

V_actual_flags=populate_consecutive_arcs(A,Prevalence_Data_List_Lower,Prevalence_Data_List_Upper,N,BT,k,number_of_viruses)

V_approx=populate_other_arcs_approx_efficient(V_actual_flags,N)
weights=extract_weights_needed_for_bellman(V_approx)

widthapprox,pathapprox=Modified_BellmanFord_narrowest(e1,e2,weights,K_allowable_changes,0,N-1)

exact_values_on_path,V_actual_flags=Determine_weights_given_pathactual(pathapprox,Prevalence_Data_List_Lower,Prevalence_Data_List_Upper,BT,k,number_of_viruses,V_actual_flags)
width=max(exact_values_on_path)

ite=0
gap=100000000000
while gap>epslion:
    V_approx=Update_V_approx_after_exact_calc_efficient(pathapprox,V_approx,exact_values_on_path,N)
    weights=extract_weights_needed_for_bellman(V_approx)
    widthapprox,pathapprox=Modified_BellmanFord_narrowest(e1,e2,weights,K_allowable_changes,0,N-1)

    exact_values_on_path,V_actual_flags=Determine_weights_given_pathactual(pathapprox,Prevalence_Data_List_Lower,Prevalence_Data_List_Upper,BT,k,number_of_viruses,V_actual_flags)
    width = max(exact_values_on_path)

    gap=(abs(width-widthapprox))/width

    print('the gap is = ', gap)
    ite=ite+1

vl=len(extract_weights_needed_for_bellman(V_actual_flags))
print('vl is, ', vl)
print(vl*100/Q)