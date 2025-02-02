import numpy as np
import os
import matplotlib.pylab as plt
import pandas as pd



def Create_gap_from_the_given_data(AllData,Karray,nbofreps):
    AllData = np.array(AllData)
    for i in range(len(AllData)):
        for j in range(len(Karray)):
            for m in range(nbofreps):
                AllData[i][j][m] = np.array(AllData[i][j][m])
                a0, a1 = np.array(AllData[i][j][m][:, 0]), np.array(AllData[i][j][m][:, 1])
                AllData[i][j][m][:, 0] = np.abs(a0 - a1)/a0
                AllData[i][j][m] = np.delete(AllData[i][j][m], 1, 1)
                AllData[i][j][m][:, 0],AllData[i][j][m][:, 1]=AllData[i][j][m][:, 0][::-1],AllData[i][j][m][:, 1][::-1]
    return AllData


def Determine_before_and_after(a,array):
    for i in range(len(array)-1):
        if a>=array[i] and a<=array[i+1]:
            return array[i],array[i+1],i,i+1
    if i==len(array)-2:
        return 10000,10000,10000,10000


def Interpolation(arrayofgaps,matrixof2cols,N):
    matrixof2cols = np.array(matrixof2cols)
    d = []
    for i in range(len(gap_array)):
        gap = arrayofgaps[i]
        before, after, indexbefore, indexafter = Determine_before_and_after(gap, matrixof2cols[:, 0])
        if before != 10000:
            percentvlbefore, percentvlafter = matrixof2cols[indexbefore][1], matrixof2cols[indexafter][1]
            interpolated = percentvlbefore - (percentvlbefore - percentvlafter) * (before - gap) / (before - after)
            d.append(interpolated)
        else:
            break
    precentifNminus1 = (N - 1) * 2 / (N * (N - 1))
    if len(d) <= len(arrayofgaps):
        for i in range(len(d), len(arrayofgaps)):
            gap = gap_array[i]
            extrapolated = d[-1] - (d[-1] - precentifNminus1) * (arrayofgaps[i - 1] - gap) / (
                        arrayofgaps[i - 1] - arrayofgaps[-1])
            d.append(extrapolated)
    Z=np.zeros((len(arrayofgaps),2),dtype=float)
    Z[:,0],Z[:,1]=arrayofgaps,d
    return Z


def Create_AllData_with_interpolations(Narray,Karray,reps,data,arrayofgaps):
    D=[]
    for i in range(len(Narray)):
        N=Narray[i]
        Z=[[100000 for p in range(reps)] for q in range(len(Karray))]
        for j in range(len(Karray)):
            for m in range(reps):
                Z[j][m]=Interpolation(arrayofgaps,data[i][j][m],N)
        D.append(Z)
    return D


def Create_AllData_with_averages_Last_step(Narray,Karray,reps,data,arrayofgaps):
    D=[]
    for i in range(len(Narray)):
        Z=[]
        for j in range(len(Karray)):
            AVGMATRIX=[np.array(data[i][j][m][:,-1]) for m in range(reps)]
            AVGColumn=np.mean([AVGMATRIX[m] for m in range(len(AVGMATRIX))],axis=0).tolist()
            dummymatrixtobeappended=np.zeros((len(arrayofgaps),2),dtype=float)
            dummymatrixtobeappended[:,0],dummymatrixtobeappended[:,1]=arrayofgaps,AVGColumn
            Z.append(dummymatrixtobeappended)
        D.append(Z)
    return D




N_array=[25,50,100]
K_allowable_changes_array=np.arange(2,11,1)
nbofreplications=1000
gap_array=np.linspace(0,1,100)


# #Read files

# #RANDOMIZED rij
Data_needed_RANDOMIZED=[]
for i in range(len(N_array)):
    with open(os.path.join('data for new runs complexity','data when N={}'.format(N_array[i])+'.npy'), 'rb') as f:
        Data_needed_RANDOMIZED.append((np.load(f,allow_pickle=True, fix_imports=True)))

Data_needed_RANDOMIZED=Create_gap_from_the_given_data(Data_needed_RANDOMIZED,K_allowable_changes_array,nbofreplications)
Data_needed_RANDOMIZED=Create_AllData_with_interpolations(N_array,K_allowable_changes_array,nbofreplications,Data_needed_RANDOMIZED,gap_array)

Data_needed_RANDOMIZED=Create_AllData_with_averages_Last_step(N_array,K_allowable_changes_array,nbofreplications,Data_needed_RANDOMIZED,gap_array)

for i in range(len(N_array)):
    for j in range(len(K_allowable_changes_array)):
        input()
        print(Data_needed_RANDOMIZED[i][j])


for i in range(len(Data_needed_RANDOMIZED)):
    Dummy=np.array(Data_needed_RANDOMIZED[i][0][:,-1].reshape(len(gap_array),1))
    for j in range(1,len(K_allowable_changes_array)):
        Dummy=np.hstack((Dummy,Data_needed_RANDOMIZED[i][j][:,-1].reshape(len(gap_array),1)))
    df=pd.DataFrame(Dummy,index=gap_array,columns=['k='+str(x) for x in K_allowable_changes_array])
    df.to_csv('caseStudyNumericalExperimentsrandomized when N={}'.format(N_array[i])+'.csv')