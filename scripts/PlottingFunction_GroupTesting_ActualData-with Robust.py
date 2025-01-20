import matplotlib.pyplot as plt
import numpy as np
import matplotlib


T = np.array([ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16,17, 18, 19])
def Plotting_function(F_of_n_equal_16,F_continuous_array, F_constrained_matrix_of_ks, FconstrainedmatrixofksRobust, path_matrix_of_ks, group_size_path_matrix_of_ks, K_allowable_changes_array,minmaxdistance,averagedistance,minmaxdistancefrom16,averagedistancefrom16,totalaveragedistancefrom16robust):
    number_for_the_subplot=1
    matplotlib.rc('text', usetex=True)
    matplotlib.rcParams.update({'font.size': 28})
    for k_allowable_changes in K_allowable_changes_array:
        plt.subplot(2, 2, number_for_the_subplot)
        a = 0
        pathK = path_matrix_of_ks[k_allowable_changes-1]
        group_size_pathK = group_size_path_matrix_of_ks[k_allowable_changes-1]
        for q in range(len(pathK)):
            t = pathK[q]
            s = group_size_pathK[a]
            a = a + 1
            if a == len(group_size_pathK):
                a = a - 1
            if q < len(pathK)-1:
                plt.axvline(x=t, ls='--', c='r')
                if   k_allowable_changes==3 and q==1:
                    plt.text(t + 0.35, 0.147*costtomultiply, "$ x^{*}$=%a" % s, fontsize=24)
                elif k_allowable_changes==3 and q==2:
                    plt.text(t + 0.35, 0.143*costtomultiply, "$ x^{*}$=%a" % s, fontsize=24)
                elif k_allowable_changes==2 and q==1:
                    plt.text(t + 0.35, 0.147*costtomultiply, "$ x^{*}$=%a" % s, fontsize=24)
                elif k_allowable_changes==1 and q==1:
                    plt.text(t + 0.35, 0.147 * costtomultiply, "$ x^{*}$=%a" % s, fontsize=24)
                elif k_allowable_changes==4 and q==1:
                    plt.text(t + 0.12, 0.147*costtomultiply, "$ x^{*}$=%a" % s, fontsize=24)
                elif k_allowable_changes==4 and q==2:
                    plt.text(t + 0.35, 0.143*costtomultiply, "$ x^{*}$=%a" % s, fontsize=24)
                elif k_allowable_changes==4 and q==3:
                    plt.text(t + 0.35, 0.18*costtomultiply, "$ x^{*}$=%a" % s, fontsize=24)
                else:
                    plt.text(t + 0.35, 0.18*costtomultiply, "$ x^{*}$=%a" % s, fontsize=24)

        props = dict(boxstyle='round', facecolor='white', pad=0.2)

        # plt.text(7.5, 0.31,
        #          r'max = '+''.join(minmaxdistance[k_allowable_changes - 1]) + '\n'+ r'avg$\;$ = '+''.join(averagedistance[k_allowable_changes - 1]),
        #              fontsize=18, verticalalignment='top', bbox=props)
        plt.text(7, 0.16*costtomultiply, r'regret:$\; \; \; \;$ \$'+"{:.1f}".format(averagedistancefrom16[k_allowable_changes - 1]) + "K" + "\n"
                 + r'minimax: \$'+"{:.1f}".format(totalaveragedistancefrom16robust[k_allowable_changes - 1]) + "K",
                 fontsize=22, verticalalignment='top', bbox=props)
        plt.text(6.45, 0.165*costtomultiply, "cost deviation from CP", fontsize=22)
        plt.xticks(np.arange(21), ('Jun', ' ', ' ', ' ', 'Jul', ' ', ' ', ' ', 'Aug', ' ', ' ', ' ', 'Sep', ' ', ' ', ' ', 'Oct', ' ', ' ', ' '))
        if k_allowable_changes > 1:
            percentage_improvement_minmax = round(100*(abs(float(minmaxdistance[k_allowable_changes - 1])- float(minmaxdistance[0])))/float((minmaxdistance[0])),2)
            percentage_improvement_avgdist = round(100*(abs(float(averagedistance[k_allowable_changes - 1]) - float(averagedistance[0]))) / float((averagedistance[0])),2)
            percentage_improvement_minmaxfrom16 = round(100*(abs(float(minmaxdistancefrom16[k_allowable_changes - 1]) - float(minmaxdistancefrom16[0]))) / float((minmaxdistancefrom16[0])), 2)
            percentage_improvement_avgdistfrom16 = round(100*(abs(float(averagedistancefrom16[k_allowable_changes - 1]) - float(averagedistancefrom16[0]))) / float((averagedistancefrom16[0])), 2)
            write_title=k_allowable_changes-1
            # plt.title("$K=%g$ " % write_title + "$(%a$" %percentage_improvement_minmax + "$\%$" + ", $%a$" %percentage_improvement_avgdist + "$\%$), "
            #           +"$(%a$" %percentage_improvement_minmaxfrom16 + "$\%$" + ", $%a$" %percentage_improvement_avgdistfrom16 + "$\%$)", fontsize=24)
            # plt.title(
            #     "$K=%g$ " % write_title + "$(%a$" % percentage_improvement_minmaxfrom16 + "$\%$" + ", $%a$" % percentage_improvement_avgdistfrom16 + "$\%$)",
            #     fontsize=24)
            plt.title("$K=%g$ " % write_title,fontsize=28)
        else:
            write_title = k_allowable_changes - 1
            plt.title("$K=%g$ " % write_title, fontsize=28)
        plt.plot(T, F_continuous_array, linewidth=0.8, label="ideal: $f^{*}_{t}$", color='k')
        plt.plot(T, F_constrained_matrix_of_ks[k_allowable_changes-1], ls='-.', color='k', linewidth=0.8, label="restricted: $\hat{f}_{t}$")
        plt.plot(T, FconstrainedmatrixofksRobust[k_allowable_changes-1],ls=":",color='k', linewidth=0.8, label="restricted: robust")
        plt.plot(T, F_of_n_equal_16, linewidth=0.8, label="current practice: $f_{CP_{t}}$", color='k',ls='--',dashes=(4,12))
        # plt.xlabel('time (weeks)')
        if k_allowable_changes in [1,3]:
            plt.ylabel(r'$\rm{Cost}$ ' + r'$(\rm{in\;} \$1\rm{K})$')
        # plt.legend(loc='upper right', fontsize=26, prop={"size":16})
        number_for_the_subplot = number_for_the_subplot + 1
    plt.show()
    return


costtomultiply=1000*13.6/len(T) #(nb of donors)
# F_of_n_star_continuous =costtomultiply*np.array([0.2558529265761047, 0.26941858214387615, 0.2884826508538203, 0.3146251613257274, 0.34930326475335916, 0.4003555937177521, 0.4490757498314597, 0.500568248033797, 0.5401095267842027, 0.5539971920096751, 0.5520548269853158, 0.527405724534145, 0.48584371332052934, 0.4366764607365723, 0.3885695349565914, 0.3427073236064069, 0.3105091263015938, 0.286223600651379, 0.26909010564213576, 0.25888608019911263])
F_of_n_star_continuous =costtomultiply*np.array([0.1396265909007297, 0.14401507089769794, 0.14933687336662094, 0.15561923045089499, 0.16280907514245047, 0.17191191422500296, 0.17952707206421215, 0.18679546270672254, 0.19195553645660968, 0.19369749438229378, 0.19345748741339897, 0.1903312656232694, 0.18478591692746726, 0.17766637980439337, 0.1699320382072531, 0.16151801197486138, 0.15468970809523164, 0.1487452793286519, 0.14391568414819145, 0.14066518089603686])


# F_of_n_equal_16=costtomultiply*np.array([0.30038527347502253, 0.3091952820269994, 0.3228771151770151, 0.34319949168452357, 0.37191757423512595, 0.41646585265225644, 0.460633555006952, 0.508506733146588, 0.5459004068638063, 0.559123817487762, 0.557268800106494, 0.5338223373311939, 0.49470842617075617, 0.4492748010451266, 0.4060016202482837, 0.36633596697942883, 0.3399100028533706, 0.3211949686072062, 0.30897098645268795, 0.30226691118503046])
F_of_n_equal_16=costtomultiply*np.array([0.18844870686390114, 0.18909721478722097, 0.1901041004216748, 0.19159914610251239, 0.19371074083062467, 0.1969837716401417, 0.2002258019816392, 0.20373640950365934, 0.2064760708576696, 0.20744437082653366, 0.2073085512924957, 0.20559140333399584, 0.20272492447264145, 0.19939232723618083, 0.1962152250212681, 0.19330043441079336, 0.19135719350242786, 0.18998032207455817, 0.18908070580482628, 0.1885872250266829])


# F_constrained_matrix_of_ks =costtomultiply*np.array([[0.27473972285977233, 0.2843286569809256, 0.2992191299716427, 0.321334472424339, 0.35258157372442334, 0.40104213296552205, 0.4490757498314597, 0.5011246529459208, 0.5417694932498147, 0.5561404061917186, 0.5541244831352099, 0.528642308894673, 0.4861243400507482, 0.4367240153411447, 0.3896600989225174, 0.34650885760620787, 0.31775495033549905, 0.29738845121217583, 0.2840845369980025, 0.27678776645015857], [0.26286234159218935, 0.2735537457173476, 0.29015387735858633, 0.3148031444925321, 0.35258157372442334, 0.40104213296552205, 0.4490757498314597, 0.5011246529459208, 0.5417694932498147, 0.5561404061917186, 0.5541244831352099, 0.528642308894673, 0.4861243400507482, 0.4367240153411447, 0.3896600989225174, 0.34650885760620787, 0.31775495033549905, 0.29738845121217583, 0.2840845369980025, 0.27678776645015857], [0.26286234159218935, 0.2735537457173476, 0.29015387735858633, 0.3148031444925321, 0.3560522443409089, 0.40321074188207984, 0.44995856588458105, 0.5006195742443911, 0.5401844968050296, 0.5541743943599293, 0.5522118937504175, 0.527405724534145, 0.48601867426754364, 0.43793701437800436, 0.3921340886619151, 0.34280379093664903, 0.3119218051532293, 0.2900441945105525, 0.275751486759086, 0.2679118246824417], [0.26286234159218935, 0.2735537457173476, 0.29015387735858633, 0.3148031444925321, 0.3560522443409089, 0.40321074188207984, 0.44995856588458105, 0.5006195742443911, 0.5401844968050296, 0.5541743943599293, 0.5522118937504175, 0.527405724534145, 0.48601867426754364, 0.43793701437800436, 0.38859260989091204, 0.3439934828961069, 0.31427157195485744, 0.2865511585637921, 0.2708221574276215, 0.26219353879736396]])
F_constrained_matrix_of_ks =costtomultiply*np.array([[0.15342564, 0.15486725, 0.15710451, 0.16042421, 0.16510837,
        0.17235838, 0.17952707, 0.18727544, 0.19331199, 0.19544339,
        0.19514449, 0.1913637 , 0.18504447, 0.17768532, 0.17065715,
        0.16419861, 0.15988715, 0.15682955, 0.15483055, 0.1537336 ],
       [0.14471657, 0.14698615, 0.15050672, 0.15572689, 0.16510837,
        0.17235838, 0.17952707, 0.18727544, 0.19331199, 0.19544339,
        0.19514449, 0.1913637 , 0.18504447, 0.17768532, 0.17065715,
        0.16419861, 0.15988715, 0.15682955, 0.15483055, 0.1537336 ],
       [0.14471657, 0.14698615, 0.15050672, 0.15572689, 0.16737472,
        0.17372386, 0.18000432, 0.18679546, 0.19208831, 0.19395755,
        0.1936954 , 0.19037985, 0.1848398 , 0.17839053, 0.17223379,
        0.16156682, 0.15564781, 0.15144759, 0.14870041, 0.1471925 ],
       [0.14471657, 0.14698615, 0.15050672, 0.15572689, 0.16737472,
        0.17372386, 0.18000432, 0.18679546, 0.19208831, 0.19395755,
        0.1936954 , 0.19037985, 0.1848398 , 0.17839053, 0.16995792,
        0.16253777, 0.15758259, 0.14898225, 0.14518951, 0.14310695]])

F_constrained_matrix_of_ks_Robust=costtomultiply*np.array([[0.1738317 , 0.17466056, 0.17594732, 0.17785765, 0.18055519,
        0.18473507, 0.18887372, 0.19335337, 0.19684793, 0.19808276,
        0.19790956, 0.19571963, 0.19206288, 0.1878099 , 0.18375374,
        0.18003108, 0.17754851, 0.17578915, 0.17463946, 0.17400875],
       [0.15408411, 0.15548969, 0.15767109, 0.160908  , 0.18055519,
        0.18473507, 0.18887372, 0.19335337, 0.19684793, 0.19808276,
        0.19790956, 0.19571963, 0.19206288, 0.1878099 , 0.18375374,
        0.18003108, 0.17754851, 0.17578915, 0.17463946, 0.17400875],
       [0.15408411, 0.15548969, 0.15767109, 0.160908  , 0.18055519,
        0.18473507, 0.18887372, 0.19335337, 0.19684793, 0.19808276,
        0.19790956, 0.19571963, 0.19206288, 0.1878099 , 0.18375374,
        0.18003108, 0.16038432, 0.15740299, 0.15545392, 0.15438437],
       [0.15408411, 0.15548969, 0.15767109, 0.160908  , 0.18055519,
        0.18473507, 0.18887372, 0.19335337, 0.19684793, 0.19808276,
        0.19790956, 0.19571963, 0.19206288, 0.1878099 , 0.17607164,
        0.17138146, 0.16825255, 0.15272457, 0.15022666, 0.1488557 ]])





path_matrix_of_ks =[np.array([ 0, 19]), np.array([ 0,  4, 19]), np.array([ 0,  4, 15, 19]), np.array([ 0,  4,  14, 17, 19])]


group_size_path_of_ks =[np.array([28]), np.array([61, 38]), np.array([61, 33, 53]), np.array([61, 33, 44, 74])]

K_allowable_changes_array =[1, 2, 3, 4]



MinMaxDistance =(costtomultiply*np.array([0.01679, 0.01577, 0.00916, 0.00903,])).tolist()
Total_Average_distance=costtomultiply*np.array([0.004360000001,  0.002970000001, 0.001750000001, 0.001280000001])
Total_Average_distance=[i*len(F_of_n_star_continuous) for i in Total_Average_distance] # we are doing this bcz those values they are averages, i.e., divided over the time horizon
MinMaxDistance_from16=(costtomultiply*np.array([0.03502, 0.04373, 0.04373, 0.04548,])).tolist()

Total_Average_distance_from16=(costtomultiply*np.array([0.02516, 0.02655, 0.02777, 0.02824])).tolist()
Total_Average_distance_from16=[i*len(F_of_n_star_continuous) for i in Total_Average_distance_from16]


Total_Average_distance_from16_Robust=(costtomultiply*np.array([0.01237, 0.01607, 0.01979, 0.02099])).tolist()
Total_Average_distance_from16_Robust=[i*len(F_of_n_star_continuous) for i in Total_Average_distance_from16_Robust]

# for m in range(len(MinMaxDistance)):
#     MinMaxDistance[m]=np.format_float_scientific(MinMaxDistance[m], exp_digits=2,precision=3)
#     Total_Average_distance[m]=np.format_float_scientific(Total_Average_distance[m], exp_digits=2, precision=3)
#     MinMaxDistance_from16[m]=np.format_float_scientific(MinMaxDistance_from16[m], exp_digits=2, precision=3)
#     Total_Average_distance_from16[m]=np.format_float_scientific(Total_Average_distance_from16[m], exp_digits=2, precision=3)

Plotting_function(F_of_n_equal_16,F_of_n_star_continuous, F_constrained_matrix_of_ks, F_constrained_matrix_of_ks_Robust, path_matrix_of_ks, group_size_path_of_ks, K_allowable_changes_array,MinMaxDistance,Total_Average_distance,MinMaxDistance_from16,Total_Average_distance_from16, Total_Average_distance_from16_Robust)