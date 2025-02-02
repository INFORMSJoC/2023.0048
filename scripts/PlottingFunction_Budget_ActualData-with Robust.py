import matplotlib.pyplot as plt
import numpy as np
import matplotlib
# from pylab import *
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]

def Plotting_function(F_continuous_array, F_constrained_matrix_of_ks,F_constrained_matrix_of_ks_robust, path_matrix_of_ks, budget_path_3d_matrix_of_ks, K_allowable_changes_array,minmaxdistance,averagedistance, averagedistanceRob):
    number_for_the_subplot=1
    matplotlib.rc('text', usetex=True)
    matplotlib.rcParams.update({'font.size': 28})
    for k_allowable_changes in K_allowable_changes_array:
        plt.subplot(2, 2, number_for_the_subplot)
        a=0
        pathK = path_matrix_of_ks[k_allowable_changes-1]
        budget_pathK= budget_path_3d_matrix_of_ks[k_allowable_changes-1]
        for q in range(len(pathK)):
            t = pathK[q]
            s=budget_pathK[a]
            s = np.array(s)
            a = a + 1
            if a == len(budget_pathK):
                a = a - 1
            if q < len(pathK)-1:
                plt.axvline(x=t, ls='--', c='r')
                if k_allowable_changes==1:
                    plt.text(t + 0.5, 0.0024*avgcostofTTIs,r"$\boldsymbol{x}^{*}=$ " "\n" r'[$ \$ $%a' %s[0] + r"" "\n" r" $ \$ $%a" % s[1] + r"" "\n" r" $ \$ $%a" % s[2] + r"" "\n" r" $ \$ $%a" % s[3] + r"" "\n" r" $ \$ $%a]" % s[4],fontsize=20)
                elif k_allowable_changes==4:
                    if q==1:
                        plt.text(t+0.5, 0.0019*avgcostofTTIs,r"$\boldsymbol{x}^{*}=$ " "\n" r'[$ \$ $%a' %s[0] + r" " "\n" r'$ \$ $%a' %s[1] + r" " "\n" r'$ \$ $%a'%s[2] + r" " "\n" r'$ \$ $%a'%s[3] + r" " "\n" r'$ \$ $%a]'%s[4],fontsize=20)
                    elif q==len(pathK)-2:
                        plt.text(t + 1.5, 0.00251*avgcostofTTIs,r"$\boldsymbol{x}^{*}=$ " "\n" r'[$ \$ $%a' % s[0] + r" " "\n" r'$ \$ $%a' % s[1] + r" " "\n" r'$ \$ $%a' % s[2] + r" " "\n" r'$ \$ $%a' % s[3] + r" " "\n" r'$ \$ $%a]' % s[4],fontsize=20)
                    elif q==len(pathK)-1:
                        plt.text(t + 0.5, 0.0019*avgcostofTTIs,r"$\boldsymbol{x}^{*}=$ " "\n" r'[$ \$ $%a' % s[0] + r" " "\n" r'$ \$ $%a' % s[1] + r" " "\n" r'$ \$ $%a' % s[2] + r" " "\n" r'$ \$ $%a' % s[3] + r" " "\n" r'$ \$ $%a]' % s[4],fontsize=20)
                    elif q==2:
                        plt.text(t + 0.5, 0.0019*avgcostofTTIs,r"$\boldsymbol{x}^{*}=$ " "\n" r'[$ \$ $%a' % s[0] + r" " "\n" r'$ \$ $%a' % s[1] + r" " "\n" r'$ \$ $%a' % s[2] + r" " "\n" r'$ \$ $%a' % s[3] + r" " "\n" r'$ \$ $%a]' % s[4],fontsize=20)
                    else:
                        plt.text(t + 0.5, 0.0021*avgcostofTTIs, r"$\boldsymbol{x}^{*}=$ " "\n" r'[$ \$ $%a' % s[0] + r" " "\n" r'$ \$ $%a' % s[1] + r" " "\n" r'$ \$ $%a' % s[2] + r" " "\n" r'$ \$ $%a' % s[3] + r" " "\n" r'$ \$ $%a]' % s[4],fontsize=20)
                elif k_allowable_changes==2:
                    if q==0:
                        plt.text(t+0.5, 0.0024*avgcostofTTIs,r"$\boldsymbol{x}^{*}=$ " "\n" r'[$ \$ $%a' %s[0] + r" " "\n" r'$ \$ $%a'%s[1] + r" " "\n" r'$ \$ $%a'%s[2] + r" " "\n" r'$ \$ $%a'%s[3] + r" " "\n" r'$ \$ $%a]'%s[4],fontsize=20)
                    else:
                        plt.text(t + 1.5, 0.0025*avgcostofTTIs, r"$\boldsymbol{x}^{*}=$ " "\n" r'[$ \$ $%a' % s[0] + r" " "\n" r'$ \$ $%a' % s[1] + r" " "\n" r'$ \$ $%a' % s[2] + r" " "\n" r'$ \$ $%a' % s[3] + r" " "\n" r'$ \$ $%a]' % s[4],fontsize=20)
                elif k_allowable_changes==3:
                    if q==0:
                        plt.text(t+0.5, 0.0021*avgcostofTTIs,r"$\boldsymbol{x}^{*}=$ " "\n" r'[$ \$ $%a' %s[0] + r" " "\n" r'$ \$ $%a'%s[1] + r" " "\n" r'$ \$ $%a'%s[2] + r" " "\n" r'$ \$ $%a'%s[3] + r" " "\n" r'$ \$ $%a]'%s[4],fontsize=20)
                    elif q==1:
                        plt.text(t + 0.5, 0.002*avgcostofTTIs,r"$\boldsymbol{x}^{*}=$ " "\n" r'[$ \$ $%a' % s[0] + r" " "\n" r'$ \$ $%a' % s[1] + r" " "\n" r'$ \$ $%a' % s[2] + r" " "\n" r'$ \$ $%a' % s[3] + r" " "\n" r'$ \$ $%a]' % s[4], fontsize=20)
                    elif q==2:
                        plt.text(t + 0.5, 0.00185*avgcostofTTIs,r"$\boldsymbol{x}^{*}=$ " "\n" r'[$ \$ $%a' % s[0] + r" " "\n" r'$ \$ $%a' % s[1] + r" " "\n" r'$ \$ $%a' % s[2] + r" " "\n" r'$ \$ $%a' % s[3] + r" " "\n" r'$ \$ $%a]' % s[4], fontsize=20)
        if k_allowable_changes==1:
            props = dict(boxstyle='round', facecolor='white', pad=0.23)
            # plt.text(22, 0.0022*1000,r'max = '+''.join(minmaxdistance[k_allowable_changes - 1]) + r" " "\n" +r'avg$\;$ = '+''.join(averagedistance[k_allowable_changes - 1]),fontsize=22, verticalalignment='top', bbox=props)
            plt.text(21, 0.34,r" " r'regret:$\;\;\;\;\;$' "\$"+"{:.3f}".format(averagedistance[k_allowable_changes - 1]) + "M" + "\n"
                     + r" " r'minimax: ' "\$"+"{:.3f}".format(averagedistanceRob[k_allowable_changes - 1]) + "M", fontsize=20, verticalalignment='top', bbox=props)
            plt.text(19, 0.35, "cost deviation from ideal", fontsize=22)
        elif k_allowable_changes==2:
            props = dict(boxstyle='round', facecolor='white', pad=0.23)
            # plt.text(20, 0.0022*1000,r'max = '+''.join(minmaxdistance[k_allowable_changes - 1]) + r" " "\n" +r'avg$\;$ = '+''.join(averagedistance[k_allowable_changes - 1]),fontsize=22, verticalalignment='top', bbox=props)
            plt.text(18, 0.33,r" " r'regret:$\;\;\;\;\;$' "\$"+"{:.3f}".format(averagedistance[k_allowable_changes - 1]) + "M"+ "\n"
                     + r" " r'minimax: ' "\$"+"{:.3f}".format(averagedistanceRob[k_allowable_changes - 1]) + "M", fontsize=20, verticalalignment='top', bbox=props)
            plt.text(16.3, 0.34, "cost deviation from ideal", fontsize=22)
        elif k_allowable_changes==3:
            props = dict(boxstyle='round', facecolor='white', pad=0.23)
            # plt.text(6, 0.0027*1000,r'max = '+''.join(minmaxdistance[k_allowable_changes - 1]) + r" " "\n" +r'avg$\;$ = '+''.join(averagedistance[k_allowable_changes - 1]),fontsize=22, verticalalignment='top', bbox=props)
            plt.text(4, 0.47,r" " r'regret:$\;\;\;\;\;$' "\$"+"{:.3f}".format(averagedistance[k_allowable_changes - 1]) + "M"+ "\n"
                     + r" " r'minimax: ' "\$"+"{:.3f}".format(averagedistanceRob[k_allowable_changes - 1]) + "M", fontsize=20, verticalalignment='top', bbox=props)
            plt.text(2, 0.48, "cost deviation from ideal", fontsize=22)
        else:
            props = dict(boxstyle='round', facecolor='white',pad=0.23)
            # plt.text(6, 0.0027*1000, r'max = '+''.join(minmaxdistance[k_allowable_changes-1])+r" " "\n" +r'avg$\;$ = '+''.join(averagedistance[k_allowable_changes-1]),fontsize=22,verticalalignment='top', bbox=props)
            plt.text(4, 0.47,r" " r'regret:$\;\;\;\;\;$' "\$"+"{:.3f}".format(averagedistance[k_allowable_changes - 1]) + "M"+ "\n"
                     + r" " r'minimax: ' "\$"+"{:.3f}".format(averagedistanceRob[k_allowable_changes - 1]) + "M", fontsize=20, verticalalignment='top', bbox=props)
            plt.text(1.2, 0.48, "cost deviation from ideal", fontsize=22)
        plt.xticks(np.arange(53), (' ', 'Jan', ' ', ' ', ' ', 'Feb', ' ', ' ', ' ', 'Mar', ' ', ' ', ' ', 'Apr', ' ', ' ', ' ',' ', 'May', ' ', ' ',' ', ' Jun', ' ', ' ', ' ', 'Jul', ' ', ' ',' ' ,'Aug', ' ', ' ', ' ', ' ', 'Sep ', ' ', ' ', ' ', 'Oct ', ' ', ' ', ' ' ,' ', 'Nov', ' ', ' ', ' ', 'Dec', ' ', ' ', ' ', ' '))
        if k_allowable_changes > 1:
            percentage_improvement_minmax = 100*(abs(float(minmaxdistance[k_allowable_changes - 1]) - float(minmaxdistance[0]))) / (float(minmaxdistance[0]))
            percentage_improvement_avgdist = 100*(abs(float(averagedistance[k_allowable_changes - 1]) - float(averagedistance[0]))) / (float(averagedistance[0]))
            # percentage_improvement_minmax_from_current_practice = 100*(abs(float(minmaxdistance_from_current_practice[k_allowable_changes - 1]) - float(minmaxdistance_from_current_practice[0]))) / (float(minmaxdistance_from_current_practice[0]))
            # percentage_improvement_avgdist_from_current_practice = 100*(abs(float(averagedistance_from_current_practice[k_allowable_changes - 1]) - float(averagedistance_from_current_practice[0]))) / (float(averagedistance_from_current_practice[0]))
            write_title=k_allowable_changes-1
            # plt.title("$K=%g$ " % write_title + "$(%a$" %round(percentage_improvement_minmax,2) + "$\%$" + ", $%a$" %round(percentage_improvement_avgdist,2) + "$\%$)", fontsize=26)
            plt.title("$K=%g$ " % write_title + "$(%a$" %round(percentage_improvement_avgdist,2) + "$\%$)", fontsize=26)

            # plt.title("$K=%g$ " % write_title + "$(%a$" % round(percentage_improvement_minmax, 2) + "$\%$" + ", $%a$" % round(percentage_improvement_avgdist, 2) + "$\%$)"
            #           + "$(%a$" % round(percentage_improvement_minmax_from_current_practice, 2) + "$\%$" + ", $%a$" % round(percentage_improvement_avgdist_from_current_practice, 2) + "$\%$)", fontsize=19)
        else:
            write_title = k_allowable_changes - 1
            plt.title("$K=%g$" % write_title, fontsize=26)
        plt.plot(T, F_continuous_array, color='k', linewidth=0.8, label="$f^{*}_{t}$ (ideal)")
        plt.plot(T, F_constrained_matrix_of_ks[k_allowable_changes-1], color='k',ls='-.', linewidth=0.8, label="restricted $\hat{f}_{t}$ (regret)")
        plt.plot(T, F_constrained_matrix_of_ks_robust[k_allowable_changes-1],color='k',ls=':', linewidth=0.8, label="restricted $\hat{f}_{t}$ (minimax)")
        # plt.plot(T, function_current_practice, color='r', ls=':', linewidth=0.8,label="$f_{CP}(t)$")
        if k_allowable_changes in [1,3]:
            # plt.legend(loc='upper right', fontsize=20)
            plt.ylabel(r'$\rm{Cost}$ ' + r'$(\rm{in\;} \$1\rm{M})$')
        # else:
            # plt.legend(loc='lower right', fontsize=20)

        # plt.yscale("log")
        number_for_the_subplot = number_for_the_subplot + 1
    plt.show()
    return


Prevalence_Data_List =[[0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007], [0.00345, 0.00345, 0.00345, 0.00345, 0.00345, 0.00345, 0.00345, 0.00345, 0.00345, 0.00345, 0.00345, 0.00345, 0.00345, 0.00345, 0.00345, 0.00345, 0.00345, 0.00345, 0.00345, 0.00345, 0.00345, 0.00345, 0.00345, 0.00345, 0.00345, 0.00345, 0.00345, 0.00345, 0.00345, 0.00345, 0.00345, 0.00345, 0.00345, 0.00345, 0.00345, 0.00345, 0.00345, 0.00345, 0.00345, 0.00345, 0.00345, 0.00345, 0.00345, 0.00345, 0.00345, 0.00345, 0.00345, 0.00345, 0.00345, 0.00345, 0.00345, 0.00345, 0.00345], [0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016], [0.000655144, 0.000590234, 0.000525325, 0.000460416, 0.000498371, 0.000536327, 0.000574282, 0.000612237, 0.000639879, 0.00066752, 0.000695162, 0.000722803, 0.000783532, 0.000844261, 0.000904989, 0.000965718, 0.001026447, 0.001258717, 0.001490988, 0.001723259, 0.001955529, 0.002429147, 0.002902765, 0.003376382, 0.00385, 0.003691907, 0.003533815, 0.003375722, 0.00321763, 0.003059537, 0.002693597, 0.002327658, 0.001961718, 0.001595778, 0.001512854, 0.001429929, 0.001347005, 0.001264081, 0.001220184, 0.001176288, 0.001132392, 0.001088495, 0.001044599, 0.000973227, 0.000901854, 0.000830481, 0.000759108, 0.000687736, 0.000616363, 0.00054499, 0.000473618, 0.000402245, 0.000330872], [3.34073e-07, 2.7044e-07, 2.06807e-07, 1.43174e-07, 2.50554e-07, 3.57935e-07, 4.65315e-07, 5.72696e-07, 6.68145e-07, 7.63594e-07, 8.59044e-07, 9.54493e-07, 1.4556e-06, 1.95671e-06, 2.45782e-06, 2.95893e-06, 3.46004e-06, 6.84252e-06, 1.0225e-05, 1.36075e-05, 1.699e-05, 5.57185e-05, 9.44471e-05, 0.000133176, 0.000171904, 0.000236523, 0.000301142, 0.000365762, 0.000430381, 0.000495, 0.000470475, 0.000445951, 0.000421426, 0.000396902, 0.00032473, 0.000252559, 0.000180387, 0.000108216, 8.92594e-05, 7.03032e-05, 5.13469e-05, 3.23907e-05, 1.34345e-05, 1.17403e-05, 1.0046e-05, 8.35181e-06, 6.65759e-06, 4.96336e-06, 3.26914e-06, 1.57491e-06, 1.57e-06, 1.57e-06, 1.57e-06]]

# avgcostofTTIs is according to the CDC data: 16000000000/26000000=615 (we have to multiply with nb of donors in the US)
avgcostofTTIs=16000000000*(13.6/len(Prevalence_Data_List[0]))/26000000

F_continuous_array = avgcostofTTIs*np.array([0.0020533000221661323, 0.0020246674772073474, 0.001993170122938167, 0.0019581183800808227, 0.001979171029192673, 0.001998882212800542, 0.0020174311481751623, 0.0020349592258259138, 0.00204716499520411, 0.002058927174143873, 0.002070267072445489, 0.0020812029558433002, 0.002104377191461207, 0.002126059216606352, 0.0021465147762097718, 0.002165787445177857, 0.002184074300001842, 0.002247986766690335, 0.0023029464613818005, 0.002351352271055341, 0.0023947552638835592, 0.0025037522773468544, 0.002601792478903882, 0.002691859443127432, 0.0027760879865071104, 0.0028261168686561744, 0.002875559183415577, 0.002924365330717299, 0.0029725771471517336, 0.0030200842355579665, 0.0029527376227678854, 0.0028799966377535124, 0.002800173494981712, 0.0027106197728262557, 0.0026218982029495358, 0.0025324583110018796, 0.0024420368942256803, 0.002350712531673501, 0.0023211074173664246, 0.0022911707128710156, 0.0022609741157507994, 0.002230254065820021, 0.0021991842765367573, 0.0021768094791286647, 0.0021530617561600278, 0.0021277858311334603, 0.0021007107388391413, 0.0020714630523903576, 0.0020395295798364987, 0.002004427004007437, 0.001967030117316828, 0.0019244037666106349, 0.001874612889725994])

function_current_practice=avgcostofTTIs*np.array([0.0034151 , 0.00335018, 0.00328526, 0.00322034, 0.00325831,
       0.00329629, 0.00333426, 0.00337223, 0.00339989, 0.00342754,
       0.0034552 , 0.00348286, 0.00354367, 0.00360447, 0.00366528,
       0.00372609, 0.0037869 , 0.0040197 , 0.0042525 , 0.0044853 ,
       0.0047181 , 0.00519781, 0.00567752, 0.00615723, 0.00663693,
       0.006489  , 0.00634107, 0.00619314, 0.00604521, 0.00589727,
       0.00552748, 0.00515768, 0.00478789, 0.00441809, 0.00432382,
       0.00422954, 0.00413527, 0.004041  , 0.00399412, 0.00394725,
       0.00390037, 0.00385349, 0.00380662, 0.00373498, 0.00366334,
       0.0035917 , 0.00352006, 0.00344842, 0.00337678, 0.00330514,
       0.00323377, 0.0031624 , 0.00309102])

F_constrained_matrix_of_ks = avgcostofTTIs*np.array([np.array([0.0022595 , 0.00224577, 0.00223205, 0.00221832, 0.0022264 ,
       0.00223449, 0.00224257, 0.00225065, 0.00225656, 0.00226246,
       0.00226836, 0.00227426, 0.00228746, 0.00230066, 0.00231386,
       0.00232706, 0.00234026, 0.00239193, 0.0024436 , 0.00249527,
       0.00254693, 0.00267789, 0.00280885, 0.00293981, 0.00307076,
       0.00308942, 0.00310809, 0.00312675, 0.00314541, 0.00316408,
       0.00306723, 0.00297039, 0.00287355, 0.00277671, 0.00270118,
       0.00262566, 0.00255013, 0.0024746 , 0.0024501 , 0.00242561,
       0.00240111, 0.00237661, 0.00235211, 0.00233571, 0.00231931,
       0.0023029 , 0.0022865 , 0.0022701 , 0.00225369, 0.00223729,
       0.00222225, 0.00220721, 0.00219216]), np.array([0.00227041, 0.00225754, 0.00224467, 0.00223181, 0.00223939,
       0.00224697, 0.00225455, 0.00226213, 0.00226766, 0.0022732 ,
       0.00227873, 0.00228427, 0.00229666, 0.00230906, 0.00232145,
       0.00233385, 0.00234624, 0.00239483, 0.00244342, 0.002492  ,
       0.00254059, 0.00266534, 0.00279009, 0.00291484, 0.00303959,
       0.00306049, 0.0030814 , 0.00310231, 0.00312322, 0.00314412,
       0.00305209, 0.00296006, 0.00286803, 0.002776  , 0.00270141,
       0.00262682, 0.00249345, 0.0023904 , 0.0023551 , 0.0023198 ,
       0.00228449, 0.00224919, 0.00221389, 0.00218562, 0.00215734,
       0.00212907, 0.0021008 , 0.00207253, 0.00204426, 0.00201598,
       0.0019894 , 0.00196282, 0.00193624]), np.array([0.00212341, 0.00211011, 0.00209681, 0.00208351, 0.00209136,
       0.0020992 , 0.00210705, 0.0021149 , 0.00212063, 0.00212636,
       0.00213209, 0.00213782, 0.00215071, 0.00216359, 0.00217648,
       0.00218936, 0.00220224, 0.00225299, 0.00230373, 0.00235447,
       0.00240522, 0.00254052, 0.00267582, 0.00282394, 0.00292212,
       0.00295335, 0.00298458, 0.0030158 , 0.00304703, 0.00307826,
       0.00300686, 0.00293546, 0.00286406, 0.00279266, 0.00271809,
       0.00262263, 0.00252718, 0.00243172, 0.00239867, 0.00236563,
       0.00233258, 0.00229953, 0.00226648, 0.0022392 , 0.00221193,
       0.00218465, 0.00215737, 0.00213009, 0.00210281, 0.00207554,
       0.00204979, 0.00202405, 0.00199831]), np.array([0.00209297, 0.00207684, 0.00206071, 0.00204458, 0.00205409,
       0.00206359, 0.00207309, 0.00208259, 0.00208952, 0.00209646,
       0.0021034 , 0.00211033, 0.00212586, 0.00214139, 0.00215692,
       0.00217245, 0.00218798, 0.00224885, 0.00230971, 0.00237058,
       0.00243144, 0.0026283 , 0.00271687, 0.00280545, 0.00289402,
       0.00292897, 0.00296392, 0.00299888, 0.00303383, 0.00306878,
       0.00302175, 0.00293854, 0.00285532, 0.00277211, 0.00269812,
       0.00260404, 0.00250348, 0.00240292, 0.00236836, 0.0023338 ,
       0.00229924, 0.00226468, 0.00223012, 0.00220219, 0.00217426,
       0.00214633, 0.0021184 , 0.00209047, 0.00206254, 0.00203461,
       0.00200832, 0.00198203, 0.00195574])])



F_constrained_matrix_of_ks_Robust=avgcostofTTIs*np.array([[0.0025847 , 0.00257811, 0.00257152, 0.00256493, 0.00256883,
        0.00257273, 0.00257663, 0.00258053, 0.00258338, 0.00258623,
        0.00258909, 0.00259194, 0.0025984 , 0.00260486, 0.00261132,
        0.00261778, 0.00262424, 0.00264995, 0.00267566, 0.00270136,
        0.00272707, 0.00280101, 0.00287496, 0.00294891, 0.00302286,
        0.00305063, 0.00307841, 0.00310618, 0.00313396, 0.00316173,
        0.00310825, 0.00305477, 0.00300128, 0.0029478 , 0.00289062,
        0.00283343, 0.00277625, 0.00271907, 0.00270182, 0.00268457,
        0.00266732, 0.00265007, 0.00263282, 0.00262448, 0.00261614,
        0.0026078 , 0.00259946, 0.00259111, 0.00258277, 0.00257443,
        0.00256723, 0.00256004, 0.00255284],
       [0.00250773, 0.00250133, 0.00249493, 0.00248854, 0.00249233,
        0.00249612, 0.00249992, 0.00250371, 0.00250649, 0.00250926,
        0.00251204, 0.00251481, 0.00252113, 0.00252745, 0.00253377,
        0.00254009, 0.00254641, 0.00257168, 0.00259695, 0.00262222,
        0.0026475 , 0.00272306, 0.00279862, 0.00287418, 0.00294974,
        0.00298306, 0.00301637, 0.00304969, 0.00313789, 0.00316207,
        0.00310247, 0.00304287, 0.00298327, 0.00292368, 0.00286595,
        0.00280823, 0.0027505 , 0.00269277, 0.00267499, 0.00265721,
        0.00263943, 0.00262166, 0.00260388, 0.0025943 , 0.00258473,
        0.00257516, 0.00256558, 0.00255601, 0.00254643, 0.00253686,
        0.00252841, 0.00251996, 0.00251151],
       [0.00250773, 0.00250133, 0.00249493, 0.00248854, 0.00249233,
        0.00249612, 0.00249992, 0.00250371, 0.00250649, 0.00250926,
        0.00251204, 0.00251481, 0.00252113, 0.00252745, 0.00253377,
        0.00254009, 0.00254641, 0.00257168, 0.00259695, 0.00262222,
        0.0026475 , 0.00272306, 0.00279862, 0.00287418, 0.00294974,
        0.00298306, 0.00301637, 0.00304969, 0.00313789, 0.00316207,
        0.00310247, 0.00304287, 0.00298327, 0.00292368, 0.00286595,
        0.00253706, 0.0024472 , 0.00235734, 0.00232901, 0.00230069,
        0.00227237, 0.00224405, 0.00221573, 0.00219881, 0.00218188,
        0.00216496, 0.00214804, 0.00213112, 0.0021142 , 0.00209728,
        0.00208204, 0.00206682, 0.00205159],
       [0.0024211 , 0.00241491, 0.00240873, 0.00240254, 0.00240622,
        0.0024099 , 0.00241357, 0.00241725, 0.00241994, 0.00242263,
        0.00242533, 0.00242802, 0.00243418, 0.00244035, 0.00244651,
        0.00245268, 0.00245884, 0.00248368, 0.00250851, 0.00253334,
        0.00255818, 0.00263612, 0.00271405, 0.00279199, 0.00286993,
        0.00291041, 0.00295089, 0.00305353, 0.00313789, 0.00316207,
        0.00310247, 0.00304287, 0.00298327, 0.00292368, 0.00286595,
        0.00253706, 0.0024472 , 0.00235734, 0.00232901, 0.00230069,
        0.00227237, 0.00224405, 0.00221573, 0.00219881, 0.00218188,
        0.00216496, 0.00214804, 0.00213112, 0.0021142 , 0.00209728,
        0.00208204, 0.00206682, 0.00205159]])


path_matrix_of_ks=[np.array([ 0, 52]), np.array([ 0, 36, 52]), np.array([ 0, 23, 34, 52]), np.array([ 0, 21, 30, 35, 52])]

budget_path_3d_matrix_of_ks= [np.array([[10. ,  9.2, 20.6,  4.1,  1.2]]), np.array([[ 9.9,  9.1, 20.5,  4.3,  1.2],
       [10.5, 10.2, 21.7,  2.6,  0. ]]), np.array([[10.2,  9.5, 21.1,  4.2,  0. ],
       [ 9.8,  8.9, 20.2,  5.2,  1. ],
       [10.4,  9.9, 21.5,  2.7,  0.5]]), np.array([[10.3,  9.7, 21.3,  3.7,  0. ],
       [ 9.7,  8.7, 20. ,  5.6,  1. ],
       [ 9.9,  9.1, 20.4,  4.6,  1. ],
       [10.5, 10.1, 21.6,  2.6,  0.2]])]

BT=45
budget_path_3d_matrix_of_ksnew=[]
for i in range(len(budget_path_3d_matrix_of_ks)):
    C=[]
    for j in range(len(budget_path_3d_matrix_of_ks[i])):
        D=[]
        for m in range(len(budget_path_3d_matrix_of_ks[i][j])):
            if m!=2:
                D.append(budget_path_3d_matrix_of_ks[i][j][m])
            else:
                D.append(round(BT-sum(budget_path_3d_matrix_of_ks[i][j])+budget_path_3d_matrix_of_ks[i][j][m],2))
        C.append(D)
    budget_path_3d_matrix_of_ksnew.append(C)


MinMaxDistance=[0.0005495678856108507, 0.0005027368502913096, 0.00031494462118598, 0.00026408477907207745]
MinMaxDistance=[i*avgcostofTTIs for i in MinMaxDistance]
Total_Average_distance =[0.00018217233894207244, 0.00012885226534257752, 0.00007079474952419054, 0.00004968434992472714] # this is the avg deviation actually
Total_Average_distance=[i*avgcostofTTIs for i in Total_Average_distance]
Total_Average_distance=[i*len(F_continuous_array) for i in Total_Average_distance]


# Below is robust case
Total_Average_distance_Rob = [0.00041355578190442106, 0.0003613414747017453, 0.00022874837428733713, 0.00018543150185251265] # this is the avg deviation actually
Total_Average_distance_Rob=[i*avgcostofTTIs for i in Total_Average_distance_Rob]
Total_Average_distance_Rob=[i*len(F_continuous_array) for i in Total_Average_distance_Rob]



# Total_Average_cost=[i*avgcostofTTIs for i in Total_Average_distance]



MinMaxDistance_from_current_practice=[0.00443, 0.00447, 0.00458, 0.00461]
Total_Average_distance_from_current_practice=[0.00254, 0.00259, 0.00265, 0.00267]

K_allowable_changes_array= [1, 2, 3, 4]
T = np.arange(0,len(Prevalence_Data_List[0]),1)
# for m in range(len(MinMaxDistance)):
#     MinMaxDistance[m]=np.format_float_scientific(MinMaxDistance[m], exp_digits=2,precision=3)
#     Total_Average_distance[m]=np.format_float_scientific(Total_Average_distance[m], exp_digits=2, precision=3)
#     # Total_Average_cost[m]=np.format_float_scientific(Total_Average_cost[m], exp_digits=1, precision=2)

Plotting_function(F_continuous_array,F_constrained_matrix_of_ks,F_constrained_matrix_of_ks_Robust,path_matrix_of_ks,budget_path_3d_matrix_of_ksnew,K_allowable_changes_array,MinMaxDistance,Total_Average_distance,Total_Average_distance_Rob)