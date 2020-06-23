
# <><><> SETUP IMPORTS <><><>

import pandas as pd
import numpy as np
import math as m
import os
from rpy2.robjects.packages import importr
import rpy2.robjects.packages as rpackages
from rpy2.robjects.vectors import StrVector
from rpy2.robjects import r
from scipy.spatial.distance import pdist,squareform
from rpy2.robjects import pandas2ri
from rpy2.robjects import Formula
import rpy2.robjects as robjects
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline

# Setup R packages used
pandas2ri.activate()

utils = rpackages.importr('utils')
utils.chooseCRANmirror(ind=1)

packnames = ('vegan','psych','lavaan')
names_to_install = []
#names_to_install = [x for packnames if not rpackages.isinstalled(x)]

for x in packnames:
    if (rpackages.isinstalled(x)==False):
        names_to_install.append(x)

if len(names_to_install) > 0:
    utils.install_packages(StrVector(names_to_install))

rvegan = importr('vegan')
rpsych = importr('psych')
rlavaan = importr('lavaan')
rsummary = r['summary']
rglm = r['glm']
rlm = r['lm']
rc = r['c']



# <><><> DEFINE FUNCTIONS <><><>

def string_find_replace(str_returned):
    replace_index = str_returned.find("Deviance Residuals")
    return str_returned[replace_index-1:]



def main():

    outDir = f"{os.path.dirname(os.path.realpath(__file__))}/output"
    # allows for cleaner execution and use of relative paths

    # Transforming the Inputs
    AE = sem_data["Archaea Evenness"] * 10
    AR = sem_data["Archaea Richness"] / 10
    BE = sem_data["Bacterial Evenness"] * 10
    BR = sem_data["Bacterial Richness"] / 100
    FE = sem_data["Fungal Evenness"] * 10
    FR = sem_data["Fungal Richness"] / 100
    PR = sem_data["Plant R"] / 10

    root_biomass = np.log(sem_data["Root_gm2"])
    litter_biomass = sem_data["Litter_gm2"] / 100
    litter_cn_ratio = sem_data["LitterCN_ratio"] / 10
    brome_shoots = sem_data["Brome_gm2"] / 100
    brome_roots = np.log(sem_data["Bromus_inermis"] + 1)
    degree_sum = (sem_data["Degree.Sum"] - min(sem_data["Degree.Sum"])) * 10
    eigen_value_sum = (sem_data["Eigenvalue.Sum"] - min(sem_data["Eigenvalue.Sum"])) * 10
    closeness_sum = (sem_data["Closeness.Sum"] - min(sem_data["Closeness.Sum"])) * 10
    betweenness_sum = (sem_data["Betweenness.Sum"] - min(sem_data["Betweenness.Sum"])) * 10

    sem_data2 = pd.concat(
        [sem_data["Horizon"], sem_data["A.horizon"], sem_data["Org_C"], sem_data["Soil_TotN"], PR, AE, AR, BE, BR, FE,
         FR, root_biomass, litter_biomass, litter_cn_ratio, brome_shoots, brome_roots, degree_sum, eigen_value_sum,
         closeness_sum, betweenness_sum, sem_data["pH"]], axis=1)

    sem_data2.rename(
        columns={'Bromus_inermis': 'BromeRoots', 'Brome_gm2': 'BromeShoots', 'Litter_gm2': 'Litter.Biomass',
                 'Root_gm2': 'Root.Biomass', 'Org_C': 'SOC', 'Soil_TotN': 'Soil.N', 'Plant R': 'PR',
                 'Archaea Evenness': 'AE', 'Archaea Richness': 'AR', 'Bacterial Evenness': 'BE',
                 'Bacterial Richness': 'BR', 'Fungal Evenness': 'FE', 'Fungal Richness': 'FR',
                 'LitterCN_ratio': 'Litter.CN.ratio'}, inplace=True)

    # <><><> PRINTING SUMMARIES <><><>

    print(string_find_replace(str(
        rsummary(rglm(Formula('Degree.Sum ~ BromeRoots + I(BromeRoots^2)'), data=sem_data2, family=r['gaussian'])))))
    print(string_find_replace(str(
        rsummary(rglm(Formula('Degree.Sum ~ BromeShoots + I(BromeShoots^2)'), data=sem_data2, family=r['gaussian'])))))
    print(string_find_replace(
        str(rsummary(rglm(Formula('BE ~ Degree.Sum + I(Degree.Sum^2)'), data=sem_data2, family=r['gaussian'])))))
    print(string_find_replace(str(
        rsummary(rglm(Formula('BE ~ Eigenvalue.Sum + I(Eigenvalue.Sum^2)'), data=sem_data2, family=r['gaussian'])))))
    print(string_find_replace(
        str(rsummary(rglm(Formula('BE ~ Closeness.Sum + I(Closeness.Sum^2)'), data=sem_data2, family=r['gaussian'])))))
    print(string_find_replace(str(
        rsummary(rglm(Formula('BE ~ Betweenness.Sum + I(Betweenness.Sum^2)'), data=sem_data2, family=r['gaussian'])))))
    print(string_find_replace(str(rsummary(rglm(Formula('FR ~ AR + I(AR^2)'), data=sem_data2, family=r['gaussian'])))))
    print(string_find_replace(
        str(rsummary(rglm(Formula('FR ~ BromeShoots + I(BromeShoots^2)'), data=sem_data2, family=r['gaussian'])))))
    print(string_find_replace(
        str(rsummary(rglm(Formula('FR ~ BromeRoots + I(BromeRoots^2)'), data=sem_data2, family=r['gaussian'])))))
    print(string_find_replace(
        str(rsummary(rglm(Formula('AR ~ BromeShoots + I(BromeShoots^2)'), data=sem_data2, family=r['gaussian'])))))
    print(string_find_replace(
        str(rsummary(rglm(Formula('AR ~ BromeRoots + I(BromeRoots^2)'), data=sem_data2, family=r['gaussian'])))))
    print(string_find_replace(str(rsummary(rglm(Formula('FE ~ AE + I(AE^2)'), data=sem_data2, family=r['gaussian'])))))
    print(string_find_replace(
        str(rsummary(rglm(Formula('FE ~ BromeShoots + I(BromeShoots^2)'), data=sem_data2, family=r['gaussian'])))))
    print(string_find_replace(
        str(rsummary(rglm(Formula('FE ~ BromeRoots + I(BromeRoots^2)'), data=sem_data2, family=r['gaussian'])))))
    print(string_find_replace(
        str(rsummary(rglm(Formula('AE ~ BromeShoots + I(BromeShoots^2)'), data=sem_data2, family=r['gaussian'])))))
    print(string_find_replace(
        str(rsummary(rglm(Formula('AE ~ BromeRoots + I(BromeRoots^2)'), data=sem_data2, family=r['gaussian'])))))
    print(string_find_replace(str(rsummary(rglm(Formula('AR ~ pH + I(pH^2)'), data=sem_data2, family=r['gaussian'])))))
    print(string_find_replace(str(rsummary(rglm(Formula('FR ~ pH + I(pH^2)'), data=sem_data2, family=r['gaussian'])))))

    print(string_find_replace(str(rsummary(rlm(Formula('Degree.Sum ~ BromeRoots'), data=sem_data2)))))
    print(string_find_replace(str(rsummary(rlm(Formula('Degree.Sum ~ BromeShoots'), data=sem_data2)))))
    print(string_find_replace(str(rsummary(rlm(Formula('BE ~ Degree.Sum'), data=sem_data2)))))
    print(string_find_replace(str(rsummary(rlm(Formula('BE ~ Eigenvalue.Sum'), data=sem_data2)))))
    print(string_find_replace(str(rsummary(rlm(Formula('BE ~ Closeness.Sum'), data=sem_data2)))))
    print(string_find_replace(str(rsummary(rlm(Formula('BE ~ Betweenness.Sum'), data=sem_data2)))))
    print(string_find_replace(str(rsummary(rlm(Formula('FR ~ AR'), data=sem_data2)))))
    print(string_find_replace(str(rsummary(rlm(Formula('FR ~ BromeShoots'), data=sem_data2)))))
    print(string_find_replace(str(rsummary(rlm(Formula('FR ~ BromeRoots'), data=sem_data2)))))
    print(string_find_replace(str(rsummary(rlm(Formula('AR ~ BromeShoots'), data=sem_data2)))))
    print(string_find_replace(str(rsummary(rlm(Formula('AR ~ BromeRoots'), data=sem_data2)))))
    print(string_find_replace(str(rsummary(rlm(Formula('FE ~ AE'), data=sem_data2)))))
    print(string_find_replace(str(rsummary(rlm(Formula('FE ~ BromeShoots'), data=sem_data2)))))
    print(string_find_replace(str(rsummary(rlm(Formula('FE ~ BromeRoots'), data=sem_data2)))))
    print(string_find_replace(str(rsummary(rlm(Formula('AE ~ BromeShoots'), data=sem_data2)))))
    print(string_find_replace(str(rsummary(rlm(Formula('AE ~ BromeRoots'), data=sem_data2)))))
    print(string_find_replace(str(rsummary(rlm(Formula('AR ~ pH'), data=sem_data2)))))
    print(string_find_replace(str(rsummary(rlm(Formula('FR ~ pH'), data=sem_data2)))))


    # <><><> PLOTTING FIGURES <><><>

    plt.scatter(sem_data2["BromeRoots"], sem_data2["Degree.Sum"], facecolors='none', edgecolors='black')
    plt.xlabel('BromeRoots')
    plt.ylabel('Degree.Sum')
    plt.show()

    plt.scatter(sem_data2["BromeShoots"], sem_data2["Degree.Sum"], facecolors='none', edgecolors='black')
    plt.xlabel('BromeShoots')
    plt.ylabel('Degree.Sum')
    plt.show()

    y = [(rlm(Formula('BE ~ Degree.Sum'), data=sem_data2))[0][0]]
    val = y[0]
    for i in range(1, m.ceil(max(sem_data2["Degree.Sum"])) + 1):
        val = val + (rlm(Formula('BE ~ Degree.Sum'), data=sem_data2))[0][1]
        y.append(val)
    plt.scatter(sem_data2["Degree.Sum"], sem_data2["BE"], facecolors='none', edgecolors='black')
    plt.plot([x for x in range(0, m.ceil(max(sem_data2["Degree.Sum"])) + 1)], y)
    plt.xlabel('Degree.Sum')
    plt.ylabel('BE')
    plt.show()

    y = [(rlm(Formula('BE ~ Eigenvalue.Sum'), data=sem_data2))[0][0]]
    val = y[0]
    for i in range(1, m.ceil(max(sem_data2["Eigenvalue.Sum"])) + 1):
        val = val + (rlm(Formula('BE ~ Eigenvalue.Sum'), data=sem_data2))[0][1]
        y.append(val)
    plt.scatter(sem_data2["Eigenvalue.Sum"], sem_data2["BE"], facecolors='none', edgecolors='black')
    plt.plot([x for x in range(0, m.ceil(max(sem_data2["Eigenvalue.Sum"])) + 1)], y)
    plt.xlabel('Eigenvalue.Sum')
    plt.ylabel('BE')
    plt.show()

    y = [(rlm(Formula('BE ~ Closeness.Sum'), data=sem_data2))[0][0]]
    val = y[0]
    for i in range(1, m.ceil(max(sem_data2["Closeness.Sum"])) + 1):
        val = val + (rlm(Formula('BE ~ Closeness.Sum'), data=sem_data2))[0][1]
        y.append(val)
    plt.scatter(sem_data2["Closeness.Sum"], sem_data2["BE"], facecolors='none', edgecolors='black')
    plt.plot([x for x in range(0, m.ceil(max(sem_data2["Closeness.Sum"])) + 1)], y)
    plt.xlabel('Closeness.Sum')
    plt.ylabel('BE')
    plt.show()

    y = [(rlm(Formula('BE ~ Betweenness.Sum'), data=sem_data2))[0][0]]
    val = y[0]
    for i in range(1, m.ceil(max(sem_data2["Betweenness.Sum"])) + 1):
        val = val + (rlm(Formula('BE ~ Betweenness.Sum'), data=sem_data2))[0][1]
        y.append(val)
    plt.scatter(sem_data2["Betweenness.Sum"], sem_data2["BE"], facecolors='none', edgecolors='black')
    plt.plot([x for x in range(0, m.ceil(max(sem_data2["Betweenness.Sum"])) + 1)], y)
    plt.xlabel('Betweenness.Sum')
    plt.ylabel('BE')
    plt.show()

    plt.scatter(sem_data2["AR"], sem_data2["FR"], facecolors='none', edgecolors='black')
    plt.xlabel('AR')
    plt.ylabel('FR')
    plt.show()

    plt.scatter(sem_data2["BromeShoots"], sem_data2["FR"], facecolors='none', edgecolors='black')
    plt.xlabel('BromeShoots')
    plt.ylabel('FR')
    plt.show()

    plt.scatter(sem_data2["BromeRoots"], sem_data2["FR"], facecolors='none', edgecolors='black')
    plt.xlabel('BromeRoots')
    plt.ylabel('FR')
    plt.show()

    y = [(rlm(Formula('AR ~ BromeShoots'), data=sem_data2))[0][0]]
    val = y[0]
    for i in range(1, m.ceil(max(sem_data2["BromeShoots"])) + 1):
        val = val + (rlm(Formula('AR ~ BromeShoots'), data=sem_data2))[0][1]
        y.append(val)
    plt.scatter(sem_data2["BromeShoots"], sem_data2["AR"], facecolors='none', edgecolors='black')
    plt.plot([x for x in range(0, m.ceil(max(sem_data2["BromeShoots"])) + 1)], y)
    plt.xlabel('BromeShoots')
    plt.ylabel('AR')
    plt.show()

    plt.scatter(sem_data2["BromeRoots"], sem_data2["AR"], facecolors='none', edgecolors='black')
    plt.xlabel('BromeRoots')
    plt.ylabel('AR')
    plt.show()

    plt.scatter(sem_data2["AE"], sem_data2["FE"], facecolors='none', edgecolors='black')
    plt.xlabel('AE')
    plt.ylabel('FE')
    plt.show()

    plt.scatter(sem_data2["BromeShoots"], sem_data2["FE"], facecolors='none', edgecolors='black')
    plt.xlabel('BromeShoots')
    plt.ylabel('FE')
    plt.show()

    plt.scatter(sem_data2["BromeRoots"], sem_data2["FE"], facecolors='none', edgecolors='black')
    plt.xlabel('BromeRoots')
    plt.ylabel('FE')
    plt.show()

    plt.scatter(sem_data2["BromeShoots"], sem_data2["AE"], facecolors='none', edgecolors='black')
    plt.xlabel('BromeShoots')
    plt.ylabel('AE')
    plt.show()

    plt.scatter(sem_data2["BromeRoots"], sem_data2["AE"], facecolors='none', edgecolors='black')
    plt.xlabel('BromeRoots')
    plt.ylabel('AE')
    plt.show()

    plt.scatter(sem_data2["pH"], sem_data2["AR"], facecolors='none', edgecolors='black')
    plt.xlabel('pH')
    plt.ylabel('AR')
    plt.show()

    plt.scatter(sem_data2["pH"], sem_data2["FR"], facecolors='none', edgecolors='black')
    plt.xlabel('pH')
    plt.ylabel('FR')
    plt.show()


    # <><><> OUTPUT RESULTS <><><>

    sem_data2a = sem_data2.loc[(sem_data2["Horizon"] == "A")]
    sem_data2b = sem_data2.loc[(sem_data2["Horizon"] == "B")]
    sem_data2a['Horizon'] = sem_data2a['Horizon'].astype('category')

    sem_data2a.rename(columns={'Litter.Biomass': 'LitterBiomass', 'Root.Biomass': 'RootBiomass', 'Soil.N': 'SoilN',
                               'Litter.CN.ratio': 'LitterCNratio', 'A.horizon': 'AHorizon', 'Degree.Sum': 'DegreeSum'},
                      inplace=True)
    sem_data2a.to_csv("sem_data2a.csv", index=False)

    sem_data2b.rename(columns={'Litter.Biomass': 'LitterBiomass', 'Root.Biomass': 'RootBiomass', 'Soil.N': 'SoilN',
                               'Litter.CN.ratio': 'LitterCNratio', 'A.horizon': 'AHorizon', 'Degree.Sum': 'DegreeSum'},
                      inplace=True)
    sem_data2b.to_csv("sem_data2b.csv", index=False)


# <><><> TEST <><><>

sem_data = pd.read_csv("Brome.SEM.data.csv")
