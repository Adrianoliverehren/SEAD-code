from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
df_ex = pd.ExcelFile("data.xlsx")
df1 = df_ex.parse("stability_data")
df2 = df_ex.parse("controlability_data")

control = {}
stability = {}

for index, row in df1.iterrows():
    control[row["Item"]] = row["Value1"]

for index, row in df2.iterrows():
    stability[row["Item"]] = row["Value"]


def find_flapped_wing_data(stability):
    #slide 38 ADSEE year 2
    delta_alpha0 = stability["delta_alpha0l_airfoil"] * np.cos(stability["sweep"]) *\
         (stability["Swf_TE"] / stability["S"])
    
    alpha0 = stability["alpha_L=0"] + delta_alpha0

    stability["C_L_alpha-wing"] = (2 * np.pi * stability["A"]) / (2 + (4 + \
        ((stability["A"] * stability["beta_stall"]) / stability["eta"])**2 * \
        (1 + (np.tan(stability["sweep"])**2) / stability["beta_stall"]**2))**0.5)

    stability["C_L_alpha-0"] = stability["C_L_alpha-wing"] * abs(alpha0 * (np.pi/180))
    
    print(stability["C_L_alpha-wing"])        
    print(delta_alpha0)
    print(stability["C_L_alpha-0"])

find_flapped_wing_data(stability)