import pandas as pd
import numpy as np
df = pd.ExcelFile("data.xlsx")
df = df.parse("data_for_weight")

vals = {}

for index, row in df.iterrows():
    vals[row["item"]] = row["value"]


vals["W_wing"] = 0.0051 * ((vals["W_dg"] * vals["N_z"]) ** 0.557) * (vals["S_w"] ** 0.649) * \
    (vals["A"] ** 0.5) * (vals["t_root"]/vals["c_root"]) ** -0.4 * \
    ((1 + vals["taper"]) ** 0.1) * ((np.cos(vals["sweep"]))** -1) * vals["S_csw"] ** 0.1

vals["W_hor_tail"] = 0.0379 * vals["K_uht"] * ((1 + vals["F_w"]/vals["B_h"]) ** -0.25) * (vals["W_dg"] ** 0.639) * \
    (vals["N_z"] ** 0.1) * (vals["S_ht"] ** 0.75) * (vals["L_t"] ** -1) * (vals["K_y"] ** 0.704) * \
        (np.cos(vals["sweep_ht"]) ** -1) * (vals["A_h"] ** 0.166) * ((1 + vals["S_e"]/vals["S_ht"]) ** 0.1)

vals["W_ver_tail"] = 0.0026 * ((1 + vals["H_t/H_v"]) ** 0.225) * (vals["W_dg"] ** 0.556) * (vals["N_z"] ** 0.536) *\
    (vals["L_t"] ** -0.5) * (vals["S_vt"] ** 0.5) * (vals["K_z"] ** 0.875) * (np.cos(vals["sweep_vt"]) ** -1) \
        * (vals["A_v"] ** 0.35) * (vals["t_root"]/vals["c_root"]) ** -0.5

vals["W_fus"] = 0.3280 * vals["K_door"] * vals["K_Lg"] * ((vals["W_dg"] * vals["N_z"]) ** 0.5) * (vals["L"] ** 0.25) *\
    (vals["S_f"] ** 0.302) * ((1 + vals["K_ws"]) ** 0.04) * ((vals["L"]/vals["D"]) ** 0.1)

vals["H_t/H_v"] #for vertical tail use this?????