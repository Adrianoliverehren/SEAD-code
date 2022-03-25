from wsgiref import validate
import pandas as pd
import numpy as np
df_ex = pd.ExcelFile("data.xlsx")
df1 = df_ex.parse("data_for_weight")
df2 = df_ex.parse("x_locations")

vals = {}
xcg_loc = {}

for index, row in df1.iterrows():
    vals[row["item"]] = row["value"]

for index, row in df2.iterrows():
    xcg_loc[row["Item"]] = [row["x_cg (from the nose) [m]"], row["x_cg/mac [-]"]]

lbs_kg = 0.45359237
ft_m = 0.3048

# print(xcg_loc)

vals["W_wing"] = 2 * 0.0051 * ((vals["W_dg"] * vals["N_z"]) ** 0.557) * (vals["S_w"] ** 0.649) * \
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

vals["W_main_gear"] = 0.0106 * vals["K_mp"] * vals["W_l"] ** 0.888 * vals["N_l"] ** 0.25 * vals["L_m"] ** 0.4 *\
    vals["N_mw"] ** 0.321 * vals["N_mss"] ** -0.5 * vals["V_stall"] ** 0.1

vals["W_nose_gear"] = 0.032 * vals["K_np"] * vals["W_l"] ** 0.646 * vals["N_l"] ** 0.2 * vals["L_n"] ** 0.5 * vals["N_nw"] ** 0.45

vals["W_nacelle"] = 0.6724 * vals["K_ng"] * vals["N_Lt"] ** 0.1 * vals["N_w"] ** 0.294 * vals["N_z"] ** 0.119 * vals["W_ec"] ** 0.611 *\
    vals["N_en"] ** 0.984 * vals["S_n"] ** 0.224

vals["W_engine_controls"] = 5 * vals["N_en"] + 0.8 * vals["L_ec"]

vals["W_starter"] = 49.19 * ((vals["N_en"] * vals["W_en"])/1000) ** 0.541

vals["W_fuel_sys"] = 2.405 * vals["V_t"] ** 0.606 * (1 + vals["V_i"] / vals["V_t"]) ** -1 * (1 + vals["V_p"] / vals["V_t"]) * vals["N_t"] ** 0.5

vals["W_flight_controls"] = 145.9 * vals["N_f"] ** 0.554 * (1 + vals["N_m"]/vals["N_f"]) ** -1 * vals["S_cs"] ** 0.2 * (vals["I_y"] * 10 ** -6) ** 0.07

vals["W_APU_installed"] = 2.2 * vals["W_APU_uninstalled"]

vals["W_instruments"] = 4.509 * vals["K_r"] * vals["K_tp"] * vals["N_c"] ** 0.541 * vals["N_en"] * (vals["L_f"] + vals["B_w"]) ** 0.5

vals["W_hydraulics"] = 0.2673 * vals["N_f"] * (vals["L_f"] + vals["B_w"]) * 0.937

vals["W_electrical"] = 7.291 * vals["R_kva"] ** 0.782 * vals["L_a"] ** 0.346 * vals["N_gen"] ** 0.1

vals["W_avionics"] = 1.73 * vals["W_uav"] ** 0.983

vals["W_furnishings"] = 0.0577 * vals["N_c"] ** 0.1 * vals["W_c"] ** 0.393 * vals["S_f"] ** 0.75

vals["W_seats"] = vals["N_p"] * vals["W_seat"]

vals["W_air_conditioning"] = 62.36 * vals["N_p"] ** 0.25 * (vals["V_pr"] / 1000) ** 0.604 * vals["W_uav"] ** 0.1

vals["W_anti-ice"] = 0.002 * vals["W_dg"]

vals["W_handling_gear"] = 3 * 10 ** -4 * vals["W_dg"]

vals["W_crew"] = vals["N_c"] * vals["W_person"]

vals["W_ec"] = vals["W_ec"] * 2

OEW_ids = ["W_wing", "W_hor_tail", "W_ver_tail", "W_fus", "W_main_gear", "W_nose_gear", "W_nacelle", \
        "W_engine_controls", "W_starter", "W_fuel_sys", "W_flight_controls", "W_APU_installed", \
        "W_instruments", "W_hydraulics", "W_electrical", "W_avionics", "W_furnishings",\
        "W_air_conditioning", "W_anti-ice", "W_handling_gear", "W_ec", "W_crew", "W_seats"]

vals["W_empty_calculated"] = 0
#print(len(OEW_ids))
for id in OEW_ids:
    vals["W_empty_calculated"] += vals[id]

#print(vals["W_empty_calculated"] * lbs_kg, vals["W_oew"] * lbs_kg)
#print((vals["W_empty_calculated"] * lbs_kg - vals["W_oew"] * lbs_kg) / (vals["W_oew"] * lbs_kg) * 100)



def find_x_cg_oew(vals, xcg_loc):
    x_weight = 0
    weight_sum = 0
    for id in OEW_ids:
        x_weight += vals[id] * xcg_loc[id][1]
        weight_sum += vals[id]

    x_cg = x_weight / weight_sum

    return x_cg

def potato_diagram(x_cg_oew, vals, xcg_loc):
    data = {}
    print(xcg_loc["W_front_cargo"][1])
    #adding cargo
    x_cg_W = x_cg_oew * vals["W_empty_calculated"]
    x_cg = (x_cg_W + xcg_loc["W_front_cargo"][1] * vals["W_front_cargo"]) /\
        (vals["W_front_cargo"] + vals["W_empty_calculated"])
    data["front cargo"] = [x_cg_oew, x_cg]
    x_cg = (x_cg_W + xcg_loc["W_aft_cargo"][1] * vals["W_aft_cargo"]) /\
        (vals["W_aft_cargo"] + vals["W_empty_calculated"])
    data["aft cargo"] = [x_cg_oew, x_cg]
    x_cg = (x_cg_W + xcg_loc["W_aft_cargo"][1] * vals["W_aft_cargo"]+ xcg_loc["W_front_cargo"][1] * vals["W_front_cargo"]) /\
        (vals["W_front_cargo"] + vals["W_aft_cargo"] + vals["W_empty_calculated"])
    data["front cargo"].append(x_cg)
    data["aft cargo"].append(x_cg)

    print(data)




    pass



x_cg_oew = find_x_cg_oew(vals, xcg_loc)
print(x_cg_oew)
potato_diagram(x_cg_oew, vals, xcg_loc)
        