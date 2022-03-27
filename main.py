from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
df_ex = pd.ExcelFile("data.xlsx")
df1 = df_ex.parse("data_for_weight")
df2 = df_ex.parse("x_locations")

vals = {}
xcg_loc = {}
vals_si = {}

for index, row in df1.iterrows():
    vals[row["item"]] = row["value"]
    vals_si[row["item"]] = row["value_SI"]

for index, row in df2.iterrows():
    xcg_loc[row["Item"]] = [row["x_cg (from the nose) [m]"], row["x_cg/mac [-]"]]

lbs_kg = 0.45359237
ft_m = 0.3048

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

def find_x_cg_oew(vals, xcg_loc):
    x_weight = 0
    xc_weight = 0
    weight_sum = 0
    for id in OEW_ids:
        xc_weight += vals[id] * xcg_loc[id][1]
        x_weight += vals[id] * xcg_loc[id][0]
        weight_sum += vals[id]

    x_cg = x_weight / weight_sum
    xc_cg = xc_weight / weight_sum

    return x_cg, xc_cg

def potato_diagram(xc_cg_oew, vals, xcg_loc):
    # front cargo
    front_cargo_data = [[],[]] #W, xc_cg
    aft_cargo_data = [[],[]] #W, xc_cg
    front_cargo_data[0].append(vals["W_empty_calculated"] * lbs_kg)
    front_cargo_data[1].append(xc_cg_oew)
    aft_cargo_data[0].append(vals["W_empty_calculated"] * lbs_kg)
    aft_cargo_data[1].append(xc_cg_oew)

    # only front cargo
    new_mass = vals["W_empty_calculated"] + vals["W_front_cargo"]
    front_cargo_data[0].append(new_mass * lbs_kg)
    new_xc_cg = ((vals["W_empty_calculated"] * xc_cg_oew + \
        vals["W_front_cargo"] * xcg_loc["W_front_cargo"][1]) / new_mass)
    front_cargo_data[1].append(new_xc_cg)

    # only rear cargo
    new_mass = vals["W_empty_calculated"] + vals["W_aft_cargo"]
    aft_cargo_data[0].append(new_mass * lbs_kg)
    new_xc_cg = ((vals["W_empty_calculated"] * xc_cg_oew + \
        vals["W_aft_cargo"] * xcg_loc["W_aft_cargo"][1]) / new_mass)
    aft_cargo_data[1].append(new_xc_cg)

    # front and rear cargo
    final_mass = vals["W_empty_calculated"] + vals["W_aft_cargo"] + vals["W_front_cargo"]
    aft_cargo_data[0].append(final_mass * lbs_kg)
    front_cargo_data[0].append(final_mass * lbs_kg)
    final_xc_cg = ((vals["W_empty_calculated"] * xc_cg_oew + \
        vals["W_aft_cargo"] * xcg_loc["W_aft_cargo"][1] + vals["W_front_cargo"] * xcg_loc["W_front_cargo"][1]) / final_mass)
    aft_cargo_data[1].append(final_xc_cg)
    front_cargo_data[1].append(final_xc_cg)

    after_cargo_cg = final_xc_cg
    after_cargo_mass = final_mass * lbs_kg
    
    def wind_front_back(after_cargo_cg, after_cargo_mass):
    #add people Window - front to back  
        new_x_cg = after_cargo_cg
        current_mass = after_cargo_mass
        data_wind_front_back = [[],[]]
            
        data_wind_front_back[0].append(current_mass)
        data_wind_front_back[1].append(new_x_cg)
            

        for i in range(25):
            current_row_pos = xcg_loc["front_seat"][1] + i * xcg_loc["d_row"][1]
            mass_before_wind_add = after_cargo_mass + 2*i * vals_si["W_person"]
            
            current_mass += 2 * vals_si["W_person"]

            new_x_cg = (current_row_pos * 2 * vals_si["W_person"] + new_x_cg * mass_before_wind_add) / (current_mass)
            
            data_wind_front_back[0].append(current_mass)
            data_wind_front_back[1].append(new_x_cg)

        return data_wind_front_back

    def wind_back_front(after_cargo_cg, after_cargo_mass):
    #add people Window - front to back  
        new_x_cg = after_cargo_cg
        current_mass = after_cargo_mass
        data_wind_back_front = [[],[]]
            
        data_wind_back_front[0].append(current_mass)
        data_wind_back_front[1].append(new_x_cg)
            

        for i in range(25):
            current_row_pos = xcg_loc["aft_seat"][1] - i * xcg_loc["d_row"][1]
            mass_before_wind_add = after_cargo_mass + 2*i * vals_si["W_person"]
            current_mass += 2 * vals_si["W_person"]
            new_x_cg = (current_row_pos * 2 * vals_si["W_person"] + new_x_cg * mass_before_wind_add)/ (current_mass)
            
            data_wind_back_front[0].append(current_mass)
            data_wind_back_front[1].append(new_x_cg)    

        return data_wind_back_front

    def aisle_front_back(after_cargo_cg, after_cargo_mass):
    #add people ailse - front to back  
        new_x_cg = after_cargo_cg
        current_mass = after_cargo_mass
        data_aisle_front_back = [[],[]]
            
        data_aisle_front_back[0].append(current_mass)
        data_aisle_front_back[1].append(new_x_cg)
            

        for i in range(25):
            current_row_pos = xcg_loc["front_seat"][1] + i * xcg_loc["d_row"][1]
            mass_before_aisle_add = after_cargo_mass + 2*i * vals_si["W_person"]
            
            current_mass += 2 * vals_si["W_person"]

            new_x_cg = (current_row_pos * 2 * vals_si["W_person"] + new_x_cg * mass_before_aisle_add) / (current_mass)
            
            data_aisle_front_back[0].append(current_mass)
            data_aisle_front_back[1].append(new_x_cg)

        return data_aisle_front_back

    def aisle_back_front(after_cargo_cg, after_cargo_mass):
    #add people aisle - front to back  
        new_x_cg = after_cargo_cg
        current_mass = after_cargo_mass
        data_aisle_back_front = [[],[]]
            
        data_aisle_back_front[0].append(current_mass)
        data_aisle_back_front[1].append(new_x_cg)
            

        for i in range(25):
            current_row_pos = xcg_loc["aft_seat"][1] - i * xcg_loc["d_row"][1]
            mass_before_aisle_add = after_cargo_mass + 2*i * vals_si["W_person"]
            current_mass += 2 * vals_si["W_person"]
            new_x_cg = (current_row_pos * 2 * vals_si["W_person"] + new_x_cg * mass_before_aisle_add)/ (current_mass)
            
            data_aisle_back_front[0].append(current_mass)
            data_aisle_back_front[1].append(new_x_cg)    

        return data_aisle_back_front

    # def add_power_juice(start_mass, start_x_cg):

    data_wind_front_back = wind_front_back(after_cargo_cg, after_cargo_mass)
    data_wind_back_front = wind_back_front(after_cargo_cg, after_cargo_mass)
    data_aisle_front_back = aisle_front_back(data_wind_back_front[1][-1], data_wind_back_front[0][-1])
    data_aisle_back_front = aisle_back_front(data_wind_back_front[1][-1], data_wind_back_front[0][-1])

    # adding fuel

    power_juice_data = [[],[]]
    power_juice_data[0].append(data_aisle_back_front[0][-1])
    power_juice_data[1].append(data_aisle_back_front[1][-1])

    new_cg = (data_aisle_back_front[1][-1] * data_aisle_back_front[0][-1] + vals_si["W_fuel_max"] * xcg_loc["W_fuel"][1]) / (data_aisle_back_front[0][-1] + vals_si["W_fuel_max"])
    new_weight = data_aisle_back_front[0][-1] + vals_si["W_fuel_max"]

    power_juice_data[0].append(new_weight)
    power_juice_data[1].append(new_cg)

    all_cgs = power_juice_data[1] + data_aisle_back_front[1] + data_aisle_front_back[1] + data_wind_back_front[1] + data_wind_front_back[1] +\
        aft_cargo_data[1] + front_cargo_data[1]

    min_cg = min(all_cgs)
    max_cg = max(all_cgs)
    min_cg = min_cg * 1.02
    max_cg = max_cg * 1.02

    print("MAX C.G. position is: " + str(max_cg))
    print("MIN C.G. position is: " + str(min_cg))

# x_seat_front
# x_seat_aft
# d_row


    #add people Window - back to front   
    

    #cargo lines added
    plt.plot(front_cargo_data[1], front_cargo_data[0], label="front cargo", marker='o')
    plt.plot(aft_cargo_data[1], aft_cargo_data[0], label="aft cargo", marker='o')
    
    # window potato added
    plt.plot(data_wind_front_back[1], data_wind_front_back[0], label="window seats front first", marker='o')
    plt.plot(data_wind_back_front[1], data_wind_back_front[0], label="window seats aft first", marker='o')

    # aisle potato added
    plt.plot(data_aisle_front_back[1], data_aisle_front_back[0], label="ailse seats front first", marker='o')
    plt.plot(data_aisle_back_front[1], data_aisle_back_front[0], label="ailse seats aft first", marker='o')
    
    # adding power juice
    plt.plot(power_juice_data[1], power_juice_data[0], label="fuel", marker='o')

    plt.vlines(min_cg, ymin=21974, ymax=vals_si["W_mtow"])
    plt.vlines(max_cg, ymin=21974, ymax=vals_si["W_mtow"])

    plt.plot
    plt.legend()
    plt.ylabel("mass [kg]")
    plt.xlabel("x_cg [mac]")
    plt.ylim(21974)
    plt.grid(True)

    plt.show()

    pass


x_cg_oew, xc_cg_oew = find_x_cg_oew(vals, xcg_loc)


potato_diagram(xc_cg_oew, vals, xcg_loc)
        