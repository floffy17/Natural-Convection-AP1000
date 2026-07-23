# ASSIGNMENT 1.2
import math
import CoolProp.CoolProp as CP

# --- GENERAL DATA ---
P_nom = 600e6  # [W_th] Nominal thermal power
decay_power_steady = 0.01  # [-] Quasi-steady state power fraction (1%)
qq = P_nom * decay_power_steady  # 1% della potenza nominale (6 MW)
gg = 9.81
delta_p_vessel_nom = 1.2  # [bar] Vessel pressure drop at nominal flow rate

delta_p_vessel_nom_Pa = delta_p_vessel_nom * 1e5  # from bar to Pa
vessel_m_nom = 3200  # [kg/s] Reference mass flow rate for vessel pressure drop
R_vessel = delta_p_vessel_nom_Pa / (vessel_m_nom ** 2)

# --- PSC (Primary System Circuit) DATA ---
psc_press = 75 * 1e5  # [Pa]  Coolant pressure
psc_outpipe_diam = 16 * 0.0254  # Pipe outside diameter (Sch 100), from inches to meters
psc_pipe_thick = 1.031 * 0.0254  # Pipe thickness (from ANSI Table Sch 100 for 16"), from inches to meters
psc_L = 8  # [m] Half pipe length
psc_pipe_rough = 2e-4  # [-] Pipe relative roughness
psc_kloss_bend = 0.45  # [-] Pressure drop coefficient for 90° bends
psc_Nbends = 4  # [-] Number of 90° bends
psc_Klossvalve = 0.12  # [-] Check valve pressure drop coefficient
psc_H_hotcold = 7  # [m] Elevation between hot and cold leg
psc_H_topbot = 3  # [m] Elevation between top and bottom of the core

psc_inpipe_diam = psc_outpipe_diam - 2 * psc_pipe_thick  # [m] Pipe inside diameter
psc_A = (math.pi * psc_inpipe_diam ** 2) / 4  # [m2] PSC area

# --- ISC (Intermediate System Circuit) DATA ---
isc_press = 70 * 1e5  # [Pa] Coolant
isc_L = 20  # [m] Half pipe length
isc_pipe_rough = 2e-4  # [-] Pipe relative roughness
isc_kloss_bend = 0.45  # [-] Pressure drop coefficient for 90° bends
isc_Nbends = 2  # [-] Number of 90° bends
isc_H = 10  # [m] Net elevation
isc_inpipe_diam = psc_inpipe_diam
isc_A = psc_A

# --- HX1 (Heat Exchanger 1) DATA ---
hx1_D_out = 19.05e-3  # [m] Tube outside diameter
hx1_tube_thick = 1.24e-3  # [m] Tube thickness
hx1_Ntubes = 897  # [-] Number of tubes
hx1_pitch = 28.5e-3  # [m] Tube pitch (triangular)
hx1_Dshell = 1.5  # [m] Shell inside diameter
hx1_Nbaffles = 2  # [-] Number of baffles
hx1_L_baffles = 1.6  # [-] Number of baffles
hx1_av_tubelength = 9.314  # [m] Average tube length
hx1_kth = 15  # [W/mK] Tube thermal conductivity
hx1_kloss_bend = 0.45  # [-] Pressure drop coefficient for 90° bends
hx1_Nbends = 2  # [-] Number of 90° bends
hx1_aveflowarea = 0.883  # header area
hx1_heatrarea = 500  # [m^2] Total heat transfer area
hx1_tube_rough = 1e-4  # -] Tube relative roughness
hx1_Ft = 0.7  # [-] LMTD correction factor

hx1_in_diam = hx1_D_out - 2 * hx1_tube_thick  # [m] Tube inside diameter  # [m]
hx1_A_tot = hx1_Ntubes * (math.pi * hx1_in_diam ** 2) / 4  # [m^2]

# shell parameters
hx1_D_eq_shell = (2 * math.sqrt(3) * hx1_pitch ** 2) / (math.pi * hx1_D_out) - hx1_D_out  # [m]
hx1_A_shell = (hx1_Dshell / hx1_pitch) * (hx1_pitch - hx1_D_out) * hx1_L_baffles  # [m^2]

hx1_A_int_lat = math.pi * hx1_av_tubelength * hx1_in_diam * hx1_Ntubes  # [m^2] Total internal lateral surface area
hx1_A_out_lat = math.pi * hx1_av_tubelength * hx1_D_out * hx1_Ntubes   # [m^2] Total external lateral surface area

# --- HX2 DATA ---
hx2_out_diam = 25.4e-3  # [m] Tube outside diameter
hx2_tube_thick = 1.24e-3  # [m] Tube thickness
hx2_Ntubes = 770  # [-] Number of tubes
hx2_av_tubelength = 7  # m] Average tube length
hx2_kth = 15  # [W/mK] Tube thermal conductivity
hx2_heatrarea = 430  # [m^2] Heat transfer area
hx2_tuberough = 1e-4  # [-] Tube relative roughness

hx2_in_diam = hx2_out_diam - 2 * hx2_tube_thick  # [m] Tube inside diameter
hx2_A_int_lat = math.pi * hx2_av_tubelength * hx2_in_diam * hx2_Ntubes   # [m^2] Total internal lateral surface area
hx2_A_out_lat = math.pi * hx2_av_tubelength * hx2_out_diam * hx2_Ntubes  # [m^2] Total external lateral surface area
hx2_A_tot = hx2_Ntubes * (math.pi * hx2_in_diam ** 2) / 4 # [m^2]

# FUNCTIONS
def estrai_parametri(temp_c, press_pa, fluido="Water"):
    temp_k = temp_c + 273.15
    dens = CP.PropsSI('D', 'T', temp_k, 'P', press_pa, fluido)
    visc = CP.PropsSI('V', 'T', temp_k, 'P', press_pa, fluido)
    pr = CP.PropsSI('Prandtl', 'T', temp_k, 'P', press_pa, fluido)
    cond = CP.PropsSI('L', 'T', temp_k, 'P', press_pa, fluido)
    cp = CP.PropsSI('C', 'T', temp_k, 'P', press_pa, fluido)
    return [dens, visc, pr, cp, cond]


def f_haaland(eps_D, Re):
    if Re < 2300: return 64 / Re
    return (1 / (-1.8 * math.log10(((eps_D) / 3.7) ** 1.11 + (6.9 / Re)))) ** 2

# INITIAL GUESSES FOR ALGORITHM
mguess = 100.0
isc_Tav_guess = 120.0
toll = 1e-10

# =========================================================================
# ISC SOLUTION
# =========================================================================

while True:
    m_old = mguess
    while True:
        T_old = isc_Tav_guess
        res = estrai_parametri(isc_Tav_guess, isc_press)
        isc_dens, isc_visc, isc_PR, isc_cp, isc_cond = res

        # Delta ISC with power balance
        deltat_isc = qq / (mguess * isc_cp)

        # h int and h out
        # Adimensional numbers HX2
        hx2_Re = (4 * mguess) / (isc_visc * math.pi * hx2_in_diam * hx2_Ntubes)
        hx2_Nu = 0.023 * (hx2_Re ** 0.8) * (isc_PR ** 0.3)
        hx2_hint = hx2_Nu * isc_cond / hx2_in_diam

        pool_T_sat = CP.PropsSI('T', 'P', 1e5, 'Q', 0, 'Water') - 273.15
        deltaT_sat = (qq / (2.257 * hx2_heatrarea)) ** (1 / 3.86)
        hx2_hout = qq / (hx2_heatrarea * deltaT_sat)

        # global HT HX2
        R_tot_hx2 = (1 / (hx2_hint * hx2_A_int_lat) +
                     math.log(hx2_out_diam / hx2_in_diam) / (2 * math.pi * hx2_kth * hx2_Ntubes * hx2_av_tubelength) +
                     1 / (hx2_hout * hx2_A_out_lat))
        hx2_U = 1 / (R_tot_hx2 * hx2_A_out_lat)

        hx2_deltaT_ml = qq / (hx2_heatrarea * hx2_U)

        isc_Tcold = pool_T_sat + deltat_isc / (math.exp(deltat_isc / hx2_deltaT_ml) - 1)
        isc_Thot = isc_Tcold + deltat_isc
        isc_Tav_new = (isc_Thot + isc_Tcold) / 2

        if abs(isc_Tav_new - T_old) / abs(T_old) < toll:
            isc_Tav_guess = isc_Tav_new
            break
        isc_Tav_guess = isc_Tav_new

    # ISC properties
    res_h = estrai_parametri(isc_Thot, isc_press)
    res_c = estrai_parametri(isc_Tcold, isc_press)
    isc_dens_hot, isc_visc_hot = res_h[0], res_h[1]
    isc_dens_cold, isc_visc_cold = res_c[0], res_c[1]

    # Reynolds numbers
    Re_h_isc = (4 * mguess) / (isc_visc_hot * math.pi * isc_inpipe_diam)
    Re_c_isc = (4 * mguess) / (isc_visc_cold * math.pi * isc_inpipe_diam)
    hx1_Re_shell = (mguess * hx1_D_eq_shell) / (isc_visc * hx1_A_shell)

    # Friction factors
    f_hot = f_haaland(isc_pipe_rough, (4 * mguess) / (isc_visc_hot * math.pi * isc_inpipe_diam))
    f_cold = f_haaland(isc_pipe_rough, (4 * mguess) / (isc_visc_cold * math.pi * isc_inpipe_diam))
    hx2_f = f_haaland(hx2_tuberough, hx2_Re)

    # Localized losses coefficients
    K_isc_elb_per_leg = isc_Nbends * isc_kloss_bend  # 2 gomiti per gamba
    K_enl_HX2 = (1 - isc_A / hx2_A_tot) ** 2
    K_cont_HX2 = 0.5 * (1 - isc_A / hx2_A_tot)
    hx1_k_shell = 8 * (0.227 / hx1_Re_shell ** 0.193) * (hx1_Dshell / hx1_D_eq_shell) * (hx1_Nbaffles + 1)

    # Distributed losses
    loss_dist_h_ISC = f_hot * isc_L / (isc_inpipe_diam) * 1 / (2 * isc_dens_hot * isc_A ** 2)
    loss_dist_c_ISC = f_cold * isc_L / (isc_inpipe_diam) * 1 / (2 * isc_dens_cold * isc_A ** 2)
    loss_dist_HX2 = hx2_f * hx2_av_tubelength / hx2_in_diam * 1 / (2 * isc_dens * hx2_A_tot ** 2)

    # Localized Losses
    loss_loc_elb_h_ISC = K_isc_elb_per_leg * 1 / (2 * isc_dens_hot * isc_A ** 2)
    loss_loc_elb_c_ISC = K_isc_elb_per_leg * 1 / (2 * isc_dens_cold * isc_A ** 2)
    loss_loc_enl_HX2 = K_enl_HX2 * 1 / (2 * isc_dens_hot * isc_A ** 2)
    loss_loc_cont_HX2 = K_cont_HX2 * 1 / (2 * isc_dens_cold * isc_A ** 2)
    loss_loc_shell = hx1_k_shell * 1 / (2 * isc_dens * hx1_A_shell ** 2)

    # flow mass rate
    loss_tot_ISC = loss_dist_h_ISC + loss_dist_c_ISC + loss_dist_HX2 + loss_loc_elb_h_ISC + loss_loc_elb_c_ISC + loss_loc_enl_HX2 + loss_loc_cont_HX2 + loss_loc_shell
    m_new = math.sqrt((isc_dens_cold - isc_dens_hot) * gg * isc_H / loss_tot_ISC)

    if abs(m_new - m_old) / abs(m_old) < toll:
        mguess = m_new
        break
    mguess = m_new

# --- PRESSURE DROPS ISC ---
dp_h_dist_ISC = loss_dist_h_ISC * mguess**2
dp_h_loc_ISC = loss_loc_elb_h_ISC * mguess**2
dp_h_tot_ISC = dp_h_dist_ISC + dp_h_loc_ISC
dp_c_dist_ISC = loss_dist_c_ISC * mguess**2
dp_c_loc_ISC = loss_loc_elb_c_ISC * mguess**2
dp_c_tot_ISC = dp_c_dist_ISC + dp_c_loc_ISC
dp_dist_HX2 = loss_dist_HX2 * mguess**2
dp_loc_HX2 = (loss_loc_enl_HX2 + loss_loc_cont_HX2) * mguess**2
dp_tot_HX2 = dp_dist_HX2 + dp_loc_HX2
dp_shell_tot_HX1 = loss_loc_shell * mguess**2
dp_tot_ISC = dp_h_tot_ISC + dp_c_tot_ISC + dp_tot_HX2 + dp_shell_tot_HX1

print('\n======================================================')
print('                     ISC Results')
print('======================================================\n')
print('--- Thermal-Hydraulic Parameters ---')
print(f'Mass flow rate:                {mguess:10.2f}  [kg/s]')
print(f'Hot leg temperature:           {isc_Thot:10.2f}  [°C]')
print(f'Cold leg temperature:          {isc_Tcold:10.2f}  [°C]\n')
print('--- Pressure Drops [Pa] ---')
print(f'Hot leg total:                 {dp_h_tot_ISC:10.2f}')
print(f'   |- Distributed:             {dp_h_dist_ISC:10.2f}')
print(f'   |- Localized:               {dp_h_loc_ISC:10.2f}\n')
print(f'Cold leg total:                {dp_c_tot_ISC:10.2f}')
print(f'   |- Distributed:             {dp_c_dist_ISC:10.2f}')
print(f'   |- Localized:               {dp_c_loc_ISC:10.2f}\n')
print(f'HX2 total:                     {dp_tot_HX2:10.2f}')
print(f'   |- Distributed:             {dp_dist_HX2:10.2f}')
print(f'   |- Localized:               {dp_loc_HX2:10.2f}\n')
print(f'HX1 shell side total:          {dp_shell_tot_HX1:10.2f}')
print('------------------------------------------------------')
print(f'TOTAL ISC PRESSURE DROP:       {dp_tot_ISC:10.2f}')
print('======================================================\n')


# =========================================================================
# PSC SOLUTION
# =========================================================================
psc_mguess = 80.0
hx1_Tav_guess = 120.0

while True:
    pm_old = psc_mguess
    while True:
        pT_old = hx1_Tav_guess
        pres = estrai_parametri(hx1_Tav_guess, psc_press)
        psc_dens, psc_visc, psc_PR, psc_cp, psc_cond = pres

        # Delta PSC with power balance
        psc_deltat = qq / (psc_mguess * psc_cp)

        # h int and h out
        # Adimensional numbers HX1
        hx1_Re = (4 * psc_mguess) / (psc_visc * math.pi * hx1_in_diam * hx1_Ntubes)
        hx1_Nu = 0.023 * (hx1_Re ** 0.8) * (psc_PR ** 0.3)
        hx1_hint = hx1_Nu * psc_cond / hx1_in_diam

        pr_shell_hx1 = isc_visc * isc_cp / isc_cond
        hx1_h_shell = 0.351 * (hx1_Re_shell ** 0.55) * isc_cond * (pr_shell_hx1 ** (1 / 3)) / hx1_D_eq_shell

        # Global HT HX1
        Rp_tot_hx1 = (1 / (hx1_hint * hx1_A_int_lat) +
                      math.log(hx1_D_out / hx1_in_diam) / (2 * math.pi * hx1_kth * hx1_Ntubes * hx1_av_tubelength) +
                      1 / (hx1_h_shell * hx1_A_out_lat))
        hx1_U = 1 / (Rp_tot_hx1 * hx1_A_out_lat)

        hx1_deltaT_ml = qq / (hx1_heatrarea * hx1_U * hx1_Ft)


        ee = math.exp((psc_deltat - deltat_isc) / hx1_deltaT_ml)
        psc_Tcold = (ee * isc_Tcold + psc_deltat - isc_Thot) / (ee - 1)
        psc_Thot = psc_Tcold + psc_deltat
        hx1_Tav_new = (psc_Thot + psc_Tcold) / 2

        if abs(hx1_Tav_new - pT_old) / abs(pT_old) < toll:
            hx1_Tav_guess = hx1_Tav_new
            break
        hx1_Tav_guess = hx1_Tav_new

    # PSC properties
    pres_h = estrai_parametri(psc_Thot, psc_press)
    pres_c = estrai_parametri(psc_Tcold, psc_press)
    psc_dens_hot, psc_visc_hot = pres_h[0], pres_h[1]
    psc_dens_cold, psc_visc_cold = pres_c[0], pres_c[1]
    psc_dens_core = (psc_dens_hot + psc_dens_cold) / 2

    # Reynolds numbers
    Re_h_psc = (4 * psc_mguess) / (psc_visc_hot * math.pi * psc_inpipe_diam)
    Re_c_psc = (4 * psc_mguess) / (psc_visc_cold * math.pi * psc_inpipe_diam)

    # Friction factors
    f_hot_p = f_haaland(psc_pipe_rough, (4 * psc_mguess) / (psc_visc_hot * math.pi * psc_inpipe_diam))
    f_cold_p = f_haaland(psc_pipe_rough, (4 * psc_mguess) / (psc_visc_cold * math.pi * psc_inpipe_diam))
    hx1_f = f_haaland(hx1_tube_rough, hx1_Re)

    # Localized losses coefficients
    K_psc_elb_per_leg = psc_Nbends * psc_kloss_bend  # 2 gomiti per gamba
    K_enl1_HX1 = (1 - psc_A / hx1_aveflowarea) ** 2
    K_cont1_HX1 = 0.5 * (1 - hx1_A_tot / hx1_aveflowarea)
    K_hx1_elb = hx1_Nbends * hx1_kloss_bend
    K_enl2_HX1 = (1 - hx1_A_tot / hx1_aveflowarea) ** 2
    K_cont2_HX1 = 0.5 * (1 - psc_A / hx1_aveflowarea)

    # Distributed losses
    loss_dist_h_PSC = f_hot_p * psc_L / (psc_inpipe_diam) * 1 / (2 * psc_dens_hot * psc_A ** 2)
    loss_dist_c_PSC = f_cold_p * psc_L / (psc_inpipe_diam) * 1 / (2 * psc_dens_cold * psc_A ** 2)
    loss_dist_HX1 = hx1_f * hx1_av_tubelength / hx1_in_diam * 1 / (2 * psc_dens * hx1_A_tot ** 2)

    # Localized losses
    #loss_loc_elb_h_PSC = K_PSC_elb * 1 / (2 * psc_dens_hot * psc_A ** 2)
    #loss_loc_elb_c_PSC = K_PSC_elb * 1 / (2 * psc_dens_cold * psc_A ** 2)
    loss_loc_elb_h_PSC = K_psc_elb_per_leg * 1 / (2 * psc_dens_hot * psc_A ** 2)
    loss_loc_elb_c_PSC = K_psc_elb_per_leg * 1 / (2 * psc_dens_cold * psc_A ** 2)
    loss_loc_elb_HX1 = K_hx1_elb * 1 / (2 * psc_dens * hx1_A_tot ** 2)
    loss_loc_enl1_HX1 = K_enl1_HX1 * 1 / (2 * psc_dens_hot * psc_A ** 2)
    loss_loc_cont1_HX1 = K_cont1_HX1 * 1 / (2 * psc_dens_hot * hx1_A_tot ** 2)
    loss_loc_enl2_HX1 = K_enl2_HX1 * 1 / (2 * psc_dens_cold * hx1_A_tot ** 2)
    loss_loc_cont2_HX1 = K_cont2_HX1 * 1 / (2 * psc_dens_cold * psc_A ** 2)
    loss_loc_valve_PSC = psc_Klossvalve * 1 / (2 * psc_dens_cold * psc_A ** 2)

    # flow mass rate
    loss_tot_PSC = loss_dist_h_PSC + loss_dist_c_PSC + loss_dist_HX1 + loss_loc_elb_h_PSC + loss_loc_elb_c_PSC + loss_loc_elb_HX1 + loss_loc_enl1_HX1 + loss_loc_cont1_HX1 + loss_loc_enl2_HX1 + loss_loc_cont2_HX1 + loss_loc_valve_PSC + R_vessel
    psc_m_new = math.sqrt((psc_dens_cold - psc_dens_hot) * gg * (psc_H_hotcold + psc_H_topbot / 2) / loss_tot_PSC)

    if abs(psc_m_new - pm_old) / abs(pm_old) < toll:
        psc_mguess = psc_m_new
        break
    psc_mguess = psc_m_new

# --- PRESSURE DROPS PSC ---
dp_h_dist_PSC = loss_dist_h_PSC * psc_mguess**2
dp_h_loc_PSC = loss_loc_elb_h_PSC * psc_mguess**2
dp_h_tot_PSC = dp_h_dist_PSC + dp_h_loc_PSC
dp_c_dist_PSC = loss_dist_c_PSC * psc_mguess**2
dp_c_loc_PSC = (loss_loc_elb_c_PSC + loss_loc_valve_PSC) * psc_mguess**2
dp_c_tot_PSC = dp_c_dist_PSC + dp_c_loc_PSC
dp_dist_HX1 = loss_dist_HX1 * psc_mguess**2
dp_loc_HX1 = (loss_loc_enl1_HX1 + loss_loc_cont1_HX1 + loss_loc_enl2_HX1 + loss_loc_cont2_HX1 + loss_loc_elb_HX1) * psc_mguess**2
dp_tot_HX1 = dp_dist_HX1 + dp_loc_HX1
dp_pressure_vessel = R_vessel * psc_mguess**2
dp_tot_PSC = dp_h_tot_PSC + dp_c_tot_PSC + dp_tot_HX1 + dp_pressure_vessel

print('\n======================================================')
print('                     PSC Results')
print('======================================================\n')
print('--- Thermal-Hydraulic Parameters ---')
print(f'Mass flow rate:                {psc_mguess:10.2f}  [kg/s]')
print(f'Hot leg temperature:           {psc_Thot:10.2f}  [°C]')
print(f'Cold leg temperature:          {psc_Tcold:10.2f}  [°C]\n')
print('--- Pressure Drops [Pa] ---')
print(f'Hot leg total:                 {dp_h_tot_PSC:10.2f}')
print(f'   |- Distributed:             {dp_h_dist_PSC:10.2f}')
print(f'   |- Localized:               {dp_h_loc_PSC:10.2f}\n')
print(f'Cold leg total:                {dp_c_tot_PSC:10.2f}')
print(f'   |- Distributed:             {dp_c_dist_PSC:10.2f}')
print(f'   |- Localized:               {dp_c_loc_PSC:10.2f}\n')
print(f'HX1 tube side total:           {dp_tot_HX1:10.2f}')
print(f'   |- Distributed:             {dp_dist_HX1:10.2f}')
print(f'   |- Localized:               {dp_loc_HX1:10.2f}\n')
print(f'Pressure Vessel total:         {dp_pressure_vessel:10.2f}')
print('------------------------------------------------------')
print(f'TOTAL PSC PRESSURE DROP:       {dp_tot_PSC:10.2f}')
print('======================================================\n')