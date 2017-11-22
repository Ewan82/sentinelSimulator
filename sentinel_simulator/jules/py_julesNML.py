# !/usr/bin/python
"""This module holds JULES namelist files. 
It has been automatically generated.
"""

from string import Template


class julesNML:
    """This is the base class for storing
    and writing JULES namelist files

    """

    def __init__(self, template, filename):
        self.t = Template(template)
        self.filename = filename
        self.mapping = {}

    def update(self, template):
        self.t = Template(template)

    def write(self):
        f = open(self.filename, "w")
        print >> f, self.t.safe_substitute(self.mapping)
        f.close()


triffid_params_txt = """
&jules_triffid
alloc_fast_io= ${jules_triffid_1_alloc_fast_io},
alloc_med_io= ${jules_triffid_1_alloc_med_io},
alloc_slow_io= ${jules_triffid_1_alloc_slow_io},
crop_io= ${jules_triffid_1_crop_io},
dpm_rpm_ratio_io= ${jules_triffid_1_dpm_rpm_ratio_io},
g_area_io= ${jules_triffid_1_g_area_io},
g_grow_io= ${jules_triffid_1_g_grow_io},
g_root_io= ${jules_triffid_1_g_root_io},
g_wood_io= ${jules_triffid_1_g_wood_io},
lai_max_io= ${jules_triffid_1_lai_max_io},
lai_min_io= ${jules_triffid_1_lai_min_io},
retran_l_io= ${jules_triffid_1_retran_l_io},
retran_r_io= ${jules_triffid_1_retran_r_io},
/

"""
triffid_params_nml = julesNML(triffid_params_txt, "triffid_params.nml")
triffid_params_nml.mapping["jules_triffid_1_lai_max_io"] = "9.00,5.00,4.00,4.00,3.00,"
triffid_params_nml.mapping["jules_triffid_1_g_grow_io"] = "15.00,20.00,20.00,20.00,20.00,"
triffid_params_nml.mapping["jules_triffid_1_dpm_rpm_ratio_io"] = "0.25,0.25,0.67,0.67,0.33,"
triffid_params_nml.mapping["jules_triffid_1_alloc_slow_io"] = "9*0.0,"
triffid_params_nml.mapping["jules_triffid_1_retran_r_io"] = "9*0.2,"
triffid_params_nml.mapping["jules_triffid_1_alloc_med_io"] = "9*0.0,"
triffid_params_nml.mapping["jules_triffid_1_alloc_fast_io"] = "9*0.0,"
triffid_params_nml.mapping["jules_triffid_1_g_area_io"] = "0.005,0.007,0.20,0.20,0.05,"
triffid_params_nml.mapping["jules_triffid_1_retran_l_io"] = "9*0.5,"
triffid_params_nml.mapping["jules_triffid_1_crop_io"] = "0,0,1,1,0,"
triffid_params_nml.mapping["jules_triffid_1_g_wood_io"] = "0.005,0.005,0.20,0.20,0.05,"
triffid_params_nml.mapping["jules_triffid_1_g_root_io"] = "0.25,0.15,0.25,0.25,0.25,"
triffid_params_nml.mapping["jules_triffid_1_lai_min_io"] = "5*1.00,"

jules_surface_txt = """
&jules_surface
beta1= ${jules_surface_1_beta1},
beta2= ${jules_surface_1_beta2},
fwe_c3= ${jules_surface_1_fwe_c3},
fwe_c4= ${jules_surface_1_fwe_c4},
hleaf= ${jules_surface_1_hleaf},
hwood= ${jules_surface_1_hwood},
iscrntdiag= ${jules_surface_1_iscrntdiag},
l_aggregate= ${jules_surface_1_l_aggregate},
l_anthrop_heat_src= ${jules_surface_1_l_anthrop_heat_src},
l_elev_land_ice= ${jules_surface_1_l_elev_land_ice},
l_elev_lw_down= ${jules_surface_1_l_elev_lw_down},
l_epot_corr= ${jules_surface_1_l_epot_corr},
l_land_ice_imp= ${jules_surface_1_l_land_ice_imp},
l_point_data= ${jules_surface_1_l_point_data},
/

"""
jules_surface_nml = julesNML(jules_surface_txt, "jules_surface.nml")
jules_surface_nml.mapping["jules_surface_1_l_epot_corr"] = ".false.,"
jules_surface_nml.mapping["jules_surface_1_l_elev_land_ice"] = ".false.,"
jules_surface_nml.mapping["jules_surface_1_hleaf"] = "5.7e4,"
jules_surface_nml.mapping["jules_surface_1_l_point_data"] = ".true.,"
jules_surface_nml.mapping["jules_surface_1_iscrntdiag"] = "0,"
jules_surface_nml.mapping["jules_surface_1_l_land_ice_imp"] = ".false.,"
jules_surface_nml.mapping["jules_surface_1_fwe_c3"] = "0.5,"
jules_surface_nml.mapping["jules_surface_1_hwood"] = "1.1e4,"
jules_surface_nml.mapping["jules_surface_1_l_elev_lw_down"] = ".false.,"
jules_surface_nml.mapping["jules_surface_1_l_aggregate"] = ".false.,"
jules_surface_nml.mapping["jules_surface_1_fwe_c4"] = "2.0e4,"
jules_surface_nml.mapping["jules_surface_1_beta2"] = "0.93,"
jules_surface_nml.mapping["jules_surface_1_beta1"] = "0.83,"
jules_surface_nml.mapping["jules_surface_1_l_anthrop_heat_src"] = ".false.,"

pft_params_txt = """
&jules_pftparm
a_wl_io= ${jules_pftparm_1_a_wl_io},
a_ws_io= ${jules_pftparm_1_a_ws_io},
aef_io= ${jules_pftparm_1_aef_io},
albsnc_max_io= ${jules_pftparm_1_albsnc_max_io},
albsnc_min_io= ${jules_pftparm_1_albsnc_min_io},
albsnf_max_io= ${jules_pftparm_1_albsnf_max_io},
albsnf_maxl_io= ${jules_pftparm_1_albsnf_maxl_io},
albsnf_maxu_io= ${jules_pftparm_1_albsnf_maxu_io},
alnir_io= ${jules_pftparm_1_alnir_io},
alnirl_io= ${jules_pftparm_1_alnirl_io},
alniru_io= ${jules_pftparm_1_alniru_io},
alpar_io= ${jules_pftparm_1_alpar_io},
alparl_io= ${jules_pftparm_1_alparl_io},
alparu_io= ${jules_pftparm_1_alparu_io},
alpha_io= ${jules_pftparm_1_alpha_io},
avg_ba_io= ${jules_pftparm_1_avg_ba_io},
b_wl_io= ${jules_pftparm_1_b_wl_io},
c3_io= ${jules_pftparm_1_c3_io},
can_struct_a_io= ${jules_pftparm_1_can_struct_a_io},
canht_ft_io= ${jules_pftparm_1_canht_ft_io},
catch0_io= ${jules_pftparm_1_catch0_io},
ccleaf_max_io= ${jules_pftparm_1_ccleaf_max_io},
ccleaf_min_io= ${jules_pftparm_1_ccleaf_min_io},
ccwood_max_io= ${jules_pftparm_1_ccwood_max_io},
ccwood_min_io= ${jules_pftparm_1_ccwood_min_io},
ci_st_io= ${jules_pftparm_1_ci_st_io},
dcatch_dlai_io= ${jules_pftparm_1_dcatch_dlai_io},
dfp_dcuo_io= ${jules_pftparm_1_dfp_dcuo_io},
dgl_dm_io= ${jules_pftparm_1_dgl_dm_io},
dgl_dt_io= ${jules_pftparm_1_dgl_dt_io},
dqcrit_io= ${jules_pftparm_1_dqcrit_io},
dz0v_dh_io= ${jules_pftparm_1_dz0v_dh_io},
emis_pft_io= ${jules_pftparm_1_emis_pft_io},
eta_sl_io= ${jules_pftparm_1_eta_sl_io},
f0_io= ${jules_pftparm_1_f0_io},
fd_io= ${jules_pftparm_1_fd_io},
fef_bc_io= ${jules_pftparm_1_fef_bc_io},
fef_ch4_io= ${jules_pftparm_1_fef_ch4_io},
fef_co2_io= ${jules_pftparm_1_fef_co2_io},
fef_co_io= ${jules_pftparm_1_fef_co_io},
fef_nox_io= ${jules_pftparm_1_fef_nox_io},
fef_oc_io= ${jules_pftparm_1_fef_oc_io},
fef_so2_io= ${jules_pftparm_1_fef_so2_io},
fl_o3_ct_io= ${jules_pftparm_1_fl_o3_ct_io},
fsmc_mod_io= ${jules_pftparm_1_fsmc_mod_io},
fsmc_of_io= ${jules_pftparm_1_fsmc_of_io},
fsmc_p0_io= ${jules_pftparm_1_fsmc_p0_io},
g_leaf_0_io= ${jules_pftparm_1_g_leaf_0_io},
glmin_io= ${jules_pftparm_1_glmin_io},
gpp_st_io= ${jules_pftparm_1_gpp_st_io},
hw_sw_io= ${jules_pftparm_1_hw_sw_io},
ief_io= ${jules_pftparm_1_ief_io},
infil_f_io= ${jules_pftparm_1_infil_f_io},
kext_io= ${jules_pftparm_1_kext_io},
kn_io= ${jules_pftparm_1_kn_io},
knl_io= ${jules_pftparm_1_knl_io},
kpar_io= ${jules_pftparm_1_kpar_io},
lai_alb_lim_io= ${jules_pftparm_1_lai_alb_lim_io},
lai_io= ${jules_pftparm_1_lai_io},
lma_io= ${jules_pftparm_1_lma_io},
mef_io= ${jules_pftparm_1_mef_io},
neff_io= ${jules_pftparm_1_neff_io},
nl0_io= ${jules_pftparm_1_nl0_io},
nmass_io= ${jules_pftparm_1_nmass_io},
nr_io= ${jules_pftparm_1_nr_io},
nr_nl_io= ${jules_pftparm_1_nr_nl_io},
ns_nl_io= ${jules_pftparm_1_ns_nl_io},
nsw_io= ${jules_pftparm_1_nsw_io},
omega_io= ${jules_pftparm_1_omega_io},
omegal_io= ${jules_pftparm_1_omegal_io},
omegau_io= ${jules_pftparm_1_omegau_io},
omnir_io= ${jules_pftparm_1_omnir_io},
omnirl_io= ${jules_pftparm_1_omnirl_io},
omniru_io= ${jules_pftparm_1_omniru_io},
orient_io= ${jules_pftparm_1_orient_io},
q10_leaf_io= ${jules_pftparm_1_q10_leaf_io},
r_grow_io= ${jules_pftparm_1_r_grow_io},
rootd_ft_io= ${jules_pftparm_1_rootd_ft_io},
sigl_io= ${jules_pftparm_1_sigl_io},
tef_io= ${jules_pftparm_1_tef_io},
tleaf_of_io= ${jules_pftparm_1_tleaf_of_io},
tlow_io= ${jules_pftparm_1_tlow_io},
tupp_io= ${jules_pftparm_1_tupp_io},
vint_io= ${jules_pftparm_1_vint_io},
vsl_io= ${jules_pftparm_1_vsl_io},
z0hm_classic_pft_io= ${jules_pftparm_1_z0hm_classic_pft_io},
z0hm_pft_io= ${jules_pftparm_1_z0hm_pft_io},
/

"""
pft_params_nml = julesNML(pft_params_txt, "pft_params.nml")
pft_params_nml.mapping["jules_pftparm_1_omegau_io"] = "0.23,0.23,0.23,0.26,0.23,0.23,0.23,0.26,0.23,"
pft_params_nml.mapping["jules_pftparm_1_dqcrit_io"] = "0.09,0.06,0.1,0.075,0.1,0.1,0.1,0.075,0.1,"
pft_params_nml.mapping["jules_pftparm_1_albsnf_maxu_io"] = "0.215,0.132,0.288,0.239,0.173,0.288,0.288,0.239,0.288,"
pft_params_nml.mapping["jules_pftparm_1_albsnc_min_io"] = "0.3,0.3,7*0.8,"
pft_params_nml.mapping["jules_pftparm_1_ccwood_min_io"] = "9*0.0,"
pft_params_nml.mapping["jules_pftparm_1_kpar_io"] = "9*0.5,"
pft_params_nml.mapping["jules_pftparm_1_fd_io"] = "0.01,0.015,0.015,0.025,0.015,0.0095,0.015,0.025,0.015,"
pft_params_nml.mapping["jules_pftparm_1_dfp_dcuo_io"] = "0.04,0.02,0.25,0.13,0.03,0.25,0.25,0.13,0.25,"
pft_params_nml.mapping["jules_pftparm_1_kext_io"] = "9*0.5,"
pft_params_nml.mapping["jules_pftparm_1_nr_io"] = "0.01726,0.00784,0.0162,0.0084,0.01726,-1,-1,-1,-1,"
pft_params_nml.mapping["jules_pftparm_1_alnirl_io"] = "0.30,0.23,7*0.39,"
pft_params_nml.mapping["jules_pftparm_1_ccwood_max_io"] = "9*0.4,"
pft_params_nml.mapping["jules_pftparm_1_fsmc_p0_io"] = "9*0.65,"
pft_params_nml.mapping["jules_pftparm_1_alniru_io"] = "0.68,0.53,7*0.87,"
pft_params_nml.mapping["jules_pftparm_1_lai_io"] = "5,4,2,4,1,2,2,4,2,"
pft_params_nml.mapping["jules_pftparm_1_vsl_io"] = "29.81,18.15,40.96,10.24,23.15,"
pft_params_nml.mapping["jules_pftparm_1_fef_oc_io"] = "4.3,9.1,9.1,3.2,9.1,9.1,9.1,3.2,9.1,"
pft_params_nml.mapping["jules_pftparm_1_fsmc_of_io"] = "9*0,"
pft_params_nml.mapping["jules_pftparm_1_omnirl_io"] = "0.50,0.30,7*0.53,"
pft_params_nml.mapping["jules_pftparm_1_dcatch_dlai_io"] = "9*0.05,"
pft_params_nml.mapping["jules_pftparm_1_glmin_io"] = "9*1.00e-6,"
pft_params_nml.mapping["jules_pftparm_1_q10_leaf_io"] = "9*2.00,"
pft_params_nml.mapping["jules_pftparm_1_fef_co_io"] = "100,106,106,64,106,106,106,64,106,"
pft_params_nml.mapping["jules_pftparm_1_nmass_io"] = "0.0210,0.0115,0.0219,0.0131,0.0219,"
pft_params_nml.mapping["jules_pftparm_1_catch0_io"] = "9*0.5,"
pft_params_nml.mapping["jules_pftparm_1_emis_pft_io"] = "9*1.00,"
pft_params_nml.mapping["jules_pftparm_1_a_ws_io"] = "10,10,1,1,10,1,1,1,1,"
pft_params_nml.mapping["jules_pftparm_1_a_wl_io"] = "0.65,0.65,0.005,0.005,0.1,0.0095,0.005,0.005,0.005,"
pft_params_nml.mapping["jules_pftparm_1_alparu_io"] = "0.15,0.11,7*0.15,"
pft_params_nml.mapping["jules_pftparm_1_omnir_io"] = "0.7,0.45,7*0.83,"
pft_params_nml.mapping["jules_pftparm_1_fsmc_mod_io"] = "9*1,"
pft_params_nml.mapping["jules_pftparm_1_dz0v_dh_io"] = "0.05,0.05,7*0.1,"
pft_params_nml.mapping["jules_pftparm_1_sigl_io"] = "0.0375,0.1,0.025,0.05,0.05,0.032,0.032,0.025,0.032,"
pft_params_nml.mapping["jules_pftparm_1_ccleaf_max_io"] = "9*1.0,"
pft_params_nml.mapping[
    "jules_pftparm_1_gpp_st_io"] = "1.29e-7,2.58e-8,2.07e-7,3.42e-7,1.68e-7,2.07e-7,2.07e-7,3.42e-7,2.07e-7,"
pft_params_nml.mapping["jules_pftparm_1_dgl_dt_io"] = "9*9,"
pft_params_nml.mapping["jules_pftparm_1_albsnf_max_io"] = "0.1,0.1,7*0.2,"
pft_params_nml.mapping[
    "jules_pftparm_1_tleaf_of_io"] = "278.15,233.15,278.15,278.15,233.15,200.15,278.15,278.15,278.15,"
pft_params_nml.mapping["jules_pftparm_1_tlow_io"] = "0,-10,0,13,0,0,0,13,0,"
pft_params_nml.mapping["jules_pftparm_1_tupp_io"] = "36,26,36,45,36,46,36,45,36,"
pft_params_nml.mapping["jules_pftparm_1_fef_nox_io"] = "2.55,3.24,3.24,2.49,3.24,3.24,3.24,2.49,3.24,"
pft_params_nml.mapping["jules_pftparm_1_kn_io"] = "9*0.78,"
pft_params_nml.mapping["jules_pftparm_1_fef_so2_io"] = "0.40,0.40,0.40,0.48,0.40,0.40,0.40,0.48,0.40,"
pft_params_nml.mapping["jules_pftparm_1_tef_io"] = "0.40,2.40,7*0.80,"
pft_params_nml.mapping["jules_pftparm_1_fl_o3_ct_io"] = "1.6,1.6,5.0,5.0,1.6,5.0,5.0,5.0,5.0,"
pft_params_nml.mapping["jules_pftparm_1_dgl_dm_io"] = "9*0,"
pft_params_nml.mapping["jules_pftparm_1_fef_co2_io"] = "1631,1576,1576,1654,1576,1576,1576,1654,1576,"
pft_params_nml.mapping["jules_pftparm_1_can_struct_a_io"] = "9*0.65,"
pft_params_nml.mapping["jules_pftparm_1_canht_ft_io"] = "19.01,16.38,0.79,1.26,1,0.5,0.79,1.26,0.79,"
pft_params_nml.mapping["jules_pftparm_1_albsnc_max_io"] = "0.15,0.15,0.6,0.6,0.4,0.6,0.6,0.6,0.6,"
pft_params_nml.mapping["jules_pftparm_1_ief_io"] = "35.0,12.0,16.0,8.0,20.0,16.0,16.0,8.0,16.0,"
pft_params_nml.mapping["jules_pftparm_1_c3_io"] = "1,1,1,0,1,0,1,0,1,"
pft_params_nml.mapping["jules_pftparm_1_avg_ba_io"] = "0.6e6,0.6e6,1.4e6,1.4e6,1.2e6,1.4e6,1.4e6,1.4e6,1.4e6,"
pft_params_nml.mapping["jules_pftparm_1_b_wl_io"] = "9*1.667,"
pft_params_nml.mapping["jules_pftparm_1_vint_io"] = "5.73,6.32,6.42,0.00,14.71,"
pft_params_nml.mapping["jules_pftparm_1_omniru_io"] = "0.90,0.65,7*0.98,"
pft_params_nml.mapping[
    "jules_pftparm_1_neff_io"] = "8.00e-4,8.00e-4,8.00e-4,4.00e-4,8.00e-4,5.70e-4,8.00e-4,4.00e-4,8.00e-4,"
pft_params_nml.mapping["jules_pftparm_1_f0_io"] = "0.875,0.875,0.9,0.8,0.9,0.4,0.9,0.8,0.9,"
pft_params_nml.mapping["jules_pftparm_1_mef_io"] = "0.60,0.90,0.60,0.90,0.57,0.60,0.60,0.90,0.60,"
pft_params_nml.mapping["jules_pftparm_1_infil_f_io"] = "4,4,7*2,"
pft_params_nml.mapping["jules_pftparm_1_r_grow_io"] = "9*0.25,"
pft_params_nml.mapping["jules_pftparm_1_g_leaf_0_io"] = "9*0.25,"
pft_params_nml.mapping["jules_pftparm_1_hw_sw_io"] = "9*0.5,"
pft_params_nml.mapping["jules_pftparm_1_rootd_ft_io"] = "3,1,7*1.7,"
pft_params_nml.mapping["jules_pftparm_1_nl0_io"] = "0.046,0.033,0.073,0.06,0.06,0.073,0.073,0.06,0.073,"
pft_params_nml.mapping["jules_pftparm_1_fef_ch4_io"] = "6.8,4.8,4.8,2.4,4.8,4.8,4.8,2.4,4.8,"
pft_params_nml.mapping["jules_pftparm_1_omega_io"] = "0.15,0.15,0.15,0.17,0.15,0.15,0.15,0.17,0.15,"
pft_params_nml.mapping["jules_pftparm_1_aef_io"] = "0.18,0.21,0.12,0.08,0.20,0.12,0.12,0.08,0.12,"
pft_params_nml.mapping["jules_pftparm_1_nr_nl_io"] = "9*0.39,"
pft_params_nml.mapping["jules_pftparm_1_nsw_io"] = "0.0072,0.0083,0.01604,0.0202,0.0072,-1,-1,-1,-1,"
pft_params_nml.mapping["jules_pftparm_1_ccleaf_min_io"] = "9*0.8,"
pft_params_nml.mapping["jules_pftparm_1_knl_io"] = "9*0.,"
pft_params_nml.mapping["jules_pftparm_1_alparl_io"] = "0.06,0.04,7*0.06,"
pft_params_nml.mapping["jules_pftparm_1_omegal_io"] = "0.10,0.10,0.10,0.12,0.10,0.10,0.10,0.12,0.10,"
pft_params_nml.mapping["jules_pftparm_1_fef_bc_io"] = "0.56,0.56,0.56,0.47,0.56,0.56,0.56,0.47,0.56,"
pft_params_nml.mapping["jules_pftparm_1_eta_sl_io"] = "9*0.01,"
pft_params_nml.mapping["jules_pftparm_1_lai_alb_lim_io"] = "9*0.5,"
pft_params_nml.mapping["jules_pftparm_1_albsnf_maxl_io"] = "0.095,0.059,0.128,0.106,0.077,0.128,0.128,0.106,0.128,"
pft_params_nml.mapping["jules_pftparm_1_lma_io"] = "0.0824,0.2263,0.0498,0.1370,0.0695,"
pft_params_nml.mapping["jules_pftparm_1_orient_io"] = "9*0,"
pft_params_nml.mapping["jules_pftparm_1_alpar_io"] = "0.1,0.07,7*0.1,"
pft_params_nml.mapping["jules_pftparm_1_alnir_io"] = "0.45,0.35,7*0.58,"
pft_params_nml.mapping["jules_pftparm_1_ci_st_io"] = "33.46,33.46,34.26,29.98,34.26,34.26,34.26,29.98,34.26,"
pft_params_nml.mapping["jules_pftparm_1_alpha_io"] = "0.08,0.08,0.12,0.06,0.08,0.12,0.12,0.06,0.12,"
pft_params_nml.mapping["jules_pftparm_1_z0hm_classic_pft_io"] = "9*0.1,"
pft_params_nml.mapping["jules_pftparm_1_ns_nl_io"] = "0.1,0.1,1,1,0.1,0.43,1,1,1,"
pft_params_nml.mapping["jules_pftparm_1_z0hm_pft_io"] = "9*0.1,"

output_txt = """
&jules_output
dump_period= ${jules_output_1_dump_period},
nprofiles= ${jules_output_1_nprofiles},
output_dir= ${jules_output_1_output_dir},
run_id= ${jules_output_1_run_id},
/

&jules_output_profile
file_period= ${jules_output_profile_1_file_period},
nvars= ${jules_output_profile_1_nvars},
output_main_run= ${jules_output_profile_1_output_main_run},
output_period= ${jules_output_profile_1_output_period},
output_spinup= ${jules_output_profile_1_output_spinup},
output_type= ${jules_output_profile_1_output_type},
profile_name= ${jules_output_profile_1_profile_name},
var= ${jules_output_profile_1_var},
var_name= ${jules_output_profile_1_var_name},
/

"""
output_nml = julesNML(output_txt, "output.nml")
output_nml.mapping["jules_output_profile_1_nvars"] = "20,"
output_nml.mapping["jules_output_profile_1_profile_name"] = "'3_hourly',"
output_nml.mapping["jules_output_profile_1_output_spinup"] = ".false.,"
output_nml.mapping["jules_output_profile_1_output_main_run"] = ".true.,"
output_nml.mapping["jules_output_profile_1_var"] = "'rainfall', 't1p5m', 'smcl', 't_soil', 'tstar', 'cropdvi', 'croprootc', " \
                                    "'cropharvc', 'cropreservec', 'croplai', 'cropcanht', 'cropyield', 'harvest_counter'," \
                                    " 'harvest_trigger', 'cropstemc', 'cropleafc', 'cropsowdate', 'fsmc', 'fsmc_gb', 'gpp_gb'"
output_nml.mapping["jules_output_profile_1_file_period"] = "0,"
output_nml.mapping["jules_output_1_run_id"] = "'wallerfing_79_12',"
output_nml.mapping["jules_output_profile_1_var_name"] = "'rainfall', 't1p5m', 'smcl', 't_soil', 'tstar', 'cropdvi', 'croprootc', " \
                                    "'cropharvc', 'cropreservec', 'croplai', 'cropcanht', 'cropyield', 'harvest_counter'," \
                                    " 'harvest_trigger', 'cropstemc', 'cropleafc', 'cropsowdate', 'fsmc', 'fsmc_gb', 'gpp_gb'"
output_nml.mapping["jules_output_profile_1_output_period"] = "10800,"
output_nml.mapping["jules_output_profile_1_output_type"] = "20*'M',"
output_nml.mapping["jules_output_1_dump_period"] = "100,"
output_nml.mapping["jules_output_1_output_dir"] = "'./output',"
output_nml.mapping["jules_output_1_nprofiles"] = "1,"

nveg_params_txt = """
&jules_nvegparm
albsnc_nvg_io= ${jules_nvegparm_1_albsnc_nvg_io},
albsnf_nvg_io= ${jules_nvegparm_1_albsnf_nvg_io},
albsnf_nvgl_io= ${jules_nvegparm_1_albsnf_nvgl_io},
albsnf_nvgu_io= ${jules_nvegparm_1_albsnf_nvgu_io},
catch_nvg_io= ${jules_nvegparm_1_catch_nvg_io},
ch_nvg_io= ${jules_nvegparm_1_ch_nvg_io},
emis_nvg_io= ${jules_nvegparm_1_emis_nvg_io},
gs_nvg_io= ${jules_nvegparm_1_gs_nvg_io},
infil_nvg_io= ${jules_nvegparm_1_infil_nvg_io},
vf_nvg_io= ${jules_nvegparm_1_vf_nvg_io},
z0_nvg_io= ${jules_nvegparm_1_z0_nvg_io},
z0hm_classic_nvg_io= ${jules_nvegparm_1_z0hm_classic_nvg_io},
z0hm_nvg_io= ${jules_nvegparm_1_z0hm_nvg_io},
/

"""
nveg_params_nml = julesNML(nveg_params_txt, "nveg_params.nml")
nveg_params_nml.mapping["jules_nvegparm_1_albsnf_nvg_io"] = "0.18,0.06,-1.00,0.75,"
nveg_params_nml.mapping["jules_nvegparm_1_infil_nvg_io"] = "0.10,0.00,0.50,0.00,"
nveg_params_nml.mapping["jules_nvegparm_1_z0hm_classic_nvg_io"] = "0.1,0.1,0.1,0.1,"
nveg_params_nml.mapping["jules_nvegparm_1_z0hm_nvg_io"] = "0.1,0.1,0.1,0.1,"
nveg_params_nml.mapping["jules_nvegparm_1_emis_nvg_io"] = "1.00,1.00,1.00,1.00,"
nveg_params_nml.mapping["jules_nvegparm_1_catch_nvg_io"] = "0.50,0.00,0.00,0.00,"
nveg_params_nml.mapping["jules_nvegparm_1_albsnc_nvg_io"] = "0.40,0.80,0.80,0.80,"
nveg_params_nml.mapping["jules_nvegparm_1_z0_nvg_io"] = "1.00,3.0e-4,3.0e-4,1.0e-4,"
nveg_params_nml.mapping["jules_nvegparm_1_albsnf_nvgu_io"] = "0.20,0.15,0.80,0.75,"
nveg_params_nml.mapping["jules_nvegparm_1_vf_nvg_io"] = "1.00,1.00,0.00,0.00,"
nveg_params_nml.mapping["jules_nvegparm_1_albsnf_nvgl_io"] = "0.16,0.06,0.03,0.75,"
nveg_params_nml.mapping["jules_nvegparm_1_ch_nvg_io"] = "0.28e6,2.11e7,0.00,0.00,"
nveg_params_nml.mapping["jules_nvegparm_1_gs_nvg_io"] = "0.00,0.00,1.0e-2,1.0e6,"

model_grid_txt = """
&jules_input_grid
	nx = ${jules_input_grid_1_nx},
	ny = ${jules_input_grid_1_ny},
	x_dim_name = ${jules_input_grid_1_x_dim_name},
	y_dim_name = ${jules_input_grid_1_y_dim_name},
	time_dim_name = ${jules_input_grid_1_time_dim_name},
	type_dim_name = ${jules_input_grid_1_type_dim_name},
	soil_dim_name = ${jules_input_grid_1_soil_dim_name},
	pft_dim_name = ${jules_input_grid_1_pft_dim_name},
/

&jules_latlon
	file = ${jules_latlon_1_file},
	lat_name = ${jules_latlon_1_lat_name},
	lon_name = ${jules_latlon_1_lon_name},
/

&jules_land_frac
	file = ${jules_land_frac_1_file},
	land_frac_name = ${jules_land_frac_1_land_frac_name},
/

&jules_model_grid
	land_only = ${jules_model_grid_1_land_only},
	use_subgrid = ${jules_model_grid_1_use_subgrid},
	latlon_region = ${jules_model_grid_1_latlon_region},
	lat_bounds = ${jules_model_grid_1_lat_bounds},
	lon_bounds = ${jules_model_grid_1_lon_bounds},
/

&jules_surf_hgt
	zero_height= ${jules_surf_hgt_1_zero_height},
/

"""
model_grid_nml = julesNML(model_grid_txt, "model_grid.nml")
model_grid_nml.mapping["jules_latlon_1_lat_name"] = " 'latitude',"
model_grid_nml.mapping["jules_input_grid_1_time_dim_name"] = " 'tstep' , "
model_grid_nml.mapping["jules_land_frac_1_land_frac_name"] = " 'lsmask' "
model_grid_nml.mapping["jules_latlon_1_file"] = " 'WFDEI-long-lat-2d-2.nc',"
model_grid_nml.mapping["jules_input_grid_1_y_dim_name"] = " 'lat' , "
model_grid_nml.mapping["jules_model_grid_1_latlon_region"] = " T,"
model_grid_nml.mapping["jules_input_grid_1_nx"] = " 720"
model_grid_nml.mapping["jules_land_frac_1_file"] = " 'WFD-EI-LandFraction2d-2.nc',"
model_grid_nml.mapping["jules_input_grid_1_ny"] = " 360"
model_grid_nml.mapping["jules_model_grid_1_lat_bounds"] = " 48.6,48.9,"
model_grid_nml.mapping["jules_model_grid_1_land_only"] = " F,"
model_grid_nml.mapping["jules_surf_hgt_1_zero_height"] = ".true.,"
model_grid_nml.mapping["jules_model_grid_1_lon_bounds"] = " 12.6,13.0"
model_grid_nml.mapping["jules_input_grid_1_x_dim_name"] = " 'lon' , "
model_grid_nml.mapping["jules_model_grid_1_use_subgrid"] = " T,"
model_grid_nml.mapping["jules_input_grid_1_soil_dim_name"] = " 'z' ,"
model_grid_nml.mapping["jules_input_grid_1_pft_dim_name"] = " 'z' ,"
model_grid_nml.mapping["jules_input_grid_1_type_dim_name"] = " 'z' , "
model_grid_nml.mapping["jules_latlon_1_lon_name"] = " 'longitude'"

crop_params_txt = """
&jules_cropparm
allo1_io= ${jules_cropparm_1_allo1_io},
allo2_io= ${jules_cropparm_1_allo2_io},
alpha1_io= ${jules_cropparm_1_alpha1_io},
alpha2_io= ${jules_cropparm_1_alpha2_io},
alpha3_io= ${jules_cropparm_1_alpha3_io},
beta1_io= ${jules_cropparm_1_beta1_io},
beta2_io= ${jules_cropparm_1_beta2_io},
beta3_io= ${jules_cropparm_1_beta3_io},
cfrac_l_io= ${jules_cropparm_1_cfrac_l_io},
cfrac_r_io= ${jules_cropparm_1_cfrac_r_io},
cfrac_s_io= ${jules_cropparm_1_cfrac_s_io},
crit_pp_io= ${jules_cropparm_1_crit_pp_io},
delta_io= ${jules_cropparm_1_delta_io},
gamma_io= ${jules_cropparm_1_gamma_io},
initial_c_dvi_io= ${jules_cropparm_1_initial_c_dvi_io},
initial_carbon_io= ${jules_cropparm_1_initial_carbon_io},
mu_io= ${jules_cropparm_1_mu_io},
nu_io= ${jules_cropparm_1_nu_io},
pp_sens_io= ${jules_cropparm_1_pp_sens_io},
remob_io= ${jules_cropparm_1_remob_io},
rt_dir_io= ${jules_cropparm_1_rt_dir_io},
sen_dvi_io= ${jules_cropparm_1_sen_dvi_io},
t_bse_io= ${jules_cropparm_1_t_bse_io},
t_max_io= ${jules_cropparm_1_t_max_io},
t_mort_io= ${jules_cropparm_1_t_mort_io},
t_opt_io= ${jules_cropparm_1_t_opt_io},
tt_emr_io= ${jules_cropparm_1_tt_emr_io},
yield_frac_io= ${jules_cropparm_1_yield_frac_io},
/

"""
crop_params_nml = julesNML(crop_params_txt, "crop_params.nml")
crop_params_nml.mapping["jules_cropparm_1_cfrac_r_io"] = "0.5,0.5,0.5,0.5,"
crop_params_nml.mapping["jules_cropparm_1_gamma_io"] = "17.3,25.9,22.5,20.9,"
crop_params_nml.mapping["jules_cropparm_1_mu_io"] = "0.02,0.0,0.0,0.0,"
crop_params_nml.mapping["jules_cropparm_1_yield_frac_io"] = "0.74,0.0,0.0,0.0,"
crop_params_nml.mapping["jules_cropparm_1_alpha2_io"] = "12.1,18.5,12.5,19.0,"
crop_params_nml.mapping["jules_cropparm_1_initial_c_dvi_io"] = "0.1,0.0,0.0,0.0,"
crop_params_nml.mapping["jules_cropparm_1_t_mort_io"] = "273.15,0.0,0.0,0.0,"
crop_params_nml.mapping["jules_cropparm_1_t_max_io"] = "315.15,313.15,315.15,315.15,"
crop_params_nml.mapping["jules_cropparm_1_delta_io"] = "-0.3307,-0.1451,-0.2587,-0.2724,"
crop_params_nml.mapping["jules_cropparm_1_cfrac_s_io"] = "0.5,0.5,0.5,0.5,"
crop_params_nml.mapping["jules_cropparm_1_t_opt_io"] = "303.15,300.15,303.15,303.15,"
crop_params_nml.mapping["jules_cropparm_1_beta3_io"] = "-14.5,-15.0,-14.0,-18.5,"
crop_params_nml.mapping["jules_cropparm_1_nu_io"] = "4.0,0.0,0.0,0.0,"
crop_params_nml.mapping["jules_cropparm_1_alpha3_io"] = "13.0,19.5,13.0,19.5,"
crop_params_nml.mapping["jules_cropparm_1_alpha1_io"] = "13.5,20.0,13.5,18.5,"
crop_params_nml.mapping["jules_cropparm_1_sen_dvi_io"] = "1.5,0.0,0.0,0.0,"
crop_params_nml.mapping["jules_cropparm_1_beta1_io"] = "-15.0,-16.5,-15.5,-19.0,"
crop_params_nml.mapping["jules_cropparm_1_cfrac_l_io"] = "0.5,0.5,0.5,0.5,"
crop_params_nml.mapping["jules_cropparm_1_pp_sens_io"] = "0.00,0.0,0.00,0.0,"
crop_params_nml.mapping["jules_cropparm_1_remob_io"] = "0.12,0.18,0.35,0.25,"
crop_params_nml.mapping["jules_cropparm_1_rt_dir_io"] = "0.0,0.0,0.0,0.0,"
crop_params_nml.mapping["jules_cropparm_1_allo1_io"] = "3.6,1.6,3.5,1.4,"
crop_params_nml.mapping["jules_cropparm_1_initial_carbon_io"] = "0.0008,0.0,0.0,0.0,"
crop_params_nml.mapping["jules_cropparm_1_allo2_io"] = "0.4,0.4,0.4,0.4,"
crop_params_nml.mapping["jules_cropparm_1_beta2_io"] = "-12.0,-14.5,-12.5,-17.0,"
crop_params_nml.mapping["jules_cropparm_1_crit_pp_io"] = "24.0,11.0,24.0,24.0,"
crop_params_nml.mapping["jules_cropparm_1_t_bse_io"] = "281.15,278.15,281.15,281.15,"
crop_params_nml.mapping["jules_cropparm_1_tt_emr_io"] = "35.0,35.0,80.0,60.0,"

urban_txt = """
"""
urban_nml = julesNML(urban_txt, "urban.nml")

jules_vegetation_txt = """
&jules_vegetation
can_model= ${jules_vegetation_1_can_model},
can_rad_mod= ${jules_vegetation_1_can_rad_mod},
ignition_method= ${jules_vegetation_1_ignition_method},
ilayers= ${jules_vegetation_1_ilayers},
l_bvoc_emis= ${jules_vegetation_1_l_bvoc_emis},
l_gleaf_fix= ${jules_vegetation_1_l_gleaf_fix},
l_ht_compete= ${jules_vegetation_1_l_ht_compete},
l_inferno= ${jules_vegetation_1_l_inferno},
l_irrig_dmd= ${jules_vegetation_1_l_irrig_dmd},
l_landuse= ${jules_vegetation_1_l_landuse},
l_leaf_n_resp_fix= ${jules_vegetation_1_l_leaf_n_resp_fix},
l_nitrogen= ${jules_vegetation_1_l_nitrogen},
l_o3_damage= ${jules_vegetation_1_l_o3_damage},
l_phenol= ${jules_vegetation_1_l_phenol},
l_prescsow= ${jules_vegetation_1_l_prescsow},
l_scale_resp_pm= ${jules_vegetation_1_l_scale_resp_pm},
l_stem_resp_fix= ${jules_vegetation_1_l_stem_resp_fix},
l_trait_phys= ${jules_vegetation_1_l_trait_phys},
l_trif_crop= ${jules_vegetation_1_l_trif_crop},
l_trif_fire= ${jules_vegetation_1_l_trif_fire},
l_triffid= ${jules_vegetation_1_l_triffid},
l_vegcan_soilfx= ${jules_vegetation_1_l_vegcan_soilfx},
/

"""
jules_vegetation_nml = julesNML(jules_vegetation_txt, "jules_vegetation.nml")
jules_vegetation_nml.mapping["jules_vegetation_1_l_stem_resp_fix"] = ".false.,"
jules_vegetation_nml.mapping["jules_vegetation_1_can_model"] = "3,"
jules_vegetation_nml.mapping["jules_vegetation_1_ilayers"] = "10,"
jules_vegetation_nml.mapping["jules_vegetation_1_l_landuse"] = ".false.,"
jules_vegetation_nml.mapping["jules_vegetation_1_l_inferno"] = ".false.,"
jules_vegetation_nml.mapping["jules_vegetation_1_l_scale_resp_pm"] = ".true.,"
jules_vegetation_nml.mapping["jules_vegetation_1_l_nitrogen"] = ".false.,"
jules_vegetation_nml.mapping["jules_vegetation_1_l_triffid"] = ".false.,"
jules_vegetation_nml.mapping["jules_vegetation_1_l_vegcan_soilfx"] = ".false.,"
jules_vegetation_nml.mapping["jules_vegetation_1_l_leaf_n_resp_fix"] = ".false.,"
jules_vegetation_nml.mapping["jules_vegetation_1_ignition_method"] = "1,"
jules_vegetation_nml.mapping["jules_vegetation_1_l_trait_phys"] = ".false.,"
jules_vegetation_nml.mapping["jules_vegetation_1_l_gleaf_fix"] = ".false.,"
jules_vegetation_nml.mapping["jules_vegetation_1_l_trif_fire"] = ".false.,"
jules_vegetation_nml.mapping["jules_vegetation_1_l_ht_compete"] = ".false.,"
jules_vegetation_nml.mapping["jules_vegetation_1_l_bvoc_emis"] = ".false.,"
jules_vegetation_nml.mapping["jules_vegetation_1_l_trif_crop"] = ".false.,"
jules_vegetation_nml.mapping["jules_vegetation_1_can_rad_mod"] = "6,"
jules_vegetation_nml.mapping["jules_vegetation_1_l_prescsow"] = ".true.,"
jules_vegetation_nml.mapping["jules_vegetation_1_l_o3_damage"] = ".false.,"
jules_vegetation_nml.mapping["jules_vegetation_1_l_irrig_dmd"] = ".false.,"
jules_vegetation_nml.mapping["jules_vegetation_1_l_phenol"] = ".false.,"

jules_soil_txt = """
&jules_soil
confrac= ${jules_soil_1_confrac},
dzsoil_io= ${jules_soil_1_dzsoil_io},
l_bedrock= ${jules_soil_1_l_bedrock},
l_dpsids_dsdz= ${jules_soil_1_l_dpsids_dsdz},
l_soil_sat_down= ${jules_soil_1_l_soil_sat_down},
l_vg_soil= ${jules_soil_1_l_vg_soil},
sm_levels= ${jules_soil_1_sm_levels},
soilhc_method= ${jules_soil_1_soilhc_method},
zsmc= ${jules_soil_1_zsmc},
zst= ${jules_soil_1_zst},
/

"""
jules_soil_nml = julesNML(jules_soil_txt, "jules_soil.nml")
jules_soil_nml.mapping["jules_soil_1_l_bedrock"] = ".false.,"
jules_soil_nml.mapping["jules_soil_1_soilhc_method"] = "2,"
jules_soil_nml.mapping["jules_soil_1_zst"] = "1.0,"
jules_soil_nml.mapping["jules_soil_1_confrac"] = "0.3,"
jules_soil_nml.mapping["jules_soil_1_dzsoil_io"] = "0.1000,0.2500,0.6500,2.0000,"
jules_soil_nml.mapping["jules_soil_1_sm_levels"] = "4,"
jules_soil_nml.mapping["jules_soil_1_zsmc"] = "1.0,"
jules_soil_nml.mapping["jules_soil_1_l_vg_soil"] = ".false.,"
jules_soil_nml.mapping["jules_soil_1_l_dpsids_dsdz"] = ".false.,"
jules_soil_nml.mapping["jules_soil_1_l_soil_sat_down"] = ".true.,"

drive_txt = """
&jules_drive
	z1_uv_in = ${jules_drive_1_z1_uv_in},
	z1_tq_in = ${jules_drive_1_z1_tq_in},
	data_start  = ${jules_drive_1_data_start},
	data_end    = ${jules_drive_1_data_end},
	data_period = ${jules_drive_1_data_period},
	read_list = ${jules_drive_1_read_list},
	nfiles    = ${jules_drive_1_nfiles},
	file = ${jules_drive_1_file},
	nvars = ${jules_drive_1_nvars},
	var      = ${jules_drive_1_var},
	var_name = ${jules_drive_1_var_name},
	tpl_name = ${jules_drive_1_tpl_name},
	interp   = ${jules_drive_1_interp},
	l_daily_disagg= ${jules_drive_1_l_daily_disagg},
	l_imogen= ${jules_drive_1_l_imogen},
	l_perturb_driving= ${jules_drive_1_l_perturb_driving},
/

"""
drive_nml = julesNML(drive_txt, "drive.nml")
drive_nml.mapping[
    "jules_drive_1_interp"] = "'nf' 'nf' 'nf' 'nf' 'i' 'i' 'i' 'i'"
drive_nml.mapping["jules_drive_1_read_list"] = " T,"
drive_nml.mapping["jules_drive_1_data_end"] = " '2013-01-01 00:00:00',"
drive_nml.mapping["jules_drive_1_data_start"] = " '1979-01-01 00:00:00',"
drive_nml.mapping["jules_drive_1_l_daily_disagg"] = ".false.,"
drive_nml.mapping["jules_drive_1_z1_tq_in"] = " 2.0,"
drive_nml.mapping["jules_drive_1_nfiles"] = " 408,"
drive_nml.mapping[
    "jules_drive_1_var"] = "'sw_down' 'lw_down' 'tot_rain' 'tot_snow' 't' 'wind' 'pstar' 'q',"
drive_nml.mapping["jules_drive_1_l_imogen"] = ".false.,"
drive_nml.mapping["jules_drive_1_data_period"] = " 10800,"
drive_nml.mapping["jules_drive_1_z1_uv_in"] = " 10.0,"
drive_nml.mapping["jules_drive_1_file"] = " 'drive-file-WFDEI-79-12.txt',"
drive_nml.mapping[
    "jules_drive_1_var_name"] = "'SWdown' 'LWdown' 'Rainf' 'Snowf' 'Tair' 'Wind' 'PSurf' 'Qair',"
drive_nml.mapping["jules_drive_1_nvars"] = " 8,"
drive_nml.mapping["jules_drive_1_l_perturb_driving"] = ".false.,"
drive_nml.mapping[
    "jules_drive_1_tpl_name"] = "'SWdown_WFDEI' 'LWdown_WFDEI' 'Rainf_WFDEI_CRU' 'Snowf_WFDEI_CRU' 'Tair_WFDEI'" \
                                " 'Wind_WFDEI' 'PSurf_WFDEI' 'Qair_WFDEI',"

prescribed_data_txt = """
&jules_prescribed
n_datasets= ${jules_prescribed_1_n_datasets},
/

"""
prescribed_data_nml = julesNML(prescribed_data_txt, "prescribed_data.nml")
prescribed_data_nml.mapping["jules_prescribed_1_n_datasets"] = "0,"

jules_radiation_txt = """
&jules_radiation
l_albedo_obs= ${jules_radiation_1_l_albedo_obs},
l_cosz= ${jules_radiation_1_l_cosz},
l_embedded_snow= ${jules_radiation_1_l_embedded_snow},
l_mask_snow_orog= ${jules_radiation_1_l_mask_snow_orog},
l_niso_direct= ${jules_radiation_1_l_niso_direct},
l_snow_albedo= ${jules_radiation_1_l_snow_albedo},
l_spec_alb_bs= ${jules_radiation_1_l_spec_alb_bs},
l_spec_albedo= ${jules_radiation_1_l_spec_albedo},
wght_alb= ${jules_radiation_1_wght_alb},
/

"""
jules_radiation_nml = julesNML(jules_radiation_txt, "jules_radiation.nml")
jules_radiation_nml.mapping["jules_radiation_1_l_spec_albedo"] = ".true.,"
jules_radiation_nml.mapping["jules_radiation_1_l_embedded_snow"] = ".false.,"
jules_radiation_nml.mapping["jules_radiation_1_l_albedo_obs"] = ".false.,"
jules_radiation_nml.mapping["jules_radiation_1_l_niso_direct"] = ".false.,"
jules_radiation_nml.mapping["jules_radiation_1_l_cosz"] = ".true.,"
jules_radiation_nml.mapping["jules_radiation_1_l_spec_alb_bs"] = ".false.,"
jules_radiation_nml.mapping["jules_radiation_1_l_snow_albedo"] = ".false.,"
jules_radiation_nml.mapping["jules_radiation_1_l_mask_snow_orog"] = ".false.,"
jules_radiation_nml.mapping["jules_radiation_1_wght_alb"] = "0.0,0.5,0.0,0.5,"

jules_rivers_txt = """
&jules_rivers
l_rivers= ${jules_rivers_1_l_rivers},
/

"""
jules_rivers_nml = julesNML(jules_rivers_txt, "jules_rivers.nml")
jules_rivers_nml.mapping["jules_rivers_1_l_rivers"] = ".false.,"

timesteps_txt = """
&jules_time
  timestep_len = ${jules_time_1_timestep_len},
  main_run_start = ${jules_time_1_main_run_start},
  main_run_end   = ${jules_time_1_main_run_end},
/

&jules_spinup
  max_spinup_cycles = ${jules_spinup_1_max_spinup_cycles},
  spinup_start = ${jules_spinup_1_spinup_start},
  spinup_end   = ${jules_spinup_1_spinup_end},
  terminate_on_spinup_fail = ${jules_spinup_1_terminate_on_spinup_fail},
  nvars = ${jules_spinup_1_nvars},
  var = ${jules_spinup_1_var},
  tolerance = ${jules_spinup_1_tolerance},
  use_percent = ${jules_spinup_1_use_percent},
/

"""
timesteps_nml = julesNML(timesteps_txt, "timesteps.nml")
timesteps_nml.mapping["jules_spinup_1_tolerance"] = " 1.0 0.2"
timesteps_nml.mapping["jules_spinup_1_use_percent"] = " T F"
timesteps_nml.mapping["jules_spinup_1_spinup_end"] = " '1983-12-31 21:00:00',"
timesteps_nml.mapping["jules_spinup_1_nvars"] = " 2 ,"
timesteps_nml.mapping["jules_time_1_main_run_start"] = " '1979-01-01 00:00:00',"
timesteps_nml.mapping["jules_spinup_1_var"] = " 'smcl' 't_soil'"
timesteps_nml.mapping["jules_spinup_1_max_spinup_cycles"] = " 0"
timesteps_nml.mapping["jules_spinup_1_terminate_on_spinup_fail"] = " F"
timesteps_nml.mapping["jules_time_1_main_run_end"] = " '2012-12-31 21:00:00'"
timesteps_nml.mapping["jules_time_1_timestep_len"] = " 1800,"
timesteps_nml.mapping["jules_spinup_1_spinup_start"] = " '1983-01-01 00:00:00',"

jules_surface_types_txt = """
&jules_surface_types
ice= ${jules_surface_types_1_ice},
lake= ${jules_surface_types_1_lake},
ncpft= ${jules_surface_types_1_ncpft},
nnvg= ${jules_surface_types_1_nnvg},
npft= ${jules_surface_types_1_npft},
soil= ${jules_surface_types_1_soil},
urban= ${jules_surface_types_1_urban},
/

"""
jules_surface_types_nml = julesNML(jules_surface_types_txt, "jules_surface_types.nml")
jules_surface_types_nml.mapping["jules_surface_types_1_ncpft"] = "4,"
jules_surface_types_nml.mapping["jules_surface_types_1_nnvg"] = "4,"
jules_surface_types_nml.mapping["jules_surface_types_1_urban"] = "10,"
jules_surface_types_nml.mapping["jules_surface_types_1_soil"] = "12,"
jules_surface_types_nml.mapping["jules_surface_types_1_npft"] = "9,"
jules_surface_types_nml.mapping["jules_surface_types_1_ice"] = "13,"
jules_surface_types_nml.mapping["jules_surface_types_1_lake"] = "11,"

imogen_txt = """
"""
imogen_nml = julesNML(imogen_txt, "imogen.nml")

jules_hydrology_txt = """
&jules_hydrology
l_pdm= ${jules_hydrology_1_l_pdm},
l_top= ${jules_hydrology_1_l_top},
l_wetland_ch4_npp= ${jules_hydrology_1_l_wetland_ch4_npp},
l_wetland_unfrozen= ${jules_hydrology_1_l_wetland_unfrozen},
nfita= ${jules_hydrology_1_nfita},
ti_max= ${jules_hydrology_1_ti_max},
ti_wetl= ${jules_hydrology_1_ti_wetl},
zw_max= ${jules_hydrology_1_zw_max},
/

"""
jules_hydrology_nml = julesNML(jules_hydrology_txt, "jules_hydrology.nml")
jules_hydrology_nml.mapping["jules_hydrology_1_l_wetland_unfrozen"] = ".false.,"
jules_hydrology_nml.mapping["jules_hydrology_1_l_pdm"] = ".false.,"
jules_hydrology_nml.mapping["jules_hydrology_1_l_top"] = ".false.,"
jules_hydrology_nml.mapping["jules_hydrology_1_nfita"] = "20,"
jules_hydrology_nml.mapping["jules_hydrology_1_zw_max"] = "6.0,"
jules_hydrology_nml.mapping["jules_hydrology_1_ti_max"] = "10.0,"
jules_hydrology_nml.mapping["jules_hydrology_1_l_wetland_ch4_npp"] = ".false.,"
jules_hydrology_nml.mapping["jules_hydrology_1_ti_wetl"] = "1.5,"

initial_conditions_txt = """
&JULES_INITIAL
file = ${JULES_INITIAL_1_file},
total_snow = ${JULES_INITIAL_1_total_snow},
dump_file = ${JULES_INITIAL_1_dump_file},
/

"""
initial_conditions_nml = julesNML(initial_conditions_txt, "initial_conditions.nml")
initial_conditions_nml.mapping["JULES_INITIAL_1_total_snow"] = " F,"
initial_conditions_nml.mapping["JULES_INITIAL_1_dump_file"] = " T,"
initial_conditions_nml.mapping["JULES_INITIAL_1_file"] = " 'output/gl4.dump.spin3.19830101.0.nc',"

fire_txt = """
&fire_switches
l_fire= ${fire_switches_1_l_fire},
/

"""
fire_nml = julesNML(fire_txt, "fire.nml")
fire_nml.mapping["fire_switches_1_l_fire"] = ".false.,"

ancillaries_txt = """
&jules_frac
file= ${jules_frac_1_file},
frac_name = ${jules_frac_1_frac_name},
/

&jules_soil_props
  const_z   = ${jules_soil_props_1_const_z},
  use_file      = ${jules_soil_props_1_file},
  nvars     = ${jules_soil_props_1_nvars},
  var       = ${jules_soil_props_1_var},
  const_val  = ${jules_soil_props_1_const_val},
/

&jules_top
const_val= ${jules_top_1_const_val},
nvars= ${jules_top_1_nvars},
use_file= ${jules_top_1_use_file},
var= ${jules_top_1_var},
/

&jules_agric
/

&jules_crop_props
nvars= ${jules_crop_props_1_nvars},
use_file= ${jules_crop_props_1_use_file},
var= ${jules_crop_props_1_var},
var_name= ${jules_crop_props_1_var_name},
const_val= ${jules_crop_props_1_const_val},
/

&jules_co2
co2_mmr= ${jules_co2_1_co2_mmr},
/

"""
ancillaries_nml = julesNML(ancillaries_txt, "ancillaries.nml")
ancillaries_nml.mapping["jules_frac_1_frac_name"] = " 'field1391'"
ancillaries_nml.mapping["jules_top_1_const_val"] = "1.0,6.0,2.0,"
ancillaries_nml.mapping["jules_top_1_nvars"] = "3,"
ancillaries_nml.mapping["jules_soil_props_1_file"] = " 9*.false.,"
ancillaries_nml.mapping["jules_top_1_var"] = "'fexp','ti_mean','ti_sig',"
ancillaries_nml.mapping["jules_crop_props_1_var"] = "'cropttveg','cropttrep','cropsowdate',"
ancillaries_nml.mapping["jules_frac_1_file"] = "'land_frac_wheat100.nc',"
ancillaries_nml.mapping["jules_soil_props_1_nvars"] = " 9,"
ancillaries_nml.mapping["jules_soil_props_1_const_val"] = " 6.631272, 0.3967309, 0.0027729999, 0.45809999, " \
                                                          "0.3283205, 0.1866142, 1185786.0, 0.2269195, 0.17,"
ancillaries_nml.mapping["jules_crop_props_1_nvars"] = "3,"
ancillaries_nml.mapping["jules_crop_props_1_var_name"] = "'',"
ancillaries_nml.mapping["jules_soil_props_1_var"] = " 'b' 'sathh' 'satcon' 'sm_sat' 'sm_crit' 'sm_wilt' 'hcap' 'hcon'" \
                                                    " 'albsoil',"
ancillaries_nml.mapping["jules_crop_props_1_const_val"] = " 510.11138916   501.136169434 110.0 "
ancillaries_nml.mapping["jules_top_1_use_file"] = ".false.,.false.,.false.,"
ancillaries_nml.mapping["jules_co2_1_co2_mmr"] = "5.94100e-4,"
ancillaries_nml.mapping["jules_soil_props_1_const_z"] = " T,"
ancillaries_nml.mapping["jules_crop_props_1_use_file"] = "F,F,F"

jules_snow_txt = """
&jules_snow
can_clump= ${jules_snow_1_can_clump},
dtland= ${jules_snow_1_dtland},
frac_snow_subl_melt= ${jules_snow_1_frac_snow_subl_melt},
graupel_options= ${jules_snow_1_graupel_options},
i_snow_cond_parm= ${jules_snow_1_i_snow_cond_parm},
kland_numerator= ${jules_snow_1_kland_numerator},
l_et_metamorph= ${jules_snow_1_l_et_metamorph},
l_snow_infilt= ${jules_snow_1_l_snow_infilt},
l_snow_nocan_hc= ${jules_snow_1_l_snow_nocan_hc},
l_snowdep_surf= ${jules_snow_1_l_snowdep_surf},
lai_alb_lim_sn= ${jules_snow_1_lai_alb_lim_sn},
maskd= ${jules_snow_1_maskd},
n_lai_exposed= ${jules_snow_1_n_lai_exposed},
nsmax= ${jules_snow_1_nsmax},
rho_snow_const= ${jules_snow_1_rho_snow_const},
snow_hcap= ${jules_snow_1_snow_hcap},
snow_hcon= ${jules_snow_1_snow_hcon},
unload_rate_cnst= ${jules_snow_1_unload_rate_cnst},
unload_rate_u= ${jules_snow_1_unload_rate_u},
/

"""
jules_snow_nml = julesNML(jules_snow_txt, "jules_snow.nml")
jules_snow_nml.mapping["jules_snow_1_kland_numerator"] = "0.3,"
jules_snow_nml.mapping["jules_snow_1_snow_hcon"] = "0.265,"
jules_snow_nml.mapping["jules_snow_1_nsmax"] = "0,"
jules_snow_nml.mapping["jules_snow_1_n_lai_exposed"] = "9*0.0,"
jules_snow_nml.mapping["jules_snow_1_maskd"] = "50.0,"
jules_snow_nml.mapping["jules_snow_1_i_snow_cond_parm"] = "0,"
jules_snow_nml.mapping["jules_snow_1_l_et_metamorph"] = ".false.,"
jules_snow_nml.mapping["jules_snow_1_l_snowdep_surf"] = ".false.,"
jules_snow_nml.mapping["jules_snow_1_unload_rate_u"] = "9*0.0,"
jules_snow_nml.mapping["jules_snow_1_rho_snow_const"] = "250.0,"
jules_snow_nml.mapping["jules_snow_1_unload_rate_cnst"] = "9*0.0,"
jules_snow_nml.mapping["jules_snow_1_dtland"] = "2.0,"
jules_snow_nml.mapping["jules_snow_1_lai_alb_lim_sn"] = "9*0.5,"
jules_snow_nml.mapping["jules_snow_1_can_clump"] = "9*0.0,"
jules_snow_nml.mapping["jules_snow_1_l_snow_infilt"] = ".false.,"
jules_snow_nml.mapping["jules_snow_1_l_snow_nocan_hc"] = ".false.,"
jules_snow_nml.mapping["jules_snow_1_graupel_options"] = "0,"
jules_snow_nml.mapping["jules_snow_1_frac_snow_subl_melt"] = "0,"
jules_snow_nml.mapping["jules_snow_1_snow_hcap"] = "0.63e6,"

jules_soil_biogeochem_txt = """
&jules_soil_biogeochem
kaps= ${jules_soil_biogeochem_1_kaps},
l_layeredc= ${jules_soil_biogeochem_1_l_layeredc},
l_q10= ${jules_soil_biogeochem_1_l_q10},
l_soil_resp_lev2= ${jules_soil_biogeochem_1_l_soil_resp_lev2},
q10_soil= ${jules_soil_biogeochem_1_q10_soil},
soil_bgc_model= ${jules_soil_biogeochem_1_soil_bgc_model},
/

"""
jules_soil_biogeochem_nml = julesNML(jules_soil_biogeochem_txt, "jules_soil_biogeochem.nml")
jules_soil_biogeochem_nml.mapping["jules_soil_biogeochem_1_l_q10"] = ".true.,"
jules_soil_biogeochem_nml.mapping["jules_soil_biogeochem_1_l_soil_resp_lev2"] = ".false.,"
jules_soil_biogeochem_nml.mapping["jules_soil_biogeochem_1_soil_bgc_model"] = "1,"
jules_soil_biogeochem_nml.mapping["jules_soil_biogeochem_1_l_layeredc"] = ".false.,"
jules_soil_biogeochem_nml.mapping["jules_soil_biogeochem_1_q10_soil"] = "2.0,"
jules_soil_biogeochem_nml.mapping["jules_soil_biogeochem_1_kaps"] = "0.5e-8,"

if __name__ == "__main__":
    # misc_params_nml.mapping["q10_leaf"]="1.5"
    # print misc_params_nml.t.safe_substitute( misc_params_nml.mapping )

    output_nml.mapping["carb_profile_name"] = "carb_cs10"
    print output_nml.t.safe_substitute(output_nml.mapping)
