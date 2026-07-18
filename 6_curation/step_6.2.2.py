###################################################################################################
# STEP 6.3.2: Annotate model with BiGG IDs
###################################################################################################

# libraries
import re
import cobra

# paths
MODEL_IN = '../models/model_3.xml'
MODEL_OUT = '../models/model_3.xml'

# stuff to force annotate from memote report
reactions_to_fix = ["13PPDH2_1","3GMPtex_1","3MBt4pp","4HPRO_LT_Et","4HTHRS_1","ABTt_1","ACACT5r_1","ACACT6r_1","ACOAD4_1","ACOAD5_1","ACPS1_1","ACS_1","ACTD_1","ADPRDP_1","AHMMPS_1","ALDD31_1","AMANK_1","AMPD2","ARGDI_1","ARGN_1","ASPCT_2","ATAH_1","BG_CELLB","BG_MADG","BG_MBDG","BPNT_1","CBPS_1","CELLBpts_2","CITt4pp_1","COt_1","DADNK_1","DGNSK_1","DHNAOT_1","DURIPP_1","ECOAH1_1","ECOAH5_1","FACOAL140_ISO","FACOAL150_ISO","FACOAL150_anteiso","FACOAL160_ISO","FACOAL170_ISO","FACOAL170_anteiso","FACOAL180_2","FCLT_2","FDMO3_1","FDMO4_1","FGLU_1","GALt2_3","GCCa_1","GCCb_1","GCCc_1","GLUDy_1","GLUKA_1","GNK_1","GTPDPK_1","HACD1_1","HEMEOS_1","HISTP_1","HMGL_2","IG3PS_1","IGPDH_1","IZPN_1","LIPATPT_2","MCTC_1","MTAM_1","NADDPp_1","NADK_1","NH4tpp_1","NMNHYD_1","NOR_syn_1","NOS1_1","NTPP10_1","NTPP11_1","NTPP1_1","NTPP2_1","NTPP8_1","NTPP9_1","OCBT_1","ORNTAC_1","ORNTA_1","PFK_4","PGLYCP_1","PLPS_1","PNTK_1","PRAIS_1","PRAMPC_1","PRFGS_1","PTHPS_1","PYDXK_1","PYNP1_1","PYNP2r_1","RBFSa_1","SHCHD2_1","SHCHF_2","THRA2","THRS_2","THZPSN_1","TMDPK_1","UAGPT2_1","UPP3MT_2","URATEt_1","UREA_1","r2465_1","sink_hemeO_c","sink_lipopb_c","sink_sheme_c","EX_23camp_e","EX_23ccmp_e","EX_23cgmp_e","EX_23cump_e","EX_23dhbzs3_e","EX_2ameph_e","EX_2ddglcn_e","EX_2hxmp_e","EX_2m35mdntha_e","EX_2pglyc_e","EX_35dnta_e","EX_3amp_e","EX_3cmp_e","EX_3gmp_e","EX_3mb_e","EX_3ump_e","EX_4abut_e","EX_4ahmmp_e","EX_4hpro_DC_e","EX_4hpro_LT_e","EX_LalaLglu_e","EX_abg4_e","EX_abt__L_e","EX_acac_e","EX_acgam_e","EX_acmana_e","EX_acon_C_e","EX_actn__R_e","EX_agm_e","EX_ala_L_asp__L_e","EX_ala_L_glu__L_e","EX_ala_L_thr__L_e","EX_ala__D_e","EX_ala__L_e","EX_ala_gln_e","EX_ala_his_e","EX_ala_leu_e","EX_alaala_e","EX_alagly_e","EX_alahis_e","EX_alaleu_e","EX_alathr_e","EX_alatrp_e","EX_arab__L_e","EX_arbt_e","EX_arg__L_e","EX_asn__L_e","EX_asp__L_e","EX_balaala_e","EX_balabala_e","EX_balagly_e","EX_balaleu_e","EX_balamd_e","EX_bhb_e","EX_but_e","EX_butso3_e","EX_bz_e","EX_carn_e","EX_cellb_e","EX_cgly_e","EX_cmcbtt_e","EX_cpgn_e","EX_cpgn_un_e","EX_cys__L_e","EX_d23hb_e","EX_dca_e","EX_dcyt_e","EX_ddca_e","EX_drib_e","EX_enter_e","EX_ethso3_e","EX_fald_e","EX_fcmcbtt_e","EX_fe3dcit_e","EX_fe3dhbzs3_e","EX_fe3mcbtt_e","EX_fe3pyovd_kt_e","EX_feenter_e","EX_feoxam_e","EX_feoxam_un_e","EX_fol_e","EX_fru_e","EX_gal_bD_e","EX_galt_e","EX_galur_e","EX_gam_e","EX_glc__D_e","EX_glc__aD_e","EX_glcn__D_e","EX_glcn_e","EX_glcr_e","EX_glcur_e","EX_gln__L_e","EX_glu__L_e","EX_gly_asn__L_e","EX_gly_asp__L_e","EX_gly_cys_e","EX_gly_gln_e","EX_gly_glu__L_e","EX_gly_leu_e","EX_gly_met_e","EX_gly_phe_e","EX_gly_pro__L_e","EX_gly_tyr_e","EX_glyald_e","EX_glygln_e","EX_glyglu_e","EX_glygly_e","EX_glyglygln_e","EX_glymet_e","EX_glyphe_e","EX_glyser_e","EX_gthox_e","EX_h2o2_e","EX_hexs_e","EX_his__L_e","EX_hisgly_e","EX_hishis_e","EX_hqn_e","EX_hxa_e","EX_hxan_e","EX_icit_e","EX_ile__L_e","EX_ind3ac_e","EX_inost_e","EX_isetac_e","EX_istfrnA_e","EX_istfrnB_e","EX_leu__L_e","EX_leuleu_e","EX_lipoate_e","EX_lys__L_e","EX_lysglugly_e","EX_madg_e","EX_mal__D_e","EX_mal__L_e","EX_malttr_e","EX_maltttr_e","EX_man_e","EX_manglyc_e","EX_mbdg_e","EX_mcbtt_e","EX_melib_e","EX_met_L_ala__L_e","EX_met__D_e","EX_met__L_e","EX_metox__R_e","EX_metox_e","EX_metsox_R__L_e","EX_metsox_S__L_e","EX_mmet_e","EX_mnl_e","EX_mso3_e","EX_mththf_e","EX_nac_e","EX_nmn_e","EX_octa_e","EX_oxa_e","EX_pacald_e","EX_phe__L_e","EX_pnto__R_e","EX_pro__L_e","EX_progly_e","EX_prohisglu_e","EX_ptrc_e","EX_pyovd_kt_e","EX_rib__D_e","EX_ribflv_e","EX_s_e","EX_salchs2_e","EX_salchs2fe_e","EX_salchs4_e","EX_salchs4fe_e","EX_salcn_e","EX_ser__D_e","EX_ser__L_e","EX_serglugly_e","EX_so3_e","EX_spmd_e","EX_stfrnA_e","EX_stfrnB_e","EX_sucr_e","EX_sulfac_e","EX_tagur_e","EX_tartr__L_e","EX_taur_e","EX_thm_e","EX_thr__L_e","EX_thymd_e","EX_tnt_e","EX_trp__L_e","EX_ttdcea_e","EX_tyr__L_e","EX_udcpdp_e","EX_udcpp_e","EX_ura_e","EX_urate_e","EX_uri_e","EX_val__L_e","EX_xan_e","EX_xyl__D_e","EX_zn2_e","Growth","PHBS_syn_1","SK_phb_c"]
metabolites_to_fix = ["glc__aD_c","glc__aD_e","glc__bD_c"]

# import model
model = cobra.io.read_sbml_model(MODEL_IN)

# Annotate reactions
rxn_fixed, rxn_missing = 0, []
for rxn_id in reactions_to_fix:
    if rxn_id not in model.reactions:
        rxn_missing.append(rxn_id)
        continue
    rxn = model.reactions.get_by_id(rxn_id)
    if rxn_id.startswith(('EX_', 'sink_', 'SK_')):
        clean_id = rxn_id
    else:
        clean_id = re.sub(r'_[a-z]$', '', rxn_id)
    if rxn.annotation is None:
        rxn.annotation = {}
    rxn.annotation["bigg.reaction"] = [clean_id]
    rxn_fixed += 1

# Annotate metabolites
met_fixed, met_missing = 0, []
for met_id in metabolites_to_fix:
    if met_id not in model.metabolites:
        met_missing.append(met_id)
        continue
    met = model.metabolites.get_by_id(met_id)
    clean_id = re.sub(r'_[a-z]$', '', met_id)
    if met.annotation is None:
        met.annotation = {}
    met.annotation["bigg.metabolite"] = [clean_id]
    met_fixed += 1

print(f"Annotated {rxn_fixed}/{len(reactions_to_fix)} reactions, "
      f"{met_fixed}/{len(metabolites_to_fix)} metabolites.")
if rxn_missing:
    print(f"Reactions not found in model ({len(rxn_missing)}): {rxn_missing}")
if met_missing:
    print(f"Metabolites not found in model ({len(met_missing)}): {met_missing}")

# export model
cobra.io.write_sbml_model(model, MODEL_OUT)
print(f"Annotated model saved to {MODEL_OUT}")

###################################################################################################
