;2017/03/10
PRO JOURNAL__20170310__DO_A_SHIFTIE_FOR_ORBS_9585_AND_9627__NOT_OBSOLETE, $
   DATE=date, $
   ORBIT=orbit, $
   USE_AVGED_FOR_SMOOTH=use_avged_for_smooth, $
   PLOT_SMOOTHED_K_COMPONENTS=plot_smoothed_k_components, $
   SAVE_PS=save_ps, $
   TO_PDF=to_pdf, $
   REMOVE_EPS=remove_eps, $
   NO_PLOTS=no_plots, $
   NOSHIFT=noShift, $
   PUBLICATION_SETTINGS=pubSettings, $
   USE_DB_FAC=use_dB_fac

  COMPILE_OPT IDL2,STRICTARRSUBS

  IF ~KEYWORD_SET(noShift) THEN BEGIN

     minShiftBack  = -20
     maxShiftBack  = 20
     stepShiftBack = 1

     minShiftFwd   = -20
     maxShiftFwd   = 20
     stepShiftFwd  = 1

     univOffset    = 0

     ;; magShift      = 5
     ;; minShiftBack  = -magShift+univOffset
     ;; maxShiftBack  = magShift+univOffset
     ;; stepShiftBack = 1

     ;; minShiftFwd   = -magShift+univOffset
     ;; maxShiftFwd   = magShift+univOffset
     ;; stepShiftFwd  = 1

     nBack         = (maxShiftBack-minShiftBack)/stepShiftBack + 1
     nFwd          = (maxShiftFwd -minShiftFwd )/stepShiftFwd + 1
     backShifts    = INDGEN(nBack)*stepShiftBack+minShiftBack
     fwdShifts     = INDGEN(nFwd )*stepShiftFwd+minShiftFwd

     tmpSuff       = STRING(FORMAT='("NB",I0,"_NF",I0,"_univ",I0)',nBack,nFwd,univOffset)

  ENDIF

  IF orbit EQ 9627 THEN BEGIN
     custom_t1 = '1999-01-27/11:32:57.542'
     custom_t2 = '1999-01-27/11:33:09'
  ENDIF

  IF N_ELEMENTS(date) EQ 0 THEN BEGIN
     date = '20170309'
  ENDIF

  JOURNAL__20170224__ALL_THE_ORBS_WE_DONE_RECENTLY______PA_N_GEORGE, $
     ORBIT=orbit, $
     /PARSE_B_AND_J_SAVEFILE, $
     /USE_TIMEBAR_TIME__FROM_FILE, $
     /FOLD_NEGFREQ_ONTO_POS, $
     /PLOT_KPERP_MAGNITUDE_FOR_KZ, $
     /PLOT_KX_VS_KY_FOR_KZ, $
     SAVE_PS=save_ps, $
     TO_PDF=to_pdf, $
     REMOVE_EPS=remove_eps, $
     NO_PLOTS=no_plots, $
     /USE_LOWRES_TIME_SERIES, $
     DATE=date, $
     PUBLICATION_SETTINGS=pubSettings, $
     BONUS_SUFF=tmpSuff, $
     USE_DB_FAC=use_dB_fac, $
     CUSTOM_T1=custom_t1, $
     CUSTOM_T2=custom_t2, $
     /USE_REPRETCAL_FILE, $
     KP__ANGLERANGE=[-90,270], $
     PRE_VIII_LAYOUT=KEYWORD_SET(save_ps), $
     PLOT_SMOOTHED_K_COMPONENTS=plot_smoothed_k_components, $
     BACKSHIFTS_FOR_AVGING=backShifts, $
     FWDSHIFTS_FOR_AVGING=fwdShifts, $
     USE_AVGED_FOR_SMOOTH=use_avged_for_smooth, $
     OUT_FREQS=out_freqs, $
     OUT_KX=out_kx, $
     OUT_KY=out_ky, $
     OUT_KZ=out_kz, $
     OUT_KP=out_kP, $
     OUT_ANGLE_KP=out_kPAngle, $
     ;; OUT_INDS=out_inds, $
     OUT_TARR=out_TArr, $
     OUT_USEDINDS=out_usedInds, $
     OUT_BX=out_Bx, $
     OUT_BY=out_By, $
     OUT_BZ=out_Bz, $
     OUT_JX=out_Jx, $
     OUT_JY=out_Jy, $
     OUT_JZ=out_Jz

END
