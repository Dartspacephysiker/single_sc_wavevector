;2017/03/10
;Wow, a time consumer:
;JOURNAL__20170310__DO_A_SHIFTIE_FOR_ORBS_9585_AND_9627,/SAVE_PS,/USE_DB_FAC,/PUBLICATION_SETTINGS,ORBIT=9585,date='20170310',/PLOT_SMOOTHED_K_COMPONENTS,/LOCK_FWDSHIFT_TO_BACKSHIFT,CUSTOM_T1='1999-01-23/14:50:55',CUSTOM_T2='1999-01-23/14:51:07',KX_SPECIALFREQS=[[0.9,2],[3.6,4.4]],/MAKE_KX_VS_KY_SPECIAL,KX_SPECIALBOUNDS=[-3.0,-0.5],/MAKE_KPANGLE_SPECIAL
;
;WINNER #1:
;JOURNAL__20170310__DO_A_SHIFTIE_FOR_ORBS_9585_AND_9627,/SAVE_PS,/USE_DB_FAC,/PUBLICATION_SETTINGS,ORBIT=9627,/PLOT_SMOOTHED_K_COMPONENTS,/LOCK_FWDSHIFT_TO_BACKSHIFT,KP__ANGLERANGE=[-45,315]
;
;WINNER #2:
;JOURNAL__20170310__DO_A_SHIFTIE_FOR_ORBS_9585_AND_9627,/SAVE_PS,/USE_DB_FAC,/PUBLICATION_SETTINGS,ORBIT=10832,date='20170310',/PLOT_SMOOTHED_K_COMPONENTS,/LOCK_FWDSHIFT_TO_BACKSHIFT,KP__ANGLERANGE=[-120,240]
;
;ET FINALEMENT (2017/03/11):
;
;JOURNAL__20170310__DO_A_SHIFTIE_FOR_ORBS_9585_AND_9627,/SAVE_PS,/USE_DB_FAC,/PUBLICATION_SETTINGS,ORBIT=10832,date='20170310',/PLOT_SMOOTHED_K_COMPONENTS,/LOCK_FWDSHIFT_TO_BACKSHIFT,KP__ANGLERANGE=[0,360],PAGE1__FREQRANGE=[0,6],PAGE2__FREQRANGE=[0,6]
;
;JOURNAL__20170310__DO_A_SHIFTIE_FOR_ORBS_9585_AND_9627,/SAVE_PS,/USE_DB_FAC,/PUBLICATION_SETTINGS,ORBIT=9627,/PLOT_SMOOTHED_K_COMPONENTS,/LOCK_FWDSHIFT_TO_BACKSHIFT,KP__ANGLERANGE=[0,360],PAGE1__FREQRANGE=[0,6],PAGE2__FREQRANGE=[0,6]
PRO JOURNAL__20170310__DO_A_SHIFTIE_FOR_ORBS_9585_AND_9627, $
   CUSTOM_T1=custom_t1, $
   CUSTOM_T2=custom_t2, $
   DATE=date, $
   ORBIT=orbit, $
   USE_AVGED_FOR_SMOOTH=use_avged_for_smooth, $
   AVG_BINSIZE=avg_binSize, $
   KX_SPECIALFREQS=kx_specialFreqs, $
   KY_SPECIALFREQS=ky_specialFreqs, $
   KPANGLE_SPECIALFREQS=kPAngle_specialFreqs, $
   KX_SPECIALBOUNDS=kx_specialBounds, $
   KY_SPECIALBOUNDS=ky_specialBounds, $
   KPANGLE_SPECIALBOUNDS=kPAngle_specialBounds, $
   MAKE_KX_VS_KY_SPECIAL=make_kx_vs_ky_special, $
   MAKE_KPANGLE_SPECIAL=make_kPAngle_special, $
   PLOT_SMOOTHED_K_COMPONENTS=plot_smoothed_k_components, $
   FREQLIMS=freqLims, $
   PAGE1__FREQRANGE=page1__freqRange, $
   PAGE2__FREQRANGE=page2__freqRange, $
   KP__ANGLERANGE=kP__angleRange, $
   SAVE_PS=save_ps, $
   TO_PDF=to_pdf, $
   REMOVE_EPS=remove_eps, $
   NO_PLOTS=no_plots, $
   LOCK_FWDSHIFT_TO_BACKSHIFT=lock_shifts, $
   NOSHIFT=noShift, $
   NO_PTSHIFT_BACK=no_ptShift_back, $
   NO_PTSHIFT_FWD=no_ptShift_fwd, $
   PUBLICATION_SETTINGS=pubSettings, $
   USE_DB_FAC=use_dB_fac

  COMPILE_OPT IDL2,STRICTARRSUBS

  IF ~KEYWORD_SET(noShift) THEN BEGIN

     ;; minShiftBack  = -20
     ;; maxShiftBack  = 20
     ;; stepShiftBack = 1

     ;; minShiftFwd   = -20
     ;; maxShiftFwd   = 20
     ;; stepShiftFwd  = 1

     univOffset    = 0

     ptShift       = 20
     stepShift     = 1
     minShiftBack  = -ptShift+univOffset
     maxShiftBack  = ptShift+univOffset
     stepShiftBack = stepShift

     minShiftFwd   = -ptShift+univOffset
     maxShiftFwd   = ptShift+univOffset
     stepShiftFwd  = stepShift

     nFwd          = (maxShiftFwd -minShiftFwd )/stepShiftFwd + 1

     IF KEYWORD_SET(no_ptShift_back) THEN BEGIN
        nBack      = 0
        backShifts = 0
        BKStr      = ''
     ENDIF ELSE BEGIN
        nBack      = (maxShiftBack-minShiftBack)/stepShiftBack + 1
        backShifts = INDGEN(nBack)*stepShiftBack+minShiftBack
        BKStr      = STRING(FORMAT='("NB",I0)',nBack)
     ENDELSE

     CASE 1 OF
        KEYWORD_SET(no_ptShift_fwd): BEGIN
           nFwd       = 0
           fwdShifts  = 0
           FWStr      = ''
        END
        KEYWORD_SET(lock_shifts): BEGIN
           nFwd       = 1
           FWStr      = '_LkF'
        END
        ELSE: BEGIN
           fwdShifts  = INDGEN(nFwd )*stepShiftFwd+minShiftFwd
           FWStr      = STRING(FORMAT='("_NF",I0)',nFwd)
        END
     ENDCASE

     univStr          = ''
     IF KEYWORD_SET(univOffset) THEN BEGIN
        univStr       = STRING(FORMAT='("_univ",I0)',univOffset)
     ENDIF

     tmpSuff          = STRING(FORMAT='(A0,A0,A0)',BKStr,FWStr,univStr)

  ENDIF

  IF N_ELEMENTS(kP__angleRange) EQ 0 THEN BEGIN
     kP__angleRange = [-90,270]
  ENDIF

  IF orbit EQ 9627 THEN BEGIN

     IF ~KEYWORD_SET(custom_t1) THEN BEGIN
        custom_t1 = '1999-01-27/11:32:56.542'
     ENDIF

     IF ~KEYWORD_SET(custom_t2) THEN BEGIN
        custom_t2 = '1999-01-27/11:33:09.000' ;actual is 08.978
     ENDIF

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
     FREQLIMS=freqLims, $
     PAGE1__FREQRANGE=page1__freqRange, $
     PAGE2__FREQRANGE=page2__freqRange, $
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
     KP__ANGLERANGE=kP__angleRange, $
     PRE_VIII_LAYOUT=KEYWORD_SET(save_ps), $
     PLOT_SMOOTHED_K_COMPONENTS=plot_smoothed_k_components, $
     BACKSHIFTS_FOR_AVGING=backShifts, $
     FWDSHIFTS_FOR_AVGING=fwdShifts, $
     LOCK_FWDSHIFT_TO_BACKSHIFT=lock_shifts, $
     USE_AVGED_FOR_SMOOTH=use_avged_for_smooth, $
     AVG_BINSIZE=avg_binSize, $
     KX_SPECIALFREQS=kx_specialFreqs, $
     KY_SPECIALFREQS=ky_specialFreqs, $
     KPANGLE_SPECIALFREQS=kPAngle_specialFreqs, $
     KX_SPECIALBOUNDS=kx_specialBounds, $
     KY_SPECIALBOUNDS=ky_specialBounds, $
     KPANGLE_SPECIALBOUNDS=kPAngle_specialBounds, $
     MAKE_KX_VS_KY_SPECIAL=make_kx_vs_ky_special, $
     MAKE_KPANGLE_SPECIAL=make_kPAngle_special, $
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
     OUT_JZ=out_Jz, $
     OUT_AVGJXBNRM=out_avgJxBNrm

  STOP

END
