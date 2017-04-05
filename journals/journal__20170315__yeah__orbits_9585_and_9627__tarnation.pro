;2017/03/15
;So Jim and I settled on orbits 9585 and 9627. Orbit 10832 could be reserved for another day.
;
;Maybe final—but you know, you never really know:
;JOURNAL__20170315__YEAH__ORBITS_9585_AND_9627__TARNATION,ORBIT=9585,/FOOTBALL_LAYOUT,/FOOTBALL_COL2TITLE,/FOOTBALL_YLOG,KP__ANGLERANGE=[-45,315],PAGE1__FREQRANGE=[0,6],PAGE2__FREQRANGE=[0,6],/SAVE_PS,/PUBLICATION_SETTINGS,CUSTOM_MULTI_SHIFTSEC=[0.D],/ITVL_MODE,INTERVAL=0
;
;
PRO JOURNAL__20170315__YEAH__ORBITS_9585_AND_9627__TARNATION, $
   CUSTOM_T1=custom_t1, $
   CUSTOM_T2=custom_t2, $
   CUSTOM_T_EXTEND=custom_t_extend, $
   CUSTOM_MULTI_T_EXTEND=custom_multi_t_extend, $
   CUSTOM_SHIFTSEC=custom_shiftSec, $
   CUSTOM_MULTI_SHIFTSEC=custom_multi_shiftSec, $
   DATE=date, $
   ORBIT=orbit, $
   PLOT_POSFREQ=plot_posFreq, $
   FOLD_NEGFREQ_ONTO_POS=fold_negFreq, $
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
   MARK_KS_BELOW_MAGERR_THRESH=mark_ks_below_magErr_thresh, $
   MARK_KS_BELOW_ERRANGLE_THRESH=mark_ks_below_errAngle_thresh, $
   MARK_KS_BELOW_PHASE_ERR_THRESH=mark_ks_below_phaseErr_thresh, $
   MARK_KS_BELOW_BOTH=mark_ks_below_both, $
   ONLY_KS_BELOW_MAGERR_THRESH=only_ks_below_magErr_thresh, $
   ;; ONLY_KS_BELOW_ERRANGLE_THRESH=only_ks_below_errAngle_thresh, $
   ;; ONLY_KS_BELOW_PHASE_ERR_THRESH=only_ks_below_phaseErr_thresh, $
   PLOT_SMOOTHED_K_COMPONENTS=plot_smoothed_k_components, $
   FREQLIMS=freqLims, $
   PAGE1__FREQRANGE=page1__freqRange, $
   PAGE2__FREQRANGE=page2__freqRange, $
   KP__ANGLERANGE=kP__angleRange, $
   WRITE_ASCIIS=write_ASCIIs, $
   SAVE_PS=save_ps, $
   FOOTBALL_LAYOUT=football_layout, $
   FOOTBALL_YLOG=football_yLog, $
   FOOTBALL_COL2TITLE=football_col2Title, $
   FOOTBALL_KMAG=football_kMag, $
   FOOTBALL_NO_ERRANGLE=football_no_errAngle, $
   FOOTBALL_NO_MAGERR=football_no_magErr, $
   TO_PDF=to_pdf, $
   REMOVE_EPS=remove_eps, $
   SHOW_PREDICTED_J=show_predicted_J, $
   NO_PLOTS=no_plots, $
   PTSHIFT=ptShift, $
   SHIFT_UNIVERSAL_OFFSET=univOffset, $
   NOSHIFT=noShift, $
   NO_PTSHIFT_BACK=no_ptShift_back, $
   NO_PTSHIFT_FWD=no_ptShift_fwd, $
   LOCK_FWDSHIFT_TO_BACKSHIFT=lock_shifts, $
   PUBLICATION_SETTINGS=pubSettings, $
   USE_DB_FAC=use_dB_fac, $
   ITVL_MODE=itvl_mode, $
   INTERVAL=interval

  COMPILE_OPT IDL2,STRICTARRSUBS

  tmpSuff = ''

  IF N_ELEMENTS(kP__angleRange) EQ 0 THEN BEGIN
     kP__angleRange = [-90,270]
  ENDIF

  IF KEYWORD_SET(itvl_mode) THEN BEGIN

     IF N_ELEMENTS(interval) EQ 0 THEN BEGIN
        interval = 0
     ENDIF

     tmpSuff += '-itvl' + STRCOMPRESS(interval,/REMOVE_ALL)

     CASE orbit OF
        9585: BEGIN

           CASE 1 OF
              interval EQ 0: BEGIN

                 ;; custom_t1 = '1999-01-23/14:50:56.0'
                 ;; custom_t2 = '1999-01-23/14:51:08.6' ;possibly #1, shifted back 0.5

                 custom_t1 = '1999-01-23/14:50:56.0'
                 custom_t2 = '1999-01-23/14:51:07.0' ;Way good, no shifting necessary

                 ;; custom_t1 = '1999-01-23/14:50:48.5'
                 ;; custom_t2 = '1999-01-23/14:51:07.0' 

                 ;; custom_t1 = '1999-01-23/14:50:53.33'
                 ;; custom_t2 = '1999-01-23/14:51:18.0'

              END
              interval EQ 1: BEGIN
                 custom_t1 = '1999-01-23/14:50:51.0'
                 custom_t2 = '1999-01-23/14:51:01.0'
              END
              interval EQ 2: BEGIN
                 custom_t1 = '1999-01-23/14:51:01.0'
                 custom_t2 = '1999-01-23/14:51:11.0'
              END
              interval EQ 3: BEGIN
                 custom_t1 = '1999-01-23/14:50:35.0'
                 ;; custom_t2 = '1999-01-23/14:50:50.0'
                 custom_t2 = '1999-01-23/14:51:13.0'
              END
           ENDCASE

        END
        9627: BEGIN

           CASE 1 OF
              interval EQ 0: BEGIN
                 ;; custom_t1 = '1999-01-27/11:32:56.5'
                 ;; custom_t2 = '1999-01-27/11:33:09.0'
                 custom_t1 = '1999-01-27/11:32:51.5'
                 custom_t2 = '1999-01-27/11:33:09.5'
              END
              interval EQ 1: BEGIN
                 custom_t1 = '1999-01-27/11:32:54.0'
                 custom_t2 = '1999-01-27/11:33:04.0'
              END
              interval EQ 2: BEGIN
                 custom_t1 = '1999-01-27/11:33:04.0'
                 custom_t2 = '1999-01-27/11:33:14.0'
              END
              interval EQ 3: BEGIN
                 custom_t1 = '1999-01-27/11:32:45.5'
                 custom_t2 = '1999-01-27/11:32:54.1'
              END
           ENDCASE

        END
        10832: BEGIN
             
           CASE 1 OF
              interval EQ 0: BEGIN
                 custom_t1 = '1999-05-18/06:50:57.1'
                 custom_t2 = '1999-05-18/06:51:07.996'
              END
              interval EQ 1: BEGIN
                 custom_t1 = '1999-05-18/06:50:45.0'
                 custom_t2 = '1999-05-18/06:50:59.5'
                 ;; custom_t1 = '1999-05-18/06:50:45.0'
                 ;; custom_t2 = '1999-05-18/06:51:04.5'
              END
              interval EQ 2: BEGIN
                 custom_t1 = '1999-05-18/06:50:59.5'
                 custom_t2 = '1999-05-18/06:51:14.0'
              END
              interval EQ 3: BEGIN
                 custom_t1 = '1999-05-18/06:50:45.0'
                 custom_t2 = '1999-05-18/06:51:15.0'
              END
           ENDCASE

        END        
        ;; 10837: BEGIN

        ;; END
     ENDCASE

  ENDIF ELSE BEGIN

     IF orbit EQ 9627 THEN BEGIN

        IF ~KEYWORD_SET(custom_t1) THEN BEGIN
           custom_t1 = '1999-01-27/11:32:56.542'
        ENDIF

        IF ~KEYWORD_SET(custom_t2) THEN BEGIN
           custom_t2 = '1999-01-27/11:33:09.000' ;actual is 08.978
        ENDIF

     ENDIF
     
  ENDELSE

  IF KEYWORD_SET(football_col2Title) THEN BEGIN

     IF (WHERE(SIZE(football_col2Title,/TYPE) EQ [1,2,3,12,13,14,15]))[0] NE -1 THEN BEGIN

        football_col2Title = STRING(FORMAT='("Orbit ",I0)',orbit)

        IF N_ELEMENTS(interval) GT 0 THEN BEGIN
           ;; football_col2Title += STRING(FORMAT='(" (itvl ",I0,")")',interval)
        ENDIF ELSE BEGIN
           ;; football_col2Title += "(original)"
        ENDELSE

     ENDIF

  ENDIF

  ;; IF N_ELEMENTS(date) EQ 0 THEN BEGIN
  ;;    date = '20170309'
  ;; ENDIF

  IF N_ELEMENTS(date) EQ 0 THEN BEGIN
     date = '20170313'
  ENDIF

  JOURNAL__20170224__ALL_THE_ORBS_WE_DONE_RECENTLY______PA_N_GEORGE, $
     ORBIT=orbit, $
     /PARSE_B_AND_J_SAVEFILE, $
     /USE_TIMEBAR_TIME__FROM_FILE, $
     PLOT_POSFREQ=plot_posFreq, $
     FOLD_NEGFREQ_ONTO_POS=fold_negFreq, $
     /PLOT_KPERP_MAGNITUDE_FOR_KZ, $
     /PLOT_KX_VS_KY_FOR_KZ, $
     FREQLIMS=freqLims, $
     PAGE1__FREQRANGE=page1__freqRange, $
     PAGE2__FREQRANGE=page2__freqRange, $
     SAVE_PS=save_ps, $
     TO_PDF=to_pdf, $
     REMOVE_EPS=remove_eps, $
     NO_PLOTS=no_plots, $
     SHOW_PREDICTED_J=show_predicted_J, $
     /USE_LOWRES_TIME_SERIES, $
     DATE=date, $
     PUBLICATION_SETTINGS=pubSettings, $
     BONUS_SUFF=tmpSuff, $
     USE_DB_FAC=use_dB_fac, $
     CUSTOM_T1=custom_t1, $
     CUSTOM_T2=custom_t2, $
     CUSTOM_T_EXTEND=custom_t_extend, $
     CUSTOM_MULTI_T_EXTEND=custom_multi_t_extend, $
     CUSTOM_SHIFTSEC=custom_shiftSec, $
     CUSTOM_MULTI_SHIFTSEC=custom_multi_shiftSec, $
     /USE_REPRETCAL_FILE, $
     KP__ANGLERANGE=kP__angleRange, $
     PRE_VIII_LAYOUT=KEYWORD_SET(save_ps) AND ~KEYWORD_SET(football_layout), $
     FOOTBALL_LAYOUT=football_layout, $
     FOOTBALL_YLOG=football_yLog, $
     FOOTBALL_COL2TITLE=football_col2Title, $
     FOOTBALL_KMAG=football_kMag, $
     FOOTBALL_NO_ERRANGLE=football_no_errAngle, $
     FOOTBALL_NO_MAGERR=football_no_magErr, $
     PLOT_SMOOTHED_K_COMPONENTS=plot_smoothed_k_components, $
     PTSHIFT=ptShift, $
     SHIFT_UNIVERSAL_OFFSET=univOffset, $
     NOSHIFT=noShift, $
     NO_PTSHIFT_BACK=no_ptShift_back, $
     NO_PTSHIFT_FWD=no_ptShift_fwd, $
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
     MARK_KS_BELOW_MAGERR_THRESH=mark_ks_below_magErr_thresh, $
     MARK_KS_BELOW_ERRANGLE_THRESH=mark_ks_below_errAngle_thresh, $
     MARK_KS_BELOW_PHASE_ERR_THRESH=mark_ks_below_phaseErr_thresh, $
     MARK_KS_BELOW_BOTH=mark_ks_below_both, $
     ONLY_KS_BELOW_MAGERR_THRESH=only_ks_below_magErr_thresh, $
     ;; ONLY_KS_BELOW_ERRANGLE_THRESH=only_ks_below_errAngle_thresh, $
     ;; ONLY_KS_BELOW_PHASE_ERR_THRESH=only_ks_below_phaseErr_thresh, $
     WRITE_ASCIIS=write_ASCIIs, $
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

  ;; STOP

END

