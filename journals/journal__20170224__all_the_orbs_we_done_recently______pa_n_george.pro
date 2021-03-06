;2017/02/24
;Here's a keeper:
;JOURNAL__20170224__ALL_THE_ORBS_WE_DONE_RECENTLY______PA_N_GEORGE,/PARSE_B_AND_J_SAVEFILE,/USE_TIMEBAR_TIME__FROM_FILE,/PLOT_KPERP_MAGNITUDE_FOR_KZ,/PLOT_KX_VS_KY_FOR_KZ,/FOLD_NEGFREQ_ONTO_POS,/USE_LOWRES_TIME_SERIES,DATE='20170224',ORBIT=10837,/SAVE_PS,/TO_PDF,/REMOVE_EPS,FFTPERCENT=25,CUSTOM_T1='1999-05-18/17:46:22.0',CUSTOM_T2='1999-05-18/17:46:40.0'
;
;Also try the above with FFTPERCENT=100. The smoothing is evident, and it looks like there's some fo' real signal happening.
;
;2017/02/27 Jim likee: JOURNAL__20170224__ALL_THE_ORBS_WE_DONE_RECENTLY______PA_N_GEORGE,ORBIT=9627,/PARSE_B_AND_J_SAVEFILE,/USE_TIMEBAR_TIME__FROM_FILE,/FOLD_NEGFREQ_ONTO_POS,/PLOT_KPERP_MAGNITUDE_FOR_KZ,/PLOT_KX_VS_KY_FOR_KZ,/SAVE_PS,/TO_PDF,/REMOVE_EPS,/USE_LOWRES_TIME_SERIES,FFTPERCENT=50,WHICH_FFTS=1,DATE='20170224',/PUBLICATION_SETTINGS,BONUS_SUFF='SDFSDF'
;
;2017/02/28 One that I like:
;JOURNAL__20170224__ALL_THE_ORBS_WE_DONE_RECENTLY______PA_N_GEORGE,ORBIT=9627,/PARSE_B_AND_J_SAVEFILE,/USE_TIMEBAR_TIME__FROM_FILE,/FOLD_NEGFREQ_ONTO_POS,/PLOT_KPERP_MAGNITUDE_FOR_KZ,/PLOT_KX_VS_KY_FOR_KZ,/SAVE_PS,/USE_LOWRES_TIME_SERIES,/PUBLICATION_SETTINGS,BONUS_SUFF='SDFSDF',DATE='20170228',/TO_PDF,/PLOT_SMOOTHED_K_COMPONENTS,KP__ANGLERANGE=[0,360],/THIRD_PAGE
;
;2017/03/08 Here's the story that I'm sticking to:
;JOURNAL__20170224__ALL_THE_ORBS_WE_DONE_RECENTLY______PA_N_GEORGE,ORBIT=9627,/PARSE_B_AND_J_SAVEFILE,/USE_TIMEBAR_TIME__FROM_FILE,/FOLD_NEGFREQ_ONTO_POS,/PLOT_KPERP_MAGNITUDE_FOR_KZ,/PLOT_KX_VS_KY_FOR_KZ,/SAVE_PS,/USE_LOWRES_TIME_SERIES,/PUBLICATION_SETTINGS,BONUS_SUFF='SDFSDF',DATE=GET_TODAY_STRING(/DO_YYYYMMDD_FMT),/TO_PDF,KP__ANGLERANGE=[0,360],/PRE_VIII_LAYOUT,CUSTOM_T1='1999-01-27/11:32:44',CUSTOM_T2='1999-01-27/11:33:09',WHICH_FFTS=1,FFTPERCENT=50
;
;2017/03/09 Nah, my story always changes (use Rep[aired]Ret[race_removed]Cal[ibrated] file):
;JOURNAL__20170224__ALL_THE_ORBS_WE_DONE_RECENTLY______PA_N_GEORGE,ORBIT=9627,/PARSE_B_AND_J_SAVEFILE,/USE_TIMEBAR_TIME__FROM_FILE,/FOLD_NEGFREQ_ONTO_POS,/PLOT_KPERP_MAGNITUDE_FOR_KZ,/PLOT_KX_VS_KY_FOR_KZ,/SAVE_PS,/USE_LOWRES_TIME_SERIES,/PUBLICATION_SETTINGS,BONUS_SUFF='SDFSDF',DATE=GET_TODAY_STRING(/DO_YYYYMMDD_FMT),/TO_PDF,KP__ANGLERANGE=[0,360],/PRE_VIII_LAYOUT,CUSTOM_T1='1999-01-27/11:32:44',CUSTOM_T2='1999-01-27/11:33:09',WHICH_FFTS=1,FFTPERCENT=50,/USE_REPRETCAL_FILE
;
;Another keeper(?):
;JOURNAL__20170309__RUN_THE_OLD_20161022_FILE__WHY_WONT_ORB_9585_COOPERATE,/USE_ALL_STREAKS,/FOLD_NEGFREQ_ONTO_POS,/USE_LOWRES_TIME_SERIES,/SAVE_PS,/TO_PDF,/USE_TIMEBAR_TIME__FROM_FILE,KP__ANGLERANGE=[0,360],/PRE_VIII_LAYOUT,/PLOT_KX_VS_KY_FOR_KZ,/PUBLICATION_SETTINGS,/PLOT_SMOOTHED_K_COMPONENTS,/USE_DB_FAC,CUSTOM_T1='1999-01-23/14:50:37',CUSTOM_T2='1999-01-23/14:50:52'
PRO JOURNAL__20170224__ALL_THE_ORBS_WE_DONE_RECENTLY______PA_N_GEORGE, $
   USE_REPRETCAL_FILE=use_RepRetCal_file,$
   SKIP_DESPIN=skip_despin, $
   ORBIT=orbit, $
   DATE=date, $
   PARSE_B_AND_J_SAVEFILE=parse_B_and_J_saveFile, $
   B_AND_J_FILE=saveFile, $
   USE_TIMEBAR_TIME__FROM_FILE=use_timeBar_time__from_file, $
   CUSTOM_T1=custom_t1, $
   CUSTOM_T2=custom_t2, $
   CUSTOM_T_EXTEND=custom_t_extend, $
   CUSTOM_MULTI_T_EXTEND=custom_multi_t_extend, $
   CUSTOM_SHIFTSEC=custom_shiftSec, $
   CUSTOM_MULTI_SHIFTSEC=custom_multi_shiftSec, $
   PTSHIFT=ptShift, $
   SHIFT_UNIVERSAL_OFFSET=univOffset, $
   NOSHIFT=noShift, $
   NO_PTSHIFT_BACK=no_ptShift_back, $
   NO_PTSHIFT_FWD=no_ptShift_fwd, $
   SHIFT_NPTS=shift_nPts, $
   BACKSHIFTS_FOR_AVGING=backShifts_for_avging, $
   FWDSHIFTS_FOR_AVGING=fwdShifts_for_avging, $
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
   ONLY_KS_BELOW_ERRANGLE_THRESH=only_ks_below_errAngle_thresh, $
   ONLY_KS_BELOW_PHASE_ERR_THRESH=only_ks_below_phaseErr_thresh, $
   EXAMPLE_MODE=example_mode, $
   PLOT_KPERP_MAGNITUDE_FOR_KZ=plot_kperp_magnitude_for_kz, $
   PLOT_KX_VS_KY_FOR_KZ=plot_kx_vs_ky_for_kz, $
   PLOT_SMOOTHED_K_COMPONENTS=plot_smoothed_ks, $
   PLOT_ABS_SMOOTHED_K_COMPONENTS=plot_abs_smoothed_ks, $
   KX_VS_KY__PLOT_SMOOTHED=kx_vs_ky__plot_smoothed, $
   KP_ANGLE__PLOT_SMOOTHED=kP_angle__plot_smoothed, $
   PLOT_POSFREQ=plot_posFreq, $
   FOLD_NEGFREQ_ONTO_POS=fold_negFreq, $
   SAVE_PS=save_ps, $
   TO_PDF=to_pdf, $
   PDF_TRANSPARENCY_LEVEL=pdf_transparency, $
   REMOVE_EPS=remove_eps, $
   NO_PLOTS=no_plots, $
   SHOW_PREDICTED_J=show_predicted_J, $
   BONUS_SUFF=bonus_suff, $
   EXTRA_SUFFIX=extra_suffix, $
   DOUBLE_CALC=double_calc, $
   HANNING=hanning, $
   USE_LOWRES_TIME_SERIES=use_lowRes_time_series, $
   USE_J_TIME_SERIES=use_J_time_series, $
   SMOOTH_J_DAT_TO_B=smooth_J_dat, $
   PRESMOOTH_MAG=presmooth_mag, $
   KSMOOTH__NSMOOTHS=smInd, $
   KSMOOTH__DOUBLENSMOOTHS=dbSmInd, $
   KSMOOTH__EDGE_TRUNCATE=kSmooth__edge_truncate, $
   KSMOOTH__EDGE_MIRROR=kSmooth__edge_mirror, $
   KSMOOTH__EDGE_WRAP=kSmooth__edge_wrap, $
   FREQLIMS=freqLims, $
   PAGE1__FREQRANGE=page1__freqRange, $
   PAGE2__FREQRANGE=page2__freqRange, $
   KP__ANGLERANGE=kP__angleRange, $
   THIRD_PAGE=third_page, $
   OVERPLOT_DOUBLY_SMOOTHED=overplot_doubly_smoothed, $
   PREPLOT_CURRENTS_AND_STOP=prePlot_currents_and_stop, $
   FITLINE__USE_ABS=fitline__use_abs, $
   FITLINE__USE_SMOOTHED=fitline__use_smoothed, $
   COMBINE_AND_AVERAGE_INTERVALS=combine_and_average_intervals, $
   USE_DB_FAC=use_dB_fac, $
   STREAKNUM=streakNum, $
   USE_ALL_STREAKS=use_all_streaks, $
   PUBLICATION_SETTINGS=pubSettings, $
   PRE_VIII_LAYOUT=PRE_VIII_layout, $
   FOOTBALL_LAYOUT=football_layout, $
   FOOTBALL_YLOG=football_yLog, $
   FOOTBALL_COL2TITLE=football_col2Title, $
   FOOTBALL_KMAG=football_kMag, $
   FOOTBALL_NO_ERRANGLE=football_no_errAngle, $
   FOOTBALL_NO_MAGERR=football_no_magErr, $
   ODDNESS_CHECK=oddness_check, $
   FFT__NEAREST_TWO_POWER=nearest_two_power, $
   FFTSIZE=FFTsize, $
   FFTPERCENT=FFTpercent,$
   WHICH_FFTS=which_FFTs, $
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

  COMPILE_OPT IDL2,STRICTARRSUBS

  IF ~KEYWORD_SET(orbit) THEN BEGIN
     PRINT,"Possibilities"
     PRINT," 10767"
     PRINT," 10832"
     PRINT," 10837"
     PRINT," 10839"
     PRINT," 9627 "
     PRINT," 10927"
     PRINT," 9585 (the original!)"
     PRINT," 6717 (Look out--FG was only sampling at 8 Hz)"
     RETURN
  ENDIF

  IF N_ELEMENTS(date) EQ 0 THEN BEGIN
     date    = '20170224'
  ENDIF

  eeb_or_ees = 'eeb'
  ieb_or_ies = 'ieb'

  saveSuff   = '-with_sc_pot'

  SUMPLOTS_AND_B_PLUS_J__GET_FILENAME, $
     ORBIT=orbit, $
     DATE=date, $
     EEB_OR_EES=eeb_or_ees, $
     IEB_OR_IES=ieb_or_ies, $
     SKIP_DESPIN=skip_despin, $
     PLOTPREF=plotPref, $
     SAVESUFF=saveSuff, $
     BONUSSUFF=bonusSuff, $
     USE_REPRETCAL_FILE=use_RepRetCal_file,$
     ANCILLARY_PLOTS=ancillary_plots, $
     OUTPLOTNAME=outPlotName, $
     SAVEFILE=saveFile

  extra_suffix = saveFile.Replace('.sav','')
  
  SINGLE_SPACECRAFT_K_MEASUREMENT_FAST, $
     PARSE_B_AND_J_SAVEFILE=parse_B_and_J_saveFile, $
     B_AND_J_FILE=saveFile, $
     USE_TIMEBAR_TIME__FROM_FILE=use_timeBar_time__from_file, $
     CUSTOM_T1=custom_t1, $
     CUSTOM_T2=custom_t2, $
     CUSTOM_T_EXTEND=custom_t_extend, $
     CUSTOM_MULTI_T_EXTEND=custom_multi_t_extend, $
     CUSTOM_SHIFTSEC=custom_shiftSec, $
     CUSTOM_MULTI_SHIFTSEC=custom_multi_shiftSec, $
     PTSHIFT=ptShift, $
     SHIFT_UNIVERSAL_OFFSET=univOffset, $
     NOSHIFT=noShift, $
     NO_PTSHIFT_BACK=no_ptShift_back, $
     NO_PTSHIFT_FWD=no_ptShift_fwd, $
     SHIFT_NPTS=shift_nPts, $
     BACKSHIFTS_FOR_AVGING=backShifts_for_avging, $
     FWDSHIFTS_FOR_AVGING=fwdShifts_for_avging, $
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
     EXAMPLE_MODE=example_mode, $
     PLOT_KPERP_MAGNITUDE_FOR_KZ=plot_kperp_magnitude_for_kz, $
     PLOT_KX_VS_KY_FOR_KZ=plot_kx_vs_ky_for_kz, $
     PLOT_SMOOTHED_K_COMPONENTS=plot_smoothed_ks, $
     PLOT_ABS_SMOOTHED_K_COMPONENTS=plot_abs_smoothed_ks, $
     KX_VS_KY__PLOT_SMOOTHED=kx_vs_ky__plot_smoothed, $
     KP_ANGLE__PLOT_SMOOTHED=kP_angle__plot_smoothed, $
     PLOT_POSFREQ=plot_posFreq, $
     FOLD_NEGFREQ_ONTO_POS=fold_negFreq, $
     SAVE_PS=save_ps, $
     TO_PDF=to_pdf, $
     PDF_TRANSPARENCY_LEVEL=pdf_transparency, $
     REMOVE_EPS=remove_eps, $
     NO_PLOTS=no_plots, $
     SHOW_PREDICTED_J=show_predicted_J, $
     BONUS_SUFF=bonus_suff, $
     EXTRA_SUFFIX=extra_suffix, $
     DOUBLE_CALC=double_calc, $
     HANNING=hanning, $
     USE_LOWRES_TIME_SERIES=use_lowRes_time_series, $
     USE_J_TIME_SERIES=use_J_time_series, $
     SMOOTH_J_DAT_TO_B=smooth_J_dat, $
     PRESMOOTH_MAG=presmooth_mag, $
     KSMOOTH__NSMOOTHS=smInd, $
     KSMOOTH__DOUBLENSMOOTHS=dbSmInd, $
     KSMOOTH__EDGE_TRUNCATE=kSmooth__edge_truncate, $
     KSMOOTH__EDGE_MIRROR=kSmooth__edge_mirror, $
     KSMOOTH__EDGE_WRAP=kSmooth__edge_wrap, $
     FREQLIMS=freqLims, $
     PAGE1__FREQRANGE=page1__freqRange, $
     PAGE2__FREQRANGE=page2__freqRange, $
     KP__ANGLERANGE=kP__angleRange, $
     THIRD_PAGE=third_page, $
     OVERPLOT_DOUBLY_SMOOTHED=overplot_doubly_smoothed, $
     PREPLOT_CURRENTS_AND_STOP=prePlot_currents_and_stop, $
     FITLINE__USE_ABS=fitline__use_abs, $
     FITLINE__USE_SMOOTHED=fitline__use_smoothed, $
     COMBINE_AND_AVERAGE_INTERVALS=combine_and_average_intervals, $
     USE_DB_FAC=use_dB_fac, $
     STREAKNUM=streakNum, $
     USE_ALL_STREAKS=use_all_streaks, $
     PUBLICATION_SETTINGS=pubSettings, $
     PRE_VIII_LAYOUT=PRE_VIII_layout, $
     FOOTBALL_LAYOUT=football_layout, $
     FOOTBALL_YLOG=football_yLog, $
     FOOTBALL_COL2TITLE=football_col2Title, $
     FOOTBALL_KMAG=football_kMag, $
     FOOTBALL_NO_ERRANGLE=football_no_errAngle, $
     FOOTBALL_NO_MAGERR=football_no_magErr, $
     ODDNESS_CHECK=oddness_check, $
     FFT__NEAREST_TWO_POWER=nearest_two_power, $
     FFTSIZE=FFTsize, $
     FFTPERCENT=FFTpercent,$
     WHICH_FFTS=which_FFTs, $
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

END
