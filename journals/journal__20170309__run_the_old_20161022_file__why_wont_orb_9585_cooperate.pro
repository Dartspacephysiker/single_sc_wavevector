;2017/03/09
;Here's the command—just uncomment the file that you want to run
;JOURNAL__20170309__RUN_THE_OLD_20161022_FILE__WHY_WONT_ORB_9585_COOPERATE,/USE_ALL_STREAKS,/FOLD_NEGFREQ_ONTO_POS,/USE_LOWRES_TIME_SERIES,/SAVE_PS,/TO_PDF,/USE_TIMEBAR_TIME__FROM_FILE,KP__ANGLERANGE=[0,360],/PRE_VIII_LAYOUT,/PLOT_KX_VS_KY_FOR_KZ,/PUBLICATION_SETTINGS,/PLOT_SMOOTHED_K_COMPONENTS
PRO JOURNAL__20170309__RUN_THE_OLD_20161022_FILE__WHY_WONT_ORB_9585_COOPERATE, $
   PARSE_B_AND_J_SAVEFILE=parse_B_and_J_saveFile, $
   B_AND_J_FILE=saveFile, $
   USE_TIMEBAR_TIME__FROM_FILE=use_timeBar_time__from_file, $
   CUSTOM_T1=custom_t1, $
   CUSTOM_T2=custom_t2, $
   EXAMPLE_MODE=example_mode, $
   PLOT_KPERP_MAGNITUDE_FOR_KZ=plot_kperp_magnitude_for_kz, $
   PLOT_KX_VS_KY_FOR_KZ=plot_kx_vs_ky_for_kz, $
   PLOT_SMOOTHED_K_COMPONENTS=plot_smoothed_ks, $
   PLOT_ABS_SMOOTHED_K_COMPONENTS=plot_abs_smoothed_ks, $
   KX_VS_KY__PLOT_SMOOTHED=kx_vs_ky__plot_smoothed, $
   KP_ANGLE____PLOT_SMOOTHED=kP_angle__plot_smoothed, $
   PLOT_POSFREQ=plot_posFreq, $
   FOLD_NEGFREQ_ONTO_POS=fold_negFreq, $
   SAVE_PS=save_ps, $
   TO_PDF=to_pdf, $
   PDF_TRANSPARENCY_LEVEL=pdf_transparency, $
   REMOVE_EPS=remove_eps, $
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
   ODDNESS_CHECK=oddness_check, $
   FFTSIZE=FFTsize, $
   FFTPERCENT=FFTpercent,$
   WHICH_FFTS=which_FFTs

  COMPILE_OPT IDL2,STRICTARRSUBS


  parse_B_and_J_saveFile = 1

  ;; saveFile               = 'Chaston_et_al_2006--B_and_J--20161022--fixed_currents.sav'
  ;; saveFile               = 'Orbit_9585--B_and_J--20161024--fixed_currents--with_sc_pot--bro.sav--alt_timebar'
  ;; bonus_suff             = '-try_to_reproduce_PRE_VIII_orig_result'

  saveFile               = 'Orbit_9585-B_and_J-20170309-eeb-ieb-with_sc_pot.sav-20170225journal__1minLoaded'
  ;; saveFile               = 'Orbit_9585-B_and_J-20170309-eeb-ieb-with_sc_pot.sav-20170225journal__50minLoaded'

  bonus_suff             = '-' + (STRSPLIT(saveFile,'__',/EXTRACT))[-1] + '_in_SDT'

  SINGLE_SPACECRAFT_K_MEASUREMENT_FAST, $
     PARSE_B_AND_J_SAVEFILE=parse_B_and_J_saveFile, $
     B_AND_J_FILE=saveFile, $
     USE_TIMEBAR_TIME__FROM_FILE=use_timeBar_time__from_file, $
     CUSTOM_T1=custom_t1, $
     CUSTOM_T2=custom_t2, $
     EXAMPLE_MODE=example_mode, $
     PLOT_KPERP_MAGNITUDE_FOR_KZ=plot_kperp_magnitude_for_kz, $
     PLOT_KX_VS_KY_FOR_KZ=plot_kx_vs_ky_for_kz, $
     PLOT_SMOOTHED_K_COMPONENTS=plot_smoothed_ks, $
     PLOT_ABS_SMOOTHED_K_COMPONENTS=plot_abs_smoothed_ks, $
     KX_VS_KY__PLOT_SMOOTHED=kx_vs_ky__plot_smoothed, $
     KP_ANGLE____PLOT_SMOOTHED=kP_angle__plot_smoothed, $
     PLOT_POSFREQ=plot_posFreq, $
     FOLD_NEGFREQ_ONTO_POS=fold_negFreq, $
     SAVE_PS=save_ps, $
     TO_PDF=to_pdf, $
     PDF_TRANSPARENCY_LEVEL=pdf_transparency, $
     REMOVE_EPS=remove_eps, $
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
     ODDNESS_CHECK=oddness_check, $
     FFTSIZE=FFTsize, $
     FFTPERCENT=FFTpercent,$
     WHICH_FFTS=which_FFTs

  
END
