;;2017/02/24
PRO JOURNAL__20170224__PRESERVE_THE_STUFF_FROM_ANALYSIS_OF_ORBIT_9585_FOR_PREVIII

  COMPILE_OPT IDL2

  parse_B_and_J_saveFile

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;The original Chaston et al. [2006] interval

  ;; saveFile = 'Chaston_et_al_2006--B_and_J.sav'
  ;; saveFile = 'Chaston_et_al_2006--B_and_J--20161022--fixed_currents.sav'
  ;; saveFile = 'Chaston_et_al_2006--B_and_J--20161022--fixed_currents_2.sav'
  ;; saveFile = 'Chaston_et_al_2006--B_and_J--20161022--fixed_currents--with_sc_pot.sav'
  ;; IF KEYWORD_SET(use_timeBar_time) THEN BEGIN
  ;;    ;; timeBar_times = ['1998-05-04/06:44:31.5','1998-05-04/06:44:56.5']
  ;;    timeBar_times = ['1998-05-04/06:44:31.5','1998-05-04/06:44:56.5']
  ;;    ;; timeBar_times = ['1998-05-04/06:44:36','1998-05-04/06:44:56']
  ;;    ;; timeBar_times = ['1998-05-04/06:44:46','1998-05-04/06:44:56']

  ;;    ;; timeBar_times = [ ['1998-05-04/06:44:21.5','1998-05-04/06:44:29.0'], $
  ;;    ;;                 ['1998-05-04/06:44:31.5','1998-05-04/06:44:56.5'] ]
  ;; ENDIF
  ;; extra_suffix = 'dude'

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;The orbit 9585
  saveFile = 'Orbit_9585--B_and_J--20161024--fixed_currents--with_sc_pot.sav'
  IF KEYWORD_SET(use_timeBar_time) THEN BEGIN
     timeBar_times = ['1999-01-23/14:50:56','1999-01-23/14:51:06']
  ENDIF
  extra_suffix = 'Orbit_9585'

  ;; saveFile = 'Orbit_9585--B_and_J--20161024--fixed_currents--with_sc_pot--wider.sav'
  ;; IF KEYWORD_SET(use_timeBar_time) THEN BEGIN
  ;;    timeBar_times = ['1999-01-23/14:50:52','1999-01-23/14:51:03']
  ;; ENDIF
  ;; extra_suffix = 'Orbit_9585--wider_tBar'

  ;; saveFile = 'Orbit_9585--B_and_J--20161024--fixed_currents--with_sc_pot--bro.sav'
  ;; IF KEYWORD_SET(use_timeBar_time) THEN BEGIN
  ;;    timeBar_times = ['1999-01-23/14:50:52','1999-01-23/14:51:03']
  ;; ENDIF
  ;; extra_suffix = 'Orbit_9585--wider_tBar--PRE_VIII'

  ;; saveFile = 'Orbit_9585--B_and_J--20161024--fixed_currents--with_sc_pot--bro.sav--alt_timebar'
  ;; IF KEYWORD_SET(use_timeBar_time) THEN BEGIN
  ;;    timeBar_times = ['1999-01-23/14:50:40','1999-01-23/14:50:51.5']

     ;; timeBar_times    = ['1999-01-23/14:50:40','1999-01-23/14:50:46.25']
     ;; extra_suffix = 'Orbit_9585--wider_alt_tBar--PRE_VIII--46p25_to_51.5'

     ;; timeBar_times = ['1999-01-23/14:50:46.25','1999-01-23/14:50:51.5']
     ;; extra_suffix = 'Orbit_9585--wider_alt_tBar--PRE_VIII--46p25_to_51.5'

     ;;Money #1
     ;; timeBar_times = ['1999-01-23/14:50:46.25','1999-01-23/14:50:51.5']
     ;; extra_suffix = 'Orbit_9585--wider_alt_tBar--PRE_VIII--46p25_to_51.5'

     ;;Possibly Money #2
     ;; timeBar_times = ['1999-01-23/14:51:02.0','1999-01-23/14:51:06.0']
     ;; extra_suffix = 'Orbit_9585--wider_alt_tBar--PRE_VIII--46p25_to_51.5'

     ;; timeBar_times = ['1999-01-23/14:50:42.0','1999-01-23/14:50:43.0']
     ;; extra_suffix = 'Orbit_9585--wider_alt_tBar--PRE_VIII--46p25_to_51.5'

     ;; timeBar_times = ['1999-01-23/14:50:52','1999-01-23/14:50:57.5']
     ;; extra_suffix = 'Orbit_9585--wider_alt_tBar--PRE_VIII--52_to_57p5'

     ;; timeBar_times = ['1999-01-23/14:50:57.5','1999-01-23/14:51:03']
     ;; extra_suffix = 'Orbit_9585--wider_alt_tBar--PRE_VIII--57p5_to_03'

     ;; timeBar_times = ['1999-01-23/14:50:44','1999-01-23/14:50:49']
     ;; extra_suffix = 'Orbit_9585--wider_alt_tBar--PRE_VIII--44_to_49'

     ;; timeBar_times = ['1999-01-23/14:50:36','1999-01-23/14:50:40']
     ;; extra_suffix = 'Orbit_9585--wider_alt_tBar--PRE_VIII--36_to_40'

     ;; timeBar_times = ['1999-01-23/14:50:40','1999-01-23/14:51:03']
  ;; ENDIF

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;The orbit 10837
  ;; saveFile = 'Orbit_10837--B_and_J--20161025--fixed_currents--with_sc_pot.sav'
  ;; saveFile = 'Orbit_10837--B_and_J--20161025--fixed_currents--with_sc_pot.sav'

  ;; IF KEYWORD_SET(use_timeBar_time) THEN BEGIN
  ;;    timeBar_times = ['1999-05-18/17:46:24','1999-05-18/17:46:39']
  ;; ENDIF
  ;; extra_suffix = 'Orbit_10837'

  ;; IF KEYWORD_SET(use_timeBar_time) THEN BEGIN
  ;;    ;; timeBar_times = ['1999-05-18/17:46:24','1999-05-18/17:46:39']
  ;;    ;; timeBar_times = ['1999-05-18/17:46:20','1999-05-18/17:46:39']
  ;;    ;; timeBar_times = ['1999-05-18/17:46:18','1999-05-18/17:46:40']
  ;;    ;; timeBar_times = ['1999-05-18/17:46:12','1999-05-18/17:46:46']
  ;;    timeBar_times = ['1999-05-18/17:46:10','1999-05-18/17:46:48']
  ;; ENDIF
  ;; extra_suffix = 'Orbit_10837--poke_around'

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;The orbit 10927
  ;; saveFile = 'Orbit_10927--B_and_J--20161025--fixed_currents--with_sc_pot.sav'
  ;; IF KEYWORD_SET(use_timeBar_time) THEN BEGIN
  ;;    timeBar_times = ['1999-05-27/00:30:30','1999-05-27/00:30:47']
  ;; ENDIF
  ;; extra_suffix = 'Orbit_10927'

  IF KEYWORD_SET(timeBar_times) THEN BEGIN

     CASE NDIMEN(timeBar_times) OF
        1: BEGIN
           CASE SIZE(timeBar_times,/TYPE) OF
              7: BEGIN
                 custom_t1  = STR_TO_TIME(timeBar_times[0])
                 custom_t2  = STR_TO_TIME(timeBar_times[1])
              END
              5: BEGIN
                 custom_t1  = timeBar_times[0]
                 custom_t2  = timeBar_times[1]
              END
              ELSE: STOP
           ENDCASE
        END
        2: BEGIN
           CASE SIZE(timeBar_times,/TYPE) OF
              7: BEGIN
                 custom_t1  = STR_TO_TIME(REFORM(timeBar_times[0,*]))
                 custom_t2  = STR_TO_TIME(REFORM(timeBar_times[1,*]))
              END
              5: BEGIN
                 custom_t1  = REFORM(timeBar_times[0,*])
                 custom_t2  = REFORM(timeBar_times[1,*])
              END
              ELSE: STOP
           ENDCASE
        END
     ENDCASE

  ENDIF

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
   PLOT_POSFREQ=plot_posFreq, $
   FOLD_NEGFREQ_ONTO_POS=fold_negFreq, $
   SAVE_PS=save_ps, $
   BONUS_SUFF=bonus_suff, $
   EXTRA_SUFFIX=extra_suffix, $
   DOUBLE_CALC=double_calc, $
   HANNING=hanning, $
   USE_J_TIME_SERIES=use_J_time_series, $
   SMOOTH_J_DAT_TO_B=smooth_J_dat, $
   PRESMOOTH_MAG=presmooth_mag, $
   KSMOOTH__NSMOOTHS=smInd, $
   KSMOOTH__DOUBLENSMOOTHS=dbSmInd, $
   KSMOOTH__EDGE_TRUNCATE=kSmooth__edge_truncate, $
   KSMOOTH__EDGE_MIRROR=kSmooth__edge_mirror, $
   KSMOOTH__EDGE_WRAP=kSmooth__edge_wrap, $
   OVERPLOT_DOUBLY_SMOOTHED=overplot_doubly_smoothed, $
   PREPLOT_CURRENTS_AND_STOP=prePlot_currents_and_stop, $
   FITLINE__USE_ABS=fitline__use_abs, $
   FITLINE__USE_SMOOTHED=fitline__use_smoothed, $
   COMBINE_AND_AVERAGE_INTERVALS=combine_and_average_intervals, $
   USE_DB_FAC=use_dB_fac, $
   STREAKNUM=streakNum, $
   USE_ALL_STREAKS=use_all_streaks, $
   PUBLICATION_SETTINGS=pubSettings, $
   ODDNESS_CHECK=oddness_check

END
