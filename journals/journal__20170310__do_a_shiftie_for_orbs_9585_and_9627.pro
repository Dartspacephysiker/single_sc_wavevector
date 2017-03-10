;2017/03/10
PRO JOURNAL__20170310__DO_A_SHIFTIE_FOR_ORBS_9585_AND_9627__OBSOLETE, $
   PLOT_SMOOTHED_K_COMPONENTS=plot_smoothed_k_components, $
   SAVE_PS=save_ps, $
   NO_PLOTS=no_plots, $
   PUBLICATION_SETTINGS=pubSettings, $
   USE_DB_FAC=use_dB_fac


  COMPILE_OPT IDL2,STRICTARRSUBS

  minShiftBack  = -21
  maxShiftBack  = 9
  stepShiftBack = 6

  minShiftFwd   = -9
  maxShiftFwd   = 21
  stepShiftFwd  = 6

  nBack         = (maxShiftBack-minShiftBack)/stepShiftBack + 1
  nFwd          = (maxShiftFwd -minShiftFwd )/stepShiftFwd + 1
  backShifts    = INDGEN(nBack)*stepShiftBack+minShiftBack
  fwdShifts     = INDGEN(nFwd )*stepShiftFwd+minShiftFwd

  junk          = MIN(ABS(backShifts),minJJInd)
  junk          = MIN(ABS(fwdShifts),minKKInd)

  count         = 0
  TArrList      = LIST()
  freqList      = LIST()
  kxList        = LIST()
  kyList        = LIST()
  kzList        = LIST()
  kPList        = LIST()
  kPAngleList   = LIST()
  indsList      = LIST()
  kzList        = LIST()
  kzList        = LIST()
  kzList        = LIST()
  FOR jj=0,nBack-1 DO BEGIN
     FOR kk=0,nFwd-1 DO BEGIN

        tmpSuff = STRING(FORMAT='("__shft_",I0,"_",I0)',backShifts[jj],fwdShifts[kk])
        PRINT,FORMAT='(I0,T10,A0)',count,tmpSuff

        JOURNAL__20170224__ALL_THE_ORBS_WE_DONE_RECENTLY______PA_N_GEORGE, $
           ORBIT=9627, $
           /PARSE_B_AND_J_SAVEFILE, $
           /USE_TIMEBAR_TIME__FROM_FILE, $
           /FOLD_NEGFREQ_ONTO_POS, $
           /PLOT_KPERP_MAGNITUDE_FOR_KZ, $
           /PLOT_KX_VS_KY_FOR_KZ, $
           SAVE_PS=save_ps, $
           TO_PDF=KEYWORD_SET(save_ps), $
           REMOVE_EPS=KEYWORD_SET(save_ps), $
           NO_PLOTS=no_plots, $
           /USE_LOWRES_TIME_SERIES, $
           DATE='20170309', $
           PUBLICATION_SETTINGS=pubSettings, $
           BONUS_SUFF=tmpSuff, $
           USE_DB_FAC=use_dB_fac, $
           CUSTOM_T1='1999-01-27/11:32:57.542', $
           CUSTOM_T2='1999-01-27/11:33:09', $
           /USE_REPRETCAL_FILE, $
           KP__ANGLERANGE=[-90,270], $
           PRE_VIII_LAYOUT=KEYWORD_SET(save_ps), $
           PLOT_SMOOTHED_K_COMPONENTS=plot_smoothed_k_components, $
           BACKSHIFTS_FOR_AVGING=backShifts_for_avging, $
           FWDSHIFTS_FOR_AVGING=fwdShifts_for_avging, $
           SHIFT_NPTS=[backShifts[jj],fwdShifts[kk]], $
           OUT_FREQS=out_freqs, $
           OUT_KX=out_kx, $
           OUT_KY=out_ky, $
           OUT_KZ=out_kz, $
           OUT_KP=out_kP, $
           OUT_KPANGLE=out_kPAngle, $
           ;; OUT_INDS=out_inds, $
           OUT_TARR=out_TArr, $
           OUT_USEDINDS=out_usedInds, $
           OUT_BX=out_Bx, $
           OUT_BY=out_By, $
           OUT_BZ=out_Bz, $
           OUT_JX=out_Jx, $
           OUT_JY=out_Jy, $
           OUT_JZ=out_Jz

        count++

        freqList.Add,out_freqs
        kxList.Add,out_kx
        kyList.Add,out_ky
        kzList.Add,out_kz
        kPList.Add,out_kP
        kPAngleList.Add,out_kPAngle
        ;; indsList.Add,out_inds
        kzList.Add,out_kz
        kzList.Add,out_kz

        IF (kk EQ minKKInd) AND (jj EQ minJJInd) THEN BEGIN

           TArr      = out_TArr
           usedInds  = out_usedInds
           Bx        = out_Bx
           By        = out_By
           Bz        = out_Bz
           Jx        = out_Jx
           Jy        = out_Jy
           Jz        = out_Jz

        ENDIF

     ENDFOR
  ENDFOR

  freqs  = LIST_TO_1DARRAY(freqList,/WARN)
  kxs    = LIST_TO_1DARRAY(kxList,/WARN)
  kys    = LIST_TO_1DARRAY(kyList,/WARN)
  kzs    = LIST_TO_1DARRAY(kzList,/WARN)
  kPs    = LIST_TO_1DARRAY(kPList,/WARN)
  kPAngles    = LIST_TO_1DARRAY(kPAngleList,/WARN)
  ;; indss    = LIST_TO_1DARRAY(indsList,/WARN)
  ;; kzs    = LIST_TO_1DARRAY(kzList,/WARN)
  ;; kzs    = LIST_TO_1DARRAY(kzList,/WARN)
  ;; kzs    = LIST_TO_1DARRAY(kzList,/WARN)

  binSz  = 0.2
  kxHist = HIST1D(freqs,kxs,BINSIZE=binSz,OBIN=oFreqsX)/count
  kyHist = HIST1D(freqs,kys,BINSIZE=binSz,OBIN=oFreqsY)/count
  kzHist = HIST1D(freqs,kzs,BINSIZE=binSz,OBIN=oFreqsY)/count
  kPHist = HIST1D(freqs,kPs,BINSIZE=binSz,OBIN=oFreqsY)/count
  kPAngleHist = HIST1D(freqs,kPAngles,BINSIZE=binSz,OBIN=oFreqsY)/count

  !P.MULTI = [0,1,2,0,0]

  kx_ysize = MAX(ABS(kxHist))
  ky_ysize = MAX(ABS(kyHist))

  SET_PLOT,'X'
  WINDOW,0,XSIZE=800,YSIZE=800

  PLOT,oFreqsX,kxHist, $
       ;; YTITLE='k!Dx!N (m!U-1!N)', $
       YTITLE='k!Dx!N (km!U-1!N)', $
       XTITLE='', $
       XRANGE=page1__freqRange, $
       YRANGE=[-kx_ysize,kx_ysize], $
       XSTYLE=1, $
       YSTYLE=1, $
       XTICKLEN=1.0, $
       YTICKLEN=1.0, $
       XGRIDSTYLE=1, $
       YGRIDSTYLE=1, $
       CHARSIZE=cs

  PLOT,oFreqsY,kYHist, $
       ;; YTITLE='k!Dy!N (m!U-1!N)', $
       YTITLE='k!Dy!N (km!U-1!N)', $
       XTITLE='', $
       XRANGE=page1__freqRange, $
       YRANGE=[-ky_ysize,ky_ysize], $
       XSTYLE=1, $
       YSTYLE=1, $
       XTICKLEN=1.0, $
       YTICKLEN=1.0, $
       XGRIDSTYLE=1, $
       YGRIDSTYLE=1, $
       CHARSIZE=cs


  STOP
END
