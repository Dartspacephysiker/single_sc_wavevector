;2017/03/10
;********************************************************************
; written by P. M. Bellan, April 2016, MS 128-95 Caltech, Pasadena CA 91125
; prefix gives file path for eps output, data input/output
; user should change to be appropriate for user's computer
PRO OUTPUT_SETUP,mode,plotDir,suff

  COMPILE_OPT idl2,strictarrsubs

  COMMON PORIG,pOrig_cThick,pOrig_thick
  

  ;; prefix='D:\CORR\MNSCRPTS\2016-JGR-spacecraft-current-wavevector\2106-Vinas-kvect-revised\'

  filename = plotDir+ $
             suff

  PRINT,STRING(FORMAT='("plotFile : ",A0)',suff)
  ;;SET UP FOR PLOTTING ON:
  ;; SCREEN,          set hardcopy = 0
  ;; PRINTER,         set hardcopy = 1
  ;; POSTSCRIPT FILE, set hardcopy = 2
  hardcopy = mode
  IF hardcopy EQ 1 THEN SET_PLOT,'printer' ;PRINTER
  IF hardcopy EQ 1 THEN DEVICE,YSIZE=25,YOFFSET=0
  IF hardcopy EQ 0 THEN BEGIN
     SET_PLOT,'X'
     WINDOW,0,XSIZE=800,YSIZE=800
  ENDIF
  IF hardcopy EQ 2 THEN BEGIN      ;;POSTSCRIPT FILE
     ;; !P.CHARSIZE  = 3
     pOrig_cThick = !P.CHARTHICK
     pOrig_thick  = !P.THICK
     
     !P.CHARTHICK = 1.2
     !P.THICK     = 1.0

     ;; @startup
     ;; POPEN,filename,XSIZE=10,YSIZE=10

     SET_PLOT,'PS'
     ;; DEVICE,FILE=filename+'.eps',/ENCAPSUL,XSIZE=10,YSIZE=10,/INCHES,YOFFSET=2,/COLOR
     DEVICE,FILE=filename+'.eps',XSIZE=5,YSIZE=5,/INCHES,YOFFSET=2,/COLOR
  ENDIF

END

;********************************************************************
; written by P. M. Bellan, April 2016, MS 128-95 Caltech, Pasadena CA 91125
; place at end of code to do housekeeping
PRO CONCLUDE_OUTPUT,mode,plotDir,suff, $
                    TO_PDF=to_pdf, $
                    PDF_TRANSPARENCY_LEVEL=pdf_transparency, $
                    REMOVE_EPS=remove_eps

  COMPILE_OPT idl2,strictarrsubs

  COMMON PORIG,pOrig_cThick,pOrig_thick

  hardcopy = mode
  IF hardcopy EQ 1 THEN DEVICE,/CLOSE
  IF hardcopy EQ 1 THEN SET_PLOT,'X'
  IF hardcopy EQ 2 THEN BEGIN
     DEVICE,/CLOSE
     ;; PCLOSE
     SET_PLOT,'X'

     !P.THICK     = pOrig_thick
     !P.CHARTHICK = pOrig_cThick

     IF KEYWORD_SET(to_pdf) THEN BEGIN
        filename = plotDir+ $
                   suff

        EPS2PDF,filename, $
                ;; /PS, $
                TRANSPARENCY_LEVEL=pdf_transparency, $
                REMOVE_EPS=remove_eps, $
                /QUIET
     ENDIF
  ENDIF
  IF hardcopy EQ 1 THEN SET_PLOT,'printer'
  IF hardcopy EQ 1 THEN DEVICE,YSIZE=25,YOFFSET=5
  IF hardcopy NE 1 THEN SET_PLOT,'X'
  ;; PRINT, 'FINISHED'
END

PRO PLOT_SPECIAL_FREQS,freq,k_specialFreqs,k_ysize,k_totSpecial_i

  COMPILE_OPT IDL2,STRICTARRSUBS

     k_totSpecial_i = !NULL
     FOR jD=0,N_ELEMENTS(k_specialFreqs)/2-1 DO BEGIN

        tmpSpFreq = k_specialFreqs[*,jD]

        ;;Special markers
        line1_freq = MAKE_ARRAY(2,VALUE=tmpSpFreq[0])
        line1_kx   = [-k_ysize,k_ysize]
        line2_freq = MAKE_ARRAY(2,VALUE=tmpSpFreq[1])
        line2_kx   = [-k_ysize,k_ysize]

        ;;Plot a special line, if there's any helping it
        special_i  = WHERE(freq GE tmpSpFreq[0] AND freq LE tmpSpFreq[1],nSpecial)
        ;; IF nSpecial GT 0 THEN BEGIN
        ;;    OPLOT,freq[special_i],kx[special_i],COLOR=specialColor
        ;; ENDIF

        ;;Now plot all those special people
        OPLOT,line1_freq,line1_kx,LINESTYLE=specialLineSty,COLOR=specialColor,THICK=specialThick
        OPLOT,line2_freq,line2_kx,LINESTYLE=specialLineSty,COLOR=specialColor,THICK=specialThick

        k_totSpecial_i = [k_totSpecial_i,special_i]
     ENDFOR



END

PRO MAKING_KX_KY_SPECIAL,kx_specialBounds,ky_specialBounds, $
                         kx_totSpecial_i,ky_totSpecial_i, $
                         OUT_KXY_SPI=kxy_spI, $
                         OUT_KPA_SPI=kPA_spI

  COMPILE_OPT IDL2,STRICTARRSUBS

  specialB_i  = !NULL
  IF KEYWORD_SET(kx_specialBounds) THEN BEGIN

     FOR jD=0,N_ELEMENTS(kx_specialBounds)/2-1 DO BEGIN

        tmpSpB     = kx_specialBounds[*,jD]
        
        special_i  = WHERE(kx GE tmpSpB[0] AND kx LE tmpSpB[1],nSpecial,/NULL)

        specialB_i = [specialB_i,special_i]
     ENDFOR

  ENDIF

  IF KEYWORD_SET(ky_specialBounds) THEN BEGIN

     FOR jD=0,N_ELEMENTS(ky_specialBounds)/2-1 DO BEGIN

        tmpSpB     = ky_specialBounds[*,jD]
        
        special_i  = WHERE(ky GE tmpSpB[0] AND ky LE tmpSpB[1],nSpecial,/NULL)

        specialB_i = [specialB_i,special_i]
     ENDFOR

  ENDIF

  CASE 1 OF
     KEYWORD_SET(kx_totSpecial_i) AND KEYWORD_SET(ky_totSpecial_i): BEGIN
        kxvky_spI = CGSETUNION(kx_totSpecial_i,ky_totSpecial_i)
     END
     KEYWORD_SET(kx_totSpecial_i): BEGIN
        kxvky_spI = kx_totSpecial_i
     END
     KEYWORD_SET(kx_totSpecial_i): BEGIN
        kxvky_spI = ky_totSpecial_i
     END
     ELSE:
  ENDCASE

  kxvky_spF_spB_i = kxvky_spI

  IF N_ELEMENTS(specialB_i) NE 0 THEN BEGIN

     specialB_i       = specialB_i[UNIQ(specialB_i,SORT(specialB_i))]
     kxvky_spF_spB_i  = CGSETINTERSECTION(kxvky_spF_spB_i,specialB_i,COUNT=nkxvky)

     IF nkxvky EQ 0 THEN STOP

  ENDIF

  kxy_spI         = kxvky_spF_spB_i

  ;;For special angles
  kPA_spI         = kxvky_spI

END

PRO GET_FREQ_INDS,freq,kx,ky,kz, $
                  kP,kPAngle, $
                  inds, $
                  BSPEC=BSpec, $
                  JSPEC=JSpec, $
                  MAGCSPEC=magCSpec, $
                  POWFREQ=powFreq, $
                  MAGERR=magErr, $
                  ERRANGLE=errAngle, $
                  OUT_FITINDS=fitInds, $
                  MAXFREQ=maxFreq, $
                  PLOT_POSFREQ=plot_posFreq, $
                  FOLD_NEGFREQ_ONTO_POS=fold_negFreq, $
                  FREQLIMS=freqLims

  COMPILE_OPT IDL2,STRICTARRSUBS

  CASE 1 OF
     KEYWORD_SET(plot_posFreq): BEGIN

        inds               = WHERE(freq GE 0.0)
        ;; fitInds            = WHERE(freq GT (0.1 < (freq[inds])[1]) AND freq LE maxFreq*1.1)

     END
     KEYWORD_SET(fold_negFreq): BEGIN
        
        ;; indNeg             = [0:((MAX(WHERE(freq LT 0.00)))-1)] 
        indNeg             = [0:((MAX(WHERE(freq LT 0.00))))] 
        indPos             = ([(WHERE(freq GE 0.00))[0]:(N_ELEMENTS(freq)-1)])[0:(N_ELEMENTS(indNeg)-1)]

        divFactor          = 2.

        adj                = freq[indPos[0]] EQ 0.00

        chopSiste          = 1

        IF N_ELEMENTS(kx) GT 0 THEN BEGIN
           kx              = (kx[indPos]-REVERSE(kx[indNeg])) / divFactor
           kx  = kx[(adj EQ 1):(-1-chopSiste)]
        ENDIF

        IF N_ELEMENTS(ky) GT 0 THEN BEGIN
           ky              = (ky[indPos]-REVERSE(ky[indNeg])) / divFactor 
           ky  = ky[(adj EQ 1):(-1-chopSiste)]
        ENDIF

        IF N_ELEMENTS(kz) GT 0 THEN BEGIN
           kz              = (kz[indPos]-REVERSE(kz[indNeg])) / divFactor
           kz  = kz[(adj EQ 1):(-1-chopSiste)]
        ENDIF

        ;;This is now handled by just calculating the power spectrum in SINGLE_SPACECRAFT_K_MEASUREMENT_FAST
        ;; IF N_ELEMENTS(BSpec) GT 0 THEN BEGIN
        ;;    BSpec           = (BSpec[indPos,*]*CONJ(REVERSE(BSpec[indNeg,*],1))) / divFactor
        ;;    BSpec = BSpec[(adj EQ 1):(-1-chopSiste),*]
        ;; ENDIF

        ;; IF N_ELEMENTS(JSpec) GT 0 THEN BEGIN
        ;;    JSpec           = (JSpec[indPos,*]*CONJ(REVERSE(JSpec[indNeg,*],1))) / divFactor
        ;;    JSpec = JSpec[(adj EQ 1):(-1-chopSiste),*]
        ;; ENDIF

        ;; IF N_ELEMENTS(magCSpec) GT 0 THEN BEGIN
        ;;    magCSpec        = (magCSpec[indPos]*CONJ(REVERSE(magCSpec[indNeg]))) / divFactor
        ;;    magCSpec = magCSpec[(adj EQ 1):(-1-chopSiste)]
        ;; ENDIF

        ;; IF N_ELEMENTS(kP) GT 0 THEN BEGIN
        ;;    kP              = (kP[indPos]+REVERSE(kP[indNeg])) / divFactor
        ;;    kP  = kP[(adj EQ 1):(-1-chopSiste)]
        ;; ENDIF

        freq               = freq[indPos]
        freq   = freq[(adj EQ 1):(-1-chopSiste)]

        inds               = INDGEN(N_ELEMENTS(freq))

        kP                 = SQRT(kx*kx+ky*ky)
        kPAngle            = ATAN(ky,kx)*!RADEG

     END
     ELSE: BEGIN

        inds       = WHERE(ABS(freq) LE 10.0)

     END
  ENDCASE

  IF KEYWORD_SET(freqLims) THEN BEGIN
     inds = WHERE(freq GE freqLims[0] AND freq LE freqLims[1],nInds)
     IF nInds LE 1 THEN STOP
  END



END

PRO PLOT_SINGLE_SPACECRAFT_K_MEASUREMENT, $
   TArr,freq, $
   Bx,By,Bz, $
   Jx,Jy,Jz, $
   kx,ky,kz, $
   kP,kPAngle, $
   ;; inds, $
   usedInds, $
   example, $
   BSPEC=BSpec, $
   JSPEC=JSpec, $
   MAGCSPEC=magCSpec, $
   POWFREQ=powFreq, $
   MAGERR=magErr, $
   ERRANGLE=errAngle, $
   PHASE_ERR=phaseErr, $
   EXAMPLE_MODE=example_mode, $
   PLOTDIR=plotDir, $
   SUFF=suff, $
   PARSE_B_AND_J_SAVEFILE=parse_B_and_J_saveFile, $
   SAVE_PS=save_ps, $
   TO_PDF=to_pdf, $
   PDF_TRANSPARENCY_LEVEL=pdf_transparency, $
   REMOVE_EPS=remove_eps, $
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
   PLOT_POSFREQ=plot_posFreq, $
   FOLD_NEGFREQ_ONTO_POS=fold_negFreq, $
   PLOT_KPERP_MAGNITUDE_FOR_KZ=plot_kperp_magnitude_for_kz, $
   PLOT_KX_VS_KY_FOR_KZ=plot_kx_vs_ky_for_kz, $
   PLOT_SMOOTHED_K_COMPONENTS=plot_smoothed_ks, $
   PLOT_ABS_SMOOTHED_K_COMPONENTS=plot_abs_smoothed_ks, $
   KX_VS_KY__PLOT_SMOOTHED=kx_vs_ky__plot_smoothed, $
   KP_ANGLE__PLOT_SMOOTHED=kP_angle__plot_smoothed, $
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
   PUBLICATION_SETTINGS=pubSettings, $
   PRE_VIII_LAYOUT=PRE_VIII_layout, $
   FOOTBALL_LAYOUT=football_layout, $
   FOOTBALL_YLOG=football_yLog, $
   FOOTBALL_COL2TITLE=football_col2Title, $
   FOOTBALL_KMAG=football_kMag, $
   FOOTBALL_NO_ERRANGLE=football_no_errAngle, $
   FOOTBALL_NO_MAGERR=football_no_magErr

  COMPILE_OPT IDL2,STRICTARRSUBS

  IF KEYWORD_SET(pubSettings) THEN BEGIN
     cs = KEYWORD_SET(football_layout) ? 0.7 : 1.1
     symSize = 0.4
  ENDIF

  ;; pSym         = 2            ;asterisk
  ;; pSym         = 1            ;plus sign
  pSym            = 7            ;x

  pSymEA          = 1            ;For the times when ISA(mark_ks_below_errAngle_thresh)
  pSymPE          = 3            ;For the times when ISA(mark_ks_below_phaseErr_thresh)
  pSymMEEA        = 2            ;For the times when ISA(mark_ks_below_both)
  
  dashSym         = [[-.5,0],[.5,0]]
  perpSym         = [[-.5,0.],[0.,0.],[0.,0.75],[0.,0.],[0.5,0.]]
  crossSym        = [[-.25,0.],[0.,0.],[0.,0.75],[0.,0.],[0.25,0.],[0.,0.],[0.,0.-.75]]

  USERSYM,perpSym               ;Just a dash
  pSymMagErr      = 8

  ;;Only for special people (see Beck)
  specialColor    = 70
  specialThick    = 3.0
  specialLineSty  = 5
  kPA_specialSty  = 5

  yStyler         = 16

  freqTitle       = 'Frequency (Hz)'
  freqTitle       = 'f!Dsp!N (Hz)'

  ;; perpSym         = STRING(120B)
  perpAll         = '!9' + STRING(120B) + '!X'

  GET_FREQ_INDS,freq,kx,ky,kz, $
                kP,kPAngle, $
                inds, $
                BSPEC=BSpec, $
                JSPEC=JSpec, $
                MAGCSPEC=magCSpec, $
                POWFREQ=powFreq, $
                MAGERR=magErr, $
                ERRANGLE=errAngle, $
                OUT_FITINDS=fitInds, $
                MAXFREQ=maxFreq, $
                PLOT_POSFREQ=plot_posFreq, $
                FOLD_NEGFREQ_ONTO_POS=fold_negFreq, $
                FREQLIMS=freqLims

  IF KEYWORD_SET(mark_ks_below_both) THEN BEGIN

     CASE N_ELEMENTS(mark_ks_below_both) OF

        1: BEGIN

           IF N_ELEMENTS(mark_ks_below_magErr_thresh) EQ 0 THEN BEGIN

              mark_ks_below_magErr_thresh = 1

           ENDIF

           IF N_ELEMENTS(mark_ks_below_errAngle_thresh) EQ 0 THEN BEGIN

              mark_ks_below_errAngle_thresh = 10

           ENDIF

           ;; IF N_ELEMENTS(mark_ks_below_phaseErr_thresh) EQ 0 THEN BEGIN

           ;;    mark_ks_below_phaseErr_thresh = 15

           ;; ENDIF

        END
        2: BEGIN

           IF N_ELEMENTS(mark_ks_below_magErr_thresh) EQ 0 THEN BEGIN

              mark_ks_below_magErr_thresh = mark_ks_below_both[0]

           ENDIF

           IF N_ELEMENTS(mark_ks_below_errAngle_thresh) EQ 0 THEN BEGIN

              mark_ks_below_errAngle_thresh = mark_ks_below_both[1]

           ENDIF

           ;; IF N_ELEMENTS(mark_ks_below_phaseErr_thresh) EQ 0 THEN BEGIN

           ;;    mark_ks_below_phaseErr_thresh = 15

           ;; ENDIF

        END
        3: BEGIN

           IF N_ELEMENTS(mark_ks_below_magErr_thresh) EQ 0 THEN BEGIN

              mark_ks_below_magErr_thresh = mark_ks_below_both[0]

           ENDIF

           IF N_ELEMENTS(mark_ks_below_errAngle_thresh) EQ 0 THEN BEGIN

              mark_ks_below_errAngle_thresh = mark_ks_below_both[1]

           ENDIF

           IF N_ELEMENTS(mark_ks_below_phaseErr_thresh) EQ 0 THEN BEGIN

              mark_ks_below_phaseErr_thresh = mark_ks_below_both[2]

           ENDIF

        END
     ENDCASE

  ENDIF

  ;;Want me to show you where it's totally rad?
  nCoolME_k = 0
  IF KEYWORD_SET(mark_ks_below_magErr_thresh) THEN BEGIN

     ;;i.e., if you are an int type, pay the boss
     IF (WHERE(SIZE(mark_ks_below_magErr_thresh,/TYPE) EQ [1,2,3,12,13,14,15]))[0] NE -1 THEN BEGIN

        CASE 1 OF
           mark_ks_below_magErr_thresh EQ 1: BEGIN
              mark_ks_below_magErr_thresh = 0.2
           END
           mark_ks_below_magErr_thresh GT 1: BEGIN
              mark_ks_below_magErr_thresh /= 100.
           END
        ENDCASE

     ENDIF

     coolME_k_ii = WHERE(magErr[inds] LE mark_ks_below_magErr_thresh,nCoolME_k)
     coolME_k_i  = inds[coolME_k_ii]
  ENDIF

  nCoolEA_k = 0
  IF KEYWORD_SET(mark_ks_below_errAngle_thresh) THEN BEGIN

        CASE 1 OF
           mark_ks_below_errAngle_thresh EQ 1: BEGIN
              mark_ks_below_errAngle_thresh = 10.D
           END
           (mark_ks_below_errAngle_thresh LE !PI): BEGIN
              PRINT,"Assuming your 'mark_ks_below_errAngle_thresh' is in radians ..."
              WAIT,2
              mark_ks_below_errAngle_thresh *= 180.D/!PI
           END
           ELSE:
        ENDCASE

     coolEA_k_ii = WHERE((errAngle[inds]*180.D/!PI) LE mark_ks_below_errAngle_thresh,nCoolEA_k)
     coolEA_k_i  = inds[coolEA_k_ii]
  ENDIF

  nCoolPE_k = 0
  IF KEYWORD_SET(mark_ks_below_phaseErr_thresh) THEN BEGIN

        CASE 1 OF
           mark_ks_below_phaseErr_thresh EQ 1: BEGIN
              mark_ks_below_phaseErr_thresh = 15.D
           END
           (mark_ks_below_phaseErr_thresh LE !PI): BEGIN
              PRINT,"Assuming your 'mark_ks_below_phaseErr_thresh' is in radians ..."
              WAIT,2
              mark_ks_below_phaseErr_thresh *= 180.D/!PI
           END
           ELSE:
        ENDCASE

        phase_zInd = 2
     coolPE_k_ii = WHERE((phaseErr[phase_zInd,inds]*180.D/!PI) LE mark_ks_below_phaseErr_thresh,nCoolPE_k)
     coolPE_k_i  = inds[coolPE_k_ii]
  ENDIF

  nCoolMEEA_k = 0
  IF KEYWORD_SET(mark_ks_below_both) THEN BEGIN

     nCoolME_k = 0
     nCoolEA_k = 0
     coolMEEA_k_i = CGSETINTERSECTION(coolME_k_i,coolEA_k_i,COUNT=nCoolMEEA_k,NORESULT=-1)

  ENDIF
  
  IF KEYWORD_SET(kP__angleRange) THEN BEGIN

     rotate_kPA = 0

     kPAngle    = NORMALIZE_ANGLE(kPAngle,MIN(kP__angleRange),MAX(kP__angleRange), $
                                  /DEGREE)

  ENDIF ELSE BEGIN
     histo   = HISTOGRAM(kPAngle,BINSIZE=90,MIN=-180,MAX=180,LOCATIONS=locs)
     rotate_kPA = (histo[0]+histo[3]) GT (histo[1]+histo[2])
     IF rotate_kPA THEN BEGIN
        kPAngle = (kPAngle + 360) MOD 360
     ENDIF
  ENDELSE


  ;;setup graphic layout
  columns   = KEYWORD_SET(football_layout) ? 2 : 3
  rows      = KEYWORD_SET(football_layout) ? 3 : 3 + KEYWORD_SET(example_mode)
  !P.MULTI  = [0, columns, rows, 0, 0]
  ;; !P.MULTI  = [0, columns, rows, 1, 0]

  ;;select output mode,
  ;;set output_mode =0 for screen, 1 for printer, 2 for postscript file
  CASE 1 OF
     KEYWORD_SET(save_ps): BEGIN
        output_mode = 2

        SET_PLOT_DIR,plotDir,/FOR_SINGLE_SC_WVEC,/ADD_TODAY
        PRINT,STRING(FORMAT='("plotDir  : ",A0)',plotDir)
     END
     ELSE: BEGIN
        output_mode = 0
     END
  ENDCASE

  ;; defFont = '!8'
  defFont = ''
  font    = defFont

  ;;setup file management, filenames for selected output mode
  IF KEYWORD_SET(save_ps) THEN BEGIN
     OUTPUT_SETUP,output_mode,plotDir,suff
  ENDIF
  
  ;;plot components of calculated k vectors
  IF example EQ 1 THEN BEGIN
     ysize = 10
  ENDIF
  IF example EQ 2 or example EQ 3 THEN BEGIN
     ysize = 2
  ENDIF
  kx_ysize = ysize
  ky_ysize = ysize
  kz_ysize = ysize

  IF KEYWORD_SET(parse_B_and_J_saveFile) THEN BEGIN
     kx_ysize = MAX(ABS(kx))
     ky_ysize = MAX(ABS(ky))
     kz_ysize = MAX(ABS(kz))
  ENDIF

  gimmeSmooths      = KEYWORD_SET(plot_smoothed_ks       ) OR KEYWORD_SET(plot_abs_smoothed_ks   ) OR $
                      KEYWORD_SET(kx_vs_ky__plot_smoothed) OR KEYWORD_SET(kP_angle__plot_smoothed)

  IF gimmeSmooths THEN BEGIN

     CASE SIZE(plot_smoothed_ks,/TYPE) OF
        8: BEGIN

           smooth_freq = plot_smoothed_ks.freq
           smooth_kx   = plot_smoothed_ks.kx
           smooth_ky   = plot_smoothed_ks.ky

           GET_FREQ_INDS,smooth_freq,smooth_kx,smooth_ky,!NULL, $
                         smooth_kP,smooth_kPAngle, $
                         smooth_inds, $
                         MAXFREQ=maxFreq, $
                         PLOT_POSFREQ=plot_posFreq, $
                         FOLD_NEGFREQ_ONTO_POS=fold_negFreq, $
                         FREQLIMS=freqLims

        END
        ELSE: BEGIN

           smooth_inds    = inds
           smooth_freq    = freq

           smooth_kx      = SMOOTH((KEYWORD_SET(plot_abs_smoothed_ks) ? ABS(kx) : kx),smInd, $
                                   EDGE_TRUNCATE=edge_truncate,EDGE_MIRROR=edge_mirror,EDGE_WRAP=edge_wrap)

           smooth_ky      = SMOOTH((KEYWORD_SET(plot_abs_smoothed_ks) ? ABS(ky) : ky),smInd, $
                                   EDGE_TRUNCATE=edge_truncate,EDGE_MIRROR=edge_mirror,EDGE_WRAP=edge_wrap)

        END
     ENDCASE


     smooth_kP         = SQRT(smooth_kx*smooth_kx+smooth_ky*smooth_ky)
     smooth_kPAngle    = ATAN(smooth_ky,smooth_kx)*!RADEG

     IF KEYWORD_SET(kP__angleRange) THEN BEGIN

        smooth_kPAngle = NORMALIZE_ANGLE(smooth_kPAngle,MIN(kP__angleRange),MAX(kP__angleRange), $
                                         /DEGREE)

        ;; those = WHERE(smooth_kPAngle LT MIN(kP__angleRange),nThose)
        ;; IF nThose GT 0 THEN BEGIN
        ;;    smooth_kPAngle[those] = (smooth_kPAngle[those] + MAX(kP__angleRange)) MOD MAX(kP__angleRange)
        ;; ENDIF

     ENDIF ELSE BEGIN
        IF rotate_kPA THEN BEGIN
           smooth_kPAngle = (smooth_kPAngle + 360) MOD 360
        ENDIF
     ENDELSE


  ENDIF

  kx_yRange = [-kx_ysize,kx_ysize]
  ky_yRange = [-ky_ysize,ky_ysize]
  k_yRange  = [MIN([kx_yRange,ky_yRange]),MAX([kx_yRange,ky_yRange])]

  muLetter = '!4' + String('154'O) + '!X'
  ;;THe following case thing handles plotting of Bx,By,Bz, Jx,Jy,Jz, and 
  CASE 1 OF
     KEYWORD_SET(football_layout): BEGIN

        bxCol = 240
        byCol = 50
        
        jMagSpecCol = 30

        thetaErrCol  = byCol

        ;;THe idea here is to kill the x axis, since it's common to Bx, By, and Jz

        ;;Row 0: plots 0,1,2
        ;;Row 1: plots 3,4,5
        ;;Row 2: plots 6,7,8

        !P.MULTI       = [0,1,1,0,0]

        colTitleSpace  = 0.09
        noTitleSpace   = 0.02
        titleSpace     = 0.08
        totForSpace    = colTitleSpace + 2*noTitleSpace + titleSpace

        rows           = 3
        panYSpace      = (1.-totForSpace)/rows

        lEdgeSpace     = 0.1
        rEdgeSpace     = 0.08
        
        lEdgeL         = lEdgeSpace
        colWid         = 0.5 - lEdgeSpace - rEdgeSpace + rEdgeSpace*0.5
        rEdgeL         = lEdgeL+colWid
        lColCtr        = MEAN([lEdgeL,rEdgeL])

        lEdgeR         = 0.5+lEdgeSpace - rEdgeSpace/2.
        rEdgeR         = lEdgeR+colWid
        rColCtr        = MEAN([lEdgeR,rEdgeR])

        ;;First column
        Bx_upper       = titleSpace+panYSpace+noTitleSpace+panYSpace+noTitleSpace+panYSpace
        Bx_lower       = titleSpace+panYSpace+noTitleSpace+panYSpace+noTitleSpace
        
        By_upper       = titleSpace+panYSpace+noTitleSpace+panYSpace
        By_lower       = titleSpace+panYSpace+noTitleSpace
        
        posBx          = [lEdgeL,Bx_lower,lEdgeL+colWid,Bx_upper]
        posBy          = [lEdgeL,By_lower,lEdgeL+colWid,By_upper]
        posJz          = [lEdgeL,titleSpace,lEdgeL+colWid,titleSpace+panYSpace]

        ;;Second column
        By_upper       = titleSpace+panYSpace+noTitleSpace+panYSpace
        By_lower       = titleSpace+panYSpace+noTitleSpace
        
        Bx_upper       = titleSpace+panYSpace+noTitleSpace+panYSpace+noTitleSpace+panYSpace
        Bx_lower       = titleSpace+panYSpace+noTitleSpace+panYSpace+noTitleSpace
        
        posSpec        = [lEdgeR,Bx_lower,lEdgeR+colWid,Bx_upper]
        poskxy         = [lEdgeR,By_lower,lEdgeR+colWid,By_upper]
        poskPA         = [lEdgeR,titleSpace,lEdgeR+colWid,titleSpace+panYSpace]

        
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        ;;Plot first column

        tPlotMe = TArr[usedInds]-TArr[usedInds[0]]

        timeTitle = font + 't since' + TIME_TO_STR(TArr[usedInds[0]]) + '(s)'
        XYOUTS,lColCtr,titleSpace/4,timeTitle, $
               /NORMAL, $
               ALIGNMENT=0.5, $
               CHARSIZE=cs

        XYOUTS,MEAN([lEdgeR,lEdgeR+colWid]),titleSpace/4,freqTitle, $
               /NORMAL, $
               ALIGNMENT=0.5, $
               CHARSIZE=cs
        
        colTitle = 'Inputs'
        XYOUTS,lColCtr,1.-colTitleSpace/2.,colTitle, $
               /NORMAL, $
               ALIGNMENT=0.5, $
               CHARSIZE=cs*1.5

        IF KEYWORD_SET(football_col2Title) THEN BEGIN
           XYOUTS,rColCtr,1.-colTitleSpace/2.,football_col2Title, $
                  /NORMAL, $
                  ALIGNMENT=0.5, $
                  CHARSIZE=cs*1.5
        ENDIF
        ;; freqBlankTickName = MAKE_ARRAY(N_ELEMENTS(freqTickV),/STRING,VALUE=' ')

        ;; !P.MULTI[0]  = 4
        PLOT,tPlotMe,Bx[usedInds], $
             ;; XTITLE='t', $
             YTITLE=font + 'B!Dx!N (nT)', $
             XSTYLE=1, $
             YSTYLE=yStyler, $
             ;; XTICKV=freqTickV, $
             ;; XTICKNAME=freqBlankTickName, $
             CHARSIZE=cs, $
             XTICKFORMAT='(A1)', $
             SYMSIZE=symSize, $
             POSITION=posBx, $
             /NOERASE

        ;; !P.MULTI[0] = 6
        PLOT,tPlotMe,By[usedInds], $
             ;; XTITLE='t since' + TIME_TO_STR(TArr[usedInds[0]]) + '(s)', $
             YTITLE=font + 'B!Dy!N (nT)', $
             XSTYLE=1, $
             YSTYLE=yStyler, $
             ;; XTICKV=freqTickV, $
             ;; XTICKNAME=freqBlankTickName, $
             XTICKFORMAT='(A1)', $
             CHARSIZE=cs, $
             SYMSIZE=symSize, $
             POSITION=posBy, $
             /NOERASE

        ;; !P.MULTI[0] = 2
        
        PLOT,tPlotMe,Jz[usedInds], $
             ;; XTITLE='t since' + TIME_TO_STR(TArr[usedInds[0]]) + '(s)', $
             ;; YTITLE='!4l!3!D0 !N' + font + 'J!Dz!N', $
             YTITLE=font + 'J!Dz!N (!4l!3A m!U-2!N)', $
             XSTYLE=1, $
             CHARSIZE=cs, $
             SYMSIZE=symSize, $
             ;; XTICKFORMAT='(A1)', $
             ;; XTICK_GET=freqTickV, $
             POSITION=posJz, $
             /NOERASE

        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        ;;Plot spectra in second column

        ;; jYRange = [0,MAX(magCSpec)]
        ;; yLog = 1
        CASE 1 OF
           KEYWORD_SET(football_yLog): BEGIN

              powPos  = WHERE((magCSpec GT 0) AND (BSpec[*,0] GT 0) AND (BSpec[*,1] GT 0))
              jYRange = [MIN(JSpec[powPos,2]) > (MIN(magCSpec) < MIN(BSpec[powPos,0]) < MIN(BSpec[powPos,1])), $
                         (1 < MAX(JSpec[powPos,2])) < (MAX(magCSpec) > MAX(BSpec[powPos,0]) > MAX(BSpec[powPos,1]))]

              currency = (MAX(ALOG10(jYRange))-MIN(ALOG10(jYRange)))+MIN(ALOG10(jYRange))

           END
           ELSE: BEGIN

              powPos  = WHERE(powFreq GT 0)
              jYRange = [0,MAX(JSpec[powPos,2]) < (MAX(magCSpec) > MAX(BSpec[powPos,0]) > MAX(BSpec[powPos,1]))]

              currency = (MAX(jYRange)-MIN(jYRange))+MIN(jYRange)

           END
        ENDCASE

        jSpecVar1 = JSpec[*,2]
        jSpecVar2 = magCSpec

        PLOT,powFreq,jSpecVar1, $
             ;; XTITLE='t', $
             YTITLE='Power (fractional)' , $
             XRANGE=page1__freqRange, $
             YRANGE=jYRange, $
             XSTYLE=1, $
             YSTYLE=yStyler, $
             YLOG=football_yLog, $
             XTICKLEN=1.0, $
             YTICKLEN=1.0, $
             XGRIDSTYLE=1, $
             YGRIDSTYLE=1, $
             XMINOR=5, $
             ;; XTICKV=freqTickV, $
             ;; XTICKNAME=freqBlankTickName, $
             CHARSIZE=cs, $
             XTICKFORMAT='(A1)', $
             YTICKFORMAT=KEYWORD_SET(football_yLog) ? 'exponentlabel' : !NULL, $
             SYMSIZE=symSize, $
             POSITION=posSpec, $
             /NOERASE, $
             /NODATA

        OPLOT,powFreq,BSpec[*,0], $
              COLOR=bxCol

        OPLOT,powFreq,BSpec[*,1], $
              COLOR=byCol

        OPLOT,powFreq,jSpecVar1

        ;;Plot mag current spec?
        ;; jLineStyle = 2 ;dotted line
        ;; jLineStyle = 1 ;dashed line??
        ;; ;; OPLOT,freq[inds],JSpecNorm[inds], $
        ;; OPLOT,powFreq,jSpecVar2, $
        ;;       LINESTYLE=jLineStyle, $
        ;;       COLOR=jMagSpecCol;; , $
              ;; MAX_VALUE=MAX(magCSpec)

        legXSymPos1 = 0.73*((MAX(powFreq)-MIN(powFreq))+MIN(powFreq))
        legXSymPos2 = 0.78*((MAX(powFreq)-MIN(powFreq))+MIN(powFreq))
        legXPos1    = 0.8*((MAX(powFreq)-MIN(powFreq))+MIN(powFreq))
        legXPos2    = 0.85*((MAX(powFreq)-MIN(powFreq))+MIN(powFreq))

        IF KEYWORD_SET(football_yLog) THEN BEGIN

           spaceFrac   = 0.6
           downLog     = 1.2
           distVec     = REVERSE(INDGEN(3)*spaceFrac - spaceFrac - downLog)

           legYSymPos  = 0.6*( distVec + currency )
           legYPos     = 0.6*( distVec + currency ) - (currency + 0.1)

           legYSymPos = 10.^(legYSymPos)
           legYPos    = 10.^(legYPos)

        ENDIF ELSE BEGIN

           spaceFrac   = 0.13
           distVec     = REVERSE(INDGEN(3)*spaceFrac - spaceFrac + 1)

           legYSymPos  = 0.93*( distVec * currency )
           legYPos     = 0.93*( distVec * currency ) - currency * 0.02

        ENDELSE

        XYOUTS,legXPos1,legYPos[0],'B!Dx!N',CHARSIZE=cs
        PLOTS,legXSymPos1,legYSymPos[0],COLOR=bxCol
        PLOTS,legXSymPos2,legYSymPos[0],COLOR=bxCol,/CONTINUE

        XYOUTS,legXPos1,legYPos[1],'B!Dy!N',CHARSIZE=cs
        PLOTS,legXSymPos1,legYSymPos[1],COLOR=byCol
        PLOTS,legXSymPos2,legYSymPos[1],COLOR=byCol,/CONTINUE

        XYOUTS,legXPos1,legYPos[2],'J!Dz!N',CHARSIZE=cs
        PLOTS,legXSymPos1,legYSymPos[2]
        PLOTS,legXSymPos2,legYSymPos[2],/CONTINUE

     END
     ELSE: BEGIN

        tPlotMe = TArr[usedInds]-TArr[usedInds[0]]
        
        PLOT,tPlotMe,Bx[usedInds], $
             XTITLE='t', $
             YTITLE=font + 'B!Dx!N' + (KEYWORD_SET(PRE_VIII_layout) ? ' (nT)' : ''), $
             XSTYLE=1, $
             YSTYLE=yStyler, $
             CHARSIZE=cs, $
             SYMSIZE=symSize

        PLOT,tPlotMe,By[usedInds], $
             XTITLE='t since' + TIME_TO_STR(TArr[usedInds[0]]) + '(s)', $
             YTITLE=font + 'B!Dy!N (nT)', $
             XSTYLE=1, $
             YSTYLE=yStyler, $
             CHARSIZE=cs, $
             SYMSIZE=symSize

        IF ~KEYWORD_SET(PRE_VIII_layout) THEN BEGIN

           PLOT,tPlotMe,Bz[usedInds], $
                XTITLE='t', $
                YTITLE=font + 'B!Dz!N', $
                XSTYLE=1, $
                YSTYLE=yStyler, $
                CHARSIZE=cs, $
                SYMSIZE=symSize

           PLOT,tPlotMe,Jx[usedInds], $
                XTITLE='t', $
                YTITLE='!4l!3!D0 !N' + font + 'J!Dx!N', $
                XSTYLE=1, $
                ;; YSTYLE=1, $
                CHARSIZE=cs, $
                SYMSIZE=symSize

           PLOT,tPlotMe,Jy[usedInds], $
                XTITLE='t since' + TIME_TO_STR(TArr[usedInds[0]]) + '(s)', $
                YTITLE='!4l!3!D0 !N' + font + 'J!Dy!N (!4l!N' + font + 'T/m)', $
                XSTYLE=1, $
                CHARSIZE=cs, $
                SYMSIZE=symSize

        ENDIF

        PLOT,tPlotMe,Jz[usedInds], $
             XTITLE='t', $
             ;; YTITLE='!4l!3!D0 !N' + font + 'J!Dz!N', $
             YTITLE=font + 'J!Dz!N (!4l!3A m!U-2!N)', $
             XSTYLE=1, $
             CHARSIZE=cs, $
             SYMSIZE=symSize

     END
  ENDCASE


  IF KEYWORD_SET(football_layout) THEN BEGIN

     PLOT,freq[inds],kx[inds], $
          ;; YTITLE='k!Dx!N (m!U-1!N)', $
          ;; COLOR=KEYWORD_SET(football_layout) ? bxCol : !NULL, $
          YTITLE=(KEYWORD_SET(football_layout) ? 'Wave number' : 'k!Dx!N') + ' (km!U-1!N)', $
          XTITLE=(KEYWORD_SET(PRE_VIII_layout) AND ~KEYWORD_SET(football_layout) ? freqTitle : ''), $
          XRANGE=page1__freqRange, $
          YRANGE=k_yRange, $
          XSTYLE=1, $
          ;; YSTYLE=3, $
          XTICKLEN=1.0, $
          YTICKLEN=1.0, $
          XGRIDSTYLE=1, $
          YGRIDSTYLE=1, $
          XMINOR=5, $
          CHARSIZE=cs, $
          XTICKFORMAT='(A1)', $
          SYMSIZE=symSize, $
          POSITION=KEYWORD_SET(football_layout) ? poskxy : !NULL, $
          NOERASE=KEYWORD_SET(football_layout), $
          /NODATA
     ;; IF KEYWORD_SET(diag) THEN BEGIN & diagInd++ & PRINT,diagInd,'  ',!P.MULTI & ENDIF

     IF ~KEYWORD_SET(only_ks_below_magErr_thresh) THEN BEGIN
        OPLOT,freq[inds],kx[inds], $
              ;; LINESTYLE=jLineStyle, $
              COLOR=bxCol
     ENDIF

        ;;Because of MARK_KS_BELOW_MAGERR_THRESH
        IF nCoolME_k GT 0 THEN BEGIN

           OPLOT,freq[coolME_k_i],kx[coolME_k_i], $
                 COLOR=bxCol, $
                 PSYM=pSym, $
                 SYMSIZE=symSize

        ENDIF

        IF nCoolEA_k GT 0 THEN BEGIN

           OPLOT,freq[coolEA_k_i],kx[coolEA_k_i], $
                 COLOR=bxCol, $
                 PSYM=pSymEA, $
                 SYMSIZE=symSize

        ENDIF

        IF nCoolPE_k GT 0 THEN BEGIN

           OPLOT,freq[coolPE_k_i],kx[coolPE_k_i], $
                 COLOR=bxCol, $
                 PSYM=pSymPE, $
                 SYMSIZE=symSize

        ENDIF

        IF nCoolMEEA_k GT 0 THEN BEGIN

           PLOTS,freq[coolMEEA_k_i],kx[coolMEEA_k_i], $
                 COLOR=bxCol, $
                 PSYM=pSymMEEA, $
                 SYMSIZE=symSize, $
                 NOCLIP=0

        ENDIF

  ENDIF ELSE BEGIN

     PLOT,freq[inds],kx[inds], $
          ;; YTITLE='k!Dx!N (m!U-1!N)', $
          COLOR=KEYWORD_SET(football_layout) ? bxCol : !NULL, $
          YTITLE=(KEYWORD_SET(football_layout) ? 'Wave number' : 'k!Dx!N') + ' (km!U-1!N)', $
          XTITLE=(KEYWORD_SET(PRE_VIII_layout) AND ~KEYWORD_SET(football_layout) ? freqTitle : ''), $
          XRANGE=page1__freqRange, $
          ;; YRANGE=kx_yRange, $
          YRANGE=kx_yRange, $
          XSTYLE=1, $
          YSTYLE=1, $
          XTICKLEN=1.0, $
          YTICKLEN=1.0, $
          XGRIDSTYLE=1, $
          YGRIDSTYLE=1, $
          XMINOR=5, $
          CHARSIZE=cs, $
          SYMSIZE=symSize, $
          POSITION=KEYWORD_SET(football_layout) ? poskxy : !NULL, $
          NOERASE=KEYWORD_SET(football_layout)
     ;; IF KEYWORD_SET(diag) THEN BEGIN & diagInd++ & PRINT,diagInd,'  ',!P.MULTI & ENDIF

  ENDELSE
  
  IF example EQ 2 THEN OPLOT,kx+0.1,LINESTYLE=2 ;dashed line

  IF KEYWORD_SET(plot_smoothed_ks) OR KEYWORD_SET(plot_abs_smoothed_ks) AND ~KEYWORD_SET(football_layout)  THEN BEGIN

     OPLOT,smooth_freq[smooth_inds], $
           smooth_kx[smooth_inds], $
           COLOR=250

  ENDIF

  IF KEYWORD_SET(overplot_doubly_smoothed) THEN BEGIN
     OPLOT,freq[inds],SMOOTH((KEYWORD_SET(plot_abs_smoothed_ks) ? ABS(kx[inds]) : kx[inds]),dbSmInd, $
                             EDGE_TRUNCATE=edge_truncate,EDGE_MIRROR=edge_mirror,EDGE_WRAP=edge_wrap),COLOR=90
  ENDIF


  IF KEYWORD_SET(kx_specialFreqs) THEN BEGIN

     PLOT_SPECIAL_FREQS,freq,kx_specialFreqs,kx_ysize,kx_totSpecial_i

  ENDIF

  CASE 1 OF
     KEYWORD_SET(football_layout): BEGIN

        IF ~KEYWORD_SET(only_ks_below_magErr_thresh) THEN BEGIN

           jLineStyle = 1       ;dotted line
           OPLOT,freq[inds],ky[inds], $
                 ;; LINESTYLE=jLineStyle, $
                 COLOR=byCol
        ENDIF

        ;;Because of MARK_KS_BELOW_MAGERR_THRESH
        IF nCoolME_k GT 0 THEN BEGIN

           OPLOT,freq[coolME_k_i],ky[coolME_k_i], $
                 COLOR=byCol, $
                 PSYM=pSym, $
                 SYMSIZE=symSize

        ENDIF

        IF nCoolEA_k GT 0 THEN BEGIN

           OPLOT,freq[coolEA_k_i],ky[coolEA_k_i], $
                 COLOR=byCol, $
                 PSYM=pSymEA, $
                 SYMSIZE=symSize

        ENDIF

        IF nCoolPE_k GT 0 THEN BEGIN

           OPLOT,freq[coolPE_k_i],ky[coolPE_k_i], $
                 COLOR=byCol, $
                 PSYM=pSymPE, $
                 SYMSIZE=symSize

        ENDIF

        IF nCoolMEEA_k GT 0 THEN BEGIN

           PLOTS,freq[coolMEEA_k_i],ky[coolMEEA_k_i], $
                 COLOR=byCol, $
                 PSYM=pSymMEEA, $
                 SYMSIZE=symSize, $
                 NOCLIP=0

        ENDIF

        spaceFrac   = 0.1
        distVec     = REVERSE(INDGEN(3)*spaceFrac - spaceFrac + 1)
        rango       = (MAX(k_yRange)-MIN(k_yRange))

        legYSymPos  = 0.925*( distVec * rango )+MIN(k_yRange)
        legYPos     = 0.920*( distVec * rango )+MIN(k_yRange)


        XYOUTS,legXPos1,legYPos[0],'k!Dx!N',CHARSIZE=cs
        PLOTS,legXSymPos1,legYSymPos[0],COLOR=bxCol
        PLOTS,legXSymPos2,legYSymPos[0],COLOR=bxCol,/CONTINUE

        XYOUTS,legXPos1,legYPos[1],'k!Dy!N',CHARSIZE=cs
        PLOTS,legXSymPos1,legYSymPos[1],COLOR=byCol
        PLOTS,legXSymPos2,legYSymPos[1],COLOR=byCol,/CONTINUE

        IF KEYWORD_SET(football_kMag) THEN BEGIN

           XYOUTS,legXPos1,legYPos[2],'k!D' + perpAll + '!N',CHARSIZE=cs 
           PLOTS,legXSymPos1,legYSymPos[2];,COLOR=byCol
           PLOTS,legXSymPos2,legYSymPos[2],/CONTINUE;,COLOR=byCol,/CONTINUE

           OPLOT,freq[inds],(SQRT(kx*kx+ky*ky))[inds];, $
                 ;; PSYM=pSym, $
                 ;; SYMSIZE=symSize

           IF nCoolME_k GT 0 THEN BEGIN

              OPLOT,freq[coolME_k_i],(SQRT(kx*kx+ky*ky))[coolME_k_i], $
                    ;; COLOR=bxCol, $
                    PSYM=pSym, $
                    SYMSIZE=symSize

           ENDIF

           IF nCoolEA_k GT 0 THEN BEGIN

              OPLOT,freq[coolEA_k_i],(SQRT(kx*kx+ky*ky))[coolEA_k_i], $
                    ;; COLOR=bxCol, $
                    PSYM=pSymEA, $
                    SYMSIZE=symSize

           ENDIF

           IF nCoolPE_k GT 0 THEN BEGIN

              OPLOT,freq[coolPE_k_i],(SQRT(kx*kx+ky*ky))[coolPE_k_i], $
                    ;; COLOR=bxCol, $
                    PSYM=pSymPE, $
                    SYMSIZE=symSize

           ENDIF

           IF nCoolMEEA_k GT 0 THEN BEGIN

              PLOTS,freq[coolMEEA_k_i],(SQRT(kx*kx+ky*ky))[coolMEEA_k_i], $
                    ;; COLOR=bxCol, $
                    PSYM=pSymMEEA, $
                    SYMSIZE=symSize, $
                    NOCLIP=0

           ENDIF

        ENDIF

     END
     ELSE: BEGIN

        PLOT,freq[inds],ky[inds], $
             ;; YTITLE='k!Dy!N (m!U-1!N)', $
             YTITLE='k!Dy!N (km!U-1!N)', $
             XTITLE=KEYWORD_SET(football_layout) ? !NULL : freqTitle, $
             XRANGE=page1__freqRange, $
             YRANGE=ky_yRange, $
             XSTYLE=1, $
             YSTYLE=1, $
             XTICKLEN=1.0, $
             YTICKLEN=1.0, $
             XGRIDSTYLE=1, $
             YGRIDSTYLE=1, $
             XMINOR=5, $
             CHARSIZE=cs, $
             SYMSIZE=symSize, $
             POSITION=KEYWORD_SET(football_layout) ? poskxy : !NULL, $
             NOERASE=KEYWORD_SET(football_layout)
        ;; IF KEYWORD_SET(diag) THEN BEGIN & diagInd++ & PRINT,diagInd,'  ',!P.MULTI & ENDIF
     END
  ENDCASE

  IF example EQ 2 THEN OPLOT,ky+0.1,LINESTYLE=2
  IF KEYWORD_SET(plot_smoothed_ks) OR KEYWORD_SET(plot_abs_smoothed_ks) AND ~KEYWORD_SET(football_layout) THEN BEGIN
     OPLOT,smooth_freq[smooth_inds], $
           smooth_ky[smooth_inds], $
           COLOR=250
  ENDIF

  IF KEYWORD_SET(overplot_doubly_smoothed) THEN BEGIN
     OPLOT,freq[inds], $
           SMOOTH((KEYWORD_SET(plot_abs_smoothed_ks) ? ABS(ky[inds]) : ky[inds]),dbSmInd, $
                  EDGE_TRUNCATE=edge_truncate,EDGE_MIRROR=edge_mirror,EDGE_WRAP=edge_wrap), $
           COLOR=90
  ENDIF

  IF KEYWORD_SET(ky_specialFreqs) THEN BEGIN

     PLOT_SPECIAL_FREQS,freq,ky_specialFreqs,ky_ysize,ky_totSpecial_i

     ;; ky_totSpecial_i = !NULL
     ;; FOR jD=0,N_ELEMENTS(ky_specialFreqs)/2-1 DO BEGIN

     ;;    tmpSpFreq = ky_specialFreqs[*,jD]

     ;;    ;;Special markers
     ;;    line1_freq = MAKE_ARRAY(2,VALUE=tmpSpFreq[0])
     ;;    line1_ky   = [-ky_ysize,ky_ysize]
     ;;    line2_freq = MAKE_ARRAY(2,VALUE=tmpSpFreq[1])
     ;;    line2_ky   = [-ky_ysize,ky_ysize]

     ;;    ;;Plot a special line, if there's any helping it
     ;;    special_i  = WHERE(freq GE tmpSpFreq[0] AND freq LE tmpSpFreq[1],nSpecial)
     ;;    ;; IF nSpecial GT 0 THEN BEGIN
     ;;    ;;    OPLOT,freq[special_i],ky[special_i],COLOR=specialColor
     ;;    ;; ENDIF

     ;;    ;;Now plot all those special people
     ;;    OPLOT,line1_freq,line1_ky,LINESTYLE=specialLineSty,COLOR=specialColor,THICK=specialThick
     ;;    OPLOT,line2_freq,line2_ky,LINESTYLE=specialLineSty,COLOR=specialColor,THICK=specialThick

     ;;    ky_totSpecial_i = [ky_totSpecial_i,special_i]
     ;; ENDFOR

  ENDIF

  IF KEYWORD_SET(make_kx_vs_ky_special) OR KEYWORD_SET(make_kPAngle_special) THEN BEGIN

     MAKING_KX_KY_SPECIAL,kx_specialBounds,ky_specialBounds, $
                          kx_totSpecial_i,ky_totSpecial_i, $
                          OUT_KXY_SPI=kxy_spI, $
                          OUT_KPA_SPI=kPA_spI

  ENDIF


  IF ~KEYWORD_SET(football_layout) THEN BEGIN
     
     IF KEYWORD_SET(plot_kx_vs_ky_for_kz) THEN BEGIN
        
        ;;Now kx vs ky

        bound = MAX(ABS([kx[inds],ky[inds]]))

        PLOT,kx[inds],ky[inds], $
             XTITLE='k!Dx!N  (km!U-1!N)', $
             YTITLE='k!Dy!N  (km!U-1!N)', $
             PSYM=pSym, $
             YRANGE=[-bound,bound], $
             XRANGE=[-bound,bound], $
             ;; XSTYLE=2, $
             ;; YSTYLE=2, $
             XTICKLEN=1.0, $
             YTICKLEN=1.0, $
             XGRIDSTYLE=1, $
             YGRIDSTYLE=1, $
             CHARSIZE=cs, $
             SYMSIZE=symSize

        IF KEYWORD_SET(kx_vs_ky__plot_smoothed) THEN BEGIN

           ;; CASE SIZE(kx_vs_ky__plot_smoothed) OF
           ;;    8: BEGIN
           ;;       OPLOT,kx_vs_ky__plot_smoothed.kx, $
           ;;             kx_vs_ky__plot_smoothed.ky, $
           ;;             COLOR=250, $
           ;;             PSYM=1
           ;;    END
           ;;    ELSE: BEGIN
           OPLOT,smooth_kx[smooth_inds],smooth_ky[smooth_inds], $
                 COLOR=250, $
                 PSYM=1
           ;;    END
           ;; ENDCASE

        ENDIF

        IF KEYWORD_SET(make_kx_vs_ky_special) THEN BEGIN

           OPLOT,kx[kxy_spI],ky[kxy_spI], $
                 PSYM=1, $
                 COLOR=specialColor, $
                 LINESTYLE=0

        ENDIF

     ENDIF ELSE BEGIN

        PLOT,freq[inds],kz[inds], $
             ;; YTITLE='k!Dz!N (m!U-1!N)', $
             YTITLE='k!Dz!N (km!U-1!N)', $
             XTITLE='', $
             XRANGE=page1__freqRange, $
             YRANGE=[-kz_ysize,kz_ysize], $
             XSTYLE=1, $
             YSTYLE=1, $
             XMINOR=5, $
             CHARSIZE=cs, $
             SYMSIZE=symSize
        IF example EQ 2 THEN OPLOT,kz+0.1,LINESTYLE=2
        IF KEYWORD_SET(plot_smoothed_ks) THEN BEGIN
           OPLOT,freq[inds],SMOOTH(kz[inds],smInd,EDGE_TRUNCATE=edge_truncate,EDGE_MIRROR=edge_mirror,EDGE_WRAP=edge_wrap),COLOR=250
        ENDIF

        IF KEYWORD_SET(overplot_doubly_smoothed) THEN BEGIN
           OPLOT,freq[inds],SMOOTH(kz[inds],dbSmInd,EDGE_TRUNCATE=edge_truncate,EDGE_MIRROR=edge_mirror,EDGE_WRAP=edge_wrap),COLOR=90
        ENDIF

     ENDELSE

  ENDIF
  

  ;; IF KEYWORD_SET(diag) THEN BEGIN & diagInd++ & PRINT,diagInd,'  ',!P.MULTI & ENDIF

  IF ~(KEYWORD_SET(PRE_VIII_layout) OR KEYWORD_SET(football_layout)) THEN BEGIN

     IF KEYWORD_SET(save_ps) THEN BEGIN
        CONCLUDE_OUTPUT,output_mode,plotDir,suff, $
                        TO_PDF=to_pdf, $
                        PDF_TRANSPARENCY_LEVEL=pdf_transparency, $
                        REMOVE_EPS=remove_eps

        ;;Don't conclude, just keep it movin'
        ;; ERASE
     ENDIF

  ENDIF

  ;;Show Hanning windowed?
  IF KEYWORD_SET(hanning) THEN BEGIN

     IF ~KEYWORD_SET(save_ps) THEN BEGIN
        WINDOW,3,XSIZE=700,YSIZE=800
     ENDIF

     ;; !P.MULTI  = [0, columns, rows, 0, 0]
     !P.MULTI[0]  = 0
     !P.MULTI[3]  = 1 ;;Next page?

     PLOT,freq[inds],kx[inds], $
          YTITLE='x component', $
          XTITLE='FFT argument', $
          CHARSIZE=cs, $
          SYMSIZE=symSize, $
          YRANGE=kx_yRange
     IF example EQ 2 THEN OPLOT,kx+0.1,LINESTYLE=2 ;dashed line

     PLOT,freq[inds],ky[inds], $
          YTITLE='y component', $
          XTITLE='FFT argument', $
          CHARSIZE=cs, $
          SYMSIZE=symSize, $
          YRANGE=ky_yRange
     IF example EQ 2 THEN OPLOT,ky+0.1,LINESTYLE=2

     PLOT,freq[inds],kz[inds], $
          YTITLE='z component', $
          XTITLE='FFT argument', $
          CHARSIZE=cs, $
          SYMSIZE=symSize, $
          YRANGE=[-kz_ysize,kz_ysize]
     IF example EQ 2 THEN OPLOT,kz+0.1,LINESTYLE=2



  ENDIF

  ;;Old or new style?
  IF KEYWORD_SET(oo_plots) THEN BEGIN
     window = WINDOW(DIMENSIONS=[700,800])

     bro = PLOT(freq[inds],kP[inds], $
                ;; XTITLE=freqTitle, $
                ;; YTITLE='|k!Dperp!N (m$^{-1}$)', $
                ;; YTITLE='|k!Dperp!N (km$^{-1}$)', $
                YTITLE='|k!D' + perpAll +  '!N (km$^{-1}$)', $
                CURRENT=window, $
                LAYOUT=[1,2,1])

     ;; bro.axes

     bro = PLOT(freq[inds],smooth_kP[inds], $
                ;; XTITLE=freqTitle, $
                ;; YTITLE='|k!Dperp!N (m$^{-1}$)', $
                ;; YTITLE='|k!Dperp!N (km$^{-1}$)', $
                YTITLE='|k!D' + perpAll +  '!N (km$^{-1}$)', $
                COLOR='RED', $
                /OVERPLOT, $
                CURRENT=window)

     ;;Kperp angle plot
     ;; kPAngle = ATAN(ky,kx)*!RADEG
     ;; kPAngle = ATAN(SMOOTH(ky,smInd),SMOOTH(kx,smInd))*!RADEG

     bro = PLOT(freq[inds],kPAngle[inds], $
                XTITLE=freqTitle, $
                ;; XTITLE='T since' + TIME_TO_STR(TArr[usedInds[0]]) + '(s)', $
                ;; YTITLE='|$\theta$(k!Dperp!N)', $
                YTITLE='|$\theta$(k!D' +perpAll + '!N)', $
                CURRENT=window, $
                LAYOUT=[1,2,2])

     ;; smkPAngle = SMOOTH(kPAngle[inds],smInd,EDGE_TRUNCATE=edge_truncate,EDGE_MIRROR=edge_mirror,EDGE_WRAP=edge_wrap)
     bro = PLOT(freq[inds],smooth_kPAngle[inds], $
                XTITLE=freqTitle, $
                ;; YTITLE='|$\theta$(k!Dperp!N)', $
                YTITLE='|$\theta$(k!D' +perpAll + '!N)', $
                COLOR='RED', $
                /OVERPLOT)
  ENDIF ELSE BEGIN


     CASE 1 OF
        KEYWORD_SET(football_layout): BEGIN

        END
        KEYWORD_SET(PRE_VIII_layout): BEGIN

        columns      = 1
        rows         = 3
        !P.MULTI[0]  = 1
        !P.MULTI[1]  = columns
        !P.MULTI[2]  = rows
        ;; !P.MULTI[3]  = 2 ;;Next page?

        page2Suff = ''

        END
        ELSE: BEGIN

           page2Suff = '-page2'

           IF KEYWORD_SET(save_ps) THEN BEGIN
              OUTPUT_SETUP,output_mode,plotDir,suff+page2Suff
           ENDIF

           columns      = 1
           rows         = 3
           !P.MULTI[0]  = 0
           !P.MULTI[1]  = columns
           !P.MULTI[2]  = rows
           !P.MULTI[3]  = 2 ;;Next page?
           ;; !P.MULTI  = [0, columns, rows, 0, 0]
           ;; !P.MULTI  = [0, columns, rows]
           ;; !P.MULTI  = [0, columns, rows, 1, 0]
           ;; !P.MULTI[0]  = 0

           IF ~KEYWORD_SET(save_ps) THEN BEGIN
              IF columns EQ 1 THEN BEGIN
                 WINDOW,2,XSIZE=700,YSIZE=800
              ENDIF ELSE BEGIN
                 WINDOW,2,XSIZE=1200,YSIZE=600
              ENDELSE
           ENDIF

           ;; k__yRange = [4e-3,1e1]
           kP__yRange = MINMAX(kP[inds])
           PLOT,freq[inds],kP[inds], $
                XTITLE=freqTitle, $
                ;; YRANGE=[4e-6,1e-2], $
                XRANGE=page2__freqRange, $
                YRANGE=kP__yRange, $
                XSTYLE=1, $
                YSTYLE=1, $
                XTICKLEN=1.0, $
                YTICKLEN=1.0, $
                XGRIDSTYLE=1, $
                YGRIDSTYLE=1, $
                ;; XLOG=1, $
                YLOG=1, $
                ;; YTITLE=font + 'k!Dperp!N (m!U-1!N)', $
                ;; YTITLE=font + 'k!Dperp!N (km!U-1!N)', $
                YTITLE=font + 'k!D' + perpAll + '!N (km!U-1!N)', $
                ;; XTICKFORMAT="(A1)", $
                CHARSIZE=cs, $
                SYMSIZE=symSize

           ;; smkP = SMOOTH(kP[inds],smInd,EDGE_TRUNCATE=edge_truncate,EDGE_MIRROR=edge_mirror,EDGE_WRAP=edge_wrap)
           IF KEYWORD_SET(plot_smoothed_ks) THEN BEGIN
              OPLOT,freq[inds],smooth_kP[inds], $
                    COLOR=250
           ENDIF

           IF KEYWORD_SET(overplot_doubly_smoothed) THEN BEGIN
              dbSmkP = SMOOTH(kP[inds],dbSmInd,EDGE_TRUNCATE=edge_truncate,EDGE_MIRROR=edge_mirror,EDGE_WRAP=edge_wrap)
              OPLOT,freq[inds],dbSmkP, $
                    COLOR=90
           ENDIF

           overplot_k_turbulence = 1
           IF KEYWORD_SET(overplot_k_turbulence) THEN BEGIN
              ;; kDoppl = (freq[inds]^(1.7)/10000)
              GET_FA_ORBIT,Tarr,/DEFINITIVE,/ALL,/TIME_ARRAY
              GET_DATA,'fa_vel',DATA=vel

              ;; speed = SQRT(vel.y[*,0]^2+vel.y[*,1]^2+vel.y[*,2]^2)*1000.0
              speed = SQRT(vel.y[*,0]^2+vel.y[*,1]^2+vel.y[*,2]^2)
              ;; avgSpeed = MEAN(speed)-6000.
              avgSpeed = MEAN(speed)

              kMajic = kx
              ;; kxTemp = kx

              ;; CASE 1 OF
              IF KEYWORD_SET(fitline__use_abs) THEN BEGIN
                 kMajic = ABS(kMajic)
              ENDIF ELSE BEGIN
                 chuck = WHERE(kMajic GT 0.0,nPosKx)
                 chuck = WHERE(kMajic LT 0.0,nNegKx)
                 IF nNegKx GT nPosKx THEN kMajic *= -1.
              ENDELSE
              IF KEYWORD_SET(fitline__use_smoothed) THEN BEGIN
                 kMajic = SMOOTH(kMajic,smInd)
              END
              ;; kMajic = SMOOTH((KEYWORD_SET(plot_abs_smoothed_ks) ? $
              ;;                  ABS(kP) : kP), $
              ;;                 smInd, $
              ;;                 EDGE_TRUNCATE=edge_truncate, $
              ;;                 EDGE_MIRROR=edge_mirror, $
              ;;                 EDGE_WRAP=edge_wrap)

              ;; wDoppl    = kMajic*(avgSpeed)/(2.*!PI)
              ;; fDoppl    = wDoppl

              minDoppl  = FLOOR(MIN(freq))
              maxDoppl  = CEIL(MAX(freq))
              fDoppl    = FINDGEN((maxDoppl-minDoppl)*20.+1)*0.05+minDoppl
              kDoppl    = fDoppl/avgSpeed

              want_kx   = 1
              IF KEYWORD_SET(want_kx) THEN BEGIN
                 ;; yArg   = ABS(kx)
                 yArg   = kx
                 yTito  = 'k!Dx!N (km!U-1!N)'
              ENDIF ELSE BEGIN
                 yArg   = ABS(ky)
                 yTito  = 'k!Dy!N (km!U-1!N)'
              ENDELSE

              k__yRange = [MIN(yArg),MAX(yArg)] 
              PLOT,freq[inds],yArg[inds], $
                   XTITLE=freqTitle, $
                   YTITLE=yTito, $
                   XRANGE=page2__freqRange, $
                   ;; YRANGE=k__yRange, $
                   XSTYLE=1, $
                   ;; YSTYLE=1, $
                   ;; XLOG=1, $
                   ;; YLOG=1, $
                   ;;NWO
                   ;; YRANGE=[(-1.)*MIN(k__yRange),MAX(k__yRange)], $
                   YRANGE=k__yRange, $
                   YLOG=0, $
                   XTICKLEN=1.0, $
                   YTICKLEN=1.0, $
                   XGRIDSTYLE=1, $
                   YGRIDSTYLE=1, $
                   CHARSIZE=cs, $
                   SYMSIZE=symSize

              ;; OPLOT,fDoppl[inds],kMajic[inds], $
              ;; OPLOT,freq[inds],kDoppl[inds], $
              OPLOT,fDoppl,kDoppl, $
                    COLOR=110

              OPLOT,fDoppl,(-1.)*kDoppl, $
                    COLOR=110

              kDopplDat = freq/avgSpeed
              posii     = WHERE((yArg[inds] GT 0) AND (yArg[inds] GT kDopplDat[inds]),nPosii)
              negii     = WHERE((yArg[inds] LT 0) AND (yArg[inds] LT (-1.)*kDopplDat[inds]),nNegii)

              dopplInds = !NULL
              IF nPosii GT 0 THEN BEGIN
                 dopplInds = [dopplInds,inds[posii]]
              ENDIF
              IF nNegii GT 0 THEN BEGIN
                 dopplInds = [dopplInds,inds[negii]]
              ENDIF

              OPLOT,freq[dopplInds],yArg[dopplInds],COLOR=230,PSYM=1

              add_Doppler_fit_string = 0
              IF KEYWORD_SET(add_Doppler_fit_string) THEN BEGIN
                 fitInds      = CGSETINTERSECTION(fitInds,WHERE(kMajic GT 0.00))
                 params       = LINFIT(ALOG10(freq[fitInds]),ALOG10(kMajic[fitInds]),YFIT=kxFitter)
                 corr         = LINCORR(ALOG10(freq[fitInds]),ALOG10(kMajic[fitInds]),T_STAT=t_stat)

                 params       = LINFIT(ALOG10(freq[fitInds]),ALOG10(ABS(kx[fitInds])),YFIT=kxFitter)
                 corr         = LINCORR(ALOG10(freq[fitInds]),ALOG10(ABS(kx[fitInds])),T_STAT=t_stat)

                 xFit         = 10.^((INDGEN(10))/ $
                                     10.*(ALOG10(MAX(freq[fitInds]))-ALOG10(MIN(freq[fitInds])))+$
                                     ALOG10(MIN(freq[fitInds])))
                 kxFit        = 10.^(params[1] * ALOG10(xFit) + params[0])
                 kxFitter     = 10.^(params[1] * ALOG10(freq[fitInds]) + params[0])

                 ;; slopeString  = STRING(FORMAT='(A-10,T15,F7.3)',"slope  =",params[1])
                 ;; corrString   = STRING(FORMAT='(A-10,T15,F7.3)',"r      =",corr[0])
                 ;; tString      = STRING(FORMAT='(A-10,T15,F7.3)',"t-test =",t_stat)
                 ;; txOutSize    = cs
                 ;; XYOUTS,0.2,0.89,slopeString,/NORMAL,CHARSIZE=txOutSize
                 ;; XYOUTS,0.2,0.86,corrString,/NORMAL,CHARSIZE=txOutSize
                 ;; XYOUTS,0.2,0.83,tString,/NORMAL,CHARSIZE=txOutSize
                 ;; XYOUTS,0.2,0.59,slopeString,/NORMAL,CHARSIZE=txOutSize
                 ;; XYOUTS,0.2,0.56,corrString,/NORMAL,CHARSIZE=txOutSize
                 ;; XYOUTS,0.2,0.53,tString,/NORMAL,CHARSIZE=txOutSize

                 ;; OPLOT,xFit,kxFit,COLOR=40

              ENDIF
           ENDIF

        END
     ENDCASE
     
     ;;Kperp angle plot     
     yMajDiv        = 45
     yMinor         = 6 * yMajDiv / 90 ;Should give spacing every 15 degrees no matter how we slice it
     IF KEYWORD_SET(kP__angleRange) THEN BEGIN

        CASE 1 OF
           (MIN(kP__angleRange) MOD 45 EQ 0): BEGIN

           END
           (MIN(kP__angleRange) MOD 60 EQ 0): BEGIN
              yMajDiv        = 60
              yMinor         = 3
           END
           ELSE: BEGIN
           END
        ENDCASE

        yARange     = kP__angleRange

        nTickV      = (MAX(kP__angleRange)-MIN(kP__angleRange))/yMajDiv+1
        yTickV      = INDGEN(nTickV)*yMajDiv+MIN(kP__angleRange) MOD 360
        ;; yTickName   = STRING(FORMAT='('+STRCOMPRESS(nTickV,/REMOVE_ALL)+'(I4))',yTickV MOD 360) ;doesn't work
        yTickName   = STRING(FORMAT='(I4)',yTickV MOD 360)

     ENDIF ELSE BEGIN

        yARange     = [-180,180]
        yTickV      = [-180,-90,0,90,180]

        IF rotate_kPA THEN BEGIN
           yARange += 180
           yTickV  += 180
        ENDIF

     ENDELSE


     PLOT,freq[inds],kPAngle[inds], $
          XTITLE=KEYWORD_SET(football_layout) ? !NULL : freqTitle, $
          ;; YTITLE='!4h!X!Dk!Dperp!N', $
          YTITLE='!4h!X!Dk!D' + perpAll + '!N', $
          XRANGE=(KEYWORD_SET(football_layout) OR KEYWORD_SET(PRE_VIII_layout))? page1__freqRange : page2__freqRange, $
          YRANGE=yARange, $
          XSTYLE=1, $
          YSTYLE=9, $
          XTICKLEN=1.0, $
          ;; YTICKLEN=1.0, $
          XGRIDSTYLE=1, $
          ;; YGRIDSTYLE=1, $
          YTICKS=N_ELEMENTS(yTickV)-1, $
          YTICKV=yTickV, $
          YTICKNAME=yTickName, $
          XMINOR=5, $
          YMINOR=yMinor, $
          CHARSIZE=cs, $
          ;; SYMSIZE=symSize, $
          POSITION=KEYWORD_SET(football_layout) ? poskPA : !NULL, $
          NOERASE=KEYWORD_SET(football_layout), $
          /NODATA

     IF KEYWORD_SET(football_layout) THEN BEGIN

        ;;And error angle
        IF ~KEYWORD_SET(football_no_errAngle) THEN BEGIN

           use_crossSym = 1
           
           IF KEYWORD_SET(use_crossSym) THEN BEGIN
              USERSYM,crossSym  ;to cross
           ENDIF
           OPLOT,freq[inds],errAngle[inds]*180.D/!PI, $
                 COLOR=thetaErrCol, $
                 PSYM=KEYWORD_SET(use_crossSym) ? 8 : !NULL

           ;;Now labels
           rango       = (MAX(yARange)-MIN(yARange))
           ;; distVec     = REVERSE(INDGEN(3)*spaceFrac - 2.*spaceFrac + 1)

           legYSymPos  = 0.925* rango + MIN(yARange)
           legYPos     = 0.900* rango + MIN(yARange)

           IF KEYWORD_SET(use_crossSym) THEN BEGIN
              ;;For cross sym
              legXPos1    = 0.68*((MAX(freq)-MIN(freq))+MIN(freq))
              legXSymPos1 = 0.66*((MAX(freq)-MIN(freq))+MIN(freq))

              PLOTS,legXSymPos1,legYSymPos[0], $
                    COLOR=thetaErrCol, $
                    PSYM=8

              ;;back to perpSym
              USERSYM,perpSym
           ENDIF ELSE BEGIN
              legXPos1    = 0.9*((MAX(freq)-MIN(freq))+MIN(freq))
              legXPos2    = 0.95*((MAX(freq)-MIN(freq))+MIN(freq))
              legXSymPos1 = 0.83*((MAX(freq)-MIN(freq))+MIN(freq))
              legXSymPos2 = 0.88*((MAX(freq)-MIN(freq))+MIN(freq))

              PLOTS,legXSymPos1,legYSymPos[0],COLOR=thetaErrCol
              PLOTS,legXSymPos2,legYSymPos[0],COLOR=thetaErrCol,/CONTINUE
           ENDELSE
           
           XYOUTS,legXPos1,legYPos[0],'!4h!X!Derr!N', $
                  CHARSIZE=cs, $
                  COLOR=thetaErrCol

        ENDIF

        AXIS,YAXIS=0, $
             YTITLE='!4h!X!Dk!D' + perpAll + '!N', $
             YRANGE=yARange, $
             XSTYLE=1, $
             YSTYLE=9, $
             YTICKS=N_ELEMENTS(yTickV)-1, $
             YTICKV=yTickV, $
             YTICKNAME=yTickName, $
             XMINOR=5, $
             YMINOR=yMinor, $
             CHARSIZE=cs, $
             /NOERASE, $
             /NODATA, $
             /SAVE

        ;;Best for last!
        IF ~KEYWORD_SET(only_ks_below_magErr_thresh) THEN BEGIN
           OPLOT,freq[inds],kPAngle[inds]
        ENDIF

        ;;Because of MARK_KS_BELOW_MAGERR_THRESH
        IF nCoolME_k GT 0 THEN BEGIN

           OPLOT,freq[coolME_k_i],kPAngle[coolME_k_i], $
                 ;; COLOR=nei, $
                 PSYM=pSym, $
                 SYMSIZE=symSize

        ENDIF

        IF nCoolEA_k GT 0 THEN BEGIN

           OPLOT,freq[coolEA_k_i],kPAngle[coolEA_k_i], $
                 ;; COLOR=nei, $
                 PSYM=pSymEA, $
                 SYMSIZE=symSize

        ENDIF

        IF nCoolPE_k GT 0 THEN BEGIN

           OPLOT,freq[coolPE_k_i],kPAngle[coolPE_k_i], $
                 ;; COLOR=nei, $
                 PSYM=pSymPE, $
                 SYMSIZE=symSize

        ENDIF

        IF nCoolMEEA_k GT 0 THEN BEGIN

           PLOTS,freq[coolMEEA_k_i],kPAngle[coolMEEA_k_i], $
                 ;; COLOR=nei, $
                 PSYM=pSymMEEA, $
                 SYMSIZE=symSize, $
                 NOCLIP=0


        ENDIF

        IF ~KEYWORD_SET(football_no_magErr) THEN BEGIN

           AXIS,YAXIS=1,YRANGE=[MIN(magErr),MAX(magErr)], $
                YSTYLE = 1, $
                YTITLE = '|Relative Error|', $
                ;; YTITLE = 'ABS(|J!Dpred!N|-|J!Dob!N|) (Normed)', $
                CHARSIZE=cs, $
                /SAVE, $
                /NOERASE, $
                COLOR=240

           ;; oldLine = !P.LINESTYLE
           ;; !P.LINESTYLE = [2,'AAAA'X]
           

           ;; jLineStyle = 4

           irregular = 1
           also_PE   = 1
           CASE 1 OF
              KEYWORD_SET(irregular): BEGIN

                 majicNivelME = KEYWORD_SET(mark_ks_below_magErr_thresh) ? mark_ks_below_magErr_thresh : 0.05
                 majicNivelPE = KEYWORD_SET(mark_ks_below_phaseErr_thresh) ? mark_ks_below_phaseErr_thresh : 45

                 CASE 1 OF
                    KEYWORD_SET(also_PE): BEGIN
                       lowErr_ii  = WHERE( ( magErr[inds] LE majicNivelME                 ) AND $
                                           ;; ( (NORMALIZE_ANGLE(phaseErr[2,inds],-!PI,!PI)*180./!PI) LE  majicNivelPE ),nLowErr, $
                                           ( ABS(phaseErr[2,inds]*180./!PI) LE  majicNivelPE ),nLowErr, $
                                          COMPLEMENT=highErr_ii, $
                                          NCOMPLEMENT=nHighErr)
                    END
                    ELSE: BEGIN
                       lowErr_ii  = WHERE(magErr[inds] LE majicNivelME,nLowErr, $
                                          COMPLEMENT=highErr_ii, $
                                          NCOMPLEMENT=nHighErr)
                    END
                 ENDCASE

                 IF nLowErr GT 1 THEN BEGIN

                    ;;If the line below is commented, you'll get a perpSym
                    USERSYM,dashSym ;Just a dash

                    OPLOT,freq[inds[lowErr_ii]],magErr[inds[lowErr_ii]], $
                          ;; LINESTYLE=jLineStyle, $
                          ;; LINESTYLE=, $
                          COLOR=240, $
                          PSYM=pSymMagErr

                 ENDIF

                 IF nHighErr GT 1 THEN BEGIN

                    USERSYM,dashSym ;Just a dash

                    OPLOT,freq[inds[highErr_ii]],magErr[inds[highErr_ii]], $
                          ;; LINESTYLE=jLineStyle, $
                          ;; LINESTYLE=, $
                          COLOR=240, $
                          PSYM=pSymMagErr

                 ENDIF

              END
              ELSE: BEGIN

                 OPLOT,freq[inds],magErr[inds], $
                       ;; LINESTYLE=jLineStyle, $
                       ;; LINESTYLE=, $
                       COLOR=240, $
                       PSYM=pSymMagErr

              END
           ENDCASE

           ;; !P.LINESTYLE = oldLine

        ENDIF

     ENDIF


     IF KEYWORD_SET(kP_angle__plot_smoothed) THEN BEGIN

        CASE SIZE(kP_angle__plot_smoothed) OF
           8: BEGIN
              OPLOT,kP_angle__plot_smoothed.freq,kP_angle__plot_smoothed.angle, $
                    COLOR=250
           END
           ELSE: BEGIN
              OPLOT,smooth_freq[smooth_inds],smooth_kPAngle[smooth_inds], $
                    COLOR=250
           END
        ENDCASE

     ENDIF

     IF KEYWORD_SET(make_kPAngle_special) THEN BEGIN

        GET_STREAKS,kPA_spI, $
                    START_I=strt_ii, $
                    STOP_I=stop_ii, $
                    MIN_STREAK_TO_KEEP=2, $
                    N_STREAKS=nStreaks, $
                    /QUIET

        FOR jj=0,nStreaks-1 DO BEGIN

           tmpI = [kPA_spI[strt_ii[jj]]:kPA_spI[stop_ii[jj]]]

           OPLOT,freq[tmpI],kPAngle[tmpI], $
                 ;; PSYM=1, $
                 LINESTYLE=kPA_specialSty, $
                 COLOR=specialColor

        ENDFOR

     ENDIF


     IF ~(KEYWORD_SET(PRE_VIII_layout) OR KEYWORD_SET(football_layout)) THEN BEGIN

        IF KEYWORD_SET(kP_angle__plot_smoothed) THEN BEGIN
           OPLOT,freq[dopplInds],smooth_kPAngle[dopplInds], $
                 COLOR=230, $
                 PSYM=1
        ENDIF

        PRINT,"Mean theta exceeding doppl: ",MEAN(kPAngle[dopplInds])
        PRINT,"Medn theta exceeding doppl: ",MEDIAN(kPAngle[dopplInds])

     ENDIF

     IF KEYWORD_SET(overplot_doubly_smoothed) THEN BEGIN
        dbSmkPAngle = SMOOTH(kPAngle,dbSmInd,EDGE_TRUNCATE=edge_truncate,EDGE_MIRROR=edge_mirror,EDGE_WRAP=edge_wrap)
        OPLOT,freq[inds],dbSmkPAngle[inds], $
              COLOR=90
     ENDIF

     CASE 1 OF
        ;; KEYWORD_SET(football_layout): BEGIN

        ;; END
        KEYWORD_SET(third_page) AND ~KEYWORD_SET(PRE_VIII_layout): BEGIN
           IF KEYWORD_SET(save_ps) THEN BEGIN
              CONCLUDE_OUTPUT,output_mode,plotDir,suff+page2Suff, $
                              ;; CONCLUDE_OUTPUT,output_mode,plotDir,suff, $
                              TO_PDF=to_pdf, $
                              PDF_TRANSPARENCY_LEVEL=pdf_transparency, $
                              REMOVE_EPS=remove_eps

              OUTPUT_SETUP,output_mode,plotDir,suff+'-page3',saveDir
           ENDIF

           columns      = 1
           rows         = 1
           !P.MULTI[0]  = 0
           !P.MULTI[1]  = columns
           !P.MULTI[2]  = rows
           !P.MULTI[3]  = 4 ;;Nextnext page?

           ;; kP__yRange = [4e-3,1e1]
           PLOT,kPAngle[inds],kP[inds], $
                ;; XTITLE='!4h!X!Dk!Dperp!N', $
                ;; YTITLE=font + 'k!Dperp!N (km!U-1!N)', $
                XTITLE='!4h!X!Dk!D' + perpAll + '!N', $
                YTITLE=font + 'k!D' + perpAll + '!N (km!U-1!N)', $
                XRANGE=yARange, $
                YRANGE=kP__yRange, $
                XSTYLE=1, $
                YSTYLE=1, $
                LINESTYLE=0, $
                PSYM=pSym, $
                XTICKLEN=1.0, $
                YTICKLEN=1.0, $
                XGRIDSTYLE=1, $
                YGRIDSTYLE=1, $
                CHARSIZE=cs, $
                SYMSIZE=symSize

           IF KEYWORD_SET(kP_angle__plot_smoothed) THEN BEGIN
              OPLOT,smooth_kPAngle[inds],smooth_kP[inds], $
                    COLOR=250, $
                    LINESTYLE=0, $
                    PSYM=pSym

           ENDIF

           IF KEYWORD_SET(save_ps) THEN BEGIN
              CONCLUDE_OUTPUT,output_mode,plotDir,suff+'-page3', $
                              ;; CONCLUDE_OUTPUT,output_mode,plotDir,suff, $
                              TO_PDF=to_pdf, $
                              PDF_TRANSPARENCY_LEVEL=pdf_transparency, $
                              REMOVE_EPS=remove_eps

              ;;Smash pages together
              IF KEYWORD_SET(to_pdf) THEN BEGIN
                 fPref = plotDir + suff

                 SPAWN,'pdfunite ' + fPref + '.pdf' + ' ' + $
                       fPref + '-page2' + '.pdf' + ' ' + $
                       fPref + '-page3' + '.pdf' + ' ' + $
                       fPref + '-all' + '.pdf'

                 SPAWN,'mv ' + fPref + '-all' + '.pdf' + ' ' + $
                       fPref + '.pdf'

                 SPAWN,'rm ' + fPref + '-page2' + '.pdf'
                 SPAWN,'rm ' + fPref + '-page3' + '.pdf'
              ENDIF

           ENDIF

        END
        ELSE: BEGIN

           IF KEYWORD_SET(save_ps) THEN BEGIN

              IF KEYWORD_SET(football_layout) THEN page2suff = ''

              CONCLUDE_OUTPUT,output_mode,plotDir,suff+page2suff, $
                              ;; CONCLUDE_OUTPUT,output_mode,plotDir,suff, $
                              TO_PDF=to_pdf, $
                              PDF_TRANSPARENCY_LEVEL=pdf_transparency, $
                              REMOVE_EPS=remove_eps

              ;;Smash pages together
              IF KEYWORD_SET(to_pdf) THEN BEGIN
                 fPref = plotDir + suff

                 CASE 1 OF
                    KEYWORD_SET(football_layout): BEGIN

                    END
                    KEYWORD_SET(PRE_VIII_layout): BEGIN

                    END
                    ELSE: BEGIN

                       SPAWN,'pdfunite ' + fPref + '.pdf' + ' ' + $
                             fPref + '-page2' + '.pdf' + ' ' + $
                             fPref + '-all' + '.pdf'

                       SPAWN,'mv ' + fPref + '-all' + '.pdf' + ' ' + $
                             fPref + '.pdf'

                       SPAWN,'rm ' + fPref + '-page2' + '.pdf'

                    END
                 ENDCASE

              ENDIF

           ENDIF

        END
     ENDCASE

  ENDELSE

END
