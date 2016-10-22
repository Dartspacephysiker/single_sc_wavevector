; This IDL code follows Section 4 in P. M. Bellan's JGR manuscript
; It gives examples of how the wavevector k(omega) can be determined
; from measurements of B(t) and J(t) by a single spacecraft
; The examples use synthetic data, but the code could be easily modified
; to use actual data.
;The code consists of:
; 1. Main program
; 2. Function, "correlation" that calculates cross-correlations
; 3. Procedure, "output_setup" that does housekeeping for setting up
;    output to screen, printer, or eps file
; 4. Procedure, "conclude_output" that does housekeeping for output at end
;    of code.
; 5. Procedure, "add_noise" that creates noise from random number generator
;    and for example 3 only adds to J and B signals
; Following IDL protocol, the function and procedures are listed before
; the main program.

;calculate cross-correlation of two functions <f(t)*g(t+tau)>
; written by P. M. Bellan, April 2016, MS 128-95 Caltech, Pasadena CA 91125
FUNCTION CORRELATION,f,g,T
  autocorr = FLTARR(T)
  FOR tau=0,T-1 do BEGIN
     FOR TT=0,T-1 do BEGIN
        TT_plus_tau= (TT+tau) mod T ; use modulo arithmetic to wrap back to beginning
        autocorr[tau] = autocorr[tau]+f[TT]*g[TT_plus_tau]
     ENDFOR
  ENDFOR
  RETURN,autocorr/T
END

;********************************************************************
; written by P. M. Bellan, April 2016, MS 128-95 Caltech, Pasadena CA 91125
; prefix gives file path for eps output, data input/output
; user should change to be appropriate for user's computer
PRO OUTPUT_SETUP,mode,plotDir,suff,saveDir
  ;; prefix='D:\CORR\MNSCRPTS\2016-JGR-spacecraft-current-wavevector\2106-Vinas-kvect-revised\'
  SET_PLOT_DIR,plotDir,/FOR_SINGLE_SC_WVEC,/ADD_TODAY
  saveDir = '/SPENCEdata/Research/Satellites/FAST/single_sc_wavevector/'
  PRINT,'plotDir: ',plotDir
  ;;
  filename = plotDir+ $
             'FAST_'+suff
  PRINT,'print plot to file ',filename
  ;;SET UP FOR PLOTTING ON:
  ;; SCREEN,          set hardcopy = 0
  ;; PRINTER,         set hardcopy = 1
  ;; POSTSCRIPT FILE, set hardcopy = 2
  hardcopy = mode
  IF hardcopy EQ 1 THEN SET_PLOT,'printer' ;PRINTER
  IF hardcopy EQ 1 THEN DEVICE,YSIZE=25,YOFFSET=0
  IF hardcopy NE 1 THEN SET_PLOT,'X'
  IF hardcopy EQ 2 THEN BEGIN      ;;POSTSCRIPT FILE
     !P.CHARTHICK = 3
     !P.THICK = 3

     ;; @startup
     ;; POPEN,filename,XSIZE=10,YSIZE=10

     SET_PLOT,'PS'
     DEVICE,FILE=filename+'.eps',/ENCAPSUL,XSIZE=10,YSIZE=10,/INCHES,YOFFSET=2,/COLOR
  ENDIF
END

;********************************************************************
; written by P. M. Bellan, April 2016, MS 128-95 Caltech, Pasadena CA 91125
; place at end of code to do housekeeping
PRO CONCLUDE_OUTPUT,mode,plotDir
  hardcopy = mode
  IF hardcopy EQ 1 THEN DEVICE,/CLOSE
  IF hardcopy EQ 1 THEN SET_PLOT,'X'
  IF hardcopy EQ 2 THEN BEGIN
     DEVICE,/CLOSE
     ;; PCLOSE
     SET_PLOT,'X'
  ENDIF
  IF hardcopy EQ 1 THEN SET_PLOT,'printer'
  IF hardcopy EQ 1 THEN DEVICE,YSIZE=25,YOFFSET=5
  IF hardcopy ne 1 THEN SET_PLOT,'X'
  PRINT, 'FINISHED'
END

;********************************************************************
PRO ADD_NOISE,Bx,By,Bz,Jx,Jy,Jz,T
  ;;create and add noise IF example = 3

  Bsqtotal = 0.0
  Jsqtotal = 0.0

  ;;calculate rms amplitude of B and J
  FOR TT=0,T-1 DO BEGIN         ; TT is time
     Bsqtotal = Bsqtotal+Bx[TT]^2+By[TT]^2+Bz[TT]^2
     Jsqtotal = Jsqtotal+Jx[TT]^2+Jy[TT]^2+Jz[TT]^2
  ENDFOR

  Brms            = SQRT(Bsqtotal/T)
  Jrms            = SQRT(Jsqtotal/T)

  ;; create noise from random number generator randomu
  Bnoise          = FLTARR(T,3)
  Jnoise          = FLTARR(T,3)
  FOR i=0,T-1 DO BEGIN
     Bnoise[i,0]  = RANDOMU(seed)
     Jnoise[i,0]  = RANDOMU(seed)
     Bnoise[i,1]  = RANDOMU(seed)
     Jnoise[i,1]  = RANDOMU(seed)
     Bnoise[i,2]  = RANDOMU(seed)
     Jnoise[i,2]  = RANDOMU(seed)
  ENDFOR

  ;;renormalize noise to have unity rms
  Bnoisetot     = 0
  Jnoisetot     = 0
  FOR TT=0,T-1 DO BEGIN         ; TT is time
     Bnoisetot  = Bnoisetot + $
                  Bnoise[TT,0]*Bnoise[TT,0] + $
                  Bnoise[TT,1]*Bnoise[TT,1] + $
                  Bnoise[TT,2]*Bnoise[TT,2]

     Jnoisetot  = Jnoisetot + $
                  Jnoise[TT,0]*Jnoise[TT,0] + $
                  Jnoise[TT,1]*Jnoise[TT,1] + $
                  Jnoise[TT,2]*Jnoise[TT,2]
  ENDFOR

  Bnoiserms     = SQRT(Bnoisetot/T)
  Jnoiserms     = SQRT(Jnoisetot/T)

  ;; below sets noise to have unity rms
  Bnoise[*,0]   = Bnoise[*,0]/Bnoiserms
  Bnoise[*,1] = Bnoise[*,1]/Bnoiserms
  Bnoise[*,2] = Bnoise[*,2]/Bnoiserms

  Jnoise[*,0] = Jnoise[*,0]/Jnoiserms
  Jnoise[*,1] = Jnoise[*,1]/Jnoiserms
  Jnoise[*,2] = Jnoise[*,2]/Jnoiserms

  ;;amount of noise set to 50% of rms signal
  Bnoisecoeff   = 0.5*Brms
  Jnoisecoeff   = 0.5*Jrms

  ;;add noise
  FOR i=0,T-1 DO BEGIN          ; add noise contributions for each time
     Bx[i]      = Bx[i] + Bnoisecoeff*Bnoise[i,0]
     By[i]      = By[i] + Bnoisecoeff*Bnoise[i,1]
     Bz[i]      = Bz[i] + Bnoisecoeff*Bnoise[i,2]
     Jx[i]      = Jx[i] + Jnoisecoeff*Jnoise[i,0]
     Jy[i]      = Jy[i] + Jnoisecoeff*Jnoise[i,1]
     Jz[i]      = Jz[i] + Jnoisecoeff*Jnoise[i,2]
  ENDFOR

  RETURN

END

PRO  SETUP_EXAMPLE,T,TArr,Bx,By,Bz,Jx,Jy,Jz,unitFactor,sPeriod,saveVar, $
                   EXAMPLE=example

  T         = 200               ; number of elements in time domain vectors
  TArr      = INDGEN(T)
  sPeriod   = 1

  ;; k_vector : first argument is for frequency, second is for Cartesian coordinate
  kvec      = FLTARR(T/2+1,3)

  omega     = INDGEN(T/2+1)*2*!Pi/T

  ;;setup k vectors for either example 1 or 2
                                ;first example
  IF example EQ 1 THEN BEGIN
     kx      = FLTARR(T/2+1)
     ky      = kx
     kz      = kx

     kx[10]  = -2
     ky[50]  = 8
     kz[70]  = 5
  ENDIF

  ;; second and third examples, continuum of k components
  ;;constructed using arbitrary complicated function
  IF example EQ 2 OR example EQ 3 THEN BEGIN

     ;; pick some arbitrary function for magnitude for k
     kmag    = SQRT(omega)

     ;;choose azimuthal spherical coordinate angle of unit k vector
     phi     = 2*!pi*INDGEN(T/2+1)/T

     ;;choose polar spherical coordinate angle of unit k vector
     theta   = !pi*INDGEN(T/2+1)/T

     ;; construct Cartesian unit vectors for k
     x_unit  = SIN(theta)*COS(phi)
     y_unit  = SIN(theta)*SIN(phi)
     z_unit  = COS(phi)

     ;;define components of k vector
     kx      = kmag*x_unit
     ky      = kmag*y_unit
     kz      = kmag*z_unit

  ENDIF

  ;;plot k vectors that have been set up
  IF example EQ 1 THEN BEGIN
     ysize  = 10
  ENDIF
  IF example EQ 2 or example EQ 3 THEN BEGIN
     ysize  = 2
  ENDIF

  PLOT,kx,XTITLE='!4x!3T/(2!4p!3)',YTITLE='!8k!Dx!N',CHARSIZE=cs,YRANGE=[-ysize,ysize]
  PLOT,ky,XTITLE='!4x!3T/(2!4p!3)',YTITLE='!8k!Dy!N',CHARSIZE=cs,YRANGE=[-ysize,ysize]
  PLOT,kz,XTITLE='!4x!3T/(2!4p!3)',YTITLE='!8k!Dz!N',CHARSIZE=cs,YRANGE=[-ysize,ysize]

  ;;store k components in a vector for use later
  FOR i=0,T/2 DO BEGIN
     kvec[i,0]= kx[i]
     kvec[i,1]= ky[i]
     kvec[i,2]= kz[i]
  ENDFOR

  ;; set up k vector
  ;; define A vector potential cosine component
  ;; first argument is frequency, second is Cartesian component
  Avec_cos = FLTARR(T/2+1,3)    ; define A vector cos component (vector potential)
  Avec_sin = FLTARR(T/2+1,3)    ; define A vector sine component (vector potential)

  ;; seed for random number generator,
  ;; get different random number each time randomu is called
  seed = 2001L

  ;; make coefficients of vector potential random,
  ;;Avec_cos[i,0] means
  ;; coefficient of the x component with cosine time dependence for frequency i
  ;;make all vector potential components random numbers
  FOR I=0,T/2 DO BEGIN
     Avec_cos[i,0]  = RANDOMU(seed)
     Avec_sin[i,0]  = RANDOMU(seed)
     Avec_cos[i,1]  = RANDOMU(seed)
     Avec_sin[i,1]  = RANDOMU(seed)
     Avec_cos[i,2]  = RANDOMU(seed)
     Avec_sin[i,2]  = RANDOMU(seed)
  ENDFOR

  ;; B vector for a single frequency, args are frequency, Cartesian direction, time
  Bvec = FLTARR(T/2+1,3,T)
  Jvec = Bvec
  FOR TT=0,T-1 DO BEGIN         ;for each time
     FOR i=0,T/2 DO BEGIN       ;for each frequency component

        ;; temp A cosine component
        Avectemp_cos = [Avec_cos[i,0], Avec_cos[i,1], Avec_cos[i,2] ]

        ;; temp A sin e component
        Avectemp_sin = [Avec_sin[i,0], Avec_sin[i,1], Avec_sin[i,2] ]

        ;; tempo k vector
        kvectemp = [kvec[i,0], kvec[i,1], kvec[i,2] ]

        ;; k vector cross product with A cosine
        kcross_Acos = CROSSP(kvectemp,Avectemp_cos)

        ;; k vector cross product with A sine
        kcross_Asin = CROSSP(kvectemp,Avectemp_sin)

        ;; put in time dependence to get total B at a single time for freq i
        Bvec[i,0,tt] = -kcross_Acos[0]*SIN(-omega[i]*tt) $
                       + kcross_Asin[0]*COS(-omega[i]*tt)

        Bvec[i,1,tt] = -kcross_Acos[1]*SIN(-omega[i]*tt)$
                       + kcross_Asin[1]*COS(-omega[i]*tt)

        Bvec[i,2,tt] = -kcross_Acos[2]*SIN(-omega[i]*tt) $
                       + kcross_Asin[2]*COS(-omega[i]*tt)

        ;; k cross (k cross A) need for current
        kcrosskcross_Acos = CROSSP(kvectemp,kcross_Acos)
        kcrosskcross_Asin = CROSSP(kvectemp,kcross_Asin)

        ;;calculate current at a single time for freq i
        Jvec[i,0,tt] = -kcrosskcross_Acos[0]* $
                       COS(-omega[i]*tt)-kcrosskcross_Asin[0]*SIN(-omega[i]*tt)

        Jvec[i,1,tt] = -kcrosskcross_Acos[1]* $
                       COS(-omega[i]*tt)-kcrosskcross_Asin[1]*SIN(-omega[i]*tt)

        Jvec[i,2,tt] = -kcrosskcross_Acos[2]* $
                       COS(-omega[i]*tt)-kcrosskcross_Asin[2]*SIN(-omega[i]*tt)
     ENDFOR
  ENDFOR

  Bx = FLTARR(T)
  By = Bx
  Bz = Bx
  Jx = Bx
  Jy = Bx
  Jz = Bx

  ;; define dimensions of B and J components
  FOR TT=0,T-1 DO BEGIN         ; TT is time
     FOR i=0,T/2 DO BEGIN       ; sum up contributions from each frequency
        Bx[TT] = Bx[TT] +Bvec[i,0,TT]
        By[TT] = By[TT] +Bvec[i,1,TT]
        Bz[TT] = Bz[TT] +Bvec[i,2,TT]
        Jx[TT] = Jx[TT] +Jvec[i,0,TT]
        Jy[TT] = Jy[TT] +Jvec[i,1,TT]

        Jz[TT] = Jz[TT] +Jvec[i,2,TT]
     ENDFOR                     ;
  ENDFOR

  ;;create noise if example EQ 3
  IF example EQ 3 THEN add_noise, Bx,By,Bz,Jx,Jy,Jz,T

  ;;Write data to file
  ;; so that code can be easily modified to input data from other sources
  B_J_file = 'Bx_By_Bz_Jx_Jy_Jz.dat'

  filename = saveDir+B_J_file

  ;;write data to file
  OPENW,1,filename
  PRINT, 'opening ', filename

  PRINTF,1,example,T,FORMAT='(I5,1x,I5,1x)'
  PRINTF,1,Bx
  PRINTF,1,By
  PRINTF,1,Bz
  PRINTF,1,Jx
  PRINTF,1,Jy
  PRINTF,1,Jz
  CLOSE,1
  PRINT,'closing ',filename

END

PRO CHUNK_SAVE_FILE,T,TArr,Bx,By,Bz,Jx,Jy,Jz,unitFactor,sPeriod,saveVar, $
                    SAVFILE_T1=s_t1, $
                    SAVFILE_T2=s_t2, $
                    USE_TIMEBAR_TIME=use_timeBar_time, $
                    USE_J_TIME_SERIES=use_J_time_series, $
                    SMOOTH_J_DAT_TO_B=smooth_J_dat, $
                    STREAKNUM=streakNum, $
                    OUT_STREAKNUM=longestInd, $
                    USE_DB_FAC=use_dB_fac


  saveDir  = '/SPENCEdata/Research/Satellites/FAST/single_sc_wavevector/saves_output_etc/'
  saveFile = 'Chaston_et_al_2006--B_and_J.sav'
  saveFile = 'Chaston_et_al_2006--B_and_J--20161022--fixed_currents.sav'

  PRINT,"Restoring " + saveFile + ' ...'
  RESTORE,saveDir+saveFile

  bUnitFactor = -9 ;'cause nT
  jUnitFactor = -6 ;'cause microA/m^2

  unitFactor  = (10.D)^(jUnitFactor)/(10.D)^(bUnitFactor)

  saveVar  = 'dB_fac_V'
  IF KEYWORD_SET(use_dB_fac) THEN saveVar = 'dB_fac'
  ;; saveVar  = 'dB_fac_V'

  ;;Align time series

  CASE STRUPCASE(saveVar) OF
     'DB_FAC_V': BEGIN
        ;;   From UCLA_MAG_DESPIN: "Field-aligned velocity-based coordinates defined as:    "
        ;;x (ind 0)-along track ((BxV)xB),
        ;;y (ind 1)-cross track (BxV),
        ;;z (ind 2)-along B" (I added "ind" marks)
        dB = dB_fac_v
     END
     'DB_FAC': BEGIN
        ;;  Field-aligned coordinates defined as:
        ;;   z-along B, y-east (BxR), x-nominally out
        dB = dB_fac
     END
  ENDCASE

  ;; TArr = dB.x
  ;; bFactor = 1.e-9 ;Get 'em out of nT
  ;; bFactor = 1.e3
  bFactor = 1.
  Bx   = dB.y[*,0] * bFactor
  By   = dB.y[*,1] * bFactor
  Bz   = dB.y[*,2] * bFactor

  mag_sRate  = 1./(dB.x  [1:-1]-dB.x  [0:-2])
  eESA_sRate = 1./(Je_z.x[1:-1]-Je_z.x[0:-2])
  iESA_sRate = 1./(Ji_z.x[1:-1]-Ji_z.x[0:-2])

  ;; this       = VALUE_CLOSEST2(je_z.x,db.x)
  ;; PRINT,ABS(je_z.x[this]-db.x)

  ;; this       = VALUE_CLOSEST2(ji_z.x,db.x)
  ;; PRINT,ABS(ji_z.x[this]-db.x)


  ;;Using j or mag time series?
  CASE 1 OF
     KEYWORD_SET(use_J_time_series): BEGIN
        frame_of_ref = 2
     END
     ;; (KEYWORD_SET(use_mag_time_series): BEGIN
     ELSE: BEGIN
        frame_of_ref = 1
     END
  ENDCASE
  CASE frame_of_ref OF
     1: BEGIN

        IF KEYWORD_SET(smooth_J_dat) THEN BEGIN

           bro = ROUND_TO_NTH_DECIMAL_PLACE(db.x[1:-1]-db.x[0:-2],-5)
           bro = bro[WHERE(ABS(bro) LT 1)]
           distFreq = HISTOGRAM(bro,MIN=MIN(bro),BINSIZE=0.00001, $
                                REVERSE_INDICES=ri, $
                                LOCATIONS=locs)
           junk = MAX(distFreq,ind)
           sPeriod = DOUBLE(locs[ind])
           even_TS = MAKE_EVENLY_SPACED_TIME_SERIES(START_T=db.x[0], $
                                                    STOP_T=db.x[-1], $
                                                    DELTA_T=sPeriod)

           je_z_interp = MAKE_ARRAY(N_ELEMENTS(even_TS),VALUE=!VALUES.F_NaN)
           ji_z_interp = MAKE_ARRAY(N_ELEMENTS(even_TS),VALUE=!VALUES.F_NaN)

           ;;The new way to smooth
           je_z_interp = SMOOTH(je_z.y,3,/NaN)
           ji_z_interp = SMOOTH(ji_z.y,3,/NaN)

           je_z_interp = DATA_CUT({x:je_z.x,y:je_z_interp},db.x)
           ji_z_interp = DATA_CUT({x:ji_z.x,y:ji_z_interp},db.x)

           ;;The old way
           ;; GET_DOUBLE_STREAKS__NTH_DECIMAL_PLACE,je_z.x,-2,NPTS=100, $
           ;;                                       GAP_TIME=0.3, $
           ;;                                       START_I=strt_i, $
           ;;                                       STOP_I=stop_i, $
           ;;                                       STREAKLENS=streakLens, $
           ;;                                       T_STREAKLENS=tStreakLens, $
           ;;                                       /PRINT_START_STOP_TIMES

           ;; jeFin = FINITE(je_z.y)
           ;; jiFin = FINITE(ji_z.y)
           ;; PRINT,"Smoothing Je and Ji ..."
           ;; FOR k=0,N_ELEMENTS(even_TS)-1 DO BEGIN
           ;;    tmpCheck = WHERE((ABS(je_z.x-even_TS[k]) LE sPeriod) AND jeFin,nCheck)
           ;;    PRINT,nCheck
           ;;    IF nCheck GT 0 THEN BEGIN
           ;;       je_z_interp[k] = MEAN(je_z.y[tmpCheck])
           ;;    ENDIF

           ;;    tmpCheck = WHERE((ABS(ji_z.x-even_TS[k]) LE sPeriod) AND jiFin,nCheck)
           ;;    IF nCheck GT 0 THEN BEGIN
           ;;       ji_z_interp[k] = MEAN(ji_z.y[tmpCheck])
           ;;    ENDIF
           ;; ENDFOR
           ;; PRINT,'Bro'

        ENDIF ELSE BEGIN

           bro = ROUND_TO_NTH_DECIMAL_PLACE(db.x[1:-1]-db.x[0:-2],-5)
           bro = bro[WHERE(ABS(bro) LT 1)]
           distFreq = HISTOGRAM(bro,MIN=MIN(bro),BINSIZE=0.00001, $
                                REVERSE_INDICES=ri, $
                                LOCATIONS=locs)
           junk = MAX(distFreq,ind)
           sPeriod = DOUBLE(locs[ind])
           even_TS = db.x

           FA_FIELDS_COMBINE,{time:dB.x,comp1:dB.y[*,0]}, $
                             {time:Je_z.x,comp1:Je_z.y}, $
                             RESULT=Je_z_interp, $
                             ;; /INTERP, $
                             /SPLINE, $
                             DELT_T=1.5, $
                             /TALK

           FA_FIELDS_COMBINE,{time:dB.x,comp1:dB.y[*,0]}, $
                             {time:Ji_z.x,comp1:Ji_z.y}, $
                             RESULT=Ji_z_interp, $
                             ;; /INTERP, $
                             /SPLINE, $
                             DELT_T=1.5, $
                             /TALK
        ENDELSE
     END
     2: BEGIN
        bro = ROUND_TO_NTH_DECIMAL_PLACE(ji_z.x[1:-1]-ji_z.x[0:-2],-5)
        bro = bro[WHERE(ABS(bro) LT 1)]
        distFreq = HISTOGRAM(bro,MIN=MIN(bro),BINSIZE=0.00001, $
                             REVERSE_INDICES=ri, $
                             LOCATIONS=locs)
        junk = MAX(distFreq,ind)
        sPeriod = DOUBLE(locs[ind])
        even_TS = MAKE_EVENLY_SPACED_TIME_SERIES(START_T=ji_z.x[0], $
                                                 STOP_T=ji_z.x[-1], $
                                                 DELTA_T=sPeriod)

        ;; even_TS = even_TS[60:-60]

        FA_FIELDS_COMBINE,{time:even_TS,comp1:even_TS}, $
                          {time:dB.x,comp1:Bx}, $
                          RESULT=Bx_interp, $
                          ;; /INTERP, $
                          /SPLINE, $
                          DELT_T=1.5, $
                          /TALK

        FA_FIELDS_COMBINE,{time:even_TS,comp1:even_TS}, $
                          {time:dB.x,comp1:By}, $
                          RESULT=By_interp, $
                          ;; /INTERP, $
                          /SPLINE, $
                          DELT_T=1.5, $
                          /TALK

        FA_FIELDS_COMBINE,{time:even_TS,comp1:even_TS}, $
                          {time:dB.x,comp1:Bz}, $
                          RESULT=Bz_interp, $
                          ;; /INTERP, $
                          /SPLINE, $
                          DELT_T=1.5, $
                          /TALK

        FA_FIELDS_COMBINE,{time:even_TS,comp1:even_TS}, $
                          {time:Je_z.x,comp1:Je_z.y}, $
                          RESULT=Je_z_interp, $
                          ;; /INTERP, $
                          /SPLINE, $
                          DELT_T=1.5, $
                          /TALK

        FA_FIELDS_COMBINE,{time:even_TS,comp1:even_TS}, $
                          {time:Ji_z.x,comp1:Ji_z.y}, $
                          RESULT=Ji_z_interp, $
                          ;; /INTERP, $
                          /SPLINE, $
                          DELT_T=1.5, $
                          /TALK

        Bx = Bx_interp
        By = By_interp
        Bz = Bz_interp

     END
  ENDCASE

  ;;Interp over gaps if the space isn't very big

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;METHOD 1

  ;; DEAL_WITH_BADNESS,je_z_interp,je_z_improv
  ;; DEAL_WITH_BADNESS,ji_z_interp,ji_z_improv

  ;; good_i = WHERE(FINITE(je_z_improv) AND FINITE(ji_z_improv),nGood)
  ;; GET_STREAKS,good_i,START_i=strt_ii,STOP_I=stop_ii

  ;; ;; good_i = (good_i[strt_ii[0]:stop_ii[0]])[0:-20]
  ;; good_i = (good_i[strt_ii[0]:stop_ii[0]])

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;METHOD 2

  je_z_improv  = je_z_interp
  ji_z_improv  = ji_z_interp

  IF KEYWORD_SET(use_timeBar_time) THEN BEGIN
     timesBarStr = ['1998-05-04/06:44:46','1998-05-04/06:44:56']
     s_t1        = STR_TO_TIME(timesBarStr[0])
     s_t2        = STR_TO_TIME(timesBarStr[1])
  ENDIF

  CASE 1 OF
     ;; KEYWORD_SET(use_timeBar_time): BEGIN
     (KEYWORD_SET(s_t1) AND KEYWORD_SET(s_t2)): BEGIN
        good_i      = WHERE(FINITE(je_z_improv) AND FINITE(ji_z_improv) AND $
                            (even_TS GE s_t1) AND (even_TS LE s_t2),nGood)
        GET_STREAKS,good_i,START_i=strt_ii,STOP_I=stop_ii,OUT_STREAKLENS=streakLens

     END
     KEYWORD_SET(s_t1): BEGIN
        good_i      = WHERE(FINITE(je_z_improv) AND FINITE(ji_z_improv) AND $
                            (even_TS GE s_t1),nGood)
        GET_STREAKS,good_i,START_i=strt_ii,STOP_I=stop_ii,OUT_STREAKLENS=streakLens

     END
     KEYWORD_SET(s_t2): BEGIN
        good_i      = WHERE(FINITE(je_z_improv) AND FINITE(ji_z_improv) AND $
                            (even_TS LE s_t2),nGood)
        GET_STREAKS,good_i,START_i=strt_ii,STOP_I=stop_ii,OUT_STREAKLENS=streakLens
     END
     ELSE: BEGIN
        good_i      = WHERE(FINITE(je_z_improv) AND FINITE(ji_z_improv),nGood)
        GET_STREAKS,good_i,START_i=strt_ii,STOP_I=stop_ii,OUT_STREAKLENS=streakLens

        ;;Did user provide a streak?
        IF N_ELEMENTS(streakNum) GT 0 THEN longestInd = streakNum

     END
  ENDCASE

  IF N_ELEMENTS(streakLens) GT 0 THEN BEGIN
     long     = MAX(streakLens,longestInd)
  ENDIF ELSE BEGIN
     longestInd = 0
  ENDELSE

  good_i      = (good_i[strt_ii[longestInd]:stop_ii[longestInd]])

  ;;END METHOD MADNESS
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  T            = good_i[-1]-good_i[0]+1
  PRINT,'N good: ',STRCOMPRESS(T,/RE)

  TArr        = even_TS[good_i]
  T           = N_ELEMENTS(TArr)
  Ji_z_improv = Ji_z_improv[good_i]
  Je_z_improv = Je_z_improv[good_i]

  Bx = Bx[good_i]
  By = By[good_i]
  Bz = Bz[good_i]

  Bx = Bx - Bx[0]
  By = By - By[0]
  Bz = Bz - Bz[0]

  ;;Sorry, Jx and Jy
  mu_0    = DOUBLE(4.0D*!PI*1e-7)
  ;; jFactor = mu_0 * 1.e-6 ;;Put it in A/m^2
  jFactor = mu_0


  Jx = MAKE_ARRAY(T,VALUE=0.) * jFactor
  Jy = MAKE_ARRAY(T,VALUE=0.) * jFactor

  Jz = (Ji_z_improv + Je_z_improv) * jFactor
  ;; Jz = (2.*Ji_z_improv + (-1.)*Je_z_improv) * jFactor ;-1 for electrons because ... you know ...


  ;; Jz = ( (-1.)*Je_z_improv) * jFactor ;-1 for electrons because ... you know ...

  ;; Jz = SMOOTH(Jz,3)

END

PRO DEAL_WITH_BADNESS,datSerie,improvSerie

  improvSerie  = datSerie

  bad_i        = WHERE(~FINITE(datSerie),nBad,COMPLEMENT=good_i,NCOMPLEMENT=nGood)
  nDat         = nGood+nBad

  IF nBad GT 0 THEN BEGIN
     GET_STREAKS,bad_i,start_i=strtB_ii,stop_i=stopB_ii
     GET_STREAKS,good_i,start_i=strtG_ii,stop_i=stopG_ii
     strtB_i  = bad_i[strtB_ii]
     stopB_i  = bad_i[stopB_ii]

     strtG_i  = bad_i[strtG_ii]
     stopG_i  = bad_i[stopG_ii]

     FOR k=0,N_ELEMENTS(strtB_i)-1 DO BEGIN

        ;;If bad streak strts at beginning of TS or ends and end of TS, no interp
        strtTmp = strtB_i[k]
        stopTmp = stopB_i[k]

        IF strtTmp EQ 0 OR stopTmp EQ nDat-1 THEN BEGIN
           PRINT,'Skipping ...'
           CONTINUE
        ENDIF

        ;;Just zero this chunk
        improvSerie[strtTmp:stopTmp] = 0.
     ENDFOR

  ENDIF

END
;********************************************************************
;MAIN PROGRAM
;********************************************************************
; written by P. M. Bellan, April 2016, MS 128-95 Caltech, Pasadena CA 91125
; this code is free to use, but please acknowledge P. M. Bellan if used

PRO SINGLE_SPACECRAFT_K_MEASUREMENT_FAST, $
   PARSE_SAVEFILE=parse_saveFile, $
   SAVFILE_T1=s_t1, $
   SAVFILE_T2=s_t2, $
   USE_TIMEBAR_TIME=use_timeBar_time, $
   EXAMPLE_MODE=example_mode, $
   PLOT_KPERP_MAGNITUDE_FOR_KZ=plot_kperp_magnitude_for_kz, $
   PLOT_KX_VS_KY_FOR_KZ=plot_kx_vs_ky_for_kz, $
   PLOT_SMOOTHED_K_COMPONENTS=plot_smoothed_ks, $
   PLOT_ABS_SMOOTHED_K_COMPONENTS=plot_abs_smoothed_ks, $
   PLOT_POSFREQ=plot_posFreq, $
   FOLD_NEGFREQ_ONTO_POS=fold_negFreq, $
   SAVE_PS=save_ps, $
   BONUS_SUFF=bonus_suff, $
   DOUBLE_CALC=double_calc, $
   HANNING=hanning, $
   USE_J_TIME_SERIES=use_J_time_series, $
   SMOOTH_J_DAT_TO_B=smooth_J_dat, $
   OVERPLOT_DOUBLY_SMOOTHED=overplot_doubly_smoothed, $
   USE_DB_FAC=use_dB_fac, $
   STREAKNUM=streakNum, $
   PUBLICATION_SETTINGS=pubSettings

   
  ;;select output mode,
  ;;set output_mode =0 for screen, 1 for printer, 2 for postscript file

  COMPILE_OPT idl2

  ;;N points to smooth k components with
  smInd   = 8
  dbSmInd = 15
  edge_truncate = 1
  edge_mirror   = 0
  edge_wrap     = 0

  CASE 1 OF
     KEYWORD_SET(save_ps): BEGIN
        output_mode = 2
     END
     ELSE: BEGIN
        output_mode = 0
     END
  ENDCASE
  ;; output_mode = 0

  ;;select one of three examples
  ;; First example : A few discrete k components
  ;; Second example: Continuum of k components following some complicated function of omega
  ;; Third example : Same as second example, except noise added

  ;; example   = 3

  suff = 'TEST'
  IF KEYWORD_SET(parse_saveFile) THEN BEGIN
     CHUNK_SAVE_FILE,T,TArr,Bx,By,Bz,Jx,Jy,Jz,unitFactor,sPeriod,saveVar, $
                     SAVFILE_T1=s_t1, $
                     SAVFILE_T2=s_t2, $
                     USE_TIMEBAR_TIME=use_timeBar_time, $
                     USE_J_TIME_SERIES=use_J_time_series, $
                     SMOOTH_J_DAT_TO_B=smooth_J_dat, $
                     STREAKNUM=streakNum, $
                     OUT_STREAKNUM=streakInd, $
                     USE_DB_FAC=use_dB_fac

     sRate = 1./(TArr[1:-1]-TArr[0:-2])

     suff = 'Chaston_et_al_2006--ionos_erosion--Bellan_method'+'--'+saveVar

  ENDIF ELSE BEGIN
     unitFactor = 1 ;Don't adjust k in this case
  ENDELSE

  IF N_ELEMENTS(streakInd) GT 0 THEN BEGIN
     suff += '--streak_' + STRCOMPRESS(streakInd,/REMOVE_ALL)
  ENDIF

  ;;Aligned to mag or j time series?
  IF KEYWORD_SET(use_J_time_series) THEN BEGIN
     suff += '--jz_ts'
  ENDIF ELSE BEGIN
     suff += '--mag_ts'
  ENDELSE

  IF KEYWORD_SET(smooth_J_dat) THEN BEGIN
     suff += '--sm_je'
  ENDIF

  IF KEYWORD_SET(plot_kperp_magnitude_for_kz) THEN BEGIN
     suff += '--kPerp'
  ENDIF

  IF KEYWORD_SET(hanning    ) THEN BEGIN
     PRINT,"You shouldn't apply a window to correlation functions ..."
     WAIT,2
     suff += '--hanning'
  ENDIF
  IF KEYWORD_SET(double_calc) THEN suff += '--double_arithmetic'
  IF KEYWORD_SET(bonus_suff ) THEN suff += bonus_suff

  example = 1

  ;;setup file management, filenames for selected output mode
  OUTPUT_SETUP,output_mode,plotDir,suff,saveDir

  ;;setup graphic layout
  columns   = 3
  rows      = 3 + KEYWORD_SET(example_mode)
  !P.MULTI  = [0, columns, rows, 0, 0]


  ;;setup frequency, and wavevector arrays if in example mode

  ;; construct array of frequencies for time domain,
  ;; zero frequency will not be used

  ;; cs        = 1.8
  IF KEYWORD_SET(pubSettings) THEN BEGIN
     cs = 1.8
  ENDIF

  IF KEYWORD_SET(example_mode) THEN BEGIN
     SETUP_EXAMPLE,T,TArr,Bx,By,Bz,Jx,Jy,Jz,unitFactor,sPeriod,saveVar, $
                   EXAMPLE=example
  ENDIF

  ;;read data from file
  IF ~KEYWORD_SET(parse_saveFile) THEN BEGIN
     filename = saveDir+B_J_file

     OPENR,1,filename
     PRINT, 'opening ', filename
     READF,1,example,T_read,FORMAT='(I5,1x,I5,1x)'
     Bx = FLTARR(T_read) & By = Bx & Bz = Bx &Jx = Bx &Jy =Bx & Jz = Bx
     readf,1,Bx
     readf,1,By
     readf,1,Bz
     readf,1,Jx
     readf,1,Jy
     readf,1,Jz
     CLOSE,1
     PRINT,'closing ',filename
  ENDIF

  JxBtotal     = 0
  norm         = 0              ; check that avg J x B =0
  FOR TT=0,T-1 DO BEGIN         ; integrate over time
     Jvectemp  = [Jx[TT], Jy[TT], Jz[TT]]
     Bvectemp  = [Bx[TT], By[TT], Bz[TT]]
     JxBtotal  = JxBtotal+CROSSP(Jvectemp,Bvectemp)

     ;; norm for denominator
     norm      = norm + SQRT(Jx[TT]^2+ Jy[TT]^2+ Jz[TT]^2)*SQRT(Bx[TT]^2+ By[TT]^2+ Bz[TT]^2)
  ENDFOR

  PRINT,'norm = ',norm
  avgJxBtotal = JxBtotal/T
  PRINT, 'avgJxBtotal/norm = ',avgJxBtotal/norm ; small (supposed to be zero)

  ;; auto-correlation of B, i.e, <B_VEC(t) dot B_VEC(t+tau)>
  BBautocorr = CORRELATION(Bx,Bx,T) + CORRELATION(By,By,T) + CORRELATION(Bz,Bz,T)

  ;;components of cross correlation <J_VEC(t) cross B_vec (t+tau)>
  JxB_xcomponent_corr  = CORRELATION(Jy,Bz,T) - CORRELATION(Jz,By,T)

  JxB_ycomponent_corr  = CORRELATION(Jz,Bx,T) - CORRELATION(Jx,Bz,T)
  JxB_zcomponent_corr  = CORRELATION(Jx,By,T) - CORRELATION(Jy,Bx,T)

  ;; Fourier transform of <B_vec(t) dot B_vec(t+tau)>
  CASE 1 OF
     KEYWORD_SET(hanning): BEGIN
        hWindow              = HANNING(N_ELEMENTS(BBautocorr),/DOUBLE)
        BB_FFT               = FFT(BBautocorr*hWindow)
     END
     ELSE: BEGIN
        BB_FFT               = FFT(BBautocorr)
     END
  ENDCASE

  ;; Fourier transform of components of <B_vec(t) dot B_vec(t+tau)>
  CASE 1 OF
     KEYWORD_SET(hanning): BEGIN
        hWindow              = HANNING(N_ELEMENTS(JxB_xcomponent_corr),/DOUBLE)
        JxB_xcomponent_FFT   = FFT(JxB_xcomponent_corr*hWindow)
        JxB_ycomponent_FFT   = FFT(JxB_ycomponent_corr*hWindow)
        JxB_zcomponent_FFT   = FFT(JxB_zcomponent_corr*hWindow)
     END
     ELSE: BEGIN
        JxB_xcomponent_FFT   = FFT(JxB_xcomponent_corr)
        JxB_ycomponent_FFT   = FFT(JxB_ycomponent_corr)
        JxB_zcomponent_FFT   = FFT(JxB_zcomponent_corr)
     END
  ENDCASE

  ;; calculate k components, put 0.001 in denom to avoid dividing zero by zero
  ;; ikx                  = -JxB_xcomponent_FFT/(BB_FFT+.001)
  ;; iky                  = -JxB_ycomponent_FFT/(BB_FFT+.001)
  ;; ikz                  = -JxB_zcomponent_FFT/(BB_FFT+.001)

  IF KEYWORD_SET(double_calc) THEN nonZ = DOUBLE(1.e-10) ELSE nonZ = .001
  ikx                  = -JxB_xcomponent_FFT/(BB_FFT+nonZ)
  iky                  = -JxB_ycomponent_FFT/(BB_FFT+nonZ)
  ikz                  = -JxB_zcomponent_FFT/(BB_FFT+nonZ)

  ;;extract imaginary part
  kx                   = IMAGINARY(ikx)
  ky                   = IMAGINARY(iky)
  kz                   = IMAGINARY(ikz)

  ;;Multiply by unit factor to get k in m^-1, if it exists
  IF unitFactor NE 1 THEN BEGIN
     kx               *= unitFactor
     ky               *= unitFactor
     kz               *= unitFactor
  ENDIF

  kP = SQRT(kx^2+ky^2)

  IF KEYWORD_SET(plot_kperp_magnitude_for_kz) THEN BEGIN
     kz = kP
  ENDIF

  ;;Frequencies
  ;; T is an integer giving the number of elements in a particular dimension
  ;; sPeriod is a floating-point number giving the sampling interval
  X = FINDGEN((T - 1)/2) + 1
  is_T_even = (T MOD 2) EQ 0
  IF (is_T_even) THEN BEGIN
     freq = [0.0, X, T/2, -T/2 + X]/(T*sPeriod)
  ENDIF ELSE BEGIN
     freq = [0.0, X, -(T/2 + 1) + X]/(T*sPeriod)
  ENDELSE


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

  IF KEYWORD_SET(parse_saveFile) THEN BEGIN
     kx_ysize = MAX(ABS(kx))
     ky_ysize = MAX(ABS(ky))
     kz_ysize = MAX(ABS(kz))
  ENDIF

  bro = 1
  IF KEYWORD_SET(bro) THEN BEGIN
     tmp = SORT(freq)
     freq = freq[tmp]
     kx   = kx[tmp]
     ky   = ky[tmp]
     kz   = kz[tmp]
     kP = kP[tmp]
  ENDIF

  muLetter = '!4' + String('154'O) + '!X'
  PLOT,TArr-TArr[0],Bx, $
       XTITLE='t', $
       YTITLE='!8B!Dx!N', $
       CHARSIZE=cs
  PLOT,TArr-TArr[0],By, $
       XTITLE='t since' + TIME_TO_STR(TArr[0]) + '(s)', $
       YTITLE='!8B!Dy!N (nT)', $
       CHARSIZE=cs
  PLOT,TArr-TArr[0],Bz, $
       XTITLE='t', $
       YTITLE='!8B!Dz!N', $
       CHARSIZE=cs

  PLOT,TArr-TArr[0],Jx, $
       XTITLE='t', $
       YTITLE='!4l!3!D0 !N!8J!Dx!N', $
       CHARSIZE=cs
  PLOT,TArr-TArr[0],Jy, $
       XTITLE='t since' + TIME_TO_STR(TArr[0]) + '(s)', $
       YTITLE='!4l!3!D0 !N!8J!Dy!N (!4l!N!8T/m)', $
       CHARSIZE=cs
  PLOT,TArr-TArr[0],Jz, $
       XTITLE='t', $
       YTITLE='!4l!3!D0 !N!8J!Dz!N', $
       CHARSIZE=cs

  CASE 1 OF
     KEYWORD_SET(plot_posFreq): BEGIN
        inds = WHERE(freq GT 0.0 AND freq LE 4.0)
     END
     KEYWORD_SET(fold_negFreq): BEGIN
        indNeg = [0:(where(freq EQ 0.00)-1)] 
        ;; indPos = ([where(freq EQ 0.00,/NULL):N_ELEMENTS(freq)-1])[0:(N_ELEMENTS(indNeg)-1)]
        indPos = ([(where(freq GT 0.00))[0]:(N_ELEMENTS(freq)-1)])[0:(N_ELEMENTS(indNeg)-1)]

        this = plot(kx[indpos]) 
        this = plot((-1.)*REVERSE(kx[indneg]),/OVERPLOT,COLOR='RED') 

        divFactor = 2.
        kx   = (kx[indPos]-REVERSE(kx[indNeg])) / divFactor 

        this = plot(kx,/OVERPLOT,COLOR='BLUE')

        ky   = (ky[indPos]-REVERSE(ky[indNeg])) / divFactor 
        kz   = (kz[indPos]-REVERSE(kz[indNeg])) / divFactor
        freq = freq[indPos]
        ;; inds = WHERE(freq GT 0.0)

        IF freq[0] EQ 0.00 THEN BEGIN
           freq = freq[1:-1]
           kx = kx[1:-1]
           ky = ky[1:-1]
           kz = kz[1:-1]
        ENDIF

        inds = INDGEN(N_ELEMENTS(freq))
     END
     ELSE: BEGIN
        ;; inds = INDGEN(N_ELEMENTS(freq))
        inds = WHERE(ABS(freq) LE 4.0)
     END
  ENDCASE

  PLOT,freq[inds],kx[inds], $
       YTITLE='k!Dx!N (m!U-1!N)', $
       XTITLE='', $
       CHARSIZE=cs, $
       YRANGE=[-kx_ysize,kx_ysize]
  IF example EQ 2 THEN OPLOT,kx+0.1,LINESTYLE=2 ;dashed line
  IF KEYWORD_SET(plot_smoothed_ks) OR KEYWORD_SET(plot_abs_smoothed_ks) THEN BEGIN
     OPLOT,freq[inds],SMOOTH((KEYWORD_SET(plot_abs_smoothed_ks) ? ABS(kx[inds]) : kx[inds]),smInd, $
                              EDGE_TRUNCATE=edge_truncate,EDGE_MIRROR=edge_mirror,EDGE_WRAP=edge_wrap),COLOR=250
  ENDIF

  IF KEYWORD_SET(overplot_doubly_smoothed) THEN BEGIN
     OPLOT,freq[inds],SMOOTH((KEYWORD_SET(plot_abs_smoothed_ks) ? ABS(kx[inds]) : kx[inds]),dbSmInd, $
                             EDGE_TRUNCATE=edge_truncate,EDGE_MIRROR=edge_mirror,EDGE_WRAP=edge_wrap),COLOR=90
  ENDIF

  PLOT,freq[inds],ky[inds], $
       YTITLE='k!Dy!N (m!U-1!N)', $
       XTITLE='Frequency (Hz)', $
       CHARSIZE=cs, $
       YRANGE=[-ky_ysize,ky_ysize]
  IF example EQ 2 THEN OPLOT,ky+0.1,LINESTYLE=2
  IF KEYWORD_SET(plot_smoothed_ks) OR KEYWORD_SET(plot_abs_smoothed_ks) THEN BEGIN
     OPLOT,freq[inds],SMOOTH((KEYWORD_SET(plot_abs_smoothed_ks) ? ABS(ky[inds]) : ky[inds]),smInd, $
                             EDGE_TRUNCATE=edge_truncate,EDGE_MIRROR=edge_mirror,EDGE_WRAP=edge_wrap),COLOR=250
  ENDIF

  IF KEYWORD_SET(overplot_doubly_smoothed) THEN BEGIN
     OPLOT,freq[inds],SMOOTH((KEYWORD_SET(plot_abs_smoothed_ks) ? ABS(ky[inds]) : ky[inds]),dbSmInd, $
                             EDGE_TRUNCATE=edge_truncate,EDGE_MIRROR=edge_mirror,EDGE_WRAP=edge_wrap),COLOR=90
  ENDIF

  IF KEYWORD_SET(plot_kx_vs_ky_for_kz) THEN BEGIN
     
     ;; smoothKx = SMOOTH(kx[inds],smInd)
     ;; smoothKy = SMOOTH(ky[inds],smInd)

     ;; smoothKx = kx[inds]
     ;; smoothKy = ky[inds]

     smoothKx = kx
     smoothKy = ky

     ;;Now kx vs ky
     bound = MAX(ABS([smoothKx,smoothKy]))
     PLOT,smoothKx,smoothKy, $
          XTITLE='k!Dx!N  (m!U-1!N)', $
          YTITLE='k!Dy!N  (m!U-1!N)', $
          PSYM=2, $
          YRANGE=[-bound,bound], $
          XRANGE=[-bound,bound], $
          CHARSIZE=cs

  ENDIF ELSE BEGIN

     PLOT,freq[inds],kz[inds], $
          YTITLE='k!Dz!N (m!U-1!N)', $
          XTITLE='', $
          CHARSIZE=cs, $
          YRANGE=[-kz_ysize,kz_ysize]
     IF example EQ 2 THEN OPLOT,kz+0.1,LINESTYLE=2
     IF KEYWORD_SET(plot_smoothed_ks) THEN BEGIN
        OPLOT,freq[inds],SMOOTH(kz[inds],smInd,EDGE_TRUNCATE=edge_truncate,EDGE_MIRROR=edge_mirror,EDGE_WRAP=edge_wrap),COLOR=250
     ENDIF

     IF KEYWORD_SET(overplot_doubly_smoothed) THEN BEGIN
        OPLOT,freq[inds],SMOOTH(kz[inds],dbSmInd,EDGE_TRUNCATE=edge_truncate,EDGE_MIRROR=edge_mirror,EDGE_WRAP=edge_wrap),COLOR=90
     ENDIF

  ENDELSE

  IF KEYWORD_SET(save_ps) THEN BEGIN
     CONCLUDE_OUTPUT,output_mode
  ENDIF

  ;;Show Hanning windowed?
  IF KEYWORD_SET(hanning) THEN BEGIN

     IF ~KEYWORD_SET(save_ps) THEN BEGIN
        WINDOW,3,XSIZE=700,YSIZE=800
     ENDIF

     ;; !P.MULTI  = [0, columns, rows, 0, 0]
     !P.MULTI[0]  = 0

     PLOT,freq[inds],kx[inds], $
          YTITLE='x component', $
          XTITLE='FFT argument', $
          CHARSIZE=cs, $
          YRANGE=[-kx_ysize,kx_ysize]
     IF example EQ 2 THEN OPLOT,kx+0.1,LINESTYLE=2 ;dashed line

     PLOT,freq[inds],ky[inds], $
          YTITLE='y component', $
          XTITLE='FFT argument', $
          CHARSIZE=cs, $
          YRANGE=[-ky_ysize,ky_ysize]
     IF example EQ 2 THEN OPLOT,ky+0.1,LINESTYLE=2

     PLOT,freq[inds],kz[inds], $
          YTITLE='z component', $
          XTITLE='FFT argument', $
          CHARSIZE=cs, $
          YRANGE=[-kz_ysize,kz_ysize]
     IF example EQ 2 THEN OPLOT,kz+0.1,LINESTYLE=2



  ENDIF

  ;;Old or new style?
  IF KEYWORD_SET(oo_plots) THEN BEGIN
     window = WINDOW(DIMENSIONS=[700,800])

     bro = PLOT(freq[inds],kP[inds], $
                ;; XTITLE='Frequency (Hz)', $
                YTITLE='|k!Dperp!N (m$^{-1}$)', $
                CURRENT=window, $
                LAYOUT=[1,2,1])

     ;; bro.axes

     smkP = SMOOTH(kP[inds],smInd,EDGE_TRUNCATE=edge_truncate,EDGE_MIRROR=edge_mirror,EDGE_WRAP=edge_wrap)
     bro = PLOT(freq[inds],smkP, $
                ;; XTITLE='Frequency (Hz)', $
                YTITLE='|k!Dperp!N (m$^{-1}$)', $
                COLOR='RED', $
                /OVERPLOT, $
                CURRENT=window)

     ;;Kperp angle plot
     kPAngle = ATAN(ky,kx)*!RADEG
     ;; kPAngle = ATAN(SMOOTH(ky,smInd),SMOOTH(kx,smInd))*!RADEG

     bro = PLOT(freq[inds],kPAngle[inds], $
                XTITLE='Frequency (Hz)', $
                ;; XTITLE='T since' + TIME_TO_STR(TArr[0]) + '(s)', $
                YTITLE='|$\theta$(k!Dperp!N)', $
                CURRENT=window, $
                LAYOUT=[1,2,2])

     smkPAngle = SMOOTH(kPAngle[inds],smInd,EDGE_TRUNCATE=edge_truncate,EDGE_MIRROR=edge_mirror,EDGE_WRAP=edge_wrap)
     bro = PLOT(freq[inds],smkPAngle, $
                XTITLE='Frequency (Hz)', $
                YTITLE='|$\theta$(k!Dperp!N)', $
                COLOR='RED', $
                /OVERPLOT)
  ENDIF ELSE BEGIN

     IF KEYWORD_SET(save_ps) THEN BEGIN
        OUTPUT_SETUP,output_mode,plotDir,suff+'--page2',saveDir
     ENDIF

     columns = 2
     rows    = 1
     !P.MULTI  = [0, columns, rows, 0, 0]
     ;; !P.MULTI[0]  = 0

     IF ~KEYWORD_SET(save_ps) THEN BEGIN
        IF columns EQ 1 THEN BEGIN
           WINDOW,2,XSIZE=700,YSIZE=800
        ENDIF ELSE BEGIN
           WINDOW,2,XSIZE=1200,YSIZE=600
        ENDELSE
     ENDIF

      PLOT,freq[inds],kP[inds], $
           XTITLE='Frequency (Hz)', $
           YRANGE=[1e-5,1e-2], $
           XLOG=1, $
           YLOG=1, $
           YTITLE='!8k!Dperp!N (m!U-1!N)', $
           ;; XTICKFORMAT="(A1)", $
           CHARSIZE=cs

      smkP = SMOOTH(kP[inds],smInd,EDGE_TRUNCATE=edge_truncate,EDGE_MIRROR=edge_mirror,EDGE_WRAP=edge_wrap)
      OPLOT,freq[inds],smkP, $
            COLOR=250

      IF KEYWORD_SET(overplot_doubly_smoothed) THEN BEGIN
         dbSmkP = SMOOTH(kP[inds],dbSmInd,EDGE_TRUNCATE=edge_truncate,EDGE_MIRROR=edge_mirror,EDGE_WRAP=edge_wrap)
         OPLOT,freq[inds],dbSmkP, $
               COLOR=90
      ENDIF

      overplot_k_turbulence = 1
      IF KEYWORD_SET(overplot_k_turbulence) THEN BEGIN
         ;; kDoppl = (freq[inds]^(1.7)/10000)
         kxTemp = SMOOTH((KEYWORD_SET(plot_abs_smoothed_ks) ? $
                          ABS(kx) : kx), $
                         smInd, $
                         EDGE_TRUNCATE=edge_truncate, $
                         EDGE_MIRROR=edge_mirror, $
                         EDGE_WRAP=edge_wrap)
         wDoppl = kxTemp*(7000)
         OPLOT,wDoppl[inds],kxTemp[inds], $
               COLOR=70
      ENDIF

     ;;Kperp angle plot
     kPAngle = ATAN(ky,kx)*!RADEG


     PLOT,freq[inds],(kPAngle[inds] + 360) MOD 360, $
          XTITLE='Frequency (Hz)', $
          YTITLE='!8|' +CGGREEK('theta',PS=save_ps) + '!Dk!Dperp!N)', $
          CHARSIZE=cs
          ;; YTITLE='|$\theta$(k!Dperp!N)'

     smkPAngle = SMOOTH((kPAngle[inds] + 360) MOD 360,smInd,EDGE_TRUNCATE=edge_truncate,EDGE_MIRROR=edge_mirror,EDGE_WRAP=edge_wrap)
     OPLOT,freq[inds],smkPAngle, $
           COLOR=250

     IF KEYWORD_SET(overplot_doubly_smoothed) THEN BEGIN
        dbSmkPAngle = SMOOTH(kPAngle[inds],dbSmInd,EDGE_TRUNCATE=edge_truncate,EDGE_MIRROR=edge_mirror,EDGE_WRAP=edge_wrap)
        OPLOT,freq[inds],dbSmkPAngle, $
              COLOR=90
     ENDIF

     IF KEYWORD_SET(save_ps) THEN BEGIN
        CONCLUDE_OUTPUT,output_mode
     ENDIF


  ENDELSE

  IF KEYWORD_SET(have_Efield) THEN BEGIN
     CHECK_E_OMEGA_B_THING,Bx,By,Bz,Ex,Ey,Ez,kx,ky,kz,freq,inds
  ENDIF
  ;; print to screen, printer, or file depending on mode=0,1,2
  ;; CONCLUDE_OUTPUT,output_mode
END

;; PRO CHECK_E_OMEGA_B_THING,Bx,By,Bz,Ex,Ey,Ez,kx,ky,kz,freq,inds

;;   ;; kxEtotal     = 0
;;   ;; norm         = 0              ; check that avg k x E =0
;;   ;; FOR TT=0,T-1 DO BEGIN         ; integrate over time
;;   ;;    kvectemp  = [kx[TT], ky[TT], kz[TT]]
;;   ;;    Evectemp  = [Ex[TT], Ey[TT], Ez[TT]]
;;   ;;    kxEtotal  = kxEtotal+CROSSP(kvectemp,Evectemp)

;;   ;;    ;; norm for denominator
;;   ;;    norm      = norm + SQRT(kx[TT]^2+ ky[TT]^2+ kz[TT]^2)*SQRT(Ex[TT]^2+ Ey[TT]^2+ Ez[TT]^2)
;;   ;; ENDFOR
  
;;   PRINT,FORMAT='(A0,T15,A0,T30,A0)',"","",""
;;   FOR k=0,N_ELEMENTS(freq)-1 DO BEGIN
     
;;   ENDFOR

;;   PRINT,'norm = ',norm
;;   avgkxEtotal = kxEtotal/T
;;   PRINT, 'avgkxEtotal/norm = ',avgkxEtotal/norm ; small (supposed to be zero)

;;   omegaB      = 

;; END