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

PRO CHUNK_SAVE_FILE,T,TArr,Bx,By,Bz,Jx,Jy,Jz,unitFactor,sPeriod,saveVar

  saveDir  = '/SPENCEdata/Research/Satellites/FAST/single_sc_wavevector/'
  saveFile = 'Chaston_et_al_2006--B_and_J.sav'

  PRINT,"Restoring " + saveFile + ' ...'
  RESTORE,saveDir+saveFile

  bUnitFactor = -9 ;'cause nT
  jUnitFactor = -6 ;'cause microA/m^2

  unitFactor  = (10.D)^(jUnitFactor)/(10.D)^(bUnitFactor)

  ;; saveVar  = 'dB_fac_v'
  saveVar  = 'dB_fac_V'

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

  TArr = dB.x
  T    = N_ELEMENTS(TArr)
  ;; bFactor = 1.e-9 ;Get 'em out of nT
  ;; bFactor = 1.e3
  bFactor = 1.
  Bx   = dB.y[*,0] * bFactor
  By   = dB.y[*,1] * bFactor
  Bz   = dB.y[*,2] * bFactor

  mag_sRate  = 1./(dB.x  [1:-1]-dB.x  [0:-2])
  eESA_sRate = 1./(Je_z.x[1:-1]-Je_z.x[0:-2])
  iESA_sRate = 1./(Ji_z.x[1:-1]-Ji_z.x[0:-2])

  this       = VALUE_CLOSEST2(je_z.x,db.x)
  PRINT,ABS(je_z.x[this]-db.x)

  this       = VALUE_CLOSEST2(ji_z.x,db.x)
  PRINT,ABS(ji_z.x[this]-db.x)


  frame_of_ref = 2
  CASE frame_of_ref OF
     1: BEGIN
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

        even_TS = even_TS[60:-60]

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

  good_i       = WHERE(FINITE(je_z_improv) AND FINITE(ji_z_improv),nGood)
  GET_STREAKS,good_i,START_i=strt_ii,STOP_I=stop_ii,OUT_STREAKLENS=streakLens
  long         = MAX(streakLens,longestInd)
  good_i       = (good_i[strt_ii[longestInd]:stop_ii[longestInd]])

  ;;END METHOD MADNESS
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  T            = good_i[-1]-good_i[0]+1
  PRINT,'N good: ',STRCOMPRESS(T,/RE)

  TArr        = Tarr[good_i]
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


  Jx = FLTARR(T) * jFactor
  Jy = FLTARR(T) * jFactor

  Jz = (2.*Ji_z_improv + (-1.)*Je_z_improv) * jFactor ;-1 for electrons because ... you know ...

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
   EXAMPLE_MODE=example_mode, $
   PLOT_KPERP_MAGNITUDE_FOR_KZ=plot_kperp_magnitude_for_kz, $
   PLOT_POSFREQ=plot_posFreq
  ;;select output mode,
  ;;set output_mode =0 for screen, 1 for printer, 2 for postscript file
  output_mode = 0
  ;;select one of three examples
  ;; First example : A few discrete k components
  ;; Second example: Continuum of k components following some complicated function of omega
  ;; Third example : Same as second example, except noise added

  ;; example   = 3

  suff = 'TEST'
  IF KEYWORD_SET(parse_saveFile) THEN BEGIN
     CHUNK_SAVE_FILE,T,TArr,Bx,By,Bz,Jx,Jy,Jz,unitFactor,sPeriod,saveVar

     sRate = 1./(TArr[1:-1]-TArr[0:-2])

     suff = 'Chaston_et_al_2006--ionos_erosion--Bellan_method'+'--'+saveVar

  ENDIF ELSE BEGIN
     unitFactor = 1 ;Don't adjust k in this case
  ENDELSE

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

  IF KEYWORD_SET(example_mode) THEN BEGIN

     T         = 200            ; number of elements in time domain vectors
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
     Avec_cos = FLTARR(T/2+1,3) ; define A vector cos component (vector potential)
     Avec_sin = FLTARR(T/2+1,3) ; define A vector sine component (vector potential)

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
     FOR TT=0,T-1 DO BEGIN      ;for each time
        FOR i=0,T/2 DO BEGIN    ;for each frequency component

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
     FOR TT=0,T-1 DO BEGIN      ; TT is time
        FOR i=0,T/2 DO BEGIN    ; sum up contributions from each frequency
           Bx[TT] = Bx[TT] +Bvec[i,0,TT]
           By[TT] = By[TT] +Bvec[i,1,TT]
           Bz[TT] = Bz[TT] +Bvec[i,2,TT]
           Jx[TT] = Jx[TT] +Jvec[i,0,TT]
           Jy[TT] = Jy[TT] +Jvec[i,1,TT]

           Jz[TT] = Jz[TT] +Jvec[i,2,TT]
        ENDFOR                  ;
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

  muLetter = '!4' + String('154'O) + '!X'
  PLOT,TArr-TArr[0],Bx, $
       XTITLE='t', $
       YTITLE='!8B!Dx!N (nT)', $
       CHARSIZE=cs
  PLOT,TArr-TArr[0],By, $
       XTITLE='t since' + TIME_TO_STR(TArr[0]) + '(s)', $
       YTITLE='!8B!Dy!N', $
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
  BB_FFT               = FFT(BBautocorr)

  ;; Fourier transform of components of <B_vec(t) dot B_vec(t+tau)>
  JxB_xcomponent_FFT   = FFT(JxB_xcomponent_corr)
  JxB_ycomponent_FFT   = FFT(JxB_ycomponent_corr)
  JxB_zcomponent_FFT   = FFT(JxB_zcomponent_corr)

  ;; calculate k components, put 0.001 in denom to avoid dividing zero by zero
  ikx                  = -JxB_xcomponent_FFT/(BB_FFT+.001)
  iky                  = -JxB_ycomponent_FFT/(BB_FFT+.001)
  ikz                  = -JxB_zcomponent_FFT/(BB_FFT+.001)

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

  CASE 1 OF
     KEYWORD_SET(plot_posFreq): BEGIN
        inds = WHERE(freq GT 0.0)
     END
     ELSE: BEGIN
        inds = INDGEN(freq)        
     END
  ENDCASE

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



  bro = PLOT(freq[inds],kP[inds], $
             XTITLE='Frequency (Hz)', $
             YTITLE='|k!Dperp!N (m$^{-1}$)')

  smInd = 10
  smkP = SMOOTH(kP[inds],smInd)
  bro = PLOT(freq[inds],smkP, $
             XTITLE='Frequency (Hz)', $
             YTITLE='|k!Dperp!N (m$^{-1}$)', $
             COLOR='RED', $
             /OVERPLOT)

  ;;Kperp angle plot
  kPAngle = ATAN(ky,kx)*!RADEG
  ;; kPAngle = ATAN(SMOOTH(ky,smInd),SMOOTH(kx,smInd))*!RADEG

  bro = PLOT(freq[inds],kPAngle[inds], $
             XTITLE='Frequency (Hz)', $
             ;; XTITLE='T since' + TIME_TO_STR(TArr[0]) + '(s)', $
             YTITLE='|$\theta$(k!Dperp!N)')

  smkPAngle = SMOOTH(kPAngle[inds],smInd)
  bro = PLOT(freq[inds],smkPAngle, $
             XTITLE='Frequency (Hz)', $
             YTITLE='|$\theta$(k!Dperp!N)', $
             COLOR='RED', $
             /OVERPLOT)

  ;; print to screen, printer, or file depending on mode=0,1,2
  CONCLUDE_OUTPUT,output_mode
END