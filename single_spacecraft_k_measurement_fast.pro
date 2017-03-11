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

  autocorr            = FLTARR(T)
  FOR tau=0,T-1 DO BEGIN
     FOR TT=0,T-1 DO BEGIN
        TT_plus_tau   = (TT+tau) MOD T ; use modulo arithmetic to wrap back to beginning
        autocorr[tau] = autocorr[tau]+f[TT]*g[TT_plus_tau]
     ENDFOR
  ENDFOR
  RETURN,autocorr/T

END

;********************************************************************
PRO ADD_NOISE,Bx,By,Bz,Jx,Jy,Jz,T
  ;;create and add noise IF example = 3

  COMPILE_OPT idl2,strictarrsubs

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

PRO SETUP_EXAMPLE,T,TArr,Bx,By,Bz,Jx,Jy,Jz,unitFactor,sPeriod,saveVar, $
                  SAVEDIR=saveDir, $
                  EXAMPLE=example

  COMPILE_OPT idl2,strictarrsubs

  T         = 200               ; number of elements in time domain vectors
  TArr      = INDGEN(T)
  sPeriod   = 1

  ;; k_vector : first argument is for frequency, second is for Cartesian coordinate
  kvec      = FLTARR(T/2+1,3)

  omega     = INDGEN(T/2+1)*2*!Pi/T

  ;;setup k vectors for either example 1 or 2
  ;;first example
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

  PLOT,kx,XTITLE='!4x!3T/(2!4p!3)',YTITLE=font + 'k!Dx!N',CHARSIZE=cs,YRANGE=[-ysize,ysize]
  PLOT,ky,XTITLE='!4x!3T/(2!4p!3)',YTITLE=font + 'k!Dy!N',CHARSIZE=cs,YRANGE=[-ysize,ysize]
  PLOT,kz,XTITLE='!4x!3T/(2!4p!3)',YTITLE=font + 'k!Dz!N',CHARSIZE=cs,YRANGE=[-ysize,ysize]

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

PRO CHECK_K_OMEGA_ODDNESS,freq,kx,ky,kz, $
                          KZ_IS_KPERP=kz_is_kPerp

  sort_i = SORT(freq)
  freqT  = freq[sort_i]
  kxT  = kx[sort_i]
  kyT  = ky[sort_i]
  kzT  = kz[sort_i]

  indNeg = WHERE(freqT LE 0.00,nNeg)
  indPos = WHERE(freqT GT 0.00,nPos)
  IF nNeg NE nPos THEN BEGIN
     IF nNeg GT nPos THEN BEGIN
        indNeg = indNeg[0:-2]
     ENDIF ELSE BEGIN
        indPos = indPos[0:-2]
     ENDELSE
  ENDIF
  
  PRINT,'avg freq oddness : ',MEAN(DOUBLE(freq[indPos])+DOUBLE(freq[indNeg]),/DOUBLE)
  PRINT,FORMAT='(A0,T20,G0.3,T30,G0.3,T40,G0.3,T50,A0)','avg k oddness    : ', $
        MEAN((DOUBLE(kxT))[indPos]+REVERSE((DOUBLE(kxT))[indNeg]),/DOUBLE), $
        MEAN((DOUBLE(kyT))[indPos]+REVERSE((DOUBLE(kyT))[indNeg]),/DOUBLE), $
        MEAN((DOUBLE(kzT))[indPos]+REVERSE((DOUBLE(kzT))[indNeg]),/DOUBLE), $
        (KEYWORD_SET(kz_is_kPerp) ? "(kz is kPerp, you know)" : '')

END

PRO CHECK_E_OMEGA_B_THING,Bx,By,Bz,Ex_sp,Ey_sp,Ez_sp,kx,ky,kz,freq ;,freq_sp,inds

  ;; omega_sp        = freq_sp * 2. * !PI
  omega_sp        = freq * 2. * !PI

  Bx_om           = FFT(Bx)
  By_om           = FFT(By)
  Bz_om           = FFT(Bz)

  Ex_om           = FFT(Ex_sp)
  Ey_om           = FFT(Ey_sp)
  Ez_om           = FFT(Ez_sp)

  omegaBtotal     = 0
  kxEtotal        = 0
  norm            = 0           ; check that avg (k x E) = 0
  diffs           = MAKE_ARRAY(3,N_ELEMENTS(freq),/DOUBLE)
  FOR k=0,N_ELEMENTS(freq)-1 DO BEGIN
     kvectemp     = [kx[k], ky[k], kz[k]]
     Evectemp     = [Ex_sp[k], Ey_sp[k], Ez_sp[k]]
     Bvectemp     = [Bx_om[k], By_om[k], Bz_om[k]]
     kxEtotal     = kxEtotal+CROSSP(kvectemp,Evectemp)
     omegaBtotal  = omegaBtotal + omega_sp[k]*Bvectemp
     diffs[*,k]   = CROSSP(kvectemp,Evectemp)-omega_sp[k]*Bvectemp

     ;; norm for denominator
     ;; norm         = norm + SQRT(kx[k]^2+ ky[k]^2+ kz[k]^2)*SQRT(Ex[k]^2+ Ey[k]^2+ Ez[k]^2)
  ENDFOR
  
  PRINT,FORMAT='(A0,T15,A0,T30,A0)',"","",""

  PRINT,'norm     : ',norm
  avgkxEtotal     = kxEtotal/T
  PRINT, 'avgkxEtotal/norm : ',avgkxEtotal/norm ; small (supposed to be zero)

END

PRO PREDICT_J,freq,kx,ky,kz,Bx,By,Bz,Jx,Jy,Jz,unitFactor

  IF N_ELEMENTS(unitFactor) EQ 0 THEN unitFactor = 1.D

  mu_0            = DOUBLE(4.0D*!PI*1e-7)
  iCmplx          = DCOMPLEX(0,1)

  Bx_om           = FFT(Bx)
  By_om           = FFT(By)
  Bz_om           = FFT(Bz)

  Jx_om           = FFT(Jx)
  Jy_om           = FFT(Jy)
  Jz_om           = FFT(Jz)

  
  nHere           = N_ELEMENTS(freq)
  
  omegaBtotal     = 0
  kxEtotal        = 0
  norm            = 0           ; check that avg (k x E) = 0
  diffs           = MAKE_ARRAY(3,nHere,/DCOMPLEX)
  JPred           = MAKE_ARRAY(3,nHere,/DCOMPLEX)
  errAngle        = MAKE_ARRAY(nHere,/DCOMPLEX)
  magErr          = MAKE_ARRAY(nHere,/DCOMPLEX)
  FOR k=0,nHere-1 DO BEGIN

     ;;Measured J
     Jvectemp     = [Jx_om[k], Jy_om[k], Jz_om[k]]

     ;;Measured B and derived k
     Bvectemp     = [Bx_om[k], By_om[k], Bz_om[k]]
     kvectemp     = [kx[k], ky[k], kz[k]] / unitFactor ;back into units native for J and B

     ;;Predicted J
     Jpred[*,k]   = iCmplx * CROSSP(kvectemp,Bvectemp) / mu_0

     diffs[*,k]   = Jpred[*,k]-Jvectemp
     tmpSum       = Jpred[*,k]+Jvectemp

     ;;Error angle
     errAngle[k]  = ACOS(DOT(Jpred[*,k],Jvectemp)/(SQRT(DOT(Jpred[*,k],CONJ(Jpred[*,k])))*SQRT(DOT(Jvectemp,CONJ(Jvectemp)))))

     ;;Magnitude error
     magErr[k]    = SQRT(DOT(diffs[*,k],CONJ(diffs[*,k])))/SQRT(DOT(tmpSum,CONJ(tmpSum)))

  ENDFOR
  
  PRINT,FORMAT='(A0,T15,A0,T30,A0)',"","",""

  ;; PRINT, 'avgkxEtotal/norm : ',avgkxEtotal/norm ; small (supposed to be zero)
  

END

PRO BELLAN_2016__BRO,T,Jx,Jy,Jz,Bx,By,Bz, $
                     freq,kx,ky,kz,kP, $
                     SPERIOD=sPeriod, $
                     UNITFACTOR=unitFactor, $
                     PLOT_KPERP_MAGNITUDE_FOR_KZ=plot_kperp_magnitude_for_kz, $
                     DOUBLE_CALC=double_calc, $
                     HANNING=hanning, $
                     EFIELD=EField, $
                     ODDNESS_CHECK=oddness_check, $
                     OUT_NORM=norm, $
                     OUT_AVGJXB=avgJxBtotal, $
                     DIAG=diag

  COMPILE_OPT idl2,strictarrsubs

  ;; diag         = 1

  JxBtotal     = 0
  norm         = 0              ; check that avg J x B =0
  FOR TT=0,T-1 DO BEGIN         ; integrate over time
     Jvectemp  = [Jx[TT], Jy[TT], Jz[TT]]
     Bvectemp  = [Bx[TT], By[TT], Bz[TT]]
     JxBtotal  = JxBtotal+CROSSP(Jvectemp,Bvectemp)

     ;; norm for denominator
     norm      = norm + SQRT(Jx[TT]^2+ Jy[TT]^2+ Jz[TT]^2)*SQRT(Bx[TT]^2+ By[TT]^2+ Bz[TT]^2)

     IF KEYWORD_SET(diag) THEN BEGIN
        PRINT,FORMAT='("avgJxBtotal/norm[",I0,"]: ",G12.5,T40,G12.5,T54,G10.5)',TT,JxBtotal/T/norm ; small (supposed to be zero)
     ENDIF

  ENDFOR

  PRINT,'norm             : ',norm
  avgJxBtotal = JxBtotal/T
  PRINT,'avgJxBtotal/norm : ',avgJxBtotal/norm ; small (supposed to be zero)

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
  ;; nonZ = COMPLEX(1.e-10,1.e-10)
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

  IF KEYWORD_SET(Efield) THEN BEGIN
     CHECK_E_OMEGA_B_THING,Bx,By,Bz,EField[0,*],EField[1,*],EField[2,*],kx,ky,kz,freq
  ENDIF

  IF KEYWORD_SET(oddness_check) THEN BEGIN
     CHECK_K_OMEGA_ODDNESS,freq,kx,ky,kz, $
                          KZ_IS_KPERP=plot_kperp_magnitude_for_kz
  ENDIF

  ;;Calculate error thing?
  predict_current = 1
  IF KEYWORD_SET(predict_current) THEN BEGIN
     PREDICT_J,freq,kx,ky,kz,Bx,By,Bz,Jx,Jy,Jz,unitFactor
  ENDIF


END

PRO DEAL_WITH_BADNESS,datSerie,improvSerie

  COMPILE_OPT idl2,strictarrsubs

  improvSerie  = datSerie

  bad_i        = WHERE(~FINITE(datSerie),nBad,COMPLEMENT=good_i,NCOMPLEMENT=nGood)
  nDat         = nGood+nBad

  IF nBad GT 0 THEN BEGIN

     GET_STREAKS,bad_i,START_I=strtB_ii,STOP_I=stopB_ii,/QUIET,/NO_PRINT_SUMMARY
     GET_STREAKS,good_i,START_I=strtG_ii,STOP_I=stopG_ii,/QUIET,/NO_PRINT_SUMMARY
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

        ;;Just zero this bit
        improvSerie[strtTmp:stopTmp] = 0.
     ENDFOR

  ENDIF

END

FUNCTION CHUNK_SAVE_FILE,T,TArr,Bx,By,Bz,Jx,Jy,Jz,unitFactor,sPeriod,saveVar, $
                         B_AND_J_FILE=saveFile, $
                         SAVEDIR=saveDir, $
                         USE_TIMEBAR_TIME__FROM_FILE=use_timeBar_time__from_file, $
                         CUSTOM_T1=custom_t1, $
                         CUSTOM_T2=custom_t2, $
                         SHIFT_NPTS=shift_nPts, $
                         USE_LOWRES_TIME_SERIES=use_lowRes_time_series, $
                         USE_J_TIME_SERIES=use_J_time_series, $
                         SMOOTH_J_DAT_TO_B=smooth_J_dat, $
                         PRESMOOTH_MAG=presmooth_mag, $
                         PREPLOT_CURRENTS_AND_STOP=prePlot_currents_and_stop, $
                         STREAKNUM=streakNum, $
                         OUT_STREAKNUM=longestInd, $
                         USE_ALL_STREAKS=use_all_streaks, $
                         USE_DB_FAC=use_dB_fac, $
                         EFIELD=EField, $
                         FFT__NEAREST_TWO_POWER=nearest_two_power, $
                         SRATES=sRates

  COMPILE_OPT idl2,strictarrsubs

  PRINT,"Restoring " + saveFile + ' ...'
  RESTORE,saveDir+saveFile

  IF KEYWORD_SET(use_timeBar_time__from_file) AND N_ELEMENTS(timeBar_times) GT 0 THEN BEGIN
     PRINT,"Using timeBar_times from file: ",timeBar_times
     CASE NDIMEN(timeBar_times) OF
        1: BEGIN
           CASE SIZE(timeBar_times,/TYPE) OF
              7: BEGIN
                 tBar_t1  = STR_TO_TIME(timeBar_times[0])
                 tBar_t2  = STR_TO_TIME(timeBar_times[1])
              END
              5: BEGIN
                 tBar_t1  = timeBar_times[0]
                 tBar_t2  = timeBar_times[1]
              END
              ELSE: STOP
           ENDCASE
        END
        2: BEGIN
           CASE SIZE(timeBar_times,/TYPE) OF
              7: BEGIN
                 tBar_t1  = STR_TO_TIME(REFORM(timeBar_times[0,*]))
                 tBar_t2  = STR_TO_TIME(REFORM(timeBar_times[1,*]))
              END
              5: BEGIN
                 tBar_t1  = REFORM(timeBar_times[0,*])
                 tBar_t2  = REFORM(timeBar_times[1,*])
              END
              ELSE: STOP
           ENDCASE
        END
     ENDCASE
  ENDIF

  ;;Just in case user is inclined to sparse provision of information
  IF N_ELEMENTS(dB) EQ 0 THEN BEGIN

     have_dB_fac   = N_ELEMENTS(dB_fac) NE 0
     have_dB_fac_v = N_ELEMENTS(dB_fac_v) NE 0

     CASE 1 OF
        (have_dB_fac AND have_dB_fac_v): BEGIN
           IF ~ARRAY_EQUAL(dB_fac.x,dB_fac_v.x) THEN BEGIN
              PRINT,"What?"
              STOP
           ENDIF

           dB = dB_fac

        END
        have_dB_fac: BEGIN

           dB = dB_fac

        END
        have_dB_fac_v: BEGIN

           dB = dB_fac_v
           
        END
        ELSE: BEGIN

           PRINT,"Got nothin!"
           STOP

        END
     ENDCASE

  ENDIF

  t1BKUP      = N_ELEMENTS(t1Zoom) GT 0 ? t1Zoom : dB.x[0]
  t2BKUP      = N_ELEMENTS(t2Zoom) GT 0 ? t2Zoom : dB.x[-1]
  
  ;;Now custom times
  IF KEYWORD_SET(custom_t1) THEN BEGIN
     PRINT,FORMAT='(A0,T25,": ",A0)','Custom T1',SIZE(custom_t1,/TYPE) EQ 7 ? custom_t1 : TIME_TO_STR(custom_t1,/MS)

     CASE NDIMEN([custom_t1]) OF
        1: BEGIN
           CASE SIZE(custom_t1,/TYPE) OF
              7: BEGIN
                 analysis_t1  = STR_TO_TIME(custom_t1)
              END
              5: BEGIN
                 analysis_t1  = custom_t1
              END
              ELSE: STOP
           ENDCASE
        END
        2: BEGIN
           CASE SIZE(custom_t1,/TYPE) OF
              7: BEGIN
                 analysis_t1  = STR_TO_TIME(REFORM(custom_t1[0,*]))
              END
              5: BEGIN
                 analysis_t1  = REFORM(custom_t1[0,*])
              END
              ELSE: STOP
           ENDCASE
        END
     ENDCASE

  ENDIF ELSE BEGIN
     analysis_t1        = KEYWORD_SET(tBar_t1) ? TEMPORARY(tBar_t1) : t1BKUP
  ENDELSE

  IF KEYWORD_SET(custom_t2) THEN BEGIN
     PRINT,FORMAT='(A0,T25,": ",A0)','Custom T2',SIZE(custom_t2,/TYPE) EQ 7 ? custom_t2 : TIME_TO_STR(custom_t2,/MS)

     CASE NDIMEN([custom_t2]) OF
        1: BEGIN
           CASE SIZE(custom_t2,/TYPE) OF
              7: BEGIN
                 analysis_t2  = STR_TO_TIME(custom_t2)
              END
              5: BEGIN
                 analysis_t2  = custom_t2
              END
              ELSE: STOP
           ENDCASE
        END
        2: BEGIN
           CASE SIZE(custom_t2,/TYPE) OF
              7: BEGIN
                 analysis_t2  = STR_TO_TIME(REFORM(custom_t2[0,*]))
              END
              5: BEGIN
                 analysis_t2  = REFORM(custom_t2[0,*])
              END
              ELSE: STOP
           ENDCASE
        END
     ENDCASE

  ENDIF ELSE BEGIN
     analysis_t2        = KEYWORD_SET(tBar_t2) ? TEMPORARY(tBar_t2) : t2BKUP
  ENDELSE

  ;; bFactor = 1.e-9 ;Get 'em out of nT
  ;; bFactor = 1.e3
  bFactor     = 1.D

  jFactor     = 1.D
  ;; jFactor = 1.e-6 ;;Put it in A/m^2
  ;; jFactor = mu_0

  EFactor     = 1.D

  bUnitFactor = (bFactor EQ 1.D) ? -9.D : -9. /ALOG10(bFactor) ;'cause nT
  jUnitFactor = (jFactor EQ 1.D) ? -6.D : -6 / ALOG10(jFactor) ;'cause microA/m^2
  EUnitFactor = (EFactor EQ 1.D) ? -3.D : -3 / ALOG10(EFactor) ;'cause mV/m

  mu_0        = DOUBLE(4.0D*!PI*1e-7)
  unitFactor  = (10.D)^(jUnitFactor)/(10.D)^(bUnitFactor)*mu_0

  ;; mu_0      = DOUBLE(4.0D*!PI*1e-7)
  ;; jFactor  *= mu_0

  ;;Now decide on magField
  saveVar  = 'fac_V'
  IF KEYWORD_SET(use_dB_fac) THEN saveVar = 'fac'
  ;; saveVar  = 'dB_fac_V'

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;Align time series
  CASE STRUPCASE(saveVar) OF
     'FAC_V': BEGIN

        ;; From UCLA_MAG_DESPIN:
        ;; "Field-aligned velocity-based coordinates defined as:    "
        ;; x (ind 0)-along track ((BxV)xB),
        ;; y (ind 1)-cross track (BxV),
        ;; z (ind 2)-along B" (I added "ind" marks)

        dB = dB_fac_v

        ;;EField too?
        ;;Presumed to be in FAC_V coordinates
        have_EField =  N_ELEMENTS(EField) GT 0

     END
     'FAC': BEGIN

        ;; From UCLA_MAG_DESPIN:
        ;; Field-aligned coordinates defined as:
        ;; z-along B, y-east (BxR), x-nominally out

        dB = dB_fac

        have_EField = 0B ;Doesn't matter if you do; you're in the wrong coordinate system
        EField      = !NULL ;Don't want to mess anyone else up

     END
  ENDCASE

  ;;Now make 'em Bellan-useable
  Bx      = dB.y[*,0] * bFactor
  By      = dB.y[*,1] * bFactor
  Bz      = dB.y[*,2] * bFactor
  Bt      = dB.x

  IF have_EField THEN BEGIN
     Ex   = EField.y[*,0] * EFactor
     Ey   = EField.y[*,1] * EFactor
     Ez   = EField.y[*,2] * EFactor
  ENDIF

  IF KEYWORD_SET(presmooth_mag) THEN BEGIN

     ;;Do one smooth
     ;; smInd  = 20
     ;; redFac = 20
     ;; Bx    = (SMOOTH(Bx,smInd))[0:-1:redFac]
     ;; By    = (SMOOTH(By,smInd))[0:-1:redFac]
     ;; Bz    = (SMOOTH(Bz,smInd))[0:-1:redFac]
     ;; Bt    = Bt[0:-1:redFac]


     ;;... Or do two smooths in a row
     smInd  = 21
     redFac = 2
     Bx    = (SMOOTH(Bx,smInd))[0:-1:redFac]
     By    = (SMOOTH(By,smInd))[0:-1:redFac]
     ;; Bx    = Bx[0:-1:redFac]
     ;; By    = By[0:-1:redFac]
     Bz    = (SMOOTH(Bz,smInd))[0:-1:redFac]
     Bt    = Bt[0:-1:redFac]

     smInd  = 11
     redFac = 2
     Bx    = (SMOOTH(Bx,smInd))[0:-1:redFac]
     By    = (SMOOTH(By,smInd))[0:-1:redFac]
     ;; Bx    = Bx[0:-1:redFac]
     ;; By    = By[0:-1:redFac]
     Bz    = (SMOOTH(Bz,smInd))[0:-1:redFac]
     Bt    = Bt[0:-1:redFac]

     ;; smInd  = 10
     ;; redFac = 2
     ;; Bx    = (SMOOTH(Bx,smInd))[0:-1:redFac]
     ;; By    = (SMOOTH(By,smInd))[0:-1:redFac]
     ;; ;; Bx    = Bx[0:-1:redFac]
     ;; ;; By    = By[0:-1:redFac]
     ;; Bz    = (SMOOTH(Bz,smInd))[0:-1:redFac]
     ;; Bt    = Bt[0:-1:redFac]

     ;; plotme = magz
     ;; bro    = (SMOOTH(plotme.y,smInd))[0:-1:smInd]
     ;; dough  = plotme.x[0:-1:smInd]
     ;; this2  = PLOT(plotme.x-plotme.x[0],plotme.y-plotme.y[0])
     ;; this   = PLOT(dough-dough[0],bro-bro[0],COLOR='RED',/OVERPLOT)
     ;; thot   = PLOT(plotme.x-plotme.x[0],SMOOTH(plotme.y,smInd)-plotme.y,COLOR='BLUE',/OVERPLOT)

     ;; PRINT,1./(dough[1]-dough[0])
     PRINT,"Final B sample frequency: ",1./(Bt[1]-Bt[0])
  ENDIF

  IF KEYWORD_SET(prePlot_currents_and_stop) THEN BEGIN

     ;;scope out mag current
     ;;Get orbit stuff
     GET_FA_ORBIT,Bt,/TIME_ARRAY ;,/all

     ;;Get speed and position for calculation of mag stuff
     GET_DATA,'fa_vel',DATA=vel
     speed                   = SQRT(vel.y[*,0]^2+vel.y[*,1]^2+vel.y[*,2]^2)*1000.0

     old_pos                 = 0.
     position                = MAKE_ARRAY(N_ELEMENTS(Bt),/DOUBLE)
     speed_mag_point         = MAKE_ARRAY(N_ELEMENTS(Bt),/DOUBLE)
     FOR j=0L,N_ELEMENTS(Bt)-2 DO BEGIN
        speed_point_ind      = MIN(ABS(vel.x-Bt[j]),ind)

        speed_mag_point[j]   = speed[ind]
        samplingperiod       = Bt[j+1] - Bt[j]

        position[j]          = old_pos + speed_mag_point[j]*samplingperiod
        old_pos              = position[j]
     ENDFOR

     ;; deltaBY                 = DERIV(position,SMOOTH(By,3))
     deltaBY                 = DERIV(position,By)
     ;; deltaBY                 = DERIV(position,By)
     ;; deltaBY                 = DERIV(position,SMOOTH(By,5))
     ;; jtemp                = ABS(1.0e-3*(deltaBx)/1.26e-6)
     ;; jtemp                = 1.0e-3*(deltaBx)/1.26e-6
     ;;in number flux units
     jtemp                   = 1.0e-3*(deltaBY)/1.26e-6
     muLetter = '!4l!X'

     t0 = analysis_t1
     this = PLOT(Bt-t0,jtemp,NAME='Magnetometer', $
                 XTITLE='Time since ' + TIME_TO_STR(analysis_t1,/MSEC), $
                 YTITLe='Current density ($\mu$A/m$^2$', $
                 XRANGE=[0,analysis_t2-analysis_t1], $
                 ;; YRANGE=MINMAX(jtemp))
                 YRANGE=[MIN(jtemp) < MIN(je_z.y) < MIN(ji_z.y), $
                         MAX(jtemp) > MAX(je_z.y) > MAX(ji_z.y)])

     this = PLOT(Je_z.x-t0,Je_z.y, $
                 ;; this = PLOT(Je_z.x-t0,Je_z.y, $
                 NAME='e!U-!N', $
                 COLOR='Red', $
                 /OVERPLOT)

     this = PLOT(Ji_z.x-t0,Ji_z.y, $
                 ;; this = PLOT(Ji_z.x-t0,SMOOTH(2*Ji_z.y+je_z.y,5)*0.6-25, $
                 NAME='i!U+!N', $
                 COLOR='Green', $
                 /OVERPLOT)
     STOP
  ENDIF

  mag_sRate  = 1./(Bt  [1:-1]-Bt  [0:-2])
  eESA_sRate = 1./(Je_z.x[1:-1]-Je_z.x[0:-2])
  iESA_sRate = 1./(Ji_z.x[1:-1]-Ji_z.x[0:-2])

  sRates     = {mag:mag_sRate, $
                eESA:eESA_sRate, $
                iESA:iESA_sRate}

  IF KEYWORD_SET(use_lowRes_time_series) THEN BEGIN

     use_J_time_series = (MEDIAN(sRates.mag)/MEDIAN(sRates.EESA)) GT 1.

  ENDIF

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

           bro = ROUND_TO_NTH_DECIMAL_PLACE(Bt[1:-1]-Bt[0:-2],-5)
           bro = bro[WHERE(ABS(bro) LT 1)]
           distFreq = HISTOGRAM(bro,MIN=MIN(bro),BINSIZE=0.00001, $
                                REVERSE_INDICES=ri, $
                                LOCATIONS=locs)
           junk = MAX(distFreq,ind)
           sPeriod = DOUBLE(locs[ind])
           even_TS = MAKE_EVENLY_SPACED_TIME_SERIES(START_T=Bt[0], $
                                                    STOP_T=Bt[-1], $
                                                    DELTA_T=sPeriod)

           je_z_interp = MAKE_ARRAY(N_ELEMENTS(even_TS),VALUE=!VALUES.F_NaN)
           ji_z_interp = MAKE_ARRAY(N_ELEMENTS(even_TS),VALUE=!VALUES.F_NaN)

           ;;The new way to smooth
           je_z_interp = SMOOTH(je_z.y,3,/NaN)
           ji_z_interp = SMOOTH(ji_z.y,3,/NaN)

           je_z_interp = DATA_CUT({x:je_z.x,y:je_z_interp},Bt)
           ji_z_interp = DATA_CUT({x:ji_z.x,y:ji_z_interp},Bt)

           ;; IF KEYWORD_SET(have_EField) THEN BEGIN
              
           ;; ENDIF

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

           IF KEYWORD_SET(have_EField) THEN BEGIN
              
           ENDIF

           bro        = ROUND_TO_NTH_DECIMAL_PLACE(Bt[1:-1]-Bt[0:-2],-5)
           bro        = bro[WHERE(ABS(bro) LT 1)]
           distFreq   = HISTOGRAM(bro,MIN=MIN(bro),BINSIZE=0.00001, $
                                REVERSE_INDICES=ri, $
                                LOCATIONS=locs)
           junk       = MAX(distFreq,ind)
           sPeriod    = DOUBLE(locs[ind])
           even_TS    = Bt

           origInterp = 0B
           spline     = 1B
           maxDelta_t = 1.5
           interp     = origInterp
           FA_FIELDS_COMBINE,{time:Bt,comp1:Bt}, $
                             {time:Je_z.x,comp1:Je_z.y}, $
                             RESULT=Je_z_interp, $
                             INTERP=interp, $
                             SPLINE=spline, $
                             DELT_T=maxDelta_t, $
                             /TALK

           interp = origInterp
           FA_FIELDS_COMBINE,{time:Bt,comp1:Bt}, $
                             {time:Ji_z.x,comp1:Ji_z.y}, $
                             RESULT=Ji_z_interp, $
                             INTERP=interp, $
                             SPLINE=spline, $
                             DELT_T=maxDelta_t, $
                             /TALK

           IF KEYWORD_SET(have_EField) THEN BEGIN

              interp = origInterp
              FA_FIELDS_COMBINE,{time:Bt,comp1:Bt}, $
                                {time:EField.x,comp1:Ex}, $
                                RESULT=Ex_interp, $
                                INTERP=interp, $
                                SPLINE=spline, $
                                DELT_T=maxDelta_t, $
                                /TALK

              interp = origInterp
              FA_FIELDS_COMBINE,{time:Bt,comp1:Bt}, $
                                {time:EField.x,comp1:Ey}, $
                                RESULT=Ey_interp, $
                                INTERP=interp, $
                                SPLINE=spline, $
                                DELT_T=maxDelta_t, $
                                /TALK

              interp = origInterp
              FA_FIELDS_COMBINE,{time:Bt,comp1:Bt}, $
                                {time:EField.x,comp1:Ez}, $
                                RESULT=Ez_interp, $
                                INTERP=interp, $
                                SPLINE=spline, $
                                DELT_T=maxDelta_t, $
                                /TALK

           ENDIF
           
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
        ;; even_TS = ji_z.x

        ;; even_TS = even_TS[60:-60]
        origInterp = 0B
        spline     = 1B
        maxDelta_t = 0.8

        interp = origInterp
        FA_FIELDS_COMBINE,{time:even_TS,comp1:even_TS}, $
                          {time:Bt,comp1:Bx}, $
                          RESULT=Bx_interp, $
                          INTERP=interp, $
                          SPLINE=spline, $
                          DELT_T=maxDelta_t, $
                          /TALK

        interp = origInterp
        FA_FIELDS_COMBINE,{time:even_TS,comp1:even_TS}, $
                          {time:Bt,comp1:By}, $
                          RESULT=By_interp, $
                          INTERP=interp, $
                          SPLINE=spline, $
                          DELT_T=maxDelta_t, $
                          /TALK

        interp = origInterp
        FA_FIELDS_COMBINE,{time:even_TS,comp1:even_TS}, $
                          {time:Bt,comp1:Bz}, $
                          RESULT=Bz_interp, $
                          INTERP=interp, $
                          SPLINE=spline, $
                          DELT_T=maxDelta_t, $
                          /TALK

        interp = origInterp
        FA_FIELDS_COMBINE,{time:even_TS,comp1:even_TS}, $
                          {time:Je_z.x,comp1:Je_z.y}, $
                          RESULT=Je_z_interp, $
                          INTERP=interp, $
                          SPLINE=spline, $
                          DELT_T=maxDelta_t, $
                          /TALK

        interp = origInterp
        FA_FIELDS_COMBINE,{time:even_TS,comp1:even_TS}, $
                          {time:Ji_z.x,comp1:Ji_z.y}, $
                          RESULT=Ji_z_interp, $
                          INTERP=interp, $
                          SPLINE=spline, $
                          DELT_T=maxDelta_t, $
                          /TALK

        IF KEYWORD_SET(have_EField) THEN BEGIN

           interp = origInterp
           FA_FIELDS_COMBINE,{time:even_TS,comp1:even_TS}, $
                             {time:Bt,comp1:Ex}, $
                             RESULT=Ex_interp, $
                             INTERP=interp, $
                             SPLINE=spline, $
                             DELT_T=maxDelta_t, $
                             /TALK

           interp = origInterp
           FA_FIELDS_COMBINE,{time:even_TS,comp1:even_TS}, $
                             {time:Bt,comp1:Ey}, $
                             RESULT=Ey_interp, $
                             INTERP=interp, $
                             SPLINE=spline, $
                             DELT_T=maxDelta_t, $
                             /TALK

           interp = origInterp
           FA_FIELDS_COMBINE,{time:even_TS,comp1:even_TS}, $
                             {time:Bt,comp1:Ez}, $
                             RESULT=Ez_interp, $
                             INTERP=interp, $
                             SPLINE=spline, $
                             DELT_T=maxDelta_t, $
                             /TALK

           Ex = Ex_interp
           Ey = Ey_interp
           Ez = Ez_interp

        ENDIF

        Bx = Bx_interp
        By = By_interp
        Bz = Bz_interp

     END
  ENDCASE

  ;;And adjust analysis_t1 and analysis_t2, based on shift_nPts
  analysis_t1 += shift_nPts[0]*sPeriod
  analysis_t2 += shift_nPts[1]*sPeriod


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

  strt_i_list  = LIST()
  stop_i_list  = LIST()
  CASE 1 OF
     (KEYWORD_SET(analysis_t1) AND KEYWORD_SET(analysis_t2)): BEGIN
        ;; IF N_ELEMENTS(analysis_t1) GT 1 THEN BEGIN

        analysis_times = [analysis_t1,analysis_t2]

        good_i      = WHERE(FINITE(je_z_improv) AND FINITE(ji_z_improv),nGood,COMPLEMENT=bad_i)
        IF nGood EQ 0 THEN STOP

        FOR k=0,N_ELEMENTS(analysis_t1)-1 DO BEGIN

           good_iTmp  = CGSETINTERSECTION(good_i, $
                                          WHERE((even_TS GE analysis_t1[k]) AND  $
                                                (even_TS LE analysis_t2[k]),nGood))
           GET_STREAKS,good_iTmp,START_i=strt_iiTmp,STOP_I=stop_iiTmp, $
                       OUT_STREAKLENS=streakLensTmp,/QUIET,/NO_PRINT_SUMMARY
           IF streakLensTmp[0] GT 1 THEN BEGIN

              IF KEYWORD_SET(nearest_two_power) THEN BEGIN

                 nHere     = N_ELEMENTS(good_iTmp)
                 n2power   = ROUND(ALOG2(nHere))
                 need      = 2^n2power - nHere
                 origStrt  = good_iTmp[MIN(strt_iiTmp)]

                 CASE 1 OF
                    need LT 0: BEGIN

                       ;;So we're trimming
                       maxGoBack = stop_iiTmp - strt_iiTmp
                       moveBack  = ABS(need) < maxGoBack
                       moveFwd   = ABS(need) - moveBack
                       stop_iiTmp -= moveBack
                       strt_iiTmp += moveFwd

                    END
                    need GT 0: BEGIN

                       ;;So we're adding
                       maxGoFwd = (N_ELEMENTS(good_i)-1) - stop_iiTmp
                       moveFwd  = ABS(need) < maxGoFwd
                       maxMoveBack = strt_iiTmp
                       moveBack    = (ABS(need) - moveFwd) < maxMoveBack
                       stop_iiTmp += moveFwd
                       strt_iiTmp -= moveBack

                    END
                    ELSE     : 
                 ENDCASE

                 ;;Check to make sure it's legit
                 IF N_ELEMENTS(WHERE((stop_iiTmp-strt_iiTmp) EQ (2^n2power-1),/NULL)) NE N_ELEMENTS(stop_iiTmp) THEN BEGIN
                    PRINT,'Wrong again'
                    STOP
                 ENDIF

                 good_iTmp = origStrt + [MIN(strt_iiTmp):MAX(stop_iiTmp)]

              ENDIF

              strt_i_list.Add,good_iTmp[strt_iiTmp]
              stop_i_list.Add,good_iTmp[stop_iiTmp]

           ENDIF

        ENDFOR

        ;; ENDIF ELSE BEGIN
        ;;    good_i      = WHERE(FINITE(je_z_improv) AND FINITE(ji_z_improv) AND $
        ;;                        (even_TS GE analysis_t1) AND (even_TS LE analysis_t2),nGood)
        ;;    GET_STREAKS,good_i,START_i=strt_ii,STOP_I=stop_ii,OUT_STREAKLENS=streakLens
        ;;    strt_i = good_i[strt_ii]
        ;;    stop_i = good_i[stop_ii]
        ;; ENDELSE

     END
     KEYWORD_SET(analysis_t1): BEGIN

        analysis_times = analysis_t1

        IF N_ELEMENTS(analysis_t1) GT 1 THEN BEGIN

           good_i      = WHERE(FINITE(je_z_improv) AND FINITE(ji_z_improv),nGood)
           IF nGood EQ 0 THEN STOP

           ;; strt_i      = !NULL
           ;; stop_i      = !NULL
           FOR k=0,N_ELEMENTS(analysis_t1)-1 DO BEGIN
              good_iTmp  = CGSETINTERSECTION(good_i, $
                                             WHERE((even_TS GE analysis_t1[k]),nGood))
              GET_STREAKS,good_iTmp,START_i=strt_iiTmp,STOP_I=stop_iiTmp, $
                          OUT_STREAKLENS=streakLensTmp,/QUIET,/NO_PRINT_SUMMARY
              IF streakLensTmp[0] GT 1 THEN BEGIN
                 ;; FOR kk=0,N_ELEMENTS(streakLensTmp)-1 DO BEGIN
                 ;; strt_i = [strt_i,good_iTmp[strt_iiTmp]]
                 ;; stop_i = [stop_i,good_iTmp[stop_iiTmp]]
                 strt_i_list.Add,good_iTmp[strt_iiTmp]
                 stop_i_list.Add,good_iTmp[stop_iiTmp]
                 ;; ENDFOR
              ENDIF
           ENDFOR

        ENDIF ELSE BEGIN
           good_i      = WHERE(FINITE(je_z_improv) AND FINITE(ji_z_improv) AND $
                               (even_TS GE analysis_t1),nGood)
           GET_STREAKS,good_i,START_i=strt_ii,STOP_I=stop_ii,OUT_STREAKLENS=streakLens,/QUIET,/NO_PRINT_SUMMARY
           IF streakLens[0] GT 1 THEN BEGIN
              ;; strt_i      = good_i[strt_ii]
              ;; stop_i      = good_i[stop_ii]
              strt_i_list.Add,good_i[strt_ii]
              stop_i_list.Add,good_i[stop_ii]
           ENDIF ELSE BEGIN
              PRINT,'No streaks'
              STOP
           ENDELSE
        ENDELSE

     END
     KEYWORD_SET(analysis_t2): BEGIN

        analysis_times = analysis_t2

        IF N_ELEMENTS(analysis_t2) GT 1 THEN BEGIN

           good_i      = WHERE(FINITE(je_z_improv) AND FINITE(ji_z_improv),nGood)
           IF nGood EQ 0 THEN STOP

           strt_i      = !NULL
           stop_i      = !NULL
           FOR k=0,N_ELEMENTS(analysis_t1)-1 DO BEGIN
              good_iTmp  = CGSETINTERSECTION(good_i, $
                                             WHERE((even_TS LE analysis_t2[k]),nGood))
              GET_STREAKS,good_iTmp,START_i=strt_iiTmp,STOP_I=stop_iiTmp, $
                          OUT_STREAKLENS=streakLensTmp,/QUIET,/NO_PRINT_SUMMARY
              IF streakLensTmp[0] GT 1 THEN BEGIN
                 ;; FOR kk=0,N_ELEMENTS(streakLensTmp)-1 DO BEGIN
                 ;; strt_i = [strt_i,good_iTmp[strt_iiTmp]]
                 ;; stop_i = [stop_i,good_iTmp[stop_iiTmp]]
                 strt_i_list.Add,good_iTmp[strt_iiTmp]
                 stop_i_list.Add,good_iTmp[stop_iiTmp]
                 ;; ENDFOR
              ENDIF
           ENDFOR

        ENDIF ELSE BEGIN
           good_i      = WHERE(FINITE(je_z_improv) AND FINITE(ji_z_improv) AND $
                               (even_TS LE analysis_t2),nGood)
           GET_STREAKS,good_i,START_i=strt_ii,STOP_I=stop_ii,OUT_STREAKLENS=streakLens,/QUIET,/NO_PRINT_SUMMARY
           IF streakLens[0] GT 1 THEN BEGIN
              strt_i_list.Add,good_i[strt_ii]
              stop_i_list.Add,good_i[stop_ii]
           ENDIF ELSE BEGIN
              PRINT,'No streaks'
              STOP
           ENDELSE
        ENDELSE

     END
     ELSE: BEGIN
        good_i      = WHERE(FINITE(je_z_improv) AND FINITE(ji_z_improv),nGood)
        GET_STREAKS,good_i,START_i=strt_ii,STOP_I=stop_ii,OUT_STREAKLENS=streakLens,/QUIET,/NO_PRINT_SUMMARY
        IF streakLens[0] GT 1 THEN BEGIN
           ;; strt_i      = good_i[strt_ii]
           ;; stop_i      = good_i[stop_ii]
           strt_i_list.Add,good_i[strt_ii]
           stop_i_list.Add,good_i[stop_ii]
        ENDIF ELSE BEGIN
           PRINT,'No streaks'
           STOP
        ENDELSE
        ;;Did user provide a streak?
        IF N_ELEMENTS(streakNum) GT 0 THEN longestInd = streakNum

     END
  ENDCASE

  IF N_ELEMENTS(strt_i_list[0]) EQ 0 OR N_ELEMENTS(stop_i_list[0]) EQ 0 THEN BEGIN
     MESSAGE,"Couldn't get any data!",/CONT
     RETURN,0
  ENDIF

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;How many streaks to use?
  CASE 1 OF
     KEYWORD_SET(use_all_streaks): BEGIN

        good_i_list = LIST()
        TArr_list   = LIST()
        T_list      = LIST()
        ;; Ji_z_list   = LIST()
        ;; Je_z_list   = LIST()
        Bx_list     = LIST()
        By_list     = LIST()
        Bz_list     = LIST()
        Jx_list     = LIST()
        Jy_list     = LIST()
        Jz_list     = LIST()

        Ex_list     = LIST()
        Ey_list     = LIST()
        Ez_list     = LIST()
        FOR k=0,N_ELEMENTS(strt_i_list)-1 DO BEGIN
           strt_iTmp = strt_i_list[k]
           stop_iTmp = stop_i_list[k]

           FOR kk=0,N_ELEMENTS(strt_iTmp)-1 DO BEGIN

              ;; good_i_list.Add,(good_i[strt_ii[k]:stop_ii[k]])
              good_i_list.Add,([strt_iTmp[kk]:stop_iTmp[kk]])

              T           = (good_i_list[-1])[-1]-(good_i_list[-1])[0]+1
              PRINT,'N good: ',STRCOMPRESS(T,/RE)

              TArr_list.Add,even_TS[good_i_list[-1]]
              T_list.   Add,N_ELEMENTS(TArr_list[-1])

              ;; Ji_z_list.Add,Ji_z_improv[good_i_list[-1]]
              ;; Je_z_list.Add,Je_z_improv[good_i_list[-1]]
              Ji_zTmp = Ji_z_improv[good_i_list[-1]]
              Je_zTmp = Je_z_improv[good_i_list[-1]]

              BxTmp = Bx[good_i_list[-1]]
              ByTmp = By[good_i_list[-1]]
              BzTmp = Bz[good_i_list[-1]]

              Bx_list.Add,BxTmp - BxTmp[0]
              By_list.Add,ByTmp - ByTmp[0]
              Bz_list.Add,BzTmp - BzTmp[0]

              ;;Sorry, Jx and Jy
              Jx_list.Add,MAKE_ARRAY(T_list[-1],VALUE=0.) * jFactor
              Jy_list.Add,MAKE_ARRAY(T_list[-1],VALUE=0.) * jFactor

              Jz_list.Add,(Ji_zTmp + Je_zTmp) * jFactor

              IF KEYWORD_SET(have_EField) THEN BEGIN

                 ExTmp = Ex[good_i_list[-1]]
                 EyTmp = Ey[good_i_list[-1]]
                 EzTmp = Ez[good_i_list[-1]]

                 Ex_list.Add,ExTmp - ExTmp[0]
                 Ey_list.Add,EyTmp - EyTmp[0]
                 Ez_list.Add,EzTmp - EzTmp[0]

              ENDIF

           ENDFOR
           
        ENDFOR

        Bx    = Bx_list
        By    = By_list
        Bz    = Bz_list
        T     = T_list
        TArr  = TArr_list
        Jx    = Jx_list
        Jy    = Jy_list
        Jz    = Jz_list

        Ex    = Ex_list
        Ey    = Ey_list
        Ez    = Ez_list

     END
     ELSE: BEGIN

        IF N_ELEMENTS(strt_i_list) GT 0 THEN BEGIN
           streakLens = stop_i_list[0]-strt_i_list[0]
           long       = MAX(streakLens,longestInd)
        ENDIF ELSE BEGIN
           longestInd = 0
        ENDELSE

        good_i      = [strt_i_list[0]:stop_i_list[0]]
        ;; good_i      = (good_i[strt_ii[longestInd]:stop_ii[longestInd]])

        T           = good_i[-1]-good_i[0]+1
        PRINT,'N good: ',STRCOMPRESS(T,/RE)

        TArr        = even_TS[good_i]
        T           = N_ELEMENTS(TArr)
        Ji_z_improv = Ji_z_improv[good_i]
        Je_z_improv = Je_z_improv[good_i]

        Bx = Bx[good_i]
        By = By[good_i]
        Bz = Bz[good_i]

        ;; Bx = Bx - Bx[0]
        ;; By = By - By[0]
        ;; Bz = Bz - Bz[0]

        IF KEYWORD_SET(have_EField) THEN BEGIN

           Ex = Ex[good_i]
           Ey = Ey[good_i]
           Ez = Ez[good_i]

        ENDIF

        ;;Sorry, Jx and Jy
        Jx = MAKE_ARRAY(T,VALUE=0.) * jFactor
        Jy = MAKE_ARRAY(T,VALUE=0.) * jFactor

        Jz = (Ji_z_improv + Je_z_improv) * jFactor

     END
  ENDCASE

  b_is_list         = SIZE(Bx,/TYPE) EQ 11
  can_convert       = b_is_list ? (NDIMEN(Bx[0]) EQ 1) : NDIMEN(Bx) EQ 1
  IF KEYWORD_SET(nearest_two_power) AND can_convert THEN BEGIN

     ;;Temporarily unpack, if need be
     IF b_is_list THEN BEGIN

        Bx    = Bx[0]
        By    = By[0]
        Bz    = Bz[0]
        T     = T[0]
        TArr  = TArr[0]
        Jx    = Jx[0]
        Jy    = Jy[0]
        Jz    = Jz[0]

        IF KEYWORD_SET(have_EField) THEN BEGIN

           Ex    = Ex[0]
           Ey    = Ey[0]
           Ez    = Ez[0]

        ENDIF

     ENDIF
     
     nPoints   = N_ELEMENTS(Bx)
     n2power   = ROUND(ALOG2(nPoints))
     newPoints = 2^n2power

     ;;Shrink or else expand
     CASE 1 OF
        newPoints GT nPoints: BEGIN

           ;;pad with zeros
           padding = MAKE_ARRAY(newPoints - nPoints,/DOUBLE,VALUE=0.D)

           Bx   = [Bx,padding]
           By   = [By,padding]
           Bz   = [Bz,padding]
           T    = newPoints
           TArr = [TArr,(DINDGEN(newPoints - nPoints)+1)*(TArr[-1]-TArr[-2])+TArr[-1]]
           Jx   = [Jx,padding]
           Jy   = [Jy,padding]
           Jz   = [Jz,padding]

           IF KEYWORD_SET(have_EField) THEN BEGIN

              Ex   = [Ex,padding]
              Ey   = [Ey,padding]
              Ez   = [Ez,padding]

           ENDIF

        END
        newPoints LT nPoints: BEGIN

           dropEm = newPoints - nPoints

           Bx   = Bx[0:dropEm]
           By   = By[0:dropEm]
           Bz   = Bz[0:dropEm]
           T    = newPoints
           TArr = TArr[0:dropEm]
           Jx   = Jx[0:dropEm]
           Jy   = Jy[0:dropEm]
           Jz   = Jz[0:dropEm]

           IF KEYWORD_SET(have_EField) THEN BEGIN

              Ex   = Ex[0:dropEm]
              Ey   = Ey[0:dropEm]
              Ez   = Ez[0:dropEm]

           ENDIF


        END
        ELSE: 
     ENDCASE

     ;;Repack, if need be
     IF b_is_list THEN BEGIN

        Bx    = LIST(Bx  ) 
        By    = LIST(By  )
        Bz    = LIST(Bz  )
        T     = LIST(T   )
        TArr  = LIST(TArr)
        Jx    = LIST(Jx  )
        Jy    = LIST(Jy  )
        Jz    = LIST(Jz  )

        IF KEYWORD_SET(have_EField) THEN BEGIN

           Ex    = LIST(Ex  ) 
           Ey    = LIST(Ey  )
           Ez    = LIST(Ez  )

        ENDIF

     ENDIF

  ENDIF

  IF KEYWORD_SET(have_EField) THEN BEGIN
     
     EField = [[Ex],[Ey],[Ez]]

  ENDIF

  RETURN,1

END

;********************************************************************
;MAIN PROGRAM
;********************************************************************
; written by P. M. Bellan, April 2016, MS 128-95 Caltech, Pasadena CA 91125
; this code is free to use, but please acknowledge P. M. Bellan if used

PRO SINGLE_SPACECRAFT_K_MEASUREMENT_FAST, $
   PARSE_B_AND_J_SAVEFILE=parse_B_and_J_saveFile, $
   B_AND_J_FILE=saveFile, $
   USE_TIMEBAR_TIME__FROM_FILE=use_timeBar_time__from_file, $
   CUSTOM_T1=custom_t1, $
   CUSTOM_T2=custom_t2, $
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
   FFT__NEAREST_TWO_POWER=nearest_two_power, $
   FFTSIZE=FFTsize, $
   FFTPERCENT=FFTpercent,$
   WHICH_FFTS=which_FFTs, $
   OUT_FREQS=out_freqs, $
   OUT_KX=out_kx, $
   OUT_KY=out_ky, $
   OUT_KZ=out_kz, $
   OUT_KP=out_kP, $
   OUT_ANGLE_KP=out_kPAngle, $
   OUT_INDS=out_inds, $
   OUT_TARR=out_TArr, $
   OUT_USEDINDS=out_usedInds, $
   OUT_BX=out_Bx, $
   OUT_BY=out_By, $
   OUT_BZ=out_Bz, $
   OUT_JX=out_Jx, $
   OUT_JY=out_Jy, $
   OUT_JZ=out_Jz, $
   OUT_AVGJXBNRM=out_avgJxBNrm


  COMPILE_OPT idl2,strictarrsubs

  doPlot    = ~KEYWORD_SET(no_plots)

  ;;diagnostics?
  ;; diag      = 1
  ;; diagInd   = 0
  
  splitFFTs = ~KEYWORD_SET(combine_and_average_intervals)

  IF N_ELEMENTS(extra_suffix) EQ 0 THEN extra_suffix = ''

  ;; IF N_ELEMENTS(double_calc) EQ 0 THEN double_calc = 1

  ;;N points to smooth k components with
  smInd         = KEYWORD_SET(smInd                 ) ? smInd                  : 5
  dbSmInd       = KEYWORD_SET(dbSmInd               ) ? dbSmInd                : 10
  edge_truncate = KEYWORD_SET(kSmooth__edge_truncate) ? kSmooth__edge_truncate : 0
  edge_mirror   = KEYWORD_SET(kSmooth__edge_mirror  ) ? kSmooth__edge_mirror   : 1
  edge_wrap     = KEYWORD_SET(kSmooth__edge_wrap    ) ? kSmooth__edge_wrap     : 0

  IF N_ELEMENTS(oddness_check) EQ 0 THEN oddness_check = 1

  ;;select one of three examples
  ;; First example : A few discrete k components
  ;; Second example: Continuum of k components following some complicated function of omega
  ;; Third example : Same as second example, except noise added

  ;; example   = 3

  saveVar  = 'dB_fac_V'
  IF KEYWORD_SET(use_dB_fac) THEN saveVar = 'dB_fac'

  suff = 'TEST'
  CASE 1 OF
     KEYWORD_SET(use_all_streaks): BEGIN
     END
     ELSE: BEGIN
        suff = 'Bellan'+'-'+ $
               extra_suffix+'-'+ $
               saveVar
     END
  ENDCASE

  saveDir  = '/SPENCEdata/Research/Satellites/FAST/single_sc_wavevector/saves_output_etc/'

  CASE 1 OF
     KEYWORD_SET(use_all_streaks): BEGIN
        suff += '-all_streaks'
     END
     ELSE: BEGIN
        IF N_ELEMENTS(streakInd) GT 0 THEN BEGIN
           suff += '-streak_' + STRCOMPRESS(streakInd,/REMOVE_ALL)
        ENDIF
     END
  ENDCASE

  ;;Aligned to mag or j time series?
  CASE 1 OF
     KEYWORD_SET(use_lowRes_time_series): BEGIN
        suff += '-lR_ts'
     END
     KEYWORD_SET(use_J_time_series): BEGIN
        suff += '-jz_ts'
     END
     ELSE: BEGIN
        suff += '-mag_ts'
     END
  ENDCASE

  IF KEYWORD_SET(use_avged_for_smooth) THEN BEGIN
     suff += '-avgd_iSm'
  ENDIF

  IF KEYWORD_SET(smooth_J_dat) THEN BEGIN
     suff += '-sm_je'
  ENDIF

  IF KEYWORD_SET(plot_kperp_magnitude_for_kz) THEN BEGIN
     ;; suff += '-kPerp'
  ENDIF

  IF KEYWORD_SET(tmpFFTsize) THEN BEGIN
     suff += '-FFTsize_' + STRCOMPRESS(FFTsize,/REMOVE_ALL)
  ENDIF

  IF KEYWORD_SET(hanning    ) THEN BEGIN
     PRINT,"You shouldn't apply a window to correlation functions ..."
     WAIT,2
     suff += '-hanning'
  ENDIF
  IF KEYWORD_SET(double_calc) THEN suff += '-double_arithmetic'
  IF KEYWORD_SET(bonus_suff ) THEN suff += (STRMATCH(bonus_suff,'-*') ? '' : '-' ) + bonus_suff

  example = 1


  CASE 1 OF
     KEYWORD_SET(shift_nPts): BEGIN
        CASE N_ELEMENTS(shift_nPts) OF
           1: BEGIN
              shift_nPts = [shift_nPts,shift_nPts]
              backShifts_for_avging = shift_nPts
              fwdShifts_for_avging  = shift_nPts
           END
           2: BEGIN
              backShifts_for_avging = shift_nPts[0]
              fwdShifts_for_avging  = shift_nPts[1]
           END
           ELSE: BEGIN
              PRINT,"Bogus! No more than two shifties"
              STOP
           END
        ENDCASE
     END
     ELSE: BEGIN
        IF N_ELEMENTS(backShifts_for_avging) EQ 0 THEN BEGIN
           backShifts_for_avging = 0
      ENDIF
        
        IF N_ELEMENTS(fwdShifts_for_avging) EQ 0 THEN BEGIN
           fwdShifts_for_avging = 0
        ENDIF
     END
  ENDCASE

  nBack = N_ELEMENTS(backShifts_for_avging)
  IF KEYWORD_SET(lock_shifts) THEN BEGIN
     nFwd  = 1
  ENDIF ELSE BEGIN
     nFwd  = N_ELEMENTS(fwdShifts_for_avging)
  ENDELSE

  IF nBack GT 1 THEN BEGIN
     junk       = MIN(ABS(backShifts_for_avging),minJJInd)
  ENDIF ELSE BEGIN
     minJJInd   = 0
  ENDELSE
  IF nFwd GT 1 THEN BEGIN
     junk       = MIN(ABS(fwdShifts_for_avging),minKKInd)
  ENDIF ELSE BEGIN
     minKKInd   = 0
  ENDELSE
  
  TArrList      = LIST()
  freqList      = LIST()
  fDiffList     = LIST()
  kxList        = LIST()
  kyList        = LIST()
  kzList        = LIST()
  ;; kPList        = LIST()
  ;; kPAngleList   = LIST()
  indsList      = LIST()
  kzList        = LIST()
  kzList        = LIST()
  kzList        = LIST()
  avgJxBNrmList = LIST()
  normList      = LIST()

  avgCount      = 0

  FOR jj=0,nBack-1 DO BEGIN
     FOR kk=0,nFwd-1 DO BEGIN

        IF KEYWORD_SET(lock_shifts) THEN BEGIN
           shift_nPts = [backShifts_for_avging[jj],backShifts_for_avging[jj]]
        ENDIF ELSE BEGIN
           shift_nPts = [backShifts_for_avging[jj],fwdShifts_for_avging[kk]]
        ENDELSE

        tmpFFTSize = KEYWORD_SET(FFTSize) ? FFTSize : !NULL
        tmpSuff = STRING(FORMAT='("__shft_",I0,"_",I0)',backShifts_for_avging[jj],fwdShifts_for_avging[kk])
        PRINT,FORMAT='("********",T10,I0,T20,A0,T40,"*********")',avgCount,tmpSuff


        IF KEYWORD_SET(parse_B_and_J_saveFile) THEN BEGIN

           IF ~CHUNK_SAVE_FILE(T,TArr,Bx,By,Bz,Jx,Jy,Jz,unitFactor,sPeriod,saveVar, $
                               B_AND_J_FILE=saveFile, $
                               SAVEDIR=saveDir, $
                               USE_TIMEBAR_TIME__FROM_FILE=use_timeBar_time__from_file, $
                               CUSTOM_T1=custom_t1, $
                               CUSTOM_T2=custom_t2, $
                               SHIFT_NPTS=shift_nPts, $
                               USE_LOWRES_TIME_SERIES=use_lowRes_time_series, $
                               USE_J_TIME_SERIES=use_J_time_series, $
                               SMOOTH_J_DAT_TO_B=smooth_J_dat, $
                               PRESMOOTH_MAG=presmooth_mag, $
                               PREPLOT_CURRENTS_AND_STOP=prePlot_currents_and_stop, $
                               STREAKNUM=streakNum, $
                               OUT_STREAKNUM=streakInd, $
                               USE_ALL_STREAKS=use_all_streaks, $
                               USE_DB_FAC=use_dB_fac, $
                               EFIELD=EField, $
                               FFT__NEAREST_TWO_POWER=nearest_two_power, $
                               SRATES=sRates) $
           THEN BEGIN
              MESSAGE,"Failed to parse!",/CONT
              RETURN
           ENDIF

           CASE 1 OF
              KEYWORD_SET(use_all_streaks): BEGIN
                 sRate = !NULL 
                 FOR k=0,N_ELEMENTS(Bx)-1 DO BEGIN
                    sRate = [sRate,1./( (TArr[k])[1:-1]-(TArr[k])[0:-2] )]
                 ENDFOR
              END
              ELSE: BEGIN
                 sRate = 1./(TArr[1:-1]-TArr[0:-2])

              END
           ENDCASE

           maxFreq = (MIN([MEDIAN(sRates.mag),MEDIAN(sRates.eESA),MEDIAN(sRates.iESA)]))/2.
           maxFreq = 4.

        ENDIF ELSE BEGIN
           unitFactor = 1       ;Don't adjust k in this case
        ENDELSE

        IF KEYWORD_SET(example_mode) THEN BEGIN
           SETUP_EXAMPLE,T,TArr,Bx,By,Bz,Jx,Jy,Jz,unitFactor,sPeriod,saveVar, $
                         SAVEDIR=saveDir, $
                         EXAMPLE=example
        ENDIF

        ;;read data from file
        IF ~KEYWORD_SET(parse_B_and_J_saveFile) THEN BEGIN
           filename = saveDir+B_J_file

           OPENR,1,filename
           PRINT,'opening ', filename
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

        ;;Now do some calcs
        defFFTsize = T[0]
        IF ( N_ELEMENTS(tmpFFTsize) EQ 0 ) AND KEYWORD_SET(splitFFTs) THEN BEGIN
           tmpFFTsize = defFFTsize ;The other way for tmpFFTsize to be set is user input, toward the beginning 
        ENDIF

        IF KEYWORD_SET(FFTpercent) THEN BEGIN
           IF FFTpercent GT 1 THEN FFTpercent /= 100.
           tmpFFTsize = FIX(T[0]*FFTpercent)
        ENDIF

        CASE 1 OF
           KEYWORD_SET(tmpFFTsize): BEGIN

              IF KEYWORD_SET(use_all_streaks) THEN BEGIN

                 ;;First get the total number of FFTs we're going to do for each streak of measurements
                 nFFTsArr = !NULL
                 FOR kk=0,N_ELEMENTS(Bx)-1 DO BEGIN
                    TTmp     = T[kk]
                    lastInd  = TTmp-1

                    nFFTsArr = [nFFTsArr,TTmp/tmpFFTSize]
                    ;; FOR k=0,nFFTs-1 DO BEGIN
                    ;;    tmpI = [(k*tmpFFTSize):( ((k+1)*tmpFFTSize-1) < lastInd )]                 
                    ;; ENDFOR
                 ENDFOR
                 nFFTs     = TOTAL(nFFTsArr)

                 fArr      = MAKE_ARRAY(nFFTs,tmpFFTSize)
                 kxArr     = MAKE_ARRAY(nFFTs,tmpFFTSize)
                 kyArr     = MAKE_ARRAY(nFFTs,tmpFFTSize)
                 kzArr     = MAKE_ARRAY(nFFTs,tmpFFTSize)
                 avgJxBArr = MAKE_ARRAY(nFFTs,3)
                 normArr   = MAKE_ARRAY(nFFTs,3)
                 ;; kpArr     = MAKE_ARRAY(nFFTs,tmpFFTSize)

                 FFTCount  = 0
                 include_i = !NULL
                 prevFFTs  = 0
                 FFTi_list = LIST()

                 ;;Remember, Bx is a LIST (not an array) and so contains sets of streaks of Bx measurements
                 ;;Here we loop over streaks of measurements
                 FOR kk=0,N_ELEMENTS(Bx)-1 DO BEGIN

                    TTmp  = T[kk]
                    TArrTmp = TArr[kk]

                    BxTmp = Bx[kk]
                    ByTmp = By[kk]
                    BzTmp = Bz[kk]

                    JxTmp = Jx[kk]
                    JyTmp = Jy[kk]
                    JzTmp = Jz[kk]
                    
                    ;; lastInd  = TTmp-1
                    lastInd  = TTmp-1
                    nArr     = !NULL
                    PRINT,"Interval: ",kk

                    ;;If we're doing multiple FFTs per streak, it happens here. Otherwise this loop is only traversed once.
                    FOR k=FFTCount,FFTCount+nFFTsArr[kk]-1 DO BEGIN
                       tmpI = [((k-prevFFTs)*tmpFFTSize):( ((k-prevFFTs+1)*tmpFFTSize-1) < lastInd )]
                       nTmp = N_ELEMENTS(tmpI)
                       nArr = [nArr,nTmp]
                       PRINT,k," ",nTmp

                       IF nTmp LT tmpFFTSize THEN CONTINUE
                       include_i = [include_i,k]

                       BELLAN_2016__BRO,nTmp,JxTmp[tmpI],JyTmp[tmpI],JzTmp[tmpI], $
                                        BxTmp[tmpI],ByTmp[tmpI],BzTmp[tmpI], $
                                        freq,kx,ky,kz,kP, $
                                        SPERIOD=sPeriod, $
                                        UNITFACTOR=unitFactor, $
                                        PLOT_KPERP_MAGNITUDE_FOR_KZ=plot_kperp_magnitude_for_kz, $
                                        DOUBLE_CALC=double_calc, $
                                        HANNING=hanning, $
                                        EFIELD=EField, $
                                        ODDNESS_CHECK=oddness_check, $
                                        OUT_NORM=norm, $
                                        OUT_AVGJXB=avgJxBtotal

                       this = ARRAY_INDICES(TRANSPOSE(fArr),(k*N_ELEMENTS(fArr[0,*])+LINDGEN(N_ELEMENTS(tmpI))))
                       ;; this = [this[1,*],this[0,*]]
                       ;; PRINT,this
                       ;; fArr [k,*]  = freq
                       ;; kxArr[k,*] = kx
                       ;; kyArr[k,*] = ky
                       ;; kzArr[k,*] = kz
                       ;; kPArr[k,*] = kP
                       ;; fArr [this[0,*],this[1,*]]  = freq
                       ;; kxArr[this[0,*],this[1,*]] = kx
                       ;; kyArr[this[0,*],this[1,*]] = ky
                       ;; kzArr[this[0,*],this[1,*]] = kz
                       ;; kPArr[this[0,*],this[1,*]] = kP

                       CASE NDIMEN(this) OF
                          1: BEGIN

                             fArr [this] = TEMPORARY(freq)
                             kxArr[this] = TEMPORARY(kx  )
                             kyArr[this] = TEMPORARY(ky  )
                             kzArr[this] = TEMPORARY(kz  )
                             ;; kPArr[this] = kP

                          END
                          2: BEGIN

                             fArr [this[1,*],this[0,*]] = TEMPORARY(freq)
                             kxArr[this[1,*],this[0,*]] = TEMPORARY(kx  )
                             kyArr[this[1,*],this[0,*]] = TEMPORARY(ky  )
                             kzArr[this[1,*],this[0,*]] = TEMPORARY(kz  )
                             ;; kPArr[this] = kP

                          END
                       ENDCASE

                       ;;avgJxB doesn't depend on FFT dimensions
                       avgJxBArr[k,*] = TEMPORARY(avgJxBtotal)
                       normArr[k,*]   = REPLICATE(TEMPORARY(norm),3)
                       
                    ENDFOR

                    FFTCount += nFFTsArr[kk]
                    prevFFTs = FFTCount
                    FFTi_list.Add,TEMPORARY(tmpI)
                 ENDFOR

                 fArr      = fArr[include_i,*]
                 kxArr     = kxArr[include_i,*]
                 kyArr     = kyArr[include_i,*]
                 kzArr     = kzArr[include_i,*]

                 freq      = MEAN(TEMPORARY(fArr) ,DIMENSION=1)
                 kx        = MEAN(TEMPORARY(kxArr),DIMENSION=1)
                 ky        = MEAN(TEMPORARY(kyArr),DIMENSION=1)
                 kz        = MEAN(TEMPORARY(kzArr),DIMENSION=1)
                 kP        = SQRT(kx*kx+ky*ky)

                 avgJxBNrm = MEAN(TEMPORARY(avgJxBArr)/TEMPORARY(normArr),DIMENSION=1)

                 usedInds  = LIST_TO_1DARRAY(FFTi_list,/WARN,/SKIP_NEG1_ELEMENTS,/SKIP_NANS)

              ENDIF ELSE BEGIN

                 eMulieribus   = T/tmpFFTSize
                 IF N_ELEMENTS(which_FFTs) GT 0 THEN BEGIN

                    IF (N_ELEMENTS(which_FFTs) GT eMulieribus) OR which_FFTs[0] GT eMulieribus THEN BEGIN
                       PRINT,"You're WRONG"
                       STOP
                    ENDIF

                    nFFTs = N_ELEMENTS(which_FFTs)

                    FFTSuff = '-FFT_itvls'
                    FOR k=0,nFFTs-1 DO FFTSuff += STRING(FORMAT='("_",I0)',which_FFTs[k])
                    suff   += FFTSuff

                 ENDIF ELSE BEGIN

                    nFFTs = eMulieribus

                 ENDELSE

                 lastInd   = T-1
                 nArr      = !NULL
                 fArr      = MAKE_ARRAY(nFFTs,tmpFFTSize)
                 kxArr     = MAKE_ARRAY(nFFTs,tmpFFTSize)
                 kyArr     = MAKE_ARRAY(nFFTs,tmpFFTSize)
                 kzArr     = MAKE_ARRAY(nFFTs,tmpFFTSize)
                 avgJxBArr = MAKE_ARRAY(nFFTs,3)
                 normArr   = MAKE_ARRAY(nFFTs,3)
                 ;; kpArr     = MAKE_ARRAY(nFFTs,tmpFFTSize)
                 FFTCount  = 0
                 FFTi_list = LIST()
                 FOR k=0,eMulieribus-1 DO BEGIN

                    IF N_ELEMENTS(which_FFTs) GT 0 THEN BEGIN
                       IF (WHERE(k EQ which_FFTs))[0] EQ -1 THEN CONTINUE
                    ENDIF

                    tmpI = [(k*tmpFFTSize):( ((k+1)*tmpFFTSize-1) < lastInd )]
                    nTmp = N_ELEMENTS(tmpI)
                    nArr = [nArr,nTmp]
                    PRINT,FORMAT='(A0,T10,A0,T20,A0,T45,A0,T72,A0)',"Itvl","nPts","Start T","Stop T","sFreq"
                    PRINT,FORMAT='(I0,T10,I0,T20,A0,T45,A0,T72,F0.3)',k,nTmp, $
                          TIME_TO_STR(TArr[tmpI[0]],/MS),TIME_TO_STR(TArr[tmpI[-1]],/MS), $
                          1.D/(TArr[1]-TArr[0])
                    BELLAN_2016__BRO,nTmp,Jx[tmpI],Jy[tmpI],Jz[tmpI],Bx[tmpI],By[tmpI],Bz[tmpI], $
                                     freq,kx,ky,kz,kP, $
                                     SPERIOD=sPeriod, $
                                     UNITFACTOR=unitFactor, $
                                     PLOT_KPERP_MAGNITUDE_FOR_KZ=plot_kperp_magnitude_for_kz, $
                                     DOUBLE_CALC=double_calc, $
                                     HANNING=hanning, $
                                     EFIELD=EField, $
                                     ODDNESS_CHECK=oddness_check, $
                                     OUT_NORM=norm, $
                                     OUT_AVGJXB=avgJxBtotal

                    fArr [FFTCount,*]     = TEMPORARY(freq       )
                    kxArr[FFTCount,*]     = TEMPORARY(kx         )
                    kyArr[FFTCount,*]     = TEMPORARY(ky         )
                    kzArr[FFTCount,*]     = TEMPORARY(kz         )

                    avgJxBArr[FFTCount,*] = TEMPORARY(avgJxBtotal)
                    normArr[FFTCount,*]   = REPLICATE(TEMPORARY(norm       ),3)

                    FFTCount++

                    FFTi_list.Add,TEMPORARY(tmpI)
                 ENDFOR

                 keepEm        = WHERE(nArr EQ MEDIAN(nArr),nKeep)
                 freq          = MEAN(fArr ,DIMENSION=1)
                 kx            = MEAN(kxArr,DIMENSION=1)
                 ky            = MEAN(kyArr,DIMENSION=1)
                 kz            = MEAN(kzArr,DIMENSION=1)
                 kP            = SQRT(kx*kx+ky*ky)

                 avgJxBNrm     = MEAN(TEMPORARY(avgJxBArr)/TEMPORARY(normArr),DIMENSION=1)

                 usedInds      = LIST_TO_1DARRAY(FFTi_list,/WARN,/SKIP_NEG1_ELEMENTS,/SKIP_NANS)

              ENDELSE

           END
           ELSE: BEGIN

              BELLAN_2016__BRO,T,Jx,Jy,Jz,Bx,By,Bz, $
                               freq,kx,ky,kz,kP, $
                               SPERIOD=sPeriod, $
                               UNITFACTOR=unitFactor, $
                               PLOT_KPERP_MAGNITUDE_FOR_KZ=plot_kperp_magnitude_for_kz, $
                               DOUBLE_CALC=double_calc, $
                               HANNING=hanning, $
                               EFIELD=EField, $
                               ODDNESS_CHECK=oddness_check, $
                               OUT_NORM=norm, $
                               OUT_AVGJXB=avgJxBtotal

              avgJxBNrm        = TEMPORARY(avgJxBtotal)/TEMPORARY(norm)
              usedInds         = LINDGEN(T)

           END
        ENDCASE

        ;;Put 'em in km^-1
        kx     *= 1000.
        ky     *= 1000.
        kz     *= 1000.
        kP     *= 1000.

        bro = 1
        IF KEYWORD_SET(bro) THEN BEGIN
           tmp  = SORT(freq)
           freq = freq[tmp]
           kx   = kx[tmp]
           ky   = ky[tmp]
           kz   = kz[tmp]
           kP   = kP[tmp]
        ENDIF

        IF KEYWORD_SET(use_all_streaks) THEN BEGIN
           Tarr = LIST_TO_1DARRAY(Tarr)
           Bx   = LIST_TO_1DARRAY(Bx)
           By   = LIST_TO_1DARRAY(By)
           Bz   = LIST_TO_1DARRAY(Bz)
           Jx   = LIST_TO_1DARRAY(Jx)
           Jy   = LIST_TO_1DARRAY(Jy)
           Jz   = LIST_TO_1DARRAY(Jz)
        END

        ;;Kperp angle
        kPAngle = ATAN(ky,kx)*!RADEG
        IF KEYWORD_SET(kP__angleRange) THEN BEGIN

           rotate_kPA = 0

           kPAngle    = NORMALIZE_ANGLE(kPAngle,MIN(kP__angleRange),MAX(kP__angleRange), $
                                        /DEGREE)

           ;; those = WHERE(kPAngle LT MIN(kP__angleRange),nThose)
           ;; IF nThose GT 0 THEN BEGIN
           ;;    kPAngle[those] = (kPAngle[those] + MAX(kP__angleRange)) MOD MAX(kP__angleRange)
           ;; ENDIF

           ;; kPAngle = UNWRAP(kPAngle,DIVISOR=360)

        ENDIF ELSE BEGIN
           histo   = HISTOGRAM(kPAngle,BINSIZE=90,MIN=-180,MAX=180,LOCATIONS=locs)
           rotate_kPA = (histo[0]+histo[3]) GT (histo[1]+histo[2])
           IF rotate_kPA THEN BEGIN
              kPAngle = (kPAngle + 360) MOD 360
           ENDIF
        ENDELSE

        ;;Rotation matrix?
        ;; angle      = 20.*!DTOR
        ;; kxprime    = kx*COS(angle) - ky*SIN(angle)
        ;; kyprime    = kx*SIN(angle) + ky*COS(angle)
        ;; kx         = TEMPORARY(kxprime)
        ;; ky         = TEMPORARY(kyprime)
        
        freqList.Add,freq
        fDiffList.Add,freq[1]-freq[0]
        kxList.Add,kx
        kyList.Add,ky
        kzList.Add,kz
        ;; kPList.Add,kP
        ;; kPAngleList.Add,kPAngle
        ;; indsList.Add,inds

        avgJxBNrmList.Add,avgJxBNrm
        
        IF (kk EQ minKKInd) AND (jj EQ minJJInd) THEN BEGIN

           ;;Get time string
           milli = 1
           CASE SIZE(Tarr,/TYPE) OF
              11: BEGIN
                 tString = '__' + ((STRSPLIT(TIME_TO_STR(TArr[0,usedInds[0]],MS=milli),'/',/EXTRACT))[1]).Replace(':','_') + '-' + $
                           ((STRSPLIT(TIME_TO_STR(TArr[0,usedInds[-1]],MS=milli),'/',/EXTRACT))[1]).Replace(':','_')
              END
              ELSE: BEGIN
                 tString = '__' + ((STRSPLIT(TIME_TO_STR(TArr[usedInds[0]],MS=milli),'/',/EXTRACT))[1]).Replace(':','_') + '-' + $
                           ((STRSPLIT(TIME_TO_STR(TArr[usedInds[-1]],MS=milli),'/',/EXTRACT))[1]).Replace(':','_')
              END
           ENDCASE
           tString = tString.Replace('.','_')
           suff   += tString

           orig_freq      = freq
           orig_kx        = kx
           orig_ky        = ky
           orig_kz        = kz
           orig_kP        = kP
           orig_kPAngle   = kPAngle
           ;; orig_inds      = inds

           send_TArr      = TArr
           send_usedInds  = usedInds
           send_Bx        = Bx
           send_By        = By
           send_Bz        = Bz
           send_Jx        = Jx
           send_Jy        = Jy
           send_Jz        = Jz

        ENDIF

        avgCount++

     ENDFOR
  ENDFOR
  
  CASE avgCount OF
     1: BEGIN
        ;;Disregard lists; we'll just take it straight
     END
     ELSE: BEGIN

        freqs     = LIST_TO_1DARRAY(freqList,/WARN)
        kxs       = LIST_TO_1DARRAY(kxList,/WARN)
        kys       = LIST_TO_1DARRAY(kyList,/WARN)
        kzs       = LIST_TO_1DARRAY(kzList,/WARN)
        ;; kPs       = LIST_TO_1DARRAY(kPList,/WARN)
        ;; kPAngles  = LIST_TO_1DARRAY(kPAngleList,/WARN)
        ;; indss  = LIST_TO_1DARRAY(indsList,/WARN)
        ;; kzs    = LIST_TO_1DARRAY(kzList,/WARN)
        ;; kzs    = LIST_TO_1DARRAY(kzList,/WARN)
        ;; kzs    = LIST_TO_1DARRAY(kzList,/WARN)

        avgJxBNrm = LIST_TO_1DARRAY(avgJxBNrmList,/WARN,/PRESERVE_DIMENSIONALITY)

        ;;Pick up rebin size from user, or else from data
        IF KEYWORD_SET(avg_binSize) THEN BEGIN

           binSz  = avg_binSize

        ENDIF ELSE BEGIN

        ;; binSz     = 0.2
        ;; binSz     = orig_freq[1]-orig_freq[0]
           binSz     = MAX(LIST_TO_1DARRAY(fDiffList,/WARN))

        ENDELSE

        kx        = HIST1D(freqs,kxs,BINSIZE=binSz,OBIN=freq)/avgCount
        ky        = HIST1D(freqs,kys,BINSIZE=binSz,OBIN=freq)/avgCount
        kz        = HIST1D(freqs,kzs,BINSIZE=binSz,OBIN=freq)/avgCount
        ;; kP        = HIST1D(freqs,kPs,BINSIZE=binSz,OBIN=freq)/avgCount
        ;; kPAngle   = HIST1D(freqs,kPAngles,BINSIZE=binSz)/avgCount
        kP        = SQRT(kx*kx+ky*ky)
        kPAngle = ATAN(ky,kx)*!RADEG

     END
  ENDCASE

  IF doPlot THEN BEGIN

     IF KEYWORD_SET(use_avged_for_smooth) THEN BEGIN

        IF KEYWORD_SET(plot_smoothed_ks) THEN BEGIN
           plot_smoothed_ks = {freq:freq, $
                               kx:kx, $
                               ky:ky}

        ENDIF

        IF KEYWORD_SET(kx_vs_ky__plot_smoothed) THEN BEGIN

           kx_vs_ky__plot_smoothed = {kx:kx, $
                                      ky:ky}

        ENDIF

        IF KEYWORD_SET(kP_angle__plot_smoothed) THEN BEGIN

           kP_angle__plot_smoothed = {freq:freq, $
                                      angle:kPAngle}

        ENDIF
        
        ;;Now give 'em what we got
        freq       =  orig_freq   
        kx         =  orig_kx     
        ky         =  orig_ky     
        kz         =  orig_kz     
        kP         =  orig_kP     
        kPAngle    =  orig_kPAngle
        ;; inds       =  orig_inds   

     ENDIF

     PLOT_SINGLE_SPACECRAFT_K_MEASUREMENT,send_TArr,freq, $
                                          send_Bx,send_By,send_Bz, $
                                          send_Jx,send_Jy,send_Jz, $
                                          kx,ky,kz, $
                                          kP,kPAngle, $
                                          ;; inds, $
                                          send_usedInds, $
                                          example, $
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
                                          PRE_VIII_LAYOUT=PRE_VIII_layout

  ENDIF

        out_inds      = TEMPORARY(inds   ) 
        out_freqs     = TEMPORARY(freq   )
        out_kx        = TEMPORARY(kx     )
        out_ky        = TEMPORARY(ky     )
        out_kz        = TEMPORARY(kz     )
        out_kP        = TEMPORARY(kP     )
        out_kPAngle   = TEMPORARY(kPAngle)

        out_avgJxBNrm = TEMPORARY(avgJxBNrm)
        slideMn_JxB   = MEAN(out_avgJxBNrm,DIMENSION=2)
        fluc          = [out_avgJxBNrm[0,*] - slideMn_JxB[0], $
                         out_avgJxBNrm[1,*] - slideMn_JxB[1], $
                         out_avgJxBNrm[2,*] - slideMn_JxB[2]]

        ;; out_inds      = (TEMPORARY(inds   ))[inds] 
        ;; out_freqs     = (TEMPORARY(freq   ))[inds]
        ;; out_kx        = (TEMPORARY(kx     ))[inds]
        ;; out_ky        = (TEMPORARY(ky     ))[inds]
        ;; out_kz        = (TEMPORARY(kz     ))[inds]
        ;; out_kP        = (TEMPORARY(kP     ))[inds]
        ;; out_kPAngle   = (TEMPORARY(kPAngle))[inds]

        out_usedInds  = TEMPORARY(usedInds)
        out_TArr      = TEMPORARY(TArr    )
        out_Bx        = TEMPORARY(Bx      )
        out_By        = TEMPORARY(By      )
        out_Bz        = TEMPORARY(Bz      )
        out_Jx        = TEMPORARY(Jx      )
        out_Jy        = TEMPORARY(Jy      )
        out_Jz        = TEMPORARY(Jz      )


END
