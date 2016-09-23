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
  return,autocorr/T
END

;********************************************************************
; written by P. M. Bellan, April 2016, MS 128-95 Caltech, Pasadena CA 91125
; prefix gives file path for eps output, data input/output
; user should change to be appropriate for user's computer
PRO OUTPUT_SETUP,mode,prefix,example
  PREFIX='D:\CORR\MNSCRPTS\2016-JGR-spacecraft-current-wavevector\2106-Vinas-kvect-revised\'
  PRINT, 'prefix used for directory',prefix
  ;;
  filename = PREFIX+ $
           'JGR-single-space-craft-k-example'+STRTRIM(STRING(example),1)+'.eps'
  PRINT, 'print plot to file ',filename
  ;;SET UP FOR PLOTTING ON:
  ;; SCREEN, SET hardcopy = 0
  ;; PRINTER, SET hardcopy = 1
  ;; POSTSCRIPT FILE, SET hardcopy = 2
  hardcopy = mode
  IF hardcopy EQ 1 THEN SET_PLOT, 'printer' ;PRINTER
  IF hardcopy EQ 1 THEN DEVICE,YSIZE=25,YOFFSET=0
  IF hardcopy NE 1 THEN SET_PLOT,'X'
  IF hardcopy EQ 2 THEN BEGIN      ;;POSTSCRIPT FILE
     SET_PLOT, 'PS'
     !P.CHARTHICK = 3
     !P.THICK = 5
     DEVICE, FILE=filename,/ENCAPSUL,XSIZE=10,YSIZE=10,/INCHES,YOFFSET=2,/COLOR
  ENDIF
END

;********************************************************************
; written by P. M. Bellan, April 2016, MS 128-95 Caltech, Pasadena CA 91125
; place at end of code to do housekeeping
PRO CONCLUDE_OUTPUT,mode,prefix
  hardcopy = mode
  IF hardcopy EQ 1 THEN DEVICE,/CLOSE
  IF hardcopy EQ 1 THEN SET_PLOT,'X'
  IF hardcopy EQ 2 THEN DEVICE,/CLOSE
  IF hardcopy EQ 2 THEN SET_PLOT,'X'
  IF hardcopy EQ 1 THEN SET_PLOT,'printer'
  IF hardcopy EQ 1 THEN device,YSIZE=25,YOFFSET=5
  IF hardcopy ne 1 THEN SET_PLOT,'X'
  PRINT, 'FINISHED'
end

;********************************************************************
PRO ADD_NOISE,Bx,By,Bz,Jx,Jy,Jz,T
  ;;create and add noise IF example = 3

  Bsqtotal = 0.0
  Jsqtotal = 0.0

  ;;calculate rms amplitude of B and J
  FOR TT=0,T-1 do BEGIN         ; TT is time
     Bsqtotal = Bsqtotal+Bx[TT]^2+By[TT]^2+Bz[TT]^2
     Jsqtotal = Jsqtotal+Jx[TT]^2+Jy[TT]^2+Jz[TT]^2
  ENDFOR
  Brms= SQRT(Bsqtotal/T)
  Jrms= SQRT(Jsqtotal/T)
  ;; create noise from random number generator randomu
  Bnoise = FLTARR(T,3)
  & Jnoise = FLTARR(T,3)
  FOR i=0,T-1 do BEGIN
     Bnoise[i,0] = RANDOMU(seed)
     & Jnoise[i,0] = RANDOMU(seed)
     Bnoise[i,1] = RANDOMU(seed)
     & Jnoise[i,1] = RANDOMU(seed)
     Bnoise[i,2] = RANDOMU(seed)
     & Jnoise[i,2] = RANDOMU(seed)
  ENDFOR
  ;;renormalize noise to have unity rms
  Bnoisetot = 0 
  Jnoisetot = 0
  FOR TT=0,T-1 do BEGIN         ; TT is time
     Bnoisetot = Bnoisetot+ Bnoise[TT,0]*Bnoise[TT,0]+Bnoise[TT,1]*Bnoise[TT,1]+Bnoise[TT,2]*Bnoise[TT,2]
     Jnoisetot = Jnoisetot+ Jnoise[TT,0]*Jnoise[TT,0]+Jnoise[TT,1]*Jnoise[TT,1]+Jnoise[TT,2]*Jnoise[TT,2]
  ENDFOR
  Bnoiserms= SQRT(Bnoisetot/T)
  Jnoiserms= SQRT(Jnoisetot/T)

  ;; below sets noise to have unity rms
  Bnoise[*,0] = Bnoise[*,0]/Bnoiserms & Bnoise[*,1] = Bnoise[*,1]/Bnoiserms & Bnoise[*,2] = Bnoise[*,2]/Bnoiserms
  Jnoise[*,0] = Jnoise[*,0]/Jnoiserms & Jnoise[*,1] = Jnoise[*,1]/Jnoiserms & Jnoise[*,2] = Jnoise[*,2]/Jnoiserms
  ;;amount of noise set to 50% of rms signal
  Bnoisecoeff = 0.5*Brms
  Jnoisecoeff = 0.5*Jrms
  ;;add noise
  FOR i=0,T-1 do BEGIN          ; add noise contributions for each time
     Bx[i] = Bx[i] +Bnoisecoeff*Bnoise[i,0]
     By[i] = By[i] +Bnoisecoeff*Bnoise[i,1]
     Bz[i] = Bz[i] +Bnoisecoeff*Bnoise[i,2]
     Jx[i] = Jx[i] +Jnoisecoeff*Jnoise[i,0]
     Jy[i] = Jy[i] +Jnoisecoeff*Jnoise[i,1]
     Jz[i] = Jz[i] +Jnoisecoeff*Jnoise[i,2]
  ENDFOR                        ;
  RETURN
END

;********************************************************************
;MAIN PROGRAM
;********************************************************************
; written by P. M. Bellan, April 2016, MS 128-95 Caltech, Pasadena CA 91125
; this code is free to use, but please acknowledge P. M. Bellan if used

  ;;select output mode,
  ;;set output_mode =0 for screen, 1 for printer, 2 for postscript file
  output_mode = 0
  ;;select one of three examples
  ;; first example, a few discrete k components
  ;;example =1
  ;;
  ;; second example,
  ;; continuum of k components following some complicated function of omega
  ;;
  ;; third example,
  ;;same as second example, except noise added
  
  example = 1
  ;;setup file management, filenames for selected output mode
  OUTPUT_SETUP,output_mode,prefix,example
  ;;setup graphic layout
  columns = 3
  rows = 4
  !P.MULTI = [0, columns, rows, 0, 0]
  T = 200
  ;; number of elements in time domain vectors
  ;;setup frequency and wavevector arrays
  ;;
  ;; construct array of frequencies for time domain,
  ;; zero frequency will not be used
  omega = INDGEN(T/2+1)*2*!Pi/T
  ;; k_vector : first argument is for frequency,
  ;;second is for Cartesian coordinate,
  kvec = FLTARR(T/2+1,3)
  ;;setup k vectors for either example 1 or 2
  ;;first example
  IF example EQ 1 THEN BEGIN
     kx = FLTARR(T/2+1) & ky = kx & kz = kx
     kx[10]=-2
     ky[50] = 8
     kz[70] = 5
  ENDIF
  ;; second and third examples, continuum of k components
  ;;constructed using arbitrary complicated function

   IF example EQ 2 or example EQ 3 THEN BEGIN
      kmag = SQRT(omega)          ; pick some arbitrary function for magnitude for k

      ;;choose azimuthal spherical coordinate angle of unit k vector
      phi= 2*!pi*INDGEN(T/2+1)/T
      ;;choose polar spherical coordinate angle of unit k vector
      theta=!pi*INDGEN(T/2+1)/T
      ;; construct Cartesian unit vectors for k
      x_unit = SIN(theta)*COS(phi)
      y_unit = SIN(theta)*SIN(phi)
      z_unit = COS(phi)
      kx = kmag*x_unit
      ky = kmag*y_unit
      kz = kmag*z_unit
   endif
   ;;define components of k vector
   ;;plot k vectors that have been set up
   cs = 2
   IF example EQ 1 THEN ysize=10
   IF example EQ 2 or example EQ 3 THEN ysize =2
   PLOT, kx ,xtitle = $
         '!4x!3T/(2!4p!3)', ytitle='!8k!Dx!N',charsize=cs, yrange=[-ysize,ysize]
   PLOT, ky ,xtitle = $
         '!4x!3T/(2!4p!3)', ytitle='!8k!Dy!N',charsize=cs, yrange=[-ysize,ysize]
   PLOT, kz ,xtitle= $
         '!4x!3T/(2!4p!3)', ytitle='!8k!Dz!N',charsize=cs, yrange=[-ysize,ysize]
   ;;store k components in a vector for use later

   for i=0, T/2 do BEGIN
      kvec[i,0]= kx[i]
      kvec[i,1]= ky[i]
      kvec[i,2]= kz[i]
   ENDFOR

   ;; set up k vector
   ;; define A vector potential cosine component
   ;; first argument is frequency, second is Cartesian component
   Avec_cos = FLTARR(T/2+1,3)     ; define A vector cos component (vector potential)
   Avec_sin = FLTARR(T/2+1,3)     ; define A vector sine component (vector potential)

   ;; seed for random number generator,
   ;; get different random number each time randomu is called
   seed = 2001L
   ;; make coefficients of vector potential random,
   ;;Avec_cos[i,0] means
   ;; coefficient of the x component with cosine time dependence for frequency i
   ;;make all vector potential components random numbers
   for i=0,T/2 do BEGIN
      Avec_cos[i,0] = RANDOMU(seed)
      & Avec_sin[i,0] = RANDOMU(seed)
      Avec_cos[i,1] = RANDOMU(seed)
      & Avec_sin[i,1] = RANDOMU(seed)
      Avec_cos[i,2] = RANDOMU(seed)
      & Avec_sin[i,2] = RANDOMU(seed)
   ENDFOR
   ;; B vector for a single frequency, args are frequency, Cartesian direction, time
   Bvec = FLTARR(T/2+1,3,T)
   Jvec = Bvec
   for TT=0,T-1 do BEGIN
      for i=0,T/2 do BEGIN
         ;;for each time
         ;;for each frequency component
         ;; temp A cosine component
         Avectemp_cos= [Avec_cos[i,0], Avec_cos[i,1], Avec_cos[i,2] ]
         ;; temp A sin e component
         Avectemp_sin= [Avec_sin[i,0], Avec_sin[i,1], Avec_sin[i,2] ]
         ;; tempo k vector
         kvectemp= [kvec[i,0], kvec[i,1], kvec[i,2] ]
         ;; k vector cross product with A cosine
         kcross_Acos = CROSSP(kvectemp,Avectemp_cos)
         ;; k vector cross product with A sine
         kcross_Asin = CROSSP(kvectemp,Avectemp_sin)
         ;; put in time dependence to get total B at a single time for freq i
         Bvec[i,0,tt] = -kcross_Acos[0]*SIN(-omega[i]*tt) $
                        + kcross_Asin[0]*COS(-omega[i]*tt)

         Bvec[i,1,tt]= -kcross_Acos[1]*SIN(-omega[i]*tt)$
                       + kcross_Asin[1]*COS(-omega[i]*tt)
         Bvec[i,2,tt]= -kcross_Acos[2]*SIN(-omega[i]*tt) $
                       + kcross_Asin[2]*COS(-omega[i]*tt)
         ;; k cross (k cross A) need for current
         kcrosskcross_Acos= CROSSP(kvectemp,kcross_Acos)
         kcrosskcross_Asin= CROSSP(kvectemp,kcross_Asin)
         ;;calculate current at a single time for freq i
         Jvec[i,0,tt]= -kcrosskcross_Acos[0]* $
                       COS(-omega[i]*tt)-kcrosskcross_Asin[0]*SIN(-omega[i]*tt)
         Jvec[i,1,tt]= -kcrosskcross_Acos[1]* $
                       COS(-omega[i]*tt)-kcrosskcross_Asin[1]*SIN(-omega[i]*tt)
         Jvec[i,2,tt]= -kcrosskcross_Acos[2]* $
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
   For TT=0,T-1 do BEGIN        ; TT is time
      for i=0,T/2 do BEGIN      ; sum up contributions from each frequency
         Bx[TT] = Bx[TT] +Bvec[i,0,TT]
         By[TT] = By[TT] +Bvec[i,1,TT]
         Bz[TT] = Bz[TT] +Bvec[i,2,TT]
         Jx[TT] = Jx[TT] +Jvec[i,0,TT]
         Jy[TT] = Jy[TT] +Jvec[i,1,TT]

         Jz[TT] = Jz[TT] +Jvec[i,2,TT]
      ENDFOR                    ;
   ENDFOR
   ;;create noise IF example = 3
   IF example EQ 3 THEN add_noise, Bx,By,Bz,Jx,Jy,Jz,T
   ;;Write data to file
   ;; so that code can be easily modified to input data from other sources
   filename = prefix+'Bx_By_Bz_Jx_Jy_Jz.dat'
   ;;write data to file
   OPENW, 1, filename & PRINT, 'opening ', filename
   printf, 1,example, T,FORMAT='(I5,1x,I5,1x)'
   printf, 1, Bx
   printf, 1, By
   printf, 1, Bz
   printf, 1, Jx
   printf, 1, Jy
   printf, 1, Jz
   CLOSE, 1 & PRINT, 'closing ', filename
   ;;read data from file
   filename = prefix+'Bx_By_Bz_Jx_Jy_Jz.dat'
   OPENr, 1, filename & PRINT, 'opening ', filename
   readf, 1, example, T_read,FORMAT='(I5,1x,I5,1x)'
   Bx = FLTARR(T_read) & By = Bx & Bz = Bx &Jx = Bx &Jy =Bx & Jz = Bx
   readf,1,Bx
   readf,1,By
   readf,1,Bz
   readf,1,Jx
   readf,1,Jy
   readf,1,Jz
   CLOSE, 1 & PRINT, 'closing ', filename

   PLOT,Bx,XTITLE='t',YTITLE='!8B!Dx!N',CHARSIZE=cs
   PLOT,By,XTITLE='t',YTITLE='!8B!Dy!N',CHARSIZE=cs
   PLOT,Bz,XTITLE='t',YTITLE='!8B!Dz!N',CHARSIZE=cs
   
   PLOT, Jx,xtitle = 't', ytitle ='!4l!3!D0 !N!8J!Dx!N',charsize=cs
   PLOT, Jy,xtitle = 't', ytitle ='!4l!3!D0 !N!8J!Dy!N',charsize=cs
   PLOT, Jz,xtitle = 't', ytitle ='!4l!3!D0 !N!8J!Dz!N',charsize=cs
   
   
   JxBtotal = 0
   norm = 0                     ; check that avg J x B =0
   For TT=0,T-1 DO BEGIN        ; integrate over time
      Jvectemp = [Jx[TT], Jy[TT], Jz[TT]]
      Bvectemp = [Bx[TT], By[TT], Bz[TT]]
      JxBtotal = JxBtotal+CROSSP(Jvectemp,Bvectemp)

      ;; norm for denominator
      norm= norm + SQRT(Jx[TT]^2+ Jy[TT]^2+ Jz[TT]^2)*SQRT(Bx[TT]^2+ By[TT]^2+ Bz[TT]^2)
   ENDFOR
   PRINT,'norm=',norm
   avgJxBtotal = JxBtotal/T
   PRINT, 'avgJxBtotal/norm=',avgJxBtotal/norm ; small (supposed to be zero)

   ;; auto-correlation of B, i.e, <B_VEC(t) dot B_VEC(t+tau)>
   BBautocorr = CORRELATION(Bx,Bx,T) + CORRELATION(By,By,T) + CORRELATION(Bz,Bz,T)

   ;;components of cross correlation <J_VEC(t) cross B_vec (t+tau)>
   JxB_xcomponent_corr = CORRELATION(Jy,Bz,T) - CORRELATION(Jz,By,T)

   JxB_ycomponent_corr = CORRELATION(Jz,Bx,T) - CORRELATION(Jx,Bz,T)
   JxB_zcomponent_corr = CORRELATION(Jx,By,T) - CORRELATION(Jy,Bx,T)

   ;; Fourier transform of <B_vec(t) dot B_vec(t+tau)>
   BB_FFT = FFT(BBautocorr)

   ;; Fourier transform of components of <B_vec(t) dot B_vec(t+tau)>
   JxB_xcomponent_FFT = FFT(JxB_xcomponent_corr)
   JxB_ycomponent_FFT = FFT(JxB_ycomponent_corr)
   JxB_zcomponent_FFT = FFT(JxB_zcomponent_corr)

   ;; calculate k components, put 0.001 in denom to avoid dividing zero by zero
   ikx=-JxB_xcomponent_FFT/(BB_FFT+.001)
   iky=-JxB_ycomponent_FFT/(BB_FFT+.001)
   ikz=-JxB_zcomponent_FFT/(BB_FFT+.001)
   kx = IMAGINARY(ikx)
   ky = IMAGINARY(iky)
   kz = IMAGINARY(ikz)
   ;;extract imaginary part

   ;;plot components of calculated k vectors
   IF example EQ 1 THEN ysize = 10
   IF example EQ 2 or example EQ 3 THEN ysize = 2
   PLOT, kx , ytitle='x component', xtitle='FFT argument',charsize=cs, yrange=[-ysize,ysize]
   IF example EQ 2 THEN oPLOT, kx+0.1, linestyle=2 ;dashed line
   PLOT, ky , ytitle='y component', xtitle='FFT argument',charsize=cs, yrange=[-ysize,ysize]
   IF example EQ 2 THEN oPLOT, ky+0.1, linestyle=2
   PLOT, kz , ytitle='z component', xtitle='FFT argument',charsize=cs, yrange=[-ysize,ysize]
   IF example EQ 2 THEN oPLOT, kz+0.1, linestyle=2

   ;; print to screen, printer, or file depending on mode=0,1,2
   CONCLUDE_OUTPUT,output_mode
end