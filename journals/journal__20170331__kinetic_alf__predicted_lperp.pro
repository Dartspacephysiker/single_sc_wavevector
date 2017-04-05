;2017/03/31
PRO JOURNAL__20170331__KINETIC_ALF__PREDICTED_LPERP

  COMPILE_OPT IDL2,STRICTARRSUBS

  l_par    = 2D3 * 1000.D ;m
  ;; l_perp   = 1D  * 1000.D ;m
  n        = 30
  omega    = 1.3 * 2.D * !PI
  altitude = 3160               ;km
  mlt      = 12.5

  this = KINETIC_ALFVEN_DISP_RELATION(L_PAR=l_par, $
                                      L_PERP=l_perp, $
                                      N=n, $
                                      OMEGA_=omega, $
                                      ALTITUDE=altitude, $
                                      MLT=mlt)

END
