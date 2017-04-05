;2017/03/28
;JOURNAL__20170328__FIGURE_OUT_FAST_SPIN_AXIS,T1='1999-01-27/11:32:10',T2='1999-01-27/11:33:10'
PRO JOURNAL__20170328__FIGURE_OUT_FAST_SPIN_AXIS, $
   T1=t1, $
   T2=t2

  COMPILE_OPT IDL2,STRICTARRSUBS

  dlist = GET_DQDS(START_TIMES=start_times,END_TIMES=end_times)
  this  = WHERE(STRMATCH(STRUPCASE(dlist),STRUPCASE('*ees*')) OR $
                STRMATCH(STRUPCASE(dlist),STRUPCASE('*eeb*')) OR $
                STRMATCH(STRUPCASE(dlist),STRUPCASE('*ies*')) OR $
                STRMATCH(STRUPCASE(dlist),STRUPCASE('*ieb*'))    $
                ,nMatch)
  IF nMatch GT 0 THEN BEGIN
     t1BKUP = start_times[this[0]]
     t2BKUP = end_times[this[0]]
  ENDIF

  IF KEYWORD_SET(t1) THEN BEGIN
     CASE SIZE(t1,/TYPE) OF
        7: BEGIN
           t1ZoomStr   = t1
        END
        5: BEGIN
           t1ZoomStr   = TIME_TO_STR(t1,/MS)
        END
        ELSE: BEGIN
           PRINT,"Bogusness?"
           STOP
        END
     ENDCASE
  ENDIF ELSE BEGIN
     t1ZoomStr         = t1BKUP
  ENDELSE

  IF KEYWORD_SET(t2) THEN BEGIN
     CASE SIZE(t2,/TYPE) OF
        7: BEGIN
           t2ZoomStr   = t2
        END
        5: BEGIN
           t2ZoomStr   = TIME_TO_STR(t2,/MS)
        END
        ELSE: BEGIN
           PRINT,"Bogusness?"
           STOP
        END
     ENDCASE
  ENDIF ELSE BEGIN
     t2ZoomStr         = t2BKUP
  ENDELSE

  t1Zoom            = STR_TO_TIME(t1ZoomStr)
  t2Zoom            = STR_TO_TIME(t2ZoomStr)

  IF GET_SDT_TIMESPAN(t1SDT,t2SDT) THEN BEGIN
     PRINT,'SDT timespan is from ',TIME_TO_STR(t1SDT),' to ',TIME_TO_STR(t2SDT)
     IF (t1SDT GT t1Zoom) OR (t2SDT LT t2Zoom) THEN BEGIN
        PRINT,"Hosed it! You've got insufficient data."
        RETURN
     ENDIF
  ENDIF ELSE BEGIN
     PRINT,' Could not get timespan! No SDT available?'
     RETURN
  ENDELSE


  PLOT_FA_ATT,[t1Zoom,t2Zoom]

  STOP

END
