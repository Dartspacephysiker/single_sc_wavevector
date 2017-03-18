;2017/03/11
PRO JOURNAL__20170311__IMF_CONDITIONS_DURING_PRE_VIII_EVENTS, $
                                  LEEWAYLOWER=leewayLower, $
                                  LEEWAYUPPER=leewayUpper

  COMPILE_OPT IDL2,STRICTARRSUBS


  ;;For the first go-round, when I thought I was going to use orb 10832
  UTCRange = [['1999-01-27/11:32:51.5','1999-01-27/11:33:09.5'], $
             ['1999-05-18/06:50:58.630','1999-05-18/06:51:09.570']]

  ;;For the second go-round, with a reprise from orb 9585 
  UTCRange = [['1999-01-23/14:50:56','1999-01-23/14:51:07'], $
              ['1999-01-27/11:32:51.5','1999-01-27/11:33:09.5']]

  PRINT_IMF_CONDITIONS_UTCRANGE,UTCRange, $
                                LEEWAYLOWER=leewayLower, $
                                LEEWAYUPPER=leewayUpper
END
