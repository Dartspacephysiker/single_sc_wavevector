;2017/03/11
PRO JOURNAL__20170311__IMF_CONDITIONS_DURING_PRE_VIII_EVENTS, $
                                  LEEWAYLOWER=leewayLower, $
                                  LEEWAYUPPER=leewayUpper

  COMPILE_OPT IDL2,STRICTARRSUBS


  UTCRange = [['1999-01-27/11:32:51.5','1999-01-27/11:33:09.5'], $
             ['1999-05-18/06:50:58.630','1999-05-18/06:51:09.570']]

  PRINT_IMF_CONDITIONS_UTCRANGE,UTCRange, $
                                LEEWAYLOWER=leewayLower, $
                                LEEWAYUPPER=leewayUpper
END
