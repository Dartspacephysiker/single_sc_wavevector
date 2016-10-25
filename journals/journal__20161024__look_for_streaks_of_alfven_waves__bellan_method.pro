;;10/24/16
PRO JOURNAL__20161024__LOOK_FOR_STREAKS_OF_ALFVEN_WAVES__BELLAN_METHOD

  COMPILE_OPT IDL2

  dbDir      = '/home/spencerh/Research/database/FAST/dartdb/saves/'

  ;; dbFile     = 'Dartdb_20151222--500-16361_inc_lower_lats--burst_1000-16361--w_Lshell--correct_pFlux--maximus.sav'

  dbFile     = 'Dartdb_20150810--1000-16361--maximus--burstmode.sav'
  dbTimeFile = 'Dartdb_20150810--1000-16361--cdbtime--burstmode.sav'
  

  RESTORE,dbDir+dbFile
  RESTORE,dbDir+dbTimeFile

  tagNames = TAG_NAMES(maximus)
  ;;Sort them
  ;; sort_i     = SORT(cdbTime)
  ;; cdbTime    = cdbTime[sort_i]


  uniq_i     = WHERE(ABS(maximus.esa_current) GE 10.)
  cdbTime    = cdbTime[uniq_i]

  maxTemp = CREATE_STRUCT(tagNames[0],(maximus.(0))[uniq_i])
  FOR k=1, N_ELEMENTS(tagNames)-1 DO BEGIN
     maxTemp = CREATE_STRUCT(maxTemp,tagNames[k],(maximus.(k))[uniq_i])
  ENDFOR

  maximus = TEMPORARY(maxTemp)


  uniq_i     = UNIQ(cdbTime,SORT(cdbTime))
  cdbTime    = cdbTime[uniq_i]

  maximus2 = CREATE_STRUCT(tagNames[0],(maximus.(0))[uniq_i])
  FOR k=1, N_ELEMENTS(tagNames)-1 DO BEGIN
     maximus2 = CREATE_STRUCT(maximus2,tagNames[k],(maximus.(k))[uniq_i])
  ENDFOR

  ;; final_i    = CGSETINTERSECTION(WHERE(cdbTime LE STR_TO_TIME('1999-11-12/22:52:14.217')), $
  ;;                                uniq_i)

  final_i    = WHERE(cdbTime LE STR_TO_TIME('1999-11-12/22:52:14.217'))
  cdbTime    = cdbTime[final_i]

  maximus = CREATE_STRUCT(tagNames[0],(maximus2.(0))[final_i])
  FOR k=1, N_ELEMENTS(tagNames)-1 DO BEGIN
     maximus = CREATE_STRUCT(maximus,tagNames[k],(maximus2.(k))[final_i])
  ENDFOR

  maximus2 = !NULL

  ;; inds = 

  ;;DOUBLE_STREAKS options
  decimal_place = -1.
  N             = 5
  gap_time      = 1.5
  print_times   = 1
  GET_DOUBLE_STREAKS__NTH_DECIMAL_PLACE,cdbTime,decimal_place, $
                                        MAXIMUS=maximus, $
                                        NPTS=n, $
                                        MIN_T_STREAKLEN=min_streakLen_t, $
                                        GAP_TIME=gap_time, $
                                        START_I=start_i, $
                                        STOP_I=stop_i, $
                                        STREAKLENS=streakLens, $
                                        T_STREAKLENS=streakLens_t, $
                                        FLOOR=floor, $
                                        CEILING=ceiling, $
                                        PRINT_START_STOP_TIMES=print_times, $
                                        /PRINT__INCLUDE_CURRENT
                                        

  ;; big_ii       = REVERSE(SORT(streakLens))
  ;; streakLens   = streakLens[big_ii]
  ;; start_i      = start_i[big_ii]
  ;; stop_i       = stop_i[big_ii]
  ;; streakLens_t = streakLens_t[big_ii]

  checkem      = [start_i[1]:stop_i[1]]

  PRINT,maximus.time[checkem]

  STOP


END
