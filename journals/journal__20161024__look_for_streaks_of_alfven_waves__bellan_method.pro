;;2016/10/24
;;Wanna see where the rubber meets the road? Check out the BATCH_JOBS repo under
;;"reproducing_figures_from_the_literature/Chaston_et_al_2006--erosion--wavevector/"
;;That's where you got it.
PRO JOURNAL__20161024__LOOK_FOR_STREAKS_OF_ALFVEN_WAVES__BELLAN_METHOD

  COMPILE_OPT IDL2

  ;;Options for current threshold
  curThresh       = 10.
  ;; curThresh       = 2.
  map_current     = 0B

  ;;Options for streaks
  min_streakLen_t = 10  ;min temporal length of streak
  decimal_place   = -1. ;don't mess, it's good
  N               = 5   ;min # points to qualify a streak
  gap_time        = 2.0 ;max allowable gap between observations

  ;;how to sort?
  sort_reverse    = 1B ;Put big guys up front

  streakLen_tSort = 0B
  streakLenSort   = 0B
  sort_by_avg_dt  = 0B
  jESA_sort       = 0B
  jMAG_sort       = 0B
  ABS_jESA_sort   = 0B
  ABS_jMAG_sort   = 0B

  dbDir           = '/home/spencerh/Research/database/FAST/dartdb/saves/'

  ;; dbFile     = 'Dartdb_20151222--500-16361_inc_lower_lats--burst_1000-16361--w_Lshell--correct_pFlux--maximus.sav'

  dbFile          = 'Dartdb_20150810--1000-16361--maximus--burstmode.sav'
  dbTimeFile      = 'Dartdb_20150810--1000-16361--cdbtime--burstmode.sav'
  dbMapFile       = 'Dartdb_20170223--1000-16361--mapRatio--burstmode.sav'
  get_mapRatio    = ~FILE_TEST(dbDir+dbMapFile)

  RESTORE,dbDir+dbFile
  RESTORE,dbDir+dbTimeFile

  tagNames       = TAG_NAMES(maximus)
  goingRate      = N_ELEMENTS(maximus.(0))

  IF goingRate NE N_ELEMENTS(cdbTime) THEN STOP

  ;;First screen for times that are actually unique
  uniq_i         = UNIQ(cdbTime,SORT(cdbTime))

  maxTemp        = CREATE_STRUCT(tagNames[0],(maximus.(0))[uniq_i])
  FOR k=1, N_ELEMENTS(tagNames)-1 DO BEGIN
     IF N_ELEMENTS(maximus.(k)) EQ goingRate THEN BEGIN
        maxTemp  = CREATE_STRUCT(maxTemp,tagNames[k],(maximus.(k))[uniq_i])
     ENDIF ELSE BEGIN
        PRINT,"Bogus: ",tagNames[k]
        STOP
     ENDELSE
  ENDFOR

  ;;Resize
  cdbTime        = cdbTime[uniq_i]
  maximus        = TEMPORARY(maxTemp)
  goingRate      = N_ELEMENTS(maximus.(0))

  IF get_mapRatio THEN BEGIN
     GET_ALT_MLT_ILAT_FROM_FAST_EPHEM,!NULL,cdbTime, $
                                      OUT_TSORTED_I=tSort_i, $
                                      OUT_ALT=alt, $
                                      OUT_MLT=mlt, $
                                      OUT_ILAT=ilat, $
                                      OUT_MAPRATIO=mapRatio, $
                                      OUT_NEVENTS=nEvents, $
                                      LOGLUN=logLun

     ;;If we have tSort_i, cdbTime had less-than-100% unique entries
     IF N_ELEMENTS(tSort_i) GT 0 THEN STOP

     SAVE,mapRatio,FILENAME=dbDir+dbMapFile
  ENDIF ELSE BEGIN
     RESTORE,dbDir+dbMapFile
  ENDELSE

  IF (goingRate NE N_ELEMENTS(cdbTime))            OR $
     (N_ELEMENTS(mapRatio) NE N_ELEMENTS(cdbTime))    $
  THEN STOP

  FASTDB__ADD_INFO_STRUCT,maximus, $
                          /FOR_ALFDB, $
                          DB_DIR=DBDir, $
                          DB_DATE=DB_date, $
                          DB_VERSION=DB_version, $
                          DB_EXTRAS=DB_extras, $
                          DB__INTO_ESPEC_FILE=DB__into_eSpec_file

  IF KEYWORD_SET(map_current) THEN BEGIN
     CORRECT_ALFVENDB_FLUXES,maximus, $
                             /MAP_ESA_CURRENT_TO_IONOS, $
                             /MAP_MAG_CURRENT_TO_IONOS, $
                             MAPRATIO=mapRatio
  ENDIF
  
  this              = WHERE(STRUPCASE(tagNames) EQ 'MODE')
  IF this[0] EQ -1 THEN STOP
  clean_these_inds = INDGEN(this)
  tagNames       = TAG_NAMES(maximus)
  
  ;;Apply current restriction
  clean_i        = BASIC_DB_CLEANER(maximus, $
                                    /CLEAN_NANS_AND_INFINITIES, $
                                    CLEAN_THESE_INDS=clean_these_inds, $
                                    /DISREGARD_SAMPLE_T)
  clean_i        = CGSETINTERSECTION(WHERE(ABS(maximus.esa_current) GE curThresh), $
                                     clean_i, $
                                     COUNT=nClean, $
                                     NORESULT=-1)
  IF nClean EQ 0 THEN STOP
  
  maxTemp        = CREATE_STRUCT(tagNames[0],(maximus.(0))[clean_i])
  FOR k=1, N_ELEMENTS(tagNames)-1 DO BEGIN
     IF N_ELEMENTS(maximus.(k)) EQ goingRate THEN BEGIN
        maxTemp  = CREATE_STRUCT(maxTemp,tagNames[k],(maximus.(k))[clean_i])
     ENDIF ELSE BEGIN
        IF STRUPCASE(tagNames[k]) EQ 'INFO' THEN BEGIN
           maxTemp = CREATE_STRUCT(maxTemp,tagNames[k],maximus.(k))
        ENDIF ELSE BEGIN
           PRINT,"Bogus: ",tagNames[k]
           STOP
        ENDELSE
     ENDELSE
  ENDFOR

  ;;Resize
  cdbTime        = cdbTime[clean_i]
  maximus        = TEMPORARY(maxTemp)
  goingRate      = N_ELEMENTS(maximus.(0))

  ;;Now cull the garbage
  final_i        = WHERE(cdbTime LE STR_TO_TIME('1999-11-12/22:52:14.217'))

  maxTemp        = CREATE_STRUCT(tagNames[0],(maximus.(0))[final_i])
  FOR k=1, N_ELEMENTS(tagNames)-1 DO BEGIN
     IF N_ELEMENTS(maximus.(k)) EQ goingRate THEN BEGIN
        maxTemp  = CREATE_STRUCT(maxTemp,tagNames[k],(maximus.(k))[final_i])
     ENDIF ELSE BEGIN
        IF STRUPCASE(tagNames[k]) EQ 'INFO' THEN BEGIN
           maxTemp = CREATE_STRUCT(maxTemp,tagNames[k],maximus.(k))
        ENDIF ELSE BEGIN
           PRINT,"Bogus: ",tagNames[k]
           STOP
        ENDELSE
     ENDELSE
  ENDFOR

  ;;Resize
  cdbTime        = cdbTime[final_i]
  maximus        = TEMPORARY(maxTemp)

  ;;DOUBLE_STREAKS options
  print_times     = 1
  GET_DOUBLE_STREAKS__NTH_DECIMAL_PLACE,cdbTime,decimal_place, $
                                        MAXIMUS=maximus, $
                                        NPTS=n, $
                                        MIN_T_STREAKLEN=min_streakLen_t, $
                                        GAP_TIME=gap_time, $
                                        START_I=start_i, $
                                        STOP_I=stop_i, $
                                        STREAKLENS=streakLens, $
                                        T_STREAKLENS=streakLens_t, $
                                        NSTREAKS=nStreaks, $
                                        FLOOR=floor, $
                                        CEILING=ceiling, $
                                        PRINT_START_STOP_TIMES=print_times, $
                                        /PRINT_MAXIMUS__INCLUDE_CURRENT, $
                                        SORT_BY_STREAKLEN=streakLenSort, $
                                        SORT_BY_T_STREAKLEN=streakLen_tSort, $
                                        SORT_BY_AVG_DT=sort_by_avg_dt, $
                                        SORT_BY_MAGNITUDE_ESA_CURRENT=ABS_jESA_sort, $
                                        SORT_BY_MAGNITUDE_MAG_CURRENT=ABS_jMAG_sort, $
                                        SORT_BY_ESA_CURRENT=jESA_sort, $
                                        SORT_BY_MAG_CURRENT=jMAG_sort, $
                                        SORT_REVERSE=sort_reverse, $
                                        NO_SORT=no_sort


  ;; big_ii       = REVERSE(SORT(streakLens))
  ;; streakLens   = streakLens[big_ii]
  ;; start_i      = start_i[big_ii]
  ;; stop_i       = stop_i[big_ii]
  ;; streakLens_t = streakLens_t[big_ii]

  checkem      = [start_i[1]:stop_i[1]]

  PRINT,maximus.time[checkem]

  STOP


END
