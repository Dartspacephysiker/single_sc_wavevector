;2017/03/31
PRO JOURNAL__20170331__READ_JIM_SPECTRUM_FILES_FOR_PRE_VIII

  COMPILE_OPT IDL2,STRICTARRSUBS

  dir         = '/SPENCEdata/Research/Satellites/FAST/single_sc_wavevector/saves_output_etc/20170331/PRE_VIII_data_from_Jim/'
  mag_files   = ['orb_9585','orb_9627']+ '-j-mag-spectrum'
  esa_files   = ['orb_9585','orb_9627']+ '-j-esa-spectrum'
  tmpltFile   = 'PRE_VII_Jimfiles_ASCII_tmplt.sav'
  outFileSuff = '-PRE_VIII_Jimfiles.sav'

  IF FILE_TEST(dir+tmpltFile) THEN BEGIN

     RESTORE,dir+tmpltFile

  ENDIF ELSE BEGIN

     tmplt   = ASCII_TEMPLATE(dir+mag_files[0])
     SAVE,tmplt,FILENAME=dir+tmpltFile     

  ENDELSE

  ;; orb_list      = LIST()
  ;; magF_list     = LIST()
  ;; mag_DB_list   = LIST()
  ;; esaF_list     = LIST()
  ;; esa_DB_list   = LIST()
  orbArr = !NULL
  magArr = !NULL
  esaArr = !NULL
  FOR k=0,N_ELEMENTS(mag_files)-1 DO BEGIN

     orbArr = [orbArr,STRMID(mag_files[k],4,4)]

     mag        = READ_ASCII(dir+mag_files[k],TEMPLATE=tmplt)
     esa        = READ_ASCII(dir+esa_files[k],TEMPLATE=tmplt)

     esaArr     = [esaArr,esa]
     magArr     = [magArr,mag]

  ENDFOR

  FOR k=0,N_ELEMENTS(orbArr)-1 DO BEGIN

     tmpOutFile = orbArr[k] + outFileSuff
     PRINT,"Saving " + tmpOutFile

     mag        = magArr[k]
     esa        = esaArr[k]

     SAVE,mag,esa,FILENAME=dir+tmpOutFile

  ENDFOR

  STOP

END
