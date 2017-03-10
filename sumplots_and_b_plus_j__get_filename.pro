;;2017/02/24
PRO SUMPLOTS_AND_B_PLUS_J__GET_FILENAME, $
   ORBIT=orbit, $
   DATE=date, $
   EEB_OR_EES=eeb_or_ees, $
   IEB_OR_IES=ieb_or_ies, $
   SKIP_DESPIN=skip_despin, $
   PLOTPREF=plotPref, $
   SAVESUFF=saveSuff, $
   BONUSSUFF=bonusSuff, $
   USE_REPRETCAL_FILE=use_RepRetCal_file,$
   ANCILLARY_PLOTS=ancillary_plots, $
   OUTPLOTNAME=outPlotName, $
   SAVEFILE=saveFile
   
  COMPILE_OPT IDL2,STRICTARRSUBS

  IF ~KEYWORD_SET(eeb_or_ees) THEN BEGIN
     eeb_or_ees     = 'eeb'
  ENDIF
  IF ~KEYWORD_SET(ieb_or_ies) THEN BEGIN
     ieb_or_ies     = 'ieb'
  ENDIF

  despunStr         = KEYWORD_SET(skip_despin) ? '-no_B_despin' : ''

  dateStr           = KEYWORD_SET(date       ) ? date           : GET_TODAY_STRING(/DO_YYYYMMDD_FMT)

  IF ~KEYWORD_SET(plotPref) THEN plotPref = ''
  IF ~KEYWORD_SET(saveSuff) THEN saveSuff = ''

  orbStr            = STRCOMPRESS(orbit,/REMOVE_ALL)
  ee_ie_string      = '-' + eeb_or_ees + '-' + ieb_or_ies

  outPlotName       = 'Orb_' + orbStr + plotPref + ee_ie_string + despunStr
  saveFile          = 'Orbit_' + orbStr + '-B_and_J-' + $
                      dateStr + $
                      ee_ie_string + despunStr + $
                      saveSuff + '.sav'


  IF N_ELEMENTS(ancillary_plots) EQ 0 THEN ancillary_plots = 1
  IF KEYWORD_SET(ancillary_plots) THEN BEGIN
     outPlotName   += '-with_ancillaries'
  ENDIF

  ;;Alternative
  IF KEYWORD_SET(bonusSuff) THEN BEGIN
     outPlotName      += bonusSuff
     saveFile         += bonusSuff
  ENDIF

  IF KEYWORD_SET(use_RepRetCal_file) THEN BEGIN
     tmpSuff        = '-' + 'RRC'
     saveFile       = saveFile.Replace('.sav',tmpSuff + '.sav')
     outPlotName   += tmpSuff
  ENDIF
  
END
