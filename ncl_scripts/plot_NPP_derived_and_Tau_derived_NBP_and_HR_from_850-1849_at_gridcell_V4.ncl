;*************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
;;;;;;;;;;;;;;;;
  fin= addfile("Timescalar_contribution_from_spectrum_anal_of_NBP_HR_Box_Model_ensemble_mean_from_850-1849_at_gridcell.nc","r")
  NPP_derived_NBP_10yr=fin->NPP_derived_NBP_10yr
  NPP_derived_NBP_100yr=fin->NPP_derived_NBP_100yr
  NPP_derived_NBP_1000yr=fin->NPP_derived_NBP_1000yr

  NPP_derived_HR_10yr=fin->NPP_derived_HR_10yr
  NPP_derived_HR_100yr=fin->NPP_derived_HR_100yr
  NPP_derived_HR_1000yr=fin->NPP_derived_HR_1000yr

  Tau_derived_NBP_10yr=fin->Tau_derived_NBP_10yr
  Tau_derived_NBP_100yr=fin->Tau_derived_NBP_100yr
  Tau_derived_NBP_1000yr=fin->Tau_derived_NBP_1000yr

  Tau_derived_HR_10yr=fin->Tau_derived_HR_10yr
  Tau_derived_HR_100yr=fin->Tau_derived_HR_100yr
  Tau_derived_HR_1000yr=fin->Tau_derived_HR_1000yr

;************************************************
  wks = gsn_open_wks("pdf","plot_NPP_and_Tau_derived_NBP_HR_from_850-1849_at_gridcell_V4") 
  gsn_define_colormap(wks,"wh-bl-gr-ye-re")

  res                      = True     ; plot mods desired
  res@gsnDraw            = False             ; don't draw yet
  res@gsnFrame           = False             ; don't advance frame yet
  res@vpWidthF              = 0.4
  res@vpHeightF             = 0.6
  res@cnFillOn             = True     ; turn on color fill
  res@cnLinesOn            = False    ; turn off contour lines
  res@cnFillOn         = True                 ; turn on color fill
  res@cnFillMode       = "RasterFill"         ; use raster fill      
;************************************************
; individual plots
;************************************************
  res@cnLevelSelectionMode = "ManualLevels"     ; manually set the contour levels 
  res@cnMinLevelValF  = 0.                   ; set the minimum contour level
  res@cnMaxLevelValF  = 1.                    ; set the maximum contour level
  res@cnLevelSpacingF         = 0.1                ; interval
  res@lbLabelStride           = 2                ; every 4th label
  res@lbLabelBarOn            = False
  res@mpLandFillColor         = "grey"          ; make land borders white
  res@mpMinLatF            = -60.
  res@tiMainFontHeightF    = 0.0125
  res@tmXBLabelFontHeightF  = 0.008
  res@tmYLLabelFontHeightF  = 0.008

  res@tiMainString            = "NBP driven by NPP only at 2-10yr timescales"
  plot1 = gsn_csm_contour_map(wks,NPP_derived_NBP_10yr,res)

  res@tiMainString            = "NBP driven by NPP only at 10-100yr timescales"
  plot2 = gsn_csm_contour_map(wks,NPP_derived_NBP_100yr,res)

  res@tiMainString            = "NBP driven by NPP only at 100-1000yr timescales"
  plot3 = gsn_csm_contour_map(wks,NPP_derived_NBP_1000yr,res)

  res@tiMainString            = "NBP driven by ~F14~k~F22~~B1~e ~B1~c ~B1~o  only at 2-10yr timescales"
  plot11 = gsn_csm_contour_map(wks,Tau_derived_NBP_10yr,res)

  res@tiMainString            = "NBP driven by ~F14~k~F22~~B1~e ~B1~c ~B1~o  only at 10-100yr timescales"
  plot22 = gsn_csm_contour_map(wks,Tau_derived_NBP_100yr,res)

  res@tiMainString            = "NBP driven by ~F14~k~F22~~B1~e ~B1~c ~B1~o  only at 100-1000yr timescales"
  plot33 = gsn_csm_contour_map(wks,Tau_derived_NBP_1000yr,res)

  res@tiMainString            = "RH driven by NPP only at 2-10yr timescales"
  plot4 = gsn_csm_contour_map(wks,NPP_derived_HR_10yr,res)

  res@tiMainString            = "RH driven by NPP only at 10-100yr timescales"
  plot5 = gsn_csm_contour_map(wks,NPP_derived_HR_100yr,res)

  res@tiMainString            = "RH driven by NPP only at 100-1000yr timescales"
  plot6 = gsn_csm_contour_map(wks,NPP_derived_HR_1000yr,res)

  res@tiMainString            = "RH driven by ~F14~k~F22~~B1~e ~B1~c ~B1~o  only at 2-10yr timescales"
  plot44 = gsn_csm_contour_map(wks,Tau_derived_HR_10yr,res)

  res@tiMainString            = "RH driven by ~F14~k~F22~~B1~e ~B1~c ~B1~o  only at 10-100yr timescales"
  plot55 = gsn_csm_contour_map(wks,Tau_derived_HR_100yr,res)

  res@tiMainString            = "RH driven by ~F14~k~F22~~B1~e ~B1~c ~B1~o  only at 100-1000yr timescales"
  plot66 = gsn_csm_contour_map(wks,Tau_derived_HR_1000yr,res)

;************************************************
; create panel
;************************************************
  resP  = True
    resP@gsnFrame          = False             ; don't advance frame yet
    resP@gsnDraw           = True            ; draw later, in maximize mod
    resP@gsnPanelLabelBar  = True              ; add common colorbar
    resP@lbOrientation     = "vertical"
    resP@lbLabelFontHeightF = 0.02
    resP@gsnPanelYWhiteSpacePercent = 5
    resP@gsnPanelXWhiteSpacePercent = 2
    resP@pmLabelBarWidthF  = 0.1             ; make thinner
    resP@pmLabelBarHeightF = 0.4
    resP@amJust   = "BottomLeft"
    resP@gsnPanelFigureStringsPerimOn = False
;    resP@gsnPanelFigureStrings = (/"(a) NBP driven by NPP only (2-10yr)","(b)NBP by Tau-derived NBP 10yr","(c) NPP-derived NBP 100yr",\
;                                    "(d) Tau-derived NBP 100yr","(e) NPP-derived NBP 1000yr","(f) Tau-derived NBP 1000yr"/)
    resP@gsnPanelFigureStrings = (/"(a)","(b)","(c)","(d)","(e)","(f)"/)
    resP@gsnPanelFigureStringsBackgroundFillColor = -1
    resP@gsnPanelFigureStringsFontHeightF = 0.02
    resP@gsnPanelCenter        = False
    gsn_panel(wks,(/plot1,plot11,plot2,plot22,plot3,plot33/),(/3,2/),resP)

    frame(wks)
;    resP@gsnPanelFigureStrings = (/"(a) NPP-derived HR 10yr","(b) Tau-derived HR 10yr","(c) NPP-derived HR 100yr",\
;                                    "(d) Tau-derived HR 100yr","(e) NPP-derived HR 1000yr","(f) Tau-derived HR 1000yr"/)

    gsn_panel(wks,(/plot4,plot44,plot5,plot55,plot6,plot66/),(/3,2/),resP)

end

