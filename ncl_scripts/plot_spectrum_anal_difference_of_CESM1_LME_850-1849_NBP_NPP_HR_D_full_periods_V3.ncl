;*************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
Startyr=850
Endyr=1849

year=ispan(Startyr,Endyr,1)
nyear=Endyr-Startyr+1

model_all="CESM1_CLM4"
model_850forcing="CESM1_CLM4_0850cntl"
model_GHG="CESM1_CLM4_GHG"
model_LULC="CESM1_CLM4_LULC_HurttPongratz"
model_ORBITAL="CESM1_CLM4_ORBITAL"
model_SSI="CESM1_CLM4_SSI_VSK_L"
model_VOLC="CESM1_CLM4_VOLC_GRA"

varname=(/"NBP","NPP","HR","D"/)
nvar=dimsizes(varname)

  PlotOutFileName="plot_spectrum_anal_difference_of_CESM1_LME_850-1849_fluxes_"
  wks = gsn_open_wks("pdf",PlotOutFileName+"NBP_NPP_HR_D_full_periods_V3")
  plot=new(nvar,graphic)
  bot_all_plot=new(nvar,graphic)
  bot_GHG_plot=new(nvar,graphic)
  bot_LULC_plot=new(nvar,graphic)
  bot_ORBITAL_plot=new(nvar,graphic)
  bot_SSI_plot=new(nvar,graphic)
  bot_VOLC_plot=new(nvar,graphic)

  Colors_line_850forcing=(/"black","black","black","black"/)
  Colors_fill_850forcing=(/"grey","grey","grey","grey"/)

  Colors_line_all=(/"black","black","black","black"/)
  Colors_fill_all=(/"grey","grey","grey","grey"/)

  Colors_line_GHG=(/"forestgreen","forestgreen","forestgreen","forestgreen"/)
  Colors_fill_GHG=(/"palegreen","palegreen","palegreen","palegreen"/)

  Colors_line_LULC=(/"blue","blue","blue","blue"/)
  Colors_fill_LULC=(/"lightblue","lightblue","lightblue","lightblue"/)

  Colors_line_ORBITAL=(/"gold","gold","gold","gold"/)
  Colors_fill_ORBITAL=(/"lightyellow","lightyellow","lightyellow","lightyellow"/)

  Colors_line_SSI=(/"purple","purple","purple","purple"/)
  Colors_fill_SSI=(/"thistle","thistle","thistle","thistle"/)

  Colors_line_VOLC=(/"red","red","red","red"/)
  Colors_fill_VOLC=(/"pink2","pink2","pink2","pink2"/)

  YValue_Max=(/30,30,8,4/)
  YValue_Min=(/0,0,0,0/)


 do i= 0, nvar-1
   print(varname(i))

FilePath="../Spectral_analysis_for_CESM-LME/"
 ; FilePath="/mnt/orico/xzhang/Last_Millennium_Carbon_Experiment/CESM1/Spectral_analysis_for_CESM-LME/" ; set to your own path

  fin_all=addfile(FilePath+"Spectrum_anal_of_"+varname(i)+"_"+model_all+"_from_"+Startyr+"-"+Endyr+".nc","r")
  spcx_data_all=fin_all->spcx_data
  frq          =fin_all->frq

  fin_850forcing=addfile(FilePath+"Spectrum_anal_of_"+varname(i)+"_"+model_850forcing+"_from_"+Startyr+"-"+Endyr+".nc","r")
  spcx_data_850forcing=fin_850forcing->spcx_data

  fin_GHG=addfile(FilePath+"Spectrum_anal_of_"+varname(i)+"_"+model_GHG+"_from_"+Startyr+"-"+Endyr+".nc","r")
  spcx_data_GHG=fin_GHG->spcx_data

  fin_LULC=addfile(FilePath+"Spectrum_anal_of_"+varname(i)+"_"+model_LULC+"_from_"+Startyr+"-"+Endyr+".nc","r")
  spcx_data_LULC=fin_LULC->spcx_data

  fin_ORBITAL=addfile(FilePath+"Spectrum_anal_of_"+varname(i)+"_"+model_ORBITAL+"_from_"+Startyr+"-"+Endyr+".nc","r")
  spcx_data_ORBITAL=fin_ORBITAL->spcx_data

  fin_SSI=addfile(FilePath+"Spectrum_anal_of_"+varname(i)+"_"+model_SSI+"_from_"+Startyr+"-"+Endyr+".nc","r")
  spcx_data_SSI=fin_SSI->spcx_data

  fin_VOLC=addfile(FilePath+"Spectrum_anal_of_"+varname(i)+"_"+model_VOLC+"_from_"+Startyr+"-"+Endyr+".nc","r")
  spcx_data_VOLC=fin_VOLC->spcx_data

  data_850forcing_avg = dim_avg_n(spcx_data_850forcing,1)
  data_850forcing_std = dim_stddev_n(spcx_data_850forcing,1)
  data_850forcing_avg = runave(data_850forcing_avg,5,1)
  data_850forcing_std = runave(data_850forcing_std,5,1)

  data_all_avg = dim_avg_n(spcx_data_all,1)
  data_all_std = dim_stddev_n(spcx_data_all,1)
  data_all_avg = runave(data_all_avg,5,1)
  data_all_std = runave(data_all_std,5,1)

  data_all_up  = data_all_avg + data_all_std
  data_all_dw  = data_all_avg - data_all_std

  data_GHG_avg = dim_avg_n(spcx_data_GHG,1)
  data_GHG_std = dim_stddev_n(spcx_data_GHG,1)
  data_GHG_avg = runave(data_GHG_avg,5,1)
  data_GHG_std = runave(data_GHG_std,5,1)

  data_GHG_up  = data_GHG_avg + data_GHG_std
  data_GHG_dw  = data_GHG_avg - data_GHG_std

  data_LULC_avg = dim_avg_n(spcx_data_LULC,1)
  data_LULC_std = dim_stddev_n(spcx_data_LULC,1)
  data_LULC_avg = runave(data_LULC_avg,5,1)
  data_LULC_std = runave(data_LULC_std,5,1)

  data_LULC_up  = data_LULC_avg + data_LULC_std
  data_LULC_dw  = data_LULC_avg - data_LULC_std

  data_ORBITAL_avg = dim_avg_n(spcx_data_ORBITAL,1)
  data_ORBITAL_std = dim_stddev_n(spcx_data_ORBITAL,1)
  data_ORBITAL_avg = runave(data_ORBITAL_avg,5,1)
  data_ORBITAL_std = runave(data_ORBITAL_std,5,1)

  data_ORBITAL_up  = data_ORBITAL_avg + data_ORBITAL_std
  data_ORBITAL_dw  = data_ORBITAL_avg - data_ORBITAL_std

  data_SSI_avg = dim_avg_n(spcx_data_SSI,1)
  data_SSI_std = dim_stddev_n(spcx_data_SSI,1)
  data_SSI_avg = runave(data_SSI_avg,5,1)
  data_SSI_std = runave(data_SSI_std,5,1)

  data_SSI_up  = data_SSI_avg + data_SSI_std
  data_SSI_dw  = data_SSI_avg - data_SSI_std

  data_VOLC_avg = dim_avg_n(spcx_data_VOLC,1)
  data_VOLC_std = dim_stddev_n(spcx_data_VOLC,1)
  data_VOLC_avg = runave(data_VOLC_avg,5,1)
  data_VOLC_std = runave(data_VOLC_std,5,1)

  data_VOLC_up  = data_VOLC_avg + data_VOLC_std
  data_VOLC_dw  = data_VOLC_avg - data_VOLC_std

;************************************************
; plotting
;************************************************
;  wks = gsn_open_wks("pdf",PlotOutFileName+varname(i))

 res =True

 res@gsnDraw              = False
 res@gsnFrame             = False

 res@gsnMaximize        = False
 res@vpKeepAspect       = True
 res@vpHeightF          = 0.4               ; change aspect ratio of plot
 res@vpWidthF           = 0.6

 res@pmLegendDisplayMode    = "Always"             ; turn on legend
 res@pmLegendSide           = "Top"                ; Change location of
 res@lgPerimOn = False
 res@tiMainFontHeightF  = 0.028                    ; change font heights
 res@tiXAxisFontHeightF = 0.024
 res@tiYAxisFontHeightF = 0.024

 res@gsnXYFillOpacities = 0.5
 
 res@xyDashPatterns    = (/1.,0.,1.,2.,3.,4.,5./)
 res@xyLineThicknesses = (/2.,2.,2.,2.,2.,2.,2./)
 
 res@xyLineColors      =  (/Colors_line_850forcing(i),Colors_line_all(i),Colors_line_GHG(i),Colors_line_LULC(i),\
                           Colors_line_ORBITAL(i),Colors_line_SSI(i),Colors_line_VOLC(i)/)

 res@xyExplicitLegendLabels = (/"Control",\
                                "Full forcings",\
                                "GHG only",\
				"Land-use only",\
				"Orbital only",\
				"Solar only",\
				"Volcanic only"/)

 res@tiMainString           = varname(i)
 if (varname(i).eq."HR") then 
 res@tiMainString           = "RH"
 end if

 res@tiYAxisString          = "Variance ((Pg C)~S1~2  yr~S1~- ~S1~2  )" 
 res@tiXAxisString          = "Period (years)"  ; xaxis

 res@pmLegendParallelPosF   = 0.45                  ; move units right
 res@pmLegendOrthogonalPosF = -0.7                ; move units down
 res@pmLegendWidthF         = 0.12                  ; Change width and
 res@pmLegendHeightF        = 0.20                   ; height of legend.
 res@lgLabelFontHeightF     = 0.018 
 
  res@trXMinF            = 1
  res@trXMaxF            = 1000
  print(YValue_Min(i))
  res@trYMinF            = YValue_Min(i)
  res@trYMaxF            = YValue_Max(i)
  res@trXLog             = True
 

  res@tmXBMajorLengthF = 0.02
  res@tmXBMinorLengthF = 0.012
  res@tmXBMajorOutwardLengthF = 0.02
  res@tmXBMinorOutwardLengthF = 0.012
  res@tmXBLabelFontHeightF = 0.022

  res@tmYLMajorLengthF = 0.02
  res@tmYLMinorLengthF = 0.012
  res@tmYLMajorOutwardLengthF = 0.02
  res@tmYLMinorOutwardLengthF = 0.012
  res@tmYLLabelFontHeightF = 0.022

  plot(i)  = gsn_csm_xy (wks,1./frq(::-1),(/data_850forcing_avg(::-1),data_all_avg(::-1),data_GHG_avg(::-1),data_LULC_avg(::-1),\
                                        data_ORBITAL_avg(::-1),data_SSI_avg(::-1),data_VOLC_avg(::-1)/),res) ; create plot

; Create a plot with the area between both curves filled in blue.
  res@pmLegendDisplayMode    = "Never" 
  delete(res@xyLineColors)

  res@gsnXYFillColors = (/Colors_fill_all(i)/)
  res@xyLineColor     = -1          ; We don't want the line, so make it transparent.
  bot_all_plot(i)  = gsn_csm_xy (wks,1./frq(::-1),(/data_all_up(::-1),data_all_dw(::-1)/),res)  ; Create filled XY plot.

  res@gsnXYFillColors = (/Colors_fill_GHG(i)/)
  res@xyLineColor     = -1          ; We don't want the line, so make it transparent.
  bot_GHG_plot(i)  = gsn_csm_xy (wks,1./frq(::-1),(/data_GHG_up(::-1),data_GHG_dw(::-1)/),res)  ; Create filled XY plot.

  res@gsnXYFillColors = (/Colors_fill_LULC(i)/)
  res@xyLineColor     = -1          ; We don't want the line, so make it transparent.
  bot_LULC_plot(i)  = gsn_csm_xy (wks,1./frq(::-1),(/data_LULC_up(::-1),data_LULC_dw(::-1)/),res)  ; Create filled XY plot.

  res@gsnXYFillColors = (/Colors_fill_ORBITAL(i)/)
  res@xyLineColor     = -1          ; We don't want the line, so make it transparent.
  bot_ORBITAL_plot(i)  = gsn_csm_xy (wks,1./frq(::-1),(/data_ORBITAL_up(::-1),data_ORBITAL_dw(::-1)/),res)  ; Create filled XY plot.

  res@gsnXYFillColors = (/Colors_fill_SSI(i)/)
  res@xyLineColor     = -1          ; We don't want the line, so make it transparent.
  bot_SSI_plot(i)  = gsn_csm_xy (wks,1./frq(::-1),(/data_SSI_up(::-1),data_SSI_dw(::-1)/),res)  ; Create filled XY plot.

  res@gsnXYFillColors = (/Colors_fill_VOLC(i)/)
  res@xyLineColor     = -1          ; We don't want the line, so make it transparent.
  bot_VOLC_plot(i)  = gsn_csm_xy (wks,1./frq(::-1),(/data_VOLC_up(::-1),data_VOLC_dw(::-1)/),res)  ; Create filled XY plot.

  overlay(bot_all_plot(i),bot_GHG_plot(i))
  overlay(bot_all_plot(i),bot_LULC_plot(i))
  overlay(bot_all_plot(i),bot_ORBITAL_plot(i))
  overlay(bot_all_plot(i),bot_SSI_plot(i))
  overlay(bot_all_plot(i),bot_VOLC_plot(i))
  overlay(bot_all_plot(i),plot(i))


;  draw(bot_plot)
;  frame(wks)

 delete(res@gsnXYFillColors)
end do

;************************************************
; create panel
;************************************************
  resP  = True
    resP@gsnFrame          = True             ; don't advance frame yet
    resP@gsnDraw           = True            ; draw later, in maximize mod
    resP@gsnPanelYWhiteSpacePercent   = 2
    resP@gsnPanelXWhiteSpacePercent   = 2
     resP@txString = "CESM1-CLM4CN Last Millennium Ensemble Project (850-1849)"
     resP@txFontHeightF = 0.022
   ;  resP@txPosYF  = 1.
     resP@txPosXF  = 0.55
    resP@amJust   = "TopLeft"
    resP@gsnPanelFigureStringsPerimOn = False
    resP@gsnPanelFigureStrings = (/"(a)","(b)","(c)","(d)"/)
    resP@gsnPanelFigureStringsBackgroundFillColor = -1
    resP@gsnPanelFigureStringsFontHeightF = 0.022
    resP@gsnPanelCenter        = False
    gsn_panel(wks,bot_all_plot(:),(/2,2/),resP)

end
