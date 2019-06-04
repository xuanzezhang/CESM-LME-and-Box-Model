;*************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
Startyr=850
Endyr=1849

year=ispan(Startyr,Endyr,1)
nyear=Endyr-Startyr+1

model="Box_Model"

runner=(/"001","002","003","004","005","006","007","008","009","010","011","012","013"/)
nrun=dimsizes(runner)


driver_string=(/" driven by NPP and varying ~F14~k~F21~~B1~e ~B1~c ~B1~o "," driven by NPP and varying ~F14~k~F21~~B1~e ~B1~c ~B1~o ",\
              " driven by NPP only"," driven by NPP only"," driven by varying ~F14~k~F21~~B1~e ~B1~c ~B1~o  only"," driven by varying ~F14~k~F21~~B1~e ~B1~c ~B1~o  only"/)
varname=(/"NBP","HR","NBP","HR","NBP","HR"/)
nvar=dimsizes(varname)

  PlotOutFileName="plot_spectrum_anal_of_Box_Model_CESM1_LME_850-1849_fluxes_V4_"
  wks = gsn_open_wks("pdf",PlotOutFileName+"NBP_HR_and_turnover")
  plot=new(nvar,graphic)
  bot_plot=new(nvar,graphic)
  Colors_line=(/"forestgreen","red","forestgreen","red","forestgreen","red"/)
  Colors_fill=(/"gray","pink2","gray","pink2","gray","pink2"/)
  YValue_Max=(/30,8,30,8,30,8/)
  YValue_Min=(/0,0,0,0,0,0/)

 do i= 0, nvar-1
   print(varname(i))

  if (i.le.3) then
  model="CESM1_CLM4"
  FilePath="/Volumes/XUANZEZHANG-RAID/xzhang/Last_Millennium_Carbon_Experiment/CESM1/Fluxes/" 
   if (i.le.3.and.i.ge.2) then
   model="Box_Model"
   FilePath="/Volumes/XUANZEZHANG-RAID/xzhang/Last_Millennium_Carbon_Experiment/CESM1/Box-Model-3-pools/" 
   end if
  else
  model="CESM1_CLM4_Box_Model"
  FilePath="/Volumes/XUANZEZHANG-RAID/xzhang/Last_Millennium_Carbon_Experiment/CESM1/Box-Model-3-pools/"
  end if  

  data=new((/nrun,nyear/),"float")

  do j= 0,nrun-1
  print(runner(j))
  FileName="Total_annual_"+varname(i)+"_"+model+"_"+runner(j)+"_from_"+Startyr+"-"+Endyr+".txt"
  data(j,:)=asciiread(FilePath+FileName,(/nyear/),"float")
  end do

  printVarSummary(data)
  data_avg = dim_avg_n(data,0)
  data_std = dim_stddev_n(data,0)

  data_up  = data_avg + data_std
  data_dw  = data_avg - data_std

;************************************************
; set function arguments
;************************************************
; detrending opt: 0=>remove mean 1=>remove mean and detrend
  d = 1
; smoothing periodogram: (0 <= sm <= ??.) should be at least 3 and odd
  sm = 5
; percent tapered: (0.0 <= pct <= 1.0) 0.10 common.
  pct = 0.010
;************************************************
; calculate spectrum
;************************************************

  spcx_data_mean = new((/500/),"float")
  spec_data_mean = specx_anal(data_avg(:),d,sm,pct)
 ; print(spec_data_mean@spcx)
  frq =spec_data_mean@frq
  spcx_data_mean(:) = spec_data_mean@spcx
  ;print(1./frq)


  spcx_data=new((/500,nrun/),"float")
  do j=0,nrun-1
  spec_data = specx_anal(data(j,:),d,sm,pct)
  frq =spec_data@frq
 ; print(spec_data@spcx)
  spcx_data(:,j)=spec_data@spcx
  delete(spec_data)
  end do

  spcx_data_avg=dim_avg_n(spcx_data,1)
  spcx_data_std=dim_stddev_n(spcx_data,1)

; spcx_data_std=runave(spcx_data_std,5,1)
  spcx_data_up=spcx_data_avg + spcx_data_std
  spcx_data_dw=spcx_data_avg - spcx_data_std
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
 res@tiMainFontHeightF  = 0.024                    ; change font heights
 res@tiXAxisFontHeightF = 0.024
 res@tiYAxisFontHeightF = 0.024
 
 res@xyDashPatterns    =(/0.,2./)
 res@xyLineThicknesses = (/2.,3./)
 res@xyLineColors      = (/Colors_line(i)/)
 res@xyExplicitLegendLabels = (/varname(i)+driver_string(i)/)
 if (varname(i).eq."HR") then
  res@xyExplicitLegendLabels = (/"RH"+driver_string(i)/)
 end if

; res@tiMainString           = varname(i)
 res@tiYAxisString          = "Variance ((Pg C)~S1~2  yr~S1~- ~S1~2  )" 
 res@tiXAxisString          = "Period (years)"  ; xaxis

 res@pmLegendParallelPosF   = 0.45                  ; move units right
 res@pmLegendOrthogonalPosF = -0.4                ; move units down
 if (i.ge.2)then
 res@pmLegendParallelPosF   = 0.36                  ; move units right
 res@pmLegendOrthogonalPosF = -0.4                ; move units down
 end if
 if (i.ge.4)then
 res@pmLegendParallelPosF   = 0.4                  ; move units right
 res@pmLegendOrthogonalPosF = -0.4                ; move units down
 end if
 res@pmLegendWidthF         = 0.08                  ; Change width and
 res@pmLegendHeightF        = 0.08                   ; height of legend.
 res@lgLabelFontHeightF     = 0.022
 
  res@trXMinF            = 1
; res@trXMaxF            = 0.5
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

  plot(i)  = gsn_csm_xy (wks,1./frq(::-1),(/spcx_data_avg(::-1)/),res) ; create plot

; Create a plot with the area between both curves filled in blue.
  delete(res@xyLineColors)
  res@gsnXYFillColors = Colors_fill(i)
  res@xyLineColor     = -1          ; We don't want the line, so make it transparent.
  res@pmLegendDisplayMode    = "Never"

  bot_plot(i)  = gsn_csm_xy (wks,1./frq(::-1),(/spcx_data_up(::-1),spcx_data_dw(::-1)/),res)  ; Create filled XY plot.

  overlay(bot_plot(i),plot(i))

;  draw(bot_plot)
;  frame(wks)

 delete(data)
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
 ;   resP@txString = "Comparison of CESM1-CLM4CN and Box-Model"
    resP@txFontHeightF = 0.022
   ;  resP@txPosYF  = 1.
    resP@txPosXF  = 0.52
    resP@amJust   = "TopLeft"
    resP@gsnPanelFigureStringsPerimOn = False
    resP@gsnPanelFigureStrings = (/"(a) CESM-LME","(b) CESM-LME","(c) Box-Model","(d) Box-Model","(e) CESM-LME minus Box-Model","(f) CESM-LME minus Box-Model"/)
    resP@gsnPanelFigureStringsBackgroundFillColor = -1
    resP@gsnPanelFigureStringsFontHeightF = 0.015
    resP@gsnPanelCenter        = False
    gsn_panel(wks,bot_plot(:),(/3,2/),resP)

end
