;*************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
Startyr1=850
Endyr1=1849

Startyr2=1850
Endyr2=2005

year=ispan(Startyr1,Endyr2,1)

nyear1=Endyr1-Startyr1+1
nyear2=Endyr2-Startyr2+1

model="CESM1_CLM4"

runner=(/"001","002","003","004","005","006","007","008","009","010","011","012","013"/)
nrun=dimsizes(runner)

varname=(/"NBP","NPP","HR","D"/)
nvar=dimsizes(varname)

  PlotOutFileName="plot_timeseries_of_CESM1_LME_850-2005_fluxes_"
  wks = gsn_open_wks("pdf",PlotOutFileName+"NBP_NPP_HR_D_V4")
  plot=new(nvar,graphic)
  bot_plot=new(nvar,graphic)
  dum=new(nvar,graphic)

  Colors_line=(/"forestgreen","blue","red","black"/)
  Colors_fill=(/"gray","lightblue","pink2","grey"/)
  YValue_Max=(/6,64,60,7/)
  YValue_Min=(/-6,52,50,2/)

 do i= 0, nvar-1
  print(varname(i))
  data1=new((/nrun,nyear1/),"float")
  data2=new((/nrun,nyear2/),"float")
;  FilePath="/mnt/orico/xzhang/Last_Millennium_Carbon_Experiment/CESM1/Fluxes/" ; set to your own path
  FilePath="../Fluxes/"

  do j= 0,nrun-1
  print(runner(j))
  FileName1="Total_annual_"+varname(i)+"_"+model+"_"+runner(j)+"_from_"+Startyr1+"-"+Endyr1+".txt"
  data1(j,:)=asciiread(FilePath+FileName1,(/nyear1/),"float")
  FileName2="Total_annual_"+varname(i)+"_"+model+"_"+runner(j)+"_from_"+Startyr2+"-"+Endyr2+".txt"
  data2(j,:)=asciiread(FilePath+FileName2,(/nyear2/),"float")
  end do

  data=new((/nrun,nyear1+nyear2/),"float")
  data(:,0:nyear1-1)=data1(:,:)
  data(:,nyear1:nyear1+nyear2-1)=data2(:,:)

  data(:,649)= 0.5*(dim_avg_n(data(:,651:655),1)+dim_avg_n(data(:,645:648),1))

  printVarSummary(data)
  data_avg = dim_avg_n(data,0)
  data_std = dim_stddev_n(data,0)

  data_up  = data_avg + data_std
  data_dw  = data_avg - data_std

;  print(avg(data_avg))
;  print(avg(data_up))
;  print(avg(data_dw))

  asciiwrite("Total_annual_"+varname(i)+"_"+model+"_ensemble_mean_from_"+Startyr1+"-"+Endyr2+".txt",data_avg)
;************************************************
; plotting
;************************************************
;  wks = gsn_open_wks("pdf",PlotOutFileName+varname(i))

 res =True

 res@gsnDraw              = False
 res@gsnFrame             = False

 res@gsnMaximize        = False
 res@vpKeepAspect       = True
 res@vpHeightF          = 0.15               ; change aspect ratio of plot
 res@vpWidthF           = 0.85

 res@trXMinF            = 850              ; set x-axis minimum
 res@trXMaxF            = 2005
 ;res@tmXBLabelsOn       = False
 res@tmXBOn             = True
 res@tmXBMode        = "Manual"
 res@tmXBTickStartF  = 850
 res@tmXBTickEndF    = 2005
 res@tmXBTickSpacingF= 100
 res@tmXBMajorOutwardLengthF = 0.01
 res@tmXBMinorOutwardLengthF = 0.00

 res@pmLegendDisplayMode    = "Never"             ; turn on legend
 res@pmLegendSide           = "Top"                ; Change location of
 res@lgPerimOn = False
 res@tiMainFontHeightF  = 0.024                    ; change font heights
 res@tiXAxisFontHeightF = 0.018
 res@tiYAxisFontHeightF = 0.018

 res@xyLineThicknesses = (/2./)
 res@xyLineColors      = Colors_line(i)
 res@xyExplicitLegendLabels = (/varname(i)/)
 res@tiYAxisString          = varname(i)+" (PgC yr~S1~- ~S1~1 )"

 if (varname(i).eq."HR") then
 res@xyExplicitLegendLabels = "RH"
 res@tiYAxisString          = "RH (PgC yr~S1~- ~S1~1 )"
 end if

 ; res@tiMainString           = "CESM1-CLM4CN Last Millennium Ensemble Project (850-1849)"
 if(i.eq.3)then
 res@tiXAxisString          = "Year"
 end if
 res@pmLegendParallelPosF   = 0.265                  ; move units right
 res@pmLegendOrthogonalPosF = -0.55                ; move units down
 res@pmLegendWidthF         = 0.08                  ; Change width and
 res@pmLegendHeightF        = 0.10                   ; height of legend.

 res@trYMinF            = YValue_Min(i)
 res@trYMaxF            = YValue_Max(i)

  plot(i)  = gsn_csm_xy (wks,year,(/data_avg/),res) ; create plot

; Create a plot with the area between both curves filled in blue.
  delete(res@xyLineColors)
  res@gsnXYFillColors = Colors_fill(i)
  res@xyLineColor     = -1          ; We don't want the line, so make it transparent.

  bot_plot(i)  = gsn_csm_xy (wks,year,(/data_up,data_dw/),res)  ; Create filled XY plot.

;  overlay(bot_plot(i),plot(i))

  x1 = (/1850,2005,2005,1850,1850/)
  y1 = (/YValue_Min(i),YValue_Min(i),YValue_Max(i),YValue_Max(i),YValue_Min(i)/)


  gres                 = True
  gres@gsFillColor     = "gray"
  gres@gsFillOpacityF  = 0.2
  gres@tfPolyDrawOrder = "Draw"  ; this can be used for polylines, polymarkers, or polygons

  dum(i) = gsn_add_polygon(wks,plot(i),x1,y1,gres)

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

    resP@txString = "The CESM Last Millennium Ensemble Project (850-2005)"
    resP@txFontHeightF = 0.022
   ;  resP@txPosYF  = 1.
    resP@txPosXF  = 0.55
    resP@tiXAxisString = "Year"
    resP@amJust   = "TopLeft"
    resP@gsnPanelFigureStringsPerimOn = False
    resP@gsnPanelFigureStrings = (/"(a)","(b)","(c)","(d)"/)
    resP@gsnPanelFigureStringsBackgroundFillColor = -1
    resP@gsnPanelFigureStringsFontHeightF = 0.015
    resP@gsnPanelCenter        = False
    gsn_panel(wks,bot_plot(:),(/nvar,1/),resP)

end
