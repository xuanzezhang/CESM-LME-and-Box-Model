;*************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

  data_Vol_in=asciiread("Global_volcanic_eruptions_501-2000_V2_data.txt",(/2000-501+1/),"float")
  data_Vol_in!0="year"
  data_Vol_in&year=ispan(501,2000,1)

  nyvol=2000-850+1
  data_Vol=data_Vol_in({year|850:2000})
  year_Vol=ispan(850,2000,1)
  Year_Vol=new((/18/),"float")
  Mass_Vol=new((/18/),"float")

  iv=0
  do iy=0,nyvol-1
    if (data_Vol(iy).ge.25) then
    Mass_Vol(iv)=data_Vol(iy)
    Year_Vol(iv)=year_Vol(iy)
    iv=iv+1
    end if
  end do

;  print(Mass_Vol)
  print(Year_Vol)  


  nrun=18  ; dimsizes(runner)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
Startyr=850
Endyr=2005

year=ispan(Startyr,Endyr,1)
nyear=Endyr-Startyr+1

model="CESM1_CLM4"

varname=(/"NBP","NPP","HR","D","TS","RAIN","FSDS"/)
varname_len=(/"NBP","NPP","RH","D","TS","RAIN","FSDS"/)
nvar=dimsizes(varname)

  PlotOutFileName="plot_scatters_and_timeseries_of_CESM1_LME_850-2005_fluxes_"
  wks = gsn_open_wks("pdf",PlotOutFileName+"NBP_NPP_HR_D_in_Volcanic_eruptions_with_T_p_value_V4")
  plot=new(nvar,graphic)
  bot_plot=new(nvar,graphic)
  plot_scatter=new(nvar,graphic)

  Colors_line=(/"forestgreen","blue","red","black"/)
  Colors_fill=(/"gray","lightblue","pink2","grey"/)
  YValue_Max=(/3,3,3,3.,1.5,50.,10./)
  YValue_Min=(/-5,-5,-5,-5,-2.5,-100.,-15./)
  YValue_Max_sc=(/2., 2., 2.,2, 1.0, 30., 10./)
  YValue_Min_sc=(/-7.,-7.,-7.,-7,-4.5,-180.,-20./)


  data_var=new((/nvar,nrun/),float)
  data_gama=new((/nvar,nrun/),float)
  var     =new(nvar,float)

 do i= 0, nvar-1
   print(varname(i))

  FilePath= "./";"/data/xzhang/Last_Millennium_Carbon_Experiment/CESM1/Plot/" ; set to your own path
  data=new((/nyear/),"float")

  FileName="Total_annual_"+varname(i)+"_"+model+"_ensemble_mean_from_"+Startyr+"-"+Endyr+".txt"
  data=asciiread(FilePath+FileName,(/nyear/),"float")
  printVarSummary(data)
  data!0="year"
  data&year=year

  data_run=new((/nrun,11/),"float")

  do irun=0, nrun-1
    data_run(irun,:)=data({year|Year_Vol(irun)-5:Year_Vol(irun)+5})  ;; 1258
  end do

  data_run_ann=dim_rmvmean_n(data_run,1)

 ; data_run(0,:)=data({year|1253:1263})  ;; 1258
 ; data_run(1,:)=data({year|1447:1457})  ;; 1452
 ; data_run(2,:)=data({year|1208:1218})  ;; 1213
 ; data_run(3,:)=data({year|1595:1605})  ;; 1560
 ; data_run(4,:)=data({year|1986:1996})  ;; 1991

  data_run_avg=dim_avg_n(data_run_ann,0)
  data_run_std = dim_stddev_n(data_run_ann,0)

  data_run_up  = data_run_avg + data_run_std
  data_run_dw  = data_run_avg - data_run_std

  diff_data_run_ann=dim_sum_n(data_run_ann(:,5:7),1)-dim_avg_n(data_run_ann(:,0:4),1)

  data_var(i,:)=diff_data_run_ann(:)
;************************************************
; Create an array to hold both the original data
; and the calculated regression line.
;************************************************
 rc    = regline(Mass_Vol, diff_data_run_ann)
 print(rc)

 R=rc*stddev(Mass_Vol)/stddev(diff_data_run_ann)

 R_2=R*R

 print(R_2)

   df   = rc@nptxy-2   ; degrees of freedom
   tval = rc@tval      ; t-statistic
   b    = 0.5          ; b must be same size as tval (and df)
   prob  = betainc(df/(df+tval^2),df/2.0,b)
   print(prob)
;;;;;;;;;;;;
 data_scatter  = new ( (/2,dimsizes(Mass_Vol)/), typeof(Mass_Vol))
 data_scatter(0,:) = diff_data_run_ann

; y = mx+b
; m is the slope:       rc      returned from regline
; b is the y intercept: rc@yave attribute of rc returned from regline

 data_scatter(1,:) = rc*(Mass_Vol-rc@xave) + rc@yave


  var(i)=rc
;************************************************
; plotting
;************************************************
;  wks = gsn_open_wks("pdf",PlotOutFileName+varname(i))

 res =True

 res@gsnDraw              = False
 res@gsnFrame             = False

 res@gsnMaximize        = False
 res@vpKeepAspect       = True
 res@vpHeightF          = 0.3               ; change aspect ratio of plot
 res@vpWidthF           = 0.7

 res@trXMinF            = -5             ; set x-axis minimum
 res@trXMaxF            = 5
 ;res@tmXBLabelsOn       = False
 res@tmXBOn             = True
 res@tmXTOn             = False
 res@tmXBMode        = "Manual"
 res@tmXBMinorOn = False
 res@tmXTMinorOn = False
 res@tmXBTickStartF  = -5
 res@tmXBTickEndF    = 5
 res@tmXBTickSpacingF= 1
 res@tmXBMajorOutwardLengthF = 0.02
 res@tmXBMinorOutwardLengthF = 0.00


 res@pmLegendDisplayMode    = "Always"             ; turn on legend
; res@pmLegendSide           = "Left"                ; Change location of
 res@lgPerimOn = False
 res@tiMainFontHeightF  = 0.024                    ; change font heights
 res@tiXAxisFontHeightF = 0.018
 res@tiYAxisFontHeightF = 0.018

 res@xyLineThicknesses = (/2.,3.,3.,3.,3./)
 res@xyLineColors      = (/"black","forestgreen","blue","red","gold"/) ;(/Colors_line(i),Colors_line(i),Colors_line(i),Colors_line(i)/)
 res@xyExplicitLegendLabels = (/"Mean of 18 volcanic eruptions","1258","1452","1815","1991 (Mount Pinatubo)"/)
 res@tiYAxisString          = "Anomaly (PgC yr~S1~- ~S1~1 )"
 if (varname(i).eq."TS") then
 res@tiYAxisString          = "Anomaly (~S1~o C)"
 end if
; res@tiMainString           = "CESM1-CLM4CN Last Millennium Ensemble Project (850-1849)"
 if (varname(i).eq."RAIN") then
 res@tiYAxisString          = "Anomaly (mm)"
 end if
 if (varname(i).eq."FSDS") then
 res@tiYAxisString          = "Anomaly (W m~S1~- ~S1~2 )"
 end if
 res@tiXAxisString          = "Event Year"

 res@pmLegendParallelPosF   = 0.32                  ; move units right
 res@pmLegendOrthogonalPosF = -0.75                 ; move units down
 res@pmLegendWidthF         = 0.11                  ; Change width and
 res@pmLegendHeightF        = 0.13
 res@lgLabelFontHeightF     = 0.018


; res@trYMinF            = -3.
; res@trYMaxF            = 3.
  res@trYMinF            = YValue_Min(i)
  res@trYMaxF            = YValue_Max(i)

  event_year=ispan(-5,5,1)

  plot(i)  = gsn_csm_xy (wks,event_year,(/data_run_avg, data_run_ann(4,:),data_run_ann(8,:),data_run_ann(15,:),data_run_ann(17,:)/),res) ; create plot

; Create a plot with the area between both curves filled in blue.
  delete(res@xyLineColors)
  res@gsnXYFillColors = "gray"; Colors_fill(i)
  res@xyLineColor     = -1          ; We don't want the line, so make it transparent.
  res@pmLegendDisplayMode    = "Never"             ; turn on legend

  bot_plot(i)  = gsn_csm_xy (wks,event_year,(/data_run_up,data_run_dw/),res)  ; Create filled XY plot.

  overlay(bot_plot(i),plot(i))

;  draw(bot_plot)
;  frame(wks)

 delete(data)
 delete(res@gsnXYFillColors)

;;;;;;;;;;;;;;
 res_sc =True

 res_sc@gsnDraw              = False
 res_sc@gsnFrame             = False

 res_sc@gsnMaximize        = False
 res_sc@vpKeepAspect       = True
 res_sc@vpHeightF          = 0.48               ; change aspect ratio of plot
 res_sc@vpWidthF           = 0.52

 res_sc@trXMinF            = 0            ; set x-axis minimum
 res_sc@trXMaxF            = 300
 res_sc@trYMinF            = YValue_Min_sc(i)
 res_sc@trYMaxF            = YValue_Max_sc(i)
 res_sc@tmXBLabelFontHeightF = 0.03
 res_sc@tmYLLabelFontHeightF = 0.03
 res_sc@pmLegendDisplayMode    = "Never"             ; turn on legend
 res_sc@lgPerimOn = False
 res_sc@tiMainFontHeightF  = 0.024                    ; change font heights
 res_sc@tiXAxisFontHeightF = 0.03
 res_sc@tiYAxisFontHeightF = 0.03
 res_sc@xyMarkLineModes     = (/"Markers","Lines"/)  ; choose which have markers
 res_sc@xyMarkers           = 16                     ; choose type of marker
 res_sc@xyMarkerColor       = "red"                  ; Marker color
 res_sc@xyMarkerSizeF       = 0.008                  ; Marker size (default 0.01)
 res_sc@xyDashPatterns      = 1                      ; solid line
 res_sc@xyLineThicknesses   = (/2,3/)                ; set second line to 2

 res_sc@tiYAxisString          = "~F33~D~F21~"+varname_len(i)+" Anomaly (PgC yr~S1~- ~S1~1 )"
 if (varname(i).eq."TS") then
 res_sc@tiYAxisString          = "~F33~D~F21~T Anomaly (~S1~o C)"
 end if

 if (varname(i).eq."FSDS") then
 res_sc@tiXAxisString          = "Volcanic Aerosols (Tg)"
 end if
 res_sc@tiXAxisString          = "Volcanic Aerosols (Tg)"
; res_sc@tiMainString           = "CESM1-CLM4CN Last Millennium Ensemble Project (850-1849)"


  plot_scatter(i)  = gsn_csm_xy (wks,Mass_Vol, data_scatter,res_sc)

  label=(/"R~S1~2 =0.27, p=0.28","R~S1~2 =0.67, p<0.001","R~S1~2 =0.63, p<0.001","R~S1~2 =0.41, p=0.004","R~S1~2 =0.56, p<0.001","R~S1~2 =0.66, p<0.001","R~S1~2 =0.68, p<0.001"/)
  tx = (/100,100,100,100,100,100,100/)
  ty = (/-6,-6,-6,-6,-3.8,-150,-16/)
  txres_sc               = True  
  txres_sc@txFontHeightF = 0.035 
  text = gsn_add_text(wks,plot_scatter(i),label(i),120,ty(i),txres_sc)


end do


  do irun=0,nrun-1
  data_gama(:,irun)=data_var(:,irun)/data_var(4,irun)
  end do

  data_gama_avg=dim_avg_n(data_gama,1)
  data_gama_std=dim_stddev_n(data_gama,1)

  print(data_gama_avg)
  print(data_gama_std)

  print(var/var(4))
;************************************************
; create panel
;************************************************

 resP1                  = True                ; mods desired
 resP1@gsnFrame         = False               ; save panel until both ready
 resP1@gsnPanelTop      = 0.965
 resP1@gsnPanelLeft     = 0.12
 resP1@gsnPanelRight    = 0.55                 ; draw from left edge to center
 resP1@amJust   = "TopLeft"
 resP1@gsnPanelFigureStringsPerimOn = False
 resP1@gsnPanelFigureStrings = (/"(a)  NBP","(b)  NPP","(c)  RH","(d)  D","(e)  T","(f)  P","(g)  R"/)
 resP1@gsnPanelFigureStringsBackgroundFillColor = -1
 resP1@gsnPanelFigureStringsFontHeightF = 0.015

 resP2                  = True                ; mods desired
 resP2@gsnFrame         = False               ; save panel until both ready
 resP2@gsnPanelTop      = 0.965  
 resP2@gsnPanelLeft     = 0.455
 resP2@gsnPanelRight    = 0.70 		      ; draw from center to right edge 
 resP2@amJust   = "TopRight"
 resP2@gsnPanelFigureStringsPerimOn = False
 resP2@gsnPanelFigureStrings = (/"NBP","NPP","RH","D","T","P","R"/)
 resP2@gsnPanelFigureStringsBackgroundFillColor = -1
 resP2@gsnPanelFigureStringsFontHeightF = 0.015

 gsn_panel(wks,bot_plot,(/7,1/),resP1)          ; create first panel
 gsn_panel(wks,plot_scatter,(/7,1/),resP2)          ; create second panel

  txres                 = True  
  txres@txFontHeightF = 0.011                     ; text font height
  txres@txJust        = "CenterLeft"              ; Default is "CenterCenter".
 
  gsn_text_ndc(wks,"Responces of Land C Fluxes to Volcanic Eruptions in the CESM-LME",0.20,0.985,txres)

 frame(wks)                                   ; draw both on same page


end
