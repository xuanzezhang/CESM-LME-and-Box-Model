;*************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
Startyr=850
Endyr=1849

year=ispan(Startyr,Endyr,1)
nyear=Endyr-Startyr+1

model="CESM1_CLM4_VOLC_GRA"

runner=(/"001","002","003","004","005"/)
nrun=dimsizes(runner)

varname=(/"NBP","NPP","HR","D"/)
nvar=dimsizes(varname)

 do i= 0, nvar-1
   print(varname(i))

  FilePath="/data/xzhang/Last_Millennium_Carbon_Experiment/CESM1/Fluxes-Forcings/" ; set to your own path
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

  nfrq=dimsizes(frq)

  spcx_data!0="frq"
  spcx_data!1="run"
  spcx_data&frq=frq
  spcx_data&run=ispan(1,nrun,1)

  OutputFileName="Spectrum_anal_of_"+varname(i)+"_"+model+"_from_"+Startyr+"-"+Endyr
  system("rm -rf "+OutputFileName+".nc")
  fout=addfile(OutputFileName+".nc","c")
  fout->spcx_data=spcx_data

  
  delete(data)

end do

end


