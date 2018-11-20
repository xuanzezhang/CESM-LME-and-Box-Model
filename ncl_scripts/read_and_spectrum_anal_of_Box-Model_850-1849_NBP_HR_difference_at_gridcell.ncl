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

varname=(/"NBP","HR"/)
nvar=dimsizes(varname)

 do i= 0, nvar-1
   print(varname(i))
   FilePath="/data/xzhang/Last_Millennium_Carbon_Experiment/CESM1/Box-Model-3-pools_at_gridcell/" ; set to your own path

  do j= 0,nrun-1
   print(runner(j))
   FileNameIn=model+"_"+runner(j)+"_from_"+Startyr+"-"+Endyr+"_of_area_annual_"

   fin=addfile(FilePath+FileNameIn+"difference.nc","r")
   data      = fin->$varname(i)$  ;;;; year, latitude, longitude

   latitude  = fin->latitude
   longitude = fin->longitude
   nlat  = dimsizes(latitude)
   nlon  = dimsizes(longitude)

   printVarSummary(data)

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

  spcx_data=new((/500,nlat,nlon/),"float")

  do ilon=0,nlon-1
   do ilat=0,nlat-1
     if (all(.not.ismissing(data(:,ilat,ilon))))
     spec_data = specx_anal(data(:,ilat,ilon),d,sm,pct)
     frq =spec_data@frq
     ; print(spec_data@spcx)
     spcx_data(:,ilat,ilon)=spec_data@spcx
     delete(spec_data)
     end if
   end do
  end do

  nfrq=dimsizes(frq)
  frq!0="frq"
  frq&frq=frq

  spcx_data!0="frq"
  spcx_data!1="latitude"
  spcx_data!2="longitude"
  spcx_data&frq=frq
  spcx_data&latitude=latitude
  spcx_data&longitude=longitude

  OutputFileName="Spectrum_anal_of_"+varname(i)+"_difference_"+model+"_from_"+Startyr+"-"+Endyr
  system("rm -rf "+OutputFileName+"_at_gridcell.nc")
  fout=addfile(OutputFileName+"_at_gridcell.nc","c")
  fout->$varname(i)$=spcx_data
  fout->frq=frq

  
  delete(data)
  delete(spcx_data)

  end do

end do

end


