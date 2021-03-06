;*************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
Startyr=850
Endyr=1849

year=ispan(Startyr,Endyr,1)
nyear=Endyr-Startyr+1

model="CESM1_CLM4_0850cntl"

nrun="001"

varname=(/"NPP"/)
nvar=dimsizes(varname)

daysofmonth=(/31,28,31,30,31,30,31,31,30,31,30,31/)

FilePath="/data/xzhang/CESM1_Last_Millennium_Ensemble_Project/NPP/" ; set to your own path

FileName="b.e11.B1850C5CN.f19_g16.0850cntl.001.clm2.h0."  ; output file name

OutFileName=model+"_"+nrun+"_from_"+Startyr+"-"+Endyr+"_of_conformed_"
OutFileNameAnn=model+"_"+nrun+"from_"+Startyr+"-"+Endyr+"_of_conformed_annual_"

OutFileName=model+"_"+nrun+"_from_"+Startyr+"-"+Endyr+"_of_area_"
OutFileNameAnn=model+"_"+nrun+"_from_"+Startyr+"-"+Endyr+"_of_area_annual_"

 do i= 0,nvar-1
   print(varname(i))

   fin=addfile(FilePath+FileName+varname(i)+".085001-184912.nc","r")
   time=fin->time
   latitude=fin->lat
   longitude=fin->lon
   latitude@units="degrees_north"
   longitude@units="degrees_east"
  ; print(latitude)
  ; print(longitude)
   data=fin->$varname(i)$  ;;; gC/m2/s
   data@_FillValue=-99999.
   data@missing_value=-99999.
   printVarSummary(data)
   ntime=dimsizes(time)
   nlat=dimsizes(latitude)
   nlon=dimsizes(longitude)

  ism=(Startyr-Startyr)*12
  iem=ntime-1

  data_cf=(/data(ism:iem,:,:)/)*1.0e-3
  data_cf!0="time"
  data_cf!1="latitude"
  data_cf!2="longitude"
  data_cf&time=ispan(1,(iem-ism+1),1)
  data_cf&latitude=latitude
  data_cf&longitude=longitude
  data_cf@units="KgC/m2/month"
  data_cf@_FillValue=-99999. ;  data@_FillValue
  data_cf@missing_value=-99999.  ; data@missing_value

  data_cf_ann=new((/nyear,nlat,nlon/),typeof(data))
  data_cf_ann!0="year"
  data_cf_ann!1="latitude"
  data_cf_ann!2="longitude"
  data_cf_ann&year=year
  data_cf_ann&latitude=latitude
  data_cf_ann&longitude=longitude
  data_cf_ann@units="KgC/m2/year"
  data_cf_ann@_FillValue=data_cf@_FillValue
  data_cf_ann@missing_value=data_cf@missing_value

  do iyear=0,nyear-1
     isy=iyear*12
     iey=isy+11
   data_cf(isy:iey,:,:)=data_cf(isy:iey,:,:)*conform(data_cf(isy:iey,:,:),daysofmonth,0)*24*3600
   data_cf_ann(iyear,:,:)=dim_sum_n(data_cf(isy:iey,:,:),0)
  end do

   system("rm -f "+OutFileName+varname(i)+".nc")
   fout=addfile(OutFileName+varname(i)+".nc","c")
   fout->$varname(i)$=data_cf

   system("rm -f "+OutFileNameAnn+varname(i)+".nc")
   fout=addfile(OutFileNameAnn+varname(i)+".nc","c")
   fout->$varname(i)$=data_cf_ann

   delete(data)
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   area=fin->area ;; unists km^2
   maskVeg=fin->landmask

   data_cf_area=data_cf
   data_cf_area=0.
   data_cf_area@units="PgC/month"  ;; from KgC/m2/month

   data_cf_ann_area=data_cf_ann
   data_cf_ann_area=0.
   data_cf_ann_area@units="PgC/year"

   area_cf=conform_dims(dimsizes(data_cf),area,(/1,2/))
   mask_cf=conform_dims(dimsizes(data_cf),maskVeg,(/1,2/))
   data_cf_area=data_cf*area_cf*1.0e6*mask_cf*1.0e-12 ;;PgC/month
   area_cf_ann=conform_dims(dimsizes(data_cf_ann),area,(/1,2/))
   mask_cf_ann=conform_dims(dimsizes(data_cf_ann),maskVeg,(/1,2/))
   data_cf_ann_area=data_cf_ann*area_cf_ann*1.0e6*mask_cf_ann*1.0e-12 ;;PgC/month


   printVarSummary(area_cf)
   printVarSummary(mask_cf)
   printVarSummary(data_cf)
   printVarSummary(data_cf_area)

   data_area_tot=dim_sum_n(data_cf_area,(/1,2/))
   print(data_area_tot)
   data_ann_area_tot=dim_sum_n(data_cf_ann_area,(/1,2/))
   print(data_ann_area_tot)

   asciiwrite("Total_"+varname(i)+"_"+model+"_"+nrun+"_from_"+Startyr+"-"+Endyr+".txt",data_area_tot)
   asciiwrite("Total_annual_"+varname(i)+"_"+model+"_"+nrun+"_from_"+Startyr+"-"+Endyr+".txt",data_ann_area_tot)

   system("rm -f "+OutFileName+varname(i)+".nc")
   fout=addfile(OutFileName+varname(i)+".nc","c")
   fout->$varname(i)$=data_cf_area
   fout->area=area
   fout->landmask=maskVeg

   system("rm -f "+OutFileNameAnn+varname(i)+".nc")
   fout=addfile(OutFileNameAnn+varname(i)+".nc","c")
   fout->$varname(i)$=data_cf_ann_area
   fout->area=area
   fout->landmask=maskVeg

   delete(data_cf)
   delete(data_cf_area)
   delete(data_cf_ann)
   delete(data_cf_ann_area)

 end do

end
