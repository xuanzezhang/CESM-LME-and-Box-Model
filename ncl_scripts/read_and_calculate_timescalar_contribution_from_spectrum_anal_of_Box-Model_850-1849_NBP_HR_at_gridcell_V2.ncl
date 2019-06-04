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

   FilePath="/data/xzhang/Last_Millennium_Carbon_Experiment/CESM1/Spectral_analysis_for_Box-Model/"
   FilePathCESM="/data/xzhang/Last_Millennium_Carbon_Experiment/CESM1/Spectral_analysis_for_CESM-LME/GridCell/"

   i=0
   j=0
   FileNameIn="Spectrum_anal_of_"+varname(i)+"_"+model+"_"+runner(j)+"_from_"+Startyr+"-"+Endyr

   fin=addfile(FilePath+FileNameIn+"_at_gridcell.nc","r")
   frq       = fin->frq
   latitude  = fin->latitude
   longitude = fin->longitude
   nfrq  = dimsizes(frq)
   nlat  = dimsizes(latitude)
   nlon  = dimsizes(longitude)
   print(frq)

   spcx_NBP=new((/nrun,nfrq,nlat,nlon/),float)
   spcx_NBP!0="run"
   spcx_NBP!1="frq"
   spcx_NBP!2="latitude"
   spcx_NBP!3="longitude"
   spcx_NBP&run=ispan(1,nrun,1)
   spcx_NBP&frq=frq
   spcx_NBP&latitude=latitude
   spcx_NBP&longitude=longitude

   spcx_HR=new((/nrun,nfrq,nlat,nlon/),float)
   spcx_HR!0="run"
   spcx_HR!1="frq"
   spcx_HR!2="latitude"
   spcx_HR!3="longitude"
   spcx_HR&run=ispan(1,nrun,1)
   spcx_HR&frq=frq
   spcx_HR&latitude=latitude
   spcx_HR&longitude=longitude

   spcx_NBP_diff=spcx_NBP
   spcx_HR_diff=spcx_HR

   spcx_NBP_cesm=spcx_NBP
   spcx_HR_cesm=spcx_HR

  do j= 0,nrun-1
   print(runner(j))
   i=0
   FileNameIn="Spectrum_anal_of_"+varname(i)+"_"+model+"_"+runner(j)+"_from_"+Startyr+"-"+Endyr
   fin=addfile(FilePath+FileNameIn+"_at_gridcell.nc","r")
   spcx_NBP(j,:,:,:)      = fin->NBP  ;;;; year, latitude, longitude

   FileNameIn="Spectrum_anal_of_"+varname(i)+"_difference_"+model+"_"+runner(j)+"_from_"+Startyr+"-"+Endyr
   fin2=addfile(FilePath+FileNameIn+"_at_gridcell.nc","r")
   spcx_NBP_diff(j,:,:,:)      = fin2->NBP  ;;;; year, latitude, longitude

   FileNameIn="Spectrum_anal_of_"+varname(i)+"_CESM1_CLM4_"+runner(j)+"_from_"+Startyr+"-"+Endyr
   fin3=addfile(FilePathCESM+FileNameIn+"_at_gridcell.nc","r")
   spcx_NBP_cesm(j,:,:,:)      = fin3->NBP  ;;;; year, latitude, longitude

   i=1
   FileNameIn="Spectrum_anal_of_"+varname(i)+"_"+model+"_"+runner(j)+"_from_"+Startyr+"-"+Endyr
   fin=addfile(FilePath+FileNameIn+"_at_gridcell.nc","r")
   spcx_HR(j,:,:,:)      = fin->HR  ;;;; year, latitude, longitude

   FileNameIn="Spectrum_anal_of_"+varname(i)+"_difference_"+model+"_"+runner(j)+"_from_"+Startyr+"-"+Endyr
   fin2=addfile(FilePath+FileNameIn+"_at_gridcell.nc","r")
   spcx_HR_diff(j,:,:,:)      = fin2->HR

   FileNameIn="Spectrum_anal_of_"+varname(i)+"_CESM1_CLM4_"+runner(j)+"_from_"+Startyr+"-"+Endyr
   fin3=addfile(FilePathCESM+FileNameIn+"_at_gridcell.nc","r")
   spcx_HR_cesm(j,:,:,:)      = fin3->HR  ;;;; year, latitude, longitude

  end do

  spcx_NBP_avg=dim_avg_n(spcx_NBP,0) 
  spcx_NBP_std=dim_stddev_n(spcx_NBP,0)  
  spcx_NBP_up=spcx_NBP_avg + spcx_NBP_std
  spcx_NBP_dw=spcx_NBP_avg - spcx_NBP_std

  spcx_HR_avg=dim_avg_n(spcx_HR,0)  
  spcx_HR_std=dim_stddev_n(spcx_HR,0)
  spcx_HR_up=spcx_HR_avg + spcx_HR_std
  spcx_HR_dw=spcx_HR_avg - spcx_HR_std

  spcx_NBP_diff_avg=dim_avg_n(spcx_NBP_diff,0)
  spcx_NBP_diff_std=dim_stddev_n(spcx_NBP_diff,0)
  spcx_NBP_diff_up=spcx_NBP_diff_avg + spcx_NBP_diff_std
  spcx_NBP_diff_dw=spcx_NBP_diff_avg - spcx_NBP_diff_std

  spcx_HR_diff_avg=dim_avg_n(spcx_HR_diff,0)
  spcx_HR_diff_std=dim_stddev_n(spcx_HR_diff,0)
  spcx_HR_diff_up=spcx_HR_diff_avg + spcx_HR_diff_std
  spcx_HR_diff_dw=spcx_HR_diff_avg - spcx_HR_diff_std

  spcx_NBP_cesm_avg=dim_avg_n(spcx_NBP_cesm,0)
  spcx_NBP_cesm_std=dim_stddev_n(spcx_NBP_cesm,0)
  spcx_NBP_cesm_up=spcx_NBP_cesm_avg + spcx_NBP_cesm_std
  spcx_NBP_cesm_dw=spcx_NBP_cesm_avg - spcx_NBP_cesm_std

  spcx_HR_cesm_avg=dim_avg_n(spcx_HR_cesm,0)
  spcx_HR_cesm_std=dim_stddev_n(spcx_HR_cesm,0)
  spcx_HR_cesm_up=spcx_HR_cesm_avg + spcx_HR_cesm_std
  spcx_HR_cesm_dw=spcx_HR_cesm_avg - spcx_HR_cesm_std

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  dfrq=frq(1)-frq(0)

  spcx_NBP_10yr=new((/nlat,nlon/),"float")
  spcx_NBP_100yr=new((/nlat,nlon/),"float")
  spcx_NBP_1000yr=new((/nlat,nlon/),"float")
  spcx_NBP_10yr=dim_sum_n(spcx_NBP_avg(100:499,:,:),0)*dfrq
  spcx_NBP_100yr=dim_sum_n(spcx_NBP_avg(10:99,:,:),0)*dfrq
  spcx_NBP_1000yr=dim_sum_n(spcx_NBP_avg(0:9,:,:),0)*dfrq

  spcx_HR_10yr=new((/nlat,nlon/),"float")
  spcx_HR_100yr=new((/nlat,nlon/),"float")
  spcx_HR_1000yr=new((/nlat,nlon/),"float")
  spcx_HR_10yr=dim_sum_n(spcx_HR_avg(100:499,:,:),0)*dfrq
  spcx_HR_100yr=dim_sum_n(spcx_HR_avg(10:99,:,:),0)*dfrq
  spcx_HR_1000yr=dim_sum_n(spcx_HR_avg(0:9,:,:),0)*dfrq

  spcx_NBP_diff_10yr=new((/nlat,nlon/),"float")
  spcx_NBP_diff_100yr=new((/nlat,nlon/),"float")
  spcx_NBP_diff_1000yr=new((/nlat,nlon/),"float")
  spcx_NBP_diff_10yr=dim_sum_n(spcx_NBP_diff_avg(100:499,:,:),0)*dfrq
  spcx_NBP_diff_100yr=dim_sum_n(spcx_NBP_diff_avg(10:99,:,:),0)*dfrq
  spcx_NBP_diff_1000yr=dim_sum_n(spcx_NBP_diff_avg(0:9,:,:),0)*dfrq

  spcx_HR_diff_10yr=new((/nlat,nlon/),"float")
  spcx_HR_diff_100yr=new((/nlat,nlon/),"float")
  spcx_HR_diff_1000yr=new((/nlat,nlon/),"float")
  spcx_HR_diff_10yr=dim_sum_n(spcx_HR_diff_avg(100:499,:,:),0)*dfrq
  spcx_HR_diff_100yr=dim_sum_n(spcx_HR_diff_avg(10:99,:,:),0)*dfrq
  spcx_HR_diff_1000yr=dim_sum_n(spcx_HR_diff_avg(0:9,:,:),0)*dfrq

  spcx_NBP_cesm_10yr=new((/nlat,nlon/),"float")
  spcx_NBP_cesm_100yr=new((/nlat,nlon/),"float")
  spcx_NBP_cesm_1000yr=new((/nlat,nlon/),"float")
  spcx_NBP_cesm_10yr=dim_sum_n(spcx_NBP_cesm_avg(100:499,:,:),0)*dfrq
  spcx_NBP_cesm_100yr=dim_sum_n(spcx_NBP_cesm_avg(10:99,:,:),0)*dfrq
  spcx_NBP_cesm_1000yr=dim_sum_n(spcx_NBP_cesm_avg(0:9,:,:),0)*dfrq

  spcx_HR_cesm_10yr=new((/nlat,nlon/),"float")
  spcx_HR_cesm_100yr=new((/nlat,nlon/),"float")
  spcx_HR_cesm_1000yr=new((/nlat,nlon/),"float")
  spcx_HR_cesm_10yr=dim_sum_n(spcx_HR_cesm_avg(100:499,:,:),0)*dfrq
  spcx_HR_cesm_100yr=dim_sum_n(spcx_HR_cesm_avg(10:99,:,:),0)*dfrq
  spcx_HR_cesm_1000yr=dim_sum_n(spcx_HR_cesm_avg(0:9,:,:),0)*dfrq

  spcx_NBP_cesm_10yr=where(spcx_NBP_cesm_10yr.gt.0, spcx_NBP_cesm_10yr,spcx_NBP_cesm_10yr@_FillValue)
  spcx_NBP_cesm_100yr=where(spcx_NBP_cesm_100yr.gt.0, spcx_NBP_cesm_100yr,spcx_NBP_cesm_100yr@_FillValue)
  spcx_NBP_cesm_1000yr=where(spcx_NBP_cesm_1000yr.gt.0, spcx_NBP_cesm_1000yr,spcx_NBP_cesm_1000yr@_FillValue)

  spcx_HR_cesm_10yr=where(spcx_HR_cesm_10yr.gt.0, spcx_HR_cesm_10yr,spcx_HR_cesm_10yr@_FillValue)
  spcx_HR_cesm_100yr=where(spcx_HR_cesm_100yr.gt.0, spcx_HR_cesm_100yr,spcx_HR_cesm_100yr@_FillValue)
  spcx_HR_cesm_1000yr=where(spcx_HR_cesm_1000yr.gt.0, spcx_HR_cesm_1000yr,spcx_HR_cesm_1000yr@_FillValue)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  NPP_derived_NBP_10yr=new((/nlat,nlon/),"float")
  NPP_derived_NBP_100yr=new((/nlat,nlon/),"float")
  NPP_derived_NBP_1000yr=new((/nlat,nlon/),"float")

  NPP_derived_HR_10yr=new((/nlat,nlon/),"float")
  NPP_derived_HR_100yr=new((/nlat,nlon/),"float")
  NPP_derived_HR_1000yr=new((/nlat,nlon/),"float")

  Tau_derived_NBP_10yr=new((/nlat,nlon/),"float")
  Tau_derived_NBP_100yr=new((/nlat,nlon/),"float")
  Tau_derived_NBP_1000yr=new((/nlat,nlon/),"float")

  Tau_derived_HR_10yr=new((/nlat,nlon/),"float")
  Tau_derived_HR_100yr=new((/nlat,nlon/),"float")
  Tau_derived_HR_1000yr=new((/nlat,nlon/),"float")

;  NPP_derived_NBP_10yr= spcx_NBP_10yr/(spcx_NBP_10yr + spcx_NBP_diff_10yr)
;  NPP_derived_NBP_100yr= spcx_NBP_100yr/(spcx_NBP_100yr + spcx_NBP_diff_100yr)
;  NPP_derived_NBP_1000yr= spcx_NBP_1000yr/(spcx_NBP_1000yr + spcx_NBP_diff_1000yr)

;  NPP_derived_HR_10yr= spcx_HR_10yr/(spcx_HR_10yr + spcx_HR_diff_10yr)
;  NPP_derived_HR_100yr= spcx_HR_100yr/(spcx_HR_100yr + spcx_HR_diff_100yr)
;  NPP_derived_HR_1000yr= spcx_HR_1000yr/(spcx_HR_1000yr + spcx_HR_diff_1000yr)

;  Tau_derived_NBP_10yr= spcx_NBP_diff_10yr/(spcx_NBP_10yr + spcx_NBP_diff_10yr)
;  Tau_derived_NBP_100yr= spcx_NBP_diff_100yr/(spcx_NBP_100yr + spcx_NBP_diff_100yr)
;  Tau_derived_NBP_1000yr= spcx_NBP_diff_1000yr/(spcx_NBP_1000yr + spcx_NBP_diff_1000yr)

;  Tau_derived_HR_10yr= spcx_HR_diff_10yr/(spcx_HR_10yr + spcx_HR_diff_10yr)
;  Tau_derived_HR_100yr= spcx_HR_diff_100yr/(spcx_HR_100yr + spcx_HR_diff_100yr)
;  Tau_derived_HR_1000yr= spcx_HR_diff_1000yr/(spcx_HR_1000yr + spcx_HR_diff_1000yr)


  NPP_derived_NBP_10yr= spcx_NBP_10yr/spcx_NBP_cesm_10yr
  NPP_derived_NBP_100yr= spcx_NBP_100yr/spcx_NBP_cesm_100yr
  NPP_derived_NBP_1000yr= spcx_NBP_1000yr/spcx_NBP_cesm_1000yr

  NPP_derived_HR_10yr= spcx_HR_10yr/spcx_HR_cesm_10yr
  NPP_derived_HR_100yr= spcx_HR_100yr/spcx_HR_cesm_100yr
  NPP_derived_HR_1000yr= spcx_HR_1000yr/spcx_HR_cesm_1000yr

  Tau_derived_NBP_10yr= spcx_NBP_diff_10yr/spcx_NBP_cesm_10yr
  Tau_derived_NBP_100yr= spcx_NBP_diff_100yr/spcx_NBP_cesm_100yr
  Tau_derived_NBP_1000yr= spcx_NBP_diff_1000yr/spcx_NBP_cesm_1000yr

  Tau_derived_HR_10yr= spcx_HR_diff_10yr/spcx_HR_cesm_10yr
  Tau_derived_HR_100yr= spcx_HR_diff_100yr/spcx_HR_cesm_100yr
  Tau_derived_HR_1000yr= spcx_HR_diff_1000yr/spcx_HR_cesm_1000yr

  NPP_derived_NBP_10yr!0="latitude"
  NPP_derived_NBP_10yr!1="longitude"
  NPP_derived_NBP_10yr&latitude=latitude
  NPP_derived_NBP_10yr&longitude=longitude

  copy_VarCoords(NPP_derived_NBP_10yr,NPP_derived_NBP_100yr)
  copy_VarCoords(NPP_derived_NBP_10yr,NPP_derived_NBP_1000yr)
  copy_VarCoords(NPP_derived_NBP_10yr,NPP_derived_HR_10yr)
  copy_VarCoords(NPP_derived_NBP_10yr,NPP_derived_HR_100yr)
  copy_VarCoords(NPP_derived_NBP_10yr,NPP_derived_HR_1000yr)

  copy_VarCoords(NPP_derived_NBP_10yr,Tau_derived_NBP_10yr)
  copy_VarCoords(NPP_derived_NBP_10yr,Tau_derived_NBP_100yr)
  copy_VarCoords(NPP_derived_NBP_10yr,Tau_derived_NBP_1000yr)
  copy_VarCoords(NPP_derived_NBP_10yr,Tau_derived_HR_10yr)
  copy_VarCoords(NPP_derived_NBP_10yr,Tau_derived_HR_100yr)
  copy_VarCoords(NPP_derived_NBP_10yr,Tau_derived_HR_1000yr)

  FileNameOut="Timescalar_contribution_from_spectrum_anal_of_NBP_HR_"+model+"_ensemble_mean_from_"+Startyr+"-"+Endyr
  system("rm -f "+FilePath+FileNameOut+"_at_gridcell_V2.nc")
  fout=addfile(FilePath+FileNameOut+"_at_gridcell_V2.nc","c")
  fout->NPP_derived_NBP_10yr=NPP_derived_NBP_10yr
  fout->NPP_derived_NBP_100yr=NPP_derived_NBP_100yr
  fout->NPP_derived_NBP_1000yr=NPP_derived_NBP_1000yr
  fout->NPP_derived_HR_10yr=NPP_derived_HR_10yr
  fout->NPP_derived_HR_100yr=NPP_derived_HR_100yr
  fout->NPP_derived_HR_1000yr=NPP_derived_HR_1000yr
  fout->Tau_derived_NBP_10yr=Tau_derived_NBP_10yr
  fout->Tau_derived_NBP_100yr=Tau_derived_NBP_100yr
  fout->Tau_derived_NBP_1000yr=Tau_derived_NBP_1000yr
  fout->Tau_derived_HR_10yr=Tau_derived_HR_10yr
  fout->Tau_derived_HR_100yr=Tau_derived_HR_100yr
  fout->Tau_derived_HR_1000yr=Tau_derived_HR_1000yr



end


