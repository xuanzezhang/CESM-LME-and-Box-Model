;*************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
Startyr=850
Endyr=1849

year=ispan(Startyr,Endyr,1)
nyear=Endyr-Startyr+1

model="CESM1_CLM4"

runner=(/"001","002","003","004","005","006","007","008","009","010","011","012","013"/)
nrun=dimsizes(runner)

varname=(/"D"/)
nvar=dimsizes(varname)

 do i= 0, nvar-1
   print(varname(i))

  FilePath="/data/xzhang/Last_Millennium_Carbon_Experiment/CESM1/Fluxes/" ; set to your own path

  do j= 0,nrun-1
  print(runner(j))

  FileNameAnn=model+"_"+runner(j)+"_from_"+Startyr+"-"+Endyr+"_of_area_annual_"
  OutFileNameAnn=model+"_"+runner(j)+"_from_"+Startyr+"-"+Endyr+"_of_area_annual_"

  fin1=addfile(FilePath+FileNameAnn+"NEP.nc","r")
  fin2=addfile(FilePath+FileNameAnn+"NBP.nc","r")
  NEP=fin1->NEP
  NBP=fin2->NBP
  area=fin1->area
  landmask=fin1->landmask

   data_cf_ann_area=NEP
   data_cf_ann_area=NEP-NBP

   data_ann_area_tot=dim_sum_n(data_cf_ann_area,(/1,2/))
   print(data_ann_area_tot)
   asciiwrite("Total_annual_"+varname(i)+"_"+model+"_"+runner(j)+"_from_"+Startyr+"-"+Endyr+".txt",data_ann_area_tot)

   data_ann_area_tropic=dim_sum_n(data_cf_ann_area(year|:,{latitude|-20:20},longitude|:),(/1,2/))
   print(data_ann_area_tropic)

   data_ann_area_NH=dim_sum_n(data_cf_ann_area(year|:,{latitude|20:90},longitude|:),(/1,2/))
   print(data_ann_area_NH)

   data_ann_area_SH=dim_sum_n(data_cf_ann_area(year|:,{latitude|-90:-20},longitude|:),(/1,2/))
   print(data_ann_area_SH)

   asciiwrite("Total_annual_"+varname(i)+"_"+model+"_"+runner(j)+"_from_"+Startyr+"-"+Endyr+"_tropic.txt",data_ann_area_tropic)

   asciiwrite("Total_annual_"+varname(i)+"_"+model+"_"+runner(j)+"_from_"+Startyr+"-"+Endyr+"_NH.txt",data_ann_area_NH)

   asciiwrite("Total_annual_"+varname(i)+"_"+model+"_"+runner(j)+"_from_"+Startyr+"-"+Endyr+"_SH.txt",data_ann_area_SH)

   system("rm -f "+OutFileNameAnn+varname(i)+".nc")
   fout=addfile(OutFileNameAnn+varname(i)+".nc","c")
   fout->$varname(i)$=data_cf_ann_area
   fout->area=area
   fout->landmask=landmask


   delete(data_cf_ann_area)
  end do
 end do

end 
