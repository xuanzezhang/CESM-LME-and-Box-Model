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

InputPath="/data/xzhang/Last_Millennium_Carbon_Experiment/CESM1/Fluxes/"

InputCpoolPath="/data/xzhang/Last_Millennium_Carbon_Experiment/CESM1/Pools/"

FilePathOut="/data/xzhang/Last_Millennium_Carbon_Experiment/CESM1/Box-Model-3-pools_at_gridcell/"

 do j= 0, 0 ;nrun-1

   FileNameIn=model+"_"+runner(j)+"_from_"+Startyr+"-"+Endyr+"_of_area_annual_" 
   FileNameOut="Box_Model_"+runner(j)+"_from_"+Startyr+"-"+Endyr+"_of_area_annual_"

   fnpp=addfile(InputPath+FileNameIn+"NPP.nc","r")
   npp_in      = fnpp->NPP  ;;;; year, latitude, longitude

   latitude  = fnpp->latitude
   longitude = fnpp->longitude
   nlat  = dimsizes(latitude)
   nlon  = dimsizes(longitude)

   fnbp=addfile(InputPath+FileNameIn+"NBP.nc","r")
   nbp_in      = fnbp->NBP

   frh=addfile(InputPath+FileNameIn+"HR.nc","r")
   rh_in       = frh->HR

   fD=addfile(InputPath+FileNameIn+"D.nc","r")
   Disturb_in = fD->D

   Clabloss_in= Disturb_in

   fCplant=addfile(InputCpoolPath+FileNameIn+"TOTVEGC.nc","r")
   Cplant_in       = fCplant->TOTVEGC

   fCsoil=addfile(InputCpoolPath+FileNameIn+"SOILC.nc","r")
   Csoil_in       = fCsoil->SOILC

   fClitter=addfile(InputCpoolPath+FileNameIn+"LITTERC.nc","r")
   Clitter_in      = fClitter->LITTERC

   fCWDC=addfile(InputCpoolPath+FileNameIn+"CWDC.nc","r")
   Ccwd_in         = fCWDC->CWDC

   Clitter_in =  Clitter_in + Ccwd_in   

   fCsoil_rh=addfile(InputPath+FileNameIn+"SOILC_HR.nc","r")
   Csoil_rh_in       = fCsoil_rh->SOILC_HR

   fClitter_rh=addfile(InputPath+FileNameIn+"LITTERC_HR.nc","r")
   Clitter_rh_in      = fClitter_rh->LITTERC_HR

   fCWDC_rh=addfile(InputPath+FileNameIn+"CWDC_HR.nc","r")
   Ccwd_rh_in         = fCWDC_rh->CWDC_HR

   Clitter_rh_in =  Clitter_rh_in + Ccwd_rh_in

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;Box Model parameters set-up;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

Cplant_in=where(Cplant_in.gt.1.e-6,Cplant_in,Cplant_in@_FillValue)
Clitter_in=where(Clitter_in.gt.1.e-6,Clitter_in,Clitter_in@_FillValue)
Csoil_in=where(Csoil_in.gt.1.e-6,Csoil_in,Csoil_in@_FillValue)


  npp_in_avg=dim_avg_n(npp_in,0) ;;;; year, latitude, longitude
  rh_in_avg=dim_avg_n(rh_in,0)
  nbp_in_avg=dim_avg_n(nbp_in,0)

  Clabloss_in_avg=dim_avg_n(Clabloss_in,0)
  Cplant_in_avg=dim_avg_n(Cplant_in,0)
  Csoil_in_avg=dim_avg_n(Csoil_in,0)
  Clitter_in_avg=dim_avg_n(Clitter_in,0)

  Clitter_rh_in_avg=dim_avg_n(Clitter_rh_in,0)
  Csoil_rh_in_avg=dim_avg_n(Csoil_rh_in,0)

  A22=npp_in ;; year, latitude, longitude
  A33=npp_in
  A21=npp_in
  A32=npp_in

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;Box Model;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 epsilon=1.
  A22=0.
  A33=0.
  A21=0.
  A32=0.

  Nyr = nyear

  Dcplant_in=Cplant_in
  Dclitter_in=Cplant_in
  Dcsoil_in=Cplant_in

  Dcplant_in=0.
  Dclitter_in=0.
  Dcsoil_in=0.

 do iyear=1,nyear-1
  Dcplant_in(iyear,:,:)=Cplant_in(iyear,:,:) - Cplant_in(iyear-1,:,:)
  Dclitter_in(iyear,:,:)=Clitter_in(iyear,:,:) - Clitter_in(iyear-1,:,:)
  Dcsoil_in(iyear,:,:)=Csoil_in(iyear,:,:) - Csoil_in(iyear-1,:,:)
 end do

;; Dcplant/Dt  = -epsilon*A21*Cplant_in_avg + npp_in_avg
;; Dclitter/nyear = epsilon*A21*Cplant_in_avg - epsilon*A32*Clitter_in_avg - epsilon*A22*Clitter_in_avg
;; Dcsoil/nyear  = epsilon*A32*Clitter_in_avg - epsilon*A33*Csoil_in_avg
;; Clitter_rh_in_avg   = epsilon*A22*Clitter_in_avg
;; Csoil_rh_in_avg     =  epsilon*A33*Csoil_in_avg

   A21 = ( npp_in - Clabloss_in - Dcplant_in ) / epsilon / Cplant_in
   A22 = (Clitter_rh_in)/ epsilon / Clitter_in
   A33 = (Csoil_rh_in) / epsilon / Csoil_in
   A32 = ( epsilon * A33 * Csoil_in + Dcsoil_in ) / epsilon / Clitter_in

  A22_avg=dim_avg_n(A22(1:nyear-1,:,:),0)
  A33_avg=dim_avg_n(A33(1:nyear-1,:,:),0)
  A21_avg=dim_avg_n(A21(1:nyear-1,:,:),0)
  A32_avg=dim_avg_n(A32(1:nyear-1,:,:),0)

  A21_avg=where(A21_avg.le.1.and.A21_avg.ge.0,A21_avg,A21_avg@_FillValue)
  A22_avg=where(A22_avg.le.1.and.A22_avg.ge.0,A22_avg,A22_avg@_FillValue)
  A32_avg=where(A32_avg.le.1.and.A32_avg.ge.0,A32_avg,A32_avg@_FillValue)
  A33_avg=where(A33_avg.le.1.and.A33_avg.ge.0,A33_avg,A33_avg@_FillValue)

   printVarSummary(A22_avg)
;   printVarSummary(A33_avg)
;   printVarSummary(A21_avg)
;   printVarSummary(A32_avg)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   Cplant_out = Cplant_in
   Clitter_out = Clitter_in
   Csoil_out = Csoil_in
   nbp_out  = nbp_in
   rh_out   = rh_in

   Cplant_out = 0.
   Clitter_out = 0.
   Csoil_out = 0.
   nbp_out = 0.
   rh_out = 0.

   nbp_out(0,:,:)=nbp_in(0,:,:)
   rh_out(0,:,:)=rh_in(0,:,:)
   Cplant_out(0,:,:) = Cplant_in(0,:,:)
   Clitter_out(0,:,:) = Clitter_in(0,:,:)
   Csoil_out(0,:,:) = Csoil_in(0,:,:)

   Cplant_out(1,:,:) = Cplant_in(1,:,:)
   Clitter_out(1,:,:) = Clitter_in(1,:,:)
   Csoil_out(1,:,:) = Csoil_in(1,:,:)

   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;;; spinup
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   do iyear = 1, nyear-1

   Cplant_out(iyear,:,:) = ((npp_in(0,:,:) - Clabloss_in(0,:,:)) - epsilon * A21_avg * Cplant_out(iyear-1,:,:)) \
                        + Cplant_out(iyear-1,:,:)
   Clitter_out(iyear,:,:) = (epsilon * A21_avg * Cplant_out(iyear-1,:,:) - epsilon * A32_avg * Clitter_out(iyear-1,:,:) \
                        - epsilon * A22_avg * Clitter_out(iyear-1,:,:))+ Clitter_out(iyear-1,:,:)
   Csoil_out(iyear,:,:) = (epsilon * A32_avg * Clitter_out(iyear-1,:,:) - epsilon * A33_avg * Csoil_out(iyear-1,:,:)) + Csoil_out(iyear-1,:,:)
   nbp_out(iyear,:,:) = ( Cplant_out(iyear,:,:) + Clitter_out(iyear,:,:) + Csoil_out(iyear,:,:) )  \
                     - ( Cplant_out(iyear-1,:,:) + Clitter_out(iyear-1,:,:) + Csoil_out(iyear-1,:,:) )
   rh_out(iyear,:,:) = ((npp_in(0,:,:) - Clabloss_in(0,:,:)) - nbp_out(iyear,:,:))
   end do

 ;  print(nbp_out)
 ;  print(rh_out)

   Cplant_out(0,:,:) = Cplant_out(nyear-1,:,:)
   Clitter_out(0,:,:) = Clitter_out(nyear-1,:,:)
   Csoil_out(0,:,:) = Csoil_out(nyear-1,:,:)

   Cplant_out(1,:,:) = Cplant_out(nyear-1,:,:)
   Clitter_out(1,:,:) = Clitter_out(nyear-1,:,:)
   Csoil_out(1,:,:) = Csoil_out(nyear-1,:,:)

   ;;;;;;;;;;;;;;;
   npp_in_avg=dim_avg_n(npp_in(0:200,:,:),0)
   Clabloss_in_avg=dim_avg_n(Clabloss_in(0:200,:,:),0)

   do iyear = 1, nyear-1

   Cplant_out(iyear,:,:) = ((npp_in_avg - Clabloss_in_avg) - epsilon * A21_avg * Cplant_out(iyear-1,:,:)) \
                        + Cplant_out(iyear-1,:,:)
   Clitter_out(iyear,:,:) = (epsilon * A21_avg * Cplant_out(iyear-1,:,:) - epsilon * A32_avg * Clitter_out(iyear-1,:,:) \
                        - epsilon * A22_avg * Clitter_out(iyear-1,:,:))+ Clitter_out(iyear-1,:,:)
   Csoil_out(iyear,:,:) = (epsilon * A32_avg * Clitter_out(iyear-1,:,:) - epsilon * A33_avg * Csoil_out(iyear-1,:,:)) + Csoil_out(iyear-1,:,:)
   nbp_out(iyear,:,:) = ( Cplant_out(iyear,:,:) + Clitter_out(iyear,:,:) + Csoil_out(iyear,:,:) )  \
                     - ( Cplant_out(iyear-1,:,:) + Clitter_out(iyear-1,:,:) + Csoil_out(iyear-1,:,:) )
   rh_out(iyear,:,:) = ((npp_in_avg - Clabloss_in_avg) - nbp_out(iyear,:,:))
   end do

   printVarSummary(nbp_out)
   printVarSummary(rh_out)

   Cplant_out(0,:,:) = Cplant_out(nyear-1,:,:)
   Clitter_out(0,:,:) = Clitter_out(nyear-1,:,:)
   Csoil_out(0,:,:) = Csoil_out(nyear-1,:,:)

   Cplant_out(1,:,:) = Cplant_out(nyear-1,:,:)
   Clitter_out(1,:,:) = Clitter_out(nyear-1,:,:)
   Csoil_out(1,:,:) = Csoil_out(nyear-1,:,:)

   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;;; running
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   do iyear = 1, nyear-1
   Cplant_out(iyear,:,:) = ((npp_in(iyear,:,:) - Clabloss_in(iyear,:,:)) - epsilon * A21_avg * Cplant_out(iyear-1,:,:)) \
                        + Cplant_out(iyear-1,:,:)
   Clitter_out(iyear,:,:) = (epsilon * A21_avg * Cplant_out(iyear-1,:,:) - epsilon * A32_avg * Clitter_out(iyear-1,:,:) \
                        - epsilon * A22_avg * Clitter_out(iyear-1,:,:)) + Clitter_out(iyear-1,:,:)
   Csoil_out(iyear,:,:) = (epsilon * A32_avg * Clitter_out(iyear-1,:,:) - epsilon * A33_avg * Csoil_out(iyear-1,:,:)) + Csoil_out(iyear-1,:,:)
   nbp_out(iyear,:,:) = ( Cplant_out(iyear,:,:) + Clitter_out(iyear,:,:) + Csoil_out(iyear,:,:) )  \
                     - ( Cplant_out(iyear-1,:,:) + Clitter_out(iyear-1,:,:) + Csoil_out(iyear-1,:,:) )
   rh_out(iyear,:,:) = ((npp_in(iyear,:,:) - Clabloss_in(iyear,:,:)) - nbp_out(iyear,:,:))
   end do

   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;


   tau_residence_in=rh_in
   tau_residence_in=0.
   tau_residence_out=rh_out
   tau_residence_out=0.

   rh_in_tau=where(rh_in.gt.1.e-6,rh_in,rh_in@_FillValue)
   rh_out_tau=where(rh_out.gt.1.e-6,rh_out,rh_out@_FillValue)
 
   tau_residence_in=Clitter_in/rh_in_tau
   tau_residence_out=Clitter_out/rh_out_tau

  ; do iyear = 0, nyear-1
  ; tau_residence_in(iyear,:,:)=(Clitter_in(iyear,:,:)+Csoil_in(iyear,:,:))/rh_in_tau(iyear,:,:)
  ; tau_residence_out(iyear,:,:)=(Clitter_out(iyear,:,:)+Csoil_out(iyear,:,:))/rh_out_tau(iyear,:,:)
  ; end do

  A21_avg!0="latitude"
  A21_avg!1="longitude"
  A21_avg&latitude=latitude
  A21_avg&longitude=longitude

  A22_avg!0="latitude"
  A22_avg!1="longitude"
  A22_avg&latitude=latitude
  A22_avg&longitude=longitude

  A32_avg!0="latitude"
  A32_avg!1="longitude"
  A32_avg&latitude=latitude
  A32_avg&longitude=longitude

  A33_avg!0="latitude"
  A33_avg!1="longitude"
  A33_avg&latitude=latitude
  A33_avg&longitude=longitude

  Cplant_out@units="PgC"
  Clitter_out@units="PgC"
  Csoil_out@units="PgC"
  tau_residence_out@units="1/year"

   system("rm -f "+FilePathOut+FileNameOut+"all.nc")
   fout=addfile(FilePathOut+FileNameOut+"all.nc","c")
   fout->NBP  = nbp_out
   fout->HR   = rh_out
   fout->NPP  = npp_in
   fout->TOTVEGC = Cplant_out
   fout->LITTERC = Clitter_out
   fout->SOILC = Csoil_out
   fout->turnover_rate  = tau_residence_out
   fout->A21=A21_avg
   fout->A22=A22_avg
   fout->A32=A32_avg
   fout->A33=A33_avg

   npp_in_global=dim_sum_n(npp_in,(/1,2/))
   nbp_in_global=dim_sum_n(nbp_in,(/1,2/))
   rh_in_global=dim_sum_n(rh_in,(/1,2/))
   nbp_out_global=dim_sum_n(nbp_out,(/1,2/))
   rh_out_global=dim_sum_n(rh_out,(/1,2/))
   tau_residence_in_global=dim_avg_n(tau_residence_in,(/1,2/))
   tau_residence_out_global=dim_avg_n(tau_residence_out,(/1,2/))

  ; print(tau_residence_in_global)
  ; print(tau_residence_out_global)
  ; print(nbp_out_global)
   print(rh_in_global)
   print(rh_out_global)  

 
   asciiwrite("Total_annual_NBP_Box_Model_"+runner(j)+"_from_"+Startyr+"-"+Endyr+"_global.txt",nbp_out_global)
   asciiwrite("Total_annual_tau_CESM1_CLM4_"+runner(j)+"_from_"+Startyr+"-"+Endyr+"_global.txt",tau_residence_in_global)
   asciiwrite("Total_annual_tau_Box_Model_"+runner(j)+"_from_"+Startyr+"-"+Endyr+"_global.txt",tau_residence_out_global)
   asciiwrite("Total_annual_HR_Box_Model_"+runner(j)+"_from_"+Startyr+"-"+Endyr+"_global.txt",rh_out_global)
   asciiwrite("Total_annual_NPP_Box_Model_"+runner(j)+"_from_"+Startyr+"-"+Endyr+"_global.txt",npp_in_global)

   asciiwrite("Total_annual_NBP_CESM1_CLM4_Box_Model_"+runner(j)+"_from_"+Startyr+"-"+Endyr+"_global.txt",nbp_in_global-nbp_out_global)
   asciiwrite("Total_annual_HR_CESM1_CLM4_Box_Model_"+runner(j)+"_from_"+Startyr+"-"+Endyr+"_global.txt",rh_in_global-rh_out_global)



   delete(nbp_in)
   delete(npp_in)
   delete(rh_in)
   delete(nbp_out)
   delete(rh_out)
   delete(tau_residence_in)
   delete(tau_residence_out)


 end do

end


