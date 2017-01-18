!> @brief Crop subroutines
!! @details
!! @author 
!! @date 

module crops

use real_precision
use pft_parameters
use site_parameters
use func
use light_methods
implicit none

contains

!**********************************************************************!
!                                                                      !
!                     seasonality :: crops                             !
!                     ---------------------                            !
!                                                                      !
! subroutine seasonality                                               !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Crops Seasonality
!! @details Following Waha et al 2012 (published version of LGJ 
!! van Bussel's thesis (2011)
!! chapter 5), calculate variation coefficients for seasonality, which may 
!! control sowing dates.  Then, get the sowing day based on precipitation or 
!! summer/winter temperature thresholds for each crop.
!
!! Following van Bussel et al 2015, Appendix S1 (van Bussel thesis 2011 fig 6.1a
!! on page 119), we also calculate the optimal photoperiod for each crop based on
!! longest and possibly shortest daylengths of the year.  For example, the number
!! of hours in the summer solstice is the optimum photoperiod for wheat, but for
!! maize it is sensitive to latitude: longest*(1-1/(longest-shortest)). Without
!! hardcoded crop types, we pass a parameter "pscale" (0 or 1) for each crop to 
!! tell us what to do: longest*(1-pscale/(longest-shortest)) 
!! 
!! Also following van Bussel et al 2015, Appendix S1 (van Bussel thesis 2011 
!! chapter 5), calculate vernalisation requirements based on temperature.
!!
!! REFERENCES:
!!
!! Waha et al 2012, Global Ecology and Biogeography 21:247-259
!! van Bussel et al 2015, Global Ecology and Biogeography (earlyview)
!! Wang and Engel 1998 (Agri Sys 58:1-24)
!!
!! tmp(12,31)        IN  Daily Temperature
!! prc(12,31)        IN  Daily Precipitation
!! cld(12)           IN  Monthly Cloud cover
!! nft               IN  Number of plant functional types
!! ssp%emnthtmp(12)      (20-year exp-mean monthly temperatures)
!! ssp%emnthprc(12)      (20-year exp-mean monthly precipitation)
!! museas                Annual average of exponential weighted *
!! seastmp               Seasonal variation coefficient of temperature
!! seasprc               Seasonal variation coefficient of precipitation
!! iseas                 Seasonality type.1 for precipitation,2 for temperature
!! mdoy(12)              First day of month in Julian days
!! mmid(12)              Midday of month in Julian days
!! hrs(12)               Daylight hours for the midday of each month
!! mnth                  Month counter
!! ssp%lat               Site latitude
!! q(12)                 Monthly photon flux density
!! pet(12)               Monthly ratio precip/(potential evapotranspiration)
!! mnthhum(12)       IN  (mean monthly humidity)
!! sowday(*)         OUT Sow day per pft
!! msow                  Sowing month
!! daysiny               Days in a year
!! mcoldest              Coldest month
!! icoldest              Midday of coldest month
!! sumtsow(*)   input Minimum temperature thresholds for sowing
!! wintsow(*)   input Maximum temperature thresholds for sowing
!! ndays(12)    input Number of days per month
!! tmp(12,31)   input Mean daily temperatures (C)
!! prc(12,31)   input Mean daily precipitation (mm)
!! pscale(*)    input Flags to reduce optimal photoperiods for each PFT
!! popt(*)     output Optimal photoperiods for each PFT
!! vernmsat(*)  input Max vernalisation requirement per PFT (days/month)
!! vernt2(*)    input Min temperature for optimal vernalisation per PFT (C)
!! vernt3(*)    input Max temperature for optimal vernalisation per PFT (C)
!! vernsat(*)  output Max vernalisation requirement per PFT (days)
!! verbose      input Logical: write out values or not
!!
!! @author LLT,EPK 
!! @date Oct 2016
!----------------------------------------------------------------------!
subroutine seasonality(tmp,prc,cld,thty_dys,nft,year)
!**********************************************************************!
implicit none
real(dp) :: tmp(12,31),prc(12,31),cld(12)
integer  :: thty_dys,nft,year

real(dp) :: museas,seastmp,seasprc,hrs(12)
integer  :: iseas,mdoy(12),mmid(12),mnth,day,ft,k,daysiny
integer  :: msow,mcoldest,icoldest,nwarm,mcold,ntot,nvern,i,j,m
integer  :: nydays,slen(12)
real(dp) :: fv,vdays,vegphu(12),repphu(12),totveg,wft,wfp,vrat
real(dp) :: totphu(12),targ,totrep
logical  :: dead
real(dp) :: qdir,qdiff,q(12),pet(12) ! internal, for getwet
!----------------------------------------------------------------------!

! Days in a year
  daysiny=0
  do mnth=1,12
    daysiny=daysiny+int(no_days(year,mnth,thty_dys))
  enddo
  
! See Page 94 of van Bussel's thesis, section 2.2.1, eqns 5.1-5.3
! or page 249 of Waha et al 2012
! Calculates variation coefficient for temperature and precipitation
  museas=sum(ssp%emnthtmp)/12.0d0+273.15
  seastmp=sqrt(1.0d0/11.0d0* &
    sum((ssp%emnthtmp(:)+273.15-museas)*(ssp%emnthtmp(:)+273.15-museas)))/museas
  museas=sum(ssp%emnthprc)/12.0d0
  seasprc=sqrt(1.0d0/11.0d0* &
  sum((ssp%emnthprc(:)-museas)*(ssp%emnthprc(:)-museas)))/museas

! Flag the seasonality type to use
  ssp%iseas=0
! Thresholds for seasonality come from van Bussel's thesis page 95
! or from Waha et al 2012 page 249

  if (seastmp.gt.0.01d0.and.seasprc.gt.0.4d0) then
! Only use one seasonality type, depending on whether 
! minimum mean monthly temperature exceeds 10C
    if (minval(ssp%emnthtmp).gt.10.0d0) then
      ssp%iseas=1 ! assume seasons controlled primarily by precipitation
    else
      ssp%iseas=2 ! assume seasons controlled primarily by temperature
    endif ! min(ssp%emnthtmp)>10C
  elseif (seastmp.gt.0.01d0) then
    ssp%iseas=2 ! seasons controlled by temperature
  elseif (seasprc.gt.0.4d0) then
    ssp%iseas=1 ! seasons controlled by precipitation
  endif ! check for combined temp and precip seasonality

! For each month, get day-of-year for the first day and mid-month day length
  mdoy(1)=1; mmid(1)=int(no_days(year,1,thty_dys)/2)
  hrs(1)=dayl(ssp%lat,15) ! Number of hours of daylight in the specified day
  call pfd(ssp%lat,15,hrs(1),cld(1),qdir,qdiff,q(1)) ! photon flux density
  do mnth=2,12
    mdoy(mnth)=mdoy(mnth-1)+no_days(year,mnth-1,thty_dys)
    mmid(mnth)=int(no_days(year,mnth,thty_dys)/2)+mdoy(mnth)-1
    hrs(mnth)=dayl(ssp%lat,mmid(mnth))
    call pfd(ssp%lat,mmid(mnth),hrs(mnth),cld(mnth),qdir,qdiff,q(mnth))
  enddo
  
  pft_tab(1:nft)%sowday=0 ! recalculate sowday each year 
  
! Follow van Bussel (section 2.2.2) or Waha et al p 249 to find sowing day
  if (ssp%iseas.eq.1) then ! If sowing depends on precipitation
         
! Run four-month sums of precip/PET ratios 
! mnth is an output variable holding the first month of the wet season
    call getwet(hrs,q,mnth,pet)
      
! Find the first day of the first wet-season month which is wet enough
! All crops will be sown on the same day in this case
    do day=1,no_days(year,mnth,thty_dys)
      if (pft_tab(3)%sowday.eq.0.and.prc(mnth,day).gt.0.1d0) &
        pft_tab(3:nft)%sowday=day-1+mdoy(mnth)
    enddo ! day=1,no_days(year,mnth,thty_dys)
  elseif (ssp%iseas.eq.2) then ! sowing depends on temperature
    pet(:)=1.0d0
    do ft=3,nft ! there will be crop-specific threshold temperatures
! Figures whether the crop will have a spring sow date
      if (pft_tab(ft)%sowthresh(2).gt.pft_tab(ft)%lethal(1)) then
        pft_tab(ft)%sowday=summerday(pft_tab(ft)%sowthresh(2),mmid)
      endif ! sowday,avtmp.ge.sumtsow
! Some crops have vernalization requirements; skip those that do not,
! usually crops have either summer or winter thresholds but not both
! (see van Bussel Table 5.1 on page 97, or Waha et al Table 1 page 250)
      if (pft_tab(ft)%sowthresh(1).ge.pft_tab(ft)%lethal(2)) cycle 
      if (pft_tab(ft)%sowthresh(1).gt.maxval(ssp%emnthtmp)) then  
        pft_tab(ft)%sowday=winterday(maxval(ssp%emnthtmp)-0.1,mmid)
      else
        pft_tab(ft)%sowday=winterday(pft_tab(ft)%sowthresh(1),mmid)
      endif
    enddo ! ft=1,nft
  endif ! seasonality
  
! Estimate vernalisation days
  do ft=3,nft
! Aseasonal climate; just sow on new years day as per van Bussel 2011 p 95 
    if (pft_tab(ft)%sowday.eq.0) pft_tab(ft)%sowday=1
    msow=1; mcoldest=1
! Get the sowing month and the coldest month
    do mnth=1,12
      if (pft_tab(ft)%sowday.ge.mdoy(mnth)) msow=mnth
      if (ssp%emnthtmp(mnth).lt.ssp%emnthtmp(mcoldest)) mcoldest=mnth
    enddo ! mnth=1,12
! Get the midday of the coldest month
    icoldest=mmid(mcoldest);
    if (ssp%iseas.lt.2) icoldest=0
    fv=1.0-pft_tab(ft)%croptype(2); vdays=0.0; dead=.FALSE.
    j=mdoy(msow)-1; vegphu(:)=0.0; repphu(:)=0.0
    nwarm=0; mcold=0; ntot=0; totveg=0.0; nvern=0

! Cycle through all 12 months starting with the sowing month
! k is the number of months added to the sowing month in each iteration
    do k=0,11 !6+5*croptype(2,ft)
      if (dead) cycle
      mnth=msow+k
      if (mnth.gt.12) mnth=mnth-12
      dead=(ssp%emnthtmp(mnth).le.pft_tab(ft)%lethal(1).or. &
        ssp%emnthtmp(mnth).ge.pft_tab(ft)%lethal(2))
      if (dead) cycle
      if (mcold.gt.0) cycle
! Calculate the mean daily heat units for this month (tmp*wft*wfp)
! wft is the fraction of the heat units allowed per day this month
! wfp is the fraction of the maximum photoperiod
      call wangengel(pft_tab(ft)%cardinal(1),pft_tab(ft)%cardinal(2), &
        pft_tab(ft)%cardinal(3),ssp%emnthtmp(mnth),pft_tab(ft)%croptype(1), &
        pft_tab(ft)%photoperiod(1),pft_tab(ft)%photoperiod(2),hrs(mnth),wft,wfp)
! For crops that require vernalization.If the crop requires vernalization but
! it has been met,meaning fv.gt.0.95,skip.
      if (fv.lt.0.95d0) then
        do i=1,no_days(year,mnth,thty_dys)
          j=j+1
          if (j.lt.pft_tab(ft)%sowday) cycle
          if (fv.lt.0.95d0) call streck(pft_tab(ft)%cardinal(4), &
            pft_tab(ft)%cardinal(5),pft_tab(ft)%cardinal(6), &
            ssp%emnthtmp(mnth),pft_tab(ft)%croptype(1),pft_tab(ft)%photoperiod(3), &
            pft_tab(ft)%photoperiod(4),dayl(ssp%lat,j),vdays,fv)
            totveg=totveg+ssp%emnthtmp(mnth)*fv*wft*wfp
            cycle
        enddo ! i=1,no_days(year,mnth,thty_dys)
        nvern=nvern+1
      else
        totveg=totveg+ssp%emnthtmp(mnth)*no_days(year,mnth,thty_dys)*wft*wfp
        call wangengel(pft_tab(ft)%cardinal(7),pft_tab(ft)%cardinal(8), &
          pft_tab(ft)%cardinal(9),ssp%emnthtmp(mnth),pft_tab(ft)%croptype(1), &
          pft_tab(ft)%photoperiod(5),pft_tab(ft)%photoperiod(6),hrs(mnth),wft,wfp)
        repphu(mnth)=repphu(mnth)+ssp%emnthtmp(mnth)*no_days(year,mnth,thty_dys)*wft*wfp
      endif
! The cumulative vegetative heat units achieved this month
      vegphu(mnth)=totveg
! If vernalisation complete, then we need to count the number of months where we
! accumulate reproductive PHU.  These months need to be contiguous.
! It will keed adding to nwarm for as long as repphu is positive.
! When pepphu becomes zero or negative,it will mean that the growing season ended
! and mcold will attain a value diff than zero which will stop PHU from accumulating
! based on a conditional above.
      if (fv.gt.0.95d0) then
        if (nwarm.gt.0.and.repphu(mnth).le.0.0d0) mcold=mnth
        if (mcold.eq.0.and.repphu(mnth).gt.0.0d0) nwarm=nwarm+1
      endif ! fv.gt.0.995d0
      if (mcold.eq.0) ntot=ntot+1
    enddo ! k=0,11
    
    pft_tab(ft)%cropgdd(1)=pft_tab(ft)%croprange(1)
    if (ssp%iseas.eq.2) then ! Temperature controlled
! Here, I apply a squared cosine function to get the GDD rather than a
! quadratic as Bondeau et al (2007) did.  The cosine provides a smoothly
! varying function which can easily be shifted depending on when the coldest
! month occurs, avoiding any need for hardcoding time windows.
      if (totveg.gt.pft_tab(ft)%cropgdd(1).and. &
        pft_tab(ft)%croprange(1).lt.pft_tab(ft)%croprange(2)) then
        nydays=daysiny
        pft_tab(ft)%cropgdd(1)=pft_tab(ft)%croprange(2)
        pft_tab(ft)%cropgdd(2)=pft_tab(ft)%croprange(4)
        vrat=pft_tab(ft)%sowday-icoldest
        pft_tab(ft)%cropgdd(1)=pft_tab(ft)%croprange(1)+ &
          (pft_tab(ft)%croprange(2)-pft_tab(ft)%croprange(1))* &!*fv
            (cos(vrat*3.14159d0/nydays)**2)
      endif ! (croprange(1,ft).lt.croprange(2,ft))
    else
! Estimate using cumulative sum for vegetative growth (vegphu), max allowed
! months for reproductive growth (3), and phenological index where senescence
! begins (cropgdd(5,ft))
! croptype(2,ft) is zero if no vernalisation required
! croprange(1,ft) is min PHU for maturity; ntot is total months accumulating PHU
! Restrict ntot to 3 if the crop does not require vernalisation
      if (pft_tab(ft)%croptype(2).eq.0 &
       .and.maxval(vegphu).gt.pft_tab(ft)%croprange(1)) ntot=min(3,ntot)
      totrep=0;
      i=0; totphu(:)=0.0d0; slen(:)=0
! Start with the month after vernalisation is complete
! nvern are the months required for vernalization
      do k=nvern+1,ntot-1
! cropgdd(1,ft) is initialised to be croprange(1,ft) but may change
! croprange(1,ft) is the min GDD required
        if (pft_tab(ft)%cropgdd(1).gt.pft_tab(ft)%croprange(1)) cycle
! Get the current calendar month, after vernalising
        mnth=msow+k
        if (mnth.gt.12) mnth=mnth-12
! targ is the target PHU for maturity, when we can harvest the plant
! targ is always greater than vegphu because we add repphu to it
        targ=vegphu(mnth)/pft_tab(ft)%cropphen(5)
        totphu(mnth)=vegphu(mnth)
! loop over reproductive (fruiting) months
        do j=k+1,ntot
! We accumulate reproductive PHU only until we achieve the PHU for harvest 
          if (totphu(mnth).ge.targ) cycle
          m=j+msow
          if (m.gt.12) m=m-12
! Add monthly reproductive PHU to total cumulative PHU
          totphu(mnth)=max(targ,totphu(mnth)+repphu(m))
          slen(mnth)=j-nvern
        enddo ! do j=k+1,ntot
        if (totphu(mnth).le.pft_tab(ft)%croprange(2).and. &
          totphu(mnth).gt.pft_tab(ft)%croprange(1)) then
          pft_tab(ft)%cropgdd(1)=totphu(mnth)
          totrep=totrep+totphu(mnth) !*4.0d0/(1+slen(mnth))
          if (slen(mnth).gt.0) i=i+1 !slen(mnth)+1
        elseif (totphu(mnth).gt.pft_tab(ft)%croprange(2)) then
          pft_tab(ft)%cropgdd(1)=pft_tab(ft)%croprange(2)
        endif
      enddo 
      if (i.gt.0) pft_tab(ft)%cropgdd(1)=totrep/i !*i/4.0d0
    endif ! iseas.eq.2, controlled by temperature or otherwise
  enddo ! ft=3,nft

end subroutine seasonality


!**********************************************************************!
!                                                                      !
!                         getwet :: crops                              !
!                     ---------------------                            !
!                                                                      !
! subroutine getwet                                                    !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief 
!! @details Calculate ratios of precipitation to evt
!! van Bussel (2011) discusses determination of start of growing period
!! in her chapter 5, section 2.2.2, and the sowing day should be the first
!! day where precip>0.1mm of that wet month.
!! Calculate 4-month sums of prec/evt and find the highest.
!! Output the starting month of the highest.
!! We depart from van Bussel in that we do not use Priestley-Taylor PET, 
!! but instead lift Penman-Montheith code from of dolyday.
!!
!! REFERENCES:
!!
!! Waha et al 2012, Global Ecology and Biogeography 21:247-259
!! van Bussel et al 2015, Global Ecology and Biogeography (earlyview)
!!
!! hrs(12)      IN    Daylight hours for the midday of each month
!! q(12)        IN    Monthly photon flux density
!! wet          OUT   Month index for start of wet season
!! pet(12)      OUT   Monthly ratio precip/(potential evapotranspiration)
!! 
!!
!! @author LLT,EPK 
!! @date Oct 2016
!----------------------------------------------------------------------!
subroutine getwet(hrs,q,wet,pet)
!**********************************************************************!    
implicit none
real(dp) :: hrs(12),q(12),pet(12)
integer  :: wet

real(dp) :: ht,windspeed,canga,t,rh,rn
real(dp) :: svp,vpd,lam,rho,s,gam,ee,eemm,petsum(12)
integer :: m    
!*----------------------------------------------------------------------*

  ht=25.1338 ! always has this value in sdgvm0.f
  windspeed= 5.0d0 ! in m/s
  canga = 0.168d0*windspeed/log((200.0d0 - &
    0.7d0*ht)/(0.1d0*ht))**2

  pet(:)=0
  do m=1,12
    if (ssp%emnthprc(m).le.0.0d0) cycle     
    t=ssp%emnthtmp(m); rh=ssp%mnthhum(m);
    rn = 0.96d0*(q(m)*1000000.0d0/4.0d0 + 208.0d0 + 6.0d0*t)
    rn = rn*0.52d0

!*----------------------------------------------------------------------*
!* Penman-Monteith equation for evapotranspiration.                     *
!* Units of ET = W/m2  (CHECK) N.B. Conductances in moles.              *
!*----------------------------------------------------------------------*
    svp = 6.108d0*exp((17.269d0*t)/(237.3d0 + t))*100.0d0
    vpd = (1.0d0 - rh/100.0d0)*svp
    lam = 2500.0d0 - 2.367d0*t
    rho = 1288.4d0 - 4.103d0*t
    s = 48.7d0*exp(0.0532d0*t)
    gam = 101325.0d0*1.012d0/(0.622d0*lam)

    ee = (s*rn + rho*1.012d0*canga*vpd)/(s + gam)
    eemm = (ee*3600.0d0*hrs(m))/(lam*1000.0d0)

    pet(m)=ssp%emnthprc(m)/eemm !Ratio

  enddo ! m=1,12

! Get the four-month sums
  petsum(1)=sum(pet(1:4))
  wet=1
  do m=2,9
    petsum(m)=sum(pet(m:m+3))
    if (petsum(m).gt.petsum(m-1)) wet=m
  enddo ! m=1,12   
  petsum(10)=sum(pet(10:12))+pet(1)
  petsum(11)=sum(pet(11:12))+sum(pet(1:2))
  petsum(11)=pet(12)+sum(pet(1:3))
  do m=10,12
    if (petsum(m).gt.petsum(m-1)) wet=m
  enddo ! m=1,12   
 
  return
end subroutine getwet


!**********************************************************************!
!                                                                      !
!                        summerday :: crops                            !
!                     ---------------------                            !
!                                                                      !
! function summerday                                                   !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Find the start day of summer as defined by a threshold temperature
!! @details van Bussel (2011, thesis) discusses determination of sowing dates
!! in her chapter 5, section 2.2.2; here we find the day number when
!! summer begins (temperature>=threshold).  We linearly interpolate
!! between the average monthly temperatures to get the mean day.
!!
!!
!! thresh       IN threshold temperature for summer
!! mmid(12)     IN Midday of month in Julian days
!!
!!
!! @author LLT,EPK 
!! @date Oct 2016
!----------------------------------------------------------------------!
integer function summerday(thresh,mmid)
!**********************************************************************!
implicit none
real(dp) :: thresh
integer  :: mmid(12)
integer  :: lasttmp,lastm,lastday,m,thistmp


  summerday=1
! If all monthly temperatures below or above threshold,return
  if (minval(ssp%emnthtmp(:)).ge.thresh) return
  if (maxval(ssp%emnthtmp(:)).lt.thresh) return
  lasttmp=0; lastm=12; lastday=-mmid(1)
  if (ssp%emnthtmp(12).ge.thresh) lasttmp=1
  do m=1,12
    thistmp=0
    if (ssp%emnthtmp(m).ge.thresh) thistmp=1
    if (thistmp.gt.lasttmp) summerday=int((thresh-ssp%emnthtmp(lastm))* &
      (mmid(m)-lastday)/(ssp%emnthtmp(m)-ssp%emnthtmp(lastm)))+lastday
    lasttmp=thistmp; lastm=m; lastday=mmid(m)
  enddo
  if (summerday.le.0) summerday=mmid(12)-summerday
      
  return 
end function summerday


!**********************************************************************!
!                                                                      !
!                       winterday :: crops                             !
!                     ---------------------                            !
!                                                                      !
! function winterday                                                   !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Find the start day of winter as defined by a threshold temperature
!! @details 
!!
!!
!!
!! thresh       IN threshold temperature for winter
!! mmid(12)     IN Midday of month in Julian days
!!
!!
!! @author LLT,EPK 
!! @date Oct 2016
!----------------------------------------------------------------------!
integer function winterday(thresh,mmid)
!**********************************************************************!
implicit none
real(dp) :: thresh
integer  :: mmid(12)
integer  :: lasttmp,lastm,lastday,m,thistmp

  winterday=1
  if (minval(ssp%emnthtmp(:)).ge.thresh) return
  if (maxval(ssp%emnthtmp(:)).lt.thresh) return
  lasttmp=0; lastm=12; lastday=-mmid(1)
  if (ssp%emnthtmp(12).le.thresh) lasttmp=1
  do m=1,12
    thistmp=0
    if (ssp%emnthtmp(m).le.thresh) thistmp=1
    if (thistmp.gt.lasttmp) winterday=int((thresh-ssp%emnthtmp(lastm))* &
      (mmid(m)-lastday)/(ssp%emnthtmp(m)-ssp%emnthtmp(lastm)))+lastday
    lasttmp=thistmp; lastm=m; lastday=mmid(m)
  enddo
  if (winterday.le.0) winterday=mmid(12)-winterday
      
  return
  end function winterday


!**********************************************************************!
!                                                                      !
!                       wangengel :: crops                             !
!                     ---------------------                            !
!                                                                      !
! subroutine wangengel                                                 !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Get temperature and photoperiod response function values for crops
!! @details These functions are from Wang and Engel 1998 (Agri Sys 58:1-24).
!! The sum of the product of the temperature and photoperiod functions gives 
!! the effective physiological days.
!!
!!
!!
!! @author LLT,EPK 
!! @date Oct 2016
!----------------------------------------------------------------------!
subroutine wangengel(tmin,topt,tmax,tmp,dtype,popt,pcrit,p,ft,fp)
!**********************************************************************!
implicit none
real(dp) :: tmin,topt,tmax ! Cardinal temperatures
real(dp) :: popt,pcrit,p ! Optimal, critical and current photoperiod
real(dp) :: tmp     ! today's mean temperature, input
real(dp) :: ft,fp   ! temperature and photoperiod functions, output
real(dp) :: alpha  ! exponent in temperature response function, internal
real(dp) :: t1,t2 ! internal
real(dp) :: dtype ! 1 for long-day, -1 for short-day, 0 for day-neutral plants

ft=0.0
! Photoperiod effect.Eq. 11,12,13 from Wang & Engel
fp=max(0.0,1-exp(-dtype*4.0*(p-pcrit)/abs(popt-pcrit)))
if (tmax.gt.tmin) then
! NB ft will be zero for tmp outside the range tmin to tmax inclusive
  if (tmp.le.tmax.and.tmp.ge.tmin) then
! alpha value from Eq.6 of Wang & Engel
    alpha=log(2.0)/log((tmax-tmin)/(topt-tmin))
    t1=(topt-tmin)**alpha
    t2=(tmp-tmin)**alpha
    ft=(2*t2*t1-(t2**2))/(t1**2)
  endif ! tmp.le.tmax.and.tmp.ge.tmin
else
endif ! vtmax.gt.vtmin

return
  
end subroutine wangengel


!**********************************************************************!
!                                                                      !
!                       streck :: crops                                !
!                     ---------------------                            !
!                                                                      !
! subroutine streck                                                    !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Get vernalisation factors for crops
!! @details  Following Streck et al 2003, Ag For Meteor 115:139-150
!! Equations 3 and 4 (temperature function), 5 (photoperiod function)
!! These functions are from Wang and Engel 1998 (Agri Sys 58:1-24).
!! In some papers including Streck et al's first one (below) the temperature
!! function has a typographical error (a minus sign is missing).
!! The sum of the product of the temperature and photoperiod functions gives 
!! the effective vernalised days VD
!! which one can use to get the vernalization response function (Eq 10)
!! or Eq 2 in Streck et al 2003 Agron J 95:155-159
!!
!! Streck et al's vernalisation function and Wang and Engel's functions
!! are used in a number of models:
!!
!! JULES-SUCROS land surface model 
!!    (van den Hoof et al 2011 Ag For Meteor 151:137-153)
!! DANUBIA crop growth model (Lenz-Wiedemann et al 2010 Ecol Model 221:314-329)
!! FROSTOL wheat model (Bergjord et al 2008 Eur J Agr 28:321-330)
!! SPACSYS C and N cycling model (Bingham and Wu 2011 Eur J Agr 34:181-189) 
!! CANDY-PLUS C, N and biomass model (Kruger et al 2013 Vadose Zone J 12)
!! Liu's pasture legume model (Liu 2007 Field Crops Res 101:331-342)
!! 
!! vtmin,vtopt,vtmax  input Cardinal temperatures for vernalisation
!! tmp                input Current mean daily temperature
!! popt,pcrit         input Optimal and critical photoperiod
!! p                  input Current photoperiod
!! dtype              input -1: short-day, 1: long-day, 0: day-neutral
!! vdays       input/output Accumulated vernalisation days
!! fv                output Current vernalisation day
!!
!!
!!
!! @author LLT,EPK 
!! @date Oct 2016
!----------------------------------------------------------------------!
subroutine streck(vtmin,vtopt,vtmax,tmp,dtype,popt,pcrit,p,vdays,fv)
!**********************************************************************!
implicit none
real(dp) :: vtmin,vtopt,vtmax
real(dp) :: tmp     ! today's mean temperature, input
real(dp) :: vdays   ! accumulated effective vernalisation days, input/output
real(dp) :: valpha  ! exponent in temperature response function, internal
real(dp) :: dtype ! 1 for long-day, -1 for short-day, 0 for day-neutral plants
real(dp) :: popt,pcrit ! photoperiod parameters, input
real(dp) :: p ! today's photoperiod, input
real(dp) :: vd5,ft,fp   ! internal
real(dp) :: fv   ! vernalisation response function, output

if (vtmax.gt.vtmin) then
! NB increment will be zero for tmp outside the range vtmin to vtmax inclusive
! therefore the heat sum vdays will not be incremented and fv will not change
  if (tmp.le.vtmax.and.tmp.ge.vtmin) then
    call wangengel(vtmin,vtopt,vtmax,tmp,dtype,popt,pcrit,p,ft,fp)
    vdays=vdays+ft*fp
    vd5=vdays**5
    fv=(vd5)/((22.5**5)+vd5)
  endif ! tmp.le.vtmax.and.tmp.ge.vtmin
else
  vdays=0.0d0; fv=1.0d0
endif ! vtmax.gt.vtmin
 
return
end subroutine streck

!**********************************************************************!
!                                                                      !
!                       crop_outputs :: crops                          !
!                     ---------------------                            !
!                                                                      !
! subroutine crop_outputs                                              !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Output files for crops
!! @details  
!!
!! @author EPK 
!! @date Dec 2016
!----------------------------------------------------------------------!
subroutine crop_outputs(stoutput,nft,sl)
!**********************************************************************!
implicit none

character :: stoutput*1000
integer :: nft,ft,sl,coun

coun=0

do ft=3,nft
  if(pft_tab(ft)%phen.NE.3) cycle
  select case (sl)
  case (0)
    coun=coun+1
    open(500+coun,file=stoutput(1:blank(stoutput))//'/'// &
      trim(pft_tab(ft)%tag)//'sow.dat')
    coun=coun+1
    open(500+coun,file=stoutput(1:blank(stoutput))//'/'// &
      trim(pft_tab(ft)%tag)//'gdd.dat')
    coun=coun+1
    open(500+coun,file=stoutput(1:blank(stoutput))//'/'// &
      'site'//'seas.dat')
  case (1)
    coun=coun+1
    close(500+coun)
    coun=coun+1
    close(500+coun)
    coun=coun+1
    close(500+coun)
  case (2)
    coun=coun+1
    WRITE(500+coun,'(f7.3,f9.3)',advance='NO') ssp%lat,ssp%lon
    coun=coun+1
    WRITE(500+coun,'(f7.3,f9.3)',advance='NO') ssp%lat,ssp%lon
    coun=coun+1
    WRITE(500+coun,'(f7.3,f9.3)',advance='NO') ssp%lat,ssp%lon
  case (3)
    coun=coun+1
    WRITE(500+coun,'(i4)',advance='NO') pft_tab(ft)%sowday
    coun=coun+1
    WRITE(500+coun,'(i5)',advance='NO') pft_tab(ft)%cropgdd(1)
    coun=coun+1
    WRITE(500+coun,'(i2)',advance='NO') ssp%iseas
  case (4)
    coun=coun+1
    WRITE(500+coun,*)
    coun=coun+1
    WRITE(500+coun,*)
    coun=coun+1
    WRITE(500+coun,*)
  end select

enddo
return
end subroutine crop_outputs

!**********************************************************************!
!                                                                      !
!                      fert_crops :: crops                             !
!                     ---------------------                            !
!                                                                      !
! subroutine fert_crops()                                              !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Fertiliser usage for crops
!! @details  
!!
!! @author LLT<EPK
!! @date Jan 2017
!----------------------------------------------------------------------!
subroutine fert_crops(nft)
!**********************************************************************!
implicit none

integer nft,ft

DO ft=3,nft
  IF(pft_tab(ft)%phen.NE.3) cycle
  pft_tab(ft)%optlai=pft_tab(ft)%cropphen(1)  
  ! cropphen gives the optimal LAI range: cropphen(1) is without fertiliser
  ! cropphen(2) is with maximum fertiliser
  ! if this crop doesn't normally get fertiliser, then cropphen(2)=cropphen(1)
  IF(pft_tab(ft)%cropphen(2).GT.pft_tab(ft)%cropphen(1)) THEN
    ! Fertiliser usage in g/m2/y,This should be read from file,here placeholder
    pft_tab(ft)%nfert=100.
    ! get optimal LAI by adding log fert usage (kg/ha) to the nonfertilised LAI
    ! kg/ha = 1000g/10000m2 so 10 x g/m2 = kg/ha = g/10m2
    ! Make sure optlai remains in the range specified by cropphen
    pft_tab(ft)%optlai=max(log10(10.d0*pft_tab(ft)%nfert),0.0d0)+ &
      pft_tab(ft)%cropphen(1) 
    pft_tab(ft)%optlai=min(pft_tab(ft)%cropphen(2),pft_tab(ft)%optlai)
  ENDIF
ENDDO


end subroutine fert_crops


end module crops

