!> @brief Methods which control the phenology, carbon allocation, and
!! the aging and updating of the active carbon stem, root, leaf
!! compartments. Along with suma which is also compartmentalised in a similar manner.
!! @details The active carbon pools of stem root and leaf, as well as suma
!! are compartmentalised into bins of an integral number of days long eg.
!! 'lai_comp_length' in the 'dims' module.
!!
!! @todo Basicly the compartments are behaving like a linked list but
!! using more memory. More memoty efficient to use linked lists.
!!
!! @section Compartments
!! Compartments are created as and when required. Each compartment
!! is a structure containing a real variable 'val' and an integer value
!! 'age'. The variable 'age' is set equal to the Julian day that the
!! compartment was created, and remains constant at that value.
!! Once created a compartment will accumulate any carbon for the following
!! number of days given by 'var_comp_length'.
module phenology_methods


use misc_values
use system_state
use pft_parameters
use site_parameters
use dims
use tuning_parameters
use misc_parameters
use crops
use light_methods
implicit none

contains

!**********************************************************************!
!                   phenology :: phenological_methods                  !
!                   ---------------------------------                  !
!                                                                      !
! Used for grasses and crops. Lai being controlled by the available    !
! npp storage.                                                         !
!                                                                      !
! SUBROUTINE phenology(yield,laiinc)                                   !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief
!! @details
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine phenology(yield,laiinc)
!**********************************************************************!
real(dp) :: yield,laiinc
integer :: co
!----------------------------------------------------------------------!

co = ssp%cohort

if (pft(co)%phen/=0) then
  if (pft(co)%phen==1) then
    call PHENOLOGY1(yield,laiinc)
  elseif (pft(co)%phen==2) then
    call PHENOLOGY2(yield,laiinc)
  elseif (pft(co)%phen==3) then
    call PHENOLOGY3(yield,laiinc)
  else
    write(*,*) 'No phenology defined for ',pft(co)%phen,pft(co)%tag
  stop
  endif
endif

end subroutine phenology





!**********************************************************************!
!                   phenology1 :: phenological_methods                 !
!                   ----------------------------------                 !
!                                                                      !
! Used for trees. Lai being controlled by a weighted memory of npp     !
! storage over the previous 5 years.                                   !
!                                                                      !
! SUBROUTINE phenology1(yield,laiinc)                                    !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief
!! @details
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine phenology1(yield,laiinc)
!**********************************************************************!
real(dp), parameter :: maxlai = 10.9
real(dp) :: rlai,lairat,laiinc,wtwp,wtfc,bb0,bbmax,bblim,sslim,bbsum, &
 yield,smtrig,tsuma,stemfr,maint,ftagh,soil2g,soilw
integer :: lai,leafls,mnth,day,ij,i,bb,gs,bbgs,sssum
integer :: bb2bbmin,bb2bbmax,bbm,ssm,sss,ss,dsbb
integer :: ftdth,harvest,chill,dschill,co
!----------------------------------------------------------------------!

co = ssp%cohort

wtwp = ssp%wilt
day  = ssp%day
mnth = ssp%mnth

wtfc     = ssp%field
ftagh    = pft(co)%crop
ftdth    = pft(co)%d2h
stemfr   = ssv(co)%stemfr
bb       = ssv(co)%bb
ss       = ssv(co)%ss
bbgs     = ssv(co)%bbgs
dsbb     = ssv(co)%dsbb
chill    = ssv(co)%chill
dschill  = ssv(co)%dschill
leafls   = pft(co)%lls
bbm      = pft(co)%bbmem
bb0      = pft(co)%bb0
bbmax    = pft(co)%bbmax
bblim    = pft(co)%bblim
ssm      = pft(co)%senm
sss      = pft(co)%sens
sslim    = pft(co)%senlim
lairat   = pft(co)%lrat

bb2bbmin = 285
bb2bbmax = 390
gs = 60

rlai = ssv(co)%lai%tot(1)

yield = 0.0
harvest = 0

ss = 0

!----------------------------------------------------------------------!
! Check for chilling.                                                  !
!----------------------------------------------------------------------!
if (chill==0) then
  bbsum = 0.0
  do i=1,20
    if (ssp%tmem(i)<-5.0) &
  bbsum = bbsum + max(-10.0,ssp%tmem(i)+5.0)
  enddo
  if (bbsum<-100) then
    chill = 1
    dschill = 1
  endif
endif
if (chill==1) then
  dschill = dschill + 1
endif
if (dschill>260) then
  chill = 0
  dschill = 0
endif

!----------------------------------------------------------------------!
! Bubburst, if no budburst set, and sufficient soil moisture, then     !
! check for  budburst.                                                 !
!----------------------------------------------------------------------!
soilw = ssv(co)%soil_h2o(1) + ssv(co)%soil_h2o(2) + &
 ssv(co)%soil_h2o(3) + ssv(co)%soil_h2o(4)
soil2g = soilw/(ssp%soil_depth*10.0)

!if (((bb==0).and.(msv%mv_soil2g>wtwp+0.25*(wtfc-wtwp))).or. &
! ((dsbb>bb2bbmax).and.(msv%mv_soil2g>wtwp+0.1*(wtfc-wtwp)))) then


if (((bb==0).and.(soil2g>wtwp+0.25*(wtfc-wtwp))).or. &
 ((dsbb>bb2bbmax).and.(soil2g>wtwp+0.1*(wtfc-wtwp)))) then

  smtrig = 0.0
  do i=1,30
    smtrig = smtrig + ssv(co)%sm_trig(i)
  enddo
  if ((smtrig>30.0).or.(dsbb>bb2bbmax)) then

!----------------------------------------------------------------------!
! Check for budburst using degree days.                                !
!----------------------------------------------------------------------!
    bbsum = 0.0
    do i=1,bbm
      if (ssp%tmem(i)>bb0)  bbsum = bbsum + min(bbmax,ssp%tmem(i)-bb0)
    enddo

! If the bbsum exceeds the lim for budburst then
    if ((real(bbsum)>=real(bblim)*exp(-0.01*real(dschill))) &
 .or.(dsbb>bb2bbmax)) then
!----------------------------------------------------------------------!
! Adjust proportion of gpp going into stem production based on suma.   !
! This is essentially the LAI control.                                 !
!----------------------------------------------------------------------!
      tsuma = ssv(co)%suma%tot
      maint = max(1.0,(real(leafls)/360.0)*1.0)
      tsuma = tsuma - msv%mv_leafmol*1.25/maint*tgp%p_opt
      if (tsuma>msv%mv_leafmol*1.25/maint*tgp%p_opt) tsuma = msv%mv_leafmol*1.25/maint*tgp%p_opt
      if (tsuma<-msv%mv_leafmol*1.25/maint*tgp%p_opt) tsuma = -msv%mv_leafmol*1.25/maint*tgp%p_opt
      stemfr = stemfr + tsuma*tgp%p_laimem*12.0

      if (stemfr<120.0) stemfr = 120.0

!----------------------------------------------------------------------!
! Budburst occurance.                                                  !
!----------------------------------------------------------------------!
      bb = (mnth-1)*30 + day
      bbgs = 0
      dsbb = 0

      if (stemfr<0.8*ssv(co)%nppstore(1)) then
        ssv(co)%nppstore(3) = ssv(co)%nppstore(1) - stemfr
      else
        if (stemfr<0.75*ssv(co)%nppstore(1)) then
          ssv(co)%nppstore(3) = ssv(co)%nppstore(1) - stemfr
        else
          ssv(co)%nppstore(3) = ssv(co)%nppstore(1)*0.25
        endif
        stemfr = stemfr*0.8
      endif
      laiinc = (ssv(co)%nppstore(1) - 0.0*ssv(co)%nppstore(3))/msv%mv_leafmol/1.25/12.0
      ssv(co)%nppstore(2) = ssv(co)%nppstore(1)
    endif
  endif
endif

!----------------------------------------------------------------------!
! Compute length of growing season, and set to zero when equal to      !
! growing season.                                                      !
!----------------------------------------------------------------------!
if (bb>0)  bbgs = bbgs + 1
if (bbgs-gs>bb2bbmin) then
  bb = 0
  bbgs = 0
endif

if (dsbb < 500) dsbb = dsbb + 1

!----------------------------------------------------------------------!
! Set LAI increase.                                                    !
!----------------------------------------------------------------------!
if ((bb>0).and.(bbgs<gs).and.(ssv(co)%nppstore(1)>0.0)) then
  laiinc = lairat*(ssv(co)%nppstore(1) - 0.0*ssv(co)%nppstore(3))/msv%mv_leafmol/1.0/12.0*2.0
  if (rlai+laiinc>maxlai)  laiinc = maxlai - rlai
  if (rlai+laiinc>11.5)  laiinc = 11.5 - rlai
  if ((rlai>0).and.(ssv(co)%nppstore(1)<0.0)) laiinc = 0.0
else
  laiinc = 0.0
endif

!----------------------------------------------------------------------!
! Senescence, if rlai is greater than zero, compute senescence.        !
!----------------------------------------------------------------------!
if (rlai>1.0e-6) then
  if (abs(ftagh)<1.0e-6) then
    if (soil2g<wtwp*0.0) then
!----------------------------------------------------------------------!
! Senescence event due to soil moisture.                               !
!----------------------------------------------------------------------!
      laiinc = -rlai
      ss = day + (mnth - 1)*30
    elseif (bbgs>100) then
!----------------------------------------------------------------------!
! Check for senescence, senescence occurs there are 'sss' days colder  *
! than 'sslim' out of the last 'ssm' days.                             !
!----------------------------------------------------------------------!
      sssum = 0
      do i=1,ssm
        if (ssp%tmem(i)<sslim)  sssum = sssum + 1
      enddo
      if (sssum>=sss) then
!----------------------------------------------------------------------!
! Senescence event due to temperature.                                 !
!----------------------------------------------------------------------!
        laiinc =-rlai
        ss = day + (mnth - 1)*30
      endif
    endif 
  else
!----------------------------------------------------------------------!
! Crop senescence
!----------------------------------------------------------------------!
    if ((bbgs>leafls).or.(bbgs>ftdth)) then
      harvest = 1
      laiinc =-rlai
      ss = day + (mnth - 1)*30
    endif
  endif
endif

!----------------------------------------------------------------------!
ssv(co)%stemfr      = stemfr
ssv(co)%bb          = bb
ssv(co)%ss          = ss
ssv(co)%bbgs        = bbgs
ssv(co)%dsbb        = dsbb
ssv(co)%chill       = chill
ssv(co)%dschill     = dschill

end subroutine phenology1





!**********************************************************************!
!                   phenology2 :: phenological_methods                 !
!                   ----------------------------------                 !
!                                                                      !
! Used for trees. Lai being controlled by a weighted memory of npp     !
! storage over the previous 5 years.                                   !
!                                                                      !
! SUBROUTINE phenology2(yield,laiinc)                                  !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief
!! @details
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine phenology2(yield,laiinc)
!**********************************************************************!
real(dp), parameter :: maxlai = 10.9
real(dp) :: rlai,lairat,wtwp,laiinc,wtfc,bb0,bbmax,bblim,sslim,bbsum,yield,smtrig, &
 tsuma,stemfr,maint,soil2g,soilw
integer :: lai,leafls,ij,i,bb,gs,bbgs,sssum,bb2bbmin,bb2bbmax,bbm,ssm, &
 sss,ss,dsbb,chill,dschill,co,day,mnth
!----------------------------------------------------------------------!

co = ssp%cohort

wtwp = ssp%wilt
day = ssp%day
mnth = ssp%mnth

wtfc     = ssp%field
stemfr   = ssv(co)%stemfr
bb       = ssv(co)%bb
ss       = ssv(co)%ss
bbgs     = ssv(co)%bbgs
dsbb     = ssv(co)%dsbb
chill    = ssv(co)%chill
dschill  = ssv(co)%dschill
leafls   = pft(co)%lls
bbm      = pft(co)%bbmem
bb0      = pft(co)%bb0
bbmax    = pft(co)%bbmax
bblim    = pft(co)%bblim
ssm      = pft(co)%senm
sss      = pft(co)%sens
sslim    = pft(co)%senlim
lairat   = pft(co)%lrat

bb2bbmin = 315
bb2bbmax = 375
gs = 30

rlai = ssv(co)%lai%tot(1)

yield = 0.0

ss = 0

!----------------------------------------------------------------------!
! Check for chilling.                                                  !
!----------------------------------------------------------------------!
! If there is no chill,sum up the bbsum over last 20 days
! which here is not the GDD.
! If it's value is less than -100,then call it chill. 
if (chill==0) then
  bbsum = 0.0
  do i=1,20
    if (ssp%tmem(i)<-5.0) &
  bbsum = bbsum + max(-10.0,ssp%tmem(i)+5.0)
  enddo
  if (bbsum<-100) then
    chill = 1
    dschill = 1
  endif
endif
! If it is chill, then add to days of chill dschill
if (chill==1) then
  dschill = dschill + 1
endif
! If days of chill exceed 260 then reset chill
if (dschill>260) then
  chill = 0
  dschill = 0
endif

!if (bb > 0) stop

!----------------------------------------------------------------------!
! Bubburst, if no budburst set, and sufficient soil moisture, then     !
! check for  budburst.                                                 !
!----------------------------------------------------------------------!
! Three checks for budburst.First it checks the soil moisture for the
! current day (soil2g) and whether budburst has already occured (bb).Then
! it checks whether we had water input in the soil over the past month
! (sm_trig).Finally,it calculates GDD (bbsum) and checks if it has
! exceeded a limit. 

soilw = ssv(co)%soil_h2o(1) + ssv(co)%soil_h2o(2) + &
 ssv(co)%soil_h2o(3) + ssv(co)%soil_h2o(4)
soil2g = soilw/(ssp%soil_depth*10.0)

! If bb=0 which means that we don't have a day of budburst yet and
! soil mositure adequate,then look for budburst occurence
if (((bb==0).and.(soil2g>wtwp+0.5*(wtfc-wtwp))).or. &
 ((dsbb>bb2bbmax).and.(soil2g>wtwp+0.1*(wtfc-wtwp)))) then

  smtrig = 0.0
  do i=1,30
    smtrig = smtrig + ssv(co)%sm_trig(i)
  enddo

  if ((smtrig>30.0).or.(dsbb>bb2bbmax)) then

!----------------------------------------------------------------------!
! Check for budburst using degree days.                                !
!----------------------------------------------------------------------!
    bbsum = 0.0
    do i=1,bbm
      if (ssp%tmem(i)>bb0)  bbsum = bbsum + min(bbmax,ssp%tmem(i)-bb0)
    enddo
    
    ! If bbsum which here is GDD has exceeded the threshold then
    ! we have budburst.
    if ((real(bbsum)>=real(bblim)*exp(-0.01*real(dschill))) &
      .or.(dsbb>bb2bbmax)) then
!----------------------------------------------------------------------!
! Adjust proportion of gpp going into stem production based on suma.   !
! This is essentially the LAI control.                                 !
!----------------------------------------------------------------------!
      ! Budburst has occured
      ! ssv(co)%suma%tot is the total leaf carbon in all comp,probably
      ! tgp%p_opt=1.5 and tgp%P_laimem=0.5
      ! leafls is the leaf life span which we have set for crops to 180!REVISE
      ! It assigns part of the suma to go to stem 
      tsuma = ssv(co)%suma%tot
      maint = max(1.0,(real(leafls)/360.0)*1.0)
      tsuma = tsuma - msv%mv_leafmol*1.25/maint*tgp%p_opt
      if (tsuma>msv%mv_leafmol*1.25/maint*tgp%p_opt) tsuma = msv%mv_leafmol*1.25/maint*tgp%p_opt
      if (tsuma<-msv%mv_leafmol*1.25/maint*tgp%p_opt) tsuma = -msv%mv_leafmol*1.25/maint*tgp%p_opt
      stemfr = stemfr + tsuma*tgp%p_laimem*12.0

      if (stemfr<120.0) stemfr = 120.0

!----------------------------------------------------------------------!
! Budburst occurance.                                                  !
!----------------------------------------------------------------------!
!     On what day did budburst occur
      bb = (mnth-1)*30 + day
      bbgs = 0
      dsbb = 0

      if (stemfr<0.75*ssv(co)%nppstore(1)) then
        ssv(co)%nppstore(3) = ssv(co)%nppstore(1) - stemfr
      else
        if (stemfr<0.75*ssv(co)%nppstore(1)) then
          ssv(co)%nppstore(3) = ssv(co)%nppstore(1) - stemfr
        else
          ssv(co)%nppstore(3) = ssv(co)%nppstore(1)*0.25
        endif
        stemfr = stemfr*0.95
      endif
! This laiinc is overwritten later,can be removed?
      laiinc = (ssv(co)%nppstore(1) - ssv(co)%nppstore(3))/msv%mv_leafmol/1.25/12.0
      ssv(co)%nppstore(2) = ssv(co)%nppstore(1)
    endif
  endif
endif

!----------------------------------------------------------------------!
! Compute length of growing season, and set to zero when equal to      !
! growing season.                                                      !
!----------------------------------------------------------------------!
! Here it is out of the loop of checking for budburst.As long as bb
! retains its value it sums up the days of the growing season.
! If the days of the growing season exceed a threshold then it sets
! bb to zero.

if (bb>0)  bbgs = bbgs + 1
if (bbgs-gs>bb2bbmin) then
  bb = 0
  bbgs = 0
endif

if (dsbb < 500) dsbb = dsbb + 1

!----------------------------------------------------------------------!
! Set LAI increase.                                                    !
!----------------------------------------------------------------------!
! Calculates laiinc (why is it calculated above as well?) and adds it
! to rlai after converting it
! Every day adds laiinc until the growing season ends (bb=0) see above,
! or the nppsore is depleted below 60?

if ((bb>0).and.(bbgs<gs).and.(ssv(co)%nppstore(1)>1.0)) then
  laiinc = lairat*(ssv(co)%nppstore(2) - ssv(co)%nppstore(3))/msv%mv_leafmol/1.25/12.0
  if (rlai+laiinc>maxlai)  laiinc = maxlai - rlai
  if (rlai+laiinc>11.5)  laiinc = 11.5 - rlai
  if ((rlai>0).and.(ssv(co)%nppstore(1)<60.0)) laiinc = 0.0
else
  laiinc = 0.0
endif

! It checks whether senescence has occured by checking both soil moisture
! and temperature.If it has laiinc is set to -lai which I guess will
! set rlai to zero later.Also keeps the julian day of senescence as ss
!----------------------------------------------------------------------!
! Senescence, if rlai is greater than zero, compute senescence.        !
!----------------------------------------------------------------------!
if (rlai>1.0e-6) then
  if (soil2g<wtwp*0.0) then
!----------------------------------------------------------------------!
! Senescence event due to soil moisture.                               !
!----------------------------------------------------------------------!
    laiinc = -rlai
    ss = day + (mnth - 1)*30
  elseif (bbgs>100) then
!----------------------------------------------------------------------!
! Check for senescence, senescence occurs there are 'sss' days colder  !
! than 'sslim' out of the last 'ssm' days.                             !
!----------------------------------------------------------------------!
    sssum = 0
    do i=1,ssm
      if (ssp%tmem(i)<sslim)  sssum = sssum + 1
    enddo
    if (sssum>=sss) then
!----------------------------------------------------------------------!
! Senescence event due to temperature.                                 !
!----------------------------------------------------------------------!
      laiinc =-rlai
      ss = day + (mnth - 1)*30
    endif
  endif
endif
!----------------------------------------------------------------------!

ssv(co)%stemfr      = stemfr
ssv(co)%bb          = bb
ssv(co)%ss          = ss
ssv(co)%bbgs        = bbgs
ssv(co)%dsbb        = dsbb
ssv(co)%chill       = chill
ssv(co)%dschill     = dschill

end subroutine phenology2

!**********************************************************************!
!                   phenology3 :: phenological_methods                 !
!                   ----------------------------------                 !
!                                                                      !
! Used for crops                                                       !
!                                                                      !
! SUBROUTINE phenology3()                                              !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief
!! @details
!! @author Lyla,EPK
!! @date Jan 2017
!----------------------------------------------------------------------!
subroutine phenology3(yield,laiinc)

real(dp), parameter :: maxlai = 10.9

real(dp) :: wtwp,wtfc,ftagh,stemfr,bb0,bbmax,bblim,sslim,lairat
real(dp) :: rlai,soilw,soil2g,bbsum,maint,tsuma,laiinc,oldphen,newopt
real(dp) :: dft,dfp,hrs,phen,oldopt,optinc,sssum,yield,gddmix(3),hi,fv
integer :: co,day,mnth,ftdth,bb,ss,bbgs,chill,dschill,leafls,bbm
integer :: ssm,sss,i

!----------------------------------------------------------------------!

  co = ssp%cohort

  wtwp = ssp%wilt
  day  = ssp%day
  mnth = ssp%mnth

  wtfc     = ssp%field
  ftagh    = pft(co)%crop
  ftdth    = pft(co)%d2h
  stemfr   = ssv(co)%stemfr
  bb       = ssv(co)%bb
  ss       = ssv(co)%ss
  bbgs     = ssv(co)%bbgs
  leafls   = pft(co)%lls
  bbm      = pft(co)%bbmem
  bb0      = pft(co)%bb0
  bbmax    = pft(co)%bbmax
  bblim    = pft(co)%bblim
  ssm      = pft(co)%senm
  sss      = pft(co)%sens
  sslim    = pft(co)%senlim
  lairat   = pft(co)%lrat

  rlai = ssv(co)%lai%tot(1)
  
  yield = 0.0
  ss=0;

  soilw = ssv(co)%soil_h2o(1) + ssv(co)%soil_h2o(2) + &
   ssv(co)%soil_h2o(3) + ssv(co)%soil_h2o(4)
  soil2g = soilw/(ssp%soil_depth*10.0)

  ! If it's the start of the year set yield to zero.We do this because yield is
  ! cumulative added as in some years we might have two harvests
  IF(mnth.EQ.1.AND.day.EQ.1) ssv(co)%yield=0.


  IF(mnth.EQ.1.AND.day.EQ.1.AND.ssv(co)%sown.EQ.1) THEN
      pft(co)%sowday(3)=pft(co)%sowday(2)
      pft(co)%cropgdd(:,3)=pft(co)%cropgdd(:,2)
  ENDIF


  IF(ssv(co)%sown.EQ.0) THEN
      IF((day+ (mnth-1)*30).EQ.pft(co)%sowday(1)) THEN
          ssv(co)%sown=1
          ssv(co)%harvest(1)=0
          pft(co)%sowday(3)=pft(co)%sowday(1)
          pft(co)%cropgdd(:,3)=pft(co)%cropgdd(:,1)
      ELSEIF((day+ (mnth-1)*30).EQ.pft(co)%sowday(2)) THEN
          ssv(co)%sown=1
          ssv(co)%harvest(1)=0
          pft(co)%sowday(3)=pft(co)%sowday(2)
          pft(co)%cropgdd(:,3)=pft(co)%cropgdd(:,2)          
      ENDIF
  ENDIF

  
  ! Add up days that the crop has been sowed,set to zero if not
  IF(ssv(co)%sown.EQ.1) THEN
      ssv(co)%sowni=ssv(co)%sowni+1
  ELSE
      ssv(co)%sowni=0
  ENDIF


  ! If nothing is going to be sowed,return
  IF(ssv(co)%sown.EQ.0) THEN 
      ssv(co)%phu=0.
      laiinc=0.
      ssv(co)%vdays=0;
      RETURN
  ENDIF
  
    
  IF (ssv(co)%sown.EQ.1.AND.(bb==0).AND.(soil2g>wtwp+0.25*(wtfc-wtwp))) THEN
      !----------------------------------------------------------------------!
      ! Check for budburst using degree days.                                !
      !----------------------------------------------------------------------!

      bbsum = 0.0
      DO i=1,ssv(co)%sowni+1
          IF(ssp%tmem(i)>bb0)  bbsum = bbsum + MIN(bbmax,ssp%tmem(i)-bb0)
      ENDDO 
 
      IF(REAL(bbsum).GE.REAL(bblim)) THEN
      !----------------------------------------------------------------------!
      ! Adjust proportion of gpp going into stem production based on suma.   !
      ! This is essentially the LAI control.                                 !
      !----------------------------------------------------------------------!
      !tsuma = ssv(co)%suma%tot
      !maint = max(1.0,(real(leafls)/360.0)*1.0)
      !tsuma = tsuma - msv%mv_leafmol*1.25/maint*tgp%p_opt
      !IF(tsuma>msv%mv_leafmol*1.25/maint*tgp%p_opt) tsuma = msv%mv_leafmol*1.25/maint*tgp%p_opt
      !IF(tsuma<-msv%mv_leafmol*1.25/maint*tgp%p_opt) tsuma = -msv%mv_leafmol*1.25/maint*tgp%p_opt
      !stemfr = stemfr + tsuma*tgp%p_laimem*12.0
      !IF(stemfr<120.0) stemfr = 120.0

      !----------------------------------------------------------------------!
      ! Budburst occurance.                                                  !
      !----------------------------------------------------------------------!
          bb = (mnth-1)*30 + day
          bbgs = 0
    
      !IF(stemfr<0.8*ssv(co)%nppstore(1)) THEN
      !  ssv(co)%nppstore(3) = ssv(co)%nppstore(1) - stemfr
      !ELSE
      !  IF(stemfr<0.75*ssv(co)%nppstore(1)) THEN
      !    ssv(co)%nppstore(3) = ssv(co)%nppstore(1) - stemfr
      !  ELSE
      !    ssv(co)%nppstore(3) = ssv(co)%nppstore(1)*0.25
      !  ENDIF
      !  stemfr = stemfr*0.8
      !ENDIF !(stemfr<0.8*ssv(co)%nppstore(1))
      !laiinc = (ssv(co)%nppstore(1) - 0.0*ssv(co)%nppstore(3))/msv%mv_leafmol/1.25/12.0
      !ssv(co)%nppstore(2) = ssv(co)%nppstore(1)
      ENDIF ! (REAL(bbsum).GE.REAL(bblim))
  ENDIF

  IF(bb>0)  bbgs = bbgs + 1

  
  ! Previous value of phenology
  oldphen=ssv(co)%phu/pft(co)%cropgdd(1,3)


  hrs=dayl(ssp%lat,(mnth-1)*30+day)

  ! Calculates phen If we have budburst and phenology is less than 1
  ! Checks whether we are in vegetative phase (LT crophen(5)), and if we 
  ! require vernalization (croptype(2)).
  ! It also checks if we are in maturity stage (GT crophen(5)).In all 3 cases
  ! it calculates the phenological heat units by multiplying temperature by
  ! fv,ft and fp which are the vernalization,development and photoperiod functions
  IF(bb.GT.0.AND.oldphen.LT.1) THEN
      dft=0.0d0; dfp=0.0d0;fv=1.; 
      IF(oldphen.LT.pft(co)%cropphen(5)) THEN
          IF(pft(co)%croptype(2).EQ.1) THEN
              CALL streck(pft(co)%cardinal(4), &
                pft(co)%cardinal(5),pft(co)%cardinal(6), &
                ssp%tmem(1),pft(co)%croptype(1),pft(co)%photoperiod(3), &
                pft(co)%photoperiod(4),hrs,ssv(co)%vdays,fv)
          ENDIF
          CALL wangengel(pft(co)%cardinal(1),pft(co)%cardinal(2),pft(co)%cardinal(3) &
            ,ssp%tmem(1),pft(co)%croptype(1),pft(co)%photoperiod(1),pft(co)%photoperiod(2) &
            ,hrs,dft,dfp)
          pft(co)%cropgdd(2,3)=pft(co)%cardinal(1)
      ELSEIF(oldphen.ge.pft(co)%cropphen(5)) THEN
          CALL wangengel(pft(co)%cardinal(7),pft(co)%cardinal(8),pft(co)%cardinal(9) &
            ,ssp%tmem(1),pft(co)%croptype(1),pft(co)%photoperiod(5),pft(co)%photoperiod(6) &
            ,hrs,dft,dfp)
          pft(co)%cropgdd(2,3)=pft(co)%cardinal(7)
      ENDIF
      ssv(co)%phu=ssv(co)%phu+ssp%tmem(1)*fv*dft*dfp
  ENDIF ! (bb.gt.0.and.oldphen.lt.1)


  ! Phenological index of maturity,PHU divided by the max PHU obtained from seasonality
  ! for the specific crop and gridcell      
  phen=ssv(co)%phu/pft(co)%cropgdd(1,3)

  !IF(ssp%year.EQ.1991) THEN
  !  WRITE(*,*)ssp%year,mnth,day,phen,rlai
  !ENDIF

  ! If the days that it is sowed exceeds a limit then set phen straight to 1
  ! so it is harvested.For most crops the limit is set to 120 but for vern crops
  ! it is set to 210 to allow a much bigger growing cycle
  IF(ssv(co)%sowni.GT.pft(co)%limdharv) phen=1.
    
  ! If plants are mature (as determined by degree-days) harvest them
  ! Just check the flags,harvest will happen in allocation sub 
  IF(phen.GE.0.95.AND.ssv(co)%harvest(1).EQ.0) THEN
      ss=day + (mnth - 1)*30
      ssv(co)%harvest(1)=1
      ssv(co)%harvest(2)=day + (mnth - 1)*30
      ssv(co)%sown=0
      ssv(co)%bb=0
      ssv(co)%sowni=0
  ELSEIF(phen.GE.pft(co)%cropphen(5)) THEN
      ! Senescence begins, start killing leaves based on PHU
      ! Here we follow Eqn 3 or 4 of Bondeau et al 2007
      
      !laiinc=((((1.-phen)/(1.-pft(co)%cropphen(5)))**pft(co)%cropphen(6))* &
      !  (1-0.9)+0.9)*pft(co)%optlai

      ! Placeholder.Removes 1% of the current LAI for senescence every day after
      ! phenology index exceeds pft(co)%cropphen(5)
      laiinc=-0.01*rlai 
      
  ELSE 
      !Leaves are still allowed to grow
      !----------------------------------------------------------------------!
      ! Set LAI increase.                                                    !
      !----------------------------------------------------------------------!
      IF (ssv(co)%nppstore(1).GT.0.0d0) THEN
          ! Optimal LAI increase (Neisch et al 2002 SWAT documentation Eqns 18.1.9)
          ! cropphen(3) and (4) are the shape parameters after Neisch et al 2002
          oldopt=oldphen/(oldphen+exp(pft(co)%cropphen(3)-pft(co)%cropphen(4)*oldphen))
          newopt=phen/(phen+exp(pft(co)%cropphen(3)-pft(co)%cropphen(4)*phen))
          optinc=(newopt-oldopt)*pft(co)%optlai
          IF(optinc*msv%mv_leafmol*12.0*1.25.LT.0.5*ssv(co)%nppstore(1)) THEN
          !IF(optinc*12.0/pft(co)%sla/25.0.LT.0.5*ssv(co)%nppstore(1)) THEN
              laiinc=optinc
          ELSE
              gddmix(1)=oldphen*pft(co)%cropgdd(1,3)
              gddmix(2)=phen*pft(co)%cropgdd(1,3)
              gddmix(3)=gddmix(1) 

              !IF(ssp%year.EQ.2001) THEN
              !  WRITE(*,*)'A: ',mnth,day,oldphen,gddmix(3),phen,gddmix(2),&
              !  optinc*msv%mv_leafmol*12.0*1.25,0.5*ssv(co)%nppstore(1)
              !ENDIF

              laiinc=0. 
              DO WHILE (gddmix(3).LE.gddmix(2).AND.0.5*ssv(co)%nppstore(1).GE.laiinc*msv%mv_leafmol*1.25*12.0)
              !DO WHILE (gddmix(3).LE.gddmix(2).AND.0.5*ssv(co)%nppstore(1).GE.laiinc*12.0/pft(co)%sla/25.0) 

                  !IF(ssp%year.EQ.2001) THEN
                  !  WRITE(*,*)'B: ',mnth,day,oldphen,gddmix(3),phen,gddmix(2),&
                  !    laiinc*msv%mv_leafmol*12.0*1.25,0.5*ssv(co)%nppstore(1)
                  !ENDIF

                  oldopt=oldphen/(oldphen+exp(pft(co)%cropphen(3)-pft(co)%cropphen(4)*oldphen))
                  phen=gddmix(3)/pft(co)%cropgdd(1,3)
                  newopt=phen/(phen+exp(pft(co)%cropphen(3)-pft(co)%cropphen(4)*phen))
                  optinc=(newopt-oldopt)*pft(co)%optlai 
                  laiinc=optinc
                  ssv(co)%phu=gddmix(3) 
                  gddmix(3)=gddmix(3)+1
              ENDDO
          ENDIF
          IF(rlai+laiinc>pft(co)%optlai) laiinc=pft(co)%optlai-rlai 
          IF((rlai>0).and.(ssv(co)%nppstore(1)<0.0)) laiinc = 0.0
      ELSE
          laiinc = 0.0
      ENDIF
  ENDIF ! crop is mature or optimal LAI attained

  !----------------------------------------------------------------------!
  ! Senescence, if rlai is greater than zero, compute senescence.        !
  !----------------------------------------------------------------------!
  IF (rlai>1.0e-6) THEN
    IF(soil2g<wtwp*0.0) THEN
      !----------------------------------------------------------------------!
      ! Senescence event due to soil moisture.                               !
      !----------------------------------------------------------------------!
      laiinc = -rlai
      ss = day + (mnth - 1)*30
    ELSE
      !----------------------------------------------------------------------!
      ! Check for senescence, senescence occurs there are 'sss' days colder  !
      ! than 'sslim' out of the last 'ssm' days.                             !
      !----------------------------------------------------------------------!
      sssum = 0
      DO i=1,ssm
!        IF (ssp%tmem(i)<sslim)  sssum = sssum + 1
        IF (ssp%tmem(i).LT.pft(co)%lethal(1).OR.ssp%tmem(i).GT.pft(co)%lethal(2)) &
           sssum = sssum + 1
      ENDDO
      IF (sssum>=sss) THEN
        !----------------------------------------------------------------------!
        ! Senescence event due to temperature.                                 !
        !----------------------------------------------------------------------!
        laiinc =-rlai
        ss = day + (mnth - 1)*30
      ENDIF
    ENDIF
  ENDIF


  ssv(co)%stemfr      = stemfr
  ssv(co)%bb          = bb
  ssv(co)%ss          = ss
  ssv(co)%bbgs        = bbgs
  ssv(co)%chill       = chill
  ssv(co)%dschill     = dschill


end subroutine phenology3


!**********************************************************************!
!                                                                      !
!                    allocation :: phenological_methods                !
!                    ----------------------------------                !
! SUBROUTINE allocation(laiinc,daygpp,resp_l,lmor_sc,resp,leaflit,&    !
! stemnpp,rootnpp,resp_s,resp_r,resp_m,check_closure)                  !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Allocate GPP between storage, stems, roots and leaves.
!! @details
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine allocation(laiinc,daygpp,resp_l,lmor_sc,resp,leaflit,&
 stemnpp,rootnpp,resp_s,resp_r,resp_m,check_closure)
!**********************************************************************!
real(dp) :: laiinc,leaflit,resp,ans,yy,stemnpp,total_carbon,old_total_carbon, &
 root_fixed,stem_fixed,lmor_sc(3600),daynpp,rootnpp,resp_s,resp_r,resp_l, &
 daygpp,resp_m,lit
integer :: i,co,k
logical :: check_closure
real(dp) :: sumrr,sumsr,sumlr,summr,yielit
save :: sumrr,sumsr,sumlr,summr
!----------------------------------------------------------------------!

co = ssp%cohort

daynpp = daygpp - resp_l

if (pft(co)%phen/=0) then

ssv(co)%nppstore(1) = ssv(co)%nppstore(1) + daygpp - resp_l

resp = resp_l

!----------------------------------------------------------------------!
! Check carbon closure.
!----------------------------------------------------------------------!
if (check_closure) then
  old_total_carbon = ssv(co)%lai%tot(1)*12.0/pft(co)%sla/25.0 + ssv(co)%nppstore(1) + &
 ssv(co)%stem%tot(1) + ssv(co)%root%tot(1) + resp + ssv(co)%bio(1) + ssv(co)%bio(2)
endif

! Subtracts carbon for new leaves of laiinc from nppstore(1) after converting
! from LAI to carbon.The reason it has the 1.25 is for the additional carbon
! cost of creating leaves.That cost is added as respiration by multiplying
! by 0.25.
!----------------------------------------------------------------------!
! Pay for new leaves.
!----------------------------------------------------------------------!
if (laiinc>0.0) then
  ssv(co)%nppstore(1) = ssv(co)%nppstore(1) - laiinc*msv%mv_leafmol*1.25*12.0
  ssv(co)%nppstore(2) = ssv(co)%nppstore(2) - laiinc*msv%mv_leafmol*1.25*12.0
  resp_m = 0.25*laiinc*msv%mv_leafmol*12.0
  resp = resp + resp_m
else
  resp_m = 0.0
endif

!----------------------------------------------------------------------!
! Age leaves by one day, kill any which have died of old age,
! adjust by laiinc (+ or -), and then sum to get 'rlai'.
!----------------------------------------------------------------------!
call LAI_ADD(laiinc,leaflit)

!----------------------------------------------------------------------!
! leaf death not through age mortality 
!----------------------------------------------------------------------!
call LAI_DIST(lmor_sc,leaflit)

! Deals with crop harvest.Needs the harvest flag and a positive lai.
! lit is the C litter that is produced after harvest and yielit is the C yield.
! The latter is completely removed from the system which is the reason why it 
! appears in the carbon closure checks,so that it wont signal a carbon breach,
! even though there is one.The sum of lit and yielit equals the crop C right 
! before it was harvested as it can be seen in the harv sub
lit=0.
yielit=0.
IF(pft(co)%phen.EQ.3.AND.ssv(co)%harvest(1).EQ.1.AND.ssv(co)%lai%tot(1).GT.0.) THEN
  CALL HARV(lit,yielit)
ENDIF
leaflit=leaflit+lit

!----------------------------------------------------------------------!
! Check carbon closure.
!----------------------------------------------------------------------!
if (check_closure) then
  ans = ssv(co)%lai%tot(1)*12.0/pft(co)%sla/25.0 + ssv(co)%nppstore(1) + &
 ssv(co)%stem%tot(1) + ssv(co)%root%tot(1) + resp + ssv(co)%bio(1) + ssv(co)%bio(2)
  if (abs(ans+leaflit+yielit-old_total_carbon) > 1.0e-3) then
    write(*,*) 'Breach of carbon after leaves:',ans+leaflit-old_total_carbon,' g/m^2.'
    write(*,*) ans,leaflit,old_total_carbon
    write(*,*) 'lai ',ssv(co)%lai
    write(*,*) 'sla ',pft(co)%sla
    write(*,*) 'nppstore ',ssv(co)%nppstore(1)
    write(*,*) '%stem ',ssv(co)%stem
    write(*,*) '%root ',ssv(co)%root
    write(*,*) 'resp ',resp
    write(*,*) 'biomass ',ssv(co)%bio(1)
    write(*,*) 'biomass ',ssv(co)%bio(2)
    stop
  endif
endif

!----------------------------------------------------------------------!
! Pay for days roots if veg exists
!----------------------------------------------------------------------!
! A fraction of nppstore goes to the roots (yy) or rootnpp.
! The fraction is constant
if ((ssv(co)%nppstore(1)>0.0).and.(daynpp>0.0)) then
  yy = ssv(co)%nppstore(1)*tgp%p_rootfr
else
  yy = 0.0
endif

ssv(co)%nppstore(1) = ssv(co)%nppstore(1) - yy
rootnpp = yy

!----------------------------------------------------------------------!
! Age roots and add todays root npp.
!----------------------------------------------------------------------!
call ROOT_ADD(yy,root_fixed)
ssv(co)%bio(2) = ssv(co)%bio(2) + root_fixed

!----------------------------------------------------------------------!
! Root resperation
!----------------------------------------------------------------------!
! Calculates root respiration.msv%mv_respref is a function of soil moisture
! and temperature and its calculated daily in SET_MISC_VALUES
call ROOT_DIST(msv%mv_respref,resp_r)
resp = resp + resp_r
rootnpp = rootnpp - resp_r

!----------------------------------------------------------------------!
! Check carbon closure.
!----------------------------------------------------------------------!
if (check_closure) then
  ans = ssv(co)%lai%tot(1)*12.0/pft(co)%sla/25.0 + ssv(co)%nppstore(1) + &
 ssv(co)%stem%tot(1) + ssv(co)%root%tot(1) + resp + ssv(co)%bio(1) + ssv(co)%bio(2)
  if (abs(ans+leaflit+yielit-old_total_carbon) > 1.0e-3) then
    write(*,*) 'Breach of carbon after roots:',ans+leaflit-old_total_carbon,' g/m^2.'
    stop
  endif
endif

!----------------------------------------------------------------------!
! Stem resperation, and NPP                                            !
!----------------------------------------------------------------------!
! A fraction of nppstore goes to the stem (yy) or stemnpp.
! The fraction is constant
if ((ssv(co)%nppstore(1)>0.0).and.(daynpp>0.0)) then
  yy = ssv(co)%nppstore(1)*tgp%p_stemfr
  ssv(co)%nppstore(1) = ssv(co)%nppstore(1) - yy
else
  yy = 0.0
endif
stemnpp = yy

!----------------------------------------------------------------------!
! Age stems and add todays stem npp.
!----------------------------------------------------------------------!
call STEM_ADD(yy,stem_fixed)
ssv(co)%bio(1) = ssv(co)%bio(1) + stem_fixed

!----------------------------------------------------------------------!
! Stem respiration.
!----------------------------------------------------------------------!
! Calculates stem respiration.msv%mv_respref is a function of soil moisture
! and temperature and its calculated daily in SET_MISC_VALUES
call STEM_DIST(msv%mv_respref,resp_s)
resp = resp + resp_s
stemnpp = stemnpp - resp_s

!----------------------------------------------------------------------!
! Check carbon closure.
!----------------------------------------------------------------------!
if (check_closure) then
  ans = ssv(co)%lai%tot(1)*12.0/pft(co)%sla/25.0 + ssv(co)%nppstore(1) + &
 ssv(co)%stem%tot(1) + ssv(co)%root%tot(1) + resp + ssv(co)%bio(1) + ssv(co)%bio(2)
  if (abs(ans+leaflit+yielit-old_total_carbon) > 1.0e-3) then
    write(*,*) 'Breach of carbon after stems:',ans+leaflit-old_total_carbon,' g/m^2.'
    stop
  endif
endif

!----------------------------------------------------------------------!

daynpp = daynpp - resp
daynpp = daygpp - resp
ssv(co)%npp = daynpp

!----------------------------------------------------------------------!
! Check carbon closure.
!----------------------------------------------------------------------!
if (check_closure) then
  total_carbon = ssv(co)%lai%tot(1)*12.0/pft(co)%sla/25.0 + ssv(co)%nppstore(1) + &
 ssv(co)%stem%tot(1) + ssv(co)%root%tot(1) + resp + ssv(co)%bio(1) + ssv(co)%bio(2)
  if (abs(total_carbon+leaflit+yielit-old_total_carbon) > 1.0e-3) then
    write(*,*) 'Breach of carbon in phenology:',ans+leaflit-old_total_carbon,' g/m^2.'
    stop
  endif
endif

endif

if (ssp%cohort == 1) then
      if ((ssp%day == 1).and.(ssp%mnth == 1)) then
        sumsr=0.0
        sumrr=0.0
        sumlr=0.0
        summr=0.0
      endif
      sumsr = sumsr + resp_s
      sumrr = sumrr + resp_r
      sumlr = sumlr + resp_l
      summr = summr + resp_m
endif

end subroutine allocation





!**********************************************************************!
!                                                                      !
!                      suma_add :: phenological_methods                !
!                      --------------------------------                !
!                                                                      !
!  SUBROUTINE suma_add(laiinc)                                         !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Update suma by adding daily increase and removing via mortality.
!! @details suma mortality is hardwired at 360 days.
!! suma_comp_length set to 30 in dims
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine suma_add(laiinc)
!**********************************************************************!
real(dp) :: laiinc
integer :: i,co,k
!----------------------------------------------------------------------!

co = ssp%cohort

!----------------------------------------------------------------------!
! Leaf death through old age.
!----------------------------------------------------------------------!
if (ssv(co)%suma%no > 0) then
  if (ssp%jday-int(ssv(co)%suma%c(1)%age) > 360) then
    k = suma_comp_length-(ssp%jday-int(ssv(co)%suma%c(1)%age)-360)+1
    if (k == 1) then
      do i=1,ssv(co)%suma%no-1
        ssv(co)%suma%c(i)%val = ssv(co)%suma%c(i+1)%val
        ssv(co)%suma%c(i)%age = ssv(co)%suma%c(i+1)%age
      enddo
      ssv(co)%suma%c(ssv(co)%suma%no)%val = 0.0
      ssv(co)%suma%c(ssv(co)%suma%no)%age = 0.0
      ssv(co)%suma%no = ssv(co)%suma%no - 1
    else
      ssv(co)%suma%c(1)%val = ssv(co)%suma%c(1)%val*(real(k-1)/real(k))
    endif
  endif
endif

!----------------------------------------------------------------------!
! Add on lai increase, create new compartment if necessary.
!----------------------------------------------------------------------!
if (laiinc>0.0) then
  if (ssv(co)%suma%no == 0) then
    ssv(co)%suma%no = 1
    ssv(co)%suma%c(ssv(co)%suma%no)%val = 0.0
    ssv(co)%suma%c(ssv(co)%suma%no)%age = ssp%jday
  elseif (ssp%jday-int(ssv(co)%suma%c(ssv(co)%suma%no)%age) >= suma_comp_length) then
    ssv(co)%suma%no = ssv(co)%suma%no + 1
    ssv(co)%suma%c(ssv(co)%suma%no)%val = 0.0
    ssv(co)%suma%c(ssv(co)%suma%no)%age = ssp%jday
  endif
  ssv(co)%suma%c(ssv(co)%suma%no)%val =  ssv(co)%suma%c(ssv(co)%suma%no)%val + laiinc
endif

ssv(co)%suma%tot = 0.0
do i=1,ssv(co)%suma%no
  ssv(co)%suma%tot = ssv(co)%suma%tot + ssv(co)%suma%c(i)%val
enddo

end subroutine suma_add
   




!**********************************************************************!
!                                                                      !
!                   stem_add  :: phenological_methods                  !
!                   ---------------------------------                  !
!                                                                      !
!  SUBROUTINE stem_add(stem_inc,stem_fixed)                            !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Update stems by adding daily increase and removing via mortality.
!! @details Stem mortality is given in the pft parameters (pft%sls).
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine stem_add(stem_inc,stem_fixed)
!**********************************************************************!
real(dp) :: stem_inc,stem_fixed
integer :: i,co,k
!----------------------------------------------------------------------!

! Check lai_add comments,this is an equivalent subroutine to lai_add but 
! for stem.Dont mind the variable names.Unlike lai_add where it deals
! with lai,this is carbon.
! The dead carbon from stem doesnt go to the litter as it does for lai_add
! but is added to bio(1) pool when it is done calling this sub.


co = ssp%cohort
stem_fixed = 0.0

!----------------------------------------------------------------------!
! Stem death through old age.
!----------------------------------------------------------------------!
if (ssv(co)%stem%no > 0) then
  if (ssp%jday-int(ssv(co)%stem%c(1,1)%age) > pft(co)%sls) then
    k = stem_comp_length-(ssp%jday-int(ssv(co)%stem%c(1,1)%age)-pft(co)%sls)+1
    if (k == 1) then
      stem_fixed = ssv(co)%stem%c(1,1)%val
      do i=1,ssv(co)%stem%no-1
        ssv(co)%stem%c(i,1)%val = ssv(co)%stem%c(i+1,1)%val
        ssv(co)%stem%c(i,1)%age = ssv(co)%stem%c(i+1,1)%age
      enddo
      ssv(co)%stem%c(ssv(co)%stem%no,1)%val = 0.0
      ssv(co)%stem%c(ssv(co)%stem%no,1)%age = 0.0
      ssv(co)%stem%no = ssv(co)%stem%no - 1
    else
      stem_fixed = ssv(co)%stem%c(1,1)%val/real(k)
      ssv(co)%stem%c(1,1)%val = ssv(co)%stem%c(1,1)%val*(real(k-1)/real(k))
    endif
  endif
endif

! Check lai_add for comments
!----------------------------------------------------------------------!
! Add on stem increase, create new compartment if necessary.
!----------------------------------------------------------------------!
if (stem_inc>0.0) then
  if (ssv(co)%stem%no == 0) then
    ssv(co)%stem%no = 1
    ssv(co)%stem%c(ssv(co)%stem%no,1)%val = 0.0
    ssv(co)%stem%c(ssv(co)%stem%no,1)%age = ssp%jday
  elseif (ssp%jday-int(ssv(co)%stem%c(ssv(co)%stem%no,1)%age) >= stem_comp_length) then
    ssv(co)%stem%no = ssv(co)%stem%no + 1
    ssv(co)%stem%c(ssv(co)%stem%no,1)%val = 0.0
    ssv(co)%stem%c(ssv(co)%stem%no,1)%age = ssp%jday
  endif
  ssv(co)%stem%c(ssv(co)%stem%no,1)%val =  ssv(co)%stem%c(ssv(co)%stem%no,1)%val + stem_inc
endif

end subroutine stem_add





!**********************************************************************!
!                                                                      !
!                    root_add :: phenology_methods                     !
!                    -----------------------------                     !
!                                                                      !
!  SUBROUTINE root_add(laiinc,leaflit)                                 !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Update roots by adding daily increase and removing via mortality.
!! @details Root mortality is given in the pft parameters (pft%rls).
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine root_add(laiinc,leaflit)
!**********************************************************************!
real(dp) :: laiinc,leaflit
integer i,co,k
!----------------------------------------------------------------------!

! Check lai_add comments,this is an equivalent subroutine to lai_add but 
! for roots.Dont mind the variable names.Unlike lai_add where it deals
! with lai,this is carbon.
! The dead carbon from roots doesnt go to the litter as it does for lai_add
! but is added to bio(2) pool when it is done calling this sub.

co = ssp%cohort
leaflit = 0.0

!----------------------------------------------------------------------!
! Root death through old age.
!----------------------------------------------------------------------!
if (ssv(co)%root%no > 0) then
  if (ssp%jday-int(ssv(co)%root%c(1,1)%age) > pft(co)%rls) then
    k = root_comp_length-(ssp%jday-int(ssv(co)%root%c(1,1)%age)-pft(co)%rls)+1
    if (k == 1) then
      leaflit = ssv(co)%root%c(1,1)%val
      do i=1,ssv(co)%root%no-1
        ssv(co)%root%c(i,1)%val = ssv(co)%root%c(i+1,1)%val
        ssv(co)%root%c(i,1)%age = ssv(co)%root%c(i+1,1)%age
      enddo
      ssv(co)%root%c(ssv(co)%root%no,1)%val = 0.0
      ssv(co)%root%c(ssv(co)%root%no,1)%age = 0.0
      ssv(co)%root%no = ssv(co)%root%no - 1
    else
      leaflit = ssv(co)%root%c(1,1)%val/real(k)
      ssv(co)%root%c(1,1)%val = ssv(co)%root%c(1,1)%val*&
 (real(k-1)/real(k))
    endif
  endif
endif

! Check lai_add for comments
!----------------------------------------------------------------------!
! Add on root increase, create new compartment if necessary.
!----------------------------------------------------------------------!
if (laiinc>0.0) then
  if (ssv(co)%root%no == 0) then
    ssv(co)%root%no = 1
    ssv(co)%root%c(ssv(co)%root%no,1)%val = 0.0
    ssv(co)%root%c(ssv(co)%root%no,1)%age = ssp%jday
  elseif (ssp%jday-int(ssv(co)%root%c(ssv(co)%root%no,1)%age) &
 >= root_comp_length) then
    ssv(co)%root%no = ssv(co)%root%no + 1
    ssv(co)%root%c(ssv(co)%root%no,1)%val = 0.0
    ssv(co)%root%c(ssv(co)%root%no,1)%age = ssp%jday
  endif
  ssv(co)%root%c(ssv(co)%root%no,1)%val =  &
 ssv(co)%root%c(ssv(co)%root%no,1)%val + laiinc
endif

end subroutine root_add





!**********************************************************************!
!                                                                      !
!                    lai_add :: phenological_methods                   !
!                    -------------------------------                   !
!                                                                      !
! SUBROUTINE lai_add(laiinc,leaflit)                                   !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief
!! @details
!! @author Mark Lomas,EPK
!! @date May,2017
!----------------------------------------------------------------------!
subroutine lai_add(laiinc,leaflit)
!**********************************************************************!
real(dp) :: laiinc,leaflit,lai_tot
integer :: i,co,k
!----------------------------------------------------------------------!

co = ssp%cohort

! For each cohort there will be carbon compartments where the leaf LAI
! is stored,several for the whole canopy.
! These compartments can be found in ssv(co)%lai where
! ssv(co)%lai%no is the number of available compartments,
! ssv(co)%lai%c(ssv(co)%lai%no)%age the julian day when the compartment
! was created and ssv(co)%lai%c(ssv(co)%lai%no)%val the lai
! in that compartment (not the carbon!)

! For old age it checks if the age of the first compartment which will always
! be the older one is greater than the leaf life span lls,usually 180 days.
! If that applies then the compartment will begin to senense which doesnt
! happen immidiately but takes  lai_comp_length days,set to 10.
! k is a countdown from lai_comp_length to 1 days, immidiately after the
! leaf reached its life span. For each day afterwards a fraction of its
! lai will be sent to the litter through the leaflit variable
! and the remainder will stay in the leaf.
! The process will continue until k==1 where the compartment will die
! ,the remainding lai goes to litter and all other compartments
! indexes move down by 1 index.Index 1 is always the oldest that dies first.

! Even though the leaflit variable holds lai,in lai_dist sub it is finally
! converted to carbon before being send to the litter.

!----------------------------------------------------------------------!
! Leaf death through old age.
!----------------------------------------------------------------------!
if (ssv(co)%lai%no > 0) then
  if (ssp%jday-int(ssv(co)%lai%c(1,1)%age) > pft(co)%lls) then
    k = lai_comp_length-(ssp%jday-int(ssv(co)%lai%c(1,1)%age)-pft(co)%lls)+1
    if (k == 1) then
      leaflit = ssv(co)%lai%c(1,1)%val
      do i=1,ssv(co)%lai%no-1
        ssv(co)%lai%c(i,1)%val = ssv(co)%lai%c(i+1,1)%val
        ssv(co)%lai%c(i,1)%age = ssv(co)%lai%c(i+1,1)%age
      enddo
      ssv(co)%lai%c(ssv(co)%lai%no,1)%val = 0.0
      ssv(co)%lai%c(ssv(co)%lai%no,1)%age = 0.0
      ssv(co)%lai%no = ssv(co)%lai%no - 1
    else
      leaflit = ssv(co)%lai%c(1,1)%val/real(k)
      ssv(co)%lai%c(1,1)%val = ssv(co)%lai%c(1,1)%val*(real(k-1)/real(k))
    endif
  endif
endif

! Checks if laiinc>0.If there is no LAI compartment it will add one. 
! If the oldest compartment has reached a certain age it will create
! another one and add LAI to the younger compartment.
! Carbon keeps getting added to each compartment.
!----------------------------------------------------------------------!
! Add on lai increase, create new compartment if necessary.
!----------------------------------------------------------------------!
if (laiinc>0.0) then
  if (ssv(co)%lai%no == 0) then
    ssv(co)%lai%no = 1
    ssv(co)%lai%c(ssv(co)%lai%no,1)%val = 0.0
    ssv(co)%lai%c(ssv(co)%lai%no,1)%age = ssp%jday
  elseif (ssp%jday-int(ssv(co)%lai%c(ssv(co)%lai%no,1)%age) >= lai_comp_length) then
    ssv(co)%lai%no = ssv(co)%lai%no + 1
    ssv(co)%lai%c(ssv(co)%lai%no,1)%val = 0.0
    ssv(co)%lai%c(ssv(co)%lai%no,1)%age = ssp%jday
  endif
  ssv(co)%lai%c(ssv(co)%lai%no,1)%val =  ssv(co)%lai%c(ssv(co)%lai%no,1)%val + laiinc
else
  if(laiinc<0.0) then
    lai_tot=ssv(co)%lai%tot(1)
    do i=1,ssv(co)%lai%no
      leaflit=leaflit-ssv(co)%lai%c(i,1)%val*laiinc/lai_tot
      ssv(co)%lai%c(i,1)%val=ssv(co)%lai%c(i,1)%val*(lai_tot+laiinc)/lai_tot
    enddo
  endif
endif


end subroutine lai_add





!**********************************************************************!
!                                                                      !
!                      stem_dist :: phenology_methods                  !
!                      ------------------------------                  !
!                                                                      !
!  SUBROUTINE stem_dist(respref,resp)                                  !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Apply stem respiration to the active stems, and sum the
!! resulting stems.
!! @details Apply stem respiration to each of the stem compartments
!! which span the age, for the current cohort.
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine stem_dist(respref,resp)
!**********************************************************************!
integer :: i,co
real(dp) :: ans,resp,respref
!**********************************************************************!

! Applies stem respiration regulated by respref.Does so for all compartments
! removing the carbon from each and adding it up to the total resp.
! Also calculates the total carbon in the stem compartments

co = ssp%cohort

resp = 0.0

do i=1,ssv(co)%stem%no
  ans = ssv(co)%stem%c(i,1)%val*respref
  resp = resp + ans
  ssv(co)%stem%c(i,1)%val = ssv(co)%stem%c(i,1)%val - ans
enddo

ssv(co)%stem%tot(1) = 0.0
do i=1,ssv(co)%stem%no
  ssv(co)%stem%tot(1) = ssv(co)%stem%tot(1) + ssv(co)%stem%c(i,1)%val
enddo

end subroutine stem_dist





!**********************************************************************!
!                                                                      !
!                   root_dist :: phenology_methods                     !
!                   ------------------------------                     !
!                                                                      !
!  SUBROUTINE ROOT_DIST(respref,resp)                                  !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Apply root respiration to the active roots, and sum the
!! resulting roots.
!! @details Apply root respiration to each of the root compartments
!! which span the age, for the current cohort.
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine ROOT_DIST(respref,resp)
!**********************************************************************!
integer :: i,co
real(dp) :: ans,resp,respref
!----------------------------------------------------------------------!

! Applies root respiration regulated by respref.Does so for all compartments
! removing the carbon from each and adding it up to the total resp.
! Also calculates the total carbon in the root compartments

co = ssp%cohort

resp = 0.0

do i=1,ssv(co)%root%no
  ans = ssv(co)%root%c(i,1)%val*respref
  resp = resp + ans
  ssv(co)%root%c(i,1)%val = ssv(co)%root%c(i,1)%val - ans
enddo

ssv(co)%root%tot(1) = 0.0
do i=1,ssv(co)%root%no
  ssv(co)%root%tot(1) = ssv(co)%root%tot(1) + ssv(co)%root%c(i,1)%val
enddo

end subroutine ROOT_DIST





!**********************************************************************!
!                                                                      !
!                     lai_dist :: phenology_methods                    !
!                     -----------------------------                    !
!                                                                      !
! SUBROUTINE LAI_DIST(lmor_sc,leaflit)                                 !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief (Currently not implemented). Apply leaf motality.
!! @details Leaf mortality is applied to the current cohort.
!! The mortality values are taken from 'lmor_sc' and are applied
!! to the lai compartments which span the age of the leaf.
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine lai_dist(lmor_sc,leaflit)
!**********************************************************************!
integer :: i,co,nlcs
real(dp) :: lmor_sc(3600),ans,leaflit
!----------------------------------------------------------------------!

co = ssp%cohort

! Applies daily leaf mortality which is different from old age.
! It is a function of the age of the compartment and its value can be found in
! lmor_sc.Here is not really implemented as the values it gives (ans) are close to
! 1 which will make the mortality (1-ans) 0.

do i=1,ssv(ssp%cohort)%lai%no
  ans = lmor_sc(int(ssp%jday-ssv(co)%lai%c(i,1)%age+1.5))*0.0+1.0
  leaflit = leaflit + ssv(co)%lai%c(i,1)%val*(1.0 - ans)
  ssv(co)%lai%c(i,1)%val = ssv(co)%lai%c(i,1)%val*ans
enddo

! Sum up the lai compartments.
ssv(co)%lai%tot(1) = 0.0
do i=1,ssv(co)%lai%no
  ssv(co)%lai%tot(1) = ssv(co)%lai%tot(1) + ssv(co)%lai%c(i,1)%val
enddo

! Leaf litter
leaflit = leaflit*12.0/pft(co)%sla/25.0

end subroutine lai_dist



end module phenology_methods
