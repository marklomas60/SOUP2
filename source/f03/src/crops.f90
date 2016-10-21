!> @brief Crop subroutines
!! @details
!! @author 
!! @date 

module crops

use real_precision
use pft_parameters
use site_parameters

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
!!     Waha et al 2012, Global Ecology and Biogeography 21:247-259
!!     van Bussel et al 2015, Global Ecology and Biogeography (earlyview)
!!
!! avtmp(*)     input (20-year exp-mean monthly temperatures)
!! avprc(*)     input (20-year exp-mean monthly precipitation)
!! mnthhum(*)   input (mean monthly humidity)
!! lat          input Latitude
!! nft          input Number of plant functional types (PFTs)
!! sumtsow(*)   input Minimum temperature thresholds for sowing
!! wintsow(*)   input Maximum temperature thresholds for sowing
!! ndays(12)    input Number of days per month
!! tmp(12,31)   input Mean daily temperatures (C)
!! prc(12,31)   input Mean daily precipitation (mm)
!! sowday(*)   output Sowing day for each PFT
!! pscale(*)    input Flags to reduce optimal photoperiods for each PFT
!! popt(*)     output Optimal photoperiods for each PFT
!! vernmsat(*)  input Max vernalisation requirement per PFT (days/month)
!! vernt2(*)    input Min temperature for optimal vernalisation per PFT (C)
!! vernt3(*)    input Max temperature for optimal vernalisation per PFT (C)
!! vernsat(*)  output Max vernalisation requirement per PFT (days)
!! verbose      input Logical: write out values or not
!!
!! @author Lyla 
!! @date 
!----------------------------------------------------------------------!
subroutine seasonality(tmp,prc,cld,nyears)
!**********************************************************************!
implicit none
real(dp) :: tmp(12,31),prc(12,31),cld(12)
integer  :: nyears

real :: sowday

real(dp) :: sowthresh(2),croprange(4),cardinal(9),phot(6),cropgdd(2)
real(dp) :: lethal(2),cropphen(6)
integer :: croptype(2)
!----------------------------------------------------------------------!

sowthresh(1)=100.
sowthresh(2)=14.
croprange(1)=1000.
croprange(2)=2000.
croprange(3)=5.
croprange(4)=5.
cardinal(1)=9.3
cardinal(2)=18.3
cardinal(3)=39.2
cardinal(4)=0
cardinal(5)=0
cardinal(6)=0
cardinal(7)=8.0
cardinal(8)=26.4
cardinal(9)=36.0
croptype(1)=-1.
croptype(2)=0
phot(1)=12.5
phot(2)=24.0
phot(3)=12.5
phot(4)=24.
phot(5)=12.5
phot(6)=24.0
lethal(1)=-1.8
lethal(2)=46
cropphen(1)=3;
cropphen(2)=7;
cropphen(3)=3.055;
cropphen(4)=13.385;
cropphen(5)=0.7;
cropphen(6)=0.5;







end subroutine seasonality

end module crops

