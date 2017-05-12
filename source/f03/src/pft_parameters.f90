module pft_parameters

use real_precision, only: dp
use dims

private :: dp, dls

public :: pft, pft_tab

type PftParameters
  character              :: tag*1000
  integer                :: itag
  real(dp)               :: mix
  integer                :: c3c4
  integer                :: phen
  real(dp)               :: crop
  integer                :: d2h
  integer                :: mort
  real(dp)               :: wden
  real(dp)               :: xyl
  real(dp)               :: pdif
  real(dp)               :: sla ! specific leaf area ()
  integer                :: lls ! Leaf life span (days)
  integer                :: sls ! Stem life span (days)
  integer                :: rls ! Root life span (days)
  real(dp)               :: lmor ! Leaf mortality
  real(dp)               :: lrat ! Leaf rate
  integer                :: bbmem ! Budburst memory
  real(dp)               :: bb0  ! Budburst lower limit
  real(dp)               :: bbmax ! Budburst upper limit
  real(dp)               :: bblim ! Budburst limit
  integer                :: senm ! Senescence memory
  integer                :: sens ! Senescence sensitivity
  real(dp)               :: senlim ! Senescence limit
  integer                :: chill ! Chilling
  integer                :: dschill
  real(dp)               :: stemx
  real(dp)               :: gr0
  real(dp)               :: grf
  real(dp)               :: ppm0 ! Plant density
  real(dp),dimension(2)  :: sowthresh
  real(dp),dimension(2)  :: lethal
  real(dp),dimension(9)  :: cardinal
  real(dp),dimension(2)  :: croptype
  real(dp),dimension(6)  :: photoperiod
  real(dp),dimension(4)  :: croprange
  real(dp),dimension(6)  :: cropphen
  real(dp),dimension(3)  :: irrig
  integer,dimension(3)   :: sowday
  integer,dimension(2,3) :: cropgdd
  real(dp),dimension(6)  :: fert
  real(dp)               :: optlai
  real(dp)               :: harvindx
  integer                :: limdharv
end type

type (PftParameters) :: pft(max_cohorts), pft_tab(max_pftps)

end module pft_parameters
