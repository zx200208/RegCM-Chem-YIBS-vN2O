!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    ICTP RegCM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

module mod_che_carbonaer
  
  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_constants
  use mod_che_common
  use mod_che_species
  use mod_che_indices

  private

  ! Parameter usefull for wet and dry deposition of carbon aerosol 
  ! densities in kg/m3

  real(rk8) , public , parameter :: rhobc   = 2000.0D0
  real(rk8) , public , parameter :: rhooc   = 1200.0D0
  real(rk8) , public , parameter :: rhobchl = 1600.0D0
  real(rk8) , public , parameter :: rhoochl = 1200.0D0

  ! S.Nakao et al. / Atmospheric Environment 68(2013) 273-277
  ! Use it carefully!!

  real(rk8) , public , parameter :: rhoaor1 = 1400.0D0  !ycq2013
  real(rk8) , public , parameter :: rhoaor2 = 1400.0D0  !ycq2013
  real(rk8) , public , parameter :: rhoaor3 = 1400.0D0  !ycq2013
  real(rk8) , public , parameter :: rhoaor4 = 1400.0D0  !ycq2013
  real(rk8) , public , parameter :: rhobor1 = 1300.0D0  !ycq2013
  real(rk8) , public , parameter :: rhobor2 = 1300.0D0  !ycq2013
  real(rk8) , public , parameter :: rhobor3 = 1300.0D0  !ycq2013
  real(rk8) , public , parameter :: rhobor4 = 1300.0D0  !ycq2013
  real(rk8) , public , parameter :: rhoisoa1 = 1200.0D0  !ycq2013
  real(rk8) , public , parameter :: rhoisoa2 = 1200.0D0  !ycq2013
  real(rk8) , public , parameter :: rhoisoa3 = 1200.0D0  !ycq2013
  real(rk8) , public , parameter :: rhoisoa4 = 1200.0D0  !ycq2013

  ! effctive dimaters ( and not radius!)  in micrometer
  ! ( should they be defined intercatively in the future ? ) 
  real(rk8) , public , parameter :: reffbc   = 0.05D0
  real(rk8) , public , parameter :: reffbchl = 0.3D0
  real(rk8) , public , parameter :: reffoc   = 0.2D0
  real(rk8) , public , parameter :: reffochl = 0.3D0

  real(rk8) , public , parameter :: reffaor1 = 0.2D0    !ycq2013
  real(rk8) , public , parameter :: reffaor2 = 0.2D0    !ycq2013
  real(rk8) , public , parameter :: reffaor3 = 0.2D0    !ycq2013
  real(rk8) , public , parameter :: reffaor4 = 0.2D0    !ycq2013
  real(rk8) , public , parameter :: reffbor1 = 0.2D0    !ycq2013
  real(rk8) , public , parameter :: reffbor2 = 0.2D0    !ycq2013
  real(rk8) , public , parameter :: reffbor3 = 0.2D0    !ycq2013
  real(rk8) , public , parameter :: reffbor4 = 0.2D0    !ycq2013
  real(rk8) , public , parameter :: reffisoa1 = 0.3D0    !ycq2013
  real(rk8) , public , parameter :: reffisoa2 = 0.3D0    !ycq2013
  real(rk8) , public , parameter :: reffisoa3 = 0.3D0    !ycq2013
  real(rk8) , public , parameter :: reffisoa4 = 0.3D0    !ycq2013

  ! aging efolding time (s), from hydrophobic to hydrophilic
  ! Cooke et al.
  real(rk8) , parameter :: chagct = 1.15D0 * 86400.0D0
  !
  ! solubility of carbon aer for rain out param of giorgi and chameides
  !
  real(rk8) , parameter :: solbc = 0.05D0
  real(rk8) , parameter :: solbchl = 0.8D0
  real(rk8) , parameter :: soloc = 0.05D0
  real(rk8) , parameter :: solochl = 0.8D0

  real(rk8) , parameter :: solaor1 = 0.05D0              !ycq2013
  real(rk8) , parameter :: solaor2 = 0.05D0              !ycq2013
  real(rk8) , parameter :: solaor3 = 0.05D0              !ycq2013
  real(rk8) , parameter :: solaor4 = 0.05D0              !ycq2013
  real(rk8) , parameter :: solbor1 = 0.6D0               !ycq2013
  real(rk8) , parameter :: solbor2 = 0.6D0               !ycq2013
  real(rk8) , parameter :: solbor3 = 0.6D0               !ycq2013
  real(rk8) , parameter :: solbor4 = 0.6D0               !ycq2013
  real(rk8) , parameter :: solisoa1 = 0.8D0              !ycq2013
  real(rk8) , parameter :: solisoa2 = 0.8D0              !ycq2013
  real(rk8) , parameter :: solisoa3 = 0.8D0              !ycq2013
  real(rk8) , parameter :: solisoa4 = 0.8D0              !ycq2013

  public :: aging_carb , solbc , solbchl , soloc , solochl ,   &   !ycq2013
            solaor1 , solaor2 , solaor3 , solaor4 , solbor1 ,  &
            solbor2 , solbor3 , solbor4 , solisoa1 , solisoa2, &
            solisoa3 , solisoa4

  ! bin size for carboneaceous aerosols
  ! ps add one dimension for sulfate too.
  real(rk8) , public , dimension(17) :: carbed     !ycq2013 from 5 to 17

  contains

    subroutine aging_carb(j)
      implicit none
      integer, intent(in) :: j
      integer(ik4) :: i , k
      real(rk8) :: agingtend1 , agingtend2
      !
      ! aging o carbon species : Conversion from hydrophobic to 
      ! hydrophilic: Carbonaceopus species time constant
      ! ( 1.15 day cooke et al.,1999)
      !
      if ( ibchb > 0 .and. ibchl > 0 ) then
        do k = 1 , kz
          do i = ici1 , ici2
            agingtend1 = -chib(j,i,k,ibchb)*(d_one-dexp(-dt/chagct))/dt
            agingtend2 = -agingtend1
            chiten(j,i,k,ibchb) = chiten(j,i,k,ibchb) + agingtend1
            chiten(j,i,k,ibchl) = chiten(j,i,k,ibchl) + agingtend2
          end do
        end do
      end if
      if ( iochb > 0  .and. iochl > 0 ) then
        do k = 1 , kz
          do i = ici1 , ici2
            agingtend1 = -chib(j,i,k,iochb)*(d_one-dexp(-dt/chagct))/dt
            agingtend2 = -agingtend1
            chiten(j,i,k,iochb) = chiten(j,i,k,iochb) + agingtend1
            chiten(j,i,k,iochl) = chiten(j,i,k,iochl) + agingtend2
          end do
        end do
      end if
    end subroutine aging_carb
!
end module mod_che_carbonaer
