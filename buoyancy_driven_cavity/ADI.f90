module ADI
  use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
  use variables_globales
  implicit none

  ! coeficiente "P", "Q"
  real(dp), dimension(nx,ny) :: P,Q

  ! Vartiables de apoyo
  real(dp) :: aephi,anphi,asphi,awphi

  public :: ADI_
  private

  contains

  ! FIXME: las temperaturas en las fronteras de dirichlet arrojan valores
  ! negativos. Las fronteras aisladas muestran valores negativos en relacion
  ! a los valores del interior del dominio

  subroutine ADI_(ap,ae,aw,an,as,b,phinew)
    real(dp), dimension(nx,ny), intent(inout) :: ap,ae,aw,an,as,b
    real(dp), dimension(nx,ny), intent(inout) :: phinew

    ! Modificacion de los coeficiences de en medio
    do j=2,ny-1
      do i=2,nx-1
        ae(i,j)=-1.0_dp*ae(i,j)
        aw(i,j)=-1.0_dp*aw(i,j)
        an(i,j)=-1.0_dp*an(i,j)
        as(i,j)=-1.0_dp*as(i,j)
      end do
    end do

    ! Barrido por LBL-X
    do j=1,ny
      P(1,j)=ae(1,j)/ap(1,j)

      Q(1,j)=(an(1,j)*phinew(1,j+1)+as(1,j)*phinew(1,j-1)+b(1,j))/ap(1,j)

      do i=2,nx

        if (j < ny) then
          anphi=an(i,j)*phinew(i,j+1)
        else
          anphi=0.0_dp
        end if

        if (j==1) then
          asphi=0.0_dp
        else
          asphi=as(i,j)*phinew(i,j-1)
        end if

        P(i,j)=ae(i,j)/(ap(i,j)-aw(i,j)*P(i-1,j))
        Q(i,j)=(anphi+asphi+ &
        b(i,j)+aw(i,j)*Q(i-1,j))/(ap(i,j)-aw(i,j)*P(i-1,j))

      end do

      phinew(nx,j)=Q(nx,j)

      do i=nx-1,1,-1
        phinew(i,j)=phinew(i+1,j)*P(i,j)+Q(i,j)
      end do
   end do

   ! Barrido por LBL-Y
    do i=1,nx
      P(i,1)=an(i,1)/ap(i,1)

      Q(i,1)=(ae(i,1)*phinew(i+1,1)+aw(i,1)*phinew(i-1,1)+b(i,1))/ap(i,1)

      do j=2,ny

        if (i < nx) then
          aephi=ae(i,j)*phinew(i+1,j)
        else
          aephi=0.0_dp
        end if

        if (i==1) then
          awphi=0.0_dp
        else
          awphi=aw(i,j)*phinew(i-1,j)
        end if

        P(i,j)=an(i,j)/(ap(i,j)-as(i,j)*P(i,j-1))
        Q(i,j)=(aephi+awphi+ &
        b(i,j)+as(i,j)*Q(i,j-1))/(ap(i,j)-as(i,j)*P(i,j-1))

      end do

      phinew(i,ny)=Q(i,ny)

      do j=ny-1,1,-1
        phinew(i,j)=phinew(i,j+1)*P(i,j)+Q(i,j)
      end do
    end do


  end subroutine ADI_

end module ADI
