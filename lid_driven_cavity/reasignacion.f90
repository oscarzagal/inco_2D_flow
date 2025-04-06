module reasignacion
  use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
  use variables_globales
  use sobrecarga_operadores
  implicit none

  public :: reasignacion_variables
  private

  contains

  subroutine reasignacion_variables()

    Pold=Pstar
    u_old=u_star
    v_old=v_star
    me_old=me_star
    mw_old=mw_star
    mn_old=mn_star
    ms_old=ms_star

    ue_n=ue
    uw_n=uw
    vn_n=vn
    vs_n=vs

    ue_i_n=ue_i
    uw_i_n=uw_i
    vn_i_n=vn_i
    vs_i_n=vs_i

  end subroutine reasignacion_variables

end module reasignacion
