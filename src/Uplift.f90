subroutine Uplift ()

  ! subroutine to apply an uplift step using the uplift function/array u

  use FastScapeContext

  implicit none

  h = h + u*dt
  b = b + u*dt

  ! Enforce a very low sealevel if the sea is very shallow and if option low_sealevel_at_shallow_sea is used.
  ! If the depth of the sea is lower than the active layer, the code will not converge.
  ! This is important for coupling with thermo-mechanical models, where this scenario can easily happen
  ! Note that you need to set sealevel at every timestep if you want to use this option.
  if (minval(h) > (sealevel - layer) .and. low_sealevel_at_shallow_sea  ) then
      sealevel = -huge(sealevel)
      write(*,*)'sealevel set to low - shallow sea: ', sealevel
  endif

  return

  end subroutine Uplift
