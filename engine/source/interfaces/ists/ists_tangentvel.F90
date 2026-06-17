!||====================================================================
!||    sts_tangentvel_global  ../engine/source/interfaces/ists/ists_tangentvel.F90
!||--- called by ------------------------------------------------------
!||    STS_CONTACT_EVAL_PAIR   ../engine/source/interfaces/ists/ists_contact_eval_pair.F90
!||--- calls ---------------------------------------------------------
!||    sts_handle_element_transition  ../engine/source/interfaces/ists/ists_elemTrans.F90
!||====================================================================
!-----------------------------------------------
!   Global convective increments with period tracking
!-----------------------------------------------
      subroutine sts_tangentvel_global(xi1, xi2, dxi1, dxi2, &
     &     EL_NR, z, q, ip, MAX_STS_SIZE_ACTUAL)

      use sts_gp_state_mod
      implicit none

!-- Input arguments
      real*8  xi1, xi2
      integer EL_NR, z, q, ip, MAX_STS_SIZE_ACTUAL

!-- Local variables
      real*8  dxi1, dxi2
      integer gp_index
      real*8  xi1_prev_local, xi2_prev_local
      real*8  tol
      INTEGER GET_GLOBAL_GP_INDEX

      tol = 1.0d-6

      dxi1 = 0.d0
      dxi2 = 0.d0

      if (EL_NR .LE. 0 .OR. EL_NR .GT. MAX_STS_SIZE_ACTUAL) return
      if (z .LE. 0 .OR. z .GT. ip .OR. q .LE. 0 .OR. q .GT. ip) return

      gp_index = GET_GLOBAL_GP_INDEX(EL_NR, z, q, ip)

      if (gp_index .LE. 0 .OR. gp_index .GT. MAX_GLOBAL_GP) then
        ! Invalid index - return zero increments
        return
      endif

      if (.NOT. GP_INITIALIZED(gp_index)) then
!       First contact: initialize global state and reset traction history
        GP_XI1_GLOBAL(gp_index)      = xi1
        GP_XI2_GLOBAL(gp_index)      = xi2
        GP_XI1_GLOBAL_PREV(gp_index) = xi1
        GP_XI2_GLOBAL_PREV(gp_index) = xi2
        GP_XI1_PERIOD(gp_index)      = 0
        GP_XI2_PERIOD(gp_index)      = 0

        GP_TTRIAL1_HIST(gp_index)    = 0.d0
        GP_TTRIAL2_HIST(gp_index)    = 0.d0

        GP_INITIALIZED(gp_index)     = .TRUE.
        dxi1 = 0.d0
        dxi2 = 0.d0
      else
!       Reconstruct previous local coordinates from global + period
        xi1_prev_local = GP_XI1_GLOBAL(gp_index) &
     &                  - 2.0d0*GP_XI1_PERIOD(gp_index)
        xi2_prev_local = GP_XI2_GLOBAL(gp_index) &
     &                  - 2.0d0*GP_XI2_PERIOD(gp_index)

!       Update period counters based on boundary crossing
        call sts_handle_element_transition(xi1_prev_local, xi1, &
     &       xi2_prev_local, xi2, &
     &       GP_XI1_PERIOD(gp_index), GP_XI2_PERIOD(gp_index), tol)

!       Update global coordinates
        GP_XI1_GLOBAL_PREV(gp_index) = GP_XI1_GLOBAL(gp_index)
        GP_XI2_GLOBAL_PREV(gp_index) = GP_XI2_GLOBAL(gp_index)

        GP_XI1_GLOBAL(gp_index) = xi1 + 2.0d0*GP_XI1_PERIOD(gp_index)
        GP_XI2_GLOBAL(gp_index) = xi2 + 2.0d0*GP_XI2_PERIOD(gp_index)

        dxi1 = GP_XI1_GLOBAL(gp_index) - GP_XI1_GLOBAL_PREV(gp_index)
        dxi2 = GP_XI2_GLOBAL(gp_index) - GP_XI2_GLOBAL_PREV(gp_index)

!       Check if dxi_global is suspiciously large (indicates lost contact)
!       If dxi > 0.2, treat as new contact and reinitialize
!       Note: Normal dxi should be very small (e.g., 1e-7 to 1e-6 per timestep)
        if (ABS(dxi1) > 0.2d0 .OR. ABS(dxi2) > 0.2d0) then
!         Reinitialize as new contact
          GP_XI1_GLOBAL(gp_index)      = xi1
          GP_XI2_GLOBAL(gp_index)      = xi2
          GP_XI1_GLOBAL_PREV(gp_index) = xi1
          GP_XI2_GLOBAL_PREV(gp_index) = xi2
          GP_XI1_PERIOD(gp_index)      = 0
          GP_XI2_PERIOD(gp_index)      = 0
          GP_TTRIAL1_HIST(gp_index)    = 0.d0
          GP_TTRIAL2_HIST(gp_index)    = 0.d0
          dxi1 = 0.d0
          dxi2 = 0.d0
        endif
      endif  ! End of if (.NOT. GP_INITIALIZED)

      return
      end

!||====================================================================
!||    sts_gp_tangential_velocity  ../engine/source/interfaces/ists/ists_tangentvel.F90
!||--- called by ------------------------------------------------------
!||    STS_CONTACT_EVAL_PAIR   ../engine/source/interfaces/ists/ists_contact_eval_pair.F90
!||====================================================================
!-----------------------------------------------
!   Relative tangential velocity at a contact GP
!   (secondary minus primary, projected onto tangent plane).
!-----------------------------------------------
      subroutine sts_gp_tangential_velocity(N_xi, N_eta, node_ids, V, &
     &     norm_contact, v_tang)

#include      "my_real.inc"
      implicit none

      real*8  N_xi(3,4), N_eta(3,4)
      integer node_ids(8)
      my_real V(3,*)
      real*8  norm_contact(3)
      real*8  v_tang(3)

      integer j
      real*8  v_prim(3), v_sec(3), v_rel(3), v_n

      v_prim = 0.d0
      v_sec  = 0.d0
      DO j = 1, 4
        v_prim(1) = v_prim(1) + N_xi(1,j)  * DBLE(V(1, node_ids(j)))
        v_prim(2) = v_prim(2) + N_xi(1,j)  * DBLE(V(2, node_ids(j)))
        v_prim(3) = v_prim(3) + N_xi(1,j)  * DBLE(V(3, node_ids(j)))
        v_sec(1)  = v_sec(1)  + N_eta(1,j) * DBLE(V(1, node_ids(j+4)))
        v_sec(2)  = v_sec(2)  + N_eta(1,j) * DBLE(V(2, node_ids(j+4)))
        v_sec(3)  = v_sec(3)  + N_eta(1,j) * DBLE(V(3, node_ids(j+4)))
      ENDDO

      v_rel(1) = v_sec(1) - v_prim(1)
      v_rel(2) = v_sec(2) - v_prim(2)
      v_rel(3) = v_sec(3) - v_prim(3)

      v_n = v_rel(1)*norm_contact(1) + v_rel(2)*norm_contact(2) &
     &    + v_rel(3)*norm_contact(3)

      v_tang(1) = v_rel(1) - v_n*norm_contact(1)
      v_tang(2) = v_rel(2) - v_n*norm_contact(2)
      v_tang(3) = v_rel(3) - v_n*norm_contact(3)

      return
      end