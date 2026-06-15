!Copyright>        OpenRadioss
!Copyright>        Copyright (C) 1986-2026 Altair Engineering Inc.
!Copyright>
!Copyright>        This program is free software: you can redistribute it and/or modify
!Copyright>        it under the terms of the GNU Affero General Public License as published by
!Copyright>        the Free Software Foundation, either version 3 of the License, or
!Copyright>        (at your option) any later version.
!Copyright>
!Copyright>        This program is distributed in the hope that it will be useful,
!Copyright>        but WITHOUT ANY WARRANTY; without even the implied warranty of
!Copyright>        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!Copyright>        GNU Affero General Public License for more details.
!Copyright>
!Copyright>        You should have received a copy of the GNU Affero General Public License
!Copyright>        along with this program.  If not, see <https://www.gnu.org/licenses/>.
!Copyright>
!Copyright>
!Copyright>        Commercial Alternative: Altair Radioss Software
!Copyright>
!Copyright>        As an alternative to this open-source version, Altair also offers Altair Radioss
!Copyright>        software under a commercial license.  Contact Altair to discuss further if the
!Copyright>        commercial version may interest you: https://www.altair.com/radioss/.
!||====================================================================
!||    ists_sts_voxel_grid_mod  ../engine/source/interfaces/ists/ists_sts_voxel_grid_mod.F90
!||--- called by ------------------------------------------------------
!||    sts_broad_phase_voxel_mod   ../engine/source/interfaces/ists/ists_broad_phase_voxel.F90
!||====================================================================
!
!   Per-interface voxel broad-phase parameters (CELL_SIZE,
!   SEARCH_PADDING, N_CELL_RADIUS). Computed once on the first
!   broad-phase call per interface index NIN.
!
      MODULE ISTS_STS_VOXEL_GRID_MOD
        USE PRECISION_MOD, ONLY : WP
        IMPLICIT NONE
        PRIVATE

        TYPE :: STS_VOXEL_GRID_STATE
          LOGICAL :: initialized = .FALSE.
          REAL(KIND=WP) :: cell_size = 0.0_WP
          REAL(KIND=WP) :: search_padding = 0.0_WP
          REAL(KIND=WP) :: pad_sq = 0.0_WP
          INTEGER :: n_cell_radius = 1
        END TYPE STS_VOXEL_GRID_STATE

        TYPE(STS_VOXEL_GRID_STATE), ALLOCATABLE, SAVE :: STS_VOXEL_GRID(:)

        PUBLIC :: ISTS_STS_VOXEL_GRID_IS_READY
        PUBLIC :: ISTS_STS_VOXEL_GRID_GET
        PUBLIC :: ISTS_STS_VOXEL_GRID_SET

      CONTAINS
!=======================================================================
!   ISTS_STS_VOXEL_GRID_ENSURE_SIZE
!   Ensure the STS_VOXEL_GRID array is large enough to store the data for the given interface index NIN.
!=======================================================================
        SUBROUTINE ISTS_STS_VOXEL_GRID_ENSURE_SIZE(NIN)
          INTEGER, INTENT(IN) :: NIN
          TYPE(STS_VOXEL_GRID_STATE), ALLOCATABLE :: TMP(:)
          INTEGER :: OLD_SIZE, I

          IF (NIN <= 0) RETURN
          IF (.NOT. ALLOCATED(STS_VOXEL_GRID)) THEN
            ALLOCATE(STS_VOXEL_GRID(NIN))
            DO I = 1, NIN
              STS_VOXEL_GRID(I)%initialized = .FALSE.
              STS_VOXEL_GRID(I)%cell_size = 0.0_WP
              STS_VOXEL_GRID(I)%search_padding = 0.0_WP
              STS_VOXEL_GRID(I)%pad_sq = 0.0_WP
              STS_VOXEL_GRID(I)%n_cell_radius = 1
            END DO
            RETURN
          END IF
          OLD_SIZE = SIZE(STS_VOXEL_GRID)
          IF (NIN <= OLD_SIZE) RETURN
          ALLOCATE(TMP(NIN))
          IF (OLD_SIZE > 0) TMP(1:OLD_SIZE) = STS_VOXEL_GRID(1:OLD_SIZE)
          DO I = OLD_SIZE + 1, NIN
            TMP(I)%initialized = .FALSE.
            TMP(I)%cell_size = 0.0_WP
            TMP(I)%search_padding = 0.0_WP
            TMP(I)%pad_sq = 0.0_WP
            TMP(I)%n_cell_radius = 1
          END DO
          DEALLOCATE(STS_VOXEL_GRID)
          CALL MOVE_ALLOC(TMP, STS_VOXEL_GRID)
        END SUBROUTINE ISTS_STS_VOXEL_GRID_ENSURE_SIZE

!=======================================================================
!   ISTS_STS_VOXEL_GRID_IS_READY
!   Check if the STS_VOXEL_GRID array is initialized for the given interface index NIN.
!=======================================================================
        LOGICAL FUNCTION ISTS_STS_VOXEL_GRID_IS_READY(NIN)
          INTEGER, INTENT(IN) :: NIN

          ISTS_STS_VOXEL_GRID_IS_READY = .FALSE.
          IF (NIN <= 0) RETURN
          IF (.NOT. ALLOCATED(STS_VOXEL_GRID)) RETURN
          IF (NIN > SIZE(STS_VOXEL_GRID)) RETURN
          ISTS_STS_VOXEL_GRID_IS_READY = STS_VOXEL_GRID(NIN)%initialized
        END FUNCTION ISTS_STS_VOXEL_GRID_IS_READY

!=======================================================================
!   ISTS_STS_VOXEL_GRID_GET
!   Get the parameters for the given interface index NIN.
!=======================================================================
        SUBROUTINE ISTS_STS_VOXEL_GRID_GET(NIN, CELL_SIZE, SEARCH_PADDING, &
     &    PAD_SQ, N_CELL_RADIUS)
          INTEGER, INTENT(IN) :: NIN
          INTEGER, INTENT(OUT) :: N_CELL_RADIUS
          REAL(KIND=WP), INTENT(OUT) :: CELL_SIZE, SEARCH_PADDING, PAD_SQ

          CELL_SIZE = STS_VOXEL_GRID(NIN)%cell_size
          SEARCH_PADDING = STS_VOXEL_GRID(NIN)%search_padding
          PAD_SQ = STS_VOXEL_GRID(NIN)%pad_sq
          N_CELL_RADIUS = STS_VOXEL_GRID(NIN)%n_cell_radius
        END SUBROUTINE ISTS_STS_VOXEL_GRID_GET

!=======================================================================
!   ISTS_STS_VOXEL_GRID_SET
!   Set the parameters for the given interface index NIN.
!=======================================================================
        SUBROUTINE ISTS_STS_VOXEL_GRID_SET(NIN, CELL_SIZE, SEARCH_PADDING, &
     &    PAD_SQ, N_CELL_RADIUS)
          INTEGER, INTENT(IN) :: NIN, N_CELL_RADIUS
          REAL(KIND=WP), INTENT(IN) :: CELL_SIZE, SEARCH_PADDING, PAD_SQ

          CALL ISTS_STS_VOXEL_GRID_ENSURE_SIZE(NIN)
          STS_VOXEL_GRID(NIN)%cell_size = CELL_SIZE
          STS_VOXEL_GRID(NIN)%search_padding = SEARCH_PADDING
          STS_VOXEL_GRID(NIN)%pad_sq = PAD_SQ
          STS_VOXEL_GRID(NIN)%n_cell_radius = N_CELL_RADIUS
          STS_VOXEL_GRID(NIN)%initialized = .TRUE.
        END SUBROUTINE ISTS_STS_VOXEL_GRID_SET

      END MODULE ISTS_STS_VOXEL_GRID_MOD
