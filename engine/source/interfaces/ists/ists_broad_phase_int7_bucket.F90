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
!||    sts_broad_phase_int7_bucket_mod   ../engine/source/interfaces/ists/ists_broad_phase_int7_bucket.F90
!||--- called by ------------------------------------------------------
!||    ists_mainf              ../engine/source/interfaces/ists/ists_mainf.F
!||--- calls      -----------------------------------------------------
!||    sts_remap_segments      ../engine/source/interfaces/ists/ists_remap_segments.F90
!||--- uses       -----------------------------------------------------
!||    intbufdef_mod           ../common_source/modules/interfaces/intbufdef_mod.F90
!||    groupdef_mod            ../engine/share/modules/groupdef_mod.F
!||====================================================================
!
!   Legacy INT7 bucket broad-phase adapter for STS contact.
!
!   Instead of running the STS-native voxel search, this path consumes
!   the candidate pairs that the legacy INT7 sorting (I7BUCE/I7TRI) has
!   already stored in INTBUF_TAB%CAND_N / CAND_E (node <-> main segment)
!   and maps them to STS segment pairs through STS_REMAP_SEGMENTS.
!
!   The output (CAND_SEC_SEG_ID, CAND_MST_SEG_ID, CONT_ELEMENT, COUNT)
!   matches STS_VOXEL_BROAD_PHASE so the downstream STS pipeline
!   (STS_CONTACT_STIFFNESS, STS_CONTACTS_ASSEMBLE) is unchanged.
!
      MODULE STS_BROAD_PHASE_INT7_BUCKET_MOD
        USE INTBUFDEF_MOD, ONLY : INTBUF_STRUCT_
        USE GROUPDEF_MOD,  ONLY : SURF_
        USE CONSTANT_MOD,  ONLY : ZERO
        USE ISTS_STS_BP_PERSIST_MOD, ONLY : ISTS_STS_BP_PERSIST_SAVE, &
     &    ISTS_STS_BP_PERSIST_TOUCH, ISTS_STS_BP_PERSIST_TRY_RESTORE
        IMPLICIT NONE
        PRIVATE
        PUBLIC :: STS_INT7_BUCKET_BROAD_PHASE
      CONTAINS
!=======================================================================
!   STS_INT7_BUCKET_BROAD_PHASE
!
!   Build STS segment pairs from the legacy INT7 candidate arrays.
!   COUNT is clamped to MAX_STS_SIZE_ACTUAL; OVERFLOW is set when the
!   storage saturated so the caller can grow capacity and retry.
!=======================================================================
        SUBROUTINE STS_INT7_BUCKET_BROAD_PHASE( &
     &      NIN, INTBUF_TAB, IGRSURF, NSURF, SEC_SURF_IDX, MST_SURF_IDX, &
     &      X, NUMNOD, NSN, NRTM, MAX_STS_SIZE_ACTUAL, &
     &      CAND_SEC_SEG_ID, CAND_MST_SEG_ID, CAND_SEC_GP_MASK, &
     &      CONT_ELEMENT, &
     &      COUNT, OVERFLOW, D_MIN)
!-----------------------------------------------
!   I m p l i c i t   T y p e s
!-----------------------------------------------
#include      "my_real.inc"
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
          INTEGER, INTENT(IN)  :: NIN
          TYPE(INTBUF_STRUCT_) :: INTBUF_TAB
          INTEGER, INTENT(IN)  :: NSURF
          TYPE(SURF_), DIMENSION(NSURF), INTENT(IN) :: IGRSURF
          INTEGER, INTENT(IN)  :: SEC_SURF_IDX, MST_SURF_IDX
          INTEGER, INTENT(IN)  :: NUMNOD
          INTEGER, INTENT(IN)  :: NSN
          INTEGER, INTENT(IN)  :: NRTM
          INTEGER, INTENT(IN)  :: MAX_STS_SIZE_ACTUAL
          my_real, INTENT(IN)  :: X(3, NUMNOD)
          INTEGER, INTENT(OUT) :: CAND_SEC_SEG_ID(MAX_STS_SIZE_ACTUAL, 5)
          INTEGER, INTENT(OUT) :: CAND_MST_SEG_ID(MAX_STS_SIZE_ACTUAL, 5)
          INTEGER, INTENT(OUT) :: CAND_SEC_GP_MASK(MAX_STS_SIZE_ACTUAL, 4)
          my_real, INTENT(OUT) :: CONT_ELEMENT(MAX_STS_SIZE_ACTUAL, 3, 8)
          INTEGER, INTENT(OUT) :: COUNT
          LOGICAL, INTENT(OUT) :: OVERFLOW
          my_real, INTENT(OUT) :: D_MIN
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
          INTEGER, PARAMETER :: PERSIST_EXTRA_LIMIT = 4
          INTEGER :: I_STOK, I, J, NSEC_BOUNDS, CAND_N_ABS, N_VALID
          INTEGER :: COUNT_FRESH, COUNT_PERSIST, PERSIST_I, FRESH_I
          INTEGER :: OVERLAP_COUNT, MISSING_COUNT, FINAL_LIMIT
          INTEGER, ALLOCATABLE :: CAND_SEC_SEG(:)
          INTEGER, ALLOCATABLE :: CAND_N_COMPACT(:), CAND_E_COMPACT(:)
          INTEGER, ALLOCATABLE :: FRESH_SEC_ID(:,:)
          INTEGER, ALLOCATABLE :: FRESH_MST_ID(:,:)
          INTEGER, ALLOCATABLE :: PERSIST_SEC_ID(:,:)
          INTEGER, ALLOCATABLE :: PERSIST_MST_ID(:,:)
          my_real, ALLOCATABLE :: FRESH_CONT_ELEMENT(:,:,:)
          my_real, ALLOCATABLE :: PERSIST_CONT_ELEMENT(:,:,:)
          LOGICAL :: PERSIST_RESTORED
          LOGICAL :: PERSIST_AUGMENTED, PERSIST_STABILIZE, PAIR_EXISTS
!-----------------------------------------------
!   S o u r c e   L i n e s
!-----------------------------------------------
          COUNT = 0
          N_VALID = 0
          OVERFLOW = .FALSE.
!         Skip is disabled for the legacy path; D_MIN is only used by the
!         voxel adaptive-skip logic, so a neutral value is returned.
          D_MIN = ZERO
!
          IF (MAX_STS_SIZE_ACTUAL <= 0) RETURN
          CAND_SEC_GP_MASK = 0
          IF (SEC_SURF_IDX <= 0 .OR. SEC_SURF_IDX > NSURF) RETURN
          IF (MST_SURF_IDX <= 0 .OR. MST_SURF_IDX > NSURF) RETURN
!
!         Number of candidate (node, main segment) pairs produced by the
!         INT7 sorting for this cycle. If the sorting was skipped this
!         cycle (e.g. distance criterion), no fresh candidates exist.
          I_STOK = INTBUF_TAB%I_STOK(1)
          IF (I_STOK <= 0) THEN
            ! Try to restore the last successful STS segment pairs from the cache.
            CALL ISTS_STS_BP_PERSIST_TRY_RESTORE( &
     &          NIN, X, NUMNOD, MAX_STS_SIZE_ACTUAL, &
     &          COUNT, CAND_SEC_SEG_ID, CAND_MST_SEG_ID, CONT_ELEMENT, &
     &          PERSIST_RESTORED)
            IF (PERSIST_RESTORED) THEN
              IF (COUNT >= MAX_STS_SIZE_ACTUAL) OVERFLOW = .TRUE.
              IF (COUNT > 0) CAND_SEC_GP_MASK(1:COUNT, 1:4) = 1
              RETURN
            END IF
            RETURN
          ENDIF
!
          NSEC_BOUNDS = NSN
          IF (INTBUF_TAB%S_NSV > 0) THEN
            NSEC_BOUNDS = INTBUF_TAB%S_NSV
          ENDIF
!
!         Keep INT7 bucket slots which reached the legacy active-force
!         state.  Testing positive fresh bucket slots here over-includes
!         candidate patches and destabilizes the STS force integration.
          N_VALID = 0
          ALLOCATE(CAND_N_COMPACT(I_STOK), CAND_E_COMPACT(I_STOK))
          DO I = 1, I_STOK
            IF (INTBUF_TAB%CAND_E(I) <= 0 .OR. &
     &          INTBUF_TAB%CAND_E(I) > NRTM) CYCLE
            IF (INTBUF_TAB%CAND_N(I) >= 0) CYCLE
            CAND_N_ABS = -INTBUF_TAB%CAND_N(I)
            IF (CAND_N_ABS > 0 .AND. &
     &              CAND_N_ABS <= NSEC_BOUNDS) THEN
              N_VALID = N_VALID + 1
              CAND_N_COMPACT(N_VALID) = CAND_N_ABS
            ELSE
              CYCLE
            ENDIF
            CAND_E_COMPACT(N_VALID) = INTBUF_TAB%CAND_E(I)
          END DO
          IF (N_VALID <= 0) THEN
            DEALLOCATE(CAND_N_COMPACT)
            DEALLOCATE(CAND_E_COMPACT)
            ! Try to restore the last successful STS segment pairs from the cache.
            CALL ISTS_STS_BP_PERSIST_TRY_RESTORE( &
     &          NIN, X, NUMNOD, MAX_STS_SIZE_ACTUAL, &
     &          COUNT, CAND_SEC_SEG_ID, CAND_MST_SEG_ID, CONT_ELEMENT, &
     &          PERSIST_RESTORED)
            IF (PERSIST_RESTORED) THEN
              IF (COUNT >= MAX_STS_SIZE_ACTUAL) OVERFLOW = .TRUE.
              IF (COUNT > 0) CAND_SEC_GP_MASK(1:COUNT, 1:4) = 1
              RETURN
            END IF
            RETURN
          ENDIF
!
          ALLOCATE(CAND_SEC_SEG(MAX_STS_SIZE_ACTUAL))
!
          CALL STS_REMAP_SEGMENTS( &
     &        INTBUF_TAB, X, NUMNOD, NRTM, NSN, CAND_SEC_SEG, &
     &        N_VALID, CAND_N_COMPACT, CAND_E_COMPACT, &
     &        INTBUF_TAB%IRECTM, CONT_ELEMENT, COUNT, &
     &        IGRSURF, CAND_SEC_SEG_ID, CAND_MST_SEG_ID, &
     &        CAND_SEC_GP_MASK, &
     &        MAX_STS_SIZE_ACTUAL, NSURF, SEC_SURF_IDX, MST_SURF_IDX)
!
          COUNT_FRESH = COUNT
          PERSIST_AUGMENTED = .FALSE.
          IF (COUNT_FRESH > 0) THEN
            ALLOCATE(PERSIST_SEC_ID(MAX_STS_SIZE_ACTUAL, 5))
            ALLOCATE(PERSIST_MST_ID(MAX_STS_SIZE_ACTUAL, 5))
            ALLOCATE(PERSIST_CONT_ELEMENT(MAX_STS_SIZE_ACTUAL, 3, 8))
            ! Try to restore the last successful STS segment pairs from the cache.
            CALL ISTS_STS_BP_PERSIST_TRY_RESTORE( &
     &        NIN, X, NUMNOD, MAX_STS_SIZE_ACTUAL, &
     &        COUNT_PERSIST, PERSIST_SEC_ID, PERSIST_MST_ID, &
     &        PERSIST_CONT_ELEMENT, PERSIST_RESTORED)

!           A partial INT7 candidate wipe is just as damaging to STS as a
!           full COUNT=0 wipe: the integrated patch loses support abruptly.
!           Compare partner identities, not only pair counts.  If the
!           fresh bucket set drops cached partners or expands abruptly,
!           rebuild from the last stable patch and add only a few new
!           fresh pairs. Projection/penetration still decides which pairs
!           actually carry force.
            IF (PERSIST_RESTORED .AND. COUNT_PERSIST > 0) THEN
              OVERLAP_COUNT = 0
              MISSING_COUNT = 0
              DO PERSIST_I = 1, COUNT_PERSIST
                PAIR_EXISTS = .FALSE.
                DO J = 1, COUNT_FRESH
                  IF (CAND_SEC_SEG_ID(J, 1) == &
     &                PERSIST_SEC_ID(PERSIST_I, 1) .AND. &
     &                CAND_MST_SEG_ID(J, 1) == &
     &                PERSIST_MST_ID(PERSIST_I, 1)) THEN
                    PAIR_EXISTS = .TRUE.
                    EXIT
                  ENDIF
                ENDDO
                IF (PAIR_EXISTS) THEN
                  OVERLAP_COUNT = OVERLAP_COUNT + 1
                ELSE
                  MISSING_COUNT = MISSING_COUNT + 1
                ENDIF
              ENDDO
              PERSIST_STABILIZE = MISSING_COUNT > 0 .OR. &
     &          COUNT_FRESH > COUNT_PERSIST + PERSIST_EXTRA_LIMIT
              IF (PERSIST_STABILIZE .AND. OVERLAP_COUNT > 0) THEN
                ALLOCATE(FRESH_SEC_ID(COUNT_FRESH, 5))
                ALLOCATE(FRESH_MST_ID(COUNT_FRESH, 5))
                ALLOCATE(FRESH_CONT_ELEMENT(COUNT_FRESH, 3, 8))
                FRESH_SEC_ID(1:COUNT_FRESH, 1:5) = &
     &            CAND_SEC_SEG_ID(1:COUNT_FRESH, 1:5)
                FRESH_MST_ID(1:COUNT_FRESH, 1:5) = &
     &            CAND_MST_SEG_ID(1:COUNT_FRESH, 1:5)
                FRESH_CONT_ELEMENT(1:COUNT_FRESH, 1:3, 1:8) = &
     &            CONT_ELEMENT(1:COUNT_FRESH, 1:3, 1:8)
                COUNT = 0
                FINAL_LIMIT = MIN(MAX_STS_SIZE_ACTUAL, &
     &            COUNT_PERSIST + PERSIST_EXTRA_LIMIT)
                DO PERSIST_I = 1, COUNT_PERSIST
                  IF (COUNT >= MAX_STS_SIZE_ACTUAL) THEN
                    OVERFLOW = .TRUE.
                    EXIT
                  ENDIF
                  COUNT = COUNT + 1
                  CAND_SEC_SEG_ID(COUNT, 1:5) = &
     &              PERSIST_SEC_ID(PERSIST_I, 1:5)
                  CAND_MST_SEG_ID(COUNT, 1:5) = &
     &              PERSIST_MST_ID(PERSIST_I, 1:5)
                  CAND_SEC_GP_MASK(COUNT, 1:4) = 1
                  CONT_ELEMENT(COUNT, 1:3, 1:8) = &
     &              PERSIST_CONT_ELEMENT(PERSIST_I, 1:3, 1:8)
                ENDDO
                DO FRESH_I = 1, COUNT_FRESH
                  IF (COUNT >= FINAL_LIMIT) EXIT
                  PAIR_EXISTS = .FALSE.
                  DO J = 1, COUNT
                    IF (FRESH_SEC_ID(FRESH_I, 1) == &
     &                  CAND_SEC_SEG_ID(J, 1) .AND. &
     &                  FRESH_MST_ID(FRESH_I, 1) == &
     &                  CAND_MST_SEG_ID(J, 1)) THEN
                      PAIR_EXISTS = .TRUE.
                      EXIT
                    ENDIF
                  ENDDO
                  IF (PAIR_EXISTS) CYCLE
                  COUNT = COUNT + 1
                  CAND_SEC_SEG_ID(COUNT, 1:5) = &
     &              FRESH_SEC_ID(FRESH_I, 1:5)
                  CAND_MST_SEG_ID(COUNT, 1:5) = &
     &              FRESH_MST_ID(FRESH_I, 1:5)
                  CAND_SEC_GP_MASK(COUNT, 1:4) = 1
                  CONT_ELEMENT(COUNT, 1:3, 1:8) = &
     &              FRESH_CONT_ELEMENT(FRESH_I, 1:3, 1:8)
                ENDDO
                DEALLOCATE(FRESH_SEC_ID)
                DEALLOCATE(FRESH_MST_ID)
                DEALLOCATE(FRESH_CONT_ELEMENT)
                PERSIST_AUGMENTED = .TRUE.
              ENDIF
            ENDIF
            DEALLOCATE(PERSIST_SEC_ID)
            DEALLOCATE(PERSIST_MST_ID)
            DEALLOCATE(PERSIST_CONT_ELEMENT)
          ENDIF
!
          DEALLOCATE(CAND_N_COMPACT)
          DEALLOCATE(CAND_E_COMPACT)
!
!         Saturated storage: signal the caller to grow capacity and retry.
          IF (COUNT >= MAX_STS_SIZE_ACTUAL) OVERFLOW = .TRUE.
!
          IF (COUNT > 0) THEN
            IF (PERSIST_AUGMENTED) THEN
              ! Touch the cached STS segment pairs to update the last save cycle number.
              CALL ISTS_STS_BP_PERSIST_TOUCH(NIN)
            ELSE
              ! Save the current STS segment pairs to the cache.
              CALL ISTS_STS_BP_PERSIST_SAVE(NIN, COUNT, CAND_SEC_SEG_ID, &
     &            CAND_MST_SEG_ID, MAX_STS_SIZE_ACTUAL)
            ENDIF
          ELSE
            ! Try to restore the last successful STS segment pairs from the cache.
            CALL ISTS_STS_BP_PERSIST_TRY_RESTORE( &
     &          NIN, X, NUMNOD, MAX_STS_SIZE_ACTUAL, &
     &          COUNT, CAND_SEC_SEG_ID, CAND_MST_SEG_ID, CONT_ELEMENT, &
     &          PERSIST_RESTORED)
            IF (PERSIST_RESTORED) THEN
              IF (COUNT >= MAX_STS_SIZE_ACTUAL) OVERFLOW = .TRUE.
              IF (COUNT > 0) CAND_SEC_GP_MASK(1:COUNT, 1:4) = 1
              DEALLOCATE(CAND_SEC_SEG)
              RETURN
            END IF
          END IF
!
          DEALLOCATE(CAND_SEC_SEG)
        END SUBROUTINE STS_INT7_BUCKET_BROAD_PHASE
      END MODULE STS_BROAD_PHASE_INT7_BUCKET_MOD
