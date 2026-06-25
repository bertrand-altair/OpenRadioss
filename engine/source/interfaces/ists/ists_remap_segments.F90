!||====================================================================
!||    sts_remap_segments  ../engine/source/interfaces/ists/ists_remap_segments.F90
!||--- called by ------------------------------------------------------
!||    sts_int7_bucket_broad_phase  ../engine/source/interfaces/ists/ists_broad_phase_int7_bucket.F90
!||--- calls ---------------------------------------------------------
!||    (none - local mapping only)
!||====================================================================
      SUBROUTINE STS_REMAP_SEGMENTS(INTBUF_TAB, X, NUMNOD, NRTM, NSN, CAND_SEC_SEG, &
     &  JLT, CAND_N_CUR, CAND_E_CUR, IRECT, CONT_ELEMENT, COUNT, &
     &  IGRSURF, CAND_SEC_SEG_ID, CAND_MST_SEG_ID, &
     &  CAND_SEC_GP_MASK, &
     &  MAX_STS_SIZE_ACTUAL, NSURF_LOCAL, SEC_SURF_ID, MST_SURF_ID)
!-----------------------------------------------
!   M o d u l e s
!----------------------------------------------- 
      USE INTBUFDEF_MOD
      USE GROUPDEF_MOD
!-----------------------------------------------
!   M o d u l e s   /   I m p l i c i t   T y p e s
!-----------------------------------------------
      use constant_mod
      implicit none
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
#include      "mvsiz_p.inc"
#include      "my_real.inc"
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      TYPE(INTBUF_STRUCT_) INTBUF_TAB
      TYPE (SURF_)   , DIMENSION(NSURF_LOCAL)   :: IGRSURF
      INTEGER JLT, NUMNOD, NRTM, NSN, CAND_N_CUR(JLT), CAND_E_CUR(JLT)
      INTEGER IRECT(4,NRTM)
      my_real X(3,NUMNOD)
      INTEGER CAND_SEC_SEG(MAX_STS_SIZE_ACTUAL)
      INTEGER CAND_MST_SEG(MAX_STS_SIZE_ACTUAL)
      INTEGER CAND_SEC_SEG_ID(MAX_STS_SIZE_ACTUAL,5)
      INTEGER CAND_MST_SEG_ID(MAX_STS_SIZE_ACTUAL,5)
      INTEGER CAND_SEC_GP_MASK(MAX_STS_SIZE_ACTUAL,4)
      my_real CONT_ELEMENT(MAX_STS_SIZE_ACTUAL,3,8)
      INTEGER COUNT, MAX_STS_SIZE_ACTUAL
      INTEGER NSURF_LOCAL, SEC_SURF_ID, MST_SURF_ID
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER I, J, K, NI
      INTEGER candidate, candidateM, sec_seg
      INTEGER SEC_SURF_IDX, NSEG, NSEC_BOUNDS
      INTEGER pair_index
      LOGICAL :: found_corner, pair_added, capacity_full
      INTEGER, ALLOCATABLE :: IGRSURF_S_TEMP(:,:)
!-----------------------------------------------
!   S o u r c e   L i n e s
!-----------------------------------------------
      SEC_SURF_IDX = SEC_SURF_ID
      NSEC_BOUNDS = NSN
      IF (INTBUF_TAB%S_NSV > 0) THEN
        NSEC_BOUNDS = INTBUF_TAB%S_NSV
      ENDIF
      IF (JLT <= 0) THEN
        COUNT = 0
        RETURN
      END IF
      IF (SEC_SURF_IDX <= 0 .OR. SEC_SURF_IDX > NSURF_LOCAL) THEN
        COUNT = 0
        RETURN
      END IF
      IF (MST_SURF_ID <= 0 .OR. MST_SURF_ID > NSURF_LOCAL) THEN
        COUNT = 0
        RETURN
      END IF
      ! Safety checks
      IF (IGRSURF(SEC_SURF_IDX)%NSEG <= 0) THEN
        COUNT = 0
        RETURN
      END IF
      IF (.NOT. ALLOCATED(IGRSURF(SEC_SURF_IDX)%NODES)) THEN
        COUNT = 0
        RETURN
      END IF

      NSEG = IGRSURF(SEC_SURF_IDX)%NSEG

      ALLOCATE(IGRSURF_S_TEMP(IGRSURF(SEC_SURF_IDX)%NSEG, 4))
      IGRSURF_S_TEMP = IGRSURF(SEC_SURF_IDX)%NODES

      CAND_SEC_SEG = -1
      CAND_MST_SEG = -1
      CAND_SEC_GP_MASK = 0
      ! Compare candidate array with IGRSURF_TEMP table to derive the respective index
      COUNT = 0

      ! Map candidate nodes to segment pairs
      DO I = 1, JLT
        candidateM = CAND_E_CUR(I)
        IF (candidateM <= 0 .OR. candidateM > NRTM) CYCLE

        found_corner = .FALSE.
        capacity_full = .FALSE.
        IF (CAND_N_CUR(I) > 0 .AND. CAND_N_CUR(I) <= NSEC_BOUNDS) THEN
          candidate = INTBUF_TAB%NSV(CAND_N_CUR(I))
          IF (candidate > 0 .AND. candidate <= NUMNOD) THEN
            sec_seg = STS_REMAP_BEST_SEC_SEG_FOR_NODE(IGRSURF_S_TEMP, &
     &        NSEG, IRECT, NRTM, candidateM, candidate, NUMNOD, X)
            IF (sec_seg > 0) THEN
              found_corner = .TRUE.
              CALL STS_REMAP_TRY_ADD_PAIR(sec_seg, candidateM, COUNT, &
     &          CAND_SEC_SEG, CAND_MST_SEG, MAX_STS_SIZE_ACTUAL, &
     &          pair_added, pair_index)
              IF (pair_added .AND. pair_index > 0) THEN
                CAND_SEC_GP_MASK(pair_index, 1:4) = 1
              ELSE
                capacity_full = .TRUE.
              END IF
            ENDIF
          ENDIF
        ENDIF
        IF (capacity_full) EXIT

!       INT7 node candidates may miss /SURF corner nodes after separation.
!       I7TRC can also invalidate CAND_N to NSN+1 while CAND_E remains a
!       useful master-segment seed. Fall back to the secondary segment whose
!       centroid is closest to the INT7 main-segment centroid.
        IF (.NOT. found_corner) THEN
          sec_seg = STS_REMAP_NEAREST_SEC_SEG(IGRSURF_S_TEMP, NSEG, &
     &        IRECT, NRTM, candidateM, NUMNOD, X)
          IF (sec_seg > 0) THEN
            CALL STS_REMAP_TRY_ADD_PAIR(sec_seg, candidateM, COUNT, &
     &        CAND_SEC_SEG, CAND_MST_SEG, MAX_STS_SIZE_ACTUAL, &
     &        pair_added, pair_index)
            IF (pair_added .AND. pair_index > 0) THEN
              CAND_SEC_GP_MASK(pair_index, 1:4) = 1
            END IF
            IF (.NOT. pair_added) EXIT
          END IF
        END IF
      END DO

      IF (COUNT <= 0) THEN
        DEALLOCATE(IGRSURF_S_TEMP)
        RETURN
      END IF

      DO I = 1, COUNT
        CAND_SEC_SEG_ID(I, 1) = CAND_SEC_SEG(I)
        CAND_SEC_SEG_ID(I, 2:5) = IGRSURF(SEC_SURF_IDX)%NODES(CAND_SEC_SEG(I), 1:4)
        CAND_MST_SEG_ID(I, 1) = CAND_MST_SEG(I)
        CAND_MST_SEG_ID(I, 2:5) = IRECT(1:4, CAND_MST_SEG(I))
      END DO

      ! Store coordinates: Primary (1-4), Secondary (5-8)
      ! PRIMARY -> FIRST (1-4)
      DO I = 1, COUNT
        J = 1
        DO K = 2, 5
          NI = CAND_MST_SEG_ID(I, K)
          CONT_ELEMENT(I, 1, J) = X(1, NI)  ! X
          CONT_ELEMENT(I, 2, J) = X(2, NI)  ! Y
          CONT_ELEMENT(I, 3, J) = X(3, NI)  ! Z
          J = J + 1
        END DO
      END DO

      ! SECONDARY -> SECOND (5-8)
      DO I = 1, COUNT
        J = 5
        DO K = 2, 5
          NI = CAND_SEC_SEG_ID(I, K)
          CONT_ELEMENT(I, 1, J) = X(1, NI)  ! X
          CONT_ELEMENT(I, 2, J) = X(2, NI)  ! Y
          CONT_ELEMENT(I, 3, J) = X(3, NI)  ! Z
          J = J + 1
        END DO
      END DO

      DEALLOCATE(IGRSURF_S_TEMP)

      RETURN
      CONTAINS

      !=======================================================================
      ! STS_REMAP_TRY_ADD_PAIR
      !
      ! Try to add a segment pair to the candidate list.
      !=======================================================================
        SUBROUTINE STS_REMAP_TRY_ADD_PAIR(SEC_SEG_IN, MST_SEG_IN, &
     &    COUNT_INOUT, CAND_SEC, CAND_MST, CAPACITY, PAIR_ADDED, &
     &    PAIR_INDEX)
          INTEGER, INTENT(IN) :: SEC_SEG_IN, MST_SEG_IN, CAPACITY
          INTEGER, INTENT(INOUT) :: COUNT_INOUT
          INTEGER, INTENT(INOUT) :: CAND_SEC(CAPACITY)
          INTEGER, INTENT(INOUT) :: CAND_MST(CAPACITY)
          LOGICAL, INTENT(OUT) :: PAIR_ADDED
          INTEGER, INTENT(OUT) :: PAIR_INDEX
          INTEGER :: K
          LOGICAL :: duplicate

          PAIR_ADDED = .FALSE.
          PAIR_INDEX = 0
          IF (SEC_SEG_IN <= 0 .OR. MST_SEG_IN <= 0) RETURN

          duplicate = .FALSE.
          DO K = 1, COUNT_INOUT
            IF (CAND_SEC(K) == SEC_SEG_IN .AND. &
     &          CAND_MST(K) == MST_SEG_IN) THEN
              duplicate = .TRUE.
              PAIR_INDEX = K
              EXIT
            END IF
          END DO
          IF (duplicate) THEN
            PAIR_ADDED = .TRUE.
            RETURN
          END IF

          COUNT_INOUT = COUNT_INOUT + 1
          IF (COUNT_INOUT > CAPACITY) THEN
            COUNT_INOUT = COUNT_INOUT - 1
            PAIR_ADDED = .FALSE.
            RETURN
          END IF

          CAND_SEC(COUNT_INOUT) = SEC_SEG_IN
          CAND_MST(COUNT_INOUT) = MST_SEG_IN
          PAIR_INDEX = COUNT_INOUT
          PAIR_ADDED = .TRUE.
        END SUBROUTINE STS_REMAP_TRY_ADD_PAIR

        !=======================================================================
        ! STS_REMAP_BEST_SEC_SEG_FOR_NODE
        !
        ! Find the secondary segment that is closest to the main segment centroid.
        !=======================================================================
        INTEGER FUNCTION STS_REMAP_BEST_SEC_SEG_FOR_NODE(IGRSURF_NODES, &
     &    NSEG_IN, IRECT_IN, NRTM_IN, MST_SEG_IN, SEC_NODE_IN, &
     &    NUMNOD_IN, X_IN)
          INTEGER, INTENT(IN) :: NSEG_IN, NRTM_IN, MST_SEG_IN, SEC_NODE_IN
          INTEGER, INTENT(IN) :: NUMNOD_IN
          INTEGER, INTENT(IN) :: IGRSURF_NODES(NSEG_IN, 4)
          INTEGER, INTENT(IN) :: IRECT_IN(4, NRTM_IN)
          my_real, INTENT(IN) :: X_IN(3, NUMNOD_IN)
          INTEGER :: J, K, NID, NC
          my_real :: XM(3), XS(3), DIST2, BEST

          STS_REMAP_BEST_SEC_SEG_FOR_NODE = 0
          IF (MST_SEG_IN <= 0 .OR. MST_SEG_IN > NRTM_IN) RETURN
          IF (SEC_NODE_IN <= 0 .OR. SEC_NODE_IN > NUMNOD_IN) RETURN
          IF (NSEG_IN <= 0) RETURN

          XM = ZERO
          NC = 0
          DO K = 1, 4
            NID = IRECT_IN(K, MST_SEG_IN)
            IF (NID <= 0 .OR. NID > NUMNOD_IN) CYCLE
            NC = NC + 1
            XM(1) = XM(1) + X_IN(1, NID)
            XM(2) = XM(2) + X_IN(2, NID)
            XM(3) = XM(3) + X_IN(3, NID)
          END DO
          IF (NC <= 0) RETURN
          XM(1) = XM(1) / NC
          XM(2) = XM(2) / NC
          XM(3) = XM(3) / NC

          BEST = HUGE(1.0D0)
          DO J = 1, NSEG_IN
            IF (.NOT. ANY(SEC_NODE_IN == IGRSURF_NODES(J, 1:4))) CYCLE
            XS = ZERO
            NC = 0
            DO K = 1, 4
              NID = IGRSURF_NODES(J, K)
              IF (NID <= 0 .OR. NID > NUMNOD_IN) CYCLE
              NC = NC + 1
              XS(1) = XS(1) + X_IN(1, NID)
              XS(2) = XS(2) + X_IN(2, NID)
              XS(3) = XS(3) + X_IN(3, NID)
            END DO
            IF (NC <= 0) CYCLE
            XS(1) = XS(1) / NC
            XS(2) = XS(2) / NC
            XS(3) = XS(3) / NC
            DIST2 = (XS(1) - XM(1))**2 + (XS(2) - XM(2))**2 + &
     &              (XS(3) - XM(3))**2
            IF (DIST2 < BEST) THEN
              BEST = DIST2
              STS_REMAP_BEST_SEC_SEG_FOR_NODE = J
            END IF
          END DO
        END FUNCTION STS_REMAP_BEST_SEC_SEG_FOR_NODE

        !=======================================================================
        ! STS_REMAP_NEAREST_SEC_SEG
        !
        ! Find the secondary segment that is closest to the main segment centroid.
        !=======================================================================
        INTEGER FUNCTION STS_REMAP_NEAREST_SEC_SEG(IGRSURF_NODES, &
     &    NSEG_IN, IRECT_IN, NRTM_IN, MST_SEG_IN, NUMNOD_IN, X_IN)
          INTEGER, INTENT(IN) :: NSEG_IN, NRTM_IN, MST_SEG_IN, NUMNOD_IN
          INTEGER, INTENT(IN) :: IGRSURF_NODES(NSEG_IN, 4)
          INTEGER, INTENT(IN) :: IRECT_IN(4, NRTM_IN)
          my_real, INTENT(IN) :: X_IN(3, NUMNOD_IN)
          INTEGER :: J, K, NID, NC
          my_real :: XM(3), XS(3), DIST2, BEST

          STS_REMAP_NEAREST_SEC_SEG = 0
          IF (MST_SEG_IN <= 0 .OR. MST_SEG_IN > NRTM_IN) RETURN
          IF (NSEG_IN <= 0) RETURN

          XM = ZERO
          NC = 0
          DO K = 1, 4
            NID = IRECT_IN(K, MST_SEG_IN)
            IF (NID <= 0 .OR. NID > NUMNOD_IN) CYCLE
            NC = NC + 1
            XM(1) = XM(1) + X_IN(1, NID)
            XM(2) = XM(2) + X_IN(2, NID)
            XM(3) = XM(3) + X_IN(3, NID)
          END DO
          IF (NC <= 0) RETURN
          XM(1) = XM(1) / NC
          XM(2) = XM(2) / NC
          XM(3) = XM(3) / NC

          BEST = HUGE(1.0D0)
          DO J = 1, NSEG_IN
            XS = ZERO
            NC = 0
            DO K = 1, 4
              NID = IGRSURF_NODES(J, K)
              IF (NID <= 0 .OR. NID > NUMNOD_IN) CYCLE
              NC = NC + 1
              XS(1) = XS(1) + X_IN(1, NID)
              XS(2) = XS(2) + X_IN(2, NID)
              XS(3) = XS(3) + X_IN(3, NID)
            END DO
            IF (NC <= 0) CYCLE
            XS(1) = XS(1) / NC
            XS(2) = XS(2) / NC
            XS(3) = XS(3) / NC
            DIST2 = (XS(1) - XM(1))**2 + (XS(2) - XM(2))**2 + &
     &              (XS(3) - XM(3))**2
            IF (DIST2 < BEST) THEN
              BEST = DIST2
              STS_REMAP_NEAREST_SEC_SEG = J
            END IF
          END DO
        END FUNCTION STS_REMAP_NEAREST_SEC_SEG

      END SUBROUTINE STS_REMAP_SEGMENTS
