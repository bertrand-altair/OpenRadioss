!||====================================================================
!||    STS_CONTACT_STIFFNESS  ../engine/source/interfaces/ists/ists_contact_stiffness.F90
!||--- called by ------------------------------------------------------
!||    i7mainf                ../engine/source/interfaces/int07/i7mainf.F
!||====================================================================
      MODULE STS_CONTACT_STIFFNESS_MOD
        USE PRECISION_MOD, ONLY : WP
        IMPLICIT NONE
        PRIVATE

        PUBLIC :: STS_CONTACT_STIFFNESS

      CONTAINS

!=======================================================================
!   STS_CONTACT_STIFFNESS
!
!   Stiffness calculation for STS.
!=======================================================================
        SUBROUTINE STS_CONTACT_STIFFNESS( &
     &      CAND_MST_SEG_ID, CAND_SEC_SEG_ID, COUNT, MAX_STS_SIZE, &
     &      IRECTM, STFM, NRTM, NSV, STFNS, NSN, NUMNOD, &
     &      IGSTI, KMIN, KMAX, STIGLO, STIF)
          INTEGER, INTENT(IN) :: COUNT, MAX_STS_SIZE
          INTEGER, INTENT(IN) :: NRTM, NSN, NUMNOD, IGSTI
          INTEGER, INTENT(IN) :: CAND_MST_SEG_ID(MAX_STS_SIZE,5)
          INTEGER, INTENT(IN) :: CAND_SEC_SEG_ID(MAX_STS_SIZE,5)
          INTEGER, INTENT(IN) :: IRECTM(:)
          INTEGER, INTENT(IN) :: NSV(:)
          REAL(KIND=WP), INTENT(IN) :: STFM(:)
          REAL(KIND=WP), INTENT(IN) :: STFNS(:)
          REAL(KIND=WP), INTENT(IN) :: KMIN, KMAX, STIGLO
          REAL(KIND=WP), INTENT(OUT) :: STIF(MAX_STS_SIZE)

          INTEGER :: I, J, IDX, MST_SEG, NODE_ID, NVALID
          INTEGER :: NSN_EFF, NRTM_EFF
          INTEGER, ALLOCATABLE :: NODE_TO_STFNS(:)
          REAL(KIND=WP) :: K_PRIMARY, K_SECONDARY, K_NODE
          REAL(KIND=WP) :: SECONDARY_FALLBACK
          REAL(KIND=WP), PARAMETER :: EPS_STIFF = 1.0E-30_WP

          STIF(1:MAX_STS_SIZE) = 0.0_WP

          ! get the effective number of nodes and segments
          NSN_EFF = MIN(NSN, SIZE(NSV), SIZE(STFNS))
          NRTM_EFF = MIN(NRTM, SIZE(STFM), SIZE(IRECTM) / 4)
          SECONDARY_FALLBACK = STS_CONTACT_AVERAGE_POSITIVE( &
     &      STFNS, NSN_EFF, ABS(STIGLO))

          ALLOCATE(NODE_TO_STFNS(MAX(1, NUMNOD)))
          NODE_TO_STFNS(:) = 0
          DO I = 1, NSN_EFF
            NODE_ID = NSV(I)
            IF (NODE_ID > 0 .AND. NODE_ID <= NUMNOD) THEN
              NODE_TO_STFNS(NODE_ID) = I
            END IF
          END DO

          DO I = 1, MIN(COUNT, MAX_STS_SIZE)
            ! find the primary segment
            ! IRECTM is prepared in starter and passed here as IRECTM/IRECT (INTBUF_TAB%IRECTM).
            MST_SEG = STS_CONTACT_FIND_PRIMARY_SEG( &
     &        CAND_MST_SEG_ID(I, 2:5), CAND_MST_SEG_ID(I, 1), &
     &        IRECTM, NRTM_EFF)

            ! get the primary segment stiffness
            K_PRIMARY = 0.0_WP
            IF (MST_SEG > 0 .AND. MST_SEG <= NRTM_EFF) THEN
              K_PRIMARY = STFM(MST_SEG)
            END IF

            K_SECONDARY = 0.0_WP
            ! get the secondary stiffness
            ! NVALID is the number of valid secondary nodes
            NVALID = 0
            DO J = 2, 5
              NODE_ID = CAND_SEC_SEG_ID(I, J)
              IF (NODE_ID <= 0 .OR. NODE_ID > NUMNOD) CYCLE

              IDX = NODE_TO_STFNS(NODE_ID)
              IF (IDX <= 0 .OR. IDX > NSN_EFF .OR. IDX > SIZE(STFNS)) CYCLE

              K_NODE = ABS(STFNS(IDX))
              IF (K_NODE <= EPS_STIFF) CYCLE

              ! Sum only valid mapped secondary stiffness values.
              K_SECONDARY = K_SECONDARY + K_NODE
              NVALID = NVALID + 1
            END DO
            IF (NVALID > 0) THEN
              ! Average over valid mapped nodes.
              K_SECONDARY = K_SECONDARY / REAL(NVALID, WP)
            ELSE
              ! Deterministic fallback when no secondary node maps to NSV.
              K_SECONDARY = SECONDARY_FALLBACK
            END IF

            ! Store the stiffness value in the STS_STIF array
            STIF(I) = STS_CONTACT_STIFFNESS_VALUE( &
     &        K_PRIMARY, K_SECONDARY, IGSTI, KMIN, KMAX)

          END DO

          DEALLOCATE(NODE_TO_STFNS)
        END SUBROUTINE STS_CONTACT_STIFFNESS

!=======================================================================
!   STS_CONTACT_FIND_PRIMARY_SEG
!=======================================================================
        INTEGER FUNCTION STS_CONTACT_FIND_PRIMARY_SEG( &
     &      PRIMARY_NODES, CAND_SEG, IRECTM, NRTM)
          INTEGER, INTENT(IN) :: PRIMARY_NODES(4)
          INTEGER, INTENT(IN) :: CAND_SEG, NRTM
          INTEGER, INTENT(IN) :: IRECTM(:)
          INTEGER :: SEG

          STS_CONTACT_FIND_PRIMARY_SEG = 0

          IF (CAND_SEG > 0 .AND. CAND_SEG <= NRTM) THEN
            ! check if the candidate segment CAND_SEG has the same nodes as the primary segment PRIMARY_NODES
            IF (STS_CONTACT_SAME_SEG_NODES(PRIMARY_NODES, IRECTM, CAND_SEG)) THEN
              STS_CONTACT_FIND_PRIMARY_SEG = CAND_SEG
              RETURN
            END IF
          END IF

          DO SEG = 1, NRTM
            IF (SEG == CAND_SEG) CYCLE
            ! check if the other segment SEG has the same nodes as the primary segment PRIMARY_NODES
            IF (STS_CONTACT_SAME_SEG_NODES(PRIMARY_NODES, IRECTM, SEG)) THEN
              STS_CONTACT_FIND_PRIMARY_SEG = SEG
              RETURN
            END IF
          END DO
          ! return 0 if no primary segment is found
          ! return the primary segment if found
        END FUNCTION STS_CONTACT_FIND_PRIMARY_SEG

!=======================================================================
!   STS_CONTACT_SAME_SEG_NODES
!=======================================================================
        LOGICAL FUNCTION STS_CONTACT_SAME_SEG_NODES(PRIMARY_NODES, IRECTM, SEG)
          INTEGER, INTENT(IN) :: PRIMARY_NODES(4)
          INTEGER, INTENT(IN) :: IRECTM(:)
          INTEGER, INTENT(IN) :: SEG
          INTEGER :: I, J, NODE_ID
          LOGICAL :: FOUND

          STS_CONTACT_SAME_SEG_NODES = .FALSE.

          DO I = 1, 4
            NODE_ID = PRIMARY_NODES(I)
            IF (NODE_ID <= 0) CYCLE
            FOUND = .FALSE.
            DO J = 1, 4
              ! check if the node NODE_ID is in the segment SEG
              IF (NODE_ID == IRECTM(4 * (SEG - 1) + J)) THEN
                FOUND = .TRUE.
                EXIT
              END IF
            END DO
            IF (.NOT. FOUND) RETURN
          END DO

          STS_CONTACT_SAME_SEG_NODES = .TRUE.
          ! return TRUE if the segment SEG has the same nodes as the primary segment PRIMARY_NODES
        END FUNCTION STS_CONTACT_SAME_SEG_NODES

!=======================================================================
!   STS_CONTACT_STIFFNESS_VALUE
!=======================================================================
        REAL(KIND=WP) FUNCTION STS_CONTACT_STIFFNESS_VALUE(K_PRIMARY, K_SECONDARY, IGSTI, KMIN, KMAX)
          REAL(KIND=WP), INTENT(IN) :: K_PRIMARY, K_SECONDARY
          REAL(KIND=WP), INTENT(IN) :: KMIN, KMAX
          INTEGER, INTENT(IN) :: IGSTI
          REAL(KIND=WP), PARAMETER :: EPS_DEN = 1.0E-30_WP

          ! IGSTI is the interface stiffness selection flag (how primary and secondary are combined).
          ! IGSTI = 1: product
          ! IGSTI = 2: average
          ! IGSTI = 3: maximum
          ! IGSTI = 4: minimum
          ! IGSTI = 5: harmonic-like blend

          IF (IGSTI <= 1) THEN
            STS_CONTACT_STIFFNESS_VALUE = K_PRIMARY * ABS(K_SECONDARY)
          ELSEIF (IGSTI == 2) THEN
            STS_CONTACT_STIFFNESS_VALUE = 0.5_WP * (K_PRIMARY + ABS(K_SECONDARY))
            STS_CONTACT_STIFFNESS_VALUE = MAX(KMIN, MIN(STS_CONTACT_STIFFNESS_VALUE, KMAX))
          ELSEIF (IGSTI == 3) THEN
            STS_CONTACT_STIFFNESS_VALUE = MAX(K_PRIMARY, ABS(K_SECONDARY))
            STS_CONTACT_STIFFNESS_VALUE = MAX(KMIN, MIN(STS_CONTACT_STIFFNESS_VALUE, KMAX))
          ELSEIF (IGSTI == 4) THEN
            STS_CONTACT_STIFFNESS_VALUE = MIN(K_PRIMARY, ABS(K_SECONDARY))
            STS_CONTACT_STIFFNESS_VALUE = MAX(KMIN, MIN(STS_CONTACT_STIFFNESS_VALUE, KMAX))
          ELSEIF (IGSTI == 5) THEN
            STS_CONTACT_STIFFNESS_VALUE = K_PRIMARY * ABS(K_SECONDARY) / MAX(EPS_DEN, K_PRIMARY + ABS(K_SECONDARY))
            STS_CONTACT_STIFFNESS_VALUE = MAX(KMIN, MIN(STS_CONTACT_STIFFNESS_VALUE, KMAX))
          ELSE
            STS_CONTACT_STIFFNESS_VALUE = K_PRIMARY * ABS(K_SECONDARY)
          END IF
        END FUNCTION STS_CONTACT_STIFFNESS_VALUE

!=======================================================================
!   STS_CONTACT_AVERAGE_POSITIVE
!=======================================================================
        REAL(KIND=WP) FUNCTION STS_CONTACT_AVERAGE_POSITIVE( &
     &      VALUES, N_EFF, FALLBACK)
          REAL(KIND=WP), INTENT(IN) :: VALUES(:)
          INTEGER, INTENT(IN) :: N_EFF
          REAL(KIND=WP), INTENT(IN) :: FALLBACK

          INTEGER :: I, NLOC, COUNT_POS
          REAL(KIND=WP) :: SUM_POS, VABS
          REAL(KIND=WP), PARAMETER :: EPS_STIFF = 1.0E-30_WP

          STS_CONTACT_AVERAGE_POSITIVE = FALLBACK
          NLOC = MIN(N_EFF, SIZE(VALUES))
          IF (NLOC <= 0) RETURN

          SUM_POS = 0.0_WP
          COUNT_POS = 0
          DO I = 1, NLOC
            VABS = ABS(VALUES(I))
            IF (VABS <= EPS_STIFF) CYCLE
            SUM_POS = SUM_POS + VABS
            COUNT_POS = COUNT_POS + 1
          END DO

          IF (COUNT_POS > 0) THEN
            STS_CONTACT_AVERAGE_POSITIVE = SUM_POS / REAL(COUNT_POS, WP)
          END IF
        END FUNCTION STS_CONTACT_AVERAGE_POSITIVE

      END MODULE STS_CONTACT_STIFFNESS_MOD
