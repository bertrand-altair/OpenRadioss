!||====================================================================
!||    STS_CONTACTS_ASSEMBLE  ../engine/source/interfaces/ists/ists_contacts_assemble.F90
!||--- called by ------------------------------------------------------
!||    i7mainf              ../engine/source/interfaces/int07/i7mainf.F
!||--- calls ---------------------------------------------------------
!||    STS_CONTACT_EVAL_PAIR    ../engine/source/interfaces/ists/ists_CONTACT_EVAL_PAIR.F90
!||====================================================================
      SUBROUTINE STS_CONTACTS_ASSEMBLE(CONT_ELEMENT, COUNT, OPTION, STS_INTERFACE_ID, NCYCLE_IN, TIME_CUR, &
     & CAND_MST_SEG_ID, CAND_SEC_SEG_ID, CAND_SEC_GP_MASK, &
     & load_arr, node_id_load, L_out, IMPACT_glob, STIF, &
     & MAX_STS_SIZE_ACTUAL, FRICC, XMU, IFPEN, QFRICT, GAP, V, MS, &
     & VISC, IVIS2, VISCFFRIC, DT2T, NELTST, ITYPTST, &
     & ECONTT_TOT, ECONVT_TOT, FN_TOT, FT_TOT)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE INTBUFDEF_MOD
!-----------------------------------------------
!   M o d u l e s   /   I m p l i c i t   T y p e s
!-----------------------------------------------
      use constant_mod
      use sts_gp_state_mod
      implicit none
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
#include      "mvsiz_p.inc"
#include      "my_real.inc"
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!      TYPE(INTBUF_STRUCT_) INTBUF_TAB(*)
      REAL*8 CONT_ELEMENT(MAX_STS_SIZE_ACTUAL,3,8)
      my_real STIF(MAX_STS_SIZE_ACTUAL)
      INTEGER COUNT, OPTION, STS_INTERFACE_ID, NCYCLE_IN
      REAL*8 TIME_CUR
      INTEGER CAND_SEC_SEG_ID(MAX_STS_SIZE_ACTUAL,5)
      INTEGER CAND_MST_SEG_ID(MAX_STS_SIZE_ACTUAL,5)
      INTEGER CAND_SEC_GP_MASK(MAX_STS_SIZE_ACTUAL,4)
      REAL*8 load_arr(MAX_STS_SIZE_ACTUAL,8,4)
      INTEGER node_id_load(MAX_STS_SIZE_ACTUAL*8)
      INTEGER L_out, IMPACT_glob, MAX_STS_SIZE_ACTUAL
      my_real FRICC(MVSIZ)
      my_real XMU(MVSIZ)
      INTEGER IFPEN(MAX_STS_SIZE_ACTUAL)     
      my_real QFRICT
      my_real GAP  ! Gap value from user input
      my_real V(3,*), MS(*), VISC, VISCFFRIC(MVSIZ), DT2T
      INTEGER IVIS2, NELTST, ITYPTST
      REAL*8 ECONTT_TOT, ECONVT_TOT
      REAL*8 FN_TOT(3), FT_TOT(3)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER I, J, K, L, IMPACT
      INTEGER selected_option, impact_gauss, impact_lobatto
      INTEGER neltst_probe, ityptst_probe
      INTEGER valid_gauss, valid_lobatto
      INTEGER LUX_STS
      REAL*8 XUPD(3,8)
      REAL*8 p_load_new(24)
      REAL*8 p_probe(24)
      REAL*8 node_stiff(8)
      REAL*8 node_stiff_probe(8)
      REAL*8 unit_gp_weight(4)
      REAL*8 p_friction(24)  ! Friction forces (separate output)
      REAL*8 p_friction_probe(24)
      REAL*8 pair_max_penetration
      REAL*8 probe_pen_gauss, probe_pen_lobatto
      REAL*8 probe_score_gauss, probe_score_lobatto
      REAL*8 min_pene_gauss, min_pene_lobatto
      REAL*8 econt_pair, econtv_pair
      REAL*8 econt_probe, econtv_probe
      my_real EFRICT_LOC
      my_real QFRICT_PROBE, DT2T_PROBE
      INTEGER node_ids(8)  ! Node IDs for velocity interpolation
      REAL*8 fx_prim, fy_prim, fz_prim, fx_sec, fy_sec, fz_sec
      REAL*8 fxf_prim, fyf_prim, fzf_prim, fxf_sec, fyf_sec, fzf_sec
      REAL*8 gap_abs, lobatto_margin
      REAL*8, ALLOCATABLE :: lobatto_gp_weight(:,:)
      LOGICAL FILE_EXISTS_STS, STS_CSV_INITIALIZED
      LOGICAL, PARAMETER :: CSV_OUTPUT_ENABLED = .FALSE.
      REAL*8, PARAMETER :: STS_MIXED_LOBATTO_GAP_MARGIN = 2.0D-2
      REAL*8, PARAMETER :: STS_MIXED_LOBATTO_REL_MARGIN = 2.0D-1
      SAVE STS_CSV_INITIALIZED
      DATA STS_CSV_INITIALIZED /.FALSE./
!-----------------------------------------------
!   I n i t i a l i z a t i o n
!-----------------------------------------------
      IMPACT_glob = 0
      ECONTT_TOT = 0.0D0
      ECONVT_TOT = 0.0D0
      FN_TOT = 0.0D0
      FT_TOT = 0.0D0
      
      ! Safety check
      IF (COUNT <= 0) THEN
        L_out = 1
        RETURN
      END IF
      
      ! Initialize counters
      K = 1
      L = 1
      unit_gp_weight = 1.0D0

      ALLOCATE(lobatto_gp_weight(4, COUNT))
      lobatto_gp_weight = 1.0D0
      IF (OPTION == 1 .OR. OPTION == 2) THEN
        lobatto_gp_weight = 0.0D0
        CALL STS_BUILD_LOBATTO_GP_WEIGHTS(COUNT, MAX_STS_SIZE_ACTUAL, &
     &    CAND_SEC_SEG_ID, CAND_MST_SEG_ID, CAND_SEC_GP_MASK, &
     &    lobatto_gp_weight)
      ENDIF

      IF (CSV_OUTPUT_ENABLED) THEN
        IF (.NOT. STS_CSV_INITIALIZED) THEN
          INQUIRE(FILE='sts_contact_forces.csv', EXIST=FILE_EXISTS_STS)
          IF (FILE_EXISTS_STS) THEN
            OPEN(NEWUNIT=LUX_STS, FILE='sts_contact_forces.csv', &
     &           STATUS='OLD', ACTION='WRITE', POSITION='APPEND')
          ELSE
            OPEN(NEWUNIT=LUX_STS, FILE='sts_contact_forces.csv', &
     &           STATUS='NEW', ACTION='WRITE')
            WRITE(LUX_STS,'(A)') &
     &        'cycle,time,interface_id,entity_id,fx,fy,fz,force_norm,fn,ft,n_pairs,max_penetration'
          ENDIF
          STS_CSV_INITIALIZED = .TRUE.
        ELSE
          OPEN(NEWUNIT=LUX_STS, FILE='sts_contact_forces.csv', &
     &         STATUS='OLD', ACTION='WRITE', POSITION='APPEND')
        ENDIF
      END IF
!-----------------------------------------------
!   M a i n   L o o p
!-----------------------------------------------
      DO I = 1, COUNT
        IMPACT = 0
        XUPD = CONT_ELEMENT(I, 1:3, 1:8)
      
        ! Get node IDs for velocity interpolation
        DO J = 1, 4
          node_ids(J)   = CAND_MST_SEG_ID(I, J+1)   ! Primary nodes
          node_ids(J+4) = CAND_SEC_SEG_ID(I, J+1)   ! Secondary nodes
        ENDDO
        XMU(1) = FRICC(MIN(I,MVSIZ)) ! Friction coefficient mu
      
        selected_option = OPTION
        impact_gauss = 0
        impact_lobatto = 0
        probe_pen_gauss = 0.0D0
        probe_pen_lobatto = 0.0D0
        probe_score_gauss = 0.0D0
        probe_score_lobatto = 0.0D0
        valid_gauss = 0
        valid_lobatto = 0
        min_pene_gauss = 0.0D0
        min_pene_lobatto = 0.0D0

        IF (OPTION == 2 .OR. OPTION == 0) THEN
          DT2T_PROBE = DT2T
          NELTST_PROBE = NELTST
          ITYPTST_PROBE = ITYPTST
          QFRICT_PROBE = 0.d0

          CALL STS_CONTACT_EVAL_PAIR(XUPD, STIF(I), p_probe, &
     &                      impact_gauss, I, node_stiff_probe, 0, &
     &                      FRICC, XMU, IFPEN, &
     &                      p_friction_probe, EFRICT_LOC, QFRICT_PROBE, &
     &                      node_ids, V, .FALSE., MAX_STS_SIZE_ACTUAL, &
     &                      GAP, unit_gp_weight, probe_pen_gauss, &
     &                      econt_probe, &
     &                      econtv_probe, MS, STS_INTERFACE_ID, VISC, &
     &                      IVIS2, VISCFFRIC(MIN(I,MVSIZ)), DT2T_PROBE, &
     &                      NELTST_PROBE, ITYPTST_PROBE, &
     &                      .FALSE., probe_score_gauss, valid_gauss, &
     &                      min_pene_gauss)

          DT2T_PROBE = DT2T
          NELTST_PROBE = NELTST
          ITYPTST_PROBE = ITYPTST
          QFRICT_PROBE = 0.d0

          CALL STS_CONTACT_EVAL_PAIR(XUPD, STIF(I), p_probe, &
     &                      impact_lobatto, I, node_stiff_probe, 1, &
     &                      FRICC, XMU, IFPEN, &
     &                      p_friction_probe, EFRICT_LOC, QFRICT_PROBE, &
     &                      node_ids, V, .FALSE., MAX_STS_SIZE_ACTUAL, &
     &                      GAP, lobatto_gp_weight(1:4,I), &
     &                      probe_pen_lobatto, &
     &                      econt_probe, &
     &                      econtv_probe, MS, STS_INTERFACE_ID, VISC, &
     &                      IVIS2, VISCFFRIC(MIN(I,MVSIZ)), DT2T_PROBE, &
     &                      NELTST_PROBE, ITYPTST_PROBE, &
     &                      .FALSE., probe_score_lobatto, &
     &                      valid_lobatto, min_pene_lobatto)

        ENDIF

        IF (OPTION == 2) THEN
          IF (impact_gauss == 0 .AND. impact_lobatto == 0) THEN
            selected_option = -1
          ELSE IF (impact_lobatto == 1 .AND. impact_gauss == 0) THEN
            selected_option = 1
          ELSE IF (impact_gauss == 1 .AND. impact_lobatto == 0) THEN
            gap_abs = MAX(DABS(DBLE(GAP)), 1.0D-30)
            lobatto_margin = MAX(STS_MIXED_LOBATTO_GAP_MARGIN*gap_abs, &
     &        STS_MIXED_LOBATTO_REL_MARGIN*probe_pen_gauss)
            IF (valid_lobatto <= 0) THEN
              selected_option = -1
            ELSE IF (min_pene_lobatto > lobatto_margin) THEN
              selected_option = -1
            ELSE
              selected_option = 0
            ENDIF
          ELSE
            gap_abs = MAX(DABS(DBLE(GAP)), 1.0D-30)
            lobatto_margin = MAX(STS_MIXED_LOBATTO_GAP_MARGIN*gap_abs, &
     &        STS_MIXED_LOBATTO_REL_MARGIN*MAX(probe_pen_gauss, &
     &        probe_pen_lobatto))
            IF (probe_pen_lobatto > probe_pen_gauss + &
     &          lobatto_margin) THEN
              selected_option = 1
            ELSE IF (probe_pen_lobatto > probe_pen_gauss .AND. &
     &          probe_score_lobatto > &
     &          (1.0D0 + STS_MIXED_LOBATTO_REL_MARGIN) * &
     &          probe_score_gauss) THEN
              selected_option = 1
            ELSE
              selected_option = 0
            ENDIF
          ENDIF
        ELSE IF (OPTION == 0) THEN
          IF (impact_gauss == 1 .AND. impact_lobatto == 0) THEN
            gap_abs = MAX(DABS(DBLE(GAP)), 1.0D-30)
            lobatto_margin = MAX(STS_MIXED_LOBATTO_GAP_MARGIN*gap_abs, &
     &        STS_MIXED_LOBATTO_REL_MARGIN*probe_pen_gauss)
            IF (valid_lobatto <= 0) THEN
              selected_option = -1
            ELSE IF (min_pene_lobatto > lobatto_margin) THEN
              selected_option = -1
            ENDIF
          ENDIF
        ENDIF

        IF (selected_option < 0) THEN
          CYCLE
        ENDIF

        ! Commit exactly one quadrature. Probe calls above are side-effect free.
        ! Normal penalty: d1 = 0.5*STIF*FAC; friction trial: d1_fric = 0.5*STIF (NTS STIF0).
        CALL STS_CONTACT_EVAL_PAIR(XUPD, STIF(I), p_load_new, IMPACT, I, &
     &                    node_stiff, selected_option, &
     &                    FRICC, XMU, IFPEN, &
     &                    p_friction, EFRICT_LOC, QFRICT, node_ids, V, &
     &                    .TRUE., MAX_STS_SIZE_ACTUAL, GAP, &
     &                    lobatto_gp_weight(1:4,I), &
     &                    pair_max_penetration, econt_pair, econtv_pair, &
     &                    MS, STS_INTERFACE_ID, VISC, IVIS2, &
     &                    VISCFFRIC(MIN(I,MVSIZ)), DT2T, NELTST, ITYPTST, &
     &                    .TRUE., probe_score_gauss, valid_gauss, &
     &                    min_pene_gauss)
      
        IF (IMPACT == 1) THEN
          IMPACT_glob = 1
          ECONTT_TOT = ECONTT_TOT + econt_pair
          ECONVT_TOT = ECONVT_TOT + econtv_pair

!         Slave-side force resultants for /TH/INTER (normal = p - p_friction).
          DO J = 5, 8
            FN_TOT(1) = FN_TOT(1) + p_load_new(3*(J-1)+1) - p_friction(3*(J-1)+1)
            FN_TOT(2) = FN_TOT(2) + p_load_new(3*(J-1)+2) - p_friction(3*(J-1)+2)
            FN_TOT(3) = FN_TOT(3) + p_load_new(3*(J-1)+3) - p_friction(3*(J-1)+3)
            FT_TOT(1) = FT_TOT(1) + p_friction(3*(J-1)+1)
            FT_TOT(2) = FT_TOT(2) + p_friction(3*(J-1)+2)
            FT_TOT(3) = FT_TOT(3) + p_friction(3*(J-1)+3)
          ENDDO

          IF (CSV_OUTPUT_ENABLED) THEN
!           Export two rows per contact pair:
!           - surface 1 (primary nodes 1..4)
!           - surface 2 (secondary nodes 5..8)
            fx_prim = 0.0D0
            fy_prim = 0.0D0
            fz_prim = 0.0D0
            fxf_prim = 0.0D0
            fyf_prim = 0.0D0
            fzf_prim = 0.0D0
            DO J = 1, 4
              fx_prim = fx_prim + p_load_new(3*(J-1)+1)
              fy_prim = fy_prim + p_load_new(3*(J-1)+2)
              fz_prim = fz_prim + p_load_new(3*(J-1)+3)
              fxf_prim = fxf_prim + p_friction(3*(J-1)+1)
              fyf_prim = fyf_prim + p_friction(3*(J-1)+2)
              fzf_prim = fzf_prim + p_friction(3*(J-1)+3)
            ENDDO

            fx_sec = 0.0D0
            fy_sec = 0.0D0
            fz_sec = 0.0D0
            fxf_sec = 0.0D0
            fyf_sec = 0.0D0
            fzf_sec = 0.0D0
            DO J = 5, 8
              fx_sec = fx_sec + p_load_new(3*(J-1)+1)
              fy_sec = fy_sec + p_load_new(3*(J-1)+2)
              fz_sec = fz_sec + p_load_new(3*(J-1)+3)
              fxf_sec = fxf_sec + p_friction(3*(J-1)+1)
              fyf_sec = fyf_sec + p_friction(3*(J-1)+2)
              fzf_sec = fzf_sec + p_friction(3*(J-1)+3)
            ENDDO

            CALL STS_CONTACT_EXPORT_CSV_PAIR(LUX_STS, NCYCLE_IN, TIME_CUR, STS_INTERFACE_ID, &
     &          CAND_MST_SEG_ID(I,1), CAND_SEC_SEG_ID(I,1), pair_max_penetration, &
     &          fx_prim, fy_prim, fz_prim, fxf_prim, fyf_prim, fzf_prim, &
     &          fx_sec, fy_sec, fz_sec, fxf_sec, fyf_sec, fzf_sec)
          END IF
      
          ! Save node IDs: Primary (1-4), Secondary (5-8)
          node_id_load(K:K+3) = CAND_MST_SEG_ID(I, 2:5)
          node_id_load(K+4:K+7) = CAND_SEC_SEG_ID(I, 2:5)
          K = K + 8
      
          ! Store forces: Primary (1-4), Secondary (5-8)
          DO J = 1, 4
            ! Primary forces
            load_arr(L, J, 1:3) = p_load_new(3*(J-1) + 1 : 3*J)
            ! Secondary forces
            load_arr(L, J + 4, 1:3) = p_load_new(12 + 3*(J-1) + 1 : 12 + 3*J)
          ENDDO
      
          ! Store stiffness info for all 8 nodes
          load_arr(L, 1:8, 4) = node_stiff(1:8)
      
          L = L + 1
          
          ! Safety check - prevent array overflow
          IF (L > MAX_STS_SIZE_ACTUAL .OR. K > MAX_STS_SIZE_ACTUAL*8) THEN
            EXIT
          END IF
        ENDIF
      ENDDO
      IF (CSV_OUTPUT_ENABLED) THEN
        CLOSE(LUX_STS)
      END IF

      DEALLOCATE(lobatto_gp_weight)
      
      L_out = L
      END SUBROUTINE STS_CONTACTS_ASSEMBLE

      SUBROUTINE STS_BUILD_LOBATTO_GP_WEIGHTS(COUNT_IN, CAPACITY, &
     & CAND_SEC_SEG_ID, CAND_MST_SEG_ID, CAND_SEC_GP_MASK, WEIGHT)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: COUNT_IN, CAPACITY
      INTEGER, INTENT(IN) :: CAND_SEC_SEG_ID(CAPACITY,5)
      INTEGER, INTENT(IN) :: CAND_MST_SEG_ID(CAPACITY,5)
      INTEGER, INTENT(IN) :: CAND_SEC_GP_MASK(CAPACITY,4)
      REAL*8, INTENT(INOUT) :: WEIGHT(4, COUNT_IN)
      INTEGER :: HASH_SIZE, I, C, SEC_SEG, MST_SEG
      INTEGER :: IDX, PROBE, CNT
      INTEGER, ALLOCATABLE :: HASH_SEC_SEG(:), HASH_CORNER(:)
      INTEGER, ALLOCATABLE :: HASH_COUNT(:)
      INTEGER(KIND=8) :: HKEY

      IF (COUNT_IN <= 0) RETURN
      HASH_SIZE = MAX(17, 8*COUNT_IN + 1)
      ALLOCATE(HASH_SEC_SEG(HASH_SIZE), HASH_CORNER(HASH_SIZE))
      ALLOCATE(HASH_COUNT(HASH_SIZE))
      HASH_SEC_SEG = 0
      HASH_CORNER = 0
      HASH_COUNT = 0

!     A Lobatto corner belongs to one secondary segment. INT7 remapping
!     can expand that same corner to several master patch segments. Keep
!     full 2x2 secondary coverage, but split the corner weight only over
!     that segment's local master-patch candidates. Do not split globally
!     by node, because adjacent secondary segments legitimately share
!     physical nodes and each segment still owns its own integrated area.
      DO I = 1, COUNT_IN
        SEC_SEG = CAND_SEC_SEG_ID(I,1)
        MST_SEG = CAND_MST_SEG_ID(I,1)
        IF (SEC_SEG <= 0 .OR. MST_SEG <= 0) CYCLE
        DO C = 1, 4
          IF (CAND_SEC_GP_MASK(I,C) == 0) CYCLE
          HKEY = INT(SEC_SEG, KIND=8) * INT(8, KIND=8) + &
     &      INT(C, KIND=8)
          IDX = INT(MOD(HKEY, INT(HASH_SIZE, KIND=8))) + 1
          DO PROBE = 1, HASH_SIZE
            IF (HASH_COUNT(IDX) == 0) THEN
              HASH_SEC_SEG(IDX) = SEC_SEG
              HASH_CORNER(IDX) = C
              EXIT
            ELSE IF (HASH_SEC_SEG(IDX) == SEC_SEG .AND. &
     &               HASH_CORNER(IDX) == C) THEN
              EXIT
            ENDIF
            IDX = IDX + 1
            IF (IDX > HASH_SIZE) IDX = 1
          ENDDO
          HASH_COUNT(IDX) = HASH_COUNT(IDX) + 1
        ENDDO
      ENDDO

      DO I = 1, COUNT_IN
        SEC_SEG = CAND_SEC_SEG_ID(I,1)
        MST_SEG = CAND_MST_SEG_ID(I,1)
        IF (SEC_SEG <= 0 .OR. MST_SEG <= 0) CYCLE
        DO C = 1, 4
          IF (CAND_SEC_GP_MASK(I,C) == 0) CYCLE
          HKEY = INT(SEC_SEG, KIND=8) * INT(8, KIND=8) + &
     &      INT(C, KIND=8)
          IDX = INT(MOD(HKEY, INT(HASH_SIZE, KIND=8))) + 1
          CNT = 1
          DO PROBE = 1, HASH_SIZE
            IF (HASH_COUNT(IDX) == 0) EXIT
            IF (HASH_SEC_SEG(IDX) == SEC_SEG .AND. &
     &          HASH_CORNER(IDX) == C) THEN
              CNT = HASH_COUNT(IDX)
              EXIT
            ENDIF
            IDX = IDX + 1
            IF (IDX > HASH_SIZE) IDX = 1
          ENDDO
          WEIGHT(C,I) = 1.0D0 / DBLE(CNT)
        ENDDO
      ENDDO

      DEALLOCATE(HASH_SEC_SEG, HASH_CORNER, HASH_COUNT)
      END SUBROUTINE STS_BUILD_LOBATTO_GP_WEIGHTS
      
      SUBROUTINE STS_CONTACT_EXPORT_CSV_PAIR(LUX_STS, NCYCLE_IN, TIME_CUR, STS_INTERFACE_ID, &
     & MST_ENTITY_ID, SEC_ENTITY_ID, pair_max_penetration, &
     & FX_PRIM, FY_PRIM, FZ_PRIM, FXF_PRIM, FYF_PRIM, FZF_PRIM, &
     & FX_SEC, FY_SEC, FZ_SEC, FXF_SEC, FYF_SEC, FZF_SEC)
      IMPLICIT NONE
      INTEGER LUX_STS, NCYCLE_IN, STS_INTERFACE_ID, MST_ENTITY_ID, SEC_ENTITY_ID
      REAL*8 TIME_CUR, pair_max_penetration
      REAL*8 FX_PRIM, FY_PRIM, FZ_PRIM, FXF_PRIM, FYF_PRIM, FZF_PRIM
      REAL*8 FX_SEC, FY_SEC, FZ_SEC, FXF_SEC, FYF_SEC, FZF_SEC
      REAL*8 force_norm, fn_mag, ft_mag

      force_norm = SQRT(FX_PRIM*FX_PRIM + FY_PRIM*FY_PRIM + FZ_PRIM*FZ_PRIM)
      ft_mag = SQRT(FXF_PRIM*FXF_PRIM + FYF_PRIM*FYF_PRIM + FZF_PRIM*FZF_PRIM)
      fn_mag = SQRT(MAX(0.0D0, (FX_PRIM-FXF_PRIM)*(FX_PRIM-FXF_PRIM) + &
     &    (FY_PRIM-FYF_PRIM)*(FY_PRIM-FYF_PRIM) + (FZ_PRIM-FZF_PRIM)*(FZ_PRIM-FZF_PRIM)))
      WRITE(LUX_STS,'(I0,'','',ES23.15,'','',I0,'','',I0)',ADVANCE='NO') &
     &  NCYCLE_IN, TIME_CUR, ABS(STS_INTERFACE_ID), MST_ENTITY_ID
      WRITE(LUX_STS,901) FX_PRIM, FY_PRIM, FZ_PRIM, force_norm, fn_mag, ft_mag, 1, &
     &    pair_max_penetration

      force_norm = SQRT(FX_SEC*FX_SEC + FY_SEC*FY_SEC + FZ_SEC*FZ_SEC)
      ft_mag = SQRT(FXF_SEC*FXF_SEC + FYF_SEC*FYF_SEC + FZF_SEC*FZF_SEC)
      fn_mag = SQRT(MAX(0.0D0, (FX_SEC-FXF_SEC)*(FX_SEC-FXF_SEC) + &
     &    (FY_SEC-FYF_SEC)*(FY_SEC-FYF_SEC) + (FZ_SEC-FZF_SEC)*(FZ_SEC-FZF_SEC)))
      WRITE(LUX_STS,'(I0,'','',ES23.15,'','',I0,'','',I0)',ADVANCE='NO') &
     &  NCYCLE_IN, TIME_CUR, -ABS(STS_INTERFACE_ID), SEC_ENTITY_ID
      WRITE(LUX_STS,901) FX_SEC, FY_SEC, FZ_SEC, force_norm, fn_mag, ft_mag, 1, &
     &    pair_max_penetration

 901  FORMAT ( ',', ES23.15, ',', ES23.15, ',', ES23.15, ',', ES23.15, ',', &
     &    ES23.15, ',', ES23.15, ',', I0, ',', ES23.15)
      END SUBROUTINE STS_CONTACT_EXPORT_CSV_PAIR
