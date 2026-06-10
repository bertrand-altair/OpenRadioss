Copyright>        OpenRadioss
Copyright>        Copyright (C) 1986-2026 Altair Engineering Inc.
Copyright>
Copyright>        This program is free software: you can redistribute it and/or modify
Copyright>        it under the terms of the GNU Affero General Public License as published by
Copyright>        the Free Software Foundation, either version 3 of the License, or
Copyright>        (at your option) any later version.
Copyright>
Copyright>        This program is distributed in the hope that it will be useful,
Copyright>        but WITHOUT ANY WARRANTY; without even the implied warranty of
Copyright>        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
Copyright>        GNU Affero General Public License for more details.
Copyright>
Copyright>        You should have received a copy of the GNU Affero General Public License
Copyright>        along with this program.  If not, see <https://www.gnu.org/licenses/>.
Copyright>
Copyright>
Copyright>        Commercial Alternative: Altair Radioss Software
Copyright>
Copyright>        As an alternative to this open-source version, Altair also offers Altair Radioss
Copyright>        software under a commercial license.  Contact Altair to discuss further if the
Copyright>        commercial version may interest you: https://www.altair.com/radioss/.
!||====================================================================
!||    i7mainf                     ../engine/source/interfaces/int07/i7mainf.F
!||--- called by ------------------------------------------------------
!||    intfop2                     ../engine/source/interfaces/interf/intfop2.F
!||--- calls      -----------------------------------------------------
!||    ancmsg                      ../engine/source/output/message/message.F
!||    arret                       ../engine/source/system/arret.F
!||    frictionparts_model_isot    ../engine/source/interfaces/int07/frictionparts_model.F
!||    frictionparts_model_ortho   ../engine/source/interfaces/int07/frictionparts_model.F
!||    i17main_pena                ../engine/source/interfaces/int17/i17main_pena.F
!||    i18dst3                     ../engine/source/interfaces/int18/i18dst3.F
!||    i18for3                     ../engine/source/interfaces/int18/i18for3.F
!||    i7ass33                     ../engine/source/interfaces/int07/i7ass3.F
!||    i7cdcor3                    ../engine/source/interfaces/int07/i7cdcor3.F
!||    i7cor3                      ../engine/source/interfaces/int07/i7cor3.F
!||    i7curvsz                    ../engine/source/interfaces/int07/i7curv.F
!||    i7dst3                      ../engine/source/interfaces/int07/i7dst3.F
!||    i7for3                      ../engine/source/interfaces/int07/i7for3.F
!||    i7norm                      ../engine/source/interfaces/int07/i7curv.F
!||    i7norme                     ../engine/source/interfaces/int07/i7rcurv.F
!||    i7normn                     ../engine/source/interfaces/int07/i7rcurv.F
!||    i7normnp                    ../engine/source/interfaces/int07/i7rcurv.F
!||    i7normp                     ../engine/source/interfaces/int07/i7curv.F
!||    i7rcurv                     ../engine/source/interfaces/int07/i7rcurv.F
!||    i7therm                     ../engine/source/interfaces/int07/i7therm.F
!||    my_barrier                  ../engine/source/system/machine.F
!||    spmd_exch_n                 ../engine/source/mpi/generic/spmd_exch_n.F
!||    spmd_i7curvsz               ../engine/source/mpi/interfaces/spmd_i7curvsz.F
!||    startime                    ../engine/source/system/timer_mod.F90
!||    stoptime                    ../engine/source/system/timer_mod.F90
!||    sum_6_float_sens            ../engine/source/system/parit.F
!||--- uses       -----------------------------------------------------
!||    ale_connectivity_mod        ../common_source/modules/ale/ale_connectivity_mod.F
!||    elbufdef_mod                ../common_source/modules/mat_elem/elbufdef_mod.F90
!||    groupdef_mod                ../common_source/modules/groupdef_mod.F
!||    h3d_mod                     ../engine/share/modules/h3d_mod.F
!||    i18dst3_mod                 ../engine/source/interfaces/int18/i18dst3.F
!||    i18for3_mod                 ../engine/source/interfaces/int18/i18for3.F
!||    intbuf_fric_mod             ../common_source/modules/interfaces/intbuf_fric_mod.F90
!||    intbufdef_mod               ../common_source/modules/interfaces/intbufdef_mod.F90
!||    message_mod                 ../engine/share/message_module/message_mod.F
!||    multi_fvm_mod               ../common_source/modules/ale/multi_fvm_mod.F90
!||    output_mod                  ../common_source/modules/output/output_mod.F90
!||    timer_mod                   ../engine/source/system/timer_mod.F90
!||====================================================================
      SUBROUTINE ISTS_MAINF(TIMERS,
     1                   IPARI          ,X              ,A              ,
     2                   STIFN          ,FCONT          ,NIN            ,
     3                   JTASK          ,NB_JLT_NEW     ,
     4                   TEMP           ,NPC            ,TF             ,
     5                   INTBUF_TAB     ,IXS            ,
     6                   H3D_DATA       ,INTBUF_FRIC_TAB,IGRSURF)
C======================================================================
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE TIMER_MOD
      USE INTBUFDEF_MOD
      USE H3D_MOD
      USE INTBUF_FRIC_MOD
      USE GROUPDEF_MOD
      USE sts_gp_state_mod
      USE Q1NP_CONTACT_DRIVER_MOD, ONLY: Q1NP_CONTACT_DRIVER_INT7,
     &  Q1NP_CONTACT_INIT_GRID_NODES
      USE STS_VOXEL_DRIVER_MOD, ONLY: STS_VOXEL_DRIVER_DATA_READY
      USE STS_BROAD_PHASE_VOXEL_MOD, ONLY: STS_VOXEL_BROAD_PHASE
      USE STS_CONTACT_STIFFNESS_MOD, ONLY: STS_CONTACT_STIFFNESS
C-----------------------------------------------
C   I m p l i c i t   T y p e s
C-----------------------------------------------
#include      "implicit_f.inc"
C-----------------------------------------------
C   G l o b a l   P a r a m e t e r s
C-----------------------------------------------
#include      "mvsiz_p.inc"
#include      "elements.inc"
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
#include      "com01_c.inc"
#include      "com04_c.inc"
#include      "com08_c.inc"
#include      "param_c.inc"
#include      "warn_c.inc"
#include      "task_c.inc"
#include      "parit_c.inc"
#include      "timeri_c.inc"
#include      "impl1_c.inc"
#include      "scr14_c.inc"
#include      "scr16_c.inc"
#include      "macro.inc"
#include      "comlock.inc"
#include      "tabsiz_c.inc"
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      TYPE(TIMER_), INTENT(INOUT) :: TIMERS
      TYPE(SURF_)   ,DIMENSION(NSURF)   :: IGRSURF
      INTEGER NIN
      INTEGER IPARI(NPARI,NINTER)
      INTEGER NB_JLT_NEW,JTASK
      INTEGER IXS(*), NPC(*)
      my_real X(3,*), A(3,*), STIFN(*), FCONT(3,*)
      my_real TEMP(*), TF(*)
      TYPE(INTBUF_STRUCT_) INTBUF_TAB
      TYPE(H3D_DATABASE) :: H3D_DATA
      TYPE(INTBUF_FRIC_STRUCT_), TARGET, DIMENSION(NINTERFRIC) ::  INTBUF_FRIC_TAB
C-----------------------------------------------
C     N E W   V a r i a b l e s
C-----------------------------------------------
      my_real, DIMENSION(:,:,:), ALLOCATABLE :: CONT_ELEMENT        !contact element array
      my_real, DIMENSION(:,:,:), ALLOCATABLE :: load_arr            !load array
      my_real, DIMENSION(:),     ALLOCATABLE :: STS_STIF            !STS pair stiffness
      INTEGER, DIMENSION(:,:),     ALLOCATABLE :: CAND_SEC_SEG_ID   !secondary segment nodes
      INTEGER, DIMENSION(:,:),     ALLOCATABLE :: CAND_MST_SEG_ID   !primary segment nodes
      INTEGER, DIMENSION(:),       ALLOCATABLE :: node_id_load      !node id array (8 nodes per element)
      INTEGER COUNT !count for the number of contact elements found
      INTEGER :: MAX_STS_SIZE_ACTUAL !actual calculated sts array size needed      
      INTEGER OPTION !option for the contact algorithm
      INTEGER IMPACT_glob !global impact flag, 0 if no impact, 1 if impact
      LOGICAL :: STS_ACTIVE !true if STS is active, false if NTS is active
      INTEGER :: IP_MAX
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER I, L, IFQ, MFROT, IGSTI,
     .        JLT_TIED, INTTH, IFORM,
     .        INTFRIC, NSETPRTS, NPARTFRIC, IORTHFRIC
      INTEGER IX1(MVSIZ), IX2(MVSIZ), IX3(MVSIZ), IX4(MVSIZ),
     .        IPARTFRICSI(MVSIZ), IPARTFRICMI(MVSIZ)
      my_real STARTT, FRIC, GAP, STOPT, VISCF, STIGLO, GAPMIN,
     .        KMIN, KMAX, GAPMAX, GAP_FOR_STS, TINT
C-----------------------------------------------
      my_real H1(MVSIZ), H2(MVSIZ), H3(MVSIZ), H4(MVSIZ),
     .        TEMPI(MVSIZ)
      INTEGER IFRIC

      my_real XFILTR_FRIC,FRIC_COEFS(MVSIZ,10),VISCFFRIC(MVSIZ),
     .        FRICC(MVSIZ), XMU(MVSIZ), QFRICT
      INTEGER, DIMENSION(:) ,POINTER  :: TABCOUPLEPARTS_FRIC
      INTEGER, DIMENSION(:) ,POINTER  :: TABPARTS_FRIC 
      INTEGER, DIMENSION(:) ,POINTER  :: ADPARTS_FRIC 
      INTEGER, DIMENSION(:) ,POINTER  :: IFRICORTH 
      my_real, DIMENSION(:) ,POINTER  :: TABCOEF_FRIC
      INTEGER,TARGET, DIMENSION(1):: TABCOUPLEPARTS_FRIC_BID
      INTEGER,TARGET, DIMENSION(1):: TABPARTS_FRIC_BID
      INTEGER,TARGET, DIMENSION(1):: ADPARTS_FRIC_BID
      INTEGER,TARGET, DIMENSION(1):: IFRICORTH_BID
      my_real,TARGET, DIMENSION(1):: TABCOEF_FRIC_BID
      INTEGER :: NSN, NTY, NOINT
      INTEGER :: NRTM

      INTEGER :: STS_FLAG
      INTEGER :: SEC_SURF_ID_STS
      INTEGER :: MST_SURF_ID_STS
      INTEGER :: NSEG_SEC, NSEG_MST
!     Secondary / master surfaces from IPARI(45/46); voxel broad-phase
!     populates candidate segment pairs up to MAX_STS_SIZE_ACTUAL.
      LOGICAL :: STS_VOXEL_OVERFLOW                       ! voxel broad-phase overflow flag
      LOGICAL, SAVE :: INT7_CONTACT_MODE_MSG_DONE = .FALSE. ! one-time flag for the contact-mode banner
      LOGICAL, SAVE :: STS_VOXEL_OVERFLOW_WARN_DONE = .FALSE.
C-----------------------------------------------
C   S o u r c e   L i n e s
C-----------------------------------------------
      STS_ACTIVE = .FALSE. ! .FALSE. = NTS (node-to-surface); .TRUE. = STS (segment-to-segment)      
      IF (IPARI(105,NIN) == 1) STS_ACTIVE = .TRUE.
C      
      NRTM = IPARI(4,NIN)
      NSN = IPARI(5,NIN)
      NTY = IPARI(7,NIN)
      IF(IPARI(33,NIN) == 1) RETURN  !precondition
      NOINT = IPARI(15,NIN)
      MFROT = IPARI(30,NIN)
      IFQ = IPARI(31,NIN)
      IGSTI = IPARI(34,NIN)
C heat interface
      INTTH = IPARI(47,NIN)
      IFORM = IPARI(48,NIN)
      STIGLO = -INTBUF_TAB%STFAC(1)
      STARTT = INTBUF_TAB%VARIABLES(3)
      STOPT = INTBUF_TAB%VARIABLES(11)
      IF(STARTT > TT) RETURN !precondition
      IF(TT > STOPT)  RETURN !precondition
C surface to surface contact i(=1 linear, =2 nurbs q1np)
      STS_FLAG = IPARI(105,NIN)    
C-----------------------------------------------
C     One-time diagnostic message: which contact branch applies (first cycle).
C-----------------------------------------------
      IF (NTY == 7 .AND. JTASK == 1 .AND. .NOT. INT7_CONTACT_MODE_MSG_DONE) THEN
        IF (ABS(TT) <= 1.0E-20 .OR. NCYCLE <= 0) THEN
          IF (STS_FLAG == 2) THEN
            WRITE(6,*) 'INT7 (NIN=', NIN, '): contact mode = Q1NP (NURBS)'
          ELSE IF (STS_ACTIVE .AND. STS_VOXEL_DRIVER_DATA_READY(IGRSURF, NSURF, NOINT,
     &                .FALSE., SEC_SURF_ID_STS, MST_SURF_ID_STS)) THEN
            WRITE(6,*) 'INT7 (NIN=', NIN, '): contact mode = STS (voxel broad-phase)'
          ELSE
            WRITE(6,*) 'INT7 (NIN=', NIN, '): contact mode = NTS (node-to-surface)'
          ENDIF
          INT7_CONTACT_MODE_MSG_DONE = .TRUE.
        ENDIF
      ENDIF
C-----------------------------------------------
      FRIC = INTBUF_TAB%VARIABLES(1)
      GAP = INTBUF_TAB%VARIABLES(2)
      GAPMIN = INTBUF_TAB%VARIABLES(13)
      VISCF= INTBUF_TAB%VARIABLES(15)
      IFRIC = 0
      IF(INTTH > 0) IFRIC =IPARI(50,NIN)
      GAPMAX = INTBUF_TAB%VARIABLES(16)
C     Scalar gap consistent with I7COR3:
      GAP_FOR_STS = MIN(GAPMAX, MAX(GAPMIN, GAP))
      KMIN   = INTBUF_TAB%VARIABLES(17)
      KMAX   = INTBUF_TAB%VARIABLES(18)
      TINT   = INTBUF_TAB%VARIABLES(22)

C--- Corresponding Friction model  
      INTFRIC=IPARI(72,NIN)
      IORTHFRIC = 0
      NSETPRTS = 0
      NPARTFRIC = 0 
      XFILTR_FRIC = ZERO
      
      IF(INTFRIC /= 0) THEN 
         TABCOUPLEPARTS_FRIC => INTBUF_FRIC_TAB(INTFRIC)%TABCOUPLEPARTS_FRIC
         TABCOEF_FRIC => INTBUF_FRIC_TAB(INTFRIC)%TABCOEF_FRIC
         TABPARTS_FRIC => INTBUF_FRIC_TAB(INTFRIC)%TABPARTS_FRIC
         ADPARTS_FRIC => INTBUF_FRIC_TAB(INTFRIC)%ADPARTS_FRIC
         XFILTR_FRIC = INTBUF_FRIC_TAB(INTFRIC)%XFILTR_FRIC
         NSETPRTS = INTBUF_FRIC_TAB(INTFRIC)%NSETPRTS
         NPARTFRIC = INTBUF_FRIC_TAB(INTFRIC)%S_TABPARTS_FRIC
         IORTHFRIC = INTBUF_FRIC_TAB(INTFRIC)%IORTHFRIC
         IFRICORTH => INTBUF_FRIC_TAB(INTFRIC)%IFRICORTH
      ELSE
         TABCOUPLEPARTS_FRIC => TABCOUPLEPARTS_FRIC_BID
         TABPARTS_FRIC => TABPARTS_FRIC_BID
         TABCOEF_FRIC => TABCOEF_FRIC_BID
         ADPARTS_FRIC => ADPARTS_FRIC_BID
         IFRICORTH => IFRICORTH_BID
         IF (IFQ/=0) XFILTR_FRIC = INTBUF_TAB%XFILTR(1)
      ENDIF

C-----------------------------------------------
C     Type 7 STS (STS_FLAG == 1): voxel broad-phase via IPARI(45/46).
C-----------------------------------------------
      IF (NTY == 7 .AND. STS_FLAG == 1) THEN
C
        MST_SURF_ID_STS = IPARI(46,NIN)
        SEC_SURF_ID_STS = IPARI(45,NIN)
C
        IF (JTASK == 1) THEN
C          Set storage capacity for the candidate arrays.
C          Upper bound: all unique (sec, mst) pairs; practical cap:
C          NSEG_sec * 512 (matches prior i7mainf STS voxel sizing).
          NSEG_SEC = IGRSURF(SEC_SURF_ID_STS)%NSEG
          NSEG_MST = IGRSURF(MST_SURF_ID_STS)%NSEG
          MAX_STS_SIZE_ACTUAL = MAX(1, MIN(NSEG_SEC * NSEG_MST,
     &                                    NSEG_SEC * 512))

C         Allocate / refresh the working arrays consumed by
C         STS_CONTACTS_ASSEMBLE and I7ASS33 (ASS_FLAG = 2).
          IF (ALLOCATED(CONT_ELEMENT))    DEALLOCATE(CONT_ELEMENT)
          IF (ALLOCATED(load_arr))        DEALLOCATE(load_arr)
          IF (ALLOCATED(STS_STIF))        DEALLOCATE(STS_STIF)
          IF (ALLOCATED(CAND_SEC_SEG_ID)) DEALLOCATE(CAND_SEC_SEG_ID)
          IF (ALLOCATED(CAND_MST_SEG_ID)) DEALLOCATE(CAND_MST_SEG_ID)
          IF (ALLOCATED(node_id_load))    DEALLOCATE(node_id_load)
          ALLOCATE(CONT_ELEMENT(MAX_STS_SIZE_ACTUAL,3,8))
          ALLOCATE(load_arr(MAX_STS_SIZE_ACTUAL,8,4))
          ALLOCATE(STS_STIF(MAX_STS_SIZE_ACTUAL))
          ALLOCATE(CAND_SEC_SEG_ID(MAX_STS_SIZE_ACTUAL,5))
          ALLOCATE(CAND_MST_SEG_ID(MAX_STS_SIZE_ACTUAL,5))
          ALLOCATE(node_id_load(MAX_STS_SIZE_ACTUAL*8))
          CONT_ELEMENT  = ZERO
          load_arr      = ZERO
          STS_STIF      = ZERO
          CAND_SEC_SEG_ID = 0
          CAND_MST_SEG_ID = 0
          node_id_load    = 0

          IP_MAX = 2
          CALL sts_gp_state_init(MAX_STS_SIZE_ACTUAL, IP_MAX)

C         Voxel broad phase populates CAND_*_SEG_ID and CONT_ELEMENT.
          IF ((IMON>0).AND.(NCYCLE > 37290)) CALL STARTIME(TIMERS,35)
          CALL STS_VOXEL_BROAD_PHASE(
     &          IGRSURF, NSURF, SEC_SURF_ID_STS, MST_SURF_ID_STS,
     &          X, NUMNOD, GAP_FOR_STS, MAX_STS_SIZE_ACTUAL,
     &          CAND_SEC_SEG_ID, CAND_MST_SEG_ID, CONT_ELEMENT,
     &          COUNT, STS_VOXEL_OVERFLOW)
          IF ((IMON>0).AND.(NCYCLE > 37290)) CALL STOPTIME(TIMERS,35)

          IF (STS_VOXEL_OVERFLOW .AND. .NOT. STS_VOXEL_OVERFLOW_WARN_DONE) THEN
            WRITE(*,'(A,I8,A,I10,A,I10)') 
     &        ' WARNING INTERFACE ', NOINT, 
     &        ': STS voxel pair overflow, kept ', COUNT, 
     &        ' / capacity ', MAX_STS_SIZE_ACTUAL
            STS_VOXEL_OVERFLOW_WARN_DONE = .TRUE.
          END IF

C         Friction model 
C         TODO: placeholder part-ID mapping; proper part-ID resolution is a separate task.
          IF (COUNT > 0) THEN
            DO I = 1, MIN(COUNT, MVSIZ)
              IPARTFRICSI(I) = 0
              IPARTFRICMI(I) = 0
            ENDDO
            IF (IORTHFRIC == 0) THEN
              JLT_TIED = 0
              CALL FRICTIONPARTS_MODEL_ISOT(
     1               INTFRIC,        MIN(COUNT,MVSIZ),  IPARTFRICSI,
     2               IPARTFRICMI,    ADPARTS_FRIC,
     3               NSETPRTS,       TABCOUPLEPARTS_FRIC,NPARTFRIC,
     4               TABPARTS_FRIC,  TABCOEF_FRIC,
     5               FRIC,           VISCF,            INTBUF_TAB%FRIC_P,
     6               FRIC_COEFS,     FRICC,
     7               VISCFFRIC,      NTY,              MFROT,
     8               IORTHFRIC,      IFRIC,
     9               JLT_TIED,       TINT,             TEMPI,
     A               NPC,            TF,
     B               TEMP,           H1,               H2,
     C               H3,             H4,
     D               IX1,            IX2,              IX3,
     E               IX4,            IFORM)
            ENDIF

C     	    Calculate Contact Stiffness for STS pairs.
            CALL STS_CONTACT_STIFFNESS(
     &        CAND_MST_SEG_ID, CAND_SEC_SEG_ID, COUNT,
     &        MAX_STS_SIZE_ACTUAL,
     &        INTBUF_TAB%IRECTM, INTBUF_TAB%STFM, NRTM,
     &        INTBUF_TAB%NSV, INTBUF_TAB%STFNS, NSN, NUMNOD,
     &        IGSTI, KMIN, KMAX, STIGLO, STS_STIF)

C           STS Gauss contact search.
            OPTION = 0
            QFRICT = 0.d0
            CALL STS_CONTACTS_ASSEMBLE(CONT_ELEMENT, COUNT, OPTION,
     &        NOINT, NCYCLE, TT,
     &        CAND_MST_SEG_ID, CAND_SEC_SEG_ID,
     &        load_arr, node_id_load, L, IMPACT_glob, STS_STIF,
     &        MAX_STS_SIZE_ACTUAL, FRICC, XMU, INTBUF_TAB%IFPEN, QFRICT, GAP_FOR_STS)

C           Fall back to Lobatto contact search if Gauss contact search reported no contact.
            IF (IMPACT_glob == 0) THEN
              OPTION = 1
              QFRICT = 0.d0
              CALL STS_CONTACTS_ASSEMBLE(CONT_ELEMENT, COUNT, OPTION,
     &          NOINT, NCYCLE, TT,
     &          CAND_MST_SEG_ID, CAND_SEC_SEG_ID,
     &          load_arr, node_id_load, L, IMPACT_glob, STS_STIF,
     &          MAX_STS_SIZE_ACTUAL, FRICC, XMU, INTBUF_TAB%IFPEN, QFRICT, GAP_FOR_STS)
            ENDIF

C           Scatter the per-pair load_arr / node_id_load directly into A / STIFN. 
            IF (IMPACT_glob /= 0) THEN
              IPARI(29,NIN) = 1
              IF (DEBUG(3) >= 1) NB_JLT_NEW = NB_JLT_NEW + L - 1
              CALL I7ASS0_STS(A, STIFN, load_arr, node_id_load, L,
     &          MAX_STS_SIZE_ACTUAL, FCONT, H3D_DATA)
            ENDIF

          ENDIF ! COUNT > 0
        ENDIF  
        CALL MY_BARRIER()
C          
C       Array deallocation
        IF (ALLOCATED(CONT_ELEMENT)) DEALLOCATE(CONT_ELEMENT)
        IF (ALLOCATED(load_arr))     DEALLOCATE(load_arr)
        IF (ALLOCATED(CAND_SEC_SEG_ID)) DEALLOCATE(CAND_SEC_SEG_ID)
        IF (ALLOCATED(CAND_MST_SEG_ID)) DEALLOCATE(CAND_MST_SEG_ID)
        IF (ALLOCATED(node_id_load))    DEALLOCATE(node_id_load)
C        
      ELSEIF ((NTY == 7).AND.(STS_FLAG == 2)) THEN
C        
        KMIN   = INTBUF_TAB%VARIABLES(17)
        KMAX   = INTBUF_TAB%VARIABLES(18)
        ! Q1NP_CONTACT_INIT_GRID_NODES:Returns IXS, NIXS, NUMELS for FCONT scatter in the force assembly.
        CALL Q1NP_CONTACT_INIT_GRID_NODES(IXS, NIXS, NUMELS)
        CALL MY_BARRIER()
        IF (JTASK == 1) THEN
          CALL Q1NP_CONTACT_DRIVER_INT7(NCYCLE, NUMNOD, X, A, STIFN,
     &                                  GAP, IGSTI, KMIN, KMAX, INTBUF_TAB%IRECTM,
     &                                  INTBUF_TAB%NSV, INTBUF_TAB%STFNS, NSN, INTBUF_TAB%STFM, NRTM,
     &                                  FCONT, INCONV == 1 .AND.
     &                                  ANIM_V(4)+OUTP_V(4)+H3D_DATA%N_VECT_CONT > 0)
        ENDIF
        CALL MY_BARRIER()
C            
      ENDIF
C
      RETURN
      END

