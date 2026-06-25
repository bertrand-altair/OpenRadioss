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
!||    ists_sts_bp_algo_mod  ../engine/source/interfaces/ists/ists_sts_bp_algo_mod.F90
!||--- called by ------------------------------------------------------
!||    ists_mainf              ../engine/source/interfaces/ists/ists_mainf.F
!||    inter_sort_07           ../engine/source/interfaces/int07/inter_sort_07.F
!||    i7main_tri              ../engine/source/interfaces/intsort/i7main_tri.F
!||====================================================================
!
!   Hardcoded broad-phase algorithm selector for STS contact.
!
!   STS_BP_ALGO_VOXEL       : use the STS-native voxel broad phase
!                             (STS_VOXEL_BROAD_PHASE).
!   STS_BP_ALGO_INT7_BUCKET : reuse the legacy INT7 bucket sorting
!                             (I7BUCE/I7TRI) candidate arrays and map
!                             them to STS segment pairs via
!                             STS_REMAP_SEGMENTS.
!
!   STS_BP_ALGO is a PARAMETER so it is patched at build time by the
!   study toolchain (patch_engine_constants.py); no input/deck change
!   is required to switch between the two paths.
!
      MODULE ISTS_STS_BP_ALGO_MOD
        IMPLICIT NONE
        PRIVATE

        INTEGER, PARAMETER, PUBLIC :: STS_BP_ALGO_VOXEL       = 0 ! This needs review
        INTEGER, PARAMETER, PUBLIC :: STS_BP_ALGO_INT7_BUCKET = 1 ! Quicker Option

        INTEGER, PARAMETER, PUBLIC :: STS_BP_ALGO = STS_BP_ALGO_INT7_BUCKET ! Set here the Option to use for the broad-phase algorithm.

      END MODULE ISTS_STS_BP_ALGO_MOD
