!copyright>        OpenRadioss
!copyright>        Copyright (C) 2026 Siemens
!copyright>
!copyright>        This program is free software: you can redistribute it and/or modify
!copyright>        it under the terms of the GNU Affero General Public License as published by
!copyright>        the Free Software Foundation, either version 3 of the License, or
!copyright>        (at your option) any later version.
!copyright>
!copyright>        This program is distributed in the hope that it will be useful,
!copyright>        but WITHOUT ANY WARRANTY; without even the implied warranty of
!copyright>        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!copyright>        GNU Affero General Public License for more details.
!copyright>
!copyright>        You should have received a copy of the GNU Affero General Public License
!copyright>        along with this program.  If not, see <https://www.gnu.org/licenses/>.
!copyright>
!copyright>
!copyright>        Commercial Alternative: Simcenter Radioss Software
!copyright>
!copyright>        As an alternative to this open-source version, Siemens also offers Simcenter(TM) Radioss(R)
!copyright>        software under a commercial license.  Contact Siemens to discuss further if the
!copyright>        commercial version may interest you:
!copyright>        https://www.siemens.com/en-us/products/simcenter/mechanical-simulation/radioss/.
!||====================================================================
!||    sph_dormant_contact_restore   ../engine/source/interfaces/interf/sph_dormant_contact_restore.F90
!||--- called by ------------------------------------------------------
!||    wrrestp                       ../engine/source/output/restart/wrrestp.F
!||--- uses       -----------------------------------------------------
!||    intbufdef_mod                 ../common_source/modules/interfaces/intbufdef_mod.F90
!||    precision_mod                 ../common_source/modules/precision_mod.F90
!||    shooting_node_mod             ../engine/share/modules/shooting_node_mod.F90
!||====================================================================
! ======================================================================================================================
!                                                   PROCEDURES
! ======================================================================================================================
!
!=======================================================================================================================
!! \brief Temporarily restore (or re-zero) the secondary stiffness STFNS of dormant SOL2SPH nodes for the restart write.
!=======================================================================================================================
!
      subroutine sph_dormant_contact_restore(set_active, ninter, shoot_struct, intbuf_tab)
! ----------------------------------------------------------------------------------------------------------------------
!                                                   MODULES
! ----------------------------------------------------------------------------------------------------------------------
        use precision_mod, only : WP
        use intbufdef_mod, only : intbuf_struct_
        use shooting_node_mod, only : shooting_node_type
! ----------------------------------------------------------------------------------------------------------------------
!                                                   IMPLICIT NONE
! ----------------------------------------------------------------------------------------------------------------------
        implicit none
! ----------------------------------------------------------------------------------------------------------------------
!                                                   ARGUMENTS
! ----------------------------------------------------------------------------------------------------------------------
        logical,                                 intent(in)    :: set_active   !< .true. = write active stiffness, .false. = re-zero
        integer,                                 intent(in)    :: ninter       !< number of interfaces
        type(shooting_node_type),                intent(inout) :: shoot_struct !< reverse map + SOL2SPH dormant-contact state
        type(intbuf_struct_), dimension(ninter), intent(inout) :: intbuf_tab   !< interface data (STFNS is updated)
! ----------------------------------------------------------------------------------------------------------------------
!                                                   LOCAL VARIABLES
! ----------------------------------------------------------------------------------------------------------------------
        integer :: e, nin, i
! ----------------------------------------------------------------------------------------------------------------------
!                                                   BODY
! ----------------------------------------------------------------------------------------------------------------------
        if (.not. allocated(shoot_struct%sph_dormant)) return

        if (set_active) then
          do e = 1, shoot_struct%size_sec_node
            if (shoot_struct%sph_dormant(e) == 1) then
              nin = shoot_struct%inter_sec_node(e)
              i   = shoot_struct%sec_node_id(e)
              intbuf_tab(nin)%stfns(i) = shoot_struct%sph_stfns_sav(e)
            end if
          end do
        else
          do e = 1, shoot_struct%size_sec_node
            if (shoot_struct%sph_dormant(e) == 1) then
              nin = shoot_struct%inter_sec_node(e)
              i   = shoot_struct%sec_node_id(e)
              intbuf_tab(nin)%stfns(i) = 0.0_WP
            end if
          end do
        end if

      end subroutine sph_dormant_contact_restore
