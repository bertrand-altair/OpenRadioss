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
!||    sph_dormant_contact_wakeup   ../engine/source/interfaces/interf/sph_dormant_contact_wakeup.F90
!||--- called by ------------------------------------------------------
!||    soltosph_on1                 ../engine/source/elements/sph/soltosph_on1.F
!||--- uses       -----------------------------------------------------
!||    intbufdef_mod                ../common_source/modules/interfaces/intbufdef_mod.F90
!||    shooting_node_mod            ../engine/share/modules/shooting_node_mod.F90
!||====================================================================
! ======================================================================================================================
!                                                   PROCEDURES
! ======================================================================================================================
!
!=======================================================================================================================
!! \brief Reactivate contact for a single SOL2SPH particle that just woke up.
!=======================================================================================================================
!
      subroutine sph_dormant_contact_wakeup(node_id, ninter, shoot_struct, intbuf_tab, newfront)
! ----------------------------------------------------------------------------------------------------------------------
!                                                   MODULES
! ----------------------------------------------------------------------------------------------------------------------
        use intbufdef_mod, only : intbuf_struct_
        use shooting_node_mod, only : shooting_node_type
! ----------------------------------------------------------------------------------------------------------------------
!                                                   IMPLICIT NONE
! ----------------------------------------------------------------------------------------------------------------------
        implicit none
! ----------------------------------------------------------------------------------------------------------------------
!                                                   ARGUMENTS
! ----------------------------------------------------------------------------------------------------------------------
        integer,                               intent(in)      :: node_id          !< reactivated SOL2SPH node (KXSP(3,NP))
        integer,                               intent(in)      :: ninter           !< number of interfaces
        integer,                               intent(inout)   :: newfront(ninter) !< force resort flag per interface
        type(shooting_node_type),              intent(inout)   :: shoot_struct     !< reverse map + SOL2SPH dormant-contact state
        type(intbuf_struct_), dimension(ninter), intent(inout) :: intbuf_tab       !< interface data (STFNS is updated)
! ----------------------------------------------------------------------------------------------------------------------
!                                                   LOCAL VARIABLES
! ----------------------------------------------------------------------------------------------------------------------
        integer :: e, nin, i
! ----------------------------------------------------------------------------------------------------------------------
!                                                   BODY
! ----------------------------------------------------------------------------------------------------------------------
        ! Restore contact on every interface where this node was deactivated by the dormant-contact logic.
        do e = shoot_struct%shift_s_node(node_id) + 1, shoot_struct%shift_s_node(node_id+1)
          if (shoot_struct%sph_dormant(e) == 1) then
            nin = shoot_struct%inter_sec_node(e)
            i   = shoot_struct%sec_node_id(e)
            intbuf_tab(nin)%stfns(i) = shoot_struct%sph_stfns_sav(e)
            shoot_struct%sph_dormant(e) = 0
            ! NEWFRONT is shared across threads (per interface); the write is idempotent (-1) but made atomic.
!$OMP ATOMIC WRITE
            newfront(nin) = -1
          end if
        end do

      end subroutine sph_dormant_contact_wakeup
