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
!||    sph_dormant_contact_init   ../engine/source/interfaces/interf/sph_dormant_contact_init.F90
!||--- called by ------------------------------------------------------
!||    init_nodal_state           ../engine/source/interfaces/interf/init_nodal_state.F
!||--- uses       -----------------------------------------------------
!||    intbufdef_mod              ../common_source/modules/interfaces/intbufdef_mod.F90
!||    precision_mod              ../common_source/modules/precision_mod.F90
!||====================================================================
! ======================================================================================================================
!                                                   PROCEDURES
! ======================================================================================================================
!
!=======================================================================================================================
!! \brief Initial deactivation of contact for SOL2SPH particles that are already dormant.
!=======================================================================================================================
!
      subroutine sph_dormant_contact_init(ninter, numnod, nisp, kxsp, nod2sp, &
                                          size_sec_node, shift_s_node, inter_sec_node, sec_node_id, &
                                          sph_stfns_sav, sph_dormant, intbuf_tab)
! ----------------------------------------------------------------------------------------------------------------------
!                                                   MODULES
! ----------------------------------------------------------------------------------------------------------------------
        use precision_mod, only : WP
        use intbufdef_mod, only : intbuf_struct_
! ----------------------------------------------------------------------------------------------------------------------
!                                                   IMPLICIT NONE
! ----------------------------------------------------------------------------------------------------------------------
        implicit none
! ----------------------------------------------------------------------------------------------------------------------
!                                                   ARGUMENTS
! ----------------------------------------------------------------------------------------------------------------------
        integer,                               intent(in)      :: ninter           !< number of interfaces
        integer,                               intent(in)      :: numnod           !< number of nodes
        integer,                               intent(in)      :: nisp             !< first dimension of KXSP
        integer,                               intent(in)      :: kxsp(nisp,*)     !< SPH connectivity (dormancy in row 2)
        integer,                               intent(in)      :: nod2sp(numnod)   !< node -> SPH particle map (0 if not SPH)
        integer,                               intent(in)      :: size_sec_node    !< size of INTER_SEC_NODE/SEC_NODE_ID
        integer,                               intent(in)      :: shift_s_node(numnod+1)        !< per-node shift into the reverse map
        integer,                               intent(in)      :: inter_sec_node(size_sec_node) !< interface id per map entry
        integer,                               intent(in)      :: sec_node_id(size_sec_node)    !< local secondary node index per map entry
        integer,                               intent(inout)   :: sph_dormant(size_sec_node)    !< 1 if entry currently forced dormant
        real(kind=WP),                         intent(inout)   :: sph_stfns_sav(size_sec_node)  !< saved active stiffness per map entry
        type(intbuf_struct_), dimension(ninter), intent(inout) :: intbuf_tab       !< interface data (STFNS is updated)
! ----------------------------------------------------------------------------------------------------------------------
!                                                   LOCAL VARIABLES
! ----------------------------------------------------------------------------------------------------------------------
        integer :: k, e, nin, i, np
! ----------------------------------------------------------------------------------------------------------------------
!                                                   BODY
! ----------------------------------------------------------------------------------------------------------------------
        do k = 1, numnod
          np = nod2sp(k)
          if (np <= 0) cycle
          if (kxsp(2,np) < 0) then
            do e = shift_s_node(k) + 1, shift_s_node(k+1)
              nin = inter_sec_node(e)
              i   = sec_node_id(e)
              if (sph_dormant(e) == 0) then
                sph_stfns_sav(e) = intbuf_tab(nin)%stfns(i)
                intbuf_tab(nin)%stfns(i) = 0.0_WP
                sph_dormant(e) = 1
              end if
            end do
          end if
        end do

      end subroutine sph_dormant_contact_init
