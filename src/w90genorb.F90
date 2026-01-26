!-*- mode: F90 -*-!
!------------------------------------------------------------!
! This file is distributed as part of the Wannier90 code and !
! under the terms of the GNU General Public License. See the !
! file `LICENSE' in the root directory of the Wannier90      !
! distribution, or http://www.gnu.org/copyleft/gpl.txt       !
!                                                            !
! The webpage of the Wannier90 code is www.wannier.org       !
!                                                            !
! The Wannier90 code is hosted on GitHub:                    !
!                                                            !
! https://github.com/wannier-developers/wannier90            !
!------------------------------------------------------------!
!                                                            !
!  w90genorb: generate orbital matrix from uIu and uHu       !
!                                                            !
!------------------------------------------------------------!

module get_orb
  !! Module to generate orb files from uIu, uHu,
  !! and also mmn matrix
  use w90_constants, only: dp, cmplx_0, cmplx_i, cmplx_1
  use w90_error
  use w90_error, only: w90_error_type
  use w90_comms, only: w90_comm_type
  use w90_types, only: kmesh_info_type, dis_manifold_type
  use w90_get_oper, only: get_gauge_overlap_matrix

  implicit none

  integer, save :: iun_mmn, iun_uHu, iun_uIu, iun_orb, iun_sIu, iun_sHu, iun_sH
  integer :: ierr

contains

  subroutine io_error(error_msg, stdout)
    !================================================
    !
    !! Abort the code giving an error message
    !
    !================================================

    implicit none

    character(len=*), intent(in) :: error_msg
    integer :: stdout

    close (stdout)
    write (*, '(1x,a)') trim(error_msg)
    write (*, '(A)') "Error: examine the output/error file for details"
    stop
  end subroutine io_error

  !================================================!
  subroutine print_usage(stdout)
    !================================================!
    !
    !! Writes the usage of the program to stdout
    !
    !================================================!

    implicit none

    integer, intent(in) :: stdout

    write (stdout, '(A)') "Usage:"
    write (stdout, '(A)') "  w90genorb.x [ARGS] [SEEDNAME]"
    write (stdout, '(A)') "where ARGS can be one of the following:"
    write (stdout, '(A)') "  -w or --wannier"
    write (stdout, '(A)') "      Output orb matrix in Wannier gauge."
    write (stdout, '(A)') "  -b or --bloch"
    write (stdout, '(A)') "      Output orb matrix in Bloch gauge."
  end subroutine print_usage

  !================================================!
  subroutine get_seedname(stdout, seedname, wan_gauge)
    !================================================!
    !
    !! Set the seedname from the command line
    !
    !================================================!
    implicit none

    integer, intent(in) :: stdout
    character(len=50), intent(inout)  :: seedname
    logical, intent(inout) :: wan_gauge
    ! .true.  : output matrix in Bloch basis (Wannier gauge)
    ! .false. : output matrix in Bloch basis (Bloch gauge)

    integer :: num_arg
    character(len=50) :: ctemp

    num_arg = command_argument_count()
    if (num_arg == 1) then
      seedname = 'wannier'
    elseif (num_arg == 2) then
      call get_command_argument(2, seedname)
    else
      call print_usage(stdout)
      call io_error('Wrong command line arguments, see logfile for usage', stdout)
    end if

    ! If on the command line the whole seedname.win was passed, I strip the last ".win"
    if (len(trim(seedname)) .ge. 5) then
      if (seedname(len(trim(seedname)) - 4 + 1:) .eq. ".win") then
        seedname = seedname(:len(trim(seedname)) - 4)
      end if
    end if

    call get_command_argument(1, ctemp)
    if ((index(ctemp, '-w') > 0) .or. (index(ctemp, '--wannier') > 0)) then
      wan_gauge = .true.
    elseif ((index(ctemp, '-b') > 0) .or. (index(ctemp, '--bloch') > 0)) then
      wan_gauge = .false.
    else
      write (stdout, '(A)') 'Wrong command line action: '//trim(ctemp)
      call print_usage(stdout)
      call io_error('Wrong command line arguments, see logfile for usage', stdout)
    end if

  end subroutine get_seedname

  subroutine get_mmn(stdout, seedname, num_bands, num_kpts, num_wann, kmesh_info, &
                     eigval, mmn, H_o, error, comm)
    !================================================!
    !
    !! read MMN matrix from seedname.mmn
    !
    !================================================!
    ! only run on root node
    implicit none

    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error
    type(kmesh_info_type), intent(in) :: kmesh_info
    
    character(len=60) :: header
    integer, intent(in) :: stdout
    character(len=50), intent(inout) :: seedname
    integer, intent(in) :: num_bands, num_kpts, num_wann
    real(kind=dp), intent(in) :: eigval(:, :)
    complex(kind=dp), allocatable, intent(inout) :: mmn(:, :, :, :)
    complex(kind=dp), allocatable, intent(inout) :: H_o(:, :, :)
    ! <u_mk|u_nk+b> m, n, nntot, k
    ! complex(kind=dp), allocatable, intent(inout) :: orb_o(:, :, :, :, :)
    complex(kind=dp), allocatable :: S_o(:, :)
    real(kind=dp) :: c_real, c_imag
    integer :: tmp_bands, tmp_kpts, tmp_nntot
    integer :: ik, ik2, nnl, nnm, nnn, nn, inn, m, n, nn1, nn2
    integer :: ncount
    logical :: nn_found

    write(stdout, '(1x,a)', advance='no') "Reading <u_k|u_k+b> from "//trim(seedname)//".mmn ..."
    open(newunit=iun_mmn, file=trim(seedname)//".mmn", &
        form='formatted', status='old', action='read')
    read(iun_mmn, *) header
    read(iun_mmn, '(3i12)') tmp_bands, tmp_kpts, tmp_nntot
    if (tmp_bands .ne. num_bands)then
      call set_error_input(error, 'Error: bands from mmn file dont match to win', comm)
      return
    endif
    if (tmp_kpts .ne. num_kpts) then
      call set_error_input(error, 'Error: kpts from mmn file dont match to win', comm)
      return
    endif
    if (tmp_nntot .ne. kmesh_info%nntot) then
      call set_error_input(error, 'Error: nntot from mmn file dont match to win', comm)
      return
    endif
    if (allocated(mmn)) then
      call set_error_input(error, 'Error: mmn matrix has been allocated before read', comm)
      return
    endif
    allocate(mmn(num_bands, num_bands, kmesh_info%nntot, num_kpts), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating mmn in get_mmn', comm)
      return
    endif
    allocate(S_o(num_bands, num_bands))
    
    mmn = cmplx_0
    ! allocate(S(num_wann, num_wann), stat=ierr)
    ! if (ierr /= 0) then
    !   call set_error_alloc(error, 'Error in allocating S in get_mmn', comm)
    !   return
    ! endif

    allocate(H_o(num_bands, num_bands, num_kpts), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating H_o in get_mmn', comm)
      return
    endif
    H_o = cmplx_0

    do ik = 1, num_kpts
      do m = 1, num_bands
        H_o(m, m, ik) = cmplx_1 * eigval(m, ik)
      enddo
    enddo

    do ncount = 1, num_kpts*kmesh_info%nntot
      !
      !Read from .mmn file the original overlap matrix
      ! S_o=<u_ik|u_ik2> between ab initio eigenstates
      !
      S_o = cmplx_0
      
      read (iun_mmn, *) ik, ik2, nnl, nnm, nnn
      do n = 1, num_bands
        do m = 1, num_bands
          read (iun_mmn, *) c_real, c_imag
          S_o(m, n) = cmplx(c_real, c_imag, kind=dp)
        enddo
      enddo
      nn = 0
      nn_found = .false.
      do inn = 1, kmesh_info%nntot
        if ((ik2 .eq. kmesh_info%nnlist(ik, inn)) .and. &
            (nnl .eq. kmesh_info%nncell(1, ik, inn)) .and. &
            (nnm .eq. kmesh_info%nncell(2, ik, inn)) .and. &
            (nnn .eq. kmesh_info%nncell(3, ik, inn))) then
          if (.not. nn_found) then
            nn_found = .true.
            nn = inn
          else
            call set_error_fatal(error, 'Error reading '//trim(seedname)//'.mmn.&
                  & More than one matching nearest neighbour found', comm)
            return
          endif
        endif
      end do
      if (nn .eq. 0) then
        write (stdout, '(/a,i8,2i5,i4,2x,3i3)') ' Error reading '//trim(seedname)//'.mmn:', &
          ncount, ik, ik2, nn, nnl, nnm, nnn
        call set_error_fatal(error, 'Neighbour not found', comm)
        return
      end if
      mmn(:, :, nn, ik) = S_o(:, :)

    enddo ! ncount over num_kpts * nntot
    
    close(iun_mmn)
    ! reading done
    write(stdout, '(1x,a/)') "Done"

    ! i <\nabla u_m | u_t >< u_t | \nabla H | u_n >
    !~i <\nabla u_m | u_t >(< u_t | Hk+b | u_nk+b > - < u_t | Hk | u_nk+b >)
    !=i <\nabla u_m | u_t >(M*H_k+b - H_k*M)
    !=i M^+ * (M*H_k+b - H_k*M)
    ! do ik = 1, num_kpts
    !   do nn2 = 1, kmesh_info%nntot
    !     S_o = mmn(:, :, nn2, ik)
    !     do nn1 = 1, kmesh_info%nntot
    !       orb_o(:, :, nn1, nn2, ik) = cmplx_i * matmul( &
    !         conjg(transpose(mmn(:, :, nn1, ik))), &
    !         ! (matmul(S_o, H_o(:, :, kmesh_info%nnlist(ik, nn2))) - &
    !         (matmul(S_o, H_o(:, :, ik)) - &
    !          matmul(H_o(:, :, ik), S_o)) &
    !       )
    !     enddo
    !   enddo
    ! enddo
    deallocate(S_o)
    ! deallocate(mmn)
  end subroutine get_mmn

  subroutine get_uHu(stdout, seedname, num_bands, num_kpts, num_wann, kmesh_info, formatted, &
                     eigval, uhu, error, comm)
    !================================================!
    !
    !! read uHu matrix from seedname.uHu
    !
    !================================================!
    ! only run on root node
    implicit none

    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error
    type(kmesh_info_type), intent(in) :: kmesh_info
    
    character(len=60) :: header
    integer, intent(in) :: stdout
    character(len=50), intent(inout) :: seedname
    integer, intent(in) :: num_bands, num_kpts, num_wann
    real(kind=dp), intent(in) :: eigval(:, :)
    logical, intent(in) :: formatted
    complex(kind=dp), allocatable, intent(inout) :: uhu(:, :, :, :, :)
    ! <\nabla u_m| H |\nabla u_n> m, n, nntot, nntot, ik
    ! complex(kind=dp), allocatable, intent(inout) :: orb_o(:, :, :, :, :)
    complex(kind=dp), allocatable :: Ho_qb1_q_qb2(:, :)
    real(kind=dp) :: c_real, c_imag
    integer :: tmp_bands, tmp_kpts, tmp_nntot
    integer :: ik, nn1, nn2, m, n

    write(stdout, '(1x,a)', advance='no') "Reading <u_k+b1|H|u_k+b2> from "//trim(seedname)//".uHu ..."
    if (formatted) then
      open(newunit=iun_uHu, file=trim(seedname)//".uHu", &
          form='formatted', status='old', action='read')
      read(iun_uHu, *) header
      read(iun_uHu, *) tmp_bands, tmp_kpts, tmp_nntot
    else ! .not. formatted
      open(newunit=iun_uHu, file=trim(seedname)//".uHu", &
          form='unformatted', status='old', action='read')
      read(iun_uHu) header
      read(iun_uHu) tmp_bands, tmp_kpts, tmp_nntot
    endif
    if (tmp_bands .ne. num_bands)then
      call set_error_input(error, 'Error: bands from uHu file dont match to win', comm)
      return
    endif
    if (tmp_kpts .ne. num_kpts) then
      call set_error_input(error, 'Error: kpts from uHu file dont match to win', comm)
      return
    endif
    if (tmp_nntot .ne. kmesh_info%nntot) then
      call set_error_input(error, 'Error: nntot from uHu file dont match to win', comm)
      return
    endif
    if (allocated(uhu)) then
      call set_error_input(error, 'Error: uhu matrix has been allocated before read', comm)
      return
    endif
    allocate(uhu(num_bands, num_bands, kmesh_info%nntot, kmesh_info%nntot, num_kpts), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating uHu in get_uHu', comm)
      return
    endif
    allocate(Ho_qb1_q_qb2(num_bands, num_bands), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating Ho_qb1_q_qb2 in get_uHu', comm)
      return
    endif
    
    do ik = 1, num_kpts
      Ho_qb1_q_qb2 = cmplx_0
      do nn2 = 1, kmesh_info%nntot
        do nn1 = 1, kmesh_info%nntot
          if (formatted) then
            do m = 1, num_bands
              do n = 1, num_bands
                read (iun_uHu, *) c_real, c_imag
                Ho_qb1_q_qb2(n, m) = cmplx(c_real, c_imag, kind=dp)
              end do
            end do
          else
            read (iun_uHu) &
              ((Ho_qb1_q_qb2(n, m), n=1, num_bands), m=1, num_bands)
          endif
          Ho_qb1_q_qb2 = transpose(Ho_qb1_q_qb2)
          uhu(:, :, nn1, nn2, ik) = Ho_qb1_q_qb2(:, :)
          ! orb_o(:, :, nn1, nn2, ik) = Ho_qb1_q_qb2
          ! orb_o(:, :, nn1, nn2, ik) = orb_o(:, :, nn1, nn2, ik) + &
          !                             cmplx_i * Ho_qb1_q_qb2
        enddo ! nn1
      enddo ! nn2
    enddo ! ik
    close(iun_uHu)
    write(stdout, '(1x,a/)') "Done"
    deallocate(Ho_qb1_q_qb2)
    ! deallocate(uhu)
  end subroutine get_uHu

  subroutine get_uIu(stdout, seedname, num_bands, num_kpts, num_wann, kmesh_info, formatted, &
                     eigval, uiu, error, comm)
    !================================================!
    !
    !! read uIu matrix from seedname.uIu
    !! Actually get < del u_n | E_m | del u_m >
    !
    !================================================!
    ! only run on root node
    implicit none

    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error
    type(kmesh_info_type), intent(in) :: kmesh_info
    
    character(len=60) :: header
    integer, intent(in) :: stdout
    character(len=50), intent(inout) :: seedname
    integer, intent(in) :: num_bands, num_kpts, num_wann
    real(kind=dp), intent(in) :: eigval(:, :)
    logical, intent(in) :: formatted
    complex(kind=dp), allocatable, intent(inout) :: uiu(:, :, :, :, :)
    ! <\nabla u_m| H \nabla|u_n> m, n, idir1, idir2, ik
    ! complex(kind=dp), allocatable, intent(inout) :: orb_o(:, :, :, :, :)
    complex(kind=dp), allocatable :: Lo_qb1_q_qb2(:, :)
    real(kind=dp) :: c_real, c_imag, temp
    integer :: tmp_bands, tmp_kpts, tmp_nntot
    integer :: ik, nn1, nn2, m, n, qb1, qb2

    write(stdout, '(1x,a)', advance='no') "Reading <u_k+b1|u_k+b2> from "//trim(seedname)//".uIu ..."
    if (formatted) then
      open(newunit=iun_uIu, file=trim(seedname)//".uIu", &
          form='formatted', status='old', action='read')
      read(iun_uIu, *) header
      read(iun_uIu, *) tmp_bands, tmp_kpts, tmp_nntot
    else ! .not. formatted
      open(newunit=iun_uIu, file=trim(seedname)//".uIu", &
          form='unformatted', status='old', action='read')
      read(iun_uIu) header
      read(iun_uIu) tmp_bands, tmp_kpts, tmp_nntot
    endif
    if (tmp_bands .ne. num_bands)then
      call set_error_input(error, 'Error: bands from uIu file dont match to win', comm)
      return
    endif
    if (tmp_kpts .ne. num_kpts) then
      call set_error_input(error, 'Error: kpts from uIu file dont match to win', comm)
      return
    endif
    if (tmp_nntot .ne. kmesh_info%nntot) then
      call set_error_input(error, 'Error: nntot from uIu file dont match to win', comm)
      return
    endif
    if (allocated(uiu)) then
      call set_error_input(error, 'Error: uiu matrix has been allocated before read', comm)
      return
    endif
    allocate(uiu(num_bands, num_bands, kmesh_info%nntot, kmesh_info%nntot, num_kpts), stat=ierr)
    uiu = cmplx_0
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating uIu in get_uIu', comm)
      return
    endif
    allocate(Lo_qb1_q_qb2(num_bands, num_bands), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating uIu and Lo_qb1_q_qb2 in get_uIu', comm)
      return
    endif
    
    do ik = 1, num_kpts
      Lo_qb1_q_qb2 = cmplx_0
      do nn2 = 1, kmesh_info%nntot
        ! qb2 = kmesh_info%nnlist(ik, nn2)
        do nn1 = 1, kmesh_info%nntot
          ! qb1 = kmesh_info%nnlist(ik, nn1)
          if (formatted) then
            do m = 1, num_bands
              do n = 1, num_bands
                read (iun_uIu, *) c_real, c_imag
                Lo_qb1_q_qb2(n, m) = cmplx(c_real, c_imag, kind=dp)
              end do
            end do
          else
            read (iun_uIu) &
              ((Lo_qb1_q_qb2(n, m), n=1, num_bands), m=1, num_bands)
          endif
          Lo_qb1_q_qb2 = transpose(Lo_qb1_q_qb2)
          uiu(:, :, nn1, nn2, ik) = Lo_qb1_q_qb2(:, :)
          ! do m = 1, num_bands
          !   do n = 1, num_bands
          !     ! Lo_qb1_q_qb2(n, m) = eigval(m, qb2) * Lo_qb1_q_qb2(n, m)
          !     Lo_qb1_q_qb2(n, m) = eigval(m, ik) * Lo_qb1_q_qb2(n, m)
          !   enddo
          ! enddo
          ! orb_o(:, :, nn1, nn2, ik) = orb_o(:, :, nn1, nn2, ik) - cmplx_i * Lo_qb1_q_qb2(:, :)
        enddo ! nn1
      enddo ! nn2
    enddo ! ik
    close(iun_uIu)
    write(stdout, '(1x,a/)') "Done"
    deallocate(Lo_qb1_q_qb2)
    ! deallocate(uiu)
  end subroutine get_uIu


  ! subroutine calc_orb(stdout, seedname, num_bands, num_kpts, num_wann, kmesh_info, dis_manifold, have_disentangled, &
  !                     eigval, v_matrix, orb_o, mmn, uhu, uiu, method, error, comm)
  !   !================================================!
  !   !
  !   !! calculate orbital matrix
  !   !
  !   !================================================!
  !   ! only run on root node
  !   use w90_io, only: io_date

  !   implicit none

  !   type(w90_comm_type), intent(in) :: comm
  !   type(w90_error_type), allocatable, intent(out) :: error
  !   type(kmesh_info_type), intent(in) :: kmesh_info
  !   type(dis_manifold_type), intent(in) :: dis_manifold
    
  !   integer, intent(in) :: stdout
  !   character(len=50), intent(inout) :: seedname
  !   integer, intent(in) :: num_bands, num_kpts, num_wann
  !   logical, intent(in) :: have_disentangled
  !   logical, intent(in) :: method
  !   real(kind=dp), intent(in) :: eigval(:, :)
  !   complex(kind=dp), allocatable, intent(in) :: v_matrix(:, :, :)
  !   complex(kind=dp), allocatable, intent(in) :: orb_o(:, :, :, :, :)
  !   complex(kind=dp), allocatable, intent(inout) :: mmn(:, :, :, :)
  !   complex(kind=dp), allocatable, intent(inout) :: uiu(:, :, :, :, :)
  !   ! <\nabla u_m| \nabla u_n> m, n, idir1, idir2, ik
  !   complex(kind=dp), allocatable, intent(inout) :: uhu(:, :, :, :, :)
  !   ! <\nabla u_m| H |\nabla u_n> m, n, idir1, idir2, ik
  !   complex(kind=dp), allocatable :: orb_ab(:, :, :, :), orb_ab_k(:, :)
  !   complex(kind=dp), allocatable :: orb_ab_h(:, :, :, :), orb_h_w_k(:, :, :), orb_h_w2o_k(:, :, :)
  !   complex(kind=dp), allocatable :: orb_ab_uiu(:, :, :, :, :), orb_uiu_w_k(:, :, :, :), orb_uiu_w2o_k(:, :, :, :)
  !   complex(kind=dp), allocatable :: orb_ab_uhu(:, :, :, :, :), orb_uhu_w_k(:, :, :, :), orb_uhu_w2o_k(:, :, :, :)
  !   complex(kind=dp), allocatable :: orb_w_k(:, :, :), orb_w2o_k(:, :, :)
  !   complex(kind=dp), allocatable :: vdag_orb_k(:, :)
  !   integer, allocatable :: num_states(:)
  !   integer :: ik, m, n, idir1, idir2, nn1, nn2, nn3, qb1, qb2

  !   integer, dimension(3), parameter :: alpha_A = (/2, 3, 1/)
  !   integer, dimension(3), parameter :: beta_A = (/3, 1, 2/)
  !   real(kind=dp), parameter :: fac = 3.674932379e-2_dp / (0.52917721092_dp)**2

  !   character(len=60) header
  !   character(len=9) cdate, ctime
  !   call io_date(cdate, ctime)
  !   header = 'Created on ' // cdate // ' at ' // ctime
    
  !   allocate(num_states(num_kpts))
  !   do ik = 1, num_kpts
  !     if (have_disentangled) then
  !       num_states(ik) = dis_manifold%ndimwin(ik)
  !     else
  !       num_states(ik) = num_wann
  !     endif
  !   enddo

  !   open (newunit=iun_orb, file=trim(seedname)//'.orb.fmt', form='formatted', status='replace', position='rewind')
  !   write(iun_orb, *) header
  !   write(iun_orb, *) num_bands, num_kpts
  !   !
  !   if (.not. method) then
  !     open (newunit=iun_sIu, file=trim(seedname)//'.sIu', form='unformatted', status='replace', position='rewind')
  !     write(iun_sIu) header
  !     write(iun_sIu) num_bands, num_kpts, kmesh_info%nntot

  !     open (newunit=iun_sHu, file=trim(seedname)//'.sHu', form='unformatted', status='replace', position='rewind')
  !     write(iun_sHu) header
  !     write(iun_sHu) num_bands, num_kpts, kmesh_info%nntot

  !     open (newunit=iun_sH, file=trim(seedname)//'.sH', form='unformatted', status='replace', position='rewind')
  !     write(iun_sH) header
  !     write(iun_sH) num_bands, num_kpts
  !   endif
  !   !
  !   if (allocated(orb_w_k)) then
  !     call set_error_input(error, 'Error: orb matrix has been allocated before calculated', comm)
  !     return
  !   endif
  !   allocate(orb_w_k(num_wann, num_wann, 3), orb_w2o_k(num_bands, num_bands, 3), stat=ierr)
  !   if (ierr /= 0) then
  !     call set_error_alloc(error, 'Error in allocating orb_w_k or orb_w2o_k in calc_orb', comm)
  !     return
  !   endif
  !   allocate(orb_ab(num_wann, num_wann, 3, 3), orb_ab_k(num_wann, num_wann), stat=ierr)
  !   if (ierr /= 0) then
  !     call set_error_alloc(error, 'Error in allocating orb_ab or orb_ab_k in calc_orb', comm)
  !     return
  !   endif

  !   if (.not. method) then
  !     allocate(orb_ab_uiu(num_wann, num_bands, 3, 3, kmesh_info%nntot), &
  !              orb_uiu_w_k(num_wann, num_bands, 3, kmesh_info%nntot), &
  !              orb_uiu_w2o_k(num_bands, num_bands, 3, kmesh_info%nntot), stat=ierr)
  !     allocate(orb_ab_uhu(num_wann, num_bands, 3, 3, kmesh_info%nntot), &
  !              orb_uhu_w_k(num_wann, num_bands, 3, kmesh_info%nntot), &
  !              orb_uhu_w2o_k(num_bands, num_bands, 3, kmesh_info%nntot), stat=ierr)
  !     allocate(orb_ab_h(num_wann, num_bands, 3, 3), &
  !              orb_h_w_k(num_wann, num_bands, 3), &
  !              orb_h_w2o_k(num_bands, num_bands, 3), stat=ierr)
  !     allocate(vdag_orb_k(num_wann, num_bands), stat=ierr)
  !     orb_ab_uhu = cmplx_0
  !     orb_ab_uiu = cmplx_0
  !     orb_ab_h = cmplx_0
  !   endif
  !   orb_w_k = cmplx_0

  !   do ik = 1, num_kpts
  !     orb_ab = cmplx_0
  !     if (.not. method) then
  !       orb_ab_uhu = cmplx_0
  !       orb_ab_uiu = cmplx_0
  !       orb_ab_h = cmplx_0
  !       do n = 1, num_bands
  !         mmn(n, :, :, ik) = eigval(n, ik) * mmn(n, :, :, ik)
  !       enddo
  !     endif
  !     do nn2 = 1, kmesh_info%nntot
  !       qb2 = kmesh_info%nnlist(ik, nn2)
  !       do nn1 = 1, kmesh_info%nntot
  !         qb1 = kmesh_info%nnlist(ik, nn1)
  !         orb_ab_k = cmplx_0
  !         call get_gauge_overlap_matrix(num_bands, num_wann, eigval, v_matrix, dis_manifold, &
  !                                       qb1, num_states(qb1), qb2, num_states(qb2), &
  !                                       orb_o(:, :, nn1, nn2, ik), have_disentangled, orb_ab_k)
  !         if (.not. method) then
  !           vdag_orb_k = matmul(conjg(transpose(v_matrix(:, :, qb1))), orb_o(:, :, nn1, nn2, ik))
  !         endif
  !         do idir2 = 1, 3
  !           do idir1 = 1, 3
  !             orb_ab(:, :, idir1, idir2) = orb_ab(:, :, idir1, idir2) + &
  !                                           kmesh_info%wb(nn1)*kmesh_info%bk(idir1, nn1, ik)* &
  !                                           kmesh_info%wb(nn2)*kmesh_info%bk(idir2, nn2, ik)*orb_ab_k(:, :)
  !             if (.not. method) then
  !               orb_ab_h(:, :, idir1, idir2) = orb_ab_h(:, :, idir1, idir2) + &
  !                                              kmesh_info%wb(nn1)*kmesh_info%bk(idir1, nn1, ik)* &
  !                                              kmesh_info%wb(nn2)*kmesh_info%bk(idir2, nn2, ik)* &
  !                                              matmul(vdag_orb_k, conjg(transpose(mmn(:, :, nn2, ik))))
  !               do nn3 = 1, kmesh_info%nntot
  !                 orb_ab_uiu(:, :, idir1, idir2, nn3) = orb_ab_uiu(:, :, idir1, idir2, nn3) + &
  !                                                       kmesh_info%wb(nn1)*kmesh_info%bk(idir1, nn1, ik)* &
  !                                                       kmesh_info%wb(nn2)*kmesh_info%bk(idir2, nn2, ik)* &
  !                                                       matmul(vdag_orb_k, uiu(:, :, nn2, nn3, ik))
  !                 orb_ab_uhu(:, :, idir1, idir2, nn3) = orb_ab_uhu(:, :, idir1, idir2, nn3) + &
  !                                                       kmesh_info%wb(nn1)*kmesh_info%bk(idir1, nn1, ik)* &
  !                                                       kmesh_info%wb(nn2)*kmesh_info%bk(idir2, nn2, ik)* &
  !                                                       matmul(vdag_orb_k, uhu(:, :, nn2, nn3, ik))
  !               enddo
  !             endif ! method
  !           enddo ! idir1
  !         enddo ! idir2
  !       enddo ! nn1
  !     enddo ! nn2
  !     orb_w2o_k = cmplx_0
  !     do idir1 = 1, 3
  !       orb_w_k(:, :, idir1) = orb_ab(:, :, alpha_A(idir1), beta_A(idir1)) - orb_ab(:, :, beta_A(idir1), alpha_A(idir1))
  !       orb_w_k(:, :, idir1) = fac * (orb_w_k(:, :, idir1) + conjg(transpose(orb_w_k(:, :, idir1))))
  !       orb_w2o_k(:, :, idir1) = matmul(v_matrix(:, :, ik), matmul(orb_w_k(:, :, idir1), conjg(transpose(v_matrix(:, :, ik)))))
  !       !
  !       if (.not. method) then
  !         orb_h_w_k(:, :, idir1) = orb_ab_h(:, :, alpha_A(idir1), beta_A(idir1)) - &
  !                                  orb_ab_h(:, :, beta_A(idir1), alpha_A(idir1))
  !         orb_h_w2o_k(:, :, idir1) = 2 * fac * matmul(v_matrix(:, :, ik), orb_h_w_k(:, :, idir1))
  !         do nn3 = 1, kmesh_info%nntot
  !           orb_uiu_w_k(:, :, idir1, nn3) = orb_ab_uiu(:, :, alpha_A(idir1), beta_A(idir1), nn3) - &
  !                                           orb_ab_uiu(:, :, beta_A(idir1), alpha_A(idir1), nn3)
  !           orb_uiu_w2o_k(:, :, idir1, nn3) = 2 * fac * matmul(v_matrix(:, :, ik), orb_uiu_w_k(:, :, idir1, nn3))
  !           orb_uhu_w_k(:, :, idir1, nn3) = orb_ab_uhu(:, :, alpha_A(idir1), beta_A(idir1), nn3) - &
  !                                           orb_ab_uhu(:, :, beta_A(idir1), alpha_A(idir1), nn3)
  !           orb_uhu_w2o_k(:, :, idir1, nn3) = 2 * fac * matmul(v_matrix(:, :, ik), orb_uhu_w_k(:, :, idir1, nn3))
  !         enddo
  !       endif
  !       !
  !     enddo ! idir1

  !     do m = 1, num_bands
  !       do n = 1, m
  !         do idir1 = 1, 3
  !           write(iun_orb, '(2es26.16)') orb_w2o_k(n, m, idir1)
  !         enddo ! idir1
  !       enddo
  !     enddo
  !     !
  !     if (.not. method) then
  !       do idir1 = 1, 3
  !         write(iun_sH) ((orb_h_w2o_k(n, m, idir1), n=1, num_bands), m=1, num_bands)
  !       enddo
  !       do nn3 = 1, kmesh_info%nntot
  !         do idir1 = 1, 3
  !           ! transpose here
  !           write(iun_sIu) ((orb_uiu_w2o_k(n, m, idir1, nn3), m=1, num_bands), n=1, num_bands)
  !           write(iun_sHu) ((orb_uiu_w2o_k(n, m, idir1, nn3), m=1, num_bands), n=1, num_bands)
  !         enddo
  !       enddo
  !     endif
  !     !
  !   enddo ! ik
  !   deallocate(orb_w_k, orb_w2o_k, orb_ab, orb_ab_k)
  !   if (.not. method) deallocate(orb_uiu_w2o_k, orb_uiu_w_k, orb_ab_uiu)
  !   if (.not. method) deallocate(orb_uhu_w2o_k, orb_uhu_w_k, orb_ab_uhu)
  !   if (.not. method) deallocate(orb_h_w2o_k, orb_h_w_k, orb_ab_h)
  !   close(iun_orb)
  !   if (.not. method) close(iun_sIu)
  !   if (.not. method) close(iun_sHu)
  !   if (.not. method) close(iun_sH)
  !   deallocate(uiu, uhu, mmn)
  ! end subroutine

    subroutine calc_orb_gh(stdout, seedname, num_bands, num_kpts, num_wann, kmesh_info, dis_manifold, have_disentangled, &
                           eigval, v_matrix, H_o, mmn, uhu, uiu, wan_gauge, error, comm)
    !================================================!
    !
    !! calculate orbital matrix
    !
    !================================================!
    ! only run on root node
    use w90_io, only: io_date

    implicit none

    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error
    type(kmesh_info_type), intent(in) :: kmesh_info
    type(dis_manifold_type), intent(in) :: dis_manifold
    
    integer, intent(in) :: stdout
    character(len=50), intent(inout) :: seedname
    integer, intent(in) :: num_bands, num_kpts, num_wann
    logical, intent(in) :: have_disentangled
    logical, intent(in) :: wan_gauge
    real(kind=dp), intent(in) :: eigval(:, :)
    complex(kind=dp), allocatable, intent(inout) :: v_matrix(:, :, :)
    complex(kind=dp), allocatable, intent(inout) :: H_o(:, :, :)
    complex(kind=dp), allocatable, intent(inout) :: mmn(:, :, :, :)
    complex(kind=dp), allocatable, intent(inout) :: uiu(:, :, :, :, :)
    ! <\nabla u_m| \nabla u_n> m, n, idir1, idir2, ik
    complex(kind=dp), allocatable, intent(inout) :: uhu(:, :, :, :, :)
    ! <\nabla u_m| H |\nabla u_n> m, n, idir1, idir2, ik
    complex(kind=dp), allocatable :: orb_g(:, :), orb_h(:, :)
    ! m, n
    complex(kind=dp), allocatable :: orb_ab(:, :, :, :)
    ! convert orb_g and orb_h to orb_ab where a,b = x, y, z
    ! m, n
    complex(kind=dp), allocatable :: orb_o(:, :, :)
    ! cross product orb_ab into orb_o
    ! m, n, idir1

    ! temp variables below here
    complex(kind=dp), allocatable :: mmn_b1(:, :), mmn_b2(:, :)
    complex(kind=dp), allocatable :: VVd(:, :)
    integer :: ik, m, n, idir1, idir2, nn1, nn2, qb1, qb2

    integer, dimension(3), parameter :: alpha_A = (/2, 3, 1/)
    integer, dimension(3), parameter :: beta_A = (/3, 1, 2/)
    real(kind=dp), parameter :: fac = 3.674932379e-2_dp / (0.52917721092_dp)**2

    character(len=60) header
    character(len=9) cdate, ctime


    call io_date(cdate, ctime)
    header = 'Created on ' // cdate // ' at ' // ctime

    open (newunit=iun_orb, file=trim(seedname)//'.orb.fmt', form='formatted', status='replace', position='rewind')
    write(iun_orb, *) header
    write(iun_orb, *) num_bands, num_kpts

    allocate(mmn_b1(num_bands, num_bands),mmn_b2(num_bands, num_bands))
    mmn_b1 = cmplx_0
    mmn_b2 = cmplx_0

    allocate(orb_g(num_bands, num_bands), &
             orb_h(num_bands, num_bands), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating orb_g or orb_h in calc_orb_gh', comm)
      return
    endif
    orb_g = cmplx_0
    orb_h = cmplx_0

    allocate(orb_ab(num_bands, num_bands, 3, 3), &
             orb_o(num_bands, num_bands, 3), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating orb_ab or orb_o in calc_orb_gh', comm)
      return
    endif

    allocate(VVd(num_bands, num_bands), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating Vdagger * V in calc_orb_gh', comm)
      return
    endif
    VVd = cmplx_0

    do ik = 1, num_kpts
      orb_g = cmplx_0
      orb_h = cmplx_0
      orb_ab = cmplx_0
      ! if (wan_gauge) VVd(:, :) = matmul(v_matrix(:, :, ik), conjg(transpose(v_matrix(:, :, ik))))
      do nn2 = 1, kmesh_info%nntot
        mmn_b2(:, :) = mmn(:, :, nn2, ik)
        do nn1 = 1, kmesh_info%nntot
          if (wan_gauge) VVd(:, :) = matmul(v_matrix(:, :, kmesh_info%nnlist(ik, nn1)), &
                            conjg(transpose(v_matrix(:, :, kmesh_info%nnlist(ik, nn1)))))
          mmn_b1(:, :) = mmn(:, :, nn1, ik)
          ! <k | k+b1> [<k+b1 | H | k+b2> - <k+b1 | k> <k | H | k> <k | k+b2>] <k+b2 | k>
          if (wan_gauge) then
            orb_g(:, :) = matmul(VVd, uhu(:, :, nn1, nn2, ik)) !- &
                          ! matmul(matmul(conjg(transpose(mmn_b1)), VVd(:, :)), &
                          ! matmul(H_o(:, :, ik), mmn_b2))
          else
            orb_g(:, :) = uhu(:, :, nn1, nn2, ik) - &
                          matmul(conjg(transpose(mmn_b1)), &
                          matmul(H_o(:, :, ik), mmn_b2))
          endif
          orb_g(:, :) = cmplx_i * matmul(mmn_b1, matmul(orb_g(:, :), conjg(transpose(mmn_b2))))

          ! <k | k+b1> [<k+b1 | k+b2> - <k+b1 | k> <k | k+b2>] <k+b2 | k> <k | H | k>
          if (wan_gauge) then
            orb_h(:, :) = matmul(VVd, matmul(uiu(:, :, nn1, nn2, ik), H_o(:, :, kmesh_info%nnlist(ik, nn2)))) ! - &
                          ! matmul(matmul(conjg(transpose(mmn_b1)), VVd(:, :)), mmn_b2)
          else
            orb_h(:, :) = uiu(:, :, nn1, nn2, ik) - &
                          matmul(conjg(transpose(mmn_b1)), mmn_b2)
            orb_h(:, :) = matmul(orb_h, H_o(:, :, ik))
          endif
          orb_h(:, :) = cmplx_i * matmul(mmn_b1, matmul(orb_h(:, :), conjg(transpose(mmn_b2))))
          ! orb_h(:, :) = matmul(orb_h, H_o(:, :, ik))
          do idir2 = 1, 3
            do idir1 = 1, 3
              orb_ab(:, :, idir1, idir2) = orb_ab(:, :, idir1, idir2) + &
                                            kmesh_info%wb(nn1)*kmesh_info%bk(idir1, nn1, ik)* &
                                            kmesh_info%wb(nn2)*kmesh_info%bk(idir2, nn2, ik)* &
                                            (orb_g(:, :) - orb_h(:, :))
            enddo ! idir1
          enddo ! idir2
        enddo ! nn1
      enddo ! nn2
      orb_o = cmplx_0
      do idir1 = 1, 3
        orb_o(:, :, idir1) = orb_ab(:, :, alpha_A(idir1), beta_A(idir1)) - orb_ab(:, :, beta_A(idir1), alpha_A(idir1))
        orb_o(:, :, idir1) = fac * (orb_o(:, :, idir1) + conjg(transpose(orb_o(:, :, idir1))))
        !
      enddo ! idir1

      do m = 1, num_bands
        do n = 1, m
          do idir1 = 1, 3
            write(iun_orb, '(2es26.16)') orb_o(n, m, idir1)
          enddo ! idir1
        enddo
      enddo
    enddo ! ik
    deallocate(orb_g, orb_h, orb_ab, orb_o)
    close(iun_orb)
    deallocate(uiu, uhu, mmn, H_o)
  end subroutine

  ! subroutine output_orb_formatted(stdout, seedname, num_bands, num_kpts, num_wann, orb, v_matrix)
  !   !================================================!
  !   !
  !   !! write orbital matrix to seedname.orb.fmt
  !   !
  !   !================================================!
  !   ! only run on root node
  !   use w90_io, only: io_date

  !   implicit none
    
  !   integer, intent(in) :: stdout
  !   character(len=50), intent(inout) :: seedname
  !   integer, intent(in) :: num_bands, num_kpts, num_wann
  !   complex(kind=dp), allocatable, intent(inout) :: orb(:, :, :, :)
  !   ! m, n, idir, ik
  !   complex(kind=dp), allocatable, intent(inout) :: v_matrix(:, :, :)
  !   ! num_bands, num_wann, ik
  !   complex(kind=dp), allocatable :: orb_o_k(:, :, :)
  !   integer :: ik, m, n, idir
  !   character(len=60) header
  !   character(len=9) cdate, ctime
  !   call io_date(cdate, ctime)
  !   header = 'Created on ' // cdate // ' at ' // ctime

  !   write (stdout, '(1x,a)', advance='no') "Writing information to formatted file "//trim(seedname)//".orb.fmt ..."

  !   open (newunit=iun_orb, file=trim(seedname)//'.orb.fmt', form='formatted', status='replace', position='rewind')
  !   allocate(orb_o_k(num_bands, num_bands, 3))

  !   write(iun_orb, *) header
  !   write(iun_orb, *) num_bands, num_kpts
  !   do ik = 1, num_kpts
  !     orb_o_k = cmplx_0
  !     do idir = 1, 3
  !       orb_o_k(:, :, idir) = matmul(v_matrix(:, :, ik), matmul(orb(:, :, idir, ik), conjg(transpose(v_matrix(:, :, ik)))))
  !     enddo
  !     do m = 1, num_bands
  !       do n = 1, m
  !         do idir = 1, 3
  !           write(iun_orb, '(2es26.16)') orb_o_k(n, m, idir)
  !         enddo ! idir
  !       enddo
  !     enddo
  !   enddo
  !   close(iun_orb)
  !   write (stdout, '(1x,a/)') "Done"

  ! end subroutine
end module get_orb


program w90genorb
  !! Program to convert spn files from formatted to unformmated
  use w90_constants, only: dp, cmplx_0, pw90_physical_constants_type
  use get_orb, only: get_seedname, get_mmn, get_uHu, get_uIu, calc_orb_gh
  use w90_error
  use w90_io
  use w90_types
  use w90_postw90_types
  use w90_postw90_readwrite
  use w90_comms, only: w90_comm_type, mpirank, mpisize
  use w90_readwrite, only: w90_readwrite_write_header, w90_readwrite_in_file, &
    w90_readwrite_read_gamma_only, w90_readwrite_read_mp_grid, w90_readwrite_read_units, &
    w90_readwrite_read_kmesh_data, w90_readwrite_read_lattice, w90_readwrite_read_kpoints, &
    w90_readwrite_read_eigvals, w90_readwrite_read_num_bands, w90_readwrite_read_chkpt, &
    w90_readwrite_clean_infile, w90_readwrite_read_final_alloc
  use w90_postw90_common, only: pw90common_wanint_data_dist
  use w90_kmesh, only: kmesh_get, kmesh_sort
  use w90_utility, only: utility_recip_lattice
  use w90_error_base, only: w90_error_type
  
#ifdef MPI
#  if !(defined(MPI08) || defined(MPI90) || defined(MPIH))
#    error "You need to define which MPI interface you are using"
#  endif
#endif

#ifdef MPI08
  use mpi_f08 ! use f08 interface if possible
#endif
#ifdef MPI90
  use mpi ! next best, use fortran90 interface
#endif

  implicit none

#ifdef MPIH
  include 'mpif.h' ! worst case, use legacy interface
#endif 
  type(pw90_physical_constants_type) :: physics

  type(atom_data_type) :: atoms
  type(w90_system_type) :: system
  type(dis_manifold_type) :: dis_window
  type(wannier_data_type) :: wann_data
  type(ws_region_type) :: ws_region
  type(kmesh_info_type) :: kmesh_info
  type(kmesh_input_type) :: kmesh_data
  type(print_output_type) :: print_output
  type(settings_type) :: settings
  type(w90_comm_type) :: comm
  type(timer_list_type) :: timer
  type(w90_error_type), allocatable :: error 
  type(kpoint_path_type) :: spec_points
  
  real(kind=dp) :: bohr
  character(len=20) :: energy_unit
  integer :: num_bands, num_kpts, num_wann
  real(kind=dp) :: real_lattice(3, 3)
  real(kind=dp) :: recip_lattice(3, 3), volume
  real(kind=dp), allocatable :: kpt_latt(:, :)
  real(kind=dp), pointer :: eigval(:, :)
  ! iband, ik, idir
  integer :: mp_grid(3)
  integer :: optimisation
  real(kind=dp), allocatable :: fermi_energy_list(:)
  logical :: gamma_only, eig_found, wan_gauge
  logical :: formatted = .false.
  logical :: effective_model = .false.
  logical :: have_disentangled

  integer :: my_node_id, num_nodes, ierr
  logical :: on_root
  integer :: stdout, stderr
  character(len=50) :: seedname
  integer :: num_exclude_bands
  integer, allocatable :: exclude_bands(:)
  ! this is a dummy that is not used, DO NOT use
  complex(kind=dp), allocatable :: m_matrix(:, :, :, :)
  ! u_matrix_opt here only for generation of v_matrix
  ! u_matrix_opt gives the num_wann dimension optimal subspace from the
  ! original bloch states
  complex(kind=dp), allocatable :: u_matrix_opt(:, :, :)
  ! optimally smooth states.
  ! m_matrix we store here, becuase it is needed for restart of wannierise
  complex(kind=dp), allocatable :: u_matrix(:, :, :)
  real(kind=dp) :: scissors_shift
  complex(kind=dp), allocatable :: v_matrix(:, :, :)
  character(len=20) :: checkpoint
  real(kind=dp) :: omega_invariant
  type(kpoint_dist_type) :: kpt_dist
  type(pw90_band_deriv_degen_type) :: pw90_ham
  type(pw90_berry_mod_type) :: berry
  type(pw90_boltzwann_type) :: boltz
  type(pw90_dos_mod_type) :: dos_data
  type(pw90_extra_io_type) :: write_data
  type(pw90_geninterp_mod_type) :: geninterp
  type(pw90_gyrotropic_type) :: gyrotropic
  type(pw90_kpath_mod_type) :: kpath
  type(pw90_kslice_mod_type) :: kslice
  type(pw90_oper_read_type) :: postw90_oper
  type(pw90_spin_hall_type) :: spin_hall
  type(pw90_spin_mod_type) :: pw90_spin
  type(wigner_seitz_type) :: ws_vec
  type(ws_distance_type) :: ws_distance
  type(pw90_calculation_type) :: pw90_calcs

  complex(kind=dp), allocatable :: H_o(:, :, :)
  complex(kind=dp), allocatable :: mmn(:, :, :, :)
  ! <u_mk|u_nk+b> m, n, nntot, k
  complex(kind=dp), allocatable :: uhu(:, :, :, :, :)
  ! <\nabla u_m| H \nabla|u_n> m, n, idir1, idir2, ik
  complex(kind=dp), allocatable :: uiu(:, :, :, :, :)
  ! <\nabla u_m| \nabla|u_n> m, n, idir1, idir2, ik
  ! complex(kind=dp), allocatable :: temp(:, :, :)


#ifdef MPI
  comm%comm = MPI_COMM_WORLD
  call mpi_init(ierr)
  if (ierr .ne. 0) then
    call set_error_fatal(error, 'MPI initialisation error', comm)
    if (allocated(error)) call print_error_halt(error, ierr, stdout, stderr, comm)
  endif
#endif

  my_node_id = mpirank(comm)
  num_nodes = mpisize(comm)
  on_root = .false.
  if (my_node_id == 0) on_root = .true.

  stdout = 6
  stderr = 0

  if (on_root) then 
    call get_seedname(stdout, seedname, wan_gauge)
    open (newunit=stdout, file=trim(seedname)//'.log')
    write(stdout, "(a,i4,a)") "Running on", num_nodes, " nodes"
    call w90_readwrite_in_file(settings, seedname, error, comm)
    call w90_postw90_readwrite_read(settings, ws_region, system, exclude_bands, print_output, &
                                  kmesh_data, kpt_latt, num_kpts, dis_window, fermi_energy_list, &
                                  atoms, num_bands, num_wann, eigval, mp_grid, real_lattice, &
                                  spec_points, pw90_calcs, postw90_oper, scissors_shift, &
                                  effective_model, pw90_spin, pw90_ham, kpath, kslice, dos_data, &
                                  berry, spin_hall, gyrotropic, geninterp, boltz, eig_found, &
                                  write_data, gamma_only, physics%bohr, optimisation, stdout, &
                                  seedname, error, comm)
    if (allocated(error)) call print_error_halt(error, ierr, stdout, stderr, comm)
    call w90_readwrite_clean_infile(settings, stdout, seedname, error, comm)
    if (allocated(error)) call print_error_halt(error, ierr, stdout, stderr, comm)

    call w90_readwrite_read_final_alloc((num_bands > num_wann), dis_window, wann_data, num_wann, &
                                        num_bands, num_kpts, error, comm)
    if (allocated(error)) call print_error_halt(error, ierr, stdout, stderr, comm)

    call w90_readwrite_read_units(settings, print_output%lenconfac, print_output%length_unit, &
                                  energy_unit, bohr, error, comm)
    if (allocated(error)) call print_error_halt(error, ierr, stdout, stderr, comm)
    ! call w90_readwrite_read_num_bands(settings, effective_model, num_bands, num_wann, error, comm)

    ! write(stdout, *) num_bands, num_wann, num_kpts
    ! call w90_readwrite_read_mp_grid(settings, effective_model, mp_grid, num_kpts, error, comm)
    ! call w90_readwrite_read_gamma_only(settings, gamma_only, num_kpts, error, comm)

    ! call w90_readwrite_read_lattice(settings, real_lattice, bohr, error, comm)
    ! call w90_readwrite_read_kmesh_data(settings, kmesh_data, error, comm)
    ! call utility_recip_lattice(real_lattice, recip_lattice, volume, error, comm)
    ! call w90_readwrite_read_kpoints(settings, effective_model, kpt_latt, num_kpts, bohr, &
    !                                 error, comm)

    ! allocate (eigval(num_bands, num_kpts)) !fixme(jj) check allocation success
    ! call w90_readwrite_read_eigvals(eig_found, eigval, num_bands, num_kpts, stdout, seedname, &
    !                                 error, comm)

    call kmesh_get(kmesh_data, kmesh_info, print_output, kpt_latt, real_lattice, &
                 num_kpts, gamma_only, stdout, timer, error, comm)
    if (allocated(error)) call print_error_halt(error, ierr, stdout, stderr, comm)
    call kmesh_sort(kmesh_info, num_kpts, error, comm)
    if (allocated(error)) call print_error_halt(error, ierr, stdout, stderr, comm)

    if (.not. effective_model) then

      ! Read files seedname.chk (overlap matrices, unitary matrices for
      ! both disentanglement and maximal localization, etc.)

      !-----------------JJ
      !if (on_root) then
      allocate (u_matrix_opt(num_bands, num_wann, num_kpts), stat=ierr)
      allocate (u_matrix(num_bands, num_wann, num_kpts), stat=ierr)
      allocate (m_matrix(num_wann, num_wann, kmesh_info%nntot, num_kpts), stat=ierr)
      !else
      !  allocate (m_matrix(0, 0, 0, 0))
      !endif
      !m_matrix = cmplx_0
      !-----------------JJ

      !if (on_root) then
      num_exclude_bands = 0
      if (allocated(exclude_bands)) num_exclude_bands = size(exclude_bands)
      call w90_readwrite_read_chkpt(dis_window, exclude_bands, kmesh_info, kpt_latt, wann_data, &
                                    m_matrix, u_matrix, u_matrix_opt, real_lattice, &
                                    omega_invariant, mp_grid, num_bands, num_exclude_bands, &
                                    num_kpts, num_wann, checkpoint, have_disentangled, .true., &
                                    seedname, stdout, error, comm)
      if (allocated(error)) call print_error_halt(error, ierr, stdout, stderr, comm)
      !endif

      ! Distribute the information in the um and chk files to the other nodes
      !
      ! Ivo: For interpolation purposes we do not need u_matrix_opt and
      !      u_matrix separately, only their product v_matrix, and this
      !      is what is distributed now
      !
      call pw90common_wanint_data_dist(num_wann, num_kpts, num_bands, u_matrix_opt, u_matrix, &
                                      dis_window, wann_data, scissors_shift, v_matrix, &
                                      system%num_valence_bands, have_disentangled, error, comm)
      if (allocated(error)) call print_error_halt(error, ierr, stdout, stderr, comm)

    end if
    ! call calc_del(stdout, num_bands, num_wann, num_kpts, v_matrix, dv, eigval, del_eig, kmesh_info, error, comm)

    call io_stopwatch_start('process mmn matrix', timer)
    call get_mmn(stdout, seedname, num_bands, num_kpts, num_wann, kmesh_info, &
                 eigval, mmn, H_o ,error, comm)
    if (allocated(error)) call print_error_halt(error, ierr, stdout, stderr, comm)
    call io_stopwatch_stop('process mmn matrix', timer)

    call io_stopwatch_start('process uHu matrix', timer)
    call get_uHu(stdout, seedname, num_bands, num_kpts, num_wann, kmesh_info, formatted, &
                 eigval, uhu ,error, comm)
    if (allocated(error)) call print_error_halt(error, ierr, stdout, stderr, comm)
    call io_stopwatch_stop('process uHu matrix', timer)

    call io_stopwatch_start('process uIu matrix', timer)
    call get_uIu(stdout, seedname, num_bands, num_kpts, num_wann, kmesh_info, formatted, &
                 eigval, uiu ,error, comm)
    if (allocated(error)) call print_error_halt(error, ierr, stdout, stderr, comm)
    call io_stopwatch_stop('process uIu matrix', timer)

    ! allocate(temp(num_bands, num_bands, 3))

    call io_stopwatch_start('calculate orb matrix', timer)
    call calc_orb_gh(stdout, seedname, num_bands, num_kpts, num_wann, kmesh_info, dis_window, have_disentangled,&
                     eigval, v_matrix, H_o, mmn, uhu, uiu, wan_gauge, error, comm)
    ! call output_orb_formatted(stdout, seedname, num_bands, num_kpts, num_wann, orb, v_matrix)
    call io_stopwatch_stop('calculate orb matrix', timer)
    call io_print_timings(timer, stdout)
    write(stdout, '(/1x,a)') "Exit"

    close (unit=stdout)
  endif

end program w90genorb
