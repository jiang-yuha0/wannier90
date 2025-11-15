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

  integer, save :: iun_mmn, iun_uHu, iun_uIu, iun_orb
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
    write (stdout, '(A)') "  -f"
    write (stdout, '(A)') "      The input uIu and uHu matrix are formatted."
    write (stdout, '(A)') "  -u"
    write (stdout, '(A)') "      The input uIu and uHu matrix are unformatted."
  end subroutine print_usage

  !================================================!
  subroutine get_seedname(stdout, seedname, formatted)
    !================================================!
    !
    !! Set the seedname from the command line
    !
    !================================================!
    implicit none

    integer, intent(in) :: stdout
    character(len=50), intent(inout)  :: seedname
    logical, intent(inout) :: formatted

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
    if (index(ctemp, '-f') > 0) then
      formatted = .true.
    elseif (index(ctemp, '-u') > 0) then
      formatted = .false.
    else
      write (stdout, '(A)') 'Wrong command line action: '//trim(ctemp)
      call print_usage(stdout)
      call io_error('Wrong command line arguments, see logfile for usage', stdout)
    end if

  end subroutine get_seedname

  subroutine get_mmn(stdout, seedname, num_bands, num_kpts, num_wann, kmesh_info, dis_manifold, have_disentangled, &
                     v_matrix, eigval, mmn, error, comm, hmmn, mhmn, dhmmn, del_H)
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
    type(dis_manifold_type), intent(in)   :: dis_manifold
    
    character(len=60) :: header
    integer, intent(in) :: stdout
    character(len=50), intent(inout) :: seedname
    integer, intent(in) :: num_bands, num_kpts, num_wann
    logical, intent(in) :: have_disentangled
    complex(kind=dp), intent(in) :: v_matrix(:, :, :)
    real(kind=dp), intent(in) :: eigval(:, :)
    complex(kind=dp), allocatable, intent(inout) :: mmn(:, :, :, :)
    ! <u_m|\nabla|u_n> m, n, idir, ik
    complex(kind=dp), allocatable :: S_o(:, :), S(:, :)
    complex(kind=dp), allocatable :: H_o(:, :, :), MH_o(:, :), HM_o(:, :)
    complex(kind=dp), allocatable :: del_H_o(:, :, :, :), del_HM_o(:, :), del_HM(:, :, :)
    complex(kind=dp), allocatable :: MH(:, :), HM(:, :)
    complex(kind=dp), allocatable, intent(inout), optional ::  hmmn(:, :, :, :), mhmn(:, :, :, :)
    complex(kind=dp), allocatable, intent(inout), optional ::  dhmmn(:, :, :, :, :)
    complex(kind=dp), allocatable, intent(inout), optional ::  del_H(:, :, :, :)
    integer, allocatable :: num_states(:)
    real(kind=dp) :: c_real, c_imag
    integer :: tmp_bands, tmp_kpts, tmp_nntot
    integer :: ik, ik2, nnl, nnm, nnn, nn, inn, m, n, idir, idir2
    integer :: ncount
    logical :: nn_found

    allocate(num_states(num_kpts))
    do ik = 1, num_kpts
      if (have_disentangled) then
        num_states(ik) = dis_manifold%ndimwin(ik)
      else
        num_states(ik) = num_wann
      endif
    enddo

    write(stdout, *) "Reading mmn ..."
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
    allocate(mmn(num_wann, num_wann, 3, num_kpts), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating mmn in get_mmn', comm)
      return
    endif
    allocate(S_o(num_bands, num_bands), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating S_o in get_mmn', comm)
      return
    endif
    allocate(S(num_wann, num_wann), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating S in get_mmn', comm)
      return
    endif
    if (present(hmmn) .or. present(mhmn) .or. present(dhmmn)) then
      allocate(H_o(num_bands, num_bands, num_kpts), HM_o(num_bands, num_bands), MH_o(num_bands, num_bands), &
              HM(num_wann, num_wann), MH(num_wann, num_wann), &
              hmmn(num_wann, num_wann, 3, num_kpts), mhmn(num_wann, num_wann, 3, num_kpts), &
              dhmmn(num_wann, num_wann, 3, 3, num_kpts), &
              del_HM_o(num_bands, num_bands), del_HM(num_wann, num_wann, 3), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error in allocating H_o, hmmn, mhmn, or dhmmn in get_mmn', comm)
        return
      endif
      hmmn = cmplx_0
      mhmn = cmplx_0
      dhmmn = cmplx_0
      H_o = cmplx_0

      do ik = 1, num_kpts
        do m = 1, num_bands
          H_o(m, m, ik) = cmplx_1 * eigval(m, ik)
        enddo
      enddo
    endif

    if (present(dhmmn)) then
      allocate(del_H_o(num_bands, num_bands, 3, ik))
      del_H_o = cmplx_0
      do ik = 1, num_kpts
        do nn = 1, kmesh_info%nntot
          do idir = 1, 3
            del_H_o(:, :, idir, ik) = del_H_o(:, :, idir, ik) + &
                          kmesh_info%wb(nn)* kmesh_info%bk(idir, nn, ik) * &
                          H_o(:, :, kmesh_info%nnlist(ik, nn))
          enddo ! idir
        enddo ! nn
      enddo
    endif
    if (present(del_H)) then
      allocate(del_H(num_wann, num_wann, 3, num_kpts))
      del_H = cmplx_0
    endif
    mmn = cmplx_0


    do ncount = 1, num_kpts*kmesh_info%nntot
      !
      !Read from .mmn file the original overlap matrix
      ! S_o=<u_ik|u_ik2> between ab initio eigenstates
      !
      S_o = cmplx_0
      S = cmplx_0
      if (present(hmmn)) then
        HM_o = cmplx_0
        HM = cmplx_0
      endif
      if (present(mhmn)) then
        MH_o = cmplx_0
        MH = cmplx_0
      endif
      if (present(dhmmn)) then
        del_HM_o = cmplx_0
        del_HM = cmplx_0
      endif
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
      call get_gauge_overlap_matrix(num_bands, num_wann, eigval, v_matrix, dis_manifold, &
                                    ik, num_states(ik), kmesh_info%nnlist(ik, nn), &
                                    num_states(kmesh_info%nnlist(ik, nn)), S_o, &
                                    have_disentangled, S)
      if (present(hmmn)) then
        HM_o = matmul(H_o(:, :, ik), S_o) ! < um | H | del un > = < um | Em | del un >
        call get_gauge_overlap_matrix(num_bands, num_wann, eigval, v_matrix, dis_manifold, &
                                      ik, num_states(ik), kmesh_info%nnlist(ik, nn), &
                                      num_states(kmesh_info%nnlist(ik, nn)), HM_o, &
                                      have_disentangled, HM)
      endif
      if (present(mhmn)) then
        MH_o = matmul(S_o, H_o(:, :, ik)) ! < um | En | del un >
        call get_gauge_overlap_matrix(num_bands, num_wann, eigval, v_matrix, dis_manifold, &
                                      ik, num_states(ik), kmesh_info%nnlist(ik, nn), &
                                      num_states(kmesh_info%nnlist(ik, nn)), MH_o, &
                                      have_disentangled, MH)
      endif
      if (present(dhmmn)) then
        do idir = 1, 3
          del_HM_o(:, :) = matmul(del_H_o(:, :, idir, ik), S_o)
          call get_gauge_overlap_matrix(num_bands, num_wann, eigval, v_matrix, dis_manifold, &
                                        ik, num_states(ik), kmesh_info%nnlist(ik, nn), &
                                        num_states(kmesh_info%nnlist(ik, nn)), del_HM_o, &
                                        have_disentangled, del_HM(:, :, idir))
        enddo
      endif
      do idir = 1, 3
        mmn(:, :, idir, ik) = mmn(:, :, idir, ik) + &
          kmesh_info%wb(nn)*kmesh_info%bk(idir, nn, ik)*S(:, :)
        if (present(hmmn)) hmmn(:, :, idir, ik) = hmmn(:, :, idir, ik) + &
            kmesh_info%wb(nn)*kmesh_info%bk(idir, nn, ik)*HM(:, :)
        if (present(mhmn)) mhmn(:, :, idir, ik) = mhmn(:, :, idir, ik) + &
          kmesh_info%wb(nn)*kmesh_info%bk(idir, nn, ik)*MH(:, :)
        if (present(dhmmn)) then
          do idir2 = 1, 3
            dhmmn(:, :, idir2, idir, ik) = dhmmn(:, :, idir2, idir, ik) + &
              kmesh_info%wb(nn)*kmesh_info%bk(idir, nn, ik)*del_HM(:, :, idir2)
          enddo
        endif
      enddo
    enddo ! ncount over num_kpts * nntot
    deallocate(S_o, S)
    if (present(hmmn) .or. present(mhmn)) then
      deallocate(H_o, HM_o, HM, MH_o, MH)
    endif
    close(iun_mmn)

    do ik = 1, num_kpts
      do idir = 1, 3
        mmn(:, :, idir, ik) = 0.5_dp * (mmn(:, :, idir, ik) - &
                        conjg(transpose(mmn(:, :, idir, ik))))
      enddo
    enddo
    if (present(del_H)) then
      do ik = 1, num_kpts
        do idir = 1, 3
          call get_gauge_overlap_matrix(num_bands, num_wann, eigval, v_matrix, dis_manifold, &
                                        ik, num_states(ik), ik, num_states(ik), &
                                        del_H_o(:, :, idir, ik), have_disentangled, del_H(:, :, idir, ik))
        enddo
      enddo
    endif
    if (present(dhmmn)) then
      deallocate(del_H_o, del_HM_o, del_HM)
    endif
  end subroutine get_mmn

  subroutine get_uHu(stdout, seedname, num_bands, num_kpts, num_wann, kmesh_info, dis_manifold, have_disentangled, formatted, &
                     v_matrix, eigval, uhu, error, comm)
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
    type(dis_manifold_type), intent(in) :: dis_manifold
    
    character(len=60) :: header
    integer, intent(in) :: stdout
    character(len=50), intent(inout) :: seedname
    integer, intent(in) :: num_bands, num_kpts, num_wann
    logical, intent(in) :: have_disentangled
    complex(kind=dp), intent(in) :: v_matrix(:, :, :)
    real(kind=dp), intent(in) :: eigval(:, :)
    logical, intent(in) :: formatted
    complex(kind=dp), allocatable, intent(inout) :: uhu(:, :, :, :, :)
    ! <\nabla u_m| H \nabla|u_n> m, n, idir1, idir2, ik
    complex(kind=dp), allocatable :: Ho_qb1_q_qb2(:, :), H_qb1_q_qb2(:, :)
    integer, allocatable :: num_states(:)
    real(kind=dp) :: c_real, c_imag
    integer :: tmp_bands, tmp_kpts, tmp_nntot
    integer :: ik, nn1, nn2, m, n, idir1, idir2, qb1, qb2

    allocate(num_states(num_kpts))
    do ik = 1, num_kpts
      if (have_disentangled) then
        num_states(ik) = dis_manifold%ndimwin(ik)
      else
        num_states(ik) = num_wann
      endif
    enddo

    write(stdout, *) "Reading uHu ..."
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
    allocate(uhu(num_wann, num_wann, 3, 3, num_kpts), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating uHu in get_uHu', comm)
      return
    endif
    uhu = cmplx_0
    allocate(Ho_qb1_q_qb2(num_bands, num_bands), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating Ho_qb1_q_qb2 in get_uHu', comm)
      return
    endif
    allocate(H_qb1_q_qb2(num_wann, num_wann), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating H_qb1_q_qb2 in get_uHu', comm)
      return
    endif
    
    do ik = 1, num_kpts
      Ho_qb1_q_qb2 = cmplx_0
      H_qb1_q_qb2 = cmplx_0
      do nn2 = 1, kmesh_info%nntot
        qb2 = kmesh_info%nnlist(ik, nn2)
        do nn1 = 1, kmesh_info%nntot
          qb1 = kmesh_info%nnlist(ik, nn1)
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
          call get_gauge_overlap_matrix(num_bands, num_wann, eigval, v_matrix, dis_manifold, &
                                        qb1, num_states(qb1), qb2, num_states(qb2), &
                                        Ho_qb1_q_qb2, have_disentangled, H_qb1_q_qb2)
          do idir2 = 1, 3
            do idir1 = 1, idir2
              uhu(:, :, idir1, idir2, ik) = uhu(:, :, idir1, idir2, ik) + &
                                            kmesh_info%wb(nn1)*kmesh_info%bk(idir1, nn1, ik)* &
                                            kmesh_info%wb(nn2)*kmesh_info%bk(idir2, nn2, ik)*H_qb1_q_qb2(:, :)
            enddo ! idir1
          enddo ! idir2
        enddo ! nn1
      enddo ! nn2
    enddo ! ik
    close(iun_uHu)
    ! End file read
    do ik = 1, num_kpts
      do idir2 = 1, 3
        do idir1 = 1, idir2
          uhu(:, :, idir2, idir1, ik) = conjg(transpose(uhu(:, :, idir1, idir2, ik)))
        enddo
      enddo
    enddo
    deallocate(Ho_qb1_q_qb2, H_qb1_q_qb2)
  end subroutine get_uHu

  subroutine get_uIu(stdout, seedname, num_bands, num_kpts, num_wann, kmesh_info, dis_manifold, have_disentangled, formatted, &
                     v_matrix, eigval, uiu, error, comm)
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
    type(dis_manifold_type), intent(in) :: dis_manifold
    
    character(len=60) :: header
    integer, intent(in) :: stdout
    character(len=50), intent(inout) :: seedname
    integer, intent(in) :: num_bands, num_kpts, num_wann
    logical, intent(in) :: have_disentangled
    complex(kind=dp), intent(in) :: v_matrix(:, :, :)
    real(kind=dp), intent(in) :: eigval(:, :)
    logical, intent(in) :: formatted
    complex(kind=dp), allocatable, intent(inout) :: uiu(:, :, :, :, :)
    ! <\nabla u_m| H \nabla|u_n> m, n, idir1, idir2, ik
    complex(kind=dp), allocatable :: Lo_qb1_q_qb2(:, :), L_qb1_q_qb2(:,:)
    integer, allocatable :: num_states(:)
    real(kind=dp) :: c_real, c_imag
    integer :: tmp_bands, tmp_kpts, tmp_nntot
    integer :: ik, nn1, nn2, m, n, idir1, idir2, qb1, qb2

    allocate(num_states(num_kpts))
    do ik = 1, num_kpts
      if (have_disentangled) then
        num_states(ik) = dis_manifold%ndimwin(ik)
      else
        num_states(ik) = num_wann
      endif
    enddo

    write(stdout, *) "Reading uIu ..."
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
    allocate(uiu(num_wann, num_wann, 3, 3, num_kpts), stat=ierr)
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
    allocate(L_qb1_q_qb2(num_wann, num_wann), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating uIu and L_qb1_q_qb2 in get_uIu', comm)
      return
    endif
    
    do ik = 1, num_kpts
      Lo_qb1_q_qb2 = cmplx_0
      L_qb1_q_qb2 = cmplx_0
      do nn2 = 1, kmesh_info%nntot
        qb2 = kmesh_info%nnlist(ik, nn2)
        do nn1 = 1, kmesh_info%nntot
          qb1 = kmesh_info%nnlist(ik, nn1)
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
          do m = 1, num_bands
            do n = 1, num_bands
              Lo_qb1_q_qb2(n, m) = eigval(m, ik) * Lo_qb1_q_qb2(n, m)
            enddo
          enddo
          call get_gauge_overlap_matrix(num_bands, num_wann, eigval, v_matrix, dis_manifold, &
                                        qb1, num_states(qb1), qb2, num_states(qb2), &
                                        Lo_qb1_q_qb2, have_disentangled, L_qb1_q_qb2)
          do idir2 = 1, 3
            do idir1 = 1, 3
              uiu(:, :, idir1, idir2, ik) = uiu(:, :, idir1, idir2, ik) + &
                                            kmesh_info%wb(nn1)*kmesh_info%bk(idir1, nn1, ik)* &
                                            kmesh_info%wb(nn2)*kmesh_info%bk(idir2, nn2, ik)*L_qb1_q_qb2(:, :)
            enddo ! idir1
          enddo ! idir2
        enddo ! nn1
      enddo ! nn2
    enddo ! ik
    close(iun_uIu)
    ! End file read
    ! do ik = 1, num_kpts
    !   do idir2 = 1, 3
    !     do idir1 = 1, idir2
    !       uiu(:, :, idir2, idir1, ik) = conjg(transpose(uiu(:, :, idir1, idir2, ik)))
    !     enddo
    !   enddo
    ! enddo
    deallocate(Lo_qb1_q_qb2, L_qb1_q_qb2)
  end subroutine get_uIu

  subroutine calc_del (stdout, num_bands, num_wann, num_kpts, v_matrix, dv, eigval, del_eig, kmesh_info, error, comm)
    !================================================!
    !
    !! calculate del v_matrix and del eigval
    !! Note that del_v * v^dagger * v = del_v
    !! We define dv = del_v * v^dagger
    !! So del (|u>V) = |del u> V + |u> dV V
    !
    !================================================!
    ! only run on root node
    implicit none

    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error
    type(kmesh_info_type), intent(in) :: kmesh_info
    
    integer, intent(in) :: stdout
    integer, intent(in) :: num_bands, num_wann, num_kpts
    complex(kind=dp), allocatable, intent(in) :: v_matrix(:, :, :)
    complex(kind=dp), allocatable, intent(inout) :: dv(:, :, :, :)
    real(kind=dp), intent(in) :: eigval(:, :)
    real(kind=dp), allocatable, intent(inout) :: del_eig(:, :, :)
    complex(kind=dp), allocatable :: del_v(:, :, :)
    integer :: ik, nn, qb, idir
    if (allocated(dv)) then
      call set_error_input(error, 'Error: dv allocated before allcated', comm)
      return
    endif
    if (allocated(del_eig)) then
      call set_error_input(error, 'Error: dv allocated before allcated', comm)
      return
    endif
    allocate(dv(num_bands, num_bands, 3, num_kpts))
    allocate(del_v(num_bands, num_wann, 3))
    allocate(del_eig(num_bands, 3, num_kpts))
    dv = cmplx_0
    del_eig = 0.0_dp
    do ik = 1, num_kpts
      del_v = cmplx_0
      do nn = 1, kmesh_info%nntot
        do idir = 1, 3
          qb = kmesh_info%nnlist(ik, nn)
          del_v(:, :, idir) = del_v(:, :, idir) + &
                        kmesh_info%wb(nn)* kmesh_info%bk(idir, nn, ik) * v_matrix(:, :, qb)
          del_eig(:, idir, ik) = del_eig(:, idir, ik) + &
                        kmesh_info%wb(nn)* kmesh_info%bk(idir, nn, ik) * eigval(:, qb)
        enddo ! idir
      enddo ! nn
      do idir = 1, 3
        dv(:, :, idir, ik) = matmul(del_v(:, :, idir), conjg(transpose(v_matrix(:, :, ik))))
      enddo
    enddo ! ik
  end subroutine

  subroutine calc_orb(stdout, num_bands, num_kpts, num_wann, eigval, del_H, v_matrix, mmn, hmmn, mhmn, dhmmn, uhu, uiu, orb, error, comm)
    !================================================!
    !
    !! calculate orbital matrix
    !
    !================================================!
    ! only run on root node
    implicit none

    type(w90_comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error
    
    integer, intent(in) :: stdout
    integer, intent(in) :: num_bands, num_kpts, num_wann
    real(kind=dp), intent(in) :: eigval(:, :)
    complex(kind=dp), allocatable, intent(inout) :: mmn(:, :, :, :)
    ! <u_m|\nabla|u_n> m, n, idir, ik
    complex(kind=dp), allocatable, intent(inout) :: hmmn(:, :, :, :), mhmn(:, :, :, :)
    ! <u_m|\nabla|u_n> m, n, idir, ik
    complex(kind=dp), allocatable, intent(inout) :: dhmmn(:, :, :, :, :)
    ! <u_m|\nabla|u_n> m, n, idir1, idir2, ik
    complex(kind=dp), allocatable, intent(inout) :: uhu(:, :, :, :, :)
    ! <\nabla u_m| H \nabla|u_n> m, n, idir1, idir2, ik
    complex(kind=dp), allocatable, intent(inout) :: uiu(:, :, :, :, :)
    ! <\nabla u_m| \nabla|u_n> m, n, idir1, idir2, ik
    complex(kind=dp), allocatable, intent(in) :: v_matrix(:, :, :)
    ! m, n, idir, ik
    complex(kind=dp), allocatable, intent(inout) :: orb(:, :, :, :)
    complex(kind=dp), allocatable :: orb_ab(:, :, :, :)
    complex(kind=dp), allocatable, intent(in) :: del_H(:, :, :, :)
    ! m, n, idir, ik
    complex(kind=dp), allocatable ::  left(:, :), right(:, :)
    integer :: ik, m, n, t, idir, a, b

    integer, dimension(3), parameter :: alpha_A = (/2, 3, 1/)
    integer, dimension(3), parameter :: beta_A = (/3, 1, 2/)
    real(kind=dp), parameter :: fac = 3.674932379e-2_dp / (0.52917721092_dp)**2

    if (allocated(orb)) then
      call set_error_input(error, 'Error: orb_o matrix has been allocated before calculated', comm)
      return
    endif
    allocate(orb(num_wann, num_wann, 3, num_kpts), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating orb_o in calc_orb', comm)
      return
    endif
    allocate(orb_ab(num_wann, num_wann, 3, 3), stat=ierr)
    orb = cmplx_0
    allocate(left(num_wann, num_wann), right(num_wann,num_wann))
    left = cmplx_0
    right = cmplx_0
    do ik = 1, num_kpts
      do b = 1, 3
        do a = 1, 3
          left = conjg(transpose(mmn(:, :, a, ik)))
          right = hmmn(:, :, b, ik) - del_H(:, :, b, ik) - mhmn(:, :, b, ik)
          orb_ab(:, :, a, b) = cmplx_i * ( &
              uhu(:, :, a, b, ik) - &
              conjg(transpose(dhmmn(:, :, b, a, ik))) - &
              uiu(:, :, a, b, ik) - &
              matmul(left, right) &
            )
        enddo ! a
      enddo ! b
      do idir = 1, 3
        orb(:, :, idir, ik) = orb_ab(:, :, alpha_A(idir), beta_A(idir)) - orb_ab(:, :, beta_A(idir), alpha_A(idir))
        orb(:, :, idir, ik) = fac * (orb(:, :, idir, ik) + conjg(transpose(orb(:, :, idir, ik))))
      enddo
    enddo ! k
    deallocate(orb_ab)
    deallocate(mmn, uhu, uiu)
  end subroutine

  subroutine output_orb_formatted(stdout, seedname, num_bands, num_kpts, num_wann, orb, v_matrix)
    !================================================!
    !
    !! write orbital matrix to seedname.orb.fmt
    !
    !================================================!
    ! only run on root node
    use w90_io, only: io_date

    implicit none
    
    integer, intent(in) :: stdout
    character(len=50), intent(inout) :: seedname
    integer, intent(in) :: num_bands, num_kpts, num_wann
    complex(kind=dp), allocatable, intent(inout) :: orb(:, :, :, :)
    ! m, n, idir, ik
    complex(kind=dp), allocatable, intent(inout) :: v_matrix(:, :, :)
    ! num_bands, num_wann, ik
    complex(kind=dp), allocatable :: orb_o_k(:, :, :)
    integer :: ik, m, n, idir
    character(len=60) header
    character(len=9) cdate, ctime
    call io_date(cdate, ctime)
    header = 'Created on ' // cdate // ' at ' // ctime

    write (stdout, '(3a)') 'Writing information to formatted file ', trim(seedname), '.orb.fmt :'

    open (newunit=iun_orb, file=trim(seedname)//'.orb.fmt', form='formatted', status='replace', position='rewind')
    allocate(orb_o_k(num_bands, num_bands, 3))

    write(iun_orb, *) header
    write(iun_orb, *) num_bands, num_kpts
    do ik = 1, num_kpts
      orb_o_k = cmplx_0
      do idir = 1, 3
        orb_o_k(:, :, idir) = matmul(v_matrix(:, :, ik), matmul(orb(:, :, idir, ik), conjg(transpose(v_matrix(:, :, ik)))))
      enddo
      do m = 1, num_bands
        do n = 1, m
          do idir = 1, 3
            write(iun_orb, '(2es26.16)') orb_o_k(n, m, idir)
          enddo
        enddo
      enddo
    enddo
    close(iun_orb)

  end subroutine
end module get_orb


program w90genorb
  !! Program to convert spn files from formatted to unformmated
  use w90_constants, only: dp, pw90_physical_constants_type
  use get_orb, only: get_seedname, get_mmn, get_uHu, get_uIu, &
    calc_del, calc_orb, output_orb_formatted
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
  real(kind=dp), allocatable :: del_eig(:, :, :)
  ! iband, ik, idir
  integer :: mp_grid(3)
  integer :: optimisation
  real(kind=dp), allocatable :: fermi_energy_list(:)
  logical :: gamma_only, eig_found, formatted
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
  complex(kind=dp), allocatable :: v_matrix(:, :, :), dv(:, :, :, :)
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


  complex(kind=dp), allocatable :: mmn(:, :, :, :)
  ! <u_m|\nabla|u_n> m, n, idir, ik
  complex(kind=dp), allocatable :: hmmn(:, :, :, :)
  ! <u_m|H|del u_n> = <u_m|E_m|del u_n> m, n, idir, ik
  complex(kind=dp), allocatable :: dhmmn(:, :, :, :, :)
  ! <u_m|H|del u_n> = <u_m|del E_m|del u_n> m, n, idir1, idir2, ik
  complex(kind=dp), allocatable :: mhmn(:, :, :, :)
  ! <u_m|E_n|del u_n> m, n, idir, ik
  complex(kind=dp), allocatable :: del_H(:, :, :, :)
  ! < u | del E | u > in Wannier Gauge
  complex(kind=dp), allocatable :: uhu(:, :, :, :, :)
  ! <\nabla u_m| H \nabla|u_n> m, n, idir1, idir2, ik
  complex(kind=dp), allocatable :: uiu(:, :, :, :, :)
  ! <\nabla u_m| \nabla|u_n> m, n, idir1, idir2, ik
  complex(kind=dp), allocatable :: orb(:, :, :, :)
  ! m, n, idir, ik
  integer :: i, j, k, idir
  complex(kind=dp), allocatable :: temp(:, :, :)


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
    call get_seedname(stdout, seedname, formatted)
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

    call get_mmn(stdout, seedname, num_bands, num_kpts, num_wann, kmesh_info, dis_window, have_disentangled, &
                 v_matrix, eigval, mmn ,error, comm, hmmn=hmmn, mhmn=mhmn, dhmmn=dhmmn, del_H=del_H)
    if (allocated(error)) call print_error_halt(error, ierr, stdout, stderr, comm)
    if (allocated(mmn)) then
      write(stdout, *) "Reading mmn ... Done"
    endif

    call get_uHu(stdout, seedname, num_bands, num_kpts, num_wann, kmesh_info, dis_window, have_disentangled, formatted, &
                 v_matrix, eigval, uhu ,error, comm)
    if (allocated(error)) call print_error_halt(error, ierr, stdout, stderr, comm)
    if (allocated(uhu)) then
      write(stdout, *) "Reading uHu ... Done"
    endif
    call get_uIu(stdout, seedname, num_bands, num_kpts, num_wann, kmesh_info, dis_window, have_disentangled, formatted, &
                 v_matrix, eigval, uiu ,error, comm)
    if (allocated(error)) call print_error_halt(error, ierr, stdout, stderr, comm)
    if (allocated(uiu)) then
      write(stdout, *) "Reading uIu ... Done"
    endif
    ! allocate(temp(num_bands, num_bands, 3))

    call calc_orb(stdout, num_bands, num_kpts, num_wann, &
                  eigval, del_H, v_matrix, mmn, hmmn, mhmn, dhmmn, uhu, uiu, orb, error, comm)
    call output_orb_formatted(stdout, seedname, num_bands, num_kpts, num_wann, orb, v_matrix)
    if (allocated(orb)) deallocate(orb)
    write(stdout, *) "TEST exiting..."
    close (unit=stdout)
  endif

end program w90genorb


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! back up deprecated codes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! subroutine calc_orb(stdout, num_bands, num_kpts, eigval, del_eig, v_matrix, dv, mmn, uhu, uiu, orb_o, error, comm)
  !   !================================================!
  !   !
  !   !! calculate orbital matrix
  !   !
  !   !================================================!
  !   ! only run on root node
  !   implicit none

  !   type(w90_comm_type), intent(in) :: comm
  !   type(w90_error_type), allocatable, intent(out) :: error
    
  !   integer, intent(in) :: stdout
  !   integer, intent(in) :: num_bands, num_kpts
  !   real(kind=dp), intent(in) :: eigval(:, :)
  !   complex(kind=dp), allocatable, intent(inout) :: mmn(:, :, :, :)
  !   ! <u_m|\nabla|u_n> m, n, idir, ik
  !   complex(kind=dp), allocatable, intent(inout) :: uhu(:, :, :, :, :)
  !   ! <\nabla u_m| H \nabla|u_n> m, n, idir1, idir2, ik
  !   complex(kind=dp), allocatable, intent(inout) :: uiu(:, :, :, :, :)
  !   ! <\nabla u_m| \nabla|u_n> m, n, idir1, idir2, ik
  !   complex(kind=dp), allocatable, intent(in) :: v_matrix(:, :, :)
  !   complex(kind=dp), allocatable, intent(in) :: dv(:, :, :, :)
  !   ! m, n, idir, ik
  !   real(kind=dp), allocatable, intent(in) :: del_eig(:, :, :)
  !   ! bands, idir, num_kpts
  !   complex(kind=dp), allocatable, intent(inout) :: orb_o(:, :, :, :)
  !   complex(kind=dp), allocatable :: orb_o_ab(:, :, :, :)
  !   complex(kind=dp), allocatable :: H_o(:, :), del_H_o(:, :, :)
  !   complex(kind=dp), allocatable :: dvh(:, :, :)
  !   complex(kind=dp), allocatable :: left(:, :), right(:, :)
  !   ! m, n, idir, ik
  !   integer :: ik, m, n, t, idir, a, b

  !   integer, dimension(3), parameter :: alpha_A = (/2, 3, 1/)
  !   integer, dimension(3), parameter :: beta_A = (/3, 1, 2/)

  !   if (allocated(orb_o)) then
  !     call set_error_input(error, 'Error: orb_o matrix has been allocated before calculated', comm)
  !     return
  !   endif
  !   allocate(orb_o(num_bands, num_bands, 3, num_kpts), stat=ierr)
  !   if (ierr /= 0) then
  !     call set_error_alloc(error, 'Error in allocating orb_o in calc_orb', comm)
  !     return
  !   endif
  !   allocate(orb_o_ab(num_bands, num_bands, 3, 3), stat=ierr)
  !   allocate(H_o(num_bands, num_bands), del_H_o(num_bands, num_bands, 3))
  !   allocate(dvh(num_bands, num_bands, 3))
  !   allocate(left(num_bands, num_bands), right(num_bands, num_bands))
  !   orb_o = cmplx_0
  !   do ik = 1, num_kpts
  !     orb_o_ab = cmplx_0
  !     H_o = cmplx_0
  !     del_H_o = cmplx_0
  !     do m = 1, num_bands
  !       H_o(m, m) = eigval(m, ik) * cmplx_1
  !       do idir = 1, 3
  !         del_H_o(m, m, idir) = del_eig(m, idir, ik) * cmplx_1
  !         dvh(:, :, idir) = conjg(transpose(dv(:, :, idir, ik)))
  !       enddo
  !     enddo
  !     do b = 1, 3
  !       do a = 1, 3
  !         left = dvh(:, :, a) + conjg(transpose(mmn(:, :, a, ik)))
  !         left = matmul(left, v_matrix(:, :, ik))

  !         right = matmul(mmn(:, :, b, ik), H_o) &
  !               - matmul(H_o, mmn(:, :, b, ik)) &
  !               + del_H_o(:, :, b)
  !         right = matmul(conjg(transpose(v_matrix(:, :, ik))), right)
  !         orb_o_ab(:, :, a, b) =  orb_o_ab(:, :, a, b) -cmplx_i * (&
  !           matmul(dvh(:, :, a), matmul(mmn(:, :, b, ik), H_o)) &
  !           + matmul(uiu(:, :, a, b, ik), H_o) &
  !           - matmul(dvh(:, :, a), matmul(H_o, mmn(:, :, b, ik))) &
  !           - uhu(:, :, a, b, ik) &
  !           + matmul(dvh(:, :, a), del_H_o(:, :, b)) &
  !           + matmul(conjg(transpose(mmn(:, :, a, ik))), del_H_o(:, :, b)) &
  !         ) + cmplx_i * matmul(left, right)
  !       enddo
  !     enddo
  !     do idir = 3, 3
  !       orb_o(:, :, idir, ik) = orb_o_ab(:, :, alpha_A(idir), beta_A(idir)) &
  !                             - orb_o_ab(:, :, beta_A(idir), alpha_A(idir))
  !     enddo
  !     ! do m = 1, num_bands
  !     !   do n = 1, num_bands
  !     !     ! idir = 1 -- x
  !     !     orb_o(n, m, 1, ik) = orb_o(n, m, 1, ik) + cmplx_i * &
  !     !       (uhu(n, m, 2, 3, ik) - uhu(n, m, 3, 2, ik) - & 
  !     !       eigval(m, ik)*(uiu(n, m, 2, 3, ik) - uiu(n, m, 3, 2, ik)))
  !     !     do t = 1, num_bands
  !     !       orb_o(n, m, 1, ik) = orb_o(n, m, 1, ik) - cmplx_i * (eigval(t, ik) - eigval(m, ik)) * &
  !     !       (conjg(mmn(t, n, 2, ik)) * mmn(t, m, 3, ik) - &
  !     !       conjg(mmn(t, n, 3, ik)) * mmn(t, m, 2, ik))
  !     !     enddo
  !     !     ! idir = 2 -- y
  !     !     orb_o(n, m, 2, ik) = orb_o(n, m, 2, ik) + cmplx_i * &
  !     !       (uhu(n, m, 3, 1, ik) - uhu(n, m, 1, 3, ik) - & 
  !     !       eigval(m, ik)  * (uiu(n, m, 3, 1, ik) - uiu(n, m, 1, 3, ik)))
  !     !     do t = 1, num_bands
  !     !       orb_o(n, m, 2, ik) = orb_o(n, m, 2, ik) - cmplx_i * &
  !     !       (conjg(mmn(t, n, 3, ik)) * (eigval(t, ik) - eigval(m, ik)) * mmn(t, m, 1, ik) - &
  !     !       conjg(mmn(t, n, 1, ik)) * (eigval(t, ik) - eigval(m, ik)) * mmn(t, m, 3, ik))
  !     !     enddo
  !     !     ! idir = 3 -- z
  !     !     orb_o(n, m, 3, ik) = orb_o(n, m, 3, ik) + cmplx_i * &
  !     !       (uhu(n, m, 1, 2, ik) - uhu(n, m, 2, 1, ik) - & 
  !     !       eigval(m, ik)  * (uiu(n, m, 1, 2, ik) - uiu(n, m, 2, 1, ik)))
  !     !     do t = 1, num_bands
  !     !       orb_o(n, m, 3, ik) = orb_o(n, m, 3, ik) - cmplx_i * &
  !     !       (conjg(mmn(t, n, 1, ik)) * (eigval(t, ik) - eigval(m, ik)) * mmn(t, m, 2, ik) - &
  !     !       conjg(mmn(t, n, 2, ik)) * (eigval(t, ik) - eigval(m, ik)) * mmn(t, m, 1, ik))
  !     !     enddo
  !     !   enddo ! n
  !     ! enddo ! m
  !     do idir = 1, 3
  !       orb_o(:, :, idir, ik) = 0.5_dp * (orb_o(:, :, idir, ik) + &
  !                         conjg(transpose(orb_o(:, :, idir, ik))))
  !     enddo
  !   enddo ! ik
  !   deallocate(orb_o_ab, H_o, del_H_o, left, right)
  !   deallocate(mmn, uhu, uiu)
  ! end subroutine
