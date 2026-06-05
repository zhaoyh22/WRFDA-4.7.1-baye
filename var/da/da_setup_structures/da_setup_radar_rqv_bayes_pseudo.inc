subroutine da_setup_radar_rqv_bayes_pseudo(grid, iv, n1, n2)

   !---------------------------------------------------------------------------
   ! Purpose: Append internal rqv-only pseudo radar observations for grid points
   !          with no radar record but strong background reflectivity.
   !---------------------------------------------------------------------------

   implicit none

   type(domain),  intent(in)    :: grid
   type(iv_type), intent(inout) :: iv
   integer,       intent(in)    :: n1, n2

   integer :: old_nlocal, old_ntotal, old_ptotal
   integer :: new_nlocal
   integer :: pseudo_count, pseudo_total, pseudo_offset, pseudo_filled
   integer :: h_stride, v_stride, pseudo_max
   integer :: i_beg, i_end, j_beg, j_end, k_beg, k_end
   integer :: i_min, i_max, j_min, j_max
   integer :: i_min_l, i_max_l, j_min_l, j_max_l
   integer :: i, j, k, n, lev, ii, jj, kk
   integer :: radius
   real    :: qgr, bg_rf, height
   logical :: found_obs
   logical, allocatable :: obs_mask(:,:,:)
#ifdef DM_PARALLEL
   integer, allocatable :: pseudo_counts(:)
   integer :: found_obs_l, found_obs_g
#endif

   if (trace_use) call da_trace_entry("da_setup_radar_rqv_bayes_pseudo")

   old_nlocal = iv%info(radar)%nlocal
   old_ntotal = iv%info(radar)%ntotal
   old_ptotal = iv%info(radar)%ptotal(iv%time) - iv%info(radar)%ptotal(iv%time-1)

   if (old_ptotal < 1) then
      if (trace_use) call da_trace_exit("da_setup_radar_rqv_bayes_pseudo")
      return
   end if

   h_stride = max(1, radar_rqv_bayes_pseudo_h_stride)
   v_stride = max(1, radar_rqv_bayes_pseudo_v_stride)
   pseudo_max = max(0, radar_rqv_bayes_pseudo_max)

   i_beg = max(its, ids)
   i_end = min(ite, ide-1, ime-1)
   j_beg = max(jts, jds)
   j_end = min(jte, jde-1, jme-1)
   k_beg = max(kts, kds)
   k_end = min(kte, kde)

   if (i_beg > i_end .or. j_beg > j_end .or. k_beg > k_end) then
      if (trace_use) call da_trace_exit("da_setup_radar_rqv_bayes_pseudo")
      return
   end if

   allocate(obs_mask(i_beg:i_end,j_beg:j_end,k_beg:k_end))
   obs_mask(:,:,:) = .false.

   i_min_l = ide
   i_max_l = ids
   j_min_l = jde
   j_max_l = jds
   found_obs = .false.

   if (n2 >= n1) then
      do n = n1, n2
         if (iv%info(radar)%levels(n) < 1) cycle

         ii = nint(iv%info(radar)%x(1,n))
         jj = nint(iv%info(radar)%y(1,n))
         if (ii < i_beg .or. ii > i_end .or. jj < j_beg .or. jj > j_end) cycle

         found_obs = .true.
         i_min_l = min(i_min_l, ii)
         i_max_l = max(i_max_l, ii)
         j_min_l = min(j_min_l, jj)
         j_max_l = max(j_max_l, jj)

         do lev = 1, iv%info(radar)%levels(n)
            kk = k_beg
            height = iv%radar(n)%height(lev)
            do k = k_beg, k_end
               if (abs(grid%xb%h(ii,jj,k)-height) < abs(grid%xb%h(ii,jj,kk)-height)) kk = k
            end do
            obs_mask(ii,jj,kk) = .true.
         end do
      end do
   end if

#ifdef DM_PARALLEL
   found_obs_l = 0
   if (found_obs) found_obs_l = 1
   call mpi_allreduce(found_obs_l, found_obs_g, 1, mpi_integer, mpi_sum, comm, ierr)
   call mpi_allreduce(i_min_l, i_min, 1, mpi_integer, mpi_min, comm, ierr)
   call mpi_allreduce(i_max_l, i_max, 1, mpi_integer, mpi_max, comm, ierr)
   call mpi_allreduce(j_min_l, j_min, 1, mpi_integer, mpi_min, comm, ierr)
   call mpi_allreduce(j_max_l, j_max, 1, mpi_integer, mpi_max, comm, ierr)

   if (found_obs_g < 1) then
      deallocate(obs_mask)
      if (trace_use) call da_trace_exit("da_setup_radar_rqv_bayes_pseudo")
      return
   end if
#else
   if (.not. found_obs) then
      deallocate(obs_mask)
      if (trace_use) call da_trace_exit("da_setup_radar_rqv_bayes_pseudo")
      return
   end if
   i_min = i_min_l
   i_max = i_max_l
   j_min = j_min_l
   j_max = j_max_l
#endif

   radius = max(0, radar_rqv_bayes_clear_radius)
   i_min = max(i_beg, i_min-radius)
   i_max = min(i_end, i_max+radius)
   j_min = max(j_beg, j_min-radius)
   j_max = min(j_end, j_max+radius)

   pseudo_count = 0

   count_j: do j = j_min, j_max
      if (mod(j-jds,h_stride) /= 0) cycle
      do i = i_min, i_max
         if (mod(i-ids,h_stride) /= 0) cycle
         do k = k_beg, k_end
            if (mod(k-kds,v_stride) /= 0) cycle
            if (obs_mask(i,j,k)) cycle
            if (radar_rqv_h_lbound >= 0.0 .and. radar_rqv_h_ubound >= 0.0) then
               if (grid%xb%h(i,j,k) < radar_rqv_h_lbound .or. &
                   grid%xb%h(i,j,k) > radar_rqv_h_ubound) cycle
            end if

            qgr = 0.0
            if (associated(grid%xb%qgr)) qgr = grid%xb%qgr(i,j,k)
            call da_setup_radar_rqv_bayes_rf(grid%xb%qrn(i,j,k), grid%xb%qsn(i,j,k), &
               qgr, grid%xb%t(i,j,k)-273.15, grid%xb%rho(i,j,k), bg_rf)
            if (bg_rf < radar_rqv_bayes_bg_thresh) cycle

            pseudo_count = pseudo_count + 1
            if (pseudo_max > 0 .and. pseudo_count >= pseudo_max) exit count_j
         end do
      end do
   end do count_j

#ifdef DM_PARALLEL
   allocate(pseudo_counts(0:num_procs-1))
   call mpi_allgather(pseudo_count, 1, mpi_integer, pseudo_counts, 1, mpi_integer, &
      comm, ierr)
   pseudo_total = sum(pseudo_counts)
   pseudo_offset = 0
   if (myproc > 0) pseudo_offset = sum(pseudo_counts(0:myproc-1))
   deallocate(pseudo_counts)
#else
   pseudo_total = pseudo_count
   pseudo_offset = 0
#endif

   if (pseudo_total < 1) then
      deallocate(obs_mask)
      if (trace_use) call da_trace_exit("da_setup_radar_rqv_bayes_pseudo")
      return
   end if

   if (pseudo_count > 0) then
      new_nlocal = old_nlocal + pseudo_count
      call da_extend_radar_rqv_bayes_pseudo(iv, new_nlocal)

      pseudo_filled = 0

      fill_j: do j = j_min, j_max
         if (mod(j-jds,h_stride) /= 0) cycle
         do i = i_min, i_max
            if (mod(i-ids,h_stride) /= 0) cycle
            do k = k_beg, k_end
               if (mod(k-kds,v_stride) /= 0) cycle
               if (obs_mask(i,j,k)) cycle
               if (radar_rqv_h_lbound >= 0.0 .and. radar_rqv_h_ubound >= 0.0) then
                  if (grid%xb%h(i,j,k) < radar_rqv_h_lbound .or. &
                      grid%xb%h(i,j,k) > radar_rqv_h_ubound) cycle
               end if

               qgr = 0.0
               if (associated(grid%xb%qgr)) qgr = grid%xb%qgr(i,j,k)
               call da_setup_radar_rqv_bayes_rf(grid%xb%qrn(i,j,k), grid%xb%qsn(i,j,k), &
                  qgr, grid%xb%t(i,j,k)-273.15, grid%xb%rho(i,j,k), bg_rf)
               if (bg_rf < radar_rqv_bayes_bg_thresh) cycle

               pseudo_filled = pseudo_filled + 1
               n = old_nlocal + pseudo_filled

               iv%info(radar)%levels(n) = 1
               iv%info(radar)%name(n) = 'RQV_BAYES_PSEUDO'
               iv%info(radar)%platform(n) = 'RQV_BAYES'
               iv%info(radar)%id(n) = 'RQV_BAYES_PSEUDO'
               iv%info(radar)%date_char(n) = analysis_date
               iv%info(radar)%lat(:,n) = grid%xb%lat(i,j)
               iv%info(radar)%lon(:,n) = grid%xb%lon(i,j)
               iv%info(radar)%elv(n) = grid%xb%terr(i,j)
               iv%info(radar)%x(:,n) = real(i)
               iv%info(radar)%y(:,n) = real(j)
               iv%info(radar)%i(:,n) = i
               iv%info(radar)%j(:,n) = j
               iv%info(radar)%dx(:,n) = 0.0
               iv%info(radar)%dxm(:,n) = 1.0
               iv%info(radar)%dy(:,n) = 0.0
               iv%info(radar)%dym(:,n) = 1.0
               iv%info(radar)%k(:,n) = 0
               iv%info(radar)%dz(:,n) = 0.0
               iv%info(radar)%dzm(:,n) = 0.0
               iv%info(radar)%zk(:,n) = missing_r
               iv%info(radar)%proc_domain(:,n) = .true.
               iv%info(radar)%thinned(:,n) = .false.
               iv%info(radar)%obs_global_index(n) = old_ntotal + pseudo_offset + pseudo_filled

               call da_init_radar_rqv_bayes_pseudo(iv%radar(n), iv%info(radar)%max_lev)
               iv%radar(n)%stn_loc%lon = grid%xb%lon(i,j)
               iv%radar(n)%stn_loc%lat = grid%xb%lat(i,j)
               iv%radar(n)%stn_loc%elv = grid%xb%terr(i,j)
               iv%radar(n)%stn_loc%x = real(i)
               iv%radar(n)%stn_loc%y = real(j)
               iv%radar(n)%stn_loc%zk = real(k)
               iv%radar(n)%height(1) = grid%xb%h(i,j,k)
               iv%radar(n)%height_qc(1) = missing

               if (pseudo_filled >= pseudo_count) exit fill_j
            end do
         end do
      end do fill_j
   end if

   iv%info(radar)%nlocal = old_nlocal + pseudo_count
   iv%info(radar)%ntotal = old_ntotal + pseudo_total
   iv%info(radar)%plocal(iv%time:max_fgat_time) = &
      iv%info(radar)%plocal(iv%time:max_fgat_time) + pseudo_count
   iv%info(radar)%ptotal(iv%time:max_fgat_time) = &
      iv%info(radar)%ptotal(iv%time:max_fgat_time) + pseudo_total

   if (rootproc .and. pseudo_total > 0) then
      write(unit=stdout,fmt='(A,I10,A,I3)') &
         ' Added internal radar rqv Bayesian pseudo obs: ', pseudo_total, &
         ' at FGAT time ', iv%time
   end if

   deallocate(obs_mask)

   if (trace_use) call da_trace_exit("da_setup_radar_rqv_bayes_pseudo")

end subroutine da_setup_radar_rqv_bayes_pseudo


subroutine da_extend_radar_rqv_bayes_pseudo(iv, new_nlocal)

   !---------------------------------------------------------------------------
   ! Purpose: Extend radar iv/info arrays before appending pseudo observations.
   !---------------------------------------------------------------------------

   implicit none

   type(iv_type), intent(inout) :: iv
   integer,       intent(in)    :: new_nlocal

   integer :: old_nlocal, max_lev
   type(radar_type), allocatable :: radar_save(:)
   character(len=40), allocatable :: name_save(:)
   character(len=12), allocatable :: platform_save(:)
   character(len=40), allocatable :: id_save(:)
   character(len=19), allocatable :: date_char_save(:)
   integer, allocatable :: levels_save(:)
   integer, allocatable :: i_save(:,:), j_save(:,:), k_save(:,:)
   integer, allocatable :: obs_global_index_save(:)
   real, allocatable :: lat_save(:,:), lon_save(:,:), elv_save(:)
   real, allocatable :: x_save(:,:), y_save(:,:), dx_save(:,:), dxm_save(:,:)
   real, allocatable :: dy_save(:,:), dym_save(:,:), dz_save(:,:), dzm_save(:,:)
   real, allocatable :: zk_save(:,:)
   logical, allocatable :: proc_domain_save(:,:), thinned_save(:,:)

   old_nlocal = iv%info(radar)%nlocal
   max_lev = max(1, iv%info(radar)%max_lev)
   iv%info(radar)%max_lev = max_lev

   if (old_nlocal < 1) then
      allocate(iv%radar(1:new_nlocal))

      allocate(iv%info(radar)%name(new_nlocal))
      allocate(iv%info(radar)%platform(new_nlocal))
      allocate(iv%info(radar)%id(new_nlocal))
      allocate(iv%info(radar)%date_char(new_nlocal))
      allocate(iv%info(radar)%levels(new_nlocal))
      allocate(iv%info(radar)%lat(max_lev,new_nlocal))
      allocate(iv%info(radar)%lon(max_lev,new_nlocal))
      allocate(iv%info(radar)%elv(new_nlocal))
      allocate(iv%info(radar)%x(kms:kme,new_nlocal))
      allocate(iv%info(radar)%y(kms:kme,new_nlocal))
      allocate(iv%info(radar)%i(kms:kme,new_nlocal))
      allocate(iv%info(radar)%j(kms:kme,new_nlocal))
      allocate(iv%info(radar)%dx(kms:kme,new_nlocal))
      allocate(iv%info(radar)%dxm(kms:kme,new_nlocal))
      allocate(iv%info(radar)%dy(kms:kme,new_nlocal))
      allocate(iv%info(radar)%dym(kms:kme,new_nlocal))
      allocate(iv%info(radar)%k(max_lev,new_nlocal))
      allocate(iv%info(radar)%dz(max_lev,new_nlocal))
      allocate(iv%info(radar)%dzm(max_lev,new_nlocal))
      allocate(iv%info(radar)%zk(max_lev,new_nlocal))
      allocate(iv%info(radar)%proc_domain(max_lev,new_nlocal))
      allocate(iv%info(radar)%obs_global_index(new_nlocal))
      allocate(iv%info(radar)%thinned(max_lev,new_nlocal))

      iv%info(radar)%name(:) = ' '
      iv%info(radar)%platform(:) = ' '
      iv%info(radar)%id(:) = ' '
      iv%info(radar)%date_char(:) = ' '
      iv%info(radar)%levels(:) = 0
      iv%info(radar)%lat(:,:) = missing_r
      iv%info(radar)%lon(:,:) = missing_r
      iv%info(radar)%elv(:) = missing_r
      iv%info(radar)%x(:,:) = missing_r
      iv%info(radar)%y(:,:) = missing_r
      iv%info(radar)%i(:,:) = 0
      iv%info(radar)%j(:,:) = 0
      iv%info(radar)%dx(:,:) = 0.0
      iv%info(radar)%dxm(:,:) = 0.0
      iv%info(radar)%dy(:,:) = 0.0
      iv%info(radar)%dym(:,:) = 0.0
      iv%info(radar)%k(:,:) = 0
      iv%info(radar)%dz(:,:) = 0.0
      iv%info(radar)%dzm(:,:) = 0.0
      iv%info(radar)%zk(:,:) = missing_r
      iv%info(radar)%proc_domain(:,:) = .false.
      iv%info(radar)%obs_global_index(:) = 0
      iv%info(radar)%thinned(:,:) = .false.

      return
   end if

   allocate(radar_save(old_nlocal))
   radar_save(1:old_nlocal) = iv%radar(1:old_nlocal)
   deallocate(iv%radar)
   allocate(iv%radar(1:new_nlocal))
   iv%radar(1:old_nlocal) = radar_save(1:old_nlocal)
   deallocate(radar_save)

   allocate(name_save(old_nlocal))
   name_save = iv%info(radar)%name
   deallocate(iv%info(radar)%name)
   allocate(iv%info(radar)%name(new_nlocal))
   iv%info(radar)%name(:) = ' '
   iv%info(radar)%name(1:old_nlocal) = name_save
   deallocate(name_save)

   allocate(platform_save(old_nlocal))
   platform_save = iv%info(radar)%platform
   deallocate(iv%info(radar)%platform)
   allocate(iv%info(radar)%platform(new_nlocal))
   iv%info(radar)%platform(:) = ' '
   iv%info(radar)%platform(1:old_nlocal) = platform_save
   deallocate(platform_save)

   allocate(id_save(old_nlocal))
   id_save = iv%info(radar)%id
   deallocate(iv%info(radar)%id)
   allocate(iv%info(radar)%id(new_nlocal))
   iv%info(radar)%id(:) = ' '
   iv%info(radar)%id(1:old_nlocal) = id_save
   deallocate(id_save)

   allocate(date_char_save(old_nlocal))
   date_char_save = iv%info(radar)%date_char
   deallocate(iv%info(radar)%date_char)
   allocate(iv%info(radar)%date_char(new_nlocal))
   iv%info(radar)%date_char(:) = ' '
   iv%info(radar)%date_char(1:old_nlocal) = date_char_save
   deallocate(date_char_save)

   allocate(levels_save(old_nlocal))
   levels_save = iv%info(radar)%levels
   deallocate(iv%info(radar)%levels)
   allocate(iv%info(radar)%levels(new_nlocal))
   iv%info(radar)%levels(:) = 0
   iv%info(radar)%levels(1:old_nlocal) = levels_save
   deallocate(levels_save)

   allocate(lat_save(max_lev,old_nlocal))
   lat_save = iv%info(radar)%lat
   deallocate(iv%info(radar)%lat)
   allocate(iv%info(radar)%lat(max_lev,new_nlocal))
   iv%info(radar)%lat(:,:) = missing_r
   iv%info(radar)%lat(:,1:old_nlocal) = lat_save
   deallocate(lat_save)

   allocate(lon_save(max_lev,old_nlocal))
   lon_save = iv%info(radar)%lon
   deallocate(iv%info(radar)%lon)
   allocate(iv%info(radar)%lon(max_lev,new_nlocal))
   iv%info(radar)%lon(:,:) = missing_r
   iv%info(radar)%lon(:,1:old_nlocal) = lon_save
   deallocate(lon_save)

   allocate(elv_save(old_nlocal))
   elv_save = iv%info(radar)%elv
   deallocate(iv%info(radar)%elv)
   allocate(iv%info(radar)%elv(new_nlocal))
   iv%info(radar)%elv(:) = missing_r
   iv%info(radar)%elv(1:old_nlocal) = elv_save
   deallocate(elv_save)

   allocate(x_save(kms:kme,old_nlocal))
   x_save = iv%info(radar)%x
   deallocate(iv%info(radar)%x)
   allocate(iv%info(radar)%x(kms:kme,new_nlocal))
   iv%info(radar)%x(:,:) = missing_r
   iv%info(radar)%x(:,1:old_nlocal) = x_save
   deallocate(x_save)

   allocate(y_save(kms:kme,old_nlocal))
   y_save = iv%info(radar)%y
   deallocate(iv%info(radar)%y)
   allocate(iv%info(radar)%y(kms:kme,new_nlocal))
   iv%info(radar)%y(:,:) = missing_r
   iv%info(radar)%y(:,1:old_nlocal) = y_save
   deallocate(y_save)

   allocate(i_save(kms:kme,old_nlocal))
   i_save = iv%info(radar)%i
   deallocate(iv%info(radar)%i)
   allocate(iv%info(radar)%i(kms:kme,new_nlocal))
   iv%info(radar)%i(:,:) = 0
   iv%info(radar)%i(:,1:old_nlocal) = i_save
   deallocate(i_save)

   allocate(j_save(kms:kme,old_nlocal))
   j_save = iv%info(radar)%j
   deallocate(iv%info(radar)%j)
   allocate(iv%info(radar)%j(kms:kme,new_nlocal))
   iv%info(radar)%j(:,:) = 0
   iv%info(radar)%j(:,1:old_nlocal) = j_save
   deallocate(j_save)

   allocate(dx_save(kms:kme,old_nlocal))
   dx_save = iv%info(radar)%dx
   deallocate(iv%info(radar)%dx)
   allocate(iv%info(radar)%dx(kms:kme,new_nlocal))
   iv%info(radar)%dx(:,:) = 0.0
   iv%info(radar)%dx(:,1:old_nlocal) = dx_save
   deallocate(dx_save)

   allocate(dxm_save(kms:kme,old_nlocal))
   dxm_save = iv%info(radar)%dxm
   deallocate(iv%info(radar)%dxm)
   allocate(iv%info(radar)%dxm(kms:kme,new_nlocal))
   iv%info(radar)%dxm(:,:) = 0.0
   iv%info(radar)%dxm(:,1:old_nlocal) = dxm_save
   deallocate(dxm_save)

   allocate(dy_save(kms:kme,old_nlocal))
   dy_save = iv%info(radar)%dy
   deallocate(iv%info(radar)%dy)
   allocate(iv%info(radar)%dy(kms:kme,new_nlocal))
   iv%info(radar)%dy(:,:) = 0.0
   iv%info(radar)%dy(:,1:old_nlocal) = dy_save
   deallocate(dy_save)

   allocate(dym_save(kms:kme,old_nlocal))
   dym_save = iv%info(radar)%dym
   deallocate(iv%info(radar)%dym)
   allocate(iv%info(radar)%dym(kms:kme,new_nlocal))
   iv%info(radar)%dym(:,:) = 0.0
   iv%info(radar)%dym(:,1:old_nlocal) = dym_save
   deallocate(dym_save)

   allocate(k_save(max_lev,old_nlocal))
   k_save = iv%info(radar)%k
   deallocate(iv%info(radar)%k)
   allocate(iv%info(radar)%k(max_lev,new_nlocal))
   iv%info(radar)%k(:,:) = 0
   iv%info(radar)%k(:,1:old_nlocal) = k_save
   deallocate(k_save)

   allocate(dz_save(max_lev,old_nlocal))
   dz_save = iv%info(radar)%dz
   deallocate(iv%info(radar)%dz)
   allocate(iv%info(radar)%dz(max_lev,new_nlocal))
   iv%info(radar)%dz(:,:) = 0.0
   iv%info(radar)%dz(:,1:old_nlocal) = dz_save
   deallocate(dz_save)

   allocate(dzm_save(max_lev,old_nlocal))
   dzm_save = iv%info(radar)%dzm
   deallocate(iv%info(radar)%dzm)
   allocate(iv%info(radar)%dzm(max_lev,new_nlocal))
   iv%info(radar)%dzm(:,:) = 0.0
   iv%info(radar)%dzm(:,1:old_nlocal) = dzm_save
   deallocate(dzm_save)

   allocate(zk_save(max_lev,old_nlocal))
   zk_save = iv%info(radar)%zk
   deallocate(iv%info(radar)%zk)
   allocate(iv%info(radar)%zk(max_lev,new_nlocal))
   iv%info(radar)%zk(:,:) = missing_r
   iv%info(radar)%zk(:,1:old_nlocal) = zk_save
   deallocate(zk_save)

   allocate(proc_domain_save(max_lev,old_nlocal))
   proc_domain_save = iv%info(radar)%proc_domain
   deallocate(iv%info(radar)%proc_domain)
   allocate(iv%info(radar)%proc_domain(max_lev,new_nlocal))
   iv%info(radar)%proc_domain(:,:) = .false.
   iv%info(radar)%proc_domain(:,1:old_nlocal) = proc_domain_save
   deallocate(proc_domain_save)

   allocate(obs_global_index_save(old_nlocal))
   obs_global_index_save = iv%info(radar)%obs_global_index
   deallocate(iv%info(radar)%obs_global_index)
   allocate(iv%info(radar)%obs_global_index(new_nlocal))
   iv%info(radar)%obs_global_index(:) = 0
   iv%info(radar)%obs_global_index(1:old_nlocal) = obs_global_index_save
   deallocate(obs_global_index_save)

   allocate(thinned_save(max_lev,old_nlocal))
   thinned_save = iv%info(radar)%thinned
   deallocate(iv%info(radar)%thinned)
   allocate(iv%info(radar)%thinned(max_lev,new_nlocal))
   iv%info(radar)%thinned(:,:) = .false.
   iv%info(radar)%thinned(:,1:old_nlocal) = thinned_save
   deallocate(thinned_save)

end subroutine da_extend_radar_rqv_bayes_pseudo


subroutine da_init_radar_rqv_bayes_pseudo(obs, max_lev)

   !---------------------------------------------------------------------------
   ! Purpose: Allocate and initialize one rqv-only pseudo radar profile.
   !---------------------------------------------------------------------------

   implicit none

   type(radar_type), intent(inout) :: obs
   integer,          intent(in)    :: max_lev

   allocate(obs%model_p(1:max_lev))
   allocate(obs%model_rho(1:max_lev))
   allocate(obs%model_qrn(1:max_lev))
   allocate(obs%model_qcl(1:max_lev))
   allocate(obs%model_qci(1:max_lev))
   allocate(obs%model_qsn(1:max_lev))
   allocate(obs%model_qgr(1:max_lev))
   allocate(obs%model_zmm(1:max_lev))
   allocate(obs%height(1:max_lev))
   allocate(obs%height_qc(1:max_lev))
   allocate(obs%rv(1:max_lev))
   allocate(obs%rf(1:max_lev))
   allocate(obs%zmm(1:max_lev))

   obs%model_p(:) = missing_r
   obs%model_rho(:) = missing_r
   obs%model_qrn(:) = missing_r
   obs%model_qcl(:) = missing_r
   obs%model_qci(:) = missing_r
   obs%model_qsn(:) = missing_r
   obs%model_qgr(:) = missing_r
   obs%model_zmm(:) = missing_r
   obs%model_ps = missing_r
   obs%height(:) = missing_r
   obs%height_qc(:) = missing

   call da_init_radar_rqv_bayes_field(obs%rv)
   call da_init_radar_rqv_bayes_field(obs%rf)
   call da_init_radar_rqv_bayes_field(obs%zmm)

   if (use_radar_rhv) then
      allocate(obs%rrn(1:max_lev))
      allocate(obs%rcl(1:max_lev))
      allocate(obs%rci(1:max_lev))
      allocate(obs%rsn(1:max_lev))
      allocate(obs%rgr(1:max_lev))
      allocate(obs%rrno(1:max_lev))
      allocate(obs%rclo(1:max_lev))
      allocate(obs%rcio(1:max_lev))
      allocate(obs%rsno(1:max_lev))
      allocate(obs%rgro(1:max_lev))
      call da_init_radar_rqv_bayes_field(obs%rrn)
      call da_init_radar_rqv_bayes_field(obs%rcl)
      call da_init_radar_rqv_bayes_field(obs%rci)
      call da_init_radar_rqv_bayes_field(obs%rsn)
      call da_init_radar_rqv_bayes_field(obs%rgr)
      obs%rrno(:) = missing_r
      obs%rclo(:) = missing_r
      obs%rcio(:) = missing_r
      obs%rsno(:) = missing_r
      obs%rgro(:) = missing_r
   end if

   if (use_radar_rqv) then
      allocate(obs%rqv(1:max_lev))
      allocate(obs%rqvo(1:max_lev))
      call da_init_radar_rqv_bayes_field(obs%rqv)
      obs%rqvo(:) = missing_r
   end if

end subroutine da_init_radar_rqv_bayes_pseudo


subroutine da_init_radar_rqv_bayes_field(field)

   !---------------------------------------------------------------------------
   ! Purpose: Initialize radar field_type arrays to missing.
   !---------------------------------------------------------------------------

   implicit none

   type(field_type), intent(inout) :: field(:)

   field(:)%inv = missing_r
   field(:)%qc = missing_data
   field(:)%error = missing_r
   field(:)%sens = missing_r
   field(:)%imp = missing_r

end subroutine da_init_radar_rqv_bayes_field


subroutine da_setup_radar_rqv_bayes_rf(ra, rs, rg, tc, rho, rf)

   !---------------------------------------------------------------------------
   ! Purpose: Local copy of radar reflectivity used during setup-time screening.
   !---------------------------------------------------------------------------

   implicit none

   real, intent(in)  :: ra, rs, rg, tc, rho
   real, intent(out) :: rf

   real :: zrr, zds, zws, zgr
   real :: zerr, zews, zeds, zegr, rze

   zrr = 3.63*1.00e+9
   zws = 4.26*1.00e+11
   zds = 9.80*1.00e+8
   zgr = 4.33*1.00e+10

   zerr = zrr*(rho*max(0.0,ra))**1.75
   zews = zws*(rho*max(0.0,rs))**1.75
   zeds = zds*(rho*max(0.0,rs))**1.75
   zegr = zgr*(rho*max(0.0,rg))**1.75

   if (tc < 0.0) then
      rze = zerr + zeds + zegr
   else
      rze = zerr + zews + zegr
   end if

   rze = max(rze, 1.0e-20)
   rf = 10.0*log10(rze)

end subroutine da_setup_radar_rqv_bayes_rf
