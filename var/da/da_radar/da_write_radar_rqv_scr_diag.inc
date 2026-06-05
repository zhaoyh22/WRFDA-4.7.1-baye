subroutine da_write_radar_rqv_scr_diag(it, grid, ob, iv, radar_bayes_rf, &
                                       scr_obs_echo, scr_assimilated)

   !-----------------------------------------------------------------------
   ! Purpose: Write SCR diagnostic masks for radar rqv Bayesian adjustment.
   !          scr_valid       : grid cells inside the radar diagnosis area.
   !          scr_candidate   : background strong echo with no/weak/missing obs.
   !          scr_assimilated : candidate cells accepted by Bayesian rqv QC.
   !-----------------------------------------------------------------------

   implicit none

   integer,       intent(in) :: it
   type(domain),  intent(in) :: grid
   type(y_type),  intent(in) :: ob
   type(iv_type), intent(in) :: iv
   real,          intent(in) :: radar_bayes_rf(ims:ime,jms:jme,kms:kme)
   logical,       intent(in) :: scr_obs_echo(ims:ime,jms:jme,kms:kme)
   real,          intent(in) :: scr_assimilated(ims:ime,jms:jme,kms:kme)

   integer :: i, j, k, n, ii, jj
   integer :: i_min, i_max, j_min, j_max
   integer :: i_min_l, i_max_l, j_min_l, j_max_l
   integer :: i_beg, i_end, j_beg, j_end, k_beg, k_end
   integer :: radius
   logical :: found_obs
   real, allocatable :: scr_valid(:,:,:)
   real, allocatable :: scr_candidate(:,:,:)
   real, allocatable :: scr_valid_g(:,:,:)
   real, allocatable :: scr_candidate_g(:,:,:)
   real, allocatable :: scr_assimilated_g(:,:,:)
#ifdef DM_PARALLEL
   integer :: found_obs_l, found_obs_g
#endif

   if (trace_use) call da_trace_entry("da_write_radar_rqv_scr_diag")

   if (it /= 1) then
      if (rootproc) write(unit=stdout,fmt='(A)') &
         ' SCR diagnostic skipped: only written for outer iteration 1'
      if (trace_use) call da_trace_exit("da_write_radar_rqv_scr_diag")
      return
   end if

#ifndef DM_PARALLEL
   if (iv%info(radar)%n2 < iv%info(radar)%n1) then
      if (rootproc) write(unit=stdout,fmt='(A)') &
         ' SCR diagnostic skipped: no radar observations are available'
      if (trace_use) call da_trace_exit("da_write_radar_rqv_scr_diag")
      return
   end if
#endif

   allocate(scr_valid(ims:ime,jms:jme,kms:kme))
   allocate(scr_candidate(ims:ime,jms:jme,kms:kme))
   allocate(scr_valid_g(ids:ide,jds:jde,kds:kde))
   allocate(scr_candidate_g(ids:ide,jds:jde,kds:kde))
   allocate(scr_assimilated_g(ids:ide,jds:jde,kds:kde))

   scr_valid(:,:,:) = 0.0
   scr_candidate(:,:,:) = 0.0
   scr_valid_g(:,:,:) = 0.0
   scr_candidate_g(:,:,:) = 0.0
   scr_assimilated_g(:,:,:) = 0.0

   i_min_l = ide
   i_max_l = ids
   j_min_l = jde
   j_max_l = jds
   found_obs = .false.

   do n = iv%info(radar)%n1, iv%info(radar)%n2
      if (iv%info(radar)%levels(n) < 1) cycle
      if (trim(iv%info(radar)%platform(n)) == 'RQV_BAYES') cycle

      ii = nint(iv%info(radar)%x(1,n))
      jj = nint(iv%info(radar)%y(1,n))
      if (ii < ids .or. ii > ide-1 .or. jj < jds .or. jj > jde-1) cycle

      found_obs = .true.
      i_min_l = min(i_min_l, ii)
      i_max_l = max(i_max_l, ii)
      j_min_l = min(j_min_l, jj)
      j_max_l = max(j_max_l, jj)
   end do

#ifdef DM_PARALLEL
   found_obs_l = 0
   if (found_obs) found_obs_l = 1
   found_obs_g = wrf_dm_sum_integer(found_obs_l)
   i_min = -wrf_dm_max_int(-i_min_l)
   i_max =  wrf_dm_max_int( i_max_l)
   j_min = -wrf_dm_max_int(-j_min_l)
   j_max =  wrf_dm_max_int( j_max_l)

   if (found_obs_g < 1) then
      if (rootproc) write(unit=stdout,fmt='(A)') &
         ' SCR diagnostic skipped: no non-RQV_BAYES radar observations found'
      deallocate(scr_valid)
      deallocate(scr_candidate)
      deallocate(scr_valid_g)
      deallocate(scr_candidate_g)
      deallocate(scr_assimilated_g)
      if (trace_use) call da_trace_exit("da_write_radar_rqv_scr_diag")
      return
   end if
#else
   if (.not. found_obs) then
      if (rootproc) write(unit=stdout,fmt='(A)') &
         ' SCR diagnostic skipped: no non-RQV_BAYES radar observations found'
      deallocate(scr_valid)
      deallocate(scr_candidate)
      deallocate(scr_valid_g)
      deallocate(scr_candidate_g)
      deallocate(scr_assimilated_g)
      if (trace_use) call da_trace_exit("da_write_radar_rqv_scr_diag")
      return
   end if
   i_min = i_min_l
   i_max = i_max_l
   j_min = j_min_l
   j_max = j_max_l
#endif

   radius = max(0, radar_rqv_bayes_clear_radius)
   i_min = max(ids, i_min-radius)
   i_max = min(ide-1, i_max+radius)
   j_min = max(jds, j_min-radius)
   j_max = min(jde-1, j_max+radius)

   i_beg = max(its, ids, i_min)
   i_end = min(ite, ide-1, ime-1, i_max)
   j_beg = max(jts, jds, j_min)
   j_end = min(jte, jde-1, jme-1, j_max)
   k_beg = max(kts, kds)
   k_end = min(kte, kds+mkz-1, kde, kme)

   if (i_beg <= i_end .and. j_beg <= j_end .and. k_beg <= k_end) then
      do k = k_beg, k_end
         do j = j_beg, j_end
            do i = i_beg, i_end
               scr_valid(i,j,k) = 1.0
               if (radar_bayes_rf(i,j,k) >= radar_rqv_bayes_bg_thresh .and. &
                   .not. scr_obs_echo(i,j,k)) then
                  scr_candidate(i,j,k) = 1.0
               end if
            end do
         end do
      end do
   end if

   call da_patch_to_global(grid, scr_valid, scr_valid_g, mkz)
   call da_patch_to_global(grid, scr_candidate, scr_candidate_g, mkz)
   call da_patch_to_global(grid, scr_assimilated, scr_assimilated_g, mkz)

   if (rootproc) then
      write(unit=stdout,fmt='(A)') ' Writing SCR diagnostic file before minimization'
      call da_write_radar_rqv_scr_diag_nc(iv%time, scr_valid_g, &
         scr_candidate_g, scr_assimilated_g)
   end if

   deallocate(scr_valid)
   deallocate(scr_candidate)
   deallocate(scr_valid_g)
   deallocate(scr_candidate_g)
   deallocate(scr_assimilated_g)

   if (trace_use) call da_trace_exit("da_write_radar_rqv_scr_diag")

end subroutine da_write_radar_rqv_scr_diag


subroutine da_write_radar_rqv_scr_diag_nc(time_index, scr_valid, &
                                          scr_candidate, scr_assimilated)

   !-----------------------------------------------------------------------
   ! Purpose: Root-process NetCDF writer for SCR diagnostic masks.
   !-----------------------------------------------------------------------

   implicit none

   include 'netcdf.inc'

   integer, intent(in) :: time_index
   real,    intent(in) :: scr_valid(ids:ide,jds:jde,kds:kde)
   real,    intent(in) :: scr_candidate(ids:ide,jds:jde,kds:kde)
   real,    intent(in) :: scr_assimilated(ids:ide,jds:jde,kds:kde)

   integer :: ncid, rcode
   integer :: dim_we, dim_sn, dim_bt
   integer :: dimids(3)
   integer :: var_valid, var_candidate, var_assimilated
   integer :: nx, ny, nz
   integer, allocatable :: mask(:,:,:)
   character(len=filename_len) :: filename
   character(len=160) :: attr

   filename = 'scr_diag.nc'

   nx = ide - ids + 1
   ny = jde - jds + 1
   nz = kde - kds + 1

   rcode = nf_create(trim(filename), nf_clobber, ncid)
   if (rcode /= nf_noerr) then
      write(unit=stdout,fmt='(A,A,A,I6)') ' Could not create ', trim(filename), &
         ', NetCDF status = ', rcode
      return
   end if

   rcode = nf_def_dim(ncid, 'west_east', nx, dim_we)
   rcode = nf_def_dim(ncid, 'south_north', ny, dim_sn)
   rcode = nf_def_dim(ncid, 'bottom_top', nz, dim_bt)

   dimids(1) = dim_we
   dimids(2) = dim_sn
   dimids(3) = dim_bt

   rcode = nf_def_var(ncid, 'scr_valid', nf_byte, 3, dimids, var_valid)
   attr = '1 where the grid cell is inside the SCR diagnosis area'
   rcode = nf_put_att_text(ncid, var_valid, 'description', len_trim(attr), attr)

   rcode = nf_def_var(ncid, 'scr_candidate', nf_byte, 3, dimids, var_candidate)
   attr = '1 where background reflectivity exceeds threshold and observation is no/weak/missing echo'
   rcode = nf_put_att_text(ncid, var_candidate, 'description', len_trim(attr), attr)

   rcode = nf_def_var(ncid, 'scr_assimilated', nf_byte, 3, dimids, var_assimilated)
   attr = '1 where a candidate SCR cell is accepted by Bayesian rqv QC and assimilated'
   rcode = nf_put_att_text(ncid, var_assimilated, 'description', len_trim(attr), attr)

   attr = trim(analysis_date)
   rcode = nf_put_att_text(ncid, nf_global, 'analysis_date', len_trim(attr), attr)
   write(attr,fmt='(I8)') time_index
   rcode = nf_put_att_text(ncid, nf_global, 'time_index', len_trim(attr), attr)
   attr = 'SCR: bg reflectivity >= radar_rqv_bayes_bg_thresh and obs no/weak/missing echo'
   rcode = nf_put_att_text(ncid, nf_global, 'scr_definition', len_trim(attr), attr)

   rcode = nf_enddef(ncid)
   if (rcode /= nf_noerr) then
      write(unit=stdout,fmt='(A,A,A,I6)') ' Could not define ', trim(filename), &
         ', NetCDF status = ', rcode
      rcode = nf_close(ncid)
      return
   end if

   allocate(mask(ids:ide,jds:jde,kds:kde))

   mask = nint(scr_valid)
   rcode = nf_put_var_int(ncid, var_valid, mask)
   mask = nint(scr_candidate)
   rcode = nf_put_var_int(ncid, var_candidate, mask)
   mask = nint(scr_assimilated)
   rcode = nf_put_var_int(ncid, var_assimilated, mask)

   deallocate(mask)

   rcode = nf_close(ncid)
   if (rcode == nf_noerr) then
      write(unit=stdout,fmt='(2A)') ' Wrote SCR diagnostic file: ', trim(filename)
   else
      write(unit=stdout,fmt='(A,A,A,I6)') ' Could not close ', trim(filename), &
         ', NetCDF status = ', rcode
   end if

end subroutine da_write_radar_rqv_scr_diag_nc
