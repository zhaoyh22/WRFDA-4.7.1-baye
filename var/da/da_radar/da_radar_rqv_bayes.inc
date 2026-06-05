subroutine da_radar_rqv_bayes (bayes_type, obs_rf, bg_rf, bg_rh, i_ob, j_ob, k_ob, &
                               bg_rf_g, bg_rh_g, pseudo_rh, bayes_ok)

   !-----------------------------------------------------------------------
   ! Purpose: Estimate pseudo RH for radar rqv assimilation using local
   !          background-reflectivity/RH neighborhoods.
   !          bayes_type=1: weak echo Bayesian adjustment.
   !          bayes_type=2: clear/missing echo false-convection suppression.
   !-----------------------------------------------------------------------

   implicit none

   integer, intent(in) :: bayes_type
   integer, intent(in) :: i_ob, j_ob, k_ob
   real,    intent(in) :: obs_rf
   real,    intent(in) :: bg_rf
   real,    intent(in) :: bg_rh
   real,    intent(in) :: bg_rf_g(ids:ide,jds:jde,kds:kde)
   real,    intent(in) :: bg_rh_g(ids:ide,jds:jde,kds:kde)
   real,    intent(out) :: pseudo_rh
   logical, intent(out) :: bayes_ok

   integer :: i, j
   integer :: i_start, i_end, j_start, j_end
   integer :: radius, count
   real    :: weight, weight_sum, rh_sum
   real    :: delta_rf, delta_rh

   pseudo_rh = 0.0
   bayes_ok = .false.

   if (k_ob < kds .or. k_ob > kde) return

   select case (bayes_type)

   case (1)
      if (radar_rqv_bayes_sigma <= 0.0) return

      radius = max(0, radar_rqv_bayes_radius)
      i_start = max(ids, i_ob - radius)
      i_end   = min(ide, i_ob + radius)
      j_start = max(jds, j_ob - radius)
      j_end   = min(jde, j_ob + radius)

      rh_sum = 0.0
      weight_sum = 0.0

      do j = j_start, j_end
         do i = i_start, i_end
            if (abs(bg_rf_g(i,j,k_ob) - missing_r) <= 1.0) cycle
            if (abs(bg_rh_g(i,j,k_ob) - missing_r) <= 1.0) cycle
            weight = exp(-((bg_rf_g(i,j,k_ob)-obs_rf)**2)/(2.0*radar_rqv_bayes_sigma**2))
            rh_sum = rh_sum + weight*bg_rh_g(i,j,k_ob)
            weight_sum = weight_sum + weight
         end do
      end do

      if (weight_sum <= 0.0) return

      pseudo_rh = rh_sum / weight_sum
      pseudo_rh = amin1(100.0, amax1(0.0, pseudo_rh))

      delta_rf = obs_rf - bg_rf
      delta_rh = pseudo_rh - bg_rh
      bayes_ok = delta_rf*delta_rh > 0.0

   case (2)
      radius = max(0, radar_rqv_bayes_clear_radius)
      i_start = max(ids, i_ob - radius)
      i_end   = min(ide, i_ob + radius)
      j_start = max(jds, j_ob - radius)
      j_end   = min(jde, j_ob + radius)

      rh_sum = 0.0
      count = 0

      do j = j_start, j_end
         do i = i_start, i_end
            if (i == i_ob .and. j == j_ob) cycle
            if (abs(bg_rf_g(i,j,k_ob) - missing_r) <= 1.0) cycle
            if (abs(bg_rh_g(i,j,k_ob) - missing_r) <= 1.0) cycle
            if (bg_rf_g(i,j,k_ob) <= radar_rqv_bayes_noise .and. &
                bg_rh_g(i,j,k_ob) < bg_rh) then
               rh_sum = rh_sum + bg_rh_g(i,j,k_ob)
               count = count + 1
            end if
         end do
      end do

      if (count < radar_rqv_bayes_min_pts) return

      pseudo_rh = rh_sum / real(count)
      pseudo_rh = amin1(100.0, amax1(radar_rqv_bayes_rh_min, pseudo_rh))
      bayes_ok = pseudo_rh < bg_rh

   end select

end subroutine da_radar_rqv_bayes
