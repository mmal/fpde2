program radauIIA_RUN

   use logger_module
   use constants_module
   use class_ode_stepper_radauIIA

   type(ode_stepper_radauIIA) :: s
   integer :: err

   integer :: i,j

   call set_log_level(FPDE_LOG_DEBUG)

   s % dim = 2

   call s%log(FPDE_LOG_INFO, "call s%init")
   call s%init(err)
   if ( err .ne. FPDE_STATUS_OK ) then
      call s%log(FPDE_LOG_ERROR, "init test failed")
      stop
   else
      call s%log(FPDE_LOG_INFO, "init test passed")
   end if

   call s%log(FPDE_LOG_INFO, "call s%bt%test")
   call s%bt%test(err)
   if ( err .ne. FPDE_STATUS_OK ) then
      call s%log(FPDE_LOG_ERROR, "Butcher tableu test failed")
      stop
   else
      call s%log(FPDE_LOG_INFO, "Butcher tableu test passed")
   end if

   print *, 'AI matrix'
   do i = 1, s%dim * s%stages
      print 11, (s % AI(i,j), j = 1, s%dim * s%stages)
   end do
11 FORMAT (2000(f12.5))

   print *, 'mu matrix'
   do i = 1, 2*s%dim
      print 11, (s % mu(i,j), j = 1,2*s%dim)
   end do

   print *, 'TT matrix'
   do i = 1, 3
      print 11, (s % TT(i,j), j = 1,3)
   end do

   print *, 'TTI matrix'
   do i = 1, 3
      print 11, (s % TTI(i,j), j = 1,3)
   end do

   print *, 'TT_cr_ID  matrix'
   do i = 1, 3*s%dim
      print 11, (s % TT_cr_ID(i,j), j = 1,3*s%dim)
   end do


   s % TT_cr_ID = matmul(s % TTI_cr_ID, s % TT_cr_ID)

   print *, 's % TT_cr_ID'
   do i = 1, 3*s%dim
      print 11, (s % TT_cr_ID(i,j), j = 1,3*s%dim)
   end do



   call s%log(FPDE_LOG_INFO, "call s%free")
   call s%free(err)
   if ( err .ne. FPDE_STATUS_OK ) then
      call s%log(FPDE_LOG_ERROR, "free test failed")
      stop
   else
      call s%log(FPDE_LOG_INFO, "free test passed")
   end if

end program radauIIA_RUN
