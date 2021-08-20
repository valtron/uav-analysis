	program new_fdm
c	This theory is from the books "Aircraft Control and Simulation" 3rd edition,
c	Stevens, Lewis, and Johnson, 2016, (some code on page 182-183), and from
c	"Introduction to Helicopter and Tiltrotor Flight Simulation" 2nd edition,
c	by Dreier. Modified by Walker.
c	This subroutine fderiv provide the derivatives of the state vector,
c	defined in the "State and control variables" intiatlization below.
c	The best description of the internal memory is looking at include_parts.h
c	and the comments at the top of subroutine fderiv below.
c	James Walker  16 March 2021
c	4/12/2021
c
c	6/7/2021  Added the A, B, and G matrices for the linearization of controls.

c	The internals of this code are in MKS.  Inputs allow for different length
c	and mass units.  However, internals assume amps, watts, ohms, for the
c	electrical end of things, and right now there are no provisions for
c	adjust those units.  Thus, include_parts includes a modified for length
c	and mass to allow chaing those in the length and mass inputs.  Also,
c	input angles are assumed to be in degrees; internally they are changed
c	to radiana.  If they are input in radians it is necessary to set that
c	conversion factor (located in include_parts.h) to 1.


c	On top of the man physics routine fderiv there is integration code (fly)
c	and a static "trim" optimization code that finds the control vaviables
c	for staic flight.

c	Future work will also answer questions like this
c	1.  How long/how much power to get to 500 ft.
c	    a) Vertically;
c	    b) Their path, after 10 ft vertically.
c	2.  What is its cruise vs. power - i.e., given a flight speed, what
c	    is the optimal configuration for least power/greatest distance.
c	    One would also get out the trim control settings for each speed.
c	    (Completed: #2 is solved by trim.f)
c	3.  Let's follow a pre-defined path.  A take off, big turn, etc.
c	4.  Can it do aerobatics?  A predefined loop?  It seems we need to
c	    let it decide how to do it.
c	5.  Beyond these, I need to find the settings for the autopilot -
c	    this is a big need that somehow I have to address.

c	We need to think about input dimensions.

	implicit none

	integer xdim, udim, i

	double precision x(75), xd(75), uc(75)  ! dummy max dimensions, since we do not yet know number of controls

	include 'include_parts.h'

c	write (6,*) '  This is the output file for the SwRI flight'
c     .       ,' dynamics model.'

c	read in from an input file
	call read_input(x,uc,aircraft,wing,propeller,battery
     .                          ,control)


c	Now we initialize the controls based on information about them.
c	call controls_init(x,uc,aircraft,wing,propeller,batter,
c     .                          ,control)


c	Based on the requested action, we have a few different paths.
c	One of them is to fly the object.

c	We now set up the dimensions for the state vector.  There is some
c	careful counting here because there are lots of possibilities,
c	And lots of links between things.

	udim = aircraft%num_controls  ! this variable is the maximum number of controls
c	xdim = 13 + udim
	xdim = 13 + aircraft%num_propellers + aircraft%num_batteries
c	Doing it this way means that we do not get output from
c	other things, such as we only get one output for ailerons
c	and flaps and we only get one motor output, which may not
c	be representative if we end up linking motors together with
c	the same controller.

c	We open the metrics output file
	open(aircraft%i_out_metrics,file='metrics.out',status='unknown')
c	We open the score output file
	open(aircraft%i_out_score,file='score.out',status='unknown')

c	We open the path file
	open(55,file='path.out',status='unknown')
c	We open the path file
	open(56,file='path2.out',status='unknown')

	if (aircraft%i_analysis_type.eq.1) then

c	   This is flying from initial conditions

	   call fly(xdim,x,xd,udim,uc,aircraft, wing, propeller,battery
     .          ,control)

	else if (aircraft%i_analysis_type.eq.2) then

c	   This is performing a trim analysis to U = x(1) forward speed,
c	   level steady flight.

	   call trim(xdim,x,xd,udim,uc,aircraft,wing,propeller,battery
     .          ,control)

	else if (aircraft%i_analysis_type.eq.3) then

c	   The autopilot and flight.
c	   This is performing a trim analysis to U = x(1) forward speed,
c	   level steady flight to give a sequence of control settings.

	   control%compute_A = .true.
	   call trim(xdim,x,xd,udim,uc,aircraft,wing,propeller,battery
     .          ,control)

c	   And now we fly.  The various settings are in the control type.

	   call fly(xdim,x,xd,udim,uc,aircraft, wing, propeller,battery
     .          ,control)


	endif

	close(aircraft%i_out_metrics)

	end

	subroutine read_input(x,uc,aircraft,wing,propeller,battery
     .                          ,control)
c	This is to read the input using the namelist file.
c	James Walker  26 March 2021
c	3/26/21

	implicit none

	include 'include_parts.h'

	double precision :: tolerance = 1.d-04     ! this is to do some checks on 1s.
	double precision ls, ms, as

	integer i, jj, k

	double precision x(75), uc(75)  ! dummy max dimensions, since we do not yet know number of controls

	double precision    Ixx,   Ixy,   Ixz,   Iyy,   Iyz,   Izz
	double precision InvIxx,InvIxy,InvIxz,InvIyy,InvIyz,InvIzz
	double precision Denom

	double precision omega, n, D, rho, thrust, power, torque, J
	double precision Ct, Cp, Kt, Kv, Rw, voltage, a, b, c, amps


	double precision alpha(-10:10), force_L(-10:10,10)
        double precision                force_D(-10:10,10)
	double precision u, qbarprime_w, C_Lw, C_Dw
	logical          iflinear
	double precision dC_Ldalpha,dC_Ddalpha

	double precision Ct11, Ct12, Ct21, Ct22
	double precision Cp11, Cp12, Cp21, Cp22
	double precision I0, motor_efficiency

	integer          jdw, jdw2, l, kpast, lpast, ierror
	double precision dxkm1, dxk, dylm1, dyl
	double precision omega0, omega1, ff, ffp, a1, b1, omega2, DeltaJ
	double precision J1, J2, Om1, Om2, Vpath
	double precision aa, bb, cc, dd, ee,   a0,   a2,   a3,  ff1, Cp2
	double precision amps1,amps2,torque_table
        double precision efficiency_ESC,MaxPower
	
	double precision xmin, xmax, ymin, ymax, zmin, zmax

	double precision uc1, temp
	
	namelist/aircraft_data/aircraft,wing,propeller,battery,control

c	This is the read of the input data, essentially all of it but
c	propeller file reads (right now).

	read (5,nml=aircraft_data)

c	Now that we've read all the data, we first do the unit conversions
	ls = aircraft% length_scale
	ms = aircraft% mass  _scale
	as = aircraft% angle _scale

	if (ls.le.0.d0) ls = 1.d0
	if (ms.le.0.d0) ms = 1.d0
	if (as.le.0.d0) as = deg_to_rad
	if (abs(ls-1.d0).lt.tolerance) ls = 1.d0
	if (abs(ms-1.d0).lt.tolerance) ms = 1.d0
	if (abs(as-1.d0).lt.tolerance) as = 1.d0
	if (abs(as-deg_to_rad).lt.tolerance) as = deg_to_rad


c
c       Aircraft properties
	aircraft%x_cm      = aircraft%x_cm      *    ls    ! (m)
	aircraft%y_cm      = aircraft%y_cm      *    ls    ! (m)
	aircraft%z_cm      = aircraft%z_cm      *    ls    ! (m)
        aircraft%Ixx       = aircraft%Ixx       * ms*ls**2 ! (kg m^2)
        aircraft%Ixy       = aircraft%Ixy       * ms*ls**2 ! (kg m^2)
	aircraft%Ixz       = aircraft%Ixz       * ms*ls**2 ! (kg m^2)
        aircraft%Iyy       = aircraft%Iyy       * ms*ls**2 ! (kg m^2)
        aircraft%Iyz       = aircraft%Iyz       * ms*ls**2 ! (kg m^2)
	aircraft%Izz       = aircraft%Izz       * ms*ls**2 ! (kg m^2)
        aircraft%X_fuseuu  = aircraft%X_fuseuu  *    ls**2 ! (m^2)
	aircraft%Y_fusevv  = aircraft%Y_fusevv  *    ls**2 ! (m^2)
	aircraft%Z_fuseww  = aircraft%Z_fuseww  *    ls**2 ! (m^2)
        aircraft%x_fuse    = aircraft%x_fuse    *    ls    ! (m)
	aircraft%y_fuse    = aircraft%y_fuse    *    ls    ! (m)
	aircraft%z_fuse    = aircraft%Z_fuse    *    ls    ! (m)
c
c	We aso inialize various items (or else there are defaults)
	aircraft%x _initial(1 ) = aircraft%x _initial(1 ) * ls  ! (m/s)
	aircraft%x _initial(2 ) = aircraft%x _initial(2 ) * ls  ! (m/s)
	aircraft%x _initial(3 ) = aircraft%x _initial(3 ) * ls  ! (m/s)
	aircraft%x _initial(4 ) = aircraft%x _initial(4 ) * as  ! (rad/s)
	aircraft%x _initial(5 ) = aircraft%x _initial(5 ) * as  ! (rad/s)
	aircraft%x _initial(6 ) = aircraft%x _initial(6 ) * as  ! (rad/s)
	aircraft%x _initial(11) = aircraft%x _initial(11) * ls  ! (m)
	aircraft%x _initial(12) = aircraft%x _initial(12) * ls  ! (m)
	aircraft%x _initial(13) = aircraft%x _initial(13) * ls  ! (m)
c	
	aircraft%xd_initial(1 ) = aircraft%xd_initial(1 ) * ls  ! (m/s^2)
	aircraft%xd_initial(2 ) = aircraft%xd_initial(2 ) * ls  ! (m/s^2)
	aircraft%xd_initial(3 ) = aircraft%xd_initial(3 ) * ls  ! (m/s^2)
	aircraft%xd_initial(4 ) = aircraft%xd_initial(4 ) * as  ! (rad/s^2)
	aircraft%xd_initial(5 ) = aircraft%xd_initial(5 ) * as  ! (rad/s^2)
	aircraft%xd_initial(6 ) = aircraft%xd_initial(6 ) * as  ! (rad/s^2)
	aircraft%xd_initial(11) = aircraft%xd_initial(11) * ls  ! (m/s)
	aircraft%xd_initial(12) = aircraft%xd_initial(12) * ls  ! (m/s)
	aircraft%xd_initial(13) = aircraft%xd_initial(13) * ls  ! (m/s)

	aircraft%Unwind    = aircraft%Unwind * ls ! (m/s) 
	aircraft%Vewind    = aircraft%Vewind * ls ! (m/s)
	aircraft%Wdwind    = aircraft%Wdwind * ls ! (m/s)
c

	do i = 1, aircraft%num_wings
	
	   wing(i)%surface_area = wing(i)%surface_area  * ls**2 !  (m^2) 
	   wing(i)%a            = wing(i)%a / as                ! (1/rad)
	   wing(i)%delta1_max   = wing(i)%delta1_max * as       ! (rad)
	   wing(i)%delta1_min   = wing(i)%delta1_min * as       ! (rad)
	   wing(i)%delta2_max   = wing(i)%delta2_max * as       ! (rad)
	   wing(i)%delta2_min   = wing(i)%delta2_min * as       ! (rad)
	   wing(i)%x            = wing(i)%x * ls                ! (m)
	   wing(i)%y            = wing(i)%y * ls                ! (m)
	   wing(i)%z            = wing(i)%z * ls                ! (m)

	enddo

	do i = 1, aircraft%num_propellers

	   propeller(i)%x       = propeller(i)%x      *    ls    ! (m) 
	   propeller(i)%y       = propeller(i)%y      *    ls    ! (m)
	   propeller(i)%z       = propeller(i)%z      *    ls    ! (m
	   propeller(i)%radius  = propeller(i)%radius *    ls    ! (m)
	   propeller(i)%Ir      = propeller(i)%Ir     * ms*ls**2 ! (kg m^2)

	enddo


c	We invert the inertia tensor, computed by REDUCE
	Ixx = aircraft%Ixx
	Ixy = aircraft%Ixy
	Ixz = aircraft%Ixz
	Iyy = aircraft%Iyy
	Iyz = aircraft%Iyz
	Izz = aircraft%Izz
	Denom = Ixx*Iyy*Izz - Ixx*Iyz**2 - Ixy**2*Izz -Ixz**2*Iyy
     .         + 2.d0*Ixy*Ixz*Iyz
	InvIxx = ( Iyy*Izz - Iyz**2 )/Denom
	InvIxy = (-Ixy*Izz + Ixz*Iyz)/Denom
	InvIxz = ( Ixy*Iyz - Ixz*Iyy)/Denom
	InvIyy = ( Ixx*Izz - Ixz**2 )/Denom
	InvIyz = (-Ixx*Iyz + Ixy*Ixz)/Denom
	InvIzz = ( Ixx*Iyy - Ixy**2 )/Denom
	aircraft%InvIxx = InvIxx
	aircraft%InvIxy = InvIxy
	aircraft%InvIxz = InvIxz
	aircraft%InvIyy = InvIyy
	aircraft%InvIyz = InvIyz
	aircraft%InvIzz = InvIzz

c	We initialize geometry for the wing normals.

	do i = 1, aircraft%num_wings

	   temp = wing(i)%nx**2 + wing(i)%ny**2 + wing(i)%nz**2
	   
	   if (temp.le.0.d0) then
c	      The default wing
	      wing(i)%nx =  0.d0
	      wing(i)%ny =  0.d0
	      wing(i)%nz = -1.d0
	      temp = 1.d0
	   endif

c	   Now we get into some more complicated coding, and these values are changed to reflect how 
c	   they are used.  The first issue is confirming that one of our two allowable cases
c	   is occurring.  We need either nx to be zero (tail) or ny to be zero (angle of incidence
c	   for the wing).  We do not allow a yaw angle (nz = 0) for the wing - no swept wings.
	   if (.not.((wing(i)%nx.eq.0.d0).or.(wing(i)%ny.eq.0.d0))) then
c	      We look for the smallest of nx or ny
	      if (abs(wing(i)%nx).lt.abs(wing(i)%ny)) then
                 wing(i)%nx = 0.d0
	      else
	         wing(i)%ny = 0.d0
	      endif
	      write (6,*) '   Warning: wing can only rotate '
     .                    ,'around X or Y axis. '
	      write (6,*) '   Now nx = ',wing(i)%nx,
     .                      ' and ny = ',wing(i)%ny
	   endif

c	   We repeat the initial step, since the normal may have now changed.
	   temp = wing(i)%nx**2 + wing(i)%ny**2 + wing(i)%nz**2
	   if (temp.le.0.d0) then
	      wing(i)%nx =  0.d0
	      wing(i)%ny =  0.d0
	      wing(i)%nz = -1.d0
	      temp = 1.d0  !  added 8/5/21
	   endif

c	   We normalize the wing normal.
	   temp = sqrt(temp)
	   wing(i)%nx = wing(i)%nx/temp
	   wing(i)%ny = wing(i)%ny/temp
	   wing(i)%nz = wing(i)%nz/temp

c	   Now there are two cases - incidence or stabilizer angle.
	   if (wing(i)%ny.eq.0.d0) then 
c	      There is an initial angle of incidence.
	      wing(i)%nx = min( 1.d0,wing(i)%nx)
	      wing(i)%nx = max(-1.d0,wing(i)%nx)
	      wing(i)%angle_of_incidence = - asin(wing(i)%nx)
c             The view here is that
c                                   ^ n = unit normal
c                        < _       |   so a positive incident angle has nx < 0 and nz < 0
c        incident angle (    - _  |                      (ny =  0)
c                  +X  <---------*--------->
c                                |
c                                |
c                                V +Z
c
c	      That is how we take care of the angle of incidence.  Because
c	      of that, we are now going to set the normal back to the regular
c	      one, but with a sign change because of how the coding is handled.
	      wing(i)%nx = 0.d0
	      wing(i)%ny = 0.d0
	      wing(i)%nz = 1.d0
	   else
c	      By the coding, we must have nx = 0, so we are rotating the stabilizer
	      wing(i)%angle_of_incidence = 0.d0
	      wing(i)%nx = 0.d0
	      wing(i)%nz = - wing(i)%nz   ! we change the sign, to conform to the coding in fderiv.
	   endif

	enddo
	      

c	We initialize some of the controller stuff

c	   We now compute the angle of attack, which includes an adjustment based
c          on the deflection of the ailerons and flaps.  We enforce the minimum 
c	   and maximum deflection of the control surfaces.
c
c	   delta_1
c     .      = min(
c     .            max(wing(i)%delta_amin, xcontrol(wing(i)%icontrol1))
c     .               ,wing(i)%delta_amax)
c	   delta_2
c     .      = min(
c     .            max(wing(i)%delta_amin, xcontrol(wing(i)%icontrol2))
c     .               ,wing(i)%delta_amax)
c
c	   delta_a = wing(i)%bias1*delta_1 + wing(i)%bias2*delta_2

	jj = 0

	do i = 1, aircraft%num_wings

c	   First we see if the channels are defined.
	   if (wing(i)%icontrol1.gt.0) wing(i)%ic1 = .true.
	   if (wing(i)%icontrol2.gt.0) wing(i)%ic2 = .true.

c	   Compute linear coefficients taking 0 <= uc <= 1 to 
c	   min and max allowable deflection, adjust for the biases,
c	   and set the initial deflection to zero.

	   if (wing(i)%ic1) then
	      wing(i)%ac1 = wing(i)%delta1_max - wing(i)%delta1_min
	      wing(i)%bc1 = wing(i)%delta1_min
	      if (wing(i)%ac1.ne.0.d0) then
	         uc(wing(i)%icontrol1) = - wing(i)%bc1/wing(i)%ac1
	      else
	         uc(wing(i)%icontrol1) = 0.d0
	         wing(i)%ac1 = 0.d0
	         wing(i)%bc1 = 0.d0
	      endif
	      wing(i)%ac1 = wing(i)%bias1 * wing(i)%ac1
	      wing(i)%bc1 = wing(i)%bias1 * wing(i)%bc1
	   endif

	   if (wing(i)%ic2) then
	      wing(i)%ac2 = wing(i)%delta2_max - wing(i)%delta2_min
	      wing(i)%bc2 = wing(i)%delta2_min
	      if (wing(i)%ac2.ne.0.d0) then
	         uc(wing(i)%icontrol2) = - wing(i)%bc2/wing(i)%ac2
	      else
	         uc(wing(i)%icontrol2) = 0.d0
	         wing(i)%ac2 = 0.d0
	         wing(i)%bc2 = 0.d0
	      endif
	      wing(i)%ac2 = wing(i)%bias2 * wing(i)%ac2
	      wing(i)%bc2 = wing(i)%bias2 * wing(i)%bc2
	   endif

c	   Add 2 or 1 control channels, depending on what is needed.  But
c	   the channels are chosen in the input.
	   jj = max(jj, wing(i)%icontrol1)
	   jj = max(jj, wing(i)%icontrol2)

c	   We now set an internal variable
	   if ((wing(i)%C_L0.eq.0.d0).and.(wing(i)%AoAatL0.ne.0.d0)) 
     .        wing(i)%C_L0 = - wing(i)%a * wing(i)%AoAatL0
	
	enddo

        aircraft%num_wing_controls = jj

c	We do batteries next because we need the voltage to set the propeller
c	control info.

c	Currently we assume battery info comes from input file, not battery
c	file  (4/12/21)

c	do i = 1, aircraft%num_batteries
c
cc	   We read the battery file.
c	   open (unit=50,file = battery(i)%fname, err=35
c     .           ,status = 'old')
c	   read (50,*) ! skipping header line
c
c	   read (50,*,err=35) battery(i)%num_cells, battery(i)%voltage
c     .           ,battery(i)%Capacity,     battery(i)%C_Peak
c     .           ,battery(i)%C_Continuous, battery(i)%Rm
c
c	   close (unit=50, err=35)
c
c	   goto 235
c
c 35	   write (6,*) ' Battery file open/reading/close error ', 
c     .                battery(i)%fname  ! currently ignored
c
c	   close (unit=50, err=235)
c
c 235	   continue
c
cc          NumCells,Voltage_V,Capacity_mAh,PeakDisRate_C,ContDisRate_C,Resistance_mO
cc          6,22.2,6000,150,75,13
c
c	enddo


	do i = 1, aircraft%num_propellers

c	   We read in the propeller data.  Currently the table
c	   has a header line stating the variables (J,Ct,Cp) and
c	   then the data points.  These are all dimensionless.
c	   J = advance ration, Ct = thrust coefficient, Cp = power coefficient
cold	   open (unit=50,file = propeller(i)%prop_fname, err=31,
cold     .           status = 'old')
cold	   read (50,*) ! skipping header line
cold	   do k = 1, 100  ! current max number of points allowed
cold	        read (50,*,err=31,end=32) propeller(i)%J(k)
cold     .               ,propeller(i)%Ct(k), propeller(i)%Cp(k)
cold	   enddo
cold 32	   continue
cold	   propeller(i)%prop_num_pts = k-1
cold	   close (unit=50, err=31)
cold
cold	   goto 231
cold
cold 31	   write (6,*) ' Propeller file open/reading/close error ', 
cold     .                propeller(i)%prop_fname
cold
cold	   close (unit=50, err=231)
cold
cold 231       continue


	   call propread(
     .       aircraft, wing, propeller, battery, control)


c
c	   Currently we are not reading the motor file - it is assumed
c	   to be in the input (4/12/21)

cc	   Now we read the motor file.  The items are
cc	   KV in rpm/volt, Kt in Nm/amp, MaxCurrent in amps
cc	   Idle current in amps, Max power in watts
cc	   Winding resistance in milli ohms.
cc
c	   open (unit=50,file = propeller(i)%motor_fname, err=33
c     .           ,status = 'old')
c	   read (50,*) ! skipping header line
c	   read (50,*,err=33) propeller(i)%KV, propeller(i)%KT,
c     .                 propeller(i)%I_max, propeller(i)%I_idle,
c     .                 propeller(i)%maxpower, propeller(i)%Rw
c
c	   close (unit=50, err=33)
c	
c	   goto 233
c
c 33	   write (6,*) ' Motor file open/reading/close error ', 
c     .                propeller(i)%motor_fname
c
c	   close (unit=50, err=233)
c
c 233       continue


c 	   To understand the relationship between the control
c	   and the motor, we need a battery voltage.
c	   The boost of the bc by 1% is due to the fact that if
c	   the voltage is exactly the value you can get a NaN
c	   coming out of our propeller routines.  This results
c	   in low (like 0.25 RPM) spin rates when uc = 0.
c	   We may need to address this at some point, if people
c	   really want to turn the motor off.
cold	   propeller(i)%bc = propeller(i)%Rw * propeller(i)%I_idle
	   propeller(i)%bc = propeller(i)%Rw * propeller(i)%I_idle
c     .                     * 1.01d0
	   propeller(i)%ac = battery(propeller(i)%ibattery)%voltage
     .                     - propeller(i)%bc 

	   if (propeller(i)%ac.le.0.d0) 
     .        write (6,*) '  Warning: battery voltage not high enough '

c	   Add the control channel.
	   jj = max(jj, propeller(i)%icontrol)

c	   We now look to see if there is initial power in the motor,
c	   because we want to set the appropriate motor speed if so
c	   (rather than have it rise up from zero, and thus enter
c	   a transient solution state), since we presume people
c	   are looking for a flight state if such a value is input.

c	   If we actually did this, this is where the coding would be.

	enddo

        aircraft%num_controls = jj

	write (6,*) '                                '
     .            ,'                     '      
     .            ,'                             '
     .            ,'          '
     .            ,' +----- Motor ----+  +---- Battery ---+'
	write (6,*) '    Motor #    omega     omega  '
     .            ,'  Voltage    Thrust  '
     .            ,'  Torque    Power    Current '
     .            ,'Efficiency'
     .            ,' Max Power  Max Cur  Peak Cur  Cont Cur'
	write (6,*)  '              (rad/s)    (RPM)  '
     .            ,'  (volts)      (N)       (Nm)   (watts)  '
     .            ,'  (amps)      (%)    '
     .            ,'(watts)    (amps)    (amps)    (amps)'


	rho = aircraft%rho   

	do i = 1, aircraft%num_propellers

c	   It is actually pretty complicated.
c	   We match the torque from the propeller file, assuming the 
c	   aircraft is at rest, with the electrical motor torque:
c
c	   Torque = Cp * rho* omega^2 * D^5 / (2pi)^3 = (Kt/Rw)*(voltage - omega/Kv)
c
c	   There are some unit issues here, as right now Kt is Newton meter / amp,
c	   Rw is ohms, and Kv is volt / rpm. 
c
c	   When we insert our maximum voltage, this is a quadratic which we 
c	   can solve to get the maximum omega.  We then compute the amperage
c	   and stuff.

	   D   = 2.d0*propeller(i)%radius
cold	   Cp  = propeller(i)%Cp(1) ! this is the approximation
cold	   Ct  = propeller(i)%Ct(1)
	   Cp  = propeller(i)%Cp(1,1) ! this is the approximation
	   Ct  = propeller(i)%Ct(1,1)
	   Kt  = propeller(i)%Kt
	   Kv  = twopi*propeller(i)%Kv/60.d0  ! unit change from rpm/volt to rad per sec/volt
	   Rw  = propeller(i)%Rw
	   I0  = propeller(i)%I_idle
	   voltage = battery(propeller(i)%ibattery)%voltage  ! max battery voltage

	   a   =   (Rw/Kt)*(Cp*rho*D**5/twopi**3)
	   b   =   1.d0/Kv
	   c   = - voltage + Rw*I0

	   omega = (-b+sqrt(b**2-4.d0*a*c))/(2.d0*a)

c	   Now we try a Newton method

	   do jdw = 1, 5  ! takes 5 iterates to converge, based on examples

	      n = omega/twopi  ! to get reveolutions per second from radians per second

	      J = Vpath / (n*D) ! J = V/nD

	      call table (propeller(i)%prop_num_pts_J, J, propeller(i)%J
     .                , dxkm1, dxk, k, kpast, ierror)
	      call table (propeller(i)%prop_num_pts_omega, omega
     .                , propeller(i)%omega
     .                , dylm1, dyl, l, lpast, ierror)
c	      I believe the first index is for J, the second index is for omega
	      J1   = propeller(i)%J    (k-1)
	      J2   = propeller(i)%J    (k  )
	      Om1  = propeller(i)%omega(l-1)
	      Om2  = propeller(i)%omega(l  )
	      Cp11 = propeller(i)%Cp(k-1,l-1)
	      Cp21 = propeller(i)%Cp(k-1,l  )
	      Cp22 = propeller(i)%Cp(k  ,l  )
	      Cp12 = propeller(i)%Cp(k  ,l-1)
	      aa = Cp11        - Cp21        + Cp22        - Cp12
	      bb =-Cp11   *J2  + Cp21   *J2  - Cp22   *J1  + Cp12    *J1
	      cc =-Cp11*Om2    + Cp21*Om1    - Cp22*Om1    + Cp12*Om2
	      dd = Cp11*Om2*J2 - Cp21*Om1*J2 + Cp22*Om1*J1 - Cp12*Om2*J1

c	      Same as above, only this a is a/Cp above
	      a   =   (Rw/Kt)*(rho*D**5/twopi**3) 
	      b   =   1.d0/Kv
	      c   = - voltage + Rw*I0

	      ee = a/((J2 - J1)*(Om2 - Om1))

	      a3 = ee*bb
	      a2 = ee*(aa*twopi*Vpath/D+dd)
	      a1 = ee* cc*twopi*Vpath/D + b
	      a0 = c

	      ff    =     a3*omega**3 +      a2*omega**2 + a1*omega + a0
	      ffp   =3.d0*a3*omega**2 + 2.d0*a2*omega    + a1
	      omega = omega - ff/ffp

	   Cp = ((omega-Om2)*(J-J2)*Cp11-(omega-Om1)*(J-J2)*Cp21
     .          +(omega-Om1)*(J-J1)*Cp22-(omega-Om2)*(J-J1)*Cp12)
     .          /((J2 - J1)*(Om2 - Om1))

	   enddo

	   n = omega/twopi 
	   J = Vpath / (n*D) ! J = V/nD
	   Ct11 = propeller(i)%Ct(k-1,l-1)
	   Ct21 = propeller(i)%Ct(k-1,l  )
	   Ct22 = propeller(i)%Ct(k  ,l  )
	   Ct12 = propeller(i)%Ct(k  ,l-1)
	   Ct = ((omega-Om2)*(J-J2)*Ct11-(omega-Om1)*(J-J2)*Ct21
     .          +(omega-Om1)*(J-J1)*Ct22-(omega-Om2)*(J-J1)*Ct12)
     .          /((J2 - J1)*(Om2 - Om1))
	   Cp = ((omega-Om2)*(J-J2)*Cp11-(omega-Om1)*(J-J2)*Cp21
     .          +(omega-Om1)*(J-J1)*Cp22-(omega-Om2)*(J-J1)*Cp12)
     .          /((J2 - J1)*(Om2 - Om1))
	   a   =   (Rw/Kt)*(Cp*rho*D**5/twopi**3)
	   b   =   1.d0/Kv
	   c   = - voltage + Rw*I0

	   omega = (-b+sqrt(b**2-4.d0*a*c))/(2.d0*a)

	   Thrust = Ct*rho*n**2*D**4      ! T = Ct*rho*n^2*D^4
	   Power  = Cp*rho*n**3*D**5      ! P = Cp*rho*n^3*D^5
           Torque = Power/omega  ! Q = P/motor omega
	
c	   write (6,*) '  Torque mechanical ', Torque

c	   write (6,*) '  Torque electrical ',(Kt/Rw)*(voltage - omega/Kv - Rw*I0)

	   amps = (voltage - omega/Kv)/Rw 

c	   Motor efficiency
	   motor_efficiency = Power/(voltage*amps)

c	   Power is electrical power
	   amps = amps/propeller(i)%efficiency_ESC  ! losses in the ESC
	   power = voltage*amps

	   write (6,10) ' Max Volt ',i,  omega, 60.d0*omega/twopi
     .         , voltage, Thrust
     .         , Torque, power, amps, motor_efficiency*100.d0
     .         , propeller(i)%maxpower, propeller(i)%I_max
     .         ,(battery(propeller(i)%ibattery)%Capacity/1000.d0)
     .          *battery(propeller(i)%ibattery)%C_Peak
     .         ,(battery(propeller(i)%ibattery)%Capacity/1000.d0)
     .          *battery(propeller(i)%ibattery)%C_Continuous

 10	   format(a10, i2, 12f10.2)

	enddo

c	Now we repeat the exercise and focus now on maximum power rating for the motor.

	do i = 1, aircraft%num_propellers

	   D   = 2.d0*propeller(i)%radius
	   Kt  = propeller(i)%Kt
	   Kv  = twopi*propeller(i)%Kv/60.d0  ! unit change from rpm/volt to rad per sec/volt
	   Rw  = propeller(i)%Rw
	   I0  = propeller(i)%I_idle
	   efficiency_ESC = propeller(i)%efficiency_ESC
	   MaxPower = propeller(i)%maxpower

c	   In this situation, we know the power that the motor can provide,
c	   which is given by I*voltage (I is current = amps).  Using the equation
c	   Voltage = omega/Kv + I*Rw and I*voltage = power we get a quadratic
c	   Voltage = 0.5*(omega/Kv+sqrt((omega/Kv)^2+4*efficiency_ESc*Power_motor*Rw))
c	   This is the equation we will solve with Torque_prop = Kt(I - I0)
c	   = (Kt/Rw)*(Voltage - omega/Kv - Rw*I0)  to determine omega.
c	   Now we need to find omega, assuming J = 0, to give the power
	   Torque_table = 0.d0
	   k = 2
	   l = 1
c	   Torque = sqrt(efficiency_ESC * MaxPower * Rw)  ! this was wrong since it is motor power we are looking at
	   Torque = sqrt(MaxPower * Rw)
	   do while (Torque_table.lt.Torque)
	      l     = l + 1
	      Om2   = propeller(i)%omega(    l  )
	      Cp21  = propeller(i)%Cp   (k-1,l  )
	      omega = Om2
	      Cp    = Cp21
	      n     = omega/twopi
	      Power = Cp*rho*n**3*D**5      ! P = Cp*rho*n^3*D^5
	      Torque_table = Power/omega    ! Q = P/motor omega
	      Torque = (Kt/Rw)*(-0.5d0*omega/Kv - Rw*I0
     .         + sqrt((0.5d0*omega/Kv)**2 + MaxPower*Rw))
c     .         + sqrt((0.5d0*omega/Kv)**2 + efficiency_ESC*MaxPower*Rw))
	   enddo
c	   Now we know we are in the right region; our goal is to find
c	   omega for the known torque.  We need to solve the equation.
	   Om1   = propeller(i)%omega(    l-1)
	   Cp11  = propeller(i)%Cp   (k-1,l-1)
	   a3    = (Cp21 - Cp11)/(Om2 - Om1)
	   a2    = (Cp11*Om2 - Cp21*Om1)/(Om2 - Om1)
	   a1    = 0.d0
	   Torque = (Kt/Rw)*(-0.5d0*omega/Kv - Rw*I0
     .         + sqrt((0.5d0*omega/Kv)**2 + MaxPower*Rw))
	   a0    = - Torque*(twopi**3)/(rho*D**5)
	   omega = 0.5d0*(Om1 + Om2)
	   do jdw = 1, 10
	      Torque = (Kt/Rw)*(-0.5d0*omega/Kv - Rw*I0
     .         + sqrt((0.5d0*omega/Kv)**2 + MaxPower*Rw))
	      a0    = - Torque*(twopi**3)/(rho*D**5)
	      ff    =     a3*omega**3 +      a2*omega**2 + a1*omega + a0
	      ffp   =3.d0*a3*omega**2 + 2.d0*a2*omega    + a1
     .              - ((twopi**3)/(rho*D**5))
     .              * (Kt/Rw)*(-0.5d0/Kv  + (0.25d0*omega/Kv**2)
     .         / sqrt((0.5d0*omega/Kv)**2 + MaxPower*Rw))
	      omega = omega - ff/ffp
              Torque_table = (a3*omega+a2)*rho*D**5*omega**2/twopi**3
	   enddo
c	   Now we know omega, so we can determine the voltage and current
	   Voltage = 0.5d0*omega/Kv
     .         + sqrt((0.5d0*omega/Kv)**2 + MaxPower*Rw)
	   amps = (Voltage - omega/Kv)/Rw

c	   To compute amps we need to make some assumptions.  We will assume 85%
c	   efficiency for the motor, then
	   motor_efficiency = Kt*(amps-I0)*omega/(amps*Voltage)
	   Ct = ((Om2 - omega)*Ct11 + (omega - Om1)*Ct21)
     .          /(Om2 - Om1)
	   Thrust = Ct*rho*n**2*D**4      ! T = Ct*rho*n^2*D^4

c	   write (6,*) ' Power ',power,torque_table*omega
c     .               ,Kt*(amps-I0)*omega,voltage*amps

	   amps = amps/efficiency_ESC
	   power = Voltage*amps
	
	   write (6,10) ' Max Power',i,  omega, 60.d0*omega/twopi
     .         , voltage, Thrust
     .         , Torque, power, amps, motor_efficiency*100.d0
     .         , propeller(i)%maxpower, propeller(i)%I_max
     .         ,(battery(propeller(i)%ibattery)%Capacity/1000.d0)
     .          *battery(propeller(i)%ibattery)%C_Peak
     .         ,(battery(propeller(i)%ibattery)%Capacity/1000.d0)
     .          *battery(propeller(i)%ibattery)%C_Continuous

	enddo

c	Now we repeat the exercise and focus now on maximum current rating for the motor.
           
	do i = 1, aircraft%num_propellers

	   D   = 2.d0*propeller(i)%radius
	   Kt  = propeller(i)%Kt
	   Kv  = twopi*propeller(i)%Kv/60.d0  ! unit change from rpm/volt to rad per sec/volt
	   Rw  = propeller(i)%Rw
	   I0  = propeller(i)%I_idle

c	   In this situation, we know the current, which means we know the torque
c	   The max current is knocked down by the ESC efficiency.
c	   amps = propeller(i)%efficiency_ESC * propeller(i)%I_max  ! wrong since this is motor amps
	   amps = propeller(i)%I_max
	   Torque = Kt*(amps - I0)  
c	   Now we need to find omega, assuming J = 0, to give the torque
	   Torque_table = 0.d0
	   k = 2
	   l = 1
	   do while (Torque_table.lt.Torque)
	      l     = l + 1
	      Om2   = propeller(i)%omega(    l  )
	      Cp21  = propeller(i)%Cp   (k-1,l  )
	      omega = Om2
	      Cp    = Cp21
	      n     = omega/twopi
	      Power = Cp*rho*n**3*D**5      ! P = Cp*rho*n^3*D^5
	      Torque_table = Power/omega    ! Q = P/motor omega
	   enddo
c	   Now we know we are in the right region; our goal is to find
c	   omega for the known torque.  It turns out to be a cubic to
c	   solve.
	   Om1   = propeller(i)%omega(    l-1)
	   Cp11  = propeller(i)%Cp   (k-1,l-1)
	   a3    = (Cp21 - Cp11)/(Om2 - Om1)
	   a2    = (Cp11*Om2 - Cp21*Om1)/(Om2 - Om1)
	   a1    = 0.d0
	   a0    = - Torque*(twopi**3)/(rho*D**5)
	   omega = 0.5d0*(Om1 + Om2)
	   do jdw = 1, 5
	      ff    =     a3*omega**3 +      a2*omega**2 + a1*omega + a0
	      ffp   =3.d0*a3*omega**2 + 2.d0*a2*omega    + a1
	      omega = omega - ff/ffp
c             Torque_table = (a3*omega+a2)*rho*D**5*omega**2/twopi**3
	   enddo
c	   Now we know omega, so we can determine the voltage
	   Voltage = Rw*amps + omega/Kv 

c	   amps1 = (voltage - omega/Kv)/Rw 
c	   amps2 = Torque/Kt + I0
c	   amps1 = amps1/propeller(i)%efficiency_ESC
c	   amps2 = amps2/propeller(i)%efficiency_ESC
c	   write (6,*) ' amps, amp1, amps2 ',amps,amps1,amps2
	   motor_efficiency = Kt*(amps-I0)*omega/(amps*Voltage)
	   Ct = ((Om2 - omega)*Ct11 + (omega - Om1)*Ct21)
     .          /(Om2 - Om1)
	   Thrust = Ct*rho*n**2*D**4      ! T = Ct*rho*n^2*D^4
	   amps = propeller(i)%I_max/propeller(i)%efficiency_ESC
	   power = Voltage*amps

	   write (6,10) ' Max Amps  ',i,  omega, 60.d0*omega/twopi
     .         , voltage, Thrust
     .         , Torque, power, amps, motor_efficiency*100.d0
     .         , propeller(i)%maxpower, propeller(i)%I_max
     .         ,(battery(propeller(i)%ibattery)%Capacity/1000.d0)
     .          *battery(propeller(i)%ibattery)%C_Peak
     .         ,(battery(propeller(i)%ibattery)%Capacity/1000.d0)
     .          *battery(propeller(i)%ibattery)%C_Continuous

	enddo


c	Now we look at printing out some motor properties, as a function
c	of voltage, to understand the behavior.
	
	if (aircraft%motordetails) then

	   write (6,*) '  '
	   write (6,*) '   This list details as uc goes from'
     .                ,' 0 to 1 (uc is last column)'
	   write (6,*) '    Motor #    omega     omega  '
     .            ,'  Voltage    Thrust  '
     .            ,'  Torque    Power    Current '
     .            ,'Efficiency      uc'
	   write (6,*)  '              (rad/s)    (RPM)  '
     .            ,'  (volts)      (N)       (Nm)   (watts)  '
     .            ,'  (amps)      (%)       (-)'

	do i = 1, aircraft%num_propellers

c	   We loop over a range of the control parameter to see what is happening
	   do jj = 0, 10

	      uc1 = dfloat(jj)/10.d0			 
	      voltage = propeller(i)%ac*uc1 + propeller(i)%bc

c	   It is actually pretty complicated.
c	   We match the torque from the propeller file, assuming the 
c	   aircraft is at rest, with the electrical motor torque:
c
c	   Torque = Cp * rho* omega^2 * D^5 / (2pi)^3 = (Kt/Rw)*(voltage - omega/Kv)
c
c	   There are some unit issues here, as right now Kt is Newton meter / amp,
c	   Rw is ohms, and Kv is volt / rpm. 
c
c	   When we insert our maximum voltage, this is a quadratic which we 
c	   can solve to get the maximum omega.  We then compute the amperage
c	   and stuff.

	   D   = 2.d0*propeller(i)%radius
cold	   Cp  = propeller(i)%Cp(1) ! this is the approximation
cold	   Ct  = propeller(i)%Ct(1)
	   Cp  = propeller(i)%Cp(1,1) ! this is the approximation
	   Ct  = propeller(i)%Ct(1,1)
	   Kt  = propeller(i)%Kt
	   Kv  = twopi*propeller(i)%Kv/60.d0  ! unit change from rpm/volt to rad per sec/volt
	   Rw  = propeller(i)%Rw
	   I0  = propeller(i)%I_idle
c	   voltage = battery(propeller(i)%ibattery)%voltage  ! max battery voltage

	   a   =   (Rw/Kt)*(Cp*rho*D**5/twopi**3)
	   b   =   1.d0/Kv
	   c   = - voltage + Rw*I0

	   omega = (-b+sqrt(b**2-4.d0*a*c))/(2.d0*a)

	   if (omega.gt.0.d0) then

c	   Now we try a Newton method

	   do jdw = 1, 5  ! takes 5 iterates to converge, based on examples

	      n = omega/twopi  ! to get reveolutions per second from radians per second

	      J = Vpath / (n*D) ! J = V/nD

	      call table (propeller(i)%prop_num_pts_J, J, propeller(i)%J
     .                , dxkm1, dxk, k, kpast, ierror)
	      call table (propeller(i)%prop_num_pts_omega, omega
     .                , propeller(i)%omega
     .                , dylm1, dyl, l, lpast, ierror)
c	      I believe the first index is for J, the second index is for omega
	      J1   = propeller(i)%J    (k-1)
	      J2   = propeller(i)%J    (k  )
	      Om1  = propeller(i)%omega(l-1)
	      Om2  = propeller(i)%omega(l  )
	      Cp11 = propeller(i)%Cp(k-1,l-1)
	      Cp21 = propeller(i)%Cp(k-1,l  )
	      Cp22 = propeller(i)%Cp(k  ,l  )
	      Cp12 = propeller(i)%Cp(k  ,l-1)
	      aa = Cp11        - Cp21        + Cp22        - Cp12
	      bb =-Cp11   *J2  + Cp21   *J2  - Cp22   *J1  + Cp12    *J1
	      cc =-Cp11*Om2    + Cp21*Om1    - Cp22*Om1    + Cp12*Om2
	      dd = Cp11*Om2*J2 - Cp21*Om1*J2 + Cp22*Om1*J1 - Cp12*Om2*J1

c	      Same as above, only this a is a/Cp above
	      a   =   (Rw/Kt)*(rho*D**5/twopi**3) 
	      b   =   1.d0/Kv
	      c   = - voltage + Rw*I0

	      ee = a/((J2 - J1)*(Om2 - Om1))

	      a3 = ee*bb
	      a2 = ee*(aa*twopi*Vpath/D+dd)
	      a1 = ee* cc*twopi*Vpath/D + b
	      a0 = c

	      ff    =     a3*omega**3 +      a2*omega**2 + a1*omega + a0
	      ffp   =3.d0*a3*omega**2 + 2.d0*a2*omega    + a1
	      omega = omega - ff/ffp

	   Cp = ((omega-Om2)*(J-J2)*Cp11-(omega-Om1)*(J-J2)*Cp21
     .          +(omega-Om1)*(J-J1)*Cp22-(omega-Om2)*(J-J1)*Cp12)
     .          /((J2 - J1)*(Om2 - Om1))

	   enddo

	   n = omega/twopi 
	   J = Vpath / (n*D) ! J = V/nD
	   Ct11 = propeller(i)%Ct(k-1,l-1)
	   Ct21 = propeller(i)%Ct(k-1,l  )
	   Ct22 = propeller(i)%Ct(k  ,l  )
	   Ct12 = propeller(i)%Ct(k  ,l-1)
	   Ct = ((omega-Om2)*(J-J2)*Ct11-(omega-Om1)*(J-J2)*Ct21
     .          +(omega-Om1)*(J-J1)*Ct22-(omega-Om2)*(J-J1)*Ct12)
     .          /((J2 - J1)*(Om2 - Om1))
	   Cp = ((omega-Om2)*(J-J2)*Cp11-(omega-Om1)*(J-J2)*Cp21
     .          +(omega-Om1)*(J-J1)*Cp22-(omega-Om2)*(J-J1)*Cp12)
     .          /((J2 - J1)*(Om2 - Om1))
	   a   =   (Rw/Kt)*(Cp*rho*D**5/twopi**3)
	   b   =   1.d0/Kv
	   c   = - voltage + Rw*I0

	   omega = (-b+sqrt(b**2-4.d0*a*c))/(2.d0*a)

	   Thrust = Ct*rho*n**2*D**4      ! T = Ct*rho*n^2*D^4
	   Power  = Cp*rho*n**3*D**5      ! P = Cp*rho*n^3*D^5
           Torque = Power/omega  ! Q = P/motor omega

	   else

	      Thrust = 0.d0
	      Power  = 0.d0
              Torque = 0.d0
	
	   endif

	
c	   write (6,*) '  Torque mechanical ', Torque

c	   write (6,*) '  Torque electrical ',(Kt/Rw)*(voltage - omega/Kv - Rw*I0)

	   amps = (voltage - omega/Kv)/Rw 

c	   Motor efficiency
	   motor_efficiency = Power/(voltage*amps)

c	   Power is electrical power
	   amps = amps/propeller(i)%efficiency_ESC  ! losses in the ESC
	   power = voltage*amps

	   write (6,10) '   Volt ',i,  omega, 60.d0*omega/twopi
     .         , voltage, Thrust
     .         , Torque, power, amps, motor_efficiency*100.d0
     .         , uc1


	   enddo ! to jj looking on uc1 values

	   write (6,*)

	enddo

	endif ! to aircraft%motordetails

	return

c	We now build a body drag table in the event that such values have not
c	been defined.  It turns out that I was getting tripped up by the unit
c	conversions, so now this algorithm with be initiated on any of the 
c	fuseuu, fusevv, or fuseww being less than -1.d-8 (and there are no checks
c	on the _fuse values themselves).

	if (aircraft%num_propellers.ge.1) then

c	if((aircraft%X_fuse+12345.d0).lt.tolerance) 
c     .              aircraft%X_fuse = aircraft%x_cm
c	if((aircraft%Y_fuse+12345.d0).lt.tolerance) 
c     .              aircraft%Y_fuse = aircraft%y_cm
c	if((aircraft%Z_fuse+12345.d0).lt.tolerance) 
c     .              aircraft%Z_fuse = aircraft%z_cm

	if (   (aircraft%X_fuseuu.le.-1.d-8)
     .     .or.(aircraft%Y_fusevv.le.-1.d-8)
     .     .or.(aircraft%Z_fuseww.le.-1.d-8)) then 

c	   We move the center of fuselage pressure to the center of mass
	   aircraft%X_fuse = aircraft%x_cm
	   aircraft%Y_fuse = aircraft%y_cm
	   aircraft%Z_fuse = aircraft%z_cm

c	   We build a box for the drag
	   xmax = propeller(1)%x
	   ymax = propeller(1)%y
	   zmax = propeller(1)%z

	   xmin = xmax
	   ymin = ymax
	   zmin = zmax

	   xmax = max(propeller(1)%x-0.07d0*propeller(1)%nx,xmax)  ! default size is 7 cm
	   ymax = max(propeller(1)%y-0.07d0*propeller(1)%ny,ymax)
	   zmax = max(propeller(1)%z-0.07d0*propeller(1)%nz,zmax)

	   xmin = min(propeller(1)%x-0.07d0*propeller(1)%nx,xmin)
	   ymin = min(propeller(1)%y-0.07d0*propeller(1)%ny,ymin)
	   zmin = min(propeller(1)%z-0.07d0*propeller(1)%nz,zmin)
	   
	   do i = 2, aircraft%num_propellers

	      xmax = max(propeller(i)%x,xmax)
	      ymax = max(propeller(i)%y,ymax)
	      zmax = max(propeller(i)%z,zmax)

	      xmin = min(propeller(i)%x,xmin)
	      ymin = min(propeller(i)%y,ymin)
	      zmin = min(propeller(i)%z,zmin)

	      xmax = max(propeller(i)%x-0.07d0*propeller(i)%nx,xmax)
	      ymax = max(propeller(i)%y-0.07d0*propeller(i)%ny,ymax)
	      zmax = max(propeller(i)%z-0.07d0*propeller(i)%nz,zmax)

	      xmin = min(propeller(i)%x-0.07d0*propeller(i)%nx,xmin)
	      ymin = min(propeller(i)%y-0.07d0*propeller(i)%ny,ymin)
	      zmin = min(propeller(i)%z-0.07d0*propeller(i)%nz,zmin)

	   enddo

	   aircraft%X_fuseuu = xmax - xmin
	   aircraft%Y_fusevv = ymax - ymin
	   aircraft%Z_fuseww = zmax - zmin
	
	   write (6,*)
	   write (6,*) ' A drag box was guessed based on motor location.'
	   write (6,*) '   X_fuse  , Y_fuse  , Z_fuse   '
     .                    ,aircraft%X_fuse, aircraft%Y_fuse
     .                    ,aircraft%Z_fuse
	   write (6,*) '   X_fuseuu, Y_fusevv, Z_fuseww '
     .                    ,aircraft%X_fuseuu, aircraft%Y_fusevv
     .                    ,aircraft%Z_fuseww

c	   Divide by a factor of 10 since it seems too big.
	   aircraft%X_fuseuu = aircraft%x_fuseuu / 10.d0
	   aircraft%Y_fusevv = aircraft%y_fusevv / 10.d0
	   aircraft%Z_fuseww = aircraft%z_fuseww / 10.d0

	endif


	endif  ! to there being at least one propeller
	      

c	We now build a small lift and drag table, so that it is possible to get a feel for 
c	where we are in terms of performance.  This table assumes there is a wing, that it is made
c	of two symmetric parts, and that it is numbered 1 and 2 (at this point).
	

	i = 1
	if (aircraft%num_wings.gt.1) then

	   write (6,*) ' Rought lift and drag forces for various ',
     .                 'angles of attack and flight speeds'

           do jj = -10, 10

c	      We do 10 increments of slope, up to stall
	      alpha(jj) = wing(i)%C_Lmax/wing(i)%a*dfloat(jj)/10.d0

	      
c	      The coefficients of lift and drag for each wing element
	     call CLDwing(alpha(jj),wing(i)%a,wing(i)%C_L0,wing(i)%C_Lmin
     .            ,wing(i)%C_Lmax,wing(i)%C_D0,wing(i)%k,wing(i)%C_Dfp
     .            ,C_Lw,C_Dw,iflinear,dC_Ldalpha,dC_Ddalpha)

c	      write (6,*) C_Lw/C_Dw, C_Lw, C_Dw  ! ratio seems to range from 10 to 20

              do k = 1, 10

	         u = dfloat(k)*5.d0  ! in 5 m/s increments
	         qbarprime_w = 0.5d0*rho*u**2

c	         Force equations  (Table 2.5-1 and page 650)
c	         times 2 to give both sides of the wing.
	         force_L(jj,k) = wing(i)%surface_area*qbarprime_w*C_Lw
     .              *2.d0
	         force_D(jj,k) = wing(i)%surface_area*qbarprime_w*C_Dw
     .              *2.d0

	      enddo
	   enddo

	   Write (6,*) ' Lift force (N)'
	   write (6,15) '  Flight speed (m/s)',(5.d0*dfloat(k),k=1,10)
	   write (6,*) ' Angle of Attack alpha'
           write (6,*)  '              (deg) '
	   do jj = -10, 10
	      write (6,20) rad_to_deg*alpha(jj),(force_L(jj,k),k=1,10)
	   enddo
	   Write (6,*) ' Drag force (N)'
	   write (6,15) '  Flight speed (m/s)',(5.d0*dfloat(k),k=1,10)
	   write (6,*) ' Angle of Attack alpha'
           write (6,*)  '              (deg) '
	   do jj = -10, 10
	      write (6,20) rad_to_deg*alpha(jj),(force_D(jj,k),k=1,10)
	   enddo

 20	   format(f20.3,10f9.2)
 15        format(a20,10f9.2)

	endif

	open(70,file='namelist.out',status='unknown')
	write (70,nml=aircraft_data)
	close (70)

	return

	end


	subroutine fly(xdim,x,xd,udim,uc,aircraft,wing,propeller,battery
     .                ,control)
c	This subroutine is a 4th order Runge-Kutta solver for the
c	flight dynamics model.
c	James Walker  26 March 2021
c	3/26/21

	implicit none

	include 'include_parts.h'

	integer xdim, udim

	integer  i, j, l, ipts
	double precision x (xdim), xd(xdim), xo(150), uc(udim)
	integer i_flight_path, i_Hackathon
	double precision Xpath(3), Xpathdot(3), UVWPQR(6)

	parameter(ipts = 100)

	double precision x1(xdim), x2(xdim), x3(xdim), x4(xdim)
	double precision k1(xdim), k2(xdim), k3(xdim), k4(xdim)
	double precision time1   , time2   , time3 

	double precision time, dt, dt_output, time_end, time_out

	double precision timeslip, temp, q0, q1, q2, q3
	double precision cospsi, sinpsi, coshalfpsi,sinhalfpsi
	double precision speed

	double precision l1_error, l_max_error, l1_error_int
	double precision max_error_time, max_error_location(3)
	double precision path_distance , max_error_velocity(3)
	double precision longitudinal_error

	double precision path_score, path_speed, ground_impact_speed
	double precision score, height

	double precision zz(ipts,16)

	double precision KK(50,12)

	double precision phi, q0phi, q1phi, q2phi, q3phi

	logical :: icomplete = .false.

	double precision bat_amps(50), bat_capacity(50)
	double precision motor_amps(50), motor_power(50), ratio
	logical :: lelectrical = .true.
	logical :: lfail = .false.

	integer imode, irad
	logical lstraight
	double precision radius

	time      = aircraft%time
	dt        = aircraft%dt
	dt_output = aircraft%dt_output
	time_end  = aircraft%time_end

	l1_error      = 0.d0
	l_max_error   = 0.d0
	l1_error_int  = 0.d0
	path_distance = 0.d0

	do i = 1, aircraft%num_propellers
	   motor_amps (i) = 0.d0
	   motor_power(i) = 0.d0
	enddo
	do i = 1, aircraft%num_batteries
	   bat_amps    (i) = 0.d0
	   bat_capacity(i) = 0.d0
	enddo

c	The wind speed, if specified
c	We are using array componentwise multiplies in F90
	x(1:13)    = aircraft%x_initial(1:13)    ! non wind speed/air density initial conditions
	uc(1:udim) = aircraft%uc_initial(1:udim) ! initial controls


c	We look to see if there is initial power in the motor,
c	because we want to set the appropriate motor speed if so
c	(rather than have it rise up from zero, and thus enter
c	a transient solution state), since we presume people

	call motor_initialize(time, xdim, x, xd, xo, udim, uc,
     .       aircraft, wing, propeller, battery, control)


 112	continue

	return

	i = 0 ! counter for additional printouts

	max_error_time     = 0.d0
	max_error_location = 0.d0
	max_error_velocity = 0.d0

	do while ((time.lt.time_end).and.(.not.icomplete))

	   call flightpath(i_Hackathon, i_flight_path, time, timeslip
     .          , speed, Xpath, Xpathdot, lstraight, radius, icomplete)
	   call autopilot(time, xdim, x, xd, xo, udim, uc,
     .                aircraft, wing, propeller, battery, control
     .              , Xpath, Xpathdot, lstraight, radius, i_flight_path)
	   call fderiv(time, xdim, x, xd, xo, udim, uc,
     .                aircraft, wing, propeller, battery, control)

c	write (6,56) ' 1 x =  ',x(1:xdim)
c	write (6,56) ' 1 xd = ',xd(1:xdim)
c	write (6,56) ' 1 uc = ',uc(1:udim)
c	write (6,56) ' DCM = ',xo(11),xo(12),xo(13)
c	write (6,56) ' DCM = ',xo(14),xo(15),xo(16)
c	write (6,56) ' DCM = ',xo(17),xo(18),xo(19)
56	format(a10,20e12.4)

	

c	   Vector adds here, for ki, xi
	   k1    = dt*xd
	   x1    = x    + 0.5d0*k1
	   time1 = time + 0.5d0*dt

	   call flightpath(i_Hackathon, i_flight_path, time1, timeslip
     .          , speed, Xpath, Xpathdot, lstraight, radius, icomplete)
	   call autopilot(time1, xdim, x1, xd, xo, udim, uc,
     .                aircraft, wing, propeller, battery, control
     .              , Xpath, Xpathdot, lstraight, radius, i_flight_path)
	   call fderiv(time1, xdim, x1, xd, xo, udim, uc,
     .                 aircraft, wing, propeller, battery, control)
c	write (6,56) ' 2 x =  ',x(1:xdim)
c	write (6,56) ' 2 xd = ',xd(1:xdim)
c	write (6,56) ' 2 uc = ',uc(1:udim)
c	write (6,56) ' DCM = ',xo(11),xo(12),xo(13)
c	write (6,56) ' DCM = ',xo(14),xo(15),xo(16)
c	write (6,56) ' DCM = ',xo(17),xo(18),xo(19)

	   k2    = dt*xd
	   x2    = x    + 0.5d0*k2
	   time2 = time + 0.5d0*dt

	   call flightpath(i_Hackathon, i_flight_path, time2, timeslip
     .          , speed, Xpath, Xpathdot, lstraight, radius, icomplete)
	   call autopilot(time2, xdim, x2, xd, xo, udim, uc,
     .                aircraft, wing, propeller, battery, control
     .              , Xpath, Xpathdot, lstraight, radius, i_flight_path)
	   call fderiv(time2, xdim, x2, xd, xo, udim, uc,
     .                 aircraft, wing, propeller, battery, control)
c	write (6,56) ' 3 x =  ',x(1:xdim)
c	write (6,56) ' 3 xd = ',xd(1:xdim)
c	write (6,56) ' 3 uc = ',uc(1:udim)
c	write (6,56) ' DCM = ',xo(11),xo(12),xo(13)
c	write (6,56) ' DCM = ',xo(14),xo(15),xo(16)
c	write (6,56) ' DCM = ',xo(17),xo(18),xo(19)

	   k3    = dt*xd
	   x3    = x    + k3
	   time3 = time + dt

	   call flightpath(i_Hackathon, i_flight_path, time3, timeslip
     .          , speed, Xpath, Xpathdot, lstraight, radius, icomplete)
	   call autopilot(time3, xdim, x3, xd, xo, udim, uc,
     .                aircraft, wing, propeller, battery, control
     .              , Xpath, Xpathdot, lstraight, radius, i_flight_path)
	   call fderiv(time3, xdim, x3, xd, xo, udim, uc,
     .                 aircraft, wing, propeller, battery, control)
c	write (6,56) ' 4 x =  ',x(1:xdim)
c	write (6,56) ' 4 xd = ',xd(1:xdim)
c	write (6,56) ' 4 uc = ',uc(1:udim)
c	write (6,56) ' DCM = ',xo(11),xo(12),xo(13)
c	write (6,56) ' DCM = ',xo(14),xo(15),xo(16)
c	write (6,56) ' DCM = ',xo(17),xo(18),xo(19)

	   k4 = dt*xd
	
	   x    = x + k1/6.d0 + k2/3.d0 + k3/3.d0 + k4/6.d0
	   time = time + dt

c	   After the official step is taken, this is where we should
c	   check to see what happened to the battery - how many amps
c	   are flowing, how much was lost.  Where should that stuff be
c	   stored?  For the motors we have spin rates and amps - 
c	   are those state variables?  Do we have more state variables
c	   that control inputs?  We want only the minimum number of
c	   control inputs because we are going to be doing various
c	   optimizations based on those control inputs.  We don't really
c	   care how many outputs there are, but it does seem there sould
c	   only be one per u.  But how do I understand what they are?

cc	   write (6,*) ' time, x, y, z ',time,x(11),x(12),x(13)
cc	   write (6,*)' Unorth, Veast, Wdown ',time,xd(11),xd(12),xd(13)
cc	   write (6,*)' Xnorth, Yeast, Zdown ',time,x (11),x (12),x (13)
cc	   write (6,*)' Wind speed ',xo(1),xo(2),xo(3)

c	      write (6,*) ' P Q R dot ',xd(4),xd(5),xd(6)

c	   We compute the error
	   path_speed = sqrt(x(1)**2 + x(2)**2 + x(3)**2)
	   path_distance = path_distance + path_speed * dt
	   call path_error(time, dt, path_speed, xdim, x, Xpath, Xpathdot
     .        , l1_error, l_max_error, l1_error_int, longitudinal_error
     .        , path_distance
     .        , max_error_time, max_error_location, max_error_velocity)

c	   Okay - here's a tricky attemp to adjust the path.  We will see
c	   if fiddling with the longitudinal expectation can help with
c	   bringing us back to the path.  Are we seeing a phugoid?
c	   timeslip = - longitudinal_error/speed


c	   We check on the current and batteries. 
	   do j = 1, aircraft%num_batteries
c	      Battery capacity units are milli amp hours
c	      the capacity of each battery is integrated in the state vector
c	      Capacity    ! (mAh = milli amp hours) Total battery charge (1 milli amp hour = 3.6 coulombs)
c	      C_Continuous         ! In regular use, max amp is (C / 1000) * capacity in mAh
c	      We are interested in time to 20% power
	      bat_capacity(j) = x(13 + aircraft%num_propellers + j)
     .                         / (battery(j)%Capacity
     .                          *3600.d0/1000.d0)
	      ratio  = xd(13 + aircraft%num_propellers + j)
     .                    /(battery(j)%C_Continuous
     .                    * battery(j)%Capacity/1000.d0)
	      bat_amps(j) = max(bat_amps(j),ratio)
	      if (.not.aircraft%ignore_electrical) then
	         if (bat_capacity(j).gt.0.8d0) goto 200
                 if (bat_amps    (j).gt.1.d0 ) goto 200
	      endif
	   enddo
	   do j = 1, aircraft%num_propellers
c	      Now we check for the current in the motors
	      ratio = xo(29 + aircraft%num_propellers + j)
     .              * propeller(j)%efficiency_ESC          ! amps through that motor
     .               / propeller(j)%I_max
	      motor_amps(j) = max(ratio,motor_amps(j))
	      ratio = xo(29 + j)
     .               * propeller(j)%efficiency_ESC ! this is the individual motor power
     .               / propeller(j)%maxpower
	      motor_power(j) = max(ratio,motor_power(j))
	      if (.not.aircraft%ignore_electrical) then
	         if (motor_amps (j).gt.1.d0) goto 200
	         if (motor_power(j).gt.1.d0) goto 200
	      endif
	   enddo


	   if (time.gt.time_out) then

	      write (6,10) time,rad_to_deg*xo(5),rad_to_deg*xo(6)
     .         ,rad_to_deg*xo(7)
     .         ,xd(11),xd(12),xd(13),x(11)
     .         ,x(12),x(13)
     .         ,xo(8),rad_to_deg*xo(9),rad_to_deg*xo(10)
     .         ,sqrt(xo(26)**2+xo(27)**2+xo(28)**2)  ! thrust
     .         ,sqrt(xo(20)**2+xo(21)**2+xo(22)**2)  ! lift
     .         ,sqrt(xo(23)**2+xo(24)**2+xo(25)**2)  ! drag
 10	      format (16f10.3)

	      time_out = time_out+dt_output

c	      Save up things for another printout
	      i = i + 1

	      i = min(i,ipts)  ! to protect against going beyond the end.
	
	      zz(i,1) = time
	      zz(i,2) = xo(26)  ! thrust
	      zz(i,3) = xo(27)
	      zz(i,4) = xo(28)
	      zz(i,5) = xo(20)  ! lift
	      zz(i,6) = xo(21)  
	      zz(i,7) = xo(22)
              zz(i,8) = xo(23)  ! drag
              zz(i,9) = xo(24)
	      zz(i,10) = xo(25)
	      zz(i,11) = time  ! temporary holder

	      zz(i,2) = xo(14)
	      zz(i,3) = xo(15)
	      zz(i,4) = xo(16)
	      zz(i,5) = xo(17)

	      zz(i,2) = x(14)
	      zz(i,3) = x(15)
	      zz(i,4) = x(16)
	      zz(i,5) = x(17)

	      zz(i,6) = xd(14)
	      zz(i,7) = xd(15)
	      zz(i,8) = xd(16)
	      zz(i,9) = xd(17)

	      write (55,55) time,x(11), xpath(1), x(12), xpath(2)
     .            , x(13), xpath(3)
     .            , x(11)-xpath(1), x(12)-xpath(2), x(13)-xpath(3)
     .            , longitudinal_error, l1_error, l_max_error
     .            , l1_error_int
 55	      format(14e12.4)
c	     Phi, theta, psi
	      write (56,55) time, x(11),x(12),x(13),l1_error,
     .              x(7), x(8), x(9), x(10), rad_to_deg*xo(5),
     .              rad_to_deg*xo(6), rad_to_deg*xo(7)
     .             ,sqrt(x(7)**2 + x(8)**2 + x(9)**2 + x(10)**2) - 1.d0

	   endif

	enddo  ! on time or completion

	lelectrical = .false.

	write (6,*)
	if (icomplete) then
	   write (6,*)  ' Calculation completed at time ',sngl(time),
     .                  ' due to course completion.'
	else
	   write (6,*)  ' Calculation completed at time ',sngl(time),
     .                  ' due to set time_end.'
	endif
 
 200	continue  ! if we come here, it was because of electrical problems

	write (6,*)

	if (lelectrical) then
	   write (6,*)  ' Calculation completed at time ',sngl(time),
     .                  ' due to an electrical issue.'
	endif
	write (6,*)
	write (6,*) ' The following information is about electrical'
     .             ,' performance.'
	if (aircraft%ignore_electrical) write (6,*)' Flight termination'
     .    ,' due to electrical issues was turned off.'
	lfail = .false.
	do j = 1, aircraft%num_batteries
	   write (6,220) j,bat_capacity(j),bat_amps(j)
 220       format('   Battery ',i2,'  fraction capacity used '
     .                  ,f10.4, ' and fraction'
     .                ,' max continuous amperage used ',f10.4)
	   if (bat_capacity(j).gt.0.8d0) lfail = .true.
	   if (bat_amps(j)    .gt.1.0d0) lfail = .true.
	enddo
	do j = 1, aircraft%num_propellers
	   write (6,221) j,motor_amps(j),motor_power(j)
 221       format('   Motor ',i2,' fraction max amps reached '
     .                ,f10.4,' and fraction max '
     .                 ,'power reached ',f10.4)
	   if (motor_amps(j) .gt.1.d0) lfail = .true.
	   if (motor_power(j).gt.1.d0) lfail = .true.
	enddo
	

c	Storing, though not sure who may need this information
	control%path_distance      = path_distance
	control%longitudinal_error = longitudinal_error
	control%l1_error           = l1_error
	control%l_max_error        = l_max_error
	control%l1_error_int       = l1_error_int
	control%max_error_time     = max_error_time
	control%max_error_location = max_error_location
	control%max_error_velocity = max_error_velocity

	path_speed = path_distance/time

c	write (6,*)
c	
c	do j = 1, i
c	   write (6,10) (zz(j,l),l=1,16)
c	enddo

	write (6,*)
	write (6,*)' Wind speed (m/s) ',sngl(xo(1)),sngl(xo(2))
     .                                 ,sngl(xo(3))
	write (6,*) ' Air density (kg/m^3) ',sngl(xo(4))

c	write (6,*) ' DCM '
c	write (6,*) xo(11),xo(12),xo(13)
c	write (6,*) xo(14),xo(15),xo(16)
c	write (6,*) xo(17),xo(18),xo(19)

c	This is the scoring section

	score = 0.d0

	write (6,*)
	write (6,*) ' Hackathon ',control%i_Hackathon
 	write (6,*) ' Path performance, flight path '
     .                ,control%i_flight_path
	write (6,*)

	write (6,*) ' Maximum distance error (measured perpendicularly)'
     .            ,' from flight path was ', sngl(l_max_error), ' m'
	if (i_flight_path.eq.1) then
c	   Straight line flight
           if (path_distance.gt.2000.d0) then
	      score = path_distance/10.d0
	      Write (6,*) ' Score is distance flown / 10 = ',sngl(score)
     .             ,' since minimum flight distance of 2 km achieved'
	      score = min(score, 400.d0)
	      Write (6,*) ' The maximum score on flight distance is 400'
     .              ,' so score is = ',sngl(score)
	      score = score - 10.d0*l_max_error
	      score = max(score, 0.d0)
	      write (6,*) ' Score is now deducted for inaccuracy'
     .            ,' by 10 * distance error ',sngl(- 10.d0*l_max_error)
	   else
	      write (6,*) ' Score is 0 since minimum flight distance '
     .             ,'of 2 km not achieved; distance was '
     .             ,sngl(path_distance),' m'
	   endif
	else if (i_flight_path.eq.2) then
c	   Straight line flight
           if (path_distance.gt.2000.d0) then
	      score = path_distance/10.d0
	      Write (6,*) ' Score is distance flown / 10 = ',sngl(score)
     .             ,' since minimum flight distance of 2 km achieved'
	      score = min(score, 400.d0)
	      Write (6,*) ' The maximum score on flight distance is 400'
     .              ,' so score is = ',sngl(score)
	      score = score - 10.d0*l_max_error
	      score = max(score, 0.d0)
	      write (6,*) ' Score is now deducted for inaccuracy'
     .            ,' by 10 * distance error ',sngl(- 10.d0*l_max_error)
	   else
	      write (6,*) ' Score is 0 since minimum flight distance '
     .             ,'of 2 km not achieved; distance was '
     .             ,sngl(path_distance),' m'
	   endif
	else if (i_flight_path.eq.3) then
c	   Circle
	   if (icomplete) then
              score = 300.d0 - 50.d0*l_max_error
	      write (6,*) ' The full circle was flown, so 300'
     .                     ,' points are initially awarded '
	      write (6,*) ' Score is now deducted for inaccuracy'
     .         ,' by 50 * '
     .         ,'maximum lateral error ', sngl(- 50.d0*l_max_error)
	      score = max(score, 0.d0)
	   else
	      write (6,*) ' Score is zero since full circle not flown'
	   endif  
	else if (i_flight_path.eq.4) then
c	   Vertical rise
	   if (i_Hackathon.eq.0) height = 200.d0
	   if (i_Hackathon.eq.1) height = 150.d0
	   if (time.gt.200.d0) then
              if (abs(x(13)+height).le.2.d0) then
                 score = min(time,400.d0)
	         write (6,*) ' Since minimum flight time of 200 s '
     .                ,'was achieved and final vehicle height was '
     .                ,'within 2 m of specified hover height, '
     .                ,'score is flight time with max 400 pts.'
	         write (6,*) ' Final height (-Z) = '
     .                      ,sngl(-x(13)),'m;  flight time = '
     .                      ,sngl(time),'s'
	         score = max(score, 0.d0)
	      else
	         write (6,*) ' Score is zero since final hover not'
     .             ,' within 2 m of specified height; height -Z = '
     .             ,sngl(-x(13)),'m'
	      endif
	   else
	      write (6,*) ' Score is zero since flight time less'
     .             ,' than 200 s; flight time = ',time,'s'
	   endif
	else if (i_flight_path.eq.5) then
c	   The racing oval
	   if (icomplete) then
              score = 200.d0 + max(0.d0, 350.d0 - time)
	      write (6,*) ' Complete oval being flown gives 200 pts.'
	      write (6,*) ' Every second less than 350 seconds gives'
     .            ,' one additional point; flight time = ',sngl(time)
	   else
	      write (6,*) ' Zero points since complete oval not flown'
	   endif
	   if (l_max_error.ge.10.d0) then
 	      score = 0.d0
	      write (6,*) ' Since distance error from specified flight'
     .                   ,' path exceeds 10 m, 0 points are awarded'
	   else
	      write (6,*) ' Stayed with 10 m of flight path (otherwise'
     .                   ,' no points would have been awarded)'
	   endif
	endif

	path_score = max(score,0.d0)

	path_score = dfloat(int(path_score + 0.5d0))

	if (aircraft%ignore_electrical.and.lfail) then
	   write (6,*) ' Score ignoring electrical issues '
     .                  ,sngl(path_score)
           path_score = 0.d0
	   write (6,*) ' However, score set to 0 because of electrical'
     .                ,' issue.'
	endif
	         
	write (6,*) ' Final score (rounded) = ',sngl(path_score)


c	The scoring write-out for the path 
	do i = 1, 2
	   if (i.eq.1) j = 6
	   if (i.eq.2) j = aircraft%i_out_metrics
	   write (j,*)
	   write (j,*) ' Measures of flight path performance '
     .      ,'(distance in meters, time is seconds, speed in '
     .      ,'meters per second)'
	   write (j,*)
	   if (path_score.gt.0.d0) then
	      write (j,*) ' Flight path was successfully traversed.'
	   else
	      write (j,*) ' Flight path was not successfully traversed.'
	   endif
	   write (j,*)
	   write (j,*) '#Metrics '
	   write (j,*) 'Flight_distance ',sngl(path_distance)
	   write (j,*) 'Time_to_traverse_path ', sngl(time)
	   write (j,*) 'Average_speed_to_traverse_path '
     .                  ,sngl(path_speed)
	   write (j,*) 'Maximimum_error_distance_during_flight '
     .                    ,sngl(l_max_error)
	   write (j,*) 'Time_of_maximum_distance_error '
     .                    ,sngl(max_error_time)
	   write (j,*) 'Location_of_maximum_distance_error '
     .                    ,sngl(max_error_location)
           write (j,*) 'Velocity_at_time_of_maximum_distance_error '
     .                    ,sngl(max_error_velocity)
           write (j,*) 'Spatial_average_distance_error '
     .                    ,sngl(l1_error_int)
           write (j,*) 'Maximum_ground_impact_speed '
     .                    ,sngl(ground_impact_speed)
           write (j,*) 'Path_traverse_score_based_on_requirements '
     .                    ,sngl(path_score)
	   write (j,*) 'Input_LQR_weights_Qp_Qv_Qav_Qang_R ', 
     .           sngl(control%Q_position),sngl(control%Q_velocity)
     .          ,sngl(control%Q_angular_velocity),sngl(control%Q_angles)
     .          ,sngl(control%R)


 99	   format(a47,5f12.6)

	enddo

	write (aircraft%i_out_score,*) '#Metrics '
	write (aircraft%i_out_score,*) 'Flight_path ',i_flight_path
	write (aircraft%i_out_score,*) 
     .          'Path_traverse_score_based_on_requirements '
     .                    ,sngl(path_score)

	return
	end




c	subroutine autopilot(time, xdim, x, xd, xo, udim, uc,
c     .                aircraft, wing, propeller, battery, control
c     .              , Xpath, Xpathdot, lstraight, radius, ipath)
cc	This is a dummy routine, to be replaced by the autopilot.
cc	Can be a user created subroutine.  Various variables are
cc	described in subroutine fderiv and include_parts.h
cc	The intent is that everything on this list is input,
cc	and the only output changes are the control array uc
cc	and (perhaps) the internal information in the type control.
c
c	implicit none
c
c	include 'include_parts.h'
c
c	integer xdim,udim
c
c	double precision time, x(xdim), xd(xdim), xo(150), uc(udim)
c
cc	Coding to update uc, etc., follows.
c
c	return
c	end




	subroutine fderiv(time, xdim, x, xd, xo, udim, uc,
     .       aircraft, wing, propeller, battery, control)

c	James D. Walker, Southwest Research Institute

c	General
c	x are aircraft state variables (input)
c	xd are aircraft state varaible derivatives (output)
c	xo are various state variables (output)
c	uc are aircraft control variables (input)
c	Others are aircraft informaation

c	Coding for the physics and geometry of the flight dynamics model.
c	Only x(1 to 13) are part of the differential update, the other
c	values are being used for state storage.

c	State and control variables
c	I'm restructuring these, for quaternions and for other basic variables.
c	
c       Body fixed frame is assumed to be in frd (x=forward, y=right, z=down)
c	That means that altitude is -z, not z, and moving upward means W < 0.
c
c	1 to 13 are state variables.  The x versions are input, the derivatives (xd) are output.
c
c	U     = x(1)         Udot     = xd(1)    body fixed frame speed and acceleration
c	V     = x(2)         Vdot     = xd(2)    body fixed frame speed and acceleration
c	W     = x(3)         Wdot     = wd(3)    body fixed frame speed and acceleration
c
c	P     = x(4)         Pdot     = xd(4)    body fixed frame rotation and acceleration (roll rate)
c	Q     = x(5)         Qdot     = xd(5)    body fixed frame rotation and acceleration (pitch rate)
c	R     = x(6)         Rdot     = xd(6)    body fixed frame rotation and acceleration (yaw rate)
c
c       q0    = x(7)         q0dot    = xd(7)    real quaternion (connecting body fixed frame to world frame)
c	q1    = x(8)         q1dot    = xd(8)    i quaternion
c	q2    = x(9)         q2dot    = xd(9)    j quaternion
c	q3    = x(10)        q3dot    = xd(10)   k quaternion
c
c	Xnorth = x(11)       Unorth   = xd(11)   North displacement and speed (world frame)
c	Yeast  = x(12)       Veast    = xd(12)   East displacement and speed (world frame)
c       Zdown  = x(13)       Wdown    = xd(13)   Down displacement and speed (world frame); down is positive
c
c	motor_omega(i) = x(13+i)   motor_omegadot= x(13+i)  These are the motor rotational motor speeds
c                                                and accelerations.  The i is the propeller index, since every
c						 propeller can spin.  In part we are storing it here so that
c						 the ODE solver that calls this routine will compute our motor
c						 response, since there is a delay from command to execution
c						 (as given by propeller(i)%tau_EXC.
c       deflections (these are states?)
c                                                Motors accelerate.  Right now we assume that ailerons, flaps,
c                                                etc. have instanteous response, so we are not solving an ODE
c						 (i.e., xd = 0 for them).
c
c	charge on battery x(13+num_propellers+num_servos+i) current xd( )  Current in amps for each batter.  We have it
c						 as a state so that the total charge from the battery is 
c						 integrated (a minus sign, so that the batter charge decreases).
c
c	Outputs in xo
c
c	Unwind = xo(1)	                   	 North wind speed in world frame (input)
c	vewind = xo(2)	                  	 East wind speed in world frame (input)
c	Wdwind = xo(3)	                  	 Down wind speed in world frame (input)
c
c	rho    = xo(4)                    	 Density of air used by this subroutine all (output)
c
c	phi    = xo(5)				 Euler about x (roll)  (output)
c	theta  = xo(6)                           Euler about y (pitch)  (output)
c	psi    = xo(7)                           Euler about z (yaw)  applied order is z then y then x (output)
c
c	Vt     = xo(8)				 Air speed (output)
c	alpha  = xo(9)				 Angle of attacka (output)
c	beta   = xo(10)				 Sideslip angle (output)
c
c       
c	In particular, the direction cosine matrix (DCM) is computed, used internally, and available as output.
c
c       dcm_xx = xo(11)                          Transformation matrix from world frame to body fixed frame
c       dcm_xy = xo(12)                          Transformation matrix from world frame to body fixed frame
c       dcm_xz = xo(13)                          Transformation matrix from world frame to body fixed frame
c       dcm_yx = xo(14)                          Transformation matrix from world frame to body fixed frame
c       dcm_yy = xo(15)                          Transformation matrix from world frame to body fixed frame
c       dcm_yz = xo(16)                          Transformation matrix from world frame to body fixed frame
c       dcm_zx = xo(17)                          Transformation matrix from world frame to body fixed frame
c       dcm_zy = xo(18)                          Transformation matrix from world frame to body fixed frame
c       dcm_zz = xo(19)                          Transformation matrix from world frame to body fixed frame
c	                                         We use this to get gravity vector in rotating frame.
c	                                         To go from body to world, take the transpose (to get world speed)
c
c	Some additional information, for better understanding of what is occuring
c
c	Fx_lift = xo(20)                         Lift forces in body fixed frame
c	Fy_lift = xo(21)
c	Fz_lift = xo(22)
c
c	Fx_drag = xo(23)		         Drag forces in body fixed frame
c	Fy_drag = xo(24)
c	Fz_drag = xo(25)
c
c       Thrust_x = xo(26)			 Total thrust
c	Thrust_y = x0(27)
c	Thrust_z = xo(28)
c
c       Power_tot = xo(29)                       Estimated total power from all propellers
c	Power(i)  = xo(29+i)                       Estimated power from each individual motor.  Primarily 
c	current at each motor (i) = xo(29+numpropellers+i)  Estimated current draw for each motor
c	Current at each battery (i) = xo(29+2*numpropellers+i)  Eastimated current for each battery
c
c	Control variables (input)
c
c	uc(icontrol)                             These are the control variables.  They are assumed to range from 0 to 1.
c						 All ranges limits are hence checked "by definition."
c
c	xdot = Ax + Bu + G			 This is the linearized control entries, for use by the autopilot.
c						 A is xdim x xdim, B is xdim x udim, and G is a xdim vector.
c						 Letters lin are used in the actual name (Alin, Blin, Glin) for linear.

	implicit none

	include 'include_parts.h'

	integer xdim,udim

	double precision time, x(xdim), xd(xdim), xo(150), uc(udim)
c	double precision Alin(xdim,xdim), Blin(xdim,udim), Glin(xdim)
c	This is a temporary sizing - the 50 will be replaced by udim when
c	the appropriate calls are all put in place.  6/9/21
	double precision Alin(12,12), Blin(12,50), Glin(12)

	double precision motor_omega(50), motor_omegadot(50)


	double precision Vt, U, V, W
	double precision alpha, beta, phi, theta, psi, P, Q, R
	double precision Fx, Fx_tot, Fy, Fy_tot, Fz, Fz_tot
	double precision Mx_tot, My_tot, Mz_tot, Mx, My, Mz
	double precision    Ixx,   Ixy,   Ixz,   Iyy,   Iyz,   Izz
	double precision InvIxx,InvIxy,InvIxz,InvIyy,InvIyz,InvIzz
	double precision Vi, Viw, VW_w
	double precision x_w, y_w, z_w, U_w, v_w, w_w
	double precision delta_1, delta_2, delta_a
	double precision alpha_w, qbarprime_w, C_Lw, C_Dw, alpha_w_d
	logical          iflinear
	double precision dC_Ldalpha,dC_Ddalpha
	double precision Force_Lift_w, Force_Drag_w
	integer          i, k, l, ierror, kpast, lpast
	double precision J, n, D, Vpath, Thrust, Power,  Torque, Cp, Ct
	double precision temp, tau, current, voltage
	double precision Udot, Vdot, Wdot
	double precision ox, oy, oz, denom, grav
	double precision rho
	double precision Unorth, Veast, Wdown
	double precision Xnorth, Yeast, Zdown
	double precision Pdot, Qdot, Rdot
	double precision q0   , q1   , q2   , q3
	double precision q0dot, q1dot, q2dot, q3dot
	double precision Unwind, Vewind, Wdwind
	double precision Uwind, Vwind, Wwind, Up, Vp, Wp
	integer          icontrol, icontrol1, icontrol2
	logical          ic1, ic2
	double precision Fx_lift, Fy_lift, Fz_lift
	double precision Fx_drag, Fy_drag, Fz_drag
	double precision dxkm1, dxk, dylm1, dyl
	double precision Thrust_x, Thrust_y, Thrust_z
	double precision a, b, c, Kt, Kv, Rw, omega, amps, Power_tot
	integer          jdw
	double precision omega0, omega1, ff, ffp, a1, b1, omega2, DeltaJ
	double precision J1, J2, O1, O2, Cp11, Cp21, Cp22, Cp12
	double precision aa, bb, cc, dd, ee,   a0,   a2,   a3
	logical       :: ip  = .false.
	double precision I0, motor_efficiency

	double precision Fx_tot_U, Fy_tot_U, Fz_tot_U
	double precision Fx_tot_V, Fy_tot_V, Fz_tot_V
	double precision FX_tot_W, Fy_tot_W, Fz_tot_W

	double precision Fx_tot_P, Fy_tot_P, Fz_tot_P
	double precision Fx_tot_Q, Fy_tot_Q, Fz_tot_Q
	double precision FX_tot_R, Fy_tot_R, Fz_tot_R

	double precision Mx_tot_U, My_tot_U, Mz_tot_U
	double precision Mx_tot_V, My_tot_V, Mz_tot_V
	double precision Mx_tot_W, My_tot_W, Mz_tot_W

	double precision Mx_tot_P, My_tot_P, Mz_tot_P
	double precision Mx_tot_Q, My_tot_Q, Mz_tot_Q
	double precision Mx_tot_R, My_tot_R, Mz_tot_R

	double precision Fx_U, Fy_U, Fz_U
	double precision Fx_V, Fy_V, Fz_V
	double precision FX_W, Fy_W, Fz_W

	double precision Fx_P, Fy_P, Fz_P
	double precision Fx_Q, Fy_Q, Fz_Q
	double precision FX_R, Fy_R, Fz_R

	double precision Mx_U, My_U, Mz_U
	double precision Mx_V, My_V, Mz_V
	double precision Mx_W, My_W, Mz_W

	double precision Mx_P, My_P, Mz_P
	double precision Mx_Q, My_Q, Mz_Q
	double precision Mx_R, My_R, Mz_R

	double precision ox_P, oy_P, oz_P
	double precision ox_Q, oy_Q, oz_Q
	double precision ox_R, oy_R, oz_R

	double precision xo_d(11:19,0:3)

	double precision dq0dq1, dq0dq2, dq0dq3
	double precision dq1dq0, dq1dq2, dq1dq3
	double precision dq2dq0, dq2dq1, dq2dq3
	double precision dq3dq0, dq3dq1, dq3dq2
	
	double precision Fx_tot_uc(50), Fy_tot_uc(50), Fz_tot_uc(50)
	double precision Mx_tot_uc(50), My_tot_uc(50), Mz_tot_uc(50)

	double precision Fx_uc1, Fy_uc1, Fz_uc1
	double precision Fx_uc2, Fy_uc2, Fz_uc2

	double precision Mx_uc1, My_uc1, Mz_uc1
	double precision Mx_uc2, My_uc2, Mz_uc2

	double precision dqbardU , dqbardV , dqbardW
	double precision dalphadU, dalphadV, dalphadW
	double precision dalphaduc1, dalphaduc2

	double precision Force_L_U, Force_L_V, Force_L_W
	double precision Force_D_U, Force_D_V, Force_D_W

	double precision Force_L_P, Force_L_Q, Force_L_R
	double precision Force_D_P, Force_D_Q, Force_D_R

	double precision Force_L_uc1, Force_L_uc2
	double precision Force_D_uc1, Force_D_uc2

	double precision dCpdJ, dCpdomega
	double precision dCtdJ, dCtdomega

	double precision dThrustdVoltage, domegadVoltage, dJdomega
	double precision dThrustduc1, dThrustdU, dThrustdV, dThrustdW

	double precision dUdP, dUdQ, dUdR, dVdP, dVdQ, dVdR
	double precision dWdP, dWdQ, dWdR

	integer qskip
	logical compute_A

	if (aircraft%debug.ne.0) ip = .true.	

	if (ip) write (6,*) ' Fderiv Time = ',time
	kpast = 0
	lpast = 0

	compute_A = control%compute_A
	if (compute_A) then
	   Alin = 0.d0
	   Blin = 0.d0
	   Glin = 0.d0
	endif


c	Read in the current states

	U      = x(1)
	V      = x(2)
	W      = x(3)

	P      = x(4)
	Q      = x(5)
	R      = x(6)

	q0     = x(7)
	q1     = x(8)
	q2     = x(9)
	q3     = x(10)

	Xnorth = x(11)
	Yeast  = x(12)
	Zdown  = x(13)

	do i= 1, aircraft%num_propellers
           motor_omega(i) = x(13 + i)
	enddo

c	Control - which quaternion to leave out of Alin, Blin, Glin
	qskip = control%qskip

c	Set atmospheric conditions.  Right now these are constants.
c	In the future they could be functions of altitude (rho) and
c	time (wind speed, gusts, etc.)
	xo(1) = aircraft%Unwind
	xo(2) = aircraft%Vewind
	xo(3) = aircraft%Wdwind
	xo(4) = aircraft%rho

	Unwind = xo(1)  ! U north world speed of wind
	Vewind = xo(2)  ! V east world speed of wind
	Wdwind = xo(3)  ! W down world speed of wind

	rho    = xo(4)  ! density

c	Some pure geometry (no physics); Direction cosine matrix (dcm) from the quaternions.
	xo(11) = q0**2 + q1**2 - q2**2 - q3**2   ! xx
	xo(12) = 2.d0*(q1*q2 + q0*q3)            ! xy
	xo(13) = 2.d0*(q1*q3 - q0*q2)            ! xz
	xo(14) = 2.d0*(q1*q2 - q0*q3)            ! yx
	xo(15) = q0**2 - q1**2 + q2**2 - q3**2   ! yy
	xo(16) = 2.d0*(q2*q3 + q0*q1)            ! yz
	xo(17) = 2.d0*(q1*q3 + q0*q2)            ! zx
	xo(18) = 2.d0*(q2*q3 - q0*q1)            ! zy
	xo(19) = q0**2 - q1**2 - q2**2 + q3**2   ! zz

c	For those interested, we compute the Euler angles
	xo(5) = atan2(xo(16),xo(19))                ! phi about x (roll)
	xo(6) = asin(-min(1.d0,max(-1.d0,xo(13))))  ! theta, about y (pitch)
	xo(7) = atan2(xo(12),xo(11))                ! psi, about z (yaw)


c	Now we get the speed of the wind in the body fixed frame
	Uwind = xo(11)*Unwind + xo(12)*Vewind + xo(13)*Wdwind
	Vwind = xo(14)*Unwind + xo(15)*Vewind + xo(16)*Wdwind
	Wwind = xo(17)*Unwind + xo(18)*Vewind + xo(19)*Wdwind

	if (ip) write (6,*) ' uwwind ',Uwind,Vwind,Wwind

c	The air speed is given by U prime - sometimes we need it, sometimes
c	we need the ground speed which is just regular U.
	Up = U - Uwind
	Vp = V - Vwind
	Wp = W - Wwind

	if (ip) write (6,*) ' upwind ',Up,Vp,Wp

c	These are the speed (Vt), angle of attack (alpha), and the sideslip (beta)
c	With the new coding, these are not being used, but they are still here for reference.
	Vt     = sqrt(Up**2+Vp**2+Wp**2)  ! total speed Vt
	alpha  = atan2(Wp,Up)             ! angle of attack alpha
	beta   = asin(Vp/Vt)              ! sideslip angle beta
	xo(8)  = Vt
	xo(9)  = alpha
	xo(10) = beta


c	This is the loop over the wing elements, which includes half main wings,
c	vertical and horizontal tails.  Each has its own flow based on where it
c	is located and the rotational motion of the vehicle.
c	Computing the induced flow vi is in the future.

	Fx_tot   = 0.d0
	Fy_tot   = 0.d0
	Fz_tot   = 0.d0

	Mx_tot   = 0.d0
	My_tot   = 0.d0
	Mz_tot   = 0.d0

	Thrust_x = 0.d0
	Thrust_y = 0.d0
	Thrust_z = 0.d0
	
	Fx_lift  = 0.d0
	Fy_lift  = 0.d0
	Fz_lift  = 0.d0

	Fx_drag  = 0.d0
	Fy_drag  = 0.d0
	Fz_drag  = 0.d0

	Power_tot = 0.d0

	if (compute_A) then
	   Fx_tot_U = 0.d0
	   Fx_tot_V = 0.d0
	   Fx_tot_W = 0.d0

	   Fy_tot_U = 0.d0
	   Fy_tot_V = 0.d0
	   Fy_tot_W = 0.d0

	   Fz_tot_U = 0.d0
	   Fz_tot_V = 0.d0
	   Fz_tot_W = 0.d0

	   Mx_tot_U = 0.d0
	   Mx_tot_V = 0.d0
	   Mx_tot_W = 0.d0

	   My_tot_U = 0.d0
	   My_tot_V = 0.d0
	   My_tot_W = 0.d0

	   Mz_tot_U = 0.d0
	   Mz_tot_V = 0.d0
	   Mz_tot_W = 0.d0

	   Fx_tot_P = 0.d0
	   Fx_tot_Q = 0.d0
	   Fx_tot_R = 0.d0

	   Fy_tot_P = 0.d0
	   Fy_tot_Q = 0.d0
	   Fy_tot_R = 0.d0

	   Fz_tot_P = 0.d0
	   Fz_tot_Q = 0.d0
	   Fz_tot_R = 0.d0

	   Mx_tot_P = 0.d0
	   Mx_tot_Q = 0.d0
	   Mx_tot_R = 0.d0

	   My_tot_P = 0.d0
	   My_tot_Q = 0.d0
	   My_tot_R = 0.d0

	   Mz_tot_P = 0.d0
	   Mz_tot_Q = 0.d0
	   Mz_tot_R = 0.d0

	   do i = 1, udim

	      Fx_tot_uc(i) = 0.d0
	      Fy_tot_uc(i) = 0.d0
	      Fz_tot_uc(i) = 0.d0

	      Mx_tot_uc(i) = 0.d0
	      My_tot_uc(i) = 0.d0
	      Mz_tot_uc(i) = 0.d0

	   enddo

	endif

c	We first include the body force and moment from drag
c	Right now the induced velocity is not being determined (from the propeller)
	Vi = 0.d0
	Viw = Up + 2.d0*vi  ! Viw is a temporary use
	Fx = - 0.5d0 * rho * aircraft%X_fuseuu * abs(Viw)*Viw
	Fy = - 0.5d0 * rho * aircraft%Y_fusevv * abs(Vp )*Vp
	Fz = - 0.5d0 * rho * aircraft%Z_fuseww * abs(Wp )*Wp
c	It should be commented that this force is not symmetric
c	in the sense that for a moving sphere, where all fuse
c	are the same, the force is not in the opposit direction
c	of motion through the air.  For that to happen the
c	force needs to be something like - f(Vt)*vec(Up).
c	Notice this force is not of that form.  It is, in
c	some sense, a low angle of attack, low slideslip
c	form of a drag force.  Maybe we will change it at some point.
c	The example of a force that works:
c	Fx = - 0.5d0 * rho * aircraft%X_fuseuu * Vt*Viw
c	Fy = - 0.5d0 * rho * aircraft%Y_fusevv * Vt*Vp
c	Fz = - 0.5d0 * rho * aircraft%Z_fuseww * Vt*Wp

	Fx_tot = Fx_tot + Fx
	Fy_tot = Fy_tot + Fy
	Fz_tot = Fz_tot + Fz

	Fx_drag = Fx
	Fy_drag = Fy
	Fz_drag = Fz

c	temp = sqrt(Fx**2+Fy**2+Fz**2)
c	write (6,*) ' F direction   ',Fx/temp, Fy/temp
c     .                               ,Fz/temp

cc	write (6,*) ' Force2',Fx_tot,Fy_tot,Fz_tot
cc	write (6,*) ' rho ',rho
cc	write (6,*) ' fuse ',aircraft%X_fuseuu,Viw
cc	write (6,*) ' fuse ',aircraft%Y_fusevv,Vp
cc	write (6,*) ' fuse ',aircraft%Z_fuseww,Wp
	
c	Now the torque, where vec M =  vec x cross vec F
	x_w = aircraft%x_fuse - aircraft%x_cm
	y_w = aircraft%y_fuse - aircraft%y_cm
	z_w = aircraft%z_fuse - aircraft%z_cm

	Mx_tot = Mx_tot + y_w*Fz - z_w*Fy
	My_tot = My_tot + z_w*Fx - x_w*Fz
	Mz_tot = Mz_tot + x_w*Fy - y_w*Fx 

	if (compute_A) then
c	   Here x y z only depends on U V W, respectively.
	   Fx_U = - rho * aircraft%X_fuseuu * sign(1.d0,Viw)*Viw

	   Fy_V = - rho * aircraft%Y_fusevv * sign(1.d0,Vp )*Vp

	   Fz_W = - rho * aircraft%Z_fuseww * sign(1.d0,Wp )*Wp

	   Fx_tot_U = Fx_tot_U + Fx_U
	   Fy_tot_V = Fy_tot_V + Fy_V
	   Fz_tot_W = Fz_tot_W + Fz_W

	   Mx_tot_V = Mx_tot_V            - z_w*Fy_V
	   Mx_tot_W = Mx_tot_W + y_w*Fz_W

	   My_tot_U = My_tot_U + z_w*Fx_U
	   My_tot_W = My_tot_W            - x_w*Fz_W

	   Mz_tot_U = Mz_tot_U            - y_w*Fx_U 
	   Mz_tot_V = Mz_tot_V + x_w*Fy_V
	endif

	if (ip) write (6,100) 'body drag',0,1, 
     .                Fx_tot, Fy_tot, Fz_tot
     .               ,Mx_tot, My_tot, Mz_tot

	do i = 1, aircraft%num_wings

c	   We view the center of mass as the center of everything.
c	   In my coding right now, I don't know Vi (induced velocity from the 
c	   propeller) or Viw (induced velocity from the wing on tail elements),
c	   so right now those effects are not included.  jdw 3/22/21
	   Vi  = 0.d0
	   Viw = 0.d0

	   x_w = wing(i)%x - aircraft%x_cm
	   y_w = wing(i)%y - aircraft%y_cm
	   z_w = wing(i)%z - aircraft%z_cm

	   ic1 = wing(i)%ic1
	   ic2 = wing(i)%ic2

c	   vec(U V W)_w = vec(U V W) + vec(P Q R) x vec(x y z)_w + stuff
c	   Notice that the rotation vector omega = vec(P Q R), so this is
c	   u_w = up + omega cross x + stuff ... for motion in the rotating frame
	   U_w = Up + Q*z_w - R*y_w + wing(i)%eta*Vi
	   V_w = Vp + R*x_w - P*z_w
	   W_w = Wp + P*y_w - Q*x_w - wing(i)%eta_w*Viw

c	   We now compute the angle of attack, which includes an adjustment based
c          on the deflection of the ailerons and flaps.  We enforce the minimum 
c	   and maximum deflection of the control surfaces.
c
c	   delta_1
c     .      = min(
c     .            max(wing(i)%delta_amin, xcontrol(wing(i)%icontrol1))
c     .               ,wing(i)%delta_amax)
c	   delta_2
c     .      = min(
c     .            max(wing(i)%delta_amin, xcontrol(wing(i)%icontrol2))
c     .               ,wing(i)%delta_amax)
c
c	   delta_a = wing(i)%bias1*delta_1 + wing(i)%bias2*delta_2

c	   Okay, so now I am thinking this looks something like, where range
c	   limits and biases are all now built in,
	   if (ic1) icontrol1  = wing(i)%icontrol1
	   if (ic2) icontrol2  = wing(i)%icontrol2
	   if (ic1) then
	      delta_1 = wing(i)%ac1*uc(icontrol1) + wing(i)%bc1
	   else
	      delta_1 = 0.d0
	   endif
	   if (ic2) then
	      delta_2 = wing(i)%ac2*uc(icontrol2) + wing(i)%bc2
	   else
	      delta_2 = 0.d0
	   endif
	   delta_a    = delta_1 + delta_2
c	write (6,*) ' u1, u2, delta_1, delta_2 ',uc(icontrol1)
c     .              ,uc(icontrol2)
c     .              , delta_1, delta_2

c	   For the specific angle of attack, it matters what the orientation of the
c	   wing segment is.  For example, for a horizontal wing,
c	      alpha_w = atan2(W_w,U_w) + wing(i)%tau_a*delta_a
c	   while for a vertical stabilizer 
c	      alpha_w = atan2(V_w,U_w) + wing(i)%tau_a*delta_a
c	   Based on this, we have

	   VW_w       = W_w*wing(i)%nz+V_w*wing(i)%ny
	   alpha_w    = atan2(VW_w,U_w)                         ! body frame angle of attck
	   alpha_w_d  = atan2(VW_w,U_w) + wing(i)%tau_a*delta_a ! augment angle of attack including flaps
     .                + wing(i)%angle_of_incidence              ! only enters through this expression

c	   We point out what this says about the orientation.  Namely, if nz = 1 then
c	   the top of the wing is up (negative z direction - sorry it is inconsistant)
c	   and positive motion for the flap is down.  If ny = 1 then the top of the 
c	   wing is to the right (frd) and the positive motion of the rudder is to the left.
c	   What does that mean for a V tail?  If we are looking from the back to the front
c	   of the aircraft, and the V is at 45 degrees,  then if we want the top of the
c	   V tail surfaces to be "up" then for the left surface (ny,nz)=(1,1)/sqrt(2)
c	   and for the right tail we have (ny,nz)=(-1,1)/sqrt(2).  Now the motion
c	   of the rudder will be positive when it moves away from the center vertical plane.

c	   The range of atan2 is -pi to pi, so we adjust in case we are outside the range.
c	   It needs to wrap around.  This puts the angle of attack between -pi and pi.
	   if (alpha_w  .gt. pi) alpha_w   = alpha_w   - twopi
	   if (alpha_w  .lt.-pi) alpha_w   = alpha_w   + twopi
	   if (alpha_w_d.gt. pi) alpha_w_d = alpha_w_d - twopi
	   if (alpha_w_d.lt.-pi) alpha_w_d = alpha_w_d + twopi

c	   We have a qbar for the wing segment
c	   qbarprime_w = 0.5d0*rho*(U_w**2+W_w**2)  ! horizontal wing
c	   qbarprime_w = 0.5d0*rho*(U_w**2+V_w**2)  ! vertical stabilizer
	   qbarprime_w = 0.5d0*rho*(U_w**2 + VW_w**2) 

c	   The coefficients of lift and drag for each wing element.
c	   This is the only place the augmented angle of attack is used.
	   call CLDwing(alpha_w_d,wing(i)%a,wing(i)%C_L0,wing(i)%C_Lmin
     .         ,wing(i)%C_Lmax,wing(i)%C_D0,wing(i)%k,wing(i)%C_Dfp
     .         ,C_Lw,C_Dw,iflinear,dC_Ldalpha,dC_Ddalpha)

c	   Force equations  (Table 2.5-1 and page 650)
	   Force_Lift_w = wing(i)%surface_area*qbarprime_w*C_Lw
	   Force_Drag_w = wing(i)%surface_area*qbarprime_w*C_Dw

c	write (6,555) 
c     .               i, icontrol1, icontrol2, uc(icontrol1)
c     .             , uc(icontrol2), delta_1, delta_2
c     .             , Force_Lift_w, Force_Drag_w

 555	format(' wing i ic1 ic2 uc1 uc2 d1 d2 FL FD ',3i3,6f10.4)

c	   After much discussion, the decision is that this force should be
c	   projected back into the body frame based on the alpha of attack from the
c	   airflow and not on the "fake" angle of attack taking into account the flaps.
c	   (Thus, the equation in Stevens needs to be revised.)

c	   Now we get to resolving forces.  Now it matters whether we are 
c	   horizontal (as in the wing and the horizontal tail) or vertical.
c	   We determine whether we are a wing or a horizontal tail based
c	   on the cname of the wing part.

c	   By definition, the lift is perpendicular to the airflow and the
c	   drag is in opposite direction to the airflow.  Notice that the standard
c	   diagram is          FD       ^ FL
c                          < - _       |
c                                - _  |   
c                                    +-----> x
c                                    | - _   ) alpha  (positive angle of attach for pitch up)
c                                    |     - > 
c                                    V z       Vrel
c 
c          so Fx =   FL*sin(alpha) - FD*cos(alpha)      
c             Fz = - FL*cos(alpha) - FD*sin(alpha)
c
c	   Force equations (Table 2.5-1 and page 650)

	   Fx =   Force_Lift_w*sin(alpha_w) - Force_Drag_w*cos(alpha_w)
	   Fz = - Force_Lift_w*cos(alpha_w) - Force_Drag_w*sin(alpha_w)  ! temporary use
	   Fy =   Fz * wing(i)%ny  ! my approach to partitioning it in the event of V tails
	   Fz =   Fz * wing(i)%nz

cc	write (6,*) ' Fx Fy Fz ',Fx,Fy,Fz

	   Fx_tot = Fx_tot + Fx
	   Fy_tot = Fy_tot + Fy
	   Fz_tot = Fz_tot + Fz

c	   Now the torques, where vec M =  vec x cross vec F
	   Mx = y_w*Fz - z_w*Fy
	   My = z_w*Fx - x_w*Fz
	   Mz = x_w*Fy - y_w*Fx 

	   Mx_tot = Mx_tot + Mx
	   My_tot = My_tot + My
	   Mz_tot = Mz_tot + Mz

	   if (ip) write (6,100) 'wing lift and drag',i ,2,
     .                Fx, Fy, Fz, Mx, My, Mz

c	   To understand the forces

	   Fx =   Force_Lift_w*sin(alpha_w)
	   Fz = - Force_Lift_w*cos(alpha_w)  ! temporary use
	   Fy =   Fz * wing(i)%ny  ! my approach to partitioning it in the event of V tails
	   Fz =   Fz * wing(i)%nz

	   Fx_lift = Fx_lift + Fx
	   Fy_lift = Fy_lift + Fy
	   Fz_lift = Fz_lift + Fz

	   Fx = - Force_Drag_w*cos(alpha_w)
	   Fz = - Force_Drag_w*sin(alpha_w)  ! temporary use
	   Fy =   Fz * wing(i)%ny  ! my approach to partitioning it in the event of V tails
	   Fz =   Fz * wing(i)%nz

	   Fx_drag = Fx_drag + Fx
	   Fy_drag = Fy_drag + Fy
	   Fz_drag = Fz_drag + Fz

	   if (compute_A) then

c	      U_w = Up + Q*z_w - R*y_w + wing(i)%eta*Vi
c	      V_w = Vp + R*x_w - P*z_w
c	      W_w = Wp + P*y_w - Q*x_w - wing(i)%eta_w*Viw
c	      VW_w    = W_w*wing(i)%nz+V_w*wing(i)%ny
c	      alpha_w   = atan2(VW_w,U_w)
c	      alpha_w_d = atan2(VW_w,U_w) + wing(i)%tau_a*delta_a
c	      qbarprime_w = 0.5d0*rho*(U_w**2 + VW_w**2) 
c	      delta_1   = wing(i)%ac1*uc(icontrol1) + wing(i)%bc1
c	      delta_2   = wing(i)%ac2*uc(icontrol2) + wing(i)%bc2
c	      delta_a   = delta_1 + delta_2
c	      so alpha_w   = arctan(VW_w/U_w)
c	      so alpha_w_d = arctan(VW_w/U_w)+wing(i)%tau_a*delta_a
c	      The only uc dependenc is through Force through alpha_w_d
c	      Force_Lift_w = wing(i)%surface_area*qbarprime_w*C_Lw
c	      Force_Drag_w = wing(i)%surface_area*qbarprime_w*C_Dw
c	      Fx =   Force_Lift_w*sin(alpha_w) - Force_Drag_w*cos(alpha_w)
c	      Fz = - Force_Lift_w*cos(alpha_w) - Force_Drag_w*sin(alpha_w)  ! temporary use
c	      Fy =   Fz * wing(i)%ny  ! my approach to partitioning it in the event of V tails
c	      Fz =   Fz * wing(i)%nz

	      dqbardU = rho *  U_w
	      dqbardV = rho * VW_w * wing(i)%ny
	      dqbardW = rho * VW_w * wing(i)%nz

c	      Notice the arctan(VW/U) is an issue for when velocity components are zero, which
c	      happens in hover trim state.

c	      d arctan (x/a) / dx = a / (a^2 + x^2), d arctan(a/x) / dx = - a / (a^2 + x^2)
	      if ((abs(U_w).lt.1.d-09).and.(abs(VW_w).lt.1.d-09)) then   ! presumably this is a small number, like zero
	         dalphadU = 0.d0
	         dalphadV = 0.d0
	         dalphadW = 0.d0
	      else
	         dalphadU = -  VW_w/(U_w**2 + VW_w**2)
	         dalphadV =   ( U_w/(U_w**2 + VW_w**2) ) * wing(i)%ny
	         dalphadW =   ( U_w/(U_w**2 + VW_w**2) ) * wing(i)%nz
	      endif
	      if (ic1) then  ! means d alpha_w_d d uc1 - not d alpha_w d uc1 which equals zero.
		 dalphaduc1 = wing(i)%ac1 * wing(i)%tau_a
	      else
		 dalphaduc1 = 0.d0
	      endif
	      if (ic2) then
	         dalphaduc2 = wing(i)%ac2 * wing(i)%tau_a
	      else
		 dalphaduc2 = 0.d0
	      endif

c	      When induced velocity is included, there will be additional terms.

c	      Now we begin building the various terms - beware the _W is a double use (not the same as above)
c	      Force_Lift_w = wing(i)%surface_area*qbarprime_w*C_Lw
c	      Force_Drag_w = wing(i)%surface_area*qbarprime_w*C_Dw
	      Force_L_U = wing(i)%surface_area * ( dqbardU*C_Lw
     .                     + qbarprime_w * dC_Ldalpha * dalphadU )
	      Force_L_V = wing(i)%surface_area * ( dqbardV*C_Lw
     .                     + qbarprime_w * dC_Ldalpha * dalphadV )
	      Force_L_W = wing(i)%surface_area * ( dqbardW*C_Lw
     .                     + qbarprime_w * dC_Ldalpha * dalphadW )
	      Force_D_U = wing(i)%surface_area * ( dqbardU*C_Dw
     .                     + qbarprime_w * dC_Ddalpha * dalphadU )
	      Force_D_V = wing(i)%surface_area * ( dqbardV*C_Dw
     .                     + qbarprime_w * dC_Ddalpha * dalphadV )
	      Force_D_W = wing(i)%surface_area * ( dqbardW*C_Dw
     .                     + qbarprime_w * dC_Ddalpha * dalphadW )

	      Fx_U =   Force_L_U*sin(alpha_w) - Force_D_U*cos(alpha_w)
     .         - (Force_Lift_w*cos(alpha_w) - Force_Drag_w*sin(alpha_w))
     .           * dalphadU
	      Fx_V =   Force_L_V*sin(alpha_w) - Force_D_V*cos(alpha_w)
     .         - (Force_Lift_w*cos(alpha_w) - Force_Drag_w*sin(alpha_w))
     .           * dalphadV
	      Fx_W =   Force_L_W*sin(alpha_w) - Force_D_W*cos(alpha_w)
     .         - (Force_Lift_w*cos(alpha_w) - Force_Drag_w*sin(alpha_w))
     .           * dalphadW

	      Fz_U = - Force_L_U*cos(alpha_w) - Force_D_U*sin(alpha_w)  ! temporary use
     .         + (Force_Lift_w*sin(alpha_w) - Force_Drag_w*cos(alpha_w))
     .           * dalphadU
	      Fz_V = - Force_L_V*cos(alpha_w) - Force_D_V*sin(alpha_w)  ! temporary use
     .         + (Force_Lift_w*sin(alpha_w) - Force_Drag_w*cos(alpha_w))
     .           * dalphadV
	      Fz_W = - Force_L_W*cos(alpha_w) - Force_D_W*sin(alpha_w)  ! temporary use
     .         + (Force_Lift_w*sin(alpha_w) - Force_Drag_w*cos(alpha_w))
     .           * dalphadW

	      Fy_U =   Fz_U * wing(i)%ny
	      Fy_V =   Fz_V * wing(i)%ny
	      Fy_W =   Fz_W * wing(i)%ny

	      Fz_U =   Fz_U * wing(i)%nz
	      Fz_V =   Fz_V * wing(i)%nz
	      Fz_W =   Fz_W * wing(i)%nz


	      if (ic1) then
	         Force_L_uc1 = wing(i)%surface_area * 
     .                          qbarprime_w * dC_Ldalpha * dalphaduc1
	         Force_D_uc1 = wing(i)%surface_area *
     .                          qbarprime_w * dC_Ddalpha * dalphaduc1
	      else
		 Force_L_uc1 = 0.d0
		 Force_D_uc1 = 0.d0
	      endif
	      if (ic2) then
	         Force_L_uc2 = wing(i)%surface_area *
     .                          qbarprime_w * dC_Ldalpha * dalphaduc2
	         Force_D_uc2 = wing(i)%surface_area *
     .                          qbarprime_w * dC_Ddalpha * dalphaduc2
	      else
		 Force_L_uc2 = 0.d0
		 Force_D_uc2 = 0.d0
	      endif

c	      Notice that dalphaduc1 means dalpha_w_d duc1 so there are no 
c	      d sin(alpha_w) = cos(alpha) dalpha terms.
	      Fx_uc1 = Force_L_uc1*sin(alpha_w)-Force_D_uc1*cos(alpha_w)
	      Fx_uc2 = Force_L_uc2*sin(alpha_w)-Force_D_uc2*cos(alpha_w)

	      Fz_uc1 =-Force_L_uc1*cos(alpha_w)-Force_D_uc1*sin(alpha_w)  ! temporary use
	      Fz_uc2 =-Force_L_uc2*cos(alpha_w)-Force_D_uc2*sin(alpha_w)  ! temporary use

	      Fy_uc1 =   Fz_uc1 * wing(i)%ny
	      Fy_uc2 =   Fz_uc2 * wing(i)%ny

	      Fz_uc1 =   Fz_uc1 * wing(i)%nz
	      Fz_uc2 =   Fz_uc2 * wing(i)%nz

c	      U_w = Up + Q*z_w - R*y_w + wing(i)%eta*Vi
c	      V_w = Vp + R*x_w - P*z_w
c	      W_w = Wp + P*y_w - Q*x_w - wing(i)%eta_w*Viw

	      Force_L_P = wing(i)%surface_area * C_Lw 
     .                      * ( - dqbardV * z_w + dqbardW * y_w)
	      Force_L_Q = wing(i)%surface_area * C_Lw
     .                      * (   dqbardU * z_w - dqbardW * x_w)
	      Force_L_R = wing(i)%surface_area * C_Lw
     .                      * ( - dqbardU * y_w + dqbardV * x_w)
	      Force_D_P = wing(i)%surface_area * C_Dw
     .                      * ( - dqbardV * z_w + dqbardW * y_w)
	      Force_D_Q = wing(i)%surface_area * C_Dw
     .                      * (   dqbardU * z_w - dqbardW * x_w)
	      Force_D_R = wing(i)%surface_area * C_Dw
     .                      * ( - dqbardU * y_w + dqbardV * x_w)

	      dUdP =   0.d0
	      dUdQ =   z_w
	      dUdR = - y_w

	      dVdP = - z_w
	      dVdQ =   0.d0
	      dVdR =   x_w

	      dWdP =   y_w
	      dWdQ = - x_w
	      dWdR =   0.d0

	      Fx_P =   Force_L_P*sin(alpha_w) - Force_D_P*cos(alpha_w)
     .         - (Force_Lift_w*cos(alpha_w) - Force_Drag_w*sin(alpha_w))
     .           * (dalphadU * dUdP + dalphadV * dVdP + dalphadW * dWdP)
	      Fx_Q =   Force_L_Q*sin(alpha_w) - Force_D_Q*cos(alpha_w)
     .         - (Force_Lift_w*cos(alpha_w) - Force_Drag_w*sin(alpha_w))
     .           * (dalphadU * dUdQ + dalphadV * dVdQ + dalphadW * dWdQ)
	      Fx_R =   Force_L_R*sin(alpha_w) - Force_D_R*cos(alpha_w)
     .         - (Force_Lift_w*cos(alpha_w) - Force_Drag_w*sin(alpha_w))
     .           * (dalphadU * dUdR + dalphadV * dVdR + dalphadW * dWdR)

	      Fz_P = - Force_L_P*cos(alpha_w) - Force_D_P*sin(alpha_w)  ! temporary use
     .         + (Force_Lift_w*sin(alpha_w) - Force_Drag_w*cos(alpha_w))
     .           * (dalphadU * dUdP + dalphadV * dVdP + dalphadW * dWdP)
	      Fz_Q = - Force_L_Q*cos(alpha_w) - Force_D_Q*sin(alpha_w)  ! temporary use
     .         + (Force_Lift_w*sin(alpha_w) - Force_Drag_w*cos(alpha_w))
     .           * (dalphadU * dUdQ + dalphadV * dVdQ + dalphadW * dWdQ)
	      Fz_R = - Force_L_R*cos(alpha_w) - Force_D_R*sin(alpha_w)  ! temporary use
     .         + (Force_Lift_w*sin(alpha_w) - Force_Drag_w*cos(alpha_w))
     .           * (dalphadU * dUdR + dalphadV * dVdR + dalphadW * dWdR)

	      Fy_P =   Fz_P * wing(i)%ny
	      Fy_Q =   Fz_Q * wing(i)%ny
	      Fy_R =   Fz_R * wing(i)%ny

	      Fz_P =   Fz_P * wing(i)%nz
	      Fz_Q =   Fz_Q * wing(i)%nz
	      Fz_R =   Fz_R * wing(i)%nz

	      Mx_U = y_w*Fz_U - z_w*Fy_U
	      My_U = z_w*Fx_U - x_w*Fz_U
	      Mz_U = x_w*Fy_U - y_w*Fx_U 

	      Mx_V = y_w*Fz_V - z_w*Fy_V
	      My_V = z_w*Fx_V - x_w*Fz_V
	      Mz_V = x_w*Fy_V - y_w*Fx_V 

	      Mx_W = y_w*Fz_W - z_w*Fy_W
	      My_W = z_w*Fx_W - x_w*Fz_W
	      Mz_W = x_w*Fy_W - y_w*Fx_W 

	      Mx_uc1 = y_w*Fz_uc1 - z_w*Fy_uc1
	      My_uc1 = z_w*Fx_uc1 - x_w*Fz_uc1
	      Mz_uc1 = x_w*Fy_uc1 - y_w*Fx_uc1 

	      Mx_uc2 = y_w*Fz_uc2 - z_w*Fy_uc2
	      My_uc2 = z_w*Fx_uc2 - x_w*Fz_uc2
	      Mz_uc2 = x_w*Fy_uc2 - y_w*Fx_uc2 

	      Mx_P = y_w*Fz_P - z_w*Fy_P
	      My_P = z_w*Fx_P - x_w*Fz_P
	      Mz_P = x_w*Fy_P - y_w*Fx_P 

	      Mx_Q = y_w*Fz_Q - z_w*Fy_Q
	      My_Q = z_w*Fx_Q - x_w*Fz_Q
	      Mz_Q = x_w*Fy_Q - y_w*Fx_Q 

	      Mx_R = y_w*Fz_R - z_w*Fy_R
	      My_R = z_w*Fx_R - x_w*Fz_R
	      Mz_R = x_w*Fy_R - y_w*Fx_R 

c	      We increment the totals
	      Fx_tot_U = Fx_tot_U + Fx_U
	      Fx_tot_V = Fx_tot_V + Fx_V
	      Fx_tot_W = Fx_tot_W + Fx_W

	      Fy_tot_U = Fy_tot_U + Fy_U
	      Fy_tot_V = Fy_tot_V + Fy_V
	      Fy_tot_W = Fy_tot_W + Fy_W

	      Fz_tot_U = Fz_tot_U + Fz_U
	      Fz_tot_V = Fz_tot_V + Fz_V
	      Fz_tot_W = Fz_tot_W + Fz_W
	
	      Fx_tot_P = Fx_tot_P + Fx_P
	      Fx_tot_Q = Fx_tot_Q + Fx_Q
	      Fx_tot_R = Fx_tot_R + Fx_R

	      Fy_tot_P = Fy_tot_P + Fy_P
	      Fy_tot_Q = Fy_tot_Q + Fy_Q
	      Fy_tot_R = Fy_tot_R + Fy_R

	      Fz_tot_P = Fz_tot_P + Fz_P
	      Fz_tot_Q = Fz_tot_Q + Fz_Q
	      Fz_tot_R = Fz_tot_R + Fz_R
	
	      if(ic1)Fx_tot_uc(icontrol1)= Fx_tot_uc(icontrol1) + Fx_uc1
	      if(ic2)Fx_tot_uc(icontrol2)= Fx_tot_uc(icontrol2) + Fx_uc2

	      if(ic1)Fy_tot_uc(icontrol1)= Fy_tot_uc(icontrol1) + Fy_uc1
	      if(ic2)Fy_tot_uc(icontrol2)= Fy_tot_uc(icontrol2) + Fy_uc2

	      if(ic1)Fz_tot_uc(icontrol1)= Fz_tot_uc(icontrol1) + Fz_uc1
	      if(ic2)Fz_tot_uc(icontrol2)= Fz_tot_uc(icontrol2) + Fz_uc2

	      Mx_tot_U = Mx_tot_U + Mx_U
	      Mx_tot_V = Mx_tot_V + Mx_V
	      Mx_tot_W = Mx_tot_W + Mx_W

	      My_tot_U = My_tot_U + My_U
	      My_tot_V = My_tot_V + My_V
	      My_tot_W = My_tot_W + My_W

	      Mz_tot_U = Mz_tot_U + Mz_U
	      Mz_tot_V = Mz_tot_V + Mz_V
	      Mz_tot_W = Mz_tot_W + Mz_W
	
	      Mx_tot_P = Mx_tot_P + Mx_P
	      Mx_tot_Q = Mx_tot_Q + Mx_Q
	      Mx_tot_R = Mx_tot_R + Mx_R

	      My_tot_P = My_tot_P + My_P
	      My_tot_Q = My_tot_Q + My_Q
	      My_tot_R = My_tot_R + My_R

	      Mz_tot_P = Mz_tot_P + Mz_P
	      Mz_tot_Q = Mz_tot_Q + Mz_Q
	      Mz_tot_R = Mz_tot_R + Mz_R
	
	      if(ic1)Mx_tot_uc(icontrol1)= Mx_tot_uc(icontrol1) + Mx_uc1
	      if(ic2)Mx_tot_uc(icontrol2)= Mx_tot_uc(icontrol2) + Mx_uc2

	      if(ic1)My_tot_uc(icontrol1)= My_tot_uc(icontrol1) + My_uc1
	      if(ic2)My_tot_uc(icontrol2)= My_tot_uc(icontrol2) + My_uc2

	      if(ic1)Mz_tot_uc(icontrol1)= Mz_tot_uc(icontrol1) + Mz_uc1
	      if(ic2)Mz_tot_uc(icontrol2)= Mz_tot_uc(icontrol2) + Mz_uc2

	   endif  ! to control_A

	enddo  ! loop on i over the wing parts
c

c	Now we do the loop over the propellers, computing thrust force and 
c	related moments, including gyroscopic moments.
c	It is assumed each propeller has its own electric motor and that it
c	is direct drive (no gearing).

c	We zero out the battery amperage
	do i = 1, aircraft%num_batteries
	   xo(29 + 2*aircraft%num_propellers + i) = 0.d0
	enddo

	do i = 1, aircraft%num_propellers

c	   We go to the table to get the information of interest.
c	   First we need the velocity in the direction of travel for the
c	   propeller so that we can then compute the advanced ratio.

	   D = 2.d0*propeller(i)%radius
	   n = motor_omega(i)/twopi  ! to get reveolutions per second from radians per second
	   Vpath = Up*propeller(i)%nx + Vp*propeller(i)%ny
     .           + Wp*propeller(i)%nz

	   if (n.gt.0.d0) then  ! the propeller is spinning

c	      J could conceivably be negative, which is not in the tables,
c	      so in that case we set it to zero.  (This could be due to me not
c	      computing an induced velocity.)
	      if (Vpath.lt.0.d0) Vpath = 0.d0

	      J = Vpath / (n*D) ! J = V/nD

cold	      call table (propeller(i)%prop_num_pts, J, propeller(i)%J
cold     .                , dxkm1, dxk, k, ierror)
cold	      Ct = dxkm1*propeller(i)%Ct(k-1) + dxk*propeller(i)%Ct(k)
cold	      Cp = dxkm1*propeller(i)%Cp(k-1) + dxk*propeller(i)%Cp(k)

	      call table (propeller(i)%prop_num_pts_J, J, propeller(i)%J
     .                , dxkm1, dxk, k, kpast, ierror)
	      call table (propeller(i)%prop_num_pts_omega,motor_omega(i)
     .                , propeller(i)%omega
     .                , dylm1, dyl, l, lpast, ierror)
c	      I believe the first index is for J, the second index is for omega
	      Cp11 = propeller(i)%Ct(k-1,l-1)
	      Cp21 = propeller(i)%Ct(k-1,l  )
	      Cp22 = propeller(i)%Ct(k  ,l  )
	      Cp12 = propeller(i)%Ct(k  ,l-1)
	      Ct   = Cp11*dxkm1*dylm1 + Cp21*dxkm1*dyl
     .             + Cp22*dxk  *dyl   + Cp12*dxk  *dylm1
	      Cp11 = propeller(i)%Cp(k-1,l-1)
	      Cp21 = propeller(i)%Cp(k-1,l  )
	      Cp22 = propeller(i)%Cp(k  ,l  )
	      Cp12 = propeller(i)%Cp(k  ,l-1)
	      Cp   = Cp11*dxkm1*dylm1 + Cp21*dxkm1*dyl
     .             + Cp22*dxk  *dyl   + Cp12*dxk  *dylm1


	      Thrust = Ct*rho*n**2*D**4      ! T = Ct*rho*n^2*D^4
	      Power  = Cp*rho*n**3*D**5      ! P = Cp*rho*n^3*D^5
              Torque = Power/motor_omega(i)  ! Q = P/motor omega

	   else

cold	      Ct = propeller(i)%Ct(1)
cold	      Cp = propeller(i)%Cp(1)

	      Ct = propeller(i)%Ct(1,1)
	      Cp = propeller(i)%Cp(1,1)

	      Thrust = 0.d0
	      Power  = 0.d0
              Torque = 0.d0

	   endif

c	   We now compute the work the motor is doing.

	   voltage = propeller(i)%ac*uc(propeller(i)%icontrol)
     .             + propeller(i)%bc ! voltage signal from the controls.

	   Kv   = twopi*propeller(i)%Kv/60.d0  ! unit change from rpm/volt to rad per sec/volt
	   Rw   = propeller(i)%Rw

	   if (voltage.gt.0.d0) then

	      amps = (voltage - motor_omega(i)/Kv)/Rw

c	      This equation is an issue if the controls are fluctuating alot,
c	      as it can result in negative current.  The underlying physics problem is
c	      probably that I am not integrating the current, which would include the 
c	      inductance of the motor, rather I am using this equation to obtain it.
c	      Perhaps in the future that will be addressed. JDW  7/26/21

cold	      amps = (Power/(propeller(i)%efficiency*voltage)) 
cold     .             + propeller(i)%I_idle  ! to change to watts

c	      Notice the voltage is not the battery voltage, it is the voltage that
c	      is being supplied by the ESC.

c	      We store and add up estimated power usage for analysis purposes
c	      There is an efficiency factor here due to the ESC
	      xo(29 + i)  = voltage*amps/propeller(i)%efficiency_ESC ! this is the individual motor power
	      xo(29 +   aircraft%num_propellers + i)
     .              = amps/propeller(i)%efficiency_ESC ! amps through that motor and ESC
	      xo(29 + 2*aircraft%num_propellers + propeller(i)%ibattery)
     .      = xo(29 + 2*aircraft%num_propellers + propeller(i)%ibattery)
     .              + amps/propeller(i)%efficiency_ESC   ! amps out of that battery
	      Power_tot = Power_tot + xo(29 + i)  ! total power

	      motor_efficiency = Power/(voltage*amps)  ! in case someone wants to know - right now not stored

	   endif


c	   Minus sign for the torque, as it is in the opposite direction
c	   of the propeller spin.

	   Mx = - Torque * propeller(i)%spin*propeller(i)%nx
	   My = - Torque * propeller(i)%spin*propeller(i)%ny
	   Mz = - Torque * propeller(i)%spin*propeller(i)%nz

	   Mx_tot = Mx_tot + Mx
	   My_tot = My_tot + My
	   Mz_tot = Mz_tot + Mz

	   if (ip) write (6,100) 'torque from prop rotation',
     .        i, 3, 0.d0, 0.d0, 0.d0, Mx, My, Mz

c	   Now we can distribute the thrust force

           Fx = Thrust * propeller(i)%nx
           Fy = Thrust * propeller(i)%ny
           Fz = Thrust * propeller(i)%nz

	   Fx_tot = Fx_tot + Fx
	   Fy_tot = Fy_tot + Fy
	   Fz_tot = Fz_tot + Fz

c	   Save up a diagnostic
	   Thrust_x = Thrust_x + Fx
	   Thrust_y = Thrust_y + Fy
	   Thrust_z = Thrust_z + Fz

c	   Now the torque from the thrust (sometimes called the diferential thrust moment Mdt)

	   x_w = propeller(i)%x - aircraft%x_cm
	   y_w = propeller(i)%y - aircraft%y_cm
	   z_w = propeller(i)%z - aircraft%z_cm

	   Mx  = y_w*Fz - z_w*Fy
	   My  = z_w*Fx - x_w*Fz
	   Mz  = x_w*Fy - y_w*Fx

	   Mx_tot = Mx_tot + Mx
	   My_tot = My_tot + My
	   Mz_tot = Mz_tot + Mz

	   if (ip) write (6,100) 'prop thrust/torque',i,4, 
     .                Fx, Fy, Fz, Mx, My, Mz

c	   Now we compute other terms involved in the moment, which have to do
c	   with the fact that we have a spinning motor and propeller, which has
c	   rotational inertia.  The first term is due to changing spin rate
c	   and the second is the gyroscopic moment.  Here we need to pay attention
c	   to the direction of rotation.  We use the center of the propeller as
c	   the spatial location.  The rotational direction is along the perpendicular.

c	   This is due to the gyroscopic behavior of the entire frame, where we 
c	   are subtracting vec (P Q R) x I motor_omega
c	   We are abusing Fi here, to avoid defining another variable.

	   temp = propeller(i)%Ir * propeller(i)%spin * motor_omega(i) ! propeller and motor angular momentum
	   Fx = temp * propeller(i)%nx
	   Fy = temp * propeller(i)%ny
	   Fz = temp * propeller(i)%nz

	   Mx = - (Q*Fz - R*Fy)
	   My = - (R*Fx - P*Fz)
	   Mz = - (P*Fy - Q*Fx)

	   Mx_tot = Mx_tot + Mx
	   My_tot = My_tot + My
	   Mz_tot = Mz_tot + Mz
	   
	   if (ip) write (6,100) 'prop gyroscopic',i,6,
     .                0.d0, 0.d0, 0.d0
     .               ,Mx, My, Mz

c	   Now we need to look at updating the spin rate, i.e., this the motor region.
c	   We are assuming the this is a speed controlled response, so that the speed
c	   is directly proportional to the voltage.  We need to find, for the current
c	   airspeed, what the target omega is based on the requested voltage.  It turns
c	   out to require the solution of a cubic, which we do through a Newton iteration.

	   Kt  = propeller(i)%Kt
	   Kv  = twopi*propeller(i)%Kv/60.d0  ! unit change from rpm/volt to rad per sec/volt
	   Rw  = propeller(i)%Rw
	   I0  = propeller(i)%I_idle
	   voltage = propeller(i)%ac*uc(propeller(i)%icontrol)
     .             + propeller(i)%bc ! voltage signal from the controls.

	   a   =   (Rw/Kt)*(Cp*rho*D**5/twopi**3)
	   b   =   1.d0/Kv
	   c   = - voltage + Rw*I0

	   omega = (-b+sqrt(b**2-4.d0*a*c))/(2.d0*a)  ! this is the requested rotation rate

c	   If the control exactly is zero, the propeller is not spinning.

	   if (uc(propeller(i)%icontrol).le.0.d0) then
	      omega = 0.d0
	      goto 200
	   endif

c	   omega0 = omega
c	   omega1 = omega

c	********** this is a test section of code *********

c	   The loop converges okay in about 10 iterations.  Can it be done faster?
c
c	   write (6,*) ' Omega ', 1, omega,cp,k,dxkm1,dxk
c	   do jdw = 2, 25
c
c	   n = omega/twopi  ! to get reveolutions per second from radians per second
c
c	      J = Vpath / (n*D) ! J = V/nD
c
c	      call table (propeller(i)%prop_num_pts, J, propeller(i)%J
c     .                , dxkm1, dxk, k, ierror)
c	      Ct = dxkm1*propeller(i)%Ct(k-1) + dxk*propeller(i)%Ct(k)
c	      Cp = dxkm1*propeller(i)%Cp(k-1) + dxk*propeller(i)%Cp(k)
c
c	   a   =   (Rw/Kt)*(Cp*rho*D**5/twopi**3)/1.d7   ! the unit change is cgs to MKS
c	   b   =   1.d0/Kv
c	   c   = - voltage
c
c	   deltaJ = propeller(i)%J(k)-propeller(i)%J(k-1)
c
c	   a1  = (a/Cp)*(propeller(i)%Cp(k-1)*propeller(i)%J(k  )
c     .             -propeller(i)%Cp(k  )*propeller(i)%J(k-1))/deltaJ
c	   b1  = b + twopi*Vpath*(a/Cp)*
c     .            (propeller(i)%Cp(k  )-propeller(i)%Cp(k-1))/(deltaJ*D)
c
c	   omega = (-b+sqrt(b**2-4.d0*a*c))/(2.d0*a)  ! this is the requested rotation rate
c	   omega2 = (-b1+sqrt(b1**2-4.d0*a1*c))/(2.d0*a1)  ! this is the requested rotation rate
c	   write (6,*) ' Omega ',jdw,omega, omega2
c
cc	   if (mod(jdw,2).eq.0) omega = 0.5d0*(omega+omega0)
cc	   omega0 = omega
c
c	   enddo

c	   Now we try a pseudo Newton method
c	   omega = omega1
c
c	   This is the original one-dimensional table version.
c
c	   do jdw = 1, 5  ! takes 5 iterates to converge, based on examples
cc
c	      n = omega/twopi  ! to get reveolutions per second from radians per second
c
c	      J = Vpath / (n*D) ! J = V/nD
c
c	      call table (propeller(i)%prop_num_pts, J, propeller(i)%J
c     .                , dxkm1, dxk, k, ierror)
c
c	      Same as above, only this a is a/Cp above
c	      a   =   (Rw/Kt)*(rho*D**5/twopi**3)   ! the unit change is cgs to MKS
c	      b   =   1.d0/Kv
c	      c   = - voltage
c	      deltaJ = propeller(i)%J(k)-propeller(i)%J(k-1)
c
c	      a1  = a*(propeller(i)%Cp(k-1)*propeller(i)%J(k  )
c     .             -propeller(i)%Cp(k  )*propeller(i)%J(k-1))/deltaJ
c	      b1  = b + twopi*Vpath*a*
c     .            (propeller(i)%Cp(k  )-propeller(i)%Cp(k-1))/(deltaJ*D)
c
c	      ff = a1*omega**2+b1*omega+c
c	      ffp = 2.d0*a1*omega + b1
c	      omega = omega - ff/ffp
c
cc	   write (6,*) ' Omega ',jdw,omega
c
c	   enddo


c	   Now we try a Newton method
c	   omega = omega1

	   do jdw = 1, 5  ! takes 5 iterates to converge, based on examples

	      n = omega/twopi  ! to get reveolutions per second from radians per second

	      J = Vpath / (n*D) ! J = V/nD

	      call table (propeller(i)%prop_num_pts_J, J, propeller(i)%J
     .                , dxkm1, dxk, k, kpast, ierror)
	      call table (propeller(i)%prop_num_pts_omega, omega
     .                , propeller(i)%omega
     .                , dylm1, dyl, l, lpast, ierror)
c	      I believe the first index is for J, the second index is for omega
	      J1   = propeller(i)%J    (k-1)
	      J2   = propeller(i)%J    (k  )
	      O1   = propeller(i)%omega(l-1)
	      O2   = propeller(i)%omega(l  )
	      Cp11 = propeller(i)%Cp(k-1,l-1)
	      Cp21 = propeller(i)%Cp(k-1,l  )
	      Cp22 = propeller(i)%Cp(k  ,l  )
	      Cp12 = propeller(i)%Cp(k  ,l-1)
	      aa =   Cp11       - Cp21       + Cp22       - Cp12
	      bb = - Cp11   *J2 + Cp21   *J2 - Cp22   *J1 + Cp12   *J1
	      cc = - Cp11*O2    + Cp21*O1    - Cp22*O1    + Cp12*O2
	      dd =   Cp11*O2*J2 - Cp21*O1*J2 + Cp22*O1*J1 - Cp12*O2*J1

c	      Same as above, only this a is a/Cp above
	      a   =   (Rw/Kt)*(rho*D**5/twopi**3) 
	      b   =   1.d0/Kv
	      c   = - voltage + Rw*I0

	      ee = a/((J2 - J1)*(O2 - O1))

	      a3 = ee*bb
	      a2 = ee*(aa*twopi*Vpath/D+dd)
	      a1 = ee* cc*twopi*Vpath/D + b
	      a0 = c

	      ff    =     a3*omega**3 +      a2*omega**2 + a1*omega + a0
	      ffp   =3.d0*a3*omega**2 + 2.d0*a2*omega    + a1
	      omega = omega - ff/ffp

cc	   write (6,*) ' Omega ',jdw,omega

	   enddo

c	   This is where we have some information regarding the relationship between
c	   the voltage and the spin, though it may be slightly out of sinc with
c	   motor_omega, as it relates to the requested spin this cycle vs. the
c	   actual spin from the state variable, which may not be quite the same,
c	   especially if the controls are erratic (hopefully not, of course).
	   if (compute_A) then
	      domegadVoltage = 1.d0/ffp
	      n = omega/twopi  ! to get reveolutions per second from radians per second
	      J = Vpath / (n*D) ! J = V/nD
	      dJdomega = - twopi*Vpath/(D*omega**2)

	      denom = (O2 - O1)*(J2 - J1)
c	      Cp = (Cp11*(omega - O2)*(J-J2) - Cp21*(omega - O1)*(J-J2)
c     .            + Cp22*(omega - O1)*(J-J1) - Cp12*(omega - O2)*(J-J1))
c     .            / denom
c	      dCpdJ = (Cp11*(omega - O2) - Cp21*(omega - O1)
c     .               + Cp22*(omega - O1) - Cp12*(omega - O2))
c     .               / denom
c	      dCpdomega = (Cp11*(J-J2) - Cp21*(J-J2)
c     .                   + Cp22*(J-J1) - Cp12*(J-J1))
c     .                   / denom

c	      Using the Cpij for Ctij
	      Cp11 = propeller(i)%Ct(k-1,l-1)
	      Cp21 = propeller(i)%Ct(k-1,l  )
	      Cp22 = propeller(i)%Ct(k  ,l  )
	      Cp12 = propeller(i)%Ct(k  ,l-1)
	      Ct = (Cp11*(omega - O2)*(J-J2) - Cp21*(omega - O1)*(J-J2)
     .            + Cp22*(omega - O1)*(J-J1) - Cp12*(omega - O2)*(J-J1))
     .            / denom
	      dCtdJ = (Cp11*(omega - O2) - Cp21*(omega - O1)
     .               + Cp22*(omega - O1) - Cp12*(omega - O2))
     .               / denom
	      dCtdomega = (Cp11*(J-J2) - Cp21*(J-J2)
     .                   + Cp22*(J-J1) - Cp12*(J-J1))
     .                   / denom

	      dThrustdVoltage =  ((rho * D**4 * omega) / (twopi**2) )
     .              * ((dCtdJ*dJdomega+dCtdomega) * omega + Ct )
     .                   * domegadVoltage

	      dThrustduc1 = dThrustdVoltage * propeller(i)%ac

c	      Vpath = Up*propeller(i)%nx + Vp*propeller(i)%ny
c       .           + Wp*propeller(i)%nz
	      dThrustdV   = rho * n * D**3 * dCtdJ

	   endif


 200	   continue ! the case where the propeller is not spinning

	   if (compute_A.and.(omega.eq.0.d0)) then

c	      Are we going to put anything here?  Or at least, anything complicated?
c	      I'm tempted to just ignore the domegadVoltage term for
c	      the zero rotation case - we should rarely be there.

	      dThrustdVoltage = 0.d0

	      dThrustduc1 = dThrustdVoltage * propeller(i)%ac

	      dThrustdV   = 0.d0

	   endif


c	******** end of the test section *****

c	   Now that we have determined the omega based on this voltage,
c	   our first order ODE to head towards the the requested omega is
	   motor_omegadot(i) = (omega - motor_omega(i))
     .                    /propeller(i)%tau_ESC


c	   This is due to the change in spin rate (subtracted since it appears on
c	   the other side of the equation, so that it becomes an applied moment).

	   temp= propeller(i)%Ir * propeller(i)%spin * motor_omegadot(i)

	   Mx = - temp * propeller(i)%nx
	   My = - temp * propeller(i)%ny
	   Mz = - temp * propeller(i)%nz

	   Mx_tot = Mx_tot + Mx
	   My_tot = My_tot + My
	   Mz_tot = Mz_tot + Mz

	   if (ip) write (6,100) 'prop acceleration',i,5, 
     .               0.d0, 0.d0, 0.d0
     .               ,Mx, My, Mz

c	   We have the collection of things for the control matrices (we need to decide how much
c	   is here).
	   if (compute_A) then

c	      First we do the most important terms, those for Blin.

c	      Minus sign for the torque, as it is in the opposite direction
c	      of the propeller spin.

	      icontrol1 = propeller(i)%icontrol

	      Mx = - dThrustduc1 * propeller(i)%spin*propeller(i)%nx
	      My = - dThrustduc1 * propeller(i)%spin*propeller(i)%ny
	      Mz = - dThrustduc1 * propeller(i)%spin*propeller(i)%nz

	      Mx_tot_uc(icontrol1) = Mx_tot_uc(icontrol1) + Mx
	      My_tot_uc(icontrol1) = My_tot_uc(icontrol1) + My
	      Mz_tot_uc(icontrol1) = Mz_tot_uc(icontrol1) + Mz

c	      Now we can distribute the thrust force

              Fx = dThrustduc1 * propeller(i)%nx
              Fy = dThrustduc1 * propeller(i)%ny
              Fz = dThrustduc1 * propeller(i)%nz

	      Fx_tot_uc(icontrol1) = Fx_tot_uc(icontrol1) + Fx
	      Fy_tot_uc(icontrol1) = Fy_tot_uc(icontrol1) + Fy
	      Fz_tot_uc(icontrol1) = Fz_tot_uc(icontrol1) + Fz

c	      Now the torque from the thrust (sometimes called the diferential thrust moment Mdt)

	      Mx = y_w*Fz - z_w*Fy
	      My = z_w*Fx - x_w*Fz
	      Mz = x_w*Fy - y_w*Fx

	      Mx_tot_uc(icontrol1) = Mx_tot_uc(icontrol1) + Mx
	      My_tot_uc(icontrol1) = My_tot_uc(icontrol1) + My
	      Mz_tot_uc(icontrol1) = Mz_tot_uc(icontrol1) + Mz


c	      Now for the parts due to U V W.  dThrustdV is from above

	      dThrustdU = dThrustdV * propeller(i)%nx
	      dThrustdW = dThrustdV * propeller(i)%nz
	      dThrustdV = dThrustdV * propeller(i)%ny

c	      Minus sign for the torque, as it is in the opposite direction
c	      of the propeller spin.

	      Mx_U = - dThrustdU * propeller(i)%spin*propeller(i)%nx
	      My_U = - dThrustdU * propeller(i)%spin*propeller(i)%ny
	      Mz_U = - dThrustdU * propeller(i)%spin*propeller(i)%nz

	      Mx_V = - dThrustdV * propeller(i)%spin*propeller(i)%nx
	      My_V = - dThrustdV * propeller(i)%spin*propeller(i)%ny
	      Mz_V = - dThrustdV * propeller(i)%spin*propeller(i)%nz

	      Mx_W = - dThrustdW * propeller(i)%spin*propeller(i)%nx
	      My_W = - dThrustdW * propeller(i)%spin*propeller(i)%ny
	      Mz_W = - dThrustdW * propeller(i)%spin*propeller(i)%nz

	      Mx_tot_U = Mx_tot_U + Mx_U
	      My_tot_U = My_tot_U + My_U
	      Mz_tot_U = Mz_tot_U + Mz_U

	      Mx_tot_V = Mx_tot_V + Mx_V
	      My_tot_V = My_tot_V + My_V
	      Mz_tot_V = Mz_tot_V + Mz_V

	      Mx_tot_W = Mx_tot_W + Mx_W
	      My_tot_W = My_tot_W + My_W
	      Mz_tot_W = Mz_tot_W + Mz_W

c	      Now we can distribute the thrust force

              Fx_U = dThrustdU * propeller(i)%nx
              Fy_U = dThrustdU * propeller(i)%ny
              Fz_U = dThrustdU * propeller(i)%nz

              Fx_V = dThrustdV * propeller(i)%nx
              Fy_V = dThrustdV * propeller(i)%ny
              Fz_V = dThrustdV * propeller(i)%nz

              Fx_W = dThrustdW * propeller(i)%nx
              Fy_W = dThrustdW * propeller(i)%ny
              Fz_W = dThrustdW * propeller(i)%nz

	      Fx_tot_U = Fx_tot_U + Fx_U
	      Fy_tot_U = Fy_tot_U + Fy_U
	      Fz_tot_U = Fz_tot_U + Fz_U

	      Fx_tot_V = Fx_tot_V + Fx_V
	      Fy_tot_V = Fy_tot_V + Fy_V
	      Fz_tot_V = Fz_tot_V + Fz_V

	      Fx_tot_W = Fx_tot_W + Fx_W
	      Fy_tot_W = Fy_tot_W + Fy_W
	      Fz_tot_W = Fz_tot_W + Fz_W

c	      Now the torque from the thrust (sometimes called the diferential thrust moment Mdt)

	      Mx_U  = y_w*Fz_U - z_w*Fy_U
	      My_U  = z_w*Fx_U - x_w*Fz_U
	      Mz_U  = x_w*Fy_U - y_w*Fx_U

	      Mx_V  = y_w*Fz_V - z_w*Fy_V
	      My_V  = z_w*Fx_V - x_w*Fz_V
	      Mz_V  = x_w*Fy_V - y_w*Fx_V

	      Mx_W  = y_w*Fz_W - z_w*Fy_W
	      My_W  = z_w*Fx_W - x_w*Fz_W
	      Mz_W  = x_w*Fy_W - y_w*Fx_W

	      Mx_tot_U = Mx_tot_U + Mx_U
	      My_tot_U = My_tot_U + My_U
	      Mz_tot_U = Mz_tot_U + Mz_U

	      Mx_tot_V = Mx_tot_V + Mx_V
	      My_tot_V = My_tot_V + My_V
	      Mz_tot_V = Mz_tot_V + Mz_V

	      Mx_tot_W = Mx_tot_W + Mx_W
	      My_tot_W = My_tot_W + My_W
	      Mz_tot_W = Mz_tot_W + Mz_W


c	      This is due to the gyroscopic behavior of the entire frame, where we 
c	      are subtracting vec (P Q R) x I motor_omega
c	      We are abusing Fi here, to avoid defining another variable.

	      temp = propeller(i)%Ir * propeller(i)%spin * motor_omega(i) ! propeller and motor angular momentum
	      Fx = temp * propeller(i)%nx
	      Fy = temp * propeller(i)%ny
	      Fz = temp * propeller(i)%nz

c	      Mx_P =   0.d0
	      My_P =   Fz
	      Mz_P = - Fy

	      Mx_Q = - Fz
c	      My_Q =   0.d0
	      Mz_Q =   Fx

	      Mx_R =   Fy
	      My_R = - Fx
c	      Mz_R =   0.d0

c	      Mx_tot_P = Mx_tot_P + Mx_P
	      My_tot_P = My_tot_P + My_P
	      Mz_tot_P = Mz_tot_P + Mz_P
	   
	      Mx_tot_Q = Mx_tot_Q + Mx_Q
c	      My_tot_Q = My_tot_Q + My_Q
	      Mz_tot_Q = Mz_tot_Q + Mz_Q
	   
	      Mx_tot_R = Mx_tot_R + Mx_R
	      My_tot_R = My_tot_R + My_R
c	      Mz_tot_R = Mz_tot_R + Mz_R

	   endif
	   
	enddo ! loop on i over the propellers


c	The acceleration from F = m(a + omega x U)

c	Gravity is in the z direction in the world frame, so that means
c	it will be in the direction given by the last column of the 
c	world frame to body frame transformation matrix.

	grav = aircraft%grav
	Udot = R*V - Q*W + grav*xo(13) + Fx_tot/aircraft%mass
	Vdot = P*W - R*U + grav*xo(16) + Fy_tot/aircraft%mass
	Wdot = Q*U - P*V + grav*xo(19) + Fz_tot/aircraft%mass

	if (ip) write (6,100) 'grav',0, 7
     .                         , grav*xo(13)*aircraft%mass
     .                         , grav*xo(16)*aircraft%mass
     .                         , grav*xo(19)*aircraft%mass


cc	write (6,*) ' dots ',Udot,Vdot,Wdot,aircraft%mass,grav
cc	write (6,*) ' Us   ',U   ,V   ,W   ,aircraft%mass,grav
cc	write (6,*) ' P, Q, R ',P, Q, R
cc	write (6,*) ' Force',Fx_tot,Fy_tot,Fz_tot
cc	write (6,*) ' Torq ',Mx_tot,My_tot,Mz_tot
c	write (6,*) ' xo( 3 6 9) ',xo(13),xo(16),xo(19)
c	temp = sqrt(Fx_tot**2+Fy_tot**2+Fz_tot**2)
c	write (6,*) ' F direction   ',Fx_tot/temp, Fy_tot/temp
c     .                               ,Fz_tot/temp
c	temp = sqrt(Up**2+Vp**2+Wp**2)
c	write (6,*) ' Up direction  ',Up/temp, Vp/temp
c     .                               ,Wp/temp
c	write (6,*) ' U  direction  ',U /temp, V /temp
c     .                               ,W /temp

c	Rotational acceleration and moments
c	See Eq. 1.7-11a and b page 39 for the description of engine effects.
c	These equation are repeated in Table 2.5-1 page 111.  Notice we are in a rotating
c	frame, because we are in the plane's body coordinates for P Q R.
c	The rotational equation is for vec(omega) = vec(P Q R),
c
c	vec(omega) dot = Iinv*(M - omega x I omega)  from M = dot(I omega) = I dot(omega) + omega x I omega
c
c	Multiplication done by REDUCE - this is omega x I omega
	Ixx = aircraft%Ixx
	Ixy = aircraft%Ixy
	Ixz = aircraft%Ixz
	Iyy = aircraft%Iyy
	Iyz = aircraft%Iyz
	Izz = aircraft%Izz
	ox  = -Ixy*P*R+Ixz*P*Q-Iyy*Q*R+Iyz*(Q*Q-R*R)+Izz*Q*R
	oy  =  Ixx*P*R+Ixy*Q*R+Ixz*(R*R-P*P)-Iyz*P*Q-Izz*P*R
	oz  = -Ixx*P*Q+Ixy*(P*P-Q*Q)-Ixz*Q*R+Iyy*P*Q+Iyz*P*R

c       This is now M - omega x I omega
	ox  = Mx_tot - ox
	oy  = My_tot - oy
	oz  = Mz_tot - oz

c	Now we hit it with the inverse to get vec omega dot
	InvIxx = aircraft%InvIxx
	InvIxy = aircraft%InvIxy
	InvIxz = aircraft%InvIxz
	InvIyy = aircraft%InvIyy
	InvIyz = aircraft%InvIyz
	InvIzz = aircraft%InvIzz
	Pdot   = InvIxx*ox + InvIxy*oy + InvIxz*oz
	Qdot   = InvIxy*ox + InvIyy*oy + InvIyz*oz
	Rdot   = InvIxz*ox + InvIyz*oy + InvIzz*oz

c	write (6,*) ' PQR dot ',Pdot, Qdot, Rdot

c	Navigation (Table 2.5-1, page 111), but with the quaternions.  
c	This computes (U V W) in the world frame and is entirely geometric (no physics).
c	We need to take the transpose.
	Unorth = xo(11)*U + xo(14)*V + xo(17)*W
	Veast  = xo(12)*U + xo(15)*V + xo(18)*W
	Wdown  = xo(13)*U + xo(16)*V + xo(19)*W


c	Quaternion update.
c	q0dot = 0.5d0*(     - P*q1 - Q*q2 - R*q3)
c	q1dot = 0.5d0*(P*q0        + R*q2 - Q*q3)
c	q2dot = 0.5d0*(Q*q0 - R*q1        + P*q3)
c	q3dot = 0.5d0*(R*q0 + Q*q1 - P*q2       )
c	The above is not adequate.  LOTS of drift.

c	We are seeing lots of drift.  Dreier suggests something a little different,
c	saying the ues of temp this way is "self correcting" (page 47)
c	He suggests multiplying by 32 (page 48).
	temp  = 32.d0*(1.d0 - q0**2 - q1**2 - q2**2 - q3**2 )
	q0dot = 0.5d0*(temp*q0 - P*q1    - Q*q2    - R*q3   )
	q1dot = 0.5d0*(P*q0    + temp*q1 + R*q2    - Q*q3   )
	q2dot = 0.5d0*(Q*q0    - R*q1    + temp*q2 + P*q3   )
	q3dot = 0.5d0*(R*q0    + Q*q1    - P*q2    + temp*q3)


c	Returning the computed derivatives of state variables.

	xd(1)  = Udot   !   body fixed frame acceleration
	xd(2)  = Vdot   !   body fixed frame acceleration
	xd(3)  = Wdot   !   body fixed frame acceleration

	xd(4)  = Pdot   !   body fixed frame rotational acceleration (roll rate)
	xd(5)  = Qdot   !   body fixed frame rotational acceleration (pitch rate)
	xd(6)  = Rdot   !   body fixed frame rotational acceleration (yaw rate)

	xd(7)  = q0dot  !   real quaternion (connecting body fixed frame to world frame) rate
	xd(8)  = q1dot  !   i quaternion rate
	xd(9)  = q2dot  !   j quaternion rate
	xd(10) = q3dot  !   k quaternion rate

	xd(11) = Unorth !   North speed (world frame)
	xd(12) = Veast  !   East  speed (world frame)
	xd(13) = Wdown  !   Down  speed (world frame); down is positive

	do i= 1, aircraft%num_propellers
           xd(13 + i) = motor_omegadot(i)
	enddo

c	Right now we have current stored in two places, but this is the place where 
c	it causally connects to the charge on the battery through time integration.
c	Still deciding whether I want it positive or negative.  Positive right now.
	do i = 1, aircraft%num_batteries
	   xd(13 + aircraft%num_propellers + i)
     .           = xo(29 + 2*aircraft%num_propellers + i)
	enddo

	xo(20) = Fx_lift
	xo(21) = Fy_lift
	xo(22) = Fz_lift

	xo(23) = Fx_drag
	xo(24) = Fy_drag
	xo(25) = Fz_drag
	
	xo(26) = Thrust_x
	xo(27) = Thrust_y
	xo(28) = Thrust_z

	xo(29) = Power_tot

	if (ip) write (6,100) 'Final Sum',0,8,
     .                Fx_tot + grav*xo(13)*aircraft%mass
     .              , Fy_tot + grav*xo(16)*aircraft%mass
     .              , Fz_tot + grav*xo(19)*aircraft%mass
     .               ,Mx_tot, My_tot, Mz_tot
 100    format(' Fx Fy Fz Mx My Mz ',a30,2i3,6f12.3)

	if (ip) write (6,49) (i,xo(i),
     .       i = 29,29+2*aircraft%num_propellers+aircraft%num_batteries)
 49	format(' xo ',6(i3,f10.3))

c	We now compute the linearized control entries

	if (compute_A) then

c	   State vector is

C	   (U V W P Q R q0 q1 q2 q3 Xnorth Yeast Zdown)

c	   The first three entries are in body frame, the quaternion links, last entry in world frame

	   Alin = 0.d0
	   Blin = 0.d0
	   Glin = 0.d0

c	   The following derivatives are needed - the complexity is due to needing to drop one
c	   of the quaternion equations, the specific one identified by qskip.  These direction
c	   cosine matrix terms are the way quaternions enter into A except for the quaternion
c	   equations themselves, hence the derivatives to compute A which is Alin = d xd / dx

c	   xo(11) = q0**2 + q1**2 - q2**2 - q3**2   ! xx
c	   xo(12) = 2.d0*(q1*q2 + q0*q3)            ! xy
c	   xo(13) = 2.d0*(q1*q3 - q0*q2)            ! xz
c	   xo(14) = 2.d0*(q1*q2 - q0*q3)            ! yx
c	   xo(15) = q0**2 - q1**2 + q2**2 - q3**2   ! yy
c	   xo(16) = 2.d0*(q2*q3 + q0*q1)            ! yz
c	   xo(17) = 2.d0*(q1*q3 + q0*q2)            ! zx
c	   xo(18) = 2.d0*(q2*q3 - q0*q1)            ! zy
c	   xo(19) = q0**2 - q1**2 - q2**2 + q3**2   ! zz

	   xo_d(11,0) =   2.d0*q0
	   xo_d(11,1) =   2.d0*q1
	   xo_d(11,2) = - 2.d0*q2
	   xo_d(11,3) = - 2.d0*q3

	   xo_d(12,0) =   2.d0*q3
	   xo_d(12,1) =   2.d0*q2
	   xo_d(12,2) =   2.d0*q1
	   xo_d(12,3) =   2.d0*q0

	   xo_d(13,0) = - 2.d0*q2
	   xo_d(13,1) =   2.d0*q3
	   xo_d(13,2) = - 2.d0*q0
	   xo_d(13,3) =   2.d0*q1

	   xo_d(14,0) = - 2.d0*q3
	   xo_d(14,1) =   2.d0*q2
	   xo_d(14,2) =   2.d0*q1
	   xo_d(14,3) = - 2.d0*q0

	   xo_d(15,0) =   2.d0*q0
	   xo_d(15,1) = - 2.d0*q1
	   xo_d(15,2) =   2.d0*q2
	   xo_d(15,3) = - 2.d0*q3

	   xo_d(16,0) =   2.d0*q1
	   xo_d(16,1) =   2.d0*q0
	   xo_d(16,2) =   2.d0*q3
	   xo_d(16,3) =   2.d0*q2

	   xo_d(17,0) =   2.d0*q2
	   xo_d(17,1) =   2.d0*q3
	   xo_d(17,2) =   2.d0*q0
	   xo_d(17,3) =   2.d0*q1

	   xo_d(18,0) = - 2.d0*q1
	   xo_d(18,1) = - 2.d0*q0
	   xo_d(18,2) =   2.d0*q3
	   xo_d(18,3) =   2.d0*q2

	   xo_d(19,0) =   2.d0*q0
	   xo_d(19,1) = - 2.d0*q1
	   xo_d(19,2) = - 2.d0*q2
	   xo_d(19,3) =   2.d0*q3

c	   Now we need to adjust the values here to reflect the fact that we are dropping one of the 
c	   quaternion equations, to get rid of the dependence.  By dropping q0, for example, all
c	   the derivatives of q0 now need to reflect the fact that they are dq0(q1,q2,q3)/dqi
c	   As an example, 
c	   d xo(11)/dq1 = partial xo(11)/partial q1 + partial xo(11)/partial q0 * partial q0/partial q1
c	                = partial xo(11)/partial q1 + partial xo(11)/partial q0 *( -q1/q0)

	   if (qskip.eq.0) then
	      dq0dq1 = -q1/q0
	      dq0dq2 = -q2/q0
	      dq0dq3 = -q3/q0
	      do i = 11, 19
	         xo_d(i,1) =   xo_d(i,1) + xo_d(i,0)*dq0dq1
	         xo_d(i,2) =   xo_d(i,2) + xo_d(i,0)*dq0dq2
	         xo_d(i,3) =   xo_d(i,3) + xo_d(i,0)*dq0dq3
	      enddo
	   else if (qskip.eq.1) then
	      dq1dq0 = -q0/q1
	      dq1dq2 = -q2/q1
	      dq1dq3 = -q3/q1
	      do i = 11, 19
	         xo_d(i,0) =   xo_d(i,0) + xo_d(i,1)*dq1dq0
	         xo_d(i,2) =   xo_d(i,2) + xo_d(i,1)*dq1dq2
	         xo_d(i,3) =   xo_d(i,3) + xo_d(i,1)*dq1dq3
	      enddo
	   else if (qskip.eq.2) then
	      dq2dq0 = -q0/q2
	      dq2dq1 = -q1/q2
	      dq2dq3 = -q3/q2
	      do i = 11, 19
	         xo_d(i,0) =   xo_d(i,0) + xo_d(i,2)*dq2dq0
	         xo_d(i,1) =   xo_d(i,1) + xo_d(i,2)*dq2dq1
	         xo_d(i,3) =   xo_d(i,3) + xo_d(i,2)*dq2dq3
	      enddo
	   else 
	      dq3dq0 = -q0/q3
	      dq3dq1 = -q1/q3
	      dq3dq2 = -q2/q3
	      do i = 11, 19
	         xo_d(i,0) =   xo_d(i,0) + xo_d(i,3)*dq3dq0
	         xo_d(i,1) =   xo_d(i,1) + xo_d(i,3)*dq3dq1
	         xo_d(i,2) =   xo_d(i,2) + xo_d(i,3)*dq3dq2
	      enddo
	   endif


c	   The first derivatives
c	   1 2 3 = U V W
c	   grav = aircraft%grav
c	   Udot = R*V - Q*W + grav*xo(13) + Fx_tot/aircraft%mass
c	   Vdot = P*W - R*U + grav*xo(16) + Fy_tot/aircraft%mass
c	   Wdot = Q*U - P*V + grav*xo(19) + Fz_tot/aircraft%mass

	   Alin( 1, 1) =       Fx_tot_U/aircraft%mass
	   Alin( 1, 2) =   R + Fx_tot_V/aircraft%mass
	   Alin( 1, 3) = - Q + Fx_tot_W/aircraft%mass
	   Alin( 1, 4) =       Fx_tot_P/aircraft%mass
	   Alin( 1, 5) = - W + Fx_tot_Q/aircraft%mass
	   Alin( 1, 6) =   V + Fx_tot_R/aircraft%mass
c           Alin( 1, 7) =       grav*xo_d(13,0)
c           Alin( 1, 8) =       grav*xo_d(13,1)
c           Alin( 1, 9) =       grav*xo_d(13,2)
c           Alin( 1,10) =       grav*xo_d(13,3)

	   Alin( 2, 1) = - R + Fy_tot_U/aircraft%mass
	   Alin( 2, 2) =       Fy_tot_V/aircraft%mass
	   Alin( 2, 3) =   P + Fy_tot_W/aircraft%mass
	   Alin( 2, 4) =   W + Fy_tot_P/aircraft%mass
	   Alin( 2, 5) =       Fy_tot_Q/aircraft%mass
	   Alin( 2, 6) = - U + Fy_tot_R/aircraft%mass
c           Alin( 2, 7) =       grav*xo_d(16,0)
c           Alin( 2, 8) =       grav*xo_d(16,1)
c           Alin( 2, 9) =       grav*xo_d(16,2)
c           Alin( 2,10) =       grav*xo_d(16,3)

	   Alin( 3, 1) =   Q + Fz_tot_U/aircraft%mass
	   Alin( 3, 2) = - P + Fz_tot_V/aircraft%mass
	   Alin( 3, 3) =       Fz_tot_W/aircraft%mass
	   Alin( 3, 4) = - V + Fz_tot_P/aircraft%mass
	   Alin( 3, 5) =   U + Fz_tot_Q/aircraft%mass
	   Alin( 3, 6) =       Fz_tot_R/aircraft%mass
c           Alin( 3, 7) =       grav*xo_d(19,0)
c           Alin( 3, 8) =       grav*xo_d(19,1)
c           Alin( 3, 9) =       grav*xo_d(19,2)
c           Alin( 3,10) =       grav*xo_d(19,3)

	   do i = 1, 3
	      if (qskip.eq.0) then
c	         0 means that the q0dot equation was dropped (7)
	         l = i
	      else if (qskip.eq.1) then
	         if (i.eq.1) l = 0
	         if (i.eq.2) l = 2
	         if (i.eq.3) l = 3
	      else if (qskip.eq.2) then
	         if (i.eq.1) l = 0
	         if (i.eq.2) l = 1
	         if (i.eq.3) l = 3
	      else
	         l = i - 1
	      endif
	      k = 6 + i
              Alin( 1, k) =    grav*xo_d(13,l)
              Alin( 2, k) =    grav*xo_d(16,l)
              Alin( 3, k) =    grav*xo_d(19,l)
	   enddo

c
c	   vec(omega) dot = Iinv*(M - omega x I omega)  from M = dot(I omega) = I dot(omega) + omega x I omega
c
c	   Multiplication done by REDUCE - this is omega x I omega
c	   Ixx = aircraft%Ixx
c	   Ixy = aircraft%Ixy
c	   Ixz = aircraft%Ixz
c	   Iyy = aircraft%Iyy
c	   Iyz = aircraft%Iyz
c	   Izz = aircraft%Izz
c	   ox  = -Ixy*P*R+Ixz*P*Q-Iyy*Q*R+Iyz*(Q*Q-R*R)+Izz*Q*R
c	   oy  =  Ixx*P*R+Ixy*Q*R+Ixz*(R*R-P*P)-Iyz*P*Q-Izz*P*R
c	   oz  = -Ixx*P*Q+Ixy*(P*P-Q*Q)-Ixz*Q*R+Iyy*P*Q+Iyz*P*R
c
c          This is now M - omega x I omega
c	   ox  = Mx_tot - ox
c	   oy  = My_tot - oy
c	   oz  = Mz_tot - oz
c
c	   Now we hit it with the inverse to get vec omega dot
c	   InvIxx = aircraft%InvIxx
c	   InvIxy = aircraft%InvIxy
c	   InvIxz = aircraft%InvIxz
c	   InvIyy = aircraft%InvIyy
c	   InvIyz = aircraft%InvIyz
c	   InvIzz = aircraft%InvIzz
c	   Pdot   = InvIxx*ox + InvIxy*oy + InvIxz*oz
c	   Qdot   = InvIxy*ox + InvIyy*oy + InvIyz*oz
c	   Rdot   = InvIxz*ox + InvIyz*oy + InvIzz*oz

	   ox_P  = -Ixy*R+Ixz*Q
	   oy_P  =  Ixx*R+Ixz*(-2.d0*P)-Iyz*Q-Izz*R
	   oz_P  = -Ixx*Q+Ixy*( 2.d0*P)+Iyy*Q+Iyz*R

	   ox_Q  =  Ixz*P-Iyy*R+Iyz*( 2.d0*Q)+Izz*R
	   oy_Q  =  Ixy*R-Iyz*P
	   oz_Q  = -Ixx*P+Ixy*(-2.d0*Q)-Ixz*R+Iyy*P

	   ox_R  = -Ixy*P-Iyy*Q+Iyz*(-2.d0*R)+Izz*Q
	   oy_R  =  Ixx*P+Ixy*Q+Ixz*( 2.d0*R)-Izz*P
	   oz_R  = -Ixz*Q+Iyz*P

	   ox_P  = Mx_tot_P - ox_P
	   oy_P  = My_tot_P - oy_P
	   oz_P  = Mz_tot_P - oz_P

	   ox_Q  = Mx_tot_Q - ox_Q
	   oy_Q  = My_tot_Q - oy_Q
	   oz_Q  = Mz_tot_Q - oz_Q

	   ox_R  = Mx_tot_R - ox_R
	   oy_R  = My_tot_R - oy_R
	   oz_R  = Mz_tot_R - oz_R


	   Alin(4,4)   = InvIxx*ox_P + InvIxy*oy_P + InvIxz*oz_P
	   Alin(4,5)   = InvIxx*ox_Q + InvIxy*oy_Q + InvIxz*oz_Q
	   Alin(4,6)   = InvIxx*ox_R + InvIxy*oy_R + InvIxz*oz_R

	   Alin(5,4)   = InvIxy*ox_P + InvIyy*oy_P + InvIyz*oz_P
	   Alin(5,5)   = InvIxy*ox_Q + InvIyy*oy_Q + InvIyz*oz_Q
	   Alin(5,6)   = InvIxy*ox_R + InvIyy*oy_R + InvIyz*oz_R

	   Alin(6,4)   = InvIxz*ox_P + InvIyz*oy_P + InvIzz*oz_P
	   Alin(6,5)   = InvIxz*ox_Q + InvIyz*oy_Q + InvIzz*oz_Q
	   Alin(6,6)   = InvIxz*ox_R + InvIyz*oy_R + InvIzz*oz_R


c	   The quaternions (xd 7 to 10) in terms of PQR (x 4 to 6) and q (x 7 to 10)
c	   Quaternion update.  This is complicated in that we are both removing
c	   a row and a column based on which one is excluded (qskip).
c	   q0dot = 0.5d0*(     - P*q1 - Q*q2 - R*q3)
c	   q1dot = 0.5d0*(P*q0        + R*q2 - Q*q3)
c	   q2dot = 0.5d0*(Q*q0 - R*q1        + P*q3)
c	   q3dot = 0.5d0*(R*q0 + Q*q1 - P*q2       )
cold	   Alin ( 7, 4) = 0.5d0*( - q1)
c	   Alin ( 7, 5) = 0.5d0*( - q2)
c	   Alin ( 7, 6) = 0.5d0*( - q3)
c	   Alin ( 8, 4) = 0.5d0*(   q0)
c	   Alin ( 8, 5) = 0.5d0*( - q3)
c	   Alin ( 8, 6) = 0.5d0*(   q2)
c	   Alin ( 9, 4) = 0.5d0*(   q3)
c	   Alin ( 9, 5) = 0.5d0*(   q0)
c	   Alin ( 9, 6) = 0.5d0*( - q1)	
c	   Alin (10, 4) = 0.5d0*( - q2)
c	   Alin (10, 5) = 0.5d0*(   q1)
c	   Alin (10, 6) = 0.5d0*(   q0)

	   Alin ( 7, 4) = 0.5d0*( - q1)
	   Alin ( 7, 5) = 0.5d0*( - q2)
	   Alin ( 7, 6) = 0.5d0*( - q3)

	   i = 8
	   if (qskip.eq.0) i = 7  ! to overwrite because q0 equation is discarded 7 8 9 = 8 9 10

	   Alin ( i, 4) = 0.5d0*(   q0)
	   Alin ( i, 5) = 0.5d0*( - q3)
	   Alin ( i, 6) = 0.5d0*(   q2)

	   i = i + 1
	   if (qskip.eq.1) i = i - 1  ! to overwrite because q1 equation is discarded 7 8 9 = 7 9 10
   
     	   Alin ( i, 4) = 0.5d0*(   q3)
	   Alin ( i, 5) = 0.5d0*(   q0)
	   Alin ( i, 6) = 0.5d0*( - q1)	

	   i = i + 1
	   if (qskip.eq.2) i = i - 1  ! to overwrite because q2 equation is discarded 7 8 9 = 7 8 10

	   Alin ( i, 4) = 0.5d0*( - q2)
	   Alin ( i, 5) = 0.5d0*(   q1)
	   Alin ( i, 6) = 0.5d0*(   q0)

c	   The last one should be qskip 3 because 7 8 9 = 7 8 9 and we skip writing 10


c	   Now it is more complicated because we are both skipping a column and doing 
c	   derivatives with respect to the qi.
c	   q0dot = 0.5d0*(     - P*q1 - Q*q2 - R*q3)
c	   q1dot = 0.5d0*(P*q0        + R*q2 - Q*q3)
c	   q2dot = 0.5d0*(Q*q0 - R*q1        + P*q3)
c	   q3dot = 0.5d0*(R*q0 + Q*q1 - P*q2       )
cold	   Alin ( 7, 7) = 0.d0
c	   Alin ( 7, 8) = 0.5d0*(  - P)
c	   Alin ( 7, 9) = 0.5d0*(  - Q)
c	   Alin ( 7,10) = 0.5d0*(  - R)
c	   Alin ( 8, 7) = 0.5d0*(    P)
c	   Alin ( 8, 8) = 0.d0
c	   Alin ( 8, 9) = 0.5d0*(    R)
c	   Alin ( 8,10) = 0.5d0*(  - Q)
c	   Alin ( 9, 7) = 0.5d0*(    Q)
c	   Alin ( 9, 8) = 0.5d0*(  - R)
c	   Alin ( 9, 9) = 0.d0
c	   Alin ( 9,10) = 0.5d0*(    P)
c	   Alin (10, 7) = 0.5d0*(    R)
c	   Alin (10, 8) = 0.5d0*(    Q)
c	   Alin (10, 9) = 0.5d0*(  - P)
c	   Alin (10,10) = 0.d0

	   if (qskip.eq.0) then
c	      q1dot = 0.5d0*(P*q0        + R*q2 - Q*q3)
c	      q2dot = 0.5d0*(Q*q0 - R*q1        + P*q3)
c	      q3dot = 0.5d0*(R*q0 + Q*q1 - P*q2       )
	      dq0dq1 = -q1/q0
	      dq0dq2 = -q2/q0
	      dq0dq3 = -q3/q0
	      Alin ( 7, 7) = 0.5d0*(       P * dq0dq1 )
	      Alin ( 7, 8) = 0.5d0*(   R + P * dq0dq2 )
	      Alin ( 7, 9) = 0.5d0*( - Q + P * dq0dq3 )

	      Alin ( 8, 7) = 0.5d0*( - R + Q * dq0dq1 )
	      Alin ( 8, 8) = 0.5d0*(       Q * dq0dq2 )
	      Alin ( 8, 9) = 0.5d0*(   P + Q * dq0dq3 )

	      Alin ( 9, 7) = 0.5d0*(   Q + R * dq0dq1 )
	      Alin ( 9, 8) = 0.5d0*( - P + R * dq0dq2 )
	      Alin ( 9, 9) = 0.5d0*(       R * dq0dq3 )

	   else if (qskip.eq.1) then
c	      q0dot = 0.5d0*(     - P*q1 - Q*q2 - R*q3)
c	      q2dot = 0.5d0*(Q*q0 - R*q1        + P*q3)
c	      q3dot = 0.5d0*(R*q0 + Q*q1 - P*q2       )
	      dq1dq0 = -q0/q1
	      dq1dq2 = -q2/q1
	      dq1dq3 = -q3/q1
	      Alin ( 7, 7) = 0.5d0*(     - P * dq1dq0 )
	      Alin ( 7, 8) = 0.5d0*( - Q - P * dq1dq2 )
	      Alin ( 7, 9) = 0.5d0*( - R - P * dq1dq3 )

	      Alin ( 8, 7) = 0.5d0*(   Q - R * dq1dq0 )
	      Alin ( 8, 8) = 0.5d0*(     - R * dq1dq2 )
	      Alin ( 8, 9) = 0.5d0*(   P - R * dq1dq3 )

	      Alin ( 9, 7) = 0.5d0*(   R + Q * dq1dq0 )
	      Alin ( 9, 8) = 0.5d0*( - P + Q * dq1dq2 )
	      Alin ( 9, 9) = 0.5d0*(       Q * dq1dq3 )

	   else if (qskip.eq.2) then
c	      q0dot = 0.5d0*(     - P*q1 - Q*q2 - R*q3)
c	      q1dot = 0.5d0*(P*q0        + R*q2 - Q*q3)
c	      q3dot = 0.5d0*(R*q0 + Q*q1 - P*q2       )
	      dq2dq0 = -q0/q2
	      dq2dq1 = -q1/q2
	      dq2dq3 = -q3/q2
	      Alin ( 7, 7) = 0.5d0*(     - Q * dq2dq0 )
	      Alin ( 7, 8) = 0.5d0*( - P - Q * dq2dq1 )
	      Alin ( 7, 9) = 0.5d0*( - R - Q * dq2dq3 )

	      Alin ( 8, 7) = 0.5d0*(   P + R * dq2dq0 )
	      Alin ( 8, 8) = 0.5d0*(       R * dq2dq1 )
	      Alin ( 8, 9) = 0.5d0*( - Q + R * dq2dq3 )

	      Alin ( 9, 7) = 0.5d0*(   R - P * dq2dq0 )
	      Alin ( 9, 8) = 0.5d0*(   Q - P * dq2dq1 )
	      Alin ( 9, 9) = 0.5d0*(     - P * dq2dq3 )

	   else 
c	      q0dot = 0.5d0*(     - P*q1 - Q*q2 - R*q3)
c	      q1dot = 0.5d0*(P*q0        + R*q2 - Q*q3)
c	      q2dot = 0.5d0*(Q*q0 - R*q1        + P*q3)
	      dq3dq0 = -q0/q3
	      dq3dq1 = -q1/q3
	      dq3dq2 = -q2/q3
	      Alin ( 7, 7) = 0.5d0*(     - R * dq3dq0 )
	      Alin ( 7, 8) = 0.5d0*( - P - R * dq3dq1 )
	      Alin ( 7, 9) = 0.5d0*( - Q - R * dq3dq2 )

	      Alin ( 8, 7) = 0.5d0*(   P - Q * dq3dq0 )
	      Alin ( 8, 8) = 0.5d0*(     - Q * dq3dq1 )
	      Alin ( 8, 9) = 0.5d0*(   R - Q * dq3dq2 )

	      Alin ( 9, 7) = 0.5d0*(   Q + P * dq3dq0 )
	      Alin ( 9, 8) = 0.5d0*( - R + P * dq3dq1 )
	      Alin ( 9, 9) = 0.5d0*(       P * dq3dq2 )

	   endif

c	   All other Al (7 8 9 10,*) are zero.  Bl (7 8 9 10,*) are zero.


c	   The world frame position (xd 11 to 13) in terms of UVW (x 1 to 3) and q (x 7 to 10)
c	   Since we dropped one of the quaternion equations, the index here runs from 10 to 12
c	   rather than 11 to 13.
c	   Unorth = xo(11)*U + xo(14)*V + xo(17)*W
c	   Veast  = xo(12)*U + xo(15)*V + xo(18)*W
c	   Wdown  = xo(13)*U + xo(16)*V + xo(19)*W
           Alin (10, 1) = xo(11)
	   Alin (10, 2) = xo(14)
	   Alin (10, 3) = xo(17)
	   Alin (11, 1) = xo(12)
	   Alin (11, 2) = xo(15)
	   Alin (11, 3) = xo(18)
	   Alin (12, 1) = xo(13)
	   Alin (12, 2) = xo(16)
	   Alin (12, 3) = xo(19)

c	   Unorth = xo(11)*U + xo(14)*V + xo(17)*W
c	   Veast  = xo(12)*U + xo(15)*V + xo(18)*W
c	   Wdown  = xo(13)*U + xo(16)*V + xo(19)*W

	   do i = 1, 3
	      if (qskip.eq.0) then
c	         0 means that the q0dot equation was dropped (7)
	         l = i
	      else if (qskip.eq.1) then
	         if (i.eq.1) l = 0
	         if (i.eq.2) l = 2
	         if (i.eq.3) l = 3
	      else if (qskip.eq.2) then
	         if (i.eq.1) l = 0
	         if (i.eq.2) l = 1
	         if (i.eq.3) l = 3
	      else
	         l = i - 1
	      endif
	      k = 6 + i
	      Alin (10, k) = xo_d(11,l)*U + xo_d(14,l)*V + xo_d(17,l)*W
	      Alin (11, k) = xo_d(12,l)*U + xo_d(15,l)*V + xo_d(18,l)*W
	      Alin (12, k) = xo_d(13,l)*U + xo_d(16,l)*V + xo_d(19,l)*W
	   enddo

c	   All other Alin(10 11 12, *) are zero.  Blin (10 11 12,*) are zero.

c	   Notice that Alin(*, 10 11 12) = 0.


c	   Now we build Blin.

c	   State vector is

C	   (U V W P Q R q0 q1 q2 q3 Xnorth Yeast Zdown)

           do i = 1, udim

c	      The acceleration from F = m(a + omega x U)
	      Blin(1,i) = Fx_tot_uc(i)/aircraft%mass
	      Blin(2,i) = Fy_tot_uc(i)/aircraft%mass
	      Blin(3,i) = Fz_tot_uc(i)/aircraft%mass

c	      The rotational acceleration term is from 
c	      vec(omega) dot = Iinv*(M - omega x I omega)  from M = dot(I omega) = I dot(omega) + omega x I omega
c	      so for B all that matters is Iinv*M
	      ox        = Mx_tot_uc(i)
	      oy        = My_tot_uc(i)
	      oz        = Mz_tot_uc(i)
	      Blin(4,i) = InvIxx*ox + InvIxy*oy + InvIxz*oz
	      Blin(5,i) = InvIxy*ox + InvIyy*oy + InvIyz*oz
	      Blin(6,i) = InvIxz*ox + InvIyz*oy + InvIzz*oz

c	      And the rest are zero.

	   enddo


c	   For now we are not going to do this.  We will assume that
c	   We only need to consider the 13 x 13 situation.
c	   This is an approximation, but probably a good one.
c	   do i = 1, aircraft%num_propellers
c	      From above we have, yielding Alin and Glin
c	      motor_omegadot(i) = (omega - motor_omega(i))
c     .                       /propeller(i)%tau_ESC
c	      Alin (13+i,13+i) = - 1.d0/propeller(i)%tau_ESC
c	      omega = motor_omegadot(i)*propeller(i)%tau_ESC+motor_omega(i)
c	      Glin (13+i) = omega/propeller(i)%tau_ESC
c	   enddo
c	   We assume all other Alin are zero (but this is not strictly true).
c	   The Blin are not zero here, as the requested voltage affects
c	   the change in motor rate.

	   if (ip) then
	      write (6,*) ' Alin = '
	      do i = 1, 12
	         write (6,1000) (Alin(i,k),k=1,12)
	      enddo
	      write (6,*) ' Blin = '
	      do i= 1, 12
	         write (6,1000) (Blin(i,k),k=1,udim)
	      enddo
	   endif
 1000	   format(12e12.5)

	   control%Alin = Alin
	   control%Blin = Blin
c	   control%Glin = Glin  ! right now this is not coded


	endif ! to comput_A


	return
	end


	subroutine CLDwing(
c			   Input
     .			   alpha,a,C_L0,C_Lmin,C_Lmax,C_D0,k,C_Dfp
c			   Output
     .                    ,C_L,C_D,iflinear,dC_Ldalpha,dC_Ddalpha)

c	This computes the wing lift and drag coefficient for the full 360 degrees.
c	Though for some reason, Stevens et al. does not seem to include any effect
c	of sideslip, so we may need to modify this in that context.

	implicit none

	include 'pi.h'

c	Input
	double precision alpha  ! angle of attack
	double precision a      ! slope of lift curve
	double precision C_L0   ! lift coefficient at zero angle of attack
	double precision C_Lmax ! maximum lift coefficient (before separation)
	double precision C_Lmin ! minimum lift coefficient
	double precision C_D0   ! zero lift drag coefficient of wing
	double precision k      ! induced drag coefficient of wing
	double precision C_Dfp  ! drag at 90 angle of attack; seems this always be 1.

c	Output
	double precision C_L ! lift coefficient
	double precision C_D ! drag coefficient
	logical iflinear     ! this is true if in the linear region, well defined region.
c		               False if not.  The point of this flag is to know if we are 
c			       in the well defined region, which may have bearing in a table search.
	double precision dC_Ldalpha    ! derivative of C_L with respect to alpha (for controls)
	double precision dC_Ddalpha    ! derivative of C_D with respect to alpha (for controls)

	double precision C_Dsep ! separation drag
	double precision alpha_L_min, alpha_L_max
	double precision dC_Dsepdalpha


c	We know the angle of attack is at least -pi.
c	These interpolations are on pg 650.
	alpha_L_min = (C_Lmin - C_L0)/a
	alpha_L_max = (C_Lmax - C_L0)/a
	iflinear  = .false.
	if (alpha.lt.-halfpi) then
           C_L    = C_Dfp*sin(2.d0*alpha)
	   C_Dsep = C_Dfp*0.5d0*(1.d0 - cos(2.d0*alpha))
	   dC_Ldalpha    = 2.d0*C_Dfp*cos(2.d0*alpha)
	   dC_Dsepdalpha =      C_Dfp*sin(2.d0*alpha)
	else if (alpha.lt.alpha_L_min) then
           C_L    = (halfpi+alpha)*C_Lmin/(halfpi+alpha_L_min)
	   C_Dsep = C_Dfp*(alpha-alpha_L_min)/(-halfpi-alpha_L_min)
	   dC_Ldalpha    = C_Lmin/( halfpi+alpha_L_min)
	   dC_Dsepdalpha = C_Dfp /(-halfpi-alpha_L_min)
	else if (alpha.lt.alpha_L_max) then
	   C_L    = a*alpha+C_L0
	   C_Dsep = 0.d0
	   dC_Ldalpha    = a
	   dC_Dsepdalpha = 0.d0
	   iflinear = .true.
	else if (alpha.lt.halfpi) then
           C_L    = (halfpi-alpha)*C_Lmax/(halfpi-alpha_L_max)
	   C_Dsep = C_Dfp*(alpha-alpha_L_max)/(halfpi-alpha_L_max)
	   dC_Ldalpha    = -C_Lmax/(halfpi-alpha_L_max)
	   dC_Dsepdalpha =  C_Dfp /(halfpi-alpha_L_max)
	else
           C_L    = C_Dfp*sin(2.d0*alpha)
	   C_Dsep = C_Dfp*0.5d0*(1.d0 - cos(2.d0*alpha))
	   dC_Ldalpha    = 2.d0*C_Dfp*cos(2.d0*alpha)
	   dC_Dsepdalpha =      C_Dfp*sin(2.d0*alpha)
	endif

	C_D = C_D0 + k*C_L**2 + C_Dsep
	dC_Ddalpha = 2.d0*k*C_L*dC_Ldalpha + dC_Dsepdalpha

c	write (6,*) ' C_L, C_D ',C_L,C_D
c	write(6,*)' a,C_Lmin,CL0,CLmax,alpha',a,C_Lmin,C_L0,C_Lmax,alpha

	return
	end


	subroutine table (n, x, x_table, dxkm1, dxk, k, kpast, ierror)
c	This is a table lookup

c	Input:
c	n          number of entries in column
c	x          x
c	x_table    x entries in table
c	           It is assumed the table is ordered in terms of increasing x and
c	           there are no repeat elements in x.
c
c       Input, Output:
c
c	Output:
c	The output allows you to calculate how ever many other colums of data you may have
c	with y = dxkm1*y_table(k-1) + dxk*y_table(k)
c	dxkm1      The k entry weighting coefficient
c	dxk        The k entry weighting coefficient
c	k          index location for interpolation
c	ierror     = 0 if no error, = 1 if below 0, = 2 if extrapolation above last point
c	           Our extrapolation is simply constant of last entry.
c
	implicit none

	double precision x, x_table(n), dxkm1, dxk
	integer          n, ierror, i, k
	integer       :: kpast    ! last call value of k
	logical       :: isuccess

c	Some of the ifs here are because the value of kpast could conceiably be from a 
c	different table, and so not related to this table.

	if ((kpast.lt.2).or.(kpast.gt.n)) goto 10  ! either first time or not from this table.

	k        = kpast
	isuccess = .false.

c	This is our quick search to see if we are near the needed value on the table 
c	based on the previous call
	if (x.ge.x_table(k-1)) then
           if (x.le.x_table(k)) then
c	      We are in the same cell as the previous call.
	      dxk      = (x - x_table(k-1))/(x_table(k) - x_table(k-1))
	      isuccess = .true.
	      ierror   = 0
	   else if ((k.eq.n).and.(x.ge.x_table(n))) then
c	      We are beyond the upper part of the table
              dxk      = 1.d0
              isuccess = .true.
	      ierror   = 2
	   else if (x.le.x_table(k+1)) then
c	      We are in the next cell up, and k < n
	      k        = k + 1
	      dxk      = (x - x_table(k-1))/(x_table(k) - x_table(k-1))
	      isuccess = .true.
	      ierror   = 0
           endif
	else if (k.eq.2) then
c          We are below the bottom of the table.
	   dxk      = 0.d0
	   isuccess = .true.
	   ierror   = 1
	else if (x.ge.x_table(k-2)) then
c	   We are in the next cell down and k > 2
	   k = k - 1
	   dxk      = (x - x_table(k-1))/(x_table(k) - x_table(k-1))
	   isuccess = .true.
	   ierror   = 0
	endif

	if (isuccess) goto 20

 10	continue

c	This is the regular, walk up the table read.

	if (x.le.x_table(1)) then
c	   Extrapolate below
	   k      = 2
	   dxk    = 0.d0
           ierror = 1
        else if (x.ge.x_table(n)) then
c	   Extrapolate above
	   k      = n
           dxk    = 1.d0
	   ierror = 2
	else
c	   We are somewhere in the table
	   k = 2
	   do while (x.gt.x_table(k))
	      k = k + 1
	   enddo
	   dxk      = (x - x_table(k-1))/(x_table(k) - x_table(k-1))
	   ierror   = 0
	endif

 20	continue

	dxkm1 = 1.d0 - dxk
	kpast = k

	return
	end

