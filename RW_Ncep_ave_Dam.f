      program NCEP
!------------this program reads the NCEP 2.5 degree or Gauss Grid data download from
!                      http://www.esrl.noaa.gov/psd/data/
!--- run in the WinXP system and the netcdf library is needed for the fortran compiler
!--- you can change the data start date and end date in line 76 as needed
!-- NOTE: when more than year data is handled, the momery maybe is unsunfficiency.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!--one time processing more than one year data is not tested!!!!!  
!
!----------------------------------------------------------------------------------
!------------ 3-D variables are u v omega RH Temp. and Hgt with 17 levels--------
!------the surface variables include latent heat/sensible heat(guass grid)
!-------------------surface u v omege RH and pressure(sigma 0.995)
!-------------------------------last modified 26/3/2013  Jinghua Chen
!---------the rainfall is pwr------------------------------------------------------- 
	include 'netcdf.inc'
	integer,parameter :: nrec=1460
	integer,parameter :: nrecr=1464
	integer,parameter :: nnr=nrecr
	integer,parameter :: lev=17 !! add surface
	integer,parameter :: n2d=14 !1lh 2sh 3pwr 4us 5vs 6 ts 7omegas 8Rhs 9Ps 10P_top, 11cpr 12DLR 13ULR 14 t2m
	integer,parameter :: n3d=6 !1u 1v 3t  4oemga 5RH 6HGT
!	real dat2d(15,10,nrecr) !(lon lat time ivar)  SH LH
!      real dat2d(15,10,17,nrecr) !(lon lat level time ivar) u,v,t w,,RH,hgt
      common/DD/ data3d(6,6,lev,nnr,n3d), !1u 1v 3t  4oemga 5RH 6HGT
     +      data2d(10,10,nnr,n2d),hgt(10,10) ! !1lh 2sh 3pwr 4us 5vs 6 ts 7omegas 8Rhs 9Ps 10P_top 11cpr 12DLR 13ULR,14t2m
	COMMON/TQ/ tvs(10,10,nnr),qvs(10,10,nnr),
     +	tv(10,10,lev,nnr),qv(10,10,lev,nnr)
	COMMON/Q12/ Q1(10,10,lev,nnr),Q2(10,10,lev,nnr),
     +	tls(10,10,lev,nnr),qls(10,10,lev,nnr)
	COMMON/DXY/ DX(10),DY
	COMMON/AVEG/HAD_Q(10,10,lev,nnr),VAD_Q(10,10,lev,nnr),
     +	TCH_Q(10,10,lev,nnr),HADT(10,10,lev,nnr),VADT(10,10,lev,nnr),
     +	TCHT(10,10,lev,nnr) 
	character(len=150) :: dir,filename,folds,file,dirout
	character(len=30) :: var3d(n3d),var2d(n2d),area(4),foldout,date
	character(len=4) :: years,yeare,year
	character(len=2) months,monthe,dys,dye
	integer IXX(4),IYY(4),IXG(4),IYG(4),nnrec
	integer X1(4),Y1(4),X2(4),Y2(4), XG1(4),YG1(4),XG2(4),YG2(4) 
	real plv(lev),dp(lev), bslatG(4),bslat(4),bslon(4)
	integer iys,iye,ims,ime,ids,ide,days(12)
!	integer ktop(10,10) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----variables with the suffix of '_xy' are average of the grid box and for output 
	real tls_xy(lev,nnr),prsf,pwr
	real qls_xy(lev,nnr),cpwr_xy(nnr)
	real u_xy(lev,nnr),the_xy(lev,nnr)
	real v_xy(lev,nnr),rh_xy(lev,nnr)
	real gz_xy(lev,nnr),Q1_xy(lev,nnr)
	real t_xy(lev,nnr),Q2_xy(lev,nnr)
	real qv_xy(lev,nnr),omega_xy(lev,nnr)
	real lh_xy(nnr),pwr_xy(nnr),psrf_xy(nnr),t2m_xy(nnr)
	real sh_xy(nnr) ! ,omega(10,10,lev,nnr)
	real DIVO(10,10,lev,nnr),OMGO(10,10,lev,nnr)
     *          ,DIV(10,10,lev,nnr),OMG(10,10,lev,nnr)
      real XY_OUT(nnr,20),XYG_OUT(nnr,10),XYN_OUT(nnr,lev,20)
	real omegas(10,10,nnr) !WSFC(10,10,nnr),
	real qlss(10,10,nnr),tlss(10,10,nnr)
	real Q1temps(10,10,nnr),Q2temps(10,10,nnr)
	real TS_XY(nnr),US_XY(nnr),VS_XY(nnr),GZS_XY,QVS_XY(nnr)
      real omegaS_XY(nnr),Q1S_XY(nnr),Q2S_XY(nnr),dyn(lev,nnr,4) ! 1-3 for vorticity 4 for divergence
	real TLSS_XY(nnr),QLSS_XY(nnr), rhs_xy(nnr),thes_xy(nnr)
!
      DIMENSION plw(lev)
!
	data rd,cp,g,rv,hlat /287.,1005.,9.81,461.,2.5e6/
	data plv/1000.,925.,850.,700.,
     +     600., 500.,400.,300.,250.,200.,150.,100.,70.,50.,30.,20.,10./
      data plw/962.5, 887.5, 775., 650., 550., 450., 350., 275., 225., 
     +      175., 125., 85., 60., 40., 25., 15., 10./   !plw(K)=(plv(k)+plv(k+1))/2 plw(17)=plv(16)/2
      data dp/75.,75.,150.,100.,100.,100.,100.,100.,50.,
     +        50.,50.,50.,30.,20.,20.,10.,10./   ! top pressure is zero donw-up
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--------------------read the start date and the end date ---------------------------
!	print*,'Please input the start year'
!	read*, iys
!      print*,'Please input the start month'
!	read*,ims
!      print*,'Please input the start day'
!	read*, ids
!	print*,'Please input the ending year'
!	read*, iye
!      print*,'Please input the ending month'
!	read*, ime
!      print*,'Please input the ending day'
!	read*, ide
!-------------------------------------------------------------------------------------
!------------- set your date that is needed here ----------------------------------------
!-- iys=start_year iye=end_yaer ims=start_month ime=end_month ids=start_day ide=end_day
!      iys=2000   ;   iye=2000
	ims=1      ;   ime=12           !
	ids=1      ;   ide=31
	do 1300 ichen=2000,2012
       iys=ichen
	 iye=ichen
	print*,iys
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	stop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----------------------check the date-----------------------------------------------
       if(iye<iys.or.iys<1979.or.iye>2012)then
	  print*,'the ending year is before the start year 
     +	  or the data has not downlaoded'
       stop
	 elseif(iye==iys.and.ime<ims)then
	  print*,'the ending month is before the start month'
       stop
       elseif(iye==iys.and.ime==ims.and.ide<ids)then
	 print*,'the ending day is before the start day'
       stop
	 endif
	!---------------------------------------------------------------------
       do i=1,12
	  days(i)=31
	 enddo
!	hour(1)='00:00'
!	hour(2)='06:00'
!	hour(3)='12:00'
!	hour(4)='18:00'
	 days(2)=28
	 days(6)=30
	 days(4)=30
	 days(9)=30
	 days(11)=30
	write(years,'(I4.4)')iys
	write(yeare,'(I4.4)')iye
	write(months,'(I2.2)')ims
	write(monthe,'(I2.2)')ime
	write(dys,'(I2.2)')ids
	write(dye,'(I2.2)')ide
	PIE=3.141592657
!!!!-------------------------set ----------------------------------------------------
! variables that named ending with character 'G' stand for the Gauss Grid  
! the following variables are for determining different regions 
! area stores the region's name, X1(:) Y1(:)store the start  gird ID,
!  X2(:) Y2(:) is end grid ID, IXX(:) IYY(:)  the number grid of different regions. 
! bslat bslon the star lat. and lon. of different regions 
! target area 3 grid X 3 grid, but when Q1 and Q2 are calculated, the loop is:
!    2 to IXX-1 and 2 to IYY-1, so 5 grid X 5 grid is needed
!-------PRD=preal River Deltr   
 	 area(1)='PRD'; X1(1)=44; X2(1)=48; Y1(1)=30; Y2(1)=26;IXX(1)=5
       IYY(1)=5; bslat(1)=17.5; XG1(1)=58; XG2(1)=64; YG1(1)=38 
    	 YG2(1)=33
       IYG(1)=6; bslatG(1)=18.095;IXG(1)=7;bslon(1)=107.5
!     MLYR=Middle & Low reaches of Yangzi River  
	area(2)='MLYR'; X1(2)=44; X2(2)=48 ;Y1(2)=27; Y2(2)=23; IXX(2)=5
	 IYY(2)=5; bslat(2)=25;XG1(2)=58;XG2(2)=64;YG1(2)=34; YG2(2)=29
       IYG(2)=6; bslatG(2)=25.7139;IXG(2)=7;bslon(2)=107.5
!     NPC=North Plain of China
	area(3)='NPC'; X1(3)=45; X2(3)=49;Y1(3)=24; Y2(3)=20; IXX(3)=5
	IYY(3)=5; bslat(3)=32.5;XG1(3)=60;XG2(3)=65;YG1(3)=30;YG2(3)=25
      IYG(3)=6; bslatG(3)=33.3328;IXG(3)=6;bslon(3)=110
!     NEC  NorthEast of China
	area(4)='NEC'; X1(4)=48; X2(4)=52; Y1(4)=21; Y2(4)=17; IXX(4)=5
      IYY(4)=5; bslat(4)=40.;XG1(4)=64;XG2(4)=69;YG1(4)=26; YG2(4)=21
       IYG(4)=6; bslatG(4)=40.9517;IXG(4)=6;bslon(4)=117.5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      var3d(1)='uwnd.'
	var3d(2)='vwnd.'
	var3d(3)='air.'
	var3d(4)='omega.'
	var3d(5)='rhum.'
	var3d(6)='hgt.'
!1lh 2sh 3pwr 4us 5vs 6 t2m 7omegas 8Rhs 
	var2d(1)='lhtfl.sfc.gauss.'
	var2d(2)='shtfl.sfc.gauss.'
      var2d(3)='prate.sfc.gauss.'
	var2d(4)='uwnd.sig995.'
	var2d(5)='vwnd.sig995.'
	var2d(6)='air.sig995.'
	var2d(7)='omega.sig995.'
	var2d(8)='rhum.sig995.'
	var2d(9)='pres.sfc.'
	var2d(10)='pres.tropp.'
	var2d(11)='cprat.sfc.gauss.'
	var2d(12)='dlwrf.sfc.gauss.'
	var2d(13)='ulwrf.sfc.gauss.'
	var2d(14)='air.2m.gauss.'

!---------------- make the output fold
	dir='Z:\DATA\LargeScale\'
      foldout=years(3:4)//months//dys//'-'//yeare(3:4)//monthe//dye
	dirout=trim(dir)//'NcepR2_Pre'
      istatus1=CHDIR(trim(dirout))
	istatus2=system("Md "//trim(foldout))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      do 1200 lo=1,4 ! start the area loop
	print*,lo,'start processing for area of ',area(lo)
	print*,'please waitting'
	   itt=0
!--------------read the surface hgt ------------------------------------
           folds='NCEP_Surface\'
           filename=trim(dir)//trim(folds)//'hgt.sfc.nc'
	 call NC_HGT(lo,hgt,filename,IXX,IYY,
     +              X1,Y1,X2,Y2)
!--------------------------------------------------------------------------------------
        do iyear=iys,iye
	      nnrec=1460
	      days(2)=28
            write(year,'(I4.4)')iyear
	      iyid=1 ! normal year nrec=1460
	     if(mod(iyear,4)==0.and.mod(iyear,100)/=0)then
	     days(2)=29 ! nrec=1464
	     nnrec=nrecr
	     elseif(mod(iyear,400)==0)then
	     days(2)=29
	     nnrec=nrecr
	     endif
	     
!	     if(iyy==0) nnrec=nrecr
           if(iyear>iys) itt=itt+nnrec ! first year don't plus
!!!!!!!!!!!!!!!!!!!!!! read and write 2D files first---------------------------------
          do i2d=1,n2d
!	    print*,i2d
	     folds='NCEP_Surface\'
	     file=trim(var2d(i2d))//year//'.nc'
	     filename=trim(dir)//trim(folds)//trim(file)
	     if(i2d<4.or.i2d>=11)then  ! 11 for convection rainfall rate
!	read the heat flux in Gauss grid
            if(i2d==14)then
               call t2msst(lo,nnrec,itt,i2d,data2d,filename,IXG,IYG,
     +              XG1,YG1,XG2,YG2,nnr,n2d)
	     else
	        call NC_2D(lo,nnrec,itt,i2d,data2d,filename,IXG,IYG,
     +              XG1,YG1,XG2,YG2,nnr,n2d)
	       endif
!
	     elseif(i2d>3.and.i2d<11)then
!---------------read surface variables in normal grid
           call NC_SRF(lo,nnrec,itt,i2d,data2d,filename,IXX,IYY,
     +              X1,Y1,X2,Y2,nnr,n2d)
!	print*,itt,'NC_3D'
	     endif
	    enddo !!!end read i2d
!!!!!!!!!!!!!!!!!!!!! End 2D files ---------------------------------------------------
!
!!---------------- read and write 3D files ---------------------------------
          do i3d=1,n3d
!	print*,i3d
	     folds='NCEPR2\'
	     file=trim(var3d(i3d))//year//'.nc'
!	      print*,data3d(2,3,1,1,1)
	     filename=trim(dir)//trim(folds)//trim(file)
!	       print*,filename
!          print*,itt,'before NC_3D'
	     call NC_3D(lo,nnrec,itt,i3d,data3d,filename,IXX,IYY,
     +              X1,Y1,X2,Y2,nnr,n3d,lev) !X1,Y1,X2,Y2)
!	      print*,data3d(2,3,1,1,1)
!       print*,itt,'NC_3D'
!	stop
	    enddo !!!end read i3d
!	stop
!       print*,'test',iyear
!!!!!!!!!!!!!!!!!!!!! End 3D files ---------------------------------------------------
        enddo !!! ending reading data
!---------------------------------------------------------------------------------------
!   the 3D variables stored in DATA3d(6,6,lev,nnr,n3d) n3d/name 1/u 2/v 3/oemga 4/t 5/RH 6/HGT
!   the 2D variables stored in DATA2d(6,6,lev,nnr,n2d) n2d/name 1/lh 2/sh 3/pwr 4/us 5/vs 6/ts 7/omegas 8/Rhs 9/Ps
!---------------------------------------------------------------------------------------
!----- the NCEP every signle file stores one year data, so the TimeID of the period we need should be checked. 
!!!!!------calculate the start TimeID---------------
!      print*,itt,'end readding'
      itm=0
	itte=0
	days(2)=28
	    if(mod(iys,4)==0.and.mod(iys,100)/=0)then
	     days(2)=29 ! nrec=1464
	     elseif(mod(iys,400)==0)then
	     days(2)=29
	     endif
	  do imm=1,ims-1
	  itm=itm+days(imm)*4
	  enddo
	 itd=(ids-1)*4
       ittb=itm+itd+1
!!!!!------calculate the ending TimeID---------------
!----------------Part 1-------------------------------
      itm=0
	itte1=0
	days(2)=28
	    if(mod(iye,4)==0.and.mod(iye,100)/=0)then
	     days(2)=29 
	     elseif(mod(iye,400)==0)then
	     days(2)=29
	     endif
	  do imm=1,ime-1
	  itm=itm+days(imm)*4
	  enddo
	 itd=ide*4
       itte1=itm+itd
!----------------Part 2-------------------------------
	 itm=0
	 itte2=0
        if(iye>iys)then
		days(2)=28
	    if(mod(iys,4)==0.and.mod(iys,100)/=0)then
	     days(2)=29 
	     elseif(mod(iys,400)==0)then
	     days(2)=29
	     endif
	     do imm=ims+1,12
	    itm=itm+days(imm)*4
	     enddo
		 itte2=itm+(days(ims)-ids+1)*4
	  endif
!----------------Part 1-------------------------------
	  itte3=0
	  if((iye-iys)>=2)then
	   do jy=iys+1,iye-1
	      njy=nrec
		 if(mod(jy,4)==0.and.mod(jy,100)/=0)then
	     njy=nrecr 
	     elseif(mod(jy,400)==0)then
	     njy=nrecr
	     endif
	    itte3=itte3+njy
	   enddo
	 endif
	itte=itte1+itte2+itte3     	 	 	        
!      print*,ittb,itte,'itt',itt    !----------------test code 
!------------open output files---------------------------------------
!-----------------surface heat flux and rainfall water
	filename=trim(dirout)//'/'//trim(foldout)//'/'//trim(area(lo))//
     +'_HF.txt'
	open(10,file=trim(filename))
	write(10,905)'From',iys,ims,ids,'to',iye,ime,ide,
     +	'Every six hour one data'
      write(10,*)'TimeID latentHeat(W/m^2) SensibleHeat(W/m^2) 
     + Precipipitable_water_for_entire_atmosphere(kg/m^2)'
!-----------------sounding ------------------------------
	filename=trim(dirout)//'/'//trim(foldout)//'/'//trim(area(lo))//
     +'_sounding.txt'
	open(11,file=trim(filename))
	write(11,905)'From',iys,ims,ids,'to',iye,ime,ide,
     +	'Every six hour one data'
      write(11,*)'Press(hPa) heigh(m) U(m/s) V(m/s) 
     + omega(pa/s) Temp(K) Theta(K) Qv(kg/kg) RH(%)'
	
!-----------large scale forcing -------------------
      filename=trim(dirout)//'/'//trim(foldout)//'/'//trim(area(lo))//
     +'_Forcing.txt'
	open(12,file=trim(filename))
	write(12,905)'From',iys,ims,ids,'to',iye,ime,ide,
     +	'Every six hour one data'
      write(12,*)'TimeID 17_levels_T_forcing(K/day)
     + 17_levels_qv_forcing(g/kg/day)'
!---------------Q1 and Q2---------------------------------------------		
	filename=trim(dirout)//'/'//trim(foldout)//'/'//trim(area(lo))//
     +'_q1q2.txt'
	open(13,file=trim(filename))
	write(13,905)'From',iys,ims,ids,'to',iye,ime,ide,
     +	'Every six hour one data'
      write(13,*)'TimeID Q1(K/day)  Q2(K/day)'
!---------------------------------------------------------------------
      filename=trim(dirout)//'/'//trim(foldout)//'/'//trim(area(lo))//
     +'_surface22.txt'
	open(1014,file=trim(filename))
	write(1014,905)'From',iys,ims,ids,'to',iye,ime,ide,
     +	'Every six hour one data'
      write(1014,'(100A)')'Press(hPa) heigh(m) U(m/s) V(m/s) 
     + omega(pa/s) Temp(K) Theta(K) Qv(kg/kg) RH(%) t2m(K)'
!----------------------------------------------------------------------
!      print*,ittb,nnrec+itt,'1099'
!	stop
!-------------------------------------------------------------------------------
!!------------start calculate the forcing and Q1 Q2 for output-----------------------
	do 1099 it=ittb,itte
!	print*,itt,'star loop 1099'
!-------------------Heat Flux-----------------------------------------------------
	     do ix=2,IXG(lo)-1
	       do iy=2,IYG(lo)-1
              lhio=lhio+DATa2d(ix,iy,it,1)  !1lh 2sh 3pwr 4us 5vs 6t2m 7omegas 8Rhs 9Ps
	        shio=shio+DATa2d(ix,iy,it,2)
			pwr=pwr+DATa2d(ix,iy,it,3)
			cpwr=cpwr+DATa2d(ix,iy,it,11)	              
	         enddo
	        enddo
	      do ix=2,IXX(lo)-1
	       do iy=2,IYY(lo)-1
			psrf=psrf+DATa2d(ix,iy,it,9)
			tmpgz=tmpgz+hgt(ix,iy)
			tmpts=tmpts+DATa2d(ix,iy,it,6)
			tmpus=tmpus+DATa2d(ix,iy,it,4)
	        tmpvs=tmpvs+DATa2d(ix,iy,it,5)
	        tmpos=tmpos+DATa2d(ix,iy,it,7)
              tmprhs=tmprhs+DATa2d(ix,iy,it,8)
			t2ms=t2ms+data2d(ix,iy,it,14)			              
	         enddo
	        enddo
           XYG=(IXG(lo)-2.0)*(IYG(lo)-2.0) ! Gauss grid
	     lh_xy(it)=lhio/XYG
	     sh_xy(it)=shio/XYG
	     XYG=(IXX(lo)-2.0)*(IYY(lo)-2.0) ! normal grid
	     pwr_xy(it)=pwr/XYG
	     cpwr_xy(it)=cpwr/XYG
	     psrf_xy(it)=psrf/XYG
	     gzs_xy=tmpgz/XYG
	     ts_xy(it)=tmpts/XYG
	     us_xy(it)=tmpus/XYG
	     vs_xy(it)=tmpvs/XYG
	     omegas_xy(it)=tmpos/XYG
	     rhs_xy(it)=tmprhs/XYG
	     t2m_xy(it)=t2ms/XYG
	     lhio=0.
	     shio=0.
	     pwr=0.
	     cpwr=0.0  ; t2ms=0.
	     psrf=0.0; tmprhs=0.
	     tmpgz=0.; tmpts=0. ; tmpus=0. ;  tmpvs=0. ; tmpos=0.
	     write(10,900)it-ittb+1,lh_xy(it),sh_xy(it),
     +		 pwr_xy(it),cpwr_xy(it)! write the surface variables into files. 
!-----------------------------------------------------------------------------------
!--------- the following three-do-loops calcultate the potential temperature, qv and omega------------
	      do 1100 ix=1,IXX(lo)
	        do 1100 iy=1,IYY(lo)
!--------calcultate surface the potential temperature, qv
!--------------------surface-------------!1lh 2sh 3pwr 4us 5vs 6t2m 7omegas 8Rhs 9Ps
	         tvs(ix,iy, it)=
     +	      DATA2D(ix,iy,it,6)*(1000*1e2/DATA2D(ix,iy,it,9))**0.286
	         den=DATA2D(ix,iy,it,9)/(rd*DATA2D(ix,iy,it,6))  ! air density 
              esat=611.*exp(hlat/rv *
     +			   (1./273.16-1./DATA2D(ix,iy,it,6)))
	        qvs00=esat/(rv*den*DATA2D(ix,iy,it,6))
			qvs(ix,iy, it)=                                ! surface qv
     +			 DATA2D(ix,iy,it,8)/100.*qvs00 
!--------------------------------------------------------------------------------
!--------------- calcultate 3D potential temperature, qv --------------------------------------------
                ilev1=1
	          do  ik=1,lev
	        tv(ix,iy, ik, it)=         !potential temperature
     +			DATA3D(ix,iy,ik,it,3)*(1000/plv(ik))**0.286
	        den=plv(ik)*1.e2/(rd*DATA3D(ix,iy,ik,it,3))
              esat=611.*exp(hlat/rv *
     +			   (1./273.16-1./DATA3D(ix,iy,ik,it,3)))
	        qvs00=esat/(rv*den*DATA3D(ix,iy,ik,it,3))
!	       if(ik>8)  DATA3D(ix,iy,ik,it,5)=0  !about 300hpa RH=0  Yanai 1998
			qv(ix,iy, ik, it)=                     !qv in kg/kg
     +			 DATA3D(ix,iy,ik,it,5)/100.*qvs00 
!

	       enddo

1100  continue   ! potential temperature and qv
!     stop   enf of calculateing the potential temperature, qv and omega for all the required gird of the need period.
!       print*,itt,'end of loop 1099'
1099  continue



!_____________________________________________________________
!----------------recalculate the omega fowllow the suggests of Yanai(1998,J.C) and Luo ()
!---------at surface  omega=-g*aird*((u/(Re*cos(lat)))*(dh/dlon) + (v/Re)*(dh/dlat))
!-------------aird=air ddensity
!---------------Re= earth radius
!----------------h=height from the NCEP DATA
!-----for normal the omega is calculated by  1/(R*cos(lat)) * (du/dlat+ d(vcoslat)/dlat)  +domega/dp =0   F_omega

!--------------------------------------------------------------------------
!----------------recalculate the omega fowllow the suggests of Yanai(1998,J.C)
!---------top omega=-(dtheta/dt+H theta)/(dtheta/dp)---------------------------
!---------surface omega=-g*aird*((u/(Re*cos(lat)))*(dh/dlon) + (v/Re)*(dh/dlat))
!----------------------------------------------------------------------------------------------------------	 
	     
!--------- the following loops calculate the forcings, Q1 and Q2
!      print*,itt,'start loop 1098'
	filename=trim(dirout)//'/'//trim(foldout)//'/'//years//
     +trim(area(lo))//'_lowlevel.txt'
	  open(14,file=trim(filename))
	filename=trim(dirout)//'/'//trim(foldout)//'/'//years//
     +trim(area(lo))//'_Middlelevel.txt'
	  open(15,file=trim(filename))
	filename=trim(dirout)//'/'//trim(foldout)//'/'//years//
     +trim(area(lo))//'_Hightlevel.txt'
	  open(16,file=trim(filename))
	filename=trim(dirout)//'/'//trim(foldout)//'/'//years//
     +trim(area(lo))//'_Abovelevel.txt'
	  open(17,file=trim(filename))
	filename=trim(dirout)//'/'//trim(foldout)//'/'//years//
     +trim(area(lo))//'_Tropsphere.txt'
	  open(18,file=trim(filename))
		filename=trim(dirout)//'/'//trim(foldout)//'/'//years//
     +trim(area(lo))//'_AllLevels.txt'
	  open(19,file=trim(filename))
	 filename=trim(dirout)//'/'//trim(foldout)//'/'//years//
     +trim(area(lo))//'_RAW.txt'
	  open(20,file=trim(filename))
      do jj=14,20
	write(jj,908)'DATE','HOUR','T_ls(k/day)','Q_ls(k/day)','U(m/s)',
     +'V(m/s)','moisture(kg/kg)','HGT(m)','AIR(K)','Adj_omega(pa/s)',
     +'RH(%)','Theta(K)','Q1(k/day)','Q2(k/day)','HADQ(K/day)',
     +'VADQ(K/day)','TCHQ(K/day)','HADT(K/day)','HADT(K/day)',
     +'TCHT(K/day)','Ori_omega(pa/s)'
      enddo
	filename=trim(dirout)//'/'//trim(foldout)//'/'//years//
     +trim(area(lo))//'_SURFACE.txt'
	  open(21,file=trim(filename))
      write(21,909)'DATE','HOUR','LH(W/m^2)','SH(W/m^2)',
     +'prate','cprate','DLR(W/m^2)','ULR(W/m^2)'
!
	filename=trim(dirout)//'/'//trim(foldout)//'/'//years//
     +trim(area(lo))//'_dyn.txt'
	  open(22,file=trim(filename))
      write(22,909)'DATE','HOUR','VorX','VorY',
     +'VorZ','DIV'
!--------------------------------------------------------------

	do 1098 it=ittb,itte
      write(11,*)it-ittb+1
	write(13,*)it-ittb+1

	 DY=2.5*111.17E3
       DO J=1,IYY(lo)
       Y=((IYY(lo)-J)*2.5+bslat(lo))/180.
       DX(J)=COS(Y*PIE)*2.5*111.17E3
       ENDDO
	DIVO=0.0
	OMGO=0.0
	DIV=0.0
	OMG=0.0
	omegas=0.0
	ikk=0
	if(lo==4) ikk=1
	 call DIVVOR(PLV,PLW,omegas,DIVO,IXX(lo),IYY(lo),bslat(lo),IT,ikk) !DIVVOR(PM,WSFC,DIVO,NXDIR,NYDIR,YMIN,IT)
       call VERVEL(PLV,PLW,omegas,DIVO,OMGO,IXX(lo),IYY(lo),IT,ikk)
	 call ADJDV(PLV,PLW,DIVO,OMGO,DIV,IXX(lo),IYY(lo),IT,ikk)
	 call VERVEL(PLV,PLW,omegas,DIV,OMG,IXX(lo),IYY(lo),IT,ikk)  ! OMG mb/s
	 call Q1Q2(PLV,PLW,OMG,IXX(lo),IYY(lo),IT,lo,ikk)   !OMG in mb/s out pa/s  
	 call dynamic(data3d,OMG,DIV,IXX(lo),IYY(lo),bslat(lo),IT,dyn) !IXX(lo),IYY(lo),bslat(lo),IT)
!----------------------------------------------------------------------------------------
       Tmcl=24*3600
       do ik=1,lev
	 do 1097 ix=2,IXX(lo)-1
	       do 1097 iy=2,IYY(lo)-1
	          Tlsio=Tlsio+tls(ix,iy, ik, it) *3600*24  !day-1
	          qlsio=Qlsio+qls(ix,iy, ik, it) *3600*24*hlat/cp !!! K/day day-1
	          uio=uio+DATA3D(ix,iy,ik,it,1)
	          vio=vio+DATA3D(ix,iy,ik,it,2)
	          qvio=qvio+qv(ix,iy, ik, it)
	          gzio=gzio+DATA3D(ix,iy,ik,it,6)
	          tio=tio+DATA3D(ix,iy,ik,it,3)
	          omegaio=omegaio+omg(ix,iy, ik, it)
	          rhio=rhio+DATA3D(ix,iy,ik,it,5)
!	          qvio=qvio+qv(ix,iy, ik, it)
	          theio=theio+tv(ix,iy, ik, it)
	          Q1io=Q1io+Q1(ix, iy, ik,it) *Tmcl
	          Q2io=Q2io+Q2(ix, iy, ik,it) *TMCL*hlat/cp     !k/day

	          tmqvs=tmqvs+qvs(ix,iy,it)
	          tmtvs=tmtvs+tvs(ix,iy,it)
	          

1097  continue
      XY=(IXX(lo)-2.0)*(IYY(lo)-2.)
      tls_xy(ik,it)=Tlsio/XY
	qls_xy(ik,it)=qlsio/XY
	u_xy(ik,it)=uio/XY
	v_xy(ik,it)=vio/XY
	gz_xy(ik,it)=gzio/XY
	t_xy(ik,it)=tio/XY
	qv_xy(ik,it)=qvio/XY
	omega_xy(ik,it)=omegaio/XY
	Q1_xy(ik,it)=Q1io/XY
	Q2_xy(ik,it)=Q2io/XY
	rh_xy(ik,it)=rhio/XY
	the_xy(ik,it)=theio/XY

	Q1s_xy(it)=Q1ios/XY
	Q2s_xy(it)=Q2ios/XY
      tlss_xy(it)=Tlsios/XY
	qlss_xy(it)=Tlsios/XY
	qvs_xy(it)=tmqvs/XY
	thes_xy(it)=tmtvs/XY
	
	
	Tlsio=0.
	qlsio=0.
	uio=0.
	vio=0.
	qvio=0.
	gzio=0.
	tio=0.0
	omegaio=0.0
	Q1io=0.0
	Q2io=0.0
	Q1ios=0.0
	Q2ios=0.0
	Tlsios=0.
	qlsios=0.
      tmqvs=0.
	rhio=0.
	qvio=0.
	theio=0.
	tmtvs=0.
	t2ms=0.
	write(11,901)plv(ik),gz_xy(ik,it),u_xy(ik,it),v_xy(ik,it),
     +omega_xy(ik,it),t_xy(ik,it),the_xy(ik,it),qv_xy(ik,it),
     +rh_xy(ik,it)
      write(13,903)plv(ik),Q1_xy(ik,it),Q2_XY(ik,it)
	enddo
      write(12,902)it-ittb+1,tlss_xy(it),(tls_xy(ikk,it),ikk=1,lev),
     +	qlss_xy(it),(qls_xy(ikk,it),ikk=1,lev)
!!!
	write(1014,901)psrf_XY(it)/100.0,gzs_xy,us_xy(it),vs_xy(it),
     +omegas_xy(it),ts_xy(it),thes_xy(it),qvs_xy(it),
     +rhs_xy(it),t2m_xy(it)
!------------------------------------------------------------------------
!	ims
!	ids

!	if
!-------------low atomosphere --------------------
        call calendar(iys,iye,ims,ime,ids,ide,it,date)
        iks=1;ike=4    ! 1000 to 700
        call NCEP_SAS(PM,PW,OMG,IXX(lo),IYY(lo),IT,iks,ike,XY_OUT)
	  write(14,906)trim(date),(XY_OUT(it,kk),kk=1,19)
!---------middle atomosphere
        iks=5;ike=8 ! 600 to 300
        call NCEP_SAS(PM,PW,OMG,IXX(lo),IYY(lo),IT,iks,ike,XY_OUT)
	  write(15,906)trim(date),(XY_OUT(it,kk),kk=1,19)
!          high atomosphere(ice cloud)
        iks=9;ike=12   ! 300 to 100
        call NCEP_SAS(PM,PW,OMG,IXX(lo),IYY(lo),IT,iks,ike,XY_OUT)
	  write(16,906)trim(date),(XY_OUT(it,kk),kk=1,19)
!----------above troposphere
        iks=13;ike=17 ! 70 to 10
        call NCEP_SAS(PM,PW,OMG,IXX(lo),IYY(lo),IT,iks,ike,XY_OUT)
	  write(17,906)trim(date),(XY_OUT(it,kk),kk=1,19)
!------------all troposphere
        iks=1;ike=12   ! 1000 to 100
        call NCEP_SAS(PM,PW,OMG,IXX(lo),IYY(lo),IT,iks,ike,XY_OUT)
	  write(18,906)trim(date),(XY_OUT(it,kk),kk=1,19)
!----------- all levels
        iks=1;ike=17  ! 1000 to 10
        call NCEP_SAS(PM,PW,OMG,IXX(lo),IYY(lo),IT,iks,ike,XY_OUT)
	  call NCEP_GUS(PM,PW,IXG(lo),IYG(lo),IT,iks,ike,XYG_OUT)
	  write(19,906)trim(date),(XY_OUT(it,kk),kk=1,19)
        write(21,907)trim(date),(XYG_OUT(it,kk),kk=1,6)
!------------ no levels ave
        call NCEP_NOSAS(PM,PW,OMG,IXX(lo),IYY(lo),IT,XYN_OUT)
	  do ikl=1,lev
	  write(20,906)trim(date),(XYN_OUT(it,ikl,kk),kk=1,19)
	  enddo
!-------------- dynamic output------------------------------------
	  do ikl=1,lev
	  write(22,910)trim(date),(dyn(ikl,it,kk),kk=1,4)
	  enddo
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
1098  continue      
      print*,'end processing for area of ',area(lo)
1200  continue

900   format(1X,I4,1X,F8.2,1X,F8.2,1X,e12.4,1X,e12.4)
901   format(1X,4(1X,F9.3),1X,e12.4,2(1X,F9.2),1X,e12.4,1X,F7.3,1X,F7.3)
902   format(1X,I4,36(1X,e12.4))
903   format(1X,F9.3,2(1X,e12.4))
905   format(1X,A4,1X,I4,2I2.2,1X,A2,1X,I4,2I2.2,1X,A14)
906   format(1X,A16,1X,19(1X,e12.4))
907   format(1X,A16,1X,6(1X,e12.4))
908   format(1X,A10,1X,A5,19(1X,A11))
909   format(1X,A10,1X,A5,6(1X,A11))
910   format(1X,A16,1X,4(1X,e12.4))
      do i=10,13
	close(i)
	enddo
      nt=1-ittb+itte
      Call Grads(years,foldout,nt)
	Call Grads_PRE(years,foldout,nt)
1300  continue
	stop


	end  ! program NCEP


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!
!------------------  the following code copy from hmbudps.f-----------------
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC  DIVVOR COMPUTES DERIVED DIVERGENCE AND VORTICITY
CC  Note the dx and dy are vector, when the HAD is handle with, carefully note the data store squence
CC  and the array  subscript.
CC  NCEP: EAST-WEAT grids, the longitude array(x+1,y) is great than array(x,y)
CC        South North: the latitude array(x,y+1) is smaller than array(x,y)
      SUBROUTINE DIVVOR(PM,PW,WSFC,DIVO,NXDIR,NYDIR,YMIN,IT,ikk)
	PARAMETER (GRID=2.50)
      integer,parameter :: nrec=1460
	integer,parameter :: nrecr=1464
	integer,parameter :: nnr=nrecr
	integer,parameter :: lev=17 !! add surface
	integer,parameter :: n2d=13 !1lh 2sh 3pwr 4us 5vs 6 ts 7omegas 8Rhs 9Ps 10P_top
	integer,parameter :: n3d=6 !!1u 1v 3t  4oemga 5RH 6HGT
	integer NXDIR,NYDIR,it,ktop
      DIMENSION WSFC(10,10,nnr),DIVO(10,10,lev,nnr),PM(lev),PW(lev)
      common/DD/ data3d(6,6,lev,nnr,n3d), !1u 1v 3t 4oemga 5RH 6HGT
     +      data2d(10,10,nnr,n2d),hgt(10,10) ! !1lh 2sh 3pwr 4us 5vs 6 ts 7omegas 8Rhs 9Ps 10P_top
      DATA PIE/0.0174532925/
CC
      DY=1./(GRID*111.17E3)
      DO 100 I=1,NXDIR
      DO 100 J=1,NYDIR
      Y=((NYDIR-J)*GRID+YMIN)/180.
      DX=1./(COS(Y*PIE)*GRID*111.17E3)
CC ******** COMPUTE TERRAIN INDUCED OMEGAS FOR LOWEST LAYER
       IF(I.EQ.1 .OR. I.EQ.NXDIR) GO TO 95
      WX=data2d(i,j,it,4)*DX/2.*(hgt(I+1,J)-hgt(I-1,J))
   95  IF(I.EQ.1) WX=data2d(i,j,it,4)*DX*(hgt(I+1,J)-hgt(I,J))
       IF(I.EQ.NXDIR) WX=data2d(i,j,it,4)*DX*(hgt(I,J)-hgt(I-1,J))
       IF(J.EQ.1 .OR. J.EQ.NYDIR) GO TO 96
      WY=data2d(i,j,it,5)*DY/2.*(hgt(I,J-1)-hgt(I,J+1))
   96  IF(J.EQ.1) WY=data2d(i,j,it,5)*DY*(hgt(I,J)-hgt(I,J+1))
       IF(J.EQ.NYDIR) WY=data2d(i,j,it,5)*DY*(hgt(I,J-1)-hgt(I,J))
      WSFC(I,J,it)=-data2d(i,j,it,9)/100*9.8/
     +  (287.*data2d(i,j,it,6))*(WX+WY)
CC UNIT: MB/S FOR WSFC
      p_top=data2d(I,J,it,10)/100.
	   pres=data2d(I,J,it,9)/100.
	   ktop=10  ! 100hPa
         do K=1,lev-1
		 if(pm(k)>p_top.and.pw(k)<P_top)then
		    ktop=k-1-ikk
		   goto 101
		 elseif(pw(k)>P_top.and.pm(k+1)<p_top)then
            ktop=k-ikk
             goto 101
		 endif
	   enddo
  101 CONTINUE
CC ******** COMPUTE DIVERGENCE AND VORTICITY
       DO 200 K=1,ktop
       IF(K.NE.1) GO TO 400
C%%%%%%%%%%%  SURFACE DIVERGENCE
CC BOUNDARY CALCULATIONS
      IF(I.EQ.1) THEN
!        DUX=(U(I+1,J,K)-U(I,J,K))*DX  
		DUX=(data3d(I+1,J,K,it,1)-data3d(I,J,K,it,1))*Dx
      ELSEIF(I.EQ.NXDIR) THEN
!        DUX=(U(I,J,K)-U(I-1,J,K))*DX
		DUX=(data3d(I,J,K,it,1)-data3d(I-1,J,K,it,1))*Dx
      ELSE
CC CENTER FINITE DIFFERENCE CALCULATIONS
!       DUX=(U(I+1,J,K)-U(I-1,J,K))*DX/2.
		DUX=(data3d(I+1,J,K,it,1)-data3d(I-1,J,K,it,1))*Dx/2.
      ENDIF
CC BOUNDARY CALCULATIONS
      IF(J.EQ.1) THEN
!        DVY=(V(I,J+1,K)-V(I,J,K))*DY
		DUY=-(data3d(I,J+1,K,it,2)-data3d(I,J,K,it,2))*DY
      ELSEIF(J.EQ.NYDIR) THEN
!        DVY=(V(I,J,K)-V(I,J-1,K))*DY
		DUY=-(data3d(I,J,K,it,2)-data3d(I,J-1,K,it,2))*DY
      ELSE
CC CENTER FINITE DIFFERENCE CALCULATIONS
!        DVY=(V(I,J+1,K)-V(I,J-1,K))*DY/2.
		DUY=-(data3d(I,J+1,K,it,2)-data3d(I,J-1,K,it,2))*DY/2.
      ENDIF
      GO TO 500
C%%%%%%%%%%%   UPPER LAYER DIVERGENCE
CC BOUNDARY CALCULATIONS
  400 IF(I.EQ.1 .OR. (data2d(I-1,J,IT,9)/100).LT.PM(K)) THEN
!        DUX=(U(I+1,J,K)-U(I,J,K))*DX
		DUX=(data3d(I+1,J,K,it,1)-data3d(I,J,K,it,1))*DX
      ELSEIF(I.EQ.NXDIR .OR. (data2d(I+1,J,IT,9)/100).LT.PM(K)) THEN
!        DUX=(U(I,J,K)-U(I-1,J,K))*DX
		DUX=(data3d(I,J,K,it,1)-data3d(I-1,J,K,it,1))*DX
      ELSE
CC CENTER FINITE DIFFERENCE CALCULATIONS
!        DUX=(U(I+1,J,K)-U(I-1,J,K))*DX/2.
		DUX=(data3d(I+1,J,K,it,1)-data3d(I-1,J,K,it,1))*DX/2.
      ENDIF
CC BOUNDARY CALCULATIONS
      IF(J.EQ.1 .OR. (data2d(I,J-1,IT,9))/100.LT.PM(K)) THEN
!        DVY=(V(I,J+1,K)-V(I,J,K))*DY
		DUY=-(data3d(I,J+1,K,it,2)-data3d(I,J,K,it,2))*DY
      ELSEIF(J.EQ.NYDIR .OR. data2d(I,J+1,IT,9).LT.PM(K)) THEN
!        DVY=(V(I,J,K)-V(I,J-1,K))*DY
		DUY=-(data3d(I,J,K,it,2)-data3d(I,J-1,K,it,2))*DY
      ELSE
CC CENTER FINITE DIFFERENCE CALCULATIONS
        DVY=-(data3d(I,J+1,K,it,2)-data3d(I,J-1,K,it,2))*DY/2.
!		DUY=(data3d(I,J+1,K,it,2)-data3d(I,J-1,K,it,2))*DY/2.
      ENDIF
C
C  COMPUTE DIVERGENCE
C
  500 DIVO(I,J,K,it)=DUX+DVY
CC  UNITS:   DIVO(1/S)
  200 CONTINUE
  100 CONTINUE
      RETURN
      END SUBROUTINE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC  VERVEL COMPUTES VERTICAL VELOCITY
      SUBROUTINE VERVEL(PM,PW,WSFC,DIVO,OMGO,NXDIR,NYDIR,IT,ikk)
	PARAMETER (GRID=2.50)
      integer,parameter :: nrec=1460
	integer,parameter :: nrecr=1464
	integer,parameter :: nnr=nrecr
	integer,parameter :: lev=17 !! add surface
	integer,parameter :: n2d=13 !1lh 2sh 3pwr 4us 5vs 6 ts 7omegas 8Rhs 9Ps 10P_top
	integer,parameter :: n3d=6 !1u 1v 3t  4oemga 5RH 6HGT
	integer NXDIR,NYDIR,it
      DIMENSION WSFC(10,10,nnr),DIVO(10,10,lev,nnr),
     + PM(lev),OMGO(10,10,lev,nnr),PW(lev)
      common/DD/ data3d(6,6,lev,nnr,n3d), !1u 1v 3oemga 4t 5RH 6HGT
     +      data2d(10,10,nnr,n2d),hgt(10,10) ! !1lh 2sh 3pwr 4us 5vs 6 ts 7omegas 8Rhs 9Ps 10P_top
	COMMON/TQ/ tvs(10,10,nnr),qvs(10,10,nnr),
     +	tv(10,10,lev,nnr),qv(10,10,lev,nnr)

CC
      
       DO 200 I=1,NXDIR
       DO 200 J=1,NYDIR
!!!-------------------the follow code get the nearest stand level of the tropopause of every gird!
         p_top=data2d(I,J,it,10)/100.
	   pres=data2d(I,J,it,9)/100.
	   ktop=10  ! 100hPa
         do K=1,lev-1
		 if(pm(k)>p_top.and.pw(k)<P_top)then
		    ktop=k-1-ikk
		   goto 101
		 elseif(pw(k)>P_top.and.pm(k+1)<p_top)then
            ktop=k-ikk
             goto 101
		 endif
	   enddo
  101 CONTINUE
!      ktop=12
!       if(ktop(i,j)>11)print*,ktop(i,j)
!--------end of ensuring the level of tropopause	   	  
      OMGO(I,J,1,it)=WSFC(I,J,it)
      OMGO(I,J,2,it)=OMGO(I,J,1,it)+DIVO(I,J,2,it)
     + *(pres-PW(2))
       IF(pres.GT.PM(2))
     *    OMGO(I,J,2,it)=OMGO(I,J,1,it)+DIVO(I,J,1,it)*(pres-PW(2))
      DO 110 K=3,ktop-1      !!!! when the k=1ev, it is the tropopause, but actually it's not!!!! so need changed!!!
       IF(K.EQ.3 .AND. pres.GT.PW(2)) THEN
      OMGO(I,J,K,it)=OMGO(I,J,1,it)+DIVO(I,J,K,it)*(pres-PW(3))
       ELSE
	   dp=pm(K-1)-pm(k)
         OMGO(I,J,K,it)=OMGO(I,J,K-1,it)+DIVO(I,J,K,it)*dp   !!!50 d of pm
       ENDIF
  110 CONTINUE
  !-----------above tropopause--------------------------------- 
	  do k=ktop,lev-1
!	         tv_t=tv(ix,iy, ik, it+1)- tv(ix,iy, ik, it)
               
	         tv_p=0.5*(tv(I,J, k+1,it)-tv(I,J, k, it))
     +              + 0.5*(tv(i,J, k, it)-tv(i,J, k-1, it))
	         dp=pm(K-1)-pm(k)
	        OMGO(I,J,K,it)=-DIVO(I,J,K,it)
     +						 	*(tv_p/dp)	
        enddo
	 IF(pres.GT.PW(2)) THEN
      OMGO(I,J,2,it)=(OMGO(I,J,1,it)+OMGO(I,J,3,it))*0.5
       ENDIF
  200 CONTINUE
  100 CONTINUE
CC  UNITS:   OMGO(MB/S)
      RETURN
      END SUBROUTINE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC  ADJDV ADJUSTS DIVERGENCE AND VERTICAL VELOCITY USING O'BRIEN (1970)
      SUBROUTINE ADJDV(PM,PW,DIVO,OMGO,DIV,NXDIR,NYDIR,IT,ikk)
      PARAMETER (GRID=2.50)
      integer,parameter :: nrec=1460
	integer,parameter :: nrecr=1464
	integer,parameter :: nnr=nrecr
	integer,parameter :: lev=17 !! add surface
	integer,parameter :: n2d=13 !1lh 2sh 3pwr 4us 5vs 6 ts 7omegas 8Rhs 9Ps 10P_top
	integer,parameter :: n3d=6 !!1u 1v 3t  4oemga 5RH 6HGT
	integer NXDIR,NYDIR,it
      DIMENSION WSFC(10,10,nnr),DIVO(10,10,lev,nnr),
     + PM(lev),OMGO(10,10,lev,nnr),PW(lev),DIV(10,10,lev,nnr)
      common/DD/ data3d(6,6,lev,nnr,n3d), !1u 1v 3oemga 4t 5RH 6HGT
     +      data2d(10,10,nnr,n2d),hgt(10,10) ! !1lh 2sh 3pwr 4us 5vs 6 ts 7omegas 8Rhs 9Ps 10P_top
CC
      DO 410 I=1,NXDIR
      DO 410 J=1,NYDIR
!!!-------------------the follow code get the nearest stand level of the tropopause of every gird!
	 p_top=data2d(I,J,it,10)/100.
	   pres=data2d(I,J,it,9)/100.
	   ktop=10  ! 100hPa
         do K=1,lev-1
		 if(pm(k)>p_top.and.pw(k)<P_top)then
		    ktop=k-1-ikk
		   goto 101
		 elseif(pw(k)>P_top.and.pm(k+1)<p_top)then
            ktop=k-ikk
             goto 101
		 endif
	   enddo
  101 CONTINUE
!     ktop=12
C     DC=(0.001-OMGO(I,J,18))/(PS(I,J)-PW(18))
      omg_top=OMGO(I,J,ktop,it)
!	omg_top=0
      DC=(omg_top-OMGO(I,J,ktop-1,it))
     +	/(data2d(I,J,it,9)/100.-PW(ktop-1))  !!!top omega
      DO 500 K=1,ktop-1
      DIV(I,J,K,it)=DIVO(I,J,K,it)+DC
  500 CONTINUE
      do k=ktop,lev
      DIV(I,J,k,it)=DIVO(I,J,k,it)
	enddo
  410 CONTINUE
      RETURN
      END SUBROUTINE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC  Q1Q2  COMPUTES Q1 AND Q2 USING S, Q, U, V, OMG
CC  Note the dx and dy are vector, when the HAD is handle with, carefully note the data store squence
CC  and the array  subscript.
CC  NCEP: EAST-WEAT grids, the longitude array(x+1,y) is great than array(x,y)
CC        South North: the latitude array(x,y+1) is smaller than array(x,y)

      SUBROUTINE Q1Q2(PM,PW,OMG,IX,IY,IT,lo,ikk)
	  integer,parameter :: nrec=1460
	  integer,parameter :: nrecr=1464
	  integer,parameter :: nnr=nrecr
	  integer,parameter :: lev=17 !! add surface
	  integer,parameter :: n2d=13 !1lh 2sh 3pwr 4us 5vs 6 ts 7omegas 8Rhs 9Ps 10P_top
	  integer,parameter :: n3d=6 !1u 1v 3oemga 4t 5RH 6HGT
	  DIMENSION PM(lev),PW(lev),OMG(10,10,lev,nnr)
      DIMENSION STOS(10,10,lev),HADS(10,10,lev),
     +	VADS(10,10,lev)
      DIMENSION STOQ(10,10,lev),HADQ(10,10,lev),
     +	VADQ(10,10,lev)
	integer NXDIR,NYDIR,it,IX,IY,M,N,lo
      DIMENSION WSFC(10,10,nnr),DIVO(10,10,nnr,lev),
     + OMGO(10,10,lev,nnr),DIV(10,10,lev,nnr)
      common/DD/ data3d(6,6,lev,nnr,n3d), !1u 1v 3oemga 4t 5RH 6HGT
     +      data2d(10,10,nnr,n2d),hgt(10,10) ! !1lh 2sh 3pwr 4us 5vs 6 ts 7omegas 8Rhs 9Ps 10P_top
      COMMON/DXY/ DX(10),DY
      COMMON/TQ/ tvs(10,10,nnr),qvs(10,10,nnr),
     +	tv(10,10,lev,nnr),qv(10,10,lev,nnr)
      COMMON/Q12/ Q1(10,10,lev,nnr),Q2(10,10,lev,nnr),
     +	tls(10,10,lev,nnr),qls(10,10,lev,nnr)
	COMMON/AVEG/HAD_Q(10,10,lev,nnr),VAD_Q(10,10,lev,nnr),
     +	TCH_Q(10,10,lev,nnr),HADT(10,10,lev,nnr),VADT(10,10,lev,nnr),
     +	TCHT(10,10,lev,nnr) 
	  
      DT=6.*3600. 
      CP=0.24
      CL=600.
      C1=86400./CP/4.18684
      C2=-0.001*86400.*CL/CP
	TIMC=3600*24.
	  OMG=OMG*100   !mb/s to pa/s
      DO 50 M=1,IX
      DO 50 N=1,IY
!      VQ1(M,N)=0.
!      VQ2(M,N)=0.
      DO 50 K=1,lev
      STOS(M,N,K)=0.
      HADS(M,N,K)=0.
      VADS(M,N,K)=0.
      STOQ(M,N,K)=0.
      HADQ(M,N,K)=0.
      VADQ(M,N,K)=0.
   50 CONTINUE
      DO 200 M=2,IX-1
      DO 200 N=2,IY-1
      p_top=data2d(M,N,it,10)/100.
	   pres=data2d(M,N,it,9)/100.
	   ktop=10  ! 100hPa
         do iK=1,lev-1
		 if(pm(ik)>p_top.and.pw(ik)<P_top)then
		    ktop=ik-1-ikk
		   goto 101
		 elseif(pw(ik)>P_top.and.pm(ik+1)<p_top)then
            ktop=ik-ikk
             goto 101
		 endif
	   enddo
  101 CONTINUE
! the calculate above troposphause  ????
      DO 150 K=2,ktop-1
	  psf=data2d(M,N,IT,9)/100.
	  STOS(M,N,K)=0.
	  STOQ(M,N,K)=0.
	  HADS(M,N,K)=0.
	  VADS(M,N,K)=0.
	  HADQ(M,N,K)=0.
	  VADQ(M,N,K)=0.
!      IF(K.EQ.2 .AND. psf.LT.PM(K))then
!	 print*,k,psf,pm(k)
!	 GO TO 900   !!!find out why goto 900
!      endif
!      IF(K.EQ.2) GO TO 900  ! the second level convert to zero
	IF(K.EQ.2 .AND. psf.LT.PM(K)) GO TO 900   !!!!! the second level not zero
	IF(K.EQ.3 .AND. psf.LT.PM(K)) GO TO 900
!!!Time boundary----------------------------
      IF(IT==1)then
	  STOS(M,N,K)=(data3d(M,N,K,IT+2,3)-data3d(M,N,K,IT,3))/(2.*DT) ! 
      STOQ(M,N,K)=(QV(M,N,K,IT+2)-QV(M,N,K,IT))/(2.*DT)
	  ELSE
!!!!!!---------Time boundary
	  if(IT==nnr)then
	STOS(M,N,K)=(data3d(M,N,K,IT,3)-data3d(M,N,K,IT-2,3))/(2.*DT) ! SN= next time tv; SP=previous time
      STOQ(M,N,K)=(QV(M,N,K,IT)-QV(M,N,K,IT-2))/(2.*DT)
	  else
      STOS(M,N,K)=(data3d(M,N,K,IT+1,3)-data3d(M,N,K,IT-1,3))/(2.*DT) ! SN= next time tv; SP=previous time
      STOQ(M,N,K)=(QV(M,N,K,IT+1)-QV(M,N,K,IT-1))/(2.*DT)
	   endif
	  ENDIF
      IF(M.EQ.2) THEN
C     IF(MM.EQ.1 .OR. PS(M-2,N).LT.PM(K)) THEN
      MM=M+2
	DXX=0.
        IF((M+2)>IX)then
	      AXS=data3d(M,N,K,IT,1)*
     +  (data3d(M+1,N,K,IT,3)-data3d(M,N,K,IT,3))/(1.*DX(N))  !Temp.
       AXQ=data3d(M,N,K,IT,1)*
     + (QV(M,N,K,IT)-QV(M,N,K,IT))/(1.*DX(N))
	  else
      AXS=0.5*data3d(M+1,N,K,IT,1)*
     +  (data3d(M+2,N,K,IT,3)-data3d(M,N,K,IT,3))/(2.*DX(N))  !Temp.
      AXQ=0.5*data3d(M+1,N,K,IT,1)*
     + (QV(M+2,N,K,IT)-QV(M,N,K,IT))/(2.*DX(N))
	    endif
       ELSEIF(M.EQ.(IX-1)) THEN
C      ELSEIF(MM.EQ.7 .OR. PS(M+2,N).LT.PM(K)) THEN
!      AXS=0.5*U(M-1,N,K)*(S(M,N,K)-S(M-2,N,K))/(2.*DX(N))
	   AXS=0.5*data3d(M-1,N,K,IT,1)*
     +  (data3d(M,N,K,IT,3)-data3d(M-2,N,K,IT,3))/(2.*DX(N)) 
!      AXQ=0.5*U(M-1,N,K)*(Q(M,N,K)-Q(M-2,N,K))/(2.*DX(N))
       AXQ=0.5*data3d(M-1,N,K,IT,1)*
     + (QV(M,N,K,IT)-QV(M-2,N,K,IT))/(2.*DX(N))
	  
       ELSE

           IF((M+2)>IX)then
	      AXS=0.5*(data3d(M,N,K,IT,1)*(data3d(M+1,N,K,It,3)
     *   -data3d(M,N,K,IT,3)) +data3d(M-2,N,K,IT,1)*
     +  (data3d(M,N,K,IT,3)-data3d(M-1,N,K,IT,3)))/(DX(N))
      AXQ=0.5*(data3d(M,N,K,IT,1)*(QV(M+1,N,K,IT)-QV(M,N,K,IT))
     *  +data3d(M-2,N,K,IT,1)*(QV(M,N,K,IT)-QV(M-1,N,K,IT)))/(DX(N))
	      else
      AXS=0.5*(data3d(M+1,N,K,IT,1)*(data3d(M+2,N,K,It,3)
     *   -data3d(M,N,K,IT,3)) +data3d(M-1,N,K,IT,1)*
     +  (data3d(M,N,K,IT,3)-data3d(M-2,N,K,IT,3)))/(2.*DX(N))
      AXQ=0.5*(data3d(M+1,N,K,IT,1)*(QV(M+2,N,K,IT)-QV(M,N,K,IT))
     *  +data3d(M-1,N,K,IT,1)*(QV(M,N,K,IT)-QV(M-2,N,K,IT)))/(2.*DX(N))
	    endif
       ENDIF
      IF(N.EQ.2) THEN
C     IF(NN.EQ.1 .OR. PS(M,N-2).LT.PM(K)) THEN
!      AYS=0.5*V(M,N+1,K)*(S(M,N+2,K)-S(M,N,K))/(2.*DY)
       NN=N+2
	 DYY=0
        IF((N+2)>IY)then
        AYS=- data3d(M,N,K,IT,2)*(data3d(M,N+1,K,IT,3)  !!! why have a '-', for N+1 and N stand the actually
     +	  -data3d(M,N,K,IT,3))/(1.*DY)
	  AYQ=-data3d(M,N,K,It,2)*(QV(M,N+1,K,IT)
     +	  -QV(M,N,K,IT))/(1.*DY)
	  else
	  AYS=-0.5*data3d(M,N+1,K,IT,2)*(data3d(M,N+2,K,IT,3)
     +	  -data3d(M,N,K,IT,3))/(2.*DY)
	  AYQ=-0.5*data3d(M,N+1,K,It,2)*(QV(M,N+2,K,IT)
     +	  -QV(M,N,K,IT))/(2.*DY)
	  endif
!      if(abs(AYS*TIMC)>30)then
!	print*,'1XXXXX'
!      print*, AYS*TIMC
!      print*, data3d(M,N+2,K,IT,3)-data3d(M,N,K,IT,3)
!	print*, data3d(M,N+1,K,IT,2)
!	print*,M,N,K,'lo',lo
!	stop
!	endif
!      AYQ=0.5*V(M,N+1,K)*(Q(M,N+2,K)-Q(M,N,K))/(2.*DY)
       ELSEIF(N.EQ.(IY-1).and.(N-2)>0) THEN
C      ELSEIF(NN.EQ.7 .OR. PS(M,N+2).LT.PM(K)) THEN
!      AYS=0.5*V(M,N-1,K)*(S(M,N,K)-S(M,N-2,K))/(2.*DY)
	  AYS=-0.5*data3d(M,N-1,K,IT,2)*(data3d(M,N,K,IT,3)
     +	  -data3d(M,N-2,K,IT,3))/(2.*DY)
!	if(abs(AYS*TIMC)>30)then
!	print*,'2XXXXX'
!      print*, AYS*TIMC
!      print*, data3d(M,N,K,IT,3)-data3d(M,N-2,K,IT,3)
!	print*, data3d(M,N-1,K,IT,2)
!	print*,M,N,K,'lo',lo
!	stop
!	endif
!      AYQ=0.5*V(M,N-1,K)*(Q(M,N,K)-Q(M,N-2,K))/(2.*DY)
	 AYQ=-0.5*data3d(M,N-1,K,IT,2)*(QV(M,N,K,IT)
     +	 -QV(M,N-2,K,IT))/(2.*DY)
       ELSE
!      AYS=0.5*(V(M,N+1,K)*(S(M,N+2,K)-S(M,N,K))
!     *        +V(M,N-1,K)*(S(M,N,K)-S(M,N-2,K)))/(2.*DY)
!      AYQ=0.5*(V(M,N+1,K)*(Q(M,N+2,K)-Q(M,N,K))
!     *        +V(M,N-1,K)*(Q(M,N,K)-Q(M,N-2,K)))/(2.*DY)
        IF((N+2)>IY)then
	  AYS=0.5*(data3d(M,N,K,IT,2)*(-data3d(M,N+1,K,It,3) + !-
     *   data3d(M,N,K,IT,3))+data3d(M,N-1,K,IT,2)*(-data3d(M,N,K,IT,3)+  ! - 
     +    data3d(M,N-1,K,IT,3)))/(1.*DY)
        AYQ=0.5*(data3d(M,N,K,IT,2)*(-QV(M,N+1,K,IT)+QV(M,N,K,IT))
     *   +data3d(M,N-1,K,IT,2)*(-QV(M,N,K,IT)+QV(M,N-1,K,IT)))/(1.*DY)
	 else
	 AYS=0.5*(data3d(M,N+1,K,IT,2)*(-data3d(M,N+2,K,It,3) + !-
     *   data3d(M,N,K,IT,3))+data3d(M,N-1,K,IT,2)*(-data3d(M,N,K,IT,3)+  !-
     +    data3d(M,N-2,K,IT,3)))/(2.*DY)
       AYQ=0.5*(data3d(M,N+1,K,IT,2)*(-QV(M,N+2,K,IT)+QV(M,N,K,IT))
     *   +data3d(M,N-1,K,IT,2)*(-QV(M,N,K,IT)+QV(M,N-2,K,IT)))/(2.*DY)
	  endif
       ENDIF
      HADS(M,N,K)=AXS+AYS
!	if(abs(HADS(M,N,K)*TIMC)>30)then
!      print*, HADS(M,N,K)*TIMC,AXS*TIMC,AYS*TIMC
!      print*,M,N,K,'lo',lo
!	print*,'*******************************************'
!	stop
!	endif

      HADQ(M,N,K)=AXQ+AYQ
          UKP1=.5*(tv(M,N,K+1,IT)-tv(M,N,K,IT))*OMG(M,N,K,IT)
          UKM1=.5*(tv(M,N,K,IT)-tv(M,N,K-1,IT))*OMG(M,N,K-1,IT)
          VKP1=.5*(QV(M,N,K+1,IT)-QV(M,N,K,IT))*OMG(M,N,K,IT)
          VKM1=.5*(QV(M,N,K,IT)-QV(M,N,K-1,IT))*OMG(M,N,K-1,IT)
          DPW=Pm(K)-Pm(K-1)
       IF(K.GE.3) GO TO 850   !!! k=5 in hmbudps is 800hPa this sentence is think that surface pressure great than 800hPa
!-------in my calculate, when k=3, the pressure=850, k=4,p=700, so changge to K GE 3	   
!       IF(K.EQ.4 .AND. psf.GE.PM(3)) GO TO 850
!         IF(K.EQ.4 .AND. psf.LT.PM(3)) THEN
!          UKM1=.5*(S(M,N,K)-S(M,N,1))*OMG(M,N,K-1)
!          VKM1=.5*(Q(M,N,K)-Q(M,N,1))*OMG(M,N,K-1)
!         ENDIF
       IF(K.EQ.3 .AND. psf.GE.PM(2)) GO TO 850
         IF(K.EQ.3 .AND. psf.LT.PW(2)) THEN
          UKM1=.5*(tv(M,N,K,It)-tv(M,N,1,IT))*OMG(M,N,1,IT)
          VKM1=.5*(QV(M,N,K,It)-QV(M,N,1,IT))*OMG(M,N,1,IT)
          DPW=PW(K)-psf
         ENDIF
         IF(K.EQ.3 .AND. psf.GE.PW(2)) THEN
          UKM1=.5*(tv(M,N,K,IT)-tv(M,N,1,IT))*OMG(M,N,K-1,It)
          VKM1=.5*(QV(M,N,K,IT)-QV(M,N,1,IT))*OMG(M,N,K-1,IT)
         ENDIF
         IF(K.EQ.2 .AND. psf.GE.PM(2)) THEN
           DPW=PW(K)-psf
         ENDIF
  850 VADS(M,N,K)=((pm(k)/1000)**0.286)*(UKP1+UKM1)/(DPW*100) !(UKP1+UKM1)/DPW 
      VADQ(M,N,K)=(VKP1+VKM1)/(DPW*100)
  900 CONTINUE
  150 CONTINUE
  200 CONTINUE
      open(999,file='F999.txt',status='unknown',position='APPEND')
      DO 11 I=1,IX
      DO 11 J=1,IY
      DO 11 K=1,lev
	  tls(I, J, k, IT)=  -HADS(I,J,K) -VADS(I,J,K) 
	  qls(I, J, k, IT)=    -HADQ(I,J,K) -VADQ(I,J,K)
      Q1(I,J,K,IT)=(STOS(I,J,K)+HADS(I,J,K)+VADS(I,J,K))
      Q2(I,J,K,IT)=-(STOQ(I,J,K)+HADQ(I,J,K)+VADQ(I,J,K))
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	HAD_Q(I,J,K,IT)=-HADQ(I,J,K)  !!! Q2 -
	VAD_Q(I,J,K,IT)=-VADQ(I,J,K)  !!!  Q2 -
     	TCH_Q(I,J,K,IT)=-STOQ(I,J,K)   !!!  Q2 -
	HADT(I,J,K,IT)=HADS(I,J,K)
	VADT(I,J,K,IT)=VADS(I,J,K) 
      TCHT(I,J,K,IT)=STOS(I,J,K)
C     Q1(I,J,K)=(HADS(I,J,K)+VADS(I,J,K))*C1
C     Q2(I,J,K)=(HADQ(I,J,K)+VADQ(I,J,K))*C2
   11 CONTINUE
CCC   UNIT:   DEG/DAY
C      DO 20 I=1,15
C      DO 20 J=1,15
C      DO 30 K=2,18
C       DPW=PW(K-1)-PW(K)
C     VQ1(I,J)=VQ1(I,J)+Q1(I,J,K)*DPW/(PW(1)-PW(18))
C     VQ2(I,J)=VQ2(I,J)+Q2(I,J,K)*DPW/(PW(1)-PW(18))
C       IF(K.EQ.2) DPW=PS(I,J)-PW(K)
C     VQ1(I,J)=VQ1(I,J)+Q1(I,J,K)*DPW*CP/CL/(9.8*24.)*100.
C     VQ2(I,J)=VQ2(I,J)+Q2(I,J,K)*DPW*CP/CL/(9.8*24.)*100.
C    UNIT:   MM/H
C     VQ1(I,J)=VQ1(I,J)+Q1(I,J,K)*DPW/(PS(I,J)-PW(18))
C     VQ2(I,J)=VQ2(I,J)+Q2(I,J,K)*DPW/(PS(I,J)-PW(18))
C      VQ1(I,J)=VQ1(I,J)+Q1(I,J,K)*DPW*CP*1000.*4.18/9.8/86400.
C      VQ2(I,J)=VQ2(I,J)+Q2(I,J,K)*DPW*CP*1000.*4.18/9.8/86400.
C    UNIT: 100W/M**2
C   40 CONTINUE
C   30 CONTINUE
C   20 CONTINUE
      RETURN
      END SUBROUTINE
C
      subroutine dynamic(d3d,OMGP,DIV,IX,IY,Blat,IT,dyn) !IXX(lo),IYY(lo),bslat(lo),IT)
	PARAMETER (GRID=2.50)
	integer,parameter :: nrec=1460
	integer,parameter :: nrecr=1464
	integer,parameter :: nnr=nrecr
	integer,parameter :: lev=17 !! add surface
	integer,parameter :: n2d=13 !1lh 2sh 3pwr 4us 5vs 6 ts 7omegas 8Rhs 9Ps 10P_top, 11cpr 12DLR 13ULR
	integer,parameter :: n3d=6 !1u 1v 3t  4oemga 5RH 6HGT
	real d3d(6,6,lev,nnr,n3d),OMG(10,10,lev,nnr)
      real OMGP(10,10,lev,nnr),blat
	real vor(6,6,lev,nnr,3),temp(4)
	real DIV(10,10,lev,nnr)
	real dyn(lev,nnr,4) ! 1-3 for vorticity 4 for divergence
	integer IX,IY,it
	DATA PIE/0.0174532925/

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------  the calculation of vorticity ---------------------------------
!  vorticity=(domega/dy-dv/dz)i+(du/dz-domega/dx)j+(dv/dx-du/dy)k
!
!
!
      DY=1./(GRID*111.17E3)
      do 101 k=1,lev
      do 101 i=1,IX
	do 101 j=1,Iy
	tempx=-9.8*461.
      OMG(i,j,k,it)=OMGP(i,j,k,it)*100./tempx
101   continue

	do k=1,lev
      if(k==1) goto 200
	do 100 i=1,IX
	do 100 j=1,IY
	Y=((IY-J)*GRID+blat)/180.
      DX=1./(COS(Y*PIE)*GRID*111.17E3)
	dz=d3d(i,j,k,it,6)-d3d(i,j,k-1,it,6)
!------------------------------------------------
      if(i==1)then
      dwy=OMG(i+2,j,k,it)-OMG(i,j,k,nnr)
      dvw=d3d(i+2,j,k,it,2)-d3d(i,j,k,it,2)
	elseif(i==IX)then
      dwy=OMG(i,j,k,it)-OMG(i-2,j,k,nnr)
      dvw=d3d(i,j,k,it,2)-d3d(i-2,j,k,it,2)
	else
      dwy=OMG(i+1,j,k,it)-OMG(i-1,j,k,nnr)
      dvw=d3d(i+1,j,k,it,2)-d3d(i-1,j,k,it,2)
	endif
	
	if(J==1)then
	dwx=OMG(i,j+2,k,it)-OMG(i,j,k,nnr)
	duw=d3d(i,j+2,k,it,1)-d3d(i,j,k,it,1)
	elseif(j==iy)then
      dwx=OMG(i,j,k,it)-OMG(i,j-2,k,nnr)
	duw=d3d(i,j,k,it,1)-d3d(i,j-2,k,it,1)
	else
	dwx=OMG(i,j+1,k,it)-OMG(i,j-1,k,nnr)
	duw=d3d(i,j+1,k,it,1)-d3d(i,j-1,k,it,1)
	endif
	dvx=d3d(i,j,k,it,2)-d3d(i,j,k-1,it,2)
      duy=d3d(i,j,k,it,1)-d3d(i,j,k-1,it,1)
	dyy=dy/2.0
	dxx=dx/2.0 
	vorX=dwx*dyy-dvx/dz
	vorW=dvw*dxx-duw*dyy
	vorY=duy/dz-dwy*dxx

      vor(i,j,k,it,1)=vorX
	vor(i,j,k,it,2)=vorY
	vor(i,j,k,it,3)=vorW
100   continue
200   continue
      enddo
      do k=1,lev
	temp=0.0
	do 103 i=2,IX-1
	do 103 j=2,IY-1
        do ip=1,3
        temp(ip)=temp(ip)+vor(i,j,k,it,ip)
	  enddo
	  temp(4)=temp(4)+div(i,j,k,it)
103   continue
        do ip=1,4
	   xy=1.0/(ix-2.0)*(iy-2.0)
         dyn(k,it,ip)=temp(ip)*XY
	  enddo
	enddo
	
	return
	end subroutine
C
      SUBROUTINE grads(year,folds,nt)
      integer, parameter :: nl=4
	integer, parameter :: lev=18
	real,allocatable :: tls(:,:,:),qls(:,:,:)
	integer  nt
	real plv(lev)
	character*100 filename,dir,filename1
	character*30 area(nl),folds
	character*4 year,month
	data plv/1003, 1000.,925., 850.,700.,
     +     600., 500.,400.,300.,250.,200.,150.,100.,70,50,30,20,10/
      allocate(tls(nt,lev,nl))
	allocate(qls(nt,lev,nl))
	!----------------------------------------------------------------
	area(1)='PRD'
	area(2)='MLYR'
	area(3)='NPC'
	area(4)='NEC'
	dir='Z:\DATA\LargeScale\NcepR2_Pre\'
!	year='2010'
!      folds='100101-100930'
	
	do in=1,nl
	filename=trim(dir)//trim(folds)//'\'//trim(area(in))//
     +'_forcing.txt'
	open(in*10,file=trim(filename))
	read(in*10,*)
	read(in*10,*)
	do it=1,nt
	 read(in*10,103)IDTIME,(tls(it,ik,in),ik=1,lev),
     +	(qls(it,ik,in),ik=1,lev)
	 enddo
	enddo
 !
      filename1=trim(dir)//'NCEP_BIN\'//year//'Forcing_Temp.dat'      
	open(99,file=filename1,form='binary')
	filename=trim(dir)//'NCEP_BIN\'//year//'Forcing_Qv.dat'      
	open(999,file=filename,form='binary')
	do it=1,nt
	  do in=1,nl
          do ik=1,lev
           write(99) tls(it,ik,in)
	     write(999)Qls(it,ik,in)
           enddo
	    enddo
	  enddo
103   format(1X,I4,36(1X,e12.4))
      open(101,file='Z:\DATA\LargeScale\NcepR2_Pre\NCEP_CTL\'
     +//year//'grads_Q.ctl')
	write(101,*)'dset',' ',trim(filename)
	write(101,*)'title Forcing'
	write(101,*)'options little_endian'
	write(101,*)'undef -9999.'
	write(101,*)'xdef  1 linear   1 1'
	write(101,*)'ydef  1 linear   1 1'
	write(101,102)'zdef',lev, 'levels', (plv(ik),ik=1,lev)
      write(101,104)'tdef', nt, 'linear 00z1jan',year,'6hr'
	write(101,105)'vars',nl
	do in=1,nl
	write(101,106)trim(area(in)),lev,'99 QLS'
	enddo
	write(101,*)'endvars'
	close(10)
!
	open(101,file='Z:\DATA\LargeScale\NcepR2_Pre\NCEP_CTL\'
     +//year//'grads_T.ctl')
	write(101,*)'dset',' ',trim(filename1)
	write(101,*)'title Forcing'
	write(101,*)'options little_endian'
	write(101,*)'undef -9999.'
	write(101,*)'xdef  1 linear   1 1'
	write(101,*)'ydef  1 linear   1 1'
	write(101,102)'zdef',lev, 'levels', (plv(ik),ik=1,lev)
      write(101,104)'tdef', nt, 'linear 00z1jan',year,'6hr'
	write(101,105)'vars',nl
	do in=1,nl
	write(101,106)trim(area(in)),lev,'99 TLS'
	enddo
	write(101,*)'endvars'
	close(10)
102   format(1X,A4,1X,I4,1X,A6,21(1X,F7.2))
104   format(1X,A4,1X,I4,1X,A14,A4,1X,A3)
105   format(1X,A4,1X,I3)
106   format(1X,A4,1X,I3,1X,A10)
      end SUBROUTINE
!!
      SUBROUTINE grads_PRE(year,folds,nt)
!	integer, parameter :: nt=1096 !!!!!
      integer, parameter :: nl=4
	integer, parameter :: lev=1
	integer nt
	real,allocatable :: lh(:,:)
	real,allocatable :: sh(:,:)
	real,allocatable :: pr(:,:)
	real,allocatable :: cpr(:,:)
	real plv(lev)
	character*100 filename(4),dir,filename1
	character*30 area(nl),folds,var(4)
	character*4 year,month
	data plv/1003 /!, 1000.,925., 850.,700.,
 !    +     600., 500.,400.,300.,250.,200.,150.,100.,70,50,30,20,10/
	allocate(lh(nt,nl))
	allocate(sh(nt,nl))
	allocate(pr(nt,nl))
	allocate(cpr(nt,nl))
	!----------------------------------------------------------------
	area(1)='PRD'
	area(2)='MLYR'
	area(3)='NPC'
	area(4)='NEC'
	var(1)='LH'
      var(2)='SH'
      var(3)='Pre'
      var(4)='Cpre'
	dir='Z:\DATA\LargeScale\NcepR2_Pre\'
!	year='2012'
!      folds='120101-120930'
	
	do in=1,nl
	filename1=trim(dir)//trim(folds)//'\'//trim(area(in))//
     +'_HF.txt'
	open(in*10,file=trim(filename1))
	read(in*10,*)
	read(in*10,*)
	read(in*10,*)
	do it=1,nt
	 read(in*10,103)IDTIME,lh(it,in),sh(it,in),pr(it,in),cpr(it,in)
	 enddo
	enddo
 !
      filename(1)=trim(dir)//'NCEP_BIN\'//year//'Surface_LH.dat'      
	  open(99,file=trim(filename(1)),form='binary')
	filename(2)=trim(dir)//'NCEP_BIN\'//year//'Surface_SH.dat'      
	  open(999,file=trim(filename(2)),form='binary')
	filename(4)=trim(dir)//'NCEP_BIN\'//year//'Surface_CPRE.dat' 
	 open(997,file=trim(filename(4)),form='binary')
      filename(3)=trim(dir)//'NCEP_BIN\'//year//'Surface_PRE.dat'      
	  open(998,file=trim(filename(3)),form='binary')
	do it=1,nt
	  do in=1,nl
           write(99) lh(it,in)
	     write(999)sh(it,in)
	     write(998)pr(it,in)*3600  ! mm/hr
		 write(997)cpr(it,in)*3600 ! to mm/hr
	    enddo
	  enddo
103   format(1X,I4,1X,F8.2,1X,F8.2,1X,e12.4,1X,e12.4)
      do i=1,4
      open(101,file='Z:\DATA\LargeScale\NcepR2_Pre\NCEP_CTL\'
     +//year//'grads_'//trim(var(i))//'.ctl')
	write(101,*)'dset',' ',trim(filename(i))
	write(101,*)'title Forcing'
	write(101,*)'options little_endian'
	write(101,*)'undef -9999.'
	write(101,*)'xdef  1 linear   1 1'
	write(101,*)'ydef  1 linear   1 1'
	write(101,102)'zdef',lev, 'levels', (plv(ik),ik=1,lev)
      write(101,104)'tdef', nt, 'linear 00z1jan',year,'6hr'
	write(101,105)'vars',nl
	do in=1,nl
	write(101,106)trim(area(in)),lev,'99',trim(var(i))
	enddo
	write(101,*)'endvars'
	close(101)
      enddo
102   format(1X,A4,1X,I4,1X,A6,21(1X,F7.2))
104   format(1X,A4,1X,I4,1X,A14,A4,1X,A3)
105   format(1X,A4,1X,I3)
106   format(1X,A4,1X,I3,1X,A3,1X,A6)
      end SUBROUTINE

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE NCEP_SAS(PM,PW,OMG,IX,IY,IT,iks,ike,XY_OUT)
	  integer,parameter :: nrec=1460
	  integer,parameter :: nrecr=1464
	  integer,parameter :: nnr=nrecr
	  integer,parameter :: lev=17 !! add surface
	  integer,parameter :: n2d=13 !1lh 2sh 3pwr 4us 5vs 6 ts 7omegas 8Rhs 9Ps 10P_top
	  integer,parameter :: n3d=6 !1u 1v 4oemga 3t 5RH 6HGT
	  DIMENSION PM(lev),PW(lev),OMG(10,10,lev,nnr)
	real temp(20),XY_OUT(nnr,20)
      DIMENSION STOS(10,10,lev),HADS(10,10,lev),
     +	VADS(10,10,lev)
      DIMENSION STOQ(10,10,lev),HADQ(10,10,lev),
     +	VADQ(10,10,lev)
	integer NXDIR,NYDIR,it,IX,IY,M,N,lo
      DIMENSION WSFC(10,10,nnr),DIVO(10,10,nnr,lev),
     + OMGO(10,10,lev,nnr),DIV(10,10,lev,nnr)
      common/DD/ data3d(6,6,lev,nnr,n3d), !1u 1v 3oemga 4t 5RH 6HGT
     +      data2d(10,10,nnr,n2d),hgt(10,10) ! !1lh 2sh 3pwr 4us 5vs 6 ts 7omegas 8Rhs 9Ps 10P_top
      COMMON/TQ/ tvs(10,10,nnr),qvs(10,10,nnr),
     +	tv(10,10,lev,nnr),qv(10,10,lev,nnr)
      COMMON/Q12/ Q1(10,10,lev,nnr),Q2(10,10,lev,nnr),
     +	tls(10,10,lev,nnr),qls(10,10,lev,nnr)
	COMMON/AVEG/HAD_Q(10,10,lev,nnr),VAD_Q(10,10,lev,nnr),
     +	TCH_Q(10,10,lev,nnr),HADT(10,10,lev,nnr),VADT(10,10,lev,nnr),
     +	TCHT(10,10,lev,nnr) 
	 
	 
	temp=0.0 
	XY_OUT=0.0
	hlat=2.5e6
	cp=1005.
!---------low atomosphere
	 do 197 ik=iks,ike
	 do 197 ix=2,IX-1
	   do 197 iy=2,IY-1
	          temp(1)= temp(1)+tls(ix,iy, ik, it) *3600*24  !day-1
	          temp(2)= temp(2)+qls(ix,iy, ik, it) *3600*24*hlat/cp !!! K/day day-1
	          temp(3)= temp(3)+DATA3D(ix,iy,ik,it,1)
	          temp(4)= temp(4)+DATA3D(ix,iy,ik,it,2)
	          temp(5)= temp(5)+qv(ix,iy, ik, it)
	          temp(6)= temp(6)+DATA3D(ix,iy,ik,it,6)/9.87  !HGT
	          temp(7)= temp(7)+DATA3D(ix,iy,ik,it,3)
	          temp(8)= temp(8)+omg(ix,iy, ik, it)
	          temp(9)= temp(9)+DATA3D(ix,iy,ik,it,5)
!	          temp(10)= temp(10)+qv(ix,iy, ik, it)
	          temp(10)= temp(10)+tv(ix,iy, ik, it)
	          temp(11)= temp(11)+Q1(ix, iy, ik,it) *3600*24   !Tmcl
	          temp(12)= temp(12)+Q2(ix, iy, ik,it) *3600*24*hlat/cp   !*TMCL*hlat/cp     !k/day
		      temp(13)= temp(13)+HAD_Q(ix,iy, ik, it) *3600*24*hlat/cp  ! day-1
	          temp(14)= temp(14)+VAD_Q(ix,iy, ik, it) *3600*24*hlat/cp   !day-1
	          temp(15)= temp(15)+TCH_Q(ix,iy, ik, it) *3600*24*hlat/cp  !!! K/day day-1
			  temp(16)= temp(16)+HADT(ix,iy, ik, it) *3600*24  !day-1
	          temp(17)= temp(17)+VADT(ix,iy, ik, it) *3600*24  !!! K/day day-1
			  temp(18)= temp(18)+TCHT(ix,iy, ik, it) *3600*24 !! K/day day-1
	          temp(19)= temp(19)+data3d(ix,iy,ik,it,4)
197   continue
	XY=(IX-2.0)*(IY-2.)
	SK=ike-iks+1.0
      XYS=XY*SK
	do ii=1,19
	XY_OUT(it,ii)=temp(ii)/XYS
      enddo

	return
	end SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE NCEP_NOSAS(PM,PW,OMG,IX,IY,IT,XYN_OUT)
	  integer,parameter :: nrec=1460
	  integer,parameter :: nrecr=1464
	  integer,parameter :: nnr=nrecr
	  integer,parameter :: lev=17 !! add surface
	  integer,parameter :: n2d=13 !1lh 2sh 3pwr 4us 5vs 6 ts 7omegas 8Rhs 9Ps 10P_top
	  integer,parameter :: n3d=6 !1u 1v 3oemga 4t 5RH 6HGT
	  DIMENSION PM(lev),PW(lev),OMG(10,10,lev,nnr)
	real temp(20),XYN_OUT(nnr,lev,20)
      DIMENSION STOS(10,10,lev),HADS(10,10,lev),
     +	VADS(10,10,lev)
      DIMENSION STOQ(10,10,lev),HADQ(10,10,lev),
     +	VADQ(10,10,lev)
	integer NXDIR,NYDIR,it,IX,IY,M,N,lo
      DIMENSION WSFC(10,10,nnr),DIVO(10,10,nnr,lev),
     + OMGO(10,10,lev,nnr),DIV(10,10,lev,nnr)
      common/DD/ data3d(6,6,lev,nnr,n3d), !1u 1v 3oemga 4t 5RH 6HGT
     +      data2d(10,10,nnr,n2d),hgt(10,10) ! !1lh 2sh 3pwr 4us 5vs 6 ts 7omegas 8Rhs 9Ps 10P_top
      COMMON/TQ/ tvs(10,10,nnr),qvs(10,10,nnr),
     +	tv(10,10,lev,nnr),qv(10,10,lev,nnr)
      COMMON/Q12/ Q1(10,10,lev,nnr),Q2(10,10,lev,nnr),
     +	tls(10,10,lev,nnr),qls(10,10,lev,nnr)
	COMMON/AVEG/HAD_Q(10,10,lev,nnr),VAD_Q(10,10,lev,nnr),
     +	TCH_Q(10,10,lev,nnr),HADT(10,10,lev,nnr),VADT(10,10,lev,nnr),
     +	TCHT(10,10,lev,nnr) 
	 
	 
	temp=0.0 
	XY_OUT=0.0
	hlat=2.5e6
	cp=1005.
!---------low atomosphere
	 do  ik=1,lev
	 do 197 ix=2,IX-1
	   do 197 iy=2,IY-1
	          temp(1)= temp(1)+tls(ix,iy, ik, it) *3600*24  !day-1
	          temp(2)= temp(2)+qls(ix,iy, ik, it) *3600*24*hlat/cp !!! K/day day-1
	          temp(3)= temp(3)+DATA3D(ix,iy,ik,it,1)
	          temp(4)= temp(4)+DATA3D(ix,iy,ik,it,2)
	          temp(5)= temp(5)+qv(ix,iy, ik, it)
	          temp(6)= temp(6)+DATA3D(ix,iy,ik,it,6)/9.87  !HGT
	          temp(7)= temp(7)+DATA3D(ix,iy,ik,it,3)
	          temp(8)= temp(8)+omg(ix,iy, ik, it)
	          temp(9)= temp(9)+DATA3D(ix,iy,ik,it,5)
!	          temp(10)= temp(10)+qv(ix,iy, ik, it)
	          temp(10)= temp(10)+tv(ix,iy, ik, it)
	          temp(11)= temp(11)+Q1(ix, iy, ik,it) *3600*24   !Tmcl
	          temp(12)= temp(12)+Q2(ix, iy, ik,it) *3600*24*hlat/cp   !*TMCL*hlat/cp     !k/day
		      temp(13)= temp(13)+HAD_Q(ix,iy, ik, it) *3600*24*hlat/cp  ! day-1
	          temp(14)= temp(14)+VAD_Q(ix,iy, ik, it) *3600*24*hlat/cp   !day-1
	          temp(15)= temp(15)+TCH_Q(ix,iy, ik, it) *3600*24*hlat/cp  !!! K/day day-1
			  temp(16)= temp(16)+HADT(ix,iy, ik, it) *3600*24  !day-1
	          temp(17)= temp(17)+VADT(ix,iy, ik, it) *3600*24  !!! K/day day-1
			  temp(18)= temp(18)+TCHT(ix,iy, ik, it) *3600*24 !! K/day day-1
	          temp(19)= temp(19)+DATA3D(ix,iy,ik,it,4)

197   continue
	XY=(IX-2.0)*(IY-2.)
		do ii=1,19
	    XYN_OUT(it,ik,ii)=temp(ii)/XY
          enddo
	temp=0.0
	enddo
!	SK=ike-iks+1.0
!      XYS=XY*SK

	return
	end SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       SUBROUTINE NCEP_GUS(PM,PW,IX,IY,IT,iks,ike,XY_OUT)
	  integer,parameter :: nrec=1460
	  integer,parameter :: nrecr=1464
	  integer,parameter :: nnr=nrecr
	  integer,parameter :: lev=17 !! add surface
	  integer,parameter :: n2d=13 !1lh 2sh 3pwr 4us 5vs 6 ts 7omegas 8Rhs 9Ps 10P_top
	  integer,parameter :: n3d=6 !1u 1v 3oemga 4t 5RH 6HGT
	  DIMENSION PM(lev),PW(lev),OMG(10,10,lev,nnr)
	real temp(20),XY_OUT(nnr,10)
      DIMENSION STOS(10,10,lev),HADS(10,10,lev),
     +	VADS(10,10,lev)
      DIMENSION STOQ(10,10,lev),HADQ(10,10,lev),
     +	VADQ(10,10,lev)
	integer NXDIR,NYDIR,it,IX,IY,M,N,lo
      DIMENSION WSFC(10,10,nnr),DIVO(10,10,nnr,lev),
     + OMGO(10,10,lev,nnr),DIV(10,10,lev,nnr)
      common/DD/ data3d(6,6,lev,nnr,n3d), !1u 1v 3oemga 4t 5RH 6HGT
     +      data2d(10,10,nnr,n2d),hgt(10,10) ! !1lh 2sh 3pwr 4us 5vs 6 ts 7omegas 8Rhs 9Ps 10P_top
      COMMON/DXY/ DX(10),DY
      COMMON/TQ/ tvs(10,10,nnr),qvs(10,10,nnr),
     +	tv(10,10,lev,nnr),qv(10,10,lev,nnr)
      COMMON/Q12/ Q1(10,10,lev,nnr),Q2(10,10,lev,nnr),
     +	tls(10,10,lev,nnr),qls(10,10,lev,nnr)
	COMMON/AVEG/HAD_Q(10,10,lev,nnr),VAD_Q(10,10,lev,nnr),
     +	TCH_Q(10,10,lev,nnr),HADT(10,10,lev,nnr),VADT(10,10,lev,nnr),
     +	TCHT(10,10,lev,nnr) 
	 
	 
	temp=0.0 
	XY_OUT=0.0
!---------low atomosphere
	 do 197 ik=iks,ike
	 do 197 ix=2,IX-1
	   do 197 iy=2,IY-1
	      do 197 iv=1,6
		     if(iv<=3)then
			 temp(iv)=temp(iv)+ data2d(ix,iy,it,iv)
			 else
			 temp(iv)=temp(iv)+ data2d(ix,iy,it,iv+7)
			 endif

197   continue
	XY=(IX-2.0)*(IY-2.)
	SK=ike-iks+1.0
      XYS=XY*SK
	do ii=1,6
	XY_OUT(it,ii)=temp(ii)/XYS
      enddo

	return
	end SUBROUTINE
!
!***********************************************************************
      subroutine calendar(iys,iye,ims,ime,ids,ide,it,cder)
      integer days(12)
	character*5 hour(5),houro
	character*30 cder
	character*4 year
	character*2 month,day


      do i=1,12
	  days(i)=31
	enddo
	hour(1)='00:00'
	hour(2)='06:00'
	hour(3)='12:00'
	hour(4)='18:00'
	 days(2)=28
	 days(6)=30
	 days(4)=30
	 days(9)=30
	 days(11)=30
	if(mod(iys,4)==0.and.mod(iys,100)/=0)then
	     days(2)=29 ! nrec=1464
	elseif(mod(iys,400)==0)then
	     days(2)=29
	endif
	write(year,'(I4)')iys
	write(month,'(I2.2)')ims
	write(day,'(I2.2)')ids
	ittb=0
	do i=1,ims-1
	ittb1=ittb1+days(i)*4
	enddo
	ittb2=(ids-1)*4
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	do i=1,ims-1
	ittb=ittb+days(i)*4
	enddo
	idays=(it-ittb)/4
	im=ims
!	print*,idays
	idd=0
	do 99 i=ims,ime
	idays1=idays-idd+1
	if(mod(it,4)==0)idays1=idays1-1
!	print*,i,idays1,days(i),idd
	if(idays1<(days(i)+1))then
!	print*,'*********************'
!	print*,i,idays1,days(i),idd
      write(day,'(I2.2)')idays1
      write(month,'(I2.2)')i
	goto 100
	else
	idd=idd+days(i)
      goto 99
	endif
99	continue
100   continue
      if(mod(it,4)==0)then
	j=4
	else
	j=mod(it,4)
	endif
	houro=hour(j)

      cder=year//'-'//month//'-'//day//' '//houro
	

      end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
!--------------Reading code-----------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine NC_2D(lo,nrec,itt,i2d,data2d,filename,IXG,IYG,
     +              XG1,YG1,XG2,YG2,nnr,n2d)
      include 'netcdf.inc'
!
!-------------------------------------------------------
!
!    Using below command line to get the head information of the NetCDF file:
!    Please execute " ncdump -h Z:\DATA\LargeScale\NCEP_Surface\lhtfl.sfc.gauss.1990.nc > output.txt " 
!-------------------------------------------------------
!
!     Define Variables.
!     Variable ids run sequentially from 1 to nvars=4 ! number of variables
      integer :: nrec,nnr,n2d ! change this to generalize
      integer*4  ncid, status    ! file control
      integer*4 recdim   ! record dimension
!-------------------------------------------------------------
 !     Below 4 variables is the data in netCDF file
      real*4         ::  lat( 94 )
      real*4         ::  lon( 192 )
      real*8 , allocatable ::  time(:)
      integer*2, allocatable  ::  f2d( :, :, : )
!	real, allocatable :: Temp2d( :, :, : )
!     above4 variables is the data in netCDF file
!-------------------------------------------------------------
      integer*4   :: start(10)
      integer*4   :: count(10)
      integer*4   :: dimids(10)! allow up to 10 dimensions
      integer*4   :: dimid, xtype
      character(len=31) :: dummy
!-------------addd--------------------------------------------------
      integer :: itt,i2d	
!	real data2d(15,10,nrecr*10,2)! once max is 10 years
	character*150 filename
	integer IXG(4),IYG(4)
	integer XG1(4),YG1(4),XG2(4),YG2(4) 
	real data2d(10,10,nnr,n2d)
	
!----------------------------------------------------------------
!  Define "scale_factor" and "add_offset" variables
       real*4     ::  scale4(1), add4(1)
!----------------------------------------------------------------
       allocate(time(nrec)) 
!	 allocate(Temp2d(192,94,nrec))
	 allocate(f2d(192,94,nrec))
!!!!-------------------------set ----------------------------------------------------
!      print*,nrec
!     PRD
!       XG1(1)=60; XG2(1)=64; YG1(1)=37; YG2(1)=34
!       IYG(1)=4; IXG(1)=5
!     MLTR 
!	 XG1(2)=60;XG2(2)=67;YG1(2)=33; YG2(2)=30
!       IYG(2)=4; IXG(2)=8
!     NPC 
!	 XG1(3)=61;XG2(3)=65;YG1(3)=29; YG2(3)=25
!       IYG(3)=5; IXG(3)=5
!     NEC
!       XG1(4)=65;XG2(4)=70;YG1(4)=25; YG2(4)=21
!       IYG(4)=5; IXG(4)=6
! Open netCDF file.
      status=nf_open(trim(filename),nf_nowrite,ncid)
          if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
 

!----------------------------------------------------
!   Retrieve data for Variable 'lat'
!   Units of 'lat' is 'degrees_north'                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
!   Long_name of 'lat' is 'Latitude'
      status=nf_inq_var(ncid,   1,dummy,xtype,ndim,dimids,natts)
          if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      do j=1,ndim
      status=nf_inq_dim(ncid,dimids(j),dummy,len)
           if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      start(j)=1 ; count(j)=len
      end do
      status=nf_get_vara_real(ncid,   1,start,count,lat)

!----------------------------------------------------
!   Retrieve data for Variable 'lon'
!   Units of 'lon' is 'degrees_east'                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
!   Long_name of 'lon' is 'Longitude'
      status=nf_inq_var(ncid,   2,dummy,xtype,ndim,dimids,natts)
          if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      do j=1,ndim
      status=nf_inq_dim(ncid,dimids(j),dummy,len)
           if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      start(j)=1 ; count(j)=len
      end do
      status=nf_get_vara_real(ncid,   2,start,count,lon)

!----------------------------------------------------
!   Retrieve data for Variable 'time'
!   Units of 'time' is 'hours since 1-1-1 00:00:0.0'                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
!   Long_name of 'time' is 'Time'
      status=nf_inq_var(ncid,   3,dummy,xtype,ndim,dimids,natts)
          if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      do j=1,ndim
      status=nf_inq_dim(ncid,dimids(j),dummy,len)
           if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      start(j)=1 ; count(j)=len
      end do
      status=nf_get_vara_double(ncid,   3,start,count,time)

!----------------------------------------------------
!   Retrieve data for Variable 'lhtfl'
!   Long_name of 'lhtfl' is '4xDaily Latent Heat Net Flux.'
!   Units of 'lhtfl' is 'W/m^2'                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
      status=nf_inq_var(ncid,   4,dummy,xtype,ndim,dimids,natts)
          if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      do j=1,ndim
      status=nf_inq_dim(ncid,dimids(j),dummy,len)
           if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      start(j)=1 ; count(j)=len
      end do
      status=nf_get_vara_int2(ncid,   4,start,count,f2d)

      scale4(1) =0.0 ; add4(1) =0.0
! Scale_factor and add_offset for variable 'lhtfl'
      status=nf_get_att_real(ncid,   4,'add_offset',add4(1))
      status=nf_get_att_real(ncid,   4,'scale_factor',scale4(1))

! add scale_factor and add_offset to get true value
! Caution: variable type of 'lhtfl' may not as the same as the type of                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
! 'scale4(1)' and 'add4(1)', you must change it youself!
!------------------------------------------------------------
!      Temp2d=f2d*scale4(1)+add4(1)

! -----------   Some useful advices  --------------
! If dimensions of a variable exceed 3, there can be an error (or warning)
! "warning LNK4084: total image size 382214144 exceeds max (268435456); image may not run"
! when link this program. The best way to resolve it: decrease dimensions
! of the variable, use "do ...  end do" cycle to get little data at one time.
! See NetCDF mannual to look for how to control the dimensions.
!------------   End suggestions   -------------

!----------------------------------------------------
!  Begin writing statements to use the data.


!     Here write your own code please!
      do it=1,nrec
	  do ix=XG1(lo),XG2(lo) ! weat to east
	    do iy=YG1(lo),YG2(lo),-1 ! south to north
	     data2d(ix-XG1(lo)+1,iy-YG2(lo)+1,it+itt,i2d)=
     +		 f2d(ix, iy, it )*scale4(1)+add4(1)
           enddo
		enddo
	 enddo	         
!----------------------------------------------------
!  End Program'
      return
      end subroutine
!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine NC_3D(lo,nrec,itt,i3d,data3d,filename,IXX,IYY,
     +              X1,Y1,X2,Y2,nnr,n3d,lev)
	include 'netcdf.inc'
!
!-------------------------------------------------------
!
!    Using below command line to get the head information of the NetCDF file:
!    Please execute " ncdump -h Z:\DATA\LargeScale\NCEPR2\air.1990.nc > output.txt " 
!-------------------------------------------------------
!
!     Define Variables.
!     Variable ids run sequentially from 1 to nvars=5 ! number of variables
      integer :: nrec,nnr,n3d,lev  ! change this to generalize
      integer*4  ncid, status    ! file control
      integer*4 recdim   ! record dimension
!-------------------------------------------------------------
 !     Below 5 variables is the data in netCDF file
      real*4         ::  level( 17 )
      real*4         ::  lat( 73 )
      real*4         ::  lon( 144 )
      real*8 , allocatable ::  time(:)
      integer*2, allocatable  ::  f3d(:, :,:,: )
!	real, allocatable     ::  temp3d(:,:,:,: )
!     above5 variables is the data in netCDF file
!-------------------------------------------------------------
      integer*4   :: start(10)
      integer*4   :: count(10)
      integer*4   :: dimids(10)! allow up to 10 dimensions
      integer*4   :: dimid, xtype
      character(len=31) :: dummy
!----------------------------------------------------------------
!  Define "scale_factor" and "add_offset" variables
       real*4     ::  scale5(1), add5(1)
!----------------------------------------------------------------
      integer :: itt,i3d	
!	real data2d(15,10,nrecr*10,2)! once max is 10 years
	character*150 filename
	integer IXX(4),IYY(4)
	integer X1(4),Y1(4),X2(4),Y2(4) 
	real data3d(6,6,lev,nnr,n3d) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       allocate(time(nrec)) 
!	 allocate(temp3d(144, 73, 17, nrec ))
       allocate(f3d(144, 73, 17, nrec ))
!-----------------------------------------------------------------
!       print*,nrec,'TTTTT',i3d
! Open netCDF file.
      status=nf_open(trim(filename),
     +	nf_nowrite,ncid)
          if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
 

!----------------------------------------------------
!   Retrieve data for Variable 'level'
!   Units of 'level' is 'millibar'                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
!   Long_name of 'level' is 'Level'
      status=nf_inq_var(ncid,   1,dummy,xtype,ndim,dimids,natts)
          if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      do j=1,ndim
      status=nf_inq_dim(ncid,dimids(j),dummy,len)
           if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      start(j)=1 ; count(j)=len
      end do
      status=nf_get_vara_real(ncid,   1,start,count,level)

!----------------------------------------------------
!   Retrieve data for Variable 'lat'
!   Units of 'lat' is 'degrees_north'                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
!   Long_name of 'lat' is 'Latitude'
      status=nf_inq_var(ncid,   2,dummy,xtype,ndim,dimids,natts)
          if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      do j=1,ndim
      status=nf_inq_dim(ncid,dimids(j),dummy,len)
           if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      start(j)=1 ; count(j)=len
      end do
      status=nf_get_vara_real(ncid,   2,start,count,lat)

!----------------------------------------------------
!   Retrieve data for Variable 'lon'
!   Units of 'lon' is 'degrees_east'                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
!   Long_name of 'lon' is 'Longitude'
      status=nf_inq_var(ncid,   3,dummy,xtype,ndim,dimids,natts)
          if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      do j=1,ndim
      status=nf_inq_dim(ncid,dimids(j),dummy,len)
           if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      start(j)=1 ; count(j)=len
      end do
      status=nf_get_vara_real(ncid,   3,start,count,lon)

!----------------------------------------------------
!   Retrieve data for Variable 'time'
!   Units of 'time' is 'hours since 1800-1-1 00:00:0.0'                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
!   Long_name of 'time' is 'Time'
      status=nf_inq_var(ncid,   4,dummy,xtype,ndim,dimids,natts)
          if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      do j=1,ndim
      status=nf_inq_dim(ncid,dimids(j),dummy,len)
           if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      start(j)=1 ; count(j)=len
      end do
      status=nf_get_vara_double(ncid,   4,start,count,time)

!----------------------------------------------------
!   Retrieve data for Variable 'air'
!   Long_name of 'air' is '6-hourly Air Temperature on Pressure Levels'
!   Units of 'air' is 'degK'                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
      status=nf_inq_var(ncid,   5,dummy,xtype,ndim,dimids,natts)
          if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      do j=1,ndim
      status=nf_inq_dim(ncid,dimids(j),dummy,len)
           if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      start(j)=1 ; count(j)=len
      end do
      status=nf_get_vara_int2(ncid,   5,start,count,f3d)

      scale5(1) =0.0 ; add5(1) =0.0
! Scale_factor and add_offset for variable 'air'
      status=nf_get_att_real(ncid,   5,'add_offset',add5(1))
      status=nf_get_att_real(ncid,   5,'scale_factor',scale5(1))

! add scale_factor and add_offset to get true value
! Caution: variable type of 'air' may not as the same as the type of                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
! 'scale5(1)' and 'add5(1)', you must change it youself!
!------------------------------------------------------------
!      Temp3d=f3d*scale5(1)+add5(1)

! -----------   Some useful advices  --------------
! If dimensions of a variable exceed 3, there can be an error (or warning)
! "warning LNK4084: total image size 382214144 exceeds max (268435456); image may not run"
! when link this program. The best way to resolve it: decrease dimensions
! of the variable, use "do ...  end do" cycle to get little data at one time.
! See NetCDF mannual to look for how to control the dimensions.
!------------   End suggestions   -------------

!----------------------------------------------------
!  Begin writing statements to use the data.


!     Here write your own code please!
       do it=1,nrec
	   do ix=X1(lo),X2(lo) ! weat to east
	     do iy=Y1(lo),Y2(lo),-1 ! south to north
	       do ik=1,17
	       data3d(ix-X1(lo)+1,iy-Y2(lo)+1,ik,it+itt,i3d)=
     +		 f3d(ix, iy, ik,it )*scale5(1)+add5(1)
!		    if(i3d==4)then
!	        print*,data3d(ix-X1(lo)+1,iy-Y2(lo)+1,ik,it+itt,i3d)
!               endif
	       enddo
		 enddo
	   enddo
	 enddo	 

!	 print*,i3d        
!----------------------------------------------------
!  End Program'
      return
      end subroutine


	subroutine NC_SRF(lo,nrec,itt,i2d,data2d,filename,IXX,IYY,
     +              X1,Y1,X2,Y2,nnr,n2d)
      include 'netcdf.inc'
!
!-------------------------------------------------------
!
!    Using below command line to get the head information of the NetCDF file:
!    Please execute " ncdump -h Z:\DATA\LargeScale\NCEP_Surface\air.sig995.1990.nc > output.txt " 
!-------------------------------------------------------
!
!     Define Variables.
!     Variable ids run sequentially from 1 to nvars=4 ! number of variables
!
      integer :: nrec,n2d  ! change this to generalize
      integer*4  ncid, status    ! file control
      integer*4 recdim   ! record dimension
!-------------------------------------------------------------
 !     Below 4 variables is the data in netCDF file
      real*4         ::  lat( 73 )
      real*4         ::  lon( 144 )
      real*8, allocatable  ::  time(:)
      integer*2, allocatable  ::  fs( :, :, : )
!     above4 variables is the data in netCDF file
!-------------------------------------------------------------
      integer*4   :: start(10)
      integer*4   :: count(10)
      integer*4   :: dimids(10)! allow up to 10 dimensions
      integer*4   :: dimid, xtype
      character(len=31) :: dummy
!----------------------------------------------------------------
!  Define "scale_factor" and "add_offset" variables
       real*4     ::  scale4(1), add4(1)
	 integer :: itt,i3d	
	real data2d(10,10,nnr,n2d)! once max is 10 years
	character*150 filename
	integer IXX(4),IYY(4)
	integer X1(4),Y1(4),X2(4),Y2(4) 
	allocate(time(nrec))
	allocate(fs(144, 73, nrec )) 
!----------------------------------------------------------------
!      X1(1)=45; X2(1)=48; Y1(1)=29; Y2(1)=27;IXX(1)=4
!       IYY(1)=3
!     MLTR 
!	 X1(2)=45; X2(2)=50 ;Y1(2)=26; Y2(2)=24; IXX(2)=6
!	 IYY(2)=3
!     NPC 
!	X1(3)=46; X2(3)=49;Y1(3)=23; Y2(3)=20; IXX(3)=4
!	IYY(3)=4
!     NEC
!      X1(4)=49; X2(4)=53; Y1(4)=20; Y2(4)=17; IXX(4)=5
!      IYY(4)=4
! Open netCDF file.
      status=nf_open(trim(filename),
     +	nf_nowrite,ncid)
          if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
 

!----------------------------------------------------
!   Retrieve data for Variable 'lat'
!   Units of 'lat' is 'degrees_north'                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
!   Long_name of 'lat' is 'Latitude'
      status=nf_inq_var(ncid,   1,dummy,xtype,ndim,dimids,natts)
          if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      do j=1,ndim
      status=nf_inq_dim(ncid,dimids(j),dummy,len)
           if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      start(j)=1 ; count(j)=len
      end do
      status=nf_get_vara_real(ncid,   1,start,count,lat)

!----------------------------------------------------
!   Retrieve data for Variable 'lon'
!   Units of 'lon' is 'degrees_east'                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
!   Long_name of 'lon' is 'Longitude'
      status=nf_inq_var(ncid,   2,dummy,xtype,ndim,dimids,natts)
          if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      do j=1,ndim
      status=nf_inq_dim(ncid,dimids(j),dummy,len)
           if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      start(j)=1 ; count(j)=len
      end do
      status=nf_get_vara_real(ncid,   2,start,count,lon)

!----------------------------------------------------
!   Retrieve data for Variable 'time'
!   Units of 'time' is 'hours since 1-1-1 00:00:0.0'                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
!   Long_name of 'time' is 'Time'
      status=nf_inq_var(ncid,   3,dummy,xtype,ndim,dimids,natts)
          if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      do j=1,ndim
      status=nf_inq_dim(ncid,dimids(j),dummy,len)
           if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      start(j)=1 ; count(j)=len
      end do
      status=nf_get_vara_double(ncid,   3,start,count,time)

!----------------------------------------------------
!   Retrieve data for Variable 'air'
!   Long_name of 'air' is '4xDaily Air temperature.'
!   Units of 'air' is 'degK'                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
      status=nf_inq_var(ncid,   4,dummy,xtype,ndim,dimids,natts)
          if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      do j=1,ndim
      status=nf_inq_dim(ncid,dimids(j),dummy,len)
           if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      start(j)=1 ; count(j)=len
      end do
      status=nf_get_vara_int2(ncid,   4,start,count,fs)

      scale4(1) =0.0 ; add4(1) =0.0
! Scale_factor and add_offset for variable 'air'
      status=nf_get_att_real(ncid,   4,'add_offset',add4(1))
      status=nf_get_att_real(ncid,   4,'scale_factor',scale4(1))

! add scale_factor and add_offset to get true value
! Caution: variable type of 'air' may not as the same as the type of                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
! 'scale4(1)' and 'add4(1)', you must change it youself!
!
!     air=air*scale4(1)+add4(1)

! -----------   Some useful advices  --------------
! If dimensions of a variable exceed 3, there can be an error (or warning)
! "warning LNK4084: total image size 382214144 exceeds max (268435456); image may not run"
! when link this program. The best way to resolve it: decrease dimensions
! of the variable, use "do ...  end do" cycle to get little data at one time.
! See NetCDF mannual to look for how to control the dimensions.
!------------   End suggestions   -------------

!----------------------------------------------------
!  Begin writing statements to use the data.


!     Here write your own code please!
       do it=1,nrec
	   do ix=X1(lo),X2(lo) ! weat to east
	     do iy=Y1(lo),Y2(lo),-1 ! south to north
	       data2d(ix-X1(lo)+1,iy-Y2(lo)+1,it+itt,i2d)=
     +		 fs(ix, iy,it )*scale4(1)+add4(1)
!	     if(i2d==6)print*,data2d(ix-X1(lo)+1,iy-Y2(lo)+1,it+itt,i2d)
!		    if(i3d==4)then
!	        print*,data3d(ix-X1(lo)+1,iy-Y2(lo)+1,ik,it+itt,i3d)
!               endif
		 enddo
	   enddo
	 enddo

!----------------------------------------------------
!  End Program
       return
       end subroutine
      
	subroutine NC_HGT(lo,hgt,filename,IXX,IYY,
     +              X1,Y1,X2,Y2)
!  Do not forget to include the -I path_to_netcdf_  includes in your compile statement Required includes.
!  Also note: need 'netcdf.lib' or 'netcdfs.lib' when link .
      include 'netcdf.inc'
!
!-------------------------------------------------------
!
!    Using below command line to get the head information of the NetCDF file:
!    Please execute " ncdump -h Z:\DATA\LargeScale\NCEP_Surface\hgt.sfc.nc > output.txt " 
!-------------------------------------------------------
!
!     Define Variables.
!     Variable ids run sequentially from 1 to nvars=4 ! number of variables
      integer,parameter :: nrec=1  ! change this to generalize
      integer*4  ncid, status    ! file control
      integer*4 recdim   ! record dimension
!-------------------------------------------------------------
 !     Below 4 variables is the data in netCDF file
      real*4         ::  lat( 73 )
      real*4         ::  lon( 144 )
      real*8         ::  time(nrec)
      integer*2      ::  heigh( 144, 73, nrec )
!     above4 variables is the data in netCDF file
!-------------------------------------------------------------
      integer*4   :: start(10)
      integer*4   :: count(10)
      integer*4   :: dimids(10)! allow up to 10 dimensions
      integer*4   :: dimid, xtype
      character(len=31) :: dummy
!----------------------------------------------------------------
!  Define "scale_factor" and "add_offset" variables
       real*4     ::  scale4(1), add4(1)
	real hgt(10,10)! once max is 10 years
	character*150 filename
	integer IXX(4),IYY(4),lo
	integer X1(4),Y1(4),X2(4),Y2(4) 
!----------------------------------------------------------------
!----------------------------------------------------------------

! Open netCDF file.
      status=nf_open(trim(filename),
     +	nf_nowrite,ncid)
          if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
 

!----------------------------------------------------
!   Retrieve data for Variable 'lat'
!   Units of 'lat' is 'degrees_north'                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
!   Long_name of 'lat' is 'Latitude'
      status=nf_inq_var(ncid,   1,dummy,xtype,ndim,dimids,natts)
          if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      do j=1,ndim
      status=nf_inq_dim(ncid,dimids(j),dummy,len)
           if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      start(j)=1 ; count(j)=len
      end do
      status=nf_get_vara_real(ncid,   1,start,count,lat)

!----------------------------------------------------
!   Retrieve data for Variable 'lon'
!   Units of 'lon' is 'degrees_east'                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
!   Long_name of 'lon' is 'Longitude'
      status=nf_inq_var(ncid,   2,dummy,xtype,ndim,dimids,natts)
          if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      do j=1,ndim
      status=nf_inq_dim(ncid,dimids(j),dummy,len)
           if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      start(j)=1 ; count(j)=len
      end do
      status=nf_get_vara_real(ncid,   2,start,count,lon)

!----------------------------------------------------
!   Retrieve data for Variable 'time'
!   Units of 'time' is 'hours since 1800-1-1 00:00:0.0'                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
!   Long_name of 'time' is 'Time'
      status=nf_inq_var(ncid,   3,dummy,xtype,ndim,dimids,natts)
          if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      do j=1,ndim
      status=nf_inq_dim(ncid,dimids(j),dummy,len)
           if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      start(j)=1 ; count(j)=len
      end do
      status=nf_get_vara_double(ncid,   3,start,count,time)

!----------------------------------------------------
!   Retrieve data for Variable 'hgt'
!   Long_name of 'hgt' is 'Geopotential Height at the Surface'
!   Units of 'hgt' is 'm'                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
      status=nf_inq_var(ncid,   4,dummy,xtype,ndim,dimids,natts)
          if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      do j=1,ndim
      status=nf_inq_dim(ncid,dimids(j),dummy,len)
           if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      start(j)=1 ; count(j)=len
      end do
      status=nf_get_vara_int2(ncid,   4,start,count,heigh)

      scale4(1) =0.0 ; add4(1) =0.0
! Scale_factor and add_offset for variable 'hgt'
      status=nf_get_att_real(ncid,   4,'add_offset',add4(1))
      status=nf_get_att_real(ncid,   4,'scale_factor',scale4(1))

! add scale_factor and add_offset to get true value
! Caution: variable type of 'hgt' may not as the same as the type of                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
! 'scale4(1)' and 'add4(1)', you must change it youself!
  !    hgt=hgt*scale4(1)+add4(1)

! -----------   Some useful advices  --------------
! If dimensions of a variable exceed 3, there can be an error (or warning)
! "warning LNK4084: total image size 382214144 exceeds max (268435456); image may not run"
! when link this program. The best way to resolve it: decrease dimensions
! of the variable, use "do ...  end do" cycle to get little data at one time.
! See NetCDF mannual to look for how to control the dimensions.
!------------   End suggestions   -------------

!----------------------------------------------------
!  Begin writing statements to use the data.


!     Here write your own code please!
        do it=1,nrec
	   do ix=X1(lo),X2(lo) ! weat to east
	     do iy=Y1(lo),Y2(lo),-1 ! south to north
	       hgt(ix-X1(lo)+1,iy-Y2(lo)+1)=
     +		 heigh(ix, iy ,it )*scale4(1)+add4(1)
!		    if(i3d==4)then
!	        print*,data3d(ix-X1(lo)+1,iy-Y2(lo)+1,ik,it+itt,i3d)
!               endif
		 enddo
	   enddo
       enddo
!----------------------------------------------------
!  End Program
       return
       end subroutine

	subroutine t2msst(lo,nrec,itt,i2d,data2d,filename,IXG,IYG,
     +              XG1,YG1,XG2,YG2,nnr,n2d)
      include 'netcdf.inc'
!
!-------------------------------------------------------
!
!    Using below command line to get the head information of the NetCDF file:
!    Please execute " ncdump -h Z:\DATA\LargeScale\NCEP_Surface\lhtfl.sfc.gauss.1990.nc > output.txt " 
!-------------------------------------------------------
!
!     Define Variables.
!     Variable ids run sequentially from 1 to nvars=4 ! number of variables
      integer :: nrec,nnr,n2d ! change this to generalize
      integer*4  ncid, status    ! file control
      integer*4 recdim   ! record dimension
!-------------------------------------------------------------
 !     Below 4 variables is the data in netCDF file
      real*4         ::  level( 1 )
      real*4         ::  lat( 94 )
      real*4         ::  lon( 192 )
      real*8 , allocatable ::  time(:)
!      integer*2, allocatable  ::  f2d( :, :, : )
	integer*2, allocatable  ::    air( :, :, :, : )
!	real, allocatable :: Temp2d( :, :, : )
!     above4 variables is the data in netCDF file
!-------------------------------------------------------------
      integer*4   :: start(10)
      integer*4   :: count(10)
      integer*4   :: dimids(10)! allow up to 10 dimensions
      integer*4   :: dimid, xtype
      character(len=31) :: dummy
!-------------addd--------------------------------------------------
      integer :: itt,i2d	
!	real data2d(15,10,nrecr*10,2)! once max is 10 years
	character*150 filename
	integer IXG(4),IYG(4)
	integer XG1(4),YG1(4),XG2(4),YG2(4) 
	real data2d(10,10,nnr,n2d)
	
!----------------------------------------------------------------
!  Define "scale_factor" and "add_offset" variables
!       real*4     ::  scale4(1), add4(1)
	 real*4     ::  scale5(1), add5(1)
!----------------------------------------------------------------
       allocate(time(nrec)) 
!	 allocate(Temp2d(192,94,nrec))
	 allocate(air(192,94,1,nrec))  !, 
!!!!-------------------------set ----------------------------------------------------
!      print*,nrec
!     PRD
!       XG1(1)=60; XG2(1)=64; YG1(1)=37; YG2(1)=34
!       IYG(1)=4; IXG(1)=5
!     MLTR 
!	 XG1(2)=60;XG2(2)=67;YG1(2)=33; YG2(2)=30
!       IYG(2)=4; IXG(2)=8
!     NPC 
!	 XG1(3)=61;XG2(3)=65;YG1(3)=29; YG2(3)=25
!       IYG(3)=5; IXG(3)=5
!     NEC
!       XG1(4)=65;XG2(4)=70;YG1(4)=25; YG2(4)=21
!       IYG(4)=5; IXG(4)=6
! Open netCDF file.
      status=nf_open(trim(filename),nf_nowrite,ncid)
C      status=nf_open('Z:\DATA\ReadNC\nctool\air.2m.gauss.1996.nc',nf_nowrite,ncid)
          if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
 

!----------------------------------------------------
!   Retrieve data for Variable 'level'
!   Units of 'level' is 'm'                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
!   Long_name of 'level' is 'Level'
      status=nf_inq_var(ncid,   1,dummy,xtype,ndim,dimids,natts)
          if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      do j=1,ndim
      status=nf_inq_dim(ncid,dimids(j),dummy,len)
           if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      start(j)=1 ; count(j)=len
      end do
      status=nf_get_vara_real(ncid,   1,start,count,level)

!----------------------------------------------------
!   Retrieve data for Variable 'lat'
!   Units of 'lat' is 'degrees_north'                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
!   Long_name of 'lat' is 'Latitude'
      status=nf_inq_var(ncid,   2,dummy,xtype,ndim,dimids,natts)
          if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      do j=1,ndim
      status=nf_inq_dim(ncid,dimids(j),dummy,len)
           if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      start(j)=1 ; count(j)=len
      end do
      status=nf_get_vara_real(ncid,   2,start,count,lat)

!----------------------------------------------------
!   Retrieve data for Variable 'lon'
!   Units of 'lon' is 'degrees_east'                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
!   Long_name of 'lon' is 'Longitude'
      status=nf_inq_var(ncid,   3,dummy,xtype,ndim,dimids,natts)
          if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      do j=1,ndim
      status=nf_inq_dim(ncid,dimids(j),dummy,len)
           if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      start(j)=1 ; count(j)=len
      end do
      status=nf_get_vara_real(ncid,   3,start,count,lon)

!----------------------------------------------------
!   Retrieve data for Variable 'time'
!   Units of 'time' is 'hours since 1800-1-1 00:00:0.0'                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
!   Long_name of 'time' is 'Time'
      status=nf_inq_var(ncid,   4,dummy,xtype,ndim,dimids,natts)
          if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      do j=1,ndim
      status=nf_inq_dim(ncid,dimids(j),dummy,len)
           if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      start(j)=1 ; count(j)=len
      end do
      status=nf_get_vara_double(ncid,   4,start,count,time)

!----------------------------------------------------
!   Retrieve data for Variable 'air'
!   Long_name of 'air' is '6-Hourly Forecast of Air temperature at 2 m'
!   Units of 'air' is 'degK'                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
      status=nf_inq_var(ncid,   5,dummy,xtype,ndim,dimids,natts)
          if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      do j=1,ndim
      status=nf_inq_dim(ncid,dimids(j),dummy,len)
           if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      start(j)=1 ; count(j)=len
      end do
      status=nf_get_vara_int2(ncid,   5,start,count,air)

      scale5(1) =0.0 ; add5(1) =0.0
! Scale_factor and add_offset for variable 'air'
      status=nf_get_att_real(ncid,   5,'add_offset',add5(1))
      status=nf_get_att_real(ncid,   5,'scale_factor',scale5(1))

! add scale_factor and add_offset to get true value
! Caution: variable type of 'air' may not as the same as the type of                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
! 'scale5(1)' and 'add5(1)', you must change it youself!

!     air=air*scale5(1)+add5(1)

! -----------   Some useful advices  --------------
! If dimensions of a variable exceed 3, there can be an error (or warning)
! "warning LNK4084: total image size 382214144 exceeds max (268435456); image may not run"
! when link this program. The best way to resolve it: decrease dimensions
! of the variable, use "do ...  end do" cycle to get little data at one time.
! See NetCDF mannual to look for how to control the dimensions.
!------------   End suggestions   -------------

!----------------------------------------------------
!  Begin writing statements to use the data.


!     Here write your own code please!

      do it=1,nrec
	  do ix=XG1(lo),XG2(lo) ! weat to east
	    do iy=YG1(lo),YG2(lo),-1 ! south to north
	     data2d(ix-XG1(lo)+1,iy-YG2(lo)+1,it+itt,i2d)=
     +		 air(ix, iy,1, it )*scale5(1)+add5(1)
           enddo
		enddo
	 enddo
!----------------------------------------------------
!  End Program
        end subroutine
