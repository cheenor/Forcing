      PROGRAM ES_Indexs
!********************************************************************************

!*******************************************************************************
      PARAMETER(NT=34,nrec=408)   !!!M time steps   N stations
      PARAMETER(MNH=NT,KS=0,KV=10,KVT=10)
      PARAMETER(ff=-9.99E+08)
      DIMENSION A(MNH,MNH),S(MNH,MNH),ER(mnh,4),V(MNH),tCF(NT,KVT),north(mnh)
!	  REAL, allocatable,dimension(:,:) :: F,evf
!	  REAL, allocatable,dimension(:) :: DF,AVF
!****  INPUT DATA ***********************************************
!      real       ::  temp( 144, 73, 17, nrec ),slptmp( 144, 73, nrec )
	  real       ::  gpctmp(144,72,12),ee
	  real       ::  uwnd( 144, 73, 17, nrec ),vwnd( 144, 73, 17,nrec )
	  real       ::  slp( 144, 73, nrec ),gpc(7,NT)

!      REAL, allocatable,dimension(:,:,:)::v850,v200,u850,u200,slp,gpc
!      REAL, allocatable,dimension(:,:,:)::v850_nor,v200_nor,   &
!    &                      	u850_nor,u200_nor,slp_nor,gpc_nor,PT
	  INTEGER :: NN,MM,KK,iiii,NX,NY,lat1,lat2,lon1,lon2
	  INTEGER :: iy ,inor
	  CHARACTER(LEN=100) :: filename,dir,dirGPCP,strGPCP
	  CHARACTER(LEN=50) :: filestr(3)
      CHARACTER(LEN=4) :: yrstr
	  character*4 :: area(7)
      integer :: lon(7,2), lat(7,2)
!!!!-------------------------set ----------------------------------------------------
  area(1)='PRD'  !!! 110 117.5 / 20 27.5
  lon(1,1)=45 ;lon(1,2)=48 
  lat(1,1)=45 ;lat(1,2)=48   
  area(2)='MLYR' !!! 110 117.5/ 27.5 35
  lon(2,1)=45 ;lon(2,2)=48 
  lat(2,1)=48 ;lat(2,2)=51 
  area(3)='NPC'  !!!! 112.5 120 /35 42.5
  lon(3,1)=46 ;lon(3,2)=49 
  lat(3,1)=51 ;lat(3,2)=54
  area(4)='NEC'  !!! 120 127.5/42.5 50
  lon(4,1)=49 ;lon(4,2)=52 
  lat(4,1)=54 ;lat(4,2)=57
    area(5)='BOB'  !!! 85-95/ 15 25
  lon(5,1)=34 ;lon(5,2)=38 
  lat(5,1)=42 ;lat(5,2)=47
    area(6)='ETP'  !!! 90-102.5/27.5 37.5 
  lon(6,1)=36 ;lon(6,2)=42 
  lat(6,1)=47 ;lat(6,2)=52
    area(7)='WTP'  !!!  77.5-90 /27.5 37.5 
  lon(7,1)=31 ;lon(7,2)=37 
  lat(7,1)=47;lat(7,2)=52

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! nc files for 850 200 u v and the slp
      filestr(1)='uwnd.mon.mean.nc'
	filestr(2)='vwnd.mon.mean.nc'
	filestr(3)='mslp.mon.mean.nc'
!-------GPCP files----------------------------------------
      dirGPCP='Z:\DATA\GPCP\MON\'
      strGPCP='gpcp_v2.2_psg.'
	dir='Z:\DATA\EOF\Data\'
!-----------uwind 
      filename=trim(dir)//trim(filestr(1))   
      call met(filename,uwnd)
!	  print*,uwnd(10,10,10,2)
!-----------vwind 
      filename=trim(dir)//trim(filestr(2))   
      call met(filename,vwnd)
!	   print*,vwnd(10,10,10,2)
!-----------slp
      filename=trim(dir)//trim(filestr(3))   
      call slprd(filename,slp)
	  print*,slp(44,16,1),'**'

	  NN=size(slp,dim=1);MM=size(slp,dim=2);KK=size(slp,dim=3)
	  print*,NN,MM,KK
      Call ZZIndex(slp,NN,MM,KK)
	  NN=size(uwnd,dim=1);MM=size(uwnd,dim=2);KK=size(uwnd,dim=3);KT=size(uwnd,dim=4)
      Call WYIndex(uwnd,NN,MM,KK,KT)
      NN=size(uwnd,dim=1);MM=size(uwnd,dim=2);KK=size(uwnd,dim=3);KT=size(uwnd,dim=4)
      Call WFIndex(uwnd,NN,MM,KK,KT)
	  NN=size(uwnd,dim=1);MM=size(uwnd,dim=2);KK=size(uwnd,dim=3);KT=size(uwnd,dim=4)
	  Call WHIndex(uwnd,vwnd,NN,MM,KK,KT)
	  NN=size(uwnd,dim=1);MM=size(uwnd,dim=2);KK=size(uwnd,dim=3);KT=size(uwnd,dim=4)
	  Call LCIndex(vwnd,NN,MM,KK,KT)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!----------GPCP---------------------------------------------------------
      do iy=1979,2012
	  do ir=1,7
        write(yrstr,'(I4)')iy
	    filename=trim(dirGPCP)//trim(strGPCP)//yrstr
	    call GPCP(filename, gpctmp)
!!!!!!! -------------test----------------------
        if(iy==2000)then
		  print*,' grid_Id       Grads     Read'
		  print*,'X45 Y50 T1',' 0.689906 ',gpctmp(45,72-50+1,1)
		  print*,'X45 Y50 T1',' 0.759595 ',gpctmp(45,72-50+1,4)
		  print*,'X45 Y50 T1',' 6.1463 ',gpctmp(45,72-50+1,7)
		endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		is=lon(ir,1);ie=lon(ir,2)
		js=lat(ir,1);je=lat(ir,2)
       do i=is,ie
	   do j=72-js+1,72-je+1,-1
           gpc(ir,iy-1978)=gpc(ir,iy-1978)+(gpctmp(i,j,6)+gpctmp(i,j,7)+gpctmp(i,j,8))/3.0
	   enddo
	   enddo
	     gpc(ir,iy-1978)=gpc(ir,iy-1978)/((ie-is+1)*(je-js+1)*1.0)
      enddo
	enddo
	!!!!!!!!!!!!
	call QIYI(gpc,area)
	!!!!!!!!!!!!!!!
	  do ir=1,7
	    filename='Z:\DATA\LargeScale\79-12\GPCP\'//trim(area(ir))//'.txt'
		write(96,*)'year rainfall'
		open(96,file=trim(filename))
     	do iy=1979,2012
         write(96,99)iy,gpc(ir,iy-1978)
	  enddo
	  close(96)
	  enddo

99 format(1X,I4,1X,F8.4)
    end
! 
Subroutine ZZIndex(slp,nx,ny,nt)
REAL, dimension(nx,ny,nt) :: slp,slpnor
! input 2.5 degree ncep R2 Monthly sea level presure(slp)
! Zhao and Zhou 2005  monsoon index
! East Asian subtropical summer monsoon index and its relationships to rainfall.
! Acta Meteor Sinica 63:933-941
!ALLOCATE(slp(nx,ny,nt)) 
!ALLOCATE(slpnor(nx,ny,nt)) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  normalization
!-------------------------------------------------------------------
NN=size(slp,dim=1);MM=size(slp,dim=2);KK=size(slp,dim=3)
!print*,NN,MM,KK,slp(44,44,3)
do i=1,nt
call normalization(slp(:,:,i),NN,MM,slpnor(:,:,i)) 
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
open(10,file='Z:\DATA\LargeScale\79-12\INdexs\ZZI.txt')
write(10,*)'Zhao and Zhou 2005 Acta Meteor Sinica 63:933-941'
write(10,*)'Year ZZI'
iyr=1979
do it=6,nt,12
  psib=0.0 ;  psub=0.0
  do iy=17,21 !! 50N-40N
  psib=psib+slpnor(45,iy,it)+slpnor(45,iy,it+1)+slpnor(45,iy,it+2)   !!! nx=44 110E
!  print*,slpnor(44,iy,it),slp(44,iy,it)
  enddo
  do iy=21,25 !! 40N-30N
  psub=psub+slpnor(65,iy,it)+slpnor(65,iy,it+1)+slpnor(65,iy,it+2)   !!! nx=64 160E
  enddo
  zzi=(psub-psib)/(5.0*3)
  write(10,99)iyr,zzi
  iyr=iyr+1
enddo
close(10)
99 format(1X,I4,1X,F8.4)
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine WYIndex(uwnd,nx,ny,nz,nt)
REAL, dimension(nx,ny,nz,nt) :: uwnd,uwndnor

!ALLOCATE(uwnd(nx,ny,nz,nt))
!ALLOCATE(uwndnor(nx,ny,nz,nt))
open(10,file='Z:\DATA\LargeScale\79-12\INdexs\WYI.txt')
write(10,*)'Webster and Yang (1992)'
write(10,*)'Year WYI'
iyr=1979
do it=6,nt,12
u850=0.0
u200=0.0
   do iy=21,33 !!! 40 -10 N
   do ix=45,57  !!!! 110-140E
   u850=u850+(uwnd(ix,iy,3,it)+uwnd(ix,iy,3,it+1) +uwnd(ix,iy,3,it+2))/3.0
   u200=u200+(uwnd(ix,iy,9,it)+uwnd(ix,iy,9,it+1) +uwnd(ix,iy,9,it+2))/3.0
   enddo
   enddo
  WYI=(u850-u200)/((32-20+1)*(56-44+1)*1.0)
   write(10,99)iyr,WYI
  iyr=iyr+1
enddo
close(10)
99 format(1X,I4,1X, F8.4)
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine WFIndex(uwnd,nx,ny,nz,nt)
REAL, dimension(nx,ny,nz,nt) :: uwnd,uwndnor

!ALLOCATE(uwnd(nx,ny,nz,nt))
!ALLOCATE(uwndnor(nx,ny,nz,nt))
open(10,file='Z:\DATA\LargeScale\79-12\INdexs\WFI.txt')
write(10,*)'Wang and Fan (1999)'
write(10,*)'Year WFI'
iyr=1979
do it=6,nt,12
u850_1=0.0
u850_2=0.0 
   do iy=31,35 !!!  15 5N
   do ix=37,53  !!!! 90 130E
   u850_1=u850_1+(uwnd(ix,iy,3,it)+uwnd(ix,iy,3,it+1) +uwnd(ix,iy,3,it+2))/3.0
   enddo
   enddo
   do iy=24,28 !!!  32.5 -22.5
   do ix=45,57  !!!! 110 140E
   u850_2=u850_2+(uwnd(ix,iy,3,it)+uwnd(ix,iy,3,it+1) +uwnd(ix,iy,3,it+2))/3.0
   enddo
   enddo
    u850_1=u850_1/((34-30+1)*(52-36+1)*1.0)
	u850_2=u850_2/((27-23+1)*(56-44+1)*1.0)
  WFI=u850_1-u850_2
   write(10,99)iyr,WFI
  iyr=iyr+1
enddo
close(10)
99 format(1X,I4,1X, F8.4)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine WHIndex(uwnd,vwnd,nx,ny,nz,nt)
REAL, dimension(nx,ny,nz,nt) :: uwnd,vwnd

!ALLOCATE(uwnd(nx,ny,nz,nt))
!ALLOCATE(vwnd(nx,ny,nz,nt))
open(10,file='Z:\DATA\LargeScale\79-12\INdexs\WHI.txt')
write(10,*)'Wang (2002)'
write(10,*)'Year WHI'
iyr=1979
do it=6,nt,12
   do iy=21,29 !!!  40,20N
   do ix=45,51  !!!! 110 125E
   u850=u850+(uwnd(ix,iy,3,it)+uwnd(ix,iy,3,it+1) +uwnd(ix,iy,3,it+2))/3.0
   v850=v850+(vwnd(ix,iy,3,it)+vwnd(ix,iy,3,it+1) +vwnd(ix,iy,3,it+2))/3.0
   enddo
   enddo
enddo
   u850=u850/((28-20+1)*(50-44+1)*(nt/12.0)*1.0)
   v850=v850/((28-20+1)*(50-44+1)*(nt/12.0)*1.0)
   mean=sqrt(u850**2+v850**2)
do it=6,nt,12
u850=0.0
v850=0.0 
   do iy=21,29 !!!  40,20N
   do ix=45,51  !!!! 110 125E
   u850=u850+(uwnd(ix,iy,3,it)+uwnd(ix,iy,3,it+1) +uwnd(ix,iy,3,it+2))/3.0
   v850=v850+(vwnd(ix,iy,3,it)+vwnd(ix,iy,3,it+1) +vwnd(ix,iy,3,it+2))/3.0
   enddo
   enddo
   u850=u850/((28-20+1)*(50-44+1)*1.0)
   v850=v850/((28-20+1)*(50-44+1)*1.0)
  WHI=sqrt(u850**2+v850**2)-mean
   write(10,99)iyr,WHI
  iyr=iyr+1
enddo


close(10)
99 format(1X,I4,1X, F8.4)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine LCIndex(vwnd,nx,ny,nz,nt)
REAL, dimension(nx,ny,nz,nt) :: vwnd
real :: lci

!ALLOCATE(vwnd(nx,ny,nz,nt))
!ALLOCATE(vwndnor(nx,ny,nz,nt))
open(10,file='Z:\DATA\LargeScale\79-12\INdexs\LCI.txt')
write(10,*)'Lu and Chan (1999)'
write(10,*)'Year LCI'
iyr=1979
! print*,vwnd(10,10,10,2),'LCI'
do it=6,nt,12
v000=0.0
   do iy=29,34 !!!  20 7.5N
   do ix=44,49  !!!! 107.5 120E
   v000=v000+(vwnd(ix,iy,1,it)+vwnd(ix,iy,1,it+1) +vwnd(ix,iy,1,it+2))/3.0
!   print*,v000
   enddo
   enddo
  LCI=v000/36.0
!  print*,lci,v000,36
   write(10,99)iyr,LCI
  iyr=iyr+1
enddo
close(10)
99 format(1X,I4,1X, F8.4)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!--------------------------------------------------------------------
      SUBROUTINE MET(filename,temp)
!!---------------------------------------------------------------------------------------
!-------------------------------------------------------
!  
!  uvmon.f
!    This file is a fortran template file designed to read the given
!     netCDF file 'uwnd.mon.mean.nc' into memory.
!  
!  History: 
!  Date       Name          Action
!-------------------------------------------------------
!  
!  ?? Oct 1993  B. Schwartz   Created.
!  24 Aug 1994  E. Boer       Modified dimension output.
!  29 Jul 1997  S. Rupert     Standardized and reinstated required include of netcdf.inc. 
!  30 Jul 1997  S. Rupert     Added usage message and command line inputs.
!  03 Apr 2003  H. Yan        Change output to fortran 90 format, modified maximum dimesions to 5.
!  30 Mar 2004  H. Yan        Added usage message, retrieve "scale_factor" and "add_offset" parameters.
!  18 May 2006  H. Yan        Fixed bugs of dimensions definition and remark lines. 
!-------------------------------------------------------
!  Do not forget to include the -I path_to_netcdf_  includes in your compile statement Required includes.
!  Also note: need 'netcdf.lib' or 'netcdfs.lib' when link .
      include 'netcdf.inc'
!
!-------------------------------------------------------
!
!    Using below command line to get the head information of the NetCDF file:
!    Please execute " ncdump -h uwnd.mon.mean.nc > output.txt " 
!-------------------------------------------------------
!
!     Define Variables.
!     Variable ids run sequentially from 1 to nvars=6 ! number of variables
      integer,parameter :: nrec=408  ! change this to generalize
      integer*4  ncid, status    ! file control
      integer*4 recdim   ! record dimension
!-------------------------------------------------------------
 !     Below 6 variables is the data in netCDF file
      real*4         ::  level( 17 )
      real*4         ::  lat( 73 )
      real*4         ::  lon( 144 )
      real*8         ::  time(nrec)
      real*8         ::  time_bnds( 2, 408 )
      integer*2    uwnd( 144, 73, 17, nrec )
!     above6 variables is the data in netCDF file
!-------------------------------------------------------------
      integer*4   :: start(10)
      integer*4   :: count(10)
      integer*4   :: dimids(10)! allow up to 10 dimensions
      integer*4   :: dimid, xtype
      character(len=31) :: dummy
	character(len=100) :: filename
!----------------------------------------------------------------
!  Define "scale_factor" and "add_offset" variables
       real*4     ::  scale6(1), add6(1)
	 real       ::  temp( 144, 73, 17, nrec )
!----------------------------------------------------------------

! Open netCDF file.
      status=nf_open(trim(filename),nf_nowrite,ncid)
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
!   Units of 'time' is 'hours since 1800-1-1 00:00:00'                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
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
!   Retrieve data for Variable 'time_bnds'
!   Long_name of 'time_bnds' is 'Time Boundaries'
      status=nf_inq_var(ncid,   5,dummy,xtype,ndim,dimids,natts)
          if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      do j=1,ndim
      status=nf_inq_dim(ncid,dimids(j),dummy,len)
           if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      start(j)=1 ; count(j)=len
      end do
      status=nf_get_vara_double(ncid,   5,start,count,time_bnds)

!----------------------------------------------------
!   Retrieve data for Variable 'uwnd'
!   Long_name of 'uwnd' is 'Monthly U-wind on Pressure Levels'
!   Units of 'uwnd' is 'm/s'                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
      status=nf_inq_var(ncid,   6,dummy,xtype,ndim,dimids,natts)
          if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      do j=1,ndim
      status=nf_inq_dim(ncid,dimids(j),dummy,len)
           if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      start(j)=1 ; count(j)=len
      end do
      status=nf_get_vara_int2(ncid,   6,start,count,uwnd)

      scale6(1) =0.0 ; add6(1) =0.0
! Scale_factor and add_offset for variable 'uwnd'
      status=nf_get_att_real(ncid,   6,'add_offset',add6(1))
      status=nf_get_att_real(ncid,   6,'scale_factor',scale6(1))

! add scale_factor and add_offset to get true value
! Caution: variable type of 'uwnd' may not as the same as the type of                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
! 'scale6(1)' and 'add6(1)', you must change it youself!
      temp=uwnd*scale6(1)+add6(1)

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


!----------------------------------------------------
!  End Program
      return
      END SUBROUTINE

	SUBROUTINE slprd(filename,slptmp)
!-------------------------------------------------------
!  
!  mslpmon.f
!    This file is a fortran template file designed to read the given
!     netCDF file 'mslp.mon.mean.nc' into memory.
!  
!  History: 
!  Date       Name          Action
!-------------------------------------------------------
!  
!  ?? Oct 1993  B. Schwartz   Created.
!  24 Aug 1994  E. Boer       Modified dimension output.
!  29 Jul 1997  S. Rupert     Standardized and reinstated required include of netcdf.inc. 
!  30 Jul 1997  S. Rupert     Added usage message and command line inputs.
!  03 Apr 2003  H. Yan        Change output to fortran 90 format, modified maximum dimesions to 5.
!  30 Mar 2004  H. Yan        Added usage message, retrieve "scale_factor" and "add_offset" parameters.
!  18 May 2006  H. Yan        Fixed bugs of dimensions definition and remark lines. 
!-------------------------------------------------------
!  Do not forget to include the -I path_to_netcdf_  includes in your compile statement Required includes.
!  Also note: need 'netcdf.lib' or 'netcdfs.lib' when link .
      include 'netcdf.inc'
!
!-------------------------------------------------------
!
!    Using below command line to get the head information of the NetCDF file:
!    Please execute " ncdump -h mslp.mon.mean.nc > output.txt " 
!-------------------------------------------------------
!
!     Define Variables.
!     Variable ids run sequentially from 1 to nvars=5 ! number of variables
      integer,parameter :: nrec=408  ! change this to generalize
      integer*4  ncid, status    ! file control
      integer*4 recdim   ! record dimension
!-------------------------------------------------------------
 !     Below 5 variables is the data in netCDF file
      real*4         ::  lat( 73 )
      real*4         ::  lon( 144 )
      integer*2      ::  mslp( 144, 73, nrec )
      real*8         ::  time(nrec)
      real*8         ::  time_bnds( 2, 408 )
!     above5 variables is the data in netCDF file
!-------------------------------------------------------------
      integer*4   :: start(10)
      integer*4   :: count(10)
       integer*4   :: dimids(10)! allow up to 10 dimensions
      integer*4   :: dimid, xtype
      character(len=31) :: dummy
	character(len=100) :: filename
!----------------------------------------------------------------
!  Define "scale_factor" and "add_offset" variables
       real*4     ::  scale3(1), add3(1)
	 real      ::  slptmp( 144, 73, nrec )
!----------------------------------------------------------------

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
!   Retrieve data for Variable 'mslp'
!   Long_name of 'mslp' is 'Monthly Mean Sea Level Pressure'
!   Units of 'mslp' is 'Pascals'                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
      status=nf_inq_var(ncid,   3,dummy,xtype,ndim,dimids,natts)
          if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      do j=1,ndim
      status=nf_inq_dim(ncid,dimids(j),dummy,len)
           if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      start(j)=1 ; count(j)=len
      end do
      status=nf_get_vara_int2(ncid,   3,start,count,mslp)

      scale3(1) =0.0 ; add3(1) =0.0
! Scale_factor and add_offset for variable 'mslp'
      status=nf_get_att_real(ncid,   3,'add_offset',add3(1))
      status=nf_get_att_real(ncid,   3,'scale_factor',scale3(1))

! add scale_factor and add_offset to get true value
! Caution: variable type of 'mslp' may not as the same as the type of                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
! 'scale3(1)' and 'add3(1)', you must change it youself!
      slptmp=mslp*scale3(1)+add3(1)

!----------------------------------------------------
!   Retrieve data for Variable 'time'
!   Units of 'time' is 'hours since 1800-1-1 00:00:00'                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
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
!   Retrieve data for Variable 'time_bnds'
!   Long_name of 'time_bnds' is 'Time Boundaries'
      status=nf_inq_var(ncid,   5,dummy,xtype,ndim,dimids,natts)
          if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      do j=1,ndim
      status=nf_inq_dim(ncid,dimids(j),dummy,len)
           if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      start(j)=1 ; count(j)=len
      end do
      status=nf_get_vara_double(ncid,   5,start,count,time_bnds)

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


!----------------------------------------------------
!  End Program
      return 
      end subroutine

      subroutine GPCP(filename, gpctmp)
      real data1(144,72)
	real gpctmp(144,72,12)
      integer i,j,month,iret,nskp,k
      character*4  year
	character*100  filename


      open( UNIT=10, FILE=trim(filename)  &
     & , ACCESS='DIRECT', FORM='UNFORMATTED', STATUS='OLD', RECL=144 &
     & , IOSTAT=iret,convert='big_endian')
       do month=1,12 
        nskp=1+(month-1)*72
        do j=1,72
        read(10,rec=j+nskp,iostat=iret) (data1(i,j),i=1,144)
        enddo
       do j=1,72
        do i=1,144
         gpctmp(i,j,month)=data1(i,j)
	  enddo
	 enddo
      enddo
     close(10)
      return
      end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111
     subroutine normalization(input,N,M,output)
      real, dimension(N,M) :: input,output
!	  real, dimension(N,M) :: mean,std
       

    
     temp=0.0
	 mean=0
	 std=0
      
  
	 temp=0.0
	  temp=sum(input(:,:))
	  mean=temp/(N*M*1.0) !!!!average	 
	 temp=0.0
	 temp2=0.0
	 do i=1,N
	 do j=1,M
	  temp=temp+(input(i,j)-mean)**2
	  temp2=temp2+input(i,j)
	 enddo
	 enddo
!	 print*,temp2,'temp2',mean(i,j)*k
	  std=sqrt(temp/((N*M-1)*1.0))
	 do i=1,N
	 do j=1,M
	  output(i,j)=(input(i,j)-mean)/std
!	  print*,output(i,j,iK),i,j,ik

!	  print*,input(i,j,iK)
	  enddo
      enddo
	return
	end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
	subroutine QIYI(gpc,area)
	real gpc(7,34)
	character*4 :: area(7)
    open(10,file='Z:\DATA\LargeScale\79-12\GPCP\qiyizhi_1std.txt')
	open(20,file='Z:\DATA\LargeScale\79-12\GPCP\qiyizhi_2std.txt')
	open(30,file='Z:\DATA\LargeScale\79-12\GPCP\qiyizhi_3std.txt')
	do ir=1,7
	avg=sum(gpc(ir,:))
	avg=avg/34.0
	temp=0.0
	 do iy=1979,2012
	 temp=temp+(avg-gpc(ir,iy-1978))**2
	 enddo
	std=sqrt(temp/34.0)
!	print*,avg, std
     do iy=1979,2012
	 rain=gpc(ir,iy-1978)
!	 print*,rain-avg-3*std,rain-avg+3*std
	 if((rain-avg-1*std>0).and.(rain-avg-2*std<0))then 
	   write(10,99)trim(area(ir)), iy, rain, avg, std, 'Wet'
      endif
     if((rain-avg+1*std<0).and.(rain-avg+2*std>0))then 
	   write(10,99)trim(area(ir)), iy, rain, avg, std, 'Dry'
      endif
	  	 if((rain-avg-2*std>0).and.(rain-avg-3*std<0))then 
	   write(20,99)trim(area(ir)), iy, rain, avg, std, 'Wet'
      endif
     if((rain-avg+2*std<0).and.(rain-avg+3*std>0))then 
	   write(20,99)trim(area(ir)), iy, rain, avg, std, 'Dry'
      endif
	   if(rain-avg-3*std>0)then 
	   write(30,99)trim(area(ir)), iy, rain, avg, std, 'Wet'
      endif
     if(rain-avg+3*std<0)then 
	   write(30,99)trim(area(ir)), iy, rain, avg, std, 'Dry'
      endif
	 enddo
	enddo
	close(10)
99 format(1X, A4,1X,I4,3(1X,F8.4),1X,A3)
    
	end subroutine
	
