      program CN5_21
!!!     Data information 
!       1961-2010_Full_daily_0.5x0.5 degree
! XDEF 142 LINEAR 69.75 0.5
! YDEF 82 LINEAR 14.75 0.5
! ZDEF 1 levels 1000
! TDEF 18262 LINEAR 01jan1961 1dy
!   UNDEF -9999.0   
! 
	implicit none
	integer,parameter :: NX=142
	integer,parameter :: NY=82
      character fnm*100,dir*100,nm(2)*80
	integer i,j,k,it,day(12)
	integer iy,iyy,nd,id,im
	real raw(2,50,12,31,NX,NY) !,rian(50,366,NX,NY)



!------- path strings--------------
	dir='Z:\DATA\CN05\CN05.2\'
      nm(2)='2400_China_Pre_1961_2010_Full_daily_05x05.dat'
	nm(1)='2400_China_Tm_1961_2010_Full_daily_05x05.dat'
!---------------------------------------------------------
      do i=1,12
	day(i)=30
	enddo
	day(1)=31;day(3)=31;day(5)=31;day(7)=31
	day(8)=31;day(10)=31;day(12)=31


	do k=1,2
	 fnm=trim(dir)//trim(nm(k))
	 open(10,file=trim(fnm),form='binary') 
	  do iyy=1961,2010
	    iy=iyy-1960
         day(2)=28
!-------------- 
	    if(mod(iyy,4)==0.and.mod(iyy,100)/=0)then
	     day(2)=29
	    elseif(mod(iyy,400)==0)then
	   	 day(2)=29
	    endif
!-----------------------------          
		do im=1,12
	       nd=day(im)
           do id=1,nd
	      do j=1,NY
	       do i=1,NX
              read(10)raw(k,iy,im,id,i,j)
	       enddo
	      enddo
	     enddo
	   enddo
	  enddo
	enddo
!----------test raw(1,1,im,id,i,j)
	print*, raw(1,1,1,1,1,1),raw(1,1,1,1,1,2)
	print*, raw(1,1,1,1,2,1),raw(1,1,1,1,2,1)
	print*, raw(1,1,1,1,1,81),raw(1,1,1,1,1,82)
!------- the result  --------------------
	call output(raw)
	
	
	
	end

	subroutine output(input)
	implicit none
	integer,parameter :: NX=142
	integer,parameter :: NY=82
	integer,parameter :: ngn=7  !! No. of regions
	real input(2,50,12,31,NX,NY)  !!!!	50 years 
	real monmean(2,50,12,NX,NY),monstd(2,12,NX,NY)
	real mymonmean(2,12,NX,NY)
	real out(2,50,12,ngn),lon,lat !!! 98-00 monthly
	real outday(2,50,12,31,ngn),outdaymy(2,12,ngn)
	real tlons,tlone,tlats,tlate
	real temp(2,ngn),XYD,temp2,mymnstd(2,12,ngn)
	real temp3(12),temp4,temp5,temp6,temp7,temp8,temp9
	real temp10
	integer iysd,indef,indef2(12),isdudef,kk2,ksp
	integer kk3
	integer iys,iye,ims,ime,ids,ide
	integer lons(ngn),lone(ngn),lats(ngn),late(ngn)
	integer day(12),kk(2,ngn)
	integer iy,im,nd,id,ig,i,j,k,iyy
	integer nt,ievn
	real undef
	real forout(19,7),prout(7,23,17),prout1(23,17)
	integer lv0,drylv1,drylv2,drylv3,wetlv1,wetlv2,wetlv3
	integer rainstrom,ikk,ikn(7,23,17),ik   !!!! area-averaged rainfall >50mm/d
	character rname(ngn)*4,vrname(23)*10

	undef=-9999.0
!-----------  
!-----------EA--------------------------------------------------
	lons(1)=100;lone(1)=140
	lats(1)=10;late(1)=55
	rname(1)='EA'
!-----------PRD--------------------------------------------------
	lons(2)=110;lone(2)=118
	lats(2)=21;late(2)=25
	rname(2)='PRD'
!-----------MLYR--------------------------------------------------
	lons(3)=110;lone(3)=122
	lats(3)=27;late(3)=33
	rname(3)='MLYR'
!-----------NPC--------------------------------------------------
	lons(4)=112;lone(4)=120
	lats(4)=34;late(4)=42
	rname(4)='NPC'
!-----------NEC--------------------------------------------------
	lons(5)=120;lone(5)=130
	lats(5)=43;late(5)=49
	rname(5)='NEC'
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------ETP--------------------------------------------------
	lons(6)=90;lone(6)=102.5
	lats(6)=27.5;late(6)=37.5 
	rname(6)='ETP'
!-----------WTP--------------------------------------------------
	lons(7)=77.5;lone(7)=90
	lats(7)=27.5;late(7)=37.5
	rname(7)='WTP'
!---------------------------------------------------------
      vrname(1)='T_forcing'
      vrname(2)='Q_forcing'
	vrname(3)='uwnd'
	vrname(4)='vwnd'
	vrname(5)='qv'
	vrname(6)='hgt'
	vrname(7)='temp'
	vrname(8)='omg_adj'
      vrname(9)='Rh'
      vrname(10)='Tv'
	vrname(11)='Q1'
      vrname(12)='Q2'
      vrname(13)='HAD_Q'
      vrname(14)='VAD_Q'
      vrname(15)='TCH_Q'
      vrname(16)='HAD_T'
	vrname(17)='VAD_T'
      vrname(18)='TCH_T'
      vrname(19)='omg_ori'
	vrname(20)='VorX'
	vrname(21)='VorY'
	vrname(22)='VorW'
	vrname(23)='Div'
      do i=1,12
	day(i)=30
	enddo
	day(1)=31;day(3)=31;day(5)=31;day(7)=31
	day(8)=31;day(10)=31;day(12)=31
!!!!!------ the output data period----------------------
	iys=1961;iye=2010
	ims=1;ime=12
	ids=1;ide=31
      

! 	real input(2,50,12,31,NX,NY)  !!!!	50 years 
      do i=1,nx ! mk1
	  do j=1,ny ! mk2
	    do k=1,2  ! temp, rian  mk3
	      temp3=0.0
	      indef2=0
            do iyy=iys,iye  ! mk4
	       iy=iyy-1960
	       day(2)=28
!-------------- 
	       if(mod(iyy,4)==0.and.mod(iyy,100)/=0)then
	         day(2)=29
	       elseif(mod(iyy,400)==0)then
	   	     day(2)=29
	       endif
!----------------------------
             do im=1,12 !mk5
	        temp2=0.0
	         indef=0
	         do id=1,day(im) ! mk6 
	           if(input(k,iy,im,id,i,j)==undef)then
	            indef=indef+1
			   else
	            temp2=temp2+input(k,iy,im,id,i,j)
                 endif
			 enddo !mk6
		      if(indef<day(im))then
			   monmean(k,iy,im,i,j)=temp2/(day(im)-indef)
	          else
                 monmean(k,iy,im,i,j)=undef
	          endif
			   temp2=0.0
			  if(monmean(k,iy,im,i,j)==undef)then
			   indef2(im)=indef2(im)+1
			  else	          
	           temp3(im)=temp3(im)+monmean(k,iy,im,i,j)
		      endif
		     enddo !mk5  month
		   enddo !mk4  years
	        do im=1,12
	          if(indef2(im)<50)then
		      mymonmean(k,im,i,j)=temp3(im)/(50.0-indef2(im))
	          else
                 mymonmean(k,im,i,j)=undef
	          endif
			enddo
	        temp3=0.0
!!!<<<<<<<<<<<<<<<<<<<  sub loop <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<              
               do im=1,12 
	          temp4=0.0
	          isdudef=0
                 do iysd=1,50
	             if(monmean(k,iysd,im,i,j)/=undef.and.
     +				mymonmean(k,im,i,j)/=undef)then
                       temp4=temp4+
     +				(monmean(k,iysd,im,i,j)-mymonmean(k,im,i,j))**2
	              else
	               isdudef=isdudef+1
	              endif                     
                  enddo
	             if(isdudef<50)then
                    monstd(k,im,i,j)=(temp4/(50.0-isdudef))**0.5
	              else
	              monstd(k,im,i,j)=undef
				 endif
                enddo
!!!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	    enddo !mk3
	   enddo !mk2
	 enddo !mk1 	 
 !---------------------------------------------------------------
 !!  results:  monstd: mean square deviation of every grid
 !             monmean : monthly mean of every grid  
 !             mmonmean : multi-year monthly mean 
 !   if these results are wanted to be dispalyed, please output these variables.          
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!!! area-averaged for different regions 
	do k=1,2  ! temp, rian 
        do iyy=iys,iye
	    iy=iyy-1960
	    day(2)=28
!-------------- 
	    if(mod(iyy,4)==0.and.mod(iyy,100)/=0)then
	     day(2)=29
	     elseif(mod(iyy,400)==0)then
	   	 day(2)=29
	    endif
!----------------------------
         do im=1,12
	    nd=day(im)
		kk=0
	    temp=0.
           do ig=1,ngn ! region
	      tlons= lons(ig);tlone= lone(ig)
	      tlats= lats(ig);tlate= late(ig)
	      do id=1,nd
	       temp7=0.0
	       ksp=0
	       do j=1,NY
	         lat=14.75+0.5*(j-1)
	         do i=1,NX
	           lon=69.75+0.5*(i-1)
			    if(input(k,iy,im,id,i,j)/=undef)then
				 if(lon>=tlons.and.lon<=tlone)then 
	              if(lat>=tlats.and.lat<=tlate)then
	               temp(k,ig)=temp(k,ig)+input(k,iy,im,id,i,j)
	               kk(k,ig)=kk(k,ig)+1
	               temp7=temp7+input(k,iy,im,id,i,j)
	               ksp=ksp+1
                    endif
				 endif
			    endif
	          enddo !NX
	        enddo !NY
	        outday(k,iyy-1960,im,id,ig)=temp7/ksp
	       enddo ! nd
	       XYD=kk(k,ig) !*nd*1.0
!	print*,kk(k,ig)
	      out(k,iyy-1960,im,ig)=temp(k,ig)/XYD
!!! ((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((    area-averaged variance 
	       temp6=0.0
	       kk2=0
	       temp9=0.0
		   kk3=0
		   do j=1,NY
	         lat=14.75+0.5*(j-1)
	         do i=1,NX
	          lon=69.75+0.5*(i-1)
	           if(monstd(k,im,i,j)/=undef)then
!!!!!!!!!!!!!!!!!!!!!!!!!!	print*,'###',lon,tlons,tlone
	            if(lon>=tlons.and.lon<=tlone)then 
	             if(lat>=tlats.and.lat<=tlate)then
	                temp6=temp6 + monstd(k,im,i,j)
	                 kk2=kk2+1
                   endif
	            endif
	           endif
!	outdaymy(2,12,ngn)
                  if(mymonmean(k,im,i,j)/=undef)then
	            if(lon>=tlons.and.lon<=tlone)then 
	             if(lat>=tlats.and.lat<=tlate)then
	                temp9=temp9 + mymonmean(k,im,i,j)
	                 kk3=kk3+1
                   endif
	            endif
	           endif
	          enddo !NX
	        enddo !NY
	         mymnstd(k,im,ig)=temp6/kk2  !!! multi-year area-averaged mean std for different regions
	         outdaymy(k,im,ig)=temp9/kk3
!!!)))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
	     enddo !ig
	  enddo !im
       enddo !iy
	enddo !K
	open(100,file='tm_mon36_61-10.txt')
	open(200,file='Pre_mon36_61-10.txt')
	open(300,file='multi-year_sd.txt')
	write(100,99)'Year','month','EA','PRD','MLYR','NPC'
     +	,'NEC','ETP','WTP'
	write(200,99)'Year','month','EA','PRD','MLYR','NPC'
     +	,'NEC','ETP','WTP'
	write(300,99)'Year','month','EA','PRD','MLYR','NPC'
     +	,'NEC','ETP','WTP'
	do iy=1,50
	  do im=1,12
	    write(100,98)iy+1960,im,(out(1,iy,im,ig),ig=1,ngn)
	    write(200,98)iy+1960,im,(out(2,iy,im,ig),ig=1,ngn)
	  enddo
	enddo
	do im=1,12 
	write(300,97)'Temperature',im,(mymnstd(1,im,ig),ig=1,ngn)
	enddo
	write(300,'(A46)')'**********************************************'
	do im=1,12 
	write(300,97)'Pericitation',im,(mymnstd(2,im,ig),ig=1,ngn)
	enddo
!!------------------------------------------------------------------------------------------
!!  results:  out: store the area-averaged temperature and pericipation for different regions  
!!            mymnstd:  multi-year area-averaged mean std for different regions 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--------------------------------------------------------------------------------------------
!***********************************************************************************************
!----------- seclected the dry and wet condition -----------------------------------------------
!-----------------------------------------------------------------------------------------------
      open(500,file='EA_temp.txt')
	open(501,file='PRD_temp.txt')
	open(502,file='MLYR_temp.txt')
	open(503,file='NPC_temp.txt')
	open(504,file='NEC_temp.txt')
	open(505,file='ETP_temp.txt')
	open(506,file='WTP_temp.txt')
	open(507,file='EA_rain.txt')
	open(508,file='PRD_rain.txt')
	open(509,file='MLYR_rain.txt')
	open(510,file='NPC_rain.txt')
	open(511,file='NEC_rain.txt')
	open(512,file='ETP_rain.txt')
	open(513,file='WTP_rain.txt')
	open(514,file='EA_rainstrom.txt')
	open(515,file='PRD_rainstrom.txt')
	open(516,file='MLYR_rainstrom.txt')
	open(517,file='NPC_rainstrom.txt')
	open(518,file='NEC_rainstrom.txt')
	open(519,file='ETP_rainstrom.txt')
	open(520,file='WTP_rainstrom.txt')
	do i=1,7
	   write(499+i,96)'Year','Month','Day', 'Temp'
     +      ,'lv0','drylv1','drylv2','drylv3','wetlv1','wetlv2','wetlv3'
	   write(506+i,96)'Year','Month','Day', 'Rain'
     +      ,'lv0','drylv1','drylv2','drylv3','wetlv1','wetlv2','wetlv3'
	   write(513+i,*)'Year',' ','Month',' ','Day', ' ','Rain'
	   open(520+i,file=trim(rname(i))//'_drylv2.txt')
	   write(520+i,*)'Year',' ','Month',' ','Day', ' ','Rain'
	   open(527+i,file=trim(rname(i))//'_wetlv2.txt')
	   write(527+i,*)'Year',' ','Month',' ','Day', ' ','Rain'
	   open(534+i,file=trim(rname(i))//'_drylv3.txt')
	   write(534+i,*)'Year',' ','Month',' ','Day', ' ','Rain'
	   open(541+i,file=trim(rname(i))//'_wetlv3.txt')
	   write(541+i,*)'Year',' ','Month',' ','Day', ' ','Rain'
	   open(548+i,file=trim(rname(i))//'_EventsDate_cloudsat.txt')
	   write(548+i,*)'Year',' ','Month',' ','Day', ' ','Rain'
!
	   open(555+i,file=trim(rname(i))//'_Intagated_Events_all.txt')
	   write(555+i,91)'Date',(vrname(k),k=1,19)
	    open(555+7+i,file=trim(rname(i))//'_Intagated_Events_low.txt')
	   write(555+7+i,91)'date',(vrname(k),k=1,19)
	   open(555+14+i,file=trim(rname(i))//'_Intagated_Events_mid.txt')
	   write(555+14+i,91)'date',(vrname(k),k=1,19)
	   open(555+21+i,file=trim(rname(i))//'_Intagated_Events_high.txt')
	   write(555+21+i,91)'date',(vrname(k),k=1,19)
	   open(555+28+i,file=trim(rname(i))//'_Intagated_Events_tro.txt')
	   write(555+28+i,91)'date',(vrname(k),k=1,19)
	   open(555+35+i,file=trim(rname(i))//'_meanprofile_Events.txt')
	   write(555+35+i,90)'level',(vrname(k),k=1,23)
	enddo
91    format(1X,A4,19(1X,A10))
90    format(1X,A5,23(1X,A10))
      do k=1,2  ! temp, rian 
       do iyy=iys,iye

	    iy=iyy-1960
	    day(2)=28
	    nt=1460
!-------------- 
	    if(mod(iyy,4)==0.and.mod(iyy,100)/=0)then
	     day(2)=29
	      nt=1464
	     elseif(mod(iyy,400)==0)then
	   	 day(2)=29
	     nt=1464
		endif
!----------------------------
         do ig=1,ngn ! region
		    ievn=0
           do im=1,12 ! month
	    nd=day(im)

	       do id=1,nd  ! day
		          lv0=0
	              drylv1=0
	              drylv2=0
	              drylv3=0 !.False.
	              wetlv1=0 !.False.
	              wetlv2=0 !.False.
	              wetlv3=0 !.False.
                    rainstrom=0 !.False.
	         temp8=outday(k,iy,im,id,ig)
     +      			 -outdaymy(k,im,ig)-mymnstd(k,im,ig)
	         temp10=outday(k,iy,im,id,ig)
     +      			 -outdaymy(k,im,ig)+mymnstd(k,im,ig)
               if(temp8>0)  wetlv1=1 !.True.
	         if(temp10<0) drylv1=1 !.True.
		     temp8=outday(k,iy,im,id,ig)
     +      			 -outdaymy(k,im,ig)-2*mymnstd(k,im,ig)
	         temp10=outday(k,iy,im,id,ig)
     +      			 -outdaymy(k,im,ig)+2*mymnstd(k,im,ig)
               if(temp8>0)  wetlv2=1 !.True.
	         if(temp10<0) drylv2=1 !.True.
			 temp8=outday(k,iy,im,id,ig)
     +      			 -outdaymy(k,im,ig)-3*mymnstd(k,im,ig)
	         temp10=outday(k,iy,im,id,ig)
     +      			 -outdaymy(k,im,ig)+3*mymnstd(k,im,ig)
	         if(temp8>0)  wetlv3=1 !.True.
	         if(temp10<0) drylv3=1 !.True.
!---------------------------------------------------------------------------
               if(k>1)then
			 if(outday(k,iy,im,id,ig)>=15.) rainstrom= 1 !.True.
			 endif			 
!!---------------output for the selected  period-------------------------------			 
                if(k==1)then
!			  if(drylv1==1.or.wetlv1==1)then
	          write(499+ig,95)iyy,im,id,outday(k,iy,im,id,ig)
     +              ,lv0,drylv1,drylv2,drylv3,wetlv1,wetlv2,wetlv3
!			  endif
	          endif
                if(k==2)then
			   if(drylv1==1.or.wetlv1==1)then
	           write(506+ig,95)iyy,im,id,outday(k,iy,im,id,ig)
     +              ,lv0,drylv1,drylv2,drylv3,wetlv1,wetlv2,wetlv3
			   endif
		        if(rainstrom==1)then
	             write(513+ig,94)iyy,im,id,outday(k,iy,im,id,ig)
			     endif
	          endif
                  if(drylv2==1.and.drylv3==0)then
	              write(520+ig,94)iyy,im,id,outday(k,iy,im,id,ig)
				endif
	            if(drylv3==1)then
	              write(534+ig,94)iyy,im,id,outday(k,iy,im,id,ig)
				endif
	            if(wetlv2==1.and.wetlv3==0)then
	              write(527+ig,94)iyy,im,id,outday(k,iy,im,id,ig)
				endif
	            if(wetlv3==1)then
	              write(541+ig,94)iyy,im,id,outday(k,iy,im,id,ig)
!	                if(iyy>2005)then
	                  if(im>4.and.im<10)then
	                   write(548+ig,94)iyy,im,id,outday(k,iy,im,id,ig)
					  endif
!	                 endif
	                if(ig>1.and.iyy>1978)then
	                  if(im>4.and.im<10)then !!!! MJJAS
	                   ievn=ievn+1
				       call RD(iyy,im,id,ig,nt,ievn
     +					   , outday(k,iy,im,id,ig),day,forout,prout1)
	             write(555+ig,93)iyy,im,id,(forout(ikk,7),ikk=1,19)
	             write(555+7+ig,93)iyy,im,id,(forout(ikk,2),ikk=1,19)
	             write(555+14+ig,93)iyy,im,id,(forout(ikk,3),ikk=1,19)
	             write(555+21+ig,93)iyy,im,id,(forout(ikk,4),ikk=1,19)
	             write(555+28+ig,93)iyy,im,id,(forout(ikk,6),ikk=1,19)
	                     prout(ig,:,:)=prout(ig,:,:)+prout1(:,:)
	                     ikn(ig,:,:)=ikn(ig,:,:)+1
                        endif
	                endif
				endif
              enddo !days
	       enddo ! month
	      enddo !region
	     enddo ! year
	     do ig=1,7
             do ik=1,17
        write(555+35+ig,92)ik,(prout(ig,ikk,ik)/ikn(ig,ikk,ik),ikk=1,23)
             enddo
		 enddo  	  
	  enddo ! k
99    format(1X,A4,1X,A5,7(1X,A4))
98    format(1X,I4,1X,I2,7(1X,F12.6))
97    format(1X,A15,1X,I2,7(1X,F12.6))
96    format(1X,A4,1X,A2,1X,A3,1X,A4,7(1X,A6))
95    format(1X,I4,1X,I2,1X,I2,1X,F12.6,7(1X,I2))
94    format(1X,I4,1X,I2,1X,I2,1X,F12.6)
93    format(1X,I4,I2.2,I2.2,19(1X,F12.6))
92    format(1X,I2,23(1X,F12.6))
	end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine RD(iyr,im,id,ig,nt,ievn,rain,day,forout,prout)
	integer,parameter :: nrec=1464
	integer,parameter :: nvr=19
	integer,parameter :: nYY=34
	character(len=100) file,fold,dir,filename
	character(len=20) area(7),levname(8)
      character(len=10) cel(1464)
	character*5 hour(1464)
      character(len=10) dayID(nrec/4)
	real for(nvr,7,1464),raw(17,nvr,1464)
	real dyn(17,4,1464),rain,temp(23)
	character(len=4) year
	real tempfor,tempro
	real forout(19,7),prout(23,17)
	character evnid*3,vrname(23)*10
      integer iyr,im,id,ig,ievn,nt,i
	integer day(12),ids,ide,it,j
	integer l,kk
      area(1)='EA'
      area(2)='PRD'  ;    area(3)='MLYR'
	area(4)='NPC'  ;    area(5)='NEC'
	area(6)='ETP'  ;    area(7)='WTP'
      levname(1)='_lowlevel.txt'
	levname(2)='_Middlelevel.txt'
	 levname(3)='_Hightlevel.txt'
	 levname(4)='_Abovelevel.txt'
	 levname(5)='_Tropsphere.txt'
	 levname(6)='_AllLevels.txt'
	levname(7)='_RAW.txt'
	levname(8)='_dyn.txt'
!!!----------------------------------      
	i=ig
!---------------------------------------------
	if(ig<6)then
	  dir='Z:\DATA\LargeScale\NcepR2_Pre\'
	else
	  dir='Z:\DATA\LargeScale\TP\NcepR2_Pre\'
	endif
      write(year,'(I4)')iyr !-1+1979
	fold=year(3:4)//'0101'//'-'//year(3:4)//'1231\'
!-------------surface---------------------------
	  file=year//trim(area(i))//'_SURFACE.txt'
	  filename=trim(dir)//trim(fold)//trim(file)
        open(10,file=trim(filename))
	  read(10,*)
	  do it=1,nt
	    read(10, 92)cel(it),hour(it),(for(k,1,it),k=1,6)
	    do k=1,6
	     if(k==3.or.k==4)then
	      for(k,1,it)=for(k,1,it)*3600. !!!convert rainfall rate to mm/hr
	     endif
	    enddo 
         enddo
	  close(10)
	 
	 do j=1,6
	  file=year//trim(area(i))//trim(levname(j))
	   filename=trim(dir)//trim(fold)//trim(file)
        open(10,file=trim(filename))
	   read(10,*)
	   do it=1,nt
	    read(10, 93)cel(it),hour(it),(for(k,j+1,it),k=1,nvr) 
          enddo
	   close(10)
	 enddo
!-----------------------------------------------------------
        file=year//trim(area(i))//trim(levname(7))
	  filename=trim(dir)//trim(fold)//trim(file)
	  open(70,file=trim(filename))
	   read(70,*)
        do it=1,nt
	   do ik=1,17 
	   read(70, 93)cel(it),hour(it),(raw(ik,k,it),k=1,nvr) 
        enddo
	 enddo
	close(70)
	file=year//trim(area(i))//trim(levname(8))
	  filename=trim(dir)//trim(fold)//trim(file)
	  open(80,file=trim(filename))
	   read(80,*)
        do it=1,nt
	   do ik=1,17 
	   read(80, 93)cel(it),hour(it),(dyn(ik,k,it),k=1,4) 
        enddo
	 enddo
	close(80)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ---------- output the events profiles --------------------------------
      vrname(1)='T_forcing'
      vrname(2)='Q_forcing'
	vrname(3)='uwnd'
	vrname(4)='vwnd'
	vrname(5)='qv'
	vrname(6)='hgt'
	vrname(7)='temp'
	vrname(8)='omg_adj'
      vrname(9)='Rh'
      vrname(10)='Tv'
	vrname(11)='Q1'
      vrname(12)='Q2'
      vrname(13)='HAD_Q'
      vrname(14)='VAD_Q'
      vrname(15)='TCH_Q'
      vrname(16)='HAD_T'
	vrname(17)='VAD_T'
      vrname(18)='TCH_T'
      vrname(19)='omg_ori'
	vrname(20)='VorX'
	vrname(21)='VorY'
	vrname(22)='VorW'
	vrname(23)='Div'
      dir='Z:\DATA\LargeScale\79-12\PHD2\Evens\'
	write(evnid,'(I3.3)')ievn
      file=trim(area(i))//'_MJJAS_'//year//evnid//'_profile.txt'
	filename=trim(dir)//trim(file)
	open(20,file=trim(filename))
	write(20,91)area(ig), iyr,im,id,rain
	write(20,90)'levels',(trim(vrname(k)),k=1,nvr+4) 
      
	ids=0
	do i=1,im-1
	 ids=ids+day(i)*4
	enddo
	ids=(id-1)*4+ids
      do ik=1,17 
	  do k=1,19
	     temp(k)=0
		 do it=ids,ids+3
	      temp(k)=temp(k)+raw(ik,k,it)
	     enddo
           temp(k)=temp(k)/4
	     prout(k,ik)=temp(k)
         enddo
!	print*,k
	   do i=1,4
         
	    temp(k)=0
	    do it=ids,ids+3 
           temp(k)=dyn(ik,i,it)+temp(k)
	    enddo
          temp(k)=temp(k)/4
	    prout(k,ik)=temp(k)
	    k=k+1
	   enddo
	 write(20,89)ik,(temp(k),k=1,23)
      enddo
	close(20)
!!!!! for()-------------------------------------
       do l=1,7
	  kk=19
	  if(l==1)then
         kk=6
	  endif
	   do k=1,kk	   
	     tempfor=0
		 do it=ids,ids+3
	      tempfor=tempfor+for(k,l,it)
	     enddo
           forout(k,l)=tempfor/4
	   enddo	  
       enddo
!	print*,k

93    format(1X,A10,1X,A5,1X,19(1X,E12.4))
92    format(1X,A10,1X,A5,1X,6(1X,E12.4))
91    format(1X,A4,1X,I4,1X,I2,1X,I2,F12.6)
90    format(1X,A6,2(1X,A10),21(1X,A7))
89    format(1X,I4,23(1X,E12.4))
!-----------------------target--------------------------------------	
!      print*,for(7,2,1,1,1),cel(1),hour(1)
	return
	end subroutine