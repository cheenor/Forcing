      program Forcing_AVE
	integer,parameter :: nrec=1464
	integer,parameter :: nvr=19
	integer,parameter :: nYY=34
      real taget(nvr,7,nrec,4,NYY),raw(17,nvr,1464,4,NYY)
	character(len=10) dayID(nrec/4,NYY)
	character(len=10) cel(1464)
	character*5 hour(1464) 
      real dyn(17,4,1464,4,NYY) ! 1-3 vor   4 DIV
      real dynd(17,4,1464/4,4,NYY)
	real dynm(17,4,12,4,NYY)
	real dyns(17,4,4,4,NYY)
	real dyny(17,4,4,NYY)
	real xdynd(17,4,1464/4,4)
	real xdynm(17,4,12,4)
	real xdyns(17,4,4,4)
	integer nt
      COMMON/DD/days(12)
	real Dayly(nvr,7,nrec/4,4,NYY), XMonth(nvr,7,12,4,NYY),
     +     Seanson(nvr,7,4,4,NYY),
     +     Yrave(nvr,7,4,NYY)
	real Dayr(17,nvr,nrec/4,4,NYY),XmonthR(17,nvr,12,4,NYY),
     +   Yrrave(17,nvr,4,NYY),Seansonr(17,nvr,4,4,NYY)
      real Xmdayly(nvr,7,nrec/4,4),XMMonth(nvr,7,12,4),
     + XMSeanson(nvr,7,4,4),Xmdayr(17,nvr,nrec/4,4),
     + XMSeasonr(17,nvr,4,4),XMMonthr(17,nvr,12,4)
      
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
	 nt=1460
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	do iyr=1,NYY
	   iys=iyr-1+1979
	   nt=1460
	   days(2)=28
         if(mod(iys,4)==0.and.mod(iys,100)/=0)then
	    days(2)=29 ! nrec=1464
	    nt=1464
	   elseif(mod(iys,400)==0)then
	   days(2)=29
	   nt=1464
	   endif
         Call RD(iyr,nt,taget,raw,dyn)
	   Call DAY(iyr,nt,taget,raw,Dayly,dayr,dyn,dynd)
	   call Monthave(iyr,nt,dayly,dayr,XMonth,XmonthR,dynd,dynm)
	   call Year(iyr,nt,taget,raw,Yrave,Yrrave,dyn,dyny)
      enddo
!	   print*,XMonth(1,2,1,1,1),'M1111111111'

         Call Seasonave(XMonth,XMonthr,Seanson,SeansonR,dynm,dyns)

!	   print*,XMonth(1,2,1,1,1),'M2222222222'

	   Call Myrave(Dayly,XMonth,Seanson,Yrave,
     + Dayr,XmonthR,Yrrave,Seansonr,
     + Xmdayly,XMMonth, XMSeanson,Xmdayr,XMSeasonr,XMMonthr,
     + dynd,dynm,dyns,xdynd,xdynm,xdyns) 

	   call output_bin(taget,raw,Dayly,XMonth,Seanson,Yrave,
     + Dayr,XmonthR,Yrrave,Seansonr,
     + Xmdayly,XMMonth, XMSeanson,Xmdayr,XMSeasonr,XMMonthr,dyn,
     + dynd,dynm,dyns,xdynd,xdynm,xdyns) 
	  call output_ACS(taget,raw,Dayly,XMonth,Seanson,Yrave,
     + Dayr,XmonthR,Yrrave,Seansonr,
     + Xmdayly,XMMonth, XMSeanson,Xmdayr,XMSeasonr,XMMonthr,dyn,
     + dynd,dynm,dyns,xdynd,xdynm,xdyns) 
!-----------------------------------------------------------------

       do i=1,12
	  days(i)=31
	 enddo
	 days(2)=28
	 days(6)=30
	 days(4)=30
	 days(9)=30
	 days(11)=30
      open(100,file='Z:\DATA\LargeScale\79-12\
     +TimeID_readme.txt')
	ntt=0
	ism=1
	imm=0
	ntt0=0
	imn=0
	do 1003 iyr=1,NYY	 
       iys=iyr-1+1979
	 nt=1460
	 days(2)=28
        if(mod(iys,4)==0.and.mod(iys,100)/=0)then
	    nt=1464
	   days(2)=29
	  elseif(mod(iys,400)==0)then
	   nt=1464
	   days(2)=29
	   endif
      if(iyr==1)then
	ntt0=0
	else
	ntt0=ntt1
	endif
      
	ntt=ntt+nt
      ntt1=ntt
	write(100,94)iys, ntt
	do im=1,12
	   imm=imm+days(im)*4
	   imn=imn+1
	   write(100,95)im, ism, imm,imn
         ism=imm+1
      enddo
1003  continue
      close(100)
94    format(1X,I4,1X,I7)
95    format(1X,3X,I2,1X,I7,1X,I7,1X,I4) 



      end
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine RD(iyr,nt,for,raw,dyn)
	integer,parameter :: nrec=1464
	integer,parameter :: nvr=19
	integer,parameter :: nYY=34
	character(len=100) file,fold,dir,filename
	character(len=20) area(4),levname(8)
      character(len=10) cel(1464)
	character*5 hour(1464)
      character(len=10) dayID(nrec/4,NYY)
	real for(nvr,7,1464,4,NYY),raw(17,nvr,1464,4,NYY)
	real dyn(17,4,1464,4,NYY)
	character(len=4) year
      integer iyr,nt,nirr

      area(1)='WTP'  ;    area(2)='ETP'
	area(3)='BOB'  ;    area(4)='NEC'
	nirr=3
      levname(1)='_lowlevel.txt'
	levname(2)='_Middlelevel.txt'
	 levname(3)='_Hightlevel.txt'
	 levname(4)='_Abovelevel.txt'
	 levname(5)='_Tropsphere.txt'
	 levname(6)='_AllLevels.txt'
	levname(7)='_RAw.txt'
	levname(8)='_dyn.txt'
	dir='Z:\DATA\LargeScale\TP\NcepR2_Pre\'
      write(year,'(I4)')iyr-1+1979
	fold=year(3:4)//'0101'//'-'//year(3:4)//'1231\'
	do i=1,nirr !!!!area
!-------------surface---------------------------
	  file=year//trim(area(i))//'_SURFACE.txt'
	  filename=trim(dir)//trim(fold)//trim(file)
        open(10,file=trim(filename))
	  read(10,*)
	  do it=1,nt
	   read(10, 98)cel(it),hour(it),(for(k,1,it,i,iyr),k=1,6)
	   do k=1,6
	   if(k==3.or.k==4)then
	   for(k,1,it,i,iyr)=for(k,1,it,i,iyr)*3600. !!!convert rainfall rate to mm/hr
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
	   read(10, 99)cel(it),hour(it),(for(k,j+1,it,i,iyr),k=1,nvr) 
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
	   read(70, 99)cel(it),hour(it),(raw(ik,k,it,i,iyr),k=1,nvr) 
        enddo
	 enddo
	file=year//trim(area(i))//trim(levname(8))
	  filename=trim(dir)//trim(fold)//trim(file)
	  open(80,file=trim(filename))
	   read(80,*)
        do it=1,nt
	   do ik=1,17 
	   read(80, 99)cel(it),hour(it),(dyn(ik,k,it,i,iyr),k=1,4) 
        enddo
	 enddo
      enddo
99    format(1X,A10,1X,A5,1X,19(1X,E12.4))
98    format(1X,A10,1X,A5,1X,6(1X,E12.4))
!-----------------------target--------------------------------------	
!      print*,for(7,2,1,1,1),cel(1),hour(1)
	return
	end subroutine

	subroutine DAY(iyr,nt,for,raw,Dayly,dayr,dyn,dynd)
	integer,parameter :: nrec=1464
	integer,parameter :: nvr=19
	integer,parameter :: nYY=34
	character(len=10) dayID(nrec/4,NYY)
	character(len=10) cel(1464)
	character*5 hour(1464)
      COMMON/DD/days(12)
	real Dayly(nvr,7,nrec/4,4,NYY)
	real Dayr(17,nvr,nrec/4,4,NYY)
	integer nt,iyr,nirr
	real dyn(17,4,1464,4,NYY)
	real dynd(17,4,1464/4,4,NYY)
	real for(nvr,7,1464,4,NYY),raw(17,nvr,1464,4,NYY)


      temp=0.0
	nirr=3
	do 1001 i=1,nirr !!!area
      do 1001 j=1,7 
	do 1001 k=1,nvr
	do 1001 it=1,nt 
	  temp=temp+for(k,j,it,i,iyr)
	if(mod(it,4)==0)then
	  itt=it/4
	  Dayly(k,j,itt,i,iyr)=temp/4.0
	  temp=0.0
	  dayID(itt,iyr)=cel(it)
      endif
1001  continue
!      print*,Dayly(7,2,1,1,1),for(7,2,1,1,1)
      temp=0.0
      do 1002 i=1,nirr !!!area
      do 1002 ik=1,17 
	do 1002 k=1,nvr
	do 1002 it=1,nt 
	  temp=temp+raw(ik,k,it,i,iyr)
	if(mod(it,4)==0)then
	  itt=it/4
	  Dayr(ik,k,itt,i,iyr)=temp/4.0
	  temp=0.0
!	  dayID(itt,iyr)=cel(it)
      endif
1002  continue
      do 1003 i=1,nirr !area
      do 1003 ik=1,17 
	do 1003 k=1,4
	do 1003 it=1,nt 
	  temp=temp+dyn(ik,k,it,i,iyr)
	if(mod(it,4)==0)then
	  itt=it/4
	  dynd(ik,k,itt,i,iyr)=temp/4.0
	  temp=0.0
!	  dayID(itt,iyr)=cel(it)
      endif
1003  continue

      return
	end subroutine

      subroutine MONTHave(iyr,nt,Dayly,Dayr,XMonth,XmonthR,dynd,dynm)
      integer,parameter :: nvr=19 
	integer,parameter :: nrec=1464 
	integer,parameter :: nYY=34
      COMMON/DD/days(12)
	real XMonth(nvr,7,12,4,NYY)
	real XmonthR(17,nvr,12,4,NYY)
	integer nt,iyr,nirr
	real Dayly(nvr,7,nrec/4,4,NYY)
	real Dayr(17,nvr,nrec/4,4,NYY)
      real dynd(17,4,1464/4,4,NYY)
	real dynm(17,4,12,4,NYY)


      temp=0.0
	itt=0
	nirr=3
!	print*,days(2),iyr
!	ntt=days(itt)*4
	ntt1=0
	im=1
	ims=1
	do 1001 i=1,nirr !!!
      do 1001 j=1,7 
	do 1001 k=1,nvr
	temp=0.0
	itt=0
	ims=1
	im=1
	nd=365
	if(nt==1464)nd=366
	do 1001 it=1,nd 
	   
	  temp=temp+dayly(k,j,it,i,iyr)
	  itt=itt+1
	  if(itt==days(ims))then
	  xt=days(ims)*1.0
	  XMonth(k,j,im,i,iyr)=temp*1./xt
!	if(k==7.and.j==2.and.iyr==4.and.im==12.and.i==1)then
!	print*,XMonth(k,j,im,i,iyr)
!      stop
!	endif
!	print*,itt,ims,xt,temp,XMonth(k,j,im,i,iyr)
!	stop
 !        print*,temp,i,j,k,it
	  temp=0.0
	  ims=ims+1
	  im=im+1
	  itt=0
        endif
	  
1001  continue
!	stop 
      temp=0.0
	do  i=1,nirr
      do  j=2,7 
	do it=1,12
      if(XMonth(7,j,it,i,iyr)==0)then
	print*, j,it,i,iyr
	stop
	endif
      enddo
	enddo
	enddo
!	ntt=days(itt)*4
	ntt1=0
	ims=1
	im=1
	do 1002 i=1,nirr
      do 1002 ik=1,17 
	do 1002 k=1,nvr
	temp=0.0
	itt=0
	ims=1
	im=1
	nd=365
	if(nt==1464)nd=366
	do 1002 it=1,nd
	  
	  temp=temp+dayr(ik,k,it,i,iyr)
	  itt=itt+1
	 if(itt==days(ims))then	 
	  xt=days(ims)*1.0 
	  XMonthr(ik,k,im,i,iyr)=temp*1.0/xt
	  temp=0.0
	  ims=ims+1
	  im=im+1
	  itt=0
        endif
	  
1002  continue

	ntt1=0
	ims=1
	im=1
	do 1003 i=1,nirr
      do 1003 ik=1,17 
	do 1003 k=1,4
	temp=0.0
	itt=0
	ims=1
	im=1
	nd=365
	if(nt==1464)nd=366
	do 1003 it=1,nd
	  
	  temp=temp+dynd(ik,k,it,i,iyr)
	  itt=itt+1
	 if(itt==days(ims))then	 
	  xt=days(ims)*1.0 
	  dynm(ik,k,im,i,iyr)=temp*1.0/xt
	  temp=0.0
	  ims=ims+1
	  im=im+1
	  itt=0
        endif
	  
1003  continue
      return
	end subroutine
!
      subroutine seasonave(XMonth,XMonthr,Seanson,SeansonR,dnym,dnys)
	  integer,parameter :: nvr=19
	integer,parameter :: nYY=34
      real days(12)
	real XMonth(nvr,7,12,4,NYY),Seanson(nvr,7,4,4,NYY)
	real XmonthR(17,nvr,12,4,NYY),SeansonR(17,nvr,4,4,NYY)
	integer nt,iyr,nirr
	real s(4)
	real dynm(17,4,12,4,NYY)
	real dyns(17,4,4,4,NYY)
       nirr=3

	 do i=1,12
	  days(i)=31
	 enddo
	 days(2)=28
	 days(6)=30
	 days(4)=30
	 days(9)=30
	 days(11)=30
       s(1)=days(12)+days(1)+days(2)
	s(2)=days(3)+days(4)+days(5)
	s(3)=days(6)+days(7)+days(8)
	s(4)=days(9)+days(10)+days(11)
	
!	print*,XMonth(1,2,1,1,1),'AAAAAAAA'
!	stop
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      do 1001 ir=1,nirr
	do 1001 k=1,nvr
	
	do  j=1,7
	idk=0
	ims=1
	temp1=0.0 
	temp=0.0
	ims=1 
	do iyr=1,NYY
	    
	  iys=iyr-1+1979
	  days(2)=28
	  s(1)=days(12)+days(1)+days(2)
        if(mod(iys,4)==0.and.mod(iys,100)/=0)then
	   S(1)=S(1)+1.0
	   days(2)=29
	  elseif(mod(iys,400)==0)then
	   S(1)=S(1)+1.
	   days(2)=29
	   endif
	 ! 2D	
	 temp1=0.0   
        if(iyr==1)then !1990
	  ims=3
	  do im=1,2
	  temp1=temp1+XMonth(k,j,im,ir,iyr)*days(im)
	  enddo
       Seanson(k,j,1,ir,iyr)=(temp1+XMonth(k,j,12,ir,NYY)*days(12))/s(1)
	if(k==7)then
	print*,temp,Seanson(k,j,2,ir,iyr),XMonth(k,j,im,ir,iyr),'S1'
	print*,s(1),XMonth(k,j,12,ir,NYY)
      endif
        temp1=0.0 
	  else
        do im=1,2
	  temp1=temp1+XMonth(k,j,im,ir,iyr)*days(im)
	  enddo    
	  Seanson(k,j,1,ir,iyr)=
     +	  (temp1+XMonth(k,j,12,ir,iyr-1)*days(12))/s(1)
	  temp1=0.0 
	  endif

	  temp=0.0
	  do im =3,5
	  temp=temp*1.0+XMonth(k,j,im,ir,iyr)*days(im)
	  enddo
	   Seanson(k,j,2,ir,iyr)=temp/s(2)
	if(k==7)then
	print*,temp,Seanson(k,j,2,ir,iyr),XMonth(k,j,im,ir,iyr),'S2'
	print*,s(2),j
      endif   
	 temp=0.0
	  do im =6,8
	  temp=temp*1.0+XMonth(k,j,im,ir,iyr)*days(im)
	  enddo
	   Seanson(k,j,3,ir,iyr)=temp/s(3)
         
	  temp=0.0
	  do im =9,11
	  temp=temp*1.0+XMonth(k,j,im,ir,iyr)*days(im)
	 enddo
	   Seanson(k,j,4,ir,iyr)=temp/s(4)
        
	enddo
      enddo
      print*, Seanson(7,2,1,1,1)
      do  ik=1,17
	idk=0
	ims=1
	temp1=0.0 
	temp=0.0
	do iyr=1,NYY
		ims=1
	  iys=iyr-1+1979
	  days(2)=28
	  s(1)=days(12)+days(1)+days(2)
        if(mod(iys,4)==0.and.mod(iys,100)/=0)then
	   S(1)=S(1)+1
	   days(2)=29
	  elseif(mod(iys,400)==0)then
	   S(1)=S(1)+1
	   days(2)=29
	   endif
	 ! 2D	  
        if(iyr==1)then !1990
	  ims=3
	  do im=1,2
	  temp1=temp1+XMonthr(ik,k,im,ir,iyr)*days(im)
	  enddo
	  Seansonr(ik,k,1,ir,iyr)=
     +	  (temp1+XMonthr(ik,k,12,ir,NYY)*days(12))/s(1)
	  temp1=0.0 
	  else
        do im=1,2
	  temp1=temp1+XMonthr(ik,k,im,ir,iyr)*days(im)
	  enddo
        Seansonr(ik,k,1,ir,iyr)=
     +	  (temp1+XMonthr(ik,k,12,ir,iyr-1)*days(12))/s(1)
	  temp1=0.0 
	  endif

       temp=0.0
	  do im =3,5
	  temp=temp*1.0+XMonthr(ik,k,im,ir,iyr)*days(im)   
	   enddo   
	   Seansonr(ik,k,2,ir,iyr)=temp/s(2)
       
        temp=0.0
	  do im =6,8
	  temp=temp*1.0+XMonthr(ik,k,im,ir,iyr)*days(im)   
	  enddo   
	   Seansonr(ik,k,3,ir,iyr)=temp/s(3)
        
	 temp=0.0
	  do im =9,11
	  temp=temp*1.0+XMonthr(ik,k,im,ir,iyr)*days(im) 
	  enddo     
	   Seansonr(ik,k,4,ir,iyr)=temp/s(4)  
	   if(Seansonr(ik,7,4,ir,iyr)>500)then
	   print*,XMonthr(ik,k,im,ir,iyr),s(4)
	stop
	   endif  
	enddo

      enddo
1001  continue
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      do 1002 ir=1,nirr
	do 1002 k=1,4
      do  ik=1,17
	idk=0
	ims=1
	temp1=0.0 
	temp=0.0
	do iyr=1,NYY
		ims=1
	  iys=iyr-1+1979
	  days(2)=28
	  s(1)=days(12)+days(1)+days(2)
        if(mod(iys,4)==0.and.mod(iys,100)/=0)then
	   S(1)=S(1)+1
	   days(2)=29
	  elseif(mod(iys,400)==0)then
	   S(1)=S(1)+1
	   days(2)=29
	   endif
	 ! 2D	  
        if(iyr==1)then !1990
	  ims=3
	  do im=1,2
	  temp1=temp1+dynm(ik,k,im,ir,iyr)*days(im)
	  enddo
	  dyns(ik,k,1,ir,iyr)=
     +	  (temp1+dynm(ik,k,12,ir,NYY)*days(12))/s(1)
	  temp1=0.0 
	  else
        do im=1,2
	  temp1=temp1+dynm(ik,k,im,ir,iyr)*days(im)
	  enddo
        dyns(ik,k,1,ir,iyr)=
     +	  (temp1+dynm(ik,k,12,ir,iyr-1)*days(12))/s(1)
	  temp1=0.0 
	  endif

       temp=0.0
	  do im =3,5
	  temp=temp*1.0+dynm(ik,k,im,ir,iyr)*days(im)   
	   enddo   
	   dyns(ik,k,2,ir,iyr)=temp/s(2)
       
        temp=0.0
	  do im =6,8
	  temp=temp*1.0+dynm(ik,k,im,ir,iyr)*days(im)   
	  enddo   
	   dyns(ik,k,3,ir,iyr)=temp/s(3)
        
	 temp=0.0
	  do im =9,11
	  temp=temp*1.0+dynm(ik,k,im,ir,iyr)*days(im) 
	  enddo     
	   dyns(ik,k,4,ir,iyr)=temp/s(4)  
	enddo

      enddo
1002  continue
  
      return
	end subroutine
!(((((((((((((((((((((((((((())))))))))))))))))))))))))))))))))))))
      subroutine Year(iyr,nt,for,raw,Yrave,Yrrave,dyn,dyny)
       integer,parameter :: nvr=19
	integer,parameter :: nYY=34
	  real for(nvr,7,1464,4,NYY),raw(17,nvr,1464,4,NYY)
      real Yrave(nvr,7,4,NYY)
      real Yrrave(17,nvr,4,NYY)
      integer iyr,nirr
	real dyn(17,4,1464,4,NYY) ! 1-3 vor   4 DIV
      real dynd(17,4,1464/4,4,NYY)
	real dynm(17,4,12,4,NYY)
	real dyns(17,4,4,4,NYY)
	real dyny(17,4,4,NYY)
      nirr=3
	do i=1,nirr
	do k=1,nvr

      do j=1,7
	temp=0.0
	do it=1,nt
	 temp=temp*1.0+for(k,j,it,i,iyr)
	enddo
      Yrave(k,j,i,iyr)=temp/nt
	temp=0.0
	enddo

	do ik=1,17
	temp=0.0
	do it=1,nt
	 temp=temp*1.0+raw(ik,k,it,i,iyr)
	enddo
      Yrrave(ik,k,i,iyr)=temp/nt
	temp=0.0
	enddo
	enddo


      do k= 1,4
	do ik=1,17
	temp=0.0
	do it=1,nt
	 temp=temp*1.0+dyn(ik,k,it,i,iyr)
	enddo
      dyny(ik,k,i,iyr)=temp/nt
	temp=0.0
	enddo
      enddo

	enddo

	return

	end subroutine
       
      subroutine  Myrave(Dayly,XMonth,Seanson,Yrave,
     + Dayr,XmonthR,Yrrave,Seasonr,
     + Xmdayly,XMMonth, XMSeanson,Xmdayr,XMSeasonr,XMMonthr,
     + dynd,dynm,dyns,xdynd,xdynm,xdyns)
      integer,parameter :: nrec=1464
      integer,parameter :: nvr=19	
	integer,parameter :: nYY=34  
	COMMON/DD/days(12)
	real Dayly(nvr,7,nrec/4,4,NYY), XMonth(nvr,7,12,4,NYY),
     +     Seanson(nvr,7,4,4,NYY),
     +     Yrave(nvr,7,4,NYY)
	real Dayr(17,nvr,nrec/4,4,NYY),XmonthR(17,nvr,12,4,NYY),
     +  Yrrave(17,nvr,4,NYY),Seasonr(17,nvr,4,4,NYY)
      real Xmdayly(nvr,7,nrec/4,4),XMMonth(nvr,7,12,4),
     + XMSeanson(nvr,7,4,4),Xmdayr(17,nvr,nrec/4,4),
     + XMSeasonr(17,nvr,4,4),XMMonthr(17,nvr,12,4)
	real dyn(17,4,1464,4,NYY) ! 1-3 vor   4 DIV
      real dynd(17,4,1464/4,4,NYY)
	real dynm(17,4,12,4,NYY)
	real dyns(17,4,4,4,NYY)
	real dyny(17,4,4,NYY)
	real xdynd(17,4,1464/4,4)
	real xdynm(17,4,12,4)
	real xdyns(17,4,4,4)
      integer nirr
!
!      print*,XMonth(1,2,1,1,1),'BBBBBBBBBBB'
      nirr=3
	do 1001 i=1,nirr
      do 1001 iv=1,nvr
        
      do j=1,7

!--------dayly---------------
	do it=1,nrec/4
	temp=0.0
      do 1002 iyr=1,NYY	 
       iys=iyr-1+1979
        if(mod(iys,4)==0.and.mod(iys,100)/=0)then
	    nt=1464
	  elseif(mod(iys,400)==0)then
	   nt=1464
	   endif
	if(it==366.and.nt/=1464)goto 1002
         temp=temp+Dayly(iv,j,it,i,iyr)
1002	continue
      Xmdayly(iv,j,it,i)=temp/(NYY*1.0)
	temp=0.0
	enddo
!------- month-----
      do im=1,12
	temp=0.0
      do iyr=1,NYY
	 temp=temp+XMonth(iv,j,im,i,iyr)
	enddo
      Xmmonth(iv,j,im,i)=temp/(NYY*1.0)
	temp=0.0
	enddo
!--------seasons
      do is=1,4
	temp=0.0
      do iyr=1,NYY
	 temp=temp+seanson(iv,j,is,i,iyr)
	enddo
      Xmseanson(iv,j,is,i)=temp/(NYY*1.0)
	temp=0.0
	enddo

	enddo
!**************************************************
      do ik=1,17

!--------dayly---------------
	do it=1,nrec/4
	temp=0.0
      do 1003 iyr=1,NYY	 
       iys=iyr-1+1979
        if(mod(iys,4)==0.and.mod(iys,100)/=0)then
	    nt=1464
	  elseif(mod(iys,400)==0)then
	   nt=1464
	   endif
	if(it==366.and.nt/=1464)goto 1003
         temp=temp+Dayr(ik,iv,it,i,iyr)
1003	continue
      Xmdayr(ik,iv,it,i)=temp/(NYY*1.0)
	temp=0.0
	enddo
!------- month-----
      do im=1,12
	temp=0.0
      do iyr=1,NYY
	 temp=temp+XMonthr(ik,iv,im,i,iyr)
	enddo
      Xmmonthr(ik,iv,im,i)=temp/(NYY*1.0)
	temp=0.0
	enddo
!--------seasons
      do is=1,4
	temp=0.0
      do iyr=1,NYY
	 temp=temp+seasonr(ik,iv,is,i,iyr)
	enddo
      Xmseasonr(ik,iv,is,i)=temp/(NYY*1.0)
	temp=0.0
	enddo

	enddo
1001  continue
!--------------------------------------------
	do 2001 i=1,nirr
      do 2001 iv=1,4

!**************************************************
      do ik=1,17

!--------dayly---------------
	do it=1,nrec/4
	temp=0.0
      do 2003 iyr=1,NYY	 
       iys=iyr-1+1979
        if(mod(iys,4)==0.and.mod(iys,100)/=0)then
	    nt=1464
	  elseif(mod(iys,400)==0)then
	   nt=1464
	   endif
	if(it==366.and.nt/=1464)goto 2003
         temp=temp+dynd(ik,iv,it,i,iyr)
2003	continue
      Xdynd(ik,iv,it,i)=temp/(NYY*1.0)
	temp=0.0
	enddo
!------- month-----
      do im=1,12
	temp=0.0
      do iyr=1,NYY
	 temp=temp+dynm(ik,iv,im,i,iyr)
	enddo
      Xdynm(ik,iv,im,i)=temp/(NYY*1.0)
	temp=0.0
	enddo
!--------seasons
      do is=1,4
	temp=0.0
      do iyr=1,NYY
	 temp=temp+dyns(ik,iv,is,i,iyr)
	enddo
      Xdyns(ik,iv,is,i)=temp/(NYY*1.0)
	temp=0.0
	enddo

	enddo
2001  continue

      return
	end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine output_bin(taget,raw,Dayly,XMonth,Seanson,Yrave,
     + Dayr,XmonthR,Yrrave,Seansonr,
     + Xmdayly,XMMonth, XMSeanson,Xmdayr,XMSeasonr,XMMonthr,dyn,
     + dynd,dynm,dyns,xdynd,xdynm,xdyns)
      integer,parameter :: nrec=1464
	integer,parameter :: nvr=19
	integer,parameter :: nYY=34
	real Dayly(nvr,7,nrec/4,4,NYY), XMonth(nvr,7,12,4,NYY),
     +     Seanson(nvr,7,4,4,NYY),
     +     Yrave(nvr,7,4,NYY)
	real Dayr(17,nvr,nrec/4,4,NYY),XmonthR(17,nvr,12,4,NYY),
     +   Yrrave(17,nvr,4,NYY),Seansonr(17,nvr,4,4,NYY)
      real Xmdayly(nvr,7,nrec/4,4),XMMonth(nvr,7,12,4),
     + XMSeanson(nvr,7,4,4),Xmdayr(17,nvr,nrec/4,4),
     + XMSeasonr(17,nvr,4,4),XMMonthr(17,nvr,12,4)
	character*100 filename,dir,filectl
	character*10 dtime
	real temp(20,20,4,1464*NYY/4),temp1(20,20,4,1464*NYY)
	real taget(nvr,7,nrec,4,NYY),raw(17,nvr,1464,4,NYY)
	real dyn(17,4,1464,4,NYY) ! 1-3 vor   4 DIV
      real dynd(17,4,1464/4,4,NYY)
	real dynm(17,4,12,4,NYY)
	real dyns(17,4,4,4,NYY)
	real dyny(17,4,4,NYY)
	real xdynd(17,4,1464/4,4)
	real xdynm(17,4,12,4)
	real xdyns(17,4,4,4)
	integer nirr
	nirr=3
!--------------------------------------------------------------
      dir='Z:\DATA\LargeScale\79-12\Grads\TP\'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      filename=trim(dir)//'dayly.dat'
      open(10,file=filename,form='binary')  
	filename=trim(dir)//'month.dat'
	open(20,file=filename,form='binary')
      dtime='1mo'
	filectl=trim(filename)//'.ctl'   
	call ctl(nirr,6,1,nvr, 12*NYY,filename,dtime,filectl)
	filename=trim(dir)//'Season.dat'
	open(30,file=filename,form='binary')
	dtime='3mo'
	filectl=trim(filename)//'.ctl'   
	call ctl(nirr,6,1,nvr,  4*NYY,filename,dtime,filectl)
	temp=0.0
	do 101 i=1,nvr
	do 101 j=2,7
	do 101 ir=1,nirr
	id=1
	do 101 iyr=1,NYY
      iys=iyr-1+1979
	nt=1460
        if(mod(iys,4)==0.and.mod(iys,100)/=0)then
	    nt=1464
	  elseif(mod(iys,400)==0)then
	   nt=1464
	   endif
	do 101 it=1,nt/4
      temp(i,j,ir,id)=Dayly(i,j,it,ir,iyr)
	id=id+1
101	continue

      id=id-1
      filename=trim(dir)//'dayly.dat'
	dtime='1dy'
	filectl=trim(filename)//'.ctl' 
	print*,filectl
	print*,trim(filename)   
	call ctl(4,6,1,nvr, id,filename,dtime,filectl)
	do 102 it=1,id
	do 102 ir=1,nvr
	do 102 iv=2,7
	do 102 k=1,nirr
	write(10)temp(ir,iv,k,it)
102   continue
	temp=0.0
!---month
      temp=0.0
!      print*,XMonth(1,2,1,1,1)
!	print*,XMonth(1,2,2,1,1)
	do 103 i=1,nvr
	do 103 j=2,7
	do 103 ir=1,nirr
	id=1
	do 103 iyr=1,NYY
	do 103 it=1,12
      temp(i,j,ir,id)=Xmonth(i,j,it,ir,iyr)
	id=id+1
103   continue
	do 104 it=1,NYY*12
	do 104 ir=1,nvr
	do 104 iv=2,7
	do 104 k=1,nirr
	write(20)temp(ir,iv,k,it)
104   continue
	temp=0.0
!---Season
      temp=0.0
	do 105 i=1,nvr
	do 105 j=2,7
	do 105 ir=1,nirr
	id=1
	do 105 iyr=1,NYY
	do 105 it=1,4
      temp(i,j,ir,id)=seanson(i,j,it,ir,iyr)
	id=id+1
105   continue
	do 106 it=1,NYY*4
	do 106 ir=1,nvr
	do 106 iv=2,7
	do 106 k=1,nirr
	write(30)temp(ir,iv,k,it)
106   continue
	close(10)
	close(20)
	close(30)
	temp=0.0
      filename=trim(dir)//'year.dat'
	open(30,file=filename,form='binary')
	dtime='1yr'
	filectl=trim(filename)//'.ctl'   
	call ctl(nirr,6,1,nvr,  NYY,filename,dtime,filectl)
	do 107 it=1,NYY
	do 107 ir=1,nvr
	do 107 iv=2,7
	do 107 k=1,nirr
	write(30)Yrave(ir,iv,k,it)
107   continue
      close(30)

 
      temp=0.0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      filename=trim(dir)//'day-levels.dat'
      open(10,file=filename,form='binary')
	
	filename=trim(dir)//'month-levels.dat'
	open(20,file=filename,form='binary')
	dtime='1mo'
	filectl=trim(filename)//'.ctl'   
	call ctl(nirr,1,17,nvr, NYY*12,filename,dtime,filectl)
	filename=trim(dir)//'season-levels.dat'
	open(30,file=filename,form='binary')
	dtime='3mo'
	filectl=trim(filename)//'.ctl'   
	call ctl(nirr,1,17,nvr, NYY*4,filename,dtime,filectl)
	temp=0.0
	do 201 i=1,nvr
	do 201 ik=1,17
	do 201 ir=1,nirr
	id=1
	do 201 iyr=1,NYY
      iys=iyr-1+1979
	nt=1460
        if(mod(iys,4)==0.and.mod(iys,100)/=0)then
	    nt=1464
	  elseif(mod(iys,400)==0)then
	   nt=1464
	   endif
	do 201 it=1,nt/4
      temp(ik,i,ir,id)=Dayr(ik,i,it,ir,iyr)
	id=id+1
201   continue
      id=id-1
      filename=trim(dir)//'day-levels.dat'
	dtime='1dy'
	filectl=trim(filename)//'.ctl'   
	call ctl(nirr,1,17,nvr, id,filename,dtime,filectl)

	do 202 it=1,id ! time
	do 202 ir=1,nvr  ! variables
	do 202 ik=1,17   ! ver
	do 202 k=1,nirr    ! areas
	write(10)temp(ik,ir,k,it)
202   continue
	temp=0.0
!---month
      temp=0.0
	do 203 i=1,nvr
	do 203 ik=1,17
	do 203 ir=1,nirr
	id=1
	do 203 iyr=1,NYY
	do 203 it=1,12
      temp(ik,i,ir,id)=Xmonthr(ik,i,it,ir,iyr)
	id=id+1
203   continue
	do 204 it=1,NYY*12
	do 204 ir=1,nvr
	do 204 ik=1,17
	do 204 k=1,nirr
	write(20)temp(ik,ir,k,it)
204 	continue 
	temp=0.0
!---Season
      temp=0.0
	do 205 i=1,nvr
	do 205 ik=1,17
	do 205 ir=1,nirr
	id=1
	do 205 iyr=1,NYY
	do 205 it=1,4
      temp(ik,i,ir,id)=seansonr(ik,i,it,ir,iyr)
	id=id+1
205   continue
	do 206 it=1,NYY*4
	do 206 ir=1,nvr
	do 206 ik=1,17
	do 206 k=1,nirr
	write(30)temp(ik,ir,k,it)
206 	continue 
	close(10)
	close(20)
	close(30)
	temp=0.0
	filename=trim(dir)//'yrrave-levels.dat'
	open(30,file=filename,form='binary')
	dtime='1yr'
	filectl=trim(filename)//'.ctl'   
	call ctl(nirr,1,17,nvr, NYY,filename,dtime,filectl)
	do 207 it=1,NYY
	do 207 ir=1,nvr
	do 207 ik=1,17
	do 207 k=1,nirr
	write(30)Yrrave(ik,ir,k,it)
207 	continue 
      close(30)
	temp=0.0
!----------------SURFACE------------------------------
      filename=trim(dir)//'dayly_S.dat'
      open(10,file=filename,form='binary')  
	
	filename=trim(dir)//'month_S.dat'
	open(20,file=filename,form='binary')
      dtime='1mo'
	filectl=trim(filename)//'.ctl'   
	call ctl(nirr,1,1,6, 12*NYY,filename,dtime,filectl)
	filename=trim(dir)//'Season_S.dat'
	open(30,file=filename,form='binary')
	dtime='3mo'
	filectl=trim(filename)//'.ctl'   
	call ctl(nirr,1,1,6,  4*NYY,filename,dtime,filectl)
	filename=trim(dir)//'year_S.dat'
	open(40,file=filename,form='binary')
	dtime='1yr'
	filectl=trim(filename)//'.ctl'   
	call ctl(nirr,1,1,6,  NYY,filename,dtime,filectl)
	      temp=0.0
      do 301 i=1,6
	do 301 ir=1,nirr
	id=1
	do 301 iyr=1,NYY
	iys=iyr-1+1979
	nt=1460
        if(mod(iys,4)==0.and.mod(iys,100)/=0)then
	    nt=1464
	  elseif(mod(iys,400)==0)then
	   nt=1464
	   endif
	do 301 it=1,nt/4
      temp(i,1,ir,id)=dayly(i,1,it,ir,iyr) !Dayly(i,j,it,ir,iyr)
	id=id+1
301   continue
      id=id-1
      filename=trim(dir)//'dayly_S.dat'
      dtime='1dy' 
	filectl=trim(filename)//'.ctl'      
	call ctl(nirr,1,1,6, id,filename,dtime,filectl)
      do 302 it=1,id
	do 302 ir=1,6
	do 302 k=1,nirr
	write(10)temp(ir,1,k,it)
302   continue
	temp=0.0
	do 303 i=1,6
	do 303 ir=1,nirr
	id=1
	do 303 iyr=1,NYY
	do 303 it=1,12
      temp(i,1,ir,id)=Xmonth(i,1,it,ir,iyr) !Dayly(i,j,it,ir,iyr)
	id=id+1
303   continue
      do 304 it=1,NYY*12
	do 304 ir=1,6
	do 304 k=1,nirr
	write(20)temp(ir,1,k,it)
304   continue
	temp=0.0
	do 305 i=1,6
	do 305 ir=1,nirr
	id=1
	do 305 iyr=1,NYY
	do 305 it=1,4
      temp(i,1,ir,id)=seanson(i,1,it,ir,iyr) !Dayly(i,j,it,ir,iyr)
	id=id+1
305   continue
      do 306 it=1,NYY*4
	do 306 ir=1,6
	do 306 k=1,nirr
	write(30)temp(ir,1,k,it)
306   continue
	temp=0.0
	close(10)
	close(20)
	close(30)
	temp=0.0
	do 307 it=1,NYY
	do 307 ir=1,nvr
	do 307 k=1,nirr
	write(40)Yrave(ir,i,k,it)
307   continue
      close(40)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ------------------- Multiyear ave output --------------------------
!  COMMON/AVE3/mdayly(18,7,nrec/4,4),MMonth(18,7,12,4),
!     +   MSeanson(18,7,4,4),mdayr(17,18,nrec/4,4),MMonthr(17,18,12,4),
!     +   MSeansonr(17,18,4,4)    
      filename=trim(dir)//'Multi_dayly.dat'
	open(10,file=filename,form='binary')  
	dtime='1dy'
	filectl=trim(filename)//'.ctl'       
	call ctl(nirr,6,1,nvr, 365,filename,dtime,filectl)
	filename=trim(dir)//'Multi_month.dat'
	open(20,file=filename,form='binary')  
	dtime='1mo'    
	filectl=trim(filename)//'.ctl'   
	call ctl(nirr,6,1,nvr, 12,filename,dtime,filectl)
	filename=trim(dir)//'Multi_season.dat'
	open(30,file=filename,form='binary')  
	dtime='3mo'    
	filectl=trim(filename)//'.ctl'   
	call ctl(nirr,6,1,nvr, 4,filename,dtime,filectl)
      do 401 it=1,365
	do 401 ir=1,nvr
	do 401 j=2,7
	do 401 k=1,nirr
	write(10)Xmdayly(ir,j,it,k)
401 	continue 
      do 402 it=1,12
	do 402 ir=1,nvr
	do 402 j=2,7
	do 402 k=1,nirr
	write(20)Xmmonth(ir,j,it,k)
402 	continue 
      do 403 it=1,4
	do 403 ir=1,nvr
	do 403 j=2,7
	do 403 k=1,nirr
	write(30)Xmseanson(ir,j,it,k)
403 	continue
      close(10)
	close(20)
	close(30)
      filename=trim(dir)//'Multi_dayly_Levs.dat'
	open(10,file=filename,form='binary')  
	dtime='1dy' 
	filectl=trim(filename)//'.ctl'      
	call ctl(nirr,1,17,nvr, 365,filename,dtime,filectl)
	filename=trim(dir)//'Multi_month_Levs.dat'
	open(20,file=filename,form='binary')  
	dtime='1mo'    
	filectl=trim(filename)//'.ctl'   
	call ctl(nirr,1,17,nvr, 12,filename,dtime,filectl)
	filename=trim(dir)//'Multi_season_Levs.dat'
	open(30,file=filename,form='binary')  
	dtime='3mo'    
	filectl=trim(filename)//'.ctl'   
	call ctl(nirr,1,17,nvr, 4,filename,dtime,filectl)
      do 501 it=1,365
	do 501 ir=1,nvr
	do 501 ik=1,17
	do 501 k=1,nirr
	write(10)Xmdayr(ik,ir,it,k)
501 	continue 
      do 502 it=1,12
	do 502 ir=1,nvr
	do 502 ik=1,17
	do 502 k=1,nirr
	write(20)Xmmonthr(ik,ir,it,k)
502 	continue 
      do 503 it=1,4
	do 503 ir=1,nvr
	do 503 ik=1,17
	do 503 k=1,nirr
	write(30)Xmseasonr(ik,ir,it,k)
503 	continue 
      close(10)
	close(20)
	close(30)
	temp=0.0
!----------------SURFACE------------------------------
      filename=trim(dir)//'Multi_dayly_S.dat'
      open(10,file=filename,form='binary')  
	dtime='1dy'    
	filectl=trim(filename)//'.ctl'
	print*,filectl  
	call ctl(nirr,1,1,6, 365,filename,dtime,filectl)
	filename=trim(dir)//'Multi_month_S.dat'
	open(20,file=filename,form='binary')
      dtime='1mo'
	filectl=trim(filename)//'.ctl'   
	call ctl(nirr,1,1,6, 12,filename,dtime,filectl)
	filename=trim(dir)//'Multi_Season_S.dat'
	open(30,file=filename,form='binary')
	dtime='3mo'
	filectl=trim(filename)//'.ctl'   
	call ctl(nirr,1,1,6,  4,filename,dtime,filectl)
      do 601 it=1,365
	do 601 ir=1,6
	do 601 k=1,nirr
	write(10)Xmdayly(ir,1,it,k)
601 	continue 
      do 602 it=1,12
	do 602 ir=1,6
	do 602 k=1,nirr
	write(20)Xmmonth(ir,1,it,k)
602 	continue
      do 603 it=1,4
	do 603 ir=1,6
	do 603 k=1,nirr
	write(30)Xmseanson(ir,1,it,k)
603 	continue 
      close(10)
	close(20)
	close(30) 
	temp1=0.0
!--------------RAW output------------------------------------
!real taget(18,7,nrec,4,23),raw(17,18,1464,4,23) temp1(20,20,4,1464*23)
      filename=trim(dir)//'RAW_dlevs.dat'
      open(10,file=filename,form='binary')  
	dtime='6hr'    
	filectl=trim(filename)//'.ctl'
	print*,filectl  
      do 701 i=1,nirr
	do 701 j=2,7
	do 701 iv=1,nvr
	id=1
	do 701 iyr=1,NYY
	iys=iyr-1+1979
	nt=1460
        if(mod(iys,4)==0.and.mod(iys,100)/=0)then
	    nt=1464
	  elseif(mod(iys,400)==0)then
	   nt=1464
	  endif
	do 701 it=1,nt
      temp1(iv,j,i,id)=taget(iv,j,it,i,iyr)
	id=id+1
 701  continue
      id=id-1
	call ctl(nirr,6,1,nvr, id,filename,dtime,filectl)
	do 702 it=1,id
      do 702 iv=1,nvr
	do 702 j=2,7
	do 702 i=1,nirr
	write(10)temp1(iv,j,i,it)
702   continue
      temp1=0.0
!  raw(17,18,1464,4,23) temp1(20,20,4,1464*23)
      filename=trim(dir)//'RAW_all_levs.dat'
      open(20,file=filename,form='binary')  
	dtime='6hr'    
	filectl=trim(filename)//'.ctl'
	print*,filectl  
      do 703 i=1,nirr
	do 703 ik=1,17
	do 703 iv=1,nvr
	id=1
	do 703 iyr=1,NYY
	iys=iyr-1+1979
	nt=1460
        if(mod(iys,4)==0.and.mod(iys,100)/=0)then
	    nt=1464
	  elseif(mod(iys,400)==0)then
	   nt=1464
	  endif
	do 703 it=1,nt
      temp1(ik,iv,i,id)=raw(ik,iv,it,i,iyr)
	id=id+1
 703  continue
      id=id-1
	call ctl(nirr,1,17,nvr, id,filename,dtime,filectl)
	do 704 it=1,id
      do 704 iv=1,nvr
	do 704 ik=1,17
	do 704 i=1,nirr
	write(20)temp1(ik,iv,i,it)
704   continue
!---------SURFACE-----------------------------------------
      temp1=0.0
!real taget(18,7,nrec,4,23),raw(17,18,1464,4,23) temp1(20,20,4,1464*23)
      filename=trim(dir)//'RAW_surface.dat'
      open(30,file=filename,form='binary')  
	dtime='6hr'    
	filectl=trim(filename)//'.ctl'
	print*,filectl  
      do 705 i=1,nirr
	do 705 iv=1,6
	id=1
	do 705 iyr=1,NYY
	iys=iyr-1+1979
	nt=1460
        if(mod(iys,4)==0.and.mod(iys,100)/=0)then
	    nt=1464
	  elseif(mod(iys,400)==0)then
	   nt=1464
	  endif
	do 705 it=1,nt
      temp1(iv,1,i,id)=taget(iv,1,it,i,iyr)
	id=id+1
 705  continue
      id=id-1
	call ctl(nirr,1,1,6, id,filename,dtime,filectl)
	do 706 it=1,id
      do 706 iv=1,6
	do 706 i=1,nirr
	write(30)temp1(iv,1,i,it)
706   continue
      close(10)
	close(20)
	close(30) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 	do 800 i=1,nirr
	do 800 ik=1,17
	do 800 ir=1,4
	id=1
	do 800 iyr=1,NYY
      iys=iyr-1+1979
	nt=1460
        if(mod(iys,4)==0.and.mod(iys,100)/=0)then
	    nt=1464
	  elseif(mod(iys,400)==0)then
	   nt=1464
	   endif
	do 800 it=1,nt/4
      temp(ik,i,ir,id)=dynd(ik,i,it,ir,iyr)
	id=id+1
800   continue
      id=id-1
      filename=trim(dir)//'day-dyn.dat'
      open(10,file=filename,form='binary')  
	dtime='1dy'
	filectl=trim(filename)//'.ctl'   
	call ctl_dyn(nirr,1,17,4, id,filename,dtime,filectl)

	do 801 it=1,id ! time
	do 801 ir=1,4  ! variables
	do 801 ik=1,17   ! ver
	do 801 k=1,nirr    ! areas
	write(10)temp(ik,ir,k,it)
801   continue
      close(10)
	temp=0.0
!  dyn  month
	do 810 i=1,nirr
	do 810 ik=1,17
	do 810 ir=1,4
	id=1
	do 810 iyr=1,NYY
	do 810 it=1,12
      temp(ik,i,ir,id)=Dayr(ik,i,it,ir,iyr)
	id=id+1
810   continue
      id=id-1
      filename=trim(dir)//'month-dyn.dat'
	open(10,file=filename,form='binary')  
	dtime='1mo'
	filectl=trim(filename)//'.ctl'   
	call ctl_dyn(nirr,1,17,4, NYY*12,filename,dtime,filectl)

	do 811 it=1,NYY*12 ! time
	do 811 ir=1,4  ! variables
	do 811 ik=1,17   ! ver
	do 811 k=1,nirr    ! areas
	write(10)temp(ik,ir,k,it)
811   continue
      close(10)
	temp=0.0
!  season  month
	do 820 i=1,nirr
	do 820 ik=1,17
	do 820 ir=1,4
	id=1
	do 820 iyr=1,NYY
	do 820 it=1,4
      temp(ik,i,ir,id)=Dayr(ik,i,it,ir,iyr)
	id=id+1
820   continue
      id=id-1
      filename=trim(dir)//'season-dyn.dat'
	open(10,file=filename,form='binary')  
	dtime='3mo'
	filectl=trim(filename)//'.ctl'   
	call ctl_dyn(nirr,1,17,4, NYY*4,filename,dtime,filectl)

	do 821 it=1,NYY*4 ! time
	do 821 ir=1,4  ! variables
	do 821 ik=1,17   ! ver
	do 821 k=1,nirr    ! areas
	write(10)temp(ik,ir,k,it)
821   continue
      close(10)
	temp=0.0
	filename=trim(dir)//'Multi_dayly_dyn.dat'
      open(10,file=filename,form='binary')  
	dtime='1dy'    
	filectl=trim(filename)//'.ctl'
	print*,filectl  
	call ctl_dyn(nirr,1,1,6, 365,filename,dtime,filectl)
	do 830 it=1,365 ! time
	do 830 ir=1,4  ! variables
	do 830 ik=1,17   ! ver
	do 830 k=1,nirr    ! areas
	write(10)xdynd(ik,ir,it,k)
830   continue
      close(10)
	temp=0.0
	filename=trim(dir)//'Multi_month_dyn.dat'
      open(10,file=filename,form='binary')  
	dtime='1dy'    
	filectl=trim(filename)//'.ctl'
	print*,filectl  
	call ctl_dyn(nirr,1,1,6, 12,filename,dtime,filectl)
	do 840 it=1,12 ! time
	do 840 ir=1,4  ! variables
	do 840 ik=1,17   ! ver
	do 840 k=1,nirr    ! areas
	write(10)xdynm(ik,ir,it,k)
840   continue
      close(10)
	temp=0.0
      filename=trim(dir)//'Multi_seanson_dyn.dat'
      open(10,file=filename,form='binary')  
	dtime='1dy'    
	filectl=trim(filename)//'.ctl'
	print*,filectl  
	call ctl_dyn(nirr,1,1,6, 4,filename,dtime,filectl)
	do 850 it=1,4 ! time
	do 850 ir=1,4  ! variables
	do 850 ik=1,17   ! ver
	do 850 k=1,nirr    ! areas
	write(10)xdyns(ik,ir,it,k)
850   continue
      close(10)
	temp=0.0
      filename=trim(dir)//'RAW_dyn.dat'
      open(10,file=filename,form='binary')  
	dtime='6hr'    
	filectl=trim(filename)//'.ctl'
	print*,filectl  
      do 860 i=1,nirr
	do 860 ik=1,17
	do 860 iv=1,4
	id=1
	do 860 iyr=1,NYY
	iys=iyr-1+1979
	nt=1460
        if(mod(iys,4)==0.and.mod(iys,100)/=0)then
	    nt=1464
	  elseif(mod(iys,400)==0)then
	   nt=1464
	  endif
	do 860 it=1,nt
      temp1(ik,iv,i,id)=dyn(ik,iv,it,i,iyr)
	id=id+1
 860  continue
      id=id-1
	call ctl_dyn(nirr,1,17,4, id,filename,dtime,filectl)
	do 861 it=1,id
      do 861 iv=1,4
	do 861 ik=1,17
	do 861 i=1,nirr
	write(20)temp1(ik,iv,i,it)
861   continue


      return
	end subroutine
      
      subroutine ctl(ix,iy,iz,iv,it,filename,dtime,filectl)
	character*100 filename,filectl
	integer,parameter :: nvr=19
	character*40 var(nvr),varn(nvr),SVAR(6),SVARN(6)
	character*10 dtime,typ
	real plv(17)
	integer days(12)
	data plv/1000.,925., 850.,700.,
     +     600., 500.,400.,300.,250.,200.,150.,100.,70,50,30,20,10/
     
      var(1)='TLS';varn(1)='Temp.forcing K/d'
	var(2)='QLS';varn(2)='Moiture.forcing K/d'
	var(3)='uwnd';varn(3)='u wind m/s'
	var(4)='vwnd';varn(4)='v.wind m/s'
	var(5)='qv';varn(5)='Moiture kg/kg'
	var(6)='hgt';varn(6)='hgt m'
	var(7)='T';varn(7)='Temp.K'
	var(8)='adjomg';varn(8)='adjust omega pa/s'
	var(9)='rh';varn(9)='RH %'
	var(10)='the';varn(10)='Theta K'
	var(11)='Q1';varn(11)='Q1 K/d'
	var(12)='Q2';varn(12)='Q2 K/d'
	var(13)='HADQ';varn(13)='HADQ K/d'
	var(14)='VADQ';varn(14)='VADQ K/d'
	var(15)='TCHQ';varn(15)='TCHQ K/d'
	var(16)='HADT';varn(16)='HADT K/d'
	var(17)='VADT';varn(17)='VADT K/d'
	var(18)='TCHT';varn(18)='TCHT K/d'
	var(19)='ncepomg';varn(19)='NCEP omega pa/s'
!=================================================
	SVAR(1)='LH';SVARN(1)='laten heat W/m^2'
	SVAR(2)='SH';SVARN(2)='sensible heat W/m^2'
	SVAR(3)='pr';SVARN(3)='raifall rate kg/m^2/s'
	SVAR(4)='cpr';SVARN(4)='convetive raifall rate kg/m^2/s'
	SVAR(5)='DLR';SVARN(5)='down longwave radiation W/m^2'
	SVAR(6)='ULR';SVARN(6)='Up longwave radiation W/m^2'
!--------------------------------------------------
      typ='.ctl'
!	print*,trim(filename)
!	filectl=trim(filename)//trim(typ)
!	print*,trim(filename)
	open(1,file=filectl)
      write(1,*)'dset',' ',trim(filename)
	write(1,*)'title Forcing_ave'
	write(1,*)'options little_endian'
	write(1,*)'undef -9999.'
      write(1,90)'xdef',ix, 'linear 1 1'
	write(1,90)'ydef',iy, 'linear 1 1'
	if(iz==17)then
	write(1,91)'zdef',iz, 'levels', (plv(ik),ik=1,iz)
	else
	write(1,911)'zdef',iz, 'levels', (plv(ik),ik=1,iz)
	endif
      write(1,92)'tdef', it, 'linear 00z1jan1979',trim(dtime)
	write(1,*)'vars',iv
	if(iv==nvr)then
	do i=1,iv
	write(1,93)trim(var(i)),iz,'99', trim(varn(i))
	enddo
	endif
	if(iv==6)then
      do i=1,iv
	write(1,93)trim(SVAR(i)),iz,'99', trim(SVARN(i))
	enddo
	endif
	write(1,*)'endvars'
90    format(1X,A4,1x,I2,1X,A10)
91    format(1X,A4,1x,I2,1X,A6,17(1X,F9.3))
911   format(1X,A4,1x,I2,1X,A6,1X,F9.3)
92    format(1X,A4,1x,I6,1X,A19,1X,A4)
93    format(1X,A6,1x,I2,1X,A2,1X,A35)
      close(1)
!-----------------------------------------

	return
	end subroutine

	subroutine ctl_dyn(ix,iy,iz,iv,it,filename,dtime,filectl)
	character*100 filename,filectl
	integer,parameter :: nvr=4
	character*40 var(nvr),varn(nvr),SVAR(6),SVARN(6)
	character*10 dtime,typ
	real plv(17)
	integer days(12)
	data plv/1000.,925., 850.,700.,
     +     600., 500.,400.,300.,250.,200.,150.,100.,70,50,30,20,10/
     
      var(1)='vorx';varn(1)='vorticity X'
	var(2)='vory';varn(2)='vorticity Y'
	var(3)='vorz';varn(3)='vorticity Z'
	var(4)='div';varn(4)='Divergence'
!	var(5)='qv';varn(5)='Moiture kg/kg'
!	var(6)='hgt';varn(6)='hgt m'
!	var(7)='T';varn(7)='Temp.K'
!	var(8)='adjomg';varn(8)='adjust omega pa/s'
!	var(9)='rh';varn(9)='RH %'
!	var(10)='the';varn(10)='Theta K'
!	var(11)='Q1';varn(11)='Q1 K/d'
!	var(12)='Q2';varn(12)='Q2 K/d'
!	var(13)='HADQ';varn(13)='HADQ K/d'
!	var(14)='VADQ';varn(14)='VADQ K/d'
!	var(15)='TCHQ';varn(15)='TCHQ K/d'
!	var(16)='HADT';varn(16)='HADT K/d'
!	var(17)='VADT';varn(17)='VADT K/d'
!	var(18)='TCHT';varn(18)='TCHT K/d'
!	var(19)='ncepomg';varn(19)='NCEP omega pa/s'
!=================================================
	SVAR(1)='LH';SVARN(1)='laten heat W/m^2'
	SVAR(2)='SH';SVARN(2)='sensible heat W/m^2'
	SVAR(3)='pr';SVARN(3)='raifall rate kg/m^2/s'
	SVAR(4)='cpr';SVARN(4)='convetive raifall rate kg/m^2/s'
	SVAR(5)='DLR';SVARN(5)='down longwave radiation W/m^2'
	SVAR(6)='ULR';SVARN(6)='Up longwave radiation W/m^2'
!--------------------------------------------------
      typ='.ctl'
!	print*,trim(filename)
!	filectl=trim(filename)//trim(typ)
!	print*,trim(filename)
	open(1,file=filectl)
      write(1,*)'dset',' ',trim(filename)
	write(1,*)'title Forcing_ave'
	write(1,*)'options little_endian'
	write(1,*)'undef -9999.'
      write(1,90)'xdef',ix, 'linear 1 1'
	write(1,90)'ydef',iy, 'linear 1 1'
	if(iz==17)then
	write(1,91)'zdef',iz, 'levels', (plv(ik),ik=1,iz)
	else
	write(1,911)'zdef',iz, 'levels', (plv(ik),ik=1,iz)
	endif
      write(1,92)'tdef', it, 'linear 00z1jan1979',trim(dtime)
	write(1,*)'vars',iv
	if(iv==nvr)then
	do i=1,iv
	write(1,93)trim(var(i)),iz,'99', trim(varn(i))
	enddo
	endif
	if(iv==6)then
      do i=1,iv
	write(1,93)trim(SVAR(i)),iz,'99', trim(SVARN(i))
	enddo
	endif
	write(1,*)'endvars'
90    format(1X,A4,1x,I2,1X,A10)
91    format(1X,A4,1x,I2,1X,A6,17(1X,F9.3))
911   format(1X,A4,1x,I2,1X,A6,1X,F9.3)
92    format(1X,A4,1x,I6,1X,A19,1X,A4)
93    format(1X,A6,1x,I2,1X,A2,1X,A35)
      close(1)
!-----------------------------------------

	return
	end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine output_ACS(taget,raw,Dayly,XMonth,Seanson,Yrave,
     + Dayr,XmonthR,Yrrave,Seansonr,
     + Xmdayly,XMMonth, XMSeanson,Xmdayr,XMSeasonr,XMMonthr,dyn,
     + dynd,dynm,dyns,xdynd,xdynm,xdyns)
      integer,parameter :: nrec=1464
	integer,parameter :: nvr=19
	integer,parameter :: nYY=34
	real Dayly(nvr,7,nrec/4,4,NYY), XMonth(nvr,7,12,4,NYY),
     +     Seanson(nvr,7,4,4,NYY),
     +     Yrave(nvr,7,4,NYY)
	real Dayr(17,nvr,nrec/4,4,NYY),XmonthR(17,nvr,12,4,NYY),
     +   Yrrave(17,nvr,4,NYY),Seansonr(17,nvr,4,4,NYY)
      real Xmdayly(nvr,7,nrec/4,4),XMMonth(nvr,7,12,4),
     + XMSeanson(nvr,7,4,4),Xmdayr(17,nvr,nrec/4,4),
     + XMSeasonr(17,nvr,4,4),XMMonthr(17,nvr,12,4)
	character*100 filename,dir,filectl
	character*10 dtime
	real temp(20,20,4,1464*NYY/4),temp1(20,20,4,1464*NYY)
	real taget(nvr,7,nrec,4,NYY),raw(17,nvr,1464,4,NYY)
	real dyn(17,4,1464,4,NYY) ! 1-3 vor   4 DIV
      real dynd(17,4,1464/4,4,NYY)
	real dynm(17,4,12,4,NYY)
	real dyns(17,4,4,4,NYY)
	real dyny(17,4,4,NYY)
	real xdynd(17,4,1464/4,4)
	real xdynm(17,4,12,4)
	real xdyns(17,4,4,4)
	integer nirr
	character(len=6) area(4)
      area(1)='WTP'  ;    area(2)='ETP'
	area(3)='BOB'  ;    area(4)='NEC'
	nirr=3
!--------------------------------------------------------------
      dir='Z:\DATA\LargeScale\79-12\ACSII\'
!!!--------JJA of every year--------------------------
!-----Q1 and Q2
      do ir=1,nirr
      filename=trim(dir)//trim(area(ir))//'_Q1Q2_79_12_JJA.txt'
      open(1,file=trim(filename))
	write(1,14)'Year', 'LQ1', 'MQ1', 'HQ1', 'AQ1', 'Q1'
     +	, 'LQ2', 'MQ2', 'HQ2', 'AQ2', 'Q2'
	filename=trim(dir)//trim(area(ir))//'Comps_79_12_JJA.txt'
      open(2,file=trim(filename))
	write(2,13)'Year', 'L_HADQ', 'M_HADQ', 'H_HADQ', 'A_HADQ', 'HADQ'
     +, 'L_VADQ', 'M_VADQ', 'H_VADQ', 'A_VADQ', 'VADQ'
     +, 'L_TADQ', 'M_TADQ', 'H_TADQ', 'A_TADQ', 'TADQ'
     +, 'L_HADT', 'M_HADT', 'H_HADT', 'A_HADT', 'HADT'
     +, 'L_VADT', 'M_VADT', 'H_VADT', 'A_VADT', 'VADT'
     +, 'L_TADT', 'M_TADT', 'H_TADT', 'A_TADT', 'TADT'
	do iy=1,NYY
	write(1,11)iy+1978
     +	,((seanson(i,j,3,ir,iy),j=2,6),i=11,12)
	write(2,12)iy+1978
     +	,((seanson(i,j,3,ir,iy),j=2,6),i=13,18)
	enddo
	enddo
      




11    format(1X,I4,10(1X,F10.5))
12    format(1X,I4,30(1X,F10.5))
13    format(1X,A4,30(1X,A6))
14    format(1X,A4,10(1X,A3))
	end subroutine