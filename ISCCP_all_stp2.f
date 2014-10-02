      program Corrl
	character dir*100,fnm*100,flnm*100
	character area(7)*4,year*4,q1q2nm(10)*3
	character compnm(30)*6,cbnm(164)*8,sr*2
	character cnm(4)*3,cldlev(3)*3
	real cldjja(27,7,9,130),gpc(7,34) !,gpc2(7,27)
      real q1(7,34,5),q2(7,34,5) !, q1c(7,27,5),q2c(7,27,5)
	real comps(7,34,30) !,compsc(7,27,30)
	real combin(27,7,130-7+5+5+30+1)  !! years areas 130-7 cloud 5+5 Q1 Q2 30 comps rain
	real tmp1(27),tmp2(27)
	character Typcld(15)*4
	character arf(1000000)*4,art(1000000)*4
	character varf(1000000)*8,vart(1000000)*8
      area(1)='PRD'  
      area(2)='MLYR' 
      area(3)='NPC' 
      area(4)='NEC' 
      area(5)='BOB'  
      area(6)='ETP'  
      area(7)='WTP'
	cnm(1)='PC'
	cnm(2)='TC'
	cnm(3)='TAU'
      cnm(4)='WP'
	cldlev(1)='LCA'
	cldlev(2)='MCA'
	cldlev(3)='HCA'
	Typcld(1)='Cu_L'
	Typcld(2)='Sc_L'
	Typcld(3)='St_L'
	Typcld(4)='Cu_I'
	Typcld(5)='Sc_I'
	Typcld(6)='St_I'
	Typcld(7)='Ac_L'
	Typcld(8)='As_L'
	Typcld(9)='Nb_L'
	Typcld(10)='Ac_I'
	Typcld(11)='As_I'
	Typcld(12)='Nb_I'
      Typcld(13)='Cr_I'
	Typcld(14)='Cs_I'
	Typcld(15)='Dc_I'

	dir='Z:\DATA\LargeScale\79-12\GPCP\'  
	do iyy=1983,2009  
	iy=iyy-1982   
	write(year,'(I4)')iyy
      fnm=trim(dir)//year//'_All.txt'
	open(999,file=trim(fnm))
	 do it=1,9
 	  do iv=1,130
          read(999,99)(cldjja(iy,ir,it,iv),ir=1,7)
	    do ir=1,7
          if(cldjja(iy,ir,it,iv)==-0.99990000E+04)
     +		cldjja(iy,ir,it,iv)=0.0
	    enddo
	  enddo
	 enddo
	close(999)
      enddo
99    format(1X,7(1X,E16.8))
      do ir=1,7
	fnm='Z:\DATA\LargeScale\79-12\GPCP\'//trim(area(ir))//'.txt'
	open(96,file=trim(fnm))
        do iy=1979,2012
         read(96,98)iyx,gpc(ir,iy-1978)
	  enddo
      close(96)
	enddo
98    format(1X,I4,1X,F8.4)
!      ip1=1983-1978
!	ip2=2009-1978
!      gpc2(:,:)=gpc(:,ip1:ip2)
!-------------------------------------------------------
	dir='Z:\DATA\LargeScale\79-12\ACSII\'
      do ir=1,7
      fnm=trim(dir)//trim(area(ir))//'_Q1Q2_79_12_JJA.txt'
      open(1,file=trim(fnm))
	read(1,13)year,(q1q2nm(i),i=1,10)
	fnm=trim(dir)//trim(area(ir))//'Comps_79_12_JJA.txt'
      open(2,file=trim(fnm))
	read(2,14)year,(compnm(i),i=1,30)
	do iy=1,34
	read(1,11)iyx
     +	,(q1(ir,iy,i),i=1,5),(q2(ir,iy,i),i=1,5)
	read(2,12)iyx
     +	, (comps(ir,iy,i),i=1,30)
	enddo
	enddo
!	ip1=1983-1978
!	ip2=2009-1978
!      q1c(:,:,:)=q1(:,ip1:ip2,:)
!	q2c(:,:,:)=q2(:,ip1:ip2,:)
!	compsc(:,:,:)=comps(:,ip1:ip2,:)	     
11    format(1X,I4,10(1X,F10.5))
12    format(1X,I4,30(1X,F10.5))
14    format(1X,A4,30(1X,A6))
13    format(1X,A4,10(1X,A3))
!------------------------------------------------------------------
! q1c q2c compsc gpc2 cldjja
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      cbnm(1)='CA'
	cbnm(2)='IR_CA'
	do i=10,19
	write(sr,'(I2)')i
	cbnm(i-7)=sr//'%_CA'
	enddo
      do i=1,4
	ii=12+(i-1)*2
      cbnm(i+ii)= cnm(i)
      cbnm(i+ii+1)= trim(cnm(i))//'_sd'
      cbnm(i+ii+2)= trim(cnm(i))//'_td'
      enddo
!!!  IR cloud types  
      k=1                     
      do i=32,40,3
	j=i-7
      cbnm(j)=cldlev(k)
	cbnm(j+1)=cldlev(k)//'_PC'
	cbnm(j+2)=cldlev(k)//'_TC'
!	print*,k,j,'AAAAAAAAAAAA'
      k=k+1
	enddo
	k=1
	do i=41,115,5
	j=i-7
       cbnm(j)=Typcld(k)
	 cbnm(j+1)=Typcld(k)//'_PC'
	 cbnm(j+2)=Typcld(k)//'_TC'
	 cbnm(j+3)=Typcld(k)//'_TAU'
	 cbnm(j+4)=Typcld(k)//'WP'
	k=k+1
      enddo
      cbnm(116-7)='TS_CK'
	cbnm(117-7)='TS_Sd'
	cbnm(118-7)='RS_CK'
	cbnm(119-7)='I-S'
	cbnm(120-7)='PS'
	cbnm(121-7)='TSA'
	cbnm(122-7)='T_740mb'
	cbnm(123-7)='T_500mb'
	cbnm(124-7)='T_375mb'
	cbnm(125-7)='PT'
	cbnm(126-7)='TT'
	cbnm(127-7)='T_50mb'
	cbnm(128-7)='PW_L'
	cbnm(129-7)='PW_H'
	cbnm(130-7)='O3'

      do ir=1,7
	do iy=1,27
	k=1
	  do ic=8,130
	   combin(iy,ir,k)=cldjja(iy,ir,9,ic)  ! k 1-123 cloud
	   k=K+1
	  enddo
!	print*,k,'**'
	   combin(iy,ir,k)=gpc(ir,iy+4)   ! 124 rain
	   cbnm(k)='Rain'
	   k=k+1
!      print*,k,'***'
	  do i=1,5
	   combin(iy,ir,k)=q1(ir,iy+4,i)  ! 125 126 127 128 129 Q1
	   cbnm(k)=trim(q1q2nm(i))
	   k=k+1
	  enddo
	   do i=1,5
	   combin(iy,ir,k)=q2(ir,iy+4,i) ! 130 131 132 133 134 Q2
	   cbnm(k)=trim(q1q2nm(i+5))
	   k=k+1
 	  enddo
        do i=1,30
	   combin(iy,ir,k)=comps(ir,iy+4,i) ! 135+29
	   cbnm(k)=trim(compnm(i))
	   k=k+1
	  enddo
 
 !     print*,k,'**********'
      enddo 
	enddo
!-----------------output--------------------------
      do ir=1,7
      open(888,file='Z:\DATA\LargeScale\79-12\ANA\'//trim(area(ir))//
     +'_ALL.txt')
      do i=1,130
	write(888,88)(combin(j,ir,i),j=1,27)
      enddo
      close(888)
	enddo
88    format(1X,27(1X,E16.8))
!!!!!!!!!!!!!! --------   K=164---------------------
      k=164
	n=27
	t005=0.381
	t001=0.487
	open(700,file='Z:\DATA\LargeScale\79-12\ANA\T001_04-05.txt')
	open(701,file='Z:\DATA\LargeScale\79-12\ANA\T001_05-06.txt')
	open(702,file='Z:\DATA\LargeScale\79-12\ANA\T001_06-07.txt')
	open(703,file='Z:\DATA\LargeScale\79-12\ANA\T001_07-08.txt')
	open(704,file='Z:\DATA\LargeScale\79-12\ANA\T001_08-10.txt')
	open(707,file='Z:\DATA\LargeScale\79-12\ANA\T001.txt')
	open(705,file='Z:\DATA\LargeScale\79-12\ANA\T005.txt')
	open(706,file='Z:\DATA\LargeScale\79-12\ANA\R_big.txt')
	m=1
      do  ir=1,7
	 do  iv=1,164
!	   print*,iv,'AAAAAAAAAAAA'
	     do iy=1,27
	     tmp1(iy)=combin(iy,ir,iv)
!	print*,tmp1(iy),'AAAAAAAAAAAAAA'
	     enddo
!		  av1=sum(tmp1)/27.
           do  ir2=1,7
	       do  iv2=1,164
		      do iy=1,27
	          tmp2(iy)=combin(iy,ir2,iv2)
!	print*,tmp2(iy),'BBBBBBBBBBBBBB'
	          enddo
!	          av2=sum(tmp2)/27.
!                x=0.0;y1=0.0;y2=0.0
!			  do i=1,27
!			  x=x+(tmp1(i)-av1)*(tmp2(i)-av2)
!	          y1=y1+(tmp1(i)-av1)**2
!	          y2=y2+(tmp2(i)-av2)**2
!			  enddo 
!	          if(y1==0.or.y2==0)then
!	           r=0
!	           else
!			  r=x/sqrt(y1*y2) 
!	           endif
			 call relation(tmp1,tmp2,n,r)
!------------ test--------------------------------------
                if(r==1.and.ir2/=ir.and.iv/=iv2)then
                print*,tmp1,ir,iv
	          print*,tmp2,ir2,iv2
	          endif
	         ttest=0.0
	         if(abs(r)>0.5)then
	         write(706,*)area(ir),iv,'**',area(ir2),iv2,'--',r 
			 endif
               if(isNaN(r))goto  777              
			 rr=r**2
	         if((1-rr)<1e-10) goto 777
			 ttest=r/sqrt(1-rr)*sqrt((n-2)*1.0)

777	continue  
      nn=0
	id1=0
	id1=index(cbnm(iv),'_sd')
	id1=id1+index(cbnm(iv),'_td')
	id1=id1+index(cbnm(iv2),'_td')
	id1=id1+index(cbnm(iv2),'_sd')
	id1=id1+index(cbnm(iv),'_Sd')
	id1=id1+index(cbnm(iv),'_Td')
	id1=id1+index(cbnm(iv2),'_Td')
	id1=id1+index(cbnm(iv2),'_Sd')
      if(m>1)then
      do im=1,m
	if(arf(im)==area(ir2).and.art(im)==area(ir).and.
     +	   varf(im)==cbnm(iv2).and.vart(im)==cbnm(iv)) id1=id1+1
	enddo
	endif
       		 
               if(abs(r)>=t001)then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               fnm='Z:\DATA\LargeScale\79-12\ANA\data\'
        		flnm=trim(adjustl(area(ir)))//'_'//trim(adjustl(cbnm(iv)))
     +//'_VS_'//trim(adjustl(area(ir2)))//'_'//trim(adjustl(cbnm(iv2)))
              
                write(707,666)area(ir),trim(cbnm(iv)),iv,
     +			   area(ir2),trim(cbnm(iv2)),iv2,r, 'T001'
	           if(abs(r)<0.5)then
                    write(700,666)area(ir),trim(cbnm(iv)),iv,
     +			   area(ir2),trim(cbnm(iv2)),iv2,r, 'T001'
	          if(m==1.or.id1==0)then
	          open(400,file=trim(fnm)//'50T001_'//trim(flnm)//'.txt')
	                write(400,*)trim(area(ir))//'_'//trim(cbnm(iv))
     +					,' ',trim(area(ir2))//'_'//trim(cbnm(iv2))
	                 do jj=1,27
	                 write(400,*)tmp1(jj),tmp2(jj)
					 enddo
	  arf(m)=area(ir);art(m)=area(ir2)
        varf(m)=cbnm(iv);vart(m)=cbnm(iv2)
	 m=m+1
	close(400)
	                 endif
	            elseif(abs(r)>=0.5.and.abs(r)<0.6)then
                    write(701,666)area(ir),trim(cbnm(iv)),iv,
     +			   area(ir2),trim(cbnm(iv2)),iv2,r, 'T001'
	if(m==1.or.id1==0)then
	          open(400,file=trim(fnm)//'60T001_'//trim(flnm)//'.txt')
	                write(400,*)trim(area(ir))//'_'//trim(cbnm(iv))
     +					,' ',trim(area(ir2))//'_'//trim(cbnm(iv2))
	                 do jj=1,27
	                 write(400,*)tmp1(jj),tmp2(jj)
					 enddo
	  arf(m)=area(ir);art(m)=area(ir2)
        varf(m)=cbnm(iv);vart(m)=cbnm(iv2)
	 m=m+1
	close(400)
	endif
	            elseif(abs(r)>=0.6.and.abs(r)<0.7)then
                    write(702,666)area(ir),trim(cbnm(iv)),iv,
     +			   area(ir2),trim(cbnm(iv2)),iv2,r, 'T001'
	if(m==1.or.id1==0)then
	          open(400,file=trim(fnm)//'70T001_'//trim(flnm)//'.txt')
	                write(400,*)trim(area(ir))//'_'//trim(cbnm(iv))
     +					,' ',trim(area(ir2))//'_'//trim(cbnm(iv2))
	                 do jj=1,27
	                 write(400,*)tmp1(jj),tmp2(jj)
					 enddo
	  arf(m)=area(ir);art(m)=area(ir2)
        varf(m)=cbnm(iv);vart(m)=cbnm(iv2)
	 m=m+1
	close(400)
	endif
	            elseif(abs(r)>=0.7.and.abs(r)<0.8)then
                    write(703,666)area(ir),trim(cbnm(iv)),iv,
     +			   area(ir2),trim(cbnm(iv2)),iv2,r, 'T001'
	if(m==1.or.id1==0)then
	           open(400,file=trim(fnm)//'80T001_'//trim(flnm)//'.txt')
	                write(400,*)trim(area(ir))//'_'//trim(cbnm(iv))
     +					,' ',trim(area(ir2))//'_'//trim(cbnm(iv2))
	                 do jj=1,27
	                 write(400,*)tmp1(jj),tmp2(jj)
					 enddo
	  arf(m)=area(ir);art(m)=area(ir2)
        varf(m)=cbnm(iv);vart(m)=cbnm(iv2)
	  m=m+1
	close(400)
	endif
	             elseif(abs(r)>=0.8)then
                    write(704,666)area(ir),trim(cbnm(iv)),iv,
     +			   area(ir2),trim(cbnm(iv2)),iv2,r, 'T001'
	                if(r<1.0.and.(n==1.or.id1==0))then
	          open(400,file=trim(fnm)//'100T001_'//trim(flnm)//'.txt')
	                write(400,*)trim(area(ir))//'_'//trim(cbnm(iv))
     +					,' ',trim(area(ir2))//'_'//trim(cbnm(iv2))
	                 do jj=1,27
	                 write(400,*)tmp1(jj),tmp2(jj)
					 enddo
	  arf(m)=area(ir);art(m)=area(ir2)
        varf(m)=cbnm(iv);vart(m)=cbnm(iv2)
	  m=m+1
	close(400)
	                 endif
	             endif
	         elseif(abs(r)>t005)then
               write(705,666)area(ir),trim(cbnm(iv)),iv,
     +			   area(ir2),trim(cbnm(iv2)),iv2,r, 'T005'
	if(m==1.or.id1==0)then
	        open(400,file=trim(fnm)//'T005_'//trim(flnm)//'.txt')
	                write(400,*)trim(area(ir))//'_'//trim(cbnm(iv))
     +					,' ',trim(area(ir2))//'_'//trim(cbnm(iv2))
	                 do jj=1,27
	                 write(400,*)tmp1(jj),tmp2(jj)
					 enddo
	  arf(m)=area(ir);art(m)=area(ir2)
        varf(m)=cbnm(iv);vart(m)=cbnm(iv2)
	 m=m+1
	close(400)
	endif
	         endif
          enddo
         enddo
       enddo
      enddo         


666   format(1X,A4,1X,A8,1X,I3,1X,A4,1X,A8,1X,I3,1X,F8.5,1X,A4)
      end



      subroutine relation(a,b,n,r)
!本程序计算两列向量的相关系数
!a,b分别是待计算的向量
!n是向量的长度，要求两列向量等长


!      integer,intent(in)::n
      real::a(27),b(27)
!      real::relation !返回的相关系数
      integer::i,j   !循环控制变量
      real::sfenzi,sfenmu1,sfenmu2,s,r !加法器
      real::amean,bmean !a,b向量的平均值

      !计算平均值
      s=0.
	ka=0
      do i=1,n
      s=s+a(i)
	if(a(i)==0.)ka=ka+1
      end do
      amean=s/n
      s=0.
      do i=1,n
      s=s+b(i)
	if(b(i)==0.)ka=ka+1
      end do
      bmean=s/n
      print*,ka
!计算相关系数
      sfenzi=0.
      sfenmu1=0.
      sfenmu2=0.
      do i=1,n
      sfenzi=sfenzi+(a(i)-amean)*(b(i)-bmean)
      sfenmu1=sfenmu1+(a(i)-amean)**2
      sfenmu2=sfenmu2+(b(i)-bmean)**2
      end do
	er=1e-16
	if(sfenmu1==0.or.sfenmu2==0.or.ka>=14)then
	r=0.
	else
      r=sfenzi/sqrt(sfenmu1*sfenmu2)
	endif
	return
      end subroutine



      