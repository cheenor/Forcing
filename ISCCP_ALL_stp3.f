      program ISCCP_STEP3
	character dir*100,fnm*100,T00*4
	character str1*8,str2*8,ar1*4,ar2*4
	character arf(1000000)*4,art(1000000)*4
	character varf(1000000)*8,vart(1000000)*8
	
	
	dir='Z:\DATA\LargeScale\79-12\ANA\'
	fnm=trim(dir)//'T001.txt'
	open(10,file=trim(fnm))
	open(11,file=trim(dir)//'T001_select.txt')
	open(12,file=trim(dir)//'40_T001_select.txt')
	open(13,file=trim(dir)//'50_T001_select.txt')
	open(14,file=trim(dir)//'60_T001_select.txt')
	open(15,file=trim(dir)//'70_T001_select.txt')
	open(16,file=trim(dir)//'80_T001_select.txt')
	k=1
	do while(istat==0)
	read(10,fmt=666,iostat=istat)ar1,str1,iv,
     +			   ar2,str2,iv2,r, T00
	id1=0
	id1=index(str1,'_sd')
	id1=id1+index(str1,'_td')
	id1=id1+index(str2,'_td')
	id1=id1+index(str2,'_sd')
	id1=id1+index(str1,'_Sd')
	id1=id1+index(str1,'_Td')
	id1=id1+index(str2,'_Td')
	id1=id1+index(str2,'_Sd')
	id1=id1+index(str1,'_sd')
	id1=id1+index(str1,'_td')
	id1=id1+index(str2,'_td')
	id1=id1+index(str2,'_sd')
	if(id1==0)then
	if(k==1)then
	write(11,666)ar1,str1,iv,ar2,str2,iv2,r,T00	 
	arf(k)=ar1;art(k)=ar2
      varf(k)=str1;vart(k)=str2
	 if(r<0.5)then
        write(12,666)ar1,str1,iv,ar2,str2,iv2,r,T00	
	 elseif(r>=0.5.and.r<0.6)then
	  write(13,666)ar1,str1,iv,ar2,str2,iv2,r,T00	
	 elseif(r>=0.6.and.r<0.7)then
	  write(14,666)ar1,str1,iv,ar2,str2,iv2,r,T00	
	 elseif(r>=0.7.and.r<0.8)then
	  write(15,666)ar1,str1,iv,ar2,str2,iv2,r,T00	
	 elseif(r>=0.8.and.r<1.0)then
        write(16,666)ar1,str1,iv,ar2,str2,iv2,r,T00	
	 endif
	endif
	if(k>1)then
	jj=0
        do j=1,k
        if(arf(j)==ar2.and.art(j)==ar1.and.
     +	   varf(j)==str2.and.vart(j)==str1) jj=jj+1
	  enddo
        if(jj==0)then
	  write(11,666)ar1,str1,iv,ar2,str2,iv2,r, T00
	  arf(k)=ar1;art(k)=ar2
        varf(k)=str1;vart(k)=str2
	  if(r<0.5)then
         write(12,666)ar1,str1,iv,ar2,str2,iv2,r,T00	
	  elseif(r>=0.5.and.r<0.6)then
	   write(13,666)ar1,str1,iv,ar2,str2,iv2,r,T00	
	  elseif(r>=0.6.and.r<0.7)then
	   write(14,666)ar1,str1,iv,ar2,str2,iv2,r,T00	
	  elseif(r>=0.7.and.r<0.8)then
	   write(15,666)ar1,str1,iv,ar2,str2,iv2,r,T00	
	  elseif(r>=0.8.and.r<1.0)then
         write(16,666)ar1,str1,iv,ar2,str2,iv2,r,T00	
	  endif	 
	 endif
	endif

	k=k+1
	endif
	enddo
	close(10)


	fnm=trim(dir)//'T005.txt'
	open(20,file=trim(fnm))
	open(21,file=trim(dir)//'T005_select.txt')
	k=1
	do while(istat2==0)
	read(20,fmt=666,iostat=istat2)ar1,str1,iv,
     +			   ar2,str2,iv2,r, T00
	id1=0
	id1=index(str1,'_sd')
	id1=id1+index(str1,'_td')
	id1=id1+index(str2,'_td')
	id1=id1+index(str2,'_sd')
	id1=id1+index(str1,'_Sd')
	id1=id1+index(str1,'_Td')
	id1=id1+index(str2,'_Td')
	id1=id1+index(str2,'_Sd')
	id1=id1+index(str1,'_sd')
	id1=id1+index(str1,'_td')
	id1=id1+index(str2,'_td')
	id1=id1+index(str2,'_sd')
	if(id1==0)then
	if(k==1)then
	write(21,666)ar1,str1,iv,ar2,str2,iv2,r,T00	 
	arf(k)=ar1;art(k)=ar2
      varf(k)=str1;vart(k)=str2
	endif
	if(k>1)then
	jj=0
        do j=1,k
        if(arf(j)==ar2.and.art(j)==ar1.and.
     +	   varf(j)==str2.and.vart(j)==str1) jj=jj+1
	  enddo
        if(jj==0)then
	  write(21,666)ar1,str1,iv,ar2,str2,iv2,r, T00
	  arf(k)=ar1;art(k)=ar2
        varf(k)=str1;vart(k)=str2
	  endif
	endif
	k=k+1
	endif
! add the following 2014/01/05   interregion and multi-regions

 
	enddo










666   format(1X,A4,1X,A8,1X,I3,1X,A4,1X,A8,1X,I3,1X,F8.5,1X,A4)



      end
