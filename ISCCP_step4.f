      program ISCCP_step4
	character dir*100,fnm*100,T00*10
	character rgn01(38081,2)*4,var01(38081,2)*8
	integer id01(38081,2),id05(37908,2)
	real r01(38081),r05(37908)
	character rgn05(37908,2)*4,var05(37908,2)*8
	character rgn(7)*4


	rgn(1)='PRD'
	rgn(2)='MLYR'
	rgn(3)='NPC'
	rgn(4)='NEC'
	rgn(5)='ETP'
	rgn(6)='WTP'
	rgn(7)='BOB'
      NN=38034  !! lines of T001
	NM=37662  !! line of T005
	dir='Z:\DATA\LargeScale\79-12\ANA\'
	fnm=trim(dir)//'T001_select.txt'
	open(10,file=trim(fnm))
	do i=1,NN
	read(10, 666) rgn01(i,1),var01(i,1),id01(i,1)
     +, rgn01(i,2),var01(i,2),id01(i,2), r01(i),T00
	enddo
	close(10)
	fnm=trim(dir)//'T005_select.txt'
	open(20,file=trim(fnm))
	do i=1,NM
	read(20,666)rgn05(i,1),var05(i,1),id05(i,1)
     +, rgn05(i,2),var05(i,2),id05(i,2), r05(i),T00
	enddo
	close(20)
!1----------  output ------------------------
      do ir=1,7
	fnm=trim(dir)//trim(rgn(ir))//'_inside.txt'
	open(20,file=trim(fnm))
	fnm=trim(dir)//trim(rgn(ir))//'_outside.txt'
	open(30,file=trim(fnm))
!!!!    0.001
	do i=1,NN
         if(trim(rgn01(i,1))==trim(rgn(ir)) )then
	      if(trim(rgn01(i,2))==trim(rgn(ir)) )then
	        write(20,665)rgn01(i,1),var01(i,1),id01(i,1)
     +,              rgn01(i,2),var01(i,2),id01(i,2), r01(i)
	      endif
             if(trim(rgn01(i,2))/=trim(rgn(ir)) )then
	        write(30,665)rgn01(i,1),var01(i,1),id01(i,1)
     +,              rgn01(i,2),var01(i,2),id01(i,2), r01(i)
	      endif
         endif
	    if(trim(rgn01(i,2))==trim(rgn(ir)) )then
	       if(trim(rgn01(i,1))/=trim(rgn(ir)) )then
	        write(30,665)rgn01(i,2),var01(i,2),id01(i,2)
     +,              rgn01(i,1),var01(i,1),id01(i,1), r01(i)
	      endif
	    endif
	end do
!!!!    0.005
	do i=1,NM
         if(trim(rgn05(i,1))==trim(rgn(ir)) )then
	      if(trim(rgn05(i,2))==trim(rgn(ir)) )then
	        write(20,665)rgn05(i,1),var05(i,1),id05(i,1)
     +,              rgn05(i,2),var05(i,2),id05(i,2), r05(i)
	      endif
             if(trim(rgn05(i,2))/=trim(rgn(ir)) )then
	        write(30,665)rgn05(i,1),var05(i,1),id05(i,1)
     +,              rgn05(i,2),var05(i,2),id05(i,2), r05(i)
	      endif
         endif
		 if(trim(rgn05(i,2))==trim(rgn(ir)) )then
	       if(trim(rgn05(i,1))/=trim(rgn(ir)) )then
	        write(30,665)rgn05(i,2),var05(i,2),id05(i,2)
     +,              rgn05(i,1),var05(i,1),id05(i,1), r05(i)
	      endif
	    endif
	end do
	close(20)
	close(30)

	enddo

665   format(1X,A4,1X,A8,1X,I3,1X,A4,1X,A8,1X,I3,1X,F8.5)
666   format(1X,A4,1X,A8,1X,I3,1X,A4,1X,A8,1X,I3,1X,F8.5,1X,A4)

      end
	
