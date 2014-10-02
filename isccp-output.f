      subroutine output(raw,iy)
      character area(7)*4,dir*100,year*4,fnm*100,tme*2
      integer lon(7,2),lat(7,2)
      real raw(12,9,38,144,72) !!!! 1-144 west to east(-180-180?) /1-72 south to north(-90-90)
      real cld_jja(7,9,38)
!!!!-------------------------set ----------------------------------------------------
      area(1)='PRD'  !!! 110 117.5 / 20 27.5
      lon(1,1)=117 ;lon(1,2)=120 
      lat(1,1)=45 ;lat(1,2)=48   
      area(2)='MLYR' !!! 110 117.5/ 27.5 35
      lon(2,1)=117 ;lon(2,2)=120 
      lat(2,1)=48 ;lat(2,2)=51 
      area(3)='NPC'  !!!! 112.5 120 /35 42.5
      lon(3,1)=118 ;lon(3,2)=121 
      lat(3,1)=51 ;lat(3,2)=54
      area(4)='NEC'  !!! 120 127.5/42.5 50
      lon(4,1)=121 ;lon(4,2)=124 
      lat(4,1)=54 ;lat(4,2)=57
      area(5)='BOB'  !!! 85-95/ 15 25
      lon(5,1)=107 ;lon(5,2)=111 
      lat(5,1)=43 ;lat(5,2)=47
      area(6)='ETP'  !!! 90-102.5/27.5 37.5 
      lon(6,1)=109 ;lon(6,2)=114 
      lat(6,1)=48 ;lat(6,2)=52
      area(7)='WTP'  !!!  77.5-90 /27.5 37.5 
      lon(7,1)=104 ;lon(7,2)=109 
      lat(7,1)=48;lat(7,2)=52
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      dir='Z:\DATA\LargeScale\79-12\ISCCP\'
	write(year,'(I4)')iy
      fnm=trim(dir)//year//'_Month.txt'
	open(999,file=trim(fnm))
	write(999,800)'Areas','CA_T','WP_T','CA_L','CA_M','CA_H'
!----------  types of low clouds
     1,'TY_1','TY1W','TY_2','TY2W','TY_3','TY3W'
     1,'TY_4','TY4W','TY_5','TY5W','TY_6','TY6W'
!----------  types of middle clouds
     1,'TY_7','TY7W','TY_8','TY8W','TY_9','TY9W'
     1,'TY10','T10W','TY11','T11W','TY12','T12W'
!----------  types of high clouds
     1,'TY13','T13W','TY14','T14W','TY15','T15W'
!----------- addtional ------------------------------
     1,'RS_C','PW_L','PW_H'
      do 100 ir=1,7
      lns=lon(ir,1);lne=lon(ir,2)
      lts=lat(ir,1);lte=lat(ir,2)
	 do 100 it=1,9 !  time steps
	  do 100 iv=1,38  !!! varaiables
	   temp=0.0
	   ims=6
	   if(iy==1983)ims=7
	    do 101 im=ims,8
	     do 101 ln=lns,lne
	      do 101 lt=lts,lte 
               temp=temp+ raw(im,it,iv,ln,lt)
101   continue
           XYT=(8-ims+1)*1.0*(lne-lns+1)*(lte-lts+1)
	     cld_jja(ir,it,iv)=temp/XYT
100   continue
      do ir=1,7
	write(999,801)area(ir),(cld(ir,9,i),i=1,38)
	fnm=trim(dir)//year//'_'//trim(area(ir))//'_Times.txt'
	open(998,file=trim(fnm))
	write(998,800)'Times','CA_T','WP_T','CA_L','CA_M','CA_H'
!----------  types of low clouds
     1,'TY_1','TY1W','TY_2','TY2W','TY_3','TY3W'
     1,'TY_4','TY4W','TY_5','TY5W','TY_6','TY6W'
!----------  types of middle clouds
     1,'TY_7','TY7W','TY_8','TY8W','TY_9','TY9W'
     1,'TY10','T10W','TY11','T11W','TY12','T12W'
!----------  types of high clouds
     1,'TY13','T13W','TY14','T14W','TY15','T15W'
!----------- addtional ------------------------------
     1,'RS_C','PW_L','PW_H'
	do it=1,8
	write(tms,'(I2.2)')(it-1)*3
      write(998,802)tme,':00',(cld(ir,it,i),i=1,38)
      enddo
	close(998)
	enddo
	close(999)
800   format(1X,A5,38(1X,A4))	
801   format(1X,A4,38(1X,F8.4))
802   format(1X,A2,A3,38(1X,F8.4)) 
      return
	end subroutine

       