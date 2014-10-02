      program ISCCP2
      character area(7)*4,TME*5,dir*100,fnm*100,year*4
	real cld_jja(27,7,9,39)

      dir='Z:\DATA\LargeScale\79-12\ISCCP\'
      do iyy=1983,2009
	iy=iyy-1982
	 write(year,'(I4)')iyy
       fnm=trim(dir)//year//'_Month.txt'
	 open(999,file=trim(fnm))
	 read(999,*)
	 do ir=1,7
	  read(999,801)area(ir),(cld_jja(iy,ir,9,i),i=1,39)
	 enddo  
       close(999)
	 do ir=1,7
	  fnm=trim(dir)//year//'_'//trim(area(ir))//'_Times.txt'
	  open(998,file=trim(fnm))
	  read(998,*)
        do it=1,8
	   read(998,802)TME,(cld_jja(iy,ir,it,i),i=1,39)
	  enddo
	  close(998)
	  enddo
      enddo
!!!----------------------------------------------------------------
      do ir=1,7
      fnm=trim(dir)//trim(area(ir))//'_Years_JJA.txt'
	open(10,file=trim(fnm))
      write(10,'(A100)')'The type 1 and 7 are Cumulus and Altocumulus
     1 with A TAU	smaller than 3.55'
	write(10,799)'Year','Cloud_Cover','WPH','CA_L','CA_M','CA_H'
     1,'TY_1','TY1W','TY_7','TY7W','TY_9','TY9W','TY12','T12W'
	1,'TY15','T15W','PW_L','PW_H','TS_K'
	 do iyy=1983,2009
	  iy=iyy-1982
	  write(10,800)iyy,(cld_jja(iy,ir,9,i),i=1,5),cld_jja(iy,ir,9,6)
     1,cld_jja(iy,ir,9,7),cld_jja(iy,ir,9,18),cld_jja(iy,ir,9,19)
	1,cld_jja(iy,ir,9,22),cld_jja(iy,ir,9,23),cld_jja(iy,ir,9,28)
	1,cld_jja(iy,ir,9,29),cld_jja(iy,ir,9,34),cld_jja(iy,ir,9,35)
	1,cld_jja(iy,ir,9,37),cld_jja(iy,ir,9,38)
     1,cld_jja(iy,ir,9,39)-273.15
	enddo
	enddo
799   format(1X,A4,1X,A11,1X,A3,16(1X,A4))
800   format(1X,I4,18(1X,F12.4))
801   format(1X,A4,39(1X,F12.4))
802   format(1X,A2,A3,39(1X,F12.4)) 
      end