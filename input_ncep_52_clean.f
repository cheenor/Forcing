      program sounding
!-------------------------------------------------------------------
!      
!
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
cc
cc data provided by Xin Lin (CSU, Dick Johnson's student) 
cc
cc we add 3 levels in the stratosphere (21, 31, 42 km) assuming
cc constant T to get deep sounding.
cc
cc number of time levels:
      parameter(ntdat=125)
      dimension tidat(ntdat),pl1(ntdat),pl2(ntdat),pl3(ntdat)
cc level of test plot
      data ktest /30/
cc below are arrays defined at clarks model levels:
      parameter(l=52,lp=l+1)
      dimension xi(lp),xis(l),z(l)
      dimension tx(2),tz(2)
      data rd,cp,g,rv,hlat /287.,1005.,9.81,461.,2.5e6/
      dimension tme1(l),the1(l),qve1(l),ue1(l),ve1(l),dtls(l),dqls(l)
      dimension out1(l),out2(l),we1(l),q1ls(l),q2ls(l)
	dimension outdy(l,4),omg_adj(l)
cc npin is number of levels
!      parameter(npin0=39,npin=42)
	parameter(npin0=17,npin=19)
      dimension press(npin),temp(npin),thet(npin),zin(npin),
     1          vap(npin),uu(npin),vv(npin),ww(npin),rh(npin)
      dimension tlsf(npin),qlsf(npin),XYN_OUT(npin,19)
	real gz(npin),the(npin)
cc
cc SST data from NMC analysis (tsst is time in days, dsst is in deg C)
      parameter(nsst=9)
	character*100 dir,filepath,path,fold,filename,filepath2
	character*4 yearstr,area(4)
	character*16 date
	integer ims(4),ids(4),ime(4),ide(4),days(12),tmid
      dimension tsst(nsst),dsst(nsst)
!--------------------------------------------------------------------
      parameter(ntm=124)
      real lhf(ntm),shf(ntm)
	data rd,cp,g,rv,hlat /287.,1005.,9.81,461.,2.5e6/
!--------------------------------------------------------------------
      character*16 celord
      real XYD(npin,19),dyn(npin,4)
!---------------set the time ----------------------------------------
      iyr=2010
	ims(1)=4  ;ime(1)=5
	ims(2)=6  ;ime(2)=7
	ims(3)=7  ;ime(3)=8
	ims(4)=8  ;ime(4)=8
	ids(1)=3  ;ide(1)=3
	ids(2)=24  ;ide(2)=23
	ids(3)=26  ;ide(3)=24
	ids(4)=1  ;ide(4)=31
	write(yearstr,'(I4)')iyr
	fold=yearstr(3:4)//'0101-'//yearstr(3:4)//'1231\'
      dir='Z:\DATA\LargeScale\NcepR2_Pre\'
	area(1)='PRD'
	area(2)='MLYR'
	area(3)='NPC'
	area(4)='NEC'
	path=trim(dir)//trim(fold)
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
	 if(mod(iyr,4)==0.and.mod(iyr,100)/=0)then
	 days(2)=29 
	 elseif(mod(iyr,400)==0)then
	 days(2)=29
	 endif
	filepath=trim(dir)//trim(fold)//'daystrt.txt'
	open(997,file=trim(filepath))
!-----------------------------------
cc
c       print*,'  ktest ??'
c       read *,ktest
cc
      tx(1)=0.
      tx(2)=40.
      tz(1)=0.
      tz(2)=0.
ccccccccccccccc code below taken out from clarks model setup:
      nzm=l-1 
cc
cc
      rat=15.
      XI(2)=0.
      DO 152 K=2,NZM
      RATZ=RAT
      DEL=100.
      nzm1=L-1
      k1=k
      XI(K+1)=XI(K)+((RATZ-1.)/FLOAT(NZM1-2)*FLOAT(K1-2)+1.)*DEL
      PRINT 501,K+1,XI(K+1)
  501 FORMAT(2X,'** GRID: K,XI:  ',I5,E16.8)
  152 CONTINUE
      xi(1)=2.*xi(2)-xi(3)
      xi(lp)=2.*xi(l)-xi(l-1)
      do k=1,l
      xis(k)=.5*(xi(k)+xi(k+1))
      z(k)=xis(k)
C	print*,z(k)
      enddo
      open(99,file='Z-Geo.txt')
	write(99,9999)(z(k),k=1,l)
9999   format(1X,52(1X,F10.3))	 
ccc sst data: convert time into hours
      do iii=1,nsst
      tsst(iii)=tsst(iii)*24.
      enddo
      nstsst=2 ! STARTING INDEX FOR INTERPOLATION FROM SST DATASET

        iwrite=0
!!!-----TOGA sounding file open--------------
!      open(10,file=
!     *'/mnt/raid50/hiba/data_toga_obs_dat/ifa_dat_new.sounding'
!     *,status='old')
!      open(20,file=
!     *'/mnt/raid50/hiba/data_toga_obs_dat/ifa_dat_new.forcing'
!     *,status='old')
!      open(30,file=
!     *'/home/wuxq/forcing_toga/toga30/2dforcing/misc.ifa'
!     *,status='old')
!      open(40,file=
!     *'/mnt/raid50/hiba/data_toga_obs_dat/ifa_dat.surface'
!     *,status='old')
!----------------------------------------------------------
      do 1014 ip=3,3   ! area loops
	tempress=0
	temptemp=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       imts=0
	 imte=0 
         do imt=1,ims(ip)-1
	     imts=imts+days(imt)
	   enddo
	   imts=imts*4+ids(ip)*4-4
	   do imt=1,ime(ip)-1
	     imte=imte+days(imt)
	   enddo	  
         imte=imte*4+ide(ip)*4
	   imt=imte-imts 
	   print*,imt   
	   write(997,*)imts ,ip
!-------------------------------------------------------------
      filepath=trim(dir)//trim(fold)//
     +yearstr//trim(area(ip))//'_RAW.txt'
      filepath2=trim(dir)//trim(fold)//
     +yearstr//trim(area(ip))//'_Surface.txt'
	     
	   call OBS(imts,ip,filepath,filepath2,imt)
      open(90,file=trim(filepath))
	read(90,*)
	do i=1,imts*17
      read(90,*)
	enddo


!----target 1 2 3 .....imts imts+1 ......imte.......1460-----------------

!!!-------NCEP imitative sounding files open---------------
      filepath=trim(dir)//trim(fold)//trim(area(ip))//'_sounding.txt'
!   Press(hPa) heigh(m) U(m/s) V(m/s)  omega(pa/s) Temp(K) Theta(K) Qv(kg/kg) RH(%)
      open(10,file=trim(filepath))
	read(10,*)
	read(10,*)
      filepath=trim(dir)//trim(fold)//trim(area(ip))//'_Forcing.txt'
!	TimeID 17_levels_T_forcing(K/day) 17_levels_qv_forcing(K/day)
      open(20,file=trim(filepath))
	read(20,*)
	read(20,*)
!      filepath=trim(dir)//trim(fold)//'2010'//trim(area(ip))//'_RAW.txt'
!	'DATE','HOUR','T_ls(k/day)','Q_ls(k/day)','U(m/s)',
!     +'V(m/s)','moisture(kg/kg)','HGT(m)','AIR(K)','Adj_omega(pa/s)',
!     +'RH(%)','Theta(K)','Q1(k/day)','Q2(k/day)','HADQ(K/day)',
!     +'VADQ(K/day)','TCHQ(K/day)','HADT(K/day)','HADT(K/day)',
!     +'TCHT(K/day)','Ori_omega(pa/s)'
!      open(20,file=trim(filepath))
!	read(20,*)

	filepath=trim(dir)//trim(fold)//trim(area(ip))//'_surface22.txt'
!	Press(hPa) heigh(m) U(m/s) V(m/s)  omega(pa/s) Temp(K) Theta(K) Qv(kg/kg) RH(%) for Surface
      open(40,file=trim(filepath))
	read(40,*)
	read(40,*)
	filepath=trim(dir)//trim(fold)//trim(area(ip))//'_HF.txt'
!	TimeID latentHeat(W/m^2) SensibleHeat(W/m^2)  Precipipitable_water_for_entire_atmosphere(kg/m^2)
      open(30,file=trim(filepath))
	read(30,*)
	read(30,*)
      read(30,*)
!!!
	filepath=trim(dir)//trim(fold)//yearstr//trim(area(ip))//'_dyn.txt'
!	TimeID latentHeat(W/m^2) SensibleHeat(W/m^2)  Precipipitable_water_for_entire_atmosphere(kg/m^2)
      open(50,file=trim(filepath))
	read(50,*)
cc
c      do ik=1,120
c      read(30,876) iy,im,id,ih,bt,sst,fsh,flh
c      enddo]
!------------skip ----------------------------------------
       do i=1,imts
       read(30,*)
	 enddo
       do i=1,imts*18
	 read(10,*)
	 enddo
       do i=1,imts
	 read(20,*)
	 enddo
	 do i=1,imts
	 read(40,*)
	 enddo
	 do i=1,imts*17
	 read(50,*)
	 enddo
!---------------------------------------------------
      do 999 itim=1,imt !ntdat  
      read(30,301) tmid,flh,fsh,pewr,cpewr
!	print*, tmid,flh,fsh,pewr,cpewr
301   format(1X,I4,1X,F8.2,1X,F8.2,1X,e12.4,1X,e12.4)
c 876  format(4i5,4f8.2)
cc read sounding data for this time level:
      read(10,*)tmid
      do k=1,npin0
	
      read(10,101) press(k),gz(k),uu(k),vv(k),ww(k),temp(k),
     +	the(k),vap(k),rh(k)
      enddo
!	print*,press(k),gz(k),uu(k),vv(k),ww(k),temp(k),
!     +	the(k),vap(k),rh(k)
c 176  format(1x,7e18.8)
cc read first level (surface)
      read(40,101) press(1),gz(1),uu(1),vv(1),ww(1),temp(1),
     +	the(1),vap(1),rh(1),sst
!	print*, press(1),gz(1),uu(1),vv(1),ww(1),temp(1),
!     +	the(1),vap(1),rh(1)
       tempress=press(1)+tempress
       temptemp=temp(1)+temptemp
101   format(1X,4(1X,F9.3),1X,e12.4,2(1X,F9.2),1X,e12.4,1X,F7.3,1X,F7.3)
!      press(1)=1008.
      ww(1)=0.
	do ik=1,17
	read(50,517)date,(dyn(ik,ir),ir=1,4)
      enddo  
517   format(1X,A16,1X,4(1X,e12.4))

cc extrapolate first level (surface)
c     press(1)=1008.
c     coe2=(press(1)-press(2))/(press(3)-press(2))
c     temp(1)=coe2*temp(3) + (1.-coe2)*temp(2)
c     vap(1)=coe2*vap(3) + (1.-coe2)*vap(2)
c     uu(1)=coe2*uu(3) + (1.-coe2)*uu(2)
c     vv(1)=coe2*vv(3) + (1.-coe2)*vv(2)

cc read ls forcing data for this time level:
      tlsf(1)=0.
      qlsf(1)=0.
       do k=npin0+1,npin
      tlsf(k)=0.
      qlsf(k)=0.
       enddo
      read(20,201)tmid,tlsss,(tlsf(K),K=1,npin0),
     + qlsfss, (qlsf(K),K=1,npin0)
201   format(1X,I4,36(1X,e12.4))

c      read(20,776) (qlsf(K),K=2,npin0)

      tlsf(1)=0.
      qlsf(1)=0.
c 776  format(1x,5e20.8)
!--------------read the q1 and q2 whcih is written by Chenjinghua 
        do ik=1,17
	   read(90,906)celord,(XYD(ik,kk),kk=1,19)
c	print*,celord
c	stop
	  enddo
906   format(1X,A16,1X,19(1X,e12.4))
convert from temperature (deg C or K) into potential temperature
      do k=1,npin0
      temp(k)=temp(k) ! +273.16
      thet(k)=temp(k)*(1.e3/press(k))**(rd/cp)
cc !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      tlsf(k)=tlsf(k)*thet(k)/temp(k) ! theta forcing in K/day
c      tlsf(k)=tlsf(k) ! temperature forcing in K/day
      qlsf(k)=qlsf(k)*cp/hlat           ! qv forcing in kg/kg/day
cc !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      enddo

compute approximated height of pressure levels:
      zin(1)=0.
      do k=2,npin0
          km=k-1
            tempk =temp(k ) * (1.+.6e-3*vap(k ))   !!!!Vap  vapor mixing  kg/kg
            tempkm=temp(km) * (1.+.6e-3*vap(km))
                delt=tempk-tempkm
                   if(delt.gt.1.e-4) then
                      tavi=alog(tempk/tempkm)/delt
                   else
                      tavi=1./tempk
                   endif
               deltz=-rd/(tavi*g) * alog(press(k)/press(km))
            zin(k)=zin(km)+deltz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	      zin(k)=gz(k)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      enddo
cc extrapolate stratosphere:
!      zin(npin0+1)=21.e3
!      zin(npin0+2)=31.e3
!      zin(npin0+3)=42.5e3
!      print*,zin(16),zin(17),zin(18)
!	zin(16)=21.e3
      zin(18)=35.e3
      zin(19)=42.5e3
!	stop
      temp00=temp(npin0)
      z00=zin(npin0)
      p00=press(npin0)
      den00=p00*1.e2/(rd*temp00)
      coe=g/(rd*temp00) 
cc get rh at k=npin0:
      rh00=rh(npin0)/100.
cc          print*,' rh at k=npin0: ',rh00

      do k=npin0-3,npin0
      den=press(k)*1.e2/(rd*temp(k))
      esat=611.*exp(hlat/rv * (1./273.16 - 1./temp(k)))
      qvs00=esat/(rv*den*temp(k))
      rh(k)=0.1*rh(k-1)
      vap(k)=rh(k)/100.*qvs00*1.e3
	vap(k)=vap(k)*1e-3  !!!!! unit  must follow the input data
      enddo

	do k=npin0+1,npin
      press(k)=p00*exp(-coe*(zin(k)-z00))
      temp(k)=temp(npin0)
      thet(k)=temp(k)*(1.e3/press(k))**(rd/cp)
      den=press(k)*1.e2/(rd*temp(k))
      esat=611.*exp(hlat/rv * (1./273.16 - 1./temp(k)))
      qvs00=esat/(rv*den*temp(k))
c     vap(k)=0.1*rh00*qvs00*1.e3
c     rh(k)=0.1*rh00*100.
      rh(k)=0.1*rh(k-1)
      vap(k)=rh(k)/100.*qvs00*1.e3
      vap(k)=vap(k)*1e-3  !!!!! unit  must follow the input data
      uu(k)=uu(npin0)
      vv(k)=vv(npin0)
      ww(k)=ww(npin0)
      enddo
ccc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	from 1000 or surface
c write initial sounding to fort.33 if itim is the starting time:
      filename=trim(dir)//trim(fold)//trim(area(ip))//'.33'
	open(33,file=trim(filename))
!      if(itim.eq.17) then
       if(itim.eq.1) then
	npin00=npin
      write(33,701) (press(k),k=1,npin00)
701    format(5x,16h  data press  / /
     1        5x,3h1  ,7(f7.2,1h,)/
     1        5x,3h1  ,7(f7.2,1h,)/
!     1        5x,3h1  ,7(f7.2,1h,)/
!     1        5x,3h1  ,7(f7.2,1h,)/
!     1        5x,3h1  ,7(f7.2,1h,)/
     1        5x,3h1  ,4(f7.2,1h,),f7.2,1h/) 
      write(33,702) (temp(k)-273.16,k=1,npin00)
702    format(5x,15h  data temp  / /
     1        5x,3h1  ,7(f7.2,1h,)/
     1        5x,3h1  ,7(f7.2,1h,)/
!     1        5x,3h1  ,7(f7.2,1h,)/
!     1        5x,3h1  ,7(f7.2,1h,)/
!     1        5x,3h1  ,7(f7.2,1h,)/
     1        5x,3h1  ,4(f7.2,1h,),f7.2,1h/) 
      write(33,703) (vap(k)*1000,k=1,npin00)
703    format(5x,14h  data vap  / /
     1        5x,3h1  ,5(e10.3,1h,)/
     1        5x,3h1  ,5(e10.3,1h,)/
     1        5x,3h1  ,5(e10.3,1h,)/
!     1        5x,3h1  ,5(e10.3,1h,)/
!     1        5x,3h1  ,5(e10.3,1h,)/
!     1        5x,3h1  ,5(e10.3,1h,)/
!     1        5x,3h1  ,5(e10.3,1h,)/
!     1        5x,3h1  ,5(e10.3,1h,)/
     1        5x,3h1  ,3(e10.3,1h,),e10.3,1h/) 
      write(33,704) (uu(k),k=1,npin00)
704    format(5x,12h  data u  / /
     1        5x,3h1  ,7(f7.2,1h,)/
     1        5x,3h1  ,7(f7.2,1h,)/
!     1        5x,3h1  ,7(f7.2,1h,)/
!     1        5x,3h1  ,7(f7.2,1h,)/
!     1        5x,3h1  ,7(f7.2,1h,)/
     1        5x,3h1  ,4(f7.2,1h,),f7.2,1h/) 
      write(33,705) (vv(k),k=1,npin00)
705    format(5x,12h  data v  / /
     1        5x,3h1  ,7(f7.2,1h,)/
     1        5x,3h1  ,7(f7.2,1h,)/
 !    1        5x,3h1  ,7(f7.2,1h,)/
 !    1        5x,3h1  ,7(f7.2,1h,)/
 !    1        5x,3h1  ,7(f7.2,1h,)/
     1        5x,3h1  ,4(f7.2,1h,),f7.2,1h/) 

ccc flag to write other data
          itims=itim
          iwrite=1

      endif
     
compute environmental profiles from sounding assuming no topography:
      iisn=1
      the1(1)=thet(iisn)
      tme1(1)=temp(iisn)
      qve1(1)=vap(iisn)      !*1.e-3
      ue1(1)=uu(1)
      ve1(1)=vv(1)
      we1(1)=ww(1)
      presst=press(iisn)
      q1ls(1)=XYD(1,11)
	q2ls(1)=XYD(1,12)
cc integrate upwards:
      filename=trim(dir)//trim(fold)//trim(area(ip))//'_uv_profiles.35'
	open(35,file=trim(filename))
	
      filename=trim(dir)//trim(fold)//trim(area(ip))//'_lsforcing.37'
	open(37,file=trim(filename))
	
      filename=trim(dir)//trim(fold)//trim(area(ip))//'_surface.39' 
	open(39,file=trim(filename))
	  filename=trim(dir)//trim(fold)//trim(area(ip))//'.49'
	open(49,file=trim(filename))
	  filename=trim(dir)//trim(fold)//trim(area(ip))//
     +'_thetaqv_profile.41'   !!!!!! unit theta(K)  qv g/kg 
	open(41,file=trim(filename))
	  filename=trim(dir)//trim(fold)//trim(area(ip))//'.43'
	open(43,file=trim(filename))
	  filename=trim(dir)//trim(fold)//trim(area(ip))//'.99'
	open(99,file=trim(filename))
!
	 filename=trim(dir)//trim(fold)//trim(area(ip))//'.dyn'
	open(998,file=trim(filename))
       filename=trim(dir)//trim(fold)//trim(area(ip))//'.dyn.417'
	open(417,file=trim(filename))
	
      do 64 k=2,l
       do kk=2,npin
        iisn=kk-1
!	   print*,zin(kk),z(k),k,kk
        if(zin(kk).ge.z(k)) go to 665
       enddo
       print*,' *** input sounding does not go high enough. stop.'
              stop 'sounding'
 665   continue 
       coe2=(z(k)-zin(iisn))/(zin(iisn+1)-zin(iisn))
       the1(k)=coe2*thet(iisn+1) + (1.-coe2)*thet(iisn)
       tme1(k)=coe2*temp(iisn+1) + (1.-coe2)*temp(iisn)
       qve1(k)=(coe2*vap(iisn+1) + (1.-coe2)*vap(iisn))!*1.e-3
       ue1(k)=coe2*uu(iisn+1) + (1.-coe2)*uu(iisn)
       ve1(k)=coe2*vv(iisn+1) + (1.-coe2)*vv(iisn)
       we1(k)=coe2*ww(iisn+1) + (1.-coe2)*ww(iisn)
       dtls(k)=coe2*tlsf(iisn+1) + (1.-coe2)*tlsf(iisn)
       dqls(k)=coe2*qlsf(iisn+1) + (1.-coe2)*qlsf(iisn)
	 q1ls(k)=coe2*XYD(iisn+1,11) + (1.-coe2)*XYD(iisn,11)
	 q2ls(k)=coe2*XYD(iisn+1,12) + (1.-coe2)*XYD(iisn,12)
	 omg_adj(k)=coe2*XYD(iisn+1,8) + (1.-coe2)*XYD(iisn,8) 
	  do idy=1,4
       outdy(k,idy)=coe2*dyn(iisn+1,idy) + (1.-coe2)*dyn(iisn,idy)
 	  enddo
        
 64   continue

cc scale and write to files:
                 if(iwrite.eq.1) then
	    itims=1
          time = float(itim-itims)*6.
          tidat(itim-itims+1)=(itim-itims)*6./24.
cc   velocity profiles for selected period (fort.35); NOTE ROTATION
          svel=10.    ! velocity scale in Clarks model
          do k=1,l
          out1(k)=-ve1(k)/svel
          out2(k)= ue1(k)/svel
          enddo
          write(35,801) time,out1,out2
801       format(10f8.3)
cc   profiles for l-s forcing terms (fort.37)
          day=24.*3600.
          out1(1)=0.
          out2(1)=0.
          do k=2,l
          out1(k)=dtls(k)/day    ! now in K/sec
          out2(k)=dqls(k)/day    ! now in kg/kg/sec
cccccccccccc
cc   set forcing to zero above 17km (17.19km at level 33)
c         if(k.ge.33) then
c         out1(k)=0.
c         out2(k)=0.
c         endif
cccccccccccccccccccc
          enddo
          write(37,802) time,out1,out2
          write(99,802) time,q1ls,q2ls
         write(998,802)time, outdy(:,1),outdy(:,2),outdy(:,3)
     +	   ,outdy(:,4),  we1
	    write(417,802)time, ww
		write(417,802)time, we1
802       format(8e12.4)

            pl1(itim-itims+1)=out1(ktest)*day
            pl2(itim-itims+1)=out2(ktest)*1.e3*day

cc   time series of ocean surface theta and qv (fort.39)
c     time2=time+itims*6.
c     do iiii=1,100
c     if(time2.gt.tsst(nstsst))   nstsst=nstsst+1 
c     enddo
ccccccc interpolate ocean temp:
c      coe2=(time2-tsst(nstsst-1))/(tsst(nstsst)-tsst(nstsst-1))
c        sst=coe2*dsst(nstsst) + (1.-coe2)*dsst(nstsst-1)
          sst=temp(1)
          ths=sst*(the1(1)+the1(2))/(tme1(1)+tme1(2))
          den=press(1)*1.e2/(rd*temp(1))
          esat=611.*exp(hlat/rv * (1./273.16 - 1./sst))
          qvss=esat/(rv*den*sst)
          write(39,803) time,ths,qvss,sst,flh,fsh   !!! what are the units?
          write(49,903) time,sst,ths,qvss
803       format(6e16.5)
903       format(4e16.5)
            pl3(itim-itims+1)=sst
cc   theta and qv profiles for selected period (fort.41)
          do k=1,l
          out1(k)=the1(k)
          out2(k)=qve1(k)
          enddo
          write(41,804) time,out1,out2
          write(43,804) time,out1,out2,tme1
804       format(8e13.5)
                     endif

999      continue
         
         write(997,*)tempress/125. ,ip
	   write(997,*)temptemp/125., ip
1014     continue  
      stop
      end


      subroutine OBS(itt,ip,fraw,fsurf,ntdat)

	character*100 fraw,fsurf,fouts
	character*16 celord(ntdat) , celors(ntdat)
      real XYD(ntdat,17,19)
      real XYS(ntdat,6)
	real Q1Q2(ntdat,2)


!	'DATE','HOUR','1 T_ls(k/day)','2 Q_ls(k/day)','3 U(m/s)',
!     +'4 V(m/s)','5 moisture(kg/kg)','6 HGT(m)','7 AIR(K)','8 Adj_omega(pa/s)',
!     +'9 RH(%)','10 Theta(K)','11 Q1(k/day)','12 Q2(k/day)','13 HADQ(K/day)',
!     +'14 VADQ(K/day)','15 TCHQ(K/day)','16 HADT(K/day)','17 HADT(K/day)',
!     +'18 TCHT(K/day)','19 Ori_omega(pa/s)'
      open(20,file=trim(fraw))
	read(20,*)
	do i=1,itt*17
      read(20,*)
	enddo
	open(30,file=trim(fsurf))
	read(30,*)
	do i=1,itt
      read(30,*)
	enddo
	do it=1,ntdat
	Q1t=0.
	Q2t=0.
        do ik=1,17
	   read(20,906)celord(it),(XYD(it,ik,kk),kk=1,19)
C	print*,celord(it)
	    Q1t=Q1t+XYD(it,ik,11)
	    Q2t=Q2t+XYD(it,ik,12)
	  enddo
	    Q1Q2(it,1)=Q1t
          Q1Q2(it,2)=Q2t
        read(30,907)celors(it),(XYS(it,kk),kk=1,6) !LH(W/m^2)   SH(W/m^2)       prate      cprate  DLR(W/m^2)  ULR(W/m^2)
	enddo
!---------------output-----------------------------------------
       ilen=len_trim(fraw)
	 fouts=fraw(1:ilen-7)//'OBS_Surface_input.txt'
       open(40,file=trim(fouts))
	 write(40,*)'Timestep  [Q1](k/d)   [Q2](k/d)  Rainfall(mm/hr)'
       fouts=fraw(1:ilen-7)//'Q1Q2_profiles_input.txt'
       open(50,file=trim(fouts))
	 write(50,*)'Timestep Q1_level1_to_17 Q2_level1_to_17'  !!!! the first level is not zero 

	 do it=1,ntdat
	 write(40,1014)it*1.0,Q1Q2(it,1),Q1Q2(it,2),XYS(it,3)*3600.0
       write(50,417)it*1.0,(XYD(it,kk,11),kk=1,17),
     +	 (XYD(it,kk,12),kk=1,17)
	 enddo 


906   format(1X,A16,1X,19(1X,e12.4))
907   format(1X,A16,1X,6(1X,e12.4))
1014  format(4e12.4)
417   format(8e12.4)
      close(20)
	close(30)
      end subroutine