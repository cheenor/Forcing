program Inteplot_V
!---------------------------------------------------------
! this program convert pressure to km in vertical direct
!
!----------------------------------------------------------
implicit none
integer ko  !!! the output levels
integer ki  !!! the input pressure levels
parameter(ko=125)
parameter(ki=17)

real zin(ki),tempk,tempkm,vap(ki),delt,tavi, deltz, gz(ki) !compute approximated height of pressure levels
integer k,km,levs(ki),iisn
real coe2
real z(ko),dz

character dirin*100,area(6)*4,dirout*100,fnm*50,path*100
character vnm(23)*10,lvls*5
real indata(23,ki),press(ki),tv(ki)
real outdata(23,ko)
real rd,cp,g,rv,hlat

data rd,cp,g,rv,hlat /287.,1005.,9.81,461.,2.5e6/
data press /1000.,925.,850.0,700.0,600.0,500.0,400.0,300.0,250.0,200.0,150.0,100.0,70.0,50.0,30.0,20.0,10.0/
integer i,j,m,n,kk

area(1)='PRD'
area(2)='MLYR'
area(3)='NPC'
area(4)='NEC'
area(5)='ETP'
area(6)='WTP'

dirin='Z:\DATA\LargeScale\79-12\PHD2\Rainfall\'
dirout='Z:\DATA\LargeScale\79-12\PHD2\Interplot\'


  Z(1)=0  !!! first level 0 km,dz=0.24km 
  dz=240  !!! m for adjusting the cloudsat data
  do i=2,ko
   z(i)=Z(i-1)+dz
  enddo
  print*, z(ko),ko
do i=1,6
  fnm=trim(area(i))//'_meanprofile_Events.txt'
  path=trim(dirin)//trim(fnm)
  open(10,file=trim(path))
  read(10,90)lvls,(vnm(j),j=1,23)
  do k=1,ki
   read(10,92)levs(k),(indata(j,k),j=1,23)
  enddo
  tv(:)=indata(10,:)
!compute approximated height of pressure levels:
    zin(1)=0.
     do k=2,ki
       km=k-1
        tempk =tv(k)  ! * (1.+.6e-3*vap(k ))   !!!!Vap  vapor mixing  kg/kg
         tempkm=tv(km)  ! * (1.+.6e-3*vap(km)) ! !!! tempk temmkm TV
          delt=tempk-tempkm
           if(delt.gt.1.e-4) then
             tavi=alog(tempk/tempkm)/delt
            else
             tavi=1./tempk
            endif
           deltz=-rd/(tavi*g) * alog(press(k)/press(km))
           zin(k)=zin(km)+deltz
		   print*, zin(k)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      enddo

!!!! interplotting -----------------------
!
!
      do 64 k=2,ko
       do kk=2,ki
        iisn=kk-1
!	   print*,zin(kk),z(k),k,kk
        if(zin(kk).ge.z(k)) go to 665
       enddo
!       print*,' *** input sounding does not go high enough. stop.'
!              stop 'sounding'
 665   continue 
       coe2=(z(k)-zin(iisn))/(zin(iisn+1)-zin(iisn))
!!! the belows are the variables needed to be interplotted --------
        do n=1,23
          outdata(n,k)=coe2*indata(n,iisn+1) + (1.-coe2)*indata(n,iisn)
	    enddo     
 64   continue
     close(10)
!!!!!　　　output the interplotted data     
    fnm=trim(area(i))//'_meanprofile_Events_inter.txt'
    path=trim(dirout)//trim(fnm)
    open(20,file=trim(path))
    write(20,90)lvls,(vnm(j),j=1,23)
    do k=1,ko
     write(20,91)z(k)/1000,(outdata(n,k),n=1,23)
	enddo


enddo  !!!　i loop

90    format(1X,A5,23(1X,A10))
91    format(1X,F12.4,23(1X,F12.6))
92    format(1X,I2,23(1X,F12.6))


end