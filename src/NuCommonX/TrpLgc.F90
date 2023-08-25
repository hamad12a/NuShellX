
module SubM

use ProjParameters

implicit none
! allocatable mvector
type lv
    integer,dimension(:),allocatable:: lword
end type lv
! allocatable vector of mvectors
type la
    type(lv),dimension(:),allocatable:: lvect
end type la
! lv test
interface operator (.and.)
    module procedure land
end interface
! lv creator
interface operator (.or.)
    module procedure lor
end interface

contains 
! allocate mvector note nspart = no_spart not words
subroutine allocatelv(lvect,nspart)
implicit none
type(lv):: lvect
integer:: nwords,nspart
!  calculate words
nwords=(nspart-1)/32+1
allocate(lvect%lword(nwords))
!  note zeros vector ie sets false
lvect%lword=0
end subroutine allocatelv

! allocate mve array nspartr=no_spart on rows, nspartc= no_spart on columns
subroutine allocatela(larray,nspartr,nspartc)
implicit none
type(la):: larray
integer:: nspartr,nspartc,nwords,i
allocate(larray%lvect(nspartc))
do i=1,nspartc
call allocatelv(larray%lvect(i),nspartr)
end do
end subroutine allocatela

pure subroutine zerolv(lvec)
implicit none
type(lv),intent(inout):: lvec
lvec%lword=0
end subroutine zerolv

pure subroutine onelv(lvec)
implicit none
type(lv),intent(inout):: lvec
lvec%lword=z'ffffffff'
end subroutine onelv

pure subroutine zerola(larr,nsc)
implicit none
type(la),intent(inout):: larr
integer,intent(in):: nsc
integer:: i
do i=1,nsc
larr%lvect(i)%lword=0
end do
end subroutine zerola


pure function lor(a,lvec)
! set bit a in lv
implicit none
type(lv),intent(in):: lvec
integer,intent(in):: a
integer:: iw,ib
type(lv):: lor
iw=(a-1)/32 + 1
ib=a - (iw-1)*32
lor=lvec
lor%lword(iw)=ibset(lor%lword(iw),ib-1)
end function lor
	
pure subroutine crproj(laa,lac,lar,nsr,nsc)
implicit none
type(la),intent(in):: laa
type(lv),intent(inout):: lac,lar
integer,intent(in)::nsc,nsr
integer:: c,r,iw,ib
logical:: nzero
lac%lword=0
lar%lword=0
do c=1,nsc
    do r=1,nsr
        iw=(r-1)/32 + 1
        ib=r - (iw-1)*32
		if (btest(laa%lvect(c)%lword(iw),ib-1)) then
		    lar%lword(iw)=ibset(lar%lword(iw),ib-1)
		end if
	end do
end do
do c=1,nsc
    nzero=.false.
    do iw=1,(nsr-1)/32 +1
        if (laa%lvect(c)%lword(iw)/=0) then
            nzero=.true.
            exit
        end if
    end do
	if (nzero) then
	    iw=(c-1)/32 + 1
	    ib=c - (iw-1)*32
	    lac%lword(iw)=ibset(lac%lword(iw),ib-1)
	end if
end do
end subroutine crproj 

subroutine lwrite(dev,nplr,nsr,nsc)
implicit none
type(la),intent(in):: nplr
integer,intent(in)::nsr,nsc,dev
integer::i,j,iw,ib
write(dev,*)
do i=1,nsr
    iw=(i-1)/32+1
    ib=i-(iw-1)*32
    do j=1,nsc
        if (btest(nplr%lvect(j)%lword(iw),ib-1) ) then
            write(dev,'(i1)',advance='no') 1
        else    
            write(dev,'(i1)',advance='no') 0
        end if
    end do
    write(dev,*)
end do
end subroutine lwrite

pure subroutine laor(nplro,nplrb,nsr,nsc)
implicit none
type(la),intent(in):: nplrb
type(la),intent(inout):: nplro
integer,intent(in)::nsr,nsc
integer::i,j,iw,ib
do i=1,nsr
    iw=(i-1)/32+1
    ib=i-(iw-1)*32
    do j=1,nsc
        if (btest(nplrb%lvect(j)%lword(iw),ib-1) ) nplro%lvect(j)%lword(iw)=&
               ibset(nplro%lvect(j)%lword(iw),ib-1)
    end do
end do
end subroutine laor

pure function land(a,lvec)
implicit none
type(lv),intent(in):: lvec
integer,intent(in):: a
integer:: iw,ib
logical:: land
iw=(a-1)/32 + 1
ib=a - (iw-1)*32
land=btest(lvec%lword(iw),ib-1)
end function land


pure subroutine condim(nsa,nsb,blocb,nJTal,nJTar,nJTbl,dar,dal,dbl,laco,laro,lbro)
implicit none
integer,intent(in):: nsa,nsb
integer,intent(inout):: dar,dbl,dal
type(lv),intent(in):: laco,laro,lbro
integer,dimension(:),intent(in):: nJTar,nJTal,nJTbl
integer,dimension(:),intent(inout):: blocb
integer:: i
! get dimension of sum of a states conditional
dar=0
do i=1,nsa
if (land(i,laco)) dar=dar+nJTar(i)
end do
dal=0
do i=1,nsa
if (land(i,laro)) dal=dal+nJTal(i)
end do
dbl=0
do i=1,nsb
blocb(i)=dbl*dar
if (land(i,lbro)) dbl=dbl+nJTbl(i)
end do
end subroutine condim

pure subroutine single_productm(bveci,posi,meV,nsa,nsb,nJTar,nJTbl,&
                     nJTbr,nplrb,trdb,posb,blocb,dar,dal,dbl,&
                     laco,lbco,laro,lbro,zero,BUFFER1,BUFFER3,ntrb,locb)
implicit none
integer,intent(in)::posi,nsa,nsb,ntrb
integer,intent(in):: dar,dbl,dal
real(kind=rkw),dimension(:),intent(inout):: BUFFER1,BUFFER3
real(kind=rkw),dimension(:),intent(in):: bveci,trdb
integer,dimension(:),intent(in):: nJTar,nJTbr,nJTbl,blocb,locb,posb
type(lv),intent(in):: laco,lbco,laro,lbro
type(la),dimension(:),intent(in):: nplrb
real(kind=rkw),dimension(:),intent(in):: meV
logical,intent(in):: zero
real(kind=rkw):: me
integer:: i,j,k,p,pp,l,ps,us,pu,pt,pm,dp,ntb
p=posi
pp=0
!loop over b partitions
do j=1,nsb
!loop over bstates
    do k=1,nJTbr(j)
!loop over a partitions
        do i=1,nsa
!check partition pair
            if (land(j,lbco).and.land(i,laco)) then
!copy one row os a states
                BUFFER1(pp+1:pp+nJTar(i))=bveci(p+1:p+nJTar(i)) 
!increment copy to pointer conditional 
                pp=pp+nJTar(i)
            end if
!increment copy from pointer unconditional
            p=p+nJTar(i)
        end do
    end do
end do
dp=dbl*dar
if (zero) BUFFER3(1:dp)=real(0,rkw)
do ntb=1,ntrb
pt=posb(locb(ntb))
!loop over partitions br of column of partition matrix 
us=0
do j=1,nsb
!loop over partitions bl of row of partition matrix
do i=1,nsb
! check partitions
ps=blocb(i)
if (land(i,nplrb(locb(ntb))%lvect(j))) then
    pu=us
! loop over states br
    do k=1,njTbr(j)
! reset bl pointer
        pm=ps
        do l=1,nJTbl(i)
            pt=pt+1
! loop over states bl increment bl pointer
! increment trdb pointer and get matrix element
            me=trdb(pt)*meV(ntb)
! multiply vector of ar states for this br and add to bl row
            BUFFER3(pm+1:pm+dar)=BUFFER3(pm+1:pm+dar)+me*BUFFER1(pu+1:pu+dar)
! increment bl pointer
            pm=pm+dar
        end do
! end loop over bl , increment br pointer
        pu=pu+dar
    end do
end if
end do
if (land(j,lbco)) us=us+dar*nJTbr(j)
end do
end do
end subroutine single_productm


pure subroutine single_product1(bveco,poso,nsa,nsb,nJTal,nJTar,nJTbl,&
                     nplra,trda,posa,dar,dal,dbl,&
                     laco,laro,lbro,BUFFER2,&
                     BUFFER3,BUFFER4)
implicit none
integer,intent(in)::poso,nsa,nsb,dar,dal,dbl,posa
real(kind=rkw),dimension(:),intent(inout):: BUFFER2,BUFFER3,BUFFER4,bveco
real(kind=rkw),dimension(:),intent(in):: trda
integer,dimension(:),intent(in):: nJTar,nJTal,nJTbl
type(lv),intent(in):: laco,laro,lbro
type(la),intent(in):: nplra
real(kind=rkw):: me
integer:: i,j,k,p,pp,l,ps,us,pu,pt,pm,dp,nta
BUFFER4(1:dal*dbl)=real(0,rkw)
pp=dar*dbl
!set  position to zero
p=0
!loop over used (bl) states
do j=1,dbl
    BUFFER2(j:pp:dbl)=BUFFER3(p+1:p+dar)
    p=p+dar
end do
pt=posa
ps=0
us=0
!loop over partitions br of column of partition matrix 
do j=1,nsa
ps=0
!loop over partitions bl of row of partition matrix
do i=1,nsa
! check partitions
if (land(i,nplra%lvect(j))) then
    pu=us
! loop over states bl
    do k=1,njTar(j)
! reset br pointer
        pm=ps
        do l=1,nJTal(i)
            pt=pt+1
! loop over states bl increment bl pointer
! increment trdb pointer and get matrix element
            me=trda(pt)
! multiply vector of ar states for this br and add to bl row
            BUFFER4(pm+1:pm+dbl)=BUFFER4(pm+1:pm+dbl)+me*BUFFER2(pu+1:pu+dbl)
! increment bl pointer
            pm=pm+dbl
        end do
! end loop over br , increment bl pointer
        pu=pu+dbl
    end do
!save bbl pointer    
end if
!end condition
if (land(i,laro)) ps=ps+dbl*nJTal(i)
end do
if (land(j,laco)) us=us+dbl*nJTar(j)
end do
!set  position to zero
p=0
pp=dal*dbl
!loop over used (bl) states
do j=1,dal
    BUFFER2(j:pp:dal)=BUFFER4(p+1:p+dbl)
    p=p+dbl
end do
p=poso
pp=0
!loop over b partitions
do j=1,nsb
!loop over bstates
    do k=1,nJTbl(j)
!loop over a partitions
        do i=1,nsa
!check partition pair
            if (land(j,lbro).and.land(i,laro)) then
!copy one row os a states
                bveco(p+1:p+nJTal(i))= bveco(p+1:p+nJTal(i)) + BUFFER2(pp+1:pp+nJTal(i)) 
!increment copy to pointer conditional 
                pp=pp+nJTal(i)
            end if
!increment copy from pointer unconditional
            p=p+nJTal(i)
        end do
    end do
end do
end subroutine single_product1

   
pure subroutine  double_product(bveco,bveci,poso,posi,posac,posbd,trdac,trdbd,&
                    nsa,nsc,nsb,nsd,nJTa,nJTc,nJTb,nJTd,&
                        nine,nplrac,nplrbd,BUFFER1,BUFFER2,BUFFER3)
implicit none
integer,intent(in)::poso,posi,posac,posbd,nsa,nsb,nsc,nsd
real(kind=rkw),dimension(:),intent(inout):: bveci
real(kind=rkw),dimension(:),intent(in):: trdac,trdbd
real(kind=rkw),dimension(:),intent(inout):: BUFFER1,BUFFER3,BUFFER2,bveco
type(la),intent(in):: nplrac,nplrbd
integer,dimension(:),intent(in):: nJTa,nJTb,nJTc,nJTd
real(kind=rkw),intent(in):: nine
real(kind=rkw):: me
integer:: i,j,da,db,ps,p,pm,us,pu,pt,k,l,dp,dg,dc,dd
da=0
do i=1,nsa
da=da+nJTa(i)
end do
db=0
do i=1,nsb
db=db+nJTb(i)
end do
dc=0
do i=1,nsc
dc=dc+nJTc(i)
end do
dd=0
do i=1,nsd
dd=dd+nJTd(i)
end do

dp=db*dc
BUFFER2(1:dp)=real(0,rkw)

pt=posbd
us=0
do j=1,nsd
ps=0
do i=1,nsb
! check partitions
if (land(i,nplrbd%lvect(j))) then
    pu=us
    do k=1,njTd(j)
        pm=ps
        do l=1,nJTb(i)
            pt=pt+1
            me=trdbd(pt)*nine
            BUFFER2(pm+1:pm+dc)=BUFFER2(pm+1:pm+dc)+&
                        me*bveci(posi+pu+1:posi+pu+dc)
            pm=pm+dc
        end do
        pu=pu+dc
    end do
end if
ps=ps+dc*nJTb(i)
end do
us=us+dc*nJTd(j)
end do

p=0
do j=1,db
    BUFFER1(j:dp:db)=BUFFER2(p+1:p+dc)
    p=p+dc
end do

dp=da*db
BUFFER3(1:dp)=real(0,rkw)

pt=posac
us=0
do j=1,nsc
ps=0
do i=1,nsa
! check partitions
if (land(i,nplrac%lvect(j))) then
    pu=us
    do k=1,njTc(j)
        pm=ps
        do l=1,nJTa(i)
            pt=pt+1
            me=trdac(pt)
            BUFFER3(pm+1:pm+db)=BUFFER3(pm+1:pm+db)+me*BUFFER1(pu+1:pu+db)
            pm=pm+db
        end do
        pu=pu+db
    end do
end if
ps=ps+db*nJTa(i)
end do
us=us+db*nJTc(j)
end do

p=0
!loop over used (a) states
do j=1,da
    BUFFER1(j:dp:da)=BUFFER3(p+1:p+db)
    p=p+db
end do
bveco(poso+1:poso+dp)=bveco(poso+1:poso+dp)+BUFFER1(1:dp)

end subroutine double_product

  
pure subroutine mulmtx(bveco,bveci,posb,posx,mtx,nsa,nsb,nJTa,nJTb,&
                            nplr,BUFFER1,BUFFER3,trans)
implicit none
integer,intent(in)::posb,posx,nsa,nsb
real(kind=rkw),dimension(:),intent(inout):: bveci
real(kind=rkw),dimension(:),intent(in):: mtx
real(kind=rkw),dimension(:),intent(inout):: BUFFER1,BUFFER3,bveco
type(la),intent(in):: nplr
logical,intent(in):: trans
integer,dimension(:),intent(in):: nJTa,nJTb
real(kind=rkw):: me
integer:: i,j,da,db,ps,p,pm,us,pu,pt,k,l,dp,dg
da=0
do i=1,nsa
da=da+nJTa(i)
end do
db=0
do i=1,nsb
db=db+nJTb(i)
end do

dp=da*db

if (trans) then
p=posb
do j=1,db
    BUFFER1(j:dp:db)=bveci(p+1:p+da)
    p=p+da
end do
end if

if (trans) BUFFER3(1:dp)=real(0,rkw)
if (trans) then
pt=posx
us=0
do j=1,nsa
ps=0
do i=1,nsa
! check partitions
if (land(i,nplr%lvect(j))) then
    pu=us
    do k=1,njTa(j)
        pm=ps
        do l=1,nJTa(i)
            pt=pt+1
            me=mtx(pt)
            BUFFER3(pu+1:pu+db)=BUFFER3(pu+1:pu+db)+me*BUFFER1(pm+1:pm+db)
            if (i/=j) BUFFER3(pm+1:pm+db)=BUFFER3(pm+1:pm+db)+me*BUFFER1(pu+1:pu+db)
            pm=pm+db
        end do
        pu=pu+db
    end do
end if
ps=ps+db*nJTa(i)
end do
us=us+db*nJTa(j)
end do
end if

if (.not.trans) then
pt=posx
us=0
do j=1,nsb
ps=0
do i=1,nsb
! check partitions
if (land(i,nplr%lvect(j))) then
    pu=us
    do k=1,njTb(j)
        pm=ps
        do l=1,nJTb(i)
            pt=pt+1
            me=mtx(pt)
            bveco(posb+pu+1:posb+pu+da)=bveco(posb+pu+1:posb+pu+da)+&
                        me*bveci(posb+pm+1:posb+pm+da)
            if (i/=j) bveco(posb+pm+1:posb+pm+da)=bveco(posb+pm+1:posb+pm+da)+&
                        me*bveci(posb+pu+1:posb+pu+da)
            pm=pm+da
        end do
        pu=pu+da
    end do
end if
ps=ps+da*nJTb(i)
end do
us=us+da*nJTb(j)
end do
end if


if (trans) then
p=0
!loop over used (a) states
do j=1,da
    BUFFER1(j:dp:da)=BUFFER3(p+1:p+db)
    p=p+db
end do
bveco(posb+1:posb+dp)=bveco(posb+1:posb+dp)+BUFFER1(1:dp)
end if

end subroutine mulmtx


end module SubM