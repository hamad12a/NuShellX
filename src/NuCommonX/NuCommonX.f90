module XParameters

integer,parameter::  isl=3,isr=4,ilm=5,ijl=1,ijr=2,inl=6,inr=8,imr=10,iml=11,ilt=12
integer,parameter::  isx=13,itl=14,itr=15,inn=16,izl=17,izr=18

integer,parameter::  trmul=4,trdiv=3,trvmul=4,trvdiv=3
real,parameter:: scale=2.0,hlim=.2

end module XParameters
module XTParameters

integer,parameter::  xsl=5,xsr=6,xsx=7,xlm=8,xjl=1,xjr=2,xlt=9,xJm=10,xTm=11,xtl=3,xtr=4,xnn=12

end module XTParameters

module OpmParameters

integer,parameter:: max_Opm=30000000


end module OpmParameters

module LnzParameters

integer,parameter:: noiter=300,max_trdb=1000,max_trda=1000

end module LnzParameters

module TrdChar

use ProjParameters
use LanczosParameters
use OrderParameters
use OperParameters
use InputNuShell
use MatrixParameters

implicit none

integer,parameter:: mask=z'ff',mask4=z'fe'

type tablebme
    integer(kind=2),dimension(:),allocatable:: b
    integer(kind=2),dimension(:),allocatable:: h
    integer(kind=2):: bn
    real(kind=rkw),dimension(:),allocatable:: me
end type tablebme

type tablear
    integer(kind=2),dimension(:),allocatable:: a
    type(tablebme),dimension(:),allocatable:: bme
end type tablear

type mat
    real(kind=rkw),dimension(:,:),allocatable:: t
end type mat

type vect
    real(kind=rkw),dimension(:),allocatable:: v
end type vect


character(len=1):: typea='?',typeb='?'
character(len=1):: typec='?',typed='?'

contains

    
             
    function no_shellsp()
    implicit none
    integer:: no_shellsp
    no_shellsp=0
    do shn=1,no_shells
    if (nucleon(shn)=='p') then
        no_shellsp=no_shellsp+1
    else
        exit
    end if
    end do
    end function no_shellsp
    
    function min_shellsa()
    implicit none
    integer:: min_shellsa
    if (typea=='p') then
        min_shellsa=1
    else
        min_shellsa=1+no_shellsp()
    end if
    end  function min_shellsa

    function min_shellsb()
    implicit none
    integer:: min_shellsb
    if (typeb=='p') then
        min_shellsb=1
    else
        min_shellsb=1+no_shellsp()
    end if
    end  function min_shellsb
    
    function max_shellsa()
    implicit none
    integer:: max_shellsa
    if (typea=='p') then
        max_shellsa=no_shellsp()
    else
        max_shellsa=no_shells
    end if
    end  function max_shellsa
    
    function max_sj2ua()
    implicit none
    integer:: max_sj2ua
    integer:: max_ashells,min_ashells,i
    max_ashells=max_shellsa()
    min_ashells=min_shellsa()
    max_sj2ua=0
    do i=min_ashells,max_ashells
        max_sj2ua=max(max_sj2ua,sj2(i))
    end do
    end function max_sj2ua
    
    function max_sj2ub()
    implicit none
    integer:: max_sj2ub
    integer:: max_bshells,min_bshells,i
    max_bshells=max_shellsb()
    min_bshells=min_shellsb()
    max_sj2ub=0
    do i=min_bshells,max_bshells
        max_sj2ub=max(max_sj2ub,sj2(i))
    end do
    end function max_sj2ub
    
    function max_shellsb()
    implicit none
    integer:: max_shellsb
    if (typeb=='p') then
        max_shellsb=no_shellsp()
    else
        max_shellsb=no_shells
    end if
    end  function max_shellsb
    
    pure function tfind(sarray,n,string)
    implicit none
    integer,intent(in):: n
    character(len=8),intent(in):: string
    integer:: tfind
    character(len=8),dimension(*),intent(in):: sarray
    integer:: in
    real:: rn,drn
    if (n==0) then
        tfind=0
        return
    end if
    if (sarray(n)(1:8)==string(1:8)) then 
        tfind=n
        return
    end if
    rn=n
    rn=rn/2.0
    drn=rn
    do
        in=nint(rn)
        if (sarray(in)(1:8)==string(1:8)) then
            tfind=in
            exit
        end if
        drn=drn/2.0
        if (abs(drn) < 0.5) then
            tfind=0
            exit
        end if
        if (lgt(string(1:8),sarray(in)(1:8))) then
            rn=rn+drn
        else
            rn=rn-drn
        end if
    end do
    end function tfind
    
    elemental function char2(i)
    implicit none
    integer,intent(in)::i
    character(len=2)::char2
    char2(1:1)=char(iand(mask,i))
    char2(2:2)=char(iand(mask,ishft(i,-8)))
    end function char2
    
    elemental function ichar2(c2)
    implicit none
    character(len=2),intent(in):: c2
    integer:: ichar2
    ichar2=ichar(c2(1:1))
    ichar2=ior(ishft(ichar(c2(2:2)),8),ichar2)
    end function ichar2
     
    function ichara(chra,n)
    implicit none
    integer,intent(in)::n
    character(len=n),intent(in):: chra
    integer,dimension(8):: ichara
    integer::i
    ichara=0
    do i=1,n
        ichara(i)=ichar(chra(i:i))
    end do
    end function ichara
       
    function chara(iarray,n)
    implicit none
    integer,intent(in)::n
    character(len=n):: chara
    integer,dimension(n),intent(in):: iarray
    integer::i
    
    do i=1,n
        chara(i:i)=char(iarray(i))
    end do
    end function chara
    
end module TrdChar    

module TrdChar3

use TrdChar

implicit none

type trx
    integer,dimension(5):: qux
end type trx

type projvec
    real(kind=rkw),dimension(:,:),allocatable:: coef
end type

contains


    pure function tfind3(sarray,n,string)
    implicit none
    integer,intent(in):: n
    character(len=12),intent(in):: string
    integer:: tfind3
    character(len=12),dimension(*),intent(in):: sarray
    integer:: in
    real:: rn,drn
    if (n==0) then
        tfind3=0
        return
    end if
    if (sarray(n)(1:12)==string(1:12)) then 
        tfind3=n
        return
    end if
    rn=n
    rn=rn/2.0
    drn=rn
    do
        in=nint(rn)
        if (sarray(in)(1:12)==string(1:12)) then
            tfind3=in
            exit
        end if
        drn=drn/2.0
        if (abs(drn) < 0.5) then
            tfind3=0
            exit
        end if
        if (lgt(string(1:12),sarray(in)(1:12))) then
            rn=rn+drn
        else
            rn=rn-drn
        end if
    end do
    end function tfind3
    

    pure function tfind4(sarray,n,string)
    implicit none
    integer,intent(in):: n
    character(len=18),intent(in):: string
    integer:: tfind4
    character(len=18),dimension(*),intent(in):: sarray
    integer:: in
    real:: rn,drn
    if (n==0) then
        tfind4=0
        return
    end if
    if (sarray(n)(1:18)==string(1:18)) then 
        tfind4=n
        return
    end if
    rn=n
    rn=rn/2.0
    drn=rn
    do
        in=nint(rn)
        if (sarray(in)(1:18)==string(1:18)) then
            tfind4=in
            exit
        end if
        drn=drn/2.0
        if (abs(drn) < 0.5) then
            tfind4=0
            exit
        end if
        if (lgt(string(1:18),sarray(in)(1:18))) then
            rn=rn+drn
        else
            rn=rn-drn
        end if
    end do
    end function tfind4
    
    pure function tfind2(sarray,n,string)
    implicit none
    integer,intent(in):: n
    character(len=8),intent(in):: string
    integer:: tfind2
    character(len=12),dimension(*),intent(in):: sarray
    integer:: in
    real:: rn,drn
    if (n==0) then
        tfind2=0
        return
    end if
    if (sarray(n)(1:8)==string(1:8)) then 
        tfind2=n
        return
    end if
    rn=n
    rn=rn/2.0
    drn=rn
    do
        in=nint(rn)
        if (sarray(in)(1:8)==string(1:8)) then
            tfind2=in
            exit
        end if
        drn=drn/2.0
        if (abs(drn) < 0.5) then
            tfind2=0
            exit
        end if
        if (lgt(string(1:8),sarray(in)(1:8))) then
            rn=rn+drn
        else
            rn=rn-drn
        end if
    end do
    end function tfind2
    
    elemental function char3(i)
    implicit none
    integer,intent(in)::i
    character(len=3)::char3
    char3(1:1)=char(iand(mask,i))
    char3(2:2)=char(iand(mask,ishft(i,-8)))
    char3(3:3)=char(iand(mask,ishft(i,-16)))
    end function char3
    
    elemental function ichar3(c3)
    implicit none
    character(len=3),intent(in):: c3
    integer:: ichar3
    ichar3=0
    ichar3=ichar(c3(1:1))
    ichar3=ior(ishft(ichar(c3(2:2)),8),ichar3)
    ichar3=ior(ishft(ichar(c3(3:3)),16),ichar3)
    end function ichar3
     
    function ichara3(chra,n)
    implicit none
    integer,intent(in)::n
    character(len=n),intent(in):: chra
    integer,dimension(12):: ichara3
    integer::i
    ichara3=0
    do i=1,n
        ichara3(i)=ichar(chra(i:i))
    end do
    end function ichara3
       
    function chara3(iarray,n)
    implicit none
    integer,intent(in)::n
    character(len=n):: chara3
    integer,dimension(n),intent(in):: iarray
    integer::i
    
    do i=1,n
        chara3(i:i)=char(iarray(i))
    end do
    end function chara3
    
end module TrdChar3    

module NextBit

use Machine
use Parameters
use Finds

implicit none

integer(kind=ik),private:: ii
integer(kind=ik),parameter:: mnn = -2**(nbits-1)
integer(kind=ik),dimension(nbits-1):: log2 = (/(2**(ii-1),ii=1,nbits-1)/)

contains

    
    pure function next_bit_lx(i)
        implicit none
        integer,intent(in):: i
        integer:: next_bit_lx
        integer(kind=ik):: ia,iik

	iik=i
        
        if (iik == mnn) then
            next_bit_lx=nbits
!            i=ibclr(i,nbits-1)
            return
        end if
        
        ia=iand(iik,-iik)
        next_bit_lx=int(ifind(log2,nbits-1,ia),4)
!        if (next_bit_lx /= 0) i=ibclr(i,next_bit_lx-1)
        
    end function next_bit_lx
    

end module NextBit

module SubM

use ProjParameters
!use Mvector
use NextBit

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

logical(kind=1),dimension(:),allocatable:: lcoef
logical(kind=1),dimension(:),allocatable::trunc,bas_log
logical(kind=1),dimension(:,:),allocatable::lbuffer1,lbuffer2
integer,dimension(:),allocatable:: para,parb,da0,db0,dc0,dd0
integer,dimension(:,:),allocatable::snJTa0,snJTb0,snJTc0,snJTd0
integer::  max_da0=0,max_db0=0
logical::  truncate

contains 

subroutine set_ab(nsa,mnJa,mxJa,deJa,nsb,mnJb,mxJb,deJb,nJTa,nJTb)
implicit none
integer,intent(in):: nsa,nsb,mnJa,mnJb,mxJa,mxJb,deJa,deJb
integer,dimension(:,0:),intent(in):: nJTa, nJTb 
integer::i,j,ja,jb
allocate(da0(mnJa/2:mxJa/2))
allocate(db0(mnJb/2:mxJb/2))
da0=0
db0=0
allocate(snJTa0(nsa,mnJa/2:mxJa/2))
allocate(snJTb0(nsb,mnJb/2:mxJb/2))
snJTa0=0
snJTb0=0
do ja=mnJa,mxJa,deJa
    do i=1,nsa
        da0(ja/2)=da0(ja/2)+nJTa(i,ja/2)
        if (i>1) snJTa0(i,ja/2)=snJTa0(i-1,ja/2)+nJTa(i-1,ja/2)
    end do
    max_da0=max(max_da0,da0(ja/2))
end do
do jb=mnJb,mxJb,deJb
    do j=1,nsb
        db0(jb/2)=db0(jb/2)+nJTb(j,jb/2)
        if (j>1) snJTb0(j,jb/2)=snJTb0(j-1,jb/2)+nJTb(j-1,jb/2)
    max_db0=max(max_db0,db0(jb/2))
    end do
end do
end subroutine set_ab
     
subroutine set_cd(nsc,mnJc,mxJc,deJc,nsd,mnJd,mxJd,deJd,nJTc,nJTd)
implicit none
integer,intent(in):: nsc,nsd,mnJc,mnJd,mxJc,mxJd,deJc,deJd
integer,dimension(:,0:),intent(in):: nJTc, nJTd 
integer::i,j,jc,jd
allocate(dc0(mnJc/2:mxJc/2))
allocate(dd0(mnJd/2:mxJd/2))
dc0=0
dd0=0
allocate(snJTc0(nsc,mnJc/2:mxJc/2))
allocate(snJTd0(nsd,mnJd/2:mxJd/2))
snJTc0=0
snJTd0=0
do jc=mnJc,mxJc,deJc
    do i=1,nsc
        dc0(jc/2)=dc0(jc/2)+nJTc(i,jc/2)
        if (i>1) snJTc0(i,jc/2)=snJTc0(i-1,jc/2)+nJTc(i-1,jc/2)
    end do
end do
do jd=mnJd,mxJd,deJd
    do j=1,nsd
        dd0(jd/2)=dd0(jd/2)+nJTd(j,jd/2)
        if (j>1) snJTd0(j,jd/2)=snJTd0(j-1,jd/2)+nJTd(j-1,jd/2)
    end do
end do
end subroutine set_cd
     
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

subroutine deallocatelv(lvect)
implicit none
type(lv):: lvect
!  calculate words
deallocate(lvect%lword)
!  note zeros vector ie sets false
end subroutine deallocatelv

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

subroutine deallocatela(larray,nspartc)
implicit none
type(la):: larray
integer:: nspartc,i
do i=1,nspartc
call deallocatelv(larray%lvect(i))
end do
deallocate(larray%lvect)
end subroutine deallocatela

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

pure subroutine next_ij(laa,i,j,nsi,nsj,word,iw,ib,done)
implicit none
type(la),intent(in):: laa
integer,intent(inout):: i,j,word,iw,ib
integer,intent(in):: nsi,nsj
integer:: max_iw
logical,intent(inout):: done
    done=.false.
    max_iw=(nsi-1)/32+1
10  if (word==0) then
        iw=iw+1
        if (iw>max_iw) then
            j=j+1
            if (j>nsj) then
                done=.true.
                return
            end if
            iw=1
        end if
        word=laa%lvect(j)%lword(iw)
        if (word==0) go to 10        
    end if
    ib=next_bit_lx(word)
    i=ib+32*(iw-1)
    if (i==0) go to 10
    word=ibclr(word,ib-1)
    return
end subroutine next_ij
        
	
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

subroutine single_productm(bveci,posi,meV,nsa,nsb,nJTar,nJTbl,&
                     nJTbr,nplrb,trdb,posb,blocb,dar,dal,dbl,&
                     laco,lbco,laro,lbro,zero,BUFFER1,BUFFER3,ntrb,locb,nt)
implicit none
integer,intent(in)::posi,nsa,nsb,ntrb,nt
integer,intent(in):: dar,dbl,dal
real(kind=rkw),dimension(:),intent(inout):: BUFFER1,BUFFER3
real(kind=rkw),dimension(:),intent(in):: bveci,trdb
integer,dimension(:),intent(in):: nJTar,nJTbr,nJTbl,blocb,locb,posb
type(lv),intent(in):: laco,lbco,laro,lbro
type(la),dimension(:),intent(in):: nplrb
real(kind=rkw),dimension(:),intent(in):: meV
logical,intent(in):: zero
logical:: done
real(kind=rkw):: me
integer,dimension(nsb):: snJTb
integer:: i,j,k,p,pp,l,ps,us,pu,pt,pm,dp,ntb,word,iw,ib
p=posi
pp=0
do j=1,nsb
    do k=1,nJTbr(j)
        do i=1,nsa
            if (land(j,lbco).and.land(i,laco).and.nJTar(i)/=0) then
                if (truncate) then
                BUFFER1(pp+1:pp+nJTar(i))=real(0,rk) 
                lbuffer1(pp+1:pp+nJTar(i),nt)=.false. 
                where(lcoef(p+1:p+nJTar(i))) BUFFER1(pp+1:pp+nJTar(i))=bveci(p+1:p+nJTar(i)) 
                where(lcoef(p+1:p+nJTar(i))) lbuffer1(pp+1:pp+nJTar(i),nt)=lbuffer1(pp+1:pp+nJTar(i),nt)&
                                            .or.trunc(p+1:p+nJTar(i)) 
                else
                BUFFER1(pp+1:pp+nJTar(i))=bveci(p+1:p+nJTar(i)) 
                end if
                pp=pp+nJTar(i)
            end if
            p=p+nJTar(i)
        end do
    end do
end do
dp=dbl*dar
if (zero) BUFFER3(1:dp)=real(0,rkw)
if (zero.and.truncate) lbuffer2(1:dp,nt)=.false.
snJTb(1)=0
do j=2,nsb
    if (land(j-1,lbco)) then
        snJTb(j)=snJTb(j-1)+dar*nJTbr(j-1)
    else
        snJTb(j)=snJTb(j-1)
    end if
end do
do ntb=1,ntrb
pt=posb(locb(ntb))
i=1
j=1
iw=1
ib=1
word=nplrb(locb(ntb))%lvect(j)%lword(iw)
do 
    call next_ij(nplrb(locb(ntb)),i,j,nsb,nsb,word,iw,ib,done)
    if (done) exit
    pu=snJTb(j)
    do k=1,njTbr(j)
        pm=blocb(i)
!DEC$ VECTOR ALWAYS
        do l=1,nJTbl(i)
            pt=pt+1
            me=trdb(pt)*meV(ntb)
            if (abs(me)>=0.0000001) then
            if (truncate)then
            where(lbuffer1(pu+1:pu+dar,nt)) BUFFER3(pm+1:pm+dar)=BUFFER3(pm+1:pm+dar)+&
                        me*BUFFER1(pu+1:pu+dar)
            where(lbuffer1(pu+1:pu+dar,nt)) lbuffer2(pm+1:pm+dar,nt)=lbuffer2(pm+1:pm+dar,nt)&
                        .or.lbuffer1(pu+1:pu+dar,nt)
            else
            BUFFER3(pm+1:pm+dar)=BUFFER3(pm+1:pm+dar)+me*BUFFER1(pu+1:pu+dar)
            end if 
            end if 
            pm=pm+dar
        end do
        pu=pu+dar
    end do
end do
end do
end subroutine single_productm


subroutine single_product1(bveco,poso,nsa,nsb,nJTal,nJTar,nJTbl,&
                     nplra,trda,posa,dar,dal,dbl,&
                     laco,laro,lbro,BUFFER2,&
                     BUFFER3,BUFFER4,nt)
implicit none
integer,intent(in)::poso,nsa,nsb,dar,dal,dbl,posa,nt
real(kind=rkw),dimension(:),intent(inout):: BUFFER2,BUFFER3,BUFFER4,bveco
real(kind=rkw),dimension(:),intent(in):: trda
integer,dimension(:),intent(in):: nJTar,nJTal,nJTbl
type(lv),intent(in):: laco,laro,lbro
type(la),intent(in):: nplra
real(kind=rkw):: me
logical::done
integer,dimension(nsa):: snJTal,snJTar
integer:: i,j,k,p,pp,l,ps,us,pu,pt,pm,dp,nta,word,iw,ib
BUFFER4(1:dal*dbl)=real(0,rkw)
pp=dar*dbl
p=0
!DEC$ VECTOR ALWAYS
do j=1,dbl
    BUFFER2(j:pp:dbl)=BUFFER3(p+1:p+dar)
    if (truncate) lbuffer1(j:pp:dbl,nt)=lbuffer2(p+1:p+dar,nt)
    p=p+dar
end do
snJTal(1)=0
do j=2,nsa
    if (land(j-1,laro)) then
        snJTal(j)=snJTal(j-1)+dbl*nJTal(j-1)
    else
        snJTal(j)=snJTal(j-1)
    end if
end do
snJTar(1)=0
do j=2,nsa
    if (land(j-1,laco)) then
        snJTar(j)=snJTar(j-1)+dbl*nJTar(j-1)
    else
        snJTar(j)=snJTar(j-1)
    end if
end do
pt=posa
i=1
j=1
iw=1
word=nplra%lvect(j)%lword(iw)
do 
    call next_ij(nplra,i,j,nsa,nsa,word,iw,ib,done)
    if (done) exit
    pu=snJTar(j)
    do k=1,nJTar(j)
        pm=snJTal(i)
!DEC$ VECTOR ALWAYS
        do l=1,nJTal(i)
            pt=pt+1
            me=trda(pt)
            if (abs(me)>=0.0000001) then
            if (truncate) then
            where(lbuffer1(pu+1:pu+dbl,nt))  BUFFER4(pm+1:pm+dbl)=BUFFER4(pm+1:pm+dbl)+ &
                                me*BUFFER2(pu+1:pu+dbl)
            else
            BUFFER4(pm+1:pm+dbl)=BUFFER4(pm+1:pm+dbl)+me*BUFFER2(pu+1:pu+dbl)
            end if
            end if
            pm=pm+dbl
        end do
        pu=pu+dbl
    end do
end do
p=0
pp=dal*dbl
!DEC$ VECTOR ALWAYS
do j=1,dal
    BUFFER2(j:pp:dal)=BUFFER4(p+1:p+dbl)
    p=p+dbl
end do
p=poso
pp=0
do j=1,nsb
    do k=1,nJTbl(j)
        do i=1,nsa
            if (land(j,lbro).and.land(i,laro).and.nJTal(i)/=0) then
                    bveco(p+1:p+nJTal(i))= bveco(p+1:p+nJTal(i)) + BUFFER2(pp+1:pp+nJTal(i)) 
                pp=pp+nJTal(i)
            end if
            p=p+nJTal(i)
        end do
    end do
end do
return
end subroutine single_product1

   
pure subroutine  double_product(bveco,bveci,poso,posi,posac,posbd,trdac,trdbd,&
                    nsa,nsc,nsb,nsd,nJTa,nJTc,nJTb,nJTd,&
                        nine,nplrac,nplrbd,BUFFER1,BUFFER2,BUFFER3,&
                        ja,jb,jc,jd)
implicit none
integer,intent(in)::poso,posi,posac,posbd,nsa,nsb,nsc,nsd,ja,jb,jc,jd
real(kind=rkw),dimension(:),intent(inout):: bveci
real(kind=rkw),dimension(:),intent(in):: trdac,trdbd
real(kind=rkw),dimension(:),intent(inout):: BUFFER1,BUFFER3,BUFFER2,bveco
type(la),intent(in):: nplrac,nplrbd
integer,dimension(:),intent(in):: nJTa,nJTb,nJTc,nJTd
real(kind=rkw),intent(in):: nine
real(kind=rkw):: me
integer:: i,j,da,db,ps,p,pm,us,pu,pt,k,l,dp,dg,dc,dd,word,iw,ib
logical:: done
da=da0(ja/2)
db=db0(jb/2)
dc=dc0(jc/2)
dd=dd0(jd/2)

dp=db*dc
BUFFER2(1:dp)=real(0,rkw)

pt=posbd
i=1
j=1
iw=1
ib=1
word=nplrbd%lvect(j)%lword(iw)
do 
    call next_ij(nplrbd,i,j,nsb,nsd,word,iw,ib,done)
    if (done) exit
    pu=dc*snJTd0(j,jd/2)
    do k=1,njTd(j)
        pm=dc*snJTb0(i,jb/2)
        do l=1,nJTb(i)
            pt=pt+1
            me=trdbd(pt)*nine
            if (me/=0.0) BUFFER2(pm+1:pm+dc)=BUFFER2(pm+1:pm+dc)+&
                        me*bveci(posi+pu+1:posi+pu+dc)
            pm=pm+dc
        end do
        pu=pu+dc
    end do
end do

p=0
do j=1,db
    BUFFER1(j:dp:db)=BUFFER2(p+1:p+dc)
    p=p+dc
end do

dp=da*db
BUFFER3(1:dp)=real(0,rkw)

pt=posac
i=1
j=1
iw=1
ib=1
word=nplrac%lvect(j)%lword(iw)
do 
    call next_ij(nplrac,i,j,nsa,nsc,word,iw,ib,done)
    if (done) exit
    pu=db*snJTc0(j,jc/2)
    do k=1,njTc(j)
        pm=db*snJTa0(i,ja/2)
        do l=1,nJTa(i)
            pt=pt+1
            me=trdac(pt)
            if (me/=0.0) BUFFER3(pm+1:pm+db)=BUFFER3(pm+1:pm+db)+me*BUFFER1(pu+1:pu+db)
            pm=pm+db
        end do
        pu=pu+db
    end do
end do

p=0
do j=1,da
    BUFFER1(j:dp:da)=BUFFER3(p+1:p+db)
    p=p+db
end do
bveco(poso+1:poso+dp)=bveco(poso+1:poso+dp)+BUFFER1(1:dp)

end subroutine double_product

  
pure subroutine mulmtx(bveco,bveci,posb,posx,mtx,nsa,nsb,nJTa,nJTb,&
                            nplr,BUFFER1,BUFFER3,trans,ja,jb)
implicit none
integer,intent(in)::posb,posx,nsa,nsb,ja,jb
real(kind=rkw),dimension(:),intent(inout):: bveci
real(kind=rkw),dimension(:),intent(in):: mtx
real(kind=rkw),dimension(:),intent(inout):: BUFFER1,BUFFER3,bveco
type(la),intent(in):: nplr
logical,intent(in):: trans
integer,dimension(:),intent(in):: nJTa,nJTb
real(kind=rkw):: me
integer:: i,j,da,db,ps,p,pm,us,pu,pt,k,l,dp,dg,word,iw,ib
logical:: done
da=da0(ja/2)
db=db0(jb/2)

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
i=1
j=1
iw=1
word=nplr%lvect(j)%lword(iw)
do 
    call next_ij(nplr,i,j,nsa,nsa,word,iw,ib,done)
    if (done) exit
    pu=db*snJTa0(j,ja/2)
    if (para(i)==para(j)) then
    do k=1,njTa(j)
        pm=db*snJTa0(i,ja/2)
        do l=1,nJTa(i)
            pt=pt+1
            me=mtx(pt)
            if (me/=0.0) BUFFER3(pu+1:pu+db)=BUFFER3(pu+1:pu+db)+me*BUFFER1(pm+1:pm+db)
            if (i/=j.and.me/=0.0) BUFFER3(pm+1:pm+db)=BUFFER3(pm+1:pm+db)+me*BUFFER1(pu+1:pu+db)
            pm=pm+db
        end do
        pu=pu+db
    end do
    end if
end do
end if

if (.not.trans) then
pt=posx
i=1
j=1
iw=1
word=nplr%lvect(j)%lword(iw)
do 
    call next_ij(nplr,i,j,nsb,nsb,word,iw,ib,done)
    if (done) exit
    pu=da*snJTb0(j,jb/2)
    if (parb(j)==parb(i)) then
    do k=1,njTb(j)
        pm=da*snJTb0(i,jb/2)
        do l=1,nJTb(i)
            pt=pt+1
            me=mtx(pt)
            if (me/=0.0) then
                        bveco(posb+pu+1:posb+pu+da)=bveco(posb+pu+1:posb+pu+da)+&
                        me*bveci(posb+pm+1:posb+pm+da)
            end if            
            if (i/=j.and.me/=0.0) then
                        bveco(posb+pm+1:posb+pm+da)=bveco(posb+pm+1:posb+pm+da)+&
                        me*bveci(posb+pu+1:posb+pu+da)
            end if            
            pm=pm+da
        end do
        pu=pu+da
    end do
    end if
end do
end if


if (trans) then
p=0
do j=1,da
    BUFFER1(j:dp:da)=BUFFER3(p+1:p+db)
    p=p+db
end do
bveco(posb+1:posb+dp)=bveco(posb+1:posb+dp)+BUFFER1(1:dp)
end if

end subroutine mulmtx


end module SubM

module TLMB
IMPLICIT NONE

REAL*8,PRIVATE:: BIN(0:99,0:99),TIN(0:99,0:99,0:99)


contains

! W D M Rae Garsington 2007. Interface to TMSH

  function moshinsky(cN,lam,sn1,sl1,sn2,sl2,m1,m2)
     implicit none
     integer,intent(in):: cN,lam,sn1,sl1,sn2,sl2,m1,m2
     real*8:: moshinsky,d
     logical,save::first=.true.
     integer:: ee,e1,e2
     
     if (first) then
        first=.false.
        call init()
     end if   
     if (2*cN+lam/=2*(sn1+sn2)+sl1+sl2) then
        moshinsky=real(0,8)
        return
     end if
     if (cN<0.or.lam<0.or.sn1<0.or.sn2<0.or.sl1<0.or.sl2<0) then
        moshinsky=real(0,8)
        return
     end if

          e1=2*sn1+sl1
          e2=2*sn2+sl2
          ee=2*cN+lam
     
     if (m1>=m2) then
          d=real(m1,8)/real(m2,8)
          moshinsky=TMB(ee,lam,0,0,e1,sl1,e2,sl2,lam,d)
     else     
          d=real(m2,8)/real(m1,8)
          moshinsky=TMB(ee,lam,0,0,e2,sl2,e1,sl1,lam,d)
     end if
     end function moshinsky
     


!** Version 1.0: May 2001 
!   Fortran 95 Version W D M rae, Garsington, 2007
!   http://knollhouse.org
!** E-mail: Gintautas_Kamuntavicius@fc.vdu.lt
!** Reference: nucl-th/0105009
!** WWW: http://www.nuclear.physics.vdu.lt 
!   Global Variables to Module BIN TIN
!** BIN-array of binomial coefficients
!** TIN-array of trinomial coefficients
!** TMB-real*8 function of HOBs:
!** Input: EE-centre of mass energy, LL-centre of mass angular moment
!**        ER-relative energy, LR-relative angular moment!
!**        E1-energy of the first particle, L1-angular moment of the first particle
!**        E2-energy of the second particle, L2-angular moment of the second particle
!**        LM-total angular moment   
!** Output: TMB-HARMONIC OSCILLATOR BRACKET

    SUBROUTINE INIT()
    IMPLICIT NONE
    REAL tim1,tim2,timd
	REAL*8 T1,skirtm
    INTEGER LAST
!c	LAST=0,1,2,3,4...
    LAST=3
!	CALL CPU_TIME(tim1)
    CALL BINOML
	CALL TRINOM
    END SUBROUTINE INIT

      SUBROUTINE BINOML
!  	  THE ARRAY OF BINOMLIAL COEFFICIENTS
!     BIN(I,J)= = I!/J!/(I-J)! 
      IMPLICIT NONE
      INTEGER I,K
    	DO I=0,99
	        BIN(I,0)=1.D0
	        BIN(I,I)=1.D0
	        DO K=1,I/2
                BIN(I,K)=DNINT(BIN(I,K-1)/DFLOAT(K)*DFLOAT(I-K+1))
	            BIN(I,I-K)=BIN(I,K)
	        END DO
        END DO
	    RETURN
      END SUBROUTINE BINOML

      INTEGER FUNCTION TRI(I,J,K)
!     TRIADIC CONDITION FOR MOMENTS I/2,J/2,K/2:
!     I+J>=K, I+K>=J, J+K>=I,
!     I/2+J/2+K/2 = INTEGER.
!     TRI=1, WHEN TRIADIC CONDITION IS FULFILLED, TRI=0 OTHERWISE	 
      IMPLICIT NONE
      INTEGER I,J,K,L
        TRI=0
    	L=I+J+K
        IF(L/2*2.NE.L) RETURN
	    L=L/2
        IF((L-I)*(L-J)*(L-K).LT.0) RETURN
        TRI=1
      RETURN
      END FUNCTION TRI

      REAL*8 FUNCTION W6(I,J,M,L,K,N)
!     RACAH W(I/2,J/2,M/2,L/2;K/2,N/2)
!     WDMR 21/03/2008       
      IMPLICIT NONE
      INTEGER I,J,K,L,M,N
      
      W6=C6J(I,J,K,L,M,N)
      IF (BTEST((I+J+M+L)/2,0)) W6=-W6
      
      RETURN
      END FUNCTION W6
      
            
      REAL*8 FUNCTION C6J(I,J,K,L,M,N)
!     6J - COEFFICIENT
!     ( I/2  J/2  K/2 )
!     ( L/2  M/2  N/2 ) 
!     [JB 65] (22.1.4)
	  IMPLICIT NONE
      INTEGER I,J,K,L,M,N,I1,I2,I3,I4,I5,IZ,JZ,KZ
 	  REAL*8 T,DZ
	  C6J=0.D0
	  IF(TRI(I,J,K)*TRI(I,M,N)*TRI(J,L,N)*TRI(K,L,M).EQ.0) RETURN
	  I1=(I+J+K)/2
      I2=(I+M+N)/2
      I3=(J+L+N)/2
      I4=(K+L+M)/2
      I5=(I+K+L+N)/2
	  T=DSQRT(DFLOAT((I1+1)*(I4+1))/DFLOAT((I2+1)*(I3+1))*&
         BIN(I2,I)*BIN(I,I2-N)*BIN(I4,L)/&
         BIN(I3,N)*BIN(L,I4-K)*BIN(I1,K)*BIN(K,I1-J)/&
         BIN(N,I3-L))/DFLOAT(I2-N+1)
	  JZ=MAX0(0,(I+L-J-M)/2)
	  DZ=1.D0
	  IF((JZ+I5)/2*2.NE.(JZ+I5)) DZ=-1.D0
      KZ=MIN0(I2-M,I3-J,I5-K)
	  DO IZ=JZ,KZ
	      C6J=C6J+DZ*T*BIN(I2-M,IZ)/BIN(I2,I5-K-IZ)*&
              BIN(I2-I,I3-J-IZ)/BIN(I4-I3+J+IZ+1,I2-N+1)
          DZ=-DZ
      END DO
	  RETURN
      END FUNCTION C6J

      REAL*8 FUNCTION C9J(J1,J2,J3,L1,L2,L3,K1,K2,K3)
!     9J COEFICIENT
!     (J1/2 J2/2 J3/2)
!     (L1/2 L2/2 L3/2)
!  	  (K1/2 K2/2 K3/2)  	 
!     [JB 65] (24.33)
      IMPLICIT NONE
      INTEGER J1,J2,J3,L1,L2,L3,K1,K2,K3,I,J,K,L
      C9J=0.D0
      L=TRI(J1,J2,J3)*TRI(L1,L2,L3)*TRI(K1,K2,K3)*&
          TRI(J1,L1,K1)*TRI(J2,L2,K2)*TRI(J3,L3,K3)
      IF(L.EQ.0) RETURN
      J=MAX0(IABS(J1-K3),IABS(J2-L3),IABS(L1-K2))
      K=MIN0(J1+K3,J2+L3,L1+K2)
      DO I=J,K,2
	      C9J=C9J+DFLOAT(I+1)*C6J(J1,J2,J3,L3,K3,I)*&
              C6J(L1,L2,L3,J2,I,K2)*C6J(K1,K2,K3,I,J1,L1)
      END DO
	  IF(J/2*2.NE.J) C9J=-C9J
 	  RETURN
      END FUNCTION C9J
      
	  REAL*8 FUNCTION KL0(I,J,K)
!  	  KLEBS-GORDAN COEFFICIENT WITH ZERO PROJECTIONS OF MOMENTA
!	  (I, J, K)
!	  (0, 0, 0)  
!	   I,J,K - MOMENTA = INTEGER NUMBERS
!	   [JB,65] (15.10)
      IMPLICIT NONE
	  REAL*8 T
      INTEGER I,J,K,L,M
      KL0=0.D0
      IF(TRI(I,J,K).EQ.0) RETURN
 	  L=(I+J+K)/2
	  M=L-K
	  T=1.D0
	  IF(M/2*2.NE.M) T=-1.D0
      KL0=T*BIN(K,L-J)*BIN(L,K)/&
           DSQRT(BIN(2*K,2*(L-J))*BIN(2*L+1,2*K+1))
      RETURN
      END FUNCTION KL0

   	  SUBROUTINE TRINOM
!	  THE ARRAY OF TRINOMIAL COEFFICIENTS
!  	  TIN(I,J,K)=I!!/J!!/K!!
	  IMPLICIT NONE
	  INTEGER I,J,K,M,N
	  TIN(0,0,0)=1.D0
	  TIN(1,1,1)=1.D0
	  DO I=2,99
	      M=I-I/2*2
	      TIN(I,I,M)=1.D0
	      TIN(I,M,I)=1.D0
	      N=M+2
	      DO J=I,N,-2
	          DO K=N,J,2
	              TIN(I,J,K)=TIN(I,J,K-2)/DFLOAT(K)
	              TIN(I,K,J)=TIN(I,J,K)
	          END DO
	          TIN(I,J-2,M)=TIN(I,J,M)*DFLOAT(J)
	          TIN(I,M,J-2)=TIN(I,J-2,M)
	      END DO
	  END DO
	  RETURN
	  END SUBROUTINE TRINOM
	
	  REAL*8 FUNCTION G(E1,L1,EA,LA,EB,LB)
	  IMPLICIT NONE
	  INTEGER E1,L1,EA,LA,EB,LB
	  G=KL0(LA,LB,L1)*DSQRT(DFLOAT((2*LA+1)*(2*LB+1))*&
        TIN(E1-L1,EA-LA,EB-LB)*TIN(E1+L1+1,EA+LA+1,EB+LB+1))
	  RETURN
	  END FUNCTION G

	  REAL*8 FUNCTION TMB(EE,LL,ER,LR,E1,L1,E2,L2,LM,D)
!     	   TALMI-MOSHINSKY BRACKET
!	       (EE,LL;ER,LR:LM/E1,L1;E2,L2:LM)D
	  IMPLICIT NONE
	  REAL*8 S,D,T
	  INTEGER EE,LL,ER,LR,E1,L1,E2,L2,LM
	  INTEGER L,M,ED,LD,EB,LB,EC,LC,EA,LA
 	  TMB=0.D0
	  IF(EE+ER.NE.E1+E2) RETURN
	  IF(TRI(2*LL,2*LR,2*LM)*TRI(2*L1,2*L2,2*LM).EQ.0) RETURN
	  T=DSQRT((D**(E1-ER))/((1.D0+D)**(E1+E2)))
	  M=MIN0(ER,E2)
	  S=1.D0
	  DO 1 ED=0,M
	      EB=ER-ED
	      EC=E2-ED
	      EA=E1-ER+ED
	      DO 2 LD=ED,0,-2 
	          DO 3 LB=EB,0,-2
	               IF(TRI(LD,LB,LR).EQ.0) GO TO 3
	               DO 4 LC=EC,0,-2
	                   IF(TRI(LD,LC,L2).EQ.0) GO TO 4
	                   DO 5 LA=EA,0,-2
	                        IF((TRI(LA,LB,L1).EQ.0).OR.(TRI(LA,LL,LC).EQ.0)) GO TO 5
	                        TMB=TMB+S*T*&
                                C9J(2*LA,2*LB,2*L1,2*LC,2*LD,2*L2,2*LL,2*LR,2*LM)*&
                                G(E1,L1,EA,LA,EB,LB)*G(E2,L2,EC,LC,ED,LD)*&
                                G(EE,LL,EA,LA,EC,LC)*G(ER,LR,EB,LB,ED,LD)
5 	                   CONTINUE
4 	               CONTINUE
3 	          CONTINUE
2 	      CONTINUE
	      S=S*(-D)
1 	  CONTINUE
	  RETURN
	  END FUNCTION TMB
	  
	  SUBROUTINE ORTNTMB(LAST,skirtm)
      IMPLICIT NONE 
      REAL*8 D,skirt,skirtm,atsak,check,tarp1,tarp2
      INTEGER LAST,LMB,EE,EE_,LL,LL_,e,e_,l,l_,l1,l2,e1,e2
      INTEGER*4 skaitl
      skaitl=0
      D=1.D0
      DO LMB=0,LAST
           DO EE=0,LAST
               LL=EE
               DO WHILE(LL.GE.0)
                    DO EE_=0,LAST
                       LL_=EE_
                       DO WHILE(LL_.GE.0)
                           DO e=0,LAST
                                l=e
                                DO WHILE(l.GE.0)
                                    IF(IABS(LL-l).GT.LMB) GOTO 42
                                    IF(LL+l.LT.LMB) GOTO 42
                                    DO e_=0,LAST
                                        l_=e_
                                        DO WHILE(l_.GE.0)
                                            IF(IABS(LL_-l_).GT.LMB) GOTO 43
                                            IF((LL_+l_).LT.LMB) GOTO 43
                                            atsak=0.D0
                                            DO e1=0,EE+e
                                                l1=e1
                                                DO WHILE(l1.GE.0)
                                                    DO e2=0,EE+e
                                                    IF((EE+e).NE.(e1+e2)) CYCLE
                                                    IF((EE_+e_).NE.(e1+e2)) CYCLE
                                                    l2=e2
                                                    DO WHILE(l2.GE.0)
                                                        IF(IABS(l1-l2).GT.LMB) GOTO 41
                                                        IF(l1+l2.LT.LMB) GOTO 41
                                                        tarp1=TMB(EE,LL,e,l,e1,l1,e2,l2,LMB,D)
                                                        tarp2=TMB(e1,l1,e2,l2,EE_,LL_,e_,l_,LMB,D)
                                                        atsak=atsak+tarp1*tarp2
                                                        skaitl=skaitl+2
41                                                      l2=l2-2
                                                    END DO
                                                END DO
                                            l1=l1-2
                                         END DO
                                     END DO
                                     IF((EE.EQ.EE_).AND.(LL.EQ.LL_).AND.&
                                            (e.EQ.e_).AND.(l.EQ.l_)) THEN
                                           check=1.0D0
                                     ELSE
                                           check=0.0D0
                                     END IF
                                     skirt=DABS(atsak-check)
                                     IF(skirtm<skirt) skirtm=skirt
43                                   l_=l_-2
                                  END DO
                               END DO
42                             l=l-2
                           END DO
                       END DO
                       LL_=LL_-2
                   END DO
               END DO
               LL=LL-2
             END DO
        END DO
      END DO
      write(6,4) skaitl
4     FORMAT(' TMB viso =',I10)
      RETURN
      END SUBROUTINE ORTNTMB
      
end module TLMB