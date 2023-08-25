module TrdChar

use ProjParameters
use LanczosParameters
use OrderParameters
use OperParameters
use InputNuShell
use MatrixParameters

implicit none

integer,parameter:: mask=z'ff'

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
    