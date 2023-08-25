module ho_shells
use Shells
use InputNuShell

implicit none
integer,parameter:: max_N=10
integer,parameter:: max_shells=66
integer,dimension(max_shells):: hcN,hn,hl,hs2,hj2,ht2
integer,dimension(:),allocatable:: maphi,mapih
contains

    subroutine init_shells()
    implicit none
    integer:: bN,n,l,j2,shells
    shells=0
    do bN=0,max_N
    do l=bN,0,-2
    do j2=2*l+1,abs(2*l-1),-2
    shells=shells+1
    end do
    end do
    end do
    if (shells>max_shells) then
        print *,' Increase max_shells to  ',shells
        stop
    end if    
    shells=0
    do bN=0,max_N
    do l=bN,0,-2
    do j2=2*l+1,abs(2*l-1),-2
    shells=shells+1
    hcN(shells)=bN
    n=(bN-l)/2
    hn(shells)=n
    hl(shells)=l
    hj2(shells)=j2
    hs2(shells)=1    
    ht2(shells)=1
    end do
    end do
    end do
    end subroutine  init_shells
    
    subroutine map_shells()
    implicit none
    integer:: h,i
    integer:: bN,n,l,j2,s2,t2
    call init_shells()
    allocate(maphi(max_shells))
    allocate(mapih(no_shells))
    mapih=0
    maphi=0
    do i=1,no_shells
    bN=cN(i)
    n=sn(i)
    l=sl(i)
    j2=sj2(i)
    do h=1,max_shells
    if (bN==hcN(h).and.n==hn(h).and.l==hl(h).and.j2==hj2(h)) then
        maphi(h)=i
        mapih(i)=h
    end if    
    end do
    end do
    end subroutine map_shells
       
end module ho_shells        

module harmonics
use ClebschGordan
use ho_shells

implicit none
complex(kind=16),parameter:: si=(0.0d0,1.0d0)
complex(kind=16),parameter:: r1=(1.0d0,0.0d0)

contains

    function siL(L)
    implicit none
    integer,intent(in)::L
    complex(kind=16):: siL
    if (mod(L,4)==0) siL=r1
    if (mod(L,4)==1) siL=si
    if (mod(L,4)==2) siL=-r1
    if (mod(L,4)==3) siL=-si
    end function siL
    
    function gradLrl(L,i,j,name)
    implicit none
    integer,intent(in)::L,i,j
    complex(kind=16):: gradLrl
    character(len=11),intent(inout)::name 
    real(kind=8):: fact
    name='grad(rLCL)'
    fact=sqrt(real(L*(2*L-1),8))
    gradLrl=-siL(L)*fact*cleb(2*hl(i),2*L,2*hl(j),0,0,0)*radial(hn(i),hl(i),hn(j),hl(j),L)
    if (i>j) gradLrl=conjg(gradLrl)
    end function gradlrl 
        
    function cCLrl(L,i,j,name)
    implicit none
    integer,intent(in)::L,i,j
    complex(kind=16):: cCLrl 
    character(len=11),intent(inout)::name 
    name=' rLCL'
    cCLrl=r1*cleb(2*hl(i),2*L,2*hl(j),0,0,0)*radial(hn(i),hl(i),hn(j),hl(j),L)
    end function cCLrl 
        
    function gradL(L,i,j,name)
    implicit none
    integer,intent(in)::L,i,j
    complex(kind=16):: gradL
    real(kind=8):: fact
    character(len=11),intent(inout)::name 
    name=' gradL'
    fact=sqrt(real(L*(2*L-1),8))
    gradL=siL(L)*fact*cleb(2*hl(i),2*L,2*hl(j),0,0,0)
    if (i>j) gradL=conjg(gradL)
    end function gradL
        
    function cCL(L,i,j,name)
    implicit none
    integer,intent(in)::L,i,j
    character(len=11),intent(inout)::name 
    complex(kind=16):: cCL
    name=' CL'
    cCL=r1*cleb(2*hl(i),2*L,2*hl(j),0,0,0)
    end function cCL
    
end module harmonics    
                        
module ASTop

use InputNushell
use ClebschGordan
use ho_shells
use harmonics
implicit none
integer,parameter::slen=80
type WRop
    character(len=slen):: id='none'
    integer:: L=-1
    integer:: S=-1
    integer:: J=-1
    integer:: T=-1
    real(kind=8):: N=0.0d0
    real(kind=8):: NDZ=0.0d0
    complex(kind=16),dimension(1:max_shells,1:max_shells):: el=0.0d0
end type WRop

type DZop
    character(len=slen)::  id='none'
    integer:: L=-1
    integer:: S=-1
    integer:: J=-1
    integer:: T=-1
    real(kind=8):: N=0.0d0
    real(kind=8):: NDZ=0.0d0
    real(kind=8),dimension(:,:),allocatable:: el
end type DZop

integer,private:: upper_N=100
integer,private:: lower_N=0
logical,private:: block=.false.

interface assignment (=)
    module procedure DZeWR
    module procedure WReDZ
    module procedure arrayeDZ
    module procedure DZearray
end interface

interface operator(+)
    module procedure add_op
    module procedure add_op_DZ
end interface    

interface operator(.x.)
    module procedure ctimes_op
end interface    

interface operator(.ix.)
    module procedure itimes_op
end interface    

interface operator(.cjg.)
    module procedure conj
end interface

interface operator(.v.)
    module procedure op_divc
end interface   

interface operator(.vi.)
    module procedure op_divi
end interface   

interface operator(-)
    module procedure sub_op
end interface    

interface operator(.dot.)
    module procedure dot_op
    module procedure dot_op_WRDZ
    module procedure dot_op_WRarray
    module procedure dot_op_arrayWR
    module procedure dot_op_DZWR
    module procedure dot_op_DZ
end interface    

contains
    
    subroutine set_bounds(low,high,blk)
    implicit none
    integer,intent(in):: low, high
    logical,intent(in):: blk
    block=blk
    if (low>=0) lower_N=low
    if (high>0) upper_N=high
    end subroutine set_bounds
   
    subroutine get_bounds(low,high,blk)
    implicit none
    integer,intent(out):: low, high
    logical,intent(out):: blk
    blk=block
    low=lower_N
    high=upper_N
    end subroutine get_bounds
   
    subroutine init(op)
    implicit none
    type(WRop),intent(inout)::op
    op%el=0.0d0
    op%L=-1
    op%S=-1
    op%J=-1
    op%T=-1
    op%N=0.0d0
    op%NDZ=0.0d0
    op%id='none'
    end subroutine init

    subroutine set(op,L,S,J,T,id,N,NDZ)
    implicit none
    type(WRop),intent(inout)::op
    integer,intent(in):: L,S,J,T
    real(kind=8),intent(in)::N,NDZ
    character(len=slen),intent(in):: id
    op%L=L
    op%S=S
    op%J=J
    op%T=T
    op%N=N
    op%NDZ=NDZ
    op%id=id
    end subroutine set
    
    subroutine get(op,L,S,J,T,id,N,NDZ)
    implicit none
    type(WRop),intent(in)::op
    integer,intent(inout):: L,S,J,T
    real(kind=8),intent(inout)::N,NDZ
    character(len=slen),intent(inout):: id
    L=op%L
    S=op%S
    J=op%J
    T=op%T
    N=op%N
    NDZ=op%NDZ
    id=op%id
    end subroutine get
    
    subroutine Ren(op)
    implicit none
    type(WRop),intent(inout)::op
    character(len=8):: chN2
    write(chN2,'(f8.4)') sqrt(op%N)
    op%N=1.0d0
    op%id='['//trim(adjustl(op%id))//'/'//trim(adjustl(chN2))//']'
    end subroutine Ren
    
    subroutine RenDZ(op)
    implicit none
    type(WRop),intent(inout)::op
    character(len=8):: chNDZ2
    write(chNDZ2,'(f8.4)') sqrt(op%NDZ)
    op%NDZ=1.0d0
    op%id='['//trim(adjustl(op%id))//'/'//trim(adjustl(chNDZ2))//']@DZ'
    end subroutine RenDZ
    
     
    function ad(i,j,name)
    implicit none
    complex(kind=16)::ad
    integer,intent(in):: i,j
    character(len=11),intent(inout):: name
    integer:: N,l
    name=' A*'
    l=hl(i)
    if (hl(j)==l+1.or.hl(j)==l-1) then
    N=2*hn(i)+l
    if (2*hn(j)+hl(j)/=N+1) then
        ad=(0.0d0,0.0d0)
        return
    end if
    ad=(0.0d0,0.0d0)
    if (hl(j)==l+1) ad=si*sqrt(real((l+1)*(N+l+3),8))    
    if (hl(j)==l-1) ad=si*sqrt(real(l*(N-l+2),8))
    return
    else
    ad=(0.0d0,0.0d0)
    return
    end if
    end function ad
    
    function a(i,j,name)
    implicit none
    complex(kind=16)::a
    integer,intent(in):: j,i
    character(len=11),intent(inout):: name
    integer:: N,l
    name=' A~'
    l=hl(j)
    if (hl(i)+1==l.or.hl(i)-1==l) then
    N=2*hn(j)+l
    if (2*hn(i)+hl(i)/=N+1) then
        a=(0.0d0,0.0d0)
        return
    end if
    a=(0.0d0,0.0d0)
    if (hl(i)==l+1) a=-si*sqrt(real((l+1)*(N+l+3),8))    
    if (hl(i)==l-1) a=-si*sqrt(real(l*(N-l+2),8))
    return
    else
    a=(0.0d0,0.0d0)
    return
    end if
    end function a
    
    function bHL(L,harm)
    implicit none
    type(WRop)::bHL
    integer,intent(in):: L
    character(len=11):: nameh
    character(len=2):: chL2
    integer:: i,j
    interface
        function harm(L,l1,l2,name)
        implicit none
        complex(kind=16):: harm
        integer,intent(in):: L,l1,l2
        character(len=11),intent(inout)::name
        end function harm
    end interface
    do i=1,max_shells
    do j=1,max_shells
    bHL%el(i,j)=sqrt(real(2*hl(i)+1,8))*harm(L/2,i,j,nameh)
    end do
    end do
    write(chL2,'(i2)') L/2
    bHL%L=L
    bHL%id=trim(adjustl(nameh))//'^'//trim(adjustl(chL2))
    call ncN0(bHL%N,bHL)
    call nDZcN(bHL%NDZ,bHL)
    end function bHL

    function bAd()
    implicit none
    type(WRop)::bAd
    integer:: i,j
    character(len=11)::name
    
    do i=1,max_shells
    do j=1,max_shells
    bAd%el(i,j)=ad(i,j,name)
    end do
    end do
    bAd%L=2
    bAd%id=trim(adjustl(name))//'^1'
    call ncN0(bAd%N,bAd)
    call nDZcN(bAd%NDZ,bAd)
    end function bAd

    function bL(no)
    implicit none
    type(WRop):: bL
    integer,intent(in):: no
    real(kind=8):: clb,fact    
    integer:: i,j,m
    character(len=2):: chno
    bL%el=(0.0d0,0.0d0)
    do i=1,max_shells
    do j=i,max_shells
    if (-2*hl(i)+no<2*hl(i)) then
    clb=cleb(2*hl(i),no,2*hl(i),-2*hl(i),no,no-2*hl(i))
    fact=Ld(no/2,hl(i),-hl(i))
    if (fact*clb/=0.0d0) then
    bL%el(i,j)=siL(no/2)*fact/clb
    end if
    end if
    end do
    end do
    do j=1,max_shells
    do i=j+1,max_shells
    bL%el(i,j)=conjg(bL%el(j,i))
    end do
    end do
    write(chno,'(i2)') no/2
    bL%id='lambda^'//trim(adjustl(chno))
    call ncN0(bL%N,bL)
    call nDZcN(bL%NDZ,bL)
    bL%L=no/2
    end function bL
    
    function bA()
    implicit none
    type(WRop)::bA
    integer:: i,j
    character(len=11)::name
    
    do i=1,max_shells
    do j=1,max_shells
    bA%el(i,j)=a(i,j,name)
    end do
    end do
    bA%L=2
    bA%id=trim(adjustl(name))//'^1'
    call ncN0(bA%N,bA)
    call nDZcN(bA%NDZ,bA)
    end function bA
    
    function bM()
    implicit none
    type(WRop):: bM
    bM%id=' I'
    bM%el=(1.0d0,0.0d0)
    bM%S=0
    bM%N=1.0d0
    bM%NDZ=1.0d0
    end function bM
    
    function bT()
    implicit none
    type(WRop):: bT
    bT%el=cmplx(-sqrt(6.0d0))
    bT%T=2
    bT%id=' Tau'
    bT%N=1.0d0
    bT%NDZ=1.0d0
    end function bT
    
    function conj(op)
    implicit none
    type(WRop)::conj
    type(WRop),intent(in)::op
    integer::i,j
    do i=1,max_shells
    do j=1,max_shells
    conj%el(i,j)=conjg(op%el(i,j))
    end do
    end do
    conj%id='['//trim(adjustl(op%id))//']*'
    conj%L=op%L
    conj%S=op%S
    conj%J=op%J
    conj%T=op%T
    conj%NDZ=op%NDZ
    conj%N=op%N
    end function conj    
        
    function anti(op)
    implicit none
    type(WRop)::anti
    type(WRop),intent(in)::op
    integer::i,j
    do i=1,max_shells
    do j=1,max_shells
    if (i<=j) then
    anti%el(i,j)=op%el(i,j)
    else if (i>j) then
    anti%el(i,j)=-op%el(i,j)
    end if
    end do
    end do
    anti%id='&['//trim(adjustl(op%id))//']&'
    anti%L=op%L
    anti%S=op%S
    anti%J=op%J
    anti%T=op%T
    anti%NDZ=op%NDZ
    anti%N=op%N
    end function anti    
        
    subroutine writeopDZ(dev,op)
    implicit none
    integer,intent(in)::dev
    type(DZop),intent(in)::op
    integer::i,j
    write(dev,*) op%id
    write(dev,'(a14,4i4)') '  2L,2S,2J,2T ',op%L,op%S,op%J,op%T
    write(dev,'(a13,2f15.3)') '  N^2, NDZ^2 ',op%N,op%NDZ
    write(dev,*)
    do i=1,no_shells
    do j=1,no_shells
    if (abs(op%el(i,j))>1.0d-7) write(dev,'(2i4,f15.9)') i,j,op%el(i,j)
    end do
    end do
    write(dev,*)
    end subroutine writeopDZ    
    
    subroutine writeDZop(dev,op)
    implicit none
    integer,intent(in)::dev
    type(WRop),intent(in)::op
    integer::i,j
    write(dev,*) op%id
    write(dev,'(a14,4i4)') '  2L,2S,2J,2T ',op%L,op%S,op%J,op%T
    write(dev,'(a13,2f15.3)') '  N^2, NDZ^2 ',op%N,op%NDZ
    write(dev,*)
    do i=1,no_shells
    do j=1,no_shells
    if (abs(op%el(mapih(i),mapih(j)))>1.0d-7) write(dev,'(2i4,2f15.9)') i,j,op%el(mapih(i),mapih(j))
    end do
    end do
    write(dev,*)
    end subroutine writeDZop    
    
    function bCPLrl(LC,LP,L)
    implicit none
    type(WRop)::bCPLrl
    integer,intent(in):: LC,LP,L
    character(len=2):: chLC2,chLP2,chL2
    real(kind=8):: rad,clb,fact
    integer:: i,j
    fact=sqrt(real(LP*(2*LP-1),8))
    do i=1,max_shells
    do j=1,max_shells
    rad=radial(hn(i),hl(i),hn(j),hl(j),(LC+LP)/2)
    clb=cleb(2*hl(i),L,2*hl(j),0,0,0)
    bCPLrl%el(i,j)=siL(LP/2)*rad*fact*clb
    if (i>j) bCPLrl%el(i,j)=conjg(bCPLrl%el(i,j))
    end do
    end do
    write(chLC2,'(i2)') LC/2
    write(chLP2,'(i2)') LP/2
    write(chL2,'(i2)') L/2
    bCPLrl%id='[cCLrl^'//trim(adjustl(chLC2))//',gradlrl^'//trim(adjustl(chLP2))//&
               ']^'//trim(adjustl(chL2))
    bCPLrl%L=L               
    call ncN0(bCPLrl%N,bCPLrl)
    call nDZcN(bCPLrl%NDZ,bCPLrl)
    end function bCPLrl
    
    function bCPL(LC,LP,L)
    implicit none
    type(WRop)::bCPL
    integer,intent(in):: LC,LP,L
    character(len=2):: chLC2,chLP2,chL2
    real(kind=8):: clb,fact
    integer:: i,j
    fact=sqrt(real(LP*(2*LP-1),8))
    do i=1,max_shells
    do j=1,max_shells
    clb=cleb(2*hl(i),L,2*hl(j),0,0,0)
    bCPL%el(i,j)=siL(LP/2)*fact*clb
    if (i>j) bCPL%el(i,j)=conjg(bCPL%el(i,j))
    end do
    end do
    write(chLC2,'(i2)') LC/2
    write(chLP2,'(i2)') LP/2
    write(chL2,'(i2)') L/2
    bCPL%id='[cCL^'//trim(adjustl(chLC2))//',gradL^'//trim(adjustl(chLP2))//&
               ']^'//trim(adjustl(chL2))
    bCPL%L=L               
    call ncN0(bCPL%N,bCPL)
    call nDZcN(bCPL%NDZ,bCPL)
    end function bCPL
    
    function bS()
    implicit none
    type(WRop)::bS
    bs%el=cmplx(-sqrt(6.0d0))    
    bs%S=2
    bs%id='Sigma'
    bs%N=1.0d0
    bs%NDZ=1.0d0    
    end function bS
  
   
    subroutine nDZcN(cN,op)
    implicit none
    type(WRop),intent(inout)::op
    real(kind=8),intent(inout):: cN
    cN=sum(op%el(mapih(:),mapih(:))*conjg(op%el(mapih(:),mapih(:))))
    if (cN/=0.0d0) op%el=op%el/sqrt(cN)
    end subroutine nDZcN
    
    function n(op)
    implicit none
    type(WRop):: n
    type(WRop),intent(in):: op
    real(kind=8):: norm2
    norm2=sum(op%el*conjg(op%el))
    if (norm2/=0.0d0) n%el=op%el/sqrt(norm2)
    n%id=op%id
    n%L=op%L
    n%S=op%S
    n%J=op%J
    n%T=op%T
    n%N=norm2
    end function n
    
    function nDZ(op)
    implicit none
    type(WRop):: nDZ
    type(WRop),intent(in):: op
    real(kind=8):: norm2
    nDZ%el=0.0d0
    norm2=sum(op%el(mapih(:),mapih(:))*conjg(op%el(mapih(:),mapih(:))))
    if (norm2/=0.0d0) nDZ%el=op%el/sqrt(norm2)
    nDZ%id=op%id
    nDZ%L=op%L
    nDZ%S=op%S
    nDZ%J=op%J
    nDZ%T=op%T
    nDZ%NDZ=norm2
    end function nDZ
    
    subroutine ncN(cN,op)
    implicit none
    type(WRop),intent(inout)::op
    real(kind=8),intent(inout):: cN
    cN=abs(sum(op%el*conjg(op%el)))
    if (cN/=0.0d0) op%el=op%el/sqrt(cN)
    end subroutine ncN
    
    subroutine ncN0(cN,op)
    implicit none
    type(WRop),intent(inout)::op
    real(kind=8),intent(inout):: cN
    cN=sum(op%el*conjg(op%el))
    end subroutine ncN0
    
    function jj(lop,sop,J2)
    implicit none
    type(WRop):: jj
    type(WRop),intent(in):: lop,sop
    integer,intent(in):: J2 
    integer:: i,j
    character(len=2):: chJ2
    real(kind=8):: fact,fact1
    complex(kind=16):: prod
    if (sop%S<0.or.lop%L<0) stop ' Error WR01 : Unassigned operator'
    jj%el=0.0d0
    write(chJ2,'(i2)') J2/2
    fact=sqrt(real((sop%S+1)*(lop%L+1),8))
    do i=1,max_shells
    do j=1,max_shells
    prod=lop%el(i,j)*sop%el(i,j)
    if (prod/=0.0d0) then
    fact1=fact*sqrt(real((hj2(i)+1)*(hj2(j)+1),8))
    jj%el(i,j)=fact1*ninej(2*hl(i),1,hj2(i),&
                   2*hl(j),1,hj2(j),lop%L,sop%S,J2)*prod
    end if
    end do
    end do
    jj%id='['//trim(adjustl(lop%id))//','//trim(adjustl(sop%id))//']^'&
                    //trim(adjustl(chJ2))
    jj%L=lop%L
    jj%S=sop%S
    jj%J=J2
    call ncN0(jj%N,jj)
    call nDZcN(jj%NDZ,jj)
    end function jj
    
    function add_op(op1,op2)
    implicit  none
    type(WRop)::add_op
    type(WRop),intent(in):: op1,op2
    real(kind=8)::x1,x2
    if (op1%L/=op2%L) stop ' Error WR02 : different L'
    if (op1%S/=op2%S) stop ' Error WR02 : different S'
    if (op1%J/=op2%J) stop ' Error WR02 : different J'
    if (op1%T/=op2%T) stop ' Error WR02 : different T'
    x1=sqrt(op1%N)
    x2=sqrt(op2%N)
    add_op%el=x1*op1%el+x2*op2%el
    add_op%id='['//trim(adjustl(op1%id))//'+'//trim(adjustl(op2%id))//']'
    add_op%L=op1%L
    add_op%S=op1%S
    add_op%J=op1%S
    call ncN0(add_op%N,add_op)
    call nDZcN(add_op%NDZ,add_op)
    
    end function add_op
    
    function add_op_DZ(op1,op2)
    implicit  none
    type(DZop)::add_op_DZ
    type(DZop),intent(in):: op1,op2
    type(WRop):: wr1,wr2
    if (op1%L/=op2%L) stop ' Error WR02 : different L'
    if (op1%S/=op2%S) stop ' Error WR02 : different S'
    if (op1%J/=op2%J) stop ' Error WR02 : different J'
    if (op1%T/=op2%T) stop ' Error WR02 : different T'
    wr1=op1
    wr2=op2
    wr1%N=op1%NDZ
    wr2%N=op2%NDZ
    add_op_DZ=wr1+wr2
    add_op_DZ%id='['//trim(adjustl(op1%id))//'+'//trim(adjustl(op2%id))//']@DZ'
    
    end function add_op_DZ

    function sub_op(op1,op2)
    implicit  none
    type(WRop)::sub_op
    type(WRop),intent(in):: op1,op2
    if (op1%L/=op2%L) stop ' Error WR02 : different L'
    if (op1%S/=op2%S) stop ' Error WR02 : different S'
    if (op1%J/=op2%J) stop ' Error WR02 : different J'
    if (op1%T/=op2%T) stop ' Error WR02 : different T'
    sub_op%el=sqrt(op1%N)*op1%el-sqrt(op2%N)*op2%el
    sub_op%id='['//trim(adjustl(op1%id))//'-'//trim(adjustl(op2%id))//']'
    sub_op%L=op1%L
    sub_op%S=op1%S
    sub_op%J=op1%S
    call ncN0(sub_op%N,sub_op)
    call nDZcN(sub_op%NDZ,sub_op)
    
    end function sub_op
    
    function sub_op_DZ(op1,op2)
    implicit  none
    type(DZop)::sub_op_DZ
    type(DZop),intent(in):: op1,op2
    type(WRop):: wr1,wr2
    if (op1%L/=op2%L) stop ' Error WR02 : different L'
    if (op1%S/=op2%S) stop ' Error WR02 : different S'
    if (op1%J/=op2%J) stop ' Error WR02 : different J'
    if (op1%T/=op2%T) stop ' Error WR02 : different T'
    wr1=op1
    wr2=op2
    wr1%N=op1%NDZ
    wr2%N=op2%NDZ
    sub_op_DZ=wr1-wr2
    sub_op_DZ%id='['//trim(adjustl(op1%id))//'-'//trim(adjustl(op1%id))//']@DZ'
    
    end function sub_op_DZ
    
    function ctimes_op(c,op)
    implicit  none
    type(WRop)::ctimes_op
    type(WRop),intent(in):: op
    real(kind=8),intent(in):: c
    character(len=8):: chc2
    write(chc2,'(f8.4)') c
    
    ctimes_op=op
    ctimes_op%N=op%N*c*c
    ctimes_op%NDZ=op%NDZ*c*c
    ctimes_op%id='['//trim(adjustl(chc2))//trim(adjustl(op%id))//']'
    
    end function ctimes_op

    function itimes_op(op)
    implicit  none
    type(WRop)::itimes_op
    type(WRop),intent(in):: op
    character(len=3):: chc2='(i)'
    
    itimes_op=op
    itimes_op%el=cmplx(0.0d0,1.0d0)*op%el
    itimes_op%id='['//trim(adjustl(chc2))//trim(adjustl(op%id))//']'
    
    end function itimes_op

    
    function op_divc(op,c)
    implicit  none
    type(WRop)::op_divc
    type(WRop),intent(in):: op
    real(kind=8),intent(in):: c
    character(len=8):: chc2
    write(chc2,'(f8.4)') c
    
    op_divc=op
    op_divc%N=op%N/c/c
    op_divc%NDZ=op%NDZ/c/c
    op_divc%id='['//trim(adjustl(op%id))//'/'//trim(adjustl(chc2))//']'
    
    end function op_divc

    function op_divi(op)
    implicit  none
    type(WRop)::op_divi
    type(WRop),intent(in):: op
    character(len=3):: chc2='(i)'
    
    op_divi=op
    op_divi%el=-cmplx(0.0d0,1.d0)*op%el
    op_divi%id='['//trim(adjustl(op%id))//'/'//trim(adjustl(chc2))//']'
    
    end function op_divi

    
    function dot_op(op1,op2)
    implicit  none
    real(kind=8)::dot_op
    type(WRop),intent(in):: op1,op2
    
    dot_op=sum(op1%el*conjg(op2%el))
    
    end function dot_op

    function dot_op_DZ(op1,op2)
    implicit  none
    real(kind=8)::dot_op_DZ
    type(DZop),intent(in):: op1,op2
    
    dot_op_DZ=sum(op1%el*op2%el)
    
    end function dot_op_DZ

    function dot_op_WRDZ(wr,dz)
    implicit  none
    complex(kind=16)::dot_op_WRDZ
    type(WRop),intent(in):: wr
    type(DZop),intent(in):: dz
    
    dot_op_WRDZ=sum(dz%el*wr%el(mapih(:),mapih(:)))
    
    end function dot_op_WRDZ

    function dot_op_WRarray(wr,array)
    implicit  none
    complex(kind=16)::dot_op_WRarray
    real(kind=8),dimension(:),intent(in):: array
    type(WRop),intent(in):: wr
    type(DZop):: dz
    allocate(dz%el(no_shells,no_shells))
    dz=array
    dot_op_WRarray=sum(dz%el*wr%el(mapih(:),mapih(:)))
    deallocate(dz%el)
    end function dot_op_WRarray

    function dot_op_DZWR(dz,wr)
    implicit  none
    complex(kind=16)::dot_op_DZWR
    type(WRop),intent(in):: wr
    type(DZop),intent(in):: dz
    
    dot_op_DZWR=sum(dz%el*wr%el(mapih(:),mapih(:)))
    
    end function dot_op_DZWR

    function dot_op_arrayWR(array,wr)
    implicit  none
    real(kind=8),dimension(:),intent(in)::array
    complex(kind=16)::dot_op_arrayWR
    type(WRop),intent(in):: wr
    type(DZop):: dz
    allocate(dz%el(no_shells,no_shells))
    dz=array
    dot_op_arrayWR=sum(dz%el*wr%el(mapih(:),mapih(:)))
    deallocate(dz%el)
    end function dot_op_arrayWR

    subroutine DZeWR(dz,wr)
    implicit none
    type(DZop),intent(inout):: dz
    type(WRop),intent(in):: wr
    integer:: i,j
    do i=1,no_shells
    do j=1,no_shells
    dz%el(i,j)=real(wr%el(mapih(i),mapih(j)))-aimag((wr%el(mapih(i),mapih(j))))
    end do
    end do
    dz%id=wr%id
    dz%L=wr%L
    dz%S=wr%S
    dz%J=wr%J
    dz%T=wr%T
    dz%N=wr%N
    dz%NDZ=wr%NDZ
    end subroutine DZeWR

    subroutine WReDZ(wr,dz)
    implicit none
    type(WRop),intent(inout):: wr
    type(DZop),intent(in):: dz
    integer:: i,j
    do i=1,no_shells
    do j=1,no_shells
    if (i>j)wr%el(mapih(i),mapih(j))=cmplx(0.0d0,dz%el(i,j))
    if (i<=j) wr%el(mapih(i),mapih(j))=cmplx(dz%el(i,j))
    end do
    end do
    wr%id=dz%id
    wr%L=dz%L
    wr%S=dz%S
    wr%J=dz%J
    wr%T=dz%T
    wr%N=dz%N
    wr%NDZ=dz%NDZ
    end subroutine WReDZ
    
    subroutine DZearray(dz,arr)
    implicit none
    type(DZop),intent(inout):: dz
    real(kind=8),dimension(size(dz%el)),intent(in):: arr
    integer:: i,j,vid
    do i=1,no_shells
    do j=1,no_shells
    vid=(i-1)*no_shells+j
    dz%el(i,j)=arr(vid)
    end do
    end do
    end subroutine DZearray
    
    subroutine arrayeDZ(arr,dz)
    implicit none
    type(DZop),intent(in):: dz
    real(kind=8),dimension(size(dz%el)),intent(inout):: arr
    integer:: i,j,vid
    do i=1,no_shells
    do j=1,no_shells
    vid=(i-1)*no_shells+j
    arr(vid)=dz%el(i,j)
    end do
    end do
    end subroutine arrayeDZ
    
end module ASTop    
