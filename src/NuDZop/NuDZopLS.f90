! This file contains the modules specific the the generic pOper
! program for NuShell, NuShell@MSU and SunShell.
module Ordered
use Parameters
use OperParameters
use Partition

implicit none

type(spartition),dimension(:),allocatable:: order_int
real(kind=ri),dimension(:),allocatable:: order_V
integer:: n_int

end module Ordered

module SU3

use TLMB
use InputNushell
use Ordered
use OperParameters
use TriDiag

contains

subroutine SU3_v(vec,dvec,n,nd)
implicit none
real(kind=ri),dimension(:),intent(inout):: vec,dvec
real(kind=ri),intent(inout):: n,nd
integer:: ih,jh,il,ir,L
real(kind=ri):: d,m,mt,md,Ni,Nj,ml
call init()
vec=0.0d0
dvec=0.0d0
do ih=1,no_shells
do jh=1,no_shells
m=0.0d0
md=0.0d0
ml=0.0d0
Ni=2*sn(ih)+sl(ih)
Nj=2*sn(jh)+sl(jh)
mt=0.0d0
if (Ni==Nj) then
do L=min(abs(sl(jh)-1),abs(sl(ih)-1)),max(abs(sl(jh)+1),abs(sl(ih)+1))
m=0.0d0
md=0.0d0
if (sl(jh)==L+1) m=-sqrt(real((sl(jh)+1)*(Nj+sl(jh)+3),ri))
if (sl(ih)==L+1) md=-sqrt(real((sl(ih)+1)*(Ni+sl(ih)+3),ri))
if (sl(ih)==L-1) md=-sqrt(real(sl(ih)*(Ni-sl(ih)+2),ri))
if (sl(jh)==L-1) m=-sqrt(real(sl(jh)*(Ni-sl(jh)+2),ri))
mt=mt+md*m*sqrt(real(6*(2*L+1),ri))*W6(2*sl(ih),2,2*sl(jh),2,2*L,4)
end do
if (sl(ih)==sl(jh).and.ih/=jh) ml=sqrt(real(sl(ih)*(sl(ih)+1),ri))
end if
il=(ih-1)*no_shells+jh
if (mt/=0.0d0) vec(il)=vec(il)+mt*&
 sqrt(real((sj2(ih)+1)*(sj2(jh)+1)*3,ri))*c9j(2*sl(ih),1,sj2(ih),2*sl(jh),1,sj2(jh),4,0,4)
if (ml/=0.0d0) dvec(il)=dvec(il)+ ml*&
 sqrt(real((sj2(ih)+1)*(sj2(jh)+1)*3,ri))*c9j(2*sl(ih),1,sj2(ih),2*sl(jh),1,sj2(jh),2,0,2)
end do  
end do
n=dot_product(vec,vec)
nd=dot_product(dvec,dvec)
if (n/=0.0) vec=vec/sqrt(n)
if (nd/=0.0) dvec=dvec/sqrt(nd)
end subroutine SU3_v

end  module SU3

module Hcm

use TLMB
use InputNushell
use Ordered
use OperParameters
use TriDiag

contains

subroutine Hcm_v(vec,dvec)
implicit none
real(kind=ri),dimension(:),intent(inout):: vec,dvec
integer:: ih,jh,il,ir
real(kind=ri):: d,m,n,md,nd
call init()

vec=0.0d0
dvec=0.0d0
do ih=1,no_shells
do jh=1,no_shells
m=0.0d0
if (sn(ih)==sn(jh).and.sl(jh)==sl(ih)+1) m=sqrt(real((sl(ih)+1)*(2*sn(ih)+2*sl(ih)+3),ri)/3.0d0)
if (sn(ih)==sn(jh).and.sl(ih)==sl(jh)+1) m=sqrt(real((sl(jh)+1)*(2*sn(jh)+2*sl(jh)+3),ri)/3.0d0)
if (sn(jh)==sn(ih)+1.and.sl(ih)==sl(jh)+1) m=sqrt(real(sl(ih)*(sn(ih)+1),ri))
if (sn(ih)==sn(jh)+1.and.sl(jh)==sl(ih)+1) m=sqrt(real(sl(jh)*(sn(jh)+1),ri))
m=-m
il=(ih-1)*no_shells+jh
if (m/=0.0d0)  vec(il)=vec(il)+m*&
 sqrt(real((sj2(ih)+1)*(sj2(jh)+1)*3,ri))*c9j(2*sl(ih),1,sj2(ih),2*sl(jh),1,sj2(jh),2,0,2)
end do
end do
n=dot_product(vec,vec)
if (n/=0.0) vec=vec/sqrt(n)
do ih=1,no_shells
do jh=1,no_shells
il=(ih-1)*no_shells+jh
dvec(il)=-vec(il)
if (2*sn(ih)+sl(ih)>2*sn(jh)+sl(jh)) dvec(il)=-dvec(il)
end do
end do

end subroutine Hcm_v

subroutine remove_Hcm(va,e,inx,nos,hcm1,hcm2,lim)
implicit none
real(kind=ri),dimension(:),intent(inout):: e
real(kind=ri),dimension(:,:),intent(inout):: va
real(kind=ri),dimension(:),intent(in):: hcm1,hcm2
real(kind=ri),intent(in):: lim
integer,dimension(:),intent(in):: inx
integer,intent(in):: nos
integer,dimension(nos*nos):: inx2
integer:: nos2,s1,s2,s3,s4,ir,ii,ic,jj
real(kind=ri):: ovl1,ovl2
nos2=nos*nos
call indexx(nos2,-abs(e),inx2)
do ic=1,nos2
ii=inx2(ic)
ovl1=0.0d0
ovl2=0.0d0
ovl1=dot_product(hcm1(:),va(:,ii))
va(:,ii)=va(:,ii)-ovl1*hcm1
ovl1=dot_product(va(:,ii),va(:,ii))
if (abs(ovl1)<lim) then
    e(ii)=0.0d0
    va(:,ii)=hcm1(:)
    print *, ' Removing Hcm vector 1'
    go to 10
end if
ovl2=dot_product(hcm2(:),va(:,ii))
va(:,ii)=va(:,ii)-ovl2*hcm2
ovl1=dot_product(va(:,ii),va(:,ii))
if (abs(ovl1)<lim) then
    e(ii)=0.0d0
    va(:,ii)=hcm2(:)
    print *, ' Removing Hcm vector 2'
    go to 10
end if
if (ovl1/=0) va(:,ii)=va(:,ii)/sqrt(ovl1)
10  do ir=1,ic-1
    jj=inx2(ir)
    if (abs(e(jj))>=0.00000000001d0) then
    ovl2=dot_product(va(:,ii),va(:,jj))
    va(:,ii)=va(:,ii)-ovl2*va(:,jj)
    end if
end do    
ovl1=dot_product(va(:,ii),va(:,ii))
if (ovl1/=0) va(:,ii)=va(:,ii)/sqrt(ovl1)
end do

end subroutine remove_Hcm

end module Hcm


module OperSupport

use InputNuShell
use OutputNuShell
use PartOrder
use Files
use OperParameters
use TrdChar
use Ordered

implicit none

integer,dimension(max_no_shells):: ishell
integer:: ni_shells,Jmax,Tmax,noint
type(spartition):: jsJT
real(kind=ri):: V
logical:: SDI=.false.
real(kind=ri),dimension(0:1):: A
real(kind=ri):: B,C,Qq,P,R,S,LS,B2,X,Qsp,Qm,Qh,Qf
character(len=1) :: xtra
character(len=80):: first_line
character(len=7)::intfile
 contains

    subroutine setup_ishell()
    implicit none
        
        integer:: i,k,j
        logical:: lid
        ishell=0
        ni_shells=1
        ishell(1)=1
        do i=2,no_shells
            lid = .false.
            do j=1,ni_shells
                k=ishell(j)   
                if (cN(i) == cN(k) .and. sn(i) == sn(k)  .and. sl(i) == sl(k) .and.   &
                sj2(i) == sj2(k) .and. st2(i) == st2(k) .and. sp(i) == sp(k)) then
                    ! identical shells
                    if ( nucleon(i) == nucleon(k) ) then
                        print *, ' Error in input. Identical shells.'
                        stop
                    end if
                    if (.not.( nucleon(i) == 'n' .and. nucleon(k) == 'p')) then
                        print *, ' Error in input. Labelling  i,p,n'
                        print *, ' i before p before n '
                        stop
                    end if
                lid=.true.    
                exit
                end if
                
            end do
            if (.not.lid) then
                    ni_shells=ni_shells+1
                    ishell(i)=ni_shells
            else
                    ishell(i)=k
            end if
         end do
 
         if (output_control>1) then
            print *, ' No of isospin shells,',ni_shells
            print *,ishell(1:no_shells)
         end if
                  
    end subroutine setup_ishell
    
    subroutine read_interaction(inExt)
    implicit none
        
        character(len=4)::inExt
        character(len=1):: commentsv
        integer:: n,i,k,j,nn
        logical:: done=.true.,more=.false.
        type(spartition):: part
        real(kind=ri),dimension(max_no):: norm

        call set_ptdim(6)
        do i=1,max_part_tree
        ipart_order(i)%part=0
        ipart_order(i)%V=real(0,ri)
        ipart_order(i)%node=0
        ipart_order(i)%lesser=0
        ipart_order(i)%greater=0
        end do
        do i=1,max_int
            order_int(i)=0
        end do
        order_v=real(0,ri)
        call output_add_int(6)
        call input_intfile(intfile)
        open(unit=10,file=intfile//inExt,err=100)
        read(10,'(a80)') first_line
        call read_comments()   
        commentsv=comment     
        norm=real(1,ri)
        
1       read(10,*,err=200,end=201) n, spe(1:no_shells)
        if (-n > max_no .and. comment == '&') then
            print *, ' Incease max_no'
            stop
        end if
        if (n < -1 .and. comment == '&') write(6,'(a19,i2,a2)',advance='no')&
                     'Normalisations n=1,',-n,'  '
        if (n < -1 .and. comment == '&') read(5,*) norm(1:-n)
        nn=1
        if (comment == '^') then
            SDI=.true.
            read(10,*,err=200,end=201) A,B,C,P,R,B2,S,LS,X
            read(10,*,err=200,end=201) Qsp,Qm,Qq,Qh,Qf
            call reset_ipart_order()
            i=0
            go to 10
        end if
        call reset_ipart_order()
        Jmax=0
        do 
            part=0
            read(10,*,err=200,end=5) (part%shell(k),k=1,6), V
            call read_comments(more)
            if (more .and. n<-1 .and. commentsv == '&') then
                nn=nn+1
                read(10,*,err=200,end=5)
                go to 4 
            end if
            if (part%shell(5) > Jmax) Jmax=part%shell(5)
            call add_part(part,V*norm(nn),.true.)
4       end do
        
5       i=0
        do
            k=get_next_ipart(done)
            i=i+1
            if (i > noint) then
                print *, ' Increase max_int in OperParameters or'
                print *, ' SET INTERACTION=I where I >>', i
                stop
            end if
            order_int(i)=ipart_order(k)%part
            order_V(i)=ipart_order(k)%V
            if (done) exit
        end do
        open(unit=30,file=Nucleus//xtra//'tint.txt',action='write',err=300) 
        write(30,*,err=400) '! Total Interaction'       
        write(30,'(i4)',err=400) i
        do j=1,i
            if (order_V(i) /= real(0,ri)) write(30,'(4i4,2i6,f9.4)',err=400) &
                  order_int(j)%shell(1:6),order_V(j)
        end do
        
        close(30)
        
10      close(10)
        
        n_int=i
        Jmax=2*Jmax
        Tmax=2
        deallocate(ipart_order)
        return
100     print *, ' Error opening ',intfile//inExt
        stop        
200     print *, ' Error reading ',intfile//inExt
        stop        
201     print *, ' EOF reading ',intfile//inExt
        stop        
300     print *, ' Error opening ',Nucleus//xtra//'tint.txt'
        stop        
400     print *, ' Error writing ',Nucleus//xtra//'tint.txt'
        stop        
   end subroutine read_interaction 
      
end module OperSupport

module Hoper

use InputNuShell
use OperSupport
use Shells
use ShellsX
use Finds
use ClebschGordan
use Clebschg
use Partition
use OperParameters
use TrdChar
use SOrder
use TriDiag
use Ordered
use Hcm
use SU3
use XParameters
use IFPORT

implicit none

character(len=15):: cute,cutm
character(len=1):: chcm,cedit
logical:: notf,lhs,rhs,npintn,npoutn
real(kind=8):: two=2.0d0
contains

    subroutine H_op(inExt)
    implicit none
    
        character(len=*):: inExt
        integer:: m1,m2,m3,m4,tz1,tz2,tz3,tz4,J,T,j1,j2,j3,j4,Tmin,Tmax,k,i,lam,max_op,tau
        integer:: sh1,sh2,sh3,sh4,iw1,iw2,iw3,iw4,ib1,ib2,ib3,ib4,p1,p2,p3,p4,lmn,lmx,taux
        real(kind=ri):: V,sumV,sumVn,sumVp,delta12,delta34,vint,vintn,vintp,phase
        real(kind=8)::me,average,four2,limit,rcute,rcutm,rqq,rll
        integer:: no_Op=0,one=1,lmid,rmid,tqlim,lamx,nz_op,lmids,rmids,cL,cS
        integer:: sigma,sigmx,l1,l2,l3,l4
        type(spartition):: parti
        character(len=8):: partm
        character(len=8):: z3
        character(len=8),dimension(:),allocatable:: opn
        real(kind=8),dimension(:),allocatable:: meV,meo,beta,Esp
        real(kind=8),dimension(:,:,:,:),allocatable:: alpha,vector
        logical(kind=1),dimension(:,:,:,:),allocatable:: lalpha
        real(kind=8),dimension(:,:,:,:,:),allocatable:: matrix,somega,omega
        logical(kind=1),dimension(:,:,:,:,:),allocatable:: lmatrix
        logical(kind=1),dimension(:),allocatable:: sh1sh2
        real(kind=8),dimension(:),allocatable:: hcmvec,hcm2c
        integer,dimension(:),allocatable:: pl,pr,plo,pro
        integer,dimension(:,:,:,:),allocatable:: indx
        logical:: warn=.false.,done,warnp=.false.
        integer:: opd,pnnp,iout,max_lmo=0
        character(len=10):: chnoop
        logical:: btest,intgj,rem_Hcm=.false.,do_edit=.false.,do_cut=.false.,add_SU3=.false.
        integer:: min_bshells,max_bshells,min_ashells,max_ashells
        limit=hlim
        call get_command_argument(5,chcm)
        call get_command_argument(6,cedit)
        if (chcm=='h'.or.chcm=='H')rem_Hcm=.true.
        if (chcm=='s'.or.chcm=='S')add_SU3=.true.
        if (cedit=='e'.or.cedit=='E') do_edit=.true.
        if (cedit=='c'.or.cedit=='C') do_cut=.true.
        if (rem_Hcm) print *, ' Hcm coupling terms will be removed.'
        if (add_SU3) print *, ' SU3 coupling terms will be added.'
        if (do_cut) then
        call get_command_argument(7,cute)
        call get_command_argument(8,cutm)
        read(cute,*) rcute
        read(cutm,*) rcutm
        print *, ' e cut sqrt((2lambda+1)(2tau+1)(2sigma+1))*',rcute 
        print *, ' v cut',rcutm        
        else
            rcute=.0000000001
            rcutm=.0000000001
        end if
        intgj=.false.
        if (btest(sj2(1)+1,0)) intgj=.true.
        
        z3=repeat(char(0),8)
        nstring=notree
        allocate(string_order(notree))
        allocate(meo(notree))
        allocate(plo(notree))
        allocate(pro(notree))
        allocate(matrix(no_shells*no_shells,no_shells*no_shells,0:max_sj2u,0:1,0:1))
        allocate(lmatrix(no_shells*no_shells,no_shells*no_shells,0:max_sj2u,0:1,0:1))
        allocate(sh1sh2(no_shells*no_shells))
        allocate(lalpha(no_shells*no_shells,0:max_sj2u,0:1,0:1))
        allocate(alpha(no_shells*no_shells,0:max_sj2u,0:1,0:1))
        allocate(vector(no_shells*no_shells,-1:1,0:1,0:1))
        allocate(indx(no_shells*no_shells,0:max_sj2u,0:1,0:1))
        allocate(somega(no_shells*no_shells,no_shells*no_shells,0:max_sj2u,0:1,0:1))
        allocate(omega(no_shells*no_shells,no_shells*no_shells,0:max_sj2u,0:1,0:1))
        allocate(hcmvec(no_shells*no_shells))
        allocate(hcm2c(no_shells*no_shells))
        if (rem_Hcm)  call Hcm_v(hcmvec,hcm2c)
        if (add_SU3)  call SU3_v(hcmvec,hcm2c,rqq,rll)
        matrix=real(0,8)
        somega=real(0,8)
        alpha=real(0,8)
        allocate(beta((no_shells*no_shells)))
        allocate(Esp(no_shells))
        Esp(1:no_shells)=spe(1:no_shells)
        open(unit=20,file=Nucleus//xtra//inExt, form='unformatted', action='write',err=100)
        open(unit=21,file=Nucleus//xtra//'.misV', form='formatted', action='write',err=200)
        open(unit=22,file=Nucleus//xtra//'oph.'//'txt', &
                                     form='formatted', action='write')
       
        !loop over particle 1 shells
        typea='p'
        typeb='n'
        if (npintn) then
        min_ashells=min_shellsa()
        max_ashells=max_shellsa()
        min_bshells=min_shellsb()
        max_bshells=max_shellsb()
        else
        min_ashells=1
        max_ashells=no_shells
        min_bshells=1
        max_bshells=no_shells
        end if         
        sigmx=2
        average=0.0
        call reset_string_order()
        max_op=0
        if (npintn) then
            taux=0
            four2=2.0
        else
            four2=4.0
            taux=2
        end if
        do sh1=min_ashells,max_ashells
!        do sh1=1,no_shells
                parti=0
                if (npintn) then
                    parti%shell(1)=sh1
                else
                    parti%shell(1)=ishell(sh1)
                end if
                partm(1:1)=char(sh1)
                j1=sj2(sh1)
                p1=sp(sh1)
                l1=2*sl(sh1)
        !loop over particle 2  shells 
                do sh2=min_bshells,max_bshells
!                do sh2=1,no_shells
                        if (shellxt.and.shtype(sh1)==shtype(sh2)) go to 30
                        if (npintn) then
                            parti%shell(2)=sh2
                        else
                            parti%shell(2)=ishell(sh2)
                        end if
                        partm(3:3)=char(sh2)
                        j2=sj2(sh2)
                        p2=sp(sh2)
                        l2=2*sl(sh2)
                        delta12=real(1,ri)
                        if (j1 ==j2 .and. ishell(sh1) == ishell(sh2) )&
                             delta12=sqrt(real(2,ri))
        !loop over particle 3 , shells 
                        do sh3=min_ashells,max_ashells
!                        do sh3=1,no_shells
                                if (npintn) then
                                    parti%shell(3)=sh3
                                else
                                    parti%shell(3)=ishell(sh3)
                                end if
                                partm(2:2)=char(sh3)
                                j3=sj2(sh3)
                                p3=sp(sh3)
                                l3=2*sl(sh3)
        !loop over particle 4 , shells 
                                do sh4=min_bshells,max_bshells
!                                do sh4=1,no_shells
                                   if (shellxt.and.shtype(sh3)==shtype(sh4)) go to 10
                                        if (npintn) then
                                            parti%shell(4)=sh4
                                        else
                                            parti%shell(4)=ishell(sh4)
                                        end if
                                        partm(4:4)=char(sh4)
                                        j4=sj2(sh4)
                                        p4=sp(sh4)
                                        l4=2*sl(sh4)
                                        delta34=real(1,ri)
!                                        if (p1*p2==p3*p4.or.p1*p2==-p3*p4) then
                                        if (j3 == j4 .and. ishell(sh3) == &
                                          ishell(sh4)) delta34=sqrt(real(2,ri))
                                          lmn=2*max(abs(sl(sh1)-sl(sh3)),abs(sl(sh2)-sl(sh4)))
                                          lmx=2*min(abs(sl(sh1)+sl(sh3)),abs(sl(sh2)+sl(sh4)))
                                          max_lmo=max(max_lmo,lmx)
                                          
                                         do  lam=lmx,lmn,-2 
                                         do  tau=0,taux,2
                                         do  sigma=0,sigmx,2
                                          partm(5:5)=char(lam)
                                          if (shellxt) partm(6:6)=char(tau)
                                          sumV=real(0,ri)
                           ! loop over J values from j1,j2
                                          do J=max(abs(j1-j2),abs(j3-j4)),min(abs(j1+j2),abs(j3+j4)),2
                                          do cS=0,2,2
                                          do cL=abs(J-cS),abs(J+cS),2
                                            Tmin=0
                                            Tmax=2
                           ! check for j1=j2 or j3=j4 and limit to J+T EVEN
                                            if (ishell(sh1) == ishell(sh2)&
                                                .or.ishell(sh3) == ishell(sh4)) then
                                                if (btest(J/2,0)) then
                                                    Tmax=0
                                                    if (intgj) Tmin=2
                                                else
                                                    Tmin=2
                                                    if (intgj) Tmax=0
                                                end if
                                            end if
                            ! loop over allowed T and check tz
                                            do T=Tmin,Tmax,2
                                                if (btest((T+cL+cS)/2,0)) then
                                                parti%shell(5)=J/2
                                                parti%shell(6)=T/2
                            ! find V
                                                if (SDI) then
                                                    call sdi_int(sh1,sh2&
                                                      ,sh3,sh4,j1,j2,j3,j4,J,T,delta12 &
                                                                 ,delta34,V)
                                                else
                                                  V=findV(parti,j1,j2,j3,j4,J,T)
                                                  if (V==0.0.and.p1*p2*p3*p4==1) warn=.true.
                                                  if (V/=0.0.and.p1*p2*p3*p4==-1) warnp=.true.
                                                end if
                                                if (V /= real(0,ri)) then
                                                    phase=real(1,rc)
                                                    if (btest((l1+l4-cL-lam)/2,0)) phase=-phase
                                                    vint=V*sqrt(real((sj2(sh1)+1)*(sj2(sh1)+1)*&
                                                    (sj2(sh1)+1)*(sj2(sh1)+1),rc))*real((cL+1)*&
                                                    (cS+1),rc)*ninej(l1,1,j1,l2,1,j2,cL,cS,J)*&
                                                    ninej(l3,1,j3,l4,1,j4,cL,cS,J)*&
                                                    real((cL+1),rc)*sqrt(real((lam+1),rc))&
                                                    *phase*racah(l1,l2,l3,l4,cL,lam)*&
                                                    delta12*delta34/4.0d0
                                                    if (btest((2-cS-sigma)/2,0)) vint=-vint
                                                    vint=vint*real((cS+1),rc)*sqrt(real((sigma+1),rc))&
                                                    *racah(1,1,1,1,cS,sigma)
                                                    if (.not.npintn) then
                                                        if (btest((2-T-tau)/2,0)) vint=-vint
                                                        vint=vint*real((T+1),rc)*sqrt(real((tau+1),rc))&
                                                        *racah(1,1,1,1,T,tau)
                                                    end if
                                                    sumV=sumV+vint
                                                end if
                                                end if
                                            end do    
                                          end do
                                          end do
                                          end do
                             !write out me
                                          me=sumV
                                          if (.not.shellxt) then
                                                partm(6:8)=z3(1:3)
                                          else
                                                partm(7:8)=z3(1:2)
                                          end if
                                          call add_string(partm,rhs)
                                          if (rhs)then
                                          max_op=max_op+1 
                                          meo(max_op)=me
                                          pro(max_op)=p2*p4
                                          plo(max_op)=p1*p3
                                          write(22,*) parti%shell(1:4),lam,tau,sigma,me
                                          end if
                                          lmid=(sh1-1)*no_shells+sh3
                                          rmid=(sh2-1)*no_shells+sh4
                                          somega(lmid,rmid,lam/2,tau/2,sigma/2)=sumV
                                          matrix(lmid,rmid,lam/2,tau/2,sigma/2)=sumV
                                          
                                         end do
                                         end do
                                         end do
10                             end do
20                     end do
30             end do
40     end do
       done=.true.
        i=0
        if (max_op/=0) then
        allocate(opn(max_op))
        allocate(meV(max_op))
        allocate(pl(max_op))
        allocate(pr(max_op))
        do
            k=get_next_string(done)
            i=i+1
            opn(i)=string_order(k)%string(1:8)
            meV(i)=meo(k)
            pl(i)=plo(k)
            pr(i)=pro(k)
            if (done) exit
        end do
        if (i/=max_op) stop ' maxop error'
        end if
    open(unit=50,file=intfile//'.e',action='write')
    open(unit=55,file=intfile//'.evc',action='write')
    open(unit=56,file=intfile//'.eve',action='write')
    lalpha=.false.
    lmatrix=.false.
    indx=0
    do lam=0,max_sj2u
    do sigma=0,1
    do tau=0,1
        call tred2(matrix(:,:,lam,tau,sigma),no_shells*no_shells,no_shells*no_shells,alpha(:,lam,tau,sigma),beta)
        tqlim=2000
        call tqli(alpha(:,lam,tau,sigma),beta,no_shells*no_shells,no_shells*no_shells,matrix(:,:,lam,tau,sigma),tqlim)
        
        if (lam==1.and.tau==0.and.sigma==0) then
            if (rem_Hcm) &
            call remove_Hcm(matrix(:,:,lam,tau,sigma),alpha(:,lam,tau,sigma),indx(:,lam,tau,sigma),no_shells,hcmvec,hcm2c,limit)
            if (add_SU3)then
                    sh1=minval(minloc(abs(alpha(:,lam,tau,sigma))))
                    print *,sh1
                    if (abs(alpha(sh1,lam,tau,sigma))<.0001) then
                         alpha(sh1,lam,tau,sigma)=sqrt(3.0d0)*rll
                         matrix(:,sh1,lam,tau,sigma)=hcm2c(:)
                    else
                        print *, ' No place to add L^2.'
                    end if
             end if
        end if                             

        if (lam==2.and.tau==0) then
            if (add_SU3)then
                    sh1=minval(minloc(abs(alpha(:,lam,tau,sigma))))
                    if (abs(alpha(sh1,lam,tau,sigma))<.0001) then
                         alpha(sh1,lam,tau,sigma)=-sqrt(5.0d0)*rqq
                         matrix(:,sh1,lam,tau,sigma)=hcmvec(:)
                    else
                        print *, ' No place to add Q.Q.'
                    end if
             end if
        end if                             

        call indexx(no_shells*no_shells,abs(alpha(:,lam,tau,sigma)),indx(:,lam,tau,sigma))      

        write (50,*) ' Lambda,tau,sigma',lam,tau,sigma
        write (55,*) ' Lambda,tau,sigma',lam,tau,sigma
        write (50,'(10f8.3)') alpha(indx(1:no_shells*no_shells,lam,tau,sigma),lam,tau,sigma)&
            /sqrt(real(2*lam+1,8))/sqrt(real(2*tau+1,8))
        do i=1,no_shells*no_shells
        if (abs(alpha(indx(i,lam,tau,sigma),lam,tau,sigma))/sqrt(real(2*lam+1,8))&
                        /sqrt(real(2*tau+1,8))>0.0000001) then
            lalpha(indx(i,lam,tau,sigma),lam,tau,sigma)=.true.
            write(55,'(a5,f14.9)') ' e =',alpha(indx(i,lam,tau,sigma),lam,tau,sigma)
            if (lalpha(indx(i,lam,tau,sigma),lam,tau,sigma)) then
            do sh1=1,no_shells
            do sh2=1,no_shells
            lmid=(sh1-1)*no_shells+sh2
            if (abs(matrix(lmid,indx(i,lam,tau,sigma),lam,tau,sigma))>0.00000001) then
                lmatrix(lmid,indx(i,lam,tau,sigma),lam,tau,sigma)=.true.
                write(55,'(2i4,5f14.9)') sh1,sh2,matrix(lmid,indx(i,lam,tau,sigma),lam,tau,sigma)
            end if
            end do  
            end do
            end if
         end if
         end do
    end do
    end do
    end do
    close(unit=55)
    matrix=0.0d0
    alpha=0.0d0
    if (do_edit) then
    print *,  ' Please edit the following file. '
    print *,  ' Do not delete or add items or lines. '
    print *,  ' You can remove terms by setting e or evec values to 0.0'   
    i=runqq('NotePad',intfile//'.evc')
    end if
    open(unit=55,file=intfile//'.evc')
    open(unit=60,file=intfile//'z.dsp',action='write')
    write(60,'(i4)') max_sj2u
    write(56,*) ' e-cut,evec-cut',rcute,rcutm
    do lam=0,max_sj2u
    write(60,'(i4)') lam
    sh1sh2=.false.
    do sigma=0,1
    do tau=0,1
        read (55,*) 
        write(56,*) ' Lambda,tau,sigma',lam,tau,sigma
        do i=1,no_shells*no_shells
        if (lalpha(indx(i,lam,tau,sigma),lam,tau,sigma)) then
            read(55,'(6x,f14.9)') alpha(indx(i,lam,tau,sigma),lam,tau,sigma)
            if (lam>0) then
            if (do_cut.and.abs(alpha(indx(i,lam,tau,sigma),lam,tau,sigma))<rcute*sqrt(real((2*lam+1)*(2*tau+1),ri)))  &
                alpha(indx(i,lam,tau,sigma),lam,tau,sigma)=0.0d0
            end if  
            write(56,'(a5,f14.9)') 'e = ',alpha(indx(i,lam,tau,sigma),lam,tau,sigma)
            if (lalpha(indx(i,lam,tau,sigma),lam,tau,sigma)) then 
            do sh1=1,no_shells
            do sh2=1,no_shells
            lmid=(sh1-1)*no_shells+sh2
            if (lmatrix(lmid,indx(i,lam,tau,sigma),lam,tau,sigma)) then
                read(55,'(2i4,f14.9)') sh3,sh4,matrix(lmid,indx(i,lam,tau,sigma),lam,tau,sigma)
                if (lam>0)  then
                if (do_cut.and.abs(matrix(lmid,indx(i,lam,tau,sigma),lam,tau,sigma))<rcutm) &
                    matrix(lmid,indx(i,lam,tau,sigma),lam,tau,sigma)=0.d0
                end if
                if (abs(matrix(lmid,indx(i,lam,tau,sigma),lam,tau,sigma))>=rcutm.and.lam>0) then
                    sh1sh2(lmid)=.true. 
                    write(56,'(2i4,f14.9)') sh1,sh2,matrix(lmid,indx(i,lam,tau,sigma),lam,tau,sigma)
                end if   
                if (lam==0) then
                    sh1sh2(lmid)=.true. 
                    write(56,'(2i4,f14.9)') sh1,sh2,matrix(lmid,indx(i,lam,tau,sigma),lam,tau,sigma)
                end if   
            end if
            end do  
            end do
            end if
         end if
         end do
         omega(:,:,lam,tau,sigma)=0.0d0
         do sh1=1,no_shells*no_shells
         omega(sh1,sh1,lam,tau,sigma)=alpha(sh1,lam,tau,sigma)
         end do
         omega(:,:,lam,tau,sigma)=matmul(matrix(:,:,lam,tau,sigma),&
            matmul(omega(:,:,lam,tau,sigma),transpose(matrix(:,:,lam,tau,sigma))))
     end do
     do sh1=1,no_shells
     do sh2=1,no_shells
     lmid=(sh1-1)+sh2
     if (sh1sh2(lmid)) then
            write(60,'(4i4)') sh1,sh2,sh1+no_shells,sh2+no_shells
     end if    
     end do   
     end do
     end do   
     write(60,'(4i4)') -1,-1,-1,-1
    end do
    close(unit=55)
    close(unit=56)
    print *,  ' Please do NOT edit the following file. '
    print *,  ' Do not delete or add items or lines. '
    print *,  ' It is the FINAL edited file. Please close it when done.'   
    i=runqq('NotePad',intfile//'.eve')
    
            
    write(20) max_op,max_lmo
    write(22,*) max_op,max_lmo
    if (max_op/=0) then       
    do i=1,max_op
    write(20) opn(i),meV(i),pl(i),pr(i)
    write(22,'(8i4,f8.4,2i4)') ichara(opn(i),6),meV(i),pl(i),pr(i)
    no_Op=no_Op+1
    end do
    end if    
    
       if (warn) then
           print *, ' WARNING some j1,j2,j3,j4,J,T ME`s were zero or missing!'
       end if 
       if (warnp) then
           print *, ' WARNING some j1,j2,j3,j4,J,T ME`s do not conserve parity!'
       end if 
     nz_op=0
     max_op=n_int 
     do no_Op=1,n_int
     sh1=order_int(no_Op)%shell(1)
     sh2=order_int(no_Op)%shell(2)
     sh3=order_int(no_Op)%shell(3)
     sh4=order_int(no_Op)%shell(4)
     if (ishell(sh1)==ishell(sh2)) delta12=1.0/sqrt(2.0)
     if (ishell(sh3)==ishell(sh4)) delta34=1.0/sqrt(2.0)
     if (ishell(sh1)/=ishell(sh2)) delta12=1.0
     if (ishell(sh3)/=ishell(sh4)) delta34=1.0
     lmid=(sh1-1)*no_shells+sh3
     rmid=(sh2-1)*no_shells+sh4
     J=2*order_int(no_Op)%shell(5)
     T=2*order_int(no_Op)%shell(6)
     me=sign(1.0,order_V(no_Op))
     if (npintn)  T=0
     sumV=0.0
     do lam=0,max_sj2u*2,2
     do sigma=0,2,2
     do tau=0,taux,2
        V=omega(lmid,rmid,lam/2,tau/2,sigma/2)
        if (sh1==sh4.and.sh2==sh3) V=somega(lmid,rmid,lam/2,tau/2,sigma/2)
        if (V /= real(0,ri)) then
           phase=real(1,rc)
           if (btest((sj2(sh1)+sj2(sh4)-J-lam)/2,0)) phase=-phase
           vint=four2*V*phase*racah(sj2(sh1),sj2(sh2),sj2(sh3),sj2(sh4),J,lam)&
                    *sqrt(real((lam+1),rc))*delta12*delta34
           if (.not.npintn) then
               if (btest((2-T-tau)/2,0)) vint=-vint
               vint=vint*racah(1,1,1,1,T,tau)*sqrt(real((tau+1),rc))
           end if
           sumV=sumV+vint
        end if
     end do   
     end do
     end do
     order_V(no_Op)=sumV
     if (npoutn) then
        if (T==2) then
        max_op=max_op+1
        parti=order_int(no_Op)
        parti%shell(1)=sh1+no_shells/2
        parti%shell(2)=sh2+no_shells/2
        parti%shell(3)=sh3+no_shells/2
        parti%shell(4)=sh4+no_shells/2
        order_int(max_op)=parti
        order_V(max_op)=sumV
        max_op=max_op+1
        parti=order_int(no_Op)
        parti%shell(1)=sh1
        parti%shell(2)=sh2+no_shells/2
        parti%shell(3)=sh3
        parti%shell(4)=sh4+no_shells/2
        order_int(max_op)=parti
        order_V(max_op)=sumV
        end if
        if (sumV/=0.0) then
        if (sh1/=sh2) then
        max_op=max_op+1
        parti=order_int(no_Op)
        parti%shell(2)=sh1+no_shells/2
        parti%shell(1)=sh2
        parti%shell(3)=sh3
        parti%shell(4)=sh4+no_shells/2
        order_int(max_op)=parti
        order_V(max_op)=sumV*(-1.0)**((sj2(sh1)+sj2(sh2)-T-J)/2)
        end if
        if (sh3/=sh4) then
        max_op=max_op+1
        parti=order_int(no_Op)
        parti%shell(1)=sh1
        parti%shell(2)=sh2+no_shells/2
        parti%shell(4)=sh3+no_shells/2
        parti%shell(3)=sh4
        order_int(max_op)=parti
        order_V(max_op)=sumV*(-1.0)**((sj2(sh3)+sj2(sh4)-T-J)/2) 
        end  if
        if (sh1/=sh2.and.sh3/=sh4) then
        max_op=max_op+1
        parti=order_int(no_Op)
        parti%shell(2)=sh1+no_shells/2
        parti%shell(1)=sh2
        parti%shell(4)=sh3+no_shells/2
        parti%shell(3)=sh4
        order_int(max_op)=parti
        order_V(max_op)=sumV*(-1.0)**((sj2(sh1)+sj2(sh2))/2)&
                            *(-1.0)**((sj2(sh3)+sj2(sh4))/2)
        end if
        if (T==0) then
            parti=order_int(no_Op)
            parti%shell(1)=sh1
            parti%shell(2)=sh2+no_shells/2
            parti%shell(3)=sh3
            parti%shell(4)=sh4+no_shells/2
            order_int(no_Op)=parti
            order_V(no_Op)=sumV
        end if
     end if 
     end if
     end do
        print * , ' No of interaction matrix elements',max_op-nz_op
        open(unit=30,file=intfile//'z.int',action='write') 
        write(30,'(a80)') first_line
        write(30,'(a42)',err=400) '! DZ Transformed and Truncated Interaction'       
        if (npoutn) then
        write(30,'(i4,50f9.4)',err=400) max_op-nz_op,Esp
        else
        write(30,'(i4,50f9.4)',err=400) max_op,Esp,Esp
        end if
        do j=1,max_op
            if (order_V(j) /= real(0,ri)) write(30,'(4i4,2i6,f9.4)') &
                  order_int(j)%shell(1:6),order_V(j)
        end do
        close(30)
       return
        
100    print *, 'Error opening ', Nucleus//xtra//inExt
       stop     
200    print *, 'Error opening ', Nucleus//xtra//'.misV'
       stop     
300    print *, 'Error writing ', Nucleus//xtra//inExt
       stop     
400    print *, 'Error writing ', Nucleus//xtra//'.misV'
       stop     
500    print *, 'Error opening ', Nucleus//xtra//'tint.txt'
       stop     
600    print *, 'Error writing ', Nucleus//xtra//'tint.txt'
       stop     
    end subroutine h_op
                                                
    function findV(parti,j1,j2,j3,j4,J,T)
    implicit none
        type(spartition),intent(in):: parti
        real(kind=ri):: findV
        integer,intent(in):: j1,j2,j3,j4,J,T
        type(spartition):: perm_part,part
        logical:: done,btest
        integer:: jtsum,findp
        real(kind=ri):: phase
        integer:: j1c,j2c,j3c,j4c
        
           done=.false.
           notf=.false. 
           lhs=.false.
           rhs=.false.
           
           part=parti    
           j1c=j1
           j2c=j2
           j3c=j3
           j4c=j4
10         findp=bfind(order_int,n_int,part)
           if (findp /= 0) then
               findV=order_V(findp)
               return
           else
               perm_part=part
               perm_part%shell(2)=part%shell(1)
               perm_part%shell(1)=part%shell(2)
               findp=bfind(order_int,n_int,perm_part)
               if (findp /= 0) then
                   jtsum=j1c+j2c-J-T
                   phase=real(1,ri)
                   if (btest(jtsum/2,0)) phase=-phase
                   findV=phase*order_V(findp)
                   lhs=.true.
                   return
               else
                   perm_part=part
                   perm_part%shell(3)=part%shell(4)
                   perm_part%shell(4)=part%shell(3)
                   findp=bfind(order_int,n_int,perm_part)
                   if (findp /= 0) then
                       jtsum=j3c+j4c-J-T
                       phase=real(1,ri)
                       if (btest(jtsum/2,0)) phase=-phase
                       findV=phase*order_V(findp)
                       rhs=.true.
                       return
                   else
                       perm_part=part
                       perm_part%shell(2)=part%shell(1)
                       perm_part%shell(1)=part%shell(2)
                       perm_part%shell(3)=part%shell(4)
                       perm_part%shell(4)=part%shell(3)
                       findp=bfind(order_int,n_int,perm_part)
                       if (findp /= 0) then
                           jtsum=j1c+j2c+j3c+j4c
                           phase=real(1,ri)
                           if (btest(jtsum/2,0)) phase=-phase
                           findV=phase*order_V(findp)
                           return
                       else
                           if (done) then
                                findV=real(0,ri)
                                notf=.true.
                                return
                           end if
                           perm_part=part
                           part%shell(1)=perm_part%shell(3)
                           part%shell(2)=perm_part%shell(4)
                           part%shell(3)=perm_part%shell(1)
                           part%shell(4)=perm_part%shell(2)
                           done=.true.
                           j1c=j3
                           j2c=j4
                           j3c=j1
                           j4c=j2
                           go to 10
                       end if
                   end if
               end if
           end if    
                                                                                         
    end function findV
    
    subroutine sdi_int(sh1,sh2,sh3,sh4,j1,j2,j3,j4,J,T,delta12,delta34,V)
        implicit none
        integer,intent(in):: sh1,sh2,sh3,sh4,j1,j2,j3,j4,J,T
        real(kind=ri),intent(in):: delta12,delta34
        real(kind=ri),intent(out):: V
        real(kind=ri):: pp1,pp2,fct,Vd,Ve
        type(spartition):: part
        logical:: found,btest
        
        part%shell(1)=sh1
        part%shell(2)=sh2
        part%shell(3)=sh3
        part%shell(4)=sh4
        part%shell(5)=J/2
        part%shell(6)=T/2
        V=find_part_V(part,found) 
        if (found) return
        V=real(0,ri)
        Vd=real(0,ri)
        Ve=real(0,ri)
! MSDI Brussard and Glaudemans, north holland. shell model book
            V=A(T/2)*(P)**(sn(sh1)+sn(sh2)+sn(sh3)+sn(sh4))*sqrt(real((j1+1)*&
               (j2+1)*(j3+1)*(j4+1),ri)/real(4*(J+1)*(J+1),ri)/delta12/delta34)&
               *( (-1.0)**((j2+j4)/2+sl(sh2)+sl(sh4))*cleb(j2,j1,J,-1,1,0)&
               *cleb(j4,j3,J,-1,1,0)*(1 - (-1)**(sl(sh1)+sl(sh2)+(J+T)/2))&
               -cleb(j2,j1,J,1,1,2)*cleb(j4,j3,J,1,1,2)*(1 + (-1)**(T/2)) )
            if (sh1 == sh3 .and. sh2 == sh4) then
                V=V+C
                if (T==2) V=V+B
                if (T==0) V=V-3*B
            end if
! P2.P2 Quadrupole zero range WDMR with formulae Brink and Satchler+Lawson
! Seems to work but could be wrong. Beware!
                pp1=real(-1,ri)**(sl(sh1)+sl(sh3))
                pp2=real(-1,ri)**(sl(sh2)+sl(sh4))                  
                fct=Qf
                if(sp(sh1)==sp(sh2).and.sp(sh3)==sp(sh2).and.sp(sh3)==sp(sh4)) fct=real(1,ri)
                if (pp1 > 0.0 .and. pp2 >0.0) then
                    Vd=Qq*fct*4*B2*B2*real(-1,ri)**((j3+j2+J)/2)*sixj(j1,j3,4,j4,j2,J)*&
                       sqrt(real((j3+1)*(j4+1),ri))*real(-1,ri)**((-j1-j2+j3+j4)/2)&
                       *cleb(j4,4,j2,1,0,1)*pp1*pp2*cleb(j3,4,j1,1,0,1)/delta12/delta34 &
                       *(R)**(sn(sh1)+sn(sh2)+sn(sh3)+sn(sh4))*radial(sn(sh1),sl(sh1),&
                       sn(sh3),sl(sh3),2)*radial(sn(sh2),sl(sh2),sn(sh4),sl(sh4),2)/real(3,ri)
                end if
! Monopole term  ri^2.rj^2
                if (pp1 > 0.0 .and. pp2 >0.0) then
                     Vd=Vd+Qm*fct*B2*B2*5*real(-1,ri)**((j3+j2+J)/2)*sixj(j1,j3,0,j4,j2,J)*&
                        sqrt(real((j3+1)*(j4+1),ri))*real(-1,ri)**((-j1-j2+j3+j4)/2)&
                        *cleb(j4,0,j2,1,0,1)*pp1*pp2*cleb(j3,0,j1,1,0,1)/delta12/delta34 &
                        *(R)**(sn(sh1)+sn(sh2)+sn(sh3)+sn(sh4))*radial(sn(sh1),sl(sh1),&
                        sn(sh3),sl(sh3),2)*radial(sn(sh2),sl(sh2),sn(sh4),sl(sh4),2)/real(3,ri)
                end if
! Hex term  ri^33.rj^3
                if (pp1 < 0.0 .and. pp2 < 0.0) then
                     Vd=Vd+Qh*B2*B2*B2*8*real(-1,ri)**((j3+j2+J)/2)*sixj(j1,j3,6,j4,j2,J)*&
                        sqrt(real((j3+1)*(j4+1),ri))*real(-1,ri)**((-j1-j2+j3+j4)/2)&
                        *cleb(j4,6,j2,1,0,1)*pp1*pp2*cleb(j3,6,j1,1,0,1)/delta12/delta34 &
                        *(R)**(sn(sh1)+sn(sh2)+sn(sh3)+sn(sh4))*radial(sn(sh1),sl(sh1),&
                        sn(sh3),sl(sh3),3)*radial(sn(sh2),sl(sh2),sn(sh4),sl(sh4),3)/real(15,ri)
                end if
! Antisymmetrise P2.P2
                pp1=real(-1,ri)**(sl(sh1)+sl(sh4))
                pp2=real(-1,ri)**(sl(sh2)+sl(sh3))
                fct=Qf
                if(sp(sh1)==sp(sh2).and.sp(sh3)==sp(sh2).and.sp(sh3)==sp(sh4)) fct=real(1,ri)
                if (pp1 > 0.0 .and. pp2 > 0.0) then
                     Ve=Qq*fct*4*B2*B2*real(-1,ri)**((j4+j2+J)/2)*sixj(j1,j4,4,j3,j2,J)*&
                        sqrt(real((j3+1)*(j4+1),ri))*real(-1,ri)**((-j1-j2+j3+j4)/2)&
                        *cleb(j3,4,j2,1,0,1)*pp1*pp2*cleb(j4,4,j1,1,0,1)/delta12/delta34 &
                        *(R)**(sn(sh1)+sn(sh2)+sn(sh3)+sn(sh4))*real(-1,ri)**((J+T-j3-j4)/2)*&
                        radial(sn(sh1),sl(sh1),sn(sh4),sl(sh4),2)&
                        *radial(sn(sh2),sl(sh2),sn(sh3),sl(sh3),2)/real(3,ri)
                end if
!monopole ri^2.rj^2
                if (pp1 > 0.0 .and. pp2 > 0.0) then
                      Ve=Ve+Qm*fct*B2*B2*5*real(-1,ri)**((j4+j2+J)/2)*sixj(j1,j4,0,j3,j2,J)*&
                      sqrt(real((j3+1)*(j4+1),ri))*real(-1,ri)**((-j1-j2+j3+j4)/2)&
                      *cleb(j3,0,j2,1,0,1)*pp1*pp2*cleb(j4,0,j1,1,0,1)/delta12/delta34 &
                      *(R)**(sn(sh1)+sn(sh2)+sn(sh3)+sn(sh4))*real(-1,ri)**((J+T-j3-j4)/2)*&
                      radial(sn(sh1),sl(sh1),sn(sh4),sl(sh4),2)&
                      *radial(sn(sh2),sl(sh2),sn(sh3),sl(sh3),2)/real(3,ri)
                end if
!hex term  ri^3.rj^3
                if (pp1 < 0.0 .and. pp2 < 0.0) then
                      Ve=Ve+Qh*8*B2*B2*B2*real(-1,ri)**((j4+j2+J)/2)*sixj(j1,j4,6,j3,j2,J)*&
                      sqrt(real((j3+1)*(j4+1),ri))*real(-1,ri)**((-j1-j2+j3+j4)/2)&
                      *cleb(j3,6,j2,1,0,1)*pp1*pp2*cleb(j4,6,j1,1,0,1)/delta12/delta34 &
                      *(R)**(sn(sh1)+sn(sh2)+sn(sh3)+sn(sh4))*real(-1,ri)**((J+T-j3-j4)/2)*&
                      radial(sn(sh1),sl(sh1),sn(sh4),sl(sh4),3)&
                      *radial(sn(sh2),sl(sh2),sn(sh3),sl(sh3),3)/real(15,ri)
                end if
                V=V+Vd+Ve
!pairing
            if (sh1 == sh2 .and. sh3 == sh4 .and. sh2 == sh3) then
                 V=V + S/real((sj2(sh1)+1),rc)
            end if
!            call add_part(part,V,.false.)


    end subroutine sdi_int
    
end module Hoper    

!  NuDZop.f90 
!
!  FUNCTIONS:
!  NuDZop      - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: NuDZop
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program NuDZop
    
    use InputNuShell
    use OutputNuShell
    use Shells
    use ShellsX
    use Mscheme
    use OperSupport
    use PartOrder  
    use Hoper
    use Files
    use Extra
      
    implicit none

    ! Variables
    
    real:: t1,t2
    
    character(len=6):: inExt
    character(len=2):: inttyp,outtyp
    character(len=11):: chint='INTERACTION'
    character(len=4):: chtree='TREE',chnode='NODE'  
    character(len=10):: chnoint,chnotree,chnonode
    integer:: res,nonodeo
    
    res=get_env_v(chint,chnoint)
    if (res /= 0) then
        read(chnoint,*) noint
    else
        noint=max_int
    end if
    allocate(order_int(noint))
    allocate(order_V(noint))
    
    res=get_env_v(chnode,chnonode)
    if (res /=0) then
        read(chnonode,*) nonode
    else
        nonode=max_node
    end if
    allocate(node(nonode))
    
    res=get_env_v(chtree,chnotree)
    if (res /=0) then
        read(chnotree,*) notree
    else
        notree=max_part_tree
    end if
    allocate(ipart_order(notree))
    
    ! Body of NuProj
    call cpu_time(t1)
    
    call output_header(6)
    call output_NuDZop(6)

    call input_nucleus(inExt,xtra)
    Nucleus=inExt    
    call output_welcome(6)    
    call input_sps(inExt)
    call output_shells(6)    
    call check_shell_order()
    call setup_shells()
    call setup_mscheme()
    call setup_lu()
    call setup_ishell()
 
    call set_ptdim(6)
    call read_interaction('.int')
    call output_int_read(6)    
    call input_type(outtyp)
!    if (inttyp=='np') then
!        npintn=.true.
!    else
        npintn=.false.
!    end if   
    if (outtyp=='np') then
        npoutn=.true.
    else    
        npoutn=.false.
    end if
    call output_end_setup(6)    
    call H_op('.oph')
    deallocate(order_int)
    deallocate(order_V)
    deallocate(node)
    call output_completed(6,'NuLnz/TRL ')   
    call cpu_time(t2)
    t1=t2-t1
    call output_time('NuDZop',t1,t1)
    print *
    
    end program NuDZop

    subroutine output_NuDZop(dev)
    
    use InputNuShell
    use OutputNuShell
    use Shells
    use ShellsX
    use Mscheme
    use OperSupport  
    use Hoper
    use Files
      
    implicit none
    integer,intent(in):: dev
    integer:: status
    character(len=14)::buffer
    call get_cmd_arg(1,buffer,status) 
    if (status /= -1) return
    
    write(dev,*) ' NuDZop Program'
    write(dev,*)
    write(dev,*) ' Compile time Parameters:'
    write(dev,*)
    write(dev,*) ' Machine Parameters:'
    write(dev,*) ' No of bits per word',nbits,' Memory Access Width Bytes', nbytes
    write(dev,*)
    write(dev,*) ' Shell Model Parameters:'
    write(dev,*) ' Maximum no of shells',max_no_shells
    write(dev,*) ' Maximum no of words with',nbits,' bits is',nwords, '    for each m-scheme state'
    write(dev,*) ' Maximum value of 2*j for a shell',max_sj2
    write(dev,*) ' Real kind for interaction', ri
    write(dev,*) ' Real kind for CG coefs', rc
    write(dev,*) ' Maaximum no interaction matrix elements',max_int
    write(dev,*) ' Maxximum size of order binary tree',max_part_tree
    write(dev,*) ' Maximm no of nodes for tree readout',max_node
    write(dev,*)

    end subroutine output_NuDZop