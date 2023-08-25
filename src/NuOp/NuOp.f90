! This file contains the modules specific the the generic pOper
! program for NuShell, NuShell@MSU and SunShell.

module OperSupport

use InputNuShell
use OutputNuShell
use PartOrder
use Files
use OperParameters
use TrdChar

implicit none

integer,dimension(max_no_shells):: ishell
integer:: ni_shells,n_int,Jmax,Tmax,noint
type(spartition):: jsJT
type(spartition),dimension(:),allocatable:: order_int
real(kind=ri):: V
real(kind=ri),dimension(:),allocatable:: order_V
logical:: SDI=.false.
real(kind=ri),dimension(0:1):: A
real(kind=ri):: B,C,Qq,P,R,S,LS,B2,X,Qsp,Qm,Qh,Qf
character(len=1) :: xtra
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
        character(len=7)::intfile
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
!        print *, ' Reading interaction file'
        call read_comments()   
        commentsv=comment     
!        if ( comment == '&') print *,' Concatented interaction'
!        if ( comment == '!') print *,' Standard single interaction'
        norm=real(1,ri)
1       read(10,*,err=200,end=201) n
        if (-n > max_no .and. comment == '&') then
            print *, ' Incease max_no'
            stop
        end if
        if (n < -1 .and. comment == '&') write(6,'(a19,i2,a2)',advance='no')&
                     'Normalisations n=1,',-n,'  '
        if (n < -1 .and. comment == '&') read(5,*) norm(1:-n)
        nn=1
        if (comment == '^') then
!            print *, ' MSDI + P2.P2 Interaction'
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

implicit none

logical:: notf,lhs,rhs,npintn
real(kind=8):: two=2.0d0
contains

    subroutine H_op(inExt)
    implicit none
    
        character(len=*):: inExt
        integer:: m1,m2,m3,m4,tz1,tz2,tz3,tz4,J,T,j1,j2,j3,j4,Tmin,Tmax,k,i,lam,max_op,tau
        integer:: sh1,sh2,sh3,sh4,iw1,iw2,iw3,iw4,ib1,ib2,ib3,ib4,p1,p2,p3,p4,lmn,lmx,taux
        real(kind=ri):: V,sumV,sumVn,sumVp,delta12,delta34,vint,vintn,vintp,phase
        real(kind=8)::me
        integer:: no_Op=0,one=1
        type(spartition):: parti
        character(len=8):: partm
        character(len=8):: z3
        character(len=8),dimension(:),allocatable:: opn
        real(kind=8),dimension(:),allocatable:: meV,meo
        integer,dimension(:),allocatable:: pl,pr,plo,pro
        logical:: warn=.false.,done,warnp=.false.
        integer:: opd,pnnp,iout,max_lmo=0
        character(len=10):: chnoop
        logical:: btest,intgj
        integer:: min_bshells,max_bshells,min_ashells,max_ashells

        intgj=.false.
        if (btest(sj2(1)+1,0)) intgj=.true.
        
        z3=repeat(char(0),8)
        nstring=notree
        allocate(string_order(notree))
        allocate(meo(notree))
        allocate(plo(notree))
        allocate(pro(notree))
        open(unit=20,file=Nucleus//xtra//inExt, form='unformatted', action='write',err=100)
        open(unit=21,file=Nucleus//xtra//'.misV', form='formatted', action='write',err=200)
        open(unit=22,file=Nucleus//xtra//'oph.'//'txt', &
                                     form='formatted', action='write')
       
        !loop over particle 1 shells
        call reset_string_order()
        max_op=0
        if (.not.shellxt) then
            taux=0
        else
            taux=2
            one=0
        end if
        do pnnp=0,one
        if (.not.shellxt) then
        if (pnnp==0) then
        
        typea='p'
        typeb='n'
        
        else
        
        typea='n'
        typeb='p'
        
        end if  
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
        do sh1=min_ashells,max_ashells
                parti=0
                if (npintn) then
                    parti%shell(1)=sh1
                else
                    parti%shell(1)=ishell(sh1)
                end if
                partm(1:1)=char(sh1)
                j1=sj2(sh1)
                p1=sp(sh1)
        !loop over particle 2  shells 
                do sh2=min_bshells,max_bshells
                        if (shellxt.and.shtype(sh1)==shtype(sh2)) go to 30
                        if (npintn) then
                            parti%shell(2)=sh2
                        else
                            parti%shell(2)=ishell(sh2)
                        end if
                        partm(3:3)=char(sh2)
                        j2=sj2(sh2)
                        p2=sp(sh2)
                        delta12=real(1,ri)
                        if (j1 ==j2 .and. ishell(sh1) == ishell(sh2) )&
                             delta12=sqrt(real(2,ri))
        !loop over particle 3 , shells 
                        do sh3=min_ashells,max_ashells
                                if (npintn) then
                                    parti%shell(3)=sh3
                                else
                                    parti%shell(3)=ishell(sh3)
                                end if
                                partm(2:2)=char(sh3)
                                j3=sj2(sh3)
                                p3=sp(sh3)
        !loop over particle 4 , shells 
                                do sh4=min_bshells,max_bshells
                                   if (shellxt.and.shtype(sh3)==shtype(sh4)) go to 10
                                   if (npintn) then
                                            parti%shell(4)=sh4
                                        else
                                            parti%shell(4)=ishell(sh4)
                                        end if
                                        partm(4:4)=char(sh4)
                                        j4=sj2(sh4)
                                        p4=sp(sh4)
                                        delta34=real(1,ri)
                                        if (p1*p2==p3*p4.or.p1*p2==-p3*p4) then
                                        if (j3 == j4 .and. ishell(sh3) == &
                                          ishell(sh4)) delta34=sqrt(real(2,ri))
                                          lmn=max(abs(j1-j3),abs(j2-j4))
                                          lmx=min(abs(j1+j3),abs(j2+j4))
                                          max_lmo=max(max_lmo,lmx)
                                         do  lam=lmx,lmn,-2 
                                         do  tau=0,taux,2
                                          partm(5:5)=char(lam)
                                          if (shellxt) partm(6:6)=char(tau)
                                          sumV=real(0,ri)
                           ! loop over J values from j1,j2
                                          do J=max(abs(j1-j2),abs(j3-j4)),min(abs(j1+j2),abs(j3+j4)),2
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
                                                    if (btest((j1+j4-J-lam)/2,0)) phase=-phase
                                                    vint=V*real((J+1),rc)*sqrt(real((lam+1),rc))&
                                                         *phase*racah(j1,j2,j3,j4,J,lam)&
                                                         *delta12*delta34/real(2,ri)
                                                    if (shellxt) then
                                                        if (btest((2-T-tau)/2,0)) vint=-vint
                                                        vint=vint*real(2,ri)*real((T+1),rc)*sqrt(real((tau+1),rc))&
                                                        *racah(1,1,1,1,T,tau)
                                                    end if
                                                    sumV=sumV+vint
                                                end if
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
                                          write(22,*) parti%shell(1:4),lam,tau,me
                                          end if
                                       end do
                                       end do
                                   end if
10                             end do
20                     end do
30             end do
40     end do
       
       end do
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
        
    write(20) max_op,max_lmo
    write(22,*) max_op,max_lmo
    if (max_op/=0) then       
    do i=1,max_op
    write(20) opn(i),meV(i),pl(i),pr(i)
    write(22,'(8i4,f8.4,2i4)') ichara(opn(i),6),meV(i),pl(i),pr(i)
    no_Op=no_Op+1
    end do
    end if

    if (.not.shellxt) then
    print *, ' Total number of PN matrix elements ',no_Op   
    else
    print *, ' Total number of ab matrix elements ',no_Op   
    end if
    
    
       if (warn) then
           print *, ' WARNING some j1,j2,j3,j4,J,T ME`s were zero or missing!'
       end if 
       if (warnp) then
           print *, ' WARNING some j1,j2,j3,j4,J,T ME`s do not conserve parity!'
       end if 
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

!  NuOp.f90 
!
!  FUNCTIONS:
!  NuOp      - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: NuOp
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program NuOp
    
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
    character(len=2):: inttyp
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
    call output_NuOp(6)

    call input_nucleus(inExt,xtra)
        
    call output_welcome(6)    
    call input(inExt)
    call output_shells(6)    
    call check_shell_order()
    call setup_shells()
    call setup_mscheme()
    call setup_lu()
    call setup_ishell()
 
    call set_ptdim(6)
    call read_interaction('.int')
    call output_int_read(6)    
    call input_type(inttyp)
    if (inttyp=='np') then
        npintn=.true.
    else
        npintn=.false.
    end if   
    call output_end_setup(6)    
    call H_op('.oph')
    deallocate(order_int)
    deallocate(order_V)
    deallocate(node)
    call output_completed(6,'NuLnz/TRL ')   
    call cpu_time(t2)
    t1=t2-t1
    call output_time('NuOp',t1,t1)
    print *
    
    end program NuOp

    subroutine output_NuOp(dev)
    
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
    
    write(dev,*) ' NuOp Program'
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

    end subroutine output_NuOp