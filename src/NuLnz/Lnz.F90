module Lnzs

use Globals
use XMltiply

implicit none

contains 

   subroutine Lnz_OP(rtbas)
   implicit none
   real,intent(in):: rtbas
        integer,parameter:: cut_orthog=100
        real(kind=rl),dimension(noiter):: alpha,beta,alpha_save,alpha_copy,beta_copy,ovo
        integer,dimension(noiter):: indx
        real(kind=rl),dimension(noiter,noiter):: eigenv
        real(kind=rl),dimension(noiter)::ovla
        integer:: n,i,imax,imin,ix,isave,np,pd,tqlim,status,idum,idum1,npre_iter=0
        real:: x,seig2,tots,totp,totn,ti,tf,tt=0.0
        real(kind=rl):: ovl,evalue,zero=0.0,eta,eta0
        integer::m,ntog=0
        integer,dimension(:),allocatable:: new
! use a fixed seed

        call random_seed(SIZE=m)
        allocate(new(m))
        new=5
        call random_seed(PUT=new(1:m))
        
        if (no_level>no_state) no_level=no_state
        eigenv=real(0,rl)
        alpha=real(0,rl)
        alpha_save=real(0,rl)
        alpha_copy=real(0,rl)
        beta_copy=real(0,rl)
        beta=real(0,rl)
        do n=1,no_state
            if (lcoef(n)) then
            call random_number(x)
            coefL(n)=0.5-x
            else
            coefL(n)=0.0
            end if
        end do
        ovl=dot(coefL,coefL,no_state)
        ovl=real(1,rl)/sqrt(ovl)
        coefR=ovl*coefL
        do i=1,npre_iter
            if (output_control>2) print *, ' Lanczos PreIteration ',i
            call cpu_time(ti)
! multiply x  must be called  first then multiply a then multiply b last         
            call multiplyx()
            call multiplya()
            call multiplyb()
            call cpu_time(tf)
            tt=tt+tf-ti
            ovl=dot(coefL,coefL,no_state)
            ovl=real(1,rl)/sqrt(ovl)
            coefR=ovl*coefL
        end do
!        if (output_control>2) print *, ' ELAPSED Time multiplication ',tt/2.0/npre_iter , ' secs'
!        print *, ' Normalised Time per pn multiplication ',tt/2.0/rtbas/npre_iter, ' secs'
        call lnzort(0)    
        beta(1)=real(0,rl)
        imin=1
        imax=min(first_iter,no_state)
10      do i=imin,imax
            if (output_control>2) print *, ' Lanczos Iteration ',i
            eigenv(i,i)=real(1,rl)
! multiply x  must be called  first then multiply a then multiply b last         
            call multiplyx()
            call multiplya()
            call multiplyb()
            if (i /= 1) coefL=coefL-beta(i)*coefT
            coefT=coefR
            alpha(i)=dot(coefR,coefL,no_state)
            if (i == no_state) go to 55
            coefL=coefL-alpha(i)*coefR
            ovl=dot(coefL,coefL,no_state)
            call lnzort(i)
            ntog=ntog+1
            ovl=dot(coefL,coefL,no_state)
            if (sqrt(abs(ovl)) > restart_limit) then
                beta(i+1)=sqrt(ovl)
                ovl=real(1,rl)/sqrt(ovl)
                coefR=ovl*coefL
            else
                print *, ' Lanczos restarted at iteration ',i
                do n=1,no_state
                    call random_number(x)
                    coefL(n)=1000.*x*coefL(n)
                end do
                call lnzort(i)
                ovl=dot(coefL,coefL,no_state)
                if (abs(ovl) < restart_limit) go to 55
                ovl=real(1,rl)/sqrt(ovl)
                coefR=ovl*coefL
                beta(i+1)=real(0,rl)
            end if         
!            coefR=coefL       
            call lnzort(-i)
        end do
        go to 56
55      imax=i
56      alpha_copy(1:imax)=alpha(1:imax)
        beta_copy(1:imax)=beta(1:imax)
        tqlim=2000
        call tqli(alpha_copy,beta_copy,imax,noiter,eigenv,tqlim)
        call indexx(imax,alpha_copy,indx)
        if (imax== no_state.or.abs(ovl) < restart_limit) go to 20
        do ix=1,min(no_state,no_level+1,imax)
            if (abs(alpha_save(ix)-alpha_copy(indx(ix))) > level_limit) then
                imin=imax+1
                imax=imin+inc_iter
                eigenv=real(0,rl)
                do i=1,imin-1
                    eigenv(i,i)=real(1,rl)
                end do
                if (imax +1> noiter) then
                    print *, ' Increase max_iter in LanczosParameters'
                    print *, ' Or SET ITERL=I where I >> ', noiter
                    imax=imin-1
                    go to 20
                end if
                alpha_save(1:min(no_state,no_level+1,imax))=&
                    alpha_copy(indx(1:min(no_state,no_level+1,imax)))
                go to 10
            end if
        end do
20      isave=min(no_state,imax,i)
        print *
        open(unit=701,file=Nucleuss//file_code(J2/2,no_cA,no_cB)//'.abz',form='unformatted',action='write')
        write(701) isave,alpha(1:isave),beta(1:isave)
        write(701) imax,no_level,tt,no_state,basd
        write(701) bas
        print '(5f10.3)', alpha_copy(indx(1:no_level))
        print *
        print *, ' No of iterations, northog',imax,ntog
   
   end subroutine Lnz_OP
   
   subroutine lnzort(i)
   implicit none
   integer,intent(in)::i
   integer::j,iol
   integer:: thrn,OMP_get_thread_num
   real(kind=rl) :: ovl
   
   if (i==0) then
        iol=no_state
        open(unit=700,file=Nucleuss//file_code(J2/2,no_cA,no_cB)//'.lnz',&
            form='unformatted',action='readwrite',access='direct',recl=iol)
        write(700,rec=1) coefR
        return
   else if (i>=nthr) then
        coefR=coefL
!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(coefX,coefR,coefL,i)
!$OMP DO ORDERED SCHEDULE(DYNAMIC,1)
        do j=1,i
            thrn=OMP_get_thread_num()+1
!$OMP ORDERED
            read(700,rec=j) coefX(:,thrn)
!$OMP END ORDERED            
            ovl=dot_product(coefR,coefX(:,thrn))
!$OMP CRITICAL
            coefL=coefL-ovl*coefX(:,thrn)
!$OMP END CRITICAL           
        end do
!$OMP END PARALLEL
   else if (i>0.and.i<nthr) then
        do j=1,i
            read(700,rec=j) coefX(:,i)
            ovl=dot_product(coefL,coefX(:,i))
            coefL=coefL-ovl*coefX(:,i)
        end do
   else if (i<0) then
        write(700,rec=-i+1) coefR
   end if     
        
   end subroutine lnzort
   
end module Lnzs