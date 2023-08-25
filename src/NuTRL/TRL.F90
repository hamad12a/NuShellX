module TRL

use XMultiply

use Globals

implicit none

contains

   subroutine Lnz_OP(rtbas)
   implicit none
   real,intent(in):: rtbas
        integer,parameter:: cut_orthog=100
        real(kind=rl),dimension(noiter):: alpha,alpha_save,alpha_copy,beta_copy,ovo
        integer,dimension(noiter):: indx
        real(kind=rl),dimension(noiter,noiter):: eigenv,matrix
        real(kind=rl),dimension(noiter)::ovla,betat
        integer:: n,i,imax,imin,ix,isave,np,pd,tqlim,status,idum,idum1,npre_iter=0
        real:: x,seig2,tots,totp,totn,ti,tf,tt=0.0
        real(kind=rl):: ovl,evalue,zero=0.0,eta,eta0
        integer::m,ntog=0,second_iter=first_iterT/2,count=1,no_keep,no_conv,nti=0
        integer(kind=2):: res
        integer,dimension(:),allocatable:: new
        integer:: thrn,OMP_get_thread_num,idumm
        logical:: betaz,cont
! use a fixed seed

        call random_seed(SIZE=m)
        allocate(new(m))
        new=5
        call random_seed(PUT=new(1:m))
        
        if (no_level>no_state) no_level=no_state
! check if trl file already exists if so use that to start
	    inquire(file=Nucleuss//file_code(J2/2,no_cA,no_cB)//'.trl',exist=cont)
	    if (cont.and.restart>0) then
             open(unit=701,file=Nucleuss//file_code(J2/2,no_cA,no_cB)//'.abz',form='unformatted',action='read')
             read(701) imax,alpha(1:imax),betat(1:imax),matrix(1:imax,1:imax)
             read(701) isave,idumm,tt,idumm,idumm,count,no_keep,no_conv,idumm,&
                        idumm,idumm
	         read(701)
 	         read(701) alpha_copy(:),betat(imax+1),indx(:)	
	         close(701)
	         alpha_save(1:no_level)=alpha_copy(indx(1:no_level))
! for small cases just start again
	         if (4*imax>no_state) go to 5
             open(unit=700,file=Nucleuss//file_code(J2/2,no_cA,no_cB)//'.trl',&
                 form='unformatted',action='readwrite',access='direct',recl=no_state)
 	         read(700,rec=imax-1) coefT
 	         read(700,rec=imax) coefR      
             imin=imax
             imax=imin+second_iter/2
             if (imax +1> noiter) then
                 print *, ' Increase max_iter in LanczosParameters'
                 print *, ' Or SET ITERL=I where I >> ', noiter
		         stop
             end if
             go to 10
	    end if     
5       eigenv=real(0,rl)
        alpha=real(0,rl)
        alpha_save=real(0,rl)
        alpha_copy=real(0,rl)
        beta_copy=real(0,rl)
        betat=real(0,rl)
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
        imin=1
        imax=min(first_iterT,no_state)
10      do i=imin,imax
            betaz=.false.
            nti=nti+1
            if (output_control>2) print *, ' Lanczos Iteration ',nti
            eigenv(i,i)=real(1,rl)
! multiply x  must be called  first then multiply a then multiply b last         
            call multiplyx()
            call multiplya()
            call multiplyb()
            alpha(i)=dot(coefR,coefL,no_state)
            if (i /= 1.and.count==1) then
                coefL=coefL-betat(i-1)*coefT
            else if (count>1.and.i>no_keep+1) then
                coefL=coefL-betat(i-1)*coefT
            else if (count>1.and.i==no_keep+1) then
                do nsc=1,no_keep
                    read(700,rec=nsc) coefT
                    coefL=coefL-betat(nsc)*coefT
                end do
            end if
            if (i == no_state) go to 15
            coefL=coefL-alpha(i)*coefR
            ovl=dot(coefL,coefL,no_state)
            if (count==1) then
                coefT=coefR
                call lnzort(i)
                ntog=ntog+1
                ovl=dot(coefL,coefL,no_state)
                if (sqrt(ovl)<restart_limit) then
                    print *, ' Lanczos restarted at iteration ',i
                    do n=1,no_state
                        call random_number(x)
                        coefL(n)=1000.*x*coefL(n)
                    end do
                    call lnzort(i)
                    ntog=ntog+1
                    betaz=.true.
                end if
            end if
            if (count>1) then
                if (i==no_keep+1) then
                    eta=alpha(i)*alpha(i) + sum(betat(1:no_keep)*betat(1:no_keep))
                else
                    eta=alpha(i)*alpha(i) + betat(i-1)*betat(i-1)
                end if
                if (ovl>=eta) then
                    coefL=coefL-dot(coefL,coefR,no_state)*coefR-&
                          dot(coefL,coefT,no_state)*coefT
                else if (sqrt(ovl)>restart_limit.and.ovl<eta) then
                    coefT=coefR
                    call lnzort(i)
                    ntog=ntog+1
                else if (sqrt(ovl)<restart_limit) then
                    print *, ' Lanczos restarted at iteration ',i
                    do n=1,no_state
                        call random_number(x)
                        coefL(n)=1000.*x*coefL(n)
                    end do
                    call lnzort(i)
                    ntog=ntog+1
                    betaz=.true.
                end if
            end if
            ovl=dot(coefL,coefL,no_state)
            if (abs(ovl) < restart_limit) go to 15
            if (betaz) then
                betat(i)=real(0,rkw)
            else
                betat(i)=sqrt(ovl)
            end if
            ovl=real(1,rl)/sqrt(ovl)
            coefR=ovl*coefL
            call lnzort(-i)
        end do
        go to 16
15      imax=i
16      if (count==1) then
            alpha_copy(1:imax)=alpha(1:imax)
            beta_copy(2:imax+1)=betat(1:imax)
            beta_copy(1)=real(0,rk)
        else if (count>1) then
            do ix=imin,imax
                matrix(ix,ix)=alpha(ix)
                matrix(ix+1,ix)=betat(ix)
                matrix(ix,ix+1)=betat(ix)
            end do
            eigenv(1:imax,1:imax)=matrix(1:imax,1:imax)
            call tred2(eigenv,imax,noiter,alpha_copy,beta_copy)
        end if    
        tqlim=2000
        call tqli(alpha_copy,beta_copy,imax,noiter,eigenv,tqlim)
        call indexx(imax,alpha_copy,indx)
        if (imax== no_state.or.abs(ovl) < restart_limit) go to 20
        do ix=1,min(no_state,no_level+restart,imax)
            if (abs(alpha_save(ix)-alpha_copy(indx(ix))) > level_limit) then
                no_conv=ix-1 
                alpha_save(1:min(no_state,no_level+restart,imax))=&
                    alpha_copy(indx(1:min(no_state,no_level+restart,imax)))
                if (count==1.or.(count>1.and.imax>no_keep+second_iter)) then
                    open(unit=800,file=Nucleuss//file_code(J2/2,no_cA,no_cB)//'.tmp',&
                    form='unformatted',action='readwrite',access='direct',recl=no_state)
                    no_keep=no_conv+min(no_level+max(6,imax/4),imax/2,imax-no_conv-2)
!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(coefX,coefL,imax,eigenv,indx,no_keep)
!$OMP DO ORDERED SCHEDULE(DYNAMIC,1)
                    do nsc=1,no_keep
                        thrn=OMP_get_thread_num()+1
                        coefX(:,thrn)=real(0,rl)
                        do lsc = 1,imax
!$OMP  CRITICAL (read)          
                        read(700,rec=lsc) coefL
                        coefX(:,thrn)=coefX(:,thrn) +eigenv(lsc,indx(nsc))*coefL(:)
!$OMP END CRITICAL (read)          
                        end do
!$OMP ORDERED          
                        write(800,rec=nsc) coefX(:,thrn)
!$OMP END ORDERED          
                    end do
!$OMP END PARALLEL                    
                    coefL=coefR
                    matrix=real(0,rkw)
                    do nsc=1,no_keep
                         matrix(nsc,nsc)=alpha_copy(indx(nsc))
                         matrix(nsc,no_keep+1)=betat(imax)*eigenv(imax,indx(nsc))                       
                         matrix(no_keep+1,nsc)=betat(imax)*eigenv(imax,indx(nsc))
                    end do
                    alpha=real(0,rkw)
                    betat=real(0,rkw)
                    alpha(1:no_keep)=alpha_copy(indx(1:no_keep))
                    betat(1:no_keep)=matrix(1:no_keep,no_keep+1)
                    imin=no_keep+1
                    imax=imin+inc_iter
                    count=count+1       
                    close(unit=700)
                    close(unit=800)
                    nsc = SYSTEM('del '//Nucleuss//file_code(J2/2,no_cA,no_cB)//'.trl > '//Nucleuss//'.del')
                    NSC = SYSTEM('ren '//Nucleuss//file_code(J2/2,no_cA,no_cB)//'.tmp'//&
                                    ' '//Nucleuss//file_code(J2/2,no_cA,no_cB)//'.trl')
                    open(unit=700,file=Nucleuss//file_code(J2/2,no_cA,no_cB)//'.trl',&
                    form='unformatted',action='readwrite',access='direct',recl=no_state)
                    call lnzort(no_keep)
                    ovl=dot(coefL,coefL,no_state)
                    ovl=real(1,rkw)/sqrt(ovl)
                    coefR=ovl*coefL
                    call lnzort(-no_keep-1)
      		        isave=min(no_state,no_keep,i)
        	        open(unit=701,file=Nucleuss//file_code(J2/2,no_cA,no_cB)//'.abz',form='unformatted',action='write')
        	        write(701) isave,alpha(1:isave),betat(1:isave),matrix(1:isave,1:isave)
        	        write(701) no_keep,no_level,tt,no_state,basd,count,no_keep,no_conv,first_iter,&
                        second_iter,noiter
        	        write(701) bas
		            write(701) alpha_copy(:),betat(imax+1),indx(:)	
		            close(701)
                else
      		        isave=min(no_state,imax,i)
        	        open(unit=701,file=Nucleuss//file_code(J2/2,no_cA,no_cB)//'.abz',form='unformatted',action='write')
        	        write(701) isave,alpha(1:isave),betat(1:isave),matrix(1:isave,1:isave)
        	        write(701) imax,no_level,tt,no_state,basd,count,no_keep,no_conv,first_iter,&
                        second_iter,noiter
        	        write(701) bas
		            write(701) alpha_copy(:),betat(imax+1),indx(:)	
		            close(701)
                    imin=imax+1
                    imax=imin+inc_iter
                    if (imax +1> noiter) then
                        print *, ' Increase max_iter in LanczosParameters'
                        print *, ' Or SET ITERL=I where I >> ', noiter
                        imax=imin-1
                        go to 20
                    end if
                end if
                go to 10
            end if
        end do
20      isave=min(no_state,imax,i)
        print *
        open(unit=701,file=Nucleuss//file_code(J2/2,no_cA,no_cB)//'.abz',form='unformatted',action='write')
        write(701) isave,alpha(1:isave),betat(1:isave),matrix(1:isave,1:isave)
        write(701) imax,no_level,tt,no_state,basd,count,no_keep,no_conv,first_iter,&
                        second_iter,noiter
        write(701) bas
        write(701) alpha_copy(:),betat(isave+1),indx(:)
	    close(701)
        print '(5f10.3)', alpha_copy(indx(1:no_level))
        print *
        print *, ' No of iterations, northog',nti,ntog
   
   end subroutine Lnz_OP
   
   subroutine lnzort(i)
   implicit none
   integer,intent(in)::i
   integer::j,iol
   integer:: thrn,OMP_get_thread_num
   real(kind=rl) :: ovl
   
   if (i==0) then
        iol=no_state
        open(unit=700,file=Nucleuss//file_code(J2/2,no_cA,no_cB)//'.trl',&
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
   
end module TRL