!!DEC$ DEFINE matrix
!!DEC$ DEFINE nulnzp
!DEC$ define nutrlp
module Globals


use Parameters
use ProjParameters
use MatrixParameters
use LanczosParameters
use OperParameters
use ClebschGordan
use InputNuShell
use TrdChar
use Files
use TriDiag
use SOrder
use DotProduct
use SubM
use LnzParameters
use XParameters
use Extra
use Mscheme


implicit none

type(tablear),dimension(:,:),allocatable:: tablebh
integer(kind=2),dimension(:,:),allocatable:: tablebhn
real(kind=rkw),dimension(:),allocatable:: coefR,trdav,trdbv,coefT,coefL,racahs,amtv,bmtv
real(kind=rkw),dimension(:,:),allocatable:: BUFFER1,BUFFER2,BUFFER3
real(kind=rkw),dimension(:,:),allocatable:: matrix,meV,coefX
character(len=8),dimension(:),allocatable::trda,trdb,rac,bas,opn
real(kind=8),dimension(:),allocatable:: opmeV
integer,dimension(:),allocatable:: posa,posb,posnv,dimnv,posma,posmb,dima,dimb
integer,dimension(:,:),allocatable:: nJTa,nJTb,trunc_bas
integer,dimension(:,:),allocatable:: posal,posbl,locb,blocb
type(la),dimension(:),allocatable:: nplra,nplrb,nplrma,nplrmb
type(spartition),dimension(:),allocatable:: sparta,spartb
integer,dimension(ptdim):: cNN
character(len=7)::opfile1,opfile2,intfile
character(len=6)::opfile16,opfile26,Nucleuss
character(len=1):: xtraa,xtrab
character(len=8):: ta,tb,pb,tp,bm,ba,rn,z1,op
character(len=4):: tropa,tropb
integer:: min_cJ2a,max_cJ2a,no_sparta,no_spart,no_sgmta,maxda,del_cJ2a
integer:: min_cJ2b,max_cJ2b,no_spartb,no_sgmtb,no_state=0,maxdb,no_bas,del_cJ2b
integer:: trdda,trddb,prddb,racd,trddat,trddbt,opd,no_baso,no_trunc,no_trunco
integer:: max_lma,max_lmr,max_lm,basd,max_lmb,max_lmo,max_lambda,del_lm
integer:: ajl,ajr,anpl,anpr,ashl,ashr,lm,idx,lmx,lmn,lmnn,lmxx
integer:: bjl,bjr,bnpl,bnpr,ajrn,ajrx,ajln,ajlx,Tparitya,Tparityb
integer:: bjln,bjlx,bshl,bshr,bjrn,bjrx,no_level,nwordsa,nwordsb
integer:: lsc,psc,atsc,btsc,bsc,nsc,hsc,basp,ptsc,rsc,restart
integer:: J2,J2u,nthr=2,ip,mna,mnb,max_dimv,av_dimv,Ntot
integer:: bseg,aseg,dimptl,dimptr,matdima,matdimb,no_cA,no_cB,kk
integer:: max_cJ2s,min_cJ2s,nmps,nep,max_spart,OMP_get_max_threads
type(lv),dimension(:),allocatable:: lvac,lvbc,lvar,lvbr
type(la),dimension(:),allocatable:: nplrao,nplrbo,nplrbop
type(lv),dimension(:),allocatable:: lbl
integer,dimension(:),allocatable:: dal,dar,dbl
logical:: done,new,proton,neutron,azero,bzero
real(kind=rkw):: me,renorm,rtbas
real(kind=rkw):: phase,rcv
real:: xa,xb,xxa,xxb

end module Globals

module XMultiply

use Globals

implicit none

contains

   subroutine multiplya()
    implicit none
    integer:: thrn,OMP_get_thread_num
    logical:: diag,trans
    if (azero) return    
!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(bas,basd,amtv,coefL,posnv,nJTa,nJTb,&
!$OMP no_sparta,no_spartb,posma,nplrma,coefR,truncate,bas_log,&
!$OMP BUFFER1,BUFFER3)
!$OMP DO SCHEDULE(DYNAMIC,1)   
    do bsc=1,basd
       if (truncate) then
          if (.not.bas_log(bsc)) go to 10
       end if
       ajl=ichar(bas(bsc)(1:1))
       bjl=ichar(bas(bsc)(2:2))
       thrn=OMP_get_thread_num()+1
       trans=.true.
       call mulmtx(coefL,coefR,posnv(bsc),posma(ajl/2),amtv, &
                      no_sparta,no_spartb,nJTa(:,ajl/2),nJTb(:,bjl/2),nplrma(ajl/2),&
                      BUFFER1(:,thrn),BUFFER3(:,thrn),trans,ajl,bjl)
10  end do
!$OMP END PARALLEL

   end subroutine multiplya
   
   subroutine multiplyb()
    implicit none
    integer:: thrn,OMP_get_thread_num
    logical:: diag,trans
    
    if (bzero) return  
!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(bas,basd,bmtv,coefL,posnv,nJTa,nJTb,&
!$OMP no_sparta,no_spartb,posmb,nplrmb,coefR,truncate,bas_log,&
!$OMP BUFFER1,BUFFER3)
!$OMP DO SCHEDULE(DYNAMIC,1)   
    do bsc=1,basd
       if (truncate) then
          if (.not.bas_log(bsc)) go to 10
       end if
       ajl=ichar(bas(bsc)(1:1))
       bjl=ichar(bas(bsc)(2:2))
       thrn=OMP_get_thread_num()+1
       trans=.false.
       call mulmtx(coefL,coefR,posnv(bsc),posmb(bjl/2),bmtv, &
                      no_sparta,no_spartb,nJTa(:,ajl/2),nJTb(:,bjl/2),nplrmb(bjl/2),&
                      BUFFER1(:,thrn),BUFFER3(:,thrn),trans,ajl,bjl)
10    end do
!$OMP END PARALLEL

   end subroutine multiplyb
   
  
   subroutine multiplyx()
    implicit none
    integer:: thrn,OMP_get_thread_num
    logical:: nzero,zero
    integer:: bjrs,lscs,ntr,loop,delj2,lp,dj2
    
    coefX=real(0,rkw)  
    delj2=(max_cJ2a-min_cJ2a)/del_cJ2a+1
    loop=basd*delj2
    thrn=1
    if (azero.or.bzero) then
        coefL=real(0,rkw) 
        return
    end if
!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(bas,racahs,tablebh,coefR,coefX,basd,truncate,&
!$OMP trdav,posa,nJTa,nJTb,lvac,lvbc,lvar,lvbr,nplra,nplrb,BUFFER1,BUFFER2,dimnv,bas_log,&
!$OMP no_sparta,no_spartb,dimb,posnv,trdbv,posb,nplrbop,nplrbo,nplrao,locb,del_cJ2b,&
!$OMP BUFFER3,blocb,tablebhn,dal,dar,dbl,mev,min_cJ2a,max_cJ2a,loop,delj2,del_cJ2a)
!$OMP DO SCHEDULE(DYNAMIC,1)
    do lp=1,loop
    bsc=(lp-1)/delj2+1
    thrn=OMP_get_thread_num()+1
    ajl=ichar(bas(bsc)(1:1))
    bjl=ichar(bas(bsc)(2:2))
    dj2=lp-(bsc-1)*delj2
    ajr=min_cJ2a-del_cJ2a+del_cJ2a*dj2
       hsc=ajr/2
       if (truncate) then
        if (.not.bas_log(bsc)) go to 75
       end if
       if (tablebhn(bsc,hsc)==0) go to 75
       do nsc=1,tablebhn(bsc,hsc)
            if (tablebh(bsc,hsc)%bme(nsc)%bn==0) go to 50
            atsc=tablebh(bsc,hsc)%a(nsc)
            call crproj(nplra(atsc),lvac(thrn),lvar(thrn),no_sparta,no_sparta)
            call zerola(nplrbo(thrn),no_spartb)
            call zerola(nplrbop(thrn),no_spartb)
            do psc=1,tablebh(bsc,hsc)%bme(nsc)%bn
                btsc=tablebh(bsc,hsc)%bme(nsc)%b(psc)
                call laor(nplrbo(thrn),nplrb(btsc),no_spartb,no_spartb)
            end do            
            call crproj(nplrbo(thrn),lvbc(thrn),lvbr(thrn),no_spartb,no_spartb)
            call condim(no_sparta,no_spartb,blocb(:,thrn),&
                    nJTa(:,ajl/2),nJTa(:,ajr/2),nJTb(:,bjl/2),&
                            dar(thrn),dal(thrn),dbl(thrn),&
                                   lvac(thrn),lvar(thrn),lvbr(thrn))
            lscs=tablebh(bsc,hsc)%bme(nsc)%h(1)
            bjrs=ichar(bas(lscs)(2:2))
            zero=.true.
            ntr=0
            do psc=1,tablebh(bsc,hsc)%bme(nsc)%bn
                btsc=tablebh(bsc,hsc)%bme(nsc)%b(psc)
                lsc=tablebh(bsc,hsc)%bme(nsc)%h(psc)
                bjr=ichar(bas(lsc)(2:2))
                me=tablebh(bsc,hsc)%bme(nsc)%me(psc)
                if (lscs/=lsc) then
                call crproj(nplrbop(thrn),lvbc(thrn),lvbr(thrn),no_spartb,no_spartb)
                call single_productm(coefR,posnv(lscs),mev(:,thrn),no_sparta,no_spartb,&
                     nJTa(:,ajr/2),nJTb(:,bjl/2),nJTb(:,bjrs/2),&
                     nplrb,trdbv,posb,blocb(:,thrn),dar(thrn),dal(thrn),&
                     dbl(thrn),lvac(thrn),lvbc(thrn),lvar(thrn),lvbr(thrn),&
                     zero,BUFFER1(:,thrn),BUFFER2(:,thrn),ntr,locb(:,thrn),thrn)
                ntr=0
                call zerola(nplrbop(thrn),no_spartb)
                bjrs=bjr
                lscs=lsc
                zero=.false.
                end if
                if (me/=0.0) then
                ntr=ntr+1
                locb(ntr,thrn)=btsc
                mev(ntr,thrn)=me
                call laor(nplrbop(thrn),nplrb(btsc),no_spartb,no_spartb)
                end if
            end do            
            if (ntr/=0) then
                call crproj(nplrbop(thrn),lvbc(thrn),lvbr(thrn),no_spartb,no_spartb)
                call single_productm(coefR,posnv(lscs),mev(:,thrn),no_sparta,no_spartb,&
                     nJTa(:,ajr/2),nJTb(:,bjl/2),nJTb(:,bjrs/2),&
                     nplrb,trdbv,posb,blocb(:,thrn),dar(thrn),dal(thrn),&
                     dbl(thrn),lvac(thrn),lvbc(thrn),lvar(thrn),lvbr(thrn),&
                     zero,BUFFER1(:,thrn),BUFFER2(:,thrn),ntr,locb(:,thrn),thrn)
            end if
            call crproj(nplrbo(thrn),lvbc(thrn),lvbr(thrn),no_spartb,no_spartb)
            call single_product1(coefX(:,thrn),posnv(bsc),no_sparta,no_spartb,nJTa(:,ajl/2),&
                 nJTa(:,ajr/2),nJTb(:,bjl/2),nplra(atsc),&
                 trdav,posa(atsc),dar(thrn),dal(thrn),dbl(thrn),&
                 lvac(thrn),lvar(thrn),lvbr(thrn),BUFFER1(:,thrn),&
                 BUFFER2(:,thrn),BUFFER3(:,thrn),thrn)
50     end do 
75  end do
    
!$OMP END PARALLEL
    coefL(:)=coefX(:,1)
    do lp=2,nthr
    coefL(:)=coefL(:)+coefX(:,lp)
    end do
    
   end subroutine multiplyx

end module XMultiply

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
        integer:: n,i,imax,imin,ix,isave,np,pd,tqlim,status,idum,idum1,npre_iter=pre_iter
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
        if (Tparity==0.and.no_level>no_bas+no_baso) no_level=no_bas+no_baso
        if (Tparity==1.and.no_level>no_bas) no_level=no_bas  
        if (Tparity==-1.and.no_level>no_baso) no_level=no_baso        
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
                 form='unformatted',action='readwrite',access='direct',recl=no_state*rkw/4)
 	         read(700,rec=imax-1) coefT
 	         read(700,rec=imax) coefR      
             imin=imax+1
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
        if (truncate) then
            where(.not.trunc) coefL=0.0
        end if
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
            if (truncate) then
            where(.not.trunc) coefL=0.0
            end if
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
            if (truncate) then
                where(.not.trunc) coefL=0.0
            end if
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
            if (Tparity==0.and.i==no_bas+no_baso) go to 15        
            if (Tparity==1.and.i==no_bas) go to 15        
            if (Tparity==-1.and.i==no_baso) go to 15        
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
                    if (truncate) then
                    do n=1,no_state
                        if (.not.trunc(n)) coefL(n)=0.0
                    end do
                    end if
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
                    if (truncate) then
                    do n=1,no_state
                        if (.not.trunc(n)) coefL(n)=0.0
                    end do
                    end if
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
            beta_copy(1)=real(0,rkw)
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
        call indexx(noiter,alpha_copy,indx)
        if (Tparity==0) ix=no_bas+no_baso
        if (Tparity==1) ix=no_bas
        if (Tparity==-1) ix=no_baso
        if (imax== ix.or.abs(ovl) < restart_limit) go to 20
        do ix=1,min(no_state,no_level+restart,imax)
            if (abs(alpha_save(ix)-alpha_copy(indx(ix))) > level_limit) then
                no_conv=ix-1 
                alpha_save(1:min(no_state,no_level+restart,imax))=&
                    alpha_copy(indx(1:min(no_state,no_level+restart,imax)))
                if (count==1.or.(count>1.and.imax>no_keep+second_iter)) then
                    open(unit=800,file=Nucleuss//file_code(J2/2,no_cA,no_cB)//'.tmp',&
                    form='unformatted',action='readwrite',access='direct',recl=no_state*rkw/4)
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
                    form='unformatted',action='readwrite',access='direct',recl=no_state*rkw/4)
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
	    do i=1,no_level
	    if (indx(i)==0) no_level=no_level-1
	    end do
       ! print '(5f10.3)', alpha_copy(indx(1:no_level))
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
        iol=no_state*rkw/4
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

module Lnzs

use Globals
use XMultiply

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
        integer:: n,i,imax,imin,ix,isave,np,pd,tqlim,status,idum,idum1,npre_iter=pre_iter
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
        if (Tparity==0.and.no_level>no_bas+no_baso) no_level=no_bas+no_baso
        if (Tparity==1.and.no_level>no_bas) no_level=no_bas  
        if (Tparity==-1.and.no_level>no_baso) no_level=no_baso        
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
        if (truncate) then
            where (.not.trunc) coefL=0.0
        end if
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
            if (truncate) then
                where (.not.trunc) coefL=0.0
            end if
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
            if (truncate) then
                where (.not.trunc) coefL=0.0
            end if
            if (i /= 1) coefL=coefL-beta(i)*coefT
            coefT=coefR
            alpha(i)=dot(coefR,coefL,no_state)
            if (i == no_state) go to 55
            if (Tparity==0.and.i==no_bas+no_baso) go to 55        
            if (Tparity==1.and.i==no_bas) go to 55        
            if (Tparity==-1.and.i==no_baso) go to 55        
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
                if (truncate) then
                do n=1,no_state
                    if (.not.trunc(n)) coefL(n)=0.0
                end do
                end if
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
        call indexx(noiter,alpha_copy,indx)
        if (Tparity==0) ix=no_bas+no_baso
        if (Tparity==1) ix=no_bas
        if (Tparity==-1) ix=no_baso
        if (imax== ix.or.abs(ovl) < restart_limit) go to 20
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
	    !do i=1,no_level
	    !if (indx(i)==0) no_level=no_level-1
	    !end do
        !print '(5f10.3)', alpha_copy(indx(1:no_level))
        !print *
        print *, ' No of iterations, northog',imax,ntog
   
   end subroutine Lnz_OP
   
   subroutine lnzort(i)
   implicit none
   integer,intent(in)::i
   integer::j,iol
   integer:: thrn,OMP_get_thread_num
   real(kind=rl) :: ovl
   
   if (i==0) then
        iol=no_state*rkw/4
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

module Lnz
    use Globals
!DEC$ IF DEFINED (nulnzp)
    use Lnzs
!DEC$ ENDIF

!DEC$ IF DEFINED (nutrlp)
    use TRL
!DEC$ ENDIF
  
contains

    subroutine Lnz_TRD(base)
    implicit none
    integer,intent(in)::base
    integer:: no_Rac=0,OMP_get_max_threads
    real::t1,t2
    nthr=OMP_get_max_threads()
    if (min_MSon/=0.or.max_MSon/=999) truncate=.true.
    cNN=0
    do nsc=1,no_shells
        cNN(nsc)=2*sn(nsc)+sl(nsc)
    end do
    allocate(locb(max_trdb,nthr))
    allocate(blocb(max_trdb,nthr))
    allocate(meV(max_trdb,nthr))
    allocate(lvac(nthr))
    allocate(lvbc(nthr))
    allocate(lvar(nthr))
    allocate(lvbr(nthr))
    allocate(nplrao(nthr))
    allocate(nplrbo(nthr))
    allocate(nplrbop(nthr))
    allocate(lbl(nthr))
    allocate(dar(nthr))
    allocate(dal(nthr))
    allocate(dbl(nthr))
    locb=0
    blocb=0
    meV=real(0,rkw)
    z1=repeat(char(0),8)
    print *, ' VERSION 4.0 '
    
    read(base+1) no_sparta,min_cJ2a,max_cJ2a,proton,neutron,xa,Tparitya,tropa,del_cJ2a
    if (proton) typea='p'
    if (neutron) typea='n'
    allocate(nJTa(no_sparta,min_cJ2a/2:max_cJ2a/2))
    allocate(sparta(no_sparta))
    read(base+1) nJTa(1:no_sparta,min_cJ2a/2:max_cJ2a/2)
    read(base+1) sparta(1:no_sparta)
    
    allocate(para(no_sparta))
    do anpl=1,no_sparta
        para(anpl)=cParity(sparta(anpl))
    end do

    do nsc=1,nthr
    call allocatelv(lvac(nsc),no_sparta)    
    call allocatelv(lvar(nsc),no_sparta)    
    end do

    read(base+2) no_spartb,min_cJ2b,max_cJ2b,proton,neutron,xb,Tparityb,tropb,del_cJ2b
    if (proton) typeb='p'
    if (neutron) typeb='n'
    if (typea==typeb) then
        print *,' Warning nucleon types are identical '
        typea='n'
        typeb='p'
    end if
    allocate(nJTb(no_spartb,min_cJ2b/2:max_cJ2b/2))
    allocate(spartb(no_spartb))
    read(base+2) nJTb(1:no_spartb,min_cJ2b/2:max_cJ2b/2)
    read(base+2) spartb(1:no_spartb)

    if (truncate) allocate(trunc_bas(no_sparta,no_spartb))
    allocate(parb(no_spartb))
    do bnpl=1,no_spartb
        parb(bnpl)=cParity(spartb(bnpl))
    end do
    
    azero=.false.
    nsc=0
    do ajl=min_cJ2a,max_cJ2a,2
    nsc=nsc+dna(no_sparta,ajl)
    end do
    if (nsc==0) azero=.true.
    
    bzero=.false.
    nsc=0
    do bjl=min_cJ2b,max_cJ2b,2
    nsc=nsc+dna(no_spartb,bjl)
    end do
    if (nsc==0) bzero=.true.
    
    do nsc=1,nthr
    call allocatelv(lvbc(nsc),no_spartb)    
    call allocatelv(lvbr(nsc),no_spartb)    
    call allocatela(nplrbo(nsc),no_spartb,no_spartb)    
    call allocatela(nplrbop(nsc),no_spartb,no_spartb)
    call allocatela(nplrao(nsc),no_sparta,no_sparta)    
    end do

    allocate(posal(no_sparta,nthr))
    allocate(posbl(no_spartb,nthr))
    do nsc=1,nthr
    call allocatelv(lbl(nsc),no_spartb)
    end do
    
    call set_ab(no_sparta,min_cJ2a,max_cJ2a,del_cJ2a,&
                 no_spartb,min_cJ2b,max_cJ2b,del_cJ2b,nJTa,nJTb)
! setup nucleon numbers    

    nmps=nmp
    if (typea=='p') then
    nmp=(no_cN+nmps)/2
    no_cB=nmp
    if (no_cB==0) bzero=.true.
    nmp=(no_cN-nmps)/2
    no_cA=nmp
    if (no_cA==0) azero=.true.
    else
    nmp=(no_cN-nmps)/2
    no_cB=nmp
    if (no_cB==0) bzero=.true.
    nmp=(no_cN+nmps)/2
    no_cA=nmp
    if (no_cA==0) azero=.true.
    end if

! end set up nucleon numbers
        
    J2=min_cJ2
    
    allocate(string_order(nstring))
    allocate(node(max_node))
    nonode=max_node

!setup basis/

    call reset_string_order()
    no_state=0
    basd=0
    do ajl=min_cJ2a,max_cJ2a,del_cJ2a
        bjln=max(abs(ajl-J2),min_cJ2b)
        bjlx=min(abs(ajl+J2),max_cJ2b)
        bjln=bjln+mod(bjln-min_cJ2b,del_cJ2b)
        do bjl=bjln,bjlx,del_cJ2b
            bm(1:1)=char(ajl)
            bm(2:2)=char(bjl)
            bm(3:8)=z1(3:8)
            call add_string(bm(1:8))
            basd=basd+1
        end do
    end do
    
    print *, ' No of matrices BAS ',basd
    
    allocate(bas(basd))
    allocate(dimnv(basd))
    allocate(posnv(basd))
    if (truncate) allocate(bas_log(basd))
    dimnv=0
    posnv=0
    dimnv=0
    done=.true.
    idx=0
    nsc=0
10  idx=get_next_string(done)
    nsc=nsc+1
    bas(nsc)=string_order(idx)%string(1:8)
    if (done) go to 20
    go to 10
20  if (basd/=nsc) then
        print *, ' basd,nsc ',basd,nsc
        stop '  basd /= nsc  '
    end if
    
    no_bas=0
    max_dimv=0
    no_state=0
    do bsc=1,basd
    ba=bas(bsc)(1:4)
    ajl=ichar(ba(1:1))
    bjl=ichar(ba(2:2))
    dimnv(bsc)= dna(no_sparta,ajl)*dnb(no_spartb,bjl)
    do anpl=1,no_sparta
    do bnpl=1,no_spartb
        if (para(anpl)*parb(bnpl)==1) no_bas=no_bas+nJTa(anpl,ajl/2)*nJTb(bnpl,bjl/2)
    end do
    end do
    posnv(bsc)=no_state
    max_dimv=max(max_dimv,dimnv(bsc))
    no_state=no_state+dimnv(bsc)                
    end do
    
    nsc=0
    no_trunc=0
    no_trunco=0
    lsc=basd
    allocate(lcoef(no_state))
    if (truncate) allocate(trunc(no_state))
    lcoef=.false.
    trunc=.false.
    do bsc=1,basd
    ba=bas(bsc)(1:4)
    ajl=ichar(ba(1:1))
    bjl=ichar(ba(2:2))
    if (truncate) then
    trunc_bas=0
    bas_log(bsc)=.true.
    end if
    do bnpl=1,no_spartb
        do btsc=1,nJTb(bnpl,bjl/2)
            do anpl=1,no_sparta
                if (truncate) then 
                    Ntot=psum(cNN*(sparta(anpl)+spartb(bnpl)))
                    if (Ntot>=min_MSon.and.Ntot<=max_MSon) trunc_bas(anpl,bnpl)=1
                end if    
                do atsc=1,nJTa(anpl,ajl/2)
                    nsc=nsc+1  
                    if (Tparity==0) then
                        lcoef(nsc)=.true.
                    else if (para(anpl)*parb(bnpl)==Tparity) then
                        lcoef(nsc)=.true.
                    end if
                    if (truncate) then
                        if (Ntot>=min_MSon.and.Ntot<=max_MSon) then
                            trunc(nsc)=.true..and.lcoef(nsc)
                            if (lcoef(nsc)) no_trunc=no_trunc+1
                            if (.not.lcoef(nsc)) no_trunco=no_trunco+1
                        end if
                    end if
                end do
            end do
       end do
    end do
    if (truncate) then
        if (sum(trunc_bas)==0) then
            bas_log(bsc)=.false.
            lsc=lsc-1
        end if
    end if
    end do
    if (truncate) print *, ' No of truncated matrices BAS ',lsc
    
    av_dimv=real(no_state)/real(basd)
    
    no_baso=no_state-no_bas
    if (truncate) then
        if (Tparity==1) then
        no_bas=no_trunc
        no_baso=no_trunco
        else if (Tparity==-1) then
        no_bas=no_trunco
        no_baso=no_trunc
        end if        
    end if 
    print *, ' No of total states ',no_bas+no_baso
    if (Tparity/=0) then
    print *, ' No of basis states with + parity ',no_bas    
    print *, ' No of basis states with - parity ',no_baso    
    end if
    if (truncate) then
        if (Tparity==1) then
        no_bas=no_trunc
        no_baso=no_trunco
        else if (Tparity==-1) then
        no_bas=no_trunco
        no_baso=no_trunc
        end if        
    end if 
    
    
! end setup basis 
    if (azero.or.bzero) go to 45

! setup racahs   
    racd=0
    call reset_string_order()
    do ajl=min_cJ2a,max_cJ2a,del_cJ2a
      rn(1:1)=char(ajl)
      do bjl=min_cJ2b,max_cJ2b,del_cJ2b
        rn(2:2)=char(bjl)
        do ajr=min_cJ2a,max_cJ2a,del_cJ2a
          rn(3:3)=char(ajr)
          do bjr=min_cJ2b,max_cJ2b,del_cJ2b
            rn(4:4)=char(bjr)
            do lm=0,max_sj2u*2,2
            rn(5:5)=char(lm)
            rn(6:8)=z1(6:8)
            rcv=racah(ajl,bjl,ajr,bjr,J2,lm)
            if (abs(rcv)>0.00000) then
                call add_string(rn(1:8))
                racd=racd+1
            end if
            end do
          end do
        end do
      end do
    end do

    if (racd/=0) then
    allocate(rac(racd))
    allocate(racahs(racd))
    racahs=real(0,rkw)
    done=.true.
    idx=0
    nsc=0
30  idx=get_next_string(done)
    nsc=nsc+1
    rac(nsc)=string_order(idx)%string(1:8)
    if (done) go to 40
    go to 30
40  if (racd/=nsc) then
        print *, ' racd,nsc ',racd,nsc
        stop '  racd /= nsc  '
    end if
    end if
 
    do nsc=1,racd
        rn=rac(nsc)(1:5)
        ajl=ichar(rn(1:1))
        bjl=ichar(rn(2:2))
        ajr=ichar(rn(3:3))
        bjr=ichar(rn(4:4))
        lm=ichar(rn(5:5))
        phase=real(1,rkw)
        if (btest((ajl+bjr-J2-lm)/2,0)) phase=-phase
        rcv=phase*racah(ajl,bjl,ajr,bjr,J2,lm)*sqrt(real(lm+1,8))
        racahs(nsc)=rcv
        no_Rac=no_Rac+1
    end do
! end setup racahs

    deallocate(string_order)
    deallocate(node)
    
! read opn         
    read(base+5) opd,max_lmo
    allocate(opn(opd))
    allocate(opmeV(opd))
    if (opd/=0) then
    do psc=1,opd
    read(base+5) opn(psc),opmeV(psc)
    if (abs(opmeV(psc))<0.0000000001) opmeV(psc)=real(0,rkw)
    opmeV(psc)=renorm*opmeV(psc)
    end do
    end if
! end read opn
! read trda
    nwordsa=((no_sparta-1)/32+1)
    read(base+6) trdda,max_lma,trddat
    allocate(trda(trdda))
    allocate(posa(trdda))
    allocate(dima(trdda))
    allocate(nplra(trdda))
    posa=0
    dima=0
    allocate(trdav(trddat))
    if (trdda/=0) then
    read(base+6) trda
    read(base+6) dima
    nsc=0
    do atsc=1,trdda
    posa(atsc)=nsc
    if (dima(atsc)/=0) then
    read(base+6) trdav(nsc+1:nsc+dima(atsc))
    do lsc=nsc+1,nsc+dima(atsc)
        if (abs(trdav(lsc))<.00000000001) trdav(lsc)=0.0
    end do
    call allocatela(nplra(atsc),no_sparta,no_sparta)
    read(base+6) (nplra(atsc)%lvect(anpl)%lword(1:nwordsa),anpl=1,no_sparta)
    nsc=nsc+dima(atsc)
    else
    call allocatela(nplra(atsc),no_sparta,no_sparta)
    end if
    end do
    end if
! end read trda 
    deallocate(dima)

! read trdb
    nwordsb=((no_spartb-1)/32+1)
    read(base+3) trddb,max_lmb,trddbt
    allocate(trdb(trddb))
    allocate(posb(trddb))
    allocate(dimb(trddb))
    allocate(nplrb(trddb))
    posb=0
    dimb=0
    allocate(trdbv(trddbt))
    if (trddb/=0) then
    read(base+3) trdb
    read(base+3) dimb
    nsc=0
    do btsc=1,trddb
    posb(btsc)=nsc
    if (dimb(btsc)/=0) then
    read(base+3) trdbv(nsc+1:nsc+dimb(btsc))
    do lsc=nsc+1,nsc+dimb(btsc)
        if (abs(trdbv(lsc))<.0000000001) trdbv(lsc)=0.0
    end do
    call allocatela(nplrb(btsc),no_spartb,no_spartb)
    read(base+3) (nplrb(btsc)%lvect(bnpl)%lword(1:nwordsb),bnpl=1,no_spartb)
    nsc=nsc+dimb(btsc)
    else
    call allocatela(nplrb(btsc),no_spartb,no_spartb)
    endif
    end do
    end if
! end read trdb 
    deallocate(dimb)

! end read prdb 
45  Nucleuss=Nucleus
    max_cJ2s=max_cJ2
    min_cJ2s=min_cJ2

    xtraa=intfile(7:7)
    xtrab=intfile(7:7)
    
! start read nn and pp matrices
! matrices B
    opfile2=adjustr(opfile2)
    opfile26(1:6)=opfile2(2:7)
    Nucleus=opfile26
!    xtrab='b'
    max_cJ2=max_cJ2b
    min_cJ2=min_cJ2b
    if (typeb=='n') then
        nmp=no_cB
    else
        nmp=-no_cB
    end if
    call open_files_read(no_spart,'.mat',400,0,-1,xtrab)
    call open_files_read(no_spart,'.sgt',500,0,-1,xtrab)
    kk=0
    allocate(posmb(min_cJ2b/2:max_cJ2b/2))
    allocate(nplrmb(min_cJ2b/2:max_cJ2b/2))
    do bjl=min_cJ2b/2,max_cJ2b/2
        call allocatela(nplrmb(bjl),no_spartb,no_spartb)
    end do
    posmb=0
    basp=0
    do bjl=min_cJ2b,max_cJ2b,2
       read (500+kk,err=50,end=50) no_sgmtb,matdimb,xxb
       if (xxb/=xb) stop ' B mat file and B proj file not consistent '
       do bseg=1,no_sgmtb
          read(400+kk) dimptl,dimptr
          basp=basp+dimptl*dimptr
       end do
       rewind(500+kk)
       rewind(400+kk)
       read(500+kk)
       read(400+kk)
       kk=kk+1
50  end do
    allocate(bmtv(basp))
    basp=0
    kk=0
    do bjl=min_cJ2b,max_cJ2b,2
       read (500+kk,err=55,end=55) no_sgmtb,matdimb,xxb
       posmb(bjl/2)=basp  
       do bseg=1,no_sgmtb
          read(400+kk) dimptl,dimptr,bnpl,bnpr,bmtv(basp+1:basp+dimptl*dimptr)
          nplrmb(bjl/2)%lvect(bnpl)=bnpr.or.nplrmb(bjl/2)%lvect(bnpl)
          basp=basp+dimptl*dimptr
       end do
       kk=kk+1
55  end do
    
    call close_files_read(400,0)
    call close_files_read(500,0)
! Matrices A
    
    opfile1=adjustr(opfile1)
    opfile16(1:6)=opfile1(2:7)
    Nucleus=opfile16
!    xtraa='a'
    max_cJ2=max_cJ2a
    min_cJ2=min_cJ2a
    if (typea=='p') then
        nmp=-no_cA
    else
        nmp=no_cA
    end if
    call open_files_read(no_spart,'.mat',400,0,-1,xtraa)
    call open_files_read(no_spart,'.sgt',500,0,-1,xtraa)
    kk=0
    allocate(posma(min_cJ2a/2:max_cJ2a/2))
    allocate(nplrma(min_cJ2a/2:max_cJ2a/2))
    do ajl=min_cJ2a/2,max_cJ2a/2
        call allocatela(nplrma(ajl),no_sparta,no_sparta)
    end do
    posma=0
    basp=0
    do ajl=min_cJ2a,max_cJ2a,2
       read (500+kk,err=60,end=60) no_sgmta,matdima,xxa
       if (xxa/=xa) stop ' a mat file and a proj file not consistent '
       do aseg=1,no_sgmta
          read(400+kk) dimptl,dimptr
          basp=basp+dimptl*dimptr
       end do
       rewind(400+kk)
       rewind(500+kk)
       read(400+kk)
       read(500+kk)
       kk=kk+1
60  end do
    allocate(amtv(basp))
    kk=0
    basp=0
    do ajl=min_cJ2a,max_cJ2a,2
       read (500+kk,err=66,end=66) no_sgmta,matdima,xxa
       posma(ajl/2)=basp  
       do aseg=1,no_sgmta
          read(400+kk) dimptl,dimptr,anpl,anpr,amtv(basp+1:basp+dimptl*dimptr)
          nplrma(ajl/2)%lvect(anpl)=anpr.or.nplrma(ajl/2)%lvect(anpl)
          basp=basp+dimptl*dimptr
       end do
       kk=kk+1
66  end do
   
    call close_files_read(400,0)
    call close_files_read(500,0)
! end setup matrices
! get max_lm
    max_lambda=del_cJ2
    if (max_lambda<=2.or.btest(max_lambda,0)) max_lambda=2*max_sj2u
    del_lm=del_cT2
    if (del_lm<2.or.btest(del_lm,0)) del_lm=2
    max_lm=min(max_lma,max_lmo,2*max_sj2u,max_lambda)
    call drdbh()   
    maxda=0 
    do ajl=min_cJ2a,max_cJ2a,del_cJ2a
    maxda=max(maxda,dna(no_sparta,ajl))
    end do
    maxdb=0 
    do bjl=min_cJ2b,max_cJ2b,del_cJ2b
    maxdb=max(maxdb,dnb(no_spartb,bjl))
    end do
    max_dimv=maxda*maxdb

    allocate(BUFFER1(max_dimv,nthr))
    allocate(BUFFER2(max_dimv,nthr))
    allocate(BUFFER3(max_dimv,nthr))
    if (truncate) allocate(lbuffer1(max_dimv,nthr))
    if (truncate) allocate(lbuffer2(max_dimv,nthr))
    BUFFER1=real(0,rkw)
    BUFFER2=real(0,rkw)
    BUFFER3=real(0,rkw)
                
! allocate and set vectors    
    allocate(coefL(no_state))
    allocate(coefX(no_state,nthr))
    allocate(coefR(no_state))
    allocate(coefT(no_state))
    coefL=real(0,rkw)
    coefR=real(0,rkw)
    coefT=real(0,rkw)
    coefX=real(0,rkw)
! end setup vectors
!DEC$ IF DEFINED (matrix)

    call matrixs()

!DEC$ ELSE    
! start lanczos iteration 
    call Lnz_OP(real(rtbas))
    
!DEC$ ENDIF
                          
   end subroutine Lnz_TRD
!DEC$ IF DEFINED (matrix)
   
   subroutine matrixs()
   implicit none
   real(kind=rkw),dimension(:,:),allocatable::matrixc,matrix1,matrix2,matrix3
   real(kind=rkw),dimension(:),allocatable::vector
   integer:: i,j
   
   allocate(matrixc(no_state,no_state))
   matrix=real(0,rkw)
   write(33,*) no_state
   write(33,*)
   write(34,*) no_state
   write(34,*)
   do i=1,no_state
   coefR=real(0,rkw)
   coefR(i)=real(1,rkw)
   call multiplyx()
   call multiplya()
   call multiplyb()
   matrixc(i,:)=coefL(:) 
   end do
   do i=1,no_state
   write(33,*) matrixc(1:no_state,i)
   write(33,*)
   end do
   matrixc=transpose(matrixc)
   do i=1,no_state
   write(34,*) matrixc(1:no_state,i)
   write(34,*)
   end do
   
   
   end subroutine matrixs
   
!DEC$ ENDIF
   
   pure function dna(n,ja)
   implicit none
   integer,intent(in):: n,ja
   integer:: dna
   if (n==0) then
        dna=0
        return
   end if
   dna=sum(nJTa(1:n,ja/2))
   end function dna
   
   pure function dnb(n,jb)
   implicit none
   integer,intent(in):: n,jb
   integer:: dnb
   if (n==0) then
        dnb=0
        return
   end if
   dnb=sum(nJTb(1:n,jb/2))
   end function dnb
   
          
    subroutine drdbh()
    implicit none
    integer:: min_ashells,max_ashells,store,nn,sts,min_bshells,max_bshells,bn,tn
    integer:: nja2,xja2,njb2,xjb2,dja2,djb2
    if (azero.or.bzero) return
    min_ashells=min_shellsa()
    max_ashells=max_shellsa()
    min_bshells=min_shellsb()
    max_bshells=max_shellsb()
    nja2=min_cJ2a/2
    xja2=max_cJ2a/2
    njb2=min_cJ2b/2
    xjb2=max_cJ2b/2
    dja2=del_cJ2a/2
    djb2=del_cJ2b/2
    allocate(tablebh(basd,nja2:xja2))
    allocate(tablebhn(basd,nja2:xja2))
    tablebhn=0
    
    do store=0,2
    if (store==1) then
        do  bsc=1,basd
        do  hsc=nja2,xja2,dja2
           if (tablebhn(bsc,hsc)/=0) then
           allocate(tablebh(bsc,hsc)%a(tablebhn(bsc,hsc)))
           allocate(tablebh(bsc,hsc)%bme(tablebhn(bsc,hsc)))
           tn=tablebhn(bsc,hsc)
           if (tn/=0) then
           do nsc=1,tn
           tablebh(bsc,hsc)%bme(nsc)%bn=0
           end do
           end if
           end if
        end do
        end do
        tablebhn=0
    end if
    
    if (store==2) then
        do  bsc=1,basd
        do  hsc=nja2,xja2,dja2
           tn=tablebhn(bsc,hsc)
           if (tn/=0) then
           do nsc=1,tn
           bn=tablebh(bsc,hsc)%bme(nsc)%bn
           if (bn/=0) then
           allocate(tablebh(bsc,hsc)%bme(nsc)%me(bn))
           allocate(tablebh(bsc,hsc)%bme(nsc)%b(bn))
           allocate(tablebh(bsc,hsc)%bme(nsc)%h(bn))
           tablebh(bsc,hsc)%bme(nsc)%bn=0
           end if
           end do
           end if
        end do
        end do
        tablebhn=0
    end if
    
    do bsc=1,basd
    ajl=ichar(bas(bsc)(1:1))
    bjl=ichar(bas(bsc)(2:2))
    ajr=min_cJ2a-del_cJ2a
    do hsc=nja2,xja2,dja2   
    ajr=ajr+del_cJ2a
    do ashl=min_ashells,max_ashells
    do ashr=min_ashells,max_ashells
    lmn=abs(sj2(ashl)-sj2(ashr))
    lmx=min(abs(sj2(ashl)+sj2(ashr)),max_lm)
    lmn=lmn+mod(lmn,del_lm)
    do lm=lmn,lmx,del_lm
        if (ajr>=abs(ajl-lm).and.ajr<=ajl+lm) then   
        ta(isl:isl)=char(ashl)
        ta(isr:isr)=char(ashr)
        ta(ijl:ijl)=char(ajl)
        ta(ijr:ijr)=char(ajr)
        ta(ilm:ilm)=char(lm)
        ta(6:8)=z1(6:8)
        atsc=tfind(trda,trdda,ta)
        if (atsc==0) go to 100
        tablebhn(bsc,hsc)=tablebhn(bsc,hsc)+1
        if (store==1) then
             tablebh(bsc,hsc)%a(tablebhn(bsc,hsc))=atsc
        end if
        if (store<1) go to 100
        bjrn=max(abs(bjl-lm),abs(ajr-J2))
        bjrx=min(abs(bjl+lm),abs(ajr+J2),max_cJ2b)
        bjrn=bjrn+mod(bjrn-min_cJ2b,del_cj2b)
        do bjr=bjrn,bjrx,del_cJ2b
        rn(1:2)=bas(bsc)(1:2)
        rn(3:3)=char(ajr)
        rn(4:4)=char(bjr)
        rn(5:5)=char(lm)
        rn(6:8)=z1(6:8)
        rsc=tfind(rac,racd,rn)
        if (rsc==0) go to 90
        bm(1:2)=rn(3:4)
        bm(3:8)=z1(3:8)
        nsc=tfind(bas,basd,bm)
        if (nsc==0) stop '  find error bh tables '
        do bshl=min_bshells,max_bshells
        do bshr=min_bshells,max_bshells
        lmn=abs(sj2(bshl)-sj2(bshr))
        lmx=min(abs(sj2(bshl)+sj2(bshr)),max_lm)
        if (lm<lmn.or.lm>lmx) go to 50
        if (bjr>=abs(bjl-lm).and.bjr<=bjl+lm) then   
            tb(isl:isl)=char(bshl)
            tb(isr:isr)=char(bshr)
            tb(ijl:ijl)=char(bjl)
            tb(ijr:ijr)=char(bjr)
            tb(ilm:ilm)=char(lm)
            tb(6:8)=z1(6:8)
            btsc=tfind(trdb,trddb,tb)
            if (btsc==0) go to 50
            op(1:2)=ta(isl:isr)
            op(3:4)=tb(isl:isr)
            op(5:5)=tb(ilm:ilm)
            op(6:8)=z1(6:8)
            psc=tfind(opn,opd,op)
            if (psc==0) go to 50
            me=opmeV(psc)
            if (me==0.0) go to 50
            tablebh(bsc,hsc)%bme(tablebhn(bsc,hsc))%bn=&
                         tablebh(bsc,hsc)%bme(tablebhn(bsc,hsc))%bn+1
            if (store==2) then
                 bn=tablebh(bsc,hsc)%bme(tablebhn(bsc,hsc))%bn
                 tablebh(bsc,hsc)%bme(tablebhn(bsc,hsc))%b(bn)=btsc
                 tablebh(bsc,hsc)%bme(tablebhn(bsc,hsc))%h(bn)=nsc
                 tablebh(bsc,hsc)%bme(tablebhn(bsc,hsc))%me(bn)=me*racahs(rsc)
            end if
        end if
50      end do
        end do
90      end do
        end if
100 end do
    end do
    end do
    end do
    end do
   
   end do
   
   deallocate(opn)
   deallocate(opmeV)
   deallocate(trda)
   deallocate(trdb)
   deallocate(rac)
   deallocate(racahs)
   
   end subroutine drdbh 
   
                            
end module Lnz

!  NuTRL.f90 
!
!  FUNCTIONS:
!  NuTRL      - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: NuTRL
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program NuTRL
    
    use Lnz
    use OperParameters
    use ProjParameters
    use InputNuShell
    use OutputNuShell
    use Mscheme
    use Files
    use Extra
    use Globals
    
    implicit none

    ! Variables
    
    real:: t1,t2,et1,et2,ms
    integer:: hr,mn,sc

    integer:: j,tmain,k,dt,k1,res,i
    
    character(len=6):: inExt
    character(len=10):: btime
    character(len=8):: bdate
    character(len=10)::chnostring,chstring='NOCSTRINGS'
    
    ! Body of NuTRL
    
    
    res=get_env_v(chstring,chnostring)
    if (res /= 0) then
        read(chnostring,*) nstring
    else
        nstring=max_part_tree
    end if
    
    
    call output_header(6)
    call input_nucleus(inExt)

    call output_welcome(6)    
    call input(inExt)

    call output_shells(6)    
    call check_shell_order()
    call setup_shells()
    if (use_isospin) then
        dt=2
    else
        dt=0
    end  if
    call setup_mscheme()
    call setup_file_ext()
    call input_intfile(opfile1,4)
    call input_intfile(opfile2,5)
    call input_intfile(intfile,6)
    
    call output_end_setup(6)        
    call cpu_time(t1)
    call date_and_time(bdate,btime)
    if (output_control > 6) open(unit=60,file=Nucleus//'trd.txt',action='write')
    if (output_control > 0) open(unit=10,file=Nucleus//'NuTRL'//'.txt',action='write')

    open(unit=801,file=opfile1//'.tri',form='unformatted',action='read')
    open(unit=806,file=opfile1//'.trd',form='unformatted',action='read')
    open(unit=802,file=opfile2//'.tri',form='unformatted',action='read')
    open(unit=803,file=opfile2//'.trd',form='unformatted',action='read')
    open(unit=805,file=intfile//'.oph',form='unformatted',action='read')
    open(unit=10,file=intfile//'.spe')
    restart=0
    read(10,*)
    read(10,*) renorm
    read(10,*) no_level
    read(10,*,end=10) restart
10  close(unit=10)
    call Lnz_TRD(800)
    close(unit=800)
    close(unit=802)
    close(unit=801)
    close(unit=803)
    read(btime,'(i2,i2,i2,f4.3)') hr,mn,sc,ms
    et1=hr*3600.+mn*60.+sc*1. +ms
    call date_and_time(bdate,btime)
    call cpu_time(t2)
    read(btime,'(i2,i2,i2,f4.3)') hr,mn,sc,ms
    et2=hr*3600.+mn*60.+sc*1. +ms
    t1=t2-t1
    et1=et2-et1
    if (et1<0.0) et1=et1+24.*3600.
 !   write(701) Nucleuss,J2,et1
    call output_completed(6,'NuVec      ')
    call output_time('NuTRL',t1,et1)
    print *
    
    end program NuTRL
    
