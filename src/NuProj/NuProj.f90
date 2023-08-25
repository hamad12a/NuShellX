! These are the modules specific tto the generic pProj for NuShell, NuShell@MSU and
! SunShell.

module ProjAccuracy

use ProjParameters

real(kind=rk):: ran_limit = .01           ! limit rep  "zero" in random vectors
real(kind=rk):: orthog_limit = .01        ! limit rep  "zero" in orthogonalisation
real(kind=rk):: stop_limit=.0001           ! norm limit for stopping Lanczos
real(kind=rk):: err_limit=.0000001         ! error limit on J^2 , T^2 for Lanczos
real(kind=rk):: ran_limrd=3.0              ! reduce ran_limit for every dim tries
real(kind=rks):: rnms_limit=.00001        ! limit for orthog process
integer:: nxtra=7                                   ! extra iterations of Lanczos if it does not terminate

contains

    subroutine read_parameters()
        implicit none
        logical::prm_exist
        
        prm_exist=.false.
        inquire(file='NuProj.prm',exist=prm_exist)
        
        if (prm_exist) then
        open(unit=13,file='NuProj.prm')
        read(13,*,err=10) ran_limit            ! increase to be more selective but not too far .01 max
        read(13,*,err=10) orthog_limit         ! increase to improve orthogonality 0.01 max
        read(13,*,err=10) stop_limit           ! limit on beta to stop Lanczos iterations if used
        read(13,*,err=10) err_limit            ! in Lanczos method J^2 and T^2 limits
        read(13,*,err=10) ran_limrd            ! every nJT tries for a vector reduces ran_limit by this
        read(13,*,err=10) rnms_limit           ! must be greater than zero. Stops sqrt of neg nummbers and Norm overflows
        read(13,*,err=10) nxtra                ! extra iterations for Lanczos if it does not terminate.
                                        ! odd numbers seem to work best and increasing can help
                                        ! difficult cases, but it is  already quite high
        close(unit=13)
        else
        open(unit=13,file='NuProj.prm',action='write')
        write(13,*,err=20) ran_limit
        write(13,*,err=20) orthog_limit
        write(13,*,err=20) stop_limit
        write(13,*,err=20) err_limit
        write(13,*,err=20) ran_limrd
        write(13,*,err=20) rnms_limit
        write(13,*,err=20) nxtra
        close(unit=13)
        end if
        return
10      print *, ' Error reading NuProj.prm'
        stop
20      print *, ' Error writing NuProj.prm'
        stop
                  
    end subroutine read_parameters
    
end module ProjAccuracy

module IntegerProject

use Shells
use Mscheme
use Mvector
use ShellsX
use ProjParameters

implicit none

real(kind=rk),dimension(0:max_sj2):: factorial,sqrt_fact,inv_fact,inv_sqrt_fact

contains

    subroutine setup_factorial()
        implicit none
        integer:: i
        
        factorial(0)=real(1,rk)
        factorial(1)=real(1,rk)
        sqrt_fact(0)=real(1,rk)
        sqrt_fact(1)=real(1,rk)
        inv_fact(0)=real(1,rk)
        inv_fact(1)=real(1,rk)
        inv_sqrt_fact(0)=real(1,rk)
        inv_sqrt_fact(1)=real(1,rk)
        
        do i=2,max_sj2
            factorial(i)=real(i,rk)*factorial(i-1)
            sqrt_fact(i)=sqrt(factorial(i))
            inv_fact(i)=real(1,rk)/factorial(i)
            inv_sqrt_fact(i)=real(1,rk)/sqrt_fact(i)
        end do
        
    end subroutine setup_factorial
        
    function norm(mvec,part)
        implicit none
        real(kind=rk):: norm,low,high
        type(spartition),intent(in)::part
        type(m_vector),intent(in):: mvec
        type(m_vector):: mvec_copy
        integer:: iw,nb,np
        integer(kind=ik):: iword
        logical:: odd
        logical,dimension(max_no_shells)::shell_logic,shell_logico,shell_logice
        
        shell_logic(1:no_shells)=.false.
        shell_logico(1:no_shells)=.false.
        shell_logice(1:no_shells)=.false.
        
        norm=real(1,rk)
        mvec_copy=mvec
        if (.not.use_isospin) then
        do shn=1,no_shells
            if ( part%shell(shn) > (sj2(shn)+1)/2) then
                mvec_copy%fn(shell_word(shn))=ieor(mvec_copy%fn(shell_word(shn))&
                                                   ,shell_mask(shn)) 
                shell_logic(shn)=.true.
            end if
        end do
        else
        do shn=1,no_shells
            iword=iand(mvec_copy%fn(shell_word(shn)),shell_masko(shn))
            np=countb(iword)
            if ( np > (sj2(shn)+1)/2) then
                mvec_copy%fn(shell_word(shn))=ieor(mvec_copy%fn(shell_word(shn))&
                                                 ,shell_masko(shn)) 
                shell_logico(shn)=.true.
            end if
            iword=iand(mvec_copy%fn(shell_word(shn)),shell_maske(shn))
            np=countb(iword)
            if ( np > (sj2(shn)+1)/2) then
                mvec_copy%fn(shell_word(shn))=ieor(mvec_copy%fn(shell_word(shn))&
                                                ,shell_maske(shn)) 
                shell_logice(shn)=.true.
            end if
        end do
        end if
        do iw=1,gwords()
            do
                if (mvec_copy%fn(iw) == 0) exit
                nb=next_bit(mvec_copy%fn(iw))
                mvec_copy%fn(iw)=ibclr(mvec_copy%fn(iw),nb-1)
                odd=.false.
                if (btest(nb,0)) odd=.true.
                shn=shn_lu(nb,iw)
                if (shell_logic(shn) .or. (shell_logico(shn) .and. odd) .or. &
                                          (shell_logice(shn) .and. .not.odd)) then
                norm=norm*inv_sqrt_fact((sj2_lu(nb,iw)+sm2_lu(nb,iw))/2)   &
                          *inv_sqrt_fact((sj2_lu(nb,iw)-sm2_lu(nb,iw))/2)
                else
                norm=norm*sqrt_fact((sj2_lu(nb,iw)+sm2_lu(nb,iw))/2)   &
                          *sqrt_fact((sj2_lu(nb,iw)-sm2_lu(nb,iw))/2)
                end if          
            end do
        end do
        
    end function norm

end module IntegerProject

module Project

! Integer Projection version
use Files
use Finds
use ProjParameters
use ShellsX
use Mscheme
use Mvector
use InputNuShell
use TriDiag
use ProjAccuracy
use DotProduct
use IntegerProject

integer:: dt,err
logical::disk
integer(kind=2)::res2
character(len=6):: inExt
type(m_vector),dimension(:),allocatable:: mvecr,mvecrJp,mvecrTp
integer,dimension(:),allocatable:: used
type(m_vector),dimension(:),allocatable:: mvec_pred
integer,dimension(max_nJT):: igood
real(kind=rk),dimension(:),allocatable,target:: coefr1,coefr2
real(kind=rk),dimension(:,:),allocatable,target:: coef
real(kind=rk),dimension(:,:),pointer:: coef0
real(kind=rk),dimension(:),allocatable:: rnorm,rnorm2
real(kind=rk),dimension(:,:),allocatable:: JJcoef,Jpcoef
real(kind=rk),dimension(:,:),allocatable:: coef_st
integer(kind=1),dimension(:,:),allocatable:: matxJp 
integer(kind=1),dimension(:,:),allocatable:: matxJm 
integer(kind=1),dimension(:,:),allocatable:: matxJpr 
integer(kind=1),dimension(:,:),allocatable:: matxJmr 
integer(kind=1),dimension(:,:),allocatable:: matxJpt 
integer(kind=1),dimension(:,:),allocatable:: matxJmt 
integer(kind=1),dimension(:,:),allocatable:: matxJptr 
integer(kind=1),dimension(:,:),allocatable:: matxJmtr 
integer:: mxopJp,mxopJm,mxopTp,mxopTm
integer,dimension(:,:),allocatable:: inxJp
integer,dimension(:,:),allocatable:: inxTp
integer,dimension(:,:),allocatable:: inxJm
integer,dimension(:,:),allocatable:: inxTm
integer,dimension(:,:),allocatable:: inxJpt
integer,dimension(:,:),allocatable:: inxTpt
integer,dimension(:,:),allocatable:: inxJmt
integer,dimension(:,:),allocatable:: inxTmt
integer,dimension(:),allocatable:: np
real(kind=rk),dimension(:,:,:),allocatable,target:: vector
integer:: Jdim,Tdim,dim2,ndisc,tz2,jz2
logical:: predict,use_real
real(kind=rk):: one_over_sqrt,orthlimit
real(kind=rk):: sqrtf(-127:127)
integer::nthr
integer,parameter:: mask=Z'00FFFFFF'
integer,parameter:: mmask=Z'000000FF'

contains

    subroutine do_projection(part,Jz,Tz,Jmax,Tmax,dim,nJT,dimJp,dimTp)
    implicit none 
    type(spartition),intent(in)::part
    integer,intent(in):: dim,Jz,Tz,Jmax,Tmax,nJT,dimJp
    integer::ii, ir,Jmin,Tmin,ngood,Tn,nstp,ngoodl,dimTp,nst,maxtNJ,maxtNT
    real(kind=rk):: rnm,randlim,high,low
    real(kind=rk),dimension(4):: errJ,errT
    real(kind=rks):: rnms
    real rndm
    logical:: lnorm,ljt2,lnmrs,logic,lrand
    integer,dimension(4):: i
    jz2=Jz
    tz2=Tz
    high=huge(rnm)/10.
    low=10.*epsilon(rnm)
    if (dim==1) then
        coef_st(1,1)=real(1,rk)
        if (disk) write(66,rec=1) coef_st
        ngood=1
        igood(1)=1
        return
    end if
    randlim=ran_limit
    orthlimit=orthog_limit
    tN=no_cI+no_cN
    if (.not.use_isospin) dimTp=0
    Jmin=Jz+2
    Tmin=Tz+2
    ngood=0
    ngoodl=0
    logic=.true.
! allocate memory for matrices
    allocate(np(max(dim,dimJp,dimTp)))
    allocate(rnorm(dim))
    allocate(rnorm2(dim))
    use_real=.false.
    do ii=1,dim
        rnorm(ii)=norm(mvecr(ii),part)
        rnorm2(ii)=rnorm(ii)*rnorm(ii)
        if (rnorm2(ii)>high .or. rnorm2(ii)<low) then
            use_real=.true.
            exit
        end if
    end do
    if (use_real) then
        rnorm(1:dim)=real(1,rk)
        rnorm2(1:dim)=real(1,rk)
    end if
    call alloc(tN,tN,dim,dimJp,dimTp,use_real)
    used(1:dim)=1
    Jdim=(Jmax-Jmin)/2+3+nxtra
    Tdim=(Tmax-Tmin)/2+3+nxtra
    if (Jdim+1 > max_lanc  .or. Tdim+1  > max_lanc) then
        print *, ' Increase  max_lanc in ProjParameters '
        stop
    end if
    dim2=max(Jdim,Tdim)
    allocate(vector(dim,4,dim2))
!setup JmJp matrix and TmTp matrix if needed
    call setup_JmJp(tN,dim,dimJp)
    if (use_isospin)     call setup_TmTp(tN,dim,dimTp)
    call setup_random(Jmax,Jmin,Tmax,Tmin,Jz,Tz,dimJp,dimTp,errJ,errT,dim,tN)
! start loop over m vectors    
        ndisc=0
    ir=0
1   do 
      ir=ir+1
      lrand=.false.
      i=0
      
      do ii=1,min(nJT-ngood,4)
2     if (maxval(used(1:dim))==0) go to 3
      call random_number(rndm)
      i(ii)=nint(rndm*dim)
      if (i(ii) < 1 .or. i(ii) > dim ) go to 2   
      if (used(i(ii))==0) go to 2
      used(i(ii))=0
! check if candidate for projection -  Random Vector Prediction
      if (abs(coefr1(i(ii))) < randlim .and. abs(coefr2(i(ii))) < randlim .and. dim/=1) go to 2
      end do
3     if (i(1)==0) go to 10
      lrand=.true.
      if (output_control > 13) write(10,*)'i,r1,r2,rlim',i, coefr1(i(1:4)),&
                                    coefr2(i(1:4)),randlim    
        lnorm=.true.
        lnmrs=.true.
        lrand=.true.
        ljt2=.true.
!project
        call Lanczos(Jmax,Tmax,dimJp,dimTp,errJ,errT,dim,i,.false.)
!check norm and J^2 error and T62 error
!if this is first vector
        if (ngood ==  0) then
                    if (output_control > 6) write(10,*) ' N, ngood, nJT, dim,rnm,errj,errt ', &
                                  i,ngood,nJT,dim,rnm,errj,errt,lnorm,ljt2
                    do ii=1,min(nJT-ngood,4)
                    if (i(ii)==0) go to 4              
                    if (errJ(ii) > err_limit .or. errT(ii) > err_limit) go to 4
                    if (abs(coef(i(ii),ii)) < orthlimit) go to 4
                    if (ngood==0) then
                        logic=.false.
                    else
                        logic=.true.
                    end if
                    igood(ngood+1)=i(ii)
!                    coef_st(1:dim,ngood+1)=coef(1:dim,ii)
                    if (ngood>0) then
                        call apply_inverse(ngood+1,dim,logic,ii)
                    else
                        if (disk) then
                            write(66,rec=1) coef(1:dim,ii)
                        else
                            coef_st(1:dim,ngood+1)=coef(1:dim,ii)
                        end if    
                    end if
                    if (logic) go to 4
                    ngood=ngood+1
                    if (ngood >= nJT) exit
                    if (ngood > ngoodl) then
                        call subtract_random(dim,ngood,ii)
                        ngoodl=ngood
                    end if
4                   end do
                    if (ngood >= nJT) exit
        else
                    if (output_control > 6) write(10,*) ' N, ngood, nJT, dim,rnm,errj,errt ', &
                                  i,ngood,nJT,dim,rnm,errj,errt,lnorm,ljt2
!accept                
                    do ii=1,min(nJT-ngood,4)
                    if (i(ii)==0) go to 5              
                    if (errJ(ii) > err_limit .or. errT(ii) > err_limit) go to 5
                    if (abs(coef(i(ii),ii)) < orthlimit) go to 5
                        logic=.true.
                        igood(ngood+1)=i(ii)
!                        coef_st(1:dim,ngood+1)= coef(1:dim,ii)
                        call apply_inverse(ngood+1,dim,logic,ii)
                        if (logic) go to 5
                        ngood=ngood+1
                        if (ngood >= nJT) exit
                        if (ngood > ngoodl) then
                            call subtract_random(dim,ngood,ii)
                            ngoodl=ngood
                        end if
5                   end do
                    if (ngood >= nJT) exit
        end if
! if  new vector has been found subtract it from random vectors
        if (output_control > 3) print '(a21,7i8)',' N(1,4),ngood,nJT,dim', i,ngood,nJT,dim
10    if (dot_product(used,used)==0) exit
      end do
    if ( ngood /= nJT) then
        randlim=randlim/ran_limrd
        orthlimit=min(orthlimit,real(2,rk)*orthlimit/ran_limrd)
        if (randlim < 0.00005) then
            print *, ' Error in Project , ngood, nJT, J*2,T*2',ngood,nJT,Jz,Tz
            stop 1111
        end if
        if (output_control > 1) print *, ' Reducing limits ',randlim
        used(1:dim)=1
        do ii=1,dim
            if (abs(coefr1(ii))<randlim .and. abs(coefr2(ii))<randlim) used(ii)=0
        end do
        do ii=1,ngood
            used(igood(ii))=0
        end do
        if (dot_product(used,used)==0) then
            call setup_random(Jmax,Jmin,Tmax,Tmin,Jz,Tz,dimJp,dimTp,errJ,errT,dim,tN)
            used(1:dim)=1
            do ii=1,ngood
            used(igood(ii))=0
            end do
        end if
        go to 1
    end if
    call dealloc()
    end subroutine do_projection
    
    subroutine apply_inverse(ngood,dim,logic,ii)
        implicit none
        integer,intent(in):: ngood,dim,ii
        logical,intent(inout):: logic
        integer:: ng1,ngu
        real(kind=rks):: rnms
        if (abs(coef(igood(ngood),ii))<orthlimit) return
        do ng1=1,ngood-1
            ngu=ng1
            if (disk) then
                read(66,rec=ng1) coef_st(1:dim,1)
                ngu=1
            end if
            rnms=coef(igood(ng1),ii)/coef_st(igood(ng1),ngu)
            if (real(1,rks)-rnms*rnms<=rnms_limit) return
            one_over_sqrt=real(1,rks)/sqrt((real(1,rks)-rnms*rnms))
            coef(1:dim,ii)=(coef(1:dim,ii)-rnms*coef_st(1:dim,ngu))*one_over_sqrt
            if (abs(coef(igood(ngood),ii))<orthlimit) return        
        end do
        if (coef(igood(ngood),ii)<0.0) coef(:,ii)=-coef(:,ii)
        if (disk) then
            write(66,rec=ngood) coef(1:dim,ii)
        else
            coef_st(1:dim,ngood)=coef(1:dim,ii)
        end if
        logic=.false.
    end subroutine apply_inverse
    
    subroutine setup_random(Jmax,Jmin,Tmax,Tmin,Jz,Tz,dimJp,dimTp,errJ,errT,dim,tN)
    implicit none
    integer,intent(in)::Jmax,Jmin,Tmax,Tmin,Jz,Tz,dimJp,dimTp,dim,tN
    real(kind=rk),dimension(4),intent(out):: errJ,errT
    integer::i,m
    integer,dimension(4)::ii
    integer,dimension(:),allocatable:: new
    real(kind=rk)::rnm
    real::rndm
! use a fixed seed
    call random_seed(SIZE=m)
    allocate(new(m))
    new=5
    call random_seed(PUT=new(1:m))
    
    ii=0
!construct  random vectors
        do i=1,dim
            call random_number(rndm)
            coef(i,1)=real(rndm - 0.5,rk)
        end do
        do i=1,dim
            call random_number(rndm)
            coef(i,2)=real(rndm - 0.5,rk)
        end do
        rnm=dot(coef(1:dim,1),coef(1:dim,1),dim)
        rnm=real(1,rk)/sqrt(rnm)
        coef(1:dim,1)=coef(1:dim,1)*rnm/rnorm(1:dim)
        rnm=dot(coef(1:dim,2),coef(1:dim,2),dim)
        rnm=real(1,rk)/sqrt(rnm)
        coef(1:dim,2)=coef(1:dim,2)*rnm/rnorm(1:dim)
!project
        ii(1)=-1
        ii(2)=-1
        call Lanczos(Jmax,Tmax,dimJp,dimTp,errJ,errT,dim,ii,.false.)
        coefr1(1:dim)=coef(1:dim,1)
        coefr2(1:dim)=coef(1:dim,2)
        rnm=dot(coefr1,coefr2,dim)
        coefr2(1:dim)=coefr2(1:dim)-rnm*coefr1(1:dim)
        rnm=sqrt(dot(coefr2,coefr2,dim))
        if (rnm/=real(0,rk)) rnm=real(1,rk)/rnm
        coefr2(1:dim)=coefr2(1:dim)*rnm
    
    end subroutine setup_random
    
    subroutine subtract_random(dim,ngood,ii)
        implicit none
        integer,intent(in):: dim,ngood,ii
        real(kind=rk):: rnm
        
            rnm=coefr1(igood(ngood))/coef(igood(ngood),ii)
            if (rnm>=real(1,rk)) then
            one_over_sqrt=real(0,rk)
            else
            one_over_sqrt=real(1,rk)/sqrt(real(1,rk)-rnm*rnm)
            end if
            coefr1(1:dim)=(coefr1(1:dim)-coef(1:dim,ii)*rnm)*one_over_sqrt
            rnm=coefr2(igood(ngood))/coef(igood(ngood),ii)
            if (rnm>=real(1,rk)) then
            one_over_sqrt=real(0,rk)
            else
            one_over_sqrt=real(1,rk)/sqrt(real(1,rk)-rnm*rnm)
            end if
            coefr2(1:dim)=(coefr2(1:dim)-coef(1:dim,ii)*rnm)*one_over_sqrt
            
    end subroutine subtract_random
        
    subroutine alloc(nopJ,nopT,dim,dim1,dim2,ureal)
        implicit none
        integer,intent(in):: nopJ,nopT,dim,dim1,dim2
        logical,intent(in):: ureal
        allocate(used(dim))
        allocate(coef(dim,4))
        allocate(coefr1(dim))
        allocate(coefr2(dim))
        allocate(JJcoef(dim,4))
        allocate(Jpcoef(max(dim1,dim2),4))
        allocate(matxJpt(nopJ,max(dim,dim1)))
        if (ureal) allocate(matxJptr(nopJ,max(dim,dim1)))
        allocate(inxJpt(nopJ,max(dim,dim1)))
        if (use_isospin) allocate(inxTpt(nopT,dim))
        allocate(matxJmt(nopJ,max(dim,dim1)))
        if (ureal) allocate(matxJmtr(nopJ,max(dim,dim1)))
        allocate(inxJmt(nopJ,max(dim,dim1)))
        if (use_isospin) allocate(inxTmt(nopT,dim))
    end subroutine alloc
    
    subroutine dealloc()
        implicit none
        deallocate(vector)
        deallocate(used)
        deallocate(coef)
        deallocate(coefr1)
        deallocate(coefr2)
        deallocate(JJcoef)
        deallocate(Jpcoef)
        deallocate(rnorm)
        deallocate(rnorm2)
        if (allocated(matxJp)) deallocate(matxJp)
        if (allocated(matxJpr)) deallocate(matxJpr)
        deallocate(inxJp)
        if (use_isospin) deallocate(inxTp)
        if (allocated(matxJm)) deallocate(matxJm)
        if (allocated(matxJmr)) deallocate(matxJmr)
        deallocate(inxJm)
        if (use_isospin) deallocate(inxTm)
        deallocate(np)
    end subroutine dealloc

    subroutine setup_JmJp(nop,dim,dimJp)
        implicit none
        integer,intent(in):: nop,dim,dimJp
        type(m_vector):: mvec_copy1,mvecJp,mvec1
        integer:: dn,nb1,iw1,np1,np2,findm,sh,mshft
        integer(kind=1):: fact1,fact2,fact3,fact4,phase1
        integer,dimension(2):: opa
        real,dimension(nwords*nbits):: rnp2
        integer,dimension(nwords*nbits):: inx,inxre
        integer(kind=1),dimension(nwords*nbits):: matre,matrer
        integer(kind=ik):: iand,ibclr,not,ieor,ishft
              
        mxopJp=0
        mxopJm=0
        
        if (use_real) then
            do dn=-127,127
                sqrtf(dn)=sqrt(real(abs(dn),rk))
            end do
        end if
        
        inxJpt=0
        inxJmt=0
        !outer loop over m_vectors, inner loops over nucleons
        np=0  
        do dn=1,dim
            np2=0      
            mvec_copy1=mvecr(dn)
            do sh=1,no_shells
                iw1=shell_word(sh)
                if (use_isospin) then
                    if (iand(mvec_copy1%fn(iw1),shell_masko(sh))==shell_masko(sh)) &
                             mvec_copy1%fn(iw1)=iand(mvec_copy1%fn(iw1),not(shell_masko(sh)))
                    if (iand(mvec_copy1%fn(iw1),shell_maske(sh))==shell_maske(sh)) &
                             mvec_copy1%fn(iw1)=iand(mvec_copy1%fn(iw1),not(shell_maske(sh)))
                else
                    if (iand(mvec_copy1%fn(iw1),shell_mask(sh))==shell_mask(sh)) &
                             mvec_copy1%fn(iw1)=iand(mvec_copy1%fn(iw1),not(shell_mask(sh)))
                end if
            end do
            mvec1=mvecr(dn)
            do iw1=1,gwords()
                opa(1)=iw1
                do 
                    if (mvec_copy1%fn(iw1) == 0) exit
                    nb1=next_bit(mvec_copy1%fn(iw1))
                    mvec_copy1%fn(iw1)=ibclr(mvec_copy1%fn(iw1),nb1-1)
                    if (sm2_lu(nb1,iw1) == sj2_lu(nb1,iw1)) go to 10
                    if (mtest(mvec1%fn(iw1),nb1+1+st2_lu(nb1,iw1))) go to 10 ! btest numbers bits from 0
                    opa(2)=nb1
                    mvecJp=opa-mvec1
                    opa(2)=opa(2)+st2_lu(nb1,iw1)+1
                    mvecJp=opa+mvecJp
                    phase1=1
                    fact1=(sj2_lu(nb1,iw1)-sm2_lu(nb1,iw1))/2
                    fact2=(sj2_lu(nb1,iw1)+sm2_lu(nb1,iw1)+2)/2
                    if (st2_lu(nb1,iw1) == 1) then 
                        if (mtest(mvec1%fn(iw1),nb1+1)) phase1=-phase1
                    end if 
                    findm=bfind(mvecrJp,dimJp,mvecJp)
                    if (findm==0) then
                        call swrite(6,mvecr(dn))
                        call swrite(6,mvecJp)
                        print *,dimJp
                        stop
                    end if
                    np(findm)=np(findm)+1
                    if (np(findm)>mxopJp) mxopJp=np(findm)
                    np2=np2+1
                    if (np2>mxopJm) mxopJm=np2
                    inxJpt(np(findm),findm)=dn
                    inxJmt(np2,dn)=findm
                    matxJpt(np(findm),findm)=fact1*phase1
                    if (use_real) matxJptr(np(findm),findm)=fact2
                    matxJmt(np2,dn)=fact2*phase1
                    if (use_real) matxJmtr(np2,dn)=fact1
10              end do
            end do                                        
        end do
        allocate(inxJp(mxopJp,dimJp))
        if (use_real) allocate(matxJp(mxopJp,dimJp))
        if (use_real) allocate(matxJpr(mxopJp,dimJp))
        allocate(inxJm(mxopJm,dim))
        if (use_real) allocate(matxJm(mxopJm,dim))
        if (use_real) allocate(matxJmr(mxopJm,dim))
        do dn=1,dimJp
            np2=0
            do sh=1,mxopJp
                np2=np2+1
                rnp2(np2)=inxJpt(np2,dn)
                if (rnp2(np2)==0.0) then
                    np2=np2-1
                    exit
                end if
            end do
            if (np2==0) go to 20
            call indexx(np2,rnp2,inx)
            do sh=1,np2
                matre(sh)=matxJpt(inx(sh),dn)
                if (use_real) matrer(sh)=matxJptr(inx(sh),dn)
                if (.not.use_real .and. inxJpt(inx(sh),dn)/=0 ) then
                    mshft=matre(sh)
                    mshft=iand(mshft,mmask)
                    inxre(sh)=ieor(inxJpt(inx(sh),dn),ishft(mshft,24))
                else
                    inxre(sh)=inxJpt(inx(sh),dn)
                end if
            end do
            if (use_real) matxJp(1:np2,dn)=matre(1:np2)
            if (use_real) matxJpr(1:np2,dn)=matrer(1:np2)
            inxJp(1:np2,dn)=inxre(1:np2)
20          if (use_real) matxJp(np2+1:mxopJp,dn)=0
            if (use_real) matxJpr(np2+1:mxopJp,dn)=0
            inxJp(np2+1:mxopJp,dn)=0
            
        end do
        do dn=1,dim
            np2=0
            do sh=1,mxopJm
                np2=np2+1
                rnp2(np2)=inxJmt(np2,dn)
                if (rnp2(np2)==0.0) then
                    np2=np2-1
                    exit
                end if
            end do
            if (np2==0) go to 30
            call indexx(np2,rnp2,inx)
            do sh=1,np2
                matre(sh)=matxJmt(inx(sh),dn)
                if (use_real) matrer(sh)=matxJmtr(inx(sh),dn)
                inxre(sh)=inxJmt(inx(sh),dn)
                if (.not.use_real .and. inxJmt(inx(sh),dn)/=0 ) then
                    mshft=matre(sh)
                    mshft=iand(mshft,mmask)
                    inxre(sh)=ieor(inxJmt(inx(sh),dn),ishft(mshft,24))
                else
                    inxre(sh)=inxJmt(inx(sh),dn)
                end if
            end do
            if (use_real) matxJm(1:np2,dn)=matre(1:np2)
            if (use_real) matxJmr(1:np2,dn)=matrer(1:np2)
            inxJm(1:np2,dn)=inxre(1:np2)
30          if (use_real) matxJm(np2+1:mxopJm,dn)=0
            if (use_real) matxJmr(np2+1:mxopJm,dn)=0
            inxJm(np2+1:mxopJm,dn)=0
        end do
        deallocate(inxJpt)
        deallocate(matxJpt)
        if (use_real) deallocate(matxJptr)
        deallocate(inxJmt)
        deallocate(matxJmt)
        if (use_real) deallocate(matxJmtr)
        return
        
    end subroutine setup_JmJp
    
    subroutine setup_TmTp(nop,dim,dimTp)
        implicit none
        integer,intent(in):: nop,dim,dimTp
        type(m_vector):: mvec_copy1,mvecTp,mvec1
        integer:: dn,nb1,iw1,np1,np2,findm,sh
        integer,dimension(2):: opa
        real,dimension(nwords*nbits):: rnp2
        integer,dimension(nwords*nbits):: inx,inxre
        integer(kind=ik):: ibclr
        
        mxopTp=0
        mxopTm=0
        
        inxTpt=0
        inxTmt=0
        !outer loop over m_vectors, inner loops over nucleons
        np=0  
        do dn=1,dim
            np2=0     
            mvec_copy1=mvecr(dn)
            mvec1=mvecr(dn)
            do iw1=1,gwords()
            opa(1)=iw1
                do
                    if (mvec_copy1%fn(iw1) == 0) exit
                    nb1=next_bit(mvec_copy1%fn(iw1))
                    mvec_copy1%fn(iw1)=ibclr(mvec_copy1%fn(iw1),nb1-1)
                    if (st2_lu(nb1,iw1)== 0) exit
                    if (stz2_lu(nb1,iw1) == st2_lu(nb1,iw1)) go to 10
                    if (mtest(mvec1%fn(iw1),nb1+1)) go to 10 ! btest numbers bits from 0
                    opa(2)=nb1
                    mvecTp=opa-mvec1
                    opa(2)=opa(2)+1
                    mvecTp=opa+mvecTp
                    findm=bfind(mvecrTp,dimTp,mvecTp)
                    if (findm==0) then
                        call swrite(6,mvecr(dn))
                        print *,dimTp
                        stop
                    end if
                    np(findm)=np(findm)+1
                    if (np(findm)>mxopTp) mxopTp=np(findm)
                    np2=np2+1
                    if (np2>mxopTm) mxopTm=np2
                    inxTpt(np(findm),findm)=dn
                    inxTmt(np2,dn)=findm
10              end do
            end do                                        
        end do
        allocate(inxTp(mxopTp,dimTp))
        allocate(inxTm(mxopTm,dim))
        do dn=1,dimTp
            np2=0
            do sh=1,mxopTp
                np2=np2+1
                rnp2(np2)=inxTpt(np2,dn)
                if (rnp2(np2)==0.0) then
                    np2=np2-1
                    exit
                end if
            end do
            if (np2==0) go to 20
            call indexx(np2,rnp2,inx)
            do sh=1,np2
                inxre(sh)=inxTpt(inx(sh),dn)
            end do
            inxTp(1:np2,dn)=inxre(1:np2)
20          inxTp(np2+1:mxopTp,dn)=0
        end do
        do dn=1,dim
            np2=0
            do sh=1,mxopTm
                np2=np2+1
                rnp2(np2)=inxTmt(np2,dn)
                if (rnp2(np2)==0.0) then
                    np2=np2-1
                    exit
                end if
            end do
            if (np2==0) go to 30
            call indexx(np2,rnp2,inx)
            do sh=1,np2
                inxre(sh)=inxTmt(inx(sh),dn)
            end do
            inxTm(1:np2,dn)=inxre(1:np2)
30          inxTm(np2+1:mxopTm,dn)=0
        end do
        deallocate(inxTpt)
        deallocate(inxTmt)
        
        return
        
    end subroutine setup_TmTp

    subroutine Lanczos(Jmax,Tmax,dimJp,dimTp,errJ,errT,dim,state,reorth)
        implicit none
        integer,intent(in)::Jmax,Tmax,dimJp,dimTp,dim
        integer,dimension(4),intent(in):: state
        logical,intent(in):: reorth
        real(kind=rk),dimension(4),intent(out):: errJ,errT
        integer,dimension(4):: endi
        real(kind=rk),dimension(max_lanc,4):: alpha,beta
        integer,dimension(max_lanc,4):: indx
        real(kind=rk),dimension(max_lanc,max_lanc,4):: eigenv
        integer:: i,ix,tqli_iter,ii,mxii
        real(kind=rk):: fact,fact1
        real(kind=rk),dimension(4):: ovl
        
        endi=0
        mxii=4
        do ii=1,mxii
        if (state(ii)<0) then
            vector(1:dim,ii,1)=coef(1:dim,ii)
        else
            if (state(ii)>0) then
            vector(1:dim,ii,1)=real(0,rk)
            vector(state(ii),ii,1)=real(1,rk)/rnorm(state(ii))
            else if (state(ii)==0) then
                mxii=ii-1
                exit
            end if
        end if
        end do
        eigenv(1:dim2,1:dim2,1:4)=real(0,rk)
        beta(1,1:mxii)=real(0,rk)
! Lanczos iteration       
        do i=1,Jdim-1
            eigenv(i,i,1:mxii)=real(1,rk)
            do ii=1,mxii
            coef0=>vector(:,:,i)
            end do
               if (use_real) then
                call PJmJpr(dim,dimJp,mxii)
               else
                call PJmJp(dim,dimJp,mxii)
               end if
            fact=real(Jmax*Jmax+2*Jmax+8,rk)/real(4,rk)
            if (jz2<0) then
                fact1=real(-jz2,rk)
            else
                fact1=real(0,rk)
            end if    
!subtract to make eienvalue desired largest negative  
!$OMP PARALLEL 
!$OMP DO PRIVATE (ii)       
            do ii=1,mxii       
            alpha(i,ii)=dot(JJcoef(1:dim,ii),coef0(1:dim,ii),rnorm2,dim)
            if (i /= 1) then
                JJcoef(1:dim,ii)=JJcoef(1:dim,ii)-alpha(i,ii)*coef0(1:dim,ii)&
                            -beta(i,ii)*vector(1:dim,ii,i-1)
            else
                JJcoef(1:dim,ii)=JJcoef(1:dim,ii)-alpha(i,ii)*coef0(1:dim,ii)
            end if
            alpha(i,ii)=alpha(i,ii)-fact
!reothogonalisation step appears not to be necessary with integer projection        
            if (reorth) call orthogonalise(i,vector,JJcoef(1:dim,ii),rnorm2,dim,ii)
            ovl(ii)=dot(JJcoef(1:dim,ii),JJcoef(1:dim,ii),rnorm2,dim)
!check for termination
            beta(i+1,ii)=sqrt(ovl(ii))
            if (abs(ovl(ii)) <= stop_limit .and. endi(ii)==0 ) then
                endi(ii)=i
            end if                
            if (ovl(ii)/=real(0,rk))ovl(ii)=real(1,rk)/sqrt(ovl(ii))
            vector(1:dim,ii,i+1)=JJcoef(1:dim,ii)*ovl(ii)
            end do
!$OMP END PARALLEL        
            if (minval(endi(1:mxii))/=0) exit
        end do
!$OMP PARALLEL 
!$OMP DO PRIVATE (ii,ix)       
        do ii=1,mxii
        if (endi(ii)==0) endi(ii)=Jdim
!tridiagonalisation        
        tqli_iter=max_tqli        
        call tqli(alpha(:,ii),beta(:,ii),endi(ii),max_lanc,eigenv(:,:,ii),tqli_iter)
        alpha(1:endi(ii),ii)=abs(alpha(1:endi(ii),ii)+fact-fact1)
!        call indexx(i,alpha(1:endi(ii),ii),indx(1:endi(ii),ii))
        indx(1,ii)=maxval(minloc(alpha(1:endi(ii),ii)))
        if (output_control > 10) write(10,*) alpha(1:endi(ii),ii),'f',fact
        
        coef(1:dim,ii)=real(0,rk)
!calculate lowest eigenvector
        do ix = 1,endi(ii)
            coef(1:dim,ii)=coef(1:dim,ii)+eigenv(ix,indx(1,ii),ii)*vector(1:dim,ii,ix)
        end do
        coef(1:dim,ii)=coef(1:dim,ii)*rnorm(1:dim)
        errJ(ii)=alpha(indx(1,ii),ii)
        errT(ii)=real(0,rk)
        end do
!$OMP END PARALLEL        
        if (.not.use_isospin) go to 100
        endi=0
        vector(1:dim,1:mxii,1)=coef(1:dim,1:mxii)
        eigenv(1:dim2,1:dim2,1:mxii)=real(0,rk)
        beta(1,1:4)=real(0,rk)
! Lanczos iteration       
        do i=1,Tdim-1
            eigenv(i,i,1:mxii)=real(1,rk)
            do ii=1,mxii
            coef0=>vector(:,:,i)
            end do
            call PTmTp(dim,dimTp,mxii)
            fact=real(Tmax*Tmax+2*Tmax+8,rk)/real(4,rk)
!subtract to make eienvalue desired largest negative  
!$OMP PARALLEL 
!$OMP DO PRIVATE (ii)       
            do ii=1,mxii          
            alpha(i,ii)=dot(JJcoef(1:dim,ii),coef0(1:dim,ii),dim)
            if (i /= 1) then
                JJcoef(1:dim,ii)=JJcoef(1:dim,ii)-alpha(i,ii)*coef0(1:dim,ii)&
                            -beta(i,ii)*vector(1:dim,ii,i-1)
            else
                JJcoef(1:dim,ii)=JJcoef(1:dim,ii)-alpha(i,ii)*coef0(1:dim,ii)
            end if
            alpha(i,ii)=alpha(i,ii)-fact
!reothogonalisation step appears not to be necessary with integer projection        
            if (reorth) call orthogonalise(i,vector,JJcoef(1:dim,ii),dim,ii)
            ovl(ii)=dot(JJcoef(1:dim,ii),JJcoef(1:dim,ii),dim)
!check for termination
            beta(i+1,ii)=sqrt(ovl(ii))
            if (abs(ovl(ii)) <= stop_limit .and. endi(ii)==0 ) then
                endi(ii)=i
            end if                
            if (ovl(ii)/=real(0,rk))ovl(ii)=real(1,rk)/sqrt(ovl(ii))
            vector(1:dim,ii,i+1)=JJcoef(1:dim,ii)*ovl(ii)
            end do
!$OMP END PARALLEL        
            if (minval(endi(1:mxii))/=0) exit
        end do
!$OMP PARALLEL 
!$OMP DO PRIVATE (ii,ix)       
        do ii=1,mxii
        if (endi(ii)==0) endi(ii)=Tdim
!tridiagonalisation        
        tqli_iter=max_tqli        
        call tqli(alpha(:,ii),beta(:,ii),endi(ii),max_lanc,eigenv(:,:,ii),tqli_iter)
        alpha(1:endi(ii),ii)=abs(alpha(1:endi(ii),ii)+fact)
!        call indexx(i,alpha(1:endi(ii),ii),indx(1:endi(ii),ii))
        indx(1,ii)=maxval(minloc(alpha(1:endi(ii),ii)))
        if (output_control > 10) write(10,*) alpha(1:endi(ii),ii),'f',fact
        coef(1:dim,ii)=real(0,rk)
!calculate lowest eigenvector
        do ix = 1,endi(ii)
            coef(1:dim,ii)=coef(1:dim,ii)+eigenv(ix,indx(1,ii),ii)*vector(1:dim,ii,ix)
        end do
        errT(ii)=alpha(indx(1,ii),ii)
        end do
!$OMP END PARALLEL        
100     return        
   end subroutine Lanczos
   
   subroutine PJmJp(dim,dimJp,mxii)
        implicit none
        integer,intent(in)::dim,dimJp,mxii
        integer ii,k
        real(kind=rk)::Jp4,JJ4
        real(kind=rk)::Jp2,JJ2
        real(kind=rk)::Jp1,JJ1
        real(kind=rk)::Jp3,JJ3
        integer(kind=1)::M
        integer::I
        integer::iand,ishft
        
        if (mxii==4) then
!apply J+   
!$OMP PARALLEL
!$OMP DO PRIVATE(ii,Jp4,Jp1,Jp2,Jp3,k,M,I)            
            do ii=1,dimJp
                Jp4=real(0,rk)
                Jp1=real(0,rk)
                Jp2=real(0,rk)
                Jp3=real(0,rk)
                do k=1,mxopJp
                    I=inxJp(k,ii)
                    if (I == 0)  exit
                    M=ishft(I,-24)
                    I=iand(mask,I)
                    Jp1=Jp1+M*coef0(I,1)
                    Jp2=Jp2+M*coef0(I,2)
                    Jp3=Jp3+M*coef0(I,3)
                    Jp4=Jp4+M*coef0(I,4)
                end do
                Jpcoef(ii,1)=Jp1
                Jpcoef(ii,2)=Jp2
                Jpcoef(ii,3)=Jp3
                Jpcoef(ii,4)=Jp4
            end do
!$OMP BARRIER
!apply J-             
!$OMP DO PRIVATE(ii,JJ4,JJ1,JJ2,JJ3,k,I,M)            
            do ii=1,dim
                JJ4=real(0,rk)
                JJ1=real(0,rk)
                JJ2=real(0,rk)
                JJ3=real(0,rk)
                do k=1,mxopJm
                    I=inxJm(k,ii)
                    if (I == 0)  exit
                    M=ishft(I,-24)
                    I=iand(mask,I)
                    JJ1=JJ1+M*Jpcoef(I,1)
                    JJ2=JJ2+M*Jpcoef(I,2)
                    JJ3=JJ3+M*Jpcoef(I,3)
                    JJ4=JJ4+M*Jpcoef(I,4)
                end do
                JJcoef(ii,1)=JJ1
                JJcoef(ii,2)=JJ2
                JJcoef(ii,3)=JJ3
                JJcoef(ii,4)=JJ4
            end do
!$OMP END PARALLEL       
            else if (mxii==3) then
!apply J+   
!$OMP PARALLEL
!$OMP DO PRIVATE(ii,Jp1,Jp2,Jp3,k,M,I)            
            do ii=1,dimJp
                Jp1=real(0,rk)
                Jp2=real(0,rk)
                Jp3=real(0,rk)
                do k=1,mxopJp
                    I=inxJp(k,ii)
                    if (I == 0)  exit
                    M=ishft(I,-24)
                    I=iand(mask,I)
                    Jp1=Jp1+M*coef0(I,1)
                    Jp2=Jp2+M*coef0(I,2)
                    Jp3=Jp3+M*coef0(I,3)
                end do
                Jpcoef(ii,1)=Jp1
                Jpcoef(ii,2)=Jp2
                Jpcoef(ii,3)=Jp3
            end do
!$OMP BARRIER
!apply J-             
!$OMP DO PRIVATE(ii,JJ1,JJ2,JJ3,k,I,M)            
            do ii=1,dim
                JJ1=real(0,rk)
                JJ2=real(0,rk)
                JJ3=real(0,rk)
                do k=1,mxopJm
                    I=inxJm(k,ii)
                    if (I == 0)  exit
                    M=ishft(I,-24)
                    I=iand(mask,I)
                    JJ1=JJ1+M*Jpcoef(I,1)
                    JJ2=JJ2+M*Jpcoef(I,2)
                    JJ3=JJ3+M*Jpcoef(I,3)
                end do
                JJcoef(ii,1)=JJ1
                JJcoef(ii,2)=JJ2
                JJcoef(ii,3)=JJ3
            end do
!$OMP END PARALLEL       
        else if (mxii==2) then
!apply J+   
!$OMP PARALLEL
!$OMP DO PRIVATE(ii,Jp1,Jp2,k,M,I)            
            do ii=1,dimJp
                Jp1=real(0,rk)
                Jp2=real(0,rk)
                do k=1,mxopJp
                    I=inxJp(k,ii)
                    if (I == 0)  exit
                    M=ishft(I,-24)
                    I=iand(mask,I)
                    Jp1=Jp1+M*coef0(I,1)
                    Jp2=Jp2+M*coef0(I,2)
                end do
                Jpcoef(ii,1)=Jp1
                Jpcoef(ii,2)=Jp2
            end do
!$OMP BARRIER
!apply J-             
!$OMP DO PRIVATE(ii,JJ1,JJ2,k,I,M)            
            do ii=1,dim
                JJ1=real(0,rk)
                JJ2=real(0,rk)
                do k=1,mxopJm
                    I=inxJm(k,ii)
                    if (I == 0)  exit
                    M=ishft(I,-24)
                    I=iand(mask,I)
                    JJ1=JJ1+M*Jpcoef(I,1)
                    JJ2=JJ2+M*Jpcoef(I,2)
                end do
                JJcoef(ii,1)=JJ1
                JJcoef(ii,2)=JJ2
            end do
!$OMP END PARALLEL       
        else if (mxii==1) then
!apply J+   
!$OMP PARALLEL
!$OMP DO PRIVATE(ii,Jp1,k,M,I)            
            do ii=1,dimJp
                Jp1=real(0,rk)
                do k=1,mxopJp
                    I=inxJp(k,ii)
                    if (I == 0)  exit
                    M=ishft(I,-24)
                    I=iand(mask,I)
                    Jp1=Jp1+M*coef0(I,1)
                end do
                Jpcoef(ii,1)=Jp1
            end do
!$OMP BARRIER
!apply J-             
!$OMP DO PRIVATE(ii,JJ1,k,I,M)            
            do ii=1,dim
                JJ1=real(0,rk)
                do k=1,mxopJm
                    I=inxJm(k,ii)
                    if (I == 0)  exit
                    M=ishft(I,-24)
                    I=iand(mask,I)
                    JJ1=JJ1+M*Jpcoef(I,1)
                end do
                JJcoef(ii,1)=JJ1
            end do
!$OMP END PARALLEL       
        end if                                 
    end subroutine PJmJp
    
   subroutine PJmJpr(dim,dimJp,mxii)
        implicit none
        integer,intent(in)::dim,dimJp,mxii
        integer ii,k
        real(kind=rk)::Jp4,JJ4
        real(kind=rk)::Jp2,JJ2
        real(kind=rk)::Jp1,JJ1
        real(kind=rk)::Jp3,JJ3
        real(kind=rk)::M
        integer::I
        
        if (mxii==4) then
!apply J+   
!$OMP PARALLEL
!$OMP DO PRIVATE(ii,Jp4,Jp1,Jp2,Jp3,k,M,I)            
            do ii=1,dimJp
                Jp4=real(0,rk)
                Jp1=real(0,rk)
                Jp2=real(0,rk)
                Jp3=real(0,rk)
                do k=1,mxopJp
                    I=inxJp(k,ii)
                    if (I == 0)  exit
                    if (matxJp(k,ii)<0) then
                        M=-sqrtf(matxJp(k,ii))*sqrtf(matxJpr(k,ii))
                    else
                        M=sqrtf(matxJp(k,ii))*sqrtf(matxJpr(k,ii))
                    end if
                    Jp1=Jp1+M*coef0(I,1)
                    Jp2=Jp2+M*coef0(I,2)
                    Jp3=Jp3+M*coef0(I,3)
                    Jp4=Jp4+M*coef0(I,4)
                end do
                Jpcoef(ii,1)=Jp1
                Jpcoef(ii,2)=Jp2
                Jpcoef(ii,3)=Jp3
                Jpcoef(ii,4)=Jp4
            end do
!$OMP BARRIER
!apply J-             
!$OMP DO PRIVATE(ii,JJ4,JJ1,JJ2,JJ3,k,I,M)            
            do ii=1,dim
                JJ4=real(0,rk)
                JJ1=real(0,rk)
                JJ2=real(0,rk)
                JJ3=real(0,rk)
                do k=1,mxopJm
                    I=inxJm(k,ii)
                    if (I == 0)  exit
                    if (matxJm(k,ii)<0) then
                        M=-sqrtf(matxJm(k,ii))*sqrtf(matxJmr(k,ii))
                    else
                        M=sqrtf(matxJm(k,ii))*sqrtf(matxJmr(k,ii))
                    end if
                    JJ1=JJ1+M*Jpcoef(I,1)
                    JJ2=JJ2+M*Jpcoef(I,2)
                    JJ3=JJ3+M*Jpcoef(I,3)
                    JJ4=JJ4+M*Jpcoef(I,4)
                end do
                JJcoef(ii,1)=JJ1
                JJcoef(ii,2)=JJ2
                JJcoef(ii,3)=JJ3
                JJcoef(ii,4)=JJ4
            end do
!$OMP END PARALLEL       
        else if (mxii==3) then
!apply J+   
!$OMP PARALLEL
!$OMP DO PRIVATE(ii,Jp1,Jp2,Jp3,k,M,I)            
            do ii=1,dimJp
                Jp1=real(0,rk)
                Jp2=real(0,rk)
                Jp3=real(0,rk)
                do k=1,mxopJp
                    I=inxJp(k,ii)
                    if (I == 0)  exit
                    if (matxJp(k,ii)<0) then
                        M=-sqrtf(matxJp(k,ii))*sqrtf(matxJpr(k,ii))
                    else
                        M=sqrtf(matxJp(k,ii))*sqrtf(matxJpr(k,ii))
                    end if
                    Jp1=Jp1+M*coef0(I,1)
                    Jp2=Jp2+M*coef0(I,2)
                    Jp3=Jp3+M*coef0(I,3)
                end do
                Jpcoef(ii,1)=Jp1
                Jpcoef(ii,2)=Jp2
                Jpcoef(ii,3)=Jp3
            end do
!$OMP BARRIER
!apply J-             
!$OMP DO PRIVATE(ii,JJ1,JJ2,JJ3,k,I,M)            
            do ii=1,dim
                JJ1=real(0,rk)
                JJ2=real(0,rk)
                JJ3=real(0,rk)
                do k=1,mxopJm
                    I=inxJm(k,ii)
                    if (I == 0)  exit
                    if (matxJm(k,ii)<0) then
                        M=-sqrtf(matxJm(k,ii))*sqrtf(matxJmr(k,ii))
                    else
                        M=sqrtf(matxJm(k,ii))*sqrtf(matxJmr(k,ii))
                    end if
                    JJ1=JJ1+M*Jpcoef(I,1)
                    JJ2=JJ2+M*Jpcoef(I,2)
                    JJ3=JJ3+M*Jpcoef(I,3)
                end do
                JJcoef(ii,1)=JJ1
                JJcoef(ii,2)=JJ2
                JJcoef(ii,3)=JJ3
            end do
!$OMP END PARALLEL       
        else if (mxii==2) then
!apply J+   
!$OMP PARALLEL
!$OMP DO PRIVATE(ii,Jp1,Jp2,k,M,I)            
            do ii=1,dimJp
                Jp1=real(0,rk)
                Jp2=real(0,rk)
                do k=1,mxopJp
                    I=inxJp(k,ii)
                    if (I == 0)  exit
                    if (matxJp(k,ii)<0) then
                        M=-sqrtf(matxJp(k,ii))*sqrtf(matxJpr(k,ii))
                    else
                        M=sqrtf(matxJp(k,ii))*sqrtf(matxJpr(k,ii))
                    end if
                    Jp1=Jp1+M*coef0(I,1)
                    Jp2=Jp2+M*coef0(I,2)
                end do
                Jpcoef(ii,1)=Jp1
                Jpcoef(ii,2)=Jp2
            end do
!$OMP BARRIER
!apply J-             
!$OMP DO PRIVATE(ii,JJ1,JJ2,k,I,M)            
            do ii=1,dim
                JJ1=real(0,rk)
                JJ2=real(0,rk)
                do k=1,mxopJm
                    I=inxJm(k,ii)
                    if (I == 0)  exit
                    if (matxJm(k,ii)<0) then
                        M=-sqrtf(matxJm(k,ii))*sqrtf(matxJmr(k,ii))
                    else
                        M=sqrtf(matxJm(k,ii))*sqrtf(matxJmr(k,ii))
                    end if
                    JJ1=JJ1+M*Jpcoef(I,1)
                    JJ2=JJ2+M*Jpcoef(I,2)
                end do
                JJcoef(ii,1)=JJ1
                JJcoef(ii,2)=JJ2
            end do
!$OMP END PARALLEL       
        else if (mxii==1) then
!apply J+   
!$OMP PARALLEL
!$OMP DO PRIVATE(ii,Jp1,k,M,I)            
            do ii=1,dimJp
                Jp1=real(0,rk)
                do k=1,mxopJp
                    I=inxJp(k,ii)
                    if (I == 0)  exit
                    if (matxJp(k,ii)<0) then
                        M=-sqrtf(matxJp(k,ii))*sqrtf(matxJpr(k,ii))
                    else
                        M=sqrtf(matxJp(k,ii))*sqrtf(matxJpr(k,ii))
                    end if
                    Jp1=Jp1+M*coef0(I,1)
                end do
                Jpcoef(ii,1)=Jp1
            end do
!$OMP BARRIER
!apply J-             
!$OMP DO PRIVATE(ii,JJ1,k,I,M)            
            do ii=1,dim
                JJ1=real(0,rk)
                do k=1,mxopJm
                    I=inxJm(k,ii)
                    if (I == 0)  exit
                    if (matxJm(k,ii)<0) then
                        M=-sqrtf(matxJm(k,ii))*sqrtf(matxJmr(k,ii))
                    else
                        M=sqrtf(matxJm(k,ii))*sqrtf(matxJmr(k,ii))
                    end if
                    JJ1=JJ1+M*Jpcoef(I,1)
                end do
                JJcoef(ii,1)=JJ1
            end do
!$OMP END PARALLEL       
        end if
                             
    end subroutine PJmJpr
    
    subroutine PTmTp(dim,dimTp,mxii)
        implicit none
        integer,intent(in)::dim,dimTp,mxii
        integer ii,k
        real(kind=rk)::Jp4,JJ4
        real(kind=rk)::Jp2,JJ2
        real(kind=rk)::Jp1,JJ1
        real(kind=rk)::Jp3,JJ3
        integer::I
        
        if (mxii==4) then
!apply T+   
!$OMP PARALLEL
!$OMP DO PRIVATE(ii,Jp4,Jp1,Jp2,Jp3,k,I)            
            do ii=1,dimTp
                Jp4=real(0,rk)
                Jp1=real(0,rk)
                Jp2=real(0,rk)
                Jp3=real(0,rk)
                do k=1,mxopTp
                    I=inxTp(k,ii)
                    if (I == 0)  exit
                    Jp1=Jp1+coef0(I,1)
                    Jp2=Jp2+coef0(I,2)
                    Jp3=Jp3+coef0(I,3)
                    Jp4=Jp4+coef0(I,4)
                end do
                Jpcoef(ii,1)=Jp1
                Jpcoef(ii,2)=Jp2
                Jpcoef(ii,3)=Jp3
                Jpcoef(ii,4)=Jp4
            end do
!$OMP BARRIER
!apply T-             
!$OMP DO PRIVATE(ii,JJ4,JJ1,JJ2,JJ3,k,I)            
            do ii=1,dim
                JJ4=real(0,rk)
                JJ1=real(0,rk)
                JJ2=real(0,rk)
                JJ3=real(0,rk)
                do k=1,mxopTm
                    I=inxTm(k,ii)
                    if (I == 0)  exit
                    JJ1=JJ1+Jpcoef(I,1)
                    JJ2=JJ2+Jpcoef(I,2)
                    JJ3=JJ3+Jpcoef(I,3)
                    JJ4=JJ4+Jpcoef(I,4)
                end do
                JJcoef(ii,1)=JJ1
                JJcoef(ii,2)=JJ2
                JJcoef(ii,3)=JJ3
                JJcoef(ii,4)=JJ4
            end do
!$OMP END PARALLEL       
        else if (mxii==3) then     
!apply T+   
!$OMP PARALLEL
!$OMP DO PRIVATE(ii,Jp1,Jp2,Jp3,k,I)            
            do ii=1,dimTp
                Jp1=real(0,rk)
                Jp2=real(0,rk)
                Jp3=real(0,rk)
                do k=1,mxopTp
                    I=inxTp(k,ii)
                    if (I == 0)  exit
                    Jp1=Jp1+coef0(I,1)
                    Jp2=Jp2+coef0(I,2)
                    Jp3=Jp3+coef0(I,3)
                end do
                Jpcoef(ii,1)=Jp1
                Jpcoef(ii,2)=Jp2
                Jpcoef(ii,3)=Jp3
            end do
!$OMP BARRIER
!apply T-             
!$OMP DO PRIVATE(ii,JJ1,JJ2,JJ3,k,I)            
            do ii=1,dim
                JJ1=real(0,rk)
                JJ2=real(0,rk)
                JJ3=real(0,rk)
                do k=1,mxopTm
                    I=inxTm(k,ii)
                    if (I == 0)  exit
                    JJ1=JJ1+Jpcoef(I,1)
                    JJ2=JJ2+Jpcoef(I,2)
                    JJ3=JJ3+Jpcoef(I,3)
                end do
                JJcoef(ii,1)=JJ1
                JJcoef(ii,2)=JJ2
                JJcoef(ii,3)=JJ3
            end do
!$OMP END PARALLEL       
        else if (mxii==2) then
!apply T+   
!$OMP PARALLEL
!$OMP DO PRIVATE(ii,Jp1,Jp2,k,I)            
            do ii=1,dimTp
                Jp1=real(0,rk)
                Jp2=real(0,rk)
                do k=1,mxopTp
                    I=inxTp(k,ii)
                    if (I == 0)  exit
                    Jp1=Jp1+coef0(I,1)
                    Jp2=Jp2+coef0(I,2)
                end do
                Jpcoef(ii,1)=Jp1
                Jpcoef(ii,2)=Jp2
            end do
!$OMP BARRIER
!apply T-             
!$OMP DO PRIVATE(ii,JJ1,JJ2,k,I)            
            do ii=1,dim
                JJ1=real(0,rk)
                JJ2=real(0,rk)
                do k=1,mxopTm
                    I=inxTm(k,ii)
                    if (I == 0)  exit
                    JJ1=JJ1+Jpcoef(I,1)
                    JJ2=JJ2+Jpcoef(I,2)
                end do
                JJcoef(ii,1)=JJ1
                JJcoef(ii,2)=JJ2
            end do
!$OMP END PARALLEL       
        else if (mxii==1) then
!apply T+   
!$OMP PARALLEL
!$OMP DO PRIVATE(ii,Jp1,k,I)            
            do ii=1,dimTp
                Jp1=real(0,rk)
                do k=1,mxopTp
                    I=inxTp(k,ii)
                    if (I == 0)  exit
                    Jp1=Jp1+coef0(I,1)
                    Jp2=Jp2+coef0(I,2)
                end do
                Jpcoef(ii,1)=Jp1
            end do
!$OMP BARRIER
!apply T-             
!$OMP DO PRIVATE(ii,JJ1,k,I)            
            do ii=1,dim
                JJ1=real(0,rk)
                do k=1,mxopTm
                    I=inxTm(k,ii)
                    if (I == 0)  exit
                    JJ1=JJ1+Jpcoef(I,1)
                    JJ2=JJ2+Jpcoef(I,2)
                end do
                JJcoef(ii,1)=JJ1
            end do
!$OMP END PARALLEL       
        end if                        
    end subroutine PTmTp
    
end module Project        

Module OutputCoef

use Parameters
use Project

implicit none

contains

    subroutine coef_write(part,pindx,ngood,j,t,k,dim,base,baser,baseb)
        implicit none
        type(spartition):: part
        integer:: pindx,ngood,j,t,dim,k,base,baser,baseb,i
        integer:: ng,ngu
        logical:: first_call=.true.
        real,save::x
        if (first_call) call random_number(x)
        first_call=.false.
        write(base+k) pindx,ngood,j,t,dim,x       
!        if (output_control > 2) write(10,*) '**********ouput**********'     
!        if (output_control > 2) write(10,*) pindx,ngood,j,t,dim,x        
        write(base+k) part
!        if (output_control > 2) write(10,*) part        
        write(baser+k) pindx,ngood,j,t,dim,x        
        write(baser+k) part
        write(baseb+k) pindx,ngood,j,t,dim,x       
        write(baseb+k) part
        if (ngood  < 1 .or. dim < 1) return
        
            write(baseb+k,err=10) (mvecr(igood(i))%fn(1:gwords()),i=1,ngood)
            do ng=1,ngood
                ngu=ng
                if (disk) then
                    read(66,rec=ng) coef_st(1:dim,1)
                    ngu=1
                end if
                write(base+k,err=20) real(coef_st(1:dim,ngu),rkw)
                if (output_control > 22) write(10,*) pindx,real(coef_st(1:dim,ngu),rkw)
                write(baser+k,err=30) real(coef_st(igood(ng:ngood),ngu),rko)
                if (output_control > 20) write(10,*) real(coef_st(igood(ng:ngood),ngu),rko)
            end do
            return
10          print *, ' Error writing jba file'
            stop            
20          print *, ' Error writing prj file'
            stop            
30          print *, ' Error writing ort file'
            stop            
        
   end subroutine coef_write
   
end module OutputCoef

!****************************************************************************
!  NuProj.f90 
!
!  FUNCTIONS:
!  NuProj      - Entry point of console application.
!
!****************************************************************************

    program NuProj
    
    use Parameters
    use InputNuShell
    use OutputNuShell
    use Shells
    use ShellsX
    use Mscheme
    use Files
    use Project
    use OutputCoef
    use ProjAccuracy
    use IntegerProject
    use Extra
    
    implicit none

    ! Variables
    
    real:: t1,t2,et1,et2,ms
    integer:: hr,mn,sc
    integer:: pindx,j,t,k,jmx,tmx,nptt,npt,dim,nJt,jj,tt,twrite,k1,length,dJ
    integer:: pindxJp,jmxJp,tmxJp,dimJp,nJtJp
    integer:: pindxTp,jmxTp,tmxTp,dimTp,nJtTp
    type(spartition):: part,partJp,partTp
    integer,dimension(max_no_cJstates*2,max_no_cTstates*2):: tno_cJTstates
    
    character(len=10):: btime,ptime
    character(len=8):: bdate,pdate
    integer:: res,noJT
    character(len=3):: chnjt='NJT'
    character(len=10):: chnojt
    integer::OMP_get_max_threads
    
    ! Body of NuProj
        
    res=get_env_v(chnjt,chnojt)
    if (res /= 0) then
        read(chnojt,*) noJT
    else
        noJT=max_nJT
    end if
    allocate(mvec_pred(noJT))
    
    call output_header(6)
    call output_NuProj(6)

    call input_nucleus(inExt) 
    call output_welcome(6)    
    call input(inExt)
    call output_shells(6)    
    call check_shell_order()
    call setup_major_shells()
    call setup_shells()
    call read_parameters()
    if (use_isospin) then
        dt=2
    else
        dt=0
    end if
    dJ=0
    if (min_cJ2<0) dJ=2
    call setup_mscheme()
    call setup_lu()
    call setup_factorial()
    call setup_file_ext()
    
    call output_end_setup(6)
    call cpu_time(t1)
    call date_and_time(bdate,btime)
    tno_cJTstates=0
    call open_files_read(nptt,'.nba',100,dt)
    call open_files_read(nptt,'.nba',1100,dt)
    call open_files_read(nptt,'.nba',1200,dt)
    call open_files(nptt,'.prj',200,dt)
    call open_files(nptt,'.ort',300,dt)
    call open_files(nptt,'.Jba',400,dt)
    if (output_control > 0) open (unit=10, file=Nucleus//'NuProj'//'.txt', action='write')

    k1=0
    k=0
    jj=-1
    do j=min_cJ2,max_cJ2+2,2
    jj=jj+2
    tt=-1
        do t=min_cT2,max_cT2+dt,2
            tt=tt+2
            do npt=1,nptt
                if (t <= max_cT2 .and. j <= max_cJ2) then
                    call read_mvec(part,pindx,k1,1,dim,jmx,tmx,nJT,100)
                    if (output_control > 0)write(10,*) ' 2*J , 2*T  ',j,t
                    if (output_control > 0)write(10,*)
                    call read_mvec(partJp,pindxJp,k1+(dt+max_cT2-min_cT2)/2+1,&
                                          2,dimJp,jmxJp,tmxJp,nJTJp,1100)
                    if (use_isospin) then
                        call read_mvec(partTp,pindxTp,k1+1,3,dimTp,jmxTp,tmxTp,nJTTp,1200)
                    end if
                    if (nJT > noJT) then
                        print *, ' Increase max_nJT in module ProjParameters or'
                        print *, ' SET NJT=I where I >>',nJT
                        stop
                    end if
                    if (nJT > 0 .and. dim > 0) then
                        allocate(coef_st(dim,nJT),stat=err)
                        disk=.false.
                        if (err/=0) then
                            disk=.true.
                            allocate(coef_st(dim,1))
                            inquire(iolength=length) coef_st
                            open(unit=66,file='coefst.tmp',form='unformatted',access='direct',recl=length)
                        end if
                        call date_and_time(pdate,ptime)
                        if (output_control > 0) write(10,*) pdate,' ',ptime, &
                                  '  starting 2*J,2*T,part no',j,t,pindx
                        if (output_control>2) print *,part%shell(1:no_shells)
                        call do_projection(part,j,t,jmx,tmx,dim,nJT,dimJp,dimTp)
                    end if
                    call coef_write(part,pindx,nJT,j,t,k,dim,200,300,400)
                    close(unit=66)
                    deallocate(mvecr)
                    deallocate(mvecrJp)
                    if (use_isospin) deallocate(mvecrTp)
                    if (nJT > 0 .and. dim > 0) deallocate(coef_st)
                    tno_cJTstates(jj,tt)=tno_cJTstates(jj,tt)+nJT
                end if
            end do
            if (t <= max_cT2 .and. j <= max_cJ2) k=k+1
            k1=k1+1
        end do
    
100 end do
    
    call close_files(200,dt)
    call close_files(300,dt)
    call close_files(400,dt)
    
    if (output_control > 0) then
        print * ,' No of states for each 2*J'
        print '(a10,10i8)', '       2*J',(j,j=min_cJ2,min_cJ2+18,2)
        
        tt=-1
        do t=min_cT2,max_cT2,2
            tt=tt+2
            twrite=t
            if (.not.use_isospin) twrite=nmp
            print '(a4,i6,10i8)', '2*T=',twrite,(tno_cJTstates(j,tt),j=1,19,2)
        end do
    end if
    
!    deallocate(mvec_pred)
    
    read(btime,'(i2,i2,i2,f4.3)') hr,mn,sc,ms
    et1=hr*3600.+mn*60.+sc*1. +ms
    call date_and_time(bdate,btime)
    call cpu_time(t2)
    read(btime,'(i2,i2,i2,f4.3)') hr,mn,sc,ms
    et2=hr*3600.+mn*60.+sc*1. +ms
    et1=et2-et1
    if (et1<0.0) et1=et1+24.*3600.
    t1=t2-t1
    call output_completed(6,'NuOper     ')    
    call output_time('NuProj',t1,et1)
    print *
    
    end program NuProj

    subroutine read_mvec(part,pindx,k,index,ibf,Jmax,Tmax,nJT,base)
    use Parameters
    use InputNuShell
    use Shells
    use ShellsX
    use Mscheme
    use Files
    use Project
    use OutputCoef
    
        implicit none
        type(spartition),intent(out):: part
        integer,intent(in):: k,base,index
        integer,intent(out)::ibf,Jmax,Tmax,pindx,nJT
        integer:: j,t,i
        
           read(base+k,err=10,end=11) pindx,nJT,j,t,ibf,Jmax,Tmax 
           if (index==1) allocate(mvecr(ibf))
           if (index==2) allocate(mvecrJp(ibf))     
           if (index==3) allocate(mvecrTp(ibf))     
           read(base+k,err=10,end=11) part
           if (ibf /= 0) then
                if (index == 1) read(base+k,err=10,end=11) (mvecr(i)%fn(1:gwords()),i=1,ibf)
                if (index == 2) read(base+k,err=10,end=11) (mvecrJp(i)%fn(1:gwords()),i=1,ibf)
                if (index == 3) read(base+k,err=10,end=11) (mvecrTp(i)%fn(1:gwords()),i=1,ibf)
           end if
           return
10         print *, ' Error reading mvectors' 
           stop      
11         print *, ' EOF reading mvectors' 
           stop      
                    
     end subroutine read_mvec
     
    subroutine output_NuProj(dev)
     
    use InputNuShell
    use OutputNuShell
    use Shells
    use ShellsX
    use Mscheme
    use Files
    use Project
    use OutputCoef
    use ProjAccuracy
     
    implicit none
    integer,intent(in)::dev
    integer:: status
    character(len=14)::buffer
    call get_cmd_arg(1,buffer,status) 
    if (status /= -1) return
     
    write(dev,*) ' NuProj Program'
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
    write(dev,*) ' Maximum no of m-vectors per partition single m',max_proj_mvec
    write(dev,*) ' Maximum no J,T values',max_no_cJstates,max_no_cTstates
    write(dev,*) ' Real kind for projection', rk
    write(dev,*) ' Real kind for output files', rkw
    write(dev,*) ' Lanczos maximum iterations',max_lanc
    write(dev,*)
    
end subroutine output_NuProj    
