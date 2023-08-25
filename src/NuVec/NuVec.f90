module Vec

use Parameters
use Partition
use Mscheme
use ProjParameters
use MatrixParameters
use LanczosParameters
use OperParameters
use InputNuShell
use Files
use TriDiag
use TrdChar
use Mscheme

implicit none
integer,parameter:: noiter=300
character(len=7)::opfile1,opfile2,intfile
character(len=6)::opfile16,opfile26,Nucleuss
character(len=1):: xtraa,xtrab
integer,dimension(:,:),allocatable:: nJTa,nJTb
type(spartition),dimension(:),allocatable:: spart,sparta,spartb
integer:: min_cJ2a,max_cJ2a,no_sparta,no_spart
integer:: min_cJ2b,max_cJ2b,no_spartb,no_state=0
real(kind=rkw),dimension(:),allocatable:: coefR,coefT
real(kind=rkw),dimension(:),allocatable:: ava,avb
real(kind=rkw),dimension(:,:,:,:),allocatable:: spartpc
character(len=8),dimension(:),allocatable:: bas,baspt
integer:: ajl,ajr,anpl,anpr,ashl,ashr,lm,idx,lmx,lmn,lmnn,lmxx
integer:: bjl,bjr,bnpl,bnpr,J2n,J2x,ajrn,ajrx,ajln,ajlx
integer:: bjln,bjlx,bshl,bshr,bjrn,bjrx,no_level
integer:: lsc,psc,atsc,btsc,bsc,osc,nsc,hsc,basp
integer:: J2,nthr=2,rsc,ip,max_dbl,max_dar,mna,mnb,apar,bpar
integer:: bseg,aseg,no_sgmtb,no_sgmta,dimptl,dimptr,matdima,matdimb,no_cA,no_cB,kk
integer:: max_cJ2s,min_cJ2s,nmps,basd,mxa,mxb,nosa,nosb
logical:: done,new,proton,neutron
real:: xa,xb,xxa,xxb,ampl
type(spartition)::part_opa,part_opb

contains

    subroutine Vec_TRD(base)
    implicit none
    integer,intent(in)::base
    integer:: no_Rac=0
    real::t1,t2
    
    read(base+1) no_sparta,min_cJ2a,max_cJ2a,proton,neutron,xa,apar
    if (proton) typea='p'
    if (neutron) typea='n'
    allocate(nJTa(no_sparta,min_cJ2a/2:max_cJ2a/2))
    allocate(sparta(no_sparta))
    read(base+1) nJTa(1:no_sparta,min_cJ2a/2:max_cJ2a/2)
    read(base+1) sparta(1:no_sparta)
    
    
    read(base+2) no_spartb,min_cJ2b,max_cJ2b,proton,neutron,xb,bpar
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

    J2=min_cJ2
    
! setup nucleon numbers    
    nmps=nmp
    if (typea=='p') then
    nmp=(no_cN+nmps)/2
    no_cB=nmp
    nmp=(no_cN-nmps)/2
    no_cA=nmp
    else
    nmp=(no_cN-nmps)/2
    no_cB=nmp
    nmp=(no_cN+nmps)/2
    no_cA=nmp
    end if
! end set up nucleon numbers

    
    Nucleuss=Nucleus
    max_cJ2s=max_cJ2
    min_cJ2s=min_cJ2

  
! start lanczos iteration 
    call Vec_OP
    
                          
   end subroutine Vec_TRD
          
   subroutine Vec_OP
   implicit none
   
        real(kind=rl),dimension(noiter):: alpha,beta
        integer,dimension(noiter):: indx
        real(kind=rl),dimension(max_no_shells):: nav
        real(kind=rl),dimension(noiter,noiter):: eigenv,matrix
        real(kind=rl),dimension(noiter)::ovla
        integer:: n,i,imax,imin,ix,isave,np,pd,tqlim,status,idum,idum1,inc_itern=2
        integer:: count,no_keep,no_conv,second_iter,firsti,noit
        real:: x,seig2,tots,totp,totn,ti,tf,tt=0.0
        real(kind=rl):: ovl,evalue,zero=0.0,paramp,nparamp,tampl
        character(len=1),dimension(noiter):: chpar
        logical:: lnzf,trlf
        
        inquire(file=Nucleuss//file_code(J2/2,no_cA,no_cB)//'.lnz',exist=lnzf)
        inquire(file=Nucleuss//file_code(J2/2,no_cA,no_cB)//'.trl',exist=trlf)
        if (.not.(trlf.or.lnzf)) then
            print *,' Lnz or Trl file not found'
            stop
        end if
        open(unit=701,file=Nucleuss//file_code(J2/2,no_cA,no_cB)//'.abz',form='unformatted',action='read')
        if (lnzf) then
            read(701) isave,alpha(1:isave),beta(1:isave)
            read(701) imax,no_level,tt,no_state,basd
            allocate(bas(basd))
            read(701) bas
        else if (trlf) then
            read(701) isave,alpha(1:isave),beta(1:isave),matrix(1:isave,1:isave)
            read(701) imax,no_level,tt,no_state,basd,count,no_keep,no_conv,firsti,&
                        second_iter,noit
            allocate(bas(basd))
            read(701) bas
        end if  
! allocate and set vectors    
        allocate(coefR(no_state))
        allocate(coefT(no_state))
        allocate(baspt(no_state))
        nsc=0
        do bsc=1,basd
            ajl=ichar(bas(bsc)(1:1))
            bjl=ichar(bas(bsc)(2:2))
            do bnpl=1,no_spartb
            do rsc=1,nJTb(bnpl,bjl/2)
            do anpl=1,no_sparta
            do lsc=1,nJTa(anpl,ajl/2)
            nsc=nsc+1
            baspt(nsc)=bas(bsc)
            baspt(nsc)(3:4)=char2(anpl)
            baspt(nsc)(5:6)=char2(bnpl)
            end do
            end do
            end do
            end do
        end do
        allocate(spartpc(no_sparta,no_spartb,min_cJ2a/2:max_cJ2a/2,min_cJ2b/2:max_cJ2b/2))
        mna=min_shellsa()
        mxa=max_shellsa()
        mnb=min_shellsb()
        mxb=max_shellsb()
        allocate(ava(no_shells))
        allocate(avb(no_shells))
! end setup vectors
        if (lnzf) open(unit=700,file=Nucleuss//file_code(J2/2,no_cA,no_cB)//'.lnz',&
            form='unformatted',action='read',access='direct',recl=no_state*rkw/4)
        if (trlf) open(unit=700,file=Nucleuss//file_code(J2/2,no_cA,no_cB)//'.trl',&
                form='unformatted',action='read',access='direct',recl=no_state*rkw/4)
        open(unit=703,file=Nucleuss//file_code(J2/2,no_cA,no_cB)//'.xvc',form='unformatted',action='write')
        open(unit=705,file=Nucleuss//file_code(J2/2,no_cA,no_cB)//'.xva',form='unformatted',action='write')
        open(unit=702,file=Nucleuss//file_code(J2/2,no_cA,no_cB)//'.lev',form='formatted',action='write')
        open(unit=704,file=Nucleuss//file_code(J2/2,no_cA,no_cB)//'.lp',form='formatted',action='write')
        call output_results(702,J2,abs(no_cA-no_cB),0)
        call output_results(704,J2,abs(no_cA-no_cB),0)
        if (no_level>no_state) no_level=no_state
        if (lnzf) then
            eigenv=real(0,rl)
            do i=1,isave
                eigenv(i,i)=real(1,rl)
            end do
            tqlim=2000
            call tqli(alpha,beta,isave,noiter,eigenv,tqlim)
            call indexx(isave,alpha,indx)
        else if (trlf) then
            if (count==1) then
            do ix=isave+1,2,-1
                beta(ix)=beta(ix-1)
            end do
            beta(1)=real(0,rk)
            eigenv=real(0,rl)
            do i=1,isave
                eigenv(i,i)=real(1,rl)
            end do
            else if (count>1) then
                eigenv(1:isave,1:isave)=matrix(1:isave,1:isave)
                call tred2(eigenv,isave,noiter,alpha,beta)
            end if    
            tqlim=2000
            call tqli(alpha,beta,isave,noiter,eigenv,tqlim)
            call indexx(isave,alpha,indx)
        end if
        write(703) no_state,min(no_level,no_state,isave)
        write(704,*)
        write(704,*) ' Type a = ',typea,'  Type b = ',typeb
        write(704,*) ' A partitions'
        do anpl=1,no_sparta 
            write(704,'(i8,18i4)')  -anpl,sparta(anpl)%shell(1:no_shells)
        end do    
        write(704,*) ' B partitions'
        do bnpl=1,no_spartb
            write(704,'(i8,18i4)')  -bnpl,spartb(bnpl)%shell(1:no_shells)
        end do    
        do i=1,min(no_level,isave)
            coefT=real(0,rl)
            do ix = 1,isave
                read(700,rec=ix) coefR
                coefT=coefT+eigenv(ix,indx(i))*coefR
            end do
            ava=real(0,rkw)
            avb=real(0,rkw)
            spartpc=real(0,rkw)
            paramp=0.0
            nparamp=0.0
            tampl=0.0
            do bsc=1,no_state
                ampl=coefT(bsc)*coefT(bsc)
                tampl=tampl+ampl
                ajl=ichar(baspt(bsc)(1:1))
                bjl=ichar(baspt(bsc)(2:2))
                anpl=ichar2(baspt(bsc)(3:4))
                bnpl=ichar2(baspt(bsc)(5:6))
                if (cParity(sparta(anpl))*cParity(spartb(bnpl))>=0) paramp=paramp + ampl
                if (cParity(sparta(anpl))*cParity(spartb(bnpl))<0) nparamp=nparamp + ampl
                spartpc(anpl,bnpl,ajl/2,bjl/2)=spartpc(anpl,bnpl,ajl/2,bjl/2)+&
                                ampl
                do nsc=1,no_shells
                    ava(nsc)=ava(nsc)+sparta(anpl)%shell(nsc)*ampl
                end do
                do nsc=1,no_shells
                    avb(nsc)=avb(nsc)+spartb(bnpl)%shell(nsc)*ampl
                end do
            end do
            chpar(i)='?'
            if (paramp/tampl>0.95) chpar(i)='+'
            if (nparamp/tampl>0.95) chpar(i)='-'
            write(704,*)
            write(704,*) 'State ',i,alpha(indx(i)),' ',chpar(i)
            write(704,*)
            write(704,*) 'Average A type nucleons  '
            write(704,'(11f7.2)') ava(1:no_shells)/tampl
            write(704,*)
            write(704,*) 'Average B type nucleons  '
            write(704,'(11f7.2)') avb(1:no_shells)/tampl
            write(704,*) ' Positive and Negative Parity '
            write(704,'(a6,f8.2,a7,f8.2,a1)') '  +ve ',paramp*100.0/tampl,&
                            '% -ve ',nparamp*100.0/tampl,'%'
            spartpc=spartpc*100./tampl
            do bjl=min_cJ2b,max_cJ2b,2
            do ajl=min_cJ2a,max_cJ2a,2
            ampl=sum(spartpc(:,:,ajl/2,bjl/2))
            if (ampl>=1.0) then
            write(704,*)
                write(704,*) '2*JA, 2*JB, amplitude ',ajl,bjl,ampl
                do bnpl=1,no_spartb
                ampl=sum(spartpc(:,bnpl,ajl/2,bjl/2))
                if (ampl>=1.0) then
                    write(704,*) 'B partition, amplitude ',bnpl,ampl
                    write(704,*) 'A partitions, amplitudes'
                    write(704,'(6(a1,i3,a1,f7.2,a1))') ('(',anpl,',',&
                           spartpc(anpl,bnpl,ajl/2,bjl/2),')',anpl=1,no_sparta)
                    do anpl=1,no_sparta
                        if (spartpc(anpl,bnpl,ajl/2,bjl/2)>1.0) then
                            if (cParity(sparta(anpl))*cParity(spartb(bnpl))>0) then
                                !chpar(i)='+'
                            else
                                !chpar(i)='-'
                            end if
                            exit
                        end if
                    end do   
                end if 
                end do
            end if      
            end do
            end do                        
            write(705) alpha(indx(i)),ava,avb,spartpc
            write(703) alpha(indx(i)),coefT
   
        end do    
        write(702+0,'(5(f10.3,a1))') &
               (alpha(indx(ix)),chpar(ix),ix=1,min(no_level,isave))
        write(704+0,'(5(f10.3,a1))') &
               (alpha(indx(ix)),chpar(ix),ix=1,min(no_level,isave))
        write(6,*)
        write(6,'(5(f10.3,a1))') &
               (alpha(indx(ix)),chpar(ix),ix=1,min(no_level,isave))
        write(6,*)
        
   end subroutine Vec_OP
   
   subroutine output_results(dev,j,t,k)
    
    use Shells
    use InputNuShell
    use OutputNuShell
    
    implicit none
    integer,intent(in):: dev,j,t,k
    integer:: twrite

    call output_headerX(dev+k)    
    write(dev+k,*) ' '//intfile
    write(dev+k,*)
    write(dev+k,*) ' ',Description(1)
    write(dev+k,*) ' ',Description(2)
    write(dev+k,*) ' ',Description(3)
    write(dev+k,*)
    twrite=t
    write(dev+k,*) ' 2*J,   2*T',j,twrite
    write(dev+k,*) 
    write(dev+k,*) ' Lowest Energy Levels'
    write(dev+k,*) 
    
    end subroutine output_results

                            
end module Vec

!  NuVec.f90 
!
!  FUNCTIONS:
!  NuVec      - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: NuVec
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program NuVec
    
    use Vec
    use OperParameters
    use ProjParameters
    use InputNuShell
    use OutputNuShell
    use ShellsX
    use Mscheme
    use Files
    use Extra
    
    implicit none

    ! Variables
    
    real:: t1,t2,et1,et2,ms
    integer:: hr,mn,sc

    integer:: j,tmain,k,dt,k1,res,i
    
    character(len=6):: inExt
    character(len=10):: btime
    character(len=8):: bdate
    character(len=10)::chnopart,chpart='PARTITIONS'
    
    ! Body of NuVec
        
    
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
    call setup_lu()
    call setup_file_ext()
    call input_intfile(opfile1,4)
    call input_intfile(opfile2,5)
    call input_intfile(intfile,6)
    call output_end_setup(6)        
    call cpu_time(t1)
    call date_and_time(bdate,btime)
    if (output_control > 6) open(unit=60,file=Nucleus//'trd.txt',action='write')
    if (output_control > 0) open(unit=10,file=Nucleus//'NuVec'//'.txt',action='write')

    open(unit=801,file=opfile1//'.tri',form='unformatted',action='read')
    open(unit=802,file=opfile2//'.tri',form='unformatted',action='read')

    call Vec_TRD(800)
    close(unit=802)
    close(unit=801)
    read(btime,'(i2,i2,i2,f4.3)') hr,mn,sc,ms
    et1=hr*3600.+mn*60.+sc*1. +ms
    call date_and_time(bdate,btime)
    call cpu_time(t2)
    read(btime,'(i2,i2,i2,f4.3)') hr,mn,sc,ms
    et2=hr*3600.+mn*60.+sc*1. +ms
    t1=t2-t1
    et1=et2-et1
    if (et1<0.0) et1=et1+24.*3600.
    
    call output_completed(6,'NuTra      ')    
    call output_time('NuVec',t1,et1)
    print *
    
    end program NuVec
    

