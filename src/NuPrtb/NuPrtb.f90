!  NuPrtb.f90 
!
!  FUNCTIONS:
!  NuPrtb      - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: NuPrtb
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program NuPrtb
    
    use TriDiag
    use ProjParameters
    use OutputNuShell
    
    implicit none

    ! Variables
    real(kind=rkw),dimension(:,:),allocatable::  matrix
    real(kind=rkw),dimension(:),allocatable::  alpha,beta,diagm
    integer,dimension(:),allocatable::  indx
    integer::  ir,tqlim,dim
    character(len=16):: filename
    logical:: diag
    ! Body of NuPrtb
    
    write(6,'(a18)',advance='no')  '  Matrix Filename '
    read(5,'(a16)') filename
    filename=adjustr(filename)
    open(unit=10,file=filename//'.ovx',form='unformatted')
    inquire(file=filename//'.ovd',exist=diag)
    read(10) dim
    allocate(matrix(dim,dim))
    allocate(alpha(dim))
    allocate(beta(dim))
    allocate(diagm(dim))
    allocate(indx(dim))
    if (diag) then
        open(unit=11,file=filename//'.ovd')
        read(11,*) diagm
        write(6,*) ' Using '//filename//'.ovd.'
        write(6,*) ' Adding diagonal matrix elements.'
    else
        write(6,*)   
        write(6,*) ' No '//filename//'.ovd text file found.'
        write(6,*) ' No additional diagonal terms added.'
        diagm=0.0
    end if           
    do ir=1,dim
        read(10) matrix(:,ir)
        matrix(ir,ir)=matrix(ir,ir)+diagm(ir)
    end do
    call tred2(matrix,dim,dim,alpha,beta)
    tqlim=2000
    call tqli(alpha,beta,dim,dim,matrix,tqlim)
    call indexx(dim,alpha,indx)      
    write(6,*)   
    write(6,*) ' Eigenvalues:'      
    write(6,*)   
    write(6,'(5f8.3)')  alpha(indx(1:dim))
    open(unit=20,file=filename//'.lp',action='write')
    call output_headerX(20)
    write(6,*)   
    write(6,*) ' Output written to ',filename//'.lp'
    write(6,*)   
    write(20,*)   
    write(20,*) ' Eigenvalues:'     
    write(20,*)   
    
    write(20,'(5f8.3)')  alpha(indx(1:dim))
    write(20,*)   
    do ir=1,dim
        write(20,*)
        write(20,*)  ' State',ir,' E =',alpha(indx(ir)),' MeV'
        write(20,*)
        write(20,'(10f8.4)') matrix(:,indx(ir))
    end do    
   

    end program NuPrtb

