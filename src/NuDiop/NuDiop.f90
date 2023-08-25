!  NuDiop.f90 
!
!  FUNCTIONS:
!  NuDiop      - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: NuDiop
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program NuDiop
    use Tridiag
    implicit none

    ! Variables
    integer,parameter:: maxshell2=800
    integer:: max_op,max_lmo,i,no_Op,lm,cl,cr,j,k,max_shl,min_shl,max_shr,min_shr
    integer:: indxl,indxr,indxcl,indxcr,tqlim,l,pt,nonz,ksave,nsp,mm
    character(len=8),dimension(:),allocatable:: opn
    real(kind=8),dimension(:),allocatable:: meV
    integer,dimension(:),allocatable:: pl,pr
    real(kind=8),dimension(maxshell2,maxshell2):: matrix,mat
    real(kind=8),dimension(maxshell2):: alpha,beta
    real(kind=8)::rcut,root2,fact,rtlm,cutm
    integer,dimension(maxshell2):: indx,mapl,mapr,sp1,sp2,sp3,sp4
    character(len=8):: filein,fileout
    character(len=15):: cut,mcut
    character(len=2),dimension(maxshell2):: chl,chr,chrc
    real,dimension(0:17)::  cutlmp,cutlmn
    logical:: cutf
    ! Body of NuDiop
    root2=sqrt(2.0)
    max_shl=0
    max_shr=0
    min_shl=999
    min_shr=999
    no_Op=0
    cutlmp=0.0
    cutlmn=0.0
    call get_command_argument(1,filein)
    call get_command_argument(2,fileout)
    fileout=adjustr(fileout)
    call get_command_argument(3,cut)
    call get_command_argument(4,mcut)
    read(cut,*) rcut
    open(unit=50,file=fileout//'.enp',action='write')
    open(unit=60,file=fileout//'.vnp',action='write')
    read(mcut,*) cutm
    do lm=1,17
    cutlmp(lm)=rcut
    cutlmn(lm)=-rcut
    end do
    inquire(file='diop.cut',exist=cutf)
    if (cutf) then
    open(unit=10,file='diop.cut')
    print *
    print *, ' Information: Using Diop.cut'
    do
    read(10,*,end=10) lm,cutlmn(lm),cutlmp(lm)
    end do    
    end if
10  filein=adjustr(filein)
    open(unit=20,file=filein//'.oph',form='unformatted')
    read(20) max_op,max_lmo
    allocate(opn(max_op))
    allocate(meV(max_op))
    allocate(pl(max_op))
    allocate(pr(max_op))
    if (max_op/=0) then       
    do i=1,max_op
    read(20) opn(i),meV(i),pl(i),pr(i)
    no_Op=no_Op+1
    end do
    end if
    cr=1
    cl=1
    chl(1)(1:2)=opn(1)(1:2)
    chr(1)(1:2)=opn(1)(3:4)
    do i=1,no_Op
    chl(cl+1)(1:2)=opn(i)(1:2)
    do j=1,cl
        cl=cl+1
        if (ichar(chl(cl)(1:1))>max_shl) max_shl=ichar(chl(cl)(1:1))
        if (ichar(chl(cl)(2:2))>max_shl) max_shl=ichar(chl(cl)(2:2))
        if (ichar(chl(cl)(1:1))<min_shl) min_shl=ichar(chl(cl)(1:1))
        if (ichar(chl(cl)(2:2))<min_shl) min_shl=ichar(chl(cl)(2:2))
        if (chl(cl)==chl(j)) then
            cl=cl-1 
            go to 20
        end if    
        if (cl+1>maxshell2) stop ' cl dimension'
        cl=cl-1
    end do
    cl=cl+1
20  chr(cr+1)(1:2)=opn(i)(3:4)
    do j=1,cr
        cr=cr+1
        if (ichar(chr(cr)(1:1))>max_shr) max_shr=ichar(chr(cr)(1:1))
        if (ichar(chr(cr)(2:2))>max_shr) max_shr=ichar(chr(cr)(2:2))
        if (ichar(chr(cr)(1:1))<min_shr) min_shr=ichar(chr(cr)(1:1))
        if (ichar(chr(cr)(2:2))<min_shr) min_shr=ichar(chr(cr)(2:2))
        if (chr(cr)==chr(j)) then
            cr=cr-1
            go to 40
        end if
        if (cr+1>maxshell2) stop ' cr dimension'
        cr=cr-1
    end do
    cr=cr+1
40  end do
    if (cr/=cl) stop ' cr/=cl not symmetric '
    if (max_shr-min_shr/=max_shl-min_shl) stop ' shells not symmetric '
    open(unit=21,file=fileout//'.dop',action='write')
    open(unit=22,file=fileout//'.lsp',action='write')
    write(21,'(i4)') max_lmo
    write(22,'(i4)') max_lmo
    do i=1,cr
        do j=1,cr
        if ((ichar(chl(i)(1:1))+max_shr/2==ichar(chr(j)(1:1))) &
            .and.(ichar(chl(i)(2:2))+max_shr/2==ichar(chr(j)(2:2))) ) then 
                        chrc(i)(1:2)=chr(j)(1:2)
        end if
        if ((ichar(chl(i)(1:1))-max_shr/2==ichar(chr(j)(1:1))) &
            .and.(ichar(chl(i)(2:2))-max_shr/2==ichar(chr(j)(2:2))) ) then 
                        chrc(i)(1:2)=chr(j)(1:2)
        end if
        end do
    end do
    do i=1,cr
        chr(i)(1:2)=chrc(i)(1:2)
    end do
    do lm=0,max_lmo,2
    do pt=-1,1,2
    nonz=0
    matrix=real(0,8)
    mat=real(0,8)
    alpha=real(0,8)
    beta=real(0,8)
    mapl=0
    mapr=0
    nsp=0
    indx=0
    sp1=0
    sp2=0
    sp3=0
    sp4=0
    rtlm=sqrt(real(lm+1,8))
        do i=1,no_Op
            if (ichar(opn(i)(5:5))/=lm) go to 100
            if (pl(i)/=pr(i)) go to 100
            if (pl(i)/=pt) go to 100
            indxl=100*ichar(opn(i)(1:1))+ ichar(opn(i)(2:2))
            indxr=100*ichar(opn(i)(3:3))+ ichar(opn(i)(4:4))
            if (indxl>indxr) go to 100
            do j=1,cr
            if (100*ichar(chl(j)(1:1))+ichar(chl(j)(2:2))==indxl) indxcl=j
            if (100*ichar(chr(j)(1:1))+ichar(chr(j)(2:2))==indxr) indxcr=j
            end do
            mat(indxcl,indxcr)=mat(indxcl,indxcr)+meV(i)
            nonz=nonz+1
100     end do
        if (nonz==0) then 
            write(21,'(2i4)') lm,pt
            write(22,'(2i4)') lm,pt
            write(22,'(i4)') 0
            write(21,'(i4)') 0
            go to 300
        end if
        k=0
        l=0
        do i=1,cr
        write(60,'(10f7.3)') mat(i,1:cr)
        end do
        do i=1,cr
            if (sum(abs(mat(1:cr,i)))/=0.0) then
                k=k+1
                mapr(k)=i
            end if
            if (sum(abs(mat(i,1:cr)))/=0.0) then
                l=l+1
                mapl(l)=i
            end if
        end do
        if (k/=l) stop '  submatrix not square '
        if (k==0) then 
            write(21,'(2i4)') lm,pt
            write(22,'(2i4)') lm,pt
            write(21,'(i4)') 0
            write(22,'(i4)') 0
            go to 300
        end if
        do i=1,l
        do j=1,k
        matrix(i,j)=mat(mapl(i),mapr(j))
        end do           
        end do 
        do i=1,k
        write(65,'(10f7.3)') matrix(i,1:k)
        end do          
        call tred2(matrix,k,maxshell2,alpha,beta)
        tqlim=2000
        call tqli(alpha,beta,k,maxshell2,matrix,tqlim)
        call indexx(k,alpha,indx)      
        ksave=k    
        write (50,*) ' Lambda, parity ',lm/2,pt
        write (50,'(10f8.3)') alpha(indx(1:k))/rtlm
        do i=1,ksave
            if (matrix(i,i)/=1.0) then
                write(55,*) lm,pt,alpha(i)/rtlm
                write(55,'(10f8.3)') matrix(1:ksave,i)
            end if
        end do
        write(21,'(2i4)') lm,pt
        write(22,'(2i4)') lm,pt
        l=0
        do i=1,ksave
            if (alpha(indx(i))/rtlm>cutlmp(lm/2).or.&
                        alpha(indx(i))/rtlm<cutlmn(lm/2)) l=l+1
        end do
        write(21,'(i4)') l
        do i=1,ksave
            if ((alpha(indx(i))/rtlm>cutlmp(lm/2).and.&
                alpha(indx(i))/rtlm>0.0).or.&
                        (alpha(indx(i))/rtlm<cutlmn(lm/2).and.&
                        alpha(indx(i))/rtlm<0.0)) then
                write(21,'(f11.7)') alpha(indx(i))
                k=0
                do j=1,ksave
                    if (abs(matrix(j,indx(i)))>cutm) then
                    if (100*ichar(chl(mapl(j))(1:1))+ichar(chl(mapl(j))(2:2))<&
                            100*ichar(chr(mapr(j))(1:1))+ichar(chr(mapr(j))(2:2))) then
                        k=k+1
                    end if
                    end if
                end do
                write(21,'(i6)') k
                do j=1,ksave
                    if (abs(matrix(j,indx(i)))>cutm) then
                    if (100*ichar(chl(mapl(j))(1:1))+ichar(chl(mapl(j))(2:2))<&
                            100*ichar(chr(mapr(j))(1:1))+ichar(chr(mapr(j))(2:2))) then
                        write(21,'(4i4,f11.7)') ichar(chl(mapl(j))(1:1)),ichar(chl(mapl(j))(2:2))&
                            ,ichar(chr(mapr(j))(1:1)),ichar(chr(mapr(j))(2:2)),matrix(j,indx(i))
                    end if
                    end if
                end do
                do j=1,ksave
                    if (abs(matrix(indx(j),indx(i)))>cutm) then
                    if (100*ichar(chl(mapl(j))(1:1))+ichar(chl(mapl(j))(2:2))<&
                            100*ichar(chr(mapr(j))(1:1))+ichar(chr(mapr(j))(2:2))) then
                            if (nsp==0) then
                                nsp=1
                                sp1(1)=ichar(chl(mapl(j))(1:1)) 
                                sp2(1)=ichar(chl(mapl(j))(2:2)) 
                                sp3(1)=ichar(chr(mapr(j))(1:1)) 
                                sp4(1)=ichar(chr(mapr(j))(2:2)) 
                            else 
                                do mm=1,nsp
                                if ( sp1(mm)==ichar(chl(mapl(j))(1:1)).and. &
                                      sp2(mm)==ichar(chl(mapl(j))(2:2))&
                                      .and. sp3(mm)==ichar(chr(mapr(j))(1:1)).and. &
                                      sp4(mm)==ichar(chr(mapr(j))(2:2)) ) go to 250
                                end do 
                                nsp=nsp+1
                                sp1(nsp)=ichar(chl(mapl(j))(1:1)) 
                                sp2(nsp)=ichar(chl(mapl(j))(2:2)) 
                                sp3(nsp)=ichar(chr(mapr(j))(1:1)) 
                                sp4(nsp)=ichar(chr(mapr(j))(2:2)) 
250                             continue     
                             end if        
                    end if
                    end if
                end do
             end if           
        end do
    write(22,'(i4)') nsp
    do mm=1,nsp
    write(22,'(4i4)') sp1(mm),sp2(mm),sp3(mm),sp4(mm) 
    end do   
300 end do
    end do
    end program NuDiop
    
