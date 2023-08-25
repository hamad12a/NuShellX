! This file contains all the modules for the generic program gOrth for
! NuShell, NuShell@MSU and SunShell

module Cholesky

use BlockParameters

implicit none

contains

      SUBROUTINE choldc(a,n,np,p)
      INTEGER n,np
      REAL(kind=rl) a(np,np),p(np)
      INTEGER i,j,k,t
      REAL(kind=rl) sum
      do i=1,n
        do j=i,n
          sum=a(i,j)
          do k=i-1,1,-1
            sum=sum-a(i,k)*a(j,k)
          end do
          if(i.eq.j)then
            if (sum<=real(0,rl)) then
                p(i)=real(0,rl)
            else
                p(i)=sqrt(sum)
            end if
          else
            if (p(i)==real(0,rl)) then
                a(j,i)=real(0,rl)
            else
                a(j,i)=sum/p(i)
            end if
          endif
        end do
      end do
      return
      END SUBROUTINE choldc
!  (C) Copr. 1986-92 Numerical Recipes Software 7$-i%.
      SUBROUTINE inverse(a,n,np,p)
      INTEGER n,np
      REAL(kind=rl) a(np,np),p(np)
      INTEGER i,j,k
      REAL(kind=rl) sum
      do i=1,n
        if (p(i)==real(0,rl)) then
            a(i,i)=real(0,rl)
        else
            a(i,i)=real(1,rl)/p(i)
        end if
        do j=i+1,n
          sum=real(0,rl)
          do k=i,j-1
            sum=sum-a(j,k)*a(k,i)
          end do
          if (p(j)==real(0,rl)) then
                a(j,i)=real(0,rl)
          else
                a(j,i)=sum/p(j)
          end if
        end do
      end do
      return
      END SUBROUTINE inverse
!  (C) Copr. 1986-92 Numerical Recipes Software 7$-i%.

end module Cholesky

module  OrthSupport

use InputNuShell
use ProjParameters
use MatrixParameters
use BlockParameters
use Partition
use Mscheme

type sgmt
    real(kind=rm),dimension(:,:),allocatable:: mat_sgmt
    integer::np
    integer::npr
    logical::zero=.false.
end type
integer:: nosgmt,noJT,nomdim,nopart
type(sgmt),dimension(:),allocatable::matort
real(kind=rko),dimension(:,:),allocatable:: orthog
integer:: mat_dim,no_spart,max_nJTu,no_sgmt
integer,dimension(:),allocatable:: dim_part,part_offset
type(spartition),dimension(:),allocatable:: spart
real(kind=rm),dimension(:,:),allocatable:: mat_seg,matrix
real(kind=rl),dimension(ptdim):: diag
real(kind=rl)::diage
real(kind=rl):: renorm
integer:: npcheck,no_level
logical:: warn=.false.,warn2=.false.,rst=.false.
character(len=1):: xtra
real::xx,yy

contains

    subroutine read_orthog(base,k,tab)
        implicit none
        integer,intent(in):: base,k
        integer:: pindx,ngood,j,t,dim,ng,i,cP
        logical,intent(in):: tab
        if (tab) then
            mat_dim=0
            max_nJTu=0
        end if
        
        read(base+k,err=10,end=11) pindx,ngood,j,t,dim,xx    
        if (ngood>noJT) then
            print *,' Increase max_nJT in ProjParameters or'
            print *,' SET NJT=I where I >>',ngood
            stop
        end if  
        npcheck=pindx
        if (output_control > 4) write(33,*)'Pindx', pindx,ngood,j,t,dim,xx        
        read(base+k,err=10,end=11) spart(pindx)
        if (output_control > 4) write(33,*) spart(pindx )%shell(1:no_shells)        
        part_offset(pindx)=mat_dim
        mat_dim=mat_dim+ngood
        if (mat_dim>nomdim) then
            print *,' Increase max_total_JT in MatrixParameters or'
            print *,' SET DIMENSION=I where I is matrix dimension which is >>',mat_dim
            stop
        end if
        dim_part(pindx)=ngood
        if (output_control > 4) write(33,*)' Offset,dim', part_offset(pindx),dim_part(pindx)   
        
        if (ngood > max_nJTu  ) max_nJTu=ngood
        if (ngood > 0 ) then
            do ng=1,ngood
                orthog(ng,1:ngood)=real(0,rkw)
                read(base+k,err=10,end=11) orthog(ng,ng:ngood)
                if (output_control > 4) write(33,*) orthog(1,1:ngood)
            end do
        end if
        if (ngood > 0)call invert(ngood)
        if (ngood > 0 ) then
            allocate(matort(pindx)%mat_sgmt(ngood,ngood))
            matort(pindx)%np=pindx
            matort(pindx)%npr=pindx
            matort(pindx)%mat_sgmt(1:ngood,1:ngood)=orthog(1:ngood,1:ngood)
            matort(pindx)%zero=.false.
        else
            matort(pindx)%np=pindx
            matort(pindx)%npr=pindx
            matort(pindx)%zero=.true.
        end if 
        return
10      print *, ' Error reading ort file'
        stop
11      print *, ' EOF reading ort file'
        stop
               
        
    end subroutine read_orthog
    
    subroutine invert(ngood)
        implicit none
        integer,intent(in):: ngood
        integer:: i,j,k
        real(kind=rko):: sum
        
        do i=1,ngood
            orthog(i,i)=1/orthog(i,i)
            do j=i+1,ngood
                sum=real(0,rk)
                do k=i,j-1
                    sum=sum-orthog(k,j)*orthog(i,k)
                end do
                orthog(i,j)=sum/orthog(j,j)
            end do
        end do 
    end subroutine invert
    
    subroutine orthog_matrix(base,k,basei)
        implicit none
        integer,intent(in)::base,k,basei
        integer:: np,i,j,l,m,npr,ii,n,dimpt,dimptr,num_seg
        logical::reset,done,zero
        real(kind=rk) diag0
        reset=.true.
        num_seg=0
        allocate(mat_seg(max_nJTu,max_nJTu))
        allocate(matrix(max_nJTu,max_nJTu))
        do
            call read_matrix_sgmt(basei,k,reset,np,npr,dimpt,dimptr,done,zero)
            reset=.false.
            if (done) exit
            if (zero) go to 100
                     mat_seg(1:dim_part(npr),1:dim_part(np))=real(0,rl)
!$OMP PARALLEL PRIVATE(m,j,l)
!$OMP DO 
                     do m=1,dim_part(np)
                         do j=1,dim_part(npr)
                             do l=1,j
                                 mat_seg(j,m)=mat_seg(j,m)+matrix(l,m) &
                                 *matort(npr)%mat_sgmt(l,j)
                             end do
                         end do
                     end do
!$OMP END PARALLEL
                    if (np==npr) then
                        diag0=0.0
                        do m=1,no_shells
                            diag0=diag0+spe(m)*spart(npr)%shell(m)
                        end do
                    end if
                    matrix(1:dim_part(npr),1:dim_part(np))= &
                                mat_seg(1:dim_part(npr),1:dim_part(np))*renorm
                    if (np==npr) then
                        do m=1,dim_part(np)
                            matrix(m,m)=matrix(m,m)+diag0
                        end do
                    end if
                    num_seg=num_seg+1
                    write(base+k) dimpt,dimptr,np,npr,matrix(1:dimptr,1:dimpt),&
                             part_offset(np),part_offset(npr),spart(np)%shell(1:no_shells)
100         end do ! over sgmt
            write(800+k) num_seg,mat_dim,xx
            deallocate(mat_seg)
            deallocate(matrix)
            do np=1,no_spart
                if (dim_part(np) > 0) deallocate(matort(np)%mat_sgmt)
            end do
            
    end subroutine orthog_matrix
    
    subroutine read_matrix_sgmt(base,k,reset,np,npr,dimpt,dimptr,done,zero)
        implicit none
        integer,intent(in):: base,k
        integer,intent(out):: dimpt,dimptr
        integer,intent(out)::npr,np
        logical,intent(in)::reset
        logical,intent(out)::done,zero
        if(reset) then
            no_sgmt=0
        end if
            zero=.false.
            done=.false.
            read(base+k,err=200,end=201)dimpt,dimptr,np,npr,zero,yy
            if (dimpt <  -1) then 
                done=.true.
                return
            end if
            if (yy/=xx) stop ' Error: .ort and .mtx files inconsistent'
            if (output_control >8) write(62,*)dimpt,dimptr,np,npr,zero
            if (dimpt < 0 .or. zero) then
                zero=.true.
                return
            end if
                if (dimpt /= dim_part(np) .or.  dimptr  /= dim_part(npr)) then
                    print *, ' Dimension error',np,npr,dimpt,dim_part(np),dimptr,dim_part(npr)
                    stop
                end if
                no_sgmt=no_sgmt+1
                read(base+k,err=200,end=201) matrix(1:dimptr,1:dimpt)
                return
200     print *, ' Error reading .mtx file '
        stop        
201     print *, ' EOF reading .mtx file '
        stop        
end subroutine read_matrix_sgmt

end module  OrthSupport
    
!  NuOrth.f90 
!
!  FUNCTIONS:
!  NuOrth      - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: NuOrth
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program NuOrth
    
    use Parameters
    use BlockParameters
    use MatrixParameters
    use ProjParameters
    use InputNuShell
    use OutputNuShell
    use Shells
    use Files
    use OrthSupport
    use Extra
    
    implicit none

    ! Variables
    
    real:: t1,t2,et1,et2,ms
    integer:: hr,mn,sc

    integer:: j,t,k,npt,dt,ns,res
    
    character(len=6):: inExt
    character(len=10):: btime
    character(len=8):: bdate
    character(len=10)::chnomdim,chnopart,chnojt,chnosgmt,chpart='PARTITIONS',chnoblock
    character(len=9)::chmdim='DIMENSION'
    character(len=8)::chsgmt='SEGMENTS'
    character(len=3)::chnjt='NJT'
    character(len=5)::chblock='BLOCK'
    
    res=get_env_v(chnjt,chnojt)
    if (res /= 0 ) then
        read(chnojt,*) noJT
    else
        noJT=max_nJT
    end if
    allocate(orthog(noJT,noJT))
    
    res=get_env_v(chmdim,chnomdim)
    if (res /= 0 ) then
        read(chnomdim,*) nomdim
    else
        nomdim=max_total_JT
    end if
    
    res=get_env_v(chsgmt,chnosgmt)
    if (res /= 0 ) then
        read(chnosgmt,*) nosgmt
    else
        nosgmt=max_sgmt
    end if
    
    res=get_env_v(chpart,chnopart)
    if (res /= 0) then
        read(chnopart,*) nopart
    else
        nopart=max_no_partitions
    end if
    allocate(part_offset(nopart))
    allocate(dim_part(nopart))
    allocate(spart(nopart))
    allocate(matort(nopart))
    
    ! Body of NuOrth
    
    call output_header(6)
    call input_nucleus(inExt,xtra)
    
    call output_welcome(6)    
    call input(inExt)
    call output_shells(6)    
    call check_shell_order()
    call setup_shells()
    if (use_isospin) then
        dt=2
    else
        dt=0
    end if
    call setup_file_ext()
    call output_end_setup(6)        
    call cpu_time(t1)
    call date_and_time(bdate,btime)
    call open_files_read(no_spart,'.ort',100,dt)
    if (no_spart>nopart) then
        print *,' Increase max_no_partitions in Parameters or'
        print *,' SET PARTITIONS=I where I >>',no_spart
        stop
    end if
    call open_files_read(no_spart,'.mtx',200,dt,-1,xtra)
    call open_files(no_spart,'.mat',400,dt,-1,xtra)
    call open_files(no_spart,'.sgt',800,dt,-1,xtra)
    if (output_control > 4) open(unit=33,file=Nucleus//'Orthog.txt', action='write')
    if (output_control > 8) open(unit=62,file=Nucleus//xtra//'matrd.txt',action='write')
    k=0
    spe=0.0
    open(unit=10, file=Nucleus//xtra//'.spe')
    read(10,*) (spe(k),k=1,no_shells)
    read(10,*) renorm
    read(10,*) no_level
    close(10)
    k=0
    
    do j=1,no_spart
        spart(j)=0
    end do

    do j=min_cJ2,max_cJ2,2
        do t=min_cT2,max_cT2,2
            if (j <= max_cJ2  .and.  t <=  max_cT2 ) then
                if (output_control > 4) write(33,*) '2*J, 2*T',j,t
                if (output_control > 8) write(62,*) '2*J, 2*T',j,t
                call read_orthog(100,k,.true.)
                if (no_spart > 1) then
                    do npt=2,no_spart
                        call read_orthog(100,k,.false.)
                    end do
                end  if
                if (mat_dim == 0) then
                    write(800+k) 0,0,xx
                    go to 50
                end if
                if (output_control > 8)write(62,*) '  2*J, 2*T',j,t
                if (output_control > 4)write(33,*) '  2*J, 2*T',j,t
                call orthog_matrix(400,k,200)
            end if
50          k=k+1
        end do
    end do
    deallocate(orthog)
    deallocate(part_offset)
    deallocate(dim_part)
    deallocate(spart)

    read(btime,'(i2,i2,i2,f4.3)') hr,mn,sc,ms
    et1=hr*3600.+mn*60.+sc*1. +ms
    call date_and_time(bdate,btime)
    call cpu_time(t2)
    read(btime,'(i2,i2,i2,f4.3)') hr,mn,sc,ms
    et2=hr*3600.+mn*60.+sc*1. +ms
    et1=et2-et1
    if (et1<0.0) et1=et1+24.*3600.
    t1=t2-t1
    call output_time('NuOrth',t1,et1)    
    print *
    
    end program NuOrth
    