! These modules are for NuShell,NushellX,SunShell and NuShell@MSU
! DO NOT COMPILE ALONE

            
module Files

use InputNuShell
use Extra

character(len=1),dimension(0:max_code):: an_code
data an_code /'0','1','2','3','4','5','6','7','8','9',  &
              '_','$','c','d','e','f','g','h','i','j',  &
              'k','l','m','n','o','p','q','r','s','t',  &
              'u','v','w','x','y','z','@','A','B','C',  &
              'D','E','F','G','H','I','J','K','L','M',  &
              'N','M','P','Q','R','S','T','U','V','W',  &
              'X','Y','Z','£','a','%','^','&','=','b',  &
              '-','+','#'/
              
contains

    function file_code(i,j,k)
    implicit none
    integer,intent(in):: i,j,k
    character(len=3):: file_code
    integer:: ii,jj,kk
    ii=i
    jj=j
    kk=k
    if (i<0) ii=i+max_code
    if (j<0) jj=j+max_code
    if (k<0) kk=k+max_code
    file_code=an_code(ii)//an_code(jj)//an_code(kk)
    
    end function file_code

    subroutine setup_file_ext()
    implicit none
    integer:: i,j,k
        
    end subroutine setup_file_ext
    
    subroutine input_nucleus(inExt,xtra)
        implicit none
        character(len=*):: inExt
        character(len=*),optional:: xtra
        integer:: status
        logical::  present
        status=-1
        call get_cmd_arg(1,InExt,status)
        inExt=adjustr(inExt)
        if (status>=1 .and. present(xtra)) call get_cmd_arg(2,xtra,status)
        if (status>=1) return   
        write(6,'(a10)',advance='no') '  Nucleus '
        read(5,'(a6)') inExt
        inExt=adjustr(inExt)
        if (present(xtra)) then
            write(6,'(a19)',advance='no') '  Interaction code '
            read(5,'(a1)') xtra
        end if
   
    end subroutine input_nucleus
    
    subroutine input_trans_one(inExt)
        implicit none
        character(len=*):: inExt
        integer:: status
        status=-1
        call get_cmd_arg(3,InExt,status)
        inExt=adjustr(inExt)
        if (status>=1) return   
        write(6,'(a9)',advance='no') '  File 1 '
        read(5,'(a13)') inExt
        inExt=adjustr(inExt)
   
    end subroutine input_trans_one
    
    subroutine input_trans_two(inExt)
        implicit none
        character(len=*):: inExt
        integer:: status
        status=-1
        call get_cmd_arg(4,InExt,status)
        inExt=adjustr(inExt)
        if (status>=1) then
            return
        else
            call get_cmd_arg(3,inext,status)
            if (status>=1) then
                inext='             '
                return
            end if    
        end if   
        write(6,'(a9)',advance='no') '  File 2 '
        read(5,'(a13)') inExt
        inExt=adjustr(inExt)
   
    end subroutine input_trans_two
    
    subroutine input_op(inExt,xtra)
        implicit none
        character(len=*):: inExt
        character(len=*),optional:: xtra
        integer:: status
        logical::  present
        status=-1
        call get_cmd_arg(3,InExt,status)
        inExt=adjustr(inExt)
        if (status>=1 .and. present(xtra)) call get_cmd_arg(4,xtra,status)
        if (status>=1) return   
        write(6,'(a11)',advance='no') '  Operator '
        read(5,'(a6)') inExt
        inExt=adjustr(inExt)
        if (present(xtra)) then
            write(6,'(a16)',advance='no') '  Filename code '
            read(5,'(a1)') xtra
        end if
   
    end subroutine input_op


    subroutine input_type(inExt)
        implicit none
        character(len=*):: inExt
        integer:: status
        status=-1
        call get_cmd_arg(4,InExt,status)
        inExt=adjustr(inExt)
        if (status>=1) return   
        write(6,'(a10)',advance='no') '  np or i '
        read(5,'(a2)') inExt
        inExt=adjustr(inExt)
   
    end subroutine input_type


    subroutine input_nucleus_final(inExt,xtra)
        implicit none
        character(len=*),intent(out):: InExt
        character(len=*),intent(out),optional:: xtra
        integer:: status
        logical::  present
        status=-1
        if (.not.present(xtra)) call get_cmd_arg(2,InExt,status)
        if (present(xtra)) call get_cmd_arg(3,InExt,status)
        inExt=adjustr(inExt)
        if (status>=1 .and. present(xtra)) call get_cmd_arg(4,xtra,status)
        if (status>=1) return   
        write(6,'(a10)',advance='no') '  Nucleus '
        read(5,'(a6)') inExt
        inExt=adjustr(inExt)
        if (present(xtra)) then
            write(6,'(a19)',advance='no') '  Interaction code '
            read(5,'(a1)') xtra
        end if
        
    end subroutine input_nucleus_final

    subroutine input_nucleus_cluster(inExt,xtra)
        implicit none
        character(len=*),intent(out):: InExt
        character(len=*),intent(out),optional:: xtra
        integer:: status
        logical::  present
        status=-1
        if (.not.present(xtra)) call get_cmd_arg(3,InExt,status)
        if (present(xtra)) call get_cmd_arg(5,InExt,status)
        inExt=adjustr(inExt)
        if (status>=1 .and. present(xtra)) call get_cmd_arg(6,xtra,status)
        if (status>=1) return   
        write(6,'(a10)',advance='no') '  Cluster '
        read(5,'(a6)') inExt
        inExt=adjustr(inExt)
        if (present(xtra)) then
            write(6,'(a19)',advance='no') '  Interaction code '
            read(5,'(a1)') xtra
        end if
        
    end subroutine input_nucleus_cluster

    subroutine input_model_space(inExt,xtra)
        implicit none
        character(len=*),intent(out):: InExt
        character(len=*),intent(out),optional:: xtra
        logical::  present
        integer:: status
        status=-1
        call get_cmd_arg(2,inExt,status)
        inExt=adjustr(inExt)
        if (status>=1 .and. present(xtra)) call get_cmd_arg(3,xtra,status)
        if (status>=1) return   
        write(6,'(a14)',advance='no') '  Model Space '
        read(5,'(a6)') inExt
        inExt=adjustr(inExt)
        if (present(xtra)) then
            write(6,'(a19)',advance='no') '  Interaction code '
            read(5,'(a1)') xtra
        end if
        
    end subroutine input_model_space

    subroutine input_intfile(intfile,no)
        implicit none
        character(len=*),intent(out):: intfile
        integer,intent(in),optional:: no
        integer::status
        logical::  present
            status=-1
            if (.not.present(no)) call get_cmd_arg(3,intfile,status)
            intfile=adjustr(intfile)
            if (status>=1) return
            if (.not.present(no))write(6,'(a19)',advance='no') '  Interaction file '
            if (present(no).and. no==1) call get_cmd_arg(1,intfile,status)
            intfile=adjustr(intfile)
            if (status>=1) return
            if (present(no).and. no==1)write(6,'(a12)',advance='no') '  MBME file '
            if (present(no).and. no==2) call get_cmd_arg(3,intfile,status)
            intfile=adjustr(intfile)
            if (status>=1) return
            if (present(no).and. no==2)write(6,'(a15)',advance='no') '  MBME op file '
            intfile=adjustr(intfile)
            if (present(no).and. no==3) call get_cmd_arg(5,intfile,status)
            intfile=adjustr(intfile)
            if (status>=1) return
            if (present(no).and. no==3)write(6,'(a19)',advance='no') '  MBME matrix file '
            if (present(no).and. no==4) call get_cmd_arg(2,intfile,status)
            intfile=adjustr(intfile)
            if (status>=1) return
            if (present(no).and. no==4)write(6,'(a19)',advance='no') '  TRDA  matrix file '
            if (present(no).and. no==5) call get_cmd_arg(3,intfile,status)
            intfile=adjustr(intfile)
            if (status>=1) return
            if (present(no).and. no==5)write(6,'(a19)',advance='no') '  TRDB matrix file '
            if (present(no).and. no==6) call get_cmd_arg(4,intfile,status)
            intfile=adjustr(intfile)
            if (status>=1) return
            if (present(no).and. no==6)write(6,'(a19)',advance='no') '  OPH  matrix file '
            if (present(no).and. no==7) call get_cmd_arg(4,intfile,status)
            intfile=adjustr(intfile)
            if (status>=1) return
            if (present(no).and. no==7)write(6,'(a19)',advance='no') '  MBME   op   file '
            if (present(no).and. no==8) call get_cmd_arg(5,intfile,status)
            intfile=adjustr(intfile)
            if (status>=1) return
            if (present(no).and. no==8)write(6,'(a19)',advance='no') '  Initial  Nucleus   '
            if (present(no).and. no==8) then
                read(5,'(a6)') intfile
                intfile=adjustr(intfile)
            else
                read(5,'(a7)') intfile
                intfile=adjustr(intfile)
            end if
    end subroutine input_intfile

    subroutine open_files(np,cExt,base,dt,kin,xtra)
    implicit none
    character(len=*):: cExt
    integer,intent(in):: np,base,dt
    integer,intent(in),optional::kin
    character(len=1),intent(in),optional::xtra
    integer i,j,k
    logical::  present
    
    integer xj,xt
    xj=0
    xt=0
    if (cExt == '.nba' .or. cExt == '.bas') then
        xj=2
        xt=dt
    end if
    
        k=0

        do i=min_cJ2,max_cJ2+xj,2
            do j=min_cT2,max_cT2+xt,2 
                 if (.not.present(xtra)) then
                 if ((present(kin) .and. k==kin) .or. .not.present(kin) .or. &
                                                    (present(kin) .and. kin <0)) then
                 open(unit=base+k, file=Nucleus//file_code(i/2,j/2,nmp/2+max_code/2)//cExt,&
                         action='write',form='unformatted')
                 write(base+k) np 
                 end if
                 end if
                 if (present(xtra)) then
                 open(unit=base+k, file=Nucleus//file_code(i/2,j/2,nmp/2+max_code/2)//xtra//cExt, &
                        action='write',form='unformatted',err=20)
                 write(base+k) np 
                 end if
                 k=k+1  
            end do
        end do
        return
10      print *, ' Error opening ',Nucleus//file_code(i/2,j/2,nmp/2+max_code/2)//cExt
        stop        
20      print *, ' Error opening ',Nucleus//file_code(i/2,j/2,nmp/2+max_code/2)//xtra//cExt
        stop        
    end subroutine open_files
 
    subroutine open_files7(name,np,cExt,base,dt,kin,xtra)
    implicit none
    character(len=*):: cExt
    character(len=*):: name
    integer,intent(in):: np,base,dt
    integer,intent(in),optional::kin
    character(len=1),intent(in),optional::xtra
    integer i,j,k
    logical::  present
    
        k=0
    
        do i=min_cJ2,max_cJ2,2
            do j=min_cT2,max_cT2,2 
                 if (.not.present(xtra)) then
                 if ((present(kin) .and. k==kin) .or. .not.present(kin) .or. &
                                                (present(kin) .and. kin <0)) then
                 open(unit=base+k, file=name//file_code(i/2,j/2,nmp/2+max_code/2)//cExt, &
                        action='write',form='unformatted', err=10)
                 write(base+k) np 
                 end if
                 end if
                 if (present(xtra)) then
                 open(unit=base+k, file=name//file_code(i/2,j/2,nmp/2+max_code/2)//xtra//cExt, &
                        action='write',form='unformatted',err=20)
                 write(base+k) np 
                 end if
                 k=k+1  
            end do
        end do
        return
10      print *, ' Error opening ',name//file_code(i/2,j/2,nmp/2+max_code/2)//cExt
        stop        
20      print *, ' Error opening ',name//file_code(i/2,j/2,nmp/2+max_code/2)//xtra//cExt
        stop        

    end subroutine open_files7
    
    subroutine open_files7_mbme(name,np,cExt,base,dt,mindmj,maxdmj,mindmt,maxdmt, &
                            min_CJ1,max_cJ1,min_cT1,max_cT1,nmp1,nmp2,kin,xtra)
    implicit none
    character(len=*):: cExt
    character(len=2):: dmdmt
    character(len=*):: name
    integer,intent(in):: np,base,dt,mindmj,mindmt,maxdmj,maxdmt,min_CJ1,&
                                max_cJ1,min_cT1,max_cT1,nmp1,nmp2
    integer,intent(in),optional::kin
    character(len=1),intent(in),optional::xtra
    integer i,j,k,delmj,delmt
    logical:: skip,present
    
        k=0
    
        do i=min_cJ2,max_cJ2,2
            do j=min_cT2,max_cT2,2 
                 do delmj=mindmj,maxdmj,2
                 skip=.false.
                 if (i+delmj > max_cJ1 .or. i+delmj < min_cJ1) skip=.true.
                 do delmt=mindmt,maxdmt,2
                 if (use_isospin) then
                 if (j+delmt > max_cT1 .or. j+delmt < min_cT1 .or. skip) go to 5
                 else
                 if (nmp2+delmt/=nmp1 .or. skip) go to 5
                 end if
                 dmdmt=an_code(delmj+18)//an_code(delmt+18)
                 if (.not.present(xtra)) then
                 if ((present(kin) .and. k==kin) .or. .not.present(kin) .or. &
                                                (present(kin) .and. kin <0)) then
                 open(unit=base+k,  &
                    file=name//file_code(i/2,j/2,nmp/2+max_code/2)//dmdmt//cExt,&
                    action='write',form='unformatted',err=10)
                 write(base+k) np 
                 end if
                 end if
                 if (present(xtra)) then
                 open(unit=base+k, &
                      file=name//file_code(i/2,j/2,nmp/2+max_code/2)//dmdmt//xtra//cExt, &
                        action='write',form='unformatted',err=20)
                 write(base+k) np 
                 end if
5                k=k+1 
                 end do 
                 end do 
            end do
        end do
        return
10      print * ,' Error opening ',name//file_code(i/2,j/2,nmp/2+max_code/2)//dmdmt//cExt
        stop        
20      print * ,' Error opening ',name//file_code(i/2,j/2,nmp/2+max_code/2)//dmdmt//xtra//cExt
        stop        

    end subroutine open_files7_mbme
    
    subroutine open_files7_read_mbme(name,np,cExt,base,dt,mindmj,maxdmj,mindmt,maxdmt, &
                                min_CJ1,max_cJ1,min_cT1,max_cT1,nmp1,nmp2,kin,xtra)
    implicit none
    character(len=*):: cExt
    character(len=2):: dmdmt
    character(len=*):: name
    integer,intent(in):: base,dt,mindmj,mindmt,maxdmj,maxdmt,min_CJ1,max_cJ1,&
                                min_cT1,max_cT1,nmp1,nmp2
    integer,intent(out):: np
    integer,intent(in),optional::kin
    character(len=1),intent(in),optional::xtra
    integer i,j,k,delmj,delmt
    logical:: skip,present
        k=0
    
        do i=min_cJ2,max_cJ2,2
            do j=min_cT2,max_cT2,2 
                 do delmj=mindmj,maxdmj,2
                 skip=.false.                
                 if (i+delmj > max_cJ1 .or. i+delmj < min_cJ1) skip=.true.
                 do delmt=mindmt,maxdmt,2
                 if (use_isospin) then
                 if (j+delmt > max_cT1 .or. j+delmt < min_cT1 .or. skip) go to 5
                 else
                 if (nmp2+delmt /= nmp1 .or. skip) go to 5
                 end if
                 dmdmt=an_code(delmj+18)//an_code(delmt+18)
                 if (.not.present(xtra)) then
                 if ((present(kin) .and. k==kin) .or. .not.present(kin) .or. &
                                                        (present(kin) .and. kin <0)) then
                 open(unit=base+k,  &
                    file=name//file_code(i/2,j/2,nmp/2+max_code/2)//dmdmt//cExt,&
                    action='read',form='unformatted',err=10)
                 read(base+k) np 
                 end if
                 end if
                 if (present(xtra)) then
                 open(unit=base+k, &
                      file=name//file_code(i/2,j/2,nmp/2+max_code/2)//dmdmt//xtra//cExt, &
                        action='read',form='unformatted',err=20)
                 read(base+k) np 
                 end if
5                k=k+1 
                 end do 
                 end do 
            end do
        end do
        return
10      print *, ' Error opening ',name//file_code(i/2,j/2,nmp/2+max_code/2)//dmdmt//cExt
        stop        
20      print *, ' Error opening ',name//file_code(i/2,j/2,nmp/2+max_code/2)//dmdmt//xtra//cExt
        stop        

    end subroutine open_files7_read_mbme
    
    subroutine open_files_form(np,cExt,base,dt,kin,xtra)
    implicit none
    character(len=*):: cExt
    integer,intent(in):: np,base,dt
    integer,intent(in),optional::kin
    character(len=1),intent(in),optional::xtra
    integer i,j,k
    logical::  present
    
        k=0
    
        do i=min_cJ2,max_cJ2,2
            do j=min_cT2,max_cT2,2 
                 if (.not.present(xtra)) then
                 if ((present(kin) .and. k==kin) .or. .not.present(kin) .or. &
                                                    (present(kin) .and. kin <0)) then
                 open(unit=base+k, file=Nucleus//file_code(i/2,j/2,nmp/2+max_code/2)//cExt,&
                         action='write',form='formatted',err=10)
                 end if
                 end if
                 if (present(xtra)) then
                 open(unit=base+k, file=Nucleus//file_code(i/2,j/2,nmp/2+max_code/2)//xtra//cExt,&
                         action='write',form='formatted',err=20)
                 end if
                 k=k+1  
            end do
        end do
        return
10      print *, ' Error opening ',Nucleus//file_code(i/2,j/2,nmp/2+max_code/2)//cExt
        stop        
20      print *, ' Error opening ',Nucleus//file_code(i/2,j/2,nmp/2+max_code/2)//xtra//cExt
        stop        
    end subroutine open_files_form
    
    
    subroutine open_files_form7(name,np,cExt,base,dt,kin,xtra)
    implicit none
    character(len=*):: cExt
    character(len=*):: name
    integer,intent(in):: np,base,dt
    integer,intent(in),optional::kin
    character(len=1),intent(in),optional::xtra
    integer i,j,k
    logical::  present
    
        k=0
    
        do i=min_cJ2,max_cJ2,2
            do j=min_cT2,max_cT2,2 
                 if (.not.present(xtra)) then
                 if ((present(kin) .and. k==kin) .or. .not.present(kin) .or. &
                                                    (present(kin) .and. kin <0)) then
                 open(unit=base+k, file=name//file_code(i/2,j/2,nmp/2+max_code/2)//cExt,&
                         action='write',form='formatted',err=10)
                 end if
                 end if
                 if (present(xtra)) then
                 open(unit=base+k, file=name//file_code(i/2,j/2,nmp/2+max_code/2)//xtra//cExt,&
                         action='write',form='formatted',err=20)
                 end if
                 k=k+1  
            end do
        end do
        return
10      print *, ' Error opening ',Nucleus//file_code(i/2,j/2,nmp/2+max_code/2)//cExt
        stop        
20      print *, ' Error opening ',Nucleus//file_code(i/2,j/2,nmp/2+max_code/2)//xtra//cExt
        stop        

    end subroutine open_files_form7
    
    subroutine open_files_form_read(np,cExt,base,dt,kin,xtra)
    implicit none
    character(len=*):: cExt
    integer,intent(in):: np,base,dt
    integer,intent(in),optional::kin
    character(len=1),intent(in),optional::xtra
    integer i,j,k
    logical::  present
    
        k=0
    
        do i=min_cJ2,max_cJ2,2
            do j=min_cT2,max_cT2,2 
                 if (.not.present(xtra)) then
                 if ((present(kin) .and. k==kin) .or. .not.present(kin) .or. &
                                                    (present(kin) .and. kin <0)) then
                 open(unit=base+k, file=Nucleus//file_code(i/2,j/2,nmp/2+max_code/2)//cExt, &
                            action='read',form='formatted',err=10)
                 end if
                 end if
                 if (present(xtra)) then
                 open(unit=base+k, file=Nucleus//file_code(i/2,j/2,nmp/2+max_code/2)//xtra//cExt,&
                             action='read',form='formatted',err=20)
                 end if
                 k=k+1  
            end do
        end do
        return
10      print *, ' Error opening ',Nucleus//file_code(i/2,j/2,nmp/2+max_code/2)//cExt
        stop        
20      print *, ' Error opening ',Nucleus//file_code(i/2,j/2,nmp/2+max_code/2)//xtra//cExt
        stop        

    end subroutine open_files_form_read
    
    subroutine open_files_read(np,cExt,base,dt,kin,xtra)
    implicit none
    character(len=*),intent(in):: cExt
    integer,intent(out):: np
    integer,intent(in):: base,dt
    integer,intent(in),optional::kin
    character(len=1),intent(in),optional::xtra
    integer i,j,k
    logical::  present
    
    integer xj,xt
    xj=0
    xt=0
    if (cExt == '.nba' .or. cExt == '.bas') then
        xj=2
        xt=dt
    end if
        
        k=0
    
        do i=min_cJ2,max_cJ2+xj,2
            do j=min_cT2,max_cT2+xt,2 
                 if (.not.present(xtra)) then
                 if ((present(kin) .and. k==kin) .or. .not.present(kin) .or. &
                                                    (present(kin) .and. kin <0)) then
                 open(unit=base+k, file=Nucleus//file_code(i/2,j/2,nmp/2+max_code/2)//cExt,&
                             action='read',form='unformatted',err=10)
                 read(base+k) np
                 end if
                 end if
                 if (present(xtra)) then
                 open(unit=base+k, file=Nucleus//file_code(i/2,j/2,nmp/2+max_code/2)//xtra//cExt,&
                            action='read',form='unformatted',err=20)
                 read(base+k) np
                 end if
                 k=k+1  
            end do
        end do
        return
10      print *, ' Error opening ',Nucleus//file_code(i/2,j/2,nmp/2+max_code/2)//cExt
        stop        
20      print *, ' Error opening ',Nucleus//file_code(i/2,j/2,nmp/2+max_code/2)//xtra//cExt
        stop        

    end subroutine open_files_read
    
    function inquire_files_exist(cExt,base,dt,kin,xtra)
    implicit none
    character(len=*),intent(in):: cExt
    integer,intent(in):: base,dt
    integer,intent(in),optional::kin
    character(len=1),intent(in),optional::xtra
    logical:: inquire_files_exist, exists
    integer i,j,k
    integer xj,xt
    logical::  present

    inquire_files_exist=.true.
    xj=0
    xt=0
    if (cExt == '.nba' .or. cExt == '.bas') then
        xj=2
        xt=dt
    end if
        
        k=0
    
        do i=min_cJ2,max_cJ2+xj,2
            do j=min_cT2,max_cT2+xt,2 
                 if (.not.present(xtra)) then
                 if ((present(kin) .and. k==kin) .or. .not.present(kin) .or.&
                                                         (present(kin) .and. kin <0)) then
                 inquire(file=Nucleus//file_code(i/2,j/2,nmp/2+max_code/2)//cExt, exist=exists)
                 inquire_files_exist=inquire_files_exist.and.exists
                 end if
                 end if
                 if (present(xtra)) then
                 inquire(file=Nucleus//file_code(i/2,j/2,nmp/2+max_code/2)//xtra//cExt, exist=exists)
                 inquire_files_exist=inquire_files_exist.and.exists
                 end if
                 k=k+1  
            end do
        end do

    end function inquire_files_exist
    
    
    subroutine open_files_read7(name,np,cExt,base,dt,kin,xtra)
    implicit none
    character(len=*),intent(in):: cExt
    character(len=*):: name
    integer,intent(out):: np
    integer,intent(in):: base,dt
    integer,intent(in),optional::kin
    character(len=1),intent(in),optional::xtra
    integer i,j,k
    logical::  present
    
    integer xj,xt
    xj=0
    xt=0
    if (cExt == '.nba' .or. cExt == '.bas') then
        xj=2
        xt=dt
    end if
        
        k=0
    
        do i=min_cJ2,max_cJ2+xj,2
            do j=min_cT2,max_cT2+xt,2 
                 if (.not.present(xtra)) then
                 if ((present(kin) .and. k==kin) .or. .not.present(kin) .or. &
                                                    (present(kin) .and. kin <0)) then
                 open(unit=base+k, file=name//file_code(i/2,j/2,nmp/2+max_code/2)//cExt,&
                             action='read',form='unformatted',err=10)
                 read(base+k) np
                 end if
                 end if
                 if (present(xtra)) then
                 open(unit=base+k, file=name//file_code(i/2,j/2,nmp/2+max_code/2)//xtra//cExt,&
                             action='read',form='unformatted',err=20)
                 read(base+k) np
                 end if
                 k=k+1  
            end do
        end do
        return
10      print *, ' Error opening ',name//file_code(i/2,j/2,nmp/2+max_code/2)//cExt
        stop        
20      print *, ' Error opening ',name//file_code(i/2,j/2,nmp/2+max_code/2)//xtra//cExt
        stop        

    end subroutine open_files_read7
    
    subroutine close_files(base,dt,kin)
    implicit none
    integer,intent(in):: base,dt
    integer,intent(in),optional::kin
    integer:: i,j,k,minus1=-1
    logical::  present
    
        k=0
    
        do i=min_cJ2,max_cJ2,2
            do j=min_cT2,max_cT2,2 
                 if ((present(kin) .and. k==kin) .or. .not.present(kin)) then
                 write(base+k) minus1,minus1,minus1,minus1,minus1,minus1,minus1,minus1
                 close(unit=base+k) 
                 end if
                 k=k+1  
            end do
        end do

    end subroutine close_files
    
    subroutine close_files_read(base,dt,kin)
    implicit none
    integer,intent(in):: base,dt
    integer,intent(in),optional::kin
    integer:: i,j,k
    logical::  present
    
        k=0
    
        do i=min_cJ2,max_cJ2,2
            do j=min_cT2,max_cT2,2 
                 if ((present(kin) .and. k==kin) .or. .not.present(kin)) then
                 close(unit=base+k) 
                 end if
                 k=k+1  
            end do
        end do

    end subroutine close_files_read
    
end module Files

module OutputNuShell

use Extra
use InputNuShell

implicit none

contains

    subroutine output_header(dev)
    implicit none
    integer,intent(in):: dev
    integer:: status
    character(len=14)::buffer
    call get_cmd_arg(1,buffer,status) 
    if (status /= -1 .and. dev==6) return
!   Intel
    include 'IHeader.FI'
!   end Intel
!   MSU
!   include 'MHeader.FI
!   end MSU
!   Sun
!   include 'SHeader.FI'
!   end Sun
    end subroutine output_header
    
    subroutine output_headerX(dev)
    implicit none
    integer,intent(in):: dev
    integer:: status
    character(len=14)::buffer
    call get_cmd_arg(1,buffer,status) 
    if (status /= -1 .and. dev==6) return
!   Intel
    include 'XHeader.FI'
!   end Intel
!   MSU
!   include 'MHeader.FI
!   end MSU
!   Sun
!   include 'SHeader.FI'
!   end Sun
    end subroutine output_headerX
    
    subroutine output_welcome(dev)
    implicit none
    integer,intent(in):: dev
    integer:: status
    character(len=14)::buffer
    call get_cmd_arg(1,buffer,status) 
    if (status /= -1) return
    
    write(dev,*) ' Welcome! Your input will now be read.'
    write(dev,*)
    write(dev,*)

    end subroutine output_welcome    
    
    subroutine output_shells(dev)
    implicit none
    integer,intent(in):: dev
    integer:: status
    character(len=14)::buffer
    call get_cmd_arg(1,buffer,status) 
    if (status /= -1) return
    
    write(dev,*) ' ',Nucleus
    write(dev,*)
    write(dev,*) ' ',Description(1)
    write(dev,*) ' ',Description(2)
    write(dev,*) ' ',Description(3)
    write(dev,*)
    write(dev,*) ' Number of ispin 1/2 nucleons',no_cI
    write(dev,*) ' Number of np formalism nucleons',no_cN
    write(dev,*) ' Parity            ',tParity
    write(dev,*) ' Number of shells  ',no_shells
    write(dev,*)
    write(dev,*) ' Shell No    N    n    l    2*j   2*tz    p    Nucleon'
    do shn=1,no_shells
        write(dev,'(i10,i5,i5,i5,i7,i7,i5,8x,a1)')  &
                shn,cN(shn),sn(shn),sl(shn),sj2(shn),st2(shn),sp(shn),nucleon(shn)
    end do
    write(dev,*)
    write(dev,*) ' Minimum and Maximum values of total (ie sum over all shells) for:'
    write(dev,*)
    write(dev,*) ' Angular Momentum x 2              ',min_cJ2,max_cJ2
    write(dev,*) ' Intermediate Isospin x 2 (if used), pn isospin ',min_cT2,max_cT2,nmp
    write(dev,*)
    write(dev,*) ' Individual shell maximum and minimum limits in each shell:'
    write(dev,*) ' Min no particles', (min_cSNo%shell(shn),shn=1,no_shells)
    write(dev,*) ' Max no particles', (max_cSNo%shell(shn),shn=1,no_shells)
    write(dev,*) ' Min no MajorSh  ', (min_MS_No%shell(shn),shn=1,no_shells)
    write(dev,*) ' Max no MajorSh  ', (max_MS_No%shell(shn),shn=1,no_shells)
    

    write(dev,*)
    write(dev,*) ' Output Control  ', output_control
    write(dev,*)
    write(dev,*) ' Your input file has been successfully read.'
    write(dev,*) ' Now setting up for calculations...'

    end subroutine output_shells    
    
    subroutine output_end_setup(dev)
    implicit none
    integer,intent(in):: dev
    integer:: status
    character(len=14)::buffer
    call get_cmd_arg(1,buffer,status) 
    if (status /= -1) return
    
    write(dev,*)
    write(dev,*) ' Setting up completed successfully.'
    write(dev,*) ' Now starting calculations...'
    write(dev,*)

    end subroutine output_end_setup
    
    subroutine output_completed(dev,next_prog)
    implicit none
    integer,intent(in)::dev
    character(len=10),intent(in):: next_prog
    integer:: status
    character(len=14)::buffer
    call get_cmd_arg(1,buffer,status) 
    if (status /= -1) return
    
    write(dev,*) 
    write(dev,*) ' Calculations successfully completed and data written to disk.'
    write(dev,*) ' Now run ',next_prog,'.'
    write(dev,*)
    
    end subroutine output_completed
    
    subroutine output_int_read(dev)
    implicit none
    integer,intent(in):: dev
    integer:: status
    character(len=14)::buffer
    call get_cmd_arg(1,buffer,status) 
    if (status /= -1) return
    
    write(dev,*) ' Your interaction file has been successfully read.'
    
    end subroutine output_int_read
    
    subroutine output_add_int(dev)
    implicit none
    integer,intent(in):: dev
    integer:: status
    character(len=14)::buffer
        call get_cmd_arg(1,buffer,status) 
        if (status /= -1) return
     
        write(dev,*)
        write(dev,*)
        write(dev,*) ' To add multiple interactions, commentenate files with'
        write(dev,*) ' at least one comment line starting with `!` between them.'
        write(dev,*) ' For the first int add a line `&` immed after comments'
        write(dev,*) ' For the first interaction set no of matrix elements '
        write(dev,*) ' to -n  where you want to add n interactions'
        write(dev,*) ' Give the normalisations when asked.'
        write(dev,*)
        write(dev,*)
        
    end subroutine output_add_int

    subroutine output_time(chr,t1,t2)
        implicit none
        character(len=*),intent(in):: chr
        real,intent(in):: t1,t2
        
        write(*,*) ' CPU Time for ',chr,' was ',t1,' seconds'
        write(*,*) ' Elapsed Time for ',chr,' was ',t2,' seconds'

     end subroutine output_time
    
end module OutputNuShell


