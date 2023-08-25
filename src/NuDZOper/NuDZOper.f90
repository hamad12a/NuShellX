! This file contains all the modules specific to generric program gShell
! for NuDZOper and SunShell only

module OutputNus

use Partition
use Extra

implicit none

character(len=6):: Nucleus
character(len=40),dimension(3):: Description
character(len=1):: comment,exclam,xtra
character(len=7),private:: filename
character(len=max_no_shells):: nucleon_line
integer:: no_cI,no_cN
integer:: tParity
integer:: min_cJ2,max_cJ2
integer:: min_MSon,max_MSon
integer:: min_cT2,max_cT2,nmp
integer:: output_control
integer:: shn,no_shells,max_sj2u
real(kind=8),dimension(max_no_shells)::  spe
real(kind=8):: A,A0,power,renorm
integer:: no_levels
integer,dimension(max_no_shells):: cN,sn,sl,sj2,st2,sp
character(len=1),dimension(max_no_shells):: nucleon
type(spartition):: min_cSNo,max_cSNo
type(spartition):: min_MS_No,max_MS_No
logical:: help

data exclam /'!'/
data filename /'NuDZOper'/

contains

    subroutine output_shell(inExt)
        implicit none
        character(len=6):: inExt
        
        open(unit=10,file=inExt//'.nus',action='write')
        
        write (10,'(a6)') Nucleus
! Intel
        write (10,'(a40)') 'NuDZOper is a nuclear shellmodel program '
        write (10,'(a40)') 'written in Fortran95 for Windows/Linux. '
        write (10,'(a40)') 'It is based on NuDZOper by W D M Rae 2007'
! end Intel
! Sun
!        write (10,'(a40)') 'SunShell is a nuclear shellmodel program'
!        write (10,'(a40)') 'written in Fortran95 for Solaris x86/Sun'
!        write (10,'(a40)') 'It is based on NuDZOper by W D M Rae 2007'
! end sun
        write (10,'(2i4)') no_cI,no_cN
        write (10,'(i4)') tParity
        write (10,'(i4)') no_shells
        write (10,'(100(a1))') nucleon(1:no_shells)
        do shn=1,no_shells
            write (10,'(6i4)') cN(shn),sn(shn),sl(shn),sj2(shn),st2(shn),sp(shn)
        end do
        write (10,'(2i4)') min_cJ2,max_cJ2
        write (10,'(3i4)') min_cT2,max_cT2,nmp
        write (10,'(2i4)') min_MSon,max_MSon
        write (10,'(100i4)') (min_cSNo%shell(shn),shn=1,no_shells)
        write (10,'(100i4)') (max_cSNo%shell(shn),shn=1,no_shells)
        write (10,'(100i4)') (min_MS_No%shell(shn),shn=1,no_shells)
        write (10,'(100i4)') (max_MS_No%shell(shn),shn=1,no_shells)
        write (10,'(i4)') output_control
        write (10,*)
               
        close(unit=10)
        
    end subroutine output_shell
       
    subroutine spe_file_shell(inExt)
        implicit none
        character(len=7):: inExt

        open(unit=10,file=inExt//'.spe',action='write')
        write (10,'(100f10.4)') spe(1:no_shells)
        write (10,'(f14.10)') renorm
        write (10,'(i4)') no_levels
        close(unit=10)
        
    end subroutine spe_file_shell
        
end module OutputNus

module InputModelSpace

use Partition
use OutputNus

implicit none

contains

    subroutine input_model_shell(inExt)
        implicit none
        character(len=6):: inExt

        open(unit=10,file=inExt//'.sps',status='old',action='read',err=20)
        
        call read_comments_model()
        read (10,*,err=30,end=31) no_shells
        if (no_shells > max_no_shells) then
            print *, ' Number of shells exceeds parameter max_no_shells in module Parameters.'
            print *, ' Increase this parameter and recompile.'
            stop
        end if
        call set_ptdim(no_shells)
        
        call read_comments_model()
        read (10,*,err=30,end=31) nucleon_line
        max_sj2u=0
        do shn=1,no_shells
            nucleon(shn)=nucleon_line(shn:shn)
            call read_comments_model()
            read(10,*,err=30,end=31) cN(shn),sn(shn),sl(shn),sj2(shn),st2(shn),sp(shn)
            if (sj2(shn) > max_sj2u) max_sj2u=sj2(shn)
        end do
        call read_comments_model()
        read (10,*,err=30,end=31) min_MSon,max_MSon
        call read_comments_model()
        read (10,*,err=30,end=31) (min_cSNo%shell(shn),shn=1,no_shells)
        call read_comments_model()
        read (10,*,err=30,end=31) (max_cSNo%shell(shn),shn=1,no_shells)
        call read_comments_model()

        read (10,*,err=30,end=31) (min_MS_No%shell(shn),shn=1,no_shells)
        call read_comments_model()
        read (10,*,err=30,end=31) (max_MS_No%shell(shn),shn=1,no_shells)
        
        do shn=1,no_shells
            if (max_cSNo%shell(shn) > (sj2(shn)+1)*(st2(shn)+1)) max_cSNo%shell(shn)=(sj2(shn)+1)*(st2(shn)+1)
            if (min_cSNo%shell(shn) > (sj2(shn)+1)*(st2(shn)+1)) min_cSNo%shell(shn)=(sj2(shn)+1)*(st2(shn)+1)
        end do
        
        close(unit=10)
        return
20      print *, ' Error opening ',inExt//'.sps'
        stop
30      print *, ' Error reading ',inExt//'.sps'
        stop
31      print *, ' EOF reading ',inExt//'.sps'
        stop
                
        
    end subroutine input_model_shell
   
    subroutine read_comments_model(more)
        implicit none
        logical,optional,intent(out):: more
        if (present(more))more=.false.
        do
            read(10,'(a1)',err=10,end=10)comment
            if (present(more)) then
                if (comment == exclam) more=.true.
            end if
            if (.not.(comment == exclam)) exit
        end do
        if (comment == '&' .or. comment == '^') return
        backspace(unit=10)
10      return
    
    end subroutine read_comments_model
    
    subroutine input_nucleus_shell(inExt,xtra)
        implicit none
        character(len=6),intent(out):: InExt
        character(len=1),intent(out),optional:: xtra
        integer:: status
        write(6,'(a10)',advance='no') '  Nucleus '
        read(5,'(a6)') inExt
        inExt=adjustr(inExt)
        if (present(xtra)) then
            write(6,'(a19)',advance='no') '  Interaction code '
            read(5,'(a1)') xtra
        end if
        
    end subroutine input_nucleus_shell

    subroutine input_nucleus_final_shell(inExt,xtra)
        implicit none
        character(len=6),intent(out):: InExt
        character(len=1),intent(out),optional:: xtra
        integer:: status
        write(6,'(a10)',advance='no') '  Nucleus '
        read(5,'(a6)') inExt
        inExt=adjustr(inExt)
        if (present(xtra)) then
            write(6,'(a19)',advance='no') '  Interaction code '
            read(5,'(a1)') xtra
        end if
        
    end subroutine input_nucleus_final_shell

    subroutine input_model_space_shell(inExt,xtra)
        implicit none
        character(len=6),intent(out):: InExt
        character(len=1),intent(out),optional:: xtra
        integer:: status
        write(6,'(a14)',advance='no') '  Model Space '
        read(5,'(a6)') inExt
        inExt=adjustr(inExt)
        if (present(xtra)) then
            write(6,'(a19)',advance='no') '  Interaction code '
            read(5,'(a1)') xtra
        end if
        
    end subroutine input_model_space_shell

    subroutine input_intfile_shell(intfile,no)
        implicit none
        character(len=7),intent(out):: intfile
        integer,intent(in),optional:: no
        integer::status
            if (.not.present(no))write(6,'(a19)',advance='no') '  Interaction file '
            if (present(no).and. no==1)write(6,'(a28)',advance='no') '  Name for transition files '
            if (present(no).and. no==2)write(6,'(a15)',advance='no') '  MBME op file '
            if (present(no).and. no==3)write(6,'(a19)',advance='no') '  MBME matrix file '
            read(5,'(a7)') intfile
            intfile=adjustr(intfile)
        
    end subroutine input_intfile_shell
    
    subroutine read_spe_shell(intfile)
        implicit none
        character(len=7),intent(in):: intfile
        character(len=1) excl
        integer:: idummy
        
        open (unit=10,file=intfile//'.int',err=20)
        read(10,*,err=30,end=31) excl,A0,A,power
        call read_comments_model()
        read (10,*,err=30,end=31) idummy,spe(1:no_shells)
        close(unit=10)
        return
20      print *,' Error opening ',intfile//'.int'
        stop
30      print *,' Error reading ',intfile//'.int'
        stop
31      print *,' EOF reading ',intfile//'.int'
        stop
        
    end subroutine read_spe_shell

end module InputModelSpace

module Levels
use extra
use InputModelSpace
use OutputNus

implicit none
!  Intel
!  end Intel
!  Sun
!  include 'system.inc'
!  end Sun
integer:: nI,nP,nN,I,errnum

contains

    subroutine Oper()
        implicit none
        character(len=6):: model
        character(len=7):: intfile
        character(len=15):: ccute,ccutm
        character(len=1):: delet,edit
        character(len=2):: npi
        integer(kind=2)::res
        real:: A00,A1,B,C
        real(kind=8):: cute,cutm
        nI=0
        nP=0
        nN=0
        min_cT2=0
        max_cT2=0
        output_control=0
        if (help) then
        write(6,*) ' If you are using a Sun model space and Nucleus filenames'
        write(6,*) ' six characters long. Interaction and transition filenames'
        write(6,*) ' must be seven characters long.'
        write(6,*) ' The interaction code is a single character of your choice.'
        end if
        xtra='x'
        call input_model_space_shell(model)
        call input_model_shell(model)
        if (nucleon(1)=='i') then
            write(6,*) ' Isospin model space '
        else
            write(6,*) ' np model space not yet supported '
            stop
        end if
        if (help) then
        write(6,*)
        write(6,*) ' If you use an Oxbash .int file add a line at top as follows:'
        write(6,*) ' ! Ai Ac pwr, eg for sd shell Ai=18, Ac=16 and pwr=0.3 for USD.'
        write(6,*) ' Remove spaces before ! ;shell order must agree with model space.'
        write(6,*) ' Here Ac is the mass of the core, which is the mass of the nucleus'
        write(6,*) ' that would be created if you put zero particles in the model space.'
        write(6,*) ' Ai is the mass of the core +2 the mass for a two-body calculation.'
        write(6,*) ' Finally pwr is the power of (Ac/(Ac+n))**pwr used to renormalise'
        write(6,*) ' the interaction, where n is the no of particles you input to NuDZOper.'
        write(6,*) ' For no renormalisation set Ac=1.0 Ai=1.0 pwr=0.0 ,since x**0.0=1.0.'
        write(6,*) ' If you use an Oxbash .int file from Windows convert it to a Unix file'
        write(6,*) ' if running SunShell on a Sun system - its all about <CR><Lf>`s etc'
        write(6,*)
        write(6,*) ' If you do not have a file of matrix elements use type none.'
        write(6,*)
        end if
        call input_intfile_shell(intfile)
        if (help) write(6,*) ' Do you want an output interaction in neutron-proton or isospin form?'
        write(6,'(a34)',advance='no') '  Output: (np) or (i) interaction '
        read(5,'(a2)') npi
        npi=adjustr(npi)
        cute=0.0
        cutm=0.0
        intfile=adjustr(intfile)
        write(6,'(a33)',advance='no') '  Do you want to remove Hcm y/n  '
        read(5,*) delet
        if (delet == 'y') then
            delet='h'
        else
        write(6,'(a30)',advance='no') '  Do you want to add SU3 y/n  '
        read(5,*) delet
        if (delet == 'y') then
            delet='s'
        end if    
        end if
        write(6,'(a38)',advance='no') '  Do you want to edit e and evecs y/n '
        read(5,*) edit
        if (edit == 'y') then
            edit='e'
        end if
        if (.not.edit=='e') then
        write(6,'(a37)',advance='no') '  Do you want to cut e and evecs y/n '
        read(5,*) edit
        if (edit == 'y') then
            edit='c'
        end if
        if (edit=='c') then
        write(6,'(a36)',advance='no') '  e cut (MeV), evec cut (amplitude) '
        read(5,*) cute,cutm
        end if
        end if 
        write(ccute,'(f15.10)') cute
        write(ccutm,'(f15.10)') cutm
        res=run('NuDZop',model//' '//xtra//' '//intfile//' '//npi//' '//delet//' '&
            //edit//' '//ccute//' '//ccutm)
        write(56,*) 'NuDZop '//model//' '//xtra//' '//intfile//' '//npi//' '//delet//' '&
            //edit//' '//ccute//' '//ccutm
        write(6,*) ' New interaction saved to ',intfile//'z.int'

        return
        
    end subroutine Oper
    

end module Levels    

module ModelSpace

use InputModelSpace

implicit none

contains

    subroutine CreateModel()
        implicit none
        logical:: use_isospin
        integer:: no_shells,cN,sn,sl,sj2,st2,sp,no_proton,no_neutron
        integer:: shn,min,max,sfl,nms,shn2
        character(len=1):: ans,ai='i',an='n',ap='p',chl
        character(len=6):: chmodel,chsh
        character(len=2):: ignore2
        character(len=3):: ignore3
        integer,dimension(50)::maxa,mina,maxna,minna,cNa,cNs
        character(len=1),dimension(10):: lch
        logical:: btest
        
        data lch /'s','p','d','f','g','h','i','j','k','l'/
        
        maxa=99
        mina=0
        maxna=99
        minna=0
        cNa=999
        cNs=0
        if (help) then
        write(6,*) ' If you are using a Sun model space and Nucleus filenames'
        write(6,*) ' six characters long. Interaction and transition filenames'
        write(6,*) ' must be seven characters long.'
        end if        
        call input_model_space_shell(chmodel)
        open(unit=11,file=chmodel//'.sps',action='write')
        chmodel=adjustr(chmodel)
        write(6,'(1x,a31)',advance='no') ' (I)sospin or p(n) calculation '
        read(5,'(a1)') ans
        use_isospin=.false.
        if (ans=='i' .or. ans=='I') use_isospin=.true.
        if (.not.use_isospin) then
            write(6,*) ' NuDZOper only works with isospin interactions at present.'
            stop
        end if    
10      write(6,'(1x,a20)',advance='no') ' Total no of shells '
        read(5,*) no_shells
        if (use_isospin) then
            write(11,'(i4)') no_shells
            do shn=1,no_shells
                write(11,'(a1)',advance='no') ai
            end do
            write(11,*)
        else
            write(6,'(1x,a21)',advance='no') ' No of proton shells '
            read(5,*) no_proton
            no_neutron=no_shells-no_proton
            if (no_neutron<0) then
                write(6,*) ' No of proton shells < Total no shells'
                go to 10
            end if    
            write(11,'(i4)') no_shells
            do shn=1,no_proton
                write(11,'(a1)',advance='no') ap
            end do
            do shn=1,no_neutron
                write(11,'(a1)',advance='no') an
            end do
            write(11,*)
        end if
        if (help) then
        write(6,*)
        write(6,*) ' If you use an Oxbash .int file get shell order from that file.'
        write(6,*) ' NuDZOper uses the convention that n starts at 0.'
        write(6,*) ' First few shells are 0s1/2 0p3/2 0p1/2 0d5/2 1s1/2 ... etc.'
        write(6,*) ' Please use this convention and lower case.'
        write(6,*)
        end if
        if (help) write(6,*) ' Enter proton shells first.'
        do shn=1,no_shells
20          write(6,'(1x,a6,i3,a1)',advance='no') ' Shell',shn,' '
            read(5,'(a6)') chsh
            chsh=adjustl(chsh)
            if (.not.chsh(6:6)=='2') then
            read(chsh,'(i1,a1,i1,a3)') sn,chl,sj2,ignore3           
            else
            read(chsh,'(i1,a1,i2,a2)') sn,chl,sj2,ignore2           
            end if
            do sfl=1,10
                if (lch(sfl)==chl) go to 30
            end do
            write(6,*) ' Spectroscopic symbol not recognised ',chl
            go to 20
30          sl=sfl-1
            sp=1
            if (btest(sl,0)) sp=-sp 
            if (use_isospin) then
                st2=1 
            else
                st2=0
            end if 
            cN=2*sn+sl
            cNa(shn)=cN
            write(11,'(6i4)') cN,sn,sl,sj2,st2,sp
        end do
        write(6,'(1x,a41)',advance='no') ' Do you want to truncate model space y/n '
        read(5,'(a1)') ans
        if (ans=='y') then
            nms=0
            do shn=1,no_shells
                min=minval(cNa)
                if (min==999) exit
                nms=nms+1
                cNs(nms)=min
                    do  shn2=1,no_shells
                         max=maxval(minloc(cNa))
                         if (cNa(max)>min) exit
                         cNa(max)=999
                    end do
            end do
            write(6,'(1x,a47)',advance='no') ' Do you want to limit the sum over N=2*n+l y/n '
            read(5,'(a1)') ans
                if (ans=='y') then
                    write(6,'(1x,a9)',advance='no') ' Min Max '
                    read(5,*) min,max
                    write(11,'(2i4)') min,max
                else
                    min=0
                    max=999
                    write(11,'(2i4)') min,max
                end if
            write(6,'(1x,a58)',advance='no') ' Do you want to limit the no of nucleons in any shell y/n '
            read(5,'(a1)') ans
40          if (ans=='y') then
                write(6,'(1x,a17)',advance='no') ' ShellNo Min Max '
                read(5,*) shn,mina(shn),maxa(shn)
            else 
                go to 50
            end if
            write(6,'(1x,a10)',advance='no') ' More y/n ' 
            read(5,'(a1)') ans
            if (ans=='y') go to 40
50          if (help) then
            write(6,*) ' A major shell (MS) is defined as shells with equal 2*n+l.'
            write(6,*)
            end if
            write(6,'(1x,a51)',advance='no') ' Do you want to limit no of nucleons in any MS y/n '
            read(5,'(a1)') ans
            if (ans=='y') then
                write(6,'(1x,a10,i4,a14)') ' There are',nms,' major shells.'
                write(6,'(1x,a4,50i5)') ' N= ',cNs(1:nms)
60              write(6,'(1x,a18)',advance='no') ' MShellNo Min Max '
                read(5,*) shn,mina(shn),maxna(shn)
                write(6,'(1x,a10)',advance='no') ' More y/n '             
                read(5,'(a1)') ans
                if (ans=='y') go to 60
            end if
            write(11,'(50i5)') mina(1:no_shells)
            write(11,'(50i5)') maxa(1:no_shells)
            write(11,'(50i5)') minna(1:no_shells)
            write(11,'(50i5)') maxna(1:no_shells)
        else
            min=0
            max=999
            write(11,'(2i4)') min,max
            write(11,'(50i5)') mina(1:no_shells)
            write(11,'(50i5)') maxa(1:no_shells)
            write(11,'(50i5)') minna(1:no_shells)
            write(11,'(50i5)') maxna(1:no_shells)
        end if
        write(11,'(i2)') 0
        write(11,*)
        close(unit=11) 
        write(6,*) ' Model Space data saved as ',chmodel//'.sps'                        
        
    end subroutine CreateModel 
    
end module ModelSpace    


!  NuDZOper.f90 
!
!  FUNCTIONS:
!  NuDZOper      - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: NuDZOper
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program NuDZOper
    
    use OutputNus
    use Levels
    use ModelSpace
    use Extra
    use OutputNuShell

    implicit none

    ! Variables
    character(len=1):: ans
    character(len=4):: buffer
    integer:: status
    
    ! Body of NuDZOper
    help=.false.
    call get_cmd_arg(1,buffer,status)
    if (status>0) then
        if (buffer=='help') then
            help=.true.
        end if
    end if
    
    call output_header(6)
    write(6,*) 
    write(6,*) ' To get help messages use switch help ie run `NuDZOper help`.'
    write(6,*) ' NuDZOper is a program which can analyse & modify interactions using'
    write(6,*) ' the formalism of Marianne Defour and Andreas Zuker, PRC 54, No 4,'
    write(6,*) ' pp1641-1659, October 1996.  It can also be used to remove the'
    write(6,*) ' linear combinations of matrix elements that couple spurious and'
    write(6,*) ' non spurious states. For more information see this reference and the'
    write(6,*) ' NuShellX paper. Currently it only works with an iso input interaction'
    if (help) then
    write(6,*)
    write(6,*) ' To run NuDZOper you need to decide on a model space and interaction.'
    write(6,*) ' If you use an Oxbash .int file check that for shells and their order.'
    write(6,*) ' Otherwise you will need to make a .int file from data in literature.'
    write(6,*) ' To do this use an Oxbash file as an example.'
    write(6,*) 
    write(6,*) ' First choose the option Model to define your model space .sps file.'
    write(6,*) ' After this or if you use an existing model space choose option Oper.'
    end if
    open(unit=56,file='DZOper.bat',action='write')
        
10  write(6,*)
    write(6,'(a24)', advance='no') '  (M)odel,(O)per,E(x)it '
    read (5,'(a1)') ans
    write(6,*)
    
    if (ans == 'm' .or. ans == 'M') call CreateModel()
    if (ans == 'O' .or. ans == 'o') call Oper()
    if (ans == 'x' .or. ans == 'X') go to 100
    
    go to 10
    
100 write(6,*) ' Command file ','DZOper.bat',' saved to disk.'
    if (help) write(6,*) ' Rename this file if you want to repeat calculations.'
    if (help) write(6,*) ' Otherwise delete it.'
    close(unit=56)

    end program NuDZOper

