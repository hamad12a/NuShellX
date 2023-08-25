! This file contains all the modules specific to generric program gShell
! for NuShell and SunShell only

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
data filename /'NuShell'/

contains

    subroutine output_shell(inExt)
        implicit none
        character(len=6):: inExt
        
        open(unit=10,file=inExt//'.nus',action='write')
        
        write (10,'(a6)') Nucleus
! Intel
        write (10,'(a40)') 'NuShell is a nuclear shellmodel program '
        write (10,'(a40)') 'written in Fortran95 for Windows/Linux. '
        write (10,'(a40)') 'It is based on NuShell by W D M Rae 2007'
! end Intel
! Sun
!        write (10,'(a40)') 'SunShell is a nuclear shellmodel program'
!        write (10,'(a40)') 'written in Fortran95 for Solaris x86/Sun'
!        write (10,'(a40)') 'It is based on NuShell by W D M Rae 2007'
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

    subroutine Level()
        implicit none
        character(len=6):: model
        character(len=7):: intfile
        character(len=1):: delet
        character(len=2):: npi
        integer(kind=2)::res
        real:: A00,A1,B,C
        
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
        call input_nucleus_shell(Nucleus,xtra)
        call input_model_space_shell(model)
        call input_model_shell(model)
        if (nucleon(1)=='i') then
            write(6,'(a17)',advance='no') '  No of nucleons '
            read(5,*) nI
            no_cI=nI
        else
            write(6,'(a16)',advance='no') '  No of protons '
            read(5,*) nP
            write(6,'(a17)',advance='no') '  No of neutrons '
            read(5,*) nN
            no_cN=nP+nN
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
        write(6,*) ' the interaction, where n is the no of particles you input to NuShell.'
        write(6,*) ' For no renormalisation set Ac=1.0 Ai=1.0 pwr=0.0 ,since x**0.0=1.0.'
        write(6,*) ' If you use an Oxbash .int file from Windows convert it to a Unix file'
        write(6,*) ' if running SunShell on a Sun system - its all about <CR><Lf>`s etc'
        write(6,*)
        write(6,*) ' If you do not have a file of matrix elements use type none.'
        write(6,*)
        end if
        call input_intfile_shell(intfile)
        if (help) write(6,*) ' Is this an interaction in neutron-proton or isospin form?'
        write(6,'(a26)',advance='no') '  (np) or (i) interaction '
        read(5,'(a2)') npi
        npi=adjustr(npi)

        intfile=adjustr(intfile)
        if (intfile(4:7)=='none' .or. intfile(4:7)=='NONE') then
            if (help) then
            write(6,*)
            write(6,*) ' The MSDI interaction requires 4 parameter A0 A1 B C.'
            write(6,*) ' Parameters depend on the A of your nucleus as 1/A**(1/3).'
            write(6,*) ' For 16O A0~A1~1.0 , B is smaller and C often zero.'
            write(6,*)
            end if
            write(6,'(1x,a27)',advance='no') ' MSDI Parameters A0 A1 B C '
            read(5,*) A00,A1,B,C
            if (help) then
            write(6,*)
            write(6,*) ' NuShell needs single particle energies for shells (MeV).'
            end if
            write(6,'(1x,a15,i4)') ' No of shells =',no_shells
            write(6,'(1x,a7)',advance='no') ' SPE`s '
            read(5,*) spe(1:no_shells) 
            renorm=1.0
            open(unit=11,file='MSDI'//'.int',action='write')
            write(11,'(a1,1x,3f8.3)') '!',1.0,1.0,1.0
            write(11,'(a8)') '^sdip2p2'
            write(11,'(i4,50f10.6)') 0,spe(1:no_shells)
            write(11,'(11f10.6)') A00,A1,B,C,-1.0,1.0,0.0,0.0,0.0,0.0,0.0
            write(11,'(5f10.6)') 0.0,0.0,0.0,0.0,0.0
            write(11,*)
            write(11,*) ' Brussaard and Glaudemans MSDI interaction see text book'
            close(unit=11)
            intfile='  MSDI'
        else           
            call read_spe_shell(intfile)
            A=A+real(no_cI+no_cN,8)
            renorm=(A0/A)**power
        end if
        write(6,'(a26)',advance='no') '  Minumum and Maximum 2*J '
        read(5,*) min_cJ2,max_cJ2
        if (nucleon(1)=='i') then
            write(6,'(a26)',advance='no') '  Minumum and Maximum 2*T '
            read(5,*) min_cT2,max_cT2
        else 
            nmp=nN-nP
        end if
        write(6,'(a15)',advance='no') '  Parity +1,-1 '
        read(5,*) tParity
        write(6,'(a15)',advance='no') '  No of levels '
        read(5,*) no_levels
        if (help) then
            write(6,*)
            write(6,*) ' It is advised that you delete existing .lev files.'
            write(6,*) ' If you want to save them copy them to another directory.'
            write(6,*)
        end if
        write(6,'(a33)',advance='no') '  Delete existing .lev files y/n '
        read(5,*) delet
        if (delet == 'y') then
            I = deletf('*.lev','sunshell.del')
        end if
        call output_shell(Nucleus)
        write(6,*) ' Data saved as ',Nucleus//'.nus'
        call spe_file_shell(Nucleus//xtra)
        write(6,*) ' SPE file created ',Nucleus//xtra//'.spe'
        write(6,*)
        res=run(pBasis,Nucleus)
        write(56,'(a7,a1,a6)') pBasis,' ',Nucleus
        if (res == -1) stop
        res=run(pProj,Nucleus)
        write(56,'(a6,a1,a6)') pProj,' ',Nucleus
        if (res == -1) stop
        res=run(pOper,Nucleus//' '//xtra//' '//intfile//' '//npi)
        write(56,'(a6,a1,a16)')pOper,' ',Nucleus//' '//xtra//' '//intfile//' '//npi
        if (res == -1) stop
        res=run(pMatrix,Nucleus//' '//xtra)
        write(56,'(a8,a1,a8)') pMatrix,' ',Nucleus//' '//xtra
        if (res == -1) stop
        res=run(pLanczos,Nucleus//' '//xtra)
        if (res == 1111) then
            res=run(pOrth,Nucleus//' '//xtra)  
            write(56,'(a9,a1,a8)') pOrth,' ',Nucleus//' '//xtra
            if (res == -1) stop
            res=run(pBlock,Nucleus//' '//xtra)  
            write(56,'(a9,a1,a8)') pBlock,' ',Nucleus//' '//xtra
            if (res == -1) stop
        else
            write(56,'(a9,a1,a8)') pLanczos,' ',Nucleus//' '//xtra
            if (res == -1) stop
        end if
        res=run(pMvec,Nucleus//' '//xtra)
        write(56,'(a6,a1,a8)') pMvec,' ',Nucleus//' '//xtra
        if (res == -1 )stop
        res=run(pPSLevel,'   ')
        write(56,'(a9,a2)') pPSLevel,' -'
        if (res == -1) stop
        write(6,*) ' Energy level diagram output as ',Nucleus//'.eps'        
        return
        
    end subroutine Level
    
    subroutine Interaction()
        implicit none
        character(len=6):: model
        character(len=7):: intfile
        character(len=1):: delet
        integer:: nI,nP,nN,I,errnum
        integer(kind=2)::res
        real:: A00,A1,B,C

        nI=0
        nP=0
        nN=0
        min_cT2=0
        max_cT2=0
        output_control=0
        call input_nucleus_shell(Nucleus,xtra)
        call input_model_space_shell(model)
        call input_model_shell(model)
        call input_intfile_shell(intfile)
        intfile=adjustr(intfile)
        if (intfile(1:4)=='none' .or. intfile(1:4)=='NONE') then
            read(5,*) A00,A1,B,C
            write(6,'(1x,a15,i4)') ' No of shells =',no_shells
            write(6,'(1x,a7)',advance='no') ' SPE`s '
            read(5,*) spe(1:no_shells) 
            renorm=1.0
            open(unit=11,file='MSDI'//'.int',action='write')
            write(11,'(a1,1x,3f8.3)') '!',1.0,1.0,1.0
            write(11,'(a8)') '^sdip2p2'
            write(11,'(i4,50f10.6)') 0,spe(1:no_shells)
            write(11,'(11f10.8)') A00,A1,B,C,-1.0,1.0,0.0,0.0,0.0,0.0,0.0
            write(11,'(5f10.8)') 0.0,0.0,0.0,0.0,0.0
            write(11,*)
            write(11,*) ' Brussaard and Glaudemans MSDI interaction see text book'
            close(unit=11)
            intfile='  MSDI'
            write(6,*) ' Interaction file saved as MSDI.int.'
        else           
            call read_spe_shell(intfile)
            A=A+real(no_cI+no_cN,8)
            renorm=(A0/A)**power
        end if
        write(6,'(a15)',advance='no') '  No of levels '
        read(5,*) no_levels
        write(6,'(a33)',advance='no') '  Delete existing .lev files y/n '
        read(5,*) delet
        call spe_file_shell(Nucleus//xtra)
        write(6,*)
        if (delet == 'y') then
            I = deletf('*.lev','sunshell.del')
        end if
        res=run(pOper,Nucleus//' '//xtra//' '//intfile)
        write(56,'(a6,a1,a16)')pOper,' ',Nucleus//' '//xtra//' '//intfile
        if (res == -1) stop
        res=run(pMatrix,Nucleus//' '//xtra)
        write(56,'(a8,a1,a8)') pMatrix,' ',Nucleus//' '//xtra
        if (res == -1) stop
        res=run(pLanczos,Nucleus//' '//xtra)
        write(56,'(a9,a1,a8)') pLanczos,' ',Nucleus//' '//xtra
        if (res == -1) stop
        res=run(pMvec,Nucleus//' '//xtra)
        write(56,'(a6,a1,a8)') pMvec,' ',Nucleus//' '//xtra
        if (res == -1 )stop
        res=run(pPSLevel,'   ')
        write(56,'(a9,a2)') pPSLevel,' -'
        if (res == -1) stop
        write(6,*) ' Energy level diagram output as ',Nucleus//'.eps'        
        return
        
    end subroutine Interaction

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
        write(6,*) ' NuShell uses the convention that n starts at 0.'
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

module Transitions

use Extra
use InputModelSpace
use OutputNus
use Levels

implicit none

real::radial=1.0

contains

    subroutine Transition()
        implicit none
        character(len=6):: Nucleusi,Nucleusf,model
        character(len=1):: xtrai,xtraf,delet,ans
        character(len=7):: mefile
        integer(kind=2):: res
        integer:: min_Ji,max_Ji,min_Ti,max_Ti,min_Jf,max_Jf,min_Tf,max_Tf
        integer:: min_Jo,maxJo,min_To,max_To
        integer:: mnDJ,mxDJ,mnDT,mxDT,mnDM,mxDM,mnDTz,mxDTz,ni,nf,i,j,amnDM,amnDTz
        logical:: OBTD,SP,TP
        real:: ep=0.0,en=0.0,glp=0.0,gln=0.0,gsp=0.0,gsn=0.0,mass=0.0,charge=0.0
        
        mnDM=999
        amnDm=999
        mxDM=-999
        mnDTz=999
        amnDTz=999
        mxDTz=-999
        
        if (help) then
        write(6,*) ' If you are using a Sun model space and Nucleus filenames'
        write(6,*) ' six characters long. Interaction and transition filenames'
        write(6,*) ' must be seven characters long.'
        end if        
        call input_model_space_shell(model)
        call input_model_shell(model)
        
        write(6,*) ' Initial States'
        call input_nucleus_shell(Nucleusi,xtrai)
        call get_JT_values(min_Ji,max_Ji,min_Ti,max_Ti,'Initial ')
        write(6,*) ' Final States'
        call input_nucleus_final_shell(Nucleusf,xtraf)
        call get_JT_values(min_Jf,max_Jf,min_Tf,max_Tf,'Final   ')
        if (help) then
        write(6,*) ' The transition filename is of your own choice'
        end if
        call input_intfile_shell(mefile,1)
        if (help) then
        write(6,*) ' For gamma decay or beta decay choose OBTD.'
        write(6,*) ' 2*J must be even for OBTD and TP, odd for SP.'
        end if
        call get_operator(OBTD,SP,TP)
        call get_JT_values(mnDJ,mxDJ,mnDT,mxDT,'Operator')
        if (OBTD) then
            if (help) then
            write(6,*)
            write(6,*)' For gamma decay NuShell require `effective charges and g-facctors`.'
            write(6,*)' These help correct for effects such as `core polarisation`.'
            write(6,*)' See textbooks or literature for more information.'
            write(6,*)' Typical input A == mass of nucleus Z == charge of nucleus.'
            write(6,*)' ep=1.5, en=0.5, glp=1.0, gln=0.0, gsp=5.586, gsn=-3.826'
            write(6,*)' Enter 0 0 0 0 0 0 for default. If you want en=0.0 enter 999.0.'
            write(6,*)
            end if
            write(6,'(1x,a5)',advance='no') ' A Z '
            read(5,*) mass,charge
            write(6,'(1x,a29)',advance='no') ' Enter ep en glp gln gsp gsn '
            read(5,*) ep,en,glp,gln,gsp,gsn
            if (ep==0.0) ep=1.5
            if (en==0.0) en=0.5
            if (en==999.0) en=0.0
            if (glp==0.0) glp=1.0
            if (gsp==0.0) gsp=5.586
            if (gsn==0.0) gsn=-3.826
        end if
        if (OBTD .or. TP) then
            mnDT=0
            mxDT=2
        end if
        if (SP) then
            mnDT=1
            mxDT=1
        end if
        write(6,'(1x,a24)',advance='no') ' Number of states ni nf '
        read(5,*) ni,nf
        if (ni==0) ni=1
        if (nf==0) nf=1
10      do i=min_Ji,max_Ji,2
            do j=min_Jf,max_Jf,2
                mnDM=min(mnDM,j-i)
                mxDM=max(mxDM,j-i)
                amnDM=min(amnDM,abs(j-i))
            end do
            if (mxDJ<abs(amnDM)) then
            write(6,*) ' The minimum |2*M| value of operator is greater than maximum 2*J'
            call get_JT_values(min_Ji,max_Ji,min_Ti,max_Ti,'Initial ')
            call get_JT_values(min_Jf,max_Jf,min_Tf,max_Tf,'Final   ')
            call get_JT_values(mnDJ,mxDJ,mnDT,mxDT,'Operator')
            mnDM=999
            amnDm=999
            mxDM=-999
            go to 10
            end if
        end do
20      do i=min_Ti,max_Ti,2
            do j=min_Tf,max_Tf,2
                mnDTz=min(mnDTz,j-i)
                mxDTz=max(mxDTz,j-i)
                amnDTz=min(amnDTz,abs(j-i))
            end do
            if (mxDTz<abs(amnDTz)) then
            write(6,*) ' The minimum |2*Tz| value of operator is greater than the maximum 2'
            call get_JT_values(min_Ji,max_Ji,min_Ti,max_Ti,'Initial ')
            call get_JT_values(min_Jf,max_Jf,min_Tf,max_Tf,'Final   ')
            call get_JT_values(mnDJ,mxDJ,mnDT,mxDT,'Operator')
            mnDTz=999
            amnDTz=999
            mxDTz=-999
            go to 20
            end if
        end do
        if (help) then
            write(6,*)
            write(6,*) ' Matrix elements depend on the sign convention for radial wavefunctions.'
            write(6,*) ' NuShell uses by default positive at origin.'
            write(6,*) ' To change this to positive  at infinity reply y to next question.'
            write(6,*)
        end if
        write(6,'(1x,a43)',advance='no') ' Change wavefunctions to + at infinity y/n '
        read(5,'(a1)') ans
        if (ans=='y' .or. ans=='Y') radial=-1.0
        open(unit=11,file=mefile//'.me',action='write')
        write(11,'(a1,1x,a6,1x,a6,1x,a7)') '!',Nucleusi,Nucleusf,mefile
        write(11,'(6f8.4,2f8.2)') ep,en,glp,gln,gsp,gsn,mass,charge
        write(11,'(10i4,f10.6)') mnDJ,mxDJ,mnDT,mxDT,mnDM,mxDM,mnDTz,mxDTz,ni,nf,radial
        if (OBTD) call obtd_out()
        if (SP)   call sp_out()
        if (TP)   call tp_out()
        close(unit=11)
        if (help) then
            write(6,*)
            write(6,*) ' It is advised that you delete existing .trs files.'
            write(6,*) ' If you want to save them copy them to another directory.'
            write(6,*)
        end if
        write(6,'(a33)',advance='no') '  Delete existing .trs files y/n '
        read(5,*) delet
        if (delet == 'y') then
            I=deletf('*.trs','sunshell.del')
        end if
        write(6,*) ' Data saved to ',mefile//'.me'
        
        res=run(pMe,mefile)
        write(56,'(a4,a1,a7)') pMe,' ',mefile
        if (res == -1 )stop
        res=run(pMbme,Nucleusi//' '//Nucleusf//' '//mefile)
        write(56,'(a6,a1,a21)') pMbme,' ',Nucleusi//' '//Nucleusf//' '//mefile
        if (res == -1 )stop
        res=run(pTramp,Nucleusi//' '//xtrai//' '//Nucleusf//' '//xtraf//' '//mefile)
        write(56,'(a7,a1,a25)') pTramp,' ', &
                  Nucleusi//' '//xtrai//' '//Nucleusf//' '//xtraf//' '//mefile
        if (res == -1 )stop
        res=run(pTrans,Nucleusi//' '//xtrai//' '//Nucleusf//' '//xtraf//' '//mefile)
        write(56,'(a7,a1,a25)') pTrans,' ', &
                  Nucleusi//' '//xtrai//' '//Nucleusf//' '//xtraf//' '//mefile
        if (res == -1 )stop
        write(6,*) ' Transition Data saved to ',mefile//'.lp'
        
    end subroutine Transition
    
    subroutine get_JT_values(mnj,mxj,mnt,mxt,text)
        implicit none
        integer, intent(out):: mnj,mxj,mnt,mxt
        character(len=8), intent(in):: text
        integer:: mmm
        
        write(6,'(2x,a8,a25)',advance='no') text,' mimimum and maximum 2*J '
        read(5,*) mnj,mxj
        mnj=abs(mnj)
        mxj=abs(mxj)
        if (mnj >mxj ) then
            mmm=mxj
            mxj=mnj
            mnj=mmm
        end if
        if (text=='Operator') return
        write(6,'(2x,a8,a25)',advance='no') text,' mimimum and maximum 2*T '
        read(5,*) mnt,mxt
        mnt=abs(mnt)
        mxt=abs(mxt)
        if (mnt >mxt ) then
            mmm=mxt
            mxt=mnt
            mnt=mmm
        end if
        
    end subroutine get_JT_values
    
    subroutine get_operator(o,s,t)
        implicit none
        logical, intent (out):: o,s,t
        character(len=1):: ans
        
        o=.false.
        s=.false.
        t=.false.
1       write(6,'(1x,a53)',advance='no') ' (O)BTD, (S)ingle a+ or (T)wo Particle a+a+ Transfer ' 
        read(5,*) ans
        if (ans=='o' .or. ans=='O') then
            o=.true.       
            return
        end if
        if (ans=='s' .or. ans=='S') then
            s=.true.       
            return
        end if
        if (ans=='t' .or. ans=='T') then
            t=.true.
            return
        end if
        go to 1 
    end subroutine get_operator 
    
    subroutine obtd_out()
        implicit none
        
        integer:: shl,shr
        
        do shl=1,no_shells
            do shr=1,no_shells
                write(11,'(8i4,f10.6)') 0,shl,shr,0,sj2(shl),1,sj2(shr),1,real(1,8)
            end do
        end do
    end subroutine obtd_out
    
    subroutine sp_out()
        implicit none
        integer:: sh 
        
        do sh=1,no_shells
                write(11,'(8i4,f10.6)') 0,sh,0,0,sj2(sh),1,0,0,real(1,8)
        end do
    end subroutine sp_out
            
    subroutine tp_out
        implicit none
        integer:: sh1,sh2,j
        
        do sh1=1,no_shells
            do sh2=sh1,no_shells
                do j=abs(sj2(sh1)-sj2(sh2)),sj2(sh1)+sj2(sh2),2
                    write(11,'(8i4,f10.6)') sh1,sh2,0,0,j,0,0,0,real(1,8)
                    write(11,'(8i4,f10.6)') sh1,sh2,0,0,j,2,0,0,real(1,8)
                end do
            end do
        end do
    end subroutine tp_out
                
end module Transitions    

!  NuShell.f90 
!
!  FUNCTIONS:
!  NuShell      - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: NuShell
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program NuShell
    
    use OutputNus
    use Levels
    use Transitions
    use ModelSpace
    use Extra
    use OutputNuShell

    implicit none

    ! Variables
    character(len=1):: ans
    character(len=4):: buffer
    integer:: status
    
    ! Body of NuShell
    help=.false.
    call get_cmd_arg(1,buffer,status)
    if (status>0) then
        if (buffer=='help') then
            help=.true.
        end if
    end if
    
    call output_header(6)
    write(6,*) 
    write(6,*) ' To get help messages use switch help ie run `NuShell help`.'
    if (help) then
    write(6,*)
    write(6,*) ' To run NuShell you need to decide on a model space and interaction.'
    write(6,*) ' If you use an Oxbash .int file check that for shells and their order.'
    write(6,*) ' Otherwise you will need to make a .int file from data in literature.'
    write(6,*) ' To do this use an Oxbash file as an example.'
    write(6,*) ' You can use the MSDI interaction that NuShell will calculate.'
    write(6,*) ' You need the four parameters A0,A1,B,C and single particle energies.'
    write(6,*) ' Again research literature. Search for MSDI , SDI or delta interaction.'
    write(6,*) ' Calculations with many shells and many particles get very big.'
    write(6,*) ' Try to keep you calculations as small as possible and use truncation.'
    write(6,*) ' Note 2*J,2*T must be odd for odd nuclei and even for even nuclei.'
    write(6,*) ' Current directory must contain interaction file (& model space file).'
    
    write(6,*) 
    write(6,*) ' First choose the option Model to define your model space .sps file.'
    write(6,*) ' After this or if you use an existing model space choose option Levels.'
    write(6,*) ' If you used option Levels use option Inter for a new interaction.'
    write(6,*) ' Choose Trans for gamma decay, beta decay or spectroscopic factors.'
    write(6,*) ' If two Nuclei are involved run Levels for initial and final nucleus.'
    end if
    open(unit=56,file='Shell.bat',action='write')
        
10  write(6,*)
    write(6,'(a34)', advance='no') '  (M)odel,(L)evels,(T)rans,E(x)it '
    read (5,'(a1)') ans
    write(6,*)
    
    if (ans == 'm' .or. ans == 'M') call CreateModel()
    if (ans == 'l' .or. ans == 'L') call Level()
!    if (ans == 'i' .or. ans == 'I') call Interaction()
    if (ans == 't' .or. ans == 'T') call Transition()
    if (ans == 'x' .or. ans == 'X') go to 100
    
    go to 10
    
100 write(6,*) ' Command file ','Shell.bat',' saved to disk.'
    if (help) write(6,*) ' Rename this file if you want to repeat calculations.'
    if (help) write(6,*) ' Otherwise delete it.'
    close(unit=56)

    end program NuShell

