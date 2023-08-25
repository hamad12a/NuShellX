! This file contains all the modules specific to generric program gShell
! for NuShellX and SunShell only

module OutputNus

use Partition
use Extra

implicit none

character(len=5):: Nucleus
character(len=40),dimension(3):: Description
character(len=1):: comment,exclam,xtra
character(len=7),private:: filename
character(len=1):: proton,neutron,ans
character(len=max_no_shells):: nucleon_line
character(len=11)::spe_file
integer:: no_cI,no_cN,no_cNe,no_cPr
integer(kind=1):: del_cJ2,del_cT2,del_cJ2p,del_cT2p,del_cJ2n,del_cT2n
integer:: tParity,tParityp,tParityn,isoi
integer(kind=1):: min_cJ2,max_cJ2,min_cJ2p,max_cJ2p,min_cJ2n,max_cJ2n
integer:: min_MSon,max_MSon,min_MSonp,max_MSonp,min_MSonn,max_MSonn
integer:: min_cT2,max_cT2,nmp,nmpp,nmpn
integer:: output_control
integer:: shn,no_shells,max_sj2u
real:: cut_energy,cut_energyp,cut_energyn
real(kind=8),dimension(max_no_shells)::  spe
real(kind=8):: A,A0,power,renorm,renormn,renormp,renormpn
integer:: no_levels
integer,dimension(max_no_shells):: cN,sn,sl,sj2,st2,sp
character(len=1),dimension(max_no_shells):: nucleon
type(spartition):: min_cSNo,max_cSNo,min_cSNop,max_cSNop,min_cSNon,max_cSNon
type(spartition):: min_MS_No,max_MS_No,min_MS_Nop,max_MS_Nop,min_MS_Non,max_MS_Non
logical:: help,iso=.true.,cut=.false.

data exclam /'!'/
data filename /'NuShellX'/

contains

    subroutine output_shell(inExt)
        implicit none
        character(len=6):: inExt
        
        
        open(unit=10,file=inExt//'.nus',action='write')
        
        write (10,'(a6)') Nucleus
! Intel
        write (10,'(a40)') 'NuShellX is a nuclear shellmodel program'
        write (10,'(a40)') 'written in Fortran95 for Windows/Linux. '
        write (10,'(a40)') 'It is based on NuShell by W D M Rae 2007'
! end Intel
! Sun
!        write (10,'(a40)') 'SunShell is a nuclear shellmodel program'
!        write (10,'(a40)') 'written in Fortran95 for Solaris x86/Sun'
!        write (10,'(a40)') 'It is based on NuShellX by W D M Rae 2007'
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
        if (del_cJ2/=0.and.del_cT2/=0) write (10,'(2i4)') del_cJ2,del_cT2
        write (10,*)
               
        close(unit=10)
        
    end subroutine output_shell
       
    subroutine output_shelln(inExt)
        implicit none
        character(len=5):: inExt
        
        
        open(unit=10,file=inExt//neutron//'.nus',action='write')
        
        write (10,'(a6)') Nucleus//neutron
! Intel
        write (10,'(a40)') 'NuShellX is a nuclear shellmodel program'
        write (10,'(a40)') 'written in Fortran95 for Windows/Linux. '
        write (10,'(a40)') 'It is based on NuShell by W D M Rae 2007'
! end Intel
! Sun
!        write (10,'(a40)') 'SunShell is a nuclear shellmodel program'
!        write (10,'(a40)') 'written in Fortran95 for Solaris x86/Sun'
!        write (10,'(a40)') 'It is based on NuShellX by W D M Rae 2007'
! end sun
        write (10,'(2i4)') -1,no_cNe
        write (10,'(i4)') tParityn
        write (10,'(i4)') no_shells
        write (10,'(100(a1))') nucleon(1:no_shells)
        do shn=1,no_shells
            write (10,'(6i4)') cN(shn),sn(shn),sl(shn),sj2(shn),st2(shn),sp(shn)
        end do
        write (10,'(2i4)') min_cJ2n,max_cJ2n
        write (10,'(3i4)') min_cT2,max_cT2,no_cNe
        write (10,'(2i4)') min_MSonn,max_MSonn
        write (10,'(100i4)') (min_cSNo%shell(shn),shn=1,no_shells)
        write (10,'(100i4)') (max_cSNo%shell(shn),shn=1,no_shells)
        write (10,'(100i4)') (min_MS_No%shell(shn),shn=1,no_shells)
        write (10,'(100i4)') (max_MS_No%shell(shn),shn=1,no_shells)
        write (10,'(i4)') output_control
        write (10,'(2i4)') del_cJ2n,del_cT2n
        
        if (cut) then
            write(10,'(a11)') spe_file
            write(10,'(f10.5)') cut_energyn
        end if
        write (10,*)
               
        close(unit=10)
        
    end subroutine output_shelln
    
    subroutine output_shellp(inExt)
        implicit none
        character(len=5):: inExt
        
        
        open(unit=10,file=inExt//proton//'.nus',action='write')
        
        write (10,'(a6)') Nucleus//proton
! Intel
        write (10,'(a40)') 'NuShellX is a nuclear shellmodel program'
        write (10,'(a40)') 'written in Fortran95 for Windows/Linux. '
        write (10,'(a40)') 'It is based on NuShell by W D M Rae 2007'
! end Intel
! Sun
!        write (10,'(a40)') 'SunShell is a nuclear shellmodel program'
!        write (10,'(a40)') 'written in Fortran95 for Solaris x86/Sun'
!        write (10,'(a40)') 'It is based on NuShellX by W D M Rae 2007'
! end sun
        write (10,'(2i4)') -1,no_cPr
        write (10,'(i4)') tParityp
        write (10,'(i4)') no_shells
        write (10,'(100(a1))') nucleon(1:no_shells)
        do shn=1,no_shells
            write (10,'(6i4)') cN(shn),sn(shn),sl(shn),sj2(shn),st2(shn),sp(shn)
        end do
        write (10,'(2i4)') min_cJ2p,max_cJ2p
        write (10,'(3i4)') min_cT2,max_cT2,-no_cPr
        write (10,'(2i4)') min_MSonp,max_MSonp
        write (10,'(100i4)') (min_cSNo%shell(shn),shn=1,no_shells)
        write (10,'(100i4)') (max_cSNo%shell(shn),shn=1,no_shells)
        write (10,'(100i4)') (min_MS_No%shell(shn),shn=1,no_shells)
        write (10,'(100i4)') (max_MS_No%shell(shn),shn=1,no_shells)
        write (10,'(i4)') output_control
        write (10,'(2i4)') del_cJ2p,del_cT2p
        
        if (cut) then
            write(10,'(a11)') spe_file
            write(10,'(f10.5)') cut_energyp
        end if
        write (10,*)
               
        close(unit=10)
        
    end subroutine output_shellp

    subroutine spe_file_shell(inExt,xtra,renorm)
        implicit none
        character(len=6):: inExt
        character(len=1):: xtra
        real(kind=8):: renorm
        
        if (xtra==' ') then
            open(unit=10,file=inExt//'.spe',action='write')
        else
            open(unit=10,file=inExt//xtra//'.spe',action='write')
        end if
        write (10,'(100f10.4)') spe(1:no_shells)
        write (10,'(f14.10)') renorm
        write (10,'(3i4)') no_levels,isoi,0
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
        character(len=5):: inExt
        integer:: non,nop
        iso=.true.
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
        non=0
        do shn=1,no_shells
            nucleon(shn)=nucleon_line(shn:shn)
            call read_comments_model()
            read(10,*,err=30,end=31) cN(shn),sn(shn),sl(shn),sj2(shn),st2(shn),sp(shn)
            if (nucleon(shn)=='n') then
                non=non+1
                if (shn-no_shells/2>0) then
                   if (sn(shn)/=sn(shn-no_shells/2).or.&
                        sl(shn)/=sl(shn-no_shells/2).or.&
                           sj2(shn)/=sj2(shn-no_shells/2)) iso=.false.
                end if
            end if
            if (sj2(shn) > max_sj2u) max_sj2u=sj2(shn)
        end do
        nop=no_shells-non
        call read_comments_model()
        read (10,*,err=30,end=31) min_MSonp,max_MSonp,min_MSonn,max_MSonn,min_MSon,max_MSon
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
    
    subroutine input_nucleus_shell(inExt,xtra,nucl)
        implicit none
        character(len=5),intent(out):: InExt
        character(len=1),intent(out),optional:: xtra,nucl
        integer:: status
        write(6,'(a10)',advance='no') '  Nucleus '
        read(5,'(a5)') inExt
        inExt=adjustr(inExt)
        if (present(xtra)) then
            if (present(nucl)) then
            write(6,'(a22)',advance='no') '  Neutron a or b code '
            else
            write(6,'(a19)',advance='no') '  Interaction code '
            end if
            read(5,'(a1)') xtra
            if (present(nucl)) nucl=xtra
        end if
        
    end subroutine input_nucleus_shell

    subroutine input_nucleus_final_shell(inExt,xtra,nucl)
        implicit none
        character(len=5),intent(out):: InExt
        character(len=1),intent(out),optional:: xtra,nucl
        integer:: status
        write(6,'(a10)',advance='no') '  Nucleus '
        read(5,'(a5)') inExt
        inExt=adjustr(inExt)
        if (present(xtra)) then
            if (present(nucl)) then
            write(6,'(a22)',advance='no') '  Neutron a or b code '
            else
            write(6,'(a19)',advance='no') '  Interaction code '
            end if
            read(5,'(a1)') xtra
            if (present(nucl)) nucl=xtra
        end if
        
    end subroutine input_nucleus_final_shell

    subroutine input_model_space_shell(inExt,xtra)
        implicit none
        character(len=5),intent(out):: InExt
        character(len=1),intent(out),optional:: xtra
        integer:: status
        write(6,'(a14)',advance='no') '  Model Space '
        read(5,'(a5)') inExt
        inExt=adjustr(inExt)
        if (present(xtra)) then
            write(6,'(a19)',advance='no') '  Interaction code '
            read(5,'(a1)') xtra
        end if
        
    end subroutine input_model_space_shell

    subroutine input_intfile_shell(intfile,no)
        implicit none
        character(len=6),intent(out):: intfile
        integer,intent(in),optional:: no
        integer::status
            if (.not.present(no))write(6,'(a19)',advance='no') '  Interaction file '
            if (present(no).and. no==1)write(6,'(a28)',advance='no') '  Name for transition files '
            if (present(no).and. no==2)write(6,'(a15)',advance='no') '  MBME op file '
            if (present(no).and. no==3)write(6,'(a19)',advance='no') '  MBME matrix file '
            read(5,*) intfile
            intfile=adjustr(intfile)
        
    end subroutine input_intfile_shell
    
    subroutine read_spe_shell(intfile)
        implicit none
        character(len=6),intent(in):: intfile
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
module TrdCharX

use OutputNus
implicit none


character(len=1):: typea='?',typeb='?'

contains

    
             
    function no_shellsp()
    implicit none
    integer:: no_shellsp,shn
    no_shellsp=0
    do shn=1,no_shells
    if (nucleon(shn)=='p') then
        no_shellsp=no_shellsp+1
    else
        exit
    end if
    end do
    end function no_shellsp
    
    function min_shellsa()
    implicit none
    integer:: min_shellsa
    if (typea=='p') then
        min_shellsa=1
    else
        min_shellsa=1+no_shellsp()
    end if
    end  function min_shellsa

    function min_shellsb()
    implicit none
    integer:: min_shellsb
    if (typeb=='p') then
        min_shellsb=1
    else
        min_shellsb=1+no_shellsp()
    end if
    end  function min_shellsb
    
    function max_shellsa()
    implicit none
    integer:: max_shellsa
    if (typea=='p') then
        max_shellsa=no_shellsp()
    else
        max_shellsa=no_shells
    end if
    end  function max_shellsa
    
    function max_shellsb()
    implicit none
    integer:: max_shellsb
    if (typeb=='p') then
        max_shellsb=no_shellsp()
    else
        max_shellsb=no_shells
    end if
    end  function max_shellsb
    

end module TrdCharX    
    
module Levels
use extra
use InputModelSpace
use TrdCharX
use TriDiag

implicit none
!  Intel
!  end Intel
!  Sun
!  include 'system.inc'
!  end Sun
integer:: nI,nP,nN,I,errnum,min_cJ2s
character(len=1),dimension(0:72):: ancode
data ancode /'0','1','2','3','4','5','6','7','8','9',  &
              '_','$','c','d','e','f','g','h','i','j',  &
              'k','l','m','n','o','p','q','r','s','t',  &
              'u','v','w','x','y','z','@','A','B','C',  &
              'D','E','F','G','H','I','J','K','L','M',  &
              'N','M','P','Q','R','S','T','U','V','W',  &
              'X','Y','Z','£','a','%','^','&','=','b',  &
              '-','+','#'/

contains

    function filecode(i,j,k)
    implicit none
    integer,intent(in):: i,j,k
    character(len=3):: filecode
    integer:: ii,jj,kk
    ii=i
    jj=j
    kk=k
    if (i<0) ii=i+72
    if (j<0) jj=j+72
    if (k<0) kk=k+72
    filecode=ancode(ii)//ancode(jj)//ancode(kk)
    
    end function filecode


    subroutine Level(expect)
        implicit none
        logical,intent(in):: expect
        character(len=5):: model
        character(len=7):: inffile
        character(len=6):: intfilep,intfilen,intfilepn,lanc
        character(len=1):: delet,BZ,dd,BZo
        character(len=2):: npi
        character(len=15):: chcutBZ,chcutBzm
        integer(kind=2)::res,nja,pja,del_J2
        real:: A00,A1,B,C,cutBZ,cutBZm
        logical:: core,llsp
        
        
        nI=0
        nP=0
        nN=0
        core=.false.
        
        min_cT2=0
        max_cT2=0
        output_control=0
!        write(6,'(a39)'),advance='no') ' Automatic nucleus naming wrt core y/n '
!        read(5,'(a1)') ans
!        if (ans=='y'.or.ans=='Y') then
!            write(6,*) ' Input core eg O16, Ca40, Pb208 '
!            core=.true.
!        end if
        if (expect) then
            write(6,'(a45)',advance='no') '  Name of .inf file for original calculation '
            read(5,'(a7)') inffile
            inffile=adjustr(inffile)
            open(unit=15,file=inffile//'.inf')
            read(15,'(a1,a5)') dd,Nucleus
            read(15,'(2a1)') dd,xtra
            xtra='@'
            Nucleus=adjustr(Nucleus)
        else
            call input_nucleus_shell(Nucleus,xtra)
        end if    
        if (expect) then 
        open(unit=123,file=Nucleus//xtra//'E.inf',action='write')
        else
        open(unit=123,file=Nucleus//xtra//'L.inf',action='write')
        end if
        write(123,*) Nucleus, '   Nucleus'
        write(123,*) xtra,'   Interaction code' 
        if (expect) then
        read(15,'(a1,a5)') dd,model
        else
        call input_model_space_shell(model)
        end if
        model=adjustr(model)
        write(123,*) model, '   Model Space'
        call input_model_shell(model)
        if (help) then
            write(6,*) ' If you want to override the automatic assignments for a and b'
            write(6,*) ' then you can input a negative no of protons and positive no'
            write(6,*) ' of neutrons to set protons to b. Negative neutons positive protons'
            write(6,*) ' sets neutrons to b. Both negative will chose the opposite of'
            write(6,*) ' automatic choice. But note that this may slow down calculations.'
        end if
        if (nucleon(1)=='i') then
            write(6,*) ' Requires a pn model space '
            stop
        else
            if (expect) then
            read(15,*) nP
            else
            write(6,'(a16)',advance='no') '  No of protons '
            read(5,*) nP
            end if              
            write(123,*) nP,'  No of protons '
            no_cPr=abs(nP)
            if (expect) then
            read(15,*) nN
            else
            write(6,'(a17)',advance='no') '  No of neutrons '
            read(5,*) nN
            end if
            write(123,*) nN,'  No of neutrons '
            no_cNe=abs(nN)
            no_cN=abs(nP)+abs(nN)
            nmpp=-abs(nP)
            nmpn=abs(nN)
            nja=1 !sum(sj2(no_shellsp()+1:no_shells))/(no_shells-no_shellsp())
            pja=1 !sum(sj2(1:no_shellsp()))/no_shellsp()
            if (nN*nja >= nP*pja.or.&
                (nN*nja==nP*pja.and.no_shells-no_shellsp()>no_shellsp())) then
                proton='b'
                neutron='a'
                write(6,*) ' NuLnz/TRL Neutron files labelled by a'
                write(6,*) ' NuLnz/TRL  Proton files labelled by b'
            else  
                proton='a'
                neutron='b'
                write(6,*) ' NuLnz Neutron files labelled by b'
                write(6,*) ' NuLnz  Proton files labelled by a'
            end if
            nN=abs(nN)
            nP=abs(np)
            iso=0
        end if
        call input_model_shell(model)
        if (expect) then 
        read(15,'(2a1)') dd,ans
        else
        write(6,'(a42)',advance='no') '  Use same interaction for n,p and np y/n '
        read(5,'(a1)') ans
        end if
        write(123,*) ans,'  Use same interaction for n,p and np y/n '
        if (help) write(6,*) ' Are the interactions in neutron-proton or isospin form?'
        if (expect) then 
        read(15,'(a1,a2)') dd,npi
        else
        write(6,'(a27)',advance='no') '  (np) or (i) interactions '
        read(5,'(a2)') npi
        end if
        npi=adjustr(npi)
        write(123,*) npi, '  (np) or (i) interactions '

        if (help) then
        write(6,*)
        write(6,*) ' If you use an Oxbash .int file add a line at top as follows:'
        write(6,*) ' ! Ai Ac pwr, eg for sd shell Ai=18, Ac=16 and pwr=0.3 for USD.'
        write(6,*) ' Remove spaces before ! ;shell order must agree with model space.'
        write(6,*) ' Here Ac is the mass of the core, which is the mass of the nucleus'
        write(6,*) ' that would be created if you put zero particles in the model space.'
        write(6,*) ' Ai is the mass of the core +2 the mass for a two-body calculation.'
        write(6,*) ' Finally pwr is the power of (Ac/(Ac+n))**pwr used to renormalise'
        write(6,*) ' the interaction, where n is the no of particles you input to NuShellX.'
        write(6,*) ' For no renormalisation set Ac=1.0 Ai=1.0 pwr=0.0 ,since x**0.0=1.0.'
        write(6,*)
        write(6,*) ' If you do not have a file of matrix elements use type none.'
        write(6,*)
        end if
        if (expect) read(15,*)    
        if  (isoi==1.or.ans=='y'.or.ans=='Y') then
        call input_intfile_shell(intfilep)
        else
        write(6,*) ' Proton Interaction ' 
        call input_intfile_shell(intfilep)
        end if
        intfilep=adjustr(intfilep)
        inquire(file=intfilep//'.lsp', exist=llsp)
        if (llsp) I = deletf('*.lsp','NushellX.del')
        write(123,*) intfilep,' Proton Interaction ' 
        if (intfilep(4:6)=='none' .or. intfilep(4:6)=='NONE') then
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
            write(6,*) ' NuShellX needs single particle energies for shells (MeV).'
            end if
            write(6,'(1x,a15,i4)') ' No of shells =',no_shells
            write(6,'(1x,a7)',advance='no') ' SPE`s '
            read(5,*) spe(1:no_shells) 
            renormp=1.0
            open(unit=11,file='MSDIp'//'.int',action='write')
            write(11,'(a1,1x,3f8.3)') '!',1.0,1.0,1.0
            write(11,'(a8)') '^sdip2p2'
            write(11,'(i4,50f10.6)') 0,spe(1:no_shells)
            write(11,'(11f10.6)') A00,A1,B,C,-1.0,1.0,0.0,0.0,0.0,0.0,0.0
            write(11,'(5f10.6)') 0.0,0.0,0.0,0.0,0.0
            write(11,*)
            write(11,*) ' Brussaard and Glaudemans MSDI interaction see text book'
            close(unit=11)
            intfilep=' MSDIp'
        else           
            call read_spe_shell(intfilep)
            A=A+real(nN+nP,8)
            renormp=(A0/A)**power
        end if
        if  (isoi==1.or.ans=='y'.or.ans=='Y') then
            intfilen=intfilep
            intfilepn=intfilep
            renormn=renormp
            renormpn=renormp
            go to 100
        end if
        write(6,*) ' Neutron Interaction ' 
        call input_intfile_shell(intfilen)
        intfilen=adjustr(intfilen)
        inquire(file=intfilen//'.lsp', exist=llsp)
        if (llsp) I = deletf('*.lsp','NushellX.del')
        if (expect) read(15,*)    
        write(123,*) intfilen,' Neutron Interaction ' 
        if (intfilen(4:6)=='none' .or. intfilen(4:6)=='NONE') then
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
            write(6,*) ' NuShellX needs single particle energies for shells (MeV).'
            end if
            write(6,'(1x,a15,i4)') ' No of shells =',no_shells
            write(6,'(1x,a7)',advance='no') ' SPE`s '
            read(5,*) spe(1:no_shells) 
            renormn=1.0
            open(unit=11,file='MSDIn'//'.int',action='write')
            write(11,'(a1,1x,3f8.3)') '!',1.0,1.0,1.0
            write(11,'(a8)') '^sdip2p2'
            write(11,'(i4,50f10.6)') 0,spe(1:no_shells)
            write(11,'(11f10.6)') A00,A1,B,C,-1.0,1.0,0.0,0.0,0.0,0.0,0.0
            write(11,'(5f10.6)') 0.0,0.0,0.0,0.0,0.0
            write(11,*)
            write(11,*) ' Brussaard and Glaudemans MSDI interaction see text book'
            close(unit=11)
            intfilen=' MSDIn'
        else           
            call read_spe_shell(intfilen)
            A=A+real(nP+nN,8)
            renormn=(A0/A)**power
        end if
        write(6,*) ' Proton-Neutron Interaction ' 
        call input_intfile_shell(intfilepn)
        intfilepn=adjustr(intfilepn)
        inquire(file=intfilepn//'.lsp', exist=llsp)
        if (llsp) I = deletf('*.lsp','NushellX.del')
        if (expect) read(15,*)    
        write(123,*) intfilepn,' Proton-Neutron Interaction ' 
        if (intfilepn(4:6)=='none' .or. intfilepn(4:6)=='NONE') then
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
            write(6,*) ' NuShellX needs single particle energies for shells (MeV).'
            end if
            write(6,'(1x,a15,i4)') ' No of shells =',no_shells
            write(6,'(1x,a7)',advance='no') ' SPE`s '
            read(5,*) spe(1:no_shells) 
            renormpn=1.0
            open(unit=11,file='MSDIpn'//'.int',action='write')
            write(11,'(a1,1x,3f8.3)') '!',1.0,1.0,1.0
            write(11,'(a8)') '^sdip2p2'
            write(11,'(i4,50f10.6)') 0,spe(1:no_shells)
            write(11,'(11f10.6)') A00,A1,B,C,-1.0,1.0,0.0,0.0,0.0,0.0,0.0
            write(11,'(5f10.6)') 0.0,0.0,0.0,0.0,0.0
            write(11,*)
            write(11,*) ' Brussaard and Glaudemans MSDI interaction see text book'
            close(unit=11)
            intfilepn='MSDIpn'
        else           
            call read_spe_shell(intfilepn)
            A=A+real(nP+nN,8)
            renormpn=(A0/A)**power
        end if
100     if (expect) then
        read(15,'(2a1)') dd,ans
        else
        write(6,'(a48)',advance='no') '  Do you want to apply sp energy truncation y/n '
        read(5,'(a1)') ans
        end if
        write(123,*)ans,'  Do you want to apply sp energy truncation y/n '
        
        if (ans=='y'.or.ans=='Y') then
            if (expect) then
            read(15,*) cut_energyp,cut_energyn
            else
            write(6,'(a45)',advance='no') '  Maximum energy cut protons, neutrons (MeV) '
            read(5,*) cut_energyp,cut_energyn
            end if
            write(123,*) cut_energyp,cut_energyn,'  Maximum energy cut protons, neutrons (MeV) '
            cut=.true.
        end if
        if (help) then
            write(6,*) '  This is an option based on ideas in PRC,54,1641(1996)'
            write(6,*) '  It is applied to pn interaction and p and n orbitalss must be the same.'
            write(6,*) '  It allows you to put a +/- cut around e=0.0 MeV eigenvalues in ref.'
        end if    
        BZ=' '
        BZo=' '
        if (expect) then
        read(15,'(2a1)') dd,ans
        if (ans=='y'.or.ans=='Y') then
            read(15,*)
            BZo='Z'
        end if    
        end if    
        if (.not.expect)  then
        write(6,'(a48)',advance='no') '  Do you want to use DufourZuker truncation y/n '
        read(5,'(a1)') ans
        write(123,*) ans,'  Do you want to use DufourZuker truncation y/n '
        if (ans=='y'.or.ans=='Y') then
            write(6,'(a44)',advance='no') '  Absolute energy cut e (MeV) , ev cut ampl '
            read(5,*) cutBZ,cutBZm
            if (cutBZ==0.0) cutBZ=.0000000001
            if (cutBZm==0.0) cutBZm=.0000000001
            write(123,*) cutBZ,cutBZm, '  Absolute energy cut e (MeV) , ev cut ampl '
            BZ='Z'
        end if
        else
            BZ='Z'
            cutBZ=.0000000001
            cutBZm=.0000000001
        end if    
        del_cJ2=0
        del_cT2=0
        if (expect) then
        if (BZo==' ') then
        read(15,'(2a1)') dd,ans
        if (ans=='y'.or.ans=='Y') then
            read(15,*)
        end if
        end if
        else       
        if (BZ==' ') then
        write(6,'(a61)',advance='no') '  Do you wish to truncate mutipolarity of pn interaction y/n '
        read(5,'(a1)') ans
        write(123,*) ans, '  Do you wish to truncate mutipolarity of pn interaction y/n '
        if (ans=='y'.or.ans=='Y') then
            write(6,'(a38)',advance='no') '  Maximum 2*lambda and delta 2*lambda '
            read(5,*) del_cJ2,del_cT2
            write(123,*) del_cJ2,del_cT2,'  Maximum 2*lambda and delta 2*lambda '
        end if
        if (del_cJ2<0.or.btest(del_cJ2,0)) del_cJ2=0
        if (del_cT2<2.or.btest(del_cT2,0)) del_cT2=2
        if (del_cJ2==0) del_cT2=0 
        end if       
        end if
        if (expect) then       
        read(15,*) min_cJ2,max_cJ2,del_J2
        else
        write(6,'(a37)',advance='no') '  Minumum and Maximum 2*J delta(2*J) '
        read(5,*) min_cJ2,max_cJ2,del_J2
        end  if
        write(123,*) min_cJ2,max_cJ2,del_J2,'  Minumum and Maximum 2*J delta(2*J) '
        if (del_J2<2) del_J2=2
        typea='p'
        typeb='n'
        del_cJ2p=2
        del_cJ2n=2
        call minmaxp(min_cJ2p,max_cJ2p)
        call minmaxn(min_cJ2n,max_cJ2n)
        if  (.not.expect) then
        write(6,*) ' If you have truncated the basis do not accept automatically.'
        write(6,*) ' For a mixed parity model space the a and b parity choice can'
        write(6,*) ' effectively act as a truncation so check carefully.'
        write(6,*) ' Maximum values of 2*J must be less than 255. ' 
        write(6,*) ' Minumum , Maximum , Delta 2*J p ',min_cJ2p,max_cJ2p,del_cJ2p
        write(6,'(a13)',advance='no') '  Accept y/n '
        read(5,'(a1)') ans
        else
        read(15,'(2a1)') dd,ans
        end if
        write(123,*) ans,'  Accept y/n '
        if (ans=='y'.or.ans=='Y') go to 120 
        if (.not.expect) then
        write(6,'(a32)',advance='no') '  Minumum, Maximum, Delta 2*J p '
        read(5,*) min_cJ2p,max_cJ2p,del_cJ2p
        else
        read(15,*) min_cJ2p,max_cJ2p,del_cJ2p
        end if
        write(123,*) min_cJ2p,max_cJ2p,del_cJ2p,'  Minumum, Maximum, Delta 2*J p '
        if (del_cJ2p<2.or.btest(del_cJ2p,0)) del_cJ2p=2
120     if (.not.expect) then
        write(6,*) ' Minumum , Maximum , Delta 2*J n ',min_cJ2n,max_cJ2n,del_cJ2n
        write(6,'(a13)',advance='no') '  Accept y/n '
        read(5,'(a1)') ans
        else
        read(15,'(2a1)') dd,ans
        end if
        write(123,*) ans,'  Accept y/n '
        if (ans=='y'.or.ans=='Y') go to 140 
        if (.not.expect) then
        write(6,'(a32)',advance='no') '  Minumum, Maximum, Delta 2*J n '
        read(5,*) min_cJ2n,max_cJ2n,del_cJ2n
        else
        read(15,*) min_cJ2n,max_cJ2n,del_cJ2n
        end if
        write(123,*) min_cJ2n,max_cJ2n,del_cJ2n,'  Minumum, Maximum, Delta 2*J n '
        if (del_cJ2n<2.or.btest(del_cJ2n,0)) del_cJ2n=2
140     del_cT2n=2
        del_cT2p=2
        nmp=nN-nP
        min_cT2=0
        max_cT2=0
        if (.not.expect) then
        if (help) write(6,*) ' Parity of 0 means both + and - parities'
        write(6,'(a19)',advance='no') '  Parity'//' '//proton//' +1,0,-1 '
        read(5,*) tParityp
        else
        read(15,*) tParityp
        end if
        write(123,*) tParityp,'  Parity'//' '//proton//' +1,0,-1 '
        if (.not.expect) then
        write(6,'(a19)',advance='no') '  Parity'//' '//neutron//' +1,0,-1 '
        read(5,*) tParityn
        else
        read(15,*) tParityn
        end if
        write(123,*) tParityn,'  Parity'//' '//neutron//' +1,0,-1 '
        if (tParityn==0.or.tParityp==0) then
            if (.not.expect) then
150         write(6,'(a34)',advance='no') '  Parity of Total Nucleus +1,0,-1 '
            read(5,*) tParity
            else
            read(15,*) tParity
            end if
            write(123,*) tParity, '  Parity of Total Nucleus +1,0,-1 '
!            if (tParity==0) then 
!                print *, ' Zero not allowed +1,-1 only '
!                go to 150
!            end if
        else
            tParity=tParityn*tParityp
        end if
        if (.not.expect) then
        write(6,'(a15)',advance='no') '  No of levels '
        read(5,*) no_levels
        else
        read(15,*) no_levels
        end if
        write(123,*) no_levels,'  No of levels '
        if (help.and..not.expect) then
            write(6,*)
            write(6,*) ' It is advised that you delete existing .lev files.'
            write(6,*) ' If you want to save them copy them to another directory.'
            write(6,*)
        end if
        if (.not.expect) then
        write(6,'(a33)',advance='no') '  Delete existing .lev files y/n '
        read(5,*) delet
        write(123,*) delet,'  Delete existing .lev files y/n '
        if (delet == 'y') then
            I = deletf('*.lev','NushellX.del')
        end if
        end if
        min_cJ2s=min_cJ2
        if (core) then
            no_cN=-no_cN
            no_cPr=-no_cPr
            no_cNe=-no_cNe
        end if
        spe_file=Nucleus//xtra//'.spe'
        call output_shell(' '//Nucleus)
        do min_cJ2=min_cJ2s,max_cJ2,del_J2
            call output_shell(Nucleus//ancode(min_cJ2/2))
        end do
        call output_shellp(Nucleus//proton)
        call output_shelln(Nucleus//neutron)
        write(6,*) ' Data saved as ',Nucleus//'.nus'
        write(6,*) ' Data saved as ',Nucleus//proton//'.nus'
        write(6,*) ' Data saved as ',Nucleus//neutron//'.nus'
        call spe_file_shell(Nucleus//xtra,' ',renormpn)
        call spe_file_shell(Nucleus//proton,xtra,renormp)
        call spe_file_shell(Nucleus//neutron,xtra,renormn)
        write(6,*) ' SPE file created ',Nucleus//xtra//'.spe'
        write(6,*) ' SPE file created ',Nucleus//proton//xtra//'.spe'
        write(6,*) ' SPE file created ',Nucleus//neutron//xtra//'.spe'
        if (help) then
            write(6,*) ' NuShellX has to run NuShell and other programs to'
            write(6,*) ' setup transition densities, matrices and tables'
            write(6,*) ' which are common to all caculations for different'
            write(6,*) ' J values of the Nucleus you are calculating.'
            write(6,*) ' If you have already done calculations for the same'
            write(6,*) ' Nucleus with identical model spaces, interactions'
            write(6,*) ' proton and neutron numbers etc etc then answer '
            write(6,*) ' no to the next question.'
        end if
        if (BZ==' ') then
        write(6,'(a59)',advance='no') '  Do you want to use (s)tandard or (t)hick restart Lanczos '
        read(5,'(a1)') ans
        write(123,*) ans,'  Do you want to use (s)tandard or (t)hick restart Lanczos '
        lanc='NuTRL '
        if (ans=='s'.or.ans=='S') lanc='NuLnz '
        else 
        if (expect) then
        lanc='NuXpct'
        else
        lanc='NuTRLZ'
        end if
        end if
        write(6,'(a40)',advance='no') '  Do you need to setup basic tables y/n '
        read(5,'(a1)') ans
        write(123,*) ans,'  Do you need to setup basic tables y/n '
        close(123)
        close(15)
        if (ans=='n'.or.ans=='N') go to 500
        res=run('NuOp',Nucleus//' '//xtra//' '//intfilepn//' '//npi)
        if (res==-1) stop
        write(56,*)'NuOp',' ',Nucleus//' '//xtra//' '//intfilepn//' '//npi
        if (BZ=='Z') then
        cutBZ=abs(cutBZ)
        cutBZm=abs(cutBZm)
        write(chcutBZ,'(f15.11)') cutBZ
        write(chcutBZm,'(f15.11)') cutBZm
        res=run('NuDiop',Nucleus//xtra//' '//intfilep//' '//chcutBZ//' '//chcutBZm)
        if (res==-1) stop
        write(56,*)'NuDiop'//' '//Nucleus//xtra//' '//intfilep//' '//chcutBZ//' '//chcutBZm
        write(6,*) ' NuDiop ran successfully'
        write(6,*)
        end if
        res=run('NuBasis',Nucleus//proton)
        write(56,*) 'NuBasis',' ',Nucleus//proton
        if (res == -1) stop
        res=run('NuProj',Nucleus//proton)
        write(56,*) 'NuProj',' ',Nucleus//proton
        if (res == -1) stop
        res=run('NuOper',Nucleus//proton//' '//xtra//' '//intfilep//' '//npi)
        write(56,*)'NuOper',' ',Nucleus//proton//' '//xtra//' '//intfilep//' '//npi
        if (res == -1) stop
        res=run('NuMatrix',Nucleus//proton//' '//xtra)
        write(56,*) 'NuMatrix',' ',Nucleus//proton//' '//xtra
        if (res == -1) stop
        res=run('NuOrth',Nucleus//proton//' '//xtra)  
        write(56,*) 'NuOrth',' ',Nucleus//proton//' '//xtra
        if (res == -1) stop
        res=run('NuOpm',Nucleus//proton//' '//Nucleus//proton//' a+a- '//proton)
        write(56,*) 'NuOpm',' ',Nucleus//proton//' '//Nucleus//proton//' a+a- '//proton
        if (res == -1) stop
        if (nP/=0.and.nN/=0) then
        res=run('NuOpd'//BZ,Nucleus//proton//' '//Nucleus//proton//' a+a- '//proton &
                    //' '//intfilep)
        write(56,*) 'NuOpd'//BZ,' ',Nucleus//proton//' '//Nucleus//proton//' a+a- '//proton &
                        //' '//intfilep
        if (res == -1) stop
        end if
        if (nP==0.or.nN==0) then
        res=run('NuOpd'//BZ,Nucleus//proton//' '//Nucleus//proton//' o '//proton &
                    //' '//intfilep)
        write(56,*) 'NuOpd'//BZ,' ',Nucleus//proton//' '//Nucleus//proton//' o '//proton &
                        //' '//intfilep
        if (res == -1) stop
        end if
        res=run('NuBasis',Nucleus//neutron)
        write(56,*) 'NuBasis',' ',Nucleus//neutron
        if (res == -1) stop
        res=run('NuProj',Nucleus//neutron)
        write(56,*) 'NuProj',' ',Nucleus//neutron
        if (res == -1) stop
        res=run('NuOper',Nucleus//neutron//' '//xtra//' '//intfilen//' '//npi)
        write(56,*)'NuOper',' ',Nucleus//neutron//' '//xtra//' '//intfilen//' '//npi
        if (res == -1) stop
        res=run('NuMatrix',Nucleus//neutron//' '//xtra)
        write(56,*) 'NuMatrix',' ',Nucleus//neutron//' '//xtra
        if (res == -1) stop
        res=run('NuOrth',Nucleus//neutron//' '//xtra)  
        write(56,*) 'NuOrth',' ',Nucleus//neutron//' '//xtra
        if (res == -1) stop
        res=run('NuOpm',Nucleus//neutron//' '//Nucleus//neutron//' a+a- '//neutron)
        write(56,*) 'NuOpm',' ',Nucleus//neutron//' '//Nucleus//neutron//' a+a- '//neutron
        if (res == -1) stop
        if (nN/=0.and.nP/=0) then
        res=run('NuOpd'//BZ,Nucleus//neutron//' '//Nucleus//neutron//' a+a- '//neutron &
                    //' '//intfilep)
        write(56,*) 'NuOpd'//BZ,' ',Nucleus//neutron//' '//Nucleus//neutron//' a+a- '//neutron &
                    //' '//intfilep
        if (res == -1) stop
        end if
        if (nN==0.or.nP==0) then
        res=run('NuOpd'//BZ,Nucleus//neutron//' '//Nucleus//neutron//' o '//neutron &
                    //' '//intfilep)
        write(56,*) 'NuOpd'//BZ,' ',Nucleus//neutron//' '//Nucleus//neutron//' o '//neutron &
                    //' '//intfilep
        if (res == -1) stop
        end if
500     do min_cJ2=min_cJ2s,max_cJ2,del_J2
        res=run(lanc,Nucleus//ancode(min_cJ2/2)//' '//Nucleus//'a '//Nucleus//'b '//Nucleus//xtra//&
                ' '//intfilep)
        write(56,*) lanc,' ',Nucleus//ancode(min_cJ2/2)//' '//Nucleus//'a '//Nucleus//'b '//Nucleus//xtra//&
                ' '//intfilep
        if (res == -1) stop
        if (.not.expect) then
        res=run('NuVec',Nucleus//ancode(min_cJ2/2)//' '//Nucleus//'a '//Nucleus//'b '//Nucleus//xtra)
        write(56,*) 'NuVec ',Nucleus//ancode(min_cJ2/2)//' '//Nucleus//'a '//Nucleus//'b '//Nucleus//xtra
        if (res == -1) stop
        end if
        end do
        if (.not.expect) then
        res=run('NuPSLevel',' ')
        write(56,*) 'NuPSLevel ',' '
        if (res == -1) stop
        end if
        write(6,*)
        if (expect) then
        write(6,*) ' Your input has been saved to '//Nucleus//xtra//'E.inf'//' Keep this file.'
        else
        write(6,*) ' Your input has been saved to '//Nucleus//xtra//'L.inf'//' Keep this file.'
        end if
        return
        
    end subroutine Level
    
    subroutine minmaxp(jmin,jmax)
        implicit none
        integer(kind=1),intent(out):: jmin,jmax
        real,dimension(max_no_shells):: rsj2
        integer,dimension(max_no_shells):: indx,sj2s
        integer:: i,c,n,s,d
        jmin=0
        if (btest(no_cPr,0)) jmin=1
        c=0
        do i = min_shellsa(),max_shellsa()
            c=c+1
            rsj2(c)=-sj2(i)
        end do
        call indexx(c,rsj2,indx)
        sj2s(1:c)=abs(rsj2(indx(1:c)))
        
        jmax=0
        d=0
        s=sj2s(1)
        do n=1,no_cPr
            d=d+1
            if (d>c) then
                d=1
                s=sj2s(d)
                jmax=jmax+s
                sj2s(d)=sj2s(d)-2
            else if (s==sj2s(d)) then
                jmax=jmax+s
                sj2s(d)=sj2s(d)-2
            else
                d=1
                s=sj2s(d)
                jmax=jmax+s
                sj2s(d)=sj2s(d)-2
            end if
        end do
        if (jmax<jmin) stop ' maxj < minj p '
        
    end subroutine minmaxp
        
    subroutine minmaxn(jmin,jmax)
        implicit none
        integer(kind=1),intent(out):: jmin,jmax
        real,dimension(max_no_shells):: rsj2
        integer,dimension(max_no_shells):: indx,sj2s
        integer:: i,c,n,s,d
        jmin=0
        if (btest(no_cNe,0)) jmin=1
        c=0
        do i = min_shellsb(),max_shellsb()
            c=c+1
            rsj2(c)=-sj2(i)
        end do
        call indexx(c,rsj2,indx)
        sj2s(1:c)=abs(rsj2(indx(1:c)))
        
        jmax=0
        d=0
        s=sj2s(1)
        do n=1,no_cNe
            d=d+1
            if (d>c) then
                d=1
                s=sj2s(d)
                jmax=jmax+s
                sj2s(d)=sj2s(d)-2
            else if (s==sj2s(d)) then
                jmax=jmax+s
                sj2s(d)=sj2s(d)-2
            else
                d=1
                s=sj2s(d)
                jmax=jmax+s
                sj2s(d)=sj2s(d)-2
            end if
        end do
        if (jmax<jmin) stop ' maxj < minj p '
        
    end subroutine minmaxn
        
end module Levels    

module ModelSpace

use InputModelSpace

implicit none

contains

    subroutine CreateModel()
        implicit none
        logical:: use_isospin
        integer:: no_shells,cN,sn,sl,sj2,st2,sp,no_proton,no_neutron
        integer:: shn,min,max,min1,max1,min2,max2,sfl,nmsp,nmsn,shn2
        character(len=1):: ans,ai='i',an='n',ap='p',chl
        character(len=6):: chsh
        character(len=5):: chmodel
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
        write(6,*) ' five characters long. Interaction and transition filenames'
        write(6,*) ' must be six characters long.'
        end if        
        call input_model_space_shell(chmodel)
        chmodel=adjustr(chmodel)
        
        open(unit=11,file=chmodel//'.sps',action='write')
        use_isospin=.false.
10      write(6,'(1x,a20)',advance='no') ' Total no of shells '
        read(5,*) no_shells
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
        if (help) then
        write(6,*)
        write(6,*) ' If you use an Oxbash .int file get shell order from that file.'
        write(6,*) ' NuShellX uses the convention that n starts at 0.'
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
            st2=0
            cN=2*sn+sl
            cNa(shn)=cN
            write(11,'(6i4)') cN,sn,sl,sj2,st2,sp
        end do 
        write(6,'(1x,a41)',advance='no') ' Do you want to truncate model space y/n '
        read(5,'(a1)') ans
        if (ans=='y') then
            nmsp=0
            nmsn=0
            do shn=1,no_proton
                min=minval(cNa(1:no_proton))
                if (min==999) exit
                nmsp=nmsp+1
                cNs(nmsp)=min
                    do  shn2=1,no_proton
                         max=maxval(minloc(cNa(1:no_proton)))
                         if (cNa(max)>min) exit
                         cNa(max)=999
                    end do
            end do
            do shn=1+no_proton,no_shells
                min=minval(cNa(1+no_proton:no_shells))
                if (min==999) exit
                nmsn=nmsn+1
                cNs(nmsn+no_proton)=min
                    do  shn2=1+no_proton,no_shells
                         max=maxval(minloc(cNa(1+no_proton:no_shells)))
                         if (cNa(max+no_proton)>min) exit
                         cNa(max+no_proton)=999
                    end do
            end do
            write(6,'(1x,a47)',advance='no') ' Do you want to limit the sum over N=2*n+l y/n '
            read(5,'(a1)') ans
                if (ans=='y') then
                    write(6,'(1x,a46)',advance='no') ' Proton Min Max, Neutron Min Max, P+N Min Max '
                    read(5,*) min,max,min1,max1,min2,max2
                    write(11,'(6i4)') min,max,min1,max1,min2,max2
                else
                    min=0
                    max=999
                    min1=0
                    max1=999
                    min2=0
                    max2=999
                    write(11,'(6i4)') min,max,min1,max1,min2,max2
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
                write(6,'(1x,a10,i4,a21)') ' There are',nmsp,' proton major shells.'
                write(6,'(1x,a4,50i5)') ' N= ',cNs(1:nmsp)
60              write(6,'(1x,a18)',advance='no') ' MShellNo Min Max '
                read(5,*) shn,minna(shn),maxna(shn)
                write(6,'(1x,a10)',advance='no') ' More y/n '             
                read(5,'(a1)') ans
                if (ans=='y') go to 60
                write(6,'(1x,a10,i4,a22)') ' There are',nmsn,' neutron major shells.'
                write(6,'(1x,a4,50i5)') ' N= ',cNs(no_proton+1:no_proton+nmsn)
61              write(6,'(1x,a18)',advance='no') ' MShellNo Min Max '
                read(5,*) shn,minna(shn+no_proton),maxna(shn+no_proton)
                write(6,'(1x,a10)',advance='no') ' More y/n '             
                read(5,'(a1)') ans
                if (ans=='y') go to 61
            end if
            write(11,'(50i5)') mina(1:no_shells)
            write(11,'(50i5)') maxa(1:no_shells)
            write(11,'(50i5)') minna(1:no_shells)
            write(11,'(50i5)') maxna(1:no_shells)
        else
            min=0
            max=999
            write(11,'(6i4)') min,max,min,max,min,max
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
        character(len=6):: model
        character(len=5):: Nucleusi,Nucleusf,modeli,modelf
        character(len=1):: xtrai,xtraf,delet,ans,ytrai,ytraf,dum
        character(len=6):: Nucleusoa1,Nucleusob1
        character(len=6):: mefile
        character(len=6):: Nucleusti,Nucleustf
        character(len=6):: Nucleusoi,Nucleusof,Nucleusoa0,Nucleusob2,Nucleusob0,Nucleusoa2
        integer(kind=2):: res,loop,lp,clus
        integer:: delNa,delNb,del_J2f,del_J2i
        integer:: J2in,J2ix,J2fn,J2fx,J2i,J2f,Nai,Naf,Nbi,Nbf,J2n,J2x,ni,nf
        character(len=4):: opa,opb,opc
        character(len=1):: opfa,opfb,xtra1,xtra2
        integer:: mnDJ,mxDJ,mnDT,mxDT,mnDM,mxDM,mnDTz,mxDTz,i,j,amnDM,amnDTz
        logical:: OBTD,SP,TP
        real:: ep,en,glp,gln,gsp,gsn,mass,charge,radial
                
        write(6,*) ' Initial State'
        if (help) then
        write(6,*) ' Here the code should be replaced by a or b code for'
        write(6,*) ' the neutrons in the initial calculation.'
        end if
        call input_nucleus_shell(Nucleusi,xtrai,dum)
        Nucleusi=adjustr(Nucleusi)
        if (xtrai=='a') then
            ytrai='b'
        else
            ytrai='a'
        end if
        write(6,'(a19)',advance='no') '  Interaction code '
        read(5,'(a1)') xtra1
        if (help) write(6,*) ' Now input 2*J min,max and the number of particles type a and b'
        write(6,'(a31)',advance='no') '  2*Jmin, 2*Jmax, 2*dJ, Na, Nb '
        read(5,*) J2in,J2ix,del_J2i,Nai,Nbi
        if (del_J2i<2) del_J2i=2
        write(6,*) ' Final State'
        if (help) then
        write(6,*) ' Here the code should be replaced by a or b code for'
        write(6,*) ' the neutrons in the final calculation.'
        end if
        call input_nucleus_final_shell(Nucleusf,xtraf,dum)
        Nucleusf=adjustr(Nucleusf)
        if (xtraf=='a') then
            ytraf='b'
        else
            ytraf='a'
        end if
        write(6,'(a19)',advance='no') '  Interaction code '
        read(5,'(a1)') xtra2
        if (help) write(6,*) ' Now input 2*J min,max and the number of particles type a and b'
        write(6,'(a31)',advance='no') '  2*Jmin, 2*Jmax, 2*dJ, Na, Nb '
        read(5,*) J2fn,J2fx,del_J2f,Naf,Nbf
        if (del_J2f<2) del_J2f=2
        delNa=Naf-Nai
        delNb=Nbf-Nbi
        if (.not.ytraf==ytrai) then
        delNa=Nbf-Nai
        delNb=Naf-Nbi
        end if
        open(unit=123,file=Nucleusi//xtra1//'L.inf')
        open(unit=124,file=Nucleusf//xtra2//'L.inf')
        read(123,*)
        read(123,*)
        read(123,'(1x,a5)') modeli
        read(124,*)
        read(124,*)
        read(124,'(1x,a5)') modelf
        if (.not.modeli==modelf) then
            write(6,*) ' Model spaces not identical ',modeli,' ',modelf
            stop
        end if
        if (help) then
        write(6,*) ' The transition filename is of your own choice a5'
        end if
        call input_intfile_shell(mefile,1)
        mefile=adjustr(mefile)
        if (help) then
        write(6,*) ' NushellX calculates transitions as a product of two matrices'
        write(6,*) ' The choice of matrices for both `a` particles and `b` particles are:'
        write(6,*) ' I , N, a+a- , a+  , a- , a+a+  , a-a- .'
        write(6,*) ' For (OBTD) gamma decay chose 1x2 or 2x1, for beta decay chose 3x4 or 4x3.'
        write(6,*) ' For (SP) single particle transfer chose 1x3 or 3x1 or 1x4 or 4x1.'
        write(6,*) ' For (TP) two particle transfer chose 1x5 or 5x1 or 1x6 or 6x1.'
        write(6,*) ' For (TBTD) two body transition densities chose 5x6 (pp) or 2x2 (ph).'
        write(6,*) ' For alpha spectrocopic factors chose 5x5 or 6x6.'
        write(6,*) ' For three particle transfer amplitudes chose 3x5 or 4x6' 
        write(6,*) ' 2*J must be even for OBTD, TBTD, alpha, and TP, odd for SP and 3 particle.'
        end if
        write(6,'(1x,a51)',advance='no') ' Operator acting on a particles in initial Nucleus '
        if (delNa==0) then 
            write(6,*)
            write(6,'(1x,a23)',advance='no') ' Input I , N or a+a- : ' 
            read(5,'(a4)') opa
        end if
        if (delNa==1) then
            write(6,*) ' a+' 
            opa='  a+'
        end if
        if (delNa==2) then
            write(6,*) ' a+a+' 
            opa='a+a+'
        end if
        if (delNa==-1.and.delNb==1) then
        write(6,*) ' If you have not changed the automatic a and b assignments then'
        write(6,*) ' this calculation will not produce a transition density. To get a'
        write(6,*) ' proper transition densiy swop initial an final states.'
        if (help) write(6,*) ' If you have changed defaults for a and b the + operator'
        if (help) write(6,*) ' must act on the appropriate particles'
        end if
        if (delNa==-1) then
            write(6,*) ' a-' 
            opa='  a-'
        end if
        if (delNa==-2) then
            write(6,*) ' a-a-' 
            opa='a-a-'
        end if
        opa=adjustr(opa)
        write(6,'(1x,a51)',advance='no') ' Operator acting on b particles in initial Nucleus '
        if (delNb==0) then 
            write(6,*)
            write(6,'(1x,a23)',advance='no') ' Input I , N or a+a- : ' 
            read(5,'(a4)') opb
        end if
        if (delNb==1) then
            write(6,*) ' a+' 
            opb='  a+'
        end if
        if (delNb==2) then
            write(6,*) ' a+a+' 
            opb='a+a+'
        end if
        if (delNb==-1) then
            write(6,*) ' a-' 
            opb='  a-'
        end if
        if (delNb==-2) then
            write(6,*) ' a-a-' 
            opb='a-a-'
        end if
        opb=adjustr(opb)
        loop=1
        clus=1
        if ((opa=='   I'.or.opa=='   i').and.opb=='a+a-') loop=2
        if ((opb=='   I'.or.opb=='   i').and.opa=='a+a-') loop=2
        if ((opa=='  a+'.or.opa=='a+a+').and.opb=='a+a+') clus=2
        if ((opb=='  a-'.or.opb=='  a-').and.opa=='a-a-') clus=2
        if ((opb=='  a+'.or.opb=='a+a+').and.opa=='a+a+') clus=2
        if ((opa=='  a-'.or.opa=='  a-').and.opb=='a-a-') clus=2
        if (loop==2) then
            mefile(1:5)=mefile(2:6)
            mefile(6:6)='x'
        end if
        do lp=1,loop
        opa=adjustr(opa)
        opb=adjustr(opb)
        if (lp==2) then
            opc=opa
            opa=opb
            opb=opc
        end if
        Nucleusoa0=adjustr(Nucleusi)//xtrai
        Nucleusob0=adjustr(Nucleusi)//ytrai
        Nucleusoa1=adjustr(Nucleusf)//xtraf
        Nucleusob1=adjustr(Nucleusf)//ytraf
        if (opa=='   I'.or.opa=='   i') then
            Nucleusoa2=adjustr(Nucleusi)//'i'
            opfa='i'
        else if (opa=='   N'.or.opa=='   n') then
            Nucleusoa2=adjustr(Nucleusi)//'x'
            opfa='x'
        else if (opa=='a+a-') then
            Nucleusoa2=adjustr(Nucleusi)//'c'
            opfa='c'
        else if (opa=='a-a-') then
            Nucleusoa2=adjustr(Nucleusi)//'s'
            opfa='s'
        else if (opa=='a+a+') then
            Nucleusoa2=adjustr(Nucleusi)//'d'
            opfa='d'
        else if (opa=='  a-') then
            Nucleusoa2=adjustr(Nucleusi)//'m'
            opfa='m'
        else if (opa=='  a+') then
            Nucleusoa2=adjustr(Nucleusi)//'p'
            opfa='p'
        else
            print *, ' Operator a not recognised ',opa
            stop
        end if    
        if (opb=='   I'.or.opb=='   i') then
            Nucleusob2=adjustr(Nucleusi)//'j'
            opfb='j'
        else if (opb=='   N'.or.opb=='   n') then
            Nucleusob2=adjustr(Nucleusi)//'y'
            opfb='y'
        else if (opb=='a+a-') then
            Nucleusob2=adjustr(Nucleusi)//'d'
            opfb='d'
        else if (opb=='a-a-') then
            Nucleusob2=adjustr(Nucleusi)//'t'
            opfb='t'
        else if (opb=='a+a+') then
            Nucleusob2=adjustr(Nucleusi)//'e'
            opfb='e'
        else if (opb=='  a-') then
            Nucleusob2=adjustr(Nucleusi)//'n'
            opfb='n'
        else if (opb=='  a+') then
            Nucleusob2=adjustr(Nucleusi)//'q'
            opfb='q'
        else
            print *, ' Operator b not recognised ',opb
            stop
        end if 
        en=0.0
        ep=0.0
        gln=0.0
        glp=0.0
        gsn=0.0
        gsp=0.0
        mass=0.0
        charge=0.0
        if (lp==2) go to 100
        if (loop==2) then 
            if (help) then
            write(6,*)
            write(6,*)' For gamma decay NuShell requires `effective charges and g-facctors`.'
            write(6,*)' These help correct for effects such as `core polarisation`.'
            write(6,*)' See textbooks or literature for more information.'
            write(6,*)' Typical input A == mass of nucleus Z == charge of initial nucleus.'
            write(6,*)' ep=1.5, en=0.5, glp=1.0, gln=0.0, gsp=5.586, gsn=-3.826'
            write(6,*)' Enter 0 0 0 0 0 0 for default. If you want en=0.0 enter 999.0.'
            write(6,*)
            end if
            write(6,'(1x,a29)',advance='no') ' Enter ep en glp gln gsp gsn '
            read(5,*) ep,en,glp,gln,gsp,gsn
            if (ep==0.0) ep=1.5
            if (en==0.0) en=0.5
            if (en==999.0) en=0.0
            if (glp==0.0) glp=1.0
            if (gsp==0.0) gsp=5.586
            if (gsn==0.0) gsn=-3.826
            write(6,'(a32)',advance='no') ' Mass(A), Charge(Z) of nucleus '
            read(5,*) mass,charge
        end if
        if (clus==2) then
             write(6,'(a32)',advance='no') ' Mass Initial and Final nuclei '
             read(5,*) mass,charge
        end  if 
        write(6,*) ' The operators are coupled to Jop.'
        call get_JT_values(J2n,J2x,'Operator')
        write(6,'(1x,a24)',advance='no') ' Number of states ni nf '
        read(5,*) ni,nf
        if (ni==0) ni=1
        if (nf==0) nf=1
        if (help) then
            write(6,*)
            write(6,*) ' Matrix elements depend on the sign convention for radial wavefunctions.'
            write(6,*) ' NuShell uses by default positive at origin.'
            write(6,*) ' To change this to positive  at infinity reply y to next question.'
            write(6,*)
        end if
        write(6,'(1x,a43)',advance='no') ' Change wavefunctions to + at infinity y/n '
        read(5,'(a1)') ans
        radial=1.0
        
        if (ans=='y' .or. ans=='Y') radial=-1.0
        open(unit=11,file=mefile//'.me',action='write')
        write(11,'(a1,1x,a6,1x,a6,1x,a7)') '!',Nucleusi,Nucleusf,mefile
        write(11,'(6f8.4,2f8.2)') ep,en,glp,gln,gsp,gsn,mass,charge
        if ((opa=='  a+'.and.opb=='  a-').or.&
            (opb=='  a+'.and.opa=='  a-')) then
            write(11,'(10i4,f10.6)') J2n,J2x,2,2,0,0,0,0,ni,nf,radial
        else
            write(11,'(10i4,f10.6)') J2n,J2x,0,0,0,0,0,0,ni,nf,radial
        end if    
        close(11)
        write(6,*) ' Data saved to ',mefile//'.me'
        if (loop==2) then
        mefile(6:6)='y'
        open(unit=11,file=mefile//'.me',action='write')
        write(11,'(a1,1x,a6,1x,a6,1x,a7)') '!',Nucleusi,Nucleusf,mefile
        write(11,'(6f8.4,2f8.2)') ep,en,glp,gln,gsp,gsn,mass,charge
        write(11,'(10i4,f10.6)') J2n,J2x,0,0,0,0,0,0,ni,nf,radial
        close(11)
        mefile(6:6)='x'
        end if
100     continue 
        res=run('NuOpm',Nucleusoa0//' '//Nucleusoa1//' '//opa//' '//opfa)
        if (res == -1 )stop
        write(56,*) 'NuOpm'//' '//Nucleusoa0//' '//Nucleusoa1//' '//opa//' '//opfa
        res=run('NuOpd',Nucleusoa0//' '//Nucleusoa1//' '//opa//' '//opfa)
        if (res == -1 )stop
        write(56,*) 'NuOpd'//' '//Nucleusoa0//' '//Nucleusoa1//' '//opa//' '//opfa
        
        res=run('NuOpm',Nucleusob0//' '//Nucleusob1//' '//opb//' '//opfb)
        if (res == -1 )stop
        write(56,*) 'NuOpm'//' '//Nucleusob0//' '//Nucleusob1//' '//opb//' '//opfb
        res=run('NuOpd',Nucleusob0//' '//Nucleusob1//' '//opb//' '//opfb)
        if (res == -1 )stop
        write(56,*) 'NuOpd'//' '//Nucleusob0//' '//Nucleusob1//' '//opb//' '//opfb
        
        if (lp==1) then
        do J2i=J2in,J2ix,del_J2i
        Nucleusoi=adjustr(Nucleusi)//ancode(J2i/2)
        do J2f=J2fn,J2fx,del_J2f
        Nucleusof=adjustr(Nucleusf)//ancode(J2f/2)
        res=run('NuTra',&
            Nucleusof//' '//Nucleusoa2//' '//Nucleusob2//' '//mefile//' '//Nucleusoi)
        if (res == -1 )stop
        write(56,*) &
        'NuTra'//' '//Nucleusof//' '//Nucleusoa2//' '//Nucleusob2//' '//mefile//' '//Nucleusoi
        if (loop==1.and.clus==1) then
        res=run('NuTrx',Nucleusoi//' '//Nucleusof//' '//mefile//filecode(J2i/2,Nai,Nbi)&
                //filecode(J2f/2,Naf,Nbf))
        if (res == -1 )stop
        write(56,*)'NuTrx'//' '//Nucleusoi//' '//Nucleusof//' '//mefile//filecode(J2i/2,Nai,Nbi)&
                //filecode(J2f/2,Naf,Nbf)
        else if (loop==1.and.clus==2) then
        res=run('NuClx',Nucleusoi//' '//Nucleusof//' '//mefile//filecode(J2i/2,Nai,Nbi)&
                //filecode(J2f/2,Naf,Nbf))
        if (res == -1 )stop
        write(56,*)'NuClx'//' '//Nucleusoi//' '//Nucleusof//' '//mefile//filecode(J2i/2,Nai,Nbi)&
                //filecode(J2f/2,Naf,Nbf)
        end if     
        end do
        end do
        else
        do J2i=J2in,J2ix,del_J2i
        Nucleusoi=adjustr(Nucleusi)//ancode(J2i/2)
        do J2f=J2fn,J2fx,del_J2f
        Nucleusof=adjustr(Nucleusf)//ancode(J2f/2)
        res=run('NuTra',&
            Nucleusof//' '//Nucleusoa2//' '//Nucleusob2//' '//mefile(1:5)//'y '//Nucleusoi)
        if (res == -1 )stop
        write(56,*) &
        'NuTra'//' '//Nucleusof//' '//Nucleusoa2//' '//Nucleusob2//' '//mefile(1:5)//'y '//Nucleusoi     
        if (loop==2) then
        res=run('NuTrx',Nucleusoi//' '//Nucleusof//' '//mefile(1:5)//'x'//filecode(J2i/2,Nai,Nbi)&
                //filecode(J2f/2,Naf,Nbf)//' '//mefile(1:5)//'y'//filecode(J2i/2,Nai,Nbi)&
                //filecode(J2f/2,Naf,Nbf))
        if (res == -1 )stop
        write(56,*) &
        'NuTrx'//' '//Nucleusoi//' '//Nucleusof//' '//mefile(1:5)//'x'//filecode(J2i/2,Nai,Nbi)&
                //filecode(J2f/2,Naf,Nbf)//' '//mefile(1:5)//'y'//filecode(J2i/2,Nai,Nbi)&
                //filecode(J2f/2,Naf,Nbf)
        end if
        end do
        end do
        end if
        end do
        return
        if (help) then
        write(6,*) ' The output files are labelled `mefile//(2*Ji)/2//Nai//Nbi//(2*Jf)/2//Naf//nbf` .'
        write(6,*) ' The test files have extension .tra and the binary files .trb .'
        write(6,*) ' You can view and print the text files with Notepad or any other '
        write(6,*) ' text editor. If you have chosen transition densities then two files'
        write(6,*) ' will have been produced for protons and neutrons with x and y in'
        write(6,*) ' the names.'
        write(6,*)
        write(6,*) ' For gamma dacay BE(L)`s and and M(L)`s run NuTrx at the command line as'
        write(6,*) ' C:\>NuTrx Nucleusi Nucleusf Trax.. Tray.. where Trax.. and Tray..' 
        write(6,*) ' are the names of the two .tra files generated by this program.'
        write(6,*) ' For one and two particle transfer and beta decay type'
        write(6,*) ' C:\>NuTrx Nucleusi Nucleusf Tra.. and for cluster amplitudes type'
        write(6,*) ' C:\>NuClx Nucleusi Nucleusf Tra.. In all cases run NuTrx and NuClx'
        write(6,*) ' for each pair of Ji and Jf requested.'
        write(6,*)
        end if
        write(6,*)
        write(6,*)
        write(6,*) ' For gamma decay run `NuTrx Nucleusi Nucleusf mbmefx... mbmefy...`'
        write(6,*) ' For all others run `NuTrx Nucleusi Nucleusf mbmef...`'
        write(6,*) ' But for clusters run `NuClx Nucleusi Nucleusf mbmef...`'
        write(6,*) ' Repeat for each Ji Jf. '
        
    end subroutine Transition
    
    subroutine get_JT_values(mnj,mxj,text)
        implicit none
        integer, intent(out):: mnj,mxj
        character(len=8), intent(in):: text
        integer:: mmm
        
        write(6,'(2x,a8,a27)',advance='no') text,' Mimimum and maximum 2*Jop '
        read(5,*) mnj,mxj
        mnj=abs(mnj)
        mxj=abs(mxj)
        if (mnj >mxj ) then
            mmm=mxj
            mxj=mnj
            mnj=mmm
        end if
        
    end subroutine get_JT_values
    
                
end module Transitions    

module DefourZuker

    use OutputNus
    use Extra
    
    implicit none
    
    contains
    
    subroutine Interaction()
        implicit none
        integer:: res
    
    
        if (help) then
            write(6,*) '  This is a new option based on ideas in PRC,54,1641(1996)'
            write(6,*) '  It allows you to put a +/- cut around e=0.0 MeV eigenvalues in ref.'
            write(6,*) '  It can be used for Hcm removal ie eliminate mixing with spurious states.'
        end if    
        write(6,'(a52)',advance='no') '  Do you want to use DufourZuker truncation/Hcm y/n '
        read(5,'(a1)') ans
        !write(123,*) ans,'  Do you want to use DufourZuker truncation/Hcm y/n '
        if (ans=='y'.or.ans=='Y') then
            write(6,*) ' You need an isospin interaction and an isospin model space.'
            write(6,*) ' You will be given option to create model space here.'
            write(6,'(a29)',advance='no') '  Do you wish to proceed y/n '
            read(5,'(a1)') ans
            !write(123,*) ans, '  Do you wish to proceed y/n '
            if (ans=='y') then
                res=run('NuDZOper',' ')
                if (res==-1) stop
            end if
        end if
    
    
    end subroutine Interaction
    
end Module DefourZuker         

!  NuShellX.f90 
!
!  FUNCTIONS:
!  NuShellX      - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: NuShellX
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program NuShellX
    
    use OutputNus
    use Levels
    use Transitions
    use ModelSpace
    use Extra
    use OutputNuShell
    use DefourZuker

    implicit none

    ! Variables
    character(len=4):: buffer
    integer:: status
    
    ! Body of NuShellX
    help=.false.
    call get_cmd_arg(1,buffer,status)
    if (status>0) then
        if (buffer=='help') then
            help=.true.
        end if
    end if
    
    call output_headerX(6)
    write(6,*) 
    write(6,*) ' To get help messages use switch help ie run `NuShellX help`.'
    if (help) then
    write(6,*)
    write(6,*) ' To run NuShellX you need to decide on a model space and interaction.'
    write(6,*) ' If you use an Oxbash .int file check that for shells and their order.'
    write(6,*) ' Otherwise you will need to make a .int file from data in literature.'
    write(6,*) ' To do this use an Oxbash file as an example.'
    write(6,*) ' You can use the MSDI interaction that NuShellX will calculate.'
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
    open(unit=56,file='ShellX.bat',action='write')
        
10  write(6,*)
    write(6,'(a49)', advance='no') '  (M)odel,(L)evels,(O)ver,(T)rans,(Z)uker,E(x)it '
    read (5,'(a1)') ans
    write(6,*)
    
    if (ans == 'm' .or. ans == 'M') call CreateModel()
    if (ans == 'l' .or. ans == 'L') call Level(.false.)
    if (ans == 'o' .or. ans == 'O') call Level(.true.)
    if (ans == 't' .or. ans == 'T') call Transition()
    if (ans == 'z' .or. ans == 'Z') call Interaction()
    if (ans == 'x' .or. ans == 'X') go to 100
    
    go to 10
    
100 write(6,*) ' Command file ','ShellX.bat',' saved to disk.'
    if (help) write(6,*) ' Rename this file if you want to repeat calculations.'
    if (help) write(6,*) ' Otherwise delete it.'
    close(unit=56)

    end program NuShellX

