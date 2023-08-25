! This file contains all the modules specific to program pPSLevel
! for NuShell,NuShell@MSU,SunShell

module PSGraphics 

use LanczosParameters

implicit none

character(len=1):: xtra

contains

    subroutine psgraphic(tot,allv,allm,indx,allj,allt,allp,tPar,Nucl,gs)    
        
        implicit none
        real(kind=rl),dimension(*):: allv,allm
        real(kind=rl):: range,allvindx1,allvindxx,Emin,Emax,gs
        character(len=1),dimension(*):: allp
        integer:: tot,tPar
        character(len=6):: Nucl
        integer,dimension(*):: indx,allj,allt
        integer:: i
        
        do i=tot,1,-1
            allvindxx=allv(indx(i))
            if (allvindxx /= real(0,rl)) exit
        end do
        Emin=real(0,rl)
        if (gs /= 0.0) then
        range = abs(gs-allvindxx)
        else
        range = abs(allv(indx(1))-allvindxx)
        end if
        read(2,*) Emin,Emax
        if (Emax > 0.0) range=abs(Emax-Emin)
        if (gs /= 0.0) then
        allvindx1=gs
        else
        allvindx1=allv(indx(1))
        end if
        do i=1,tot
            if (allv(i) /= real(0,rl)) then
                allm(i)=(allv(i)-allvindx1)
                allv(i)=-(allv(i)-allvindx1-Emin)/range*940 - 30
            end if
        end do       
!  Build as a Graphics ap.
         open(unit=45,file=Nucl//'.eps',action='write')
         write(45,'(a23)') '%!PS-Adobe-2.0 EPSF-2.0'
         write(45,'(a31)') '%%BoundingBox:   1  592  1  838'
         write(45,'(a23)') '%%Creator: W. D. M. Rae'
         write(45,'(a13)') '%%EndComments'
         CALL psdrawlines(tot,allv,allm,allj,allt,allp,indx,Nucl,allvindx1,tPar)
 
    end subroutine psgraphic
    
    function h(x)
    implicit none
    real,intent(in)::x
    real::h
    h=593.0*x/1000.0
    end function h
    
    function v(y)
    implicit none
    real,intent(in)::y
    real::v
    v=838.0*y/1000.0
    end function v

!  DRAWLINES - This subroutine draws a box and
!  energy level lines.
 SUBROUTINE psdrawlines(tot,allv,allm,allj,allt,allp,indx,nucl,gs ,tpar)
    integer  tot,tpar
    character(len=6) nucl
    character(len=19) text
    character(len=1) par
    real (kind=rl) gs
    integer, dimension (*) :: indx
    real(kind=rl),dimension (*) :: allv,allm
    integer,dimension (*) :: allj,allt
    character(len=1),dimension(*):: allp
    real int2allv,rone,rtwo,rthree,r0,rrhs
    logical rhs,three,two,one
    integer i
    integer:: results
! 
!  Draw the box.
    if (tpar > 0) par='+'
    if (tpar < 0) par='-'
    write(45,*) '1.0 setlinewidth'
    write(45,'(a8,2f8.2,a8,2f8.2,a14)')'newpath ',h(2.0),v(2.0),' moveto ',h(999.0),v(2.0),' lineto stroke'
    write(45,'(a8,2f8.2,a8,2f8.2,a14)')'newpath ',h(999.0),v(2.0),' moveto ',h(999.0),v(999.0),' lineto stroke'
    write(45,'(a8,2f8.2,a8,2f8.2,a14)')'newpath ',h(2.0),v(999.0),' moveto ',h(999.0),v(999.0),' lineto stroke'
    write(45,'(a8,2f8.2,a8,2f8.2,a14)')'newpath ',h(2.0),v(2.0),' moveto ',h(2.0),v(999.0),' lineto stroke'
!  This sets the new origin to 0 for x and 1000 for y. 
!  Draw the lines and annotate.
    write(45,'(a24,f8.2,a,2f8.2,a7)') '/Helvetica    findfont  ',32.5,' scalefont setfont',&
                            h(10.0),v(700.0),' moveto'
    write(45,'(f7.2,a7)') 0.0,' rotate' 
    write(45,'(a1,a6,a5)') '(',adjustl(Nucl),')show'                           
    write(45,'(a24,f8.2,a,2f8.2,a7)') '/Helvetica    findfont  ',18.5,' scalefont setfont',&
                            h(10.0),v(500.0),' moveto'
    write(45,'(f7.2,a7)') 0.0,' rotate' 
    write(45,'(a1,a5,a5)') '(','Int '//xtra,')show'                           
    write(text,'(f8.3,a6)') gs,'MeV gs'
    write(45,'(a24,f8.2,a,2f8.2,a7)') '/Helvetica    findfont  ',13.5,' scalefont setfont',&
                            h(10.0),v(300.0),' moveto'
    write(45,'(f7.2,a7)') 0.0,' rotate' 
    write(45,'(a1,a19,a5)') '(',text,')show'                           
    write(text,'(a19)')'Ex(MeV)    J     T '
    write(45,'(a24,f8.2,a,2f8.2,a7)') '/Helvetica    findfont  ',8.5,' scalefont setfont',&
                            h(490.0),v(15.0),' moveto'
    write(45,'(f7.2,a7)') 0.0,' rotate' 
    write(45,'(a1,a19,a5)') '(',text,')show'                           
    rhs=.false.
    one=.false.
    two=.false.
    three=.false.
    r0=0.0; rone=0.0; rtwo=0.0;
    rthree=0.0; rrhs=0.0;
    do i=1,tot
        if (i/=1) then
            if (abs(allv(indx(i))-r0) < 7.) then
                if (abs(allv(indx(i))-rrhs) > 7.) then
                    rhs=.true.
                    go to 10
                end if
                if (abs(allv(indx(i))-rone) > 7.) then
                    one=.true.
                    go to 10
                end if
                if (abs(allv(indx(i))-rtwo) > 7.) then
                    two=.true.
                    go to 10
                end if
                if (abs(allv(indx(i))-rthree) > 7.) then
                    three=.true.
                    go to 10
                end if
            end if    
        end if
10      if (allv(indx(i)) == real(0,rl)) go to 50
        if (allp(indx(i)) == ' ') go to 50
! no more levels
        if (btest(allj(indx(i)),0)) then
            write(text,'(f8.3,1x,i2,a2,a1,1x,i2,a2)')&
                     allm(indx(i)),allj(indx(i)),'/2',allp(indx(i)),abs(allt(indx(i))),'/2'
        else
            write(text,'(f8.3,1x,i2,a2,a1,1x,i2,a2)')&
                 allm(indx(i)),int(allj(indx(i))/2),'  ',allp(indx(i)),int(abs(allt(indx(i))/2)),'  '
        end if
        int2allv=allv(indx(i))
! stagger text for closely spaced levels
        if (indx(i+1)==0)  indx(i+1)=indx(i)
        if (rhs ) then
            write(45,'(a24,f8.2,a,2f8.2,a7)') '/Helvetica    findfont  ',8.5,' scalefont setfont',&
                            h(830.0),v(-int2allv),' moveto'
            write(45,'(f7.2,a7)') 0.0,' rotate' 
            write(45,'(a1,a19,a5)') '(',text,')show'                           
            rhs=.false.
            rrhs=allv(indx(i))
        else     if (three) then
            write(45,'(a24,f8.2,a,2f8.2,a7)') '/Helvetica    findfont  ',8.5,' scalefont setfont',&
                            h(200.0),v(-int2allv),' moveto'
            write(45,'(f7.2,a7)') 0.0,' rotate' 
            write(45,'(a1,a19,a5)') '(',text,')show'                           
            three=.false.
            rthree=allv(indx(i))
        else     if (two) then
            write(45,'(a24,f8.2,a,2f8.2,a7)') '/Helvetica    findfont  ',8.5,' scalefont setfont',&
                            h(300.0),v(-int2allv),' moveto'
            write(45,'(f7.2,a7)') 0.0,' rotate' 
            write(45,'(a1,a19,a5)') '(',text,')show'                           
            two=.false.
            rtwo=allv(indx(i))
        else     if (one) then
            write(45,'(a24,f8.2,a,2f8.2,a7)') '/Helvetica    findfont  ',8.5,' scalefont setfont',&
                            h(400.0),v(-int2allv),' moveto'
            write(45,'(f7.2,a7)') 0.0,' rotate' 
            write(45,'(a1,a19,a5)') '(',text,')show'                           
            one=.false.
            rone=allv(indx(i))
        else
            write(45,'(a24,f8.2,a,2f8.2,a7)') '/Helvetica    findfont  ',8.5,' scalefont setfont',&
                            h(500.0),v(-int2allv),' moveto'
            write(45,'(f7.2,a7)') 0.0,' rotate' 
            write(45,'(a1,a19,a5)') '(',text,')show'                           
            r0=allv(indx(i))
        end if
!draw level
        write(45,'(a8,2f8.2,a8,2f8.2,a14)')'newpath ',h(600.0),v(-int2allv),' moveto ',&
                    h(800.0),v(-int2allv),' lineto stroke'
50   end do
     write(text,'(a19)')'Ex(MeV)    J     T '
     write(45,'(a24,f8.2,a,2f8.2,a7)') '/Helvetica    findfont  ',8.5,' scalefont setfont',&
                            h(10.0),v(15.0),' moveto'
     write(45,'(f7.2,a7)') 0.0,' rotate' 
     write(45,'(a1,a52,a5)')&
             '(',' NuShell V5.0 R2.015 (c) W D M Rae, Garsington, 2007',')show'                           
     write(45,'(a8)') 'showpage'
     
 END SUBROUTINE psdrawlines

end Module PSGraphics

!  NuPSLevel.f90 
!
!  FUNCTIONS:
!  NuPSLevel      - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: NuPSLevel
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program NuPSLevel
    
    use LanczosParameters
    use InputNuShell
    use OutputNuShell
    use Files
    use PSGraphics
    use TriDiag   
    use Extra
    

    implicit none

    ! Variables
    integer,parameter:: max_tot=2000000
    real(kind=rl):: gs
    real(kind=rl),dimension(max_tot):: all_level,all_mev
    integer,dimension(max_tot):: all_levelj,all_levelt,indx
    character(len=1),dimension(max_tot):: all_levelp
    integer:: tot_level,npt
    logical:: done
    integer:: j,t,dim,twrite
    character(len=5):: file
    character(len=6):: inExt
    character(len=1):: blank
    character(len=12):: dummy2
    character(len=14):: filename
    ! Body of NuPSLevel
    tot_level=0
    call output_header(6)    
    print *, ' NuPSLevel Program'
    file(1:5)='*.lev'
2   call get_next_file(file,filename,done)
    if (done) go to 20
    open(unit=500,file=filename)
    do dim=1,6
        read(500,*)
    end do
    read(500,'(a1,a1,a6,a1)') blank,blank,Nucleus,xtra
    do dim=1,5
        read(500,*)
    end do
    read(500,'(a12,i12,i12)') dummy2,j,t
    do dim=1,3
        read(500,*)
    end do
    twrite=t
    do
    read(500,'(5(f10.3,a1))',end=10,err=10) &
               (all_level(npt),all_levelp(npt),npt=tot_level+1,tot_level+5)
    do npt=tot_level+1,tot_level+5
        if (all_level(npt) /= real(0,rl)) then
            tot_level=tot_level+1
        end if
        all_levelj(npt)=j
        all_levelt(npt)=twrite
    end do
    
    if (tot_level+5 > max_tot) stop ' Increase max_tot'
    end do
10  go to 2
    
20  call pssetup_prm()
    open(unit=2,file='NuPSLevel.prm')
    read(2,*) gs
    
    call indexx(tot_level,all_level,indx)
    call psgraphic(tot_level,all_level,all_mev,indx,all_levelj,all_levelt,all_levelp,tParity,Nucleus,gs)   
    
    end program NuPSLevel
    
    subroutine pssetup_prm()
        use LanczosParameters
        implicit none
        real(kind=rl):: gs=0.0
        real(kind=rl):: Emin=0.0
        real(kind=rl):: Emax=-1.0
        logical:: exst
        
        inquire(file='NuPSLevel.prm',exist=exst)
        if (exst) return
        open(unit=2,file='NuPSLevel.prm',action='write')
        write(2,'(f10.6)') gs
        write(2,'(2f10.6)') Emin,Emax
        close(unit=2)
        
    end subroutine pssetup_prm
        