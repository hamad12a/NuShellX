! This file contains all the modules specific to program NuLevel for NuShell

module Graphics 

use LanczosParameters
use IFQWIN

implicit none

character(len=1):: xtra

contains

    subroutine graphic(tot,allv,allm,indx,allj,allt,allp,tPar,Nucl,gs)    
        USE IFQWIN
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
         CALL graphicsmode( )
         CALL drawlines(tot,allv,allm,allj,allt,allp,indx,Nucl,allvindx1,tPar)
 
    end subroutine graphic

  SUBROUTINE graphicsmode()
    USE IFQWIN
    LOGICAL             modestatus
    INTEGER(2)          maxx, maxy
    TYPE (windowconfig) myscreen
    COMMON              maxx, maxy
!  Set highest resolution graphics mode.
    myscreen.numxpixels=-1
    myscreen.numypixels=-1
    myscreen.numtextcols=-1
    myscreen.numtextrows=-1
    myscreen.numcolors=-1
    myscreen.fontsize=-1
    myscreen.title = "NuShell"C ! blank
    modestatus=SETWINDOWCONFIG(myscreen)
!  Determine the maximum dimensions.
    modestatus=GETWINDOWCONFIG(myscreen)
    maxx=myscreen.numxpixels - 1
    maxy=myscreen.numypixels - 1
  END SUBROUTINE graphicsmode
  
! NEWX - This function finds new x-coordinates.
  FUNCTION newx( xcoord )
    INTEGER(2) newx
    INTEGER(2) xcoord, maxx, maxy
    REAL(4) tempx
    COMMON maxx, maxy
    tempx = maxx / 1000.0
    tempx = xcoord * tempx + 0.5
    newx = tempx
  END FUNCTION newx
! NEWY - This function finds new y-coordinates.
!
  FUNCTION newy( ycoord )
    INTEGER(2) newy
    INTEGER(2) ycoord, maxx, maxy
    REAL(4) tempy
    COMMON maxx, maxy
    tempy = maxy / 1000.0
    tempy = ycoord * tempy + 0.5
    newy = tempy
  END FUNCTION newy
    
!  DRAWLINES - This subroutine draws a box and
!  energy level lines.
 SUBROUTINE drawlines(tot,allv,allm,allj,allt,allp,indx,nucl,gs ,tpar)
    USE IFQWIN
    integer  tot,tpar
    character(len=6) nucl
    character(len=19) text
    character(len=1) par
    real (kind=rl) gs
    integer, dimension (*) :: indx
    real(kind=rl),dimension (*) :: allv,allm
    integer,dimension (*) :: allj,allt
    character(len=1),dimension(*):: allp
    integer(2) int2allv
    real rone,rtwo,rthree,r0,rrhs
    logical rhs,three,two,one
    integer i
    INTEGER(2) fontnum, numfonts
    INTEGER(2)      status, maxx, maxy
    TYPE (xycoord)  xy
    TYPE (wxycoord) wxylr,wxyul
    COMMON          maxx, maxy
    integer:: results
! 
!  Draw the box.
    if (tpar > 0) par='+'
    if (tpar < 0) par='-'
    status = RECTANGLE( $GBORDER, INT2(0), INT2(0), maxx, maxy )
    CALL SETVIEWORG( INT2(0), newy( INT2( 1000 ) ), xy )     
!  This sets the new origin to 0 for x and 1000 for y. 
!  Draw the lines and annotate.
    numfonts = INITIALIZEFONTS ( )
!  Set typeface to Arial, character height to 72,
!  character width to 40, and italic
    fontnum = SETFONT ('t''Arial''h50w20i')
    CALL MOVETO (newx(INT2(10)), newy(INT2(-700)), xy)
! Title  1
    CALL OUTGTEXT(adjustl(Nucl))
    fontnum = SETFONT ('t''Arial''h30w15i')
    CALL MOVETO (newx(INT2(10)), newy(INT2(-500)), xy)
! Title
    CALL OUTGTEXT('Int '//xtra)
    write(text,'(f8.3,a6)') gs,'MeV gs'
    fontnum = SETFONT ('t''Arial''h18w10i')
    CALL MOVETO (newx(INT2(10)), newy(INT2(-300)), xy)
! gs binding energy
    CALL OUTGTEXT(text)
    fontnum = SETFONT ('t''Arial''h12w06i')
    write(text,'(a19)')'Ex(MeV)    J     T '
    CALL MOVETO (newx(INT2(490)), newy(INT2(-15)), xy)
! Ex 2*J 2*T
    rhs=.false.
    one=.false.
    two=.false.
    three=.false.
    r0=0.0; rone=0.0; rtwo=0.0;
    rthree=0.0; rrhs=0.0;
    CALL OUTGTEXT(text)
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
                 allm(indx(i)),allj(indx(i))/2,'  ',allp(indx(i)),abs(allt(indx(i))/2),'  '
        end if
        int2allv=INT2(allv(indx(i)))
        if (indx(i+1)==0)  indx(i+1)=indx(i)
! stagger text for closely spaced levels
        if (rhs ) then
            CALL MOVETO( newx(INT2(830)), newy(int2allv), xy )
            rhs=.false.
            rrhs=allv(indx(i))
        else     if (three) then
            CALL MOVETO( newx(INT2(200)), newy(int2allv), xy )
            three=.false.
            rthree=allv(indx(i))
        else     if (two) then
            CALL MOVETO( newx(INT2(300)), newy(int2allv), xy )
            two=.false.
            rtwo=allv(indx(i))
        else     if (one) then
            CALL MOVETO( newx(INT2(400)), newy(int2allv), xy )
            one=.false.
            rone=allv(indx(i))
        else
            CALL MOVETO( newx(INT2(500)), newy(int2allv), xy )
            r0=allv(indx(i))
        end if
! out Ex 2J,pi,2T
        CALL OUTGTEXT(text)
!draw level
        CALL MOVETO( newx(INT2(600)), newy(int2allv), xy )
        status = LINETO( newx(INT2( 800 )), newy(int2allv))
50   end do
     CALL MOVETO (newx(INT2(10)), newy(INT2(-15)), xy)
     CALL OUTGTEXT(' NuShell V5.0 R2.015 (c) W D M Rae, Garsington, 2007')
     CALL GETWINDOWCOORD (newx(int2(0)), newy(int2(-1000)), wxyul)
     CALL GETWINDOWCOORD (newx(int2(1000)), newy(int2(0)), wxylr)
     results = SAVEIMAGE_W(nucl//'.bmp',wxyul%wx,wxyul%wy,wxylr%wx,wxylr%wy)

 END SUBROUTINE drawlines

end Module Graphics

!  NuLevel.f90 
!
!  FUNCTIONS:
!  NuLevel      - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: NuLevel
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program NuLevel
    
    use Parameters
    use LanczosParameters
    use InputNuShell
    use OutputNuShell
    use Shells
    use Files
    use Graphics
    use TriDiag   
    USE IFPORT
    USE IFCORE
 
    implicit none

    ! Variables
    integer,parameter:: max_tot=2000000
    real(kind=rl):: gs
    real(kind=rl),dimension(max_tot):: all_level,all_mev
    integer,dimension(max_tot):: all_levelj,all_levelt,indx
    character(len=1),dimension(max_tot):: all_levelp
    integer:: tot_level,npt

    integer:: j,t,dim,twrite

    character(len=6):: inExt
    character(len=1):: blank
    character(len=12):: dummy2
    character(len=14):: filename
    CHARACTER(80) file
    INTEGER(KIND=INT_PTR_KIND( )) handle
    TYPE(FILE$INFO) info
    INTEGER(4) length
    ! Body of NuLevel
    tot_level=0
    call output_header(6)    
    print *, ' NuLevel Program'

    file(1:5)='*.lev'
    handle = FILE$FIRST
2   length = GETFILEINFOQQ(file, info, handle)
    IF (handle .EQ. FILE$LAST) GO TO 20
    IF ((handle .EQ. FILE$ERROR)) THEN
        SELECT CASE (GETLASTERRORQQ( ))
           CASE (ERR$NOMEM)
             WRITE (*,*) 'Out of memory'
           CASE (ERR$NOENT)
             WRITE (*,*) ' File not found'
             STOP
           CASE DEFAULT
             WRITE (*,*) 'Invalid file or path name'
             STOP
        END SELECT
    END IF
    if (length > 14) then
        print *, ' Filename too long',length,info%name
        stop
    end if 
    filename=info%name(1:length)
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
    
20  call setup_prm()
    open(unit=2,file='NuLevel.prm')
    read(2,*) gs
    
    call indexx(tot_level,all_level,indx)
    call graphic(tot_level,all_level,all_mev,indx,all_levelj,all_levelt,all_levelp,tParity,Nucleus,gs)   
    
    end program NuLevel
    
    subroutine setup_prm()
        use LanczosParameters
        implicit none
        real(kind=rl):: gs=0.0
        real(kind=rl):: Emin=0.0
        real(kind=rl):: Emax=-1.0
        logical:: exst
        
        inquire(file='NuLevel.prm',exist=exst)
        if (exst) return
        open(unit=2,file='NuLevel.prm',action='write')
        write(2,'(f10.6)') gs
        write(2,'(2f10.6)') Emin,Emax
        close(unit=2)
        
    end subroutine setup_prm
        