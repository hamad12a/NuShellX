!  NuMake.f90 
!
!  FUNCTIONS:
!  NuMake      - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: NuMake
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program NuMake
    
    USE IFPORT
    USE IFCORE

    implicit none

    ! Variables
       character(len=25):: filename,project
       CHARACTER(80) file
       INTEGER(KIND=INT_PTR_KIND( )) handle
       TYPE(FILE$INFO) info
       INTEGER(4) length
       character(len=1):: au
       character(len=40):: aa
       character(len=40),dimension(100):: a
       character(len=70):: f
       integer,dimension(100):: jmax
       integer:: i,j,iq,imax,jjmax,ip
       file(1:12)='*Order.build'
       handle = FILE$FIRST
       length = GETFILEINFOQQ(file, info, handle)
       IF (handle .EQ. FILE$LAST) GO TO 500
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
       if (length > 25) then
         print *, ' Filename too long',length,info%name
         stop
       end if 
       filename=info%name(1:length)
       filename=adjustr(filename)
       open(unit=10,file=filename)
       open(unit=20,file='make.bat')
       filename=adjustl(filename)
       do i=1,24
         if (filename(i:i)=='O') exit
       end do
       ip=i-1
       project(1:ip)=filename(1:ip)

       i = 0
       do 
       read(10,'(a40)') aa
       if(aa.eq.' ') exit
       i = i + 1
       a(i) = aa
       do j = 40,1,-1
         if(.not.(aa(j:j)==' ')) exit
       end do
       jmax(i) = j
       end do
       iq = 0
       read(10,*,end=300)
       read(10,'(a1)',end=300) au
       if(au=='/') iq = 1
       read(10,'(a1)',end=300) au
       if(au=='/') iq = 2
300    imax = i
       write(20,'(a79)',advance='no') 'ifort /Ob2 /O3 /assume:buffered_io /heap-arrays /Qprec-sqrt /Qftz /Qipo /QaxWT '
       if(iq == 1 .or. iq == 2) write(20,'(a  )',advance='no') '/Qopenmp'
       do i = 1,imax
       aa = a(i)
       if (aa(1:2)=='..') then
         aa(1:2)='n:'
         write(20,'(1x,<jmax(i)>a1)',advance='no') (aa(j:j),j=1,jmax(i))
       else
         write(20,'(1x,a3,<ip>a1,a1,<jmax(i)>a1)',advance='no') 'n:\', &
                            (project(j:j),j=1,ip),'\',(aa(j:j),j=1,jmax(i))
       end if
       end do
       jjmax = jmax(imax)-4
       write(20,'(a6,<jjmax>a1)',advance='no') ' /exe:',(aa(j:j),j=1,jjmax)
       if(iq.eq.2) write(20,'(a23)', advance='no') ' /link /STACK:100000000'
       write(20,*)
       close(10)
       close(20)


    ! Body of NuMake

500  end program NuMake

