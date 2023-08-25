module TrdChar3

use TrdChar

implicit none



type trx
    integer,dimension(5):: qux
end type trx




contains


    pure function tfind3(sarray,n,string)
    implicit none
    integer,intent(in):: n
    character(len=12),intent(in):: string
    integer:: tfind3
    character(len=12),dimension(*),intent(in):: sarray
    integer:: in
    real:: rn,drn
    if (n==0) then
        tfind3=0
        return
    end if
    if (sarray(n)(1:12)==string(1:12)) then 
        tfind3=n
        return
    end if
    rn=n
    rn=rn/2.0
    drn=rn
    do
        in=nint(rn)
        if (sarray(in)(1:12)==string(1:12)) then
            tfind3=in
            exit
        end if
        drn=drn/2.0
        if (abs(drn) < 0.5) then
            tfind3=0
            exit
        end if
        if (lgt(string(1:12),sarray(in)(1:12))) then
            rn=rn+drn
        else
            rn=rn-drn
        end if
    end do
    end function tfind3
    
    elemental function char3(i)
    implicit none
    integer,intent(in)::i
    character(len=3)::char3
    char3(1:1)=char(iand(mask,i))
    char3(2:2)=char(iand(mask,ishft(i,-8)))
    char3(3:3)=char(iand(mask,ishft(i,-16)))
    end function char3
    
    elemental function ichar3(c3)
    implicit none
    character(len=3),intent(in):: c3
    integer:: ichar3
    ichar3=0
    ichar3=ichar(c3(1:1))
    ichar3=ior(ishft(ichar(c3(2:2)),8),ichar3)
    ichar3=ior(ishft(ichar(c3(3:3)),16),ichar3)
    end function ichar3
     
    function ichara3(chra,n)
    implicit none
    integer,intent(in)::n
    character(len=n),intent(in):: chra
    integer,dimension(12):: ichara3
    integer::i
    ichara3=0
    do i=1,n
        ichara3(i)=ichar(chra(i:i))
    end do
    end function ichara3
       
    function chara3(iarray,n)
    implicit none
    integer,intent(in)::n
    character(len=n):: chara3
    integer,dimension(n),intent(in):: iarray
    integer::i
    
    do i=1,n
        chara3(i:i)=char(iarray(i))
    end do
    end function chara3
    
end module TrdChar3    
    