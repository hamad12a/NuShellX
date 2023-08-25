module TrdCharX

use InputNuShell

implicit none


character(len=1):: typea='?',typeb='?'

contains

    
             
    function no_shellsp()
    implicit none
    integer:: no_shellsp
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
    