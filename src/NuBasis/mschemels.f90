module Mscheme

use Parameters
use Shells
use Partition
use Machine

implicit none

!integer,parameter:: max_sj2p1 = max_sj2+1
integer,parameter:: max_sj2p1 = max_sj2
integer,parameter:: max_sj2p1d2 =max_sj2p1/2
integer,parameter:: m1dim =((max_sj2p1-1)/nbytes+1)*nbytes
integer,parameter:: m1d2dim =((max_sj2p1d2-1)/nbytes+1)*nbytes


!integer(kind=2),dimension(m1d2dim,m1dim):: sm2v
integer(kind=2),dimension(0:m1d2dim,0:m1dim):: sm2v
!integer(kind=2),dimension(m1d2dim,m1d2dim):: sm2vp
integer(kind=2),dimension(0:m1d2dim,0:m1d2dim):: sm2vp

contains

    subroutine setup_mscheme()
        implicit none
        integer:: i,j,k,m2
        
        do i=0,max_sj2p1d2
!        do i=1,max_sj2p1d2
            j=0
            k=0
            do m2=-sj2f(i,0),sj2f(i,0),2  ! any other calls to this elsewhere need 0 adding
!            do m2=-sj2f(i),sj2f(i),2
                j=j+1
                sm2v(i,j)=m2
                if (m2 >= 0) then           ! not sure here
!                if (m2 > 0) then
                    k=k+1
                    sm2vp(i,k)=m2
                end if
            end do
        end do
        
    end subroutine setup_mscheme
    
    function sj2f(i,n)
        implicit none
        integer:: sj2f
        integer,intent(in):: i
        integer,optional,intent(in):: n
        
        if (present(n)) then
            if (btest(n,0)) then
                sj2f=i*2-1
            else
                sj2f=i*2
            end if
        else
            sj2f=i*2-1
        end if
        
    end function sj2f
    
    function isj2(j2)
        implicit none
        integer:: isj2
        integer,intent(in):: j2
        
        if (btest(j2,0)) then
            isj2=(j2+1)/2
        else
            isj2=j2/2
        end if
        
    end function isj2
    
    function cParity(part)
        implicit none
        integer:: cParity
        type(spartition),intent(in):: part
        integer:: i
        
        cParity=1
        do i=1,no_shells
            if (part%shell(i) > 0) then
                cParity=cParity*(sp(i))**abs(part%shell(i))
            end if
        end do
        
        if (cParity ==0) then
            print *, part%shell(1:no_shells)
            print *, 'Parity',cParity
            stop
        end if
        
    end function cParity
    
    function M_graph(mpart)
        implicit none
        logical:: M_graph
        type(spartition),intent(in):: mpart
        
        if (pintegral(mpart) < 0) then
            M_graph=.false.
            return
        end if
        
        M_graph=.true.
        
    end function M_graph
        
               
end module Mscheme