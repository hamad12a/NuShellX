module SOrder4

! Note stringOrder and other Order progs share node array
! Do not use more than one without care in same program

use OrderParameters
use PartOrder

implicit none

type string_tree4
    character(len=18):: string
    integer:: node
    integer:: greater
    integer:: lesser
end type string_tree4



interface operator (<)
    module procedure sls4
end interface

interface operator (>)
    module procedure sgs4
end interface

    
type(string_tree4),dimension(:),allocatable:: string_order4
integer,private,save:: next_string=1
integer:: nstring

contains

    elemental function sls4(a,b)
    implicit none
    character(len=18),intent(in):: a,b
    logical:: sls4
    sls4=llt(a(1:18),b(1:18))
    end function sls4
    
    elemental function sgs4(a,b)
    implicit none
    character(len=18),intent(in):: a,b
    logical:: sgs4
    sgs4=lgt(a(1:18),b(1:18))
    end function sgs4
    

    function read_next_string4()
        implicit none
        integer:: read_next_string4
        
        read_next_string4=next_string
        
    end function read_next_string4
        
    subroutine reset_string_order4()
        implicit none
        
        next_string=1
        string_order4(1)%string=' '
        string_order4(1)%node=1
        string_order4(1)%greater=0
        string_order4(1)%lesser=0
        
    end subroutine reset_string_order4
    
    subroutine add_string4(string,new,iequ)
        implicit none
        character(len=18),intent(in):: string
        logical,optional,intent(out):: new
        integer,optional,intent(out):: iequ
        integer:: i
        if (present(new)) new=.true.
        if (present(iequ)) iequ=0
        !store string in a binary tree
        
        if (next_string == 1) then
            string_order4(1)%string=string
            next_string=next_string+1
            if (next_string > nstring) call string_message4()
            return
        end if
        
        i=1
        
10      if (string == string_order4(i)%string)  then
            if (present(new)) new=.false.
            if (present(iequ)) iequ=i
            return
        end if
        if (string > string_order4(i)%string) then
            if (string_order4(i)%greater == 0) then
                    string_order4(i)%greater=next_string
                    string_order4(next_string)%node=i
                    string_order4(next_string)%string=string
                    string_order4(next_string)%greater=0
                    string_order4(next_string)%lesser=0
                    next_string=next_string+1
                    if (next_string > nstring) call string_message4()
                    return
            else
                i=string_order4(i)%greater
                go to 10
            end if
        end if
        if (string < string_order4(i)%string) then
            if (string_order4(i)%lesser == 0) then
                    string_order4(i)%lesser=next_string
                    string_order4(next_string)%node=i
                    string_order4(next_string)%string=string
                    string_order4(next_string)%greater=0
                    string_order4(next_string)%lesser=0
                    next_string=next_string+1
                    if (next_string > nstring) call string_message4()
                    return
            else
                i=string_order4(i)%lesser
                go to 10
            end if
         end if
         
         ! should never reach here
         
         print *, ' Error 2 in module StringOrder'                      !ERROR2
         
    end subroutine add_string4
    
    
    function find_string4(string,found)
        implicit none
        character(len=18),intent(in):: string
        logical,intent(out)::found
        integer:: find_string4
        integer:: i
        
        !find string in a binary tree
        found=.false.
        find_string4=0
        if (next_string == 1) then
            return
        end if
        
        i=1
        
10      if (string == string_order4(i)%string)  then
            found=.true.
            find_string4=i
            return
        end if
        if (string > string_order4(i)%string) then
            if (string_order4(i)%greater == 0) then
                    return
            else
                i=string_order4(i)%greater
                go to 10
            end if
        end if
        if (string < string_order4(i)%string) then
            if (string_order4(i)%lesser == 0) then
                    return
            else
                i=string_order4(i)%lesser
                go to 10
            end if
         end if
         
         ! should never reach here
         
         print *, ' Error 4 in module StringOrder'                      !ERROR4
         
    end function find_string4
    
                        
    function get_next_string4(done)
        implicit none
        integer:: get_next_string4
        logical,intent(inout):: done
        integer,save:: i,gl,n,nd
        
        !read out strings's in increasing order of sort variables
        
        if (done) then
            i=1
            n=1
            nd=1
            node=0
            gl=-1
            done=.false.
        end if
        if (next_string==1) then
            done=.true.
            get_next_string4=0
            return
        end if
        if (gl == -1) go to 20
        if (gl == +1) go to 30
        
20      if (string_order4(i)%lesser /=0 ) then 
            i=string_order4(i)%lesser                  !go down tree
            go to 20
        else 
            get_next_string4=i                             !bottom of tree, return value
            n=n+1
            if (n == next_string) done=.true.
            gl=+1                                   !ready to go back up
            return
        end if
        
30      if (string_order4(i)%greater /=0 ) then
            i=string_order4(i)%greater                 !found upward branch
            nd=nd+1
            if (nd > nonode) then
                print *, ' Increase parameter nonode'
                stop
            end if
            node(nd)=string_order4(i)%node             !push current node on stack
            if (gl == -1) then
                get_next_string4=string_order4(i)%node        !last movement was down, return value
                return
            else
                gl=-1                               !get set to go down branch
                go to 20
            end if
        else
40          i=string_order4(i)%node
            if (i == node(nd)) then                 !if node equals to of sstack
                nd=nd-1                             !pop stack
                if (nd < 1) then
                    done=.true.                     !reached root node, done
                    get_next_string4=0
                    return
                end if
                go to 40                            !go back up
            else 
                go to 50
            end if
50          n=n+1                                   !return value at top of branch
            get_next_string4=i
            if (n == next_string) done=.true.
            gl=+1                                   !prepare to go up again or pop
            return
        end if
        
        done=.true.                                 !should never reach here NO ERROR TRAP!!!  
              
       print *, ' Error 4 in module StringOrder'                      !ERROR4
       stop
            
    end function get_next_string4
                        
    
    subroutine string_message4()
        implicit none
        
        print * ,' Increase parameter nstring'
        stop
        
    end subroutine string_message4
                
end module SOrder4
