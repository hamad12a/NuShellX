module SOrder

! Note stringOrder and other Order progs share node array
! Do not use more than one without care in same program

use OrderParameters
use PartOrder

implicit none

type string_tree
    character(len=8):: string
    integer:: node
    integer:: greater
    integer:: lesser
end type string_tree


interface operator (<)
    module procedure sls
end interface

interface operator (>)
    module procedure sgs
end interface

    
type(string_tree),dimension(:),allocatable:: string_order
integer,private,save:: next_string=1
integer:: nstring

contains

    elemental function sls(a,b)
    implicit none
    character(len=8),intent(in):: a,b
    logical:: sls
    sls=llt(a(1:8),b(1:8))
    end function sls
    
    elemental function sgs(a,b)
    implicit none
    character(len=8),intent(in):: a,b
    logical:: sgs
    sgs=lgt(a(1:8),b(1:8))
    end function sgs
    
    function read_next_string()
        implicit none
        integer:: read_next_string
        
        read_next_string=next_string
        
    end function read_next_string
        
    subroutine reset_string_order()
        implicit none
        
        next_string=1
        string_order(1)%string=' '
        string_order(1)%node=1
        string_order(1)%greater=0
        string_order(1)%lesser=0
        
    end subroutine reset_string_order
    
    subroutine add_string(string,new,iequ)
        implicit none
        character(len=8),intent(in):: string
        logical,optional,intent(out):: new
        integer,optional,intent(out):: iequ
        integer:: i
        if (present(new)) new=.true.
        if (present(iequ)) iequ=0
        !store string in a binary tree
        
        if (next_string == 1) then
            string_order(1)%string=string
            next_string=next_string+1
            if (next_string > nstring) call string_message()
            return
        end if
        
        i=1
        
10      if (string == string_order(i)%string)  then
            if (present(new)) new=.false.
            if (present(iequ)) iequ=i
            return
        end if
        if (string > string_order(i)%string) then
            if (string_order(i)%greater == 0) then
                    string_order(i)%greater=next_string
                    string_order(next_string)%node=i
                    string_order(next_string)%string=string
                    string_order(next_string)%greater=0
                    string_order(next_string)%lesser=0
                    next_string=next_string+1
                    if (next_string > nstring) call string_message()
                    return
            else
                i=string_order(i)%greater
                go to 10
            end if
        end if
        if (string < string_order(i)%string) then
            if (string_order(i)%lesser == 0) then
                    string_order(i)%lesser=next_string
                    string_order(next_string)%node=i
                    string_order(next_string)%string=string
                    string_order(next_string)%greater=0
                    string_order(next_string)%lesser=0
                    next_string=next_string+1
                    if (next_string > nstring) call string_message()
                    return
            else
                i=string_order(i)%lesser
                go to 10
            end if
         end if
         
         ! should never reach here
         
         print *, ' Error 2 in module StringOrder'                      !ERROR2
         
    end subroutine add_string
    
    
    function find_string(string,found)
        implicit none
        character(len=8),intent(in):: string
        logical,intent(out)::found
        integer:: find_string
        integer:: i
        
        !find string in a binary tree
        found=.false.
        find_string=0
        if (next_string == 1) then
            return
        end if
        
        i=1
        
10      if (string == string_order(i)%string)  then
            found=.true.
            find_string=i
            return
        end if
        if (string > string_order(i)%string) then
            if (string_order(i)%greater == 0) then
                    return
            else
                i=string_order(i)%greater
                go to 10
            end if
        end if
        if (string < string_order(i)%string) then
            if (string_order(i)%lesser == 0) then
                    return
            else
                i=string_order(i)%lesser
                go to 10
            end if
         end if
         
         ! should never reach here
         
         print *, ' Error 4 in module StringOrder'                      !ERROR4
         
    end function find_string
    
                        
    function get_next_string(done)
        implicit none
        integer:: get_next_string
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
            get_next_string=0
            return
        end if
        if (gl == -1) go to 20
        if (gl == +1) go to 30
        
20      if (string_order(i)%lesser /=0 ) then 
            i=string_order(i)%lesser                  !go down tree
            go to 20
        else 
            get_next_string=i                             !bottom of tree, return value
            n=n+1
            if (n == next_string) done=.true.
            gl=+1                                   !ready to go back up
            return
        end if
        
30      if (string_order(i)%greater /=0 ) then
            i=string_order(i)%greater                 !found upward branch
            nd=nd+1
            if (nd > nonode) then
                print *, ' Increase parameter nonode'
                stop
            end if
            node(nd)=string_order(i)%node             !push current node on stack
            if (gl == -1) then
                get_next_string=string_order(i)%node        !last movement was down, return value
                return
            else
                gl=-1                               !get set to go down branch
                go to 20
            end if
        else
40          i=string_order(i)%node
            if (i == node(nd)) then                 !if node equals to of sstack
                nd=nd-1                             !pop stack
                if (nd < 1) then
                    done=.true.                     !reached root node, done
                    get_next_string=0
                    return
                end if
                go to 40                            !go back up
            else 
                go to 50
            end if
50          n=n+1                                   !return value at top of branch
            get_next_string=i
            if (n == next_string) done=.true.
            gl=+1                                   !prepare to go up again or pop
            return
        end if
        
        done=.true.                                 !should never reach here NO ERROR TRAP!!!  
              
       print *, ' Error 4 in module StringOrder'                      !ERROR4
       stop
            
    end function get_next_string
                        
    
    subroutine string_message()
        implicit none
        
        print * ,' Increase parameter nstring'
        stop
        
    end subroutine string_message        
                
end module SOrder
