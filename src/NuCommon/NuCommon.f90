! Common modules for NuShell, NuShell@MSU and SunShell.
! NuShell and NuShell@MSU are Intel x86 base and can compile
! under Intel Fortran 9.1 or later on Windows and Linux in 32 bit
! or 64 bit platforms. SunShell compiles with Sun Studio 12 on Solaris
! and Linux operating sytems using  Sun SPARC or x86 hardware.
! Again both 64 bit or 32 bit platforms can be used with no changes.
! Compiler specific routines are collected here in the module Extra.
! Modules DotProduct , Clebsch and TriDiag have different include
! since Sun do not emulate real*16. These are to be found between
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ! and
! -------------------------------------------------------------- !
! markers below. Changing between Sun and Intel is just a matter of
! commenting out beween Intel and end Intel and uncommenting between
! Sun and end Sun and visa versa. 
! DO NOT COMPILE ALONE

module Machine

implicit none

! integer word size bits 

integer,parameter:: nbits = 64 !  32bits or 64bits usually

! memory access/boundaries bytes=8bits

integer,parameter:: nbytes = 4 !  4bytes or 8bytes usually

!For EMT64 this can be left at 4.

end module Machine

module Parameters

use Machine

implicit none

integer,parameter:: max_no_shells = 30 ! maximum number of neutron and proton shells nbytes*integer
integer,parameter:: nwords = 5        ! no of words per m-vector
integer,parameter:: ik = nbits/8       ! integer kind for m_vector words
integer,parameter:: max_sj2 = 17       ! maximum value of sj2
integer,parameter:: max_no_partitions  = 10000 ! maximum no of shell partitions
integer,parameter:: max_combine = 10000   ! maximum combinations >= (2*(max_sj2+1))!/(((max_sj2+1))!)^2
integer,parameter:: max_allmvector = 600000 ! maximum no of m_vectors per partition all m,tz
integer,parameter:: max_no_cJstates = 1009   ! maximum no of J states
integer,parameter:: max_no_cTstates = 200    ! maximum no of T  states for each J
integer,parameter:: max_code = 72         ! maximum no of alphaa numeric characters 0-9,a-z
integer,parameter:: max_proj_mvec = 200000  ! maximum no of m-vectors per partition in NuProj
integer,parameter:: oxbash_lw =nwords        ! oxbash long word dimension
integer,parameter:: max_threads=1000                ! maximum no of threads

end module Parameters

module ProjParameters

use Parameters

implicit none

integer,parameter:: rks = 8                      ! real kind for orthogonalisation 8 or 16
integer,parameter:: rk = 8                       ! real kind for projection routines and storage
integer,parameter:: rkw = 4                      ! real kind for angular momentum projection writing 4 or 8
integer,parameter:: rko = 8                      ! real kind for orthog writing 4 or 8
integer,parameter:: max_nJT=10000                ! maximum value of nJT
integer,parameter:: max_lanc=1000                ! maximum dimension for Lanczos (Jmax-Jmin)+2
integer,parameter:: max_tqli=2000                ! maximum iterations for tqli
integer,parameter:: dim_par=3000                 ! use parallel above this dimension
end module ProjParameters

module OperParameters

implicit none

integer,parameter:: rc = 8                  !real kind for CG
integer,parameter:: ri = 8                  !real kind for interaction
real(kind=ri),parameter:: op_limit=real(.000001,ri)   !me below this value are ignored 
integer,parameter:: max_int =100000           !maximum enties in interaction file
logical,parameter:: npint=.false.            !use oxbash pn interaction files for pn, 
                                            !otherwise when false use isospin int for pn
integer,parameter:: max_no=10               !max number of interactions to add


end module

module MatrixParameters

integer,parameter:: rm=4                            !real kind for matrix
integer,parameter:: max_total_JT  = 100000             !max total jscheme dimension
integer,parameter:: max_op  = 4000000                 !max total mscheme t.b. me's
integer,parameter:: max_index =100000               !max total me's pert partition
integer,parameter:: np_cntl  =1                   ! =1 npr >= np; -1 npr <= np; 0 npr=all
logical,parameter:: hpart=.false.                   ! use positive h partitions  (set np_cntl /= 0)

end module MatrixParameters

module LanczosCommonParameters

implicit none

integer,parameter:: rl=4                                    ! real kind for lanczos
real(kind=rl),parameter:: level_limit=real( .0005,rl)       ! converged limit
real(kind=rl),parameter:: level_limitJ=real( .0005,rl)       ! converged limit jennings
real(kind=rl),parameter:: restart_limit= real(.0000001,rl)  ! restart  limit
integer,parameter:: max_sgmt=100000                         ! maximum no matrix  segments
real(kind=rl),parameter::  min_matrix=real(.000001,rl)      ! matrix element below this limit set to zero
real(kind=rl),parameter::  differ_limit = real(.002,rl)     ! allowable difference in <eig|H|eig> and alpha

end module LanczosCommonParameters

module LanczosParameters

use LanczosCommonParameters

implicit none

integer,parameter:: pre_iter=2                             ! no of pre-iterations
integer,parameter:: max_iter=390                            ! maximumm no of iterations
integer,parameter:: first_iterT=80                           ! first TRL try for iterations
integer,parameter:: first_iter=90                           ! first try for iterations
integer,parameter:: first_iterJ=1                          ! first try for iterations jennings
integer,parameter:: inc_iter=3                             ! increment iterations

end module LanczosParameters

module BlockParameters

use LanczosCommonParameters

implicit none

integer,parameter:: pre_iter=0                             ! no of pre-iterations blocks
integer,parameter:: max_iter=72                            ! maximumm no of iterations
integer,parameter:: first_iter=18                           ! first try for iterations
integer,parameter:: inc_iter=3                             ! increment iterations blocks
integer,parameter:: def_block=9   				                ! blocksize for NuBlkLncz

end module BlockParameters

module MvecParameters

implicit none

integer,parameter:: max_no_level=50            ! maximum no_level

end module MvecParameters

module OrderParameters

implicit none

integer,parameter:: ro = 8                      !real kind for matrix elements
integer,parameter:: max_part_tree = 2000000    !maximum dimension for ordering tree
integer,parameter:: max_node=1000000            !max nodes for binary lookuo 
integer,parameter:: max_me=10000               ! for NuMe 

end module OrderParameters

module MbmeParameters

integer,parameter:: rm=4                            !real kind for mbme
integer,parameter:: max_indexm =1000                !max total me's per dnr
integer,parameter:: max_total_mvr=6000000           ! max no of rhs mvectors

end module MbmeParameters

module MoverParameters

use LanczosCommonParameters
use MvecParameters

implicit none

integer,parameter:: max_level=max_no_level        !copy fro Mvec Parameters

end module MoverParameters

module TrampParameters

end module TrampParameters

module ClusterParameters

integer,parameter:: rm=4                            !real kind for mbme
integer,parameter:: max_indexm =10000                !max total me's per dnr
integer,parameter:: max_total_mvr=6000000           ! max no of rhs mvectors
integer,parameter:: max_total_mvc=2000000           ! max no of cluster mvectors

end module ClusterParameters

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!   Intel
    include 'ICommon.FI'
!   end Intel
!   Sun
!    include 'SCommon.FI'
!   end Sun
! --------------------------------------------------------------------- !

module Partition

! Routines for partition manipulations

use Parameters

implicit none

integer,parameter:: ptdim=max(8,(((max_no_shells-1)/nbytes)+1)*nbytes)
integer,private:: use_dim

type spartition
    integer(kind=2),dimension(ptdim):: shell =0
end type spartition

interface operator(==)
    module procedure ptn_lequals_int
    module procedure ptn_lequals_ptn
end interface

interface operator(/=)
    module procedure ptn_nequal_int
    module procedure ptn_nequal_ptn
end interface

interface assignment (=)
    module procedure ptn_aequals_array
    module procedure array_aequals_ptn
    module procedure ptn_aequals_int
end interface

interface operator (<)
    module procedure ptn_lt_int
    module procedure ptn_lt_ptn
    module procedure ptn_lt_array
end interface

interface operator (>)
    module procedure ptn_gt_int
    module procedure ptn_gt_ptn
    module procedure ptn_gt_array
end interface

interface operator (<=)
    module procedure ptn_lteq_int
    module procedure ptn_lteq_ptn
end interface

interface operator (>=)
    module procedure ptn_gteq_int
    module procedure ptn_gteq_ptn
end interface

interface operator (+)
    module procedure ptn_add
end interface

interface operator (-)
    module procedure ptn_sub
    module procedure ptn_minus
end interface

interface operator (*)
    module procedure array_x_ptn
    module procedure rarray_x_ptn
end interface 

interface operator (**)
    module procedure array_xx_ptn
end interface 

contains
    pure function gptdim()
        implicit none
        integer:: gptdim
        
        gptdim=use_dim

    end function gptdim
            
    subroutine set_ptdim(i)
        implicit none
        integer,intent(in):: i
        
        use_dim  = i
        
    end subroutine set_ptdim
    
    pure subroutine ptn_aequals_array(part,array)
        implicit none
        type(spartition),intent(out):: part
        integer,dimension(ptdim),intent(in):: array
        
        part%shell(1:use_dim)=array(1:use_dim)
        
    end subroutine ptn_aequals_array

    pure subroutine ptn_aequals_int(part,i)
        implicit none
        type(spartition),intent(out):: part
        integer,intent(in):: i
        
        part%shell(1:use_dim)=i
        
    end subroutine ptn_aequals_int

    pure subroutine array_aequals_ptn(array,part)
        implicit none
        type(spartition),intent(in):: part
        integer,dimension(ptdim),intent(out):: array
        
        array(1:use_dim)=part%shell(1:use_dim)
        
    end subroutine array_aequals_ptn

    pure function ptn_lequals_int(part,i)
        implicit none
        logical:: ptn_lequals_int
        type(spartition),intent(in)::part
        integer,intent(in):: i
        integer:: one
       
        do one=1,use_dim
        if (part%shell(one) /= i) then
            ptn_lequals_int=.false.
            return
        end if
        end do
        ptn_lequals_int=.true.
        
    end function ptn_lequals_int

    pure function ptn_lequals_ptn(parta,partb)
        implicit none
        logical:: ptn_lequals_ptn
        type(spartition),intent(in)::parta,partb
        integer:: one
        
        do one=1,use_dim
        if (parta%shell(one) /= partb%shell(one)) then
            ptn_lequals_ptn=.false.
            return
        end if
        end do
        ptn_lequals_ptn=.true.
        
    end function ptn_lequals_ptn

    pure function ptn_nequal_int(part,i)
        implicit none
        logical:: ptn_nequal_int
        type(spartition),intent(in)::part
        integer,intent(in):: i
       
        ptn_nequal_int=.not.ptn_lequals_int(part,i)
                
    end function ptn_nequal_int

    pure function ptn_nequal_ptn(parta,partb)
        implicit none
        logical:: ptn_nequal_ptn
        type(spartition),intent(in)::parta,partb
        
        ptn_nequal_ptn=.not.ptn_lequals_ptn(parta,partb)
        
    end function ptn_nequal_ptn
    
    pure function hop(part)
        implicit none
        logical:: hop
        type(spartition),intent(in):: part
        integer:: one
        
        do one=use_dim,1,-1
        if (part%shell(one) > 0) then
            hop=.true.
            return
        end if
        if (part%shell(one) < 0) then
            hop=.false.
            return
        end if
        end do
        hop=.true.
    
    end function hop

    
    pure function ptn_lt_int(part,i)
        implicit none
        logical:: ptn_lt_int
        type(spartition),intent(in):: part
        integer,intent(in):: i
        integer:: one
        
        do one=use_dim,1,-1
        if (part%shell(one) < i) then
            ptn_lt_int=.true.
            return
        end if
        end do
        ptn_lt_int=.false.
    
    end function ptn_lt_int

    pure function ptn_lteq_int(part,i)
        implicit none
        logical:: ptn_lteq_int
        type(spartition),intent(in):: part
        integer,intent(in):: i
        
        ptn_lteq_int=.not.(part > i)
    
    end function ptn_lteq_int

    pure function ptn_gteq_int(part,i)
        implicit none
        logical:: ptn_gteq_int
        type(spartition),intent(in):: part
        integer,intent(in):: i
        
        ptn_gteq_int=.not.(part < i)
    
    end function ptn_gteq_int

    pure function ptn_lt_ptn(parta,partb)
        implicit none
        logical:: ptn_lt_ptn
        type(spartition),intent(in):: parta,partb
        integer:: i
        do i=use_dim,1,-1
            if (parta%shell(i) > partb%shell(i)) then
                ptn_lt_ptn=.false.
                return
            else if (parta%shell(i) < partb%shell(i)) then
                    ptn_lt_ptn=.true.
                    return
            end if
        end do
        
        ptn_lt_ptn=.false.
        
    end function ptn_lt_ptn

    pure function ptn_lteq_ptn(parta,partb)
        implicit none
        logical:: ptn_lteq_ptn
        type(spartition),intent(in):: parta,partb

        ptn_lteq_ptn=.not.(parta > partb)
                
    end function ptn_lteq_ptn

    pure function ptn_gteq_ptn(parta,partb)
        implicit none
        logical:: ptn_gteq_ptn
        type(spartition),intent(in):: parta,partb

        ptn_gteq_ptn=.not.(parta < partb)
                
    end function ptn_gteq_ptn

    pure function ptn_gt_int(part,i)
        implicit none
        logical:: ptn_gt_int
        type(spartition),intent(in):: part
        integer,intent(in):: i
        integer:: one
        
        do one=use_dim,1,-1
        if (part%shell(one) > i) then
            ptn_gt_int=.true.
            return
        end if
        end do
        ptn_gt_int=.false.
    
    end function ptn_gt_int

    pure function ptn_gt_ptn(parta,partb)
        implicit none
        logical:: ptn_gt_ptn
        type(spartition),intent(in):: parta,partb
        integer:: i
        do i=use_dim,1,-1
            if (parta%shell(i) < partb%shell(i)) then
                ptn_gt_ptn=.false.
                return
            else if (parta%shell(i) > partb%shell(i)) then
                    ptn_gt_ptn=.true.
                    return
            end if
        end do
        
        ptn_gt_ptn=.false.
        
    end function ptn_gt_ptn
    
    pure function ptn_gt_array(parta,array)
        implicit none
        logical:: ptn_gt_array
        type(spartition),intent(in):: parta
        
        integer,dimension(max_no_shells),intent(in):: array
        integer:: i
        do i=use_dim,1,-1

            if (parta%shell(i) > array(i)) then
                    ptn_gt_array=.true.
                    return
            end if
        end do
        
        ptn_gt_array=.false.
        
    end function ptn_gt_array
    
    pure function ptn_lt_array(parta,array)
        implicit none
        logical:: ptn_lt_array
        type(spartition),intent(in):: parta
        
        integer,dimension(max_no_shells),intent(in):: array
        integer:: i
        do i=use_dim,1,-1

            if (parta%shell(i) < array(i)) then
                    ptn_lt_array=.true.
                    return
            end if
        end do
        
        ptn_lt_array=.false.
        
    end function ptn_lt_array
    
    pure function ptn_add(parta,partb)
        implicit none
        type(spartition):: ptn_add
        type(spartition),intent(in):: parta,partb
        
        ptn_add%shell(1:use_dim)=parta%shell(1:use_dim)+partb%shell(1:use_dim)
        
    end function ptn_add

    pure function ptn_sub(parta,partb)
        implicit none
        type(spartition):: ptn_sub
        type(spartition),intent(in):: parta,partb
        
        ptn_sub%shell(1:use_dim)=parta%shell(1:use_dim)-partb%shell(1:use_dim)
        
    end function ptn_sub

    pure function ptn_minus(parta)
        implicit none
        type(spartition):: ptn_minus
        type(spartition),intent(in):: parta
        
        ptn_minus%shell(1:use_dim)=-parta%shell(1:use_dim)
        
    end function ptn_minus
    
    pure function psum(part)
        implicit none
        integer:: psum
        type(spartition),intent(in):: part
        integer:: i
        psum=0
        do i=1,use_dim
        psum=psum+part%shell(i)
        end do
        
    end function psum

    pure function apsum(part)
        implicit none
        integer:: apsum
        type(spartition),intent(in):: part
        integer:: i
        apsum=0
        do i=1,use_dim
        apsum=apsum+abs(part%shell(i))
        end do
        
    end function apsum

    pure function pprod(part)
        implicit none
        integer:: pprod
        type(spartition),intent(in):: part
        integer:: i
        pprod=1
        do i=1,use_dim
        if (part%shell(i) /=0) pprod=pprod*part%shell(i)
        end do
        
    end function pprod

    pure function pintegral(part)
        implicit none
        type(spartition):: pintegral
        type(spartition),intent(in):: part
        integer::i
        
        pintegral%shell(1)=part%shell(1)
        if (use_dim == 1) return
        do i=2,use_dim
        pintegral%shell(i)=part%shell(i)+pintegral%shell(i-1)
        end do
        
    end function pintegral
    
    pure function array_x_ptn(array,part)
        implicit none
        type(spartition):: array_x_ptn
        integer,dimension(ptdim),intent(in):: array
        type(spartition),intent(in):: part
        
        array_x_ptn%shell(1:use_dim)=array(1:use_dim)*part%shell(1:use_dim)
        
    end function array_x_ptn
       
    pure function array_xx_ptn(array,part)
        implicit none
        type(spartition):: array_xx_ptn
        integer,dimension(ptdim),intent(in):: array
        type(spartition),intent(in):: part
        
        array_xx_ptn%shell(1:use_dim)=array(1:use_dim)**part%shell(1:use_dim)
        
    end function array_xx_ptn
   
    pure function rarray_x_ptn(rarray,part)
        implicit none
        real,dimension(ptdim):: rarray_x_ptn
        real,dimension(ptdim),intent(in):: rarray
        type(spartition),intent(in):: part
        
        rarray_x_ptn(1:use_dim)=rarray(1:use_dim)*part%shell(1:use_dim)
        
   end function rarray_x_ptn
        
end module Partition

module OpPartition

! Routines for partition manipulations

use Parameters

implicit none

integer,parameter:: optdim=4
integer,private:: use_dim=optdim

type opartition
    integer(kind=2),dimension(optdim):: op =0
end type opartition

interface operator(==)
    module procedure optn_lequals_int
    module procedure optn_lequals_optn
end interface

interface operator(/=)
    module procedure optn_nequal_int
    module procedure optn_nequal_optn
end interface

interface assignment (=)
    module procedure optn_aequals_array
    module procedure array_aequals_optn
    module procedure optn_aequals_int
end interface

interface operator (<)
    module procedure optn_lt_int
    module procedure optn_lt_optn
    module procedure optn_lt_array
end interface

interface operator (>)
    module procedure optn_gt_int
    module procedure optn_gt_optn
    module procedure optn_gt_array
end interface

interface operator (<=)
    module procedure optn_lteq_int
    module procedure optn_lteq_optn
end interface

interface operator (>=)
    module procedure optn_gteq_int
    module procedure optn_gteq_optn
end interface

interface operator (+)
    module procedure optn_add
end interface

interface operator (-)
    module procedure optn_sub
    module procedure optn_minus
end interface

interface operator (*)
    module procedure array_x_optn
    module procedure rarray_x_optn
end interface 

contains
    pure function goptdim()
        implicit none
        integer:: goptdim
        
        goptdim=use_dim

    end function goptdim
            
    subroutine set_optdim(i)
        implicit none
        integer,intent(in):: i
        
        use_dim  = i
        
    end subroutine set_optdim
    
    pure subroutine optn_aequals_array(part,array)
        implicit none
        type(opartition),intent(out):: part
        integer,dimension(optdim),intent(in):: array
        
        part%op(1:use_dim)=array(1:use_dim)
        
    end subroutine optn_aequals_array

    pure subroutine optn_aequals_int(part,i)
        implicit none
        type(opartition),intent(out):: part
        integer,intent(in):: i
        
        part%op(1:use_dim)=i
        
    end subroutine optn_aequals_int

    pure subroutine array_aequals_optn(array,part)
        implicit none
        type(opartition),intent(in):: part
        integer,dimension(optdim),intent(out):: array
        
        array(1:use_dim)=part%op(1:use_dim)
        
    end subroutine array_aequals_optn

    pure function optn_lequals_int(part,i)
        implicit none
        logical:: optn_lequals_int
        type(opartition),intent(in)::part
        integer,intent(in):: i
        integer:: one
       
        do one=1,use_dim
        if (part%op(one) /= i) then
            optn_lequals_int=.false.
            return
        end if
        end do
        optn_lequals_int=.true.
        
    end function optn_lequals_int

    pure function optn_lequals_optn(parta,partb)
        implicit none
        logical:: optn_lequals_optn
        type(opartition),intent(in)::parta,partb
        integer:: one
        
        do one=1,use_dim
        if (parta%op(one) /= partb%op(one)) then
            optn_lequals_optn=.false.
            return
        end if
        end do
        optn_lequals_optn=.true.
        
    end function optn_lequals_optn

    pure function optn_nequal_int(part,i)
        implicit none
        logical:: optn_nequal_int
        type(opartition),intent(in)::part
        integer,intent(in):: i
       
        optn_nequal_int=.not.optn_lequals_int(part,i)
                
    end function optn_nequal_int

    pure function optn_nequal_optn(parta,partb)
        implicit none
        logical:: optn_nequal_optn
        type(opartition),intent(in)::parta,partb
        
        optn_nequal_optn=.not.optn_lequals_optn(parta,partb)
        
    end function optn_nequal_optn
    
    pure function optn_lt_int(part,i)
        implicit none
        logical:: optn_lt_int
        type(opartition),intent(in):: part
        integer,intent(in):: i
        integer:: one
        
        do one=1,use_dim
        if (part%op(one) > i) then
            optn_lt_int=.false.
            return
        end if
        end do
        optn_lt_int=.true.
    
    end function optn_lt_int

    pure function optn_lteq_int(part,i)
        implicit none
        logical:: optn_lteq_int
        type(opartition),intent(in):: part
        integer,intent(in):: i
        
        optn_lteq_int=.not.(part > i)
    
    end function optn_lteq_int

    pure function optn_gteq_int(part,i)
        implicit none
        logical:: optn_gteq_int
        type(opartition),intent(in):: part
        integer,intent(in):: i
        
        optn_gteq_int=.not.(part < i)
    
    end function optn_gteq_int

    pure function optn_lt_optn(parta,partb)
        implicit none
        logical:: optn_lt_optn
        type(opartition),intent(in):: parta,partb
        integer:: i
        do i=use_dim,1,-1
            if (parta%op(i) > partb%op(i)) then
                optn_lt_optn=.false.
                return
            else if (parta%op(i) < partb%op(i)) then
                    optn_lt_optn=.true.
                    return
            end if
        end do
        
        optn_lt_optn=.false.
        
    end function optn_lt_optn

    pure function optn_lteq_optn(parta,partb)
        implicit none
        logical:: optn_lteq_optn
        type(opartition),intent(in):: parta,partb

        optn_lteq_optn=.not.(parta > partb)
                
    end function optn_lteq_optn

    pure function optn_gteq_optn(parta,partb)
        implicit none
        logical:: optn_gteq_optn
        type(opartition),intent(in):: parta,partb

        optn_gteq_optn=.not.(parta < partb)
                
    end function optn_gteq_optn

    pure function optn_gt_int(part,i)
        implicit none
        logical:: optn_gt_int
        type(opartition),intent(in):: part
        integer,intent(in):: i
        integer:: one
        
        do one=1,use_dim
        if (part%op(one) < i) then
            optn_gt_int=.false.
            return
        end if
        end do
        optn_gt_int=.true.
    
    end function optn_gt_int

    pure function optn_gt_optn(parta,partb)
        implicit none
        logical:: optn_gt_optn
        type(opartition),intent(in):: parta,partb
        integer:: i
        do i=use_dim,1,-1
            if (parta%op(i) < partb%op(i)) then
                optn_gt_optn=.false.
                return
            else if (parta%op(i) > partb%op(i)) then
                    optn_gt_optn=.true.
                    return
            end if
        end do
        
        optn_gt_optn=.false.
        
    end function optn_gt_optn
    
    pure function optn_gt_array(parta,array)
        implicit none
        logical:: optn_gt_array
        type(opartition),intent(in):: parta
        
        integer,dimension(optdim),intent(in):: array
        integer:: i
        do i=use_dim,1,-1

            if (parta%op(i) > array(i)) then
                    optn_gt_array=.true.
                    return
            end if
        end do
        
        optn_gt_array=.false.
        
    end function optn_gt_array
    
    pure function optn_lt_array(parta,array)
        implicit none
        logical:: optn_lt_array
        type(opartition),intent(in):: parta
        
        integer,dimension(optdim),intent(in):: array
        integer:: i
        do i=use_dim,1,-1

            if (parta%op(i) < array(i)) then
                    optn_lt_array=.true.
                    return
            end if
        end do
        
        optn_lt_array=.false.
        
    end function optn_lt_array
    
    pure function optn_add(parta,partb)
        implicit none
        type(opartition):: optn_add
        type(opartition),intent(in):: parta,partb
        
        optn_add%op(1:use_dim)=parta%op(1:use_dim)+partb%op(1:use_dim)
        
    end function optn_add

    pure function optn_sub(parta,partb)
        implicit none
        type(opartition):: optn_sub
        type(opartition),intent(in):: parta,partb
        
        optn_sub%op(1:use_dim)=parta%op(1:use_dim)-partb%op(1:use_dim)
        
    end function optn_sub

    pure function optn_minus(parta)
        implicit none
        type(opartition):: optn_minus
        type(opartition),intent(in):: parta
        
        optn_minus%op(1:use_dim)=-parta%op(1:use_dim)
        
    end function optn_minus
    
    pure function opsum(part)
        implicit none
        integer:: opsum
        type(opartition),intent(in):: part
        integer:: i
        opsum=0
        do i=1,use_dim
        opsum=opsum+part%op(i)
        end do
        
    end function opsum

    pure function opintegral(part)
        implicit none
        type(opartition):: opintegral
        type(opartition),intent(in):: part
        integer::i
        
        opintegral%op(1)=part%op(1)
        if (use_dim == 1) return
        do i=2,use_dim
        opintegral%op(i)=part%op(i)+opintegral%op(i-1)
        end do
        
    end function opintegral
    
    pure function array_x_optn(array,part)
        implicit none
        type(opartition):: array_x_optn
        integer,dimension(optdim),intent(in):: array
        type(opartition),intent(in):: part
        
        array_x_optn%op(1:use_dim)=array(1:use_dim)*part%op(1:use_dim)
        
    end function array_x_optn
   
    pure function rarray_x_optn(rarray,part)
        implicit none
        real,dimension(optdim):: rarray_x_optn
        real,dimension(optdim),intent(in):: rarray
        type(opartition),intent(in):: part
        
        rarray_x_optn(1:use_dim)=rarray(1:use_dim)*part%op(1:use_dim)
        
    end function rarray_x_optn
        
end module OpPartition

module Unsigned

! Routines to emulate unsigned integers

use Machine
use Parameters

implicit none

interface operator (.ult.)
    module procedure unsigned_lt
end interface

interface operator (.ugt.)
    module procedure unsigned_gt
end interface

interface operator (.ule.)
    module procedure unsigned_le
end interface

interface operator (.uge.)
    module procedure unsigned_ge
end interface

contains

    pure function unsigned_lt(i,j)
        implicit none
        logical:: unsigned_lt
        integer(kind=ik),intent(in):: i,j
        
        if (i >= 0) then
            if (j >= 0) then
                unsigned_lt=(i < j)
                return
            else
                unsigned_lt=.true.
                return
            end if
        else
            if (j >= 0) then
                unsigned_lt=.false.
                return
            else
                unsigned_lt=(i < j)
                return
            end if
        end if
        
        unsigned_lt=.false.
        
    end function unsigned_lt
        
    function unsigned_gt(i,j)
        implicit none
        logical:: unsigned_gt
        integer(kind=ik),intent(in):: i,j
        
        unsigned_gt=unsigned_lt(j,i)
        
    end function unsigned_gt
    
    function unsigned_le(i,j)
        implicit none
        logical:: unsigned_le
        integer(kind=ik),intent(in):: i,j
        
        unsigned_le=.not.unsigned_lt(j,i)
        
    end function unsigned_le
        
    function unsigned_ge(i,j)
        implicit none
        logical:: unsigned_ge
        integer(kind=ik),intent(in):: i,j
        
        unsigned_ge=.not.unsigned_lt(i,j)
        
    end function unsigned_ge
    
end module Unsigned        

module Mvector
! Routines for m-vector manipulations

use Unsigned

implicit none

integer:: use_words
integer(kind=8),parameter:: PAIRSEL=Z'AAAAAAAAAAAAAAAA'
integer(kind=8),parameter:: QUADSEL=Z'CCCCCCCCCCCCCCCC'
integer(kind=8),parameter:: OCTSEL=Z'F0F0F0F0F0F0F0F0'
integer(kind=8),parameter:: HEXSEL=Z'FF00FF00FF00FF00'
integer(kind=8),parameter:: WORDEL=Z'FFFF0000FFFF0000'
integer(kind=8),parameter:: LONGEL=Z'FFFFFFFF00000000'
integer(kind=8),parameter:: PAIRSOL=Z'5555555555555555'
integer(kind=8),parameter:: QUADSOL=Z'3333333333333333'
integer(kind=8),parameter:: OCTSOL=Z'0F0F0F0F0F0F0F0F'
integer(kind=8),parameter:: HEXSOL=Z'00FF00FF00FF00FF'
integer(kind=8),parameter:: WORDOL=Z'0000FFFF0000FFFF'
integer(kind=8),parameter:: LONGOL=Z'00000000FFFFFFFF'
integer(kind=4),parameter:: PAIRSE=Z'AAAAAAAA'
integer(kind=4),parameter:: QUADSE=Z'CCCCCCCC'
integer(kind=4),parameter:: OCTSE=Z'F0F0F0F0'
integer(kind=4),parameter:: HEXSE=Z'FF00FF00'
integer(kind=4),parameter:: WORDE=Z'FFFF0000'
integer(kind=4),parameter:: PAIRSO=Z'55555555'
integer(kind=4),parameter:: QUADSO=Z'33333333'
integer(kind=4),parameter:: OCTSO=Z'0F0F0F0F'
integer(kind=4),parameter:: HEXSO=Z'00FF00FF'
integer(kind=4),parameter:: WORDO=Z'0000FFFF'

type m_vector
    integer(kind=ik),dimension(nwords):: fn =0
end type m_vector

interface countb
    module procedure count32,count64
end interface countb

interface next_bit
    module procedure next_bit32,next_bit64
end interface next_bit

interface operator (==)
    module procedure mvec_lequals_mvec
    module procedure mvec_lequals_int
end interface

interface operator (/=)
    module procedure mvec_nequals_mvec
    module procedure mvec_nequals_int
end interface

interface assignment (=)
    module procedure mvec_aequals_array
    module procedure array_aequals_mvec
    module procedure mvec_aequals_int
end interface

interface operator (+)
    module procedure create
    module procedure qwick_create
end interface

interface operator (-)
    module procedure destroy
    module procedure qwick_destroy
end interface

interface operator (<)
    module procedure mvec_lt_mvec
end interface

interface operator (>)
    module procedure mvec_gt_mvec
end interface

interface operator (<=)
    module procedure mvec_le_mvec
end interface

interface operator (>=)
    module procedure mvec_ge_mvec
end interface

interface operator (.sub.)
    module procedure mvec_subset_mvec
end interface

interface operator (.and.)
    module procedure mvec_test
end interface

interface operator (.or.)
    module procedure creator
end interface

contains

    pure function gwords()
        implicit none
        integer:: gwords
        
        gwords=use_words
        
    end function gwords
    
    subroutine set_words(i)
        implicit none
        integer,intent(in):: i
        
        use_words=i
        
    end subroutine set_words
    
	pure function count32(iword)
		implicit none
		
		integer:: count32,cnt
		integer(kind=4),intent(in):: iword
		
		cnt=iword
		cnt=ISHFT(IAND(cnt,PAIRSE),-1)+IAND(cnt,PAIRSO)
		cnt=ISHFT(IAND(cnt,QUADSE),-2)+IAND(cnt,QUADSO)
		cnt=ISHFT(IAND(cnt,OCTSE), -4)+IAND(cnt,OCTSO)
		cnt=ISHFT(IAND(cnt,HEXSE), -8)+IAND(cnt,HEXSO)
		cnt=ISHFT(IAND(cnt,WORDE),-16)+IAND(cnt,WORDO)

		count32=cnt
		
	end function count32
	
	pure function next_bit32(iwordin)
	    implicit none
	    
	    integer:: next_bit32
		integer(kind=4),intent(in):: iwordin
		integer(kind=4)::iword
		
        iword=iwordin
		if (IAND(WORDO,iword)/=0) then
		    iword=IAND(WORDO,iword)
		    next_bit32=0
		else 
		    iword=IAND(WORDE,iword)
		    next_bit32=16
		end if
		if (IAND(HEXSO,iword)/=0) then
		    iword=IAND(HEXSO,iword)
		else 
		    iword=IAND(HEXSE,iword)
		    next_bit32=next_bit32+8
		end if
		if (IAND(OCTSO,iword)/=0) then
		    iword=IAND(OCTSO,iword)
		else 
		    iword=IAND(OCTSE,iword)
		    next_bit32=next_bit32+4
		end if
		if (IAND(QUADSO,iword)/=0) then
		    iword=IAND(QUADSO,iword)
		else 
		    iword=IAND(QUADSE,iword)
		    next_bit32=next_bit32+2
		end if
		if (IAND(PAIRSO,iword)==0) then
		    next_bit32=next_bit32+1
		end if
		next_bit32=next_bit32+1
	    return
	    
    end function next_bit32
		
	pure function next_bit64(iwordin)
	    implicit none
	    
	    integer:: next_bit64
		integer(kind=8),intent(in):: iwordin
		integer(kind=8)::iword
		
        iword=iwordin
		if (IAND(LONGOL,iword)/=0) then
		    iword=IAND(LONGOL,iword)
		    next_bit64=0
		else 
		    iword=IAND(LONGEL,iword)
		    next_bit64=32
		end if
		if (IAND(WORDOL,iword)/=0) then
		    iword=IAND(WORDOL,iword)
		else 
		    iword=IAND(WORDEL,iword)
		    next_bit64=next_bit64+16
		end if
		if (IAND(HEXSOL,iword)/=0) then
		    iword=IAND(HEXSOL,iword)
		else 
		    iword=IAND(HEXSEL,iword)
		    next_bit64=next_bit64+8
		end if
		if (IAND(OCTSOL,iword)/=0) then
		    iword=IAND(OCTSOL,iword)
		else 
		    iword=IAND(OCTSEL,iword)
		    next_bit64=next_bit64+4
		end if
		if (IAND(QUADSOL,iword)/=0) then
		    iword=IAND(QUADSOL,iword)
		else 
		    iword=IAND(QUADSEL,iword)
		    next_bit64=next_bit64+2
		end if
		if (IAND(PAIRSOL,iword)==0) then
		    next_bit64=next_bit64+1
		end if
		next_bit64=next_bit64+1
	    return
	    
    end function next_bit64
		
	pure function count64(iword)
		implicit none
		integer:: count64
		integer(kind=8),intent(in):: iword
        integer(kind=8):: cnt
        			
		cnt=iword
		cnt=ISHFT(IAND(cnt,PAIRSEL),-1)+IAND(cnt,PAIRSOL)
		cnt=ISHFT(IAND(cnt,QUADSEL),-2)+IAND(cnt,QUADSOL)
		cnt=ISHFT(IAND(cnt,OCTSEL), -4)+IAND(cnt,OCTSOL)
		cnt=ISHFT(IAND(cnt,HEXSEL), -8)+IAND(cnt,HEXSOL)
		cnt=ISHFT(IAND(cnt,WORDEL),-16)+IAND(cnt,WORDOL)
		cnt=ISHFT(IAND(cnt,LONGEL),-32)+IAND(cnt,LONGOL)

		count64=cnt

    end function count64	

    pure function mvec_lequals_mvec(veca,vecb)
        implicit none
        logical:: mvec_lequals_mvec
        type(m_vector),intent(in):: veca,vecb
        integer:: one
        
        do one=1,use_words
        if (veca%fn(one) /= vecb%fn(one)) then    
            mvec_lequals_mvec=.false.
            return
        end if
        end do
        mvec_lequals_mvec=.true.

    end function mvec_lequals_mvec

    pure function mvec_nequals_mvec(veca,vecb)
        implicit none
        logical:: mvec_nequals_mvec
        type(m_vector),intent(in):: veca,vecb
        
        mvec_nequals_mvec=.not.(mvec_lequals_mvec(veca,vecb))
        
    end function mvec_nequals_mvec
        
    pure function mvec_lequals_int(veca,i)
        implicit none
        logical:: mvec_lequals_int
        type(m_vector),intent(in):: veca
        integer(kind=ik),intent(in):: i
        integer:: one
        
        do one=1,use_words
        if (veca%fn(one) /= i) then    
            mvec_lequals_int=.false.
            return
        end if
        end do
        mvec_lequals_int=.true.
        
    end function mvec_lequals_int

    pure function mvec_nequals_int(veca,i)
        implicit none
        logical:: mvec_nequals_int
        type(m_vector),intent(in):: veca
        integer(kind=ik),intent(in):: i
        
        mvec_nequals_int=.not.(mvec_lequals_int(veca,i))
        
    end function mvec_nequals_int

	pure subroutine array_aequals_mvec(array,vec)
		implicit none
		type(m_vector),intent(in):: vec
		integer,dimension(nbits*nwords),intent(out):: array
		integer:: n,iw,ib
    		n=0
    		array=0
			do iw=1,use_words
				do ib=0,nbits-1
					n=n+1
					if  (btest(vec%fn(iw),ib)) then
						array(n)=1
					else	
						array(n)=0
					end if
				end do
	    	end do
	    	
	end subroutine array_aequals_mvec
				
	pure subroutine mvec_aequals_array(vec,array)
		implicit none
		type(m_vector),intent(out):: vec
		integer,intent(in),dimension(nbits*nwords):: array
		integer:: iw,ib,n
				 
            vec%fn(1:use_words)=0
    		n=0
			do iw=1,use_words
				do ib=0,nbits-1
					n=n+1
					if (array(n) == 1) vec%fn(iw)=ibset(vec%fn(iw),ib)	
				end do
	    	end do
				
	end subroutine mvec_aequals_array

	pure subroutine mvec_aequals_int(vec,i)
		implicit none
		type(m_vector),intent(out):: vec
		integer(kind=ik),intent(in):: i
				 
        vec%fn(1:use_words)=i
				
	end subroutine mvec_aequals_int
	
	subroutine mwrite(dev,mvec)
	    implicit none
	    integer,intent(in):: dev
	    type(m_vector),intent(in):: mvec
	    integer,dimension(nbits*nwords):: array
	    integer:: iw,ib,n
	    
	    array=mvec

	    n=0
	    write(dev,'(a1)',advance='no') ' '
	    do iw=1,use_words
	        do ib=0,nbits-1
	            n=n+1
	            write(dev,'(i1)',advance='no') array(n)
	            if ((mod(n,64) == 0).and.(n /= nbits*nwords)) write(dev,'(a1,/,a1)',advance='no') '&',' '
	        end do
	    end do
	    write(dev,*)
	    
	end subroutine mwrite
	
	pure function cNop(mvec)
	    implicit none
	    integer:: cNop
	    type(m_vector),intent(in):: mvec
	    integer:: i
	    integer:: wsum
	    wsum=0
	    do i=1,use_words
	        wsum=wsum+countb(mvec%fn(i))
	    end do
	    
	    cNop=wsum
	    
    end function cNop
    
    pure function create(a,mvec)
		implicit none
		type(m_vector),intent(in):: mvec
		integer,intent(in):: a
		integer:: iw,ib
		type(m_vector):: create

			iw=(a-1)/nbits + 1
			ib=a - (iw-1)*nbits
		
			if (btest(mvec%fn(iw),ib-1)) then
				create=int(0,ik)
				return
			end if
                
            create=mvec
			create%fn(iw)=ibset(create%fn(iw),ib-1)
			
	end function create

    pure function creator(a,mvec)
		implicit none
		type(m_vector),intent(in):: mvec
		integer,intent(in):: a
		integer:: iw,ib
		type(m_vector):: creator

			iw=(a-1)/nbits + 1
			ib=a - (iw-1)*nbits
		
			if (btest(mvec%fn(iw),ib-1)) then
				creator=mvec
				return
			end if
                
            creator=mvec
			creator%fn(iw)=ibset(creator%fn(iw),ib-1)
			
	end function creator
	
    pure function mvec_x_mvec(mveca,da1,da,mvecb,db1,db)
		implicit none
		type(m_vector),dimension(:),intent(in):: mveca,mvecb
		integer,intent(in)::da,db,da1,db1
		integer:: kva,iwb,ibb,jvb,iwa,iba,iva,jwb,jbb,kvb
		type(m_vector),dimension(2)::mvec_x_mvec
		
		do jvb=1,2
		    mvec_x_mvec(jvb)%fn(1:use_words)=0
		end do
		do jvb=1,db
        do kvb=1,db1
            iwb=(kvb-1)/nbits + 1
            ibb=kvb - (iwb-1)*nbits
		    if (btest(mvecb(jvb)%fn(iwb),ibb-1)) then
		    mvec_x_mvec(2)%fn(iwb)=ibset(mvec_x_mvec(2)%fn(iwb),ibb-1)
		    end if
	    end do
	    end do
	    do iva=1,da
            iwa=(iva-1)/nbits + 1
            iba=iva - (iwa-1)*nbits
	        if (mveca(iva)/=int(0,ik)) then
		    mvec_x_mvec(1)%fn(iwa)=ibset(mvec_x_mvec(1)%fn(iwa),iba-1)
		    end if
		end do
			
	end function mvec_x_mvec


    pure function qwick_create(a,mvec)
		implicit none
		type(m_vector),intent(in):: mvec
		integer,dimension(*),intent(in):: a
		type(m_vector):: qwick_create

            qwick_create=mvec
			qwick_create%fn(a(1))=ibset(qwick_create%fn(a(1)),a(2)-1)
			
	end function qwick_create

    pure function qwick_destroy(a,mvec)
		implicit none
		type(m_vector),intent(in):: mvec
		integer,dimension(*),intent(in):: a
		type(m_vector):: qwick_destroy

            qwick_destroy=mvec
			qwick_destroy%fn(a(1))=ibclr(qwick_destroy%fn(a(1)),a(2)-1)
			
	end function qwick_destroy

    pure function mvec_test(a,mvec)
		implicit none
		type(m_vector),intent(in):: mvec
		integer,intent(in):: a
		integer:: iw,ib
		logical:: mvec_test

			iw=(a-1)/nbits + 1
			ib=a - (iw-1)*nbits
		
			mvec_test=btest(mvec%fn(iw),ib-1)

	end function mvec_test

    pure function destroy(a,mvec)
		implicit none
		type(m_vector),intent(in):: mvec
		integer,intent(in):: a
		integer:: iw,ib
		type(m_vector):: destroy

			iw=(a-1)/nbits + 1
			ib=a - (iw-1)*nbits
		
			if (.not.(btest(mvec%fn(iw),ib-1))) then
				destroy=int(0,ik)
				return
			end if
                
            destroy=mvec
			destroy%fn(iw)=ibclr(destroy%fn(iw),ib-1)
			
	end function destroy

    pure function qphase(a,mvec)
		implicit none
		type(m_vector),intent(in):: mvec
		integer,dimension(2),intent(in):: a
		integer:: iw,ib,i
        real(kind=4):: qphase
        integer(kind=ik):: mfn
        integer:: wsum
        
            wsum=0
			iw=a(1)
			ib=a(2)
            if (iw > 1) then	
                do i=1,iw-1
	                wsum=wsum+countb(mvec%fn(i))
	            end do
	        end if
	        mfn=0
	        call mvbits(mvec%fn(iw),0,ib-1,mfn,0)
			
		    if (btest(countb(mfn)+wsum,0)) then
				qphase=-1.0
		    else
		        qphase=1.0
			end if
						
	end function qphase
	
    pure function mphase(iw,mvec)
		implicit none
		type(m_vector),intent(in):: mvec
		integer,intent(in):: iw
		integer:: i
        real(kind=4):: mphase
        integer:: wsum
        
            wsum=0
            if (iw > 1) then	
                do i=1,iw-1
	                wsum=wsum+countb(mvec%fn(i))
	            end do
			
		        if (btest(wsum,0)) then
				    mphase=-1.0
		        else
		            mphase=1.0
			    end if
			else
			    mphase=1.0
	        end if
						
	end function mphase
	
    pure function wphase(ib,mvecfn)
		implicit none
		integer(kind=ik),intent(in):: mvecfn
		integer,intent(in):: ib
        real(kind=4):: wphase
        integer(kind=ik):: mfn
        
	        mfn=0
	        call mvbits(mvecfn,0,ib-1,mfn,0)
			
		    if (btest(countb(mfn),0)) then
				wphase=-1.0
		    else
		        wphase=1.0
			end if
						
	end function wphase
	
    pure function gphase(a,mvec)
		implicit none
		type(m_vector),intent(in):: mvec
		integer,intent(in):: a
		integer:: iw,ib,i
        real(kind=4):: gphase
        integer(kind=ik):: mfn
        integer:: wsum
        
            wsum=0
			iw=(a-1)/nbits + 1
			ib=a - (iw-1)*nbits
            if (iw > 1) then	
                do i=1,iw-1
	                wsum=wsum+countb(mvec%fn(i))
	            end do
	        end if
	        mfn=0
	        call mvbits(mvec%fn(iw),0,ib-1,mfn,0)
			
		    if (btest(countb(mfn)+wsum,0)) then
				gphase=-1.0
		    else
		        gphase=1.0
			end if
						
	end function gphase
	
	pure function mvec_lt_mvec(mveca,mvecb)
	    implicit none
	    logical:: mvec_lt_mvec
	    type(m_vector),intent(in):: mveca,mvecb
	    integer:: i
	    do i=use_words,1,-1
	        if (unsigned_lt(mvecb%fn(i),mveca%fn(i))) then
	            mvec_lt_mvec=.false.
	            return
	        else if (unsigned_lt(mveca%fn(i),mvecb%fn(i))) then
	            mvec_lt_mvec=.true.
	            return
	        end if
	    end do
	    
	    mvec_lt_mvec=.false.
	    
    end function mvec_lt_mvec
    
    pure function mvec_le_mvec(mveca,mvecb)
	    implicit none
	    logical:: mvec_le_mvec
	    type(m_vector),intent(in):: mveca,mvecb
            
        mvec_le_mvec=.not.(mvec_gt_mvec(mveca,mvecb))
    
    end function mvec_le_mvec

    pure function mvec_ge_mvec(mveca,mvecb)
	    implicit none
	    logical:: mvec_ge_mvec
	    type(m_vector),intent(in):: mveca,mvecb
            
        mvec_ge_mvec=.not.(mvec_lt_mvec(mveca,mvecb))
    
    end function mvec_ge_mvec

	pure function mvec_gt_mvec(mveca,mvecb)
	    implicit none
	    logical:: mvec_gt_mvec
	    type(m_vector),intent(in):: mveca,mvecb
	    integer:: i
	    do i=use_words,1,-1
	        if (unsigned_lt(mveca%fn(i),mvecb%fn(i))) then
	            mvec_gt_mvec=.false.
	            return
	        else if (unsigned_lt(mvecb%fn(i),mveca%fn(i))) then
	            mvec_gt_mvec=.true.
	            return
	        end if
	    end do
	    
	    mvec_gt_mvec=.false.
	    
    end function mvec_gt_mvec
    
    pure function mvec_subset_mvec(mveca,mvecb)
        implicit none
        logical:: mvec_subset_mvec
        type(m_vector),intent(in):: mveca,mvecb
        integer,dimension(nwords):: one
        
        one=0
        where (iand(mveca%fn,mvecb%fn) == mveca%fn) one=1
        if (sum(one) == nwords) then
            mvec_subset_mvec=.true.
        else
            mvec_subset_mvec=.false.
        end if
        
    end function mvec_subset_mvec
    
    pure function mand(mveca,mvecb)
        implicit none
        type(m_vector):: mand
        type(m_vector),intent(in):: mveca,mvecb
        integer:: i
        
         do i=1,use_words
         mand%fn(i)=iand(mveca%fn(i),mvecb%fn(i))
         end do
        
    end function mand

    pure function mor(mveca,mvecb)
        implicit none
        type(m_vector):: mor
        type(m_vector),intent(in):: mveca,mvecb
        integer:: i
        
         do i=1,use_words
         mor%fn(i)=ior(mveca%fn(i),mvecb%fn(i))
         end do
        
    end function mor


    pure function meor(mveca,mvecb)
        implicit none
        type(m_vector):: meor
        type(m_vector),intent(in):: mveca,mvecb
        integer:: i
        
         do i=1,use_words
         meor%fn(i)=ieor(mveca%fn(i),mvecb%fn(i))
         end do
        
    end function meor
    
    pure function ph(mvec,word_mask)
        implicit none
        type(m_vector):: ph
        type(m_vector),intent(in):: mvec
        integer(kind=ik),dimension(nwords),intent(in):: word_mask
        integer::n
    
        do n=1,use_words
            ph%fn(n)=iand(word_mask(n),not(mvec%fn(n)))
        end do
    end function ph  
   
    pure function sph(mvec,shell_mask,shell_word)
        implicit none
        type(m_vector):: sph
        type(m_vector),intent(in):: mvec
        integer(kind=ik),intent(in):: shell_mask
        integer,intent(in)::shell_word
    
        sph%fn(shell_word)=iand(shell_mask,not(mvec%fn(shell_word)))
        
    end function sph  
   
    pure function phs(mvec,shell_mask,shell_word)
        implicit none
        type(m_vector):: phs
        type(m_vector),intent(in):: mvec
        integer(kind=ik),intent(in):: shell_mask
        integer,intent(in)::shell_word
        
        phs%fn(1:use_words)=mvec%fn(1:use_words)
        phs%fn(shell_word)=ieor(shell_mask,phs%fn(shell_word))
        
    end function phs  
   
    pure function mtest(i,n)
        implicit none
        integer,intent(in):: n
        logical:: mtest
        integer(kind=ik),intent(in):: i
   
        mtest=btest(i,n-1)
    
    end function mtest
   
end module Mvector

module Shells

use Partition
use Mvector

implicit none

integer:: no_shells,no_nucl_types,shn,wdn,tzn,nMS,ptn,ptni,max_sj2u
logical:: use_isospin
integer,dimension(max_no_shells):: cN,sn,sl,sj2,st2,sp
integer,dimension(ptdim):: cNc
character(len=1),dimension(ptdim):: nucleon
integer,dimension(max_no_shells):: shell_word,shell_start,shell_2jp1
integer,dimension(ptdim):: cNMS,max_MSI,max_MSN
integer(kind=ik),dimension(max_no_shells):: shell_mask,shell_masko,shell_maske,nshell_mask
integer(kind=ik),dimension(nwords):: word_mask
integer,dimension(0:1):: no_shells_t2
type(spartition),dimension(ptdim):: cMS_mask
data st2 /max_no_shells*0/

contains

    subroutine check_shell_order()
        implicit none
        integer,dimension(max_no_shells):: one
        
        no_nucl_types=0
        no_shells_t2=0
        use_isospin=.false.
        one=0
        if (st2(1) == 0) then
            no_shells_t2(1)=0
            no_shells_t2(0)=no_shells
            no_nucl_types=1
            use_isospin=.false.
            where(st2 == 0) one=1
            if (sum(one) /= max_no_shells) then
                print * , ' First shell has tz=0, subsequent shells have tz>0. Reorder shells.'
                call shell_message()
                stop
            end if
        end if
        if (st2(1) == 1) then
            no_nucl_types=1
            use_isospin=.true.
            do shn=2,no_shells
                if (st2(shn) /= 1) exit
            end do
            no_shells_t2(1)=shn-1
            no_shells_t2(0)=no_shells-shn+1
            if (shn-1 /= no_shells) no_nucl_types=2

            where(st2 == 0) one=1
            if (sum(one) /= max_no_shells-shn+1) then
                print *, ' First shell has tz=1/2. This is repeated up to shell',shn-1,'.'
                print *, ' Subsequent shells do not have only tz=0. Reorder shells.'
                call shell_message()
                stop
            end if
        end if

        if (no_nucl_types > nwords) then
            print *, ' The parameter NWORDS in module Parameters must at least equal no of nucleon types.'
            print *, ' Nucleon types have tz=0, tz=+1/2.'
            call shell_message()
            stop
        end if
                
   end subroutine check_shell_order
   
   subroutine shell_message
        implicit none

        print * , ' Shells must be ordered tz=1/2 first then tz=0.'
        
    end subroutine shell_message
   
    subroutine setup_major_shells()
        implicit none
        integer,dimension(max_no_shells):: cN_copy
        integer:: i
        
        cNc(1:no_shells)=cN(1:no_shells)
        cN_copy=cN
        if (no_shells < max_no_shells) then
            do i=no_shells+1,max_no_shells
                 cN_copy(i)=999
            end do
        end if       
        nMS=0
        do
            if (minval(cN_copy) == 999) exit
                nMS=nMS+1
                cNMS(nMS)=minval(cN_copy)
            where(cN_copy == cNMS(nMS)) cN_copy=999
        end do
        
        do i=1,nMS
            where( cNc == cNMS(i) ) cMS_mask(i)%shell = 1
        end do
    end subroutine setup_major_shells
        
    subroutine setup_shells(logic)
        implicit none
        logical,intent(in),optional::logic
        integer:: iw,ib
        
        iw=1
        ib=1
        shn=1
        shell_mask=0
        shell_masko=0
        shell_maske=0
        word_mask=0
        
        if (no_shells_t2(1) > 0) then
            do shn=shn,no_shells_t2(1)
                shell_start(shn)=ib
                shell_word(shn)=iw
                do ib=ib,ib+2*(sj2(shn)+1)-1
                    shell_mask(shn)=ibset(shell_mask(shn),ib-1)
                    if (btest(ib,0)) then
                        shell_masko(shn)=ibset(shell_masko(shn),ib-1)
                    else
                        shell_maske(shn)=ibset(shell_maske(shn),ib-1)
                    end if
                    word_mask(iw)=ibset(word_mask(iw),ib-1)
                end do
                if (shn+1  > no_shells) go to 50
                if (ib + 2*(sj2(shn+1)+1)-1 > nbits) then
                    iw=iw+1
                    if (iw > nwords) call word_message()
                    ib=1
                end if
            end do
            if (shn-1 < no_shells) then
            if (ib /= 1) then
                iw=iw+1
                if (iw > nwords) call word_message()
                ib=1
            end if
            end if
        end if
        
        if (no_shells_t2(0) > 0) then
            do shn=shn,no_shells_t2(0)+no_shells_t2(1)
                shell_start(shn)=ib
                shell_word(shn)=iw
                do ib=ib,ib+sj2(shn)
                    shell_mask(shn)=ibset(shell_mask(shn),ib-1)
                    word_mask(iw)=ibset(word_mask(iw),ib-1)
                end do
                if (shn+1  > no_shells) go to 50
                if (ib + sj2(shn+1) > nbits) then
                    iw=iw+1
                    if (iw > nwords) call word_message()
                    ib=1
                end if
            end do
            if (shn-1 < no_shells) then     
            if (ib /= 1) then
                iw=iw+1
                if (iw > nwords) call word_message()
                ib=1
            end if
            end if
        
        end if
50      do shn=1,no_shells
            nshell_mask(shn)=not(shell_mask(shn))
        end do
        
        if (present(logic) .and. logic ) print * ,' No of words used ',iw
        call set_words(iw)
        
    end subroutine setup_shells    
                
    subroutine word_message()
        implicit none
        
        print *, ' The parameter NWORDS may need to be increased in module Parameters.'    
        print *, ' Different nucleon types are stored in different words.'
        call shell_message()
        print *, ' Shells are not spread over two words.'
        print *, ' Reordering shells to make full use of each word might help.'
        stop
        
    end subroutine word_message

    subroutine swrite(dev,mvec)
	    implicit none
	    integer,intent(in):: dev
	    type(m_vector),intent(in):: mvec
	    integer,dimension(nbits*nwords):: array
	    integer:: sm2,stz2,ib,n,m
	    
	    array=mvec
	    m=0
	    write(dev,'(a1)',advance='no') ' '
	    do shn=1,no_shells
	       ib=shell_start(shn)+(shell_word(shn)-1)*nbits
	       n=ib
	       do sm2=-sj2(shn),sj2(shn),2 
	       do stz2=-st2(shn),st2(shn),2

	            m=m+1
	            write(dev,'(i1)',advance='no') array(n)
	            n=n+1
	            if ((mod(m,64) == 0).and.(n /= nbits*nwords)) write(dev,'(a1,/,a1)',advance='no') '&',' '
	       end do
	       end do
	       write(dev,'(a1)',advance='no') ':'
	       m=m+1
	       if ((mod(m,64) == 0).and.(n /= nbits*nwords)) write(dev,'(a1,/,a1)',advance='no') '&',' '
	    end do
	    write(dev,*)
	    
	end subroutine swrite
	
	function shell_offset(shn)
	implicit none
	integer:: shell_offset,shn
	
	shell_offset=shell_start(shn) +(shell_word(shn)-1)*nbits
	
	end function shell_offset
	    
end module Shells

module Finds

use Mvector
use OpPartition
use Partition

interface bfind
    module procedure ifind,mfind,pfind,pifind,sfind,opfind
end interface

contains

    pure function ifind(iarray,n,i,offset)
        implicit none
        integer,intent(in):: n
        integer(kind=ik),intent(in):: i
        integer,optional,intent(in):: offset
        integer:: ifind
        integer(kind=ik),dimension(*),intent(in):: iarray
        integer:: in,start
        real:: rn,drn,sign
        
        if (present(offset)) then
            start=offset
        else
            start=0
        end if
        
        if (iarray(n+start) > iarray(1+start)) then
            sign=1.0
        else
            sign=-1.0
        end if
        
        if (i == iarray(n+start)) then 
            ifind=n+start
            return
        end if
        rn=n
        rn=rn/2.0
        drn=rn
        
        do
            in=nint(rn)
            if (i == iarray(in+start)) then
                ifind=in+start
                exit
            end if
            drn=drn/2.0
            if (drn < 0.5) then
                ifind=0
                exit
            end if
            if (i > iarray(in+start)) then
                rn=rn+drn*sign
            else
                rn=rn-drn*sign
            end if
        end do
            
    end function ifind
                
    function pmfind(marray,n,mvec,offset)
        implicit none
        integer,intent(in):: n
        type(m_vector),intent(in):: mvec
        integer,optional,intent(in):: offset
        integer:: pmfind
        integer,dimension(max_threads)::findp,off,nn
        type(m_vector),dimension(:),intent(in):: marray
        integer:: in,start,thrn,nthr,ndth
        real:: rn,drn,sign
        integer:: OMP_get_max_threads,OMP_get_thread_num
        
        if (n*gwords()<100000) then
            if (present(offset)) then
                pmfind=mfind(marray,n,mvec,offset)
            else
                pmfind=mfind(marray,n,mvec)
            end if
            return
        else
!$OMP PARALLEL
!        nthr=OMP_get_max_threads()         !un comment this line for parallel
!$OMP END PARALLEL
        if (present(offset)) then
            start=offset
        else
            start=0
        end if
        if (nthr>max_threads) then
            print *, ' Increase max_threads in Parameters'
            stop
        end if
        findp(1:nthr)=0  
        ndth=n/nthr
        if (nthr*ndth<n) ndth=ndth+1
        do in=1,nthr
            off(in)=start+(in-1)*ndth
            nn(in)=ndth
            if (in==nthr) nn(in)=n-off(in)
        end do
!$OMP PARALLEL
!$OMP DO PRIVATE(in,thrn)
        do in=1,nthr
            findp(in)=mfind(marray,nn(in),mvec,off(in))
        end do        
!$OMP END DO
!$OMP END PARALLEL
        pmfind=maxval(findp(1:nthr))
        end if
            
    end function pmfind

    pure function mfind(marray,n,mvec,offset)
        implicit none
        integer,intent(in):: n
        type(m_vector),intent(in):: mvec
        integer,optional,intent(in):: offset
        integer:: mfind
        type(m_vector),dimension(:),intent(in):: marray
        integer:: in,start
        real:: rn,drn,sign
                
        if (present(offset)) then
            start=offset
        else
            start=0
        end if
        
        if (marray(n+start) > marray(1+start)) then
            sign=1.0
        else
            sign=-1.0
        end if
                
        if (mvec == marray(n+start)) then 
            mfind=n+start
            return
        end if
        rn=n
        rn=rn/2.0
        drn=rn*sign
        
        do
            in=nint(rn)+start
            if (mvec == marray(in)) then
                mfind=in
                exit
            end if
            drn=drn/2.0
            if (abs(drn) < 0.5) then
                mfind=0
                exit
            end if
            if (mvec > marray(in)) then
                rn=rn+drn
            else
                rn=rn-drn
            end if
        end do
            
    end function mfind
                
    pure function pfind(parray,n,part,offset)
        implicit none
        integer,intent(in):: n
        type(spartition),intent(in):: part
        integer,optional,intent(in):: offset
        integer:: pfind
        type(spartition),dimension(*),intent(in):: parray
        integer:: in,start
        real:: rn,drn,sign
                
        if (present(offset)) then
            start=offset
        else
            start=0
        end if
                
        if (parray(n+start) > parray(1+start)) then
            sign=1.0
        else
            sign=-1.0
        end if
                
        if (part == parray(n+start)) then 
            pfind=n+start
            return
        end if
        rn=n
        rn=rn/2.0
        drn=rn
        
        do
            in=nint(rn)
            if (part == parray(in+start)) then
                pfind=in+start
                exit
            end if
            drn=drn/2.0
            if (drn < 0.5) then
                pfind=0
                exit
            end if
            if (part > parray(in+start)) rn=rn+drn*sign
            if (part < parray(in+start)) rn=rn-drn*sign
        end do
            
    end function pfind
    
    pure function opfind(parray,n,part,offset)
        implicit none
        integer,intent(in):: n
        type(opartition),intent(in):: part
        integer,optional,intent(in):: offset
        integer:: opfind
        type(opartition),dimension(*),intent(in):: parray
        integer:: in,start
        real:: rn,drn,sign
                
        if (present(offset)) then
            start=offset
        else
            start=0
        end if
                
        if (parray(n+start) > parray(1+start)) then
            sign=1.0
        else
            sign=-1.0
        end if
                
        if (part == parray(n+start)) then 
            opfind=n+start
            return
        end if
        rn=n
        rn=rn/2.0
        drn=rn
        
        do
            in=nint(rn)
            if (part == parray(in+start)) then
                opfind=in+start
                exit
            end if
            drn=drn/2.0
            if (drn < 0.5) then
                opfind=0
                exit
            end if
            if (part > parray(in+start)) rn=rn+drn*sign
            if (part < parray(in+start)) rn=rn-drn*sign
        end do
            
    end function opfind
    
    pure function pifind(parray,indx,n,part)
        implicit none
        integer,intent(in):: n
        type(spartition),intent(in):: part
        integer:: pifind
        type(spartition),dimension(*),intent(in):: parray
        integer,dimension(*),intent(in):: indx
        integer:: in
        real:: rn,drn,sign
        
        if (parray(indx(n)) > parray(indx(1))) then
            sign=1.0
        else
            sign=-1.0
        end if
                
        if (part == parray(indx(n))) then 
            pifind=n
            return
        end if
        rn=n
        rn=rn/2.0
        drn=rn
        
        do
            in=nint(rn)
            if (part == parray(indx(in))) then
                pifind=in
                exit
            end if
            drn=drn/2.0
            if (drn < 0.5) then
                pifind=0
                exit
            end if
            if (part > parray(indx(in))) rn=rn+drn*sign
            if (part < parray(indx(in))) rn=rn-drn*sign
        end do
            
    end function pifind
    
    function sfind(marray,n,mvec,svec)
        implicit none
        integer,intent(in):: n
        type(m_vector),intent(in):: mvec,svec
        integer:: sfind
        type(m_vector),dimension(n),intent(in):: marray
        integer:: in
        real:: rn,drn,sign
        
        if (marray(n) > marray(1)) then
            sign=1.0
        else
            sign=-1.0
        end if
        
        if (svec.sub.marray(n)) then 
            sfind=n
            return
        end if
        rn=n
        rn=rn/2.0
        drn=rn
        
        do
            in=nint(rn)
            if (svec.sub.marray(in)) then
                sfind=in
                exit
            end if
            drn=drn/2.0
            if (drn < 0.5) then
                sfind=0
                exit
            end if
            if (mvec > marray(in)) rn=rn+drn*sign
            if (mvec < marray(in)) rn=rn-drn*sign
        end do
            
    end function sfind
                
end module Finds                

module NZ

character(len=1),dimension(0:9):: digit
character(len=2),dimension(0:100):: zcode
data zcode /' n',' H','He','Li','Be',' B',' C',' N',' O',' F',  &
            'Ne','Na','Mg','Al','Si',' P',' S','Cl','Ar',' K',  &
            'Ca','Sc','Ti',' V','Cr','Mn','Fe','Co','Ni','Cu',  &
            'Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr',' Y',  &
            'Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In',  &
            'Sn','Sb',' I','Te','Xe','Cs','Ba','La','Ce','Pr',  &
            'Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm',  &
            'Yb','Lu','Hf','Ta',' W','Re','Os','Ir','Pt','Au',  &
            'Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac',  &
            'Th','Pa',' U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm'/
data digit /'0','1','2','3','4','5','6','7','8','9'/
            
contains 

    subroutine get_NZ(Nucleus,N,Z)
    implicit none
    integer,intent(inout):: N,Z
    character(len=*),intent(inout):: Nucleus
    character(len=3):: No
    character(len=2):: Chc
    integer:: i,j,k,l
    integer,dimension(0:2):: d
    
    Nucleus=adjustr(Nucleus)
    l=len(Nucleus)
    do i=0,3
        do j=0,9
            if (Nucleus(l-i:l-i)==digit(j)) exit
        end do
        if (j>9.or.j<0) exit
    end do
    if (i<0.or.i>2) stop ' error in Nucleus name digits'
    i=i-1
    No='   '
    No(3-i:3)=Nucleus(l-i:l)
    Chc(1:2)=Nucleus(l-i-2:l-i-1)
    read(No,'(i3)') N
    do k=0,100
        if (Chc==zcode(k)) exit
    end do
    if (k<0.or.k>100) stop ' error in Nucleus name chemical code'
    Z=k
    N=N-Z
    end subroutine get_NZ
    
    subroutine set_NZ(Nucleus,N,Z) 
    implicit none
    integer,intent(in):: N,Z
    character(len=*),intent(inout):: Nucleus
    integer:: l,i,j
    
    l=len(Nucleus)
    Nucleus=repeat(' ',l)
    if (l<5) stop ' character variable Nucleus too short min 5 ch'
    i=l-5
    if (N+Z>99) then
        write(Nucleus(l-2:l),'(i3)') N+Z      
        Nucleus(l-4:l-3)=zcode(Z)(1:2)
    else if (N+Z>9.and.N+Z<100) then
        write(Nucleus(l-1:l),'(i2)') N+Z      
        Nucleus(l-3:l-2)=zcode(Z)(1:2)
    else if (N+Z<10) then
        write(Nucleus(l:l),'(i1)') N+Z      
        Nucleus(l-2:l-1)=zcode(Z)(1:2)
    end if
    end subroutine set_NZ
    
end module NZ            

module InputNuShell

use Shells
use Partition
use NZ

implicit none

character(len=6):: Nucleus,QShell
character(len=40),dimension(3):: Description
character(len=1):: comment,exclam
character(len=7),private:: filename
character(len=11):: spe_file
character(len=max_no_shells),private:: nucleon_line
integer:: no_cI,no_cN,no_cIs
integer:: tParity
integer:: min_cJ2,max_cJ2,del_cJ2,min_cNu,max_cNu
integer:: min_MSon,max_MSon
integer:: min_cT2,max_cT2,del_cT2,nmp
integer:: output_control
type(spartition):: min_cSNo,max_cSNo
type(spartition):: min_MS_No,max_MS_No
integer,dimension(max_no_shells):: typ,maps
real,dimension(ptdim):: spe
real:: cut_energy
character(len=1),dimension(max_no_shells):: shtype
character(len=1):: typex
logical:: shellx=.false.,shellxt=.false.,cut,qspin=.false.
integer:: bigN,bigZ,shnc

data exclam /'!'/
data filename /'NuShell'/

contains

    subroutine input(inExt)
        implicit none
        character(len=*),intent(in):: inExt

        open(unit=10,file=inExt//'.nus',status='old',action='read',err=100)
        read(10,'(a6)') QShell
        if (QShell=='!QSpin') then
            qspin=.true.
        else
            qspin=.false.
            backspace(unit=10)
        end if    
        call read_comments()
        read (10,'(a6)',err=200,end=201) Nucleus
        Nucleus=adjustr(Nucleus)
        call read_comments()
        read (10,'(a40)',err=200,end=201) Description(1)
        read (10,'(a40)',err=200,end=201) Description(2)
        read (10,'(a40)',err=200,end=201) Description(3)
        call read_comments()
        read (10,*,err=200,end=201) no_cI,no_cN
        no_cIs=no_cI
        if  (no_cI < 0) then
            shellx=.true.
            no_cIs=abs(no_cI)
            no_cI=0
        else
            shellx=.false.
        end  if
        if  (no_cN < 0) then
            shellxt=.true.
            no_cN=abs(no_cN)
            no_cI=no_cIs
            if (no_cI==0) typex='i'
            if (no_cI==1) typex='a'
            if (no_cI==2) typex='b'
            if (no_cI>2) typex='u'
            no_cI=no_cN
            no_cN=0
        else
            shellxt=.false.
        end  if
        call read_comments()
        read (10,*,err=200,end=201) tParity
        call read_comments()
        read (10,*,err=200,end=201) no_shells
        if (no_shells > max_no_shells) then
            print *, ' Number of shells exceeds parameter max_no_shells in module Parameters.'
            print *, ' Increase this parameter and recompile.'
            stop
        end if
        call set_ptdim(no_shells)
        
        call read_comments()
        read (10,*,err=200,end=201) nucleon_line
        max_sj2u=0
        typ=0
        shnc=0
        do shn=1,no_shells
            nucleon(shn)=nucleon_line(shn:shn)
            call read_comments()
            if (shellxt) then
            read(10,*,err=200,end=201) cN(shn),sn(shn),sl(shn),sj2(shn),st2(shn),sp(shn),&
                                       typ(shn)
            else
            read(10,*,err=200,end=201) cN(shn),sn(shn),sl(shn),sj2(shn),st2(shn),sp(shn)
            end if                       
            if (typ(shn)==0) then
                shtype(shn)='i'
            else if (typ(shn)==1) then
                shtype(shn)='a'
            else if (typ(shn)==2) then
                shtype(shn)='b'
            else 
                shtype(shn)='u'
            end if    
            if (shellxt.and.typex=='a'.and.shtype(shn)/='a') then
                shnc=shnc+1
                maps(shnc)=shn
            end if    
            if (shellxt.and.typex=='b'.and.shtype(shn)/='b') then
                shnc=shnc+1
                maps(shnc)=shn
            end if    
            if (sj2(shn) > max_sj2u) max_sj2u=sj2(shn)
        end do
        call read_comments()
        del_cJ2=2
        read (10,*,err=200,end=201) min_cJ2,max_cJ2
        if (max_cJ2/2 > max_code) then
            print * ,' File labelling problem, max_cJ2 too large'
            stop
        end if
        if (max_cJ2-min_cJ2+5 > max_no_cJstates*2) then
             print *, ' Increase parameter max_cJstates in module Parameters.'
             stop
        end if
        if (qspin) then
        call read_comments()
        read (10,*,err=200,end=201) min_cNu,max_cNu
        if (max_cNu>no_cI+no_cN) max_cNu=no_cI+no_cN
        if (max_cNu > max_code) then
            print * ,' File labelling problem, max_cNu too large'
            stop
        end if
        if (max_cNu-min_cNu+5 > (max_sj2+1)*2) then
             print *, ' Increase parameter max_sj2 in module Parameters.'
             stop
        end if
        end if
        call read_comments()
        del_cT2=2
        read (10,*,err=200,end=201) min_cT2,max_cT2,nmp
        if (max_cT2/2 > max_code) then
            print * ,' File labelling problem, max_cT2 too large'
            stop
        end if
        if (abs(nmp/2) > max_code/2) then
            print * ,' File labelling problem, nmp too large'
            stop
        end if
        if (max_cT2-min_cT2+5 > max_no_cTstates*2) then
             print *, ' Increase parameter max_cTstates in module Parameters.'
             stop
        end if
        call read_comments()
        read (10,*,err=200,end=201) min_MSon,max_MSon
        call read_comments()
        read (10,*,err=200,end=201) (min_cSNo%shell(shn),shn=1,no_shells)
        if (shellxt) then
            do shn=1,shnc
                min_cSNo%shell(maps(shn))=0
            end do
        end if
        call read_comments
        read (10,*,err=200,end=201) (max_cSNo%shell(shn),shn=1,no_shells)
        if (shellxt) then
            do shn=1,shnc
                max_cSNo%shell(maps(shn))=0
            end do
        end if
        call read_comments

        read (10,*,err=200,end=201) (min_MS_No%shell(shn),shn=1,no_shells)
        call read_comments
        read (10,*,err=200,end=201) (max_MS_No%shell(shn),shn=1,no_shells)
        
        
        do shn=1,no_shells
            if (max_cSNo%shell(shn) > (sj2(shn)+1)*(st2(shn)+1)) max_cSNo%shell(shn)=(sj2(shn)+1)*(st2(shn)+1)
            if (min_cSNo%shell(shn) > (sj2(shn)+1)*(st2(shn)+1)) min_cSNo%shell(shn)=(sj2(shn)+1)*(st2(shn)+1)
        end do
        call read_comments()
        read (10,*,err=200,end=201) output_control
        cut=.false.
        read (10,*,end=80) del_cJ2,del_cT2
        if (del_cJ2<2.or.btest(del_cJ2,0)) del_cJ2=2
        if (del_cT2<2.or.btest(del_cT2,0)) del_cT2=2
        read (10,'(a11)',end=99) spe_file
        read (10,*,end=99) cut_energy
        if (len_trim(spe_file)/=0)  then
            close(unit=10)
            cut=.true.
            spe_file=adjustr(spe_file)
            open(unit=10,file=spe_file)
            read(10,*) spe(1:no_shells)
            close(unit=10)
            return
         end if
! ---------------------------
!    if (shellxt) then
!        call get_NZ(Nucleus,bigN,bigZ)
!        bigN=bigN+(no_cN+nmp)/2
!        bigZ=bigZ+(no_cN-nmp)/2
!        call set_NZ(Nucleus,bigN,bigZ)
!    end if
    
    
! ----------------------------        
        go to 99
 80     del_cJ2=2
        del_cT2=2            
 99     close(unit=10)
        return
100     print *, ' Error opening ',inExt//'.nus'        
        stop
200     print *, ' Error reading ',inExt//'.nus'
        stop        
201     print *, ' EOF reading ',inExt//'.nus'
        stop        
    end subroutine input
   
    subroutine input_sps(inExt)
        implicit none
        character(len=*),intent(in):: inExt

        open(unit=10,file=inExt//'.sps',status='old',action='read',err=100)
        call read_comments()
        read (10,*,err=200,end=201) no_shells
        if (no_shells > max_no_shells) then
            print *, ' Number of shells exceeds parameter max_no_shells in module Parameters.'
            print *, ' Increase this parameter and recompile.'
            stop
        end if
        call set_ptdim(no_shells)
        
        call read_comments()
        read (10,*,err=200,end=201) nucleon_line
        max_sj2u=0
        typ=0
        shnc=0
        do shn=1,no_shells
            nucleon(shn)=nucleon_line(shn:shn)
            call read_comments()
            if (shellxt) then
            read(10,*,err=200,end=201) cN(shn),sn(shn),sl(shn),sj2(shn),st2(shn),sp(shn),&
                                       typ(shn)
            else
            read(10,*,err=200,end=201) cN(shn),sn(shn),sl(shn),sj2(shn),st2(shn),sp(shn)
            end if                       
            if (typ(shn)==0) then
                shtype(shn)='i'
            else if (typ(shn)==1) then
                shtype(shn)='a'
            else if (typ(shn)==2) then
                shtype(shn)='b'
            else 
                shtype(shn)='u'
            end if    
            if (shellxt.and.typex=='a'.and.shtype(shn)/='a') then
                shnc=shnc+1
                maps(shnc)=shn
            end if    
            if (shellxt.and.typex=='b'.and.shtype(shn)/='b') then
                shnc=shnc+1
                maps(shnc)=shn
            end if    
            if (sj2(shn) > max_sj2u) max_sj2u=sj2(shn)
        end do
        call read_comments()
        read (10,*,err=200,end=201) min_MSon,max_MSon
        call read_comments()
        read (10,*,err=200,end=201) (min_cSNo%shell(shn),shn=1,no_shells)
        if (shellxt) then
            do shn=1,shnc
                min_cSNo%shell(maps(shn))=0
            end do
        end if
        call read_comments
        read (10,*,err=200,end=201) (max_cSNo%shell(shn),shn=1,no_shells)
        if (shellxt) then
            do shn=1,shnc
                max_cSNo%shell(maps(shn))=0
            end do
        end if
        call read_comments

        read (10,*,err=200,end=201) (min_MS_No%shell(shn),shn=1,no_shells)
        call read_comments
        read (10,*,err=200,end=201) (max_MS_No%shell(shn),shn=1,no_shells)
        
        
        do shn=1,no_shells
            if (max_cSNo%shell(shn) > (sj2(shn)+1)*(st2(shn)+1)) max_cSNo%shell(shn)=(sj2(shn)+1)*(st2(shn)+1)
            if (min_cSNo%shell(shn) > (sj2(shn)+1)*(st2(shn)+1)) min_cSNo%shell(shn)=(sj2(shn)+1)*(st2(shn)+1)
        end do
        call read_comments()
        read (10,*,err=200,end=201) output_control         
 99     close(unit=10)
        return
100     print *, ' Error opening ',inExt//'.nus'        
        stop
200     print *, ' Error reading ',inExt//'.nus'
        stop        
201     print *, ' EOF reading ',inExt//'.nus'
        stop        
    end subroutine input_sps
   
   
    subroutine read_comments(more)
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
    
    end subroutine read_comments

end module InputNuShell

module Mscheme

use Parameters
use Shells
use Partition
use Machine

implicit none

integer,parameter:: max_sj2p1 = ((max_sj2+1)/2)*2
integer,parameter:: max_sj2p1d2 =max_sj2p1/2
integer,parameter:: m1dim =((max_sj2p1-1)/nbytes+1)*nbytes
integer,parameter:: m1d2dim =((max_sj2p1d2-1)/nbytes+1)*nbytes


integer(kind=2),dimension(0:m1d2dim,0:m1dim):: sm2v
integer(kind=2),dimension(0:m1d2dim,0:m1d2dim):: sm2vp

contains

    subroutine setup_mscheme()
        implicit none
        integer:: i,j,k,m2,st
        logical:: btest
        
        st=0
        if (btest(max_sj2u,0)) st=1
        
        do i=st,max_sj2p1d2
            j=0
            k=0
            do m2=-sj2f(i,st),sj2f(i,st),2
                j=j+1
                sm2v(i,j)=m2
                if (m2 >= 0) then
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
                           
end module Mscheme

module ShellsX

use Shells
use Mscheme
use Mvector
use Partition

implicit none
integer,parameter:: rx=8
integer,private:: ii

integer,dimension(0:nbits,0:nwords):: sj2_lu,sm2_lu,shn_lu,st2_lu,stz2_lu,sp_lu
type(m_vector),dimension(-max_sj2:max_sj2,-1:1):: mvec_jz_tz
real(kind=rx),dimension(0:nbits,0:nwords):: fact1_lu,fact2_lu

contains

    subroutine setup_lu()
        implicit none
        integer:: ib,iw,j
        
        integer,dimension(2):: opa
    
        ib=0
        iw=0
        sj2_lu=0
        sm2_lu=0
        shn_lu=0
        st2_lu=0
        stz2_lu=0
        fact1_lu=real(0,rx)
        fact2_lu=real(0,rx)
        sp_lu=0
        do j=-max_sj2,max_sj2,2
        mvec_jz_tz(j,-1)%fn=int(0,ik)
        mvec_jz_tz(j,1)%fn=int(0,ik)
        end do
        
        do shn=1,no_shells
            iw=shell_word(shn)
            opa(1)=iw
            ib=shell_start(shn)
            do j=1,sj2(shn)+1
                sj2_lu(ib,iw)=sj2(shn)
                st2_lu(ib,iw)=st2(shn)
                shn_lu(ib,iw)=shn
                sp_lu(ib,iw)=sp(shn)
                if (st2(shn) == 1) then
                   
                        sm2_lu(ib,iw)=sm2v(isj2(sj2(shn)),j)
                        stz2_lu(ib,iw)=-st2(shn)
                        opa(2)=ib
                        mvec_jz_tz(sm2_lu(ib,iw),stz2_lu(ib,iw))=   &
                                 opa+mvec_jz_tz(sm2_lu(ib,iw),stz2_lu(ib,iw))
                        fact1_lu(ib,iw)=sqrt(real((sj2_lu(ib,iw)-sm2_lu(ib,iw)),rx)*   &
                            real((sj2_lu(ib,iw)+sm2_lu(ib,iw)+2),rx)/real(4,rx))
                        fact2_lu(ib,iw)=sqrt(real((sj2_lu(ib,iw)+sm2_lu(ib,iw)),rx)*   &
                            real((sj2_lu(ib,iw)-sm2_lu(ib,iw)+2),rx)/real(4,rx))
                        ib=ib+1
                        sj2_lu(ib,iw)=sj2(shn)
                        st2_lu(ib,iw)=st2(shn)
                        stz2_lu(ib,iw)=st2(shn)
                        shn_lu(ib,iw)=shn
                        sp_lu(ib,iw)=sp(shn)
                        sm2_lu(ib,iw)=sm2v(isj2(sj2(shn)),j)
                        opa(2)=ib
                        mvec_jz_tz(sm2_lu(ib,iw),stz2_lu(ib,iw))=   &
                                 opa+mvec_jz_tz(sm2_lu(ib,iw),stz2_lu(ib,iw))
                        fact1_lu(ib,iw)=sqrt(real((sj2_lu(ib,iw)-sm2_lu(ib,iw)),rx)*   &
                            real((sj2_lu(ib,iw)+sm2_lu(ib,iw)+2),rx)/real(4,rx))
                        fact2_lu(ib,iw)=sqrt(real((sj2_lu(ib,iw)+sm2_lu(ib,iw)),rx)*   &
                            real((sj2_lu(ib,iw)-sm2_lu(ib,iw)+2),rx)/real(4,rx))
                        ib=ib+1
                   
                else
                    sm2_lu(ib,iw)=sm2v(isj2(sj2(shn)),j)
                    if (nucleon(shn) == 'n') stz2_lu(ib,iw) = 1
                    if (nucleon(shn) == 'p') stz2_lu(ib,iw) = -1
                    opa(2)=ib
                    mvec_jz_tz(sm2_lu(ib,iw),stz2_lu(ib,iw))=   &
                              opa+mvec_jz_tz(sm2_lu(ib,iw),stz2_lu(ib,iw))
                    fact1_lu(ib,iw)=sqrt(real((sj2_lu(ib,iw)-sm2_lu(ib,iw)),rx)*   &
                            real((sj2_lu(ib,iw)+sm2_lu(ib,iw)+2),rx)/real(4,rx))
                    fact2_lu(ib,iw)=sqrt(real((sj2_lu(ib,iw)+sm2_lu(ib,iw)),rx)*   &
                            real((sj2_lu(ib,iw)-sm2_lu(ib,iw)+2),rx)/real(4,rx))
                    ib=ib+1
                end if
                
            end do
        end do
        
    end subroutine setup_lu
    
    function mpart(mvec)
        implicit none
        type(spartition):: mpart
        type(m_vector),intent(in):: mvec
        type(m_vector):: mvec_copy
        integer:: iw,nb
        
        mpart=0
        mvec_copy=mvec
        do iw=1,nwords
            do
                if (mvec_copy%fn(iw) == 0) exit
                nb=next_bit(mvec_copy%fn(iw))
                mpart%shell(shn_lu(nb,iw))=mpart%shell(shn_lu(nb,iw))+sm2_lu(nb,iw)
            end do
        end do
        
    end function mpart
    
    function tpart(mvec)
        implicit none
        type(spartition):: tpart
        type(m_vector),intent(in):: mvec
        type(m_vector):: mvec_copy
        integer:: iw,nb
        
        tpart=0
        mvec_copy=mvec
        do iw=1,nwords
            do
                if (mvec_copy%fn(iw) == 0) exit
                nb=next_bit(mvec_copy%fn(iw))
                tpart%shell(shn_lu(nb,iw))=tpart%shell(shn_lu(nb,iw))+stz2_lu(nb,iw)
            end do
        end do
        
    end function tpart

    function ptcle(mvec,np)
        implicit none
        integer(kind=2),dimension(2,nwords*nbits):: ptcle
        type(m_vector),intent(in):: mvec
        integer,intent(out):: np
        type(m_vector):: mvec_copy
        integer:: iw,nb
        
        ptcle=0
        np=0
        mvec_copy=mvec
        do iw=1,nwords
            do
                if (mvec_copy%fn(iw) == 0) exit
                nb=next_bit(mvec_copy%fn(iw))
                np=np+1
                ptcle(1,np)=nb
                ptcle(2,np)=iw
            end do
        end do
        
    end function ptcle

end module ShellsX

module ClebschGordan
! Library of angular momentum coupling coefficient routines in fortran 90
! Paul Stevenson, Oxford University/Oak Ridge National Laboratory.
! spaul@mail.phy.ornl.gov
use OperParameters

implicit none 
integer,parameter::rr=rc
logical,private,save:: not_done=.true.
contains

  function cleb(j1,j2,j,m1,m2,m)
    implicit none
    ! calculate a clebsch-gordan coefficient < j1/2 m1/2 j2/2 m2/2 | j/2 m/2 >
    ! arguments are integer and twice the true value. 

    real(kind=rc)    :: cleb,factor,sum
    integer,intent(in):: j1,m1,j2,m2,j,m
    integer::par,z,zmin,zmax

    ! some checks for validity (let's just return zero for bogus arguments)

    if (2*(j1/2)-int(2*(j1/2.0)) /= 2*(abs(m1)/2)-int(2*(abs(m1)/2.0)) .or. &
         2*(j2/2)-int(2*(j2/2.0)) /= 2*(abs(m2)/2)-int(2*(abs(m2)/2.0)) .or. &
         2*(j/2)-int(2*(j/2.0)) /= 2*(abs(m)/2)-int(2*(abs(m)/2.0)) .or. &
         j1<0 .or. j2<0 .or. j<0 .or. abs(m1)>j1 .or. abs(m2)>j2 .or.&
         abs(m)>j .or. j1+j2<j .or. abs(j1-j2)>j .or. m1+m2/=m) then
       cleb= real(0,rc)
    else
    
       factor = real(0,rc)
       factor = binm(j1,(j1+j2-j)/2) / binm((j1+j2+j+2)/2,(j1+j2-j)/2)
       factor = factor * binm(j2,(j1+j2-j)/2) / binm(j1,(j1-m1)/2)
       factor = factor / binm(j2,(j2-m2)/2) / binm(j,(j-m)/2)
       factor = sqrt(factor)
       
       zmin = max(0,j2+(j1-m1)/2-(j1+j2+j)/2,j1+(j2+m2)/2-(j1+j2+j)/2)
       zmax = min((j1+j2-j)/2,(j1-m1)/2,(j2+m2)/2)
       
       sum= real(0,rc)
       do z = zmin,zmax
          par=1
          if(2*(z/2)-int(2*(z/2.0)) /= 0) par=-1
          sum=sum+par*binm((j1+j2-j)/2,z)*binm((j1-j2+j)/2,(j1-m1)/2-z)*&
               binm((-j1+j2+j)/2,(j2+m2)/2-z)
       end do
       
       cleb = factor*sum
    end if
  end function cleb
  
  function threej(j1,j2,j,m1,m2,m)
    implicit none
    integer,intent(in):: j1,m1,j2,m2,j,m
    real(kind=rc):: threej,phase
    
    phase=real(1,rc)
    if (btest((j2+m-j1)/2,0)) phase=-phase
    threej=cleb(j1,j2,j,m1,m2,-m)*phase/sqrt(real(j+1,rc))
    
  end function threej

  function radial(n1,l1,n2,l2,l)
!bab this is the radial me:   <n1,l1| r**l |n2,l2>
!bab in units of b**l
      implicit none
      integer,intent(in)::n1,n2,l1,l2,l
      real(kind=rr)::xn1,xn2,xl1,xl2
      real(kind=rr):: radial,xx,x2,x3,xi,xj,x6p,x6,x7,x8,par,s,t,xmm
      integer:: i,j,im,jm
      real(kind=rr),parameter:: pi = 3.14159265358979
      xn1=n1
      xl1=l1
      xn2=n2
      xl2=l2
      xx = l
      x2 = real(2,rr)*xl1+real(2,rr)*xn1+real(1,rr)
      x3 = real(2,rr)*xl2+real(2,rr)*xn2+real(1,rr)
      im = xn1+real(1,rr)
      jm = xn2 + real(1,rr)
      s = real(0,rr)
      do i = 1,im
        do j = 1,jm
           xi = i-1
           xj = j-1
           if(btest(int(xl1+xl2+xx),0)) x6p = xi+xj+(xl1+xl2+xx+real(1,rr))/real(2,rr)
           x6 = xl1+xl2+real(2,rr)*xi+real(2,rr)*xj+xx+real(1,rr)
           x7 = real(2,rr)*xl1+real(2,rr)*xi+real(1,rr)
           x8 = real(2,rr)*xl2+real(2,rr)*xj+real(1,rr)
           par=real(1,rr)
           if (btest(int(xi+xj),0)) par=-par
           t = par/(factorial(int(xi))*factorial(int(xj)))
           if(.not.btest(int(xl1+xl2+xx),0)) go to 31
30         t = t*factorial(int(x6p))
           xmm = (x6+real(1,rr))/real(2,rr)
           go to 32
31         t = t*dfactorial(int(x6))
32         t = t/(dfactorial(int(x7))*dfactorial(int(x8))*factorial(int(xn1-xi))*factorial(int(xn2-xj)))
           if(.not.btest(int(xl1+xl2+xx),0)) go to 13
40         t = t*(real(2,rr)**xmm)/sqrt(pi)
13         s = s + t
        end do
      end do
      s = s*sqrt(dfactorial(int(x2))*dfactorial(int(x3)))
      s = s*sqrt(factorial(int(xn1))*factorial(int(xn2)))
      s = s/sqrt(real(2,rr)**(xn1+xn2+xx))
      radial = s
      return
  end function radial
  
  function ZN(cN,m1,m2)
     implicit none
     real(kind=rc):: ZN
     integer,intent(in):: cN
     real(kind=rc),intent(in)::m1,m2
     
           ZN=(m1/(m1-m2))**(real(cN)/2.0)
           
  end function ZN       
  
    
  function oradial(n,l,np,lp,p)
    implicit none
    integer, intent(in):: n,l,np,lp,p
    real(kind=rc):: oradial
    integer:: sm,ss,s
    
    if (p==0) then
        if (n==np .and. l==lp)  then
            oradial=real(1,rc)
        else
            oradial=real(0,rc)
        end if
        return
    end if

    if (btest(l+lp+p,0)) then
        oradial=real(0,rc)
        return
    end if

    sm=min(n,np)
    ss=max(np-(l-lp+p)/2,n-(lp-l+p)/2,0)
    oradial=real(0,rc)
    do s=ss,sm
        oradial=oradial+dfactorial(l+lp+p+2*s+1)/ &
               ((reaL(2,rc)**s)*factorial(s)*factorial(n-s)*factorial(np-s)* &
                factorial(s+(l-lp+p)/2-np)*factorial(s+(lp-l+p)/2-n))
    end do
    if (oradial==real(0,rc)) then
        return
    end if
    oradial=oradial*factorial((l-lp+p)/2)*factorial((lp-l+p)/2)*real(-1,rc)**(n-np)
    oradial=oradial*sqrt(factorial(n)*factorial(np)*(real(2,rc)**(n+np-p))/ &
           (dfactorial(2*n+2*l+1)*dfactorial(2*np+2*lp+1)))
  
  end function oradial
  
  function racah(j1,j2,j3,j4,J,K)
    implicit none
    integer, intent(in):: j1,j2,j3,j4,J,K
    real(kind=rc):: racah
    real(kind=rc):: phase
    
    phase=real(1,rc)
    if (btest((j1+j2+j3+j4)/2,0)) phase=-phase
    racah=sixj(j1,j2,J,j4,j3,K)*phase
    
  end function racah
    
  
  function sixj(a,b,c,d,e,f)
    implicit none
    integer, intent(in) :: a,b,c,d,e,f
    real(kind=rc) :: sixj
    integer :: nlo, nhi, n
    real(kind=rc) :: outfactors, sum, sumterm
    ! calculates a Wigner 6-j symbol. Argument a-f are integer and are
    ! twice the true value of the 6-j's arguments, in the form
    ! { a b c }
    ! { d e f }
    ! Calculated using binmial coefficients to allow for (reasonably) high
    ! arguments.
    ! First check for consistency of arguments:
    sixj= real(0,rc)
    if(mod(a+b,2)/=mod(c,2)) return
    if(mod(c+d,2)/=mod(e,2)) return
    if(mod(a+e,2)/=mod(f,2)) return
    if(mod(b+d,2)/=mod(f,2)) return
    if(abs(a-b)>c .or. a+b<c) return
    if(abs(c-d)>e .or. c+d<e) return
    if(abs(a-e)>f .or. a+e<f) return
    if(abs(b-d)>f .or. b+d<f) return

    outfactors = angdelta(a,e,f)/angdelta(a,b,c)
    outfactors = outfactors * angdelta(b,d,f)*angdelta(c,d,e)

    nlo = max( (a+b+c)/2, (c+d+e)/2, (b+d+f)/2, (a+e+f)/2 )
    nhi = min( (a+b+d+e)/2, (b+c+e+f)/2, (a+c+d+f)/2)

    sum= real(0,rc)
    do n=nlo,nhi
       sumterm = real((-1)**n,rc)
       sumterm = sumterm * binm(n+1,n-(a+b+c)/2)
       sumterm = sumterm * binm((a+b-c)/2,n-(c+d+e)/2)
       sumterm = sumterm * binm((a-b+c)/2,n-(b+d+f)/2)
       sumterm = sumterm * binm((b-a+c)/2,n-(a+e+f)/2)
       sum=sum+sumterm
    end do

    sixj = sum * outfactors

  end function sixj

  function angdelta(a,b,c)
    implicit none
    integer :: a,b,c
    real(kind=rc)    :: angdelta, scr1
    ! calculate the function delta as defined in varshalovich et al. for
    ! use in 6-j symbol:
    scr1= factorial((a+b-c)/2)
    scr1=scr1/factorial((a+b+c)/2+1)
    scr1=scr1*factorial((a-b+c)/2)
    scr1=scr1*factorial((-a+b+c)/2)
    angdelta=sqrt(scr1)
  end function angdelta

  function ninej(a,b,c,d,e,f,g,h,i)
    implicit none
    integer :: a,b,c,d,e,f,g,h,i
    real(kind=rc)    :: ninej, sum
    integer :: xlo, xhi
    integer :: x
    ! calculate a 9-j symbol. The arguments are given as integers twice the
    ! value of the true arguments in the form
    ! { a b c }
    ! { d e f }
    ! { g h i }
    ninej= real(0,rc)
    ! first check for bogus arguments (and return zero if so)
    if(abs(a-b)>c .or. a+b<c) return
    if(abs(d-e)>f .or. d+e<f) return
    if(abs(g-h)>i .or. g+h<i) return
    if(abs(a-d)>g .or. a+d<g) return
    if(abs(b-e)>h .or. b+e<h) return
    if(abs(c-f)>i .or. c+f<i) return
    
    xlo = max(abs(b-f),abs(a-i),abs(h-d))
    xhi = min(b+f,a+i,h+d)
    
    sum= real(0,rc)
    do x=xlo,xhi,2
       sum=sum+(-1)**x*(x+1)*sixj(a,b,c,f,i,x)*sixj(d,e,f,b,x,h)*&
            sixj(g,h,i,x,a,d)
    end do
    ninej=sum

  end function ninej

  recursive function factorial(n) result(res)
    implicit none
    integer :: n
    real(kind=rr) :: res

    if (n==0 .or. n==1) then
       res= real(1,rc)
    else
       res=real(n,rc)*factorial(n-1)
    end if
  end function factorial

  recursive function dfactorial(n) result(res)
    implicit none
    integer :: n
    real(kind=rr) :: res

    if (n==0 .or. n==1 .or. n < 0) then
       res= real(1,rc)
    else
       res=real(n,rc)*dfactorial(n-2)
    end if
  end function dfactorial
  
  recursive function Ld(n,L,M) result(res)
    implicit none
    real(kind=rc) :: res
    integer :: L,M,n
    
    if (n==0) then
        res = real(1,rc)
        return
    end if
    if (n<0) then
        res = real(0,rc)
        return
    end if
    res  = sqrt(real((L+M+n)*(L-M-n+1),rc)/2.0d0)*Ld(n-1,L,M)
    
  end function Ld      
   
    
  
  recursive function binm(n,r) result(res)
    implicit none
    integer :: n,r
    real(kind=rc) :: res

    if(n==r .or. r==0) then
       res =  real(1,rc)
    else if (r==1) then
       res = real(n,rc)
    else
       res = real(n,rc)/real(n-r,rc)*binm(n-1,r)
    end if
  end function binm
end module ClebschGordan

module PartOrder

! Note PartOrder and OpOrder share node array
! Do not use both in same program

use OrderParameters
use Partition

implicit none

type part_tree
    type(spartition):: part
    real(kind=ro):: V
    integer:: node
    integer:: greater
    integer:: lesser
end type part_tree
    
type(part_tree),dimension(:),allocatable:: ipart_order
integer,dimension(:),allocatable:: node
integer:: notree,nonode
integer,private,save:: next_ipart=1


contains

    function read_next_ipart()
        implicit none
        integer:: read_next_ipart
        
        read_next_ipart=next_ipart
        
    end function read_next_ipart
        
    subroutine reset_ipart_order()
        implicit none
        
        next_ipart=1
        ipart_order(1)%part=0
        ipart_order(1)%V=real(0,ro)
        ipart_order(1)%node=1
        ipart_order(1)%greater=0
        ipart_order(1)%lesser=0
        
    end subroutine reset_ipart_order
    
    subroutine add_part(ipart,V,add,new)
        implicit none
        type(spartition),intent(in):: ipart
        real(kind=ro):: V
        logical,intent(in),optional::add
        logical,intent(out),optional::new
        integer:: i
        
        if (present(new)) new=.true.
        
        !store ipart in a binary tree
        
        if (next_ipart == 1) then
            ipart_order(1)%part=ipart
            ipart_order(1)%V=V
            next_ipart=next_ipart+1
            if (next_ipart > max_part_tree) call ipart_message()
            return
        end if
        
        i=1
        
10      if (ipart == ipart_order(i)%part)  then
            if (present(add).and. add) then
            ipart_order(i)%V=ipart_order(i)%V+V
            if (present(new)) new=.false.
            else if (present(add).and. .not.add) then
            ipart_order(i)%V=V
            if (present(new)) new=.false.
            else if (.not.present(add)) then
            stop  ' Identical matrix elements in me file'
            end if
            return
        end if
        if (ipart > ipart_order(i)%part) then
            if (ipart_order(i)%greater == 0) then
                    ipart_order(i)%greater=next_ipart
                    ipart_order(next_ipart)%node=i
                    ipart_order(next_ipart)%part=ipart
                    ipart_order(next_ipart)%V=V
                    ipart_order(next_ipart)%greater=0
                    ipart_order(next_ipart)%lesser=0
                    next_ipart=next_ipart+1
                    if (next_ipart > max_part_tree) call ipart_message()
                    return
            else
                i=ipart_order(i)%greater
                go to 10
            end if
        end if
        if (ipart < ipart_order(i)%part) then
            if (ipart_order(i)%lesser == 0) then
                    ipart_order(i)%lesser=next_ipart
                    ipart_order(next_ipart)%node=i
                    ipart_order(next_ipart)%part=ipart
                    ipart_order(next_ipart)%V=V
                    ipart_order(next_ipart)%greater=0
                    ipart_order(next_ipart)%lesser=0
                    next_ipart=next_ipart+1
                    if (next_ipart > max_part_tree) call ipart_message()
                    return
            else
                i=ipart_order(i)%lesser
                go to 10
            end if
         end if
         
         ! should never reach here
         
         print *, ' Error 2 in module Order'                      !ERROR2
         
    end subroutine add_part
    
    
    function find_part_V(ipart,found)
        implicit none
        type(spartition),intent(in):: ipart
        logical,intent(out)::found
        real(kind=ro):: find_part_V
        integer:: i
        
        !find ipart in a binary tree
        found=.false.
        find_part_V=real(0,ro)
        if (next_ipart == 1) then
            return
        end if
        
        i=1
        
10      if (ipart == ipart_order(i)%part)  then
            found=.true.
            find_part_V=ipart_order(i)%V
            return
        end if
        if (ipart > ipart_order(i)%part) then
            if (ipart_order(i)%greater == 0) then
                    return
            else
                i=ipart_order(i)%greater
                go to 10
            end if
        end if
        if (ipart < ipart_order(i)%part) then
            if (ipart_order(i)%lesser == 0) then
                    return
            else
                i=ipart_order(i)%lesser
                go to 10
            end if
         end if
         
         ! should never reach here
         
         print *, ' Error 4 in module Order'                      !ERROR4
         
    end function find_part_V
    
                        
    function get_next_ipart(done)
        implicit none
        integer:: get_next_ipart
        logical,intent(inout):: done
        integer,save:: i,gl,n,nd
        
        !read out svec's in increasing order of sort variables
        
        if (done) then
            i=1
            n=1
            nd=1
            node=0
            gl=-1
            done=.false.
        end if
        
        if (gl == -1) go to 20
        if (gl == +1) go to 30
        
20      if (ipart_order(i)%lesser /=0 ) then 
            i=ipart_order(i)%lesser                  !go down tree
            go to 20
        else 
            get_next_ipart=i                             !bottom of tree, return value
            n=n+1
            if (n == next_ipart) done=.true.
            gl=+1                                   !ready to go back up
            return
        end if
        
30      if (ipart_order(i)%greater /=0 ) then
            i=ipart_order(i)%greater                 !found upward branch
            nd=nd+1
            if (nd > max_node) then
                print *, ' Increase parameter max_node in module OrderParameters or'
                print *, ' SET NODE=I where I >>',nd
                stop
            end if
            node(nd)=ipart_order(i)%node             !push current node on stack
            if (gl == -1) then
                get_next_ipart=ipart_order(i)%node        !last movement was down, return value
                return
            else
                gl=-1                               !get set to go down branch
                go to 20
            end if
        else
40          i=ipart_order(i)%node
            if (i == node(nd)) then                 !if node equals to of sstack
                nd=nd-1                             !pop stack
                if (nd < 1) then
                    done=.true.                     !reached root node, done
                    get_next_ipart=0
                    return
                end if
                go to 40                            !go back up
            else 
                go to 50
            end if
50          n=n+1                                   !return value at top of branch
            get_next_ipart=i
            if (n == next_ipart) done=.true.
            gl=+1                                   !prepare to go up again or pop
            return
        end if
        
        done=.true.                                 !should never reach here NO ERROR TRAP!!!  
              
       print *, ' Error 4 in module Order'                      !ERROR4
       stop
            
    end function get_next_ipart
                        
    
    subroutine ipart_message()
        implicit none
        
        print *, ' If running NuOper then :'
        print * ,' Increase parameter max_part_tree in module OrDerParameters or'
        print * ,' SET TREE=I where I >>',notree+1
        print *, ' If running NuMe then :'
        print * ,' Increase parameter max_me in module OrDerParameters or'
        print * ,' SET NOME=I where I >>',notree+1
        stop
        
    end subroutine ipart_message        
                
end module PartOrder

module OpOrder

! Note PartOrder and OpOrder share node array
! Do not use both in same program

use OrderParameters
use OpPartition
use PartOrder

implicit none

type opart_tree
    type(opartition):: part
    integer:: node
    integer:: greater
    integer:: lesser
end type opart_tree
    
type(opart_tree),dimension(:),allocatable:: opart_order
integer,private,save:: next_opart=1
integer:: nopart

contains

    function read_next_opart()
        implicit none
        integer:: read_next_opart
        
        read_next_opart=next_opart
        
    end function read_next_opart
        
    subroutine reset_opart_order()
        implicit none
        
        next_opart=1
        opart_order(1)%part=0
        opart_order(1)%node=1
        opart_order(1)%greater=0
        opart_order(1)%lesser=0
        
    end subroutine reset_opart_order
    
    subroutine add_op(opart)
        implicit none
        type(opartition),intent(in):: opart
        integer:: i
        
        !store opart in a binary tree
        
        if (next_opart == 1) then
            opart_order(1)%part=opart
            next_opart=next_opart+1
            if (next_opart > nopart) call opart_message()
            return
        end if
        
        i=1
        
10      if (opart == opart_order(i)%part)  then
            return
        end if
        if (opart > opart_order(i)%part) then
            if (opart_order(i)%greater == 0) then
                    opart_order(i)%greater=next_opart
                    opart_order(next_opart)%node=i
                    opart_order(next_opart)%part=opart
                    opart_order(next_opart)%greater=0
                    opart_order(next_opart)%lesser=0
                    next_opart=next_opart+1
                    if (next_opart > nopart) call opart_message()
                    return
            else
                i=opart_order(i)%greater
                go to 10
            end if
        end if
        if (opart < opart_order(i)%part) then
            if (opart_order(i)%lesser == 0) then
                    opart_order(i)%lesser=next_opart
                    opart_order(next_opart)%node=i
                    opart_order(next_opart)%part=opart
                    opart_order(next_opart)%greater=0
                    opart_order(next_opart)%lesser=0
                    next_opart=next_opart+1
                    if (next_opart > nopart) call opart_message()
                    return
            else
                i=opart_order(i)%lesser
                go to 10
            end if
         end if
         
         ! should never reach here
         
         print *, ' Error 2 in module Order'                      !ERROR2
         
    end subroutine add_op
    
    
    function find_op(opart,found)
        implicit none
        type(opartition),intent(in):: opart
        logical,intent(out)::found
        integer:: find_op
        integer:: i
        
        !find opart in a binary tree
        found=.false.
        find_op=0
        if (next_opart == 1) then
            return
        end if
        
        i=1
        
10      if (opart == opart_order(i)%part)  then
            found=.true.
            find_op=i
            return
        end if
        if (opart > opart_order(i)%part) then
            if (opart_order(i)%greater == 0) then
                    return
            else
                i=opart_order(i)%greater
                go to 10
            end if
        end if
        if (opart < opart_order(i)%part) then
            if (opart_order(i)%lesser == 0) then
                    return
            else
                i=opart_order(i)%lesser
                go to 10
            end if
         end if
         
         ! should never reach here
         
         print *, ' Error 4 in module Order'                      !ERROR4
         
    end function find_op
    
                        
    function get_next_opart(done)
        implicit none
        integer:: get_next_opart
        logical,intent(inout):: done
        integer,save:: i,gl,n,nd
        integer,dimension(max_node):: node
        
        !read out svec's in increasing order of sort variables
        
        if (done) then
            i=1
            n=1
            nd=1
            node=0
            gl=-1
            done=.false.
        end if
        if (next_opart==1) then
            done=.true.
            get_next_opart=0
            return
        end if
        if (gl == -1) go to 20
        if (gl == +1) go to 30
        
20      if (opart_order(i)%lesser /=0 ) then 
            i=opart_order(i)%lesser                  !go down tree
            go to 20
        else 
            get_next_opart=i                             !bottom of tree, return value
            n=n+1
            if (n == next_opart) done=.true.
            gl=+1                                   !ready to go back up
            return
        end if
        
30      if (opart_order(i)%greater /=0 ) then
            i=opart_order(i)%greater                 !found upward branch
            nd=nd+1
            if (nd > nonode) then
                print *, ' Increase parameter nonode'
                stop
            end if
            node(nd)=opart_order(i)%node             !push current node on stack
            if (gl == -1) then
                get_next_opart=opart_order(i)%node        !last movement was down, return value
                return
            else
                gl=-1                               !get set to go down branch
                go to 20
            end if
        else
40          i=opart_order(i)%node
            if (i == node(nd)) then                 !if node equals to of sstack
                nd=nd-1                             !pop stack
                if (nd < 1) then
                    done=.true.                     !reached root node, done
                    get_next_opart=0
                    return
                end if
                go to 40                            !go back up
            else 
                go to 50
            end if
50          n=n+1                                   !return value at top of branch
            get_next_opart=i
            if (n == next_opart) done=.true.
            gl=+1                                   !prepare to go up again or pop
            return
        end if
        
        done=.true.                                 !should never reach here NO ERROR TRAP!!!  
              
       print *, ' Error 4 in module Order'                      !ERROR4
       stop
            
    end function get_next_opart
                        
    
    subroutine opart_message()
        implicit none
        
        print * ,' Increase parameter nopart'
        stop
        
    end subroutine opart_message        
                
end module OpOrder

module PartMvec

use Mvector
use Partition
use Shells

implicit none

interface assignment (=)
    module procedure part_aequals_mvec
end interface


contains

	subroutine part_aequals_mvec(part,mvec)
	    implicit none
	    type(spartition),intent(out):: part
	    type(m_vector),intent(in):: mvec
	    integer:: i
	    integer(kind=ik),dimension(max_no_shells):: mfn
	    
	    
	    do i=1,no_shells
	        mfn(i)=iand(mvec%fn(shell_word(i)),shell_mask(i))
	        part%shell(i)=countb(mfn(i))
	    end do
	    
	end subroutine part_aequals_mvec

end module PartMvec
