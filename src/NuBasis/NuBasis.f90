! NuShell, NuShell@MSU, SunShell main source for NuBasis/sunbasis
! All modules in this file are specific to the (generic) pBasis
! program. This file contains the .main function. 

module Combine

contains

subroutine comb ( m, n, l, iarray )
!
!*******************************************************************************
!
!! COMB returns the L-th combination of N things out of M.
!
!
!  Discussion:
!
!    The combinations are ordered lexically.
!
!    Lexical order can be illustrated for the general case of N and M as
!    follows:
!
!    1:       1,     2,     3,     ..., N-2, N-1, N
!    2:       1,     2,     3,     ..., N-2, N-1, N+1
!    3:       1,     2,     3,     ..., N-2, N-1, N+2
!    ...
!    M-N+1:   1,     2,     3,     ..., N-2, N-1, M
!    M-N+2:   1,     2,     3,     ..., N-2, N,   N+1
!    M-N+3:   1,     2,     3,     ..., N-2, N,   N+2
!    ...
!    LAST-2:  M-N,   M-N+1, M-N+3, ..., M-2, M-1, M
!    LAST-1:  M-N,   M-N+2, M-N+3, ..., M-2, M-1, M
!    LAST:    M-N+1, M-N+2, M-N+3, ..., M-2, M-1, M
!
!    There are a total of M!/(N!*(M-N)!) combinations of M
!    things taken N at a time.
!
!  Reference:
!
!    B P Buckles and M Lybanon,
!    Algorithm 515,
!    Generation of a Vector from the Lexicographical Index,
!    ACM Transactions on Mathematical Software,
!    Volume 3, Number 2, pages 180-182, June 1977.
!
!  Modified:
!
!    15 April 1999
!
!  Parameters:
!
!    Input, integer M, the size of the set.
!
!    Input, integer N, the number of things in the combination.
!    N must be greater than 0, and no greater than M.
!
!    Input, integer L, the lexicographical index of combination
!    sought.  L must be at least 1, and no greater than M!/(N!*(M-N)!).
!
!    Output, integer IARRAY(N), array containing the combination set.
!
  implicit none
!
  integer n
!
  integer i
  integer iarray(n)
  integer j
  integer k
  integer l
  integer m
!
!  Initialize lower bound index at zero.
!
  k = 0
!
!  Loop to select elements in ascending order.
!
  do i = 1, n - 1
!
!  Set lower bound element number for next element value.
!
    iarray(i) = 0

    if ( i /= 1 ) then
      iarray(i) = iarray(i-1)
    end if
!
!  Check each element value.
!
    do

      iarray(i) = iarray(i) + 1
      call combin2 ( m-iarray(i), n-i, j )
      k = k + j

      if ( k >= l ) then
        exit
      end if

    end do

    k = k - j

  end do

  iarray(n) = iarray(n-1) + l - k

  return
end subroutine comb

subroutine comb_next ( done, n, k, iarray )
!
!*******************************************************************************
!
!! COMB_NEXT computes combinations of K things out of N.
!
!
!  Discussion:
!
!    The combinations are computed one at a time, in lexicographical order.
!
!  Reference:
!
!    Charles Mifsud,
!    Combination in Lexicographic Order,
!    ACM algorithm 154,
!    March 1963.
!
!  Modified:
!
!    15 April 1999
!
!  Parameters:
!
!    Input/output, logical DONE, indicator.
!    On input, if this is the first call, the user should set
!    DONE to FALSE.  On each subsequent call, the value of
!    DONE should be left at its output value from the previous
!    call.
!
!    On output, if DONE is TRUE, then a new combination was
!    computed and returned.  If DONE is FALSE, then the list
!    of combinations was finished on the previous call.
!
!    Input, integer N, the total number of things.
!
!    Input, integer K, the number of things in each combination.
!
!    Output, integer IARRAY(K), contains the list of elements in
!    the current combination.
!
  implicit none
!
  integer k
!
  logical done
  integer i
  integer iarray(k)
  integer j
  integer n
!
  if ( done ) then

    call ivec_identity ( k, iarray )

    if ( k > 1 ) then
      done = .false.
    else
      done = .true.
    end if

  else

    if ( iarray(k) < n ) then
      iarray(k) = iarray(k) + 1
      return
    end if

    do i = k, 2, -1

      if ( iarray(i-1) < n-k+i-1 ) then

        iarray(i-1) = iarray(i-1) + 1

        do j = i, k
          iarray(j) = iarray(i-1) + j - ( i-1 )
        end do

        return

      end if

    end do

    done = .true.

  end if

  return
end subroutine comb_next

subroutine combin2 ( n, k, icnk )
!
!*******************************************************************************
!
!! COMBIN2 computes the binomial coefficient C(N,K).
!
!
!  Discussion:
!
!    The value is calculated in such a way as to avoid overflow and
!    roundoff.  The calculation is done in integer arithmetic.
!
!  Formula:
!
!    ICNK = C(N,K) = N! / ( K! * (N-K)! )
!
!  Reference:
!
!    M L Wolfson and H V Wright,
!    Combinatorial of M Things Taken N at a Time,
!    ACM algorithm 160,
!    Communications of the ACM,
!    April, 1963.
!
!  Modified:
!
!    17 January 1999
!
!  Parameters:
!
!    Input, integer N, K, are the values of N and K.
!
!    Output, integer ICNK, the number of combinations of N
!    things taken K at a time.
!
  implicit none
!
  integer i
  integer icnk
  integer k
  integer mn
  integer mx
  integer n
!
  mn = min ( k, n-k )

  if ( mn < 0 ) then

    icnk = 0

  else if ( mn == 0 ) then

    icnk = 1

  else

    mx = max ( k, n-k )
    icnk = mx + 1

    do i = 2, mn
      icnk = ( icnk * ( mx + i ) ) / i
    end do

  end if

  return
end subroutine combin2

subroutine comp_next ( n, k, iarray, more )
!
!*******************************************************************************
!
!! COMP_NEXT computes the compositions of the integer N into K parts.
!
!
!  Discussion:
!
!    The routine computes one composition on each call until there are no more.
!    For instance, one composition of 6 into 3 parts is
!    3+2+1, another would be 6+0+0.
!
!  Example:
!
!    The 28 compositions of 6 into three parts are:
!
!      6 0 0,  5 1 0,  5 0 1,  4 2 0,  4 1 1,  4 0 2,
!      3 3 0,  3 2 1,  3 1 2,  3 0 3,  2 4 0,  2 3 1,
!      2 2 2,  2 1 3,  2 0 4,  1 5 0,  1 4 1,  1 3 2,
!      1 2 3,  1 1 4,  1 0 5,  0 6 0,  0 5 1,  0 4 2,
!      0 3 3,  0 2 4,  0 1 5,  0 0 6.
!
!  Reference:
!
!    A Nijenhuis and H Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Modified:
!
!    15 April 1999
!
!  Parameters:
!
!    Input, integer N, the integer whose compositions are desired.
!
!    Input, integer K, the number of parts in the composition.
!
!    Output, integer IARRAY(K).  IARRAY(I) is the I-th part of
!    the composition.
!
!    Input/output, logical MORE.
!    Set MORE = .FALSE. on first call.  It will be reset to .TRUE. on return
!    with a new composition.  Each new call returns another composition until
!    MORE is set to .FALSE. when the last composition has been computed
!    and returned.
!
  implicit none
!
  integer k
!
  integer iarray(k)
  integer, save :: ih = 0
  integer, save :: it = 0
  logical more
  integer n
!
  if ( .not. more ) then

    it = n
    ih = 0
    iarray(1) = n
    iarray(2:k) = 0

  else

    if ( it > 1 ) then
      ih = 0
    end if

    ih = ih + 1
    it = iarray(ih)
    iarray(ih) = 0
    iarray(1) = it - 1
    iarray(ih+1) = iarray(ih+1) + 1

  end if

  more = iarray(k) /= n

  return
end subroutine comp_next

subroutine ivec_identity ( n, a )
!
!*******************************************************************************
!
!! IVEC_IDENTITY sets an integer vector to the identity vector A(I)=I.
!
!
!  Modified:
!
!    09 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of elements of A.
!
!    Output, integer A(N), the array to be initialized.
!
  implicit none
!
  integer n
!
  integer a(n)
  integer i
!
  do i = 1, n
    a(i) = i
  end do

  return
end subroutine ivec_identity

subroutine pos_next ( done, n, k, a)

    implicit none

    logical done

    integer n , k

    integer a(k)

    integer, save :: p

    if (k /=1 ) return

    if (done) p = 0

    done = .false.

    p = p + 1
    if ( p > n ) done = .true.

    a(1) = p

    return
    
end subroutine pos_next

subroutine perm_next ( n, p, more, even )
!
!*******************************************************************************
!
!! PERM_NEXT computes all of the permutations of N objects, one at a time.
!
!
!  Discussion:
!
!    The routine is initialized by calling with MORE = TRUE, in which case
!    it returns the identity permutation.
!
!    If the routine is called with MORE = FALSE, then the successor of the
!    input permutation is computed.
!
!  Reference:
!
!    A Nijenhuis and H Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Modified:
!
!    30 March 2001
!
!  Parameters:
!
!    Input, integer N, the number of objects being permuted.
!
!    Input/output, integer P(N), the permutation, in standard index form.
!
!    On the first call, the input value is unimportant.
!    On subsequent calls, the input value should be the same
!    as the output value from the previous call.  In other words, the
!    user should just leave P alone.
!
!    On output, contains the "next" permutation.
!
!    Input/output, logical MORE.
!
!    Set MORE = .FALSE. before the first call.
!
!    MORE will be reset to .TRUE. and a permutation will be returned.
!
!    Each new call produces a new permutation until
!    MORE is returned .FALSE.
!
!    Input/output, logical EVEN.
!
!    The input value of EVEN should simply be its output value from the
!    previous call; (the input value on the first call doesn't matter.)
!
!    On output, EVEN is .TRUE. if the output permutation is even, that is,
!    involves an even number of transpositions.
!
  implicit none
!
  integer n
!
  logical even
  integer i
  integer i1
  integer ia
  integer id
  integer is
  integer j
  integer l
  integer m
  logical more
  integer p(*)
!
  if ( .not. more ) then

    call ivec_identity ( n, p )
    more = .true.
    even = .true.

    if ( n == 1 ) then
      more = .false.
      return
    end if

    if ( p(n) /= 1 .or. p(1) /= 2 + mod ( n, 2 ) ) then
      return
    end if

    do i = 1, n-3
      if ( p(i+1) /= p(i)+1 ) then
        return
      end if
    end do

    more = .false.

  else

    if ( n == 1 ) then
      p(1) = 0
      more = .false.
      return
    end if

    if ( even ) then

      ia = p(1)
      p(1) = p(2)
      p(2) = ia
      even = .false.

      if ( p(n) /= 1 .or. p(1) /= 2 + mod ( n, 2 ) ) then
        return
      end if

      do i = 1, n-3
        if ( p(i+1) /= p(i)+1 ) then
          return
        end if
      end do

      more = .false.
      return

    else

      more = .false.

      is = 0

      do i1 = 2, n

        ia = p(i1)
        i = i1 - 1
        id = 0

        do j = 1, i
          if ( p(j) > ia ) then
            id = id + 1
          end if
        end do

        is = id + is
        if ( id /= i * mod ( is, 2 ) ) then
          more = .true.
          exit
        end if

      end do

      if ( .not. more ) then
        p(1) = 0
        return
      end if

    end if

    m = mod ( is+1, 2 ) * ( n + 1 )

    do j = 1, i

      if ( sign ( 1, p(j)-ia ) /= sign ( 1, p(j)-m ) ) then
        m = p(j)
        l = j
      end if

    end do

    p(l) = ia
    p(i1) = m
    even = .true.

  end if

  return
end subroutine perm_next

end module combine

module Counter

use Parameters
use Partition
use Shells

implicit none

integer(kind=2),dimension(:,:,:),allocatable:: combin
integer(kind=2),dimension(:,:),allocatable:: combM2
integer(kind=2),dimension(:,:),allocatable:: combT2

type(spartition),private,save:: max_counter,min_counter

contains

    function set_counter(counter,part,cM,cT,current_counter)
        implicit none
        type(spartition),intent(in):: counter,part
        type(spartition),intent(out):: current_counter
        type(spartition):: set_counter
        integer,intent(inout):: cM,cT
        cM=0
        cT=0
        min_counter=1
        max_counter=counter
        current_counter=min_counter
        do  shn=1,no_shells
            if (part%shell(shn)/=0) then
                cM=cM+combM2(current_counter%shell(shn),shn)
                cT=cT+combT2(current_counter%shell(shn),shn)
            end if
        end do
        set_counter=current_counter
        
    end function set_counter
    
    function reset_counter(part,cM,cT,current_counter)
        implicit none
        type(spartition),intent(in):: part
        type(spartition),intent(out):: current_counter
        type(spartition):: reset_counter
        integer,intent(inout):: cM,cT
        cM=0
        cT=0
        current_counter=min_counter
        do  shn=1,no_shells
            if (part%shell(shn)/=0) then
                cM=cM+combM2(current_counter%shell(shn),shn)
                cT=cT+combT2(current_counter%shell(shn),shn)
            end if
        end do
        reset_counter=current_counter
        
    end function reset_counter
    
    pure subroutine next(next_counter,part,cM,cT,current_counter)
        implicit none
        type(spartition),intent(in):: part
        type(spartition),intent(inout):: current_counter
        type(spartition),intent(out):: next_counter
        integer,intent(inout):: cM,cT
        integer:: digit
                
        do digit=1,no_shells
            if (part%shell(digit)/=0) then
                cM=cM-combM2(current_counter%shell(digit),digit)
                cT=cT-combT2(current_counter%shell(digit),digit)
            end  if
            current_counter%shell(digit)=current_counter%shell(digit)+1
            if (current_counter%shell(digit) <= max_counter%shell(digit)) then
                next_counter=current_counter
                if (part%shell(digit)/=0) then
                    cM=cM+combM2(current_counter%shell(digit),digit)
                    cT=cT+combT2(current_counter%shell(digit),digit)
                end if
                return
            else
                if (digit == no_shells) then
                    current_counter=0
                    next_counter=current_counter
                    cM=0
                    cT=0
                    return
                end if
                    
                current_counter%shell(digit)=min_counter%shell(digit)
                if (part%shell(digit)/=0) then
                    cM=cM+combM2(current_counter%shell(digit),digit)
                    cT=cT+combT2(current_counter%shell(digit),digit)
                end if
            end if
        end do
        
        current_counter=-1
        next_counter=-1
   
   end subroutine next
   
end module Counter

module SetupPartitions

use Combine
use InputNuShell
use Mscheme
use Finds

implicit none

interface assignment (=)
    module procedure array_aequals_chara
end interface


type(spartition),dimension(:),allocatable:: spart
type(spartition),private:: part
integer:: no_spart,nopart
integer,dimension(ptdim):: sttz2
real:: tot_energy

contains

    subroutine setup_isospin()
        implicit none
        
        sttz2=0
        sttz2=nucleon
        
    end subroutine setup_isospin      

    subroutine array_aequals_chara(array,chara)
        implicit none
        integer,dimension(ptdim),intent(out):: array
        character(len=1),dimension(ptdim),intent(in):: chara
        integer:: i
        
        array=0
        do i=1,max_no_shells
            if (chara(i) == 'p') array(i)=-1
            if (chara(i) == 'n') array(i)=+1
        end do
        
    end subroutine array_aequals_chara

    subroutine setup_partitions()
        implicit none
        integer:: n,k,m,sum_cN,np
        integer,dimension(ptdim):: a
        type(spartition),dimension(max_no_shells):: sMS
        logical:: more
        
        more=.false.
        n=no_cI+no_cN
        k=no_shells
        no_spart=0
        call setup_isospin()
        if (no_cN==abs(nmp)) then
            if (shellx) then
            do m=1,no_shells
                if (sign(1,sttz2(m))/=sign(1,nmp).or.no_cN==0) then
                    min_cSNo%shell(m)=0
                    max_cSNo%shell(m)=0
                end if
            end do
            end if
        end if
        
! NuShell SunShell        
        include 'ReadPartitions.FI'
! end
! NuShell@MSU
!        include 'Brown.FI'
! end                
        do

            call comp_next ( n, k, a, more )
            if (sum(a*sttz2) /= nmp) go to 10
            part=a
! check Shell linits min max 
            if (min_cSNo  > a) go to 10
            if (max_cSNo  < a) go to 10              
! check parity
            if (max_MSon<0) then
                do m=1,no_shells
                    if (sp(m)<0) then
                        if (btest(a(m),0)) go to 10
                    end if
                end do
            end if
            if (tParity /=0) then
                if (tParity /= cParity(part)) go to 10
            end if
! check Major shell limits sum over N (energy condition) and no in each MS min,max
            sum_cN=0
            if (cut) then
                tot_energy=0.0
                tot_energy=sum(a*spe)
                if (tot_energy>cut_energy ) go to 10
            end if
            do m=1,nMS
                sMS(m)=a*cMs_mask(m)
                sum_cN=sum_cN+psum(sMS(m))*cNMS(m)
                if (psum(sMS(m)) < min_MS_No%shell(m)) go to 10
                if (psum(sMS(m)) > max_MS_No%shell(m)) go to 10
            end do
            if (sum_cN < min_MSon) go to 10
            if (sum_cn > abs(max_MSon)) go to 10
! passed all limit tests
            no_spart=no_spart+1
            if (no_spart >  nopart) then
                print *, ' Increase parameter max_no_partitions in module Parameters,or'
                print *, ' SET PARTITIONS=I where I >>',no_spart
                stop
            end if
            spart(no_spart)=part
                    
10          if ( .not. more )  then
                exit
            end if

        end do
        
50      if (output_control > 1) then
            open(unit=40,file='partitions.txt',status='unknown',action='write',err=60)
            write (40,*,err=59) 'Partition Information'
            write (40,*,err=59)
            write (40,*,err=59) 'No of Shell Partitions',no_spart

            write (40,*,err=59) 'Shell partitions:'
            do ptn=1,no_spart
                write(40,'(a2,i5,a2)', advance='no',err=59) 'No',ptn,' >'
                write(40,'(20i3)',err=59) (spart(ptn)%shell(shn),shn=1,no_shells)
            end do
59          close(unit=40)
        end if
60      return                   
    end subroutine setup_partitions
    
    subroutine write_partitions()
        implicit none
        integer:: i
        
        open (unit=34, file=Nucleus//'.paa', status='unknown', action='write',err=10)
        write (34,'(i7)',err=9) no_spart
        do i=1,no_spart
            write (34,'(50i5)',err=9) spart(i)%shell(1:no_shells)
        end do
9       close (unit=34)
10      return        
    end subroutine write_partitions  
    
end module SetupPartitions

module SetupMvectors

use InputNuShell
use SetupPartitions
use Combine
use Counter

implicit none

type(m_vector),dimension(:),allocatable:: mvec
type(m_vector),dimension(:),allocatable:: vec
integer(kind=2),dimension(:),allocatable:: cM2_mvec,cTz2_mvec
integer,dimension(max_no_cJstates*2+6,max_no_cTstates*2+6):: no_cMTz2  
type(spartition):: maxJ,maxT
integer:: max_cM2,max_cTz2,max_J,Max_T,diff_J,diff_T,dJ
integer:: no_mvec,noprojmvec,nocombins,noallmvec,res,max_nJT=0,max_ibf=0
type(spartition):: max_count,count,current_counter
integer,dimension(max_no_cJstates*2,max_no_cTstates*2):: tot_ibf

contains

    subroutine setup_mvectors(part,pindx,dt)
        implicit none
        type(spartition),intent(in):: part
        type(m_vector):: mvecm
        integer,intent(in):: pindx,dt
        integer:: cM2,nptcle,cTz2
        integer,dimension(2)::opa
        logical:: nowr
        
1       max_cM2=0       ! maximum M*2 for partition
        max_cTz2=0      ! maximum T*2 for partition
        max_J=0         ! maximum total J*2 for partition
        max_T=0         ! maximum total T*2 for partition
        maxJ=0          ! maximum total J*2 per shell
        maxT=0          ! maximum total T*2 per shell
        max_count=0     ! maximum no of cominations per shell
        no_cMTz2=0      ! no of states with particular M*2 Tz*2
        no_mvec=0       ! no of mvectors
        
        if ((output_control > 3).and.(pindx == 1)) then
            nowr=.true.
            open(unit=42,file='Mvectors.txt',status='unknown',action='write',err=2)
            write (42,*) ' Mvectors for 2*M =',min_cJ2,'  2*Tz =',min_cT2
            go to 3
2           nowr=.false.
            stop            
        end if
3       if (output_control > 5) print *, ' pindx,no_spart',pindx,no_spart
        if (output_control > 3) then
                if (nowr) write(42,'(10i5)') part
        end if
        
        ! get all possible combinations for each shell
        call get_combinations(part,.false.)        
        ! calculate max_J
        ! calculate max_T
        do shn=1,no_shells
            max_J=max_J+maxJ%shell(shn)
            max_T=max_T+maxT%shell(shn)
        end do
        diff_J=max_J-min_cJ2    ! difference between minimum J*2 we require and maximum J*2 possible
        diff_T=max_T-min_cT2    ! difference between minimum T*2 we require and maximum T*2 possible
        
        ! get combinations again eliminating those not haveing large enough M2, Tz2
        if (dJ==0) call get_combinations(part,.true.)
        
        !Start counter. Counter counts through the different cominations in each shell.
        !it works like a digital counter increasing the first digit (shell1) from 1 to
        !max_count%shell(1)=ncomb(shell1), then resets shell1 digit to 1 and increments the
        !next digit (shell2) etc. It also computes cM2 and cTz2 the total M*2 and Tz*2 for this
        !set of digits.
        count=set_counter(max_count,part,cM2,cTz2,current_counter)
        do
        
        ! check M2 , Tz2, parity
            if (cM2 < min_cJ2-dJ) go to 10
            if (cTz2 < min_cT2) go to 10
            if (cM2 > max_cM2) max_cM2=cM2
            if (cTz2 > max_cTz2) max_cTz2=cTz2
            if (cTz2 > max_cT2+dt+2) go to 10
            if (cM2 > max_cJ2+2+2 ) go to 10
            if (cM2-min_cJ2+dJ+1 > 2*max_no_cJstates) then
                print *, ' Increase max_no_cJtates'
                stop
            end if
            if (cTz2-min_cT2+1 > 2*max_no_cTstates) then
                print *, ' Increase max_no_cTtates'
                stop
            end if
            if (cM2>=min_cJ2-dJ .and. cTz2>=min_cT2) & 
            no_cMTz2(cM2-min_cJ2+dJ+1,cTz2-min_cT2+1)=no_cMTz2(cM2-min_cJ2+dJ+1,cTz2-min_cT2+1)+1
            if (cTz2 == max_cT2+dt .and. cM2 == max_cJ2+2 .and. use_isospin) go to 10
            if (cTz2 > max_cT2+dt) go to 10
            if (cM2 > max_cJ2+2 ) go to 10

        ! M2 Tz2 accepted
        
            no_mvec=no_mvec+1
            if (no_mvec+1 > noallmvec) then
                print *, ' Increase parameter max_allmvector in module Parameters or'
                print *, ' SET TOTALMVECS=I where I >>',no_mvec+1
                stop
            end if
            cM2_mvec(no_mvec)=cM2
            cTz2_mvec(no_mvec)=cTz2
            cM2_mvec(no_mvec+1)=-1
            cTz2_mvec(no_mvec+1)=-1
            
        ! Create m-vector
            mvecm=int(0,ik)
            do shn=1,no_shells
                if (part%shell(shn) > 0) then
                    opa(1)=shell_word(shn)
                    do nptcle=1,part%shell(shn)
                        opa(2)=shell_offset(shn)+combin(nptcle,count%shell(shn),shn)-1
                        mvecm=opa+mvecm
                    end do
                end if
            end do
            mvec(no_mvec)=mvecm
        
            if ((output_control > 3).and.(cM2 == min_cJ2).and.(cTz2 == min_cT2) ) then
                if (nowr) call swrite(42,mvec(no_mvec))
            end if
                                    
10          call next(count,part,cM2,cTz2,current_counter)   !increment counter
            if (count <= 0) exit
        end do
        
        if ((output_control > 3).and.(pindx == no_spart)) then
            if (nowr) close(unit=42)                    
        end if
        if (output_control >  8) write(15,*) 'no_mvec,max_cm2,max_ctz2',&
                              no_mvec,max_cm2,max_ctz2 
    end subroutine setup_mvectors
    
    subroutine write_mvectors(part,pindx,base,baset,dt)
        implicit none
        type(spartition),intent(in):: part
        integer,intent(in):: dt,pindx,base,baset
        integer:: j,t,ibf,nmv,nJT,i,s
        integer,save:: k,jj,tt
        
        k=0
        jj=-1+dJ
        do j=min_cJ2,max_cJ2+2,2
        jj=jj+2
        tt=-1
        s=1
        if (j<0) s=-1
        do t=min_cT2,max_cT2+dt,2   
        tt=tt+2
        if (output_control > 10) write(15,*) ' j,t,pindx',j,t,pindx
        if (pindx==1) write(baset+k,'(i8)') no_spart
        
            ibf=0
            
            if (no_cMTz2(j-min_cJ2+dJ+1,t-min_cT2+1) /= 0 ) then
                do nmv=no_mvec,1,-1
                    if ((cM2_mvec(nmv) == j) .and. (cTz2_mvec(nmv) == t) & 
                        )then
                        ibf=ibf+1
                        if (ibf > noprojmvec) then
                            print *, ' Increase max_proj_mvec  in Parammeters or'
                            print *, ' SET MVECTORS=I where I >>',ibf
                            stop
                        end if 
                        if (output_control > 10) call swrite(15,mvec(nmv))
                        vec(ibf)= mvec(nmv)
                    end if
                    
                end do
                
           end if
           max_ibf=max(ibf,max_ibf)
           nJT=no_cMTz2(jj,tt)-no_cMTz2(jj+s*2,tt)  &
                        -no_cMTz2(jj,tt+2)+no_cMTz2(jj+s*2,tt+2)
           if (t<=max_cT2 .and. j<=max_cJ2 .and. nJT>max_nJT) max_nJT=nJT                        
           write(base+k,err=100) pindx,nJT,j,t,ibf,max_cM2,max_cTz2       
           if (output_control > 8) write(15,*) no_mvec,pindx,nJT,j,t,ibf,max_cM2,max_cTz2       
           write(base+k,err=100) part
           write(baset+k,'(100i4)',advance='no') part%shell(1:no_shells)
           write(baset+k,'(2i9)') nJT,ibf
           if (output_control > 8) write(15,*) part
           if (ibf /=0)  write(base+k,err=100) (vec(i)%fn(1:gwords()),i=1,ibf)
           tot_ibf(jj,tt)=tot_ibf(jj,tt)+ibf
           k=k+1
       end do
       end do
       return
100    print *, ' Error writing mvectors'
       stop
                    
    end subroutine write_mvectors     
    
    function cmvpt(mvec)
        implicit none
        type(spartition):: cmvpt
        type(m_vector):: mvec
        integer(kind=ik)::word
        integer::shn
        
        do shn=1,no_shells
            word=iand(shell_mask(shn),mvec%fn(shell_word(shn)))
            cmvpt%shell(shn)=countb(word)
        end do
    end  function cmvpt
    
    subroutine get_combinations(part,option)
        implicit none
        logical, intent(in):: option
        integer:: n,k,ncomb,nptcle,sM2,sTz2
        logical:: done,even_odd
        integer,dimension(2*(max_sj2+1)):: a,b
        type(spartition),intent(in):: part
        integer:: part_pos
        
        do shn=1,no_shells
            k=part%shell(shn)
            n=(st2(shn)+1)*(sj2(shn)+1)
            done=.true.
            ncomb=0
            a=0
            b=0
            if (k /= 0) then
                do
                    if (k > 1) then 
                        call comb_next ( done, n, k, a )
                        where (a /= 0) b=n-a+1
                    else
                        call pos_next (done, n, k, a )
                        where (a /= 0) b=n-a+1
                    end if
                    if ( done ) then
                        exit
                    end if
                    sM2=0
                    sTz2=0
                    if (st2(shn) == 0) then
                        do nptcle=1,k
                            sM2=sM2+sm2v(isj2(sj2(shn)),b(nptcle))
                            sTz2=sTz2+0
                        end do
                    else
                        do nptcle=1,k
                            even_odd=btest(b(nptcle),0)
                            if (even_odd) then
                                part_pos=b(nptcle)/2+1
                                sTz2=sTz2-1
                            else
                                part_pos=b(nptcle)/2
                                sTz2=sTz2+1
                            end if
                            sM2=sM2+sm2v(isj2(sj2(shn)),part_pos)
                        end do
                    end if
                    if (.not.option) then
                        ncomb=ncomb+1
                        if (ncomb.gt.nocombins) then
                            print *, ' Increase Parameter max_combine in module Parameters or'
                            print *, ' SET COMBINATIONS=I where I >> ', ncomb
                            stop
                        end if
                        combin(:,ncomb,shn)=b
                        combM2(ncomb,shn)=sM2   ! M*2 for this comination/shell
                        if (sM2>maxJ%shell(shn)) maxJ%shell(shn)=sM2
                        combT2(ncomb,shn)=sTz2  ! T*2 for this combination/shell
                        if (sTz2>maxT%shell(shn)) maxT%shell(shn)=sTz2
                        if (output_control > 12) write(15,*) ncomb,shn,k,combin(1:k,ncomb,shn)&
                                        ,combM2(ncomb,shn),combT2(ncomb,shn)
                    else
                    !check that M*2 is greater than or equal the minumum that might be used
                    !check that Tz*2 is greater than or equal the minumum that might be used
                        if (sM2 >= maxJ%shell(shn)-diff_J .and. sTz2 >= maxT%shell(shn)-diff_T) then
                        !we have a new valid combination/shell
                            ncomb=ncomb+1
                            combin(:,ncomb,shn)=b
                            combM2(ncomb,shn)=sM2
                            combT2(ncomb,shn)=sTz2
                            if (output_control > 12) write(15,*) ncomb,shn,k,combin(1:k,ncomb,shn)&
                                        ,combM2(ncomb,shn),combT2(ncomb,shn)
                        end if
                    endif
                end do
            end if
            if (ncomb == 0) ncomb=1
            max_count%shell(shn)=ncomb
        end do
        
    end subroutine get_combinations
    
end module SetupMvectors

!  NuBasis.f90 
!
!  FUNCTIONS:
!  NuBasis      - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: NuBasis
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program NuBasis
    
    use InputNuShell
    use OutputNuShell
    use SetupPartitions
    use SetupMvectors
    use Files
    use Extra
    use NZ
    
    implicit none

    ! Variables
    
    real:: t1,t2,et1,et2,ms
    integer:: hr,mn,sc
    integer:: pindx,s,j,t,jj,tt,dt,twrite,max_dim=0,error,max_tot_ibf=0
    integer,dimension(max_no_cJstates*2,max_no_cTstates*2):: tno_cJTstates,total_dim
    character(len=6):: inExt
    character(len=10):: btime
    character(len=8):: bdate
        character(len=10):: chnocombins,chnoallmvec,chnoprojmvec,challmvec='TOTALMVECS'
        character(len=10):: chnodim,chnonjt,chnoibf
        character(len=12):: chcombins='COMBINATIONS'
        character(len=8):: chprojmvec='MVECTORS'
        character(len=10):: chpart='PARTITIONS',chnopart
        
        total_dim=0
        
        res=get_env_v(chpart,chnopart)
        if (res /= 0) then
            read(chnopart,*) nopart
        else
            nopart=max_no_partitions
        end if
        allocate(spart(nopart))
        
        res=get_env_v(chcombins,chnocombins)
        if (res /=0 ) then
            read(chnocombins,*) nocombins
        else
            nocombins=max_combine
        end if
        allocate(combin(2*(max_sj2+1),nocombins,max_no_shells))
        allocate(combM2(nocombins,max_no_shells))
        allocate(combT2(nocombins,max_no_shells))
        
        res=get_env_v(challmvec,chnoallmvec)
        if (res /= 0) then
            read(chnoallmvec,*) noallmvec
        else
            noallmvec=max_allmvector
        end if
        allocate(mvec(noallmvec))
        allocate(cM2_mvec(noallmvec))
        allocate(cTz2_mvec(noallmvec))
        
        res=get_env_v(chprojmvec,chnoprojmvec)
        if (res /= 0) then
            read(chnoprojmvec,*) noprojmvec
        else
            noprojmvec=max_proj_mvec
        end if
        allocate(vec(noprojmvec))
    
    ! Body of NuBasis
    
    call output_header(6)
    call output_NuBasis(6)
    call input_nucleus(inExt)
    
    call output_welcome(6)    
    call input(inExt)
    call output_shells(6)    
    call check_shell_order()
    call setup_major_shells()
    call setup_shells()
    if (use_isospin) then
        dt=2
    else
        dt=0
    end if
    dJ=0
    if (min_cJ2<0) dJ=2
    call setup_mscheme()
    call setup_file_ext()
    call output_end_setup(6)
    
    call cpu_time(t1)
    call date_and_time(bdate,btime)
    call setup_partitions()
    if (output_control  > 0) then
        write(6,*)
        write (6,*) ' No of Shell Partitions',no_spart
        write(8,*)
    end if
    tno_cJTstates=0
    call open_files(no_spart,'.nba',100,dt)
    call open_files_form(no_spart,'.npt',200,dt)
    if (output_control > 1) open(unit=15,file=Nucleus//'NuBasis.txt',action='write')
    
    do pindx=1,no_spart
        call setup_mvectors(spart(pindx),pindx,dt)
        jj=-1+dJ
        do j=min_cJ2,max_cJ2+2,2
            jj=jj+2
            s=1
            if (j<0) s=-1
            tt=-1
            do t=min_cT2,max_cT2+dt,2
                tt=tt+2
                tno_cJTstates(jj,tt)=tno_cJTstates(jj,tt)+no_cMTz2(jj,tt)-no_cMTz2(jj+s*2,tt)  &
                                    -no_cMTz2(jj,tt+2)+no_cMTz2(jj+s*2,tt+2)
            end do
        end do
        call write_mvectors(spart(pindx),pindx,100,200,dt)
    end do
    max_tot_ibf=maxval(tot_ibf)
    
    if (output_control > 0) then
        print * ,' No of states for each 2*J'
        print '(a10,10i8)', '       2*J',(j,j=min_cJ2,min_cJ2+18-dJ,2)
    end if
    open(unit=1,file=Nucleus//'.JT' ,action='write')   
        tt=-1
        do t=min_cT2,max_cT2+dt,2
            twrite=t
            if (.not.use_isospin) twrite=nmp
            tt=tt+2
            if (t<=max_cT2) max_dim=max(max_dim,maxval(tno_cJTstates(1:max_cJ2-min_cJ2+dJ+1,tt)))
            if (output_control > 0) print '(a4,i6,10i8)', '2*T=',twrite,(tno_cJTstates(j,tt),j=1+dJ,19+dJ,2)
            write(1,'(100i7)')(tno_cJTstates(j,tt),j=1,max_cJ2-min_cJ2+2,2)
        end do
        
    deallocate(spart)
    deallocate(combin)
    deallocate(combM2)
    deallocate(combT2)
    deallocate(cM2_mvec)
    deallocate(cTz2_mvec)
    deallocate(vec)
    write(chnoallmvec,'(i10)') max_tot_ibf
    write(chnodim,'(i10)') max_dim
    write(chnonjt,'(i10)') max_njt
    write(chnoibf,'(i10)') max_ibf
    write(chnopart,'(i10)') no_spart
    chnoallmvec=adjustl(chnoallmvec)
    chnodim=adjustl(chnodim)
    chnonjt=adjustl(chnonjt)
    chnopart=adjustl(chnopart)
    chnoibf=adjustl(chnoibf)
    open (unit=1,file='SetNuShell.bat',action='write',err=10)
        write(1,'(a25)') 'SET TOTALMVECS='//chnoallmvec
        write(1,'(a24)') 'SET DIMENSION='//chnodim
        write(1,'(a18)') 'SET NJT='//chnonjt
        write(1,'(a23)') 'SET MVECTORS='//chnoibf
        write(1,'(a25)') 'SET PARTITIONS='//chnopart
    close (unit=1)
10  read(btime,'(i2,i2,i2,f4.3)') hr,mn,sc,ms
    et1=hr*3600.+mn*60.+sc*1. +ms
    call date_and_time(bdate,btime)
    call cpu_time(t2)
    read(btime,'(i2,i2,i2,f4.3)') hr,mn,sc,ms
    et2=hr*3600.+mn*60.+sc*1. +ms
    et1=et2-et1
    if (et1 < 0.0) et1=et1+24.*3600.
    t1=t2-t1
    call output_completed(6,'NuProj     ')
    call output_time('NuBasis',t1,et1)
    print *

    
    end program NuBasis

    subroutine output_NuBasis(dev)
    
    use InputNuShell
    use OutputNuShell
    use SetupPartitions
    use SetupMvectors
    use Files

    implicit none
    integer,intent(in):: dev
    integer:: status
    character(len=14)::buffer
    call get_cmd_arg(1,buffer,status) 
    if (status /= -1) return
    
    write(dev,*) ' NuBasis Program'
    write(dev,*)
    write(dev,*) ' Compile time Parameters:'
    write(dev,*)
    write(dev,*) ' Machine Parameters:'
    write(dev,*) ' No of bits per word',nbits,' Memory Access Width Bytes', nbytes
    write(dev,*)
    write(dev,*) ' Shell Model Parameters:'
    write(dev,*) ' Maximum no of shells',max_no_shells
    write(dev,*) ' Maximum no of words with',nbits,' bits is',nwords, '    for each m-scheme state'
    write(dev,*) ' Maximum no of oxbash words with',nbits,' bits is',nwords, '    for each m-scheme state'
    write(dev,*) ' Maximum value of 2*j for a shell',max_sj2
    write(dev,*) ' Maximum no of shell and isospin shell partitions',max_no_partitions
    write(dev,*) ' Maximum no of arrangements of nucleons in each shell (max_combine)',max_combine
    write(dev,*) ' Maximum no of m-vectors per partition all m',max_allmvector
    write(dev,*) ' Maximum no of m-vectors per partition single m',max_proj_mvec
    write(dev,*) ' Maximum no J,T values',max_no_cJstates,max_no_cTstates
    write(dev,*)

    end subroutine output_NuBasis