! This file contains all the module for the generic pMatrix
! program for NuShell, NuShell@MSU and SunShell.

module Matrix

use ProjParameters
use MatrixParameters
use ShellsX
use Mscheme
use Mvector
use InputNuShell
use OperParameters
use Finds
use Partition
use OpPartition

implicit none
integer:: nomdim,noprojmvec,noindex,noop,nopart,max_sgmt
integer,dimension(nwords*nbits,nwords*nbits)::  op_start,op_dim
type(m_vector),dimension(:),allocatable:: mvec_lhs
type(m_vector),dimension(:),allocatable:: mvec_rhs
integer:: dim_lhs,dim_rhs,no_spart,no_me,index,use_wordm
real(kind=rkw),dimension(:,:),allocatable:: coef_lhs
integer,dimension(:,:),allocatable::index_lhs
real(kind=rkw),dimension(:,:),allocatable:: coef_rd
real(kind=rm),dimension(:,:),allocatable:: mat_sgmt
integer,dimension(:),allocatable::  part_lu
integer,dimension(:),allocatable::  part_offset
integer,dimension(:),allocatable:: dim_lu
integer,dimension(:),allocatable:: dim_part
type(spartition),dimension(:),allocatable:: spart
type(opartition),dimension(:),allocatable:: me_part
real(kind=ri),dimension(:),allocatable:: me
character(len=1):: xtra
real,save:: x,y
logical:: defin

contains

    subroutine mat_op(baseb,basei,baser,baseo,k,k1)
        implicit none
        integer,intent(in):: baseb,basei,baser,baseo,k,k1
        type(m_vector):: mvec_copy1,mvec_copy2,mvecam,mvecam2,mvec1
        type(m_vector):: mvec_copy3,mvec_copy4,mvecap,mvecap2
        integer:: dnr,nb1,nb2,iw1,iw2,findm,np,npr,npst,npsp
        integer:: nb3,nb4,iw3,iw4,j1,j2,j3,j4,jmn12,jmx12,jmn34,jmx34
        integer:: m1,m2,m3,m4,tz1,tz2,tz3,tz4,n,i
        integer:: sh1,sh2,sh3,sh4,shw1,shw2,shw3,shw4,iwmn,iwmx
        real(kind=ri):: Vm,dot
        real(kind=4):: phase1,ph12,ph123,ph1234,pphase1,pphase2,pphase3,pphase4
        integer:: op1,op2,op3,op4,use_wordm,opn1,opn2,opn3,opn4
        integer,dimension(2)::aop1,aop2,aop3,aop4
        type(spartition):: part_lhs,part_rhs,part_op,op_op1,op_op2,op_op3,op_op4
        logical:: warn=.false.,zero,main_warn=.false.
        character(len=10):: ctime,ptime
        character(len=8):: date,pdate
        integer(kind=ik):: shmask1,shmask2,shmask3,shmask4
        INTEGER:: OMP_get_thread_num,thrn,OMP_get_max_threads,nthr
        integer(kind=ik):: ibclr,iand,not
	
        if (hpart .and. np_cntl  == 0 ) then
            print *, ' hpart is true  and np_cntl = 0, for  hpart true np_cntl /=0 '
            stop
        end if
        use_wordm=gwords()    
        call read_rhs_basis(baser,k)

        nthr=OMP_get_max_threads()
        allocate(coef_lhs(noindex,nthr))
        allocate(index_lhs(noindex,nthr))
        !outer loop over m_vectors, inner loops over nucleons
        do np=1,no_spart
        if (output_control > 0) then
            call date_and_time(date,ctime)
            write(10,*) date,' ',ctime,' npart left',np
        end if
        call read_lhs_basis(baseb,k1)
        if (allocated(coef_rd)) deallocate(coef_rd)
        allocate(coef_rd(dim_lhs,dim_part(np)))
        call read_coef_rd(basei,k)
        part_lhs=spart(np)
        if (dim_part(np) < 1) then
            deallocate(coef_rd)
            go to 50
        end if
        if (np_cntl > 0 .and. .not.hpart) then
            npst=np
            npsp=no_spart
        end if
        if (np_cntl < 0 .and. .not.hpart) then
            npst=1
            npsp=np
        end if
        if (np_cntl == 0 .or. hpart) then
            npst=1
            npsp=no_spart
        end if
        do npr=npst,npsp
        
        if (output_control > 3) then
            call date_and_time(pdate,ptime)
            write(10,*) pdate,' ',ptime,' npart right',npr
        end if
        if (dim_part(npr) <  1) go to 45
        zero=.true.
        part_rhs=spart(npr)
        if (cParity(part_rhs) /= cParity(part_lhs)) go to 45
        part_op=part_lhs-part_rhs
        if (hpart .and. hop(part_op)) go to 45
        if (part_op < -2 .or. part_op > 2 .or. analyse_part_op(part_op,sh1,sh2,sh3,sh4)) then
            if (output_control > 6) then
                write(60,*) ' Eliminated partitions'
                write(60,*)  part_lhs
                write(60,*)  part_rhs
                write(60,*)  part_op
            end if
            go to 45
        end if
        zero=.false.
        allocate(mat_sgmt(dim_part(npr),dim_part(np)))
        defin=.false.
        if (sh1/=0 .and. sh2/=0 .and. sh3/=0 .and. sh4/=0) then
            defin=.true.
            shmask1=shell_mask(sh1)
            shw1=shell_word(sh1)
            shmask3=shell_mask(sh3)
            shw3=shell_word(sh3)
            shmask2=shell_mask(sh2)
            shw2=shell_word(sh2)
            shmask4=shell_mask(sh4)
            shw4=shell_word(sh4)
        end if
!$OMP PARALLEL DEFAULT(PRIVATE) &
!$OMP SHARED(part_offset,dim_part,mvec_rhs,mvec_lhs,coef_rd,&
!$OMP mat_sgmt,mvec_jz_tz,word_mask,shn_lu,sm2_lu,sj2_lu,&
!$OMP stz2_lu,shmask1,shmask2,shmask3,shmask4,shw1,shw2,&
!$OMP shw3,shw4,defin,nshell_mask,use_wordm,npr,np,part_op,&
!$OMP op_dim,op_start,me,me_part,no_me,&
!$OMP dim_lhs,max_sj2u,noindex,index_lhs,coef_lhs) 
!$OMP DO
        do dnr=1+part_offset(npr),part_offset(npr)+dim_part(npr)
            thrn=OMP_get_thread_num()+1
            index=0
            mvec_copy1=mvec_rhs(dnr)
            mvec1=mvec_rhs(dnr)
            iwmn=1
            iwmx=use_wordm
            if (defin) then
                iwmx=shw1
                iwmn=shw1
            end if
            do iw1=iwmn,iwmx
            aop1(1)=iw1
            opn1=(iw1-1)*nbits
            if (defin) mvec_copy1%fn(iw1)=iand(shmask1,mvec_copy1%fn(iw1))
            pphase1=mphase(iw1,mvec1)
                do
                    if (mvec_copy1%fn(iw1) == 0) exit
                    nb1=next_bit(mvec_copy1%fn(iw1))
                    mvec_copy1%fn(iw1)=ibclr(mvec_copy1%fn(iw1),nb1-1)
                    op_op1=0
                    op_op1%shell(shn_lu(nb1,iw1))=-1
                    if ( test1(part_op-op_op1)) then 
                        mvec_copy1%fn(iw1)=iand(mvec_copy1%fn(iw1), &
                                                   nshell_mask(shn_lu(nb1,iw1)))
                        go to 35
                    end if
                    m1=sm2_lu(nb1,iw1)
                    tz1=stz2_lu(nb1,iw1)
                    j1=sj2_lu(nb1,iw1)
                    op1=opn1+nb1
                    aop1(2)=nb1
                    mvecam=aop1-mvec1
                    phase1=wphase(nb1,mvec1%fn(iw1))*pphase1
                    if (defin) then
                        mvec_copy2=mvecam
                    else
                        mvec_copy2=mvec_copy1
                    end if
                    iwmn=iw1
                    iwmx=use_wordm
                    if (defin) then
                        iwmn=shw2
                        iwmx=shw2
                    end if
                    do iw2=iwmn,iwmx
                    aop2(1)=iw2
                    opn2=(iw2-1)*nbits
                    if (defin) mvec_copy2%fn(iw2)=iand(shmask2,mvec_copy2%fn(iw2))
                    pphase2=mphase(iw2,mvecam)
                        do
                            if (mvec_copy2%fn(iw2) == 0) exit
                            nb2=next_bit(mvec_copy2%fn(iw2))
                            mvec_copy2%fn(iw2)=ibclr(mvec_copy2%fn(iw2),nb2-1)
                            op2=opn2+nb2
                            if (op_dim(op1,op2) == 0) go to 30
                            op_op2=op_op1
                            op_op2%shell(shn_lu(nb2,iw2))=op_op2%shell(shn_lu(nb2,iw2))-1
                            if ( test2(part_op-op_op2)) then
                                mvec_copy2%fn(iw2)=iand(mvec_copy2%fn(iw2), &
                                                   nshell_mask(shn_lu(nb2,iw2)))
                                go to 30
                            end if
                            m2=sm2_lu(nb2,iw2)
                            tz2=stz2_lu(nb2,iw2)
                            j2=sj2_lu(nb2,iw2)
                            jmx12=j1+j2
                            jmn12=abs(j1-j2)
                            aop2(2)=nb2
                            mvecam2=aop2-mvecam
                            ph12=phase1*wphase(nb2,mvecam%fn(iw2))*pphase2
                            mvec_copy3=ph(mvecam2,word_mask)
                            iwmn=1
                            iwmx=use_wordm
                            if (defin) then
                                iwmn=shw3
                                iwmx=shw3
                            end if
                            do iw3=iwmn,iwmx
                            aop3(1)=iw3
                            opn3=(iw3-1)*nbits
                            if (defin) mvec_copy3%fn(iw3)=iand(shmask3,mvec_copy3%fn(iw3))
                            pphase3=mphase(iw3,mvecam2)
                                do
                                    if (mvec_copy3%fn(iw3) == 0) exit
                                    nb3=next_bit(mvec_copy3%fn(iw3))
                                    mvec_copy3%fn(iw3)=ibclr(mvec_copy3%fn(iw3),nb3-1)
                                    op_op3=op_op2
                                    op_op3%shell(shn_lu(nb3,iw3))=op_op3%shell(shn_lu(nb3,iw3))+1
                                    if ( test3(part_op-op_op3)) then
                                        mvec_copy3%fn(iw3)=iand(mvec_copy3%fn(iw3), &
                                                           nshell_mask(shn_lu(nb3,iw3)))
                                        go to 20
                                    end if
                                    m3=sm2_lu(nb3,iw3)
                                    tz3=stz2_lu(nb3,iw3)
                                    j3=sj2_lu(nb3,iw3)
                                    m4=m1+m2-m3
                                    tz4=tz1+tz2-tz3
                                    if (abs(tz4) > 1) go to 20
                                    if (abs(m4) > max_sj2u) go to 20
                                    op3=opn3+nb3
                                    aop3(2)=nb3
                                    mvecap=aop3+mvecam2
                                    ph123=ph12*wphase(nb3,mvecam2%fn(iw3))*pphase3
                                    if (defin) then
                                        mvec_copy4=ph(mvecap,word_mask)
                                    else
                                        mvec_copy4=mvec_copy3
                                    end if
                                    mvec_copy4=mand(mvec_jz_tz(m4,tz4),mvec_copy4)
                                    if (mvec_copy4 == int(0,ik)) go to 20
                                    iwmn=iw3
                                    iwmx=use_wordm
                                    if (defin) then
                                        iwmn=shw4
                                        iwmx=shw4
                                    end if
                                    do iw4=iwmn,iwmx
                                    aop4(1)=iw4
                                    opn4=(iw4-1)*nbits
                                    if (defin) mvec_copy4%fn(iw4)=iand(shmask4,mvec_copy4%fn(iw4))
                                    pphase4=mphase(iw4,mvecap)
                                        do
                                            if (mvec_copy4%fn(iw4) == 0) exit
                                            nb4=next_bit(mvec_copy4%fn(iw4))
                                            mvec_copy4%fn(iw4)=ibclr(mvec_copy4%fn(iw4),nb4-1)
                                            op4=opn4+nb4
                                            if (op_dim(op3,op4) == 0) go to 15
                                            op_op4=op_op3
                                            op_op4%shell(shn_lu(nb4,iw4))=op_op4%shell(shn_lu(nb4,iw4))+1
                                            if (part_op-op_op4 /= 0) then
                                                mvec_copy4%fn(iw4)=iand(mvec_copy4%fn(iw4), &
                                                           nshell_mask(shn_lu(nb4,iw4)))
                                                go to 15
                                            end if
                                            j4=sj2_lu(nb4,iw4)
                                            jmx34=j3+j4
                                            jmn34=abs(j3-j4)
                                            if (max(jmn12,jmn34) > min(jmx12,jmx34)) go to 15
                                                aop4(2)=nb4
                                                mvecap2=aop4+mvecap
                                                ph1234=ph123*wphase(nb4,mvecap%fn(iw4))*pphase4
                                                Vm=findVm(op1,op2,op3,op4)
                                                if (Vm  == real(0,ri)) then
                                                    go to 15
                                                end if
                                                findm=bfind(mvec_lhs,dim_lhs,mvecap2)
                                                if (findm /=0 ) then
                                                    index=index+1
                                                    if (index > max_index) then
                                                        print *, ' Increase max_index in MatrixParams or'
                                                        print *, ' SET INDEXX=I where I >>',index
                                                        stop
                                                    end if
                                                    coef_lhs(index,thrn)=Vm*ph1234
                                                    index_lhs(index,thrn)=findm
                                                else
                                                     print *,part_rhs%shell(1:no_shells)
                                                     print *,part_lhs%shell(1:no_shells)
                                                     print *,part_op%shell(1:no_shells)
                                                     print *,dim_rhs,dim_lhs
                                                     print *,index,dnr,Vm,ph1234
                                                     print *,k,k1
                                                     call swrite(6,mvec_rhs(dnr))
                                                     call swrite(6,mvecap2)
                                                     do shn=1,dim_rhs
                                                     call swrite(6,mvec_rhs(shn))
                                                     end do
                                                     do shn=1,dim_lhs
                                                     call swrite(6,mvec_lhs(shn))
                                                     end do
                                                     stop
                                                end if
15                                       end do
                                     end do
20                               end do
                            end do                                        
30                      end do
                    end do
35              end do
            end do
! do overlaps 
40       do n=1, dim_part(np)
            dot=real(0,rm)
                if (index /=0) then
                do i=1,index
                    dot=dot+coef_rd(index_lhs(i,thrn),n)*coef_lhs(i,thrn)
                end do
                end if
            mat_sgmt(dnr-part_offset(npr),n)=dot
         end do

         end do  !dnr
!$OMP END PARALLEL         
        max_sgmt=max_sgmt+1
        write(baseo+k,err=100)dim_part(np),dim_part(npr),np, npr,zero,x
        write(baseo+k,err=100)mat_sgmt(1:dim_part(npr),1:dim_part(np))
        if (output_control > 6) then
            write(60,*)dim_part(np),dim_part(npr),np,npr, zero,x
            write(60,*)mat_sgmt(1:dim_part(npr),1:dim_part(np))
        end if
         
         deallocate(mat_sgmt)
45       end do         !npr
         write(baseo+k,err=100) -1,-1,-1,-1,.true.,x         
         if (output_control > 6) write(60,*) -1,-1,-1,-1,.true.,x         

         deallocate(coef_rd)
50       end do            !np
         write(baseo+k,err=100) -2,-2,-2,-2,.true.,x         
         if (output_control > 6) write(60,*) -2,-2,-2,-2,.true.,x         
         if (main_warn) then
            print *, 'Warning some t.b. mcheme matrix elements were not found.'
         end if
         deallocate(coef_lhs)
         deallocate(index_lhs)
         return
100      print *, ' Error writing matrix'
         stop         

    end subroutine mat_op
    
    subroutine read_lhs_basis(base,k)
        implicit none
        integer,intent(in):: base,k
        integer:: nbuf,pindx,nJT,j,t,Jmax,Tmax,i
                
            read(base+k,err=10,end=11) pindx,nJT,j,t,nbuf,Jmax,Tmax      
            read(base+k,err=10,end=11) ! part
            if (nbuf > max_proj_mvec) then
                print * ,' Increase max_proj_mvec in Module Parameters or'
                print * ,' SET MVECTORS=I  where I >> ',nbuf
                stop
            end if
            if (nbuf > 0) read(base+k,err=10,end=11) (mvec_lhs(i)%fn(1:use_wordm),i=1,nbuf)

            dim_lhs=nbuf
            return
10          print *, ' Error reading mvectors'
            stop            
11          print *, ' EOF reading mvectors'
            stop            
    end subroutine read_lhs_basis
    
    subroutine read_coef_rd(basei,k)
        implicit none
        integer,intent(in):: basei,k
        integer:: pindx,nJT,j,t,dim,ngd
        
            read(basei+k,err=10,end=11) pindx,nJT,j,t,dim,y
            if (y/=x) then
                print *, ' Error: .prj and .jba files not consistent',x,y
                stop
            end if       
            read(basei+k,err=10,end=11)    !part
            if (nJT >0) then
            do ngd=1,nJT
                read(basei+k,err=10,end=11) coef_rd(1:dim,ngd)
            end do
            end if
            return
10          print *, ' Error reading prj file '
            stop            
11          print *, ' EOF reading prj file '
            stop            
            
    end subroutine read_coef_rd     
    
    function analyse_part_op(part,sh1,sh2,sh3,sh4)
        implicit none
        logical:: analyse_part_op
        type(spartition),intent(in)::part
        integer,intent(out)::sh1,sh2,sh3,sh4
        integer:: ones,twos,mones,mtwos,i,ps

        analyse_part_op=.false.
        sh1=0
        sh2=0
        sh3=0
        sh4=0
        
        if (part == 0) then
             return
        end if
        ones=0
        twos=0
        mones=0
        mtwos=0
        
        do i=1,no_shells
            ps=part%shell(i)
            if (ps/=0) then
            if (ps == 1) then
                ones=ones+1
                if (ones==1) sh3=i
                if (ones==2) sh4=i
            end if
            if (ps == -1) then
                mones=mones+1
                if (mones==1) sh1=i
                if (mones==2) sh2=i
            end if    
            if (ps == 2) then
                twos=twos+1
                if (twos==1) then
                    sh3=i
                    sh4=i
                end if
            end if
            if (ps == -2) then
                mtwos=mtwos+1
                if (mtwos==1) then
                    sh1=i
                    sh2=i
                end if
            end if
            end if
        end do
        if (ones > 2 .or. mones  > 2 .or. twos > 1 .or. mtwos  >  1) then
            analyse_part_op=.true.
            return
        end if
        if (ones - mones + 2*twos - 2*mtwos  /=  0) then
            analyse_part_op=.true.
            return
        end if
                
        if  (ones == 1 .and. mones == 1 .and. twos == 0 .and. mtwos == 0) then
            return
        end if
        if  (ones == 2 .and. mones == 2 .and. twos == 0 .and. mtwos == 0) then
            return
        end if
        if  (ones == 0 .and. mones == 2 .and. twos == 1 .and. mtwos == 0) then
            return
        end if
        if  (ones == 2 .and. mones == 0 .and. twos == 0 .and. mtwos == 1) then
            return
        end if
        if  (ones == 0 .and. mones == 0 .and. twos == 1 .and. mtwos == 1) then
            return
        end if
        analyse_part_op=.true.
    end function analyse_part_op 

    pure function test1(part)
        implicit none
        logical:: test1
        type(spartition),intent(in)::part
        integer:: ones,twos,mones,i,ps
        
        test1=.false.
        if (defin) return
        
        test1=.true.

        if (part == 0) return

        ones=0
        twos=0
        mones=0

        do i=1,no_shells
            ps=part%shell(i)
            if (ps/=0) then
            if (ps > 2) return  
            if (ps <= -2) return
            if (ps == 1) then
                ones=ones+1
                if (ones > 2) return
            end if
            if (ps == -1) then
                mones=mones+1
                if (mones > 1) return
            end if  
            if (ps == 2) then
                twos=twos+1
                if (twos > 1) return
            end if
            end if
        end do

        test1=.false.
        
    end function test1 


    pure function test2(part)
        implicit none
        logical:: test2
        type(spartition),intent(in)::part
        integer:: ones,twos,i,ps

        test2=.false.
        if (defin) return
        
        test2=.true.
        if (part == 0) return
        
        ones=0
        twos=0
        
        do i=1,no_shells
            ps=part%shell(i)
            if (ps/=0) then
            if (ps <= -1) return
            if (ps > 2) return
            if (ps == 1) then
                ones=ones+1
                if (ones > 2) return
            end if
            if (ps == 2) then
                twos=twos+1
                if (twos > 1) return
            end if
            end if
        end do

        test2=.false.
        
    end function test2 

    pure function test3(part)
        implicit none
        logical:: test3
        type(spartition),intent(in)::part
        integer:: ones,i,ps

        test3=.false.
        if (defin) return
        
        test3=.true.
        
        if (part == 0) return
        
        ones=0
        
        do i=1,no_shells
            ps=part%shell(i)
            if (ps/=0) then
            if (ps <= -1) return
            if (ps >= 2) return
            if (ps == 1) then
                ones=ones+1
                if (ones > 1) return
            end if
            end if
        end do
  
        test3=.false.
              
    end function test3 

    pure function  findVm(a1,a2,a3,a4)
        implicit none
        real(kind=ri)::findVm
        integer,intent(in):: a1,a2,a3,a4
        type(opartition):: part
        integer:: findp
        if (op_dim(a3,a4) == 0) then
            findVm=real(0,rm)
            return
        end if
        part%op(4)=a4
        part%op(3)=a3
        part%op(2)=a2
        part%op(1)=a1
        findp=bfind(me_part,op_dim(a3,a4),part,op_start(a3,a4)-1)
        if (findp /=0) then
            findVm=-me(findp)
            return
        else
            findVm=real(0,rm)
            return
        end if
        
    end function findVm
    
    subroutine read_rhs_basis(base,k)
        implicit none
        integer,intent(in):: base,k
        integer::pindx,ngood,j,t,dim,np,i
        use_wordm=gwords()
        dim_rhs=0
        do np=1,no_spart        
            read(base+k,err=10,end=11) pindx,ngood,j,t,dim,x
            read(base+k,err=10,end=11) spart(np)
            if (dim_lhs+ngood > nomdim) then
                print * ,' Increase max_total_JT in module MatrixParameters or'
                print * ,' SET DIMENSION=I where I=dimension of matrix and is >>',dim_lhs+ngood
                stop
            end if
            part_offset(pindx)=dim_rhs
            dim_part(pindx)=ngood
            if (ngood > 0) then
                read(base+k,err=10,end=11) (mvec_rhs(i)%fn(1:use_wordm),i=dim_rhs+1,ngood+dim_rhs)
                part_lu(dim_rhs+1:ngood+dim_rhs)=np
                dim_lu(dim_rhs+1:ngood+dim_rhs)=ngood
            end if
            dim_rhs=dim_rhs+ngood
        end do
        close(unit=base+k)
        return
10      print *, ' Error reading nba file'
        stop        
11      print *, ' EOF reading nba file'
        stop        
    end subroutine 
    
end module Matrix

module MatrixSupport

use Parameters
use ProjParameters
use MatrixParameters
use Matrix

implicit none

contains

    subroutine read_interaction(inExt)
        implicit none
        character(len=*):: inExt
        integer:: opst
        op_start=0
        op_dim=0
        no_me=1
        open(unit=20, file=Nucleus//xtra//inExt, form='unformatted',err=10)
        do
            read(20,err=20,end=21) me_part(no_me),me(no_me)
            if (me_part(no_me)%op(1) < 0) exit
            if (no_me > 1) then
                if (me_part(no_me) <= me_part(no_me-1)) then
                    print *, 'me order error'
                    print *, me_part(no_me-1)
                    print *, me_part(no_me) 
                end if
            end if
            opst=op_start(me_part(no_me)%op(3),me_part(no_me)%op(4))
            if (opst == 0 ) then
                op_start(me_part(no_me)%op(3),me_part(no_me)%op(4))=no_me
                op_dim(me_part(no_me)%op(3),me_part(no_me)%op(4))=1
            else
                op_dim(me_part(no_me)%op(3),me_part(no_me)%op(4))=  &
                op_dim(me_part(no_me)%op(3),me_part(no_me)%op(4))+1
            end if             
            no_me=no_me+1
            if (no_me > noop) then
                print * ,' Increase max_op in module MatrixParametere or'
                print * ,' SET NOP=I wheree I >>',no_me
                stop
            end if
        end do
        no_me=no_me-1
        close(unit=20)
        return
10      print *, ' Error opening ', Nucleus//xtra//inExt
        stop       
20      print *, ' Error reading ', Nucleus//xtra//inExt
        stop       
21      print *, ' EOF reading ', Nucleus//xtra//inExt
        stop       
                
    end subroutine read_interaction
    
end module MatrixSupport

!  NuMatrix.f90 
!
!  FUNCTIONS:
!  Numatrix      - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: NuMatrix
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program NuMatrix
    
    use MatrixParameters
    use OperParameters
    use ProjParameters
    use InputNuShell
    use OutputNuShell
    use ShellsX
    use Mscheme
    use Files
    use MatrixSupport
    use Extra
    
    implicit none

    ! Variables
    
    real:: t1,t2,et1,et2,ms
    integer:: hr,mn,sc

    integer:: j,t,k,dt,k1,res,maxsgmt=0
    
    character(len=6):: inExt
    character(len=10):: btime
    character(len=8):: bdate
    character(len=10)::chnoop,chnomdim,chnoprojmvec,chnoindex,chnopart,chpart='PARTITIONS',chnosgmt
    character(len=9)::chmdim='DIMENSION'
    character(len=8)::chprojmvec='MVECTORS'
    character(len=6)::chindex='INDEXX'
    character(len=3)::chop='NOP'
    
    ! Body of NuMatrix
    
    res=get_env_v(chop,chnoop)
    if (res /= 0 ) then
        read(chnoop,*) noop
        noop=noop+2
    else
        noop=max_op-2
    end if
    allocate(me(noop))
    allocate(me_part(noop))
    
    res=get_env_v(chmdim,chnomdim)
    if (res /= 0 ) then
        read(chnomdim,*) nomdim
    else
        nomdim=max_total_JT
    end if
    allocate(part_lu(nomdim))
    allocate(dim_lu(nomdim))
    allocate(mvec_rhs(nomdim))
    
    res=get_env_v(chprojmvec,chnoprojmvec)
    if (res /= 0 ) then
        read(chnoprojmvec,*) noprojmvec
    else
        noprojmvec=max_proj_mvec
    end if
    allocate(mvec_lhs(noprojmvec))
    
    res=get_env_v(chpart,chnopart)
    if (res /= 0) then
        read(chnopart,*) nopart
    else
        nopart=max_no_partitions
    end if
    allocate(part_offset(nopart))
    allocate(dim_part(nopart))
    allocate(spart(nopart))
    
    res=get_env_v(chindex,chnoindex)
    if (res /= 0 ) then
        read(chnoindex,*) noindex
    else
        noindex=max_index
    end if
    
    call output_header(6)
    call output_NuMatrix(6)
    
    call input_nucleus(inExt,xtra)

    call output_welcome(6)    
    call input(inExt)

    call output_shells(6)    
    call check_shell_order()
    call setup_shells()
    if (use_isospin) then
        dt=2
    else
        dt=0
    end  if
    call setup_mscheme()
    call setup_lu()
    call setup_file_ext()
    call read_interaction('.op')
    call output_end_setup(6)        
    call cpu_time(t1)
    call date_and_time(bdate,btime)
    if (output_control > 6) open(unit=60,file=Nucleus//xtra//'mat.txt',action='write')
    open(unit=99,file=Nucleus//xtra//'misVm.txt',action='write')
    call open_files_read(no_spart,'.nba',100,dt)
    if (no_spart>nopart) then
        print *, ' Increase max_no_partitions in Parameters or'
        print *, ' SET PARTITIONS=I where I >>',no_spart
        stop
    end if
    call open_files_read(no_spart,'.prj',200,dt)
    call open_files_read(no_spart,'.Jba',300,dt)
    call open_files(no_spart,'.mtx',500,dt,-1,xtra)
    if (output_control > 0) open(unit=10,file=Nucleus//'NuMatrix'//xtra//'.txt',action='write')
    k=0
    k1=0
    do j=min_cJ2,max_cJ2+2,2
        do t=min_cT2,max_cT2+dt,2
            if (j <= max_cJ2  .and. t <= max_cT2) then
                if (output_control > 0) write(10,*) ' 2*J, 2*T',j,t
                if (output_control > 6) write(60,*) ' 2*J, 2*T',j,t
                max_sgmt=0
                call mat_op(100,200,300,500,k,k1)
                maxsgmt=max(maxsgmt,max_sgmt)
            end if
            if (j <= max_cJ2  .and. t <= max_cT2) k=k+1
            k1=k1+1
        end do
    end do
    deallocate(me)
    deallocate(me_part)
    deallocate(part_lu)
    deallocate(dim_lu)
    deallocate(mvec_rhs)
    deallocate(mvec_lhs)
    deallocate(spart)
    deallocate(part_offset)
    deallocate(dim_part)
    write(chnosgmt,'(i10)') maxsgmt
    chnosgmt=adjustl(chnosgmt)
    open (unit=1, file='SetLanczos.bat', action='write')
        write(1,'(a23)') 'SET SEGMENTS='//chnosgmt
    close (unit=1) 
    read(btime,'(i2,i2,i2,f4.3)') hr,mn,sc,ms
    et1=hr*3600.+mn*60.+sc*1. +ms
    call date_and_time(bdate,btime)
    call cpu_time(t2)
    read(btime,'(i2,i2,i2,f4.3)') hr,mn,sc,ms
    et2=hr*3600.+mn*60.+sc*1. +ms
    t1=t2-t1
    et1=et2-et1
    if (et1<0.0) et1=et1+24.*3600.
    call output_completed(6,'NuLanczos ')    
    call output_time('NuMatrix',t1,et1)
    print *
    
    end program NuMatrix
    
    subroutine output_NuMatrix(dev)
    use MatrixParameters
    use OperParameters
    use ProjParameters
    use InputNuShell
    use OutputNuShell
    use ShellsX
    use Mscheme
    use Files
    use MatrixSupport
    
    implicit none
    integer,intent(in)::dev
    integer:: status
    character(len=14)::buffer
    call get_cmd_arg(1,buffer,status) 
    if (status /= -1) return
    
    write(dev,*) ' NuMatrix Program'
    write(dev,*)
    write(dev,*) ' Compile time Parameters:'
    write(dev,*)
    write(dev,*) ' Machine Parameters:'
    write(dev,*) ' No of bits per word',nbits,' Memory Access Width Bytes', nbytes
    write(dev,*)
    write(dev,*) ' Shell Model Parameters:'
    write(dev,*) ' Maximum no of shells',max_no_shells
    write(dev,*) ' Maximum no of words with',nbits,' bits is',nwords, '    for each m-scheme state'
    write(dev,*) ' Maximum value of 2*j for a shell',max_sj2
    write(dev,*) ' Maximum no of shell and isospin shell partitions',max_no_partitions
    write(dev,*) ' Maximum no of arrangements of nucleons in each shell (max_combine)',max_combine
    write(dev,*) ' Maximum no of m-vectors per partition all m',max_allmvector
    write(dev,*) ' Maximum no of m-vectors per partition single m',max_proj_mvec
    write(dev,*) ' Maximum no J,T values',max_no_cJstates,max_no_cTstates
    write(dev,*) ' Maximum Jscheme dimenson',max_total_JT
    write(dev,*) ' Maximum t.b. operator dimension',max_op 
    write(dev,*) ' Real kind for projection', rk
    write(dev,*) ' Real kind for matrix', rm
    write(dev,*) ' Real kind for t.b. operator', ri
    write(dev,*) ' Symmetry control', np_cntl
    write(dev,*) ' Hpart', hpart
    write(dev,*)

end subroutine output_NuMatrix