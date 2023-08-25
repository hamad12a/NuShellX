! This file conation al the modules specific to geneergic program pOpm for
! npShell.
module Opm

use ProjParameters
use MbmeParameters
use OpmParameters
use ShellsX
use Shells
use Mscheme
use Mvector
use InputNuShell
use OperParameters
use Finds
use Partition
use OpPartition
use TrdChar3
use XParameters

implicit none

integer:: nopart,notmvr,noOpm,noprojmvec,no_disk,indexo=0
integer:: j,t,j1,it1,nn,nmp1,nmp2,twrite1,twrite,dt,index=0
integer:: mindj,maxdj,mindt,maxdt,mindmj,maxdmj,mindmt,maxdmt,delmj,delmt
type(m_vector),dimension(:),allocatable:: mvec_lhs
type(m_vector),dimension(:),allocatable:: mvec_rhs
type(spartition),dimension(:),allocatable::spartr
integer,dimension(:),allocatable:: part_offsetr,dim_partr
integer:: dim_lhs,dim_rhs,no_spart,no_spartr,use_wordso
integer,dimension(:),allocatable:: part_offset
integer,dimension(:),allocatable:: dim_part
type(spartition),dimension(:),allocatable:: spart
character(len=1):: xtra
logical,save:: shop,proton=.false.,neutron=.false.,first=.true.
character(len=7)::opfile
character(len=6)::z1=repeat(char(0),6)
character(len=12),dimension(:),allocatable:: Opmidx 
type(trx),dimension(:),allocatable:: Opmidxx
character(len=6),dimension(:),allocatable:: Opm3idx
contains

    subroutine  Apm_op(baseb,baseo,k,k1,k2)
        implicit none
        integer,intent(in):: baseb,baseo,k,k1,k2
        type(m_vector):: mvec_copy2,mvecam2
        type(m_vector):: mvec_copy3,mvecap
        integer:: dnr,nb2,iw2,findm,np,npr
        integer:: nb3,iw3,sh2,sh3,icode
        integer:: m2,m3,tz2,tz3,iwm2,iwm3 
        real(kind=4):: phase2,phase3
        integer,dimension(2)::aop2,aop3
        type(spartition):: part_lhs,part_rhs,part_op,op_op2,op_op3,fin_op
        character(len=10):: ctime,ptime
        character(len=8):: date,pdate
        integer:: i
        integer(kind=ik):: iand,not,ibclr
        use_wordso=gwords()
        do np=1,no_spart
        if (output_control > 0) then
            call date_and_time(date,ctime)
            write(10,*) date,' ',ctime,' npart left',np
        end if
        call oread_lhs_basis(baseb,k1)
        if (dim_lhs==0) go to 50
        part_lhs=spart(np)
        do npr=1,no_spartr
        if (output_control > 3) then
            call date_and_time(pdate,ptime)
            write(10,*) pdate,' ',ptime,' npart right',npr
        end if
        if (dim_partr(npr) ==0) go to 45
        part_rhs=spartr(npr)
        part_op=part_lhs-part_rhs
        if (oanalyse_part_opap(part_op)) go to 45
        iwm3=use_wordso
        iwm2=use_wordso
        do dnr=1+part_offsetr(npr),part_offsetr(npr)+dim_partr(npr)
                    mvecam2=mvec_rhs(dnr)
                            m2=0
                            tz2=0
                            phase2=1.0
                            mvec_copy3=ph(mvecam2,word_mask)
                            m3=delmj
                            if (abs(m3)>max_sj2u) go to 40
                            tz3=delmt
                            if (abs(tz3)>1) go to 40
                            mvec_copy3=mand(mvec_jz_tz(m3,tz3),mvec_copy3)
                            if (mvec_copy3 == int(0,ik)) go to 40
                            do iw3=1,iwm3
                                aop3(1)=iw3
                                do
                                    if (mvec_copy3%fn(iw3) == 0) exit
                                    nb3=next_bit(mvec_copy3%fn(iw3))
                                    mvec_copy3%fn(iw3)=ibclr(mvec_copy3%fn(iw3),nb3-1)
                                    op_op3=0
                                    op_op3%shell(shn_lu(nb3,iw3))=op_op3%shell(shn_lu(nb3,iw3))+1
                                    aop3(2)=nb3
                                    sh3=shn_lu(nb3,iw3)
                                    if (shellxt.and.shtype(sh3)/=typex) go to 15
                                    mvecap=aop3+mvecam2
                                    phase3=qphase(aop3,mvecam2)
3                                   fin_op=part_op-op_op3
                                    if (fin_op /= 0) then
                                        go to 15
                                    end if
                                    findm=bfind(mvec_lhs,dim_lhs,mvecap)
                                    if (findm /=0 ) then
                                        index=index+1
                                        if (first) then
                                        first=.false.
                                        if (.not.shellxt) then
                                            if (sh3<=no_shellsp()) then
                                                proton=.true.
                                                print *, ' Proton APM'
                                            else
                                                neutron=.true.
                                                print *, ' Neutron APM'
                                            end if
                                        else
                                        if (shtype(sh3)=='a') then
                                            neutron=.true.
                                            print *, ' Type a APM'
                                        else
                                            proton=.true.
                                            print *, ' Type b APM'
                                        end if
                                        end if
                                        
                                        end if
                                        if (index > noOpm) then
                                            print *, ' Increase max_Opm in OpmParams or'
                                            print *, ' SET OpmDIM=I where I >>',index+1
                                            stop
                                        end if
                                        phase2=phase2*phase3
                                        call store(sh3,0,j1,j,np,npr,m3,0,findm,dnr&
                                                        ,phase2,delmj,delmt)
                                        if (shellxt) call store3(0,it1,t,nn,tz3,0)                              
                                    else
                                        print *, 'Mvec not found'
                                        print *, 'Final mvec'
                                        call swrite(6,mvecap)
                                        print *, 'Initial mvec'
                                        call swrite(6,mvec_rhs(dnr))
                                        stop
                                    end if
15                              end do
                            end do                                        

40       end do
45       end do         !npr
50       end do            !np
         return
    end subroutine Apm_op
    
    subroutine Amm_op(baseb,baseo,k,k1,k2)
        implicit none
        integer,intent(in):: baseb,baseo,k,k1,k2
        type(m_vector):: mvec_copy2,mvecam2
        type(m_vector):: mvec_copy3,mvecap
        integer:: dnr,nb2,iw2,findm,np,npr
        integer:: nb3,iw3,sh2,sh3,icode
        integer:: m2,m3,tz2,tz3,iwm2,iwm3 
        real(kind=4):: phase2,phase3
        integer,dimension(2)::aop2,aop3
        type(spartition):: part_lhs,part_rhs,part_op,op_op2,op_op3,fin_op
        character(len=10):: ctime,ptime
        character(len=8):: date,pdate
        integer:: i
        integer(kind=ik):: iand,not,ibclr
        use_wordso=gwords()
        do np=1,no_spart
        if (output_control > 0) then
            call date_and_time(date,ctime)
            write(10,*) date,' ',ctime,' npart left',np
        end if
        call oread_lhs_basis(baseb,k1)
        if (dim_lhs==0) go to 50
        part_lhs=spart(np)
        do npr=1,no_spartr
        if (output_control > 3) then
            call date_and_time(pdate,ptime)
            write(10,*) pdate,' ',ptime,' npart right',npr
        end if
        if (dim_partr(npr) ==0) go to 45
        part_rhs=spartr(npr)
        part_op=part_lhs-part_rhs
        if (oanalyse_part_opam(part_op)) go to 45
        iwm3=use_wordso
        iwm2=use_wordso
        do dnr=1+part_offsetr(npr),part_offsetr(npr)+dim_partr(npr)
                    mvec_copy2=mvec_rhs(dnr)
                    do iw2=1,iwm2
                        aop2(1)=iw2
                        do
                            if (mvec_copy2%fn(iw2) == 0) exit
                            nb2=next_bit(mvec_copy2%fn(iw2))
                            mvec_copy2%fn(iw2)=ibclr(mvec_copy2%fn(iw2),nb2-1)
                            m2=sm2_lu(nb2,iw2)
                            tz2=stz2_lu(nb2,iw2)
                            sh2=shn_lu(nb2,iw2)
                            if (shellxt.and.typex/=shtype(sh2)) go to 30
                            if (first) then
                                first=.false.
                                if(.not.shellxt)then
                                if (sh2<=no_shellsp()) then
                                    proton=.true.
                                    print *, ' Proton AMM'
                                else
                                    neutron=.true.
                                    print *, ' Neutron AMM'
                                end if
                                        else
                                        if (shtype(sh2)=='a') then
                                            neutron=.true.
                                            print *, ' Type a AMM'
                                        else
                                            proton=.true.
                                            print *, ' Type b AMM'
                                        end if
                                        end if
                            end if
                            op_op2=0
                            op_op2%shell(shn_lu(nb2,iw2))=op_op2%shell(shn_lu(nb2,iw2))-1
                            aop2(2)=nb2
                            mvecam2=aop2-mvec_rhs(dnr)
                            phase2=qphase(aop2,mvec_rhs(dnr))
                            m3=m2+delmj
                            if (abs(m3)/=0) go to 30
                            tz3=tz2+delmt
                            if (abs(tz3)/=0) go to 30
                                    phase3=1.0
3                                   fin_op=part_op-op_op2
                                    if (fin_op /= 0) then
                                        go to 30
                                    end if
                                    findm=bfind(mvec_lhs,dim_lhs,mvecam2)
                                    if (findm /=0 ) then
                                        index=index+1
                                        if (index > noOpm) then
                                            print *, ' Increase max_Opm in OpmParams or'
                                            print *, ' SET OPMDIM=I where I >>',index+1
                                            stop
                                        end if
                                        phase2=phase2*phase3
                                        if (btest((sj2(sh2)-m2)/2,0)) phase2=-phase2
                                        call store(0,sh2,j1,j,np,npr,0,m2,findm,dnr&
                                                        ,phase2,delmj,delmt)
                                        if (shellxt) call store3(0,it1,t,nn,0,tz2)                              
                                    else
                                        print *, 'Mvec not found'
                                        print *, 'Final mvec'
                                        call swrite(6,mvecap)
                                        print *, 'Initial mvec'
                                        call swrite(6,mvec_rhs(dnr))
                                        stop
                                    end if
30                      end do
                    end do

40       end do
45       end do         !npr
50       end do            !np
         return
    end subroutine Amm_op
    
    subroutine Idm_op(baseb,baseo,k,k1,k2)
        implicit none
        integer,intent(in):: baseb,baseo,k,k1,k2
        type(m_vector):: mvec_copy2,mvecam2
        type(m_vector):: mvec_copy3,mvecap
        integer:: dnr,nb2,iw2,findm,np,npr
        integer:: nb3,iw3,sh2,sh3,icode
        integer:: m2,m3,tz2,tz3,iwm2,iwm3 
        real(kind=4):: phase2,phase3
        integer,dimension(2)::aop2,aop3
        type(spartition):: part_lhs,part_rhs,part_op,op_op2,op_op3,fin_op
        character(len=10):: ctime,ptime
        character(len=8):: date,pdate
        integer:: i
        integer(kind=ik):: iand,not,ibclr
        use_wordso=gwords()
        do np=1,no_spart
        if (output_control > 0) then
            call date_and_time(date,ctime)
            write(10,*) date,' ',ctime,' npart left',np
        end if
        call oread_lhs_basis(baseb,k1)
        if (dim_lhs==0) go to 50
        part_lhs=spart(np)
        do npr=np,np
        if (output_control > 3) then
            call date_and_time(pdate,ptime)
            write(10,*) pdate,' ',ptime,' npart right',npr
        end if
        if (dim_partr(npr) ==0) go to 45
        part_rhs=spartr(npr)
        part_op=part_lhs-part_rhs
        if  (part_op/=0)  go to 45
        iwm3=use_wordso
        iwm2=use_wordso
        do dnr=1+part_offsetr(npr),part_offsetr(npr)+dim_partr(npr)
                    mvec_copy2=mvec_rhs(dnr)
                    do iw2=1,iwm2
                        do
                            if (mvec_copy2%fn(iw2) == 0) exit
                            nb2=next_bit(mvec_copy2%fn(iw2))
                            sh2=shn_lu(nb2,iw2)
                            if (first) then
                                first=.false.
                                if (.not.shellxt) then
                                if (sh2<=no_shellsp()) then
                                    proton=.true.
                                    print *, ' Proton IDM'
                                else
                                    neutron=.true.
                                    print *, ' Neutron IDM'
                                end if
                                        else
                                        if (shtype(sh2)=='a') then
                                            neutron=.true.
                                            print *, ' Type a IDM'
                                        else
                                            proton=.true.
                                            print *, ' Type b IDM'
                                        end if
                                        end if
                            end if
                            exit
                         end do
                      end do
                      sh2=0
                      sh3=0
                      m2=0
                      m3=0 
                                    findm=bfind(mvec_lhs,dim_lhs,mvec_copy2)
                                    if (findm /=0 ) then
                                        index=index+1
                                        if (index > noOpm) then
                                            print *, ' Increase max_Opm in OpmParams or'
                                            print *, ' SET  OPMDIM=I where I >>',index+1
                                            stop
                                        end if
                                        phase2=1.0
                                        call store(sh3,sh2,j1,j,np,npr,m3,m2,findm,dnr&
                                                        ,phase2,delmj,delmt)
                                        if (shellxt) call store3(0,it1,t,nn,tz3,tz2)                              
                                    else
                                        print *, 'Mvec not found'
                                        print *, 'Final mvec'
                                        call swrite(6,mvec_copy2)
                                        print *, 'Initial mvec'
                                        call swrite(6,mvec_rhs(dnr))
                                        stop
                                    end if

40       end do
45       end do         !npr
50       end do            !np
         return
    end subroutine Idm_op
    
    subroutine Trm_op(baseb,baseo,k,k1,k2)
        implicit none
        integer,intent(in):: baseb,baseo,k,k1,k2
        type(m_vector):: mvec_copy2,mvecam2
        type(m_vector):: mvec_copy3,mvecap
        integer:: dnr,nb2,iw2,findm,np,npr
        integer:: nb3,iw3,sh2,sh3,icode
        integer:: m2,m3,tz2,tz3,iwm2,iwm3 
        real(kind=4):: phase2,phase3
        integer,dimension(2)::aop2,aop3
        type(spartition):: part_lhs,part_rhs,part_op,op_op2,op_op3,fin_op
        character(len=10):: ctime,ptime
        character(len=8):: date,pdate
        integer:: i
        integer(kind=ik):: iand,not,ibclr
        use_wordso=gwords()
        do np=1,no_spart
        if (output_control > 0) then
            call date_and_time(date,ctime)
            write(10,*) date,' ',ctime,' npart left',np
        end if
        call oread_lhs_basis(baseb,k1)
        if (dim_lhs==0) go to 50
        part_lhs=spart(np)
        do npr=1,no_spartr
        if (output_control > 3) then
            call date_and_time(pdate,ptime)
            write(10,*) pdate,' ',ptime,' npart right',npr
        end if
        if (dim_partr(npr) ==0) go to 45
        part_rhs=spartr(npr)
        part_op=part_lhs-part_rhs
        if (oanalyse_part_optd(part_op)) go to 45
        iwm3=use_wordso
        iwm2=use_wordso
        do dnr=1+part_offsetr(npr),part_offsetr(npr)+dim_partr(npr)
                    mvec_copy2=mvec_rhs(dnr)
                    do iw2=1,iwm2
                        aop2(1)=iw2
                        do
                            if (mvec_copy2%fn(iw2) == 0) exit
                            nb2=next_bit(mvec_copy2%fn(iw2))
                            mvec_copy2%fn(iw2)=ibclr(mvec_copy2%fn(iw2),nb2-1)
                            m2=sm2_lu(nb2,iw2)
                            tz2=stz2_lu(nb2,iw2)
                            sh2=shn_lu(nb2,iw2)
                            if (first) then
                                first=.false.
                                if (.not.shellxt) then
                                if (sh2<=no_shellsp()) then
                                    proton=.true.
                                    print *, ' Proton TRM'
                                else
                                    neutron=.true.
                                    print *, ' Neutron TRM'
                                end if
                                        else
                                        if (shtype(sh2)=='a') then
                                            neutron=.true.
                                            print *, ' Type a TRM'
                                        else
                                            proton=.true.
                                            print *, ' Type b TRM'
                                        end if
                                        end if
                            end if
                            op_op2=0
                            op_op2%shell(shn_lu(nb2,iw2))=op_op2%shell(shn_lu(nb2,iw2))-1
                            aop2(2)=nb2
                            mvecam2=aop2-mvec_rhs(dnr)
                            phase2=qphase(aop2,mvec_rhs(dnr))
                            mvec_copy3=ph(mvecam2,word_mask)
                            m3=m2+delmj
                            if (abs(m3)>max_sj2u) go to 30
                            tz3=tz2+delmt
                            if (abs(tz3)>1) go to 30
                            mvec_copy3=mand(mvec_jz_tz(m3,tz3),mvec_copy3)
                            if (mvec_copy3 == int(0,ik)) go to 30
                            do iw3=1,iwm3
                                aop3(1)=iw3
                                do
                                    if (mvec_copy3%fn(iw3) == 0) exit
                                    nb3=next_bit(mvec_copy3%fn(iw3))
                                    mvec_copy3%fn(iw3)=ibclr(mvec_copy3%fn(iw3),nb3-1)
                                    op_op3=op_op2
                                    op_op3%shell(shn_lu(nb3,iw3))=op_op3%shell(shn_lu(nb3,iw3))+1
                                    aop3(2)=nb3
                                    sh3=shn_lu(nb3,iw3)
                                    if (shellxt.and.shtype(sh3)/=shtype(sh2)) go to 15
                                    mvecap=aop3+mvecam2
                                    phase3=qphase(aop3,mvecam2)
3                                   fin_op=part_op-op_op3
                                    if (fin_op /= 0) then
                                        go to 15
                                    end if
                                    findm=bfind(mvec_lhs,dim_lhs,mvecap)
                                    if (findm /=0 ) then
                                        index=index+1
                                        if (index > noOpm) then
                                            print *, ' Increase max_Opm in OpmParams or'
                                            print *, ' SET OPMDIM=I where I >>',index+1
                                            stop
                                        end if
                                        phase2=phase2*phase3
                                        if (btest((sj2(sh2)-m2)/2,0)) phase2=-phase2
                                        call store(sh3,sh2,j1,j,np,npr,m3,m2,findm,dnr&
                                                        ,phase2,delmj,delmt)
                                        if (shellxt) call store3(0,it1,t,nn,tz3,tz2)                              
                                    else
                                        print *, 'Mvec not found'
                                        print *, 'Final mvec'
                                        call swrite(6,mvecap)
                                        print *, 'Initial mvec'
                                        call swrite(6,mvec_rhs(dnr))
                                        stop
                                    end if
15                              end do
                            end do                                        
30                      end do
                    end do

40       end do
45       end do         !npr
50       end do            !np
         return
    end subroutine Trm_op
       
    subroutine Tmm_op(baseb,baseo,k,k1,k2)
        implicit none
        integer,intent(in):: baseb,baseo,k,k1,k2
        type(m_vector):: mvec_copy2,mvecam2
        type(m_vector):: mvec_copy3,mvecap
        integer:: dnr,nb2,iw2,findm,np,npr
        integer:: nb3,iw3,sh2,sh3,icode
        integer:: m2,m3,tz2,tz3,iwm2,iwm3 
        real(kind=4):: phase2,phase3
        integer,dimension(2)::aop2,aop3
        type(spartition):: part_lhs,part_rhs,part_op,op_op2,op_op3,fin_op
        character(len=10):: ctime,ptime
        character(len=8):: date,pdate
        integer:: i
        integer(kind=ik):: iand,not,ibclr
        use_wordso=gwords()
        do np=1,no_spart
        if (output_control > 0) then
            call date_and_time(date,ctime)
            write(10,*) date,' ',ctime,' npart left',np
        end if
        call oread_lhs_basis(baseb,k1)
        if (dim_lhs==0) go to 50
        part_lhs=spart(np)
        do npr=1,no_spartr
        if (output_control > 3) then
            call date_and_time(pdate,ptime)
            write(10,*) pdate,' ',ptime,' npart right',npr
        end if
        if (dim_partr(npr) ==0) go to 45
        part_rhs=spartr(npr)
        part_op=part_lhs-part_rhs
        if (oanalyse_part_optm(part_op)) go to 45
        iwm3=use_wordso
        iwm2=use_wordso
        do dnr=1+part_offsetr(npr),part_offsetr(npr)+dim_partr(npr)
                    mvec_copy2=mvec_rhs(dnr)
                    do iw2=1,iwm2
                        aop2(1)=iw2
                        do
                            if (mvec_copy2%fn(iw2) == 0) exit
                            nb2=next_bit(mvec_copy2%fn(iw2))
                            mvec_copy2%fn(iw2)=ibclr(mvec_copy2%fn(iw2),nb2-1)
                            m2=sm2_lu(nb2,iw2)
                            tz2=stz2_lu(nb2,iw2)
                            sh2=shn_lu(nb2,iw2)
                            op_op2=0
                            op_op2%shell(shn_lu(nb2,iw2))=op_op2%shell(shn_lu(nb2,iw2))-1
                            aop2(2)=nb2
                            mvecam2=aop2-mvec_rhs(dnr)
                            phase2=qphase(aop2,mvec_rhs(dnr))
                            mvec_copy3=mvecam2
                            m3=-m2-delmj
                            if (abs(m3)>max_sj2u) go to 30
                            tz3=-tz2-delmt
                            if (abs(tz3)>1) go to 30
                            mvec_copy3=mand(mvec_jz_tz(m3,tz3),mvec_copy3)
                            if (mvec_copy3 == int(0,ik)) go to 30
                            do iw3=1,iwm3
                                aop3(1)=iw3
                                do
                                    if (mvec_copy3%fn(iw3) == 0) exit
                                    nb3=next_bit(mvec_copy3%fn(iw3))
                                    mvec_copy3%fn(iw3)=ibclr(mvec_copy3%fn(iw3),nb3-1)
                                    op_op3=op_op2
                                    op_op3%shell(shn_lu(nb3,iw3))=op_op3%shell(shn_lu(nb3,iw3))-1
                                    aop3(2)=nb3
                                    sh3=shn_lu(nb3,iw3)
                                    if (shellxt.and.shtype(sh3)/=shtype(sh2)) go to 15
                                    mvecap=aop3-mvecam2
                                    phase3=qphase(aop3,mvecam2)
3                                   fin_op=part_op-op_op3
                                    if (fin_op /= 0) then
                                        go to 15
                                    end if
                                    findm=bfind(mvec_lhs,dim_lhs,mvecap)
                                    if (first) then
                                    first=.false.
                                        if (.not.shellxt) then
                                        if (sh2<=no_shellsp()) then
                                        proton=.true.
                                        print *, ' Proton TMM'
                                    else
                                        neutron=.true.
                                        print *, ' Neutron TMM'
                                    end if
                                        else
                                        if (shtype(sh2)=='a') then
                                            neutron=.true.
                                            print *, ' Type a TPM'
                                        else
                                            proton=.true.
                                            print *, ' Type b TPM'
                                        end if
                                        end if
                                    end if
                                   if (findm /=0 ) then
                                        index=index+1
                                        if (index > noOpm) then
                                            print *, ' Increase max_Opm in OpmParams or'
                                            print *, ' SET OPMDIM=I where I >>',index+1
                                            stop
                                        end if
                                        phase2=phase2*phase3
                                        if (btest((sj2(sh2)+sj2(sh3)-m2-m3)/2,0)) phase2=-phase2
                                        call store(sh3,sh2,j1,j,np,npr,m3,m2,findm,dnr&
                                                        ,phase2,delmj,delmt)
                                        if (shellxt) call store3(0,it1,t,nn,tz3,tz2)                              
                                    else
                                        print *, 'Mvec not found'
                                        print *, 'Final mvec'
                                        call swrite(6,mvecap)
                                        print *, 'Initial mvec'
                                        call swrite(6,mvec_rhs(dnr))
                                        stop
                                    end if
15                              end do
                            end do                                        
30                      end do
                    end do

40       end do
45       end do         !npr
50       end do            !np
         return
    end subroutine Tmm_op
    
    
    subroutine Tpm_op(baseb,baseo,k,k1,k2)
        implicit none
        integer,intent(in):: baseb,baseo,k,k1,k2
        type(m_vector):: mvec_copy2,mvecam2
        type(m_vector):: mvec_copy3,mvecap
        integer:: dnr,nb2,iw2,findm,np,npr
        integer:: nb3,iw3,sh2,sh3,icode
        integer:: m2,m3,tz2,tz3,iwm2,iwm3 
        real(kind=4):: phase2,phase3
        integer,dimension(2)::aop2,aop3
        type(spartition):: part_lhs,part_rhs,part_op,op_op2,op_op3,fin_op
        character(len=10):: ctime,ptime
        character(len=8):: date,pdate
        integer:: i
        integer(kind=ik):: iand,not,ibclr
        use_wordso=gwords()
        do np=1,no_spart
        if (output_control > 0) then
            call date_and_time(date,ctime)
            write(10,*) date,' ',ctime,' npart left',np
        end if
        call oread_lhs_basis(baseb,k1)
        if (dim_lhs==0) go to 50
        part_lhs=spart(np)
        do npr=1,no_spartr
        if (output_control > 3) then
            call date_and_time(pdate,ptime)
            write(10,*) pdate,' ',ptime,' npart right',npr
        end if
        if (dim_partr(npr) ==0) go to 45
        part_rhs=spartr(npr)
        part_op=part_lhs-part_rhs
        if (oanalyse_part_optp(part_op)) go to 45
        iwm3=use_wordso
        iwm2=use_wordso
        do dnr=1+part_offsetr(npr),part_offsetr(npr)+dim_partr(npr)
                    mvec_copy2=ph(mvec_rhs(dnr),word_mask)
                    do iw2=1,iwm2
                        aop2(1)=iw2
                        do
                            if (mvec_copy2%fn(iw2) == 0) exit
                            nb2=next_bit(mvec_copy2%fn(iw2))
                            mvec_copy2%fn(iw2)=ibclr(mvec_copy2%fn(iw2),nb2-1)
                            m2=sm2_lu(nb2,iw2)
                            tz2=stz2_lu(nb2,iw2)
                            sh2=shn_lu(nb2,iw2)
                            op_op2=0
                            op_op2%shell(shn_lu(nb2,iw2))=op_op2%shell(shn_lu(nb2,iw2))+1
                            aop2(2)=nb2
                            mvecam2=aop2+mvec_rhs(dnr)
                            phase2=qphase(aop2,mvec_rhs(dnr))
                            mvec_copy3=ph(mvecam2,word_mask)
                            m3=-m2+delmj
                            if (abs(m3)>max_sj2u) go to 30
                            tz3=-tz2+delmt
                            if (abs(tz3)>1) go to 30
                            mvec_copy3=mand(mvec_jz_tz(m3,tz3),mvec_copy3)
                            if (mvec_copy3 == int(0,ik)) go to 30
                            do iw3=1,iwm3
                                aop3(1)=iw3
                                do
                                    if (mvec_copy3%fn(iw3) == 0) exit
                                    nb3=next_bit(mvec_copy3%fn(iw3))
                                    mvec_copy3%fn(iw3)=ibclr(mvec_copy3%fn(iw3),nb3-1)
                                    op_op3=op_op2
                                    op_op3%shell(shn_lu(nb3,iw3))=op_op3%shell(shn_lu(nb3,iw3))+1
                                    aop3(2)=nb3
                                    sh3=shn_lu(nb3,iw3)
                                    if (shellxt.and.shtype(sh3)/=shtype(sh2)) go to 15
                                    mvecap=aop3+mvecam2
                                    phase3=qphase(aop3,mvecam2)
3                                   fin_op=part_op-op_op3
                                    if (fin_op /= 0) then
                                        go to 15
                                    end if
                                    findm=bfind(mvec_lhs,dim_lhs,mvecap)
                                    if (first) then
                                        if (.not.shellxt) then
                                            if (sh2<=no_shellsp()) then
                                                proton=.true.
                                                print *, ' Proton TPM'
                                            else
                                            neutron=.true.
                                            print *, ' Neutron TPM'
                                            end if
                                        else
                                            if (shtype(sh2)=='a') then
                                                neutron=.true.
                                                print *, ' Type a TPM'
                                            else
                                                proton=.true.
                                                print *, ' Type b TPM'
                                            end if
                                        end if
                                    first=.false.
                                    end if
                                    if (findm /=0 ) then
                                        index=index+1
                                        if (index > noOpm) then
                                            print *, ' Increase max_Opm in OpmParams or'
                                            print *, ' SET OPMDIM=I where I >>',index+1
                                            stop
                                        end if
                                        phase2=phase2*phase3
                                        call store(sh3,sh2,j1,j,np,npr,m3,m2,findm,dnr&
                                                        ,phase2,delmj,delmt)
                                        if (shellxt) call store3(0,it1,t,nn,tz3,tz2)                              
                                    else
                                        print *, 'Mvec not found'
                                        print *, 'Final mvec'
                                        call swrite(6,mvecap)
                                        print *, 'Initial mvec'
                                        call swrite(6,mvec_rhs(dnr))
                                        stop
                                    end if
15                              end do
                            end do                                        
30                      end do
                    end do

40       end do
45       end do         !npr
50       end do            !np
         return
    end subroutine Tpm_op
    
    subroutine Tpm2_op(baseb,baseo,k,k1,k2)
        implicit none
        integer,intent(in):: baseb,baseo,k,k1,k2
        type(m_vector):: mvec_copy1,mvec_copy2,mvecam,mvecam2
        type(m_vector):: mvec_copy3,mvec_copy4,mvecap,mvecap2
        integer:: dnr,nb1,nb2,iw1,iw2,findm,np,npr
        integer:: nb3,nb4,iw3,iw4
        integer:: m1,m2,m3,m4,tz1,tz2,tz3,tz4,iwm1,iwm2,iwm3,iwm4 
        real(kind=4):: phase1,phase2,phase4,phase3,ph12,ph123
        integer:: sh1,sh2,sh3,sh4,use_wordst
        integer,dimension(2)::aop1,aop2,aop3,aop4
        type(spartition):: part_lhs,part_rhs,part_op,op_op1,op_op2,op_op3,op_op4,fin_op
        logical:: warnmvec=.false.,warnmxi=.false.
        character(len=10):: ctime,ptime
        character(len=8):: date,pdate
        integer:: tindex,i,shmin1,shmin3,shl
        integer::thrn,OMP_get_thread_num
        integer(kind=ik):: ibclr,iadd,ieor,ior,not

        use_wordst=gwords()
        use_wordso=use_wordst
        tindex=0
        do np=1,no_spart
        if (output_control > 0) then
            call date_and_time(date,ctime)
            write(10,*) date,' ',ctime,' npart left',np
        end if
        call oread_lhs_basis(baseb,k1)
        if (dim_lhs==0) go to 50
        part_lhs=spart(np)
        do npr=1,no_spartr
        if (output_control > 3) then
            call date_and_time(pdate,ptime)
            write(10,*) pdate,' ',ptime,' npart right',npr
        end if
        if (dim_partr(npr) ==0) go to 45
        part_rhs=spartr(npr)
        part_op=part_lhs-part_rhs
        if (tptm2_analyse_part_op(part_op)) go to 45
!$OMP PARALLEL DEFAULT(PRIVATE) &
!$OMP SHARED(part_offsetr,dim_partr,mvec_rhs,mvec_lhs,mvec_jz_tz,word_mask,index,&
!$OMP shn_lu,sm2_lu,sj2_lu,stz2_lu,use_wordst,npr,np,part_op,noOpm,Opmidx,Opmidxx,Opm3idx,&
!$OMP dim_lhs,max_sj2u,delmj,delmt,baseo,k2,first,shtype,neutron,proton,j1,j,it1,t,nn,&
!$OMP shellxt) 
        iwm1=use_wordst
        iwm3=use_wordst
        iwm2=use_wordst
        iwm4=use_wordst
!$OMP DO ORDERED SCHEDULE(DYNAMIC)
        do dnr=1+part_offsetr(npr),part_offsetr(npr)+dim_partr(npr)
            thrn=OMP_get_thread_num()+1
            mvec_copy1=mvec_rhs(dnr)
            if (mvec_copy1==int(0,ik)) go to 40
!loop over rhmost annihilation operator words then bits   a+a-a-         
            do iw1=1,iwm1
                aop1(1)=iw1
                do
                    if (mvec_copy1%fn(iw1) == 0) exit
!identify next bit                    
                    nb1=next_bit(mvec_copy1%fn(iw1))
                    mvec_copy1%fn(iw1)=ibclr(mvec_copy1%fn(iw1),nb1-1)
!lots of stuff to help speed up and keep track
                    op_op1=0
                    sh1=shn_lu(nb1,iw1)
                    op_op1%shell(sh1)=-1
                    m1=sm2_lu(nb1,iw1)
                    tz1=stz2_lu(nb1,iw1)
                    aop1(2)=nb1
!apply annihilation operator                    
                    mvecam=aop1-mvec_rhs(dnr)
                    mvec_copy2=mvecam
                    phase1=qphase(aop1,mvec_rhs(dnr))
!speedup trick                    
!loop over lh annihilation operator  words then bits a+a+(a)a                  
                    do iw2=1,iwm2
                        aop2(1)=iw2
                        do
                            if (mvec_copy2%fn(iw2) == 0) exit
!identify next bit                            
                            nb2=next_bit(mvec_copy2%fn(iw2))
                            mvec_copy2%fn(iw2)=ibclr(mvec_copy2%fn(iw2),nb2-1)
!gain lots of work to keep track and speed up                            
                            sh2=shn_lu(nb2,iw2)
                            if (shellxt.and.shtype(sh1)/=shtype(sh2)) go to 30
                            op_op2=op_op1
                            op_op2%shell(sh2)=op_op2%shell(sh2)-1
                            m2=sm2_lu(nb2,iw2)
                            tz2=stz2_lu(nb2,iw2)
                            aop2(2)=nb2
                            m3=m1+m2+delmj
                            tz3=tz1+tz2+delmt
                            if (abs(tz3) > 1) go to 30
                            if (abs(m3) > max_sj2u) go to 30
!apply annihilation op                            
                            mvecam2=aop2-mvecam
                            phase2=qphase(aop2,mvecam)
                            mvec_copy3=ph(mvecam2,word_mask)
                            mvec_copy3=mand(mvec_jz_tz(m3,tz3),mvec_copy3)
2                           ph12=phase1*phase2
!loop over  rh creation operator words then bits a+(a+)aa                          
3                           do iw3=1,iwm3
                                aop3(1)=iw3
                                do
                                    if (mvec_copy3%fn(iw3) == 0) exit
!identify next bit                                    
                                    nb3=next_bit(mvec_copy3%fn(iw3))
                                    mvec_copy3%fn(iw3)=ibclr(mvec_copy3%fn(iw3),nb3-1)
! then all the usual things                                    
                                    sh3=shn_lu(nb3,iw3)
                                    if (shellxt.and.shtype(sh3)/=shtype(sh1)) go to 20
                                    op_op3=op_op2
                                    op_op3%shell(sh3)=op_op3%shell(sh3)+1
                                    m3=sm2_lu(nb3,iw3)
                                    tz3=stz2_lu(nb3,iw3)
                                    aop3(2)=nb3
                                    mvecap=aop3+mvecam2
                                    phase3=qphase(aop3,mvecam2)
                                    ph123=ph12*phase3
!final partition check                                            
4                                           fin_op=part_op-op_op3
                                            if (fin_op /= 0) then
                                                go to 20
                                            end if
                                            findm=bfind(mvec_lhs,dim_lhs,mvecap)
!$OMP ORDERED
                                            if (findm /=0 ) then
                                                index=index+1
                                                if (index > noOpm) then
                                                    print *, ' Increase max_indexm in MbmeParams or'
                                                    print *, ' SET OPMDIM=I where I >>',index+1
                                                    stop
                                                end if
                                            if (first) then
                                            first=.false.
                                            if (shtype(sh2)=='a') then
                                                    neutron=.true.
                                                    print *, ' Type a PMM'
                                            else
                                                    proton=.true.
                                                    print *, ' Type b PMM'
                                            end if
                                            end if
                                                call store(sh2,sh1,j1,j,np,npr,m2,m1,findm,dnr&
                                                        ,ph123,delmj,delmt)
                                                if (shellxt) call store3(sh3,it1,t,nn,tz2,tz1)                              
                                            else
                                                print *, 'Mvec not found'
                                                print *, 'Final mvec'
                                                call swrite(6,mvecap)
                                                print *, 'Initial mvec'
                                                call swrite(6,mvec_rhs(dnr))
                                                stop
                                            end if
!$OMP END ORDERED
20                               end do
                            end do                                        
30                      end do
                    end do
35              end do
            end do
! write out me's 
40          end do
!$OMP END PARALLEL
45       end do         !npr
50       end do            !np

    end subroutine Tpm2_op

    subroutine Tp2m_op(baseb,baseo,k,k1,k2)
        implicit none
        integer,intent(in):: baseb,baseo,k,k1,k2
        type(m_vector):: mvec_copy1,mvec_copy2,mvecam,mvecam2
        type(m_vector):: mvec_copy3,mvec_copy4,mvecap,mvecap2
        integer:: dnr,nb1,nb2,iw1,iw2,findm,np,npr
        integer:: nb3,nb4,iw3,iw4
        integer:: m1,m2,m3,m4,tz1,tz2,tz3,tz4,iwm1,iwm2,iwm3,iwm4 
        real(kind=4):: phase1,phase2,phase4,phase3,ph12,ph123
        integer:: sh1,sh2,sh3,sh4,use_wordst
        integer,dimension(2)::aop1,aop2,aop3,aop4
        type(spartition):: part_lhs,part_rhs,part_op,op_op1,op_op2,op_op3,op_op4,fin_op
        logical:: warnmvec=.false.,warnmxi=.false.
        character(len=10):: ctime,ptime
        character(len=8):: date,pdate
        integer:: tindex,i,shmin1,shmin3,shl
        integer::thrn,OMP_get_thread_num
        integer(kind=ik):: ibclr,iadd,ieor,ior,not
        use_wordst=gwords()
        use_wordso=use_wordst
        tindex=0
        do np=1,no_spart
        if (output_control > 0) then
            call date_and_time(date,ctime)
            write(10,*) date,' ',ctime,' npart left',np
        end if
        call oread_lhs_basis(baseb,k1)
        if (dim_lhs==0) go to 50
        part_lhs=spart(np)
        do npr=1,no_spartr
        if (output_control > 3) then
            call date_and_time(pdate,ptime)
            write(10,*) pdate,' ',ptime,' npart right',npr
        end if
        if (dim_partr(npr) ==0) go to 45
        part_rhs=spartr(npr)
        part_op=part_lhs-part_rhs
        if (tp2tm_analyse_part_op(part_op)) go to 45
!$OMP PARALLEL DEFAULT(PRIVATE) &
!$OMP SHARED(part_offsetr,dim_partr,mvec_rhs,mvec_lhs,mvec_jz_tz,word_mask,index,&
!$OMP shn_lu,sm2_lu,stz2_lu,use_wordst,npr,np,part_op,noOpm,Opmidx,Opmidxx,Opm3idx,&
!$OMP dim_lhs,max_sj2u,delmj,delmt,baseo,k2,first,shtype,neutron,proton,j1,j,it1,t,nn,&
!$OMP shellxt) 
        iwm1=use_wordst
        iwm3=use_wordst
        iwm2=use_wordst
        iwm4=use_wordst
!$OMP DO ORDERED SCHEDULE(DYNAMIC)
        do dnr=1+part_offsetr(npr),part_offsetr(npr)+dim_partr(npr)
            thrn=OMP_get_thread_num()+1
            mvec_copy2=mvec_rhs(dnr)
            if (mvec_copy2==int(0,ik)) go to 40
                    m1=0
                    tz1=0
                    phase1=1.0
!apply annihilation operator                    
                    mvecam=mvec_rhs(dnr)
!loop over lh annihilation operator  words then bits a+a+(a)a                  
                    do iw2=1,iwm2
                        aop2(1)=iw2
                        do
                            if (mvec_copy2%fn(iw2) == 0) exit
!identify next bit                            
                            nb2=next_bit(mvec_copy2%fn(iw2))
                            mvec_copy2%fn(iw2)=ibclr(mvec_copy2%fn(iw2),nb2-1)
                            sh2=shn_lu(nb2,iw2)
                            op_op2=0
                            op_op2%shell(sh2)=op_op2%shell(sh2)-1
                            m2=sm2_lu(nb2,iw2)
                            tz2=stz2_lu(nb2,iw2)
                            aop2(2)=nb2
!apply annihilation op                            
                            mvecam2=aop2-mvecam
                            mvec_copy3=ph(mvecam2,word_mask)
                            phase2=qphase(aop2,mvecam)
2                           ph12=phase1*phase2
!loop over  rh creation operator words then bits a+(a+)aa                          
                            do iw3=1,iwm3
                                aop3(1)=iw3
                                do
                                    if (mvec_copy3%fn(iw3) == 0) exit
!identify next bit                                    
                                    nb3=next_bit(mvec_copy3%fn(iw3))
                                    mvec_copy3%fn(iw3)=ibclr(mvec_copy3%fn(iw3),nb3-1)
! then all the usual things                                    
                                    sh3=shn_lu(nb3,iw3)
                                    if (shellxt.and.shtype(sh3)/=shtype(sh2)) go to 20
                                    op_op3=op_op2
                                    op_op3%shell(sh3)=op_op3%shell(sh3)+1
                                    m3=sm2_lu(nb3,iw3)
                                    tz3=stz2_lu(nb3,iw3)
                                    m4=m2-m3+delmj
                                    tz4=tz2-tz3+delmt
                                    if (abs(tz4) /= 1) go to 20
                                    if (abs(m4) > max_sj2u) go to 20
                                    aop3(2)=nb3
                                    mvecap=aop3+mvecam2
                                    phase3=qphase(aop3,mvecam2)
3                                   ph123=ph12*phase3
                                    mvec_copy4=ph(mvecap,word_mask)
!application of another trick and with mvector with bits only with m4 and tz4                                       
                                    mvec_copy4=mand(mvec_jz_tz(m4,tz4),mvec_copy4)
!final loop over lh creation operator words then bits (a+)a+aa                                    
                                    do iw4=1,iwm4
                                        aop4(1)=iw4
                                        do
                                            if (mvec_copy4%fn(iw4) == 0) exit
                                            nb4=next_bit(mvec_copy4%fn(iw4))
                                            mvec_copy4%fn(iw4)=ibclr(mvec_copy4%fn(iw4),nb4-1)
!identify last bit                                            
                                            sh4=shn_lu(nb4,iw4)
                                            if (shellxt.and.shtype(sh4)/=shtype(sh2)) go to 15
                                            aop4(2)=nb4
                                            op_op4=op_op3
                                            op_op4%shell(sh4)=op_op4%shell(sh4)+1
                                            mvecap2=aop4+mvecap
                                            if (mvecap2==int(0,ik))go to 15
                                            phase4=qphase(aop4,mvecap)
!final partition check                                            
4                                           fin_op=part_op-op_op4
                                            if (fin_op /= 0) then
                                                go to 15
                                            end if
                                            findm=bfind(mvec_lhs,dim_lhs,mvecap2)
!$OMP ORDERED
                                            if (findm /=0 ) then
                                                index=index+1
                                                if (index > noOpm) then
                                                    print *, ' Increase max_indexm in MbmeParams or'
                                                    print *, ' SET OPMDIM=I where I >>',index+1
                                                    stop
                                                end if
                                                if (first) then
                                                first=.false.
                                                if (shtype(sh2)=='a') then
                                                    neutron=.true.
                                                    print *, ' Type a PPM'
                                                else
                                                    proton=.true.
                                                    print *, ' Type b PPM'
                                                end if
                                                end if
                                                call store(sh3,sh2,j1,j,np,npr,m3,m2,findm,dnr&
                                                        ,ph123*phase4,delmj,delmt)
                                                if (shellxt) call store3(sh4,it1,t,nn,tz3,tz2)                              
                               
                                            else
                                                print *, 'Mvec not found'
                                                print *, 'Final mvec'
                                                call swrite(6,mvecap2)
                                                print *, 'Initial mvec'
                                                call swrite(6,mvec_rhs(dnr))
                                                stop
                                            end if
!$OMP END ORDERED
15                                       end do
                                     end do
20                               end do
                            end do                                        
30                      end do
                    end do
40       end do
!$OMP END PARALLEL
45       end do         !npr
50       end do            !np

    end subroutine Tp2m_op


    subroutine store(sh3,sh2,j1,j,np,npr,m3,m2,findm,dnr,phase2,delmj,delmt)
    implicit none
    integer, intent(in)::sh3,sh2,j1,j,np,npr,m3,m2,findm,dnr,delmj,delmt
    integer:: i
    real(kind=4),intent(in):: phase2
                                        Opmidx(index)(isl:isl)=char(sh3)
                                        Opmidx(index)(isr:isr)=char(sh2)
                                        Opmidx(index)(ilm:ilm)=char(128+delmj)
                                        Opmidx(index)(ijl:ijl)=char(abs(j1))
                                        Opmidx(index)(ijr:ijr)=char(abs(j))
                                        Opmidx(index)(inl:inl+1)=char2(np)
                                        Opmidx(index)(inr:inr+1)=char2(npr)
                                        Opmidx(index)(iml:iml)=char(128+m3)
                                        Opmidx(index)(imr:imr)=char(128+m2)
                                        Opmidx(index)(ilt:ilt)=char(128+delmt)
                                        Opmidxx(index)%qux(1)=findm
                                        Opmidxx(index)%qux(2)=dnr-part_offsetr(npr)
                                        Opmidxx(index)%qux(3)=int(phase2)
                                        Opmidxx(index)%qux(4)=j1
                                        Opmidxx(index)%qux(5)=j
                                        Opm3idx(index)(1:4)=char(0)//char(0)//char(0)//char(0)
    if (index==max_Opm) then
    do i=1,index
    write(2) Opmidx(i)(1:12),Opmidxx(i)  !,Opm3idx(i)(1:6)
    end do
    no_disk=no_disk+1
    index=0   
    end if                           
    end subroutine store
    
    subroutine store3(shx,itl,itr,nn,tzl,tzr)
    implicit none
    
    integer, intent(in):: shx,itl,itr,nn,tzl,tzr
    
                                         Opm3idx(index)(1:1)=char(shx)
                                         Opm3idx(index)(2:2)=char(abs(itl))
                                         Opm3idx(index)(3:3)=char(abs(itr))
                                         Opm3idx(index)(4:4)=char(nn)
                                         Opm3idx(index)(5:5)=char(128+tzl)
                                         Opm3idx(index)(6:6)=char(128+tzr)
    end subroutine store3
    
    subroutine oread_lhs_basis(base,kk)
        implicit none
        integer,intent(in):: base,kk
        integer:: nbuf,pindx,nJT,jj,tt,Jmax,Tmax,i
                
            read(base+kk,err=10,end=11) pindx,nJT,jj,tt,nbuf,Jmax,Tmax      
            read(base+kk,err=10,end=11) spart(pindx)
            if (nbuf > noprojmvec) then
                print * ,' Increase max_proj_mvec in Parameters or'
                print * ,' SET MVECTORS=I where I >>',nbuf
                stop
            end if
            if (nbuf > 0 ) read(base+kk,err=10,end=11) (mvec_lhs(i)%fn(1:use_wordso),i=1,nbuf)
            dim_lhs=nbuf
            dim_part(pindx)=nbuf
            if (nJT==0) then
            dim_part(pindx)=0
            dim_lhs=0
            end if
            return
10          print *, ' Error reading lhs mvectors '
            stop
11          print *, ' EOF reading lhs mvectors '
            stop
    end subroutine oread_lhs_basis

    pure function oanalyse_part_opap(part)
        implicit none
        logical:: oanalyse_part_opap
        type(spartition),intent(in)::part
        integer:: ones,twos,mones,mtwos,i,ps

        oanalyse_part_opap=.false.
        
        if (part > 1 .or. part < 0) then
            oanalyse_part_opap=.true.
            return
        end if
         
        ones=0
        twos=0
        mones=0
        mtwos=0
        
        do i=1,no_shells
            ps=part%shell(i)
            if (ps == 1) then
                ones=ones+1
            end if
            if (ps == -1) then
                mones=mones+1
            end if    
            if (ps == 2) then
                twos=twos+1
            end if
            if (ps == -2) then
                mtwos=mtwos+1
            end if
        end do
        if (ones > 1 .or. mones  > 0 .or. twos > 0 .or. mtwos  >  0) then
            oanalyse_part_opap=.true.
            return
        end if
                
        if  (ones == 1 .and. mones == 0 .and. twos == 0 .and. mtwos == 0) then
            return
        end if
        
        oanalyse_part_opap=.true.
        
    end function oanalyse_part_opap

    pure function oanalyse_part_opam(part)
        implicit none
        logical:: oanalyse_part_opam
        type(spartition),intent(in)::part
        integer:: ones,twos,mones,mtwos,i,ps

        oanalyse_part_opam=.false.
        
        if (part > 0 .or. part < -1) then
            oanalyse_part_opam=.true.
            return
        end if
         
        ones=0
        twos=0
        mones=0
        mtwos=0
        
        do i=1,no_shells
            ps=part%shell(i)
            if (ps == 1) then
                ones=ones+1
            end if

            if (ps == -1) then
                mones=mones+1
            end if    
            if (ps == 2) then
                twos=twos+1
            end if
            if (ps == -2) then
                mtwos=mtwos+1
            end if
        end do

        if (ones > 0 .or. mones  > 1 .or. twos > 0 .or. mtwos  >  0) then
            oanalyse_part_opam=.true.
            return
        end if
                
        if  (ones == 0 .and. mones == 1 .and. twos == 0 .and. mtwos == 0) then
            return
        end if
        
        oanalyse_part_opam=.true.
        
    end function oanalyse_part_opam 

    pure function oanalyse_part_optd(part)
        implicit none
        logical:: oanalyse_part_optd
        type(spartition),intent(in)::part
        integer:: ones,twos,mones,mtwos,i,ps

        oanalyse_part_optd=.false.
        
        if (part > 1 .or. part < -1) then
            oanalyse_part_optd=.true.
            return
        end if
         
        ones=0
        twos=0
        mones=0
        mtwos=0
        
        do i=1,no_shells
            ps=part%shell(i)
            if (ps == 1) then
                ones=ones+1
            end if
            if (ps == -1) then
                mones=mones+1
            end if    
            if (ps == 2) then
                twos=twos+1
            end if
            if (ps == -2) then
                mtwos=mtwos+1
            end if
        end do
        if (ones > 2 .or. mones  > 2 .or. twos > 1 .or. mtwos  >  1) then
            oanalyse_part_optd=.true.
            return
        end if
                
        if  (ones == 1 .and. mones == 1 .and. twos == 0 .and. mtwos == 0) then
            return
        end if
        if  (ones == 0 .and. mones == 0 .and. twos == 0 .and. mtwos == 0) then
            return
        end if
        
        oanalyse_part_optd=.true.
        
    end function oanalyse_part_optd 
    
    pure function oanalyse_part_optm(part)
        implicit none
        logical:: oanalyse_part_optm
        type(spartition),intent(in)::part
        integer:: ones,twos,mones,mtwos,i,ps

        oanalyse_part_optm=.false.
        
        if (part > 0 .or. part < -2) then
            oanalyse_part_optm=.true.
            return
        end if
         
        ones=0
        twos=0
        mones=0
        mtwos=0
        
        do i=1,no_shells
            ps=part%shell(i)
            if (ps == 1) then
                ones=ones+1
            end if
            if (ps == -1) then
                mones=mones+1
            end if    
            if (ps == 2) then
                twos=twos+1
            end if
            if (ps == -2) then
                mtwos=mtwos+1
            end if
        end do
        if (ones > 0 .or. mones  > 2 .or. twos > 0 .or. mtwos  >  1) then
            oanalyse_part_optm=.true.
            return
        end if
                
        if  (ones == 0 .and. mones == 2 .and. twos == 0 .and. mtwos == 0) then
            return
        end if
        if  (ones == 0 .and. mones == 0 .and. twos == 0 .and. mtwos == 1) then
            return
        end if
        
        oanalyse_part_optm=.true.
        
    end function oanalyse_part_optm
    
    pure function oanalyse_part_optp(part)
        implicit none
        logical:: oanalyse_part_optp
        type(spartition),intent(in)::part
        integer:: ones,twos,mones,mtwos,i,ps

        oanalyse_part_optp=.false.
        
        if (part > 2 .or. part < 0) then
            oanalyse_part_optp=.true.
            return
        end if
         
        ones=0
        twos=0
        mones=0
        mtwos=0
        
        do i=1,no_shells
            ps=part%shell(i)
            if (ps == 1) then
                ones=ones+1
            end if
            if (ps == -1) then
                mones=mones+1
            end if    
            if (ps == 2) then
                twos=twos+1
            end if
            if (ps == -2) then
                mtwos=mtwos+1
            end if
        end do
        if (ones > 2 .or. mones  > 0 .or. twos > 1 .or. mtwos  >  0) then
            oanalyse_part_optp=.true.
            return
        end if
                
        if  (ones == 2 .and. mones == 0 .and. twos == 0 .and. mtwos == 0) then
            return
        end if
        if  (ones == 0 .and. mones == 0 .and. twos == 1 .and. mtwos == 0) then
            return
        end if
        
        oanalyse_part_optp=.true.
        
    end function oanalyse_part_optp
    
    pure function tptm2_analyse_part_op(part)
        implicit none
        logical:: tptm2_analyse_part_op
        type(spartition),intent(in)::part
        integer:: ones,twos,mones,mtwos,i,ps

        tptm2_analyse_part_op=.false.
        
        if (part > 1 .or. part < -2) then
            tptm2_analyse_part_op=.true.
            return
        end if
         
        ones=0
        twos=0
        mones=0
        mtwos=0
        
        do i=1,no_shells
            ps=part%shell(i)
            if (ps == 1) then
                ones=ones+1
            end if
            if (ps == -1) then
                mones=mones+1
            end if    
            if (ps == 2) then
                twos=twos+1
            end if
            if (ps == -2) then
                mtwos=mtwos+1
            end if
        end do
        if (ones > 1 .or. mones  > 2 .or. twos > 0 .or. mtwos  >  1) then
            tptm2_analyse_part_op=.true.
            return
        end if
        
        if (ones ==1 .and. mones ==2 .and. twos ==0 .and. mtwos ==0) return
        if (ones ==1 .and. mones ==0 .and. twos ==0 .and. mtwos ==1) return
        if (ones ==0 .and. mones ==1 .and. twos ==0 .and. mtwos ==0) return
                
        tptm2_analyse_part_op=.true.
        
    end function tptm2_analyse_part_op 
    
    pure function tp2tm_analyse_part_op(part)
        implicit none
        logical:: tp2tm_analyse_part_op
        type(spartition),intent(in)::part
        integer:: ones,twos,mones,mtwos,i,ps

        tp2tm_analyse_part_op=.false.
        
        if (part > 2 .or. part < -1) then
            tp2tm_analyse_part_op=.true.
            return
        end if
         
        ones=0
        twos=0
        mones=0
        mtwos=0
        
        do i=1,no_shells
            ps=part%shell(i)
            if (ps == 1) then
                ones=ones+1
            end if
            if (ps == -1) then
                mones=mones+1
            end if    
            if (ps == 2) then
                twos=twos+1
            end if
            if (ps == -2) then
                mtwos=mtwos+1
            end if
        end do
        if (ones > 2 .or. mones  > 1 .or. twos > 1 .or. mtwos  >  0) then
            tp2tm_analyse_part_op=.true.
            return
        end if
        
        if (ones ==2 .and. mones ==1 .and. twos ==0 .and. mtwos ==0) return
        if (ones ==0 .and. mones ==1 .and. twos ==1 .and. mtwos ==0) return
        if (ones ==1 .and. mones ==0 .and. twos ==0 .and. mtwos ==0) return
                
        tp2tm_analyse_part_op=.true.
        
    end function tp2tm_analyse_part_op 
    
    
end module Opm
        
module OpmSupport

use Parameters
use ProjParameters
use MbmeParameters
use Opm
use Files
use Partition
use OpPartition
use InputNuShell

implicit none

contains
    
    subroutine oread_rhs_base(base,k)
        implicit none
        integer,intent(in):: base,k
        integer::pindx,jj,tt,np,i,ibf,Jmax,Tmax,nJT
        dim_rhs=0
        do np=1,no_spartr        
           read(base+k,err=20,end=21) pindx,nJT,jj,tt,ibf,Jmax,Tmax      
           read(base+k,err=20,end=21) spartr(pindx)
            if (dim_lhs+ibf > notmvr) then
                print * ,' Increase max_total_mvr in MbmeParameters or'
                print * ,' SET TOTALMVECS=I where I >>',dim_lhs+ibf
                stop
            end if
            part_offsetr(pindx)=dim_rhs
            dim_partr(pindx)=0
            if (ibf /= 0) read(base+k,err=20,end=21) (mvec_rhs(i)%fn(1:gwords()),i=1+dim_rhs,ibf+dim_rhs)
            if (nJT == 0) go to 10
            part_offsetr(pindx)=dim_rhs
            dim_partr(pindx)=ibf
            dim_rhs=dim_rhs+ibf
10      end do
        return
20      print *, ' Error reading rhs mvectors '
        stop
21      print *, ' EOF reading rhs mvectors '
        stop
    end subroutine oread_rhs_base


end module OpmSupport

!  NuOpm.f90 
!
!  FUNCTIONS:
!  NuOpm      - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: NuOpm
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program NuOpm
    
    use Parameters
    use MbmeParameters
    use OpmParameters
    use OperParameters
    use ProjParameters
    use InputNuShell
    use OutputNuShell
    use Shells
    use ShellsX
    use Mscheme
    use Files
    use Opm
    use OpmSupport
    use Extra

    implicit none

    ! Variables
    
    real:: t1,t2,et1,et2,ms
    integer:: hr,mn,sc

    integer:: k,npt,k1,k2,i
    integer::nthr,OMP_get_max_threads    
    character(len=6):: inExt1,inExt2,Nucleus1,Nucleus2
    character(len=10):: btime
    character(len=8):: bdate
    
    integer:: max_cJ21,max_cJ22,min_cJ21,min_cJ22,max_cT21,max_cT22,min_cT21,min_cT22
    integer:: status,res
    character(len=14)::buffer
    character(len=4):: op
    character(len=1):: opf
    character(len=10):: chpart='PARTITIONS',chnopart,chnoOpm,chnotmvr,chtmvr='TOTALMVECS',chnoprojmvec
    character(len=6)::chOpm='OPMDIM'
    character(len=8)::chprojmvec='MVECTORS'
    
        res=get_env_v(chpart,chnopart)
        if (res /= 0) then
            read(chnopart,*) nopart
        else
            nopart=max_no_partitions
        end if
        allocate(spart(nopart))
        allocate(spartr(nopart))
        allocate(part_offset(nopart))
        allocate(part_offsetr(nopart))
        allocate(dim_part(nopart))
        allocate(dim_partr(nopart))
        
        res=get_env_v(chtmvr,chnotmvr)
        if (res /= 0) then
            read(chnotmvr,*) notmvr
        else
            notmvr=max_total_mvr
        end if
        allocate(mvec_rhs(notmvr))
        
        res=get_env_v(chprojmvec,chnoprojmvec)
        if (res /= 0) then
            read(chnoprojmvec,*) noprojmvec
        else
            noprojmvec=max_proj_mvec
        end if
        allocate(mvec_lhs(noprojmvec))
        
        res=get_env_v(chOpm,chnoOpm)
        if (res /= 0) then
            read(chnoOpm,*) noOpm
        else
            noOpm=max_Opm
        end if
        allocate(Opmidx(noOpm))
        allocate(Opm3idx(noOpm))
        allocate(Opmidxx(noOpm))
    
    ! Body of NuOpm
    call get_cmd_arg(1,buffer,status) 

    call output_header(6)  
    if (status == -1 ) then  
        print *, ' NuOpm Program'
        print *
        print *, ' Initial Nucleus'
    end if
    call input_nucleus(inExt2)
    if (status == -1 ) print *, ' Final Nucleus'
    call output_welcome(6)    
    call input_nucleus_final(inExt1)
    call input(inExt1)
    
    Nucleus1=Nucleus
    max_cJ21=max_cJ2
    max_cT21=max_cT2
    min_cJ21=min_cJ2
    min_cT21=min_cT2
    nmp1=nmp
    call output_shells(6)    
    call input(inExt2)
    Nucleus2=Nucleus
    max_cJ22=max_cJ2
    max_cT22=max_cT2
    min_cJ22=min_cJ2
    min_cT22=min_cT2
    nmp2=nmp
    nn=no_cI
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
    call output_end_setup(6)        
    call cpu_time(t1)
    call date_and_time(bdate,btime)
    call input_op(op,opf)
    op=adjustr(op)
    if (op=='   I'.or.op=='   i') then
    maxdmj=0
    mindmj=-maxdmj
    maxdmt=0
    mindmt=-maxdmt
    else if (op=='  a+'.or.op=='  a-')then 
    maxdmj=max_sj2u
    mindmj=-maxdmj
    maxdmt=1
    mindmt=-1
    else if (op=='a+a+'.or.op=='a-a-') then
    maxdmj=2*max_sj2u
    mindmj=-maxdmj
    maxdmt=2
    mindmt=-2
    else if (op=='a+a-') then
    maxdmj=2*max_sj2u
    mindmj=-maxdmj
    maxdmt=2
    mindmt=-2
    else if (shellxt.and.op=='a+--') then
    maxdmj=max_sj2u
    mindmj=-maxdmj
    maxdmt=1
    mindmt=-1
    else if (shellxt.and.op=='a++-') then
    maxdmj=max_sj2u
    mindmj=-maxdmj
    maxdmt=1
    mindmt=-1
    else if (op=='   N'.or.op=='   n') then
    maxdmj=0
    mindmj=-maxdmj
    maxdmt=0
    mindmt=-maxdmt
    end if
    
    if (output_control > 6) open(unit=60,file='op.txt',action='write')
    Nucleus=Nucleus1
    max_cJ2=max_cJ21
    max_cT2=max_cT21
    min_cJ2=min_cJ21
    min_cT2=min_cT21
    nmp=nmp1

    call open_files_read(no_spart,'.nba',100,dt)
    if (no_spart > nopart) then
        print *, ' Increase max_no_partitions in Parameeters or'
        print *, ' SET PARTITIONS=I where I >>',no_spart
        stop
    end if
    Nucleus=Nucleus2
    max_cJ2=max_cJ22
    max_cT2=max_cT22
    min_cJ2=min_cJ22
    min_cT2=min_cT22
    nmp=nmp2
    call open_files_read(no_spartr,'.nba',300,dt)
    if (no_spartr > nopart) then
        print *, ' Increase max_no_partitions in Parameeters or'
        print *, ' SET PARTITIONS=I where I >>',no_spartr
        stop
    end if
    Nucleus2=adjustr(Nucleus2)
    opfile(1:1)=' '
    opfile(2:6)=Nucleus2(1:5)
    opfile(7:7)=opf
    open(unit=2,file=opfile//'.mpt' ,form='unformatted', action='readwrite')
    
    if (output_control > 0) open(unit=10,file='NuOpm'//'.txt',action='write')
    k=0
    k2=0
    do j=min_cJ22,max_cJ22+2,2
        do t=min_cT22,max_cT2+dt,2
            if (j <= max_cJ22  .and. t <= max_cT22) then
                if (output_control > 0) write(10,*) ' 2*J2, 2*T2',j,t
                if (output_control > 6) write(60,*) ' 2*J2, 2*T2',j,t
                call oread_rhs_base(300,k)
                if (shellxt)use_isospin=.true.
                do delmj=mindmj,maxdmj,2
                do delmt=mindmt,maxdmt,2
                j1=j+delmj
                if (use_isospin) then
                    it1=t+delmt
                else
                    it1=t
                end if
                if (output_control > 0) write(10,*) ' 2*J1, 2*T1',j1,it1
                if (output_control > 6) write(60,*) ' 2*J1, 2*T1',j1,it1
                if (j1 <= max_cJ21 .and. j1>= 0 .and. it1 <= max_cT21 .and. it1>=0) then
                    k1=((j1-min_cJ21)/2)*((max_cT21+dt-min_cT21)/2+1)+((it1-min_cT21)/2)
                    twrite=t
                    if (.not.use_isospin) twrite=abs(nmp2)
                    twrite1=it1
                    if (.not.use_isospin) twrite1=abs(nmp1)
                    if (op=='   I'.or.op=='   i') call Idm_op(100,500,k,k1,k2)
                    if (op=='  a+') call Apm_op(100,500,k,k1,k2)
                    if (op=='  a-') call Amm_op(100,500,k,k1,k2)
                    if (op=='a-a-') call Tmm_op(100,500,k,k1,k2)
                    if (op=='a+a+') call Tpm_op(100,500,k,k1,k2)
                    if (shellxt.and.op=='a+--') call Tpm2_op(100,500,k,k1,k2)
                    if (shellxt.and.op=='a++-') call Tp2m_op(100,500,k,k1,k2)
                    if (op=='a+a-'.or.op=='   N'.or.op=='   n')&
                                 call Trm_op(100,500,k,k1,k2)
                    rewind(100+k1)
                    read(100+k1)
                end if
                k2=k2+1
                end do
                end do
            end if
            k=k+1
        end do
100 end do
    indexo=index
    if (no_disk/=0) then
        indexo=indexo+no_disk*max_Opm
    end if
    
    print * ,' No of OPM matrix elements ',indexo,proton,neutron,op
    open(unit=1,file=opfile//'.trm' ,form='unformatted', action='write')
    if (output_control>4) open(unit=11,file=opfile//'.mrt', action='write')
    write(1) indexo,proton,neutron,op
    if (output_control>4) write(11,*) indexo,proton,neutron,op
    do i=1,index
    write(1) Opmidx(i)(1:12),Opmidxx(i)  !,Opm3idx(i)(1:6)
    if (output_control>4) write(11,*) ichara3(Opmidx(i),12),Opmidxx(i)  !,ichara3(Opm3idx(i),6)
    end do
    if (no_disk/=0) then
    rewind(unit=2)
    do k=1,no_disk
    do i=1,max_Opm
    read(2) Opmidx(i)(1:12),Opmidxx(i)  !,Opm3idx(i)(1:6)
    write(1) Opmidx(i)(1:12),Opmidxx(i)  !,Opm3idx(i)(1:6)
    if (output_control>4) write(11,*) ichara3(Opmidx(i),12),Opmidxx(i)  !,ichara3(Opm3idx(i),6)
    end do
    end do    
    end if
    close(unit=1)
    
    deallocate(spart)
    deallocate(spartr)
    deallocate(part_offset)
    deallocate(part_offsetr)
    deallocate(dim_part)
    deallocate(dim_partr)
    deallocate(mvec_rhs)
    deallocate(mvec_lhs)
    deallocate(Opmidx)
    deallocate(Opm3idx)
    deallocate(Opmidxx)
    
    read(btime,'(i2,i2,i2,f4.3)') hr,mn,sc,ms
    et1=hr*3600.+mn*60.+sc*1. +ms
    call date_and_time(bdate,btime)
    call cpu_time(t2)
    read(btime,'(i2,i2,i2,f4.3)') hr,mn,sc,ms
    et2=hr*3600.+mn*60.+sc*1. +ms
    t1=t2-t1
    et1=et2-et1
    if (et1<0.0) et1=et1+24.*3600.
    if (status==-1) call output_completed(6,'NuOpd      ')
    call output_time('NuOpm',t1,et1)
    print *
    
    end program NuOpm
    

