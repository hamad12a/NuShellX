
module TransSupport

use Parameters
use LanczosParameters
use Partition
use OperParameters

implicit none
type(spartition),dimension(:),allocatable::meop
real,dimension(:),allocatable::me
real(kind=rl),dimension(:,:,:),allocatable::trampa
real(kind=rl),dimension(:),allocatable:: rvalue,lvalue
integer:: nome,j,twrite,j1,twrite1,delja,delta,delmja,delmta
integer:: deljb,deltb,delmjb,delmtb,nomex
logical:: ob,tc,pc
integer:: mindj,maxdj,mindt,maxdt,mindmj,maxdmj,mindmt,maxdmt,ni,nf
character(len=6):: Nucleus1,Nucleus2
character(len=7)::mbmemf
real(kind=rc):: ep,en,glp,gln,gsp,gsn,mass,charge,factor
character(len=13):: inExt1,inExt2

contains
    
    subroutine read_ordered_me(onetwo)
        implicit none
        integer,intent(in):: onetwo
        integer:: i,pdim
        if (mbmemf=='       ') return
        call set_ptdim(8)
        open(unit=53, file=mbmemf//'.meo',err=25)
        read(53,*,err=30) ep,en,glp,gln,gsp,gsn,mass,charge
        read(53,*,err=30) mindj,maxdj,mindt,maxdt,mindmj,maxdmj,mindmt,maxdmt,ni,nf,factor
        if (onetwo==1) then
            i=0
        else
            i=nome
        end if
        do 
            if (i+1 > nomex) then
                print *, ' Increase max_me in OrderParameters or'
                print *, ' SET NOME=I where I >>',i+1
                stop
            end if
            i=i+1
            read(53,*,err=30, end=20) meop(i)%shell(1:8),me(i)
            if(i==1) then
                if (meop(i)%shell(1) == 0 .and. meop(i)%shell(2) /= 0 .and. &
                    meop(i)%shell(3) /= 0 .and. meop(i)%shell(4) == 0) ob=.true.
                if (meop(i)%shell(1) == 0 .and. meop(i)%shell(2) /= 0 .and. &
                    meop(i)%shell(3) == 0 .and. meop(i)%shell(4) == 0) pc=.true.
                if (meop(i)%shell(1) /= 0 .and. meop(i)%shell(2) /= 0 .and. &
                    meop(i)%shell(3) == 0 .and. meop(i)%shell(4) == 0) tc=.true.
            end if
        end do
20      nome=i-1
        close(unit=53)
        return
25      print *, ' Error opening ',mbmemf//'.meo'
        stop  
30      print *, ' Error reading ',mbmemf//'.meo'
        stop  
    end subroutine read_ordered_me
       
end module TransSupport 

module CalcTrans

use TransSupport
use ClebschGordan
use OperParameters
use Shells
use InputNushell
use TLMB
use Parameters

implicit none

real(kind=rc),parameter:: nu0=0.0241126

contains

    subroutine Trans(isospin,one,two)
        implicit none
        logical,intent(in)::isospin,one,two
        real(kind=rc):: isofact0,isofact1,hbarw,nu,nu12m1,sumE,sumM,sumF,sumG,optz,omtz,betat
        real(kind=rc):: Elp,Eln,Mps,Mns,Mpl,Mnl,F1,G1,three,Q2,MD,par,radial,nulm,nulm1,S2,eep,een
        real(kind=rc):: C0,C1,CP,factc,factd,sqrtm2
        real(kind=rc):: nine1,nine2,mosh1,mosh2,mosh3,ampl,tampl,mosh4,xampl
        real(kind=rc),dimension(:),allocatable:: amplv,mvectv
        integer:: li,lf,nm,loop,i,comp,sh1,sh2,sh3,sh4,j12,j34,l12,l34,cn12,cn34
        integer:: cn1234,l1234,della,m12,m34,phase1,phase2
        integer,dimension(2):: delj,delmj
        character(len=1):: chsign
        character(len=1),dimension(10):: lch
        real(kind=rc),dimension(2):: facta, factb
        integer,dimension(:),allocatable:: csj2,csn,csl
        
        data lch /'s','p','d','f','g','h','i','j','k','l'/
  
        allocate(csj2(0:no_shells))
        allocate(csn(0:no_shells))
        allocate(csl(0:no_shells))
        allocate(amplv(nome))
        allocate(mvectv(nome))
        csj2=0
        csn=0
        csl=0
        csj2(1:no_shells)=sj2(1:no_shells)
        csn(1:no_shells)=sn(1:no_shells)
        csl(1:no_shells)=sl(1:no_shells)
       
        sqrtm2=real(1)/sqrt(real(2))
        three=sqrt(real(3))
        
        write(500,*) ' The relative spectroscopic factors given here can be directly'
        write(500,*) ' compared with G(SU(3))^2 from the SU(3) Shell Model. If you'
        write(500,*) ' divide the spectroscopic factors given here by G(SU(3))^2 and'
        write(500,*) ' multiply by 100 you get the % fraction of the SU(3) limit.'
        write(500,*) ' See Anyas-Weiss et al PHYSICS REPORTS (Section C of Physics '
        write(500,*) ' Letters)12,no.3 (1974) 201—272 for a formula for G(SU(3)).'
        
        delja=abs(delja)
        write(500,*) ' Initial  Nucleus ',Nucleus2
        write(500,*) ' Final    Nucleus ',Nucleus1
        write(500,*) ' Initial  2*J,2*T ',j,twrite
        write(500,*) ' Final    2*J,2*T ',j1,twrite1
        write(500,*) ' Spec Amp 2*J,2*T ',delja,abs(twrite-twrite1)
        write(500,*)
            
            do li=1,ni
                do lf=1,nf
                    mvectv=real(0)
                    amplv=real(0)
                    write(500,'(5x,a22,f8.3,a4)') ' Initial level energy ',rvalue(li),' MeV'
                    write(500,'(5x,a22,f8.3,a4)') ' Final   level energy ',lvalue(lf),' MeV'
                    write(500,*)
                    write(500,*) 's1 s2 s3 s4 na la nb lb  N  L  amplitude  amplitude^2  Sum'
                    tampl=real(0,rc)
                    ampl=real(0,rc)
                        do nm=1,nome
                        if (trampa(li,lf,nm)==0.0) go to 300
                        sh1=meop(nm)%shell(1)
                        sh2=meop(nm)%shell(2)
                        sh3=meop(nm)%shell(3)
                        sh4=meop(nm)%shell(4)
                        j12=meop(nm)%shell(5)
                        j34=meop(nm)%shell(7)
                        m12=2
                        m34=2
                        l12=j12
                        l34=j34
                        phase1=1
                        phase2=1
                        if (sh3==sh4) factd=factd*sqrtm2 
                        if (sh1==0.or.sh2==0) then
                            l12=2*max(csl(sh1),csl(sh2))
                            m12=1
                            phase1=(-1)**((1-max(csj2(sh1),csj2(sh2))+2*max(csl(sh1),csl(sh2)))/2)
                        end if    
                        if (sh3==0.or.sh4==0) then
                            l34=2*max(csl(sh3),csl(sh4))
                            m34=1
                            phase2=(-1)**((1-max(csj2(sh3),csj2(sh4))+2*max(csl(sh3),csl(sh4)))/2)
                        end if
                        nine1=sqrt(real((csj2(sh1)+1)*(csj2(sh2)+1)*(j12+1),rc))*&
                            ninej(mod(csj2(sh1),2),2*csl(sh1),csj2(sh1), &
                                  mod(csj2(sh2),2),2*csl(sh2),csj2(sh2), &
                        mod(mod(csj2(sh1),2)+mod(csj2(sh2),2),2),l12,j12)
                        factc=real(1)
                        if (sh1==sh2) factc=factc*sqrtm2                                      
                        nine2=sqrt(real((csj2(sh3)+1)*(csj2(sh4)+1)*(j34+1),rc))*&
                            ninej(mod(csj2(sh3),2),2*csl(sh3),csj2(sh3), &
                                  mod(csj2(sh4),2),2*csl(sh4),csj2(sh4), &
                        mod(mod(csj2(sh3),2)+mod(csj2(sh4),2),2),l34,j34)
                        factd=real(1)
                        if (sh3==sh4) factd=factd*sqrtm2                                      
                        cn12=(2*csn(sh1)+csl(sh1)+2*csn(sh2)+csl(sh2)-l12/2)/2
                        mosh1=moshinsky(cn12,l12/2,csn(sh1),csl(sh1),csn(sh2),csl(sh2),1,1)
                        cn34=(2*csn(sh3)+csl(sh3)+2*csn(sh4)+csl(sh4)-l34/2)/2
                        mosh2=moshinsky(cn34,l34/2,csn(sh3),csl(sh3),csn(sh4),csl(sh4),1,1)
                        della=2*(delja/2)
                        if (btest(delja,0).and.btest((l12+l34-della)/2,0)) then
                            della=della+2 
                        end if
                        cn1234=(2*cn12+l12/2+2*cn34+l34/2-della/2)/2
                        mosh3=moshinsky(cn1234,della/2,cn12,l12/2,cn34,l34/2,m12,m34)
                        ampl=factc*factd*nine1*nine2*mosh1*mosh2*mosh3*phase1*phase2
                        if (factor==1.0) then
                            ampl=ampl*(real(-1,rc)**(csn(sh1)+csn(sh2)+csn(sh3)+csn(sh4)))
                        end if
                        ampl=ampl*real(2)
                        if (ampl/=0.0) then
                        if ((sh1==0.or.sh2==0.or.sh1<=sh2).and.&
                            (sh3==0.or.sh4==0.or.sh3<=sh4)) then
                            amplv(nm)=ampl
                            mvectv(nm)=trampa(li,lf,nm)
                           
                            ampl=trampa(li,lf,nm)*ampl
                            tampl=tampl+ampl
                            write(500,'(10i3,2x,3(e10.4,2x))') sh1,sh2,sh3,&
                                sh4,cn12,l12/2,cn34,l34/2,&
                                cn1234,della/2,ampl,ampl**2,tampl
                        end if
                        end if
300                 end do
                    tampl=dot_product(amplv,mvectv)**2
                    write(500,*)
                    write(500,'(a52,2x,1(e10.4,2x))') ' Relative Spectroscopic Factor (independent of N,L) ',tampl
                    write(500,*)
                end do
            end do
     
    end subroutine trans
    
end module CalcTrans

!  NuClx.f90 
!
!  FUNCTIONS:
!  NuClx      - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: NuClx
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program NuClx
    
    use LanczosParameters
    use InputNuShell
    use OutputNuShell
    use MvecParameters
    use Mscheme
    use Files
    use OrderParameters
    use TransSupport
    use CalcTrans
    use Extra

    implicit none

    ! Variables
    
    real:: t1,t2,et1,et2,mns
    integer:: hr,mn,sc
    integer:: k,k1,k2,ll,indx,nolevel,inx,delt,lld,delj,delmj,delmt,dt,t,it1
    integer:: min_DJ,max_DJ,min_DT,max_DT,ndummy,no_levul,no_levur
    logical:: first,one,two,done
    
    real(kind=8)::clebsch
    
    character(len=10):: btime
    character(len=8):: bdate
    character(len=1):: xtra1,xtra2
    
    integer:: status,res,signl
    character(len=14)::buffer
    character(len=10)::chnolevel,chnome
    character(len=6)::chlevel='LEVELS',in1,in2
    character(len=4)::chme='NOME'

    character(len=15):: filename
    character(len=5):: file

        res=get_env_v(chlevel,chnolevel)
        if (res /= 0) then
            read(chnolevel,*) nolevel
        else
            nolevel=max_no_level
        end if
                
        res=get_env_v(chme,chnome)
        if (res /= 0) then
            read(chnome,*) nomex
        else
            nomex=max_me
        end if
        
    allocate(meop(nomex))
    allocate(me(nomex))
    
    call get_cmd_arg(1,buffer,status) 
    ! Body of NuClx

    call output_headerX(6) 
    if (status == -1) then   
        print *, ' NuClx Program'
        print *
        print *, ' Initial Nucleus '
        print *, ' Final   Nucleus '
        print *, ' Transition File 1'
    end if
    call input_Nucleus(in2)
    call input_nucleus_final(in1)   
    call input_trans_one(inExt1)
   
    call output_welcome(6)    
    call input(in1)
    
    call output_shells(6)    
    call input(in2)

    call output_shells(6)    
    call check_shell_order()
    call setup_shells()
    if (use_isospin) then
        dt=2
    else
        dt=0
    end  if
    call setup_mscheme()
    call setup_file_ext()
    mbmemf=inExt1(1:7)
    call read_ordered_me(1)
    min_DJ=mindj
    Min_DT=mindt
    Max_DJ=maxdj
    Max_Dt=maxdt        
    call output_end_setup(6)
    call cpu_time(t1)
    call date_and_time(bdate,btime)
    open(unit=500,file=inExt1//'.lp',action='write',err=55)
    call output_headerX(500)
    write(500,*)
    open(unit=601,err=400,file=inExt1//'.trb',form='unformatted')
10  read(601,end=40) j,t,j1,it1,twrite,twrite1,ni,nf
    read(601,end=40) min_DJ,max_DJ,min_DT,max_DT,delmj,delmt
    allocate(trampa(ni,nf,nome))
    allocate(lvalue(nf))
    allocate(rvalue(ni))
    do
    delja=0
    deljb=0
    delta=0
    deltb=0
    delmta=0
    delmtb=0
    delmja=0
    delmjb=0
    nome=0
    read(601,end=40) signl 
    if (signl==0) go to 30
    read(601,end=40) Nucleus2
    read(601,end=40) Nucleus1
    read(601,end=40) j,twrite
    read(601,end=40) j1,twrite1
    read(601,end=40) delja,delta,delmja,delmta,use_isospin
    delj=delja
    read(601,end=40) no_levul,no_levur,lvalue(1:no_levul),rvalue(1:no_levur)
    trampa(1:ni,1:nf,1:nome)=real(0,rl)
    do 
        read(601,err=40,end=40) indx,inx
        if (indx==-1) exit
        read(601,end=40)meop(inx)%shell(1:8)
        do ll=1,nf
            read(601,end=40) trampa(1:ni,ll,inx)
        end do
    end do
30  end do
40  nome=inx
    if (nome>0) call trans(use_isospin,one,two)
    close(unit=601)
    go to 60
50  print *, ' Error reading ',inext1
    go to 60 
55  print *, ' Error opening ',mbmemf//'.lp'
60  deallocate(trampa)
    deallocate(meop)
    deallocate(me)
    deallocate(rvalue)
    deallocate(lvalue)
    if (nome==0) go to 600
    read(btime,'(i2,i2,i2,f4.3)') hr,mn,sc,mns
    et1=hr*3600.+mn*60.+sc*1. +mns
    call date_and_time(bdate,btime)
    call cpu_time(t2)
    read(btime,'(i2,i2,i2,f4.3)') hr,mn,sc,mns
    et2=hr*3600.+mn*60.+sc*1. +mns
    et1=et2-et1
    t1=t2-t1
    if (et1<0.0) et1=et1+24.*3600.
    print *, ' Output file',inExt1//'.lp'
    if (output_control>0) then
        print *
        print *, ' Calculations successfully completed and data written to disk.'
        print *
    end if
    call output_time('NuClx',t1,et1)
    stop    
    print *
400 print *, ' File not found ',inExt1//'.trb'
    stop
500 print *, ' File not found ',inExt2//'.trb'
    stop
600 print *, ' No of trd`s = 0'
    stop
    
    end program NuClx
    