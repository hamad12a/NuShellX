! This file contains a\ll the modules specific to the generic program pTrans
! for NuShell and SunShell

module HarmonicOscillator

use Operparameters

implicit none

contains

      function hbw(a,z)
      implicit none
      real(kind=rc),intent(in)::a,z
      real(kind=rc)::hbw, bb ,con, hbc, amu, xmass
      real(kind=rc),dimension(1000):: b(1000)
      integer:: i
      
      data b /0.,0.,1.67,1.37,0.,1.88,1.74,0.,         &
       1.763,1.611,1.611,1.628,1.628,1.645,1.678,      &
       1.769,1.763,1.821,1.833,1.869,1.845,1.822,1.810,&
       1.813,1.793,1.802,1.804,1.827,1.825,1.835,1.848,&
       1.881,1.881,1.881,1.921,1.938,1.921,1.948,1.950,&
       1.963,960*0./
       
       if (a==real(0,rc)) then
            hbw=real(0,rc)
            return
       end if
       amu = 931.50
       hbc = 197.32
       xmass = ((a-z)*1.00728+z*1.00866)/a
       i = int(a)
       con = hbc*hbc/(amu*xmass)
       if (b(i) > real(0,rc)) then
       i = a
       bb = b(i)
       hbw = con/(bb*bb)
       return
       else
       hbw=45.*(a**(-1./3.))-25.*(a**(-2./3.))
       end if
       return
      end function hbw
      
end module HarmonicOscillator   

module TransSupport

use Parameters
use LanczosParameters
use Partition
use OperParameters

implicit none
type(spartition),dimension(:),allocatable::meop
real,dimension(:),allocatable::me
real(kind=rl),dimension(:,:,:),allocatable::trampa,trampb
real(kind=rl),dimension(:),allocatable:: rvalue,lvalue
integer:: nome,j,twrite,j1,twrite1,delja,delta,delmja,delmta
integer:: deljb,deltb,delmjb,delmtb,nomex
logical:: ob,tc,pc,iso
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
        iso=.false.
        if (mindt==2.or.maxdt==2) iso=.true.
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
            read(53,*,err=30, end=20) meop(i)%shell(1:8)
            me(i)=real(1)
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

module SingleParticle

use InputNuShell
use Shells
use OperParameters
use ClebschGordan

interface ML
    module procedure ML1,MLG
end interface
interface MS
    module procedure MS1,MSG
end interface

real(kind=rc),parameter:: pi4=3.14159265358979*real(4,rc)

contains 

! formulae based on B A Brown lectures 2005, MSU

    function E(lambda,shl,shr,efc)
        implicit none
        integer,intent(in):: lambda,shl,shr
        real(kind=rc),intent(in):: efc
        real(kind=rc)::phase1,phase2,three,ri,fact2,sqrtf,E
        logical:: btest
        
        phase1=real(1,rc)
        if (btest((sj2(shl)+1)/2,0)) phase1=-phase1
        phase2=real(1,rc)
        if (btest(sl(shr)+sl(shl)+lambda,0)) phase2=-phase2
        fact2=(real(1,rc)+phase2)/real(2,rc)
        if (fact2==real(0,rc)) then
            E=fact2
            return
        end if
        sqrtf=sqrt(real((sj2(shl)+1)*(sj2(shr)+1)*(2*lambda+1),rc)/pi4)
        three=threej(sj2(shl),2*lambda,sj2(shr),1,0,-1)
        ri=radial(sn(shl),sl(shl),sn(shr),sl(shr),lambda)
        E=phase1*fact2*sqrtf*three*ri*efc
        
    end function E
        
    function ML1(shl,shr,gl)
        implicit none
        integer,intent(in):: shl,shr
        real(kind=rc),intent(in):: gl
        real(kind=rc):: pif,phase,sqrtj,six,sqrtl,ML1
        logical:: btest
        
        ML1=real(0,rc)
        if (sn(shl)/=sn(shr)) return
        if (sl(shl)/=sl(shr)) return
        
        pif=sqrt(real(3,rc)/pi4)
        phase=real(1,rc)
        if (btest((2*sl(shl)+sj2(shr)+3)/2,0)) phase=-phase
        sqrtj=sqrt(real((sj2(shl)+1)*(sj2(shr)+1),rc))
        sqrtl=sqrt(real(sl(shl)*(sl(shl)+1)*(2*sl(shl)+1),rc))
        six=sixj(2*sl(shl),2*sl(shr),2,sj2(shr),sj2(shl),1)
        ML1=pif*phase*sqrtj*six*sqrtl*gl
        
    end function ML1
    
    function MS1(shl,shr,gs)
        implicit none
        integer,intent(in):: shl,shr
        real(kind=rc),intent(in):: gs
        real(kind=rc):: pif,phase,sqrtj,six,sqrts,MS1,sf
        logical:: btest
        
        MS1=real(0,rc)
        if (sn(shl)/=sn(shr)) return
        if (sl(shl)/=sl(shr)) return
        
        sf=sqrt(real(3,rc)/real(2,rc))
        pif=sqrt(real(3,rc)/pi4)
        phase=real(1,rc)
        if (btest((2*sl(shl)+sj2(shl)+3)/2,0)) phase=-phase
        sqrtj=sqrt(real((sj2(shl)+1)*(sj2(shr)+1),rc))
        sqrts=sqrt(real(3,rc)/real(2,rc))
        six=sixj(1,1,2,sj2(shr),sj2(shl),sl(shl))
        MS1=pif*phase*sqrtj*six*sqrts*sf*gs
        
    end function MS1
    
    function MLG(lambda,shl,shr,gl)
        implicit none
        integer,intent(in):: lambda,shl,shr
        real(kind=rc),intent(in):: gl
        real(kind=rc):: MLG,sqrtlm,phase1,sqrtj,six1,ri,phase2,sqrtll1,six2
        real(kind=rc):: phase3,sqrtll2,three,two
        logical:: btest
        
        two=real(2,rc)
        sqrtlm=sqrt(real(lambda*(2*lambda+1),rc))/real(lambda+1,rc)
        phase1=real(1,rc)
        if (btest((2*sl(shl)+1+sj2(shr)+2*lambda)/2,0)) phase1=-phase1
        sqrtj=sqrt(real((sj2(shl)+1)*(sj2(shr)+1),rc))
        six1=sixj(2*sl(shl),2*sl(shr),2*lambda,sj2(shr),sj2(shl),1)
        ri=radial(sn(shl),sl(shl),sn(shr),sl(shr),lambda-1)
        phase2=real(1,rc)
        if (btest(lambda+sl(shl)+sl(shr),0)) phase2=-phase2
        sqrtll1=sqrt(real((2*lambda+1)*sl(shr)*(sl(shr)+1)*(2*sl(shr)+1),rc))
        six2=sixj(2*lambda-2,2,2*lambda,2*sl(shr),2*sl(shl),2*sl(shr))
        phase3=(1,rc)
        if (btest(sl(shl),0)) phase3=-phase3
        sqrtll2=sqrt(real((2*sl(shl)+1)*(2*sl(shr)+1)*(2*lambda-1),rc)/pi4)
        three=threej(2*sl(shl),2*lambda-2,2*sl(shr),0,0,0)
        MLG=sqrtlm*phase1*sqrtj*six1*phase2*sqrtll1*six2*phase3*sqrtll2* &
                   two*three*ri*gl
                   
    end function MLG

    function MSG(lambda,shl,shr,gs)
        implicit none
        integer,intent(in):: lambda,shl,shr
        real(kind=rc),intent(in):: gs
        real(kind=rc):: MSG,sqrtlm,sqrtj,nine,phase,sqrtll,three,sf,ri
        logical:: btest
        
        sf=sqrt(real(3,rc)/real(2,rc))
        sqrtlm=sqrt(real(lambda*(2*lambda+1),rc))
        sqrtj=sqrt(real((sj2(shl)+1)*(sj2(shr)+1)*(2*lambda+1),rc))
        ri=radial(sn(shl),sl(shl),sn(shr),sl(shr),lambda-1)
        phase=(1,rc)
        if (btest(sl(shl),0)) phase=-phase
        sqrtll=sqrt(real((2*sl(shl)+1)*(2*sl(shr)+1)*(2*lambda-1),rc)/pi4)
        three=threej(2*sl(shl),2*lambda-2,2*sl(shr),0,0,0)
        nine=ninej(2*sl(shl),1,sj2(shl),2*sl(shr),1,sj2(shr),2*lambda-2,2,2*lambda)
        MSG=sqrtlm*sqrtj*nine*phase*sqrtll*three*sf*ri*gs
        
    end function MSG
    
    function BF(shl,shr)
        implicit none
        integer,intent(in):: shl,shr
        real(kind=rc):: BF
        logical:: btest
               
        if (sn(shl)/=sn(shr)) return
        if (sl(shl)/=sl(shr)) return
        if (sj2(shl)/=sj2(shr)) return
            BF=real(1,rc)
        
    end function BF
    
    function BGT(shl,shr,obdp)
        implicit none
        integer, intent(in):: shl,shr
        real(kind=rc):: BGT,sf,six,sqrtj,phase,two
        logical,intent(in):: obdp
        
        BGT=real(0,rc)
        if (sn(shl)/=sn(shr)) return
        if (sl(shl)/=sl(shr)) return
        two=real(2,rc)*sqrt(real(2,rc))
        sf=sqrt(real(3,rc)/real(2,rc))
        phase=real(1,rc)
        if (obdp.and.btest((sj2(shl)+sj2(shr)-2)/2,0)) phase=-phase
        if (btest((2*sl(shl)+sj2(shr)+3)/2,0)) phase=-phase
        sqrtj=sqrt(real((sj2(shl)+1)*(sj2(shr)+1),rc))
        six=sixj(1,1,2,sj2(shr),sj2(shl),2*sl(shl))
        BGT=phase*sqrtj*six*sf*two
        
    end function BGT
    
end module SingleParticle

module CalcTrans

use TransSupport
use ClebschGordan
use SingleParticle
use OperParameters
use Shells
use HarmonicOscillator
use Parameters

implicit none

real(kind=rc),parameter:: nu0=0.0241126
logical::obdp

contains

    subroutine Trans(isospin,one,two)
        implicit none
        logical,intent(in)::isospin,one,two
        real(kind=rc):: isofact0,isofact1,hbarw,nu,nu12m1,sumE,sumM,sumF,sumG,optz,omtz,betat
        real(kind=rc):: Elp,Eln,Mps,Mns,Mpl,Mnl,F1,G1,three,Q2,MD,par,radial,nulm,nulm1,S2,eep,een
        real(kind=rc):: C0,C1,CP,factc,factd
        integer:: li,lf,nm,loop,i,comp
        integer,dimension(2):: delj,delmj
        character(len=1):: chsign
        character(len=1),dimension(10):: lch
        real(kind=rc),dimension(2):: facta, factb
        
        data lch /'s','p','d','f','g','h','i','j','k','l'/
        
        if (ob) then
        hbarw=hbw(mass,charge)
        write(500,*)
        write(500,*) ' hbarw',real(hbarw),' MeV'
        write(500,*)
        write(500,'(1x,a23,6f8.3)') ' ep,en,glp,gln,gsp,gsn ',ep,en,glp,gln,gsp,gsn
        write(500,*)
        end if
        eep=(mass-charge)/mass
        een=-charge/mass
        if (one .and. two) then
             if (delja==deljb .and. delmja==delmjb) then
                loop=1
                facta(1)=reaL(1,rc)
                factb(1)=real(1,rc)
                delj(1)=delja
                delmj(1)=delmja
             else
                loop=2
                facta(1)=reaL(1,rc)
                factb(1)=real(0,rc)
                delj(1)=delja
                delmj(1)=delmja
                facta(2)=reaL(0,rc)
                factb(2)=real(1,rc)
                delj(2)=deljb
                delmj(2)=delmjb
             end if
        end if
        if (one  .and. .not.two) then
                loop=1
                facta(1)=reaL(1,rc)
                factb(1)=real(0,rc)
                delj(1)=delja
                delmj(1)=delmja
        end if
        if (two .and. .not.one) then
                loop=1
                facta(1)=reaL(0,rc)
                factb(1)=real(1,rc)
                delj(1)=deljb
                delmj(1)=delmjb
        end if   
        do i=1,loop
        write(500,*) ' Initial  Nucleus ',Nucleus2
        write(500,*) ' Final    Nucleus ',Nucleus1
        write(500,*) ' Initial  2*J,2*T ',j,twrite
        write(500,*) ' Final    2*J,2*T ',j1,twrite1
        write(500,*) ' Operator 2*J,2*T ',delj(i),delta,deltb
        write(500,*) ' Operator 2*M,2*Tz',delmj(i),delmta,delmtb
        if (ob) then                                  
            par=real(1,rc)
            if (btest((twrite1-delmtb)/2,0)) par=-par
            three=threej(twrite1,delta,twrite,-twrite1,delmta,twrite)/real(2,rc)
            isofact0=real(1,rc)
            if (isospin) isofact0=par*sqrt(real(2,rc))*three
            if (three==real(0,rc)) isofact0=real(0,rc)
            three=threej(twrite1,deltb,twrite,-twrite1,delmtb,twrite)/real(2,rc)
            isofact1=sqrt(real(3,rc))
            if (isospin) isofact1=par*sqrt(real(6,rc))*three
            if (three==real(0,rc)) isofact1=real(0,rc)
            nu=nu0*hbarw
            nu12m1=real(1,rc)/sqrt(nu)
            write(500,*)  
            if (delj(i)==0) nulm=real(1,rc)
            if (delj(i)==0) nulm1=real(0,rc)
            if (delj(i)==2) nulm=nu12m1                                       
            if (delj(i)==2) nulm1=real(1,rc)
            if (delj(i)>2)  nulm=(nu12m1)**(real(delj(i)/2,rc))     
            if (delj(i)>2)  nulm1=(nu12m1)**(real(delj(i)/2-1,rc))     
            do li=1,ni
                do lf=1,nf
                    write(500,'(5x,a22,f8.3,a4)') ' Initial level energy ',rvalue(li),' MeV'
                    write(500,'(5x,a22,f8.3,a4)') ' Final   level energy ',lvalue(lf),' MeV'
                    sumE=real(0,rc)
                    sumM=real(0,rc)
                    sumF=real(0,rc)
                    sumG=real(0,rc)
                    if (isospin) then 
                        factc=real(1,rc)
                        factd=real(1,rc)
                    else 
                        comp=0
                        do nm=1,nome
                            if (abs(trampa(li,lf,nm))<0.00001) comp=comp+1
                        end do
                        factc=real(1,rc)
                        factd=real(0,rc)
                        if (comp==nome) then
                            factc=real(0,rc)
                            factd=real(1,rc)
                        end if
                    end if 
                    do nm=1,nome
                        Mpl=real(0,rc)
                        Mnl=real(0,rc)
                        Mps=real(0,rc)
                        Mns=real(0,rc)
                        Elp=real(0,rc)
                        Eln=real(0,rc)
                        G1=real(0,rc)
                        F1=real(0,rc)
                        optz=(facta(i)*factc*trampa(li,lf,nm)*isofact0+  &
                              factb(i)*factd*trampb(li,lf,nm)*isofact1)
                        omtz=(facta(i)*factc*trampa(li,lf,nm)*isofact0-  &
                              factb(i)*factd*trampb(li,lf,nm)*isofact1)
                        betat=factb(i)*trampb(li,lf,nm)*isofact1
                        radial=factor
                        radial=radial**(sn(int(meop(nm)%shell(2)))+sn(int(meop(nm)%shell(3))))
                        if (delj(i)/=2) then
                            if (nucleon(int(meop(nm)%shell(2),4))=='i' .or. &
                               (nucleon(int(meop(nm)%shell(2),4))=='p' .and. &
                                nucleon(int(meop(nm)%shell(3),4))=='p')) then 
                            Elp=E(delj(i)/2,int(meop(nm)%shell(2),4),int(meop(nm)%shell(3),4),ep)&
                                    *nulm*radial
                            end if
                        else
                            if (nucleon(int(meop(nm)%shell(2),4))=='i' .or. &
                               (nucleon(int(meop(nm)%shell(2),4))=='p' .and. &
                                nucleon(int(meop(nm)%shell(3),4))=='p')) then 
                            Elp=E(delj(i)/2,int(meop(nm)%shell(2),4),int(meop(nm)%shell(3),4),eep)&
                                    *nulm*radial
                            end if
                        end if
                        if (delj(i)==0) then
                            if (nucleon(int(meop(nm)%shell(2),4))=='i' .or. &
                               (nucleon(int(meop(nm)%shell(2),4))=='p' .and. &
                                nucleon(int(meop(nm)%shell(3),4))=='n') .or. &
                               (nucleon(int(meop(nm)%shell(2),4))=='n' .and. &
                                nucleon(int(meop(nm)%shell(3),4))=='p')) then 
                            F1=BF(int(meop(nm)%shell(2),4),int(meop(nm)%shell(3),4))
                            end if
                        end if
                        if (delj(i)==2) then
                            if (nucleon(int(meop(nm)%shell(2),4))=='i' .or. &
                               (nucleon(int(meop(nm)%shell(2),4))=='p' .and. &
                                nucleon(int(meop(nm)%shell(3),4))=='n') .or. &
                               (nucleon(int(meop(nm)%shell(2),4))=='n' .and. &
                                nucleon(int(meop(nm)%shell(3),4))=='p')) then 
                            G1=BGT(int(meop(nm)%shell(2),4),int(meop(nm)%shell(3),4),obdp)
                            end if
                        end if
                        if (delj(i)/=2) then
                            if (nucleon(int(meop(nm)%shell(2),4))=='i' .or. &
                               (nucleon(int(meop(nm)%shell(2),4))=='n' .and. &
                                nucleon(int(meop(nm)%shell(3),4))=='n')) then 
                            Eln=E(delj(i)/2,int(meop(nm)%shell(2),4),int(meop(nm)%shell(3),4),en) &
                                    *nulm*radial
                            end if
                        else
                            if (nucleon(int(meop(nm)%shell(2),4))=='i' .or. &
                               (nucleon(int(meop(nm)%shell(2),4))=='n' .and. &
                                nucleon(int(meop(nm)%shell(3),4))=='n')) then 
                            Eln=E(delj(i)/2,int(meop(nm)%shell(2),4),int(meop(nm)%shell(3),4),een) &
                                     *nulm*radial
                            end if
                        end if
                        sumE=sumE+Elp*omtz+Eln*optz
                        sumG=sumG+G1*betat
                        sumF=sumF+F1*betat
                        
                        if (nucleon(int(meop(nm)%shell(2),4))=='i' .or. &
                               (nucleon(int(meop(nm)%shell(2),4))=='p' .and. &
                                nucleon(int(meop(nm)%shell(3),4))=='p')) then 
                        Mpl=ML(delj(i)/2,int(meop(nm)%shell(2),4),int(meop(nm)%shell(3),4),glp)&
                                    *nulm1
                        end if
                        if (delj(i)>2) Mpl=Mpl*radial
                        if (nucleon(int(meop(nm)%shell(2),4))=='i' .or. &
                               (nucleon(int(meop(nm)%shell(2),4))=='n' .and. &
                                nucleon(int(meop(nm)%shell(3),4))=='n')) then 
                        Mnl=ML(delj(i)/2,int(meop(nm)%shell(2),4),int(meop(nm)%shell(3),4),gln)&
                                    *nulm1
                        end if
                        if (delj(i)>2) Mnl=Mnl*radial
                        if (nucleon(int(meop(nm)%shell(2),4))=='i' .or. &
                               (nucleon(int(meop(nm)%shell(2),4))=='p' .and. &
                                nucleon(int(meop(nm)%shell(3),4))=='p')) then 
                        Mps=MS(delj(i)/2,int(meop(nm)%shell(2),4),int(meop(nm)%shell(3),4),gsp)&
                                    *nulm1
                        end if
                        if (delj(i)>2) Mps=Mps*radial
                        if (nucleon(int(meop(nm)%shell(2),4))=='i' .or. &
                               (nucleon(int(meop(nm)%shell(2),4))=='n' .and. &
                                nucleon(int(meop(nm)%shell(3),4))=='n')) then 
                        Mns=MS(delj(i)/2,int(meop(nm)%shell(2),4),int(meop(nm)%shell(3),4),gsn)&
                                    *nulm1
                        end if
                        if (delj(i)>2) Mns=Mns*radial
                        sumM=sumM+(Mpl+Mps)*omtz+(Mnl+Mns)*optz
                    end do
                    if (sumE /= real(0,rc)) then
                    chsign='+'
                    if (int(sign(1.0,sumE))<0) chsign='-'            
                    write(500,'(10x,a4,i1,a2,e10.4,a3,a3,i1,a6,a1)') &
                            ' BE(',delj(i)/2,') ',sumE*sumE/real(j+1,rc),' e2',' fm',delj(i) &
                                 , ' sign ', chsign
                    if (j1==j .and. twrite==twrite1 .and. li==lf .and. Nucleus1==Nucleus2 &
                        .and. lvalue(lf)==rvalue(li) .and. delj(i)/2==2) then
                        three=threej(j,delj(i),j,-j,0,j)*sqrt(real(2,rc)*pi4/real(5,rc))
                        Q2=sumE*three
                        write(500,'(10x,a19,e10.4,a6)') ' Quadrupole Moment ',Q2,' e fm2'
                    end if
                    end if
                    if (sumM /= real(0,rc)) then
                    chsign='+'
                    if (int(sign(1.0,sumM))<0) chsign='-'            
                    write(500,'(10x,a3,i1,a2,e10.4,a4,a3,i1,a6,a1)') &
                       ' M(',delj(i)/2,') ',sumM*sumM/real(j+1,rc),' uN2',' fm',delj(i)-2 &
                           , ' sign ',chsign
                    if (j1==j .and. twrite==twrite1 .and. li==lf .and. Nucleus1==Nucleus2 &
                        .and. lvalue(lf)==rvalue(li) .and. delj(i)/2==1) then
                        three=threej(j,delj(i),j,-j,0,j)*sqrt(real(2,rc)*pi4/real(3,rc))
                        MD=sumM*three
                        write(500,'(10x,a17,e10.4,a3)') ' Magnetic Moment ',MD,' uN'
                    end if  
                    end if          
                    if (sumG/=real(0,rc)) then
                        chsign='+'
                        if (int(sign(1.0,sumG))<0) chsign='-'            
                        write(500,'(10x,a6,e10.3,a6,a1)') ' B(GT)',sumG*sumG/real(j+1,rc) &
                                 , ' sign ', chsign
                    end if    
                    if (sumF/=real(0,rc)) then
                        chsign='+'
                        if (int(sign(1.0,sumF))<0) chsign='-'            
                        write(500,'(10x,a6,e10.3,a6,a1)') ' BF',sumF*sumF/real(j+1,rc) &
                                 , ' sign ', chsign
                    end if    
                    write(500,*)
                end do
            end do
        end if
        if (pc) then
            C0=real(1,rc)
            three=real(1,rc)
            isofact0=sqrt(real((delj(i)+1),rc))/sqrt(real((j1+1),rc)) &
                                *sqrt(real((delta+1),rc))
            if (isospin) then
                C0=cleb(twrite,delta,twrite1,twrite,delmta,twrite1)
                three=C0/sqrt(real((twrite1+1),rc))
                isofact0=three*isofact0
            end if
            if (three==real(0,rc)) isofact0=real(0,rc)
            do li=1,ni
                do lf=1,nf
                    write(500,'(5x,a22,f8.3,a4)') ' Initial level energy ',rvalue(li),' MeV'
                    write(500,'(5x,a22,f8.3,a4)') ' Final   level energy ',lvalue(lf),' MeV'
                    do nm=1,nome
                        S2=real(0,rc)
                        if (int(meop(nm)%shell(5),4)==delj(i)) then
                            if (int(meop(nm)%shell(6),4)==delta) then
                                optz=facta(i)*trampa(li,lf,nm)
                                S2=(isofact0*optz)
                            end if
                            if (S2/=real(0,rc)) then
                            chsign='+'
                            if (int(sign(1.0,S2))<0) chsign='-'     
                            write(500,'(10x,a29,i1,a1,i2,a3,a1,a2,e10.4,a6,a1)') &
                                ' Spectroscopic factor C2S   (',&
                                sn(int(meop(nm)%shell(2),4)),&
                                lch(sl(int(meop(nm)%shell(2),4))+1), &
                                sj2(int(meop(nm)%shell(2),4)),'/2)', &
                                nucleon(int(meop(nm)%shell(2),4)),'  ', & 
                                S2*S2,' sign ',chsign 
                            write(500,'(10x,a3,e10.4,a3,e10.4)')' C ',C0, ' S ',S2*S2/C0**2
                            end if
                        end if
                    end do
                    write(500,*)
                end do
            end do
        end if
        if (tc) then
            C0=real(1,rc)
            three=real(1,rc)
            isofact0=sqrt(real((delj(i)+1),rc))/sqrt(real((j1+1),rc)) &
                        *sqrt(real((delta+1),rc))
            if (isospin) then
                C0=cleb(twrite,delta,twrite1,twrite,delmta,twrite1)
                three=C0/sqrt(real((twrite1+1),rc))
                isofact0=three*isofact0
            end if
            if (three==real(0,rc)) isofact0=real(0,rc)
            C1=real(1,rc)
            three=real(1,rc)
            isofact1=sqrt(real((delj(i)+1),rc))/sqrt(real((j1+1),rc))* &
                        sqrt(real((deltb+1),rc))
            if (isospin) then
                C1=cleb(twrite,deltb,twrite1,twrite,delmtb,twrite1)
                three=C1/sqrt(real((twrite1+1),rc))
                isofact1=three*isofact1
            end if
            if (three==real(0,rc)) isofact1=real(0,rc)   
            do li=1,ni
                do lf=1,nf
                    write(500,'(5x,a22,f8.3,a4)') ' Initial level energy ',rvalue(li),' MeV'
                    write(500,'(5x,a22,f8.3,a4)') ' Final   level energy ',lvalue(lf),' MeV'
                    do nm=1,nome
                        S2=real(0,rc)
                        if (int(meop(nm)%shell(5),4)==delj(i) ) then
                            if (int(meop(nm)%shell(6),4)==delta) then
                                optz=facta(i)*trampa(li,lf,nm)
                                S2=(isofact0*optz)
                                CP=C0
                            end if
                            if (S2==real(0,rc)) then
                            if (int(meop(nm)%shell(6),4)==deltb) then
                                omtz=factb(i)*trampb(li,lf,nm)
                                S2=(isofact1*omtz)
                                CP=C1
                            end if
                            end if
                            if (S2/=real(0,rc)) then
                            chsign='+'
                            if (int(sign(1.0,S2))<0) chsign='-'     
                            write(500,'(10x,a29,i1,a1,i2,a3,a1,a3,i1,a1,i2,a3,a1,a2,e10.4,a6,a1)') &
                                ' Spectroscopic factor C2S   (',&
                                sn(int(meop(nm)%shell(1),4)),&
                                lch(sl(int(meop(nm)%shell(1),4))+1), &
                                sj2(int(meop(nm)%shell(1),4)),'/2)',&
                                nucleon(int(meop(nm)%shell(1),4)),'  (', & 
                                sn(int(meop(nm)%shell(2),4)),&
                                lch(sl(int(meop(nm)%shell(2),4))+1), &
                                sj2(int(meop(nm)%shell(2),4)),'/2)',&
                                nucleon(int(meop(nm)%shell(2),4)),'  ', & 
                                S2*S2,' sign ',chsign 
                            write(500,'(10x,a2,i1,a1,e10.4,a3,e10.4)') &
                                    ' C',meop(nm)%shell(6)/2,' ',CP,' S ',S2*S2/CP**2
                            end if
                        end if
                    end do
                    write(500,*)
                end do
            end do
         end if
         end do
     
    end subroutine trans
    
end module CalcTrans

!  NuTrx.f90 
!
!  FUNCTIONS:
!  NuTrx      - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: NuTrx
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program NuTrx
    
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
    integer:: k,k1,k2,ll,indx,nolevel,inx,inx2,delt,lld,delj,delmj,delmt,dt,t,it1
    integer:: min_DJ,max_DJ,min_DT,max_DT,ndummy,no_levul,no_levur,notrd
    logical:: first,one,two,done
    
    real(kind=8)::clebsch
    
    character(len=10):: btime
    character(len=8):: bdate
    character(len=1):: xtra1,xtra2
    
    integer:: status,res,signl,nome1
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
    ! Body of NuTrx
    notrd=0
    call output_headerX(6) 
    if (status == -1) then   
        print *, ' NuTrx Program'
        print *
        print *, ' Initial Nucleus '
        print *, ' Final   Nucleus '
        print *, ' Transition File 1'
        print *, ' Transition File 2'
    end if
    call input_Nucleus(in2)
    call input_nucleus_final(in1)   
    call input_trans_one(inExt1)
    call input_trans_two(inExt2)
   
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
    iso=.false.
    mbmemf=inExt1(1:7)
    nome1=0
    call read_ordered_me(1)
    nome1=nome
    mbmemf=inExt2(1:7)
    call read_ordered_me(2)
    min_DJ=mindj
    Min_DT=mindt
    Max_DJ=maxdj
    Max_Dt=maxdt        
    call output_end_setup(6)
    call cpu_time(t1)
    call date_and_time(bdate,btime)
    mbmemf=inExt1(1:7)
    open(unit=500,file=inExt1//'.lp',action='write',err=55)
    call output_headerX(500)
    write(500,*)
    first=.true.
    open(unit=601,err=400,file=inExt1//'.trb',form='unformatted')
    one=.true.
    two=.false. 
    if (.not.inExt2=='             ') then
        two=.true.
        open(unit=602,err=500,file=inExt2//'.trb',form='unformatted')
    end if
10  read(601,end=45) j,t,j1,it1,twrite,twrite1,ni,nf
    read(601,end=45) min_DJ,max_DJ,min_DT,max_DT,delmj,delmt
    if (first) then
        allocate(trampa(ni,nf,nome))
        allocate(trampb(ni,nf,nome))
        allocate(lvalue(nf))
        allocate(rvalue(ni))
        first=.false.
    end if    
    if (two) then
        read(602,end=45)
        read(602,end=45)
    end if
    do
    read(601,end=45) signl 
    if (two) read(602,end=45)
    if (signl==0) go to 42
    read(601,end=40) Nucleus2
    if (two) read(602,end=40)
    read(601,end=40) Nucleus1
    if (two) read(602,end=40)
    read(601,end=40) j,twrite
    if (two) read(602,end=40)
    read(601,end=40) j1,twrite1
    if (two) read(602,end=40)
    read(601,end=40) delja,delta,delmja,delmta,use_isospin,obdp
    deltb=0
    delj=abs(delja)
    if (two) read(602,end=40)
    read(601,end=40) no_levul,no_levur,lvalue(1:no_levul),rvalue(1:no_levur)
    if (two) read(602,end=40)
    if ((.not.two.and.(abs(min_DT)==2.or.abs(max_dT)==2)).or.iso) then
	one=.false.
	iso=.true.
	deljb=delja
	deltb=2
	delmjb=delmja
	delmtb=twrite1-twrite
    end if	
    trampa(1:ni,1:nf,1:nome)=real(0,rl)
    trampb(1:ni,1:nf,1:nome)=real(0,rl)
    nome=0
    do 
        read(601,err=30,end=30) indx,inx
        if (indx==-1) exit
        read(601,end=30)meop(inx)%shell(1:8)
        
        do ll=1,nf
            if (one)read(601,end=30) trampa(1:ni,ll,inx)
    	    if (iso)read(601,end=30) trampb(1:ni,ll,inx)
        end do
    end do
30  if (one.and.two) then
    do 
        read(602,err=40,end=40) indx,inx2
        if (indx==-1) exit
        read(602,end=40)meop(inx2+inx)%shell(1:8)
        
        do ll=1,nf
            read(602,end=40) trampa(1:ni,ll,inx2+inx)
        end do
    end do
    end if
40  nome=inx+inx2
    notrd=notrd+nome
    if (one.and.two.and.nome>0) call trans(use_isospin,.true.,.false.)
    if (iso.and.nome>0) call trans(use_isospin,.false.,.true.)
    if (.not.two.and.one.and.nome>0) call trans(use_isospin,.true.,.false.)
42  end do
45  close(unit=601)
    if (two.and.deltb==0) close(602)
    go to 60
50  print *, ' Error reading ',inext1
    go to 60 
55  print *, ' Error opening ',inext1//'.lp'
60  deallocate(trampa)
    deallocate(trampb)
    deallocate(meop)
    deallocate(me)
    deallocate(rvalue)
    deallocate(lvalue)
    if (notrd==0) go to 600
    read(btime,'(i2,i2,i2,f4.3)') hr,mn,sc,mns
    et1=hr*3600.+mn*60.+sc*1. +mns
    call date_and_time(bdate,btime)
    call cpu_time(t2)
    read(btime,'(i2,i2,i2,f4.3)') hr,mn,sc,mns
    et2=hr*3600.+mn*60.+sc*1. +mns
    et1=et2-et1
    t1=t2-t1
    if (et1<0.0) et1=et1+24.*3600.
    if (output_control>0) then
        print *
        print *, ' Calculations successfully completed and data written to disk',mbmemf//'.lp.'
        print *
    else 
        print *, ' Output file ',inExt1//'.lp'
    end if
    call output_time('NuTrx',t1,et1)    
    print *
    stop
400 print *, ' File not found ',inExt1//'.trb'
    stop
500 print *, ' File not found ',inExt2//'.trb'
    stop
600 print *, ' No of trd`s = 0'
    stop
    end program NuTrx
    