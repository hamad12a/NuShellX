module Tra

use Parameters
use ProjParameters
use MatrixParameters
use LanczosParameters
use OperParameters
use ClebschGordan
use InputNuShell
use TrdChar
use Files
use TriDiag
use SOrder
use DotProduct
use SubM

implicit none
integer,parameter:: noiter=300,max_trdb=1000,max_trda=1000
integer,parameter::  isl=3,isr=4,ilm=5,ijl=1,ijr=2,inl=6,inr=8,imr=10,iml=11,ilt=12
type(tablear),dimension(:,:),allocatable:: tablebh
integer(kind=2),dimension(:,:),allocatable:: tablebhn
real(kind=rkw),dimension(:),allocatable:: coefR,trdav,trdbv,coefL,ninejs,coefT
real(kind=rkw),dimension(:,:),allocatable:: BUFFER1,BUFFER2,BUFFER3
real(kind=rkw),dimension(:,:),allocatable:: coefX
character(len=8),dimension(:),allocatable::trda,trdb,prdb,basR,basL,nine
integer,dimension(:),allocatable:: posa,posb,posnvL,dimnvL,posnvR,dimnVR,dima,dimb,csj2
integer,dimension(:,:),allocatable:: posal,posbl
type(la),dimension(:),allocatable:: nplra,nplrb,nplrma,nplrmb
integer:: no_cA,no_cB,no_cC,no_cD,nwordsa,nwordsb
real(kind=rkw),dimension(:,:,:,:),allocatable:: tramp
integer:: delj,delt,mindj,maxdj,mindt,maxdt,mindmj,maxdmj,mindmt,maxdmt,ni,nf
character(len=7)::opfile1,opfile2,mbmemf
character(len=6)::opfile16,opfile26,NucleusR,NucleusL,Nucleuss
character(len=1):: xtraa,xtrab,xtrac,xtrad
integer,dimension(:,:),allocatable:: nJTa,nJTb,nJTc,nJTd
character(len=8),dimension(:),allocatable:: meop
character(len=2),dimension(:),allocatable:: meopp
integer:: min_cJ2a,max_cJ2a,no_sparta,no_spart,maxda,del_cJ2a
integer:: min_cJ2b,max_cJ2b,no_spartb,no_stateL=0,maxdb,del_cJ2b
integer:: min_cJ2c,max_cJ2c,no_spartc,maxdc,del_cJ2c
integer:: min_cJ2d,max_cJ2d,no_spartd,no_stateR=0,maxdd,del_cJ2d
real(kind=rkw),dimension(:),allocatable:: rvalue,lvalue
character(len=8):: ta,tb,pb,tp,bm,ba,rn,z1
character(len=4):: opa,opb,opc,opd
integer:: trdda,trddb,prddb,nined,trddat,trddbt,max_dimv
integer:: max_lma,max_lmb,max_lmd,max_lmc,basdR,basdL,max_lm
integer:: min_lma,min_lmb,min_lmd,min_lmc,min_lm,tpa,tpb,tpc,tpd
integer:: ajl,ajr,anpl,anpr,ashl,ashr,lma,lmb,lm,idx,lmx,lmn,lmnn,lmxx,lmp
integer:: bjl,bjr,bnpl,bnpr,J2R,J2L,ajrn,ajrx,ajln,ajlx
integer:: bjln,bjlx,bshl,bshr,bjrn,bjrx,no_level,no_levelR,no_levelL
integer:: cjl,cjr,cnpl,cnpr,cshl,cshr,cjrn,cjrx,cjln,cjlx
integer:: djl,djr,dnpl,dnpr,dshl,dshr,djrn,djrx,djln,djlx
integer:: lsc,psc,atsc,btsc,bsc,osc,nsc,hsc,basp,ptsc,csc,dsc,rsc
integer:: J2,J2u,nthr=2,max_dar,max_dbl,max_dcr,max_ddl
integer:: min_ashells,max_ashells,min_bshells,max_bshells,nspa,nspb
integer:: max_cJ2s,min_cJ2s,nmps,nep,no_levul,no_levur,nome,ll,lr
integer:: mial,maal,miar,maar,mibl,mabl,mibr,mabr
logical:: done,new,proton,neutron,obdp,ijphase=.false.
real(kind=rkw):: me,renorm,rtbas
real(kind=rkw):: phase,ninejv
real:: xa,xb,xxa,xxb,xc,xd,obtdf

contains

    subroutine Tra_TRD(base,dase)
    implicit none
    integer,intent(in)::base,dase
    integer:: no_Rac=0,OMP_get_max_threads
    real::t1,t2
    nthr=OMP_get_max_threads()

    allocate(csj2(0:no_shells))
    csj2=0
    z1=repeat(char(0),8)
    do nsc=1,no_shells
        csj2(nsc)=sj2(nsc)
    end do
    print *, ' VERSION 4.0 '
    
    read(dase+1) no_sparta,min_cJ2a,max_cJ2a,proton,neutron,xa,tpa,opa,del_cJ2a
    if (proton) typea='p'
    if (neutron) typea='n'
    allocate(nJTa(no_sparta,min_cJ2a/2:max_cJ2a/2))
    read(dase+1) nJTa(1:no_sparta,min_cJ2a/2:max_cJ2a/2)

    read(dase+2) no_spartb,min_cJ2b,max_cJ2b,proton,neutron,xb,tpb,opb,del_cJ2b
    if (proton) typeb='p'
    if (neutron) typeb='n'
    if (typea==typeb) then
        print *,' Warning nucleon types are identical '
        typea='n'
        typeb='p'
    end if
    allocate(nJTb(no_spartb,min_cJ2b/2:max_cJ2b/2))
    read(dase+2) nJTb(1:no_spartb,min_cJ2b/2:max_cJ2b/2)

    read(base+1) no_spartc,min_cJ2c,max_cJ2c,proton,neutron,xc,tpc,opc,del_cJ2c
    if (proton) typec='p'
    if (neutron) typec='n'
    allocate(nJTc(no_spartc,min_cJ2c/2:max_cJ2c/2))
    read(base+1) nJTc(1:no_spartc,min_cJ2c/2:max_cJ2c/2)

    read(base+2) no_spartd,min_cJ2d,max_cJ2d,proton,neutron,xd,tpd,opd,del_cJ2d
    if (proton) typed='p'
    if (neutron) typed='n'
    if (typec==typed) then
        print *,' Warning nucleon types are identical '
        typec='n'
        typed='p'
    end if
    allocate(nJTd(no_spartd,min_cJ2d/2:max_cJ2d/2))
    read(base+2) nJTd(1:no_spartd,min_cJ2d/2:max_cJ2d/2)
    
    call set_ab(no_sparta,min_cJ2a,max_cJ2a,del_cJ2a,no_spartb,min_cJ2b,max_cJ2b,del_cJ2b,nJTa,nJTb)
    call set_cd(no_spartc,min_cJ2c,max_cJ2c,del_cJ2c,no_spartd,min_cJ2d,max_cJ2d,del_cJ2d,nJTc,nJTd)
    
    read(base+6) trdda,max_lma,trddat,nwordsa,nspa
    read(base+3) trddb,max_lmb,trddbt,nwordsb,nspb
    min_lma=0
    min_lmb=0
    min_ashells=min_shellsa()
    max_ashells=max_shellsa()
    min_bshells=min_shellsb()
    max_bshells=max_shellsb()
    if (opa=='   i'.or.opa=='   I') then
    min_ashells=0
    max_ashells=0
    end if
    if (opb=='   i'.or.opb=='   I') then
    min_bshells=0
    max_bshells=0
    end if
    
    obtdf=real(1)
    delt=0
    if (opa=='a+a+'.or.opa=='a-a-') delt=2
    if (opa=='a+a+'.or.opb=='a-a-') delt=2
    if (opa=='a+a+'.and.opb=='a+a+') then
        delt=0
        obtdf=real(3)
    end if    
    if (opa=='a-a-'.and.opb=='a-a-') then
        obtdf=real(3)
        delt=0
    end if
    
    mial=min_ashells
    maal=max_ashells
    miar=min_ashells
    maar=max_ashells
    
    mibl=min_bshells
    mabl=max_bshells
    mibr=min_bshells
    mabr=max_bshells
    
    if (opa=='  a+'.and.opb=='  a-') then
    obtdf=real(2,rkw)/sqrt(real(3,rkw))
    delt=2
    obdp=.true.
    end if
    if (opa=='  a-'.and.opb=='  a+') then
    obtdf=real(2,rkw)/sqrt(real(3,rkw))
    delt=2
    obdp=.true.
    print *, ' WARNING: if you have not changed the default a an b allocations'
    print *, '          this will not produce normal transition densities !!!!'
    end if
    if (opa=='  a+'.and.opb=='  a+') then
    obtdf=real(2,rkw)
    delt=0
    end if
    if (opa=='  a-'.and.opb=='  a-') then
    delt=0
    obtdf=real(2,rkw)
    end if
    
    if (opa=='  a+'.and.no_cA==no_cC+1) then
        miar=0
        maar=0
    max_lma=max_sj2u
    min_lma=1
    delt=1
    end if
    if (opa=='  a-'.and.no_cA==no_cC-1) then
        mial=0
        maal=0
    max_lma=max_sj2u
    min_lma=1
    delt=1
    end if
    
    if (opb=='  a+'.and.no_cB==no_cD+1) then
        mibr=0
        mabr=0
    max_lmb=max_sj2u
    min_lmb=1
    delt=1
    end if
    if (opb=='  a-'.and.no_cB==no_cD-1) then
        mibl=0
        mabl=0
    max_lmb=max_sj2u
    min_lmb=1
    delt=1
    end if
    
5   nsc=0
    do ashl=mial,maal
    do ashr=miar,maar
    lmn=max(abs(csj2(ashl)-csj2(ashr)),min_lma)
    lmx=min(abs(csj2(ashl)+csj2(ashr)),max_lma)
    do lm=lmn,lmx,2
        do bshl=mibl,mabl
        do bshr=mibr,mabr
        lmnn=max(abs(csj2(bshl)-csj2(bshr)),min_lmb)
        lmxx=min(abs(csj2(bshl)+csj2(bshr)),max_lmb)
        do lmp=lmnn,lmxx,2
        nsc=nsc+1
        end do
        end do
        end do
    end do        
    end do        
    end do        
    nome=nsc
    
    allocate(tramp(ni,nf,nome,mindj/2:maxdj/2))
    tramp=real(0,rkw)
    allocate(meop(nome))
    allocate(meopp(nome))
    
    
    print *, ' Operator combination ',opa,' ',opb
    
    allocate(string_order(max_part_tree))
    allocate(node(max_node))
    nstring=max_part_tree
    nonode=max_node

!setup basis/

    call reset_string_order()
    no_stateL=0
    basdL=0
    do ajl=min_cJ2a,max_cJ2a,del_cJ2a
        bjln=max(abs(ajl-J2L),min_cJ2b)
        bjlx=min(abs(ajl+J2L),max_cJ2b)
        bjln=bjln+mod(bjln-min_cJ2b,del_cJ2b)
        do bjl=bjln,bjlx,2
            bm(1:1)=char(ajl)
            bm(2:2)=char(bjl)
            bm(3:8)=z1(3:8)
            call add_string(bm(1:8))
            basdL=basdL+1
        end do
    end do
    
    print *, ' No of matrices BAS L ',basdL
    allocate(basL(basdL))
    allocate(dimnvL(basdL))
    allocate(posnvL(basdL))
    dimnvL=0
    posnvL=0
    dimnvL=0
    done=.true.
    idx=0
    nsc=0
10  idx=get_next_string(done)
    nsc=nsc+1
    basL(nsc)=string_order(idx)%string(1:8)
    if (done) go to 20
    go to 10
20  if (basdL/=nsc) then
        print *, ' basdL,nsc ',basdL,nsc
        stop '  basdL /= nsc  '
    end if
    
    max_dimv=0
    no_stateL=0
    do bsc=1,basdL
    ba=basL(bsc)(1:2)
    ajl=ichar(ba(1:1))
    bjl=ichar(ba(2:2))
    dimnvL(bsc)= dna(no_sparta,ajl)*dnb(no_spartb,bjl)
    posnvL(bsc)=no_stateL
    max_dimv=max(max_dimv,dimnvL(bsc))
    no_stateL=no_stateL+dimnvL(bsc)                
    end do
    
    print *, ' No of basis states L ',no_stateL    
    
    call reset_string_order()
    no_stateR=0
    basdR=0
    do cjl=min_cJ2c,max_cJ2c,del_cJ2c
        djln=max(abs(cjl-J2R),min_cJ2d)
        djlx=min(abs(cjl+J2R),max_cJ2d)
        djln=djln+mod(djln-min_cJ2d,del_cJ2d)
        do djl=djln,djlx,del_cJ2d
            bm(1:1)=char(cjl)
            bm(2:2)=char(djl)
            bm(3:8)=z1(3:8)
            call add_string(bm(1:8))
            basdR=basdR+1
        end do
    end do
    
    print *, ' No of matrices BAS R ',basdR
    allocate(basR(basdR))
    allocate(dimnvR(basdR))
    allocate(posnvR(basdR))
    dimnvR=0
    posnvR=0
    dimnvR=0
    done=.true.
    idx=0
    nsc=0
30  idx=get_next_string(done)
    nsc=nsc+1
    basR(nsc)=string_order(idx)%string(1:8)
    if (done) go to 40
    go to 30
40  if (basdR/=nsc) then
        print *, ' basdR,nsc ',basdR,nsc
        stop '  basdR /= nsc  '
    end if
    
    no_stateR=0
    do bsc=1,basdR
    ba=basR(bsc)(1:2)
    cjl=ichar(ba(1:1))
    djl=ichar(ba(2:2))
    dimnvR(bsc)= dnc(no_spartc,cjl)*dnd(no_spartd,djl)
    posnvR(bsc)=no_stateR
    max_dimv=max(max_dimv,dimnvR(bsc))
    no_stateR=no_stateR+dimnvR(bsc)                
    end do

    print *, ' No of basis states R ',no_stateR    
    
    call reset_string_order()
! end setup basis 
! setup nineJ's   
    do ajl=min_cJ2a,max_cJ2a,del_cJ2a
      rn(1:1)=char(ajl)
      do bjl=min_cJ2b,max_cJ2b,del_cJ2b
        rn(2:2)=char(bjl)
        do cjr=min_cJ2c,max_cJ2c,del_cJ2c
          rn(3:3)=char(cjr)
          do djr=min_cJ2d,max_cJ2d,del_cJ2d
            rn(4:4)=char(djr)
            do lma=min_lma,max_lma,2
              rn(5:5)=char(lma)
                 do lmb=min_lmb,max_lmb,2
                   rn(6:6)=char(lmb)
                   lmn=mindj
                   lmx=maxdj
                   do lm=lmn,lmx,2
                        rn(7:7)=char(lm)
                        ninejv=ninej(cjr,lma,ajl,djr,lmb,bjl,J2R,lm,J2L)
                        if (abs(ninejv)>0.00000) then
                           rn(8:8)=char(0)
                           call add_string(rn(1:8),new)
                           if (new) nined=nined+1
                        end if
                     end do          
                 end do          
            end do
          end do
        end do
      end do
    end do
    if (nined/=0) then
    allocate(nine(nined))
    allocate(ninejs(nined))
    ninejs=real(0,rkw)
    done=.true.
    idx=0
    nsc=0
50  idx=get_next_string(done)
    nsc=nsc+1
    nine(nsc)=string_order(idx)%string(1:8)
    if (done) go to 60
    go to 50
60  if (nined/=nsc) then
        print *, ' nined,nsc ',nined,nsc
        stop '  nined /= nsc  '
    end if
    end if
    do nsc=1,nined
        rn=nine(nsc)(1:8)
        ajl=ichar(rn(1:1))
        bjl=ichar(rn(2:2))
        cjr=ichar(rn(3:3))
        djr=ichar(rn(4:4))
        lma=ichar(rn(5:5))
        lmb=ichar(rn(6:6))
        lm=ichar(rn(7:7))
        ninejv=ninej(cjr,lma,ajl,djr,lmb,bjl,J2R,lm,J2L)*&
                  sqrt(real((J2L+1)*(lma+1)*(lmb+1)*(J2R+1),rkw))
        if (no_cA/=no_cC.or.no_cB/=no_cD) then
            if (btest((ajl+bjl+lma+cjr+djr+lmb+J2L+J2R+lm),0))&
                         ninejv=-ninejv
            if (btest((J2R-J2L+lma-lmb)/2,0).and.ijphase)& 
                         ninejv=-ninejv
        end if
        ninejs(nsc)=ninejv*obtdf
        no_Rac=no_Rac+1
    end do
! end setup racahs


    deallocate(string_order)
    deallocate(node)

! read trda
!    nwordsa=((no_sparta-1)/32+1)
    allocate(trda(trdda))
    allocate(posa(trdda))
    allocate(dima(trdda))
    allocate(nplra(trdda))
    posa=0
    dima=0
    allocate(trdav(trddat))
    if (trdda/=0) then
    read(  base+6) trda
    read(base+6) dima
    nsc=0
    do atsc=1,trdda
    posa(atsc)=nsc
    if (dima(atsc)/=0) read(base+6) trdav(nsc+1:nsc+dima(atsc))
    do lsc=nsc+1,nsc+dima(atsc)
        if (abs(trdav(lsc))<.0000001) trdav(lsc)=0.0
    end do
    call allocatela(nplra(atsc),no_sparta,no_spartc)
    if (dima(atsc)/=0) read(base+6) (nplra(atsc)%lvect(anpl)%lword(1:nwordsa),anpl=1,nspa)
    nsc=nsc+dima(atsc)
    end do
    end if
! end read trda 
    deallocate(dima)

! read trdb
!    nwordsb=((no_spartb-1)/32+1)
    allocate(trdb(trddb))
    allocate(posb(trddb))
    allocate(dimb(trddb))
    allocate(nplrb(trddb))
    posb=0
    dimb=0
    allocate(trdbv(trddbt))
    if (trddb/=0) then
    read(base+3) trdb
    read(base+3) dimb
    nsc=0
    do btsc=1,trddb
    posb(btsc)=nsc
    if (dimb(btsc)/=0) read(base+3) trdbv(nsc+1:nsc+dimb(btsc))
    do lsc=nsc+1,nsc+dimb(btsc)
        if (abs(trdbv(lsc))<.0000001) trdbv(lsc)=0.0
    end do
    call allocatela(nplrb(btsc),no_spartb,no_spartd)
    if (dimb(btsc)/=0) read(base+3) (nplrb(btsc)%lvect(bnpl)%lword(1:nwordsb),bnpl=1,nspb)
    nsc=nsc+dimb(btsc)
    end do
    end if
! end read trdb 
    deallocate(dimb)
! read prdb  

! get max_lm
    
    maxda=0 
    do ajl=min_cJ2a,max_cJ2a,del_cJ2a
    maxda=max(maxda,dna(no_sparta,ajl))
    end do
    maxdb=0 
    do bjl=min_cJ2b,max_cJ2b,del_cJ2b
    maxdb=max(maxdb,dnb(no_spartb,bjl))
    end do
    max_dimv=maxda*maxdb
    maxdc=0 
    do cjl=min_cJ2c,max_cJ2c,del_cJ2c
    maxdc=max(maxdc,dnc(no_spartc,cjl))
    end do
    maxdd=0 
    do djl=min_cJ2d,max_cJ2d,del_cJ2d
    maxdd=max(maxdd,dnd(no_spartd,djl))
    end do
    max_dimv=max(max_dimv,maxdc*maxdd)
    
    
    allocate(BUFFER1(max_dimv,nthr))
    allocate(BUFFER2(max_dimv,nthr))
    allocate(BUFFER3(max_dimv,nthr))
    BUFFER1=real(0,rkw)
    BUFFER2=real(0,rkw)
    BUFFER3=real(0,rkw)
    allocate(coefL(no_stateL))
    allocate(coefR(no_stateR))
    allocate(coefT(no_stateL))
    allocate(coefX(no_stateL,nthr))
    
    call read_vecT(0)                
    call read_vecR(0)
    allocate(rvalue(ni))
    allocate(lvalue(nf))
                    
    do lr=1,ni 
    call read_vecR(lr)       
    do ll=1,nf
    call read_vecT(ll)        
    do delj=mindj,maxdj,2
    nome=0
    do ashl=mial,maal
    do ashr=miar,maar
    lmn=max(abs(csj2(ashl)-csj2(ashr)),min_lma)
    lmx=min(abs(csj2(ashl)+csj2(ashr)),max_lma)
    do lma=lmn,lmx,2
        do bshl=mibl,mabl
        do bshr=mibr,mabr
        lmnn=max(abs(csj2(bshl)-csj2(bshr)),min_lmb)
        lmxx=min(abs(csj2(bshl)+csj2(bshr)),max_lmb)
        do lmb=lmnn,lmxx,2
        nome=nome+1
        if (lr==1.and.ll==1.and.delj==mindj) then
        meop(nome)(1:1)=char(ashl)
        meop(nome)(2:2)=char(ashr)
        meop(nome)(3:3)=char(bshl)
        meop(nome)(4:4)=char(bshr)
        meop(nome)(5:5)=char(lma)
        meop(nome)(6:6)=char(lmb)
        meop(nome)(7:7)=char(J2R)
        meop(nome)(8:8)=char(J2L)
        end if
            if (lmb>=abs(lma-delj).and.lmb<=lma+delj) then
            call drdbh(delj,ashl,ashr,bshl,bshr,lma,lmb)
            call multiplyx()
            tramp(lr,ll,nome,delj/2)=tramp(lr,ll,nome,delj/2)+dot_product(coefL,coefT)
            end if
        end do           
        end do           
        end do
   end do                
   end do                
   end do 
   end do               
   end do               
   end do  
   
                          
   end subroutine Tra_TRD
   
   pure function dna(n,ja)
   implicit none
   integer,intent(in):: n,ja
   integer:: dna
   if (n==0) then
        dna=0
        return
   end if
   dna=sum(nJTa(1:n,ja/2))
   end function dna
   
   pure function dnb(n,jb)
   implicit none
   integer,intent(in):: n,jb
   integer:: dnb
   if (n==0) then
        dnb=0
        return
   end if
   dnb=sum(nJTb(1:n,jb/2))
   end function dnb
   
   
   pure function dnc(n,jc)
   implicit none
   integer,intent(in):: n,jc
   integer:: dnc
   if (n==0) then
        dnc=0
        return
   end if
   dnc=sum(nJTc(1:n,jc/2))
   end function dnc
   
   pure function dnd(n,jd)
   implicit none
   integer,intent(in):: n,jd
   integer:: dnd
   if (n==0) then
        dnd=0
        return
   end if
   dnd=sum(nJTd(1:n,jd/2))
   end function dnd
  
   subroutine multiplyx()
    implicit none
    integer:: thrn,OMP_get_thread_num,loop,lp
    
    coefX=real(0,rkw)  
    loop=basdL*basdR
    thrn=1
!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(basL,basR,tablebh,coefR,coefX,basdL,&
!$OMP trdav,posa,nJTa,nJTb,nplra,nplrb,BUFFER1,BUFFER2,posnvL,basdR,&
!$OMP no_sparta,no_spartb,dimb,posnvR,trdbv,posb,loop,&
!$OMP BUFFER3,tablebhn,no_spartc,no_spartd,&
!$OMP nJTc,nJTd)
!$OMP DO SCHEDULE(DYNAMIC,1)
    do bsc=1,basdL
!    bsc=(lp-1)/basdR+1
    thrn=OMP_get_thread_num()+1
    ajl=ichar(basL(bsc)(1:1))
    bjl=ichar(basL(bsc)(2:2))
!   hsc=lp-(bsc-1)*basdR
    do hsc=1,basdR
       if (tablebhn(bsc,hsc)==0) go to 75
       cjr=ichar(basR(hsc)(1:1))
       djr=ichar(basR(hsc)(2:2))
       do nsc=1,tablebhn(bsc,hsc)
            atsc=tablebh(bsc,hsc)%a(nsc)
            do psc=1,tablebh(bsc,hsc)%bme(nsc)%bn
                btsc=tablebh(bsc,hsc)%bme(nsc)%b(psc)
                ninejv=tablebh(bsc,hsc)%bme(nsc)%me(psc)
                call double_product(coefX(:,thrn),coefR,posnvL(bsc),posnvR(hsc),&
                       posa(atsc),posb(btsc),trdav,trdbv,no_sparta,no_spartc,&
                       no_spartb,no_spartd,nJTa(:,ajl/2),nJTc(:,cjr/2),nJTb(:,bjl/2),&
                       nJTd(:,djr/2),ninejv,nplra(atsc),nplrb(btsc),&
                       BUFFER1(:,thrn),BUFFER2(:,thrn),BUFFER3(:,thrn),&
                       ajl,bjl,cjr,djr)
            end do            
50     end do 
75  end do
    end do
    
!$OMP END PARALLEL
    coefL(:)=coefX(:,1)
    do lp=2,nthr
    coefL(:)=coefL(:)+coefX(:,lp)
    end do
    
   end subroutine multiplyx
      
    subroutine drdbh(delj,ashl,ashr,bshl,bshr,lma,lmb)
    implicit none
    integer,intent(in):: delj,ashl,ashr,bshl,bshr,lma,lmb
    integer:: store
    integer,save:: njc2,xjc2,bn=1,tn=1
    logical,save:: first=.true.
    bn=max(basdL,basdR)
    tn=bn
    if (first) then    
    allocate(tablebh(basdL,basdR))
    allocate(tablebhn(basdL,basdR))
    tablebhn=0
    else
        do  bsc=1,basdL
        do  hsc=1,basdR
           do nsc=1,tn
           tablebh(bsc,hsc)%bme(nsc)%me=real(0,rkw)
           tablebh(bsc,hsc)%bme(nsc)%b=0
           tablebh(bsc,hsc)%bme(nsc)%h=0
           tablebh(bsc,hsc)%bme(nsc)%bn=0
           end do
        end do
        end do
        do  bsc=1,basdL
        do  hsc=1,basdR
           tablebh(bsc,hsc)%a=0
        end do
        end do
        tablebhn=0
    end if
    
    if (first) then
        do  bsc=1,basdL
        do  hsc=1,basdR
           allocate(tablebh(bsc,hsc)%a(tn))
           allocate(tablebh(bsc,hsc)%bme(tn))
           do nsc=1,tn
           tablebh(bsc,hsc)%bme(nsc)%bn=0
           end do
        end do
        end do
    
        do  bsc=1,basdL
        do  hsc=1,basdR
           do nsc=1,tn
           allocate(tablebh(bsc,hsc)%bme(nsc)%me(bn))
           allocate(tablebh(bsc,hsc)%bme(nsc)%b(bn))
           allocate(tablebh(bsc,hsc)%bme(nsc)%h(bn))
           tablebh(bsc,hsc)%bme(nsc)%me=real(0,rkw)
           tablebh(bsc,hsc)%bme(nsc)%h=0
           tablebh(bsc,hsc)%bme(nsc)%b=0
           end do
        end do
        end do
        first=.false.
    end if


    if (abs(lma-lmb)>delj.or.lmb+lma<delj) return
    
    do bsc=1,basdL
    ajl=ichar(basL(bsc)(1:1))
    bjl=ichar(basL(bsc)(2:2))
    do hsc=1,basdR
    cjr=ichar(basR(hsc)(1:1))
    djr=ichar(basR(hsc)(2:2))
    
        if (cjr>=abs(ajl-lma).and.cjr<=ajl+lma) then   
        if (djr>=abs(bjl-lmb).and.djr<=bjl+lmb) then   
        ta(isl:isl)=char(ashl)
        ta(isr:isr)=char(ashr)
        ta(ijl:ijl)=char(ajl)
        ta(ijr:ijr)=char(cjr)
        ta(ilm:ilm)=char(lma)
        ta(6:8)=z1(6:8)
        atsc=tfind(trda,trdda,ta)
        if (atsc==0) go to 100
        
        tablebhn(bsc,hsc)=tablebhn(bsc,hsc)+1
        tablebh(bsc,hsc)%a(tablebhn(bsc,hsc))=atsc
        
        rn(1:2)=basL(bsc)(1:2)
        rn(3:4)=basR(hsc)(1:2)
        rn(5:5)=char(lma)
        rn(6:6)=char(lmb)
        rn(7:7)=char(delj)
        rsc=tfind(nine,nined,rn)
        if (rsc==0) go to 100
        
            tb(isl:isl)=char(bshl)
            tb(isr:isr)=char(bshr)
            tb(ijl:ijl)=char(bjl)
            tb(ijr:ijr)=char(djr)
            tb(ilm:ilm)=char(lmb)
            tb(6:8)=z1(6:8)
            btsc=tfind(trdb,trddb,tb)
            if (btsc==0) go to 100
            
            tablebh(bsc,hsc)%bme(tablebhn(bsc,hsc))%bn=&
                         tablebh(bsc,hsc)%bme(tablebhn(bsc,hsc))%bn+1
            
                bn=tablebh(bsc,hsc)%bme(tablebhn(bsc,hsc))%bn
                tablebh(bsc,hsc)%bme(tablebhn(bsc,hsc))%b(bn)=btsc
                if (.not.ijphase) &
                    tablebh(bsc,hsc)%bme(tablebhn(bsc,hsc))%me(bn)=ninejs(rsc)
        endif
        endif
100 end do
    end do
   
   
   end subroutine drdbh 

    subroutine read_vecT(lsc) 
       implicit none
       integer,intent(in):: lsc
       integer:: no_state
       logical:: xst
       
       if (lsc==0) then
       inquire(file=NucleusL//file_code(J2L/2,no_cA,no_cB)//'.xvc',exist=xst)
       if (xst)then
       open(unit=703,file=NucleusL//file_code(J2L/2,no_cA,no_cB)//'.xvc',form='unformatted',action='read')
       else
       inquire(file=NucleusL//file_code(J2L/2,no_cB,no_cA)//'.xvc',exist=xst)
       if (xst) then
       print *, ' Warning ',NucleusL//file_code(J2L/2,no_cA,no_cB)//'.xvc not found'
       print *, ' Using   ',NucleusL//file_code(J2L/2,no_cB,no_cA)//'.xvc instead'
       open(unit=703,file=NucleusL//file_code(J2L/2,no_cB,no_cA)//'.xvc',form='unformatted',action='read')
       else
       stop '  Final vector file not found.'
       end if
       end if
       read(703) no_state,no_levelL
       if (no_state/=no_stateL) then
        print *, no_state, no_stateL
        print *, no_sparta,no_spartb
        stop '  error no_stateL '
       end if
       no_levul=min(nf,no_levelL)
       nf=no_levul
       return
       end if
       if (lsc==1) then
       rewind(703)
       read(703)
       end if
       read(703) lvalue(lsc),coefT
!       if (dot_product(coefT,coefT)<0.95)  stop ' coefT norm error'
    end subroutine read_vecT   
     
    subroutine read_vecR(rsc) 
       implicit none
       integer,intent(in):: rsc
       integer:: no_state
       logical:: xst
       
       if (rsc==0) then
       inquire(file=NucleusR//file_code(J2R/2,no_cC,no_cD)//'.xvc',exist=xst)
       if (xst)then
       open(unit=704,file=NucleusR//file_code(J2R/2,no_cC,no_cD)//'.xvc',form='unformatted',action='read')
       else
       inquire(file=NucleusR//file_code(J2R/2,no_cD,no_cC)//'.xvc',exist=xst)
       if (xst) then
       print *, ' Warning ',NucleusR//file_code(J2R/2,no_cC,no_cD)//'.xvc not found'
       print *, ' Using   ',NucleusR//file_code(J2R/2,no_cD,no_cC)//'.xvc instead'
       open(unit=704,file=NucleusR//file_code(J2R/2,no_cD,no_cC)//'.xvc',form='unformatted',action='read')
       else
       stop '  Initial vector file not found.'
       end if
       end if
       
       read(704) no_state,no_levelR
       if (no_state/=no_stateR) then
        print *, no_state,no_stateR
        print *, no_spartc,no_spartd
        stop '  error no_stateR '
       end if 
       no_levur=min(ni,no_levelR)
       ni=no_levur
       return
       end if
       if (rsc==1) then
       rewind(704)
       read(704)
       end if
       read(704) rvalue(rsc),coefR
!       if (dot_product(coefR,coefR)<0.95)  stop ' coefR norm error'
    end subroutine read_vecR
    
    subroutine reorder()
    implicit none
    
    if (opa=='a+a-'.and.(opb=='   I'.or.opb=='   i')) then
        do nsc=1,nome
            meopp(nsc)(1:1)=meop(nsc)(8:8)
            meopp(nsc)(2:2)=char(0)
            meop(nsc)(2:3)=meop(nsc)(1:2)
            meop(nsc)(1:1)=char(0)
            meop(nsc)(5:5)=char(sj2(ichar(meop(nsc)(2:2))))
            meop(nsc)(7:7)=char(sj2(ichar(meop(nsc)(3:3))))
            meop(nsc)(6:6)=char(1)
            meop(nsc)(8:8)=char(1)
        end do
    end if
    if (opb=='a+a-'.and.(opa=='   I'.or.opa=='   i')) then
        do nsc=1,nome
            meopp(nsc)(1:1)=meop(nsc)(8:8)
            meopp(nsc)(2:2)=char(0)
            meop(nsc)(2:3)=meop(nsc)(3:4)
            meop(nsc)(4:4)=char(0)
            meop(nsc)(5:5)=char(sj2(ichar(meop(nsc)(2:2))))
            meop(nsc)(7:7)=char(sj2(ichar(meop(nsc)(3:3))))
            meop(nsc)(6:6)=char(1)
            meop(nsc)(8:8)=char(1)
        end do
    end if
    if (opa=='  a+'.and.opb=='  a-') then
        do nsc=1,nome
            meopp(nsc)(1:1)=meop(nsc)(8:8)
            meopp(nsc)(2:2)=char(2)
            meop(nsc)(2:2)=meop(nsc)(1:1)
            meop(nsc)(3:3)=meop(nsc)(4:4)
            meop(nsc)(1:1)=char(0)
            meop(nsc)(4:4)=char(0)
            meop(nsc)(5:5)=char(sj2(ichar(meop(nsc)(2:2))))
            meop(nsc)(7:7)=char(sj2(ichar(meop(nsc)(3:3))))
            meop(nsc)(6:6)=char(1)
            meop(nsc)(8:8)=char(1)
        end do
    end if
    if (opb=='  a+'.and.opa=='  a-') then
        do nsc=1,nome
            meopp(nsc)(1:1)=meop(nsc)(8:8)
            meopp(nsc)(2:2)=char(2)
            meop(nsc)(1:1)=char(0)
            meop(nsc)(4:4)=meop(nsc)(2:2)
            meop(nsc)(2:2)=meop(nsc)(3:3)
            meop(nsc)(3:3)=meop(nsc)(4:4)
            meop(nsc)(4:4)=char(0)
            meop(nsc)(5:5)=char(sj2(ichar(meop(nsc)(3:3))))
            meop(nsc)(7:7)=char(sj2(ichar(meop(nsc)(2:2))))
            meop(nsc)(6:6)=char(1)
            meop(nsc)(8:8)=char(1)
        end do
    end if
    if (opa=='  a+'.and.opb=='  a+') then
        do nsc=1,nome
            meopp(nsc)(1:1)=meop(nsc)(8:8)
            meop(nsc)(2:2)=meop(nsc)(3:3)
            meop(nsc)(3:3)=char(0)
            meop(nsc)(4:4)=char(0)
            meop(nsc)(5:5)=meop(nsc)(8:8)
            meop(nsc)(7:7)=char(0)
            meop(nsc)(6:6)=char(0)
            meop(nsc)(8:8)=char(0)
        end do
    end if
    if (opa=='  a-'.and.opb=='  a-') then
        do nsc=1,nome
            meopp(nsc)(1:1)=meop(nsc)(8:8)
            meop(nsc)(3:3)=meop(nsc)(2:2)
            meopp(nsc)(2:2)=char(0)
            meopp(nsc)(1:1)=char(0)
            meop(nsc)(7:7)=meop(nsc)(8:8)
            meop(nsc)(5:5)=char(0)
            meop(nsc)(6:6)=char(0)
            meop(nsc)(8:8)=char(0)
        end do
    end if
    if (opb=='  a+'.and.(opa=='   I'.or.opa=='   i')) then
        do nsc=1,nome
            meopp(nsc)(1:1)=meop(nsc)(8:8)
            meopp(nsc)(2:2)=char(1)
            meop(nsc)(2:2)=meop(nsc)(3:3)
            meop(nsc)(3:3)=char(0)
            meop(nsc)(7:7)=char(0)
            meop(nsc)(5:5)=char(sj2(ichar(meop(nsc)(2:2))))
            meop(nsc)(6:6)=char(1)
            meop(nsc)(8:8)=char(0)
        end do
    end if
    if (opa=='  a+'.and.(opb=='   I'.or.opb=='   i')) then
        do nsc=1,nome
            meopp(nsc)(1:1)=meop(nsc)(8:8)
            meopp(nsc)(2:2)=char(0)
            meop(nsc)(2:2)=meop(nsc)(1:1)
            meop(nsc)(1:1)=char(0)
            meop(nsc)(5:5)=char(sj2(ichar(meop(nsc)(2:2))))
            meop(nsc)(7:7)=char(0)
            meop(nsc)(6:6)=char(1)
            meop(nsc)(8:8)=char(0)
        end do
    end if
    if (opb=='  a-'.and.(opa=='   I'.or.opa=='   i')) then
        do nsc=1,nome
            meopp(nsc)(1:1)=meop(nsc)(8:8)
            meopp(nsc)(2:2)=char(1)
            meop(nsc)(3:3)=meop(nsc)(4:4)
            meop(nsc)(4:4)=char(0)
            meop(nsc)(5:5)=char(0)
            meop(nsc)(7:7)=char(sj2(ichar(meop(nsc)(3:3))))
            meop(nsc)(6:6)=char(0)
            meop(nsc)(8:8)=char(1)
        end do
    end if
    if (opa=='  a-'.and.(opb=='   I'.or.opb=='   i')) then
        do nsc=1,nome
            meopp(nsc)(1:1)=meop(nsc)(8:8)
            meop(nsc)(3:3)=meop(nsc)(2:2)
            meopp(nsc)(2:2)=char(0)
            meop(nsc)(1:1)=char(0)
            meop(nsc)(7:7)=char(sj2(ichar(meop(nsc)(3:3))))
            meop(nsc)(7:7)=char(0)
            meop(nsc)(6:6)=char(0)
            meop(nsc)(8:8)=char(1)
        end do
    end if
    if (opa=='a+a+'.and.(opb=='   I'.or.opb=='   i')) then
        do nsc=1,nome
            meopp(nsc)(1:1)=meop(nsc)(8:8)
            meopp(nsc)(2:2)=char(2)
            meop(nsc)(5:5)=meop(nsc)(5:5)
            meop(nsc)(6:6)=char(2)
            meop(nsc)(8:8)=char(0)
        end do
    end if
    if (opb=='a+a+'.and.(opa=='   I'.or.opa=='   i')) then
        do nsc=1,nome
            meopp(nsc)(1:1)=meop(nsc)(8:8)
            meopp(nsc)(2:2)=char(2)
            meop(nsc)(1:1)=meop(nsc)(3:3)
            meop(nsc)(2:2)=meop(nsc)(4:4)
            meop(nsc)(3:3)=char(0)
            meop(nsc)(4:4)=char(0)
            meop(nsc)(5:5)=meop(nsc)(6:6)
            meop(nsc)(6:6)=char(2)
            meop(nsc)(8:8)=char(0)
        end do
    end if
    if (opa=='a+a+'.and.opb=='a-a-') then
        do nsc=1,nome
            meopp(nsc)(1:1)=meop(nsc)(8:8)
            meopp(nsc)(2:2)=char(0)
            meop(nsc)(7:7)=meop(nsc)(6:6)
            meop(nsc)(6:6)=char(2)
            meop(nsc)(8:8)=char(2)
        end do
    end if
    if (opa=='a+a+'.and.opb=='a-a-') then
        do nsc=1,nome
            meopp(nsc)(1:1)=meop(nsc)(8:8)
            meopp(nsc)(2:2)=char(0)
            meop(nsc)(7:7)=meop(nsc)(6:6)
            meop(nsc)(6:6)=char(2)
            meop(nsc)(8:8)=char(2)
        end do
    end if
    if (opa=='a+a+'.and.opb=='a+a+') then
        do nsc=1,nome
            meopp(nsc)(1:1)=meop(nsc)(8:8)
            meopp(nsc)(2:2)=char(0)
            meop(nsc)(7:7)=meop(nsc)(6:6)
            meop(nsc)(6:6)=char(2)
            meop(nsc)(8:8)=char(2)
        end do
    end if
    if (opa=='a-a-'.and.opb=='a-a-') then
        do nsc=1,nome
            meopp(nsc)(1:1)=meop(nsc)(8:8)
            meopp(nsc)(2:2)=char(0)
            meop(nsc)(7:7)=meop(nsc)(6:6)
            meop(nsc)(6:6)=char(2)
            meop(nsc)(8:8)=char(2)
        end do
    end if
    if (opa=='a+a-'.and.opb=='a+a-') then
        do nsc=1,nome
            meopp(nsc)(1:1)=meop(nsc)(8:8)
            meopp(nsc)(2:2)=char(0)
            meop(nsc)(7:7)=meop(nsc)(6:6)
            meop(nsc)(6:6)=char(0)
            meop(nsc)(8:8)=char(0)
        end do
    end if
    if (opa=='  a+'.and.opb=='a+a-') then
        do nsc=1,nome
            meopp(nsc)(1:1)=meop(nsc)(8:8)
            meopp(nsc)(2:2)=char(1)
            meop(nsc)(7:7)=meop(nsc)(6:6)
            meop(nsc)(6:6)=char(1)
            meop(nsc)(8:8)=char(0)
        end do
    end if
    if (opb=='  a+'.and.opa=='a+a-') then
        do nsc=1,nome
            meopp(nsc)(1:1)=meop(nsc)(8:8)
            meopp(nsc)(2:2)=char(1)
            meop(nsc)(7:7)=meop(nsc)(6:6)
            meop(nsc)(6:6)=char(0)
            meop(nsc)(8:8)=char(1)
        end do
    end if
    if (opa=='  a+'.and.opb=='a+a+') then
        do nsc=1,nome
            meopp(nsc)(1:1)=meop(nsc)(8:8)
            meopp(nsc)(2:2)=char(1)
            meop(nsc)(7:7)=meop(nsc)(6:6)
            meop(nsc)(6:6)=char(1)
            meop(nsc)(8:8)=char(2)
        end do
    end if
    if (opb=='  a+'.and.opa=='a+a+') then
        do nsc=1,nome
            meopp(nsc)(1:1)=meop(nsc)(8:8)
            meopp(nsc)(2:2)=char(1)
            meop(nsc)(7:7)=meop(nsc)(6:6)
            meop(nsc)(6:6)=char(2)
            meop(nsc)(8:8)=char(1)
        end do
    end if
    if (opa=='  a-'.and.opb=='a+a-') then
        do nsc=1,nome
            meopp(nsc)(1:1)=meop(nsc)(8:8)
            meopp(nsc)(2:2)=char(1)
            meop(nsc)(7:7)=meop(nsc)(6:6)
            meop(nsc)(6:6)=char(1)
            meop(nsc)(8:8)=char(0)
        end do
    end if
    if (opb=='  a-'.and.opa=='a+a-') then
        do nsc=1,nome
            meopp(nsc)(1:1)=meop(nsc)(8:8)
            meopp(nsc)(2:2)=char(1)
            meop(nsc)(7:7)=meop(nsc)(6:6)
            meop(nsc)(6:6)=char(0)
            meop(nsc)(8:8)=char(1)
        end do
    end if
    if (opa=='  a-'.and.opb=='a-a-') then
        do nsc=1,nome
            meopp(nsc)(1:1)=meop(nsc)(8:8)
            meopp(nsc)(2:2)=char(1)
            meop(nsc)(7:7)=meop(nsc)(6:6)
            meop(nsc)(6:6)=char(1)
            meop(nsc)(8:8)=char(2)
        end do
    end if
    if (opb=='  a-'.and.opa=='a-a-') then
        do nsc=1,nome
            meopp(nsc)(1:1)=meop(nsc)(8:8)
            meopp(nsc)(2:2)=char(1)
            meop(nsc)(7:7)=meop(nsc)(6:6)
            meop(nsc)(6:6)=char(2)
            meop(nsc)(8:8)=char(1)
        end do
    end if
end subroutine reorder
    
    
        
end module Tra  



!  NuTra.f90 
!
!  FUNCTIONS:
!  NuTra      - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: NuTra
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program NuTra
    
    use Tra
    use LanczosParameters
    use ProjParameters
    use InputNuShell
    use OutputNuShell
    use OrderParameters
    use MvecParameters
    use Mscheme
    use Files
    use Extra
    use TrdChar
    implicit none

    ! Variables
    
    real:: t1,t2,et1,et2,ms
    integer:: hr,mn,sc
    integer:: indx,err,nis,niv,llbab,ifbab
    integer:: j,j1,twrite,twrite1,k=0,nmp1,nmp2 ,delmj   
    character(len=6):: inExt1,inExt2,Nucleus1,Nucleus2
    character(len=10):: btime
    character(len=8):: bdate
    character(len=1):: xtra1,xtra2
    integer:: max_cJ21,max_cJ22,min_cJ21,min_cJ22,max_cT21,max_cT22,min_cT21,min_cT22
    integer:: status,res
    character(len=14)::buffer
    real::  ep,en,glp,gln,gsp,gsn,mass,charge,rdial
    integer,dimension(:),allocatable:: meindx
    
    call output_header(6)
    call input_nucleus(inExt1)

    call output_welcome(6)    
    call input(inExt1)

    call output_shells(6)    
    call check_shell_order()
    call setup_shells()

    call setup_mscheme()
    call setup_file_ext()
    call input_intfile(opfile1,4)
    call input_intfile(opfile2,5)
    call input_intfile(mbmemf,7)
    call input_intfile(inExt2,8)
    call output_end_setup(6)        
    call cpu_time(t1)
    call date_and_time(bdate,btime)
    open(unit=801,file=opfile1//'.tri',form='unformatted',action='read')
    open(unit=806,file=opfile1//'.trd',form='unformatted',action='read')
    open(unit=802,file=opfile2//'.tri',form='unformatted',action='read')
    open(unit=803,file=opfile2//'.trd',form='unformatted',action='read')
    open(unit=901,file=opfile1//'.trf',form='unformatted',action='read')
    open(unit=902,file=opfile2//'.trf',form='unformatted',action='read')
        open(unit=10, file=mbmemf//'.me')
        call read_comments()
        read(10,*) ep,en,glp,gln,gsp,gsn,mass,charge
        read(10,*) mindj,maxdj,mindt,maxdt,mindmj,maxdmj,mindmt,maxdmt,ni,nf,rdial
        close(10)
        open(unit=12,file=mbmemf//'.meo',action='write')
        write(12,'(8f9.4)') ep,en,glp,gln,gsp,gsn,mass,charge
        write(12,'(10i5,f9.4)') mindj,maxdj,mindt,maxdt,mindmj,maxdmj,mindmt,maxdmt,ni,nf,rdial
    call output_end_setup(6)
    call cpu_time(t1)
    call date_and_time(bdate,btime)

    J2L=min_cJ2
    min_cJ22=min_cJ2
    nmp2=nmp
    NucleusL=Nucleus
    Nucleus2=Nucleus
    j=min_cJ2
    read(901) no_sparta,min_cJ2a,max_cJ2a,proton,neutron,xa
    if (proton) typea='p'
    if (neutron) typea='n'
    rewind(901)
    read(902) no_spartb,min_cJ2b,max_cJ2b,proton,neutron,xb
    if (proton) typeb='p'
    if (neutron) typeb='n'
    if (typea==typeb) then
        print *,' Warning nucleon types are identical '
        typea='n'
        typeb='p'
    end if
    rewind(902)
! setup nucleon numbers         
    nmps=nmp
    if (typea=='p') then
    nmp=(no_cN+nmps)/2
    no_cB=nmp
    nmp=(no_cN-nmps)/2
    no_cA=nmp
    else
    nmp=(no_cN-nmps)/2
    no_cB=nmp
    nmp=(no_cN+nmps)/2
    no_cA=nmp
    end if
! end set up nucleon numbers

    call input(inExt2)
    
    min_cJ21=min_cj2
    nmp1=nmp
    j1=min_cJ2
    J2R=min_cJ2
    NucleusR=Nucleus
    Nucleus1=Nucleus
    read(801) no_spartc,min_cJ2c,max_cJ2c,proton,neutron,xc
    if (proton) typec='p'
    if (neutron) typec='n'
    rewind(801)
    read(802) no_spartd,min_cJ2d,max_cJ2d,proton,neutron,xd
    if (proton) typed='p'
    if (neutron) typed='n'
    if (typec==typed) then
        print *,' Warning nucleon types are identical '
        typec='n'
        typed='p'
    end if
    rewind(802)
! setup nucleon numbers    
    nmps=nmp
    if (typec=='p') then
    nmp=(no_cN+nmps)/2
    no_cD=nmp
    nmp=(no_cN-nmps)/2
    no_cC=nmp
    else
    nmp=(no_cN-nmps)/2
    no_cD=nmp
    nmp=(no_cN+nmps)/2
    no_cC=nmp
    end if
! end set up nucleon numbers
    delmj=-j1+j    
    twrite=abs(nmp2)
    twrite1=abs(nmp1)
    call Tra_TRD(800,900)
    open(unit=500,&
    file=mbmemf//file_code(J2R/2,no_cC,no_cD)//file_code(J2L/2,no_cA,no_cB)//'.tra',&
                    action='write')
    open(unit=600,&
    file=mbmemf//file_code(J2R/2,no_cC,no_cD)//file_code(J2L/2,no_cA,no_cB)//'.trb',&
                    form='unformatted', action='write')
    open(unit=700,&
    file=mbmemf//file_code(J2R/2,no_cC,no_cD)//file_code(J2L/2,no_cA,no_cB)//'.DIN',&
                    form='formatted', action='write')
    write(700,*) ' Special Dens Input File for Alex Brown '
    write(700,*)              
    write(600+k,err=200) j,twrite,j1,twrite1,twrite,twrite1,ni,nf
    write(600+k,err=200) mindj,maxdj,mindt,maxdt,delmj,-twrite1+twrite
    call reorder()
    allocate(meindx(nome))
    lsc=0
    do nsc=1,nome
        if (sum(abs(tramp(:,:,nsc,:)))>0.0001) then
            write(12,'(8i5,f9.4)') ichara(meop(nsc),8),1.0
            lsc=lsc+1
            meindx(nsc)=lsc
        end if
    end do
    close(unit=12)
    do delj=mindj,maxdj,2
        if (sum(abs(tramp(:,:,:,delj/2)))>0.0001) then
            write(600+k,err=200) 1
        else 
            write(600+k,err=200) 0
            go to 50
        end if    
        call output_headerX(500+k)
        write(500+k,*,err=100) ' NuTra(mp)X Program'
        write(500+k,*,err=100)
        write(500+k,*,err=100) ' Initial Nucleus ',NucleusR
        write(600+k,err=200) NucleusR
        write(500+k,*,err=100) ' Final   Nucleus ',NucleusL
        write(600+k,err=200) NucleusL
        write(500+k,*,err=100)
        write(500+k,*,err=100) ' Initial 2*J, 2*T ',j1,twrite1
        write(600+k,err=200) j1,twrite1
        write(500+k,*,err=100) ' Final   2*J, 2*T ',j,twrite
        write(600+k,err=200) j,twrite
        write(500+k,*,err=100)
        write(500+k,*,err=100) ' Mbme Name ', mbmemf
        write(500+k,*,err=100)
        write(500+k,*,err=100) ' 2*J , 2*T for operator',delj,delt
        write(500+k,*,err=100) ' 2*M , 2*Tz for operator',0,0
        write(600+k,err=200) delj,abs(twrite1-twrite),delmj,-twrite1+twrite,.false.,obdp
        write(500+k,*,err=100)
        write(500+k,*,err=100) ' Reduced Matrix Elements'
        write(500+k,*,err=100)
        write(500+k,*,err=100) ' The states of the Initial Nucleus run horizontally.'
        write(500+k,*,err=100) ' The states of the Final   Nucleus run vertically.'
        write(500+k,*,err=100)
        write(600+k,err=200) no_levul,no_levur,lvalue(1:no_levul),rvalue(1:no_levur)
        nsc=0
        do indx=1,nome
            if (sum(abs(tramp(:,:,indx,delj/2)))>0.0001) then
            nsc=nsc+1
            write(600+k) meindx(indx),nsc
            write (500+k,'(a42,10i3)',err=100) ' sh1,sh2,sh3,sh4,lm12,t12,lm34,t34,lm,t : ',&
                      ichara(meop(indx),8),delj,ichar(meopp(indx)(2:2))
            write(600+k) int(ichara(meop(indx),8),2),delj,ichar(meopp(indx)(2:2))
            write(500+k,'(a5,10(2x,i4,1hi,1x))',err=100)'    ',(lr,lr=1,ni)
                do ll=1,nf
                     write(500+k,'(i4,a1,10f8.4)',err=100)ll,'f',tramp(1:ni,ll,indx,delj/2)
                     write(600+k,err=200) tramp(1:ni,ll,indx,delj/2)
                end do    
                do llbab=1,ni
                do ll=1,nf
                     write(700+k,'(a6,2i3,a6,2i3,10i3,2i4,1x,f8.4)',err=100) &
                      NucleusR,j1,twrite1,NucleusL,j,twrite, &
                      ichara(meop(indx),8),delj,ichar(meopp(indx)(2:2)), &
                     llbab,ll,tramp(llbab,ll,indx,delj/2)
                end do  
                end do     
                write(500+k,*,err=100)
            end if    
        end  do
        write(600+k) -1
50  end do
    
        
    read(btime,'(i2,i2,i2,f4.3)') hr,mn,sc,ms
    et1=hr*3600.+mn*60.+sc*1. +ms
    call date_and_time(bdate,btime)
    call cpu_time(t2)
    read(btime,'(i2,i2,i2,f4.3)') hr,mn,sc,ms
    et2=hr*3600.+mn*60.+sc*1. +ms
    t1=t2-t1
    et1=et2-et1
    if (et1<0.0) et1=et1+24.*3600.
    if (status==-1) call output_completed(6,'NuTrs     ') 
    call output_time('NuTra',t1,et1)    
    print *
    stop
100 print *, ' Error writing trd file'
    stop
200 print *, ' Error writing trs file'
    stop
    
    end program NuTra
    
