module Opd

use Parameters
use Partition
use ProjParameters
use Clebschg
use InputNuShell
use Finds
use TrdChar3
use SOrder3
use SubM
use XParameters

implicit none
type(trx):: trxn
integer:: max_cJ21,max_cJ22,min_cJ21,min_cJ22,max_cT21,max_cT22,min_cT21,min_cT22
integer:: del_cJ21,del_cJ22,del_cT21,del_cT22
character(len=6):: inExt1,inExt2,Nucleus1,Nucleus2
character(len=7)::opfile,intfile
character(len=6)::nsfile
character(len=12):: tm,tn
character(len=12):: ta,tb
character(len=4)::op,opin
type(projvec),dimension(:,:),allocatable:: proj_array,proj_arrayr    
integer::no_Opm,indexx,no_spart,nz1=0,no_Opd=0,min_cJ2u,ic,icp,no_spartr
integer:: maxmru,minmru
type(spartition),dimension(:),allocatable:: spart,spartr
integer,dimension(:,:),allocatable:: nJTa,nJTar
character(len=12),dimension(:),allocatable:: trda,trma
character(len=8),dimension(:),allocatable:: trdao
type(mat),dimension(:),allocatable:: trmam
type(vect),dimension(:),allocatable:: trdav
integer,dimension(:),allocatable:: trddv,csj2         
integer,dimension(:,:,:),allocatable:: posn,dimn
logical(kind=1),dimension(:,:,:),allocatable:: lsp
type(la),dimension(:),allocatable::nplr
character(len=12):: z1
integer::  trdd,trmd,trddu,max_lm,Tparity1,Tparity2,totlms,totlm,nelm,nel
logical:: proton,neutron,special=.false.,llsp
real:: x1,x2


contains

    subroutine read_proj(base,k)
    implicit none
    integer,intent(in):: base,k
    integer:: pindx,nJT,j,t,dim,ngd,np,sh
    do np=1,no_spart
        read(base+k,err=10,end=11) pindx,nJT,j,t,dim,x1
        allocate(proj_array(np,j/2)%coef(dim,nJT))
        nJTa(np,j/2)=nJT
        read(base+k,err=10,end=11) spart(np)
        if (nJT >0) then
            do ngd=1,nJT
                read(base+k,err=10,end=11) proj_array(np,j/2)%coef(1:dim,ngd)
            end do
        end if
        
    end do    
    return
10  print *, ' Error reading prj file '
    stop            
11  print *, ' EOF reading prj file '
    stop            
    
    end subroutine read_proj
    
    subroutine read_projr(base,k)
    implicit none
    integer,intent(in):: base,k
    integer:: pindx,nJT,j,t,dim,ngd,np,sh
    do np=1,no_spartr
        read(base+k,err=10,end=11) pindx,nJT,j,t,dim,x2
        allocate(proj_arrayr(np,j/2)%coef(dim,nJT))
        nJTar(np,j/2)=nJT
        read(base+k,err=10,end=11) spartr(np)
        if (nJT >0) then
            do ngd=1,nJT
                read(base+k,err=10,end=11) proj_arrayr(np,j/2)%coef(1:dim,ngd)
            end do
        end if
        
    end do    
    return
10  print *, ' Error reading prjr file '
    stop            
11  print *, ' EOF reading prjr file '
    stop            
    
    end subroutine read_projr

    subroutine Opd_Opm(baseo)
    implicit none
    integer,intent(in)::baseo
    integer:: jl,jr,npl,npr,nJL,nJR,shl,shr,lm,idx,lmx,lmn,mr,ml,phs,ipl,ipr,lmm
    integer:: bsc,rsc,dsc,lsc,nsc,jrn,jrx,sjl,sjr,mru,nmiss=0,max_sj2ud,words,icm
    real(kind=rkw)::clebsch1,clebsch2,me,c1c2,sqrt2m1,sqrt3m1
    type(spartition)::part_op
    logical:: done,new,ftypea,ftypeb,warn=.false.,ltrd,lnum,single,two,idop
    integer:: min_shells,max_shells,maxlmm=0,lmt,sh1,sh2,sh3,sh4
! open TRM file and deduce types
 
    allocate(lsp(0:max_sj2u,no_shells,no_shells))
    lsp=.false.
    allocate(csj2(0:no_shells))
    csj2=0
    do lsc=1,no_shells
        csj2(lsc)=sj2(lsc)
    end do
    sqrt2m1=real(1,rkw)/sqrt(real(2,rkw))
    sqrt3m1=real(1,rkw)/sqrt(real(3,rkw))
   
    z1=repeat(char(0),12)
    ftypea=.false.
    ftypeb=.true.
    open(unit=2,file=opfile//'.mst',action='write')
    open(unit=1,file=opfile//'.trm',form='unformatted')
    read(1) indexx,proton,neutron,op
    if (opin=='   O'.or.opin=='   o') return
    if (.not.(opin==op)) then
        print *, ' Operators not consistent. Input,file ',opin,' ',op
        stop
    end if
    idop=.false.
    two=.false.
    single=.false.
    lnum=.false.
    ltrd=.false.
    if (op=='   I'.or.op=='   i') idop=.true.
    if (op=='   N'.or.op=='   n') lnum=.true.
    if (op=='  a+'.or.op=='  a-') single=.true.
    if (op=='a+a+'.or.op=='a-a-') two=.true.
    if (op=='a+a-') ltrd=.true.
    if (proton.and.neutron) stop ' Error p.and.n '
    if (proton) print *,' Proton OPD '
    if (neutron) print *,' Neutron OPD '
    nsfile=adjustr(Nucleus1)
    if (nsfile(6:6)=='a') then
        ftypea=.true.
        ftypeb=.false.
    else if (.not.(nsfile(6:6)=='b')) then
        stop ' Error type not equal to a or b'
    end if
    if (proton.and.ftypea) then
        typea='p'
        typeb='n'
    else if(neutron.and.ftypeb) then
        typeb='n'
        typea='p'
    else if(proton.and.ftypeb) then
        typeb='p'
        typea='n'
    else if(neutron.and.ftypea) then
        typea='n'
        typeb='p'
    end if
    if (ftypea) then
        min_shells=min_shellsa()    
        max_shells=max_shellsa()
        max_sj2ud=min(max_sj2ua(),max_sj2ub())    
        max_sj2u=max_sj2ua()
    else    
        min_shells=min_shellsb()    
        max_shells=max_shellsb()    
        max_sj2ud=min(max_sj2ub(),max_sj2ub())    
        max_sj2u=max_sj2ub()
    end if
    
    if (ltrd) then
    call get_command_argument(5,intfile) 
    intfile=adjustr(intfile)
    llsp=.false.   
!    inquire(file=intfile//'.lsp',exist=llsp)
    if (llsp) then
    open(unit=43,file=intfile//'.lsp')
    read(43,*) totlms
    totlms=totlms+1
    do totlm=1,totlms
    read(43,*) lm
    do 
    read(43,*) sh1,sh2,sh3,sh4
    if (sh1<0) exit
    if ((ftypea.and.typea=='p').or.(ftypeb.and.typeb=='p')) then
        lsp(lm,sh1,sh2)=.true.
    else if ((ftypeb.and.typeb=='n').or.(ftypea.and.typea=='n')) then
        lsp(lm,sh3,sh4)=.true.
    end if
    end do
    end do
    close(unit=43)
    end if
    end if
    
    print *,' Nucleon types ab ',typea,typeb
    nsc=0
! allocate for ordering
    allocate(node(max_node))
!    nstring=max_part_tree
    nonode=max_node
    allocate(string_order3(nstring))
    call reset_string_order3()
    trmd=0
    max_lm=0
    if (single) maxlmm=max_sj2u
    if (two) maxlmm=2*max_sj2u
    if (idop.or.lnum) maxlmm=0
    if (ltrd) maxlmm=2*max_sj2u
    
    do idx=1,indexx
        read(1) tm(1:12),trxn
        jl=ichar(tm(ijl:ijl))    
        jr=ichar(tm(ijr:ijr))
        mr=ichar(tm(imr:imr))
        shr=ichar(tm(isr:isr))
        shl=ichar(tm(isl:isl))
        lmm=ichar(tm(ilm:ilm))-128
        max_lm=max(max_lm,abs(lmm))
            call add_string3(tm(1:12),new)
            if (new) trmd=trmd+1
    end do
    
    if (single) print * , ' One Particle Op ',op
    if (two) print * , ' Two Particle Op ',op
    if (ltrd) print * , ' Transition Density Op ',op
    if (lnum) print * , ' Number Op ',op
    if (idop) print * , ' Identity Op ',op
    
    if (trmd/=0) then
    allocate(trma(trmd))
    allocate(trmam(trmd))
    done=.true.
    idx=0
    lsc=0
160 idx=get_next_string3(done)
    lsc=lsc+1
    trma(lsc)(1:12)=string_order3(idx)%string(1:12)
    if (done) go to 170
    go to 160
170 if (trmd/=lsc) then
        print *, ' trmd,lsc',trmd,lsc
        stop '  trmd /= lsc '
    end if
    end if
    print *, ' Opmd = ',trmd 
    
    call reset_string_order3()
    trdd=0 
    if (trmd/=0) then
    do nsc=1,trmd
        tm(1:12)=trma(nsc)(1:12)
        jl=ichar(tm(ijl:ijl))    
        jr=ichar(tm(ijr:ijr))
        lmm=ichar(tm(ilm:ilm))-128
        shl=ichar(tm(isl:isl))
        shr=ichar(tm(isr:isr))
        npl=ichar2(tm(inl:inl+1))
        npr=ichar2(tm(inr:inr+1))
        lmn=abs(csj2(shr)-csj2(shl))
        lmx=abs(csj2(shr)+csj2(shl))
        do lm=lmn,lmx,2
            if ((llsp.and.lsp(lm/2,shl,shr).and.ltrd).or..not.ltrd.or..not.llsp) then
            ta(1:5)=tm(1:5)
            ta(ilm:ilm)=char(lm)
            ta(6:12)=z1(6:12)
            call add_string3(ta(1:12),new)
            if (new) then
                trdd=trdd+1
            end if
            end if
        end do
    end do  
    end if
    
   if (trdd/=0) then   
   allocate(trda(trdd))
   allocate(trdav(trdd))
   allocate(trddv(trdd))
   allocate(posn(no_spart,no_spartr,trdd))
   allocate(dimn(no_spart,no_spartr,trdd))
   allocate(nplr(trdd))
   do nsc =1,trdd
        call allocatela(nplr(nsc),no_spart,no_spartr)
   end do
   posn=0
   dimn=0
   done=.true.
   idx=0
   lsc=0
260 idx=get_next_string3(done)
   lsc=lsc+1
   trda(lsc)(1:12)=string_order3(idx)%string(1:12)
   if (done) go to 270
   go to 260
270 if (trdd/=lsc) then
        print *, ' trdd,lsc ',trdd,lsc
        stop '  trdd /= lsc '
    end if
    end if
! allocate memory for trdm
    deallocate(string_order3)
    deallocate(node)
    do nsc=1,trmd
        tm(1:12)=trma(nsc)(1:12)
        jl=ichar(tm(ijl:ijl))    
        jr=ichar(tm(ijr:ijr))
        shl=ichar(tm(isl:isl))
        shr=ichar(tm(isr:isr))
        npl=ichar2(tm(inl:inl+1))
        npr=ichar2(tm(inr:inr+1))
        lmm=ichar(tm(ilm:ilm))-128
        lmn=abs(csj2(shr)-csj2(shl))
        lmx=abs(csj2(shr)+csj2(shl))
        do lm=lmn,lmx,2
            if ((llsp.and.lsp(lm/2,shl,shr).and.ltrd).or..not.ltrd.or..not.llsp) then
            ta(1:5)=tm(1:5)
            ta(ilm:ilm)=char(lm)
            ta(6:12)=z1(6:12)
            lsc=tfind2(trda,trdd,ta(1:8))
            if (lsc/=0.and.nJTar(npr,jr/2)/=0) then
                nplr(lsc)%lvect(npr)=npl.or.nplr(lsc)%lvect(npr)
            end if
            end if
        end do
    end do 
      
    
    trddu=trdd  
    do nsc=1,trdd
        ta(1:12)=trda(nsc)(1:12)
        jl=ichar(ta(ijl:ijl))  
        jr=ichar(ta(ijr:ijr))
        lsc=0
        do npr=1,no_spartr
        do npl=1,no_spart
            if (npl.and.nplr(nsc)%lvect(npr)) then
               posn(npl,npr,nsc)=lsc
               dimn(npl,npr,nsc)=nJTa(npl,jl/2)*nJTar(npr,jr/2)
               lsc=lsc+dimn(npl,npr,nsc)
            end if
        end do       
        end do
        trddv(nsc)=lsc
    end do           
    
    if (trdd/=0) then
    do nsc=1,trdd
        ta(1:12)=trda(nsc)(1:12)
        jl=ichar(ta(ijl:ijl))  
        jr=ichar(ta(ijr:ijr))
        lm=ichar(ta(ilm:ilm))
        if (trddv(nsc)/=0) then
            no_Opd=no_Opd+trddv(nsc)
        end if
    end do
    end if
    print *, ' Opdd = ',trddu
    

    if (trmd/=0) then   
    do nsc=1,trmd
        tm(1:12)=trma(nsc)(1:12)
        npl=ichar2(tm(inl:inl+1))
        npr=ichar2(tm(inr:inr+1))
        jl=ichar(tm(ijl:ijl))    
        jr=ichar(tm(ijr:ijr))
        allocate(trmam(nsc)%t(nJTa(npl,jl/2),nJTar(npr,jr/2)))
        trmam(nsc)%t=real(0,rkw)
    end do 
    end if
    rewind(1)
    read(1) indexx

    do idx=1,indexx 
        read(1) tm(1:12),trxn
        jl=ichar(tm(ijl:ijl))    
        jr=ichar(tm(ijr:ijr))
        lmm=ichar(tm(ilm:ilm))-128
        lsc=tfind3(trma,trmd,tm(1:12))
        if (lsc==0) then
            write(2,*) ' Cannot find trma for idx ',idx,' & for tm =',ichara3(tm,12)
            warn=.true.
            nmiss=nmiss+1 
            go to 54
        end if
        phs=trxn%qux(3)
        ipl=trxn%qux(1)
        ipr=trxn%qux(2)
        sjl=trxn%qux(4)
        sjr=trxn%qux(5)
        if (sjr<0) go to 54
        if (sjl<0) go to 54
        npl=ichar2(tm(inl:inl+1))
        npr=ichar2(tm(inr:inr+1))
        do nJR=1,nJTar(npr,sjr/2)
        do nJL=1,nJTa(npl,sjl/2)
            trmam(lsc)%t(nJL,nJR)=trmam(lsc)%t(nJL,nJR)+&
                    real(phs)*proj_array(npl,sjl/2)%coef(ipl,nJL)*&
                        proj_arrayr(npr,sjr/2)%coef(ipr,nJR)
        end do
        end do
54  end do
    close(unit=1)
    
    if (trmd/=0) then
    do nsc=1,trmd
    ta(1:12)=trma(nsc)(1:12)
    jl=ichar(ta(ijl:ijl))    
    jr=ichar(ta(ijr:ijr))
    shl=ichar(ta(isl:isl))
    shr=ichar(ta(isr:isr))
    lmm=ichar(ta(ilm:ilm))-128
    lmt=ichar(ta(ilt:ilt))-128
    mr=ichar(ta(imr:imr))-128
    ml=ichar(ta(iml:iml))-128
    npl=ichar2(ta(inl:inl+1))
    npr=ichar2(ta(inr:inr+1))
    lmn=abs(csj2(shr)-csj2(shl))
    lmx=abs(csj2(shr)+csj2(shl))
        do lm=lmn,lmx,2
        if ((llsp.and.lsp(lm/2,shl,shr).and.ltrd).or..not.ltrd.or..not.llsp) then
        clebsch1=clebg(jl,jr,lm,jl,-jr,jl-jr)
        if (abs(clebsch1)>0.00000) then
            tn(1:5)=ta(1:5)
            tn(ilm:ilm)=char(lm)
            tn(6:12)=z1(6:12)
            lsc=tfind2(trda,trdd,tn(1:8))
            if (lsc/=0) then
            if (trddv(lsc)/=0) then
                if (single) clebsch2=clebg(csj2(shr),csj2(shl),lm,mr,ml,ml+mr)*sqrt2m1
                if (ltrd.or.lnum) clebsch2=clebg(csj2(shl),csj2(shr),lm,ml,-mr,ml-mr)
                if (lnum) clebsch2= clebsch2*sqrt(real(csj2(shl)+1,8))
                if (two) then
                    clebsch2=clebg(csj2(shr),csj2(shl),lm,mr,ml,mr+ml)*sqrt3m1
                    if (shr==shl.and.btest(lm/2,0)) then
                        clebsch2=real(0)
                    else if (shr==shl.and..not.btest(lm/2,0)) then
                        clebsch2=clebsch2*sqrt2m1
                    end if
                end if
                if (idop) clebsch2=real(1,rkw)
                if (abs(clebsch2)>0.00000) then
                if (.not.allocated(trdav(lsc)%v)) then
                    allocate(trdav(lsc)%v(trddv(lsc)))
                    trdav(lsc)%v=real(0,rkw)
                end if                    
                    c1c2=clebsch2/clebsch1
                    call avs(c1c2*trmam(nsc)%t,int(posn(npl,npr,lsc)),&
                            nJTa(npl,jl/2),nJTar(npr,jr/2),trdav(lsc)%v)
                end if    
            end if
            end if
        end if
        end if
        end do
    end do         
    end if

90  if (allocated(trma)) deallocate(trma)
    allocate(trdao(trdd))
    
    do nsc=1,trdd
        ta(1:12)=trda(nsc)(1:12)
        trdao(nsc)(1:8)=ta(1:8)
        if (.not.allocated(trdav(nsc)%v)) then 
            trddv(nsc)=0
            go to 100
        end if  
        if (output_control>5) then
            jl=ichar(ta(ijl:ijl))    
            jr=ichar(ta(ijr:ijr))
            write(333,'(i9,8i5)') nsc,ichara(trda(nsc),8)
            do npr=1,no_spartr
            do npl=1,no_spart
                c1c2=sum(trdav(nsc)%v&
                (posn(npl,npr,nsc)+1:posn(npl,npr,nsc)+dimn(npl,npr,nsc)))
                if (abs(c1c2)>=.0001) write(333,'(i9,2i5,f9.4)') nsc,npl,npr,c1c2
            end do
            end do
        end if
100 end do 
    
    words=(no_spart-1)/32 + 1
    write(baseo) trdd,2*max_sj2ud,no_Opd,words,no_spartr,no_spart,op
    if (trdd/=0) then
        write(baseo) trdao
        write(baseo) trddv
        do nsc=1,trdd
        if (trddv(nsc)/=0) then
            write(baseo) trdav(nsc)%v
            if (output_control>5) write(333,'(i9,8i5)') nsc,ichara(trda(nsc),8)
            if (output_control>4) write(333,'(8f9.4)') trdav(nsc)%v
            write(baseo) (nplr(nsc)%lvect(npr)%lword(1:words),npr=1,no_spartr)
            if (output_control>5) call lwrite(333,nplr(nsc),no_spart,no_spartr)
        end if
        end do
    end if
    
    if (warn) then
        print * ,nmiss ,'OPM`s not found'
        print * , ' See ',opfile//'.mst', ' for details '
    end if
    
   end subroutine Opd_Opm
   
   pure subroutine avs(atd,pos,diml,dimr,av)
        implicit none
        integer,intent(in):: pos,diml,dimr
        real(kind=rkw),dimension(:,:),intent(in):: atd
        real(kind=rkw),dimension(:),intent(inout):: av
        integer:: j,mv,nv,lv
        lv=pos
        mv=diml
        nv=lv-mv
        do j=1,dimr
            nv=nv+mv
            av(1+nv:mv+nv)=av(1+nv:mv+nv)+atd(1:mv,j)
        end do
   end subroutine avs
                            
end module Opd

!  NuOpd.f90 
!
!  FUNCTIONS:
!  NuOpd      - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: NuOpd
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program NuOpd
    
    use Opd
    use OperParameters
    use ProjParameters
    use InputNuShell
    use OutputNuShell
    use ShellsX
    use Mscheme
    use Files
    use Extra
    
    implicit none

    ! Variables
    
    real:: t1,t2,et1,et2,ms
    integer:: hr,mn,sc

    integer:: j,t,k,dt,k1,res,i
    
    character(len=6):: inExt,Nucl
    character(len=10):: btime
    character(len=8):: bdate
    character(len=10)::chnomdim,chnostring,chstring='NOCSTRINGS'
    integer:: status,nmp1,nmp2
    character(len=14)::buffer
    character(len=1):: opf
        
    ! Body of NuOpd
    
    
    res=get_env_v(chstring,chnostring)
    if (res /= 0) then
        read(chnostring,*) nstring
    else
        nstring=max_part_tree
    end if
    call get_cmd_arg(1,buffer,status) 

    call output_header(6)  
    if (status == -1 ) then  
        print *, ' NuOpd Program'
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
    del_cJ21=del_cJ2
    del_cT21=del_cT2
    Tparity1=Tparity
    nmp1=nmp
    call output_shells(6)    
    call input(inExt2)
    Nucleus2=Nucleus
    max_cJ22=max_cJ2
    max_cT22=max_cT2
    min_cJ22=min_cJ2
    min_cT22=min_cT2
    del_cJ22=del_cJ2
    del_cT22=del_cT2
    Tparity2=Tparity
    nmp2=nmp
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
    call cpu_time(t1)
    call date_and_time(bdate,btime)
    call input_op(opin,opf)
    Nucl=adjustr(Nucleus2)
    opfile(1:1)=' '
    opfile(2:6)=Nucl(1:5)
    opfile(7:7)=opf
    Nucleus=Nucleus1
    max_cJ2=max_cJ21
    max_cT2=max_cT21
    min_cJ2=min_cJ21
    min_cT2=min_cT21
    nmp=nmp1
    call open_files_read(no_spart,'.prj',200,dt)
    allocate(proj_array(no_spart,min_cJ2/2:max_cJ2/2))
    allocate(nJTa(no_spart,min_cJ2/2:max_cJ2/2))
    allocate(spart(no_spart))
    k=0
    k1=0
    do j=min_cJ2,max_cJ2+2,2
        do t=min_cT2,max_cT2+dt,2
            if (j <= max_cJ2  .and. t <= max_cT2) then
                call read_proj(200,k)
            end if
            if (j <= max_cJ2  .and. t <= max_cT2) k=k+1
            k1=k1+1
        end do
    end  do    
    call close_files_read(200,dt)
    Nucleus=Nucleus2
    max_cJ2=max_cJ22
    max_cT2=max_cT22
    min_cJ2=min_cJ22
    min_cT2=min_cT22
    nmp=nmp2
    call open_files_read(no_spartr,'.prj',200,dt)
    allocate(proj_arrayr(no_spartr,min_cJ2/2:max_cJ2/2))
    allocate(nJTar(no_spartr,min_cJ2/2:max_cJ2/2))
    allocate(spartr(no_spartr))
    k=0
    k1=0
    do j=min_cJ2,max_cJ2+2,2
        do t=min_cT2,max_cT2+dt,2
            if (j <= max_cJ2  .and. t <= max_cT2) then
                call read_projr(200,k)
            end if
            if (j <= max_cJ2  .and. t <= max_cT2) k=k+1
            k1=k1+1
        end do
    end do
    call close_files_read(200,dt)

    open(unit=800,file=opfile//'.trd',form='unformatted',action='write')
    open(unit=333,file=opfile//'.drt',action='write')
    open(unit=801,file=opfile//'.tri',form='unformatted',action='write')
    open(unit=802,file=opfile//'.trf',form='unformatted',action='write')
    call Opd_Opm(800)
    close(unit=800)
    write(801) no_spartr,min_cJ22,max_cJ22,proton,neutron,x2,Tparity2,op,del_cJ22
    write(801) nJTar(1:no_spartr,min_cJ22/2:max_cJ22/2)
    write(801) spartr(1:no_spartr)
    close(unit=801)
    write(802) no_spart,min_cJ21,max_cJ21,proton,neutron,x1,Tparity1,op,del_cJ21
    write(802) nJTa(1:no_spart,min_cJ21/2:max_cJ21/2)
    write(802) spart(1:no_spart)
    close(unit=802)
!    print *, ' Del_J2a, Del_J2b ', del_cJ21,del_cJ22
    print *, ' No of OPD matrix elements ',no_Opd,proton,neutron
    read(btime,'(i2,i2,i2,f4.3)') hr,mn,sc,ms
    et1=hr*3600.+mn*60.+sc*1. +ms
    call date_and_time(bdate,btime)
    call cpu_time(t2)
    read(btime,'(i2,i2,i2,f4.3)') hr,mn,sc,ms
    et2=hr*3600.+mn*60.+sc*1. +ms
    t1=t2-t1
    et1=et2-et1
    if (et1<0.0) et1=et1+24.*3600.
    call output_completed(6,'NuPrd/NuTra')    
    call output_time('NuOpd',t1,et1)
    print *
    
    end program NuOpd
    
