        program miniproject
        IMPLICIT NONE
        include 'mpif.h'
                
        DOUBLE PRECISION a,b,c,d,x0,y0,hx,hy,x,y
        DOUBLE PRECISION t0,dt,tf,Lx,Ly,Q0,Q,Ss,kxy,dum1,dum2
        INTEGER dnx,dny,nx,ny,iproc,ndims,left,right,bottom,top,reorder
        INTEGER nprocs,nprocx,nprocy,myid,ierr,comm2d,diff,ipnt1,ipnt2
        INTEGER,dimension(:),allocatable :: dims,isoperiodic
        INTEGER,dimension(:),allocatable :: nxproc,nyproc,nxtot,nytot
        INTEGER proccoord(2),iprocx,iprocy,status,dum
        INTEGER,dimension(:,:),allocatable ::proccoordall,proccordtomyid

        DOUBLE PRECISION,dimension(:),allocatable :: head,head0
        DOUBLE PRECISION,dimension(:),allocatable :: kx,ky,k0,sendvals
        DOUBLE PRECISION,dimension(:,:),allocatable :: coord

        INTEGER i,j
        DOUBLE PRECISION stime,etime,time,pi,MFLOPS,MFLOPSt,L1,L1norm
        DOUBLE PRECISION, dimension (:), allocatable::lnodrval,rnodrval
        DOUBLE PRECISION, dimension (:), allocatable::bnodrval,tnodrval
        DOUBLE PRECISION, dimension (:), allocatable::lnodsval,rnodsval
        DOUBLE PRECISION, dimension (:), allocatable::bnodsval,tnodsval
        integer, dimension (:), allocatable :: lnodr,rnodr,bnodr,tnodr
        integer, dimension (:), allocatable :: lnods,rnods,bnods,tnods
        
        DOUBLE PRECISION, dimension (:), allocatable :: dumrnods
        integer, dimension (:), allocatable :: bnodin,tnodin
        integer lnodrcnt,rnodrcnt,bnodrcnt,tnodrcnt
        integer lnodscnt,rnodscnt,bnodscnt,tnodscnt
        integer bnodincnt,tnodincnt
        integer inod0,inod1,inod2,inod3,inod4
        integer tag1,tag2,tag3,tag4
        character(len=100) :: outname        
                
        call MPI_init(ierr)
        call MPI_Comm_rank(MPI_COMM_WORLD,myid,ierr )
        call MPI_Comm_size(MPI_COMM_WORLD,nprocs,ierr)
                
        a=-1e-3; b=-0.2e-2; c=0.4e-5; d=2.d0
        Lx=1000; Ly=500; Ss=0.1; Q0=0.001
        nx=126; ny=63
        t0=0; dt=1e-3; tf=1000
        hx=Lx/(nx-1);   hy=Ly/(ny-1)
        
        if(nprocs==1) then
! Sequential programing
        call singleproc(a,b,c,d,t0,dt,tf,Lx,Ly,Q0,Ss,nx,ny)
        else
! Parallel programming
! 2D domain decomposition of processors

        nprocx=nprocs; nprocy=1;
        diff=abs(nprocx-nprocy)
        do iproc=1,nprocs/2
        if(mod(nprocs,iproc)==0) then
                if(abs(nprocs/iproc-iproc)<diff) then
                    nprocx=iproc; nprocy=nprocs/iproc
                    diff=abs(nprocs/iproc-iproc)
                endif
        endif
        enddo
        
        if(nprocx<nprocy) then
        diff=nprocx; nprocx=nprocy; nprocy=diff
        endif
        
        if(myid==nprocs-1) then
        write(*,*) "Number of Processors in x-direction:",nprocx
        write(*,*) "Number of Processors in y-direction:",nprocy
        endif
        
        ndims=2
        allocate(dims(1:ndims),isoperiodic(1:ndims))
        
        isoperiodic(1)=0; isoperiodic(2)=0
        reorder=1
        dims(1)=nprocx;   dims(2)=nprocy
! End of 2D Domain Decomposition        
        
      	call MPI_Cart_create(MPI_COMM_WORLD,ndims,dims,isoperiodic,
     &        reorder,comm2d,ierr)
        call MPI_Comm_rank(comm2d,myid,ierr)
        call MPI_cart_shift(comm2d,0,1,left,right,ierr)
        call MPI_cart_shift(comm2d,1,1,bottom,top,ierr)
        
! FInding the coordinates of all the processors
        call MPI_CART_COORDS(comm2d,myid,2,proccoord,ierr)
!        Write(*,*) myid,"Processor coordiantes",proccoord
        
        allocate(proccoordall(1:2,1:nprocs)) 
        call MPI_ALLGATHER(proccoord,2,MPI_INTEGER,proccoordall,2,
     &   MPI_INTEGER,comm2d,ierr)
     
!         write(*,*) myid,proccoordall
        allocate(proccordtomyid(0:nprocx-1,0:nprocy-1))
        proccordtomyid=0
         
         do iproc=1,nprocs
         ipnt1=proccoordall(1,iproc)
         ipnt2=proccoordall(2,iproc)

         proccordtomyid(ipnt1,ipnt2)=iproc-1
         
         enddo
         
!         write(*,*) myid,left,right,bottom,top
!         write(*,*) "processor coordinates to id:",proccordtomyid

!----------------------------------------------------------
! Dividing the mesh into multiple processors
!----------------------------------------------------------
        if(myid==0) then
        allocate(nxproc(0:nprocx-1),nyproc(0:nprocy-1))

        nxproc=0; nyproc=0
        dnx=nx/nprocx;   dny=ny/nprocy

! gird points along x-direction
        do iprocx=0,nprocx-2
        nxproc(iprocx)=dnx
        enddo
        nxproc(nprocx-1)=nx-sum(nxproc)

! gird points along y-direction
        do iprocy=0,nprocy-2
        nyproc(iprocy)=dny
        enddo
        nyproc(nprocy-1)=ny-sum(nyproc)

        do iprocy=nprocy-1,0,-1
        do iprocx=nprocx-1,0,-1
        iproc=iprocy*nprocx+iprocx

        if(iproc.ne.0) then
        call MPI_SSend(nxproc(iprocx),1,MPI_INTEGER
     &  ,proccordtomyid(iprocx,iprocy),1,comm2d,ierr)
     
        call MPI_SSend(nyproc(iprocy),1,MPI_INTEGER
     &  ,proccordtomyid(iprocx,iprocy),2,comm2d,ierr)

        call MPI_SSend(sum(nxproc(0:iprocx-1))*hx,1,MPI_DOUBLE_PRECISION
     &  ,proccordtomyid(iprocx,iprocy),3,comm2d,ierr)

        call MPI_SSend(sum(nyproc(0:iprocy-1))*hy,1,MPI_DOUBLE_PRECISION
     &  ,proccordtomyid(iprocx,iprocy),4,comm2d,ierr)

        else
        nx=nxproc(iprocx);  ny=nyproc(iprocy)
        x0=0; y0=0
        endif
        
        enddo
        enddo
        
        write(*,*) "---------------------------------------------------"
        write(*,*) "myid:",myid,"nx:",nx,"ny:",ny
        write(*,*) "myid:",myid,"x0:",x0,"y0:",y0
        write(*,*)"myid:",myid,"xend:",x0+hx*(nx-1),"yend:",y0+hy*(ny-1)
        write(*,*) "---------------------------------------------------"

        else
       
        call MPI_RECV(nx,1,MPI_INTEGER,0,1,comm2d,status,ierr)
        call MPI_RECV(ny,1,MPI_INTEGER,0,2,comm2d,status,ierr)
        call MPI_RECV(x0,1,MPI_DOUBLE_PRECISION,0,3,comm2d,status,ierr)
        call MPI_RECV(y0,1,MPI_DOUBLE_PRECISION,0,4,comm2d,status,ierr)
        
        write(*,*) "---------------------------------------------------"
        write(*,*) "myid:",myid,"nx:",nx,"ny:",ny
        write(*,*) "myid:",myid,"x0:",x0,"y0:",y0
        write(*,*)"myid:",myid,"xend:",x0+hx*(nx-1),"yend:",y0+hy*(ny-1)
        write(*,*) "---------------------------------------------------"
        endif
!----------------------------------------------------------
! End of Dividing the mesh into multiple processors
!----------------------------------------------------------

!----------------------------------------------------------
! Solving on multiple processors
!----------------------------------------------------------
        
        allocate(head(1:(nx+2)*(ny+2)),coord(1:(nx+2)*(ny+2),1:2))
        allocate(kx(1:(nx+1)*(ny+2)),ky(1:(nx+2)*(ny+1)))
        allocate(k0(1:(nx+2)*(ny+2)),head0(1:(nx+2)*(ny+2)))

        
        pi=3.1415926535897932384626433832795028841971693993
! Calculating the coordinates of different mesh node points        
        coord=0.d0
        x0=x0-hx; y0=y0-hy
        coord=0.d0
        do i=1,nx+2
        do j=1,ny+2
        inod0=(ny+2)*(i-1)+j
        
        coord(inod0,1)=x0+dble(i-1)*hx; coord(inod0,2)=y0+dble(j-1)*hy
        enddo
        enddo
!        write(*,*) myid,coord(:,1)
! Calculating the values of hydralic conductivity along x
        kx=0.d0
        do i=1,nx+1
        do j=1,ny+2
            inod0=(i-1)*(ny+2)+j
            x=0.5*coord((i-1)*(ny+2)+j,1)+0.5*coord(i*(ny+2)+j,1)
            y=0.5*coord((i-1)*(ny+2)+j,2)+0.5*coord(i*(ny+2)+j,2)
            kx(inod0)=a*x+b*y+c*x*y+d
            kx(inod0)=kx(inod0)*0.5/hx**2
        enddo
        enddo

! Calculating the values of hydralic conductivity along y
        ky=0.d0
        do i=1,nx+2
        do j=1,ny+1
            inod0=(i-1)*(ny+1)+j
            x=0.5*coord((i-1)*(ny+2)+j,1)+0.5*coord((i-1)*(ny+2)+j+1,1)
            y=0.5*coord((i-1)*(ny+2)+j,2)+0.5*coord((i-1)*(ny+2)+j+1,2)
            ky(inod0)=a*x+b*y+c*x*y+d
            ky(inod0)=ky(inod0)*0.5/hy**2
        enddo
        enddo

! Calculating the values of hydralic conductivity at each node
        k0=0.d0
        do i=2,nx+1
        do j=2,ny+1
            inod0=(i-1)*(ny+2)+j
            k0(inod0)=kx((i-2)*(ny+2)+j)+kx((i-1)*(ny+2)+j)
            k0(inod0)=k0(inod0)+ky((i-1)*(ny+1)+j-1)+ky((i-1)*(ny+1)+j)
        enddo
        enddo
        
! Left and right boundary node

        lnodrcnt=ny; allocate(lnodr(1:ny),lnodrval(1:ny))
        rnodrcnt=ny; allocate(rnodr(1:ny),rnodrval(1:ny))

        lnodscnt=ny; allocate(lnods(1:ny),lnodsval(1:ny))
        rnodscnt=ny; allocate(rnods(1:ny),rnodsval(1:ny))
        
        do j=1,ny
        lnodr(j)=j+1
        rnodr(j)=(nx+1)*(ny+2)+j+1

        lnods(j)=j+1+(ny+2)
        rnods(j)=(nx)*(ny+2)+j+1
        enddo

! Bottom and top boundary nodes

        bnodrcnt=nx;  allocate(bnodr(1:nx),bnodrval(1:nx))
        tnodrcnt=nx;  allocate(tnodr(1:nx),tnodrval(1:nx))
        
        bnodscnt=nx;  allocate(bnods(1:nx),bnodsval(1:nx))
        tnodscnt=nx;  allocate(tnods(1:nx),tnodsval(1:nx))

        tnodincnt=nx; allocate(tnodin(1:tnodincnt))
        bnodincnt=nx; allocate(bnodin(1:bnodincnt))

        do i=1,nx
        bnodr(i)=i*(ny+2)+1
        tnodr(i)=(i+1)*(ny+2)

        bnods(i)=i*(ny+2)+2
        tnods(i)=(i+1)*(ny+2)-1

        bnodin(i)=i*(ny+2)+3
        tnodin(i)=(i+1)*(ny+2)-2
        enddo

! intial conditions
        head=0.d0
        head=10.d0-coord(:,1)/200.d0
        head0=head
        time=t0

        call MPI_cart_shift(comm2d,0,1,left,right,ierr)
        call MPI_cart_shift(comm2d,1,1,bottom,top,ierr)
!        write(*,*) myid,left,right,bottom,top

        stime=MPI_Wtime()
        do while(time<tf)
        
        L1=sum(abs(head0))

        call MPI_REDUCE(L1,L1norm,1,MPI_DOUBLE_PRECISION
     &  ,MPI_SUM,0,comm2d,ierr)

!        if(myid==0) write(*,*) myid,time,L2norm
        
        if(proccoord(1).ne.nprocx-1)   rnodsval=head0(rnods)
        if(proccoord(1).ne.0)          lnodsval=head0(lnods)
        if(proccoord(2).ne.nprocy-1)   tnodsval=head0(tnods)
        if(proccoord(2).ne.0)          bnodsval=head0(bnods)
        
        call MPI_Sendrecv(rnodsval,ny,MPI_DOUBLE_PRECISION
     &  ,right,1,lnodrval,ny,MPI_DOUBLE_PRECISION,left,1
     &  ,comm2d,status,ierr)
     
        call MPI_Sendrecv(lnodsval,ny,MPI_DOUBLE_PRECISION
     &  ,left,2,rnodrval,ny,MPI_DOUBLE_PRECISION,right,2
     &  ,comm2d,status,ierr)
     
        call MPI_Sendrecv(tnodsval,nx,MPI_DOUBLE_PRECISION
     &  ,top,3,bnodrval,nx,MPI_DOUBLE_PRECISION,bottom,3
     &  ,comm2d,status,ierr)
        
        call MPI_Sendrecv(bnodsval,nx,MPI_DOUBLE_PRECISION
     &  ,bottom,4,tnodrval,nx,MPI_DOUBLE_PRECISION,top,4
     &  ,comm2d,status,ierr)


        ! boundary conditions
        if(proccoord(1)==nprocx-1) then
             head0(rnodr)=5.d0
        else
             head0(rnodr)=rnodrval
        endif

        if(proccoord(1)==0)         then
             head0(lnodr)=10.d0
        else
             head0(lnodr)=lnodrval
        endif

        if(proccoord(2)==nprocy-1)  then
           head0(tnodr)=head0(tnodin)
        else
           head0(tnodr)=tnodrval
        endif

        if(proccoord(2)==0)         then
           head0(bnodr)=head0(bnodin)
        else
           head0(bnodr)=bnodrval
        endif
        
        do i=2,nx+1
        do j=2,ny+1
        inod0=(ny+2)*(i-1)+j
        inod1=(ny+2)*(i-2)+j
        inod2=(ny+2)*(i)+j
        inod3=inod0-1
        inod4=inod0+1
        
        head(inod0)=kx((i-2)*(ny+2)+j)*head0(inod1)**2
        head(inod0)=head(inod0)+kx((i-1)*(ny+2)+j)*head0(inod2)**2
        head(inod0)=head(inod0)+ky((i-1)*(ny+1)+j-1)*head0(inod3)**2
        head(inod0)=head(inod0)+ky((i-1)*(ny+1)+j)*head0(inod4)**2
        
        head(inod0)=head(inod0)-k0(inod0)*head0(inod0)**2
        head(inod0)=head(inod0)+Q0*(1+sin(pi*time/300.d0))
        head(inod0)=dt*head(inod0)/Ss+head0(inod0)
        enddo
        enddo

        head0=head
        time=time+dt

        end do        
        
        etime=MPI_Wtime()
        etime=etime-stime
                
        MFLOPS= 22*ny*nx*(tf/dt+1.d0)/etime*1e-6
        call MPI_REDUCE(MFLOPS,MFLOPSt,1,MPI_DOUBLE_PRECISION
     &  ,MPI_SUM,0,comm2d,ierr)

        if(myid==0) then
        open(100,file="comp_time.dat",ACCESS="Append")        
        Write(*,*) "Comp. time:",etime,"MFLOPS:",MFLOPSt
        write(*,*) "Phy . time:",time,"L1Norm:",L1norm/(nx*ny*nprocs)
        write(100,*) nx*ny*nprocs,int(time),nprocs,etime,int(MFLOPSt)
     &   ,L1norm/(nx*ny*nprocs)
        close(100)
        endif

!----------------------------------------------------------------
! Processor 0 receiving the data from other processors and writing to 
! a single file
!----------------------------------------------------------------
        allocate(sendvals(1:6))   
        sendvals=0.d0

        if(myid.ne.0) then

! other processors sending the data to first processor

         call MPI_Send(nx,1,MPI_INTEGER,0,(myid-1)*4+5,comm2d,ierr)
         call MPI_Send(ny,1,MPI_INTEGER,0,(myid-1)*4+6,comm2d,ierr)

         call MPI_Send(coord,2*(nx+2)*(ny+2),MPI_DOUBLE_PRECISION,0,
     &   (myid-1)*4+7,comm2d,ierr)
     
         call MPI_Send(head0,(nx+2)*(ny+2),MPI_DOUBLE_PRECISION,0,
     &   (myid-1)*4+8,comm2d,ierr)

        else
! Processor 0 is writing the data

        write (outname,"('solution_parallel',I3.3,'_',I3.3,'.dat')") 
     &  nprocx,nprocy   
        open(26,file=outname,action="write")
 
        do j=2,ny+1
        do i=2,nx+1
        inod0=(ny+2)*(i-1)+j
        write(26,*) coord(inod0,1),coord(inod0,2),head0(inod0)
        enddo
        enddo

        deallocate(head0,coord);        
        write(*,*) "Writing data using processor:",myid
        write(*,*) "Saving into the file:",outname
        
        do iproc=1,nprocs-1
!        write(26,*) "Processor",iproc
        call MPI_RECV(nx,1,MPI_INTEGER,iproc,
     &   (iproc-1)*4+5,comm2d,status,ierr) 

        call MPI_RECV(ny,1,MPI_INTEGER,iproc,
     &   (iproc-1)*4+6,comm2d,status,ierr) 

        allocate(head0(1:(nx+2)*(ny+2)),coord(1:(nx+2)*(ny+2),1:2))

        call MPI_RECV(coord,2*(nx+2)*(ny+2),MPI_DOUBLE_PRECISION,iproc,
     &   (iproc-1)*4+7,comm2d,status,ierr)
     
        call MPI_RECV(head0,(nx+2)*(ny+2),MPI_DOUBLE_PRECISION,iproc,
     &   (iproc-1)*4+8,comm2d,status,ierr)
    
        do j=2,ny+1
        do i=2,nx+1
        inod0=(ny+2)*(i-1)+j
        write(26,*) coord(inod0,1),coord(inod0,2),head0(inod0)
        enddo
        enddo     
        
        deallocate(head0,coord);        
        enddo
        endif
        
                

        
        endif

        call MPI_Finalize ( ierr )
        
        end

!---------------------------------------------------------------------
! Subroutine for solving the problem sequentially
!---------------------------------------------------------------------

        
        subroutine singleproc(a,b,c,d,t0,dt,tf,Lx,Ly,Q0,Ss,nx,ny)
!---------------------------------------------------------------------
! Subroutine for solving the problem sequentially
!---------------------------------------------------------------------
        
        IMPLICIT NONE
        DOUBLE PRECISION a,b,c,d
        DOUBLE PRECISION t0,tf,dt,Lx,Ly,Q0,Q,Ss
        DOUBLE PRECISION,dimension(:),allocatable :: head,head0
        DOUBLE PRECISION,dimension(:),allocatable :: kx,ky,k0
        DOUBLE PRECISION,dimension(:,:),allocatable :: coord
        INTEGER nx,ny,i,j
        
        DOUBLE PRECISION hx,hy,x0,y0,stime,etime,time,pi,MFLOPS
        integer, dimension (:), allocatable :: lnodr,rnodr,bnodr,tnodr
        integer, dimension (:), allocatable :: bnodin,tnodin
        integer lnodrcnt,rnodrcnt,bnodrcnt,tnodrcnt,bnodincnt,tnodincnt
        integer inod0,inod1,inod2,inod3,inod4
        
        allocate(head(1:(nx+2)*(ny+2)),coord(1:(nx+2)*(ny+2),1:2))
        allocate(kx(1:(nx+1)*(ny+2)),ky(1:(nx+2)*(ny+1)))
        allocate(k0(1:(nx+2)*(ny+2)),head0(1:(nx+2)*(ny+2)))
        
        pi=3.1415926535897932384626433832795028841971693993
! Calculating the coordinates of different mesh node points        
        hx=Lx/(nx-1); hy=Ly/(ny-1); coord=0.d0
        do i=1,nx+2
        do j=1,ny+2
        inod0=(ny+2)*(i-1)+j
        
        coord(inod0,1)=-hx+dble(i-1)*hx; coord(inod0,2)=-hy+dble(j-1)*hy
        enddo
        enddo

! Calculating the values of hydralic conductivity along x
        kx=0.d0
        do i=1,nx+1
        do j=1,ny+1
            inod0=(i-1)*(ny+2)+j
            x0=0.5*coord((i-1)*(ny+2)+j,1)+0.5*coord(i*(ny+2)+j,1)
            y0=0.5*coord((i-1)*(ny+2)+j,2)+0.5*coord(i*(ny+2)+j,2)
            kx(inod0)=a*x0+b*y0+c*x0*y0+d
            kx(inod0)=kx(inod0)*0.5/hx**2
        enddo
        enddo

! Calculating the values of hydralic conductivity along y
        ky=0.d0
        do i=1,nx+1
        do j=1,ny+1
            inod0=(i-1)*(ny+1)+j
            x0=0.5*coord((i-1)*(ny+2)+j,1)+0.5*coord((i-1)*(ny+2)+j+1,1)
            y0=0.5*coord((i-1)*(ny+2)+j,2)+0.5*coord((i-1)*(ny+2)+j+1,2)
            ky(inod0)=a*x0+b*y0+c*x0*y0+d
            ky(inod0)=ky(inod0)*0.5/hy**2
        enddo
        enddo

! Calculating the values of hydralic conductivity at each node
        k0=0.d0
        do i=2,nx+1
        do j=2,ny+1
            inod0=(i-1)*(ny+2)+j
            k0(inod0)=kx((i-2)*(ny+2)+j)+kx((i-1)*(ny+2)+j)
            k0(inod0)=k0(inod0)+ky((i-1)*(ny+1)+j-1)+ky((i-1)*(ny+1)+j)
        enddo
        enddo

! Left and right boundary node
        lnodrcnt=ny; allocate(lnodr(1:lnodrcnt))
        rnodrcnt=ny; allocate(rnodr(1:rnodrcnt))
        
        do j=1,ny
        lnodr(j)=j+1
        rnodr(j)=(nx+1)*(ny+2)+j+1
        enddo

! Bottom and top boundary nodes
        bnodrcnt=nx; allocate(bnodr(1:bnodrcnt))
        tnodrcnt=nx; allocate(tnodr(1:tnodrcnt))

        bnodincnt=nx; allocate(bnodin(1:bnodincnt))
        tnodincnt=nx; allocate(tnodin(1:tnodincnt))

        do i=1,nx
        bnodr(i)=i*(ny+2)+1
        tnodr(i)=(i+1)*(ny+2)

        bnodin(i)=i*(ny+2)+3
        tnodin(i)=(i+1)*(ny+2)-2
        enddo

! intial conditions
        head=0.d0
        head(:)=10.d0-coord(:,1)/200.d0
        head0=head
        time=t0
        
        call cpu_time(stime)
        
        do while(time<tf)

!        write(*,*) time, sum(abs(head))

c        ! boundary conditions
        head0(lnodr)=10.d0;          head0(rnodr)=5.d0
        head0(tnodr)=head0(tnodin);    head0(bnodr)=head0(bnodin)
!        write(*,*) time, sum(abs(head))
        
        do i=2,nx+1
        do j=2,ny+1
        inod0=(ny+2)*(i-1)+j
        inod1=(ny+2)*(i-2)+j
        inod2=(ny+2)*(i)+j
        inod3=inod0-1
        inod4=inod0+1
        
        head(inod0)=kx((i-2)*(ny+2)+j)*head0(inod1)**2
        head(inod0)=head(inod0)+kx((i-1)*(ny+2)+j)*head0(inod2)**2
        head(inod0)=head(inod0)+ky((i-1)*(ny+1)+j-1)*head0(inod3)**2
        head(inod0)=head(inod0)+ky((i-1)*(ny+1)+j)*head0(inod4)**2
        
        head(inod0)=head(inod0)-k0(inod0)*head0(inod0)**2
        head(inod0)=head(inod0)+Q0*(1+sin(pi*time/300.d0))
        head(inod0)=dt*head(inod0)/Ss+head0(inod0)
        enddo
        enddo
        
        head0=head
        time=time+dt

        end do        
        
        call cpu_time(etime)
        etime=etime-stime
                
        MFLOPS= 22*ny*nx*(tf/dt+1.d0)/etime*1e-6
        Write(*,*) "Comp. time:",etime,"MFLOPS:",MFLOPS,
     &   sum(abs(head0))/(nx*ny)

        open(100,file="comp_time.dat",ACCESS="Append")        
        write(100,*) nx*ny,int(time),1,etime,int(MFLOPS),
     &  sum(abs(head0))/(nx*ny)
        close(100)
                
        open(25,file="mini_project_results.dat")
        write(*,*) "done"
        do j=2,ny+1
        do i=2,nx+1
        inod0=(ny+2)*(i-1)+j
        write(25,*) coord(inod0,1),coord(inod0,2),head(inod0)
        enddo
        enddo
        close(25)
        

        end

 
