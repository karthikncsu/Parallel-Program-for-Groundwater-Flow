clc
clear all
close all

% Reading the results from FORTRAN code
fileID1=fopen('solution_parallel008_008.dat','r');

dum=fscanf(fileID1,['%f','%f','%f'],[3 Inf]);
[n Ntot]=size(dum);
fclose(fileID1)

syms x
eqn = (x+1)*(x+2) ==2*Ntot;
solx = solve(eqn,x);

Nx=double((solx(2)))+1;
Ny=double((solx(2)))/2+1;
xFOR=dum(1,:);
 yFOR=dum(2,:);
headFOR=dum(3,:);
dx=1000/(Nx-1);
dy=500/(Ny-1);
loc=zeros(Nx*Ny,2);
for iy=1:Ny
    for ix=1:Nx
        loc((iy-1)*Nx+ix,1)=round(xFOR((iy-1)*Nx+ix)/dx)+1;
        loc((iy-1)*Nx+ix,2)=round(yFOR((iy-1)*Nx+ix)/dy)+1;
    end
end

xFOR1=linspace(0,1000,Nx);
yFOR1=linspace(0,500,Ny);
[XX, YY]=meshgrid(xFOR1,yFOR1);

    head=zeros(Nx,Ny);
    for iy=1:Ny
    for ix=1:Nx
        indx=loc((iy-1)*Nx+ix,1);
        indy=loc((iy-1)*Nx+ix,2);
%         [indx, indy]
        head(indx,indy)=headFOR((iy-1)*Nx+ix);
    end
    end
   
   contourf(XX',YY',head)
   colorbar

   figure
   title('Hydraulic head along x-direction at  y=250 m')
   xlabel('x(m)')
   ylabel('head(m)')
   hold on
   plot(xFOR1,head(:,(Ny-1)/2+1),'Linewidth',2)
   
    figure
   title('Hydraulic head along y-direction at  x=500 m')
   xlabel('y(m)')
   ylabel('head(m)')
   hold on
   plot(yFOR1,head((Nx-1)/2+1,:),'Linewidth',2)
   
   a=-1e-3; b=-0.2e-2; c=0.4e-5; d=2.0
   ZZ=a*XX+b*YY+c*XX.*YY+d;
   
   figure
   title('Hydraulic Conductivity')
   xlabel('x(m)')
   ylabel('y(m)')
   hold on
   contourf(XX',YY',ZZ')
   colorbar
   
   