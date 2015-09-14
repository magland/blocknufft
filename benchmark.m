close all;

if 1
N1=200; N2=200; N3=200;
Nr=100; Ntheta=100; Nphi=100;
M=Nr*Ntheta*Nphi;
fprintf('M = %d\n',M);
r=(0.5:1:Nr-0.5)/Nr;
theta=(0.5:1:Ntheta-0.5)/Ntheta*2*pi;
phi=(0.5:1:Nphi-0.5)/Nphi*pi-pi/2;
[R,THETA,PHI]=ndgrid(r,theta,phi);
x0=R.*cos(THETA).*cos(PHI); x0=x0(:)*pi+pi;
y0=R.*sin(THETA).*cos(PHI); y0=y0(:)*pi+pi;
z0=R.*sin(PHI); z0=z0(:)*pi+pi;
xyz=cat(2,x0,y0,z0);
d=rand(M,1)*2-1;
eps=1e-5;
K1=50000; K2=50000; K3=50000;
num_threads=1;
%incr=25;
%figure; plot3(x0(25:incr:end),y0(25:incr:end),z0(25:incr:end),'.');
end

if 0
N1=200; N2=200; N3=200;
M=10^6;
fprintf('M = %d\n',M);
xyz=rand(M,3)*2*pi;
d=rand(M,1)*2-1;
eps=1e-5;
K1=50000; K2=50000; K3=50000;
num_threads=1;
end

if 1
disp('***** nufft3d1f90 *****');
tic;
[A1,ierr]=nufft3d1f90(xyz(:,1),xyz(:,2),xyz(:,3),d,0,eps,N1,N2,N3);
toc
%figure; imagesc(squeeze(real(A1(:,:,6)))); colormap('gray'); drawnow;
writemda(A1,'A1.mda');
%writemda(spread,'A1_spread.mda');
%writemda(fftn(spread),'A1_spread_fft.mda');
end

if 1
disp('***** New implementation, single thread, blocking off *****');
tic
[A2,spread]=blocknufft3d(N1,N2,N3,xyz,d,eps,K1,K2,K3,num_threads);
toc
%figure; imagesc(squeeze(real(A2(:,:,6)))); colormap('gray'); drawnow;
writemda(A2,'A2.mda');
end

fprintf('Max difference in images: %.10f\n',max(abs(A1(:)-A2(:))));

if 1
disp('***** New implementation, single thread, blocking on *****');
tic
K1=100; K2=100; K3=100;
A3=blocknufft3d(N1,N2,N3,xyz,d,eps,K1,K2,K3,num_threads);
toc
%figure; imagesc(squeeze(abs(A3(1,:,:)))); colormap('gray'); drawnow;
end;

if 1
disp('***** New implementation, six threads, blocking on *****');
tic
num_threads=6;
A4=blocknufft3d(N1,N2,N3,xyz,d,eps,K1,K2,K3,num_threads);
toc
%figure; imagesc(squeeze(abs(A4(1,:,:)))); colormap('gray'); drawnow;
end;


