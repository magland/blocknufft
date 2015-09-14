close all;

N1=10; N2=10; N3=10;
M=1e4;
xyz=rand(M,3)*2*pi;
d=rand(M,1)*2-1; %xyz(1,:)=0.5; d(:)=0; d(1)=1;
eps=1e-5;
K1=50000; K2=50000; K3=50000;
num_threads=1;

% N1=10; N2=10; N3=10;
% M=1;
% xyz=zeros(M,3)*2*pi;
% d=ones(M,1);
% eps=1e-5;
% K1=50000; K2=50000; K3=50000;
% num_threads=1;

if 1
disp('***** nufft3d1f90 *****');
tic;
[A1,ierr]=nufft3d1f90(xyz(:,1),xyz(:,2),xyz(:,3),d,0,eps,N1,N2,N3);
disp(max(abs(A1(:))));
A1=A1/max(abs(A1(:)));
toc
%figure; imagesc(squeeze(real(A1(:,:,6)))); colormap('gray'); drawnow;
writemda(A1,'A1.mda');
%writemda(spread,'A1_spread.mda');
%writemda(fftn(spread),'A1_spread_fft.mda');
end

if 1
disp('***** New implementation, single thread, blocking off *****');
tic
[A2]=blocknufft3d(N1,N2,N3,xyz,d,eps,K1,K2,K3,num_threads);
toc
disp(max(abs(A2(:))));
A2=A2/max(abs(A2(:)));
%figure; imagesc(squeeze(real(A2(:,:,6)))); colormap('gray'); drawnow;
writemda(A2,'A2.mda');
end

%A1=ifftshift(A1);
fprintf('Max difference in images: %.10f\n',max(abs(A1(:)-A2(:))));
fprintf('Max error in A1: %.10f\n',max(abs(A1(:)-1/M)));
fprintf('Max error in A2: %.10f\n',max(abs(A2(:)-1/M)));

if 1
disp('***** New implementation, single thread, blocking on *****');
tic
K1=50; K2=50; K3=50;
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


