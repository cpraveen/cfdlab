%% Taken from
%%  Damien BIAU, The first-digit frequencies in data of turbulent flows
%%  http://arxiv.org/abs/1509.06740v1
%%——————————————————————- %% simulation of Taylor-Green Vortex %%——————————————————————- clear all,
Re = 1600;
N = 2^8;
dt = 5e-3;
Tend=20; 
%——————————————————————-
%% Fourier pseudo-spectral method
x = 2*pi/N*[0:N-1]';
kx = [0:N/2-1 0 -N/2+1:-1]';
for m=1:N
   for n=1:N,
      IKX(:,m,n)=i*kx; 
      IKY(m,:,n)=i*kx; 
      IKZ(m,n,:)=i*kx;
   end
end

K2=-(IKX.^2+IKY.^2+IKZ.^2);
K2p=K2; 
K2p(1,1,1)=1; 
K2p(N/2+1,:,:)=1; 
K2p(:,N/2+1,:)=1; 
K2p(:,:,N/2+1)=1; 
Z=ones(N,N,N); 
Z(N/2+1,:,:)=0;
Z = 1 - (sqrt(K2) > round(N/3)+1);
Z(:,N/2+1,:)=0; 
Z(:,:,N/2+1)=0;
E = exp(-dt*K2/Re); 
E2 = exp(-dt/2*K2/Re);
%% Newcomb-Benford’s law 
NBL=log10(1+1./[1:9]); 
SNB(1)=-sum(NBL.*log10(NBL));
%% initial condition
[Xs, Ys, Zs] = meshgrid(x,x,x); %% 3D grid
uf = fftn( 2/sqrt(3)*sin( 2*pi/3)*sin(Xs).*cos(Ys).*cos(Zs) ); 
vf = fftn( 2/sqrt(3)*sin(-2*pi/3)*cos(Xs).*sin(Ys).*cos(Zs) ); 
wf = zeros(N,N,N);

clear Xs Ys Zs

for k=1:round(Tend/dt) %%===========================
u=real(ifftn(uf)); 
v=real(ifftn(vf)); 
w=real(ifftn(wf)); 
nltu=Z.*fftn(v.*real(ifftn(IKX.*vf-IKY.*uf ))-w.*real(ifftn(IKZ.*uf-IKX.*wf ))); 
nltv=Z.*fftn(w.*real(ifftn(IKY.*wf-IKZ.*vf ))-u.*real(ifftn(IKX.*vf-IKY.*uf ))); 
nltw=Z.*fftn(u.*real(ifftn(IKZ.*uf-IKX.*wf ))-v.*real(ifftn(IKY.*wf-IKZ.*vf ))); 
pf=-(IKX.*nltu+IKY.*nltv+IKZ.*nltw)./K2p; 
pf(1,1,1)=0; nltu=nltu-IKX.*pf; 
ufa = E2.*(uf + dt/2*nltu);
nltv=nltv-IKY.*pf; 
vfa = E2.*(vf + dt/2*nltv);
nltw=nltw-IKZ.*pf; 
wfa = E2.*(wf + dt/2*nltw);
u=real(ifftn(ufa)); 
v=real(ifftn(vfa)); 
w=real(ifftn(wfa)); 
nltua=Z.*fftn(v.*real(ifftn(IKX.*vfa-IKY.*ufa))-w.*real(ifftn(IKZ.*ufa-IKX.*wfa))); 
nltva=Z.*fftn(w.*real(ifftn(IKY.*wfa-IKZ.*vfa))-u.*real(ifftn(IKX.*vfa-IKY.*ufa))); 
nltwa=Z.*fftn(u.*real(ifftn(IKZ.*ufa-IKX.*wfa))-v.*real(ifftn(IKY.*wfa-IKZ.*vfa))); 
pf=-(IKX.*nltua+IKY.*nltva+IKZ.*nltwa)./K2p; 
pf(1,1,1)=0;


nltua=nltua-IKX.*pf; 
ufb = E2.*(uf + dt/2*nltua); 
nltva=nltva-IKY.*pf; 
vfb = E2.*(vf + dt/2*nltva); 
nltwa=nltwa-IKZ.*pf; 
wfb = E2.*(wf + dt/2*nltwa);

u=real(ifftn(ufb)); 
v=real(ifftn(vfb)); 
w=real(ifftn(wfb)); 

nltub=Z.*fftn(v.*real(ifftn(IKX.*vfb-IKY.*ufb))-w.*real(ifftn(IKZ.*ufb-IKX.*wfb))); 
nltvb=Z.*fftn(w.*real(ifftn(IKY.*wfb-IKZ.*vfb))-u.*real(ifftn(IKX.*vfb-IKY.*ufb))); 
nltwb=Z.*fftn(u.*real(ifftn(IKZ.*ufb-IKX.*wfb))-v.*real(ifftn(IKY.*wfb-IKZ.*vfb))); 
pf=-(IKX.*nltub+IKY.*nltvb+IKZ.*nltwb)./K2p; 
pf(1,1,1)=0; 
nltub=nltub-IKX.*pf; 
ufb = E.*(ufa + dt*nltub);
nltvb=nltvb-IKY.*pf; 
vfb = E.*(vfa + dt*nltvb); 
nltwb=nltwb-IKZ.*pf; 
wfb = E.*(wfa + dt*nltwb);
u=real(ifftn(ufb)); 
v=real(ifftn(vfb)); 
w=real(ifftn(wfb)); 
nltuc=Z.*fftn(v.*real(ifftn(IKX.*vfb-IKY.*ufb))-w.*real(ifftn(IKZ.*ufb-IKX.*wfb))); 
nltvc=Z.*fftn(w.*real(ifftn(IKY.*wfb-IKZ.*vfb))-u.*real(ifftn(IKX.*vfb-IKY.*ufb))); 
nltwc=Z.*fftn(u.*real(ifftn(IKZ.*ufb-IKX.*wfb))-v.*real(ifftn(IKY.*wfb-IKZ.*vfb))); 
pf=-(IKX.*nltuc+IKY.*nltvc+IKZ.*nltwc)./K2p; 
pf(1,1,1)=0;
nltuc=nltuc-IKX.*pf;
nltvc=nltvc-IKY.*pf; 
nltwc=nltwc-IKZ.*pf;
uf = E.*(uf + dt/6*(nltu + 2*(nltua+nltub) + nltuc)); 
vf = E.*(vf + dt/6*(nltv + 2*(nltva+nltvb) + nltvc)); 
wf = E.*(wf + dt/6*(nltw + 2*(nltwa+nltwb) + nltwc));
%% plot result
time(k) = k*dt;
EPS=(2*(real(ifftn(IKX.*uf)).^2+real(ifftn(IKY.*vf)).^2+ ... 
    real(ifftn(IKZ.*wf)).^2)+real(ifftn(IKY.*uf + IKX.*vf)).^2+ ...
    real(ifftn(IKZ.*vf + IKY.*wf)).^2+real(ifftn(IKX.*wf + IKZ.*uf)).^2)/Re;
Diss(k) = mean(mean(mean(EPS)));
X = reshape(EPS,N^3,1);
fd = floor(X./(10.^floor(log10(X)))); % extract the first digit
stat(:,k) = histc(fd,1:9)'/length(X); % compute the probability 
H(k)=-sum(stat(:,k).*log10(stat(:,k)));
div = max(max(max(abs( real(ifftn(IKX.*uf+IKY.*vf+IKZ.*wf)) ))));
disp([' time=',num2str(time(k)),' divergence=',num2str(div),' dissipation=',num2str(Diss(k))]) 
subplot(1,2,1), plot(time,Diss), xlabel t, ylabel dissipation
subplot(1,2,2), plot(time,H,[0 time(k)],SNB*[1 1],'k'), xlabel t, ylabel H
drawnow,
end 
% TIME STEPPING =================================
