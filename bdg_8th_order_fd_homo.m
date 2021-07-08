Nx = 128
Xmax = 10

dx = 2 * Xmax / Nx;

x = -Xmax + [0:Nx-1]*dx;

D2x = -1/560*(diag(ones(1,Nx-4),4) + diag(ones(1,Nx-4),-4) + diag(ones(1,4),Nx-4) + diag(ones(1,4),-(Nx-4)))...
    + 8/315*(diag(ones(1,Nx-3),3) + diag(ones(1,Nx-3),-3)  + diag(ones(1,3),Nx-3) + diag(ones(1,3),-(Nx-3)))...
    + -1/5*(diag(ones(1,Nx-2),2) + diag(ones(1,Nx-2),-2)   + diag(ones(1,2),Nx-2) + diag(ones(1,2),-(Nx-2)))...
    + 8/5*(diag(ones(1,Nx-1),1) + diag(ones(1,Nx-1),-1)    + diag(ones(1,1),Nx-1) + diag(ones(1,1),-(Nx-1)) )...
    + -205/72*(diag(ones(1,Nx-0),0));
% D2x(1,end-4) = -1/560;
% D2x(1,end-3) = 8/315;
% D2x(1,end-2) = -1/5;
% D2x(1,end-1) = -205/72;
% D2x(end,4) = -1/560;
% D2x(end,33) = 8/315;
% D2x(end,2) = -1/5;
% D2x(end,1) = -205/72;

D2x = D2x/dx^2;

f = x;%ones(size(x));%exp(-x.^2);

dfdx2= D2x*f(:);

figure, plot(x,f,x,dfdx2)

%% GPE Basis
g = 1 
u = ones(size(x));

H0  = -0.5 * D2x;
Hgp = H0 + g * diag(abs(u(:)).^2);
[V0,D0] = eig(Hgp);
D0 = diag(D0); D0(1)

%%
Vgpe = V0(:,2:end);
L_bdg = H0 + 2 * g * diag(abs(u(:)).^2) - 1 * diag(ones(size(x)));
L_bdg = ctranspose(Vgpe) * L_bdg * Vgpe;
H_bdg_uv = g*diag(u(:).^2);
H_bdg_uv = ctranspose(Vgpe) * H_bdg_uv * Vgpe;
H_bdg = [L_bdg -H_bdg_uv; -ctranspose(H_bdg_uv)   L_bdg  ];
tz   = (diag([ones(1,size(L_bdg,1)) ones(1,size(L_bdg,1))*-1]));

[V,D] = eig(H_bdg,full(tz));
D = diag(D);




