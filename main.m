%Start script
%%
%modeling Euler-Bernoulli beam in MFE

%Parametros
N = 5;%Nelementos
E = 7e10;%modulo de young
Es=E+5*E/100.*randn(100,1);
Le = 5; %comprimento da barra
le = Le/N; %comprimento do elemento
rho = 2700;%densidade da barra
h = 0.5;%comprimento
b = 0.5;%largura
Ae = b*h;%area
I = (b*h*(b^2 + h^2))/24;%inercia

%Obtendo as matrizez de massa e rigidez
[K3,M3] = beam_ef_creating(N,rho,E,Ae,I,Le);
%amortecimento proporcional
alfa=1e-2;
beta=1e-5;
C3=alfa*M3 + beta*K3;
[m,n]=size(M3);
%% random field parameters
%% initial
d_KL  = 4;           % number of terms in the expansion
ns = 50;
mu    = 1;           % mean of the field
sigma = 1;         % standard deviation of the field
l_c   = 1;           % correlation length x
% domain of definition
a    = 0.5;
b    = 1;
D    = b-a;
n    = 3e2;
xnod = linspace(-1/2,1/2,n);
uo = 0.001;

cov_kernel = @(x1,x2) sigma^2* exp(-abs(x1-x2)/l_c); %função de covariança em simbólico
[eigval, eigvec, eigfun, wn,alphaVec] = KL_analytical(d_KL, l_c, sigma, xnod, 'true'); %função de covariança em numérico obtém { autovalores, autovetores, wn e autofunção}

%%
%monte carlo

for p = 1: ns
    
    [K_s] = beam_ef_creating(N,rho,Es(p),Ae,I,Le); %cria a matriz de rigidez em relação ao parâmetro incerto 
    K_ie = zeros(m);
    for i = 1:d_KL
        K_ie = K_ie + sqrt(alphaVec(i).*eigval(i))*K_s; %matriz estocástica
    end
    Kt = K3 + K_ie;  %soma das matrizes, matriz randômica e matriz elementar  
    
    [U,t]=time_response(M3,Kt,C3,1, uo);%Deflexão na extremidade da viga
    
    figure(3)
    plot(t,U(:,end-1));
    xlabel('Tempo[s]','Interpreter','latex','fontsize',16)
    ylabel('$u_{n}$[m]','Interpreter','latex','fontsize',16)
    hold on
end





