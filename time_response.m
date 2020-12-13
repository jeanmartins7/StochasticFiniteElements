function [U,t]=time_response(M,K,C,tf, uo)

U0=[0 0 0 0 0 0 0 0 0 0 0 0 -0.0000 0 0.000 0.000 0.000 0.000 uo  0.000];

np=1024;
gdl=length(M);
A = [zeros(gdl) eye(gdl);-inv(M)*K -inv(M)*C];
B = [zeros(gdl);inv(M)];
Ct=eye(2*(gdl));

%Aquisicao dos sinais de deslocamento
sys1=ss(A,B,Ct,0);
Ua=zeros(np+1,(gdl));

tf=0.5;
T=[0:tf/np:tf];
[U,t]=lsim(sys1,Ua,T,U0);
end
