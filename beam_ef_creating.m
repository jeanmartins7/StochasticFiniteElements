function [K,M] = beam_ef_creating(N,rho,E,A,I,Le)
L=Le/N;
% Stiffness Matrix
k=(E*I/(L^3))*[12, 6*L, -12, 6*L; 6*L, 4*L^2, -6*L, 2*L^2; -12, -6*L, 12, -6*L; 6*L, 2*L^2, -6*L, 4*L^2];
k1=(E*I/(L^3))*[24, 0, -12, 6*L; 0, 8*L^2, -6*L, 2*L^2; -12, -6*L, 24, 0; 6*L, 2*L^2, 0, 8*L^2];
% Mass Matrix
m=(rho*A*L/420)*[156, 22*L, 54, -13*L; 22*L, 4*L^2, 13*L, -3*L^2; 54, 13*L, 156, -22*L; -13*L, -3*L^2, -22*L, 4*L^2];
m1=(rho*A*L/420)*[312, 0, 54, -13*L; 0, 8*L^2, 13*L, -3*L^2; 54, 13*L, 312, 0; -13*L, -3*L^2, 0, 8*L^2];
for i=1:4
    for j=1:4
            K(i,j)=k(i,j);
            K(2*N-2+i,2*N-2+j)=k(i,j);
    end
    
end
for n=1:(N-2)
    for i=1:4
          for j=1:4
          K(i+2*n,j+2*n)=k1(i,j);
          end
    end
end
for i=1:4
    for j=1:4
            M(i,j)=m(i,j);
            M(2*N-2+i,2*N-2+j)=m(i,j);
    end
end
for n=1:(N-2)
    for i=1:4
          for j=1:4
          M(i+2*n,j+2*n)=m1(i,j);
          end
    end
end
% Boundary conditions (cantilever beam)
    K(1,:)=[]; 
    K(1,:)=[]; % Second row
    K(:,1)=[];
    K(:,1)=[]; % Second column
    M(1,:)=[];
    M(1,:)=[]; % Second row
    M(:,1)=[];
    M(:,1)=[]; % Second column
    
end