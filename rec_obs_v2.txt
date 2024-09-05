% Exemplo de Projeto de um observador de estados de ordem mínima

clc
clear
clf

C1=100e-9;
C2=680e-9;
R1=36e3;
R2=18e3;

A = [0 1/C1;-1/(R1*R2*C2) -(R1+R2)/(R1*R2*C2)];
B = [0;1/(R1*R2*C2)];
C = [1 0];
D = [0];

sistema=ss(A,B,C,D);

% resposta ao degrau
[Y,t,X] = step(sistema);
% separando as saídas
figure(1)
subplot(2,2,1)
plot(t,X(:,1))
title('V1')
hold on
subplot(2,2,2)
plot(t,X(:,2))
title('Ic')
hold on
hold on
subplot(2,2,3)
plot(t,Y)
title('y')
hold on

% observador de ordem mínima
Abb = (-1/(R1*C2))+(-1/(R2*C2));
Aab = -1/(R1*R2*C2);
L   = -10;
Ke  = acker(Abb',Aab',L)'


Aaa=0;
Aba=1/C1;
Ba=0;
Bb=1/(R1*R2*C2);

Achapeu=Abb-Ke*Aab;
Bchapeu=Achapeu*Ke+Aba-Ke*Aaa;
Fchapeu=Bb-Ke*Ba;

Cchapeu=[0;1];
Dchapeu=[1;Ke];

% simulação com equações recursivas
T=0.0001;
k=0:1:round(.12/T);
u=ones(1,length(k));
% condições iniciais
x1(1)=0;
x2(1)=0;

x1_ponto(1)=0;
x2_ponto(1)=0;

y(1)=C(1)*x1(1)+C(2)*x2(1)

eta1(1)=0;

eta1_ponto(1)=Fchapeu(1)*u(1);

%%
for j=2:length(k)
% sistema original
    % Equações dos integradores
    x1(j)=T*x1_ponto(j-1)+x1(j-1);
    x2(j)=T*x2_ponto(j-1)+x2(j-1);    
  
    % Equação diferencial de estados:  Xponto=A*X+B*U
    x1_ponto(j)=A(1,1)*x1(j)+A(1,2)*x2(j)+B(1)*u(j);  
    x2_ponto(j)=A(2,1)*x1(j)+A(2,2)*x2(j)+B(2)*u(j); 
      
    % Equação de Saída: Y=C*X+D*U
    y(j)=C(1)*x1(j)+C(2)*x2(j)+D*u(j);
    
    % observador de ordem mínima
    % Equações dos integradores do observador
    eta1(j)=T*eta1_ponto(j-1)+eta1(j-1);
    
    % Equação diferencial de estados do observador
    eta1_ponto(j)=Achapeu*eta1(j)+Bchapeu*y(j)+Fchapeu*u(j);
    
    % equação de saída do observador de ordem mínima
    xtil1(j)=Cchapeu(1)*eta1(j)+Dchapeu(1)*y(j);
    xtil2(j)=Cchapeu(2)*eta1(j)+Dchapeu(2)*y(j);
    
end

figure(2)
% plotando as variáveis originais
subplot(2,2,1)
plot(k*T,x1,'-')
subplot(2,2,2)
plot(k*T,x2,'-')
subplot(2,2,3)
plot(k*T,y,'-')

% plotando as variáveis do observador
subplot(2,2,1)
plot(k*T,xtil1,'-')
subplot(2,2,2)
plot(k*T,xtil2,'-')






hold off