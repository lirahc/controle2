% Representação do Sistema no Espaço de Estados
clc;
clear;
close all;
format long;

C1=100e-9;
C2=680e-9;
R1=36e3;
R2=18e3;

A = [0 1/C1;-1/(R1*R2*C2) -(R1+R2)/(R1*R2*C2)];
B = [0;1/(R1*R2*C2)];
C = [1 0];
D = [0];

% Resposta ao degrau da função original
sistema=ss(A,B,C,D)

figure(1)
step(sistema)

[Y,t,X] = step(sistema);
% Separando as partes
figure(2)
subplot(2,2,1)
plot(t,X(:,1))
hold on
title('Vc1')
xlabel('Tempo (s)')
ylabel('Vc1 (V)')
subplot(2,2,2)
plot(t,X(:,2))
hold on
title('Ic1')
xlabel('Tempo (s)')
ylabel('Ic1 (A)')
subplot(2,2,3)
plot(t,Y)
hold on
title('y')
xlabel('Tempo (s)')
ylabel('y (V)')

% Especificações do projeto
Mp = 0.12
ts5 = 0.018                       

% Calculo do zeta e wn
syms zeta wn
zeta = solve(Mp == exp(-pi*(zeta/sqrt(1-zeta^2))), zeta);
zeta = eval(zeta(1))
wn = 3/(zeta*ts5)

% Polos de malha fechada
s1=-zeta*wn+j*wn*sqrt(1-zeta^2)
s2=-zeta*wn-j*wn*sqrt(1-zeta^2)
s3= 5*real(s1)   

% Matriz de Controlabilidade
M=ctrb(A,B)
% Teste de controlabilidade
rank(M)  

% Projeto do controlador / servosistema
A_chapeu = [A zeros(2,1);-C 0];
B_chapeu = [B;0];

% Matriz de ganhos do controlador K_chapeu = [k1 k2 -ki]
K_chapeu=acker(A_chapeu,B_chapeu,[s1 s2 s3])

% Verificando
eig(A_chapeu-B_chapeu*K_chapeu);

K = [K_chapeu(1) K_chapeu(2)];
Ki = -K_chapeu(3);
AA = [A-B*K B*Ki;-C 0];
BB = [0;0;1];
CC = [1 0 0];
DD = 0;

% Resposta ao degrau do controlador e original
figure(3)
step(sistema)
hold on
step(AA,BB,CC,DD)
title('Resposta ao Degrau da Planta Original e Controlador')
legend('Planta Original', 'Controlador')
hold off

%% observador de ordem mínima
Abb=[A(2,2)];
Aab=[A(1,2)];
L=5*real(s1);
Ke=acker(Abb',Aab',L)'

Aaa=A(1,1);
Aba=A(2,1);
Ba=0;
Bb=B(2,1);

Achapeu=Abb-Ke*Aab;
Bchapeu=Achapeu*Ke+Aba-Ke*Aaa;
Fchapeu=Bb-Ke*Ba;
Cchapeu=[0;1];
Dchapeu=[1;Ke];

AAA = [A zeros(2,1); Bchapeu*C Achapeu]
BBB = [B;Fchapeu]
CCC = [Dchapeu*C Cchapeu]
DDD = zeros(length(A),1)

figure(4)
hold on
t=0:0.00001:0.12;
obs = step(AAA,BBB,CCC,DDD,1,t);
plot(t,obs(:,1),'ko')
sori = step(A,B,C,D,1,t);
plot(t,sori,'r','LineWidth',2)
title('Resposta ao Degrau da Planta Original e Observador de Ordem Minima')
legend('Observador','Original');
hold off

%% Equação recursiva
T = 1e-4;
k = 0:1:(0.30/T); 
r = ones(1,length(k)/2); 
r = [r 1.5.*r];

% Condições iniciais
x1(1)=0;
x2(1)=0;
x1_ponto(1)=0;
x2_ponto(1)=0;
y(1)=C(1)*x1(1)+C(2)*x2(1);
u(1) = 0;          
eta1(1)=0;
eta1_ponto(1)=Fchapeu(1)*u(1);   
x1til(1) = 0; 
x2til(1) = 0; 
E(1) = 0;            
y_til(1) = 0;  
e(1) = 0;         
v(1) = 0;        
E_ponto(1) = r(1)-y(1);

for j=2:length(k)

    x1(j)= T*x1_ponto(j-1)+x1(j-1);
    x2(j)= T*x2_ponto(j-1)+x2(j-1);    
    E(j) = T*E_ponto(j-1)+E(j-1);    
    v(j) = K(1)*x1(j)+K(2)*x2(j);
    u(j) = Ki*E(j)-v(j);
    x1_ponto(j)=A(1,1)*x1(j)+A(1,2)*x2(j)+B(1)*u(j);  
    x2_ponto(j)=A(2,1)*x1(j)+A(2,2)*x2(j)+B(2)*u(j); 
    y(j)=C(1)*x1(j)+C(2)*x2(j)+D*u(j);
    E_ponto(j) = r(j)-y(j);
    eta1(j)=T*eta1_ponto(j-1)+eta1(j-1);
    eta1_ponto(j)=Achapeu(1,1)*eta1(j)+Bchapeu(1)*y(j)+Fchapeu(1)*u(j);
    xtil1(j)=Cchapeu(1,1)*eta1(j)+Dchapeu(1)*y(j);
    xtil2(j)=Cchapeu(2,1)*eta1(j)+Dchapeu(2)*y(j);    
    y_til(j)=C(1)*xtil1(j);
    e(j) = y(j)-y_til(j);
end

% Variáveis originais e do observador
figure(5)
subplot(2,1,1)
plot(k*T,x1,'*')
title('x1(Vc1)')
hold on
subplot(2,1,2)
plot(k*T,x2,'*')
title('x2(Ic1)')
hold on
subplot(2,1,1)
plot(k*T,xtil1,'o')
subplot(2,1,2)
plot(k*T,xtil2,'o')

hold off

% Saida do Sistema
figure(6)
sist=ss(AA,BB,CC,DD)
opt = stepDataOptions('InputOffset',1,'StepAmplitude',.5);
t=0.0:T:0.15;
sist_=step(sist,t,opt);
plot(t+.15,sist_(:,1))
hold on
plot(k*T,y,'r*')
title('Resposta ao Degrau da Planta Controlada')
legend('Planta Original', 'Planta Controlada')
xlabel('Tempo (s)')
xlim([.15 .30])
hold off

% Ação de controle
figure(7);
plot(k*T, u, k*T, y);
%plot(k*T, y);
title('Ação de Controle do Sistema com Observador de Ordem Mínima');
xlabel('Tempo (s)');
ylabel('Ação de Controle');
title('Forma de onda do sinal de saída')
legend('Ação de Controle', 'Controlador')

%xlim([0.15 0.30])
grid on;
