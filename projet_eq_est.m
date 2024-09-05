% Thiago Santos de Lira

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

[NUM,DEN] = ss2tf(A,B,C,D)

% Resposta ao degrau da função original
sistema=ss(A,B,C,D);

[Y,t,X] = step(sistema);
% Separando as partes
figure()
subplot(2,2,1)
plot(t,X(:,1))
hold on
title('x1')
grid on;
subplot(2,2,2)
plot(t,X(:,2))
hold on
title('x2')
grid on;
subplot(2,2,3)
plot(t,Y)
hold on
title('y')
grid on;

clear u
T=0.0001;
k=0:.12/T;

u=ones(1,length(k));

% Condições Iniciais
x1(1)=0;  %  para k = 0
x2(1)=0;  %  para k = 0
x1_ponto(1)=0;  %  para k = 0
x2_ponto(1)=0;  %  para k = 0
y(1)=0;  %  para k = 0

for j=2:length(k)
    % Equações dos integradores
    x1(j)=T*x1_ponto(j-1)+x1(j-1);
    x2(j)=T*x2_ponto(j-1)+x2(j-1);    
      
    x1_ponto(j)=A(1,1)*x1(j)+A(1,2)*x2(j);  
    x2_ponto(j)=A(2,1)*x1(j)+A(2,2)*x2(j)+B(2)*u(j);  
    
    y(j)=C(1)*x1(j);
end

subplot(2,2,1)
plot(k*T, x1, 'yx')  % 'y' define a cor amarelo
ylabel('x1')
grid on;

subplot(2,2,2)
plot(k*T, x2, 'bx')  % 'b' define a cor azul
ylabel('x2')
grid on;

subplot(2,2,3)
plot(k*T, y, 'gx')   % 'm' define a cor verde
ylabel('y')
grid on;
hold off


% Informacoes de Projeto
Mp = 0.12
ts5 = 0.018                     

% Calculo do zeta
syms zeta
zeta = solve(Mp == exp(-pi*(zeta/sqrt(1-zeta^2))), zeta);
zeta = eval(zeta(1))

% Calculo da frequencia natutral nao amortecida (Wn), dos polos dominantes
% de malha fechada.
syms wn
wn = 3/(zeta*ts5)

% Calculo da frequencia natural amortecida (Wd).
wd = wn*sqrt(1-zeta^2)

s1=-zeta*wn+i*wn*sqrt(1-zeta^2)
s2=-zeta*wn-i*wn*sqrt(1-zeta^2)
s3=5*real(s1)   % Autovalor adicional para o projeto, valor mais alto que a parte real de s1 para maior estabilidade.

% Autovalores (polos de malha fechada desejados)
u1=s1;
u2=s2;
u3=s3;


% Matriz de Controlabilidade
M=ctrb(A,B)
% Teste de controlabilidade
rank(M)  

% Projeto do controlador / servosistema
A_chapeu = [A zeros(2,1);-C 0];
B_chapeu = [B;0];

% Matriz de ganhos do controlador K_chapeu = [k1 k2 -ki]
K_chapeu=acker(A_chapeu,B_chapeu,[u1 u2 u3])

% Verificando
eig(A_chapeu-B_chapeu*K_chapeu)

K = [K_chapeu(1) K_chapeu(2)];
Ki = -K_chapeu(3);
AA = [A-B*K B*Ki;-C 0];
BB = [0;0;1];
CC = [1 0 0];
DD = 0;

% Resposta ao degrau do controlador e original
figure()
step(sistema)
hold on
step(AA,BB,CC,DD)
title('Resposta ao Degrau da Planta Original X Controlador')
legend('Planta Original', 'Controlador')
hold off


% observador de ordem mínima
Abb=[A(2,2)];
Aab=[A(1,2)];
L=[s3];
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

figure()
t=0:0.00001:0.18;
AAA = [A zeros(2,1); Bchapeu*C Achapeu]
BBB = [B;Fchapeu]
CCC = [Dchapeu*C Cchapeu]
DDD = zeros(length(A),1)

step(sistema)
hold on
so = step(AAA,BBB,CCC,DDD,1,t);
plot(t,so(:,1),'m')
title('Resposta ao Degrau da Planta Original X Observador de Ordem minima')
legend('Planta Original', 'Observador de Ordem minima')
hold off



figure()
ax = gca
ax.YLim = [0 1.4];
hold on
so = step(AAA,BBB,CCC,DDD,1,t);
plot(t,so(:,1),'ko')
step(AA,BB,CC,DD,1,t);
sori = step(A,B,C,D,1,t);
plot(t,sori,'r','LineWidth',1)
title('Resposta ao Degrau da Planta Original X Observador X Controlador')
legend('Observador','Controlado','Original');
hold off

T = 1e-4;
k = 0:1:(0.30/T); 
r = ones(1,length(k)/2); 
r = [r 1.5.*r];

% condições iniciais
x1(1)=0;
x2(1)=0;

x1_ponto(1)=0;
x2_ponto(1)=0;

y(1)=C(1)*x1(1)+C(2)*x2(1);

eta1(1)=0;
eta1_ponto(1)=Fchapeu(1)*u(1);
   
x1til(1) = 0; 
x2til(1) = 0; 
E(1) = 0;            % Saida do integrador 
y_til(1) = 0;  
e(1) = 0;          % Erro
v(1) = 0; 
u(1) = 0;          % Sinal de controle          
E_ponto(1) = r(1)-y(1);

for j=2:length(k)
% sistema original
    % Equações dos integradores
    x1(j)=T*x1_ponto(j-1)+x1(j-1);
    x2(j)=T*x2_ponto(j-1)+x2(j-1);
    
    E(j) =T*E_ponto(j-1)+E(j-1);
    
    v(j) = K(1)*x1(j)+K(2)*x2(j);
    u(j) = Ki*E(j)-v(j);
    
    % Equação diferencial de estados:  Xponto=A*X+B*U
    x1_ponto(j)=A(1,1)*x1(j)+A(1,2)*x2(j)+B(1)*u(j);  
    x2_ponto(j)=A(2,1)*x1(j)+A(2,2)*x2(j)+B(2)*u(j); 

    % Equação de Saída: Y=C*X+D*U
    y(j)=C(1)*x1(j)+C(2)*x2(j)+D*u(j);
    E_ponto(j) = r(j)-y(j);
    % observador de ordem mínima
    % Equações dos integradores do observador
    eta1(j)=T*eta1_ponto(j-1)+eta1(j-1);
    % Equação diferencial de estados do observador
    eta1_ponto(j)=Achapeu(1,1)*eta1(j)+Bchapeu(1)*y(j)+Fchapeu(1)*u(j);
    % equação de saída do observador de ordem mínima
    xtil1(j)=Cchapeu(1,1)*eta1(j)+Dchapeu(1)*y(j);
    xtil2(j)=Cchapeu(2,1)*eta1(j)+Dchapeu(2)*y(j);
    
    y_til(j)=C(1)*xtil1(j);
    e(j) = y(j)-y_til(j);
end
figure()
% plotando as variáveis originais
subplot(2,2,1)
plot(k*T, x1, 'o', 'LineWidth', 0.5, 'Color', 'b')  % 'b' para azul
hold on
subplot(2,2,2)
plot(k*T, x2, 'o', 'LineWidth', 0.5, 'Color', 'b')
hold on
subplot(2,2,3)
plot(k*T, y, 'o', 'LineWidth', 0.5, 'Color', 'b')
hold on
% plotando as variáveis do observador
subplot(2,2,1)
plot(k*T, xtil1, 'o', 'LineWidth', 0.5, 'Color', 'r')  % 'r' para red
subplot(2,2,2)
plot(k*T, xtil2, 'o', 'LineWidth', 0.5, 'Color', 'r')
hold off

figure()
step(sistema)
hold on
step(AA,BB,CC,DD)
plot(k*T,y, 'o', 'LineWidth', 0.5)
title('Resposta ao Degrau da Planta Original X Controlador')
legend('Planta Original', 'Controlador')
hold off

% Plot da ação de controle
figure();
plot(k*T, u, 'LineWidth', 2);
title('Ação de Controle do Sistema com Observador de Ordem Mínima');
xlabel('Tempo (s)');
ylabel('Ação de Controle');
grid on;