clear
clc
% 3º Trabalho de Dinamica dos Fluidos Computacional
% EXERCÍCIO 1, ITEM 1.1 e ITEM 1.2
% Algoritmo para calculo da distribuição de velocidades e temperaturas
% ao longo da coordenada radial de um cilindro de seção anular.
% A) PRE-PROCESSAMENTO---------------------------------------------------%
% A.1-Dados do problema:

% Parâmetros do duto:

RR=4;                       % razão de raios
rint=1;                   % raio interno
rext=RR*rint;               % raio externo
deltax=1;                       % comprimento axial
A=15*pi()*(rint)^2*deltax;      % area da seção transversal
L=rext-rint;                    % comprimento entre as paredes do duto
Pm=2*pi()*(rint+rext);          % perímetro molhado
Paq=2*pi()*rext;                % perímetro aquecido
Dh=4*A/Pm;                      % diametro hidraulico
De=2*rext;                      % diametro externo
Di=2*rint;                      % diametro interno

% Parâmetros físicos:

visc=1;                     % viscosidade
dP=-2;                          % queda de pressão entre os pontos
rho=1000;                       % massa específica da água
k=1;                            % coeficiente de condução
qw=20;                          % fluxo prescrito na parede externa
Tw=300;                         % temperatura da parede externa

% A.1.2-Volume de controle e número de pontos 

nvc=20;             
N=nvc+2;    

% A.1.3-Criação da malha

% Vetor da coordenada das faces
alpha=1; % coeficiente constante que depende da uniformidade da malha
for i=1:N+1
    if (i==1) || (i==2)
        r_face(i)=rint;
    elseif (i>=3) && (i<N+1)
        r_face(i)=rint+L*((i-2)/(N-2))^alpha;
    elseif i==N+1
        r_face(i)=r_face(N);
    end
end
% Vetor da distancia entre as faces
for i=1:N 
    dist_face(i)=r_face(i+1)-r_face(i);
end
% Vetor coordenada dos pontos
for i=1:N
    r_ponto(i)=(((r_face(i+1)+r_face(i)))/2); 
    
end
r_ponto(N)=(r_face(N));
% Vetor distancia entre os pontos internos
for i=1:N-1 
    dist_p(i)=r_ponto(i+1)-r_ponto(i); 
end

% ITEM 1.1: TDMA PARA DISTRIBUIÇÃO DE VELOCIDADES ------------------------%

% Discretização dimensional-cálculo dos coeficientes  

   for i=1:N
       if i ==1
           % condição no contorno, velocidade nula na parede.
           a(i)=1;
           b(i)=0;
           c(i)=0;
           d(i)=0;
       elseif (i>=2) && (i<N)
           b(i)=visc*r_face(i+1)/dist_p(i);
           c(i)=visc*r_face(i)/dist_p(i-1);
           d(i)=-dP*r_ponto(i)*(r_face(i+1)-r_face(i));
           a(i)=b(i)+c(i);
       else
           % condição no contorno, velocidade nula na parede.
           b(i)=0;
           c(i)=0;
           a(i)=1;  
           d(i)=0;
       end
   end 
   
% TDMA 
        p(1)=b(1)/a(1); 
        q(1)=d(1)/a(1);
for i=2:N
        t(i)=a(i)-c(i)*p(i-1);
        p(i)=b(i)/t(i);
        q(i)=(d(i)+c(i)*q(i-1))/t(i);
end
  
% Velocidade nos pontos

u(N)=q(N); 
for i=N-1:-1:1
    u(i)=q(i)+p(i)*u(i+1);    
end

%------------------------FIM TDMA VELOCIDADE----------------------------%
% B) POS-PROCESSAMENTO (ITEM 1.1)----------------------------------------%

% Perfil de velocidade obtida do TDMA

figure (1)
hold on
plot(u,r_ponto, 'b-o'), grid;
title ('Perfil de velocidade')
xlabel ('u')
ylabel ('r')

% Balanço global

Jin=-visc*r_ponto(1)*(u(2)-u(1))/(r_ponto(2)-r_ponto(1));
Jout=-visc*r_ponto(N)*(u(N)-u(N-1))/(r_ponto(N)-r_ponto(N-1));
soma=0;
for i=1:N
    S(i)=-dP*r_ponto(i)*(r_face(i+1)-r_face(i));
    soma=soma+S(i);
end
Bal=Jin-Jout+soma;

% Velocidade média numérica

sum_u=0;
area=0;
for i=1:N-1
   d_area=r_ponto(i)*dist_face(i);
   area=area+d_area;
   sum_u=sum_u+u(i)*d_area;
end
u_medio=sum_u/area;

% Fator de atrito 
f=(-dP)*Dh/(0.5*rho*u_medio^2);
% Numero de Reynolds
Re=rho*u_medio*Dh/visc;
% Produto numero de Reynolds e fator de atrito
Prod=f*Re;

% Velocidade adimensional 

for i=1:N
    ua(i)=(visc*u(i))/((rint^2)*(-dP));  
end

% Velocidade adimensional exata

for i=1:N
    ue(i)=0.25+(0.25*15/log(4))*log(r_ponto(i)/rint)-0.25*(r_ponto(i)/rint)^2;
end

% Perfil de velocidades adimensionais numerica e analitica 

figure (2)
hold on
plot(ua,r_ponto, 'b-o'), grid;
plot(ue,r_ponto, 'g'), grid;
legend('ua=velocidade analítica','ue=velocidade exata')
title ('Comparação do perfil de velocidade numerica e exata para 5 volumes')
xlabel ('u')
ylabel ('r')

% Erros

sum_erro=0;
for i=1:N
    erro(i)=abs(ua(i)-ue(i));   % erro absoluto
    sum_erro=sum_erro+erro(i);
    erro_med=sum_erro/N;        % erro médio
end

% Perfil do erro versus coordenada radial

figure (3)
hold on
plot(r_ponto/rint,erro, 'b-o'), grid
title ('Distribuição dos erros absolutos')
xlabel ('r/rint')
ylabel ('Erro absoluto')

% ITEM 1.2: TDMA PARA DISTRIBUIÇÃO DE TEMPERATURAS----------------------%
% Discretização Adimensional   

   for i=1:N
       if i ==1
% condição de contorno na parede externa com fluxo prescrito
           at(i)=1;
           bt(i)=1;
           ct(i)=0;
           dt(i)=0;
           
       elseif (i>=2) && (i<N)
           bt(i)=k*r_face(i+1)/dist_p(i);
           ct(i)=k*r_face(i)/dist_p(i-1);
           dt(i)=-u(i)/u_medio*qw*Paq/A*(r_ponto(i)*(r_face(i+1)-r_face(i)));
           at(i)=bt(i)+ct(i);
           
       else
% condição de contorno na parede externa com temperatura prescrita           
           bt(i)=0;
           ct(i)=0;
           at(i)=1;  
           dt(i)=Tw;
       end
   end
   
for i=1:N
    if i==1
        pt(i)=bt(i)/at(i); 
        qt(i)=dt(i)/at(i);
        
    else
        tt(i)=at(i)-ct(i)*pt(i-1);
        pt(i)=bt(i)/tt(i);
        qt(i)=(dt(i)+ct(i)*qt(i-1))/tt(i);
        
    end
end
%------------------------------TDMA TEMPERATURA---------------------------%
% Distribuição de temperaturas
T(N)=qt(N);
for i = N-1:-1:1
    T(i)=pt(i)*T(i+1)+qt(i);
end

%-----------------------FIM TDMA TEMPERATURA------------------------------%
% C) POS-PROCESSAMENTO (ITEM 1.2)--------------------------------------%
% Temperatura adimensionalizada versus raio do ponto pelo raio interno
figure (4)
hold on
plot((r_ponto/rint),((Tw-T)/(qw*Dh/k)), 'b-o'), grid
title ('Distribuição de temperatura adimensional')
xlabel ('r/rint')
ylabel ('\theta')

% Temperatura com dimensão versus raio do ponto
figure (5)
hold on
plot(r_ponto,T, 'b-o'), grid
title ('Distribuição de temperatura dimensional')
xlabel ('r')
ylabel ('T')

% Temperatura média de mistura e determinação do número de Nusselt

T_mist=0;
area2=0;
for i=1:N-1
   d_area2=r_ponto(i)*(r_face(i+1)-r_face(i));
   area2=area2+d_area2;
   T_mist=T_mist+T(i)*u(i)*d_area2;         %temperatura de mistura   
end
const=1/(u_medio*area2);
T_medmist=const*T_mist;

% Coeficiente médio convectivo

hmedio=qw/(Tw-T_medmist);

% Numero de Nusselt

Nu=hmedio*(Dh)/k;



