clear all;
clc;

%============= Lista 5 de Dinamica dos Fluidos Computacional=============%

% Aluna: Carolina Maria Nunes Bezerra

% P2: Obtenha a solução numérica da equação de conservação genérica
...em regime permanente, com rho, u e gamma constantes utilizando os seguintes
...esquemas:
% a) upwind b) power-law, c) diferença central d) TVD-Van Leer

%============================ PRE-PROCESSAMENTO =========================%

% Dados do problema------------------------------------------------------%
L=1;      % comprimento em x
Area=1;   % area de face
rho=1;    % massa específica do fluído 
phi0=0;   % condicao de contorno para phi em x=0;
phi1=1;   % condicao de contorno para phi em x=L;
P1=10;    % numero de Peclet
P2=500;   % numero de Peclet
Sp=0;     % fonte nula devida S ser constante
% Considerando Peclet P=10 e fonte adimensional S=0:
u_p1=10;                        % velocidade do fluido
S1=0;
Gamma=rho*u_p1*L/P1;            % coeficiente de difusao 
Sc1=S1*Gamma*(phi1-phi0)/L^2;   % termo fonte

% Considerando Peclet P=500 e fonte adimensional S=10:
u_p2=500;                       % velocidade do fluido
S2=10;
Gamma=rho*u_p2*L/P2;            % coeficiente de difusao 
Sc2=S2*Gamma*(phi1-phi0)/L^2;   % termo fonte

% --------------------------- MÉTODO B ----------------------------------%
% Volume de controle e numero de pontos:
nvc=10;             
N=nvc+2;    

% Geração da malha-------------------------------------------------------%

% Vetor da coordenada das faces
alpha=1; 
x_face(1)=0;
for i=2:N
        x_face(i)=L*((i-2)/(N-2))^alpha;   
end

% vetor distancia entre as faces
x_face(N+1)=x_face(N);
for i=1:N 
    dist_face(i)=x_face(i+1)-x_face(i);
end

% vetor coordenada dos pontos
for i=1:N-1
    x_ponto(i)=((x_face(i+1)+x_face(i)))/2; 
end
x_ponto(N)=x_face(N);

% vetor distancia entre os pontos
for i=1:N-1 
    dist_p(i)=x_ponto(i+1)-x_ponto(i); 
end
%=======================================================================%
%============================== PROCESSAMENTO ==========================%
%=======================================================================%

% primeira condicao P=10
  
  Fe1=rho*u_p1*Area;
  Fw1=rho*u_p1*Area;
  
% segunda condicao P=500

  Fe2=rho*u_p2*Area;
  Fw2=rho*u_p2*Area;

%================== ESQUEMA UPWIND 1 (Peclet=10 e Sc=0)=================%
%-------------------------- Calculo dos Coeficientes -------------------%      
% funcao A(|P|) para esquema upwind
  A_upe=1;     
  A_upw=1;
     
% condicao de contorno i=1
  ap_up1(1)=1;
  ae_up1(1)=0;
  aw_up1(1)=0;
  d_up1(1)=phi0;
  
% pontos internos
for i=2:N-1
  
  % coeficiente D é igual para todos os esquemas. Não irá ser repetido
  % dentro dos demais loops.
  De(i)=Gamma*Area/(dist_p(i));
  Dw(i)=Gamma*Area/(dist_p(i-1));
           
  ae_up1(i)=De(i)*A_upe+max(-Fe1,0);
  aw_up1(i)=Dw(i)*A_upw+max(Fw1,0);
  ap_up1(i)=ae_up1(i)+aw_up1(i)-Sp*Area*dist_face(i);
           
  d_up1(i)=Sc1*Area*dist_face(i);
  
end
  
  % condicao de contorno i=N
  ae_up1(N)=0;
  aw_up1(N)=0;
  ap_up1(N)=1;
  d_up1(N)=phi1;

  for i=1:N
      a_up1(i)=ap_up1(i);
      b_up1(i)=ae_up1(i);
      c_up1(i)=aw_up1(i);
  end
  
%----------------------------- TDMA UP1 ---------------------------------%
% Calculo de phi nos pontos

p_up1(1)=b_up1(1)/a_up1(1); 
q_up1(1)=d_up1(1)/a_up1(1);

for i=2:N 
  t_up1(i)=a_up1(i)-c_up1(i)*p_up1(i-1);
  p_up1(i)=b_up1(i)/t_up1(i);
  q_up1(i)=(d_up1(i)+c_up1(i)*q_up1(i-1))/t_up1(i);
end 
      
phi_up1(N)=q_up1(N);
for i=N-1:-1:1
  phi_up1(i)=q_up1(i)+p_up1(i)*phi_up1(i+1);           
end
%----------------------------- FIM TDMA UP1 -----------------------------%

%================== ESQUEMA UPWIND 2 (Peclet=10 e Sc=10)=================%
%-------------------------- Calculo dos Coeficientes -------------------%            

% condicao de contorno i=1
  ap_up2(1)=1;
  ae_up2(1)=0;
  aw_up2(1)=0;
  d_up2(1)=phi0;
  
% pontos internos
for i=2:N-1
                 
  ae_up2(i)=De(i)*A_upe+max(-Fe1,0);
  aw_up2(i)=Dw(i)*A_upw+max(Fw1,0);
  ap_up2(i)=ae_up2(i)+aw_up2(i)-Sp*Area*dist_face(i);
           
  d_up2(i)=Sc2*Area*dist_face(i);
  
end
  
  % condicao de contorno i=N
  ae_up2(N)=0;
  aw_up2(N)=0;
  ap_up2(N)=1;
  d_up2(N)=phi1;

  for i=1:N
      a_up2(i)=ap_up2(i);
      b_up2(i)=ae_up2(i);
      c_up2(i)=aw_up2(i);
  end
  
%----------------------------- TDMA UP2 ---------------------------------%
% Calculo de phi nos pontos

p_up2(1)=b_up2(1)/a_up2(1); 
q_up2(1)=d_up2(1)/a_up2(1);

for i=2:N 
  t_up2(i)=a_up2(i)-c_up2(i)*p_up2(i-1);
  p_up2(i)=b_up2(i)/t_up2(i);
  q_up2(i)=(d_up2(i)+c_up2(i)*q_up2(i-1))/t_up2(i);
end 
      
phi_up2(N)=q_up2(N);
for i=N-1:-1:1
  phi_up2(i)=q_up2(i)+p_up2(i)*phi_up2(i+1);           
end
%----------------------------- FIM TDMA UP2 -----------------------------%

%================= ESQUEMA UPWIND 3 (Peclet=500 e S=0) =================%

%-------------------------- Calculo dos Coeficientes---------------------%      

% condicao de contorno i=1
  ap_up3(1)=1;
  ae_up3(1)=0;
  aw_up3(1)=0;
  d_up3(1)=phi0;
  
% pontos internos
for i=2:N-1
           
  ae_up3(i)=De(i)*A_upe+max(-Fe2,0);
  aw_up3(i)=Dw(i)*A_upw+max(Fw2,0);
  ap_up3(i)=ae_up3(i)+aw_up3(i)-Sp*Area*dist_face(i);
           
  d_up3(i)=Sc1*Area*dist_face(i);
  
end
  
% condicao de contorno i=N
  ae_up3(N)=0;
  aw_up3(N)=0;
  ap_up3(N)=1;
  d_up3(N)=phi1;  

  for i=1:N
      a_up3(i)=ap_up3(i);
      b_up3(i)=ae_up3(i);
      c_up3(i)=aw_up3(i);
  end
  
%----------------------------- TDMA UP3 -----------------------------------%
% Calculo de phi nos pontos

p_up3(1)=b_up3(1)/a_up3(1); 
q_up3(1)=d_up3(1)/a_up3(1);

for i=2:N 
  t_up3(i)=a_up3(i)-c_up3(i)*p_up3(i-1);
  p_up3(i)=b_up3(i)/t_up3(i);
  q_up3(i)=(d_up3(i)+c_up3(i)*q_up3(i-1))/t_up3(i);
end 
      
phi_up3(N)=q_up3(N);
for i=N-1:-1:1
  phi_up3(i)=q_up3(i)+p_up3(i)*phi_up3(i+1);           
end
%----------------------------- FIM TDMA UP3 -----------------------------%

%================= ESQUEMA UPWIND 4 (Peclet=500 e S=10) =================%

%-------------------------- Calculo dos Coeficientes---------------------%      
      
% condicao de contorno i=1
  ap_up4(1)=1;
  ae_up4(1)=0;
  aw_up4(1)=0;
  d_up4(1)=phi0;
  
% pontos internos
for i=2:N-1
           
  ae_up4(i)=De(i)*A_upe+max(-Fe2,0);
  aw_up4(i)=Dw(i)*A_upw+max(Fw2,0);
  ap_up4(i)=ae_up4(i)+aw_up4(i)-Sp*Area*dist_face(i);
           
  d_up4(i)=Sc2*Area*dist_face(i);
  
end
  
% condicao de contorno i=N
  ae_up4(N)=0;
  aw_up4(N)=0;
  ap_up4(N)=1;
  d_up4(N)=phi1;  

  for i=1:N
      a_up4(i)=ap_up4(i);
      b_up4(i)=ae_up4(i);
      c_up4(i)=aw_up4(i);
  end
  
%----------------------------- TDMA UP4 -----------------------------------%
% Calculo de phi nos pontos

p_up4(1)=b_up4(1)/a_up4(1); 
q_up4(1)=d_up4(1)/a_up4(1);

for i=2:N 
  t_up4(i)=a_up4(i)-c_up4(i)*p_up4(i-1);
  p_up4(i)=b_up4(i)/t_up4(i);
  q_up4(i)=(d_up4(i)+c_up4(i)*q_up4(i-1))/t_up4(i);
end 
      
phi_up4(N)=q_up4(N);
for i=N-1:-1:1
  phi_up4(i)=q_up4(i)+p_up4(i)*phi_up4(i+1);           
end
%----------------------------- FIM TDMA UP4 -----------------------------%

%===================== POWER-LAW 1 (Peclet=10 e Sc=0)====================%

%-------------------------- Calculo dos Coeficientes---------------------%
      
% condicao de contorno i=1
  ap_pl1(1)=1;
  ae_pl1(1)=0;
  aw_pl1(1)=0;
  d_pl1(1)=phi0;
  
% pontos internos
for i=2:N-1
  % o coeficiente de Peclet1 servirá para os demais esquemas.             
  Pe1(i)=Fe1/De(i);
  Pw1(i)=Fw1/Dw(i);
  
  Ae_pl1(i)=max(0,(1-0.5*abs(Pe1(i)))^5); 
  Aw_pl1(i)=max(0,(1-0.5*abs(Pw1(i)))^5);
           
  ae_pl1(i)=De(i)*Ae_pl1(i)+max(-Fe1,0);
  aw_pl1(i)=Dw(i)*Aw_pl1(i)+max(Fw1,0);
  ap_pl1(i)=ae_pl1(i)+aw_pl1(i)-Sp*Area*dist_face(i);
           
  d_pl1(i)=Sc1*Area*dist_face(i);
  
end
  
% condicao de contorno i=N
  ap_pl1(N)=1;
  ae_pl1(N)=0;
  aw_pl1(N)=0;
  d_pl1(N)=phi1;  

  for i=1:N
      a_pl1(i)=ap_pl1(i);
      b_pl1(i)=ae_pl1(i);
      c_pl1(i)=aw_pl1(i);
  end
 
%----------------------------- TDMA PL1 -----------------------------------%
% Calculo de phi_pl nos pontos

p_pl1(1)=b_pl1(1)/a_pl1(1); 
q_pl1(1)=d_pl1(1)/a_pl1(1);

for i=2:N 
  t_pl1(i)=a_pl1(i)-c_pl1(i)*p_pl1(i-1);
  p_pl1(i)=b_pl1(i)/t_pl1(i);
  q_pl1(i)=(d_pl1(i)+c_pl1(i)*q_pl1(i-1))/t_pl1(i);
end 
      
phi_pl1(N)=q_pl1(N);
for i=N-1:-1:1
  phi_pl1(i)=q_pl1(i)+p_pl1(i)*phi_pl1(i+1);           
end
%----------------------------- FIM TDMA PL1 -----------------------------%

%===================== POWER-LAW 2 (Peclet=10 e Sc=10)===================%

%-------------------------- Calculo dos Coeficientes---------------------%
    
% condicao de contorno i=1
  ap_pl2(1)=1;
  ae_pl2(1)=0;
  aw_pl2(1)=0;
  d_pl2(1)=phi0;
  
% pontos internos
for i=2:N-1
                    
  ae_pl2(i)=De(i)*Ae_pl1(i)+max(-Fe1,0);
  aw_pl2(i)=Dw(i)*Aw_pl1(i)+max(Fw1,0);
  ap_pl2(i)=ae_pl2(i)+aw_pl2(i)-Sp*Area*dist_face(i);
           
  d_pl2(i)=Sc2*Area*dist_face(i);
  
end
  
% condicao de contorno i=N
  ap_pl2(N)=1;
  ae_pl2(N)=0;
  aw_pl2(N)=0;
  d_pl2(N)=phi1;  

  for i=1:N
      a_pl2(i)=ap_pl2(i);
      b_pl2(i)=ae_pl2(i);
      c_pl2(i)=aw_pl2(i);
  end
 
%----------------------------- TDMA PL2 -----------------------------------%
% Calculo de phi_pl nos pontos

p_pl2(1)=b_pl2(1)/a_pl2(1); 
q_pl2(1)=d_pl2(1)/a_pl2(1);

for i=2:N 
  t_pl2(i)=a_pl2(i)-c_pl2(i)*p_pl2(i-1);
  p_pl2(i)=b_pl2(i)/t_pl2(i);
  q_pl2(i)=(d_pl2(i)+c_pl2(i)*q_pl2(i-1))/t_pl2(i);
end 
      
phi_pl2(N)=q_pl2(N);
for i=N-1:-1:1
  phi_pl2(i)=q_pl2(i)+p_pl2(i)*phi_pl2(i+1);           
end
%----------------------------- FIM TDMA PL2 -----------------------------%

%===================== POWER-LAW 3 (Peclet=500 e Sc=0)==================%

%-------------------------- Calculo dos Coeficientes---------------------%
 
% condicao de contorno i=1
  ap_pl3(1)=1;
  ae_pl3(1)=0;
  aw_pl3(1)=0;
  d_pl3(1)=phi0;
  
% pontos internos
for i=2:N-1
  % o coeficiente de Peclet2 servirá para os demais esquemas.           
  Pe2(i)=Fe2/De(i);
  Pw2(i)=Fw2/Dw(i);
  
  Ae_pl2(i)=max(0,(1-0.5*abs(Pe2(i)))^5);
  Aw_pl2(i)=max(0,(1-0.5*abs(Pw2(i)))^5);
           
  ae_pl3(i)=De(i)*Ae_pl2(i)+max(-Fe2,0);
  aw_pl3(i)=Dw(i)*Aw_pl2(i)+max(Fw2,0);
  ap_pl3(i)=ae_pl3(i)+aw_pl3(i)-Sp*Area*dist_face(i);
           
  d_pl3(i)=Sc1*Area*dist_face(i);
  
end
  
% condicao de contorno i=N
  ap_pl3(N)=1;
  ae_pl3(N)=0;
  aw_pl3(N)=0;
  d_pl3(N)=phi1;  

  for i=1:N
      a_pl3(i)=ap_pl3(i);
      b_pl3(i)=ae_pl3(i);
      c_pl3(i)=aw_pl3(i);
  end
 
%----------------------------- TDMA PL3 ---------------------------------%
% Calculo de phi_pl3 nos pontos

p_pl3(1)=b_pl3(1)/a_pl3(1); 
q_pl3(1)=d_pl3(1)/a_pl3(1);

for i=2:N 
  t_pl3(i)=a_pl3(i)-c_pl3(i)*p_pl3(i-1);
  p_pl3(i)=b_pl3(i)/t_pl3(i);
  q_pl3(i)=(d_pl3(i)+c_pl3(i)*q_pl3(i-1))/t_pl3(i);
end 
      
phi_pl3(N)=q_pl3(N);
for i=N-1:-1:1
  phi_pl3(i)=q_pl3(i)+p_pl3(i)*phi_pl3(i+1);           
end
%----------------------------- FIM TDMA PL3 -----------------------------%


%===================== POWER-LAW 4 (Peclet=500 e Sc=10)==================%

%-------------------------- Calculo dos Coeficientes---------------------%

% condicao de contorno i=1
  ap_pl4(1)=1;
  ae_pl4(1)=0;
  aw_pl4(1)=0;
  d_pl4(1)=phi0;
  
% pontos internos
for i=2:N-1
                   
  ae_pl4(i)=De(i)*Ae_pl2(i)+max(-Fe2,0);
  aw_pl4(i)=Dw(i)*Aw_pl2(i)+max(Fw2,0);
  ap_pl4(i)=ae_pl4(i)+aw_pl4(i)-Sp*Area*dist_face(i);
           
  d_pl4(i)=Sc2*Area*dist_face(i);
  
end
  
% condicao de contorno i=N
  ap_pl4(N)=1;
  ae_pl4(N)=0;
  aw_pl4(N)=0;
  d_pl4(N)=phi1;  

  for i=1:N
      a_pl4(i)=ap_pl4(i);
      b_pl4(i)=ae_pl4(i);
      c_pl4(i)=aw_pl4(i);
  end
 
%----------------------------- TDMA PL4 ---------------------------------%
% Calculo de phi_pl4 nos pontos

p_pl4(1)=b_pl4(1)/a_pl4(1); 
q_pl4(1)=d_pl4(1)/a_pl4(1);

for i=2:N 
  t_pl4(i)=a_pl4(i)-c_pl4(i)*p_pl4(i-1);
  p_pl4(i)=b_pl4(i)/t_pl4(i);
  q_pl4(i)=(d_pl4(i)+c_pl4(i)*q_pl4(i-1))/t_pl4(i);
end 
      
phi_pl4(N)=q_pl4(N);
for i=N-1:-1:1
  phi_pl4(i)=q_pl4(i)+p_pl4(i)*phi_pl4(i+1);           
end
%----------------------------- FIM TDMA PL4 -----------------------------%

%=============== DIFERENÇA CENTRAL 1 (Peclet=10 e Sc=0)==================%
%-------------------------- Calculo dos Coeficientes---------------------%

% condicao de contorno i=1
  ap_dc1(1)=1;
  ae_dc1(1)=0;
  aw_dc1(1)=0;
  d_dc1(1)=phi0;
  
% pontos internos
for i=2:N-1
           
  Ae_dc1(i)=1-(abs(Pe1(i)))/2;
  Aw_dc1(i)=1-(abs(Pw1(i)))/2;
           
  ae_dc1(i)=De(i)*Ae_dc1(i)+max(-Fe1,0);
  aw_dc1(i)=Dw(i)*Aw_dc1(i)+max(Fw1,0);
  ap_dc1(i)=ae_dc1(i)+aw_dc1(i)-Sp*Area*dist_face(i);
           
  d_dc1(i)=Sc1*Area*dist_face(i);
  
end
  
% condicao de contorno i=N
  ap_dc1(N)=1;
  ae_dc1(N)=0;
  aw_dc1(N)=0;
  d_dc1(N)=phi1;  

  for i=1:N
      a_dc1(i)=ap_dc1(i);
      b_dc1(i)=ae_dc1(i);
      c_dc1(i)=aw_dc1(i);
  end
 
%----------------------------- TDMA DC1 -----------------------------------%
% Calculo de phi_dc1 nos pontos

p_dc1(1)=b_dc1(1)/a_dc1(1); 
q_dc1(1)=d_dc1(1)/a_dc1(1);

for i=2:N 
  t_dc1(i)=a_dc1(i)-c_dc1(i)*p_dc1(i-1);
  p_dc1(i)=b_dc1(i)/t_dc1(i);
  q_dc1(i)=(d_dc1(i)+c_dc1(i)*q_dc1(i-1))/t_dc1(i);
end 
      
phi_dc1(N)=q_dc1(N);
for i=N-1:-1:1
  phi_dc1(i)=q_dc1(i)+p_dc1(i)*phi_dc1(i+1);           
end
%----------------------------- FIM TDMA DC1 -------------------------------%

%=============== DIFERENÇA CENTRAL 2 (Peclet=10 e Sc=10)==================%
%-------------------------- Calculo dos Coeficientes---------------------%

% condicao de contorno i=1
  ap_dc2(1)=1;
  ae_dc2(1)=0;
  aw_dc2(1)=0;
  d_dc2(1)=phi0;
  
% pontos internos
for i=2:N-1
           
  ae_dc2(i)=De(i)*Ae_dc1(i)+max(-Fe1,0);
  aw_dc2(i)=Dw(i)*Aw_dc1(i)+max(Fw1,0);
  ap_dc2(i)=ae_dc2(i)+aw_dc2(i)-Sp*Area*dist_face(i);
           
  d_dc2(i)=Sc2*Area*dist_face(i);
  
end
  
% condicao de contorno i=N
  ap_dc2(N)=1;
  ae_dc2(N)=0;
  aw_dc2(N)=0;
  d_dc2(N)=phi1;  

  for i=1:N
      a_dc2(i)=ap_dc2(i);
      b_dc2(i)=ae_dc2(i);
      c_dc2(i)=aw_dc2(i);
  end
 
%----------------------------- TDMA DC2 -----------------------------------%
% Calculo de phi_dc1 nos pontos

p_dc2(1)=b_dc2(1)/a_dc2(1); 
q_dc2(1)=d_dc2(1)/a_dc2(1);

for i=2:N 
  t_dc2(i)=a_dc2(i)-c_dc2(i)*p_dc2(i-1);
  p_dc2(i)=b_dc2(i)/t_dc2(i);
  q_dc2(i)=(d_dc2(i)+c_dc2(i)*q_dc2(i-1))/t_dc2(i);
end 
      
phi_dc2(N)=q_dc2(N);
for i=N-1:-1:1
  phi_dc2(i)=q_dc2(i)+p_dc2(i)*phi_dc2(i+1);           
end
%----------------------------- FIM TDMA DC2 -----------------------------%
%=============== DIFERENÇA CENTRAL 3 (Peclet=500 e Sc=0)================%
%-------------------------- Calculo dos Coeficientes---------------------%

% condicao de contorno i=1
  ap_dc3(1)=1; 
  ae_dc3(1)=0;
  aw_dc3(1)=0;
  d_dc3(1)=phi0;
  
% pontos internos
for i=2:N-1
             
  Ae_dc2(i)=1-(abs(Pe2(i)))/2;
  Aw_dc2(i)=1-(abs(Pw2(i)))/2;

  ae_dc3(i)=De(i)*Ae_dc2(i)+max(-Fe2,0);
  aw_dc3(i)=Dw(i)*Aw_dc2(i)+max(Fw2,0);
  ap_dc3(i)=ae_dc3(i)+aw_dc3(i)-Sp*Area*dist_face(i);
           
  d_dc3(i)=Sc1*Area*dist_face(i);
  
end
  
% condicao de contorno i=N
  ap_dc3(N)=1;
  ae_dc3(N)=0;
  aw_dc3(N)=0;
  d_dc3(N)=phi1;  

  for i=1:N
      a_dc3(i)=ap_dc3(i);
      b_dc3(i)=ae_dc3(i);
      c_dc3(i)=aw_dc3(i);
  end
 
%----------------------------- TDMA DC3 -----------------------------------%
% Calculo de phi_pl3 nos pontos

p_dc3(1)=b_dc3(1)/a_dc3(1); 
q_dc3(1)=d_dc3(1)/a_dc3(1);

for i=2:N 
  t_dc3(i)=a_dc3(i)-c_dc3(i)*p_dc3(i-1);
  p_dc3(i)=b_dc3(i)/t_dc3(i);
  q_dc3(i)=(d_dc3(i)+c_dc3(i)*q_dc3(i-1))/t_dc3(i);
end 
      
phi_dc3(N)=q_dc3(N);
for i=N-1:-1:1
  phi_dc3(i)=q_dc3(i)+p_dc3(i)*phi_dc3(i+1);           
end
%----------------------------- FIM TDMA DC3 -------------------------------%

%=============== DIFERENÇA CENTRAL 4 (Peclet=500 e Sc=10)================%
%-------------------------- Calculo dos Coeficientes---------------------%

% condicao de contorno i=1
  ap_dc4(1)=1;
  ae_dc4(1)=0;
  aw_dc4(1)=0;
  d_dc4(1)=phi0;
  
% pontos internos
for i=2:N-1
           
  ae_dc4(i)=De(i)*Ae_dc2(i)+max(-Fe2,0);
  aw_dc4(i)=Dw(i)*Aw_dc2(i)+max(Fw2,0);
  ap_dc4(i)=ae_dc4(i)+aw_dc4(i)-Sp*Area*dist_face(i);
           
  d_dc4(i)=Sc2*Area*dist_face(i);
  
end
  
% condicao de contorno i=N
  ap_dc4(N)=1;
  ae_dc4(N)=0;
  aw_dc4(N)=0;
  d_dc4(N)=phi1;  

  for i=1:N
      a_dc4(i)=ap_dc4(i);
      b_dc4(i)=ae_dc4(i);
      c_dc4(i)=aw_dc4(i);
  end
 
%----------------------------- TDMA DC4 -----------------------------------%
% Calculo de phi_pl4 nos pontos

p_dc4(1)=b_dc4(1)/a_dc4(1); 
q_dc4(1)=d_dc4(1)/a_dc4(1);

for i=2:N 
  t_dc4(i)=a_dc4(i)-c_dc4(i)*p_dc4(i-1);
  p_dc4(i)=b_dc4(i)/t_dc4(i);
  q_dc4(i)=(d_dc4(i)+c_dc4(i)*q_dc4(i-1))/t_dc4(i);
end 
      
phi_dc4(N)=q_dc4(N);
for i=N-1:-1:1
  phi_dc4(i)=q_dc4(i)+p_dc4(i)*phi_dc4(i+1);           
end
%----------------------------- FIM TDMA DC4 ----------------------------%

%===================== PÓS-PROCESSAMENTO ==============================%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SOLUÇÃO ANALÍTICA %%%%%%%%%%%%%%%%%%%%%%%%%%
% para P=10 e S*=0
for i=1:N
    
    P1S1_1(i)=exp(P1*x_ponto(i)/L);
    P1S1_2(i)=phi0-((1/(1-exp(P1)))*((phi0*((S1/P1)-exp(P1))+phi1*(1-(S1/P1)))));
    P1S1_3(i)=(1/(1-exp(P1)))*((phi0*((S1/P1)-exp(P1))+phi1*(1-(S1/P1))));
    P1S1_4(i)=-(S1*x_ponto(i)/(L*P1))*phi0+(S1*x_ponto(i)/(L*P1))*phi1;
    
    phiP1S1(i)=P1S1_1(i)*P1S1_2(i)+P1S1_3(i)+P1S1_4(i)
end

% para P=10 e S*=10
for i=1:N
    
    P1S2_1(i)=exp(P1*x_ponto(i)/L);
    P1S2_2(i)=phi0-((1/(1-exp(P1)))*((phi0*((S2/P1)-exp(P1))+phi1*(1-(S2/P1)))));
    P1S2_3(i)=(1/(1-exp(P1)))*((phi0*((S2/P1)-exp(P1))+phi1*(1-(S2/P1))));
    P1S2_4(i)=-(S2*x_ponto(i)/(L*P1))*phi0+(S2*x_ponto(i)/(L*P1))*phi1;
    
    phiP1S2(i)=P1S2_1(i)*P1S2_2(i)+P1S2_3(i)+P1S2_4(i);
end

% para P=500 e S*=0
for i=1:N
    
    P2S1_1(i)=exp(P2*x_ponto(i)/L);
    P2S1_2(i)=phi0-((1/(1-exp(P2)))*((phi0*((S1/P2)-exp(P2))+phi1*(1-(S1/P2)))));
    P2S1_3(i)=(1/(1-exp(P2)))*((phi0*((S1/P2)-exp(P2))+phi1*(1-(S1/P2))));
    P2S1_4(i)=-(S1*x_ponto(i)/(L*P2))*phi0+(S1*x_ponto(i)/(L*P2))*phi1;
    
    phiP2S1(i)=P2S1_1(i)*P2S1_2(i)+P2S1_3(i)+P2S1_4(i);
end

% para P=500 e S*=10
for i=1:N
    
    P2S2_1(i)=exp(P2*x_ponto(i)/L);
    P2S2_2(i)=phi0-((1/(1-exp(P2)))*((phi0*((S2/P2)-exp(P2))+phi1*(1-(S2/P2)))));
    P2S2_3(i)=(1/(1-exp(P2)))*((phi0*((S2/P2)-exp(P2))+phi1*(1-(S2/P2))));
    P2S2_4(i)=-(S2*x_ponto(i)/(L*P2))*phi0+(S2*x_ponto(i)/(L*P2))*phi1;
    
    phiP2S2(i)=P2S2_1(i)*P2S2_2(i)+P2S2_3(i)+P2S2_4(i);
end

%%%%%% UPWIND %%%%%%%
figure (1)
hold on
plot(x_ponto,phi_up1, 'r-.')
plot(x_ponto,phiP1S1, 'bs-')
title ('Solucao Upwind versus Exata para P=10 e S*=0')
legend('Solução numérica','Solução exata')
xlabel ('x')
ylabel ('phiUP e phiExata')
grid on

figure (2)
hold on
plot(x_ponto,phi_up2, 'r-.')
plot(x_ponto,phiP1S2, 'bs-')
title ('Solucao Upwind versus Exata para P=10 e S*=10')
legend('Solução numérica','Solução exata')
xlabel ('x')
ylabel ('phiUP e phiExata')
grid on

figure (3)
hold on
plot(x_ponto,phi_up3, 'r-.')
plot(x_ponto,phiP2S1, 'bs-')
title ('Solucao Upwind versus Exata para P=500 e S*=0')
legend('Solução numérica','Solução exata')
xlabel ('x')
ylabel ('phiUP e phiExata')
grid on

figure (4)
hold on
plot(x_ponto,phi_up4, 'r-.'), grid
plot(x_ponto,phiP2S2, 'bs-'), grid
title ('Solucao Upwind versus Exata para P=500 e S*=10')
legend('Solução numérica','Solução exata')
xlabel ('x')
ylabel ('phiUP e phiExata')
grid on

%%%%%% POWER-LAW %%%%%%%
figure (5)
hold on
plot(x_ponto,phi_pl1, 'r-.'), grid
plot(x_ponto,phiP1S1, 'bs-'), grid
title ('Solucao power-law versus exata para P=10 e S*=0')
legend('Solução numérica','Solução exata')
xlabel ('x')
ylabel ('phiPL e phiExata')
grid on

figure (6)
hold on
plot(x_ponto,phi_pl2, 'r-.'), grid
plot(x_ponto,phiP1S2, 'bs-'), grid
title ('Solucao power-law versus exata para P=10 e S*=10')
legend('Solução numérica','Solução exata')
xlabel ('x')
ylabel ('phiPL2 e phiExata')
grid on

figure (7)
hold on
plot(x_ponto,phi_pl3, 'r-.'), grid
plot(x_ponto,phiP2S1, 'bs-'), grid
title ('Solucao power-law versus exata para P=500 e S*=0')
legend('Solução numérica','Solução exata')
xlabel ('x')
ylabel ('phiPL e phiExata')
grid on

figure (8)
hold on
plot(x_ponto,phi_pl4, 'r-.'), grid
plot(x_ponto,phiP2S2, 'bs-'), grid
title ('Solucao power-law versus exata para P=500 e S*=10')
legend('Solução numérica','Solução exata')
xlabel ('x')
ylabel ('phiPL e phiExata')
grid on

%%%%%% DIFERENÇA CENTRAL %%%%%%%
figure (9)
hold on
plot(x_ponto,phi_dc1, 'r-.'), grid
plot(x_ponto,phiP1S1, 'bs-'), grid
title ('Solucao diferença central versus exata para P=10 e S*=0')
legend('Solução numérica','Solução exata')
xlabel ('x')
ylabel ('phiDC1 e phiExata')
grid on

figure (10)
hold on
plot(x_ponto,phi_dc2, 'r-.'), grid
plot(x_ponto,phiP1S2, 'bs-'), grid
title ('Solucao diferença central versus exata para P=10 e S*=10')
legend('Solução numérica','Solução exata')
xlabel ('x')
ylabel ('phiDC e phiExata')
grid on

figure (11)
hold on
plot(x_ponto,phi_dc3, 'r-.'), grid
plot(x_ponto,phiP2S1, 'bs-'), grid
title ('Solucao diferença central versus exata para P=500 e S*=0')
legend('Solução numérica','Solução exata')
xlabel ('x')
ylabel ('phiDC e phiExata')
grid on

figure (12)
hold on
plot(x_ponto,phi_dc4, 'r-.'), grid
plot(x_ponto,phiP2S2, 'bs-'), grid
title ('Solucao diferença central versus exata para P=500 e S*=10')
legend('Solução numérica','Solução exata')
xlabel ('x')
ylabel ('phiDC e phiExata')
grid on




