% 2º Trabalho de Dinamica dos Fluidos Computacional
% Algoritmo para calculo da distribuição de temperaturas em uma aleta
% finita de seção constante

% 1-PRE-PROCESSAMENTO------------------------------------------------------
% 1.1-Dados do problema
% 1.1.a-Parametros adimensionais (para considerar na discretização adimensionalizada):
A_st1=1;            % area da seção transversal
L1=1;               % comprimento da aleta
k1=1;               % coefeciente de transferencia de calor por condução
gamma1=1;           % coeficiente de difusão 
betaq1=3;           % coeficiente beta ao quadrado                             
biot1=20;           % número de Biot
tetab1=20;          % temperatura da base

% 1.1.b-Parametros variados (para considerar na discretização dimensionalizada):
A_st2=1;                
P2=1;                    
L2=3/20;                 
k2=1;                % coefeciente de transferencia de calor por condução
h2=400/3;            % coeficiente de transferencia de calor por convecção
betaq2=h2*P2*(L2^2)/(k2*A_st2);
beta2=sqrt(betaq2);        
biot2=h2*L2/k2; 

A_st3=25;                
P3=0.75;                    
L3=5;                 
k3=1;               
h3=4;            
betaq3=h3*P3*(L3^2)/(k3*A_st3);
beta3=sqrt(betaq3);        
biot3=h3*L3/k3; 

tetab2=20;
tetab3=20;
tetainf=5;          % temperatura no infinito
% 1.1.c-Volume de controle e número de pontos 
nvc=10;             % numero de volumes de controle (10, 20 ou 40) 
N=nvc+2;            % numero de pontos
%-------------------------------------------------------------------------%
% 1.2-Criação da malha
% 1.2.a-Vetor da coordenada das faces
alpha=1; % coeficiente constante que depende da uniformidade da malha
for i=1:N+1
    if (i==1) && (i==2)
        x_face1(i)=0;
    elseif (i>=3) && (i<N+1)
        x_face1(i)=L1*((i-2)/(N-2))^alpha; % coordenada das faces
    elseif i==N+1
        x_face1(i)=x_face1(N);
    end
end
for i=1:N+1
    if (i==1) && (i==2)
        x_face2(i)=0;
    elseif (i>=3) && (i<N+1)
        x_face2(i)=L2*((i-2)/(N-2))^alpha; 
    elseif i==N+1
        x_face2(i)=x_face2(N);
    end
end
for i=1:N+1
    if (i==1) && (i==2)
        x_face3(i)=0;
    elseif (i>=3) && (i<N+1)
        x_face3(i)=L3*((i-2)/(N-2))^alpha; 
    elseif i==N+1
        x_face3(i)=x_face3(N);
    end
end
% 1.2.b-vetor distancia entre as faces (tamanho volume)
for i=1:N %(N=nvc+1 e nvc=N-1, nro pontos=nro faces, dist_faces=nro volumes)
    dist_face1(i)=x_face1(i+1)-x_face1(i);
end
for i=1:N 
    dist_face2(i)=x_face2(i+1)-x_face2(i);
end
for i=1:N 
    dist_face3(i)=x_face3(i+1)-x_face3(i);
end
% 1.2.c-vetor coordenada dos pontos
for i=1:N
    x_ponto1(i)=((x_face1(i+1)+x_face1(i)))/2; % coordenada dos pontos
    x_ponto1(N)=x_face1(N);
end
for i=1:N
    x_ponto2(i)=((x_face2(i+1)+x_face2(i)))/2; % coordenada dos pontos
    x_ponto2(N)=x_face2(N);
end
for i=1:N
    x_ponto3(i)=((x_face3(i+1)+x_face3(i)))/2; % coordenada dos pontos
    x_ponto3(N)=x_face3(N);
end
% 1.2.d-vetor distancia entre os pontos
for i=1:N-1 %quantidade desse dist corresponde a quantidade de volumes
    dist_p1(i)=x_ponto1(i+1)-x_ponto1(i); %distancia entre os pontos internos
end
for i=1:N-1 %quantidade desse dist corresponde a quantidade de volumes
    dist_p2(i)=x_ponto2(i+1)-x_ponto2(i); %distancia entre os pontos internos
end
for i=1:N-1 %quantidade desse dist corresponde a quantidade de volumes
    dist_p3(i)=x_ponto3(i+1)-x_ponto3(i); %distancia entre os pontos internos
end
%-------------------------------------------------------------------------%
% PROCESSAMENTO
% 2-Calculo dos coeficientes 
% 2.1-discretização adimensional
for i=1:N
       if i ==1
           a1(i)=1;
           b1(i)=0;
           c1(i)=0;
           d1(i)=tetab1;
       elseif (i>=2) && (i<N)
           b1(i)=gamma1*A_st1/dist_p1(i);
           c1(i)=gamma1*A_st1/dist_p1(i-1);
           a1(i)=b1(i)+c1(i)+betaq1*dist_face1(i);
           d1(i)=0;
       else
           d1(i)=0;
           b1(i)=0;
           c1(i)=gamma1*A_st1/dist_p1(i-1); % condição no contorno
           a1(i)=b1(i)+c1(i)+biot1;         % condição no contorno
       end
   end
% 2.2-discretização dimensional   
   for i=1:N
       if i ==1
           a2(i)=1;
           b2(i)=0;
           c2(i)=0;
           d2(i)=tetab2;
       elseif (i>=2) && (i<N)
           b2(i)=k2*A_st2/dist_p2(i);
           c2(i)=k2*A_st2/dist_p2(i-1);
           a2(i)=b2(i)+c2(i)+h2*P2*dist_face2(i);
           d2(i)=h2*P2*tetainf*dist_face2(i);
       else
           b2(i)=0;
           d2(i)=h2*A_st2*tetainf;
           c2(i)=k2*A_st2/dist_p2(i-1); % condição no contorno
           a2(i)=b2(i)+c2(i)+h2*A_st2;  % condição no contorno
       end
   end
      for i=1:N
       if i ==1
           a3(i)=1;
           b3(i)=0;
           c3(i)=0;
           d3(i)=tetab3;
       elseif (i>=2) && (i<N)
           b3(i)=k3*A_st3/dist_p3(i);
           c3(i)=k3*A_st3/dist_p3(i-1);
           a3(i)=b3(i)+c3(i)+h3*P3*dist_face3(i);
           d3(i)=h3*P3*tetainf*dist_face3(i);
       else
           b3(i)=0;
           d3(i)=h3*A_st3*tetainf;
           c3(i)=k3*A_st3/dist_p3(i-1);
           a3(i)=b3(i)+c3(i)+h3*A_st3;
       end
   end
% 2.3-procedimento de Thomas
% Caso adimensional
for i=1:N
    if i==1
        p1(i)=b1(i)/a1(i); 
        q1(1)=d1(i)/a1(i);
    else
        t1(i)=a1(i)-c1(i)*p1(i-1);
        p1(i)=b1(i)/t1(i);
        q1(i)=(d1(i)+c1(i)*q1(i-1))/t1(i);
    end
end
for i=N:-1:1
    if i==N
        x1(i)=q1(i);
    elseif (i<=N-1) && (i>=2)
        x1(i)=q1(i)+p1(i)*x1(i+1);
    else
        x1(i)=tetab1;
    end
end
% Caso dimensional
for i=1:N
    if i==1
        p2(i)=b2(i)/a2(i); 
        q2(1)=d2(i)/a2(i);
    else
        t2(i)=a2(i)-c2(i)*p2(i-1);
        p2(i)=b2(i)/t2(i);
        q2(i)=(d2(i)+c2(i)*q2(i-1))/t2(i);
    end
end
for i=N:-1:1
    if i==N
        x2(i)=q2(i);
    elseif (i<=N-1) && (i>=2)
        x2(i)=q2(i)+p2(i)*x2(i+1);
    else
        x2(i)=tetab2;
    end
end
for i=1:N
    if i==1
        p3(i)=b3(i)/a3(i); 
        q3(1)=d3(i)/a3(i);
    else
        t3(i)=a3(i)-c3(i)*p3(i-1);
        p3(i)=b3(i)/t3(i);
        q3(i)=(d3(i)+c3(i)*q3(i-1))/t3(i);
    end
end
for i=N:-1:1
    if i==N
        x3(i)=q3(i);
    elseif (i<=N-1) && (i>=2)
        x3(i)=q3(i)+p3(i)*x3(i+1);
    else
        x3(i)=tetab3;
    end
end
for i=1:N
    Tadm1(i)=x1(i)/tetab1;
    Tadm2(i)=(x2(i)-tetainf)/(tetab2-tetainf);
    Tadm3(i)=(x3(i)-tetainf)/(tetab3-tetainf); 
end
% Temperaturas adimensionais
figure (1)
hold on
% perfil adimensional (L1) ao longo da aleta
plot(x_ponto1/L1,Tadm1, 'g'), grid
% distribuição de temperaturas adimensionais
plot(x_ponto1/L1,Tadm1, 'g'), grid
plot(x_ponto2/L2,Tadm2, 'b-o'), grid
plot(x_ponto3/L3,Tadm3, 'c*'), grid

legend('L1=1','L2=0.15', 'L3=5')
title ('Distribuição de temperatura adimensional')
xlabel ('Comprimento da aleta adimensional')
ylabel ('Temperatura adimensional')

%-------------------------------------------------------------------------%
% 3-POS-PROCESSAMENTO
% 3.1-Solução exata
m=sqrt(h2*P2/(k2*A_st2));
den2=cosh(m*L2)+(h2/(m*L2))*sinh(m*L2);
den3=cosh(m*L3)+(h3/(m*L3))*sinh(m*L3);
for i=1:N
    Ex2(i)=(cosh(m*L2*(1-x_ponto2(i)/L2)+(h2/(m*L2))*sinh(m*L2*(1-x_ponto2(i)/L2))))/den2;
    Ex3(i)=(cosh(m*L3*(1-x_ponto3(i)/L3)+(h3/(m*L3))*sinh(m*L3*(1-x_ponto3(i)/L3))))/den3;
end

% 3.2-Erros normalizados
for i=1:N
    erron2(i)=abs((Tadm2(i)-Ex2(i))/(Ex2(1)-Ex2(N)));
    erron3(i)=abs((Tadm3(i)-Ex3(i))/(Ex3(1)-Ex3(N)));
end
figure (2)
hold on

plot(x_ponto2/L2,erron2, 'b-o'), grid
plot(x_ponto3/L3,erron3, 'c*'), grid
title ('Distribuição dos erros normalizados')
xlabel ('Comprimento da aleta adimensional')
ylabel ('Erro normalizado')
csvwrite('Erro normalizado2_10vc.csv',[(x_ponto2/L2), erron2]);
csvwrite('Erro normalizado3_10vc.csv',[(x_ponto3/L3), erron3]);

% 3.3-Perda de calor adimensional na base da aleta
M2=(sqrt(h2*P2*k2*A_st2))*tetab2;
M3=(sqrt(h3*P3*k3*A_st3))*tetab3;
Mad2=M2*(L2/(k2*A_st2*tetab2)); %coeficiente M adimensionalizado. Madm=raiz(Betaq)
Mad3=M3*(L3/(k3*A_st3*tetab3));
q_a2=Mad2*((sinh(m*L2)+(h2/(m*L2))*cosh(m*L2)))/(cosh(m*L2)+(h2/(m*L2))*sinh(m*L2));
q_a3=Mad3*((sinh(m*L3)+(h3/(m*L3))*cosh(m*L3)))/(cosh(m*L3)+(h3/(m*L3))*sinh(m*L3));

% 3.4-Balanço global
betaq=3;
for i=1:N
    x_adm2(i)=(x_ponto2(i)/L2);
    x_adm3(i)=(x_ponto3(i)/L3);
end

Jin2=-(Tadm2(2)-Tadm2(1))/(x_adm2(2)-x_adm2(1));
Jout2=-(Tadm2(N)-Tadm2(N-1))/(x_adm2(N)-x_adm2(N-1));
soma2=0;
for i=1:N
    S2(i)=-betaq*Tadm2(i)*(dist_face2(i))/L2;
    soma2=soma2+S2(i);
end
Bal2=Jin2-Jout2+soma2;

Jin3=-(Tadm3(2)-Tadm3(1))/(x_adm3(2)-x_adm3(1));
Jout3=-(Tadm3(N)-Tadm3(N-1))/(x_adm3(N)-x_adm3(N-1));
soma3=0;
for i=1:N
    S3(i)=-betaq*Tadm3(i)*(dist_face3(i))/L3;
    soma3=soma3+S3(i);
end
Bal3=Jin3-Jout3+soma3;


% %-------------------------------------------------------------------------%
% 
% % Exportação dos dados para o excel
% 
% % 10 volumes de controle
% csvwrite('trabcf10volTadm1.csv',[(x_ponto1/L1)', Tadm1']); 
% csvwrite('trabcf10volTadm2.csv',[(x_ponto2/L2)', Tadm2']); 
% csvwrite('trabcf10volTadm3.csv',[(x_ponto3/L3)', Tadm3']); 
% 
% csvwrite('trabcf10volx_ponto1.csv',[x_ponto1]);
% csvwrite('trabcf10volx_ponto2.csv',[x_ponto2]); 
% csvwrite('trabcf10volx_ponto3.csv',[x_ponto3]);
% 
% csvwrite('trabcf10volx_ponto1adm.csv',[x_ponto1/L1]);
% csvwrite('trabcf10volx_ponto2adm.csv',[x_ponto2/L2]); 
% csvwrite('trabcf10volx_ponto3adm.csv',[x_ponto3/L3]); 
% 
% csvwrite('trabcf10Temperatura1.csv',[x1]);
% csvwrite('trabcf10Temperatura2.csv',[x2]);
% csvwrite('trabcf10Temperatura3.csv',[x3]);
% 
% csvwrite('trabcf10dist_face1.csv',[dist_face1]);
% csvwrite('trabcf10dist_face2.csv',[dist_face2]);
% csvwrite('trabcf10dist_face3.csv',[dist_face3]);
% 
% % % 20 volumes de controle
% csvwrite('trabcf20volTadm1.csv',[(x_ponto1/L1)', Tadm1']); 
% csvwrite('trabcf20volTadm2.csv',[(x_ponto2/L2)', Tadm2']); 
% csvwrite('trabcf20volTadm3.csv',[(x_ponto3/L3)', Tadm3']); 
% 
% csvwrite('trabcf20volx_ponto1.csv',[x_ponto1]);
% csvwrite('trabcf20volx_ponto2.csv',[x_ponto2]); 
% csvwrite('trabcf20volx_ponto2.csv',[x_ponto3]); 
% 
% csvwrite('trabcf20Temperatura1.csv',[x1]);
% csvwrite('trabcf20Temperatura2.csv',[x2]);
% csvwrite('trabcf20Temperatura3.csv',[x3]);
% 
% csvwrite('trabcf20dist_face1.csv',[dist_face1]);
% csvwrite('trabcf20dist_face2.csv',[dist_face2]);
% csvwrite('trabcf20dist_face3.csv',[dist_face3]);
% % 
% % % 40 volumes de controle
% csvwrite('trabcf40volTadm1.csv',[(x_ponto1/L1)', Tadm1']); 
% csvwrite('trabcf40volTadm2.csv',[(x_ponto2/L2)', Tadm2']); 
% csvwrite('trabcf40volTadm3.csv',[(x_ponto3/L3)', Tadm3']); 
% 
% csvwrite('trabcf40volx_ponto1.csv',[x_ponto1]);
% csvwrite('trabcf40volx_ponto2.csv',[x_ponto2]); 
% csvwrite('trabcf40volx_ponto2.csv',[x_ponto3]); 
% 
% csvwrite('trabcf40Temperatura1.csv',[x1]);
% csvwrite('trabcf40Temperatura2.csv',[x2]);
% csvwrite('trabcf40Temperatura3.csv',[x3]);
% 
% csvwrite('trabcf40dist_face1.csv',[dist_face1]);
% csvwrite('trabcf40dist_face2.csv',[dist_face2]);
% csvwrite('trabcf40dist_face3.csv',[dist_face3]);

