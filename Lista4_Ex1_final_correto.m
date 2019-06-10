clear all;
clc;

%=============4� Lista de Dinamica dos Fluidos Computacional=============%

% Aluna: Carolina Maria Nunes Bezerra

% Algoritmo para resolver uma equa��o de difus�o bidimensional em
% coordenadas cartesianas, regime permamente, utilizando o m�todo A. Testar
% o algoritmo utilizando a solu��o manufaturada.

%============================1 PRE-PROCESSAMENTO=========================%

% 1.1-Dados do problema

Lx=1;   % comprimento unit�rio em x
Ly=1;   % comprimento unit�rio em y
Nx=41;  % numero de pontos em x
Ny=41;  % numero de pontos em y

% Cria��o da malha (M�TODO A)--------------------------------------------%

% Vetor coordenada dos pontos

alpha=1;
% Na dire��o x
for i=1:Nx
    x_ponto(i)=Lx*((i-1)/(Nx-1))^alpha;
end
% Na dire��o y
for j=1:Ny
    y_ponto(j)=Ly*((j-1)/(Ny-1))^alpha;
end

% vetor distancia entre os pontos

% Na dire��o x
for i=1:Nx-1 
    dist_p_x(i)=x_ponto(i+1)-x_ponto(i); 
end
% Na dire��o y
for j=1:Ny-1 
    dist_p_y(j)=y_ponto(j+1)-y_ponto(j); 
end

% Vetor da coordenada das faces

% Na dire��o x
for i=1:Nx-1
    x_face(i)=(x_ponto(i)+x_ponto(i+1))/2;
end
% Na dire��o y
for j=1:Ny-1
    y_face(j)=(y_ponto(j)+y_ponto(j+1))/2;
end

% Vetor distancia entre as faces (tamanho volume)

% Na dire��o x
dist_face_x(1)=x_face(1);
for i=2:Nx-1
    dist_face_x(i)=x_face(i)-x_face(i-1);
end
dist_face_x(Nx)=x_ponto(Nx)-x_face(Nx-1);
% Na dire��o y
dist_face_y(1)=y_face(1);
for j=2:Ny-1 
    dist_face_y(j)=y_face(j)-y_face(j-1);
end
dist_face_y(Ny)=y_ponto(Ny)-y_face(Ny-1);

%==========================2 PROCESSAMENTO===============================%

% Calculo da fun��o manufaturada |Phi| para o vetor [d] nas fronteiras 
% (condi��o dada)

for i=1:Nx
    for j=1:Ny
        Phi(i,j)=x_ponto(i)^4*y_ponto(j)^5;
    end
end

% Calculo dos coeficientes ----------------------------------------------%

for i=1:Nx  
    for j=1:Ny 
        
    % Escolha para Sp ser nulo
        Sp(i,j)=0;
    % Coeficientes de difus�o
        gamma_e(i)=1;
        gamma_w(i)=1;
        gamma_n(j)=1;
        gamma_s(j)=1;
        
    % Condi��o para os pontos do contorno (1,j), (i,1), (Nx,j), (i,Ny) 
    if i==1||i==Nx||j==1||j==Ny
        ae(i,j)=0;    % leste
        aw(i,j)=0;    % oeste
        an(i,j)=0;    % norte
        as(i,j)=0;    % sul
        ap(i,j)=1;    % ponto i (principal)
        Sc(i,j)=0;
        d(i,j)=Phi(i,j);
    
        a(i,j)=ap(i,j);
        b(i,j)=ae(i,j);
        c(i,j)=aw(i,j);
        e(i,j)=an(i,j);
        f(i,j)=as(i,j);
        
    % Pontos internos 
    % Fonte Sc a partir da solu��o fun��o manufaturada 
    else

        Sc(i,j)=-(12*x_ponto(i)^2*y_ponto(j)^5+20*x_ponto(i)^4*y_ponto(j)^3);          
        
        ae(i,j)=gamma_e(i)*dist_face_y(j)/dist_p_x(i);
        aw(i,j)=gamma_w(i)*dist_face_y(j)/dist_p_x(i-1);
        an(i,j)=gamma_n(j)*dist_face_x(i)/dist_p_y(j);
        as(i,j)=gamma_s(j)*dist_face_x(i)/dist_p_y(j-1);
        ap(i,j)=ae(i,j)+aw(i,j)+an(i,j)+as(i,j)-Sp(i,j)*dist_face_x(i)*dist_face_y(i);
        d(i,j)=Sc(i,j)*dist_face_x(i)*dist_face_y(i);         
  
        a(i,j)=ap(i,j);
        b(i,j)=ae(i,j);
        c(i,j)=aw(i,j);
        e(i,j)=an(i,j);
        f(i,j)=as(i,j);
        
    end
    end
end

%========================================================================%
%===========================TESTE DOS M�TODOS============================%
%========================================================================%

tolerancia=10^-8;   % toler�ncia como crit�rio para convergencia

%================PROCEDIMENTO ITERATIVO DE GAUSS SEIDEL==================%

ResGS=1;            % Res�duo de Guass Seidel(inicializa��o)
ItGS=0;             % itera��es de Guass Seidel(inicializa��o)

% Calculo para os valores iniciais da fun��o nos pontos
for i=1:Nx
    for j=1:Ny 
    % Para pontos da fronteira (1,j), (i,1), (Nx,j), (i,Ny) o valor de 
    % PhiGS � a fun��o manufaturada Phi
        if i==1||i==Nx||j==1||j==Ny
           PhiGS(i,j)=Phi(i,j);
        else
           PhiGS(i,j)=0; % chute inicial 
        end   
    end
    
end

% Itera��o usando o m�todo de Gauss Seidel
while ResGS>tolerancia
    
    for i=2:Nx-1
        for j=2:Ny-1
            PhiGS(i,j)=(c(i,j)*PhiGS(i-1,j)+b(i,j)*PhiGS(i+1,j)+f(i,j)*PhiGS(i,j-1)+e(i,j)*PhiGS(i,j+1)+d(i,j))/a(i,j);
        end
    end
    
% C�lculo do res�duo
    for i=2:Nx-1
        for j=2:Ny-1
            RGS(i,j)=abs(a(i,j)*PhiGS(i,j)-(b(i,j)*PhiGS(i+1,j)+c(i,j)*PhiGS(i-1,j)+e(i,j)*PhiGS(i,j+1)+f(i,j)*PhiGS(i,j-1)+d(i,j)));
        end
    end
    
 ItGS=ItGS+1
 % Maior valor da matriz dos res�duos
 ResGS=max(max(RGS))  
end

%========================TDMA LINHA POR LINHA===========================%
% d(i,j) ser� o termo somando os coeficientes norte e sul mais o termo
% fonte Ss manufaturado.
for i=1:Nx
    for j=1:Ny 
    % Para pontos da fronteira (1,j), (i,1), (Nx,j), (i,Ny) o valor de 
    % PhiT � a fun��o manufaturada Phi
        if i==1||i==Nx||j==1||j==Ny
           PhiT(i,j)=Phi(i,j);
        else
           PhiT(i,j)=0; % chute inicial 
        end   
    end
    
end

% Itera��o usando o m�todo de TDMA Linha por Linha
ResT=1;     % Res�duo do TDMA (inicializa��o)
ItT=0;      % Itera��es do TDMA(inicializa��o)

% Como os coeficientes escolhidos para o b* foram o norte e sul, o TDMA ir�
% caminhar coluna por coluna a partir dos [an] e [as] calculados anteriormente
% e para cada coluna os valores na dire��o i ser�o calculados de uma vez
% por itera��o.
% Na primeira itera��o os valores de PhiT recebem os valores do loop fora 
% do while, onde calculou-se a solu��o PhiT a partir da fun��o manufaturada
% no contorno e chutado o valor de zero nos pontos internos. A seguir,
% esses valores ser�o atualizados com o TDMA at� atingir a convergencia
% desejada para o res�duo.

while ResT>tolerancia 
    for j=2:Ny-1
        
        bf(1,j)=e(1,j)*PhiT(1,j+1)+f(1,j)*PhiT(1,j-1)+d(1,j);
        P(1)=b(1,j)/a(1,j); 
        Q(1)=bf(1,j)/a(1,j);
                
        for i=2:Nx
            bf(i,j)=e(i,j)*PhiT(i,j+1)+f(i,j)*PhiT(i,j-1)+d(i,j);
            T(i)=a(i,j)-c(i,j)*P(i-1);
            P(i)=b(i,j)/T(i);
            Q(i)=(bf(i,j)+c(i,j)*Q(i-1))/T(i);
        end
        PhiT(Nx,j)=Q(Nx);

        for i=Nx-1:-1:1
            PhiT(i,j)=P(i)*PhiT(i+1,j)+Q(i);
        end   
    end
   % C�lculo do res�duo 
    for i=2:Nx-1
        for j=2:Ny-1
            RT(i,j)=abs(a(i,j)*PhiT(i,j)-(b(i,j)*PhiT(i+1,j)+c(i,j)*PhiT(i-1,j)+e(i,j)*PhiT(i,j+1)+f(i,j)*PhiT(i,j-1)+d(i,j)));
        end
    end
    
ItT=ItT+1
% Maior valor da matriz dos res�duos
ResT=max(max(RT))
end

%=========================================================================%
%=========COMPARA��O DA SOLU��O DOS M�TODOS COM A SOLU��O EXATA===========%
%=========================================================================%

% Solu��o exata
for i=1:Nx
    for j=1:Ny
        Phi_exata(i,j)=x_ponto(i)^4*y_ponto(j)^5;    
    end
end

% Erro m�ximo
for i=1:Nx
    for j=1:Ny
        EGS(i,j)=abs(PhiGS(i,j)-Phi_exata(i,j));
        ET(i,j)=abs(PhiT(i,j)-Phi_exata(i,j));
    end
end
EGS_max=max(max(EGS));
ET_max=max(max(ET));

% �ndice da ordem para o erro de discretiza��o
i=3;
IGS=EGS_max/(dist_face_x(i))^2; 
IT=ET_max/(dist_face_x(i))^2;

% Superf�cie da solu��o exata
figure (1)
[ x , y ] = meshgrid( x_ponto , y_ponto) ;
surf(x,y,Phi_exata) ;
axis([0 1 0 1 0 max(max(Phi_exata))]);
title('Superf�cie da solu��o exata');
xlabel('y');
ylabel('x');
zlabel('Phi_exata');

% Superf�cie da solu��o num�rica de Guass Seidel
figure (2)
[ x , y ] = meshgrid( x_ponto , y_ponto) ;
surf(x,y,PhiGS) ;
axis([0 1 0 1 0 max(max(PhiGS))]);
title('Superf�cie da solu��o num�rica de Gauss Seidel');
xlabel('y');
ylabel('x');
zlabel('PhiGS');

% Superf�cie da solu��o num�rica do TDMA Linha por Linha
figure (3)
[ x , y ] = meshgrid( x_ponto , y_ponto) ;
surf(x,y,PhiT) ;
axis([0 1 0 1 0 max(max(PhiT))]);
title('Superf�cie da solu��o num�rica do TDMA Linha por Linha');
xlabel('y');
ylabel('x');
zlabel('PhiT');









