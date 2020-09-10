

%% %%%%%%%%% Donn?es %%%%%%%%%%%

% G?om?trie
L=1.25;
w=0.08;
h=0.02;

% Mat?riau
rho=2690;
E=6.45e10;
v=0.39;

% Support
k=10800;

% Importation des donn?es experimentales
[frequences]=load('Project_2019_freq.txt');
[modes]=load('Project_2019_modes.txt');

%% %%%%%%%%% Initialisation %%%%%%%%%%%

n=10;

M=zeros(2,2);
K=zeros(2,2);


x_abs = -0.625:0.001:0.625;
X = [-0.59 -0.5 -0.4 -0.3 -0.2 -0.1 0 0.1 0.2 0.3 0.4 0.5 0.59];


mode_Raleigh_Ritz=zeros(10,length(X));

% Initialisation de l'inconnue
syms x

M_RR=zeros(10,10);
K_RR=zeros(10,10);

p=zeros(length(X),5);

%% %%%%%%%%%%%%%%%%%%%%%% 2 DDL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Masses
m_poutre=L*w*h*rho;
m1=1;
m2=1;
jm1=2.112895833e-3; % JM1 + (61*10e-3/2 + 20*10e-3/2)^2
jm2=2.112895833e-3;
% Inertie flexionnelle
I=w*h^3/12;

% Moment d'inertie de la masse 
d=0.3; % distance masse-x7
 

% Moment d'inertie de la poutre


% Matrice de masse
%M=[m_poutre+m2 0 ; 0 0.77827];
M=[7.38 -0.2; -0.2 1.0447];

% Distance des ressorts ? x7
d1=-0.59;
d2=0.59;

% Matrice de raideur
%K=[2*k (d1+d2)*k ; (d1+d2)*k (d1^2+d2^2)*k];
K=[2*k 0; 0 2*k*(0.59)^2];

% D?termination des pulsation propres et modes propres
[vecteur_propre valeur_propre]=eigs(K,M);

% Normalisation et calculs des fr?quences
for i=1:2
   f_2DDL(i)=sqrt(valeur_propre(i,i))/(2*pi);
   mode_2DDL(2,i) = vecteur_propre(2,i)/vecteur_propre(2,i);
   mode_2DDL(1,i) = vecteur_propre(1,i)/vecteur_propre(2,i);
end

 f_2DDL
 mode_2DDL



%% %%%%%%%%%%%%%%%%%%%%%% RALEIGH-RITZ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calcul masse lineique
m_lin=m_poutre/L;


% D?termination de la matrice de masse
for i=0:9
    for j=0:9
        M_RR(i+1,j+1)=int(m_lin*(x/L)^i*(x/L)^j,x,-L/2,L/2)+(0.3/L)^i*(0.3/L)^j +(-0.5/L)^i*(-0.5/L)^j + jm1*( (-0.5)^(i-1)/ (L^i) )*( (-0.3)^(j-1)/ (L^j) ) + jm2*( (-0.5)^(i-1)/ (L^i) )*( (-0.3)^(j-1)/ (L^j) )  ; 
        %rotations et translation des masses
    end
end

% D?termination de la matrice de raideur
for i=0:9
    for j=0:9
        K_RR(i+1,j+1)=int(E*I*diff(diff((x/L)^i,x),x)*diff(diff((x/L)^j,x),x),x,-L/2,L/2)+k*(-0.59/L)^i*(-0.59/L)^j+k*(0.59/L)^i*(0.59/L)^j;
    end
end

% D?termination des pulsation propres et modes propres
[vecteur_propre valeur_propre]=eig(K_RR,M_RR);

% Calcul de la valeur des fonctions d'approximation

for i = 1:length(x_abs)
    for j=0:9
    w(j+1,i)= (x_abs(i)/L)^j;
    end
end

mode_Raleigh_Ritz = zeros(10,length(x_abs));
% Calcul des modes 
for i=1:10
     f_RR(i)=sqrt(valeur_propre(i,i))/(2*pi);
    for j=1:10
     mode_Raleigh_Ritz(j,:)=mode_Raleigh_Ritz(j,:)+vecteur_propre(i,j)*w(i,:);
    end   
end

f_RR =fliplr(f_RR) ;
f_RR= f_RR(1:5) ;
f_RR

    

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% ANALYSE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Erreur 2DDL
for i=1:2
    erreur_2DDL(i)=abs( (f_2DDL(i)-frequences(i))/frequences(i) );
end

erreur_2DDL

% Erreur Raleigh-Ritz
for i=1:5
    erreur_RR(i)=abs( (f_RR(i)-frequences(i))/frequences(i) );
end

erreur_RR

%% Poutre
poutre=zeros(1,length(X));

%% Graphiques

% D?termination et normalisation des modes de Raleigh Ritz

    %%mode 1
    mode_exp_graphe_1=modes(:,1)/modes(7,1);
    p_1=interp1(X,mode_exp_graphe_1,x_abs,'spline');
    mode_Raleigh_Ritz_graphe_1=mode_Raleigh_Ritz(11-2,:)/mode_Raleigh_Ritz(11-2,find(x_abs==0));
    
    %%mode 2
    mode_exp_graphe_2=modes(:,2)/modes(7,2);
    p_2=interp1(X,mode_exp_graphe_2,x_abs,'spline');
    mode_Raleigh_Ritz_graphe_2=mode_Raleigh_Ritz(11-1,:)/mode_Raleigh_Ritz(11-1,find(x_abs==0));
    
    
    % Affichage du graphique mode 1 et 2
    figure;
    
    subplot(2,1,1)
    hold on;
    %plot(X,poutre,'color',[0.8 0.8 0.8],'linewidth',6)
    plot(X,poutre,'k','linewidth',3)
    plot(x_abs,p_1,'b')
    plot(x_abs,mode_Raleigh_Ritz_graphe_1,'r')
    set(gca, 'fontsize', 10);
    title(['Mode n° ' num2str(1)],'fontsize',12);
    legend('Poutre non-déformée','Mode expérimental','Mode théorique');
    xlabel('x [m]','fontsize',12);
    ylabel('allure de déformée y [-]','fontsize',11);
    grid on;
    
    subplot(2,1,2)
    hold on;
    %plot(X,poutre,'color',[0.8 0.8 0.8],'linewidth',6)
    plot(X,poutre,'k','linewidth',3)
    plot(x_abs,p_2,'b')
    plot(x_abs,mode_Raleigh_Ritz_graphe_2,'r')
    set(gca, 'fontsize', 10);
    title(['Mode n° ' num2str(2)],'fontsize',12);
    legend('Poutre non-déformée','Mode expérimental','Mode théorique');
    xlabel('x [m]','fontsize',12);
    ylabel('allure de déformée y [-]','fontsize',11);
    grid on;

    
    
    %%mode 3
    mode_exp_graphe_3=modes(:,3)/modes(7,3);
    p_3=interp1(X,mode_exp_graphe_3,x_abs,'spline');
    mode_Raleigh_Ritz_graphe_3=mode_Raleigh_Ritz(11-3,:)/mode_Raleigh_Ritz(11-3,find(x_abs==0));
    
    
    % Affichage du graphique mode 3
    figure;
    subplot(2,1,1)
    hold on;
    %plot(X,poutre,'color',[0.8 0.8 0.8],'linewidth',6)
    plot(X,poutre,'k','linewidth',3)
    plot(x_abs,p_3,'b')
    plot(x_abs,mode_Raleigh_Ritz_graphe_3,'r')
    set(gca, 'fontsize', 10);
    title(['Mode n° ' num2str(3)],'fontsize',12);
    legend('Poutre non-déformée','Mode expérimental','Mode théorique');
    xlabel('x [m]','fontsize',12);
    ylabel('allure de déformée y [-]','fontsize',11);
    grid on;
    
    
    %%mode 4
    mode_exp_graphe_4=modes(:,4)/modes(7,4);
    p_4=interp1(X,mode_exp_graphe_4,x_abs,'spline');
    mode_Raleigh_Ritz_graphe_4=mode_Raleigh_Ritz(11-4,:)/mode_Raleigh_Ritz(11-4,find(x_abs==0));
    
    %%mode 5
    mode_exp_graphe_5=modes(:,5)/modes(7,5);
    p_5=interp1(X,mode_exp_graphe_5,x_abs,'spline');
    mode_Raleigh_Ritz_graphe_5=mode_Raleigh_Ritz(11-5,:)/mode_Raleigh_Ritz(11-5,find(x_abs==0));
    
    % Affichage du graphique mode 4 et 5
    %figure;
    
    %subplot(2,1,1)
    subplot(2,1,2)
    hold on;
    %plot(X,poutre,'color',[0.8 0.8 0.8],'linewidth',6)
    plot(X,poutre,'k','linewidth',3)
    plot(x_abs,p_4,'b')
    plot(x_abs,mode_Raleigh_Ritz_graphe_4,'r')
    set(gca, 'fontsize', 8);
    title(['Mode n° ' num2str(4)],'fontsize',12);
    legend('Poutre non-déformée','Mode expérimental','Mode théorique');
    xlabel('x [m]','fontsize',12);
    ylabel('allure de déformée y [-]','fontsize',11);
    grid on;
    
    figure;
    %subplot(2,1,2)
    hold on;
    %plot(X,poutre,'color',[0.8 0.8 0.8],'linewidth',6)
    plot(X,poutre,'k','linewidth',3)
    plot(x_abs,p_5,'b')
    plot(x_abs,mode_Raleigh_Ritz_graphe_5,'r')
    set(gca, 'fontsize', 10);
    title(['Mode n° ' num2str(5)],'fontsize',12);
    legend('Poutre non-déformée','Mode expérimental','Mode théorique');
    xlabel('x [m]','fontsize',12);
    ylabel('allure de déformée y [-]','fontsize',11);
    grid on;
