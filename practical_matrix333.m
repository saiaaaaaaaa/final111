clc,clear,close all ;
step=0.001;  
s_time=50;

time=zeros(1, s_time/step);

 alpha11=zeros(1, s_time/step);
 alpha21=zeros(1, s_time/step);
 alpha_1_1_2=zeros(1, s_time/step);
 alpha_1_2_2=zeros(1, s_time/step);

 d_omega_1_2=zeros(1, s_time/step);
 d_omega_2_2=zeros(1, s_time/step);
 dy_1_d=zeros(1, s_time/step);
 dy_2_d=zeros(1, s_time/step);
 fu_k_a_1=zeros(1, s_time/step);
 fu_k_a_2=zeros(1, s_time/step);
 fu_k_c_1=zeros(1, s_time/step);
 fu_k_c_2=zeros(1, s_time/step);

 g11=zeros(1, s_time/step);
 g12=zeros(1, s_time/step);
 g21=zeros(1, s_time/step);
 g22=zeros(1, s_time/step);

 k_a_1=zeros(1, s_time/step);
 k_a_2=zeros(1, s_time/step);
 k_b11=zeros(1, s_time/step);
 k_b11diff=zeros(1, s_time/step);
 k_b12=zeros(1, s_time/step);
 k_b12diff=zeros(1, s_time/step);
 k_b21=zeros(1, s_time/step);
 k_b21diff=zeros(1, s_time/step);
 k_b22=zeros(1, s_time/step);
 k_b22diff=zeros(1, s_time/step);
 k_c_1=zeros(1, s_time/step);
 k_c_2=zeros(1, s_time/step);

 kensai_1_1=zeros(1, s_time/step);
 kensai_1_2=zeros(1, s_time/step);
 kensai_2_1=zeros(1, s_time/step);
 kensai_2_2=zeros(1, s_time/step);
 M_1=zeros(1, s_time/step);
 M_2=zeros(1, s_time/step);

 u_1_1_2=zeros(1, s_time/step);
 u_1_2_2=zeros(1, s_time/step);
 v_1_1_2=zeros(1, s_time/step);
 v_1_2_2=zeros(1, s_time/step);
 var1=zeros(1, s_time/step);
 var2=zeros(1, s_time/step);

 y_1_d=zeros(1, s_time/step);
 y_2_d=zeros(1, s_time/step);
 z_1_1=zeros(1, s_time/step);
 z_1_2=zeros(1, s_time/step);
 z_2_1=zeros(1, s_time/step);
 z_2_2=zeros(1, s_time/step);
 
 

 omega_1_2_2=zeros(1, s_time/step+1);
 omega_2_1=zeros(1, s_time/step+1);
 omega_2_2=zeros(1, s_time/step+1);
 omega_2_2_2=zeros(1, s_time/step+1);

 Theta_1_1_1=zeros(1, s_time/step+1);
 Theta_1_1_2=zeros(1, s_time/step+1);
 Theta_1_2_1=zeros(1, s_time/step+1);
 Theta_1_2_2=zeros(1, s_time/step+1);

 varrho_1_1_1=zeros(1, s_time/step+1);
 varrho_1_1_2=zeros(1, s_time/step+1);
 varrho_1_2_1=zeros(1, s_time/step+1);
 varrho_1_2_2=zeros(1, s_time/step+1);

 x_j1_i1=zeros(1, s_time/step+1);
 x_j1_i2=zeros(1, s_time/step+1);
 x_j2_i1=zeros(1, s_time/step+1);
 x_j2_i2=zeros(1, s_time/step+1);



 func_Theta_1_1_1=zeros(1, s_time/step);
 func_Theta_1_1_2=zeros(1, s_time/step);
 func_Theta_1_2_1=zeros(1, s_time/step);
 func_Theta_1_2_2=zeros(1, s_time/step);

 func_varrho_1_1_1=zeros(1, s_time/step);
 func_varrho_1_1_2=zeros(1, s_time/step);
 func_varrho_1_2_1=zeros(1, s_time/step);
 func_varrho_1_2_2=zeros(1, s_time/step);

 func_x_i1_j1=zeros(1, s_time/step);
 func_x_i2_j1=zeros(1, s_time/step);
 func_x_i1_j2=zeros(1, s_time/step);
 func_x_i2_j2=zeros(1, s_time/step);
 

 

r_1_1=zeros(1, s_time/step+1);
r_2_1=zeros(1, s_time/step+1);
r_1_2=zeros(1, s_time/step+1);
r_2_2=zeros(1, s_time/step+1);









x_j1_i1(1)=0.2; 
x_j1_i2(1)=1.5;
x_j2_i1(1)=0.2;  
x_j2_i2(1)=1.6;

k_b11(1)=3;
k_b12(1)=3;
k_b21(1)=3;
k_b22(1)=3;

Theta_1_1_1(1)=0.5;
Theta_1_1_2(1)=0.6;
Theta_1_2_1(1)=0.1;
Theta_1_2_2(1)=0.5;

varrho_1_1_1(1)=0.2;
varrho_1_1_2(1)=0.1;
varrho_1_2_1(1)=0.5;
varrho_1_2_2(1)=0.2;



iota_1=162;
iota_2=175;
gamma_1=0.95;

r_1_1(1)=0;
r_2_1(1)=0;
r_1_2(1)=0;
r_2_2(1)=0;
omega_1_2_2(1)=0;
omega_2_2_2(1)=0;


k11=2.5;
k12=2.3;


k21=2.2;
k22=2.6;

beta_1_1_1=0.8;
beta_1_1_2=0.2;
beta_1_2_1=0.8;
beta_1_2_2=0.2;

gamma_1_1_1=0.002;
gamma_1_1_2=0.005;
gamma_1_2_1=0.002;
gamma_1_2_2=0.005;

zeta_1_1_1=0.5;
zeta_1_1_2=0.6;
zeta_1_2_1=0.5;
zeta_1_2_2=0.3;

varsigma_1_1_1=18;
varsigma_1_1_2=15;
varsigma_1_2_1=16;  
varsigma_1_2_2=12;

Y_1_1=0.1;
Y_2_1=0.1;

Y0_1_2=1.1;
Y0_2_2=1.1;

% Self-triggered parameter
s11=0.15;
s12=0.3;
s13=12;
s14=0.15;

s21=0.15;
s22=0.3;
s23=9;
s24=0.1;

varpi_1=0.12;
varpi_2=0.1;

rou=3;
h=1;

numa_1=0;
numa_2=0;

count1_1=0;
count2_2=0;

D_1=0.1;
D_2=0.1;

t_next_1=1;
t_next_2=1;




rt1=[];
rr1=[];
rt11=[];
rr11=[];

rt2=[];
rr2=[];
rt22=[];
rr22=[];
a=0.25;

m_1=66.88;
m_2=56.81;
k1=0.08;
k2=0.08;

l_1=3.83;
l_2=2.27;
g=9.8;


for t=1:500000
    time(t)=(t-1)*step;
    
    y_1_d(t)=1*sin(1.5*t*step); 
    y_2_d(t)=0.8*sin(2*t*step);  
    dy_1_d(t)=1*1.5*cos(1.5*t*step)*exp(sin(1.5*t*step)); 
    dy_2_d(t)=-2*0.8*cos(2*t*step)*exp(sin(2*t*step)); 

    z_1_1(t)=x_j1_i1(t)-y_1_d(t);
    kensai_1_1(t)=z_1_1(t)-r_1_1(t);
    z_2_1(t)=x_j2_i1(t)-y_2_d(t);
    kensai_2_1(t)=z_2_1(t)-r_2_1(t);
    
    k_b11(t)=sin(1.1*t*step)+3.5;
    k_b12(t)=0.5*sin(0.9*t*step)+5;
    k_b21(t)=0.5*sin(1.1*t*step)+4.5;
    k_b22(t)=0.6*sin(0.5*t*step)+6;
    
    k_b11diff(t)=1.1*cos(1.1*t*step);
    k_b12diff(t)=0.45*cos(0.9*t*step);
    k_b21diff(t)=cos(t*step);
    k_b22diff(t)=0.3*cos(0.5*t*step); 
     
    k_a_1(t)=k_b11(t)+0.5*Y_1_1;
    k_a_2(t)=k_b21(t)+0.2*Y_2_1;
   
    k_c_1(t)=k_b12(t)+Y0_1_2;
    k_c_2(t)=k_b22(t)+Y0_2_2;
    
    fu_k_a_1(t)= (-(k_b11(t)+0.5*Y_1_1));
    fu_k_a_2(t)= (-(k_b21(t)+0.2*Y_2_1));
    
    fu_k_c_1(t)= (-(k_b12(t)+Y0_1_2));
    fu_k_c_2(t)= (-(k_b22(t)+Y0_2_2));
     
    Gaussian_a1_1=exp(-(x_j1_i1(t)-2)^2/2^2)*exp(-(y_1_d(t)-2)^2/2^2)*exp(-(dy_1_d(t)-2)^2/2^2);
    Gaussian_a1_2=exp(-(x_j1_i1(t)-1)^2/2^2)*exp(-(y_1_d(t)-1)^2/2^2)*exp(-(dy_1_d(t)-1)^2/2^2); 
    Gaussian_a1_3=exp(-(x_j1_i1(t)-0)^2/2^2)*exp(-(y_1_d(t)-0)^2/2^2)*exp(-(dy_1_d(t)-0)^2/2^2); 
    Gaussian_a1_4=exp(-(x_j1_i1(t)+1)^2/2^2)*exp(-(y_1_d(t)+1)^2/2^2)*exp(-(dy_1_d(t)+1)^2/2^2);
    Gaussian_a1_5=exp(-(x_j1_i1(t)+2)^2/2^2)*exp(-(y_1_d(t)+2)^2/2^2)*exp(-(dy_1_d(t)+2)^2/2^2);
    Gaussian_a1_sum=Gaussian_a1_1+Gaussian_a1_2+Gaussian_a1_3+Gaussian_a1_4+Gaussian_a1_5;   
    
    kappa_1_1=[Gaussian_a1_1/Gaussian_a1_sum,Gaussian_a1_2/Gaussian_a1_sum,Gaussian_a1_3/Gaussian_a1_sum,Gaussian_a1_4/Gaussian_a1_sum,Gaussian_a1_5/Gaussian_a1_sum]';
    
    Gaussian_c1_1=exp(-(x_j2_i1(t)-2)^2/2^2)*exp(-(y_2_d(t)-2)^2/2^2)*exp(-(dy_2_d(t)-2)^2/2^2); 
    Gaussian_c1_2=exp(-(x_j2_i1(t)-1)^2/2^2)*exp(-(y_2_d(t)-1)^2/2^2)*exp(-(dy_2_d(t)-1)^2/2^2); 
    Gaussian_c1_3=exp(-(x_j2_i1(t)-0)^2/2^2)*exp(-(y_2_d(t)-0)^2/2^2)*exp(-(dy_2_d(t)-0)^2/2^2); 
    Gaussian_c1_4=exp(-(x_j2_i1(t)+1)^2/2^2)*exp(-(y_2_d(t)+1)^2/2^2)*exp(-(dy_2_d(t)+1)^2/2^2);
    Gaussian_c1_5=exp(-(x_j2_i1(t)+2)^2/2^2)*exp(-(y_2_d(t)+2)^2/2^2)*exp(-(dy_2_d(t)+2)^2/2^2); 
    Gaussian_c1_sum=Gaussian_c1_1+Gaussian_c1_2+Gaussian_c1_3+Gaussian_c1_4+Gaussian_c1_5;   
    
    kappa_2_1=[Gaussian_c1_1/Gaussian_c1_sum,Gaussian_c1_2/Gaussian_c1_sum,Gaussian_c1_3/Gaussian_c1_sum,Gaussian_c1_4/Gaussian_c1_sum,Gaussian_c1_5/Gaussian_c1_sum]';
 
    g11(t)=1+0.01*sin(x_j1_i1(t));
    g12(t)=1+0.01*sin(x_j1_i2(t))*x_j1_i1(t); 
    g21(t)=1+0.01*sin(x_j2_i1(t));
    g22(t)=1+0.01*sin(x_j2_i2(t))*x_j2_i1(t);    
    alpha11(t)=(1/(-g11(t))*(((k11*(k_b11(t)^2)*sin((pi*(kensai_1_1(t))^2)/(2*(k_b11(t))^2))*cos((pi*(kensai_1_1(t))^2)/(2*(k_b11(t))^2)))/(pi*kensai_1_1(t)))-abs(k_b11diff(t)/k_b11(t))*(kensai_1_1(t))+(-1/g11(t))*(Theta_1_1_1(t)/(2*(varsigma_1_1_1)^2))*((norm(kappa_1_1))^2)*((sec((pi*(kensai_1_1(t))^2)/(2*(k_b11(t))^2)))^2)*(kensai_1_1(t))+0.5*(((sec((pi*(kensai_1_1(t))^2)/(2*(k_b11(t))^2)))^2)*(kensai_1_1(t)))+(1/2)*(((sec((pi*(kensai_1_1(t))^2)/(2*(k_b11(t))^2)))^2)*(kensai_1_1(t)))+varrho_1_1_1(t)));
    alpha21(t)=(1/(-g21(t))*(((k21*(k_b21(t)^2)*sin((pi*(kensai_2_1(t))^2)/(2*(k_b21(t))^2))*cos((pi*(kensai_2_1(t))^2)/(2*(k_b21(t))^2)))/(pi*kensai_2_1(t)))-abs(k_b21diff(t)/k_b21(t))*(kensai_2_1(t))+(-1/g21(t))*(Theta_1_2_1(t)/(2*(varsigma_1_2_1)^2))*((norm(kappa_2_1))^2)*((sec((pi*(kensai_2_1(t))^2)/(2*(k_b21(t))^2)))^2)*(kensai_2_1(t))+0.5*(((sec((pi*(kensai_2_1(t))^2)/(2*(k_b21(t))^2)))^2)*(kensai_2_1(t)))+(1/2)*(((sec((pi*(kensai_2_1(t))^2)/(2*(k_b21(t))^2)))^2)*(kensai_2_1(t)))+varrho_1_2_1(t)));
 
    omega_2_1(1)=alpha11(1);
    omega_2_1(t+1)=omega_2_1(t)+(iota_1 *omega_1_2_2(t) )*step;
    omega_1_2_2(t+1)=omega_1_2_2(t)+(-2*gamma_1 *iota_1*omega_1_2_2(t)-iota_1 *(omega_2_1(t)-alpha11(t)))*step;
    d_omega_1_2(t)=iota_1*omega_1_2_2(t);
    z_1_2(t)=x_j1_i2(t)-omega_2_1(t);
    kensai_1_2(t)=z_1_2(t)-r_1_2(t);
    
    omega_2_2(1)=alpha21(1);
    omega_2_2(t+1)=omega_2_2(t)+(iota_2*omega_2_2_2(t) )*step;
    omega_2_2_2(t+1)=omega_2_2_2(t)+(-2*gamma_1 *iota_2*omega_2_2_2(t)-iota_2 *(omega_2_2(t)-alpha21(t)))*step;
    d_omega_2_2(t)=iota_2*omega_2_2_2(t);
    z_2_2(t)=x_j2_i2(t)-omega_2_2(t);
    kensai_2_2(t)=z_2_2(t)-r_2_2(t);
    
    Gaussian_b1_1=exp(-(x_j1_i1(t)-2)^2/2^2)*exp(-(x_j1_i2(t)-2)^2/2^2)*exp(-(d_omega_1_2(t)-2)^2/2^2); 
    Gaussian_b1_2=exp(-(x_j1_i1(t)-1)^2/2^2)*exp(-(x_j1_i2(t)-1)^2/2^2)*exp(-(d_omega_1_2(t)-1)^2/2^2); 
    Gaussian_b1_3=exp(-(x_j1_i1(t)-0)^2/2^2)*exp(-(x_j1_i2(t)-0)^2/2^2)*exp(-(d_omega_1_2(t)-0)^2/2^2); 
    Gaussian_b1_4=exp(-(x_j1_i1(t)+1)^2/2^2)*exp(-(x_j1_i2(t)+1)^2/2^2)*exp(-(d_omega_1_2(t)+1)^2/2^2);
    Gaussian_b1_5=exp(-(x_j1_i1(t)+2)^2/2^2)*exp(-(x_j1_i2(t)+2)^2/2^2)*exp(-(d_omega_1_2(t)+2)^2/2^2); 
    Gaussian_b1_sum=Gaussian_b1_1+Gaussian_b1_2+Gaussian_b1_3+Gaussian_b1_4+Gaussian_b1_5;   
    
    kappa_1_2=[Gaussian_b1_1/Gaussian_b1_sum,Gaussian_b1_2/Gaussian_b1_sum,Gaussian_b1_3/Gaussian_b1_sum,Gaussian_b1_4/Gaussian_b1_sum,Gaussian_b1_5/Gaussian_b1_sum]';
    
    Gaussian_d1_1=exp(-(x_j2_i1(t)-2)^2/2^2)*exp(-(x_j2_i2(t)-2)^2/2^2)*exp(-(d_omega_2_2(t)-2)^2/2^2); 
    Gaussian_d1_2=exp(-(x_j2_i1(t)-1)^2/2^2)*exp(-(x_j2_i2(t)-1)^2/2^2)*exp(-(d_omega_2_2(t)-1)^2/2^2); 
    Gaussian_d1_3=exp(-(x_j2_i1(t)-0)^2/2^2)*exp(-(x_j2_i2(t)-0)^2/2^2)*exp(-(d_omega_2_2(t)-0)^2/2^2); 
    Gaussian_d1_4=exp(-(x_j2_i1(t)+1)^2/2^2)*exp(-(x_j2_i2(t)+1)^2/2^2)*exp(-(d_omega_2_2(t)+1)^2/2^2);
    Gaussian_d1_5=exp(-(x_j2_i1(t)+2)^2/2^2)*exp(-(x_j2_i2(t)+2)^2/2^2)*exp(-(d_omega_2_2(t)+2)^2/2^2); 
    Gaussian_d1_sum=Gaussian_d1_1+Gaussian_d1_2+Gaussian_d1_3+Gaussian_d1_4+Gaussian_d1_5;      
    
    kappa_2_2=[Gaussian_d1_1/Gaussian_d1_sum,Gaussian_d1_2/Gaussian_d1_sum,Gaussian_d1_3/Gaussian_d1_sum,Gaussian_d1_4/Gaussian_d1_sum,Gaussian_d1_5/Gaussian_d1_sum]';
        
    Theta_1_1_1(t+1)=Theta_1_1_1(t)+(((((norm(kappa_1_1))^2))/(2*(varsigma_1_1_1)^2))*((((sec((pi*(kensai_1_1(t))^2)/(2*(k_b11(t))^2)))^2)^2*(kensai_1_1(t)))^2)-beta_1_1_1*Theta_1_1_1(t))*step;
    Theta_1_1_2(t+1)=Theta_1_1_2(t)+(((((norm(kappa_1_2))^2))/(2*(varsigma_1_1_2)^2))*((((sec((pi*(kensai_1_2(t))^2)/(2*(k_b12(t))^2)))^2)^2*(kensai_1_2(t)))^2)-beta_1_1_2*Theta_1_1_2(t))*step;
    Theta_1_2_1(t+1)=Theta_1_2_1(t)+(((((norm(kappa_2_1))^2))/(2*(varsigma_1_2_1)^2))*((((sec((pi*(kensai_2_1(t))^2)/(2*(k_b21(t))^2)))^2)^2*(kensai_2_1(t)))^2)-beta_1_2_1*Theta_1_2_1(t))*step;
    Theta_1_2_2(t+1)=Theta_1_2_2(t)+(((((norm(kappa_2_2))^2))/(2*(varsigma_1_2_2)^2))*((((sec((pi*(kensai_2_2(t))^2)/(2*(k_b22(t))^2)))^2)^2*(kensai_2_2(t)))^2)-beta_1_2_2*Theta_1_2_2(t))*step;
    
    varrho_1_1_1(t+1)=varrho_1_1_1(t)+((gamma_1_1_1*((sec((pi*(kensai_1_1(t))^2)/(2*(k_b11(t))^2)))^2)*(kensai_1_1(t)))-zeta_1_1_1*varrho_1_1_1(t))*step;
    varrho_1_1_2(t+1)=varrho_1_1_2(t)+((gamma_1_1_2*((sec((pi*(kensai_1_2(t))^2)/(2*(k_b12(t))^2)))^2)*(kensai_1_2(t)))-zeta_1_1_2*varrho_1_1_2(t))*step;
    varrho_1_2_1(t+1)=varrho_1_2_1(t)+((gamma_1_2_1*((sec((pi*(kensai_2_1(t))^2)/(2*(k_b21(t))^2)))^2)*(kensai_2_1(t)))-zeta_1_2_1*varrho_1_2_1(t))*step;
    varrho_1_2_2(t+1)=varrho_1_2_2(t)+((gamma_1_2_2*((sec((pi*(kensai_2_2(t))^2)/(2*(k_b22(t))^2)))^2)*(kensai_2_2(t)))-zeta_1_2_2*varrho_1_2_2(t))*step;

    alpha_1_1_2(t)=-(1/(g12(t)))*(((k12*(k_b12(t)^2)*sin((pi*(kensai_1_2(t))^2)/(2*(k_b12(t))^2))*cos((pi*(kensai_1_2(t))^2)/(2*(k_b12(t))^2)))/(pi*kensai_1_2(t)))+(-abs(k_b12diff(t)/k_b12(t))*(kensai_1_2(t))+(Theta_1_1_2(t)/(2*(varsigma_1_1_2)^2))*((norm(kappa_1_2))^2)*((sec((pi*(kensai_1_2(t))^2)/(2*(k_b12(t))^2)))^2)*(kensai_1_2(t))+0.5*(((sec((pi*(kensai_1_2(t))^2)/(2*(k_b12(t))^2)))^2)*(kensai_1_2(t)))+varrho_1_1_2(t)));
    alpha_1_2_2(t)=-(1/(g22(t)))*(((k22*(k_b22(t)^2)*sin((pi*(kensai_2_2(t))^2)/(2*(k_b22(t))^2))*cos((pi*(kensai_2_2(t))^2)/(2*(k_b22(t))^2)))/(pi*kensai_2_2(t)))+(-abs(k_b22diff(t)/k_b22(t))*(kensai_2_2(t))+(Theta_1_2_2(t)/(2*(varsigma_1_2_2)^2))*((norm(kappa_2_2))^2)*((sec((pi*(kensai_2_2(t))^2)/(2*(k_b22(t))^2)))^2)*(kensai_2_2(t))+0.5*(((sec((pi*(kensai_2_2(t))^2)/(2*(k_b22(t))^2)))^2)*(kensai_2_2(t)))+varrho_1_2_2(t)));

    var1(t)=((sec((pi*(kensai_1_2(t))^2)/(2*(k_b12(t))^2)))^2)*(kensai_1_2(t));
    var2(t)=((sec((pi*(kensai_2_2(t))^2)/(2*(k_b22(t))^2)))^2)*(kensai_2_2(t));
    v_1_1_2(t)=-(1+varpi_1)*(alpha_1_1_2(t)*tanh((g12(t)*var1(t)*alpha_1_1_2(t))/rou)+h*tanh((g12(t)*var1(t)*h)/rou));
    v_1_2_2(t)=-(1+varpi_2)*(alpha_1_2_2(t)*tanh((g22(t)*var2(t)*alpha_1_2_2(t))/rou)+h*tanh((g22(t)*var2(t)*h)/rou));

M_1(t)=s11+(s12-s11)*tanh(1/(s13-s14*abs(z_1_1(t))));
if(t>=t_next_1)
    numa_1=numa_1+1;
    u_1_1_2(t)=v_1_1_2(t);
     rt1=[rt1;t*step];
     rr1=[rr1;t*step-count1_1];
     count1_1=t*step;
     rt11=[rt11;t*step];
     rr11=[rr11;1];
     count1_1=t*step;
    if(t_next_1==1)
        bb_1=D_1;
    else
        tt_1=t;
        hh_1=(v_1_1_2(tt_1)-v_1_1_2(tt_1-1))/step;
        if(hh_1>D_1)
            bb_1=hh_1;
        else
            bb_1=D_1;
        end
    end
    t_next_1=t;
    t_next_1=t_next_1+(varpi_1*abs(u_1_1_2(t))+M_1(t))/(bb_1);
else
    u_1_1_2(t)=u_1_1_2(t-1);
end

M_2(t)=s21+(s22-s21)*tanh(1/(s23-s24*abs(z_2_1(t))));
if(t>=t_next_2)
     numa_2=numa_2+1;
     u_1_2_2(t)=v_1_2_2(t);
     rt2=[rt2;t*step];
     rr2=[rr2;t*step-count2_2];
     count2_2=t*step;
     rt22=[rt22;t*step];
     rr22=[rr22;2];
     count2_2=t*step;
    if(t_next_2==1)
        bb_2=D_2;
    else
        tt_2=t;
        hh_2=(v_1_2_2(tt_2)-v_1_2_2(tt_2-1))/step;
        if(hh_2>D_2)
            bb_2=hh_2;
        else
            bb_2=D_2;
        end
    end
    t_next_2=t;
    t_next_2=t_next_2+(varpi_2*abs(u_1_2_2(t))+M_2(t))/(bb_2);
else
    u_1_2_2(t)=u_1_2_2(t-1);
end

    r_1_1(t+1)=r_1_1(t)+g11(t)*r_1_2(t)+g11(t)*(omega_2_1(t)-alpha11(t))*step;
    r_1_2(t+1)=r_1_2(t)+0*step;
    r_2_1(t+1)=r_2_1(t)+g21(t)*r_2_2(t)+g21(t)*(omega_2_2(t)-alpha21(t))*step;
    r_2_2(t+1)=r_2_2(t)+0*step;

    x_j1_i1(t+1)=x_j1_i1(t)+(x_j1_i2(t)+0.1*cos(2*x_j1_i1(t)))*step; 
    x_j1_i2(t+1)=x_j1_i2(t)+(1/(m_1*((l_1)^2))*u_1_1_2(t)-g/l_1*sin(x_j1_i1(t))-k1/m_1*x_j1_i2(t)+0.1*sin((2*x_j1_i2(t))*x_j1_i1(t)))*step; 
    x_j2_i1(t+1)=x_j2_i1(t)+(x_j2_i2(t)+0.1*cos(x_j2_i1(t)))*step; 
    x_j2_i2(t+1)=x_j2_i2(t)+(1/(m_2*((l_2)^2))*u_1_2_2(t)-g/l_2*sin(x_j2_i1(t))-k2/m_2*x_j2_i2(t)+0.1*sin((2*x_j2_i2(t))*x_j2_i1(t)))*step;
   
    func_x_i1_j1(t)=x_j1_i1(t);
    func_x_i2_j1(t)=x_j2_i1(t);
    func_x_i1_j2(t)=x_j1_i2(t);
    func_x_i2_j2(t)=x_j2_i2(t);
    func_Theta_1_1_1(t)=Theta_1_1_1(t);
    func_Theta_1_1_2(t)=Theta_1_1_2(t);
    func_Theta_1_2_1(t)=Theta_1_2_1(t);
    func_Theta_1_2_2(t)=Theta_1_2_2(t);
    
    func_varrho_1_1_1(t)=varrho_1_1_1(t);
    func_varrho_1_1_2(t)=varrho_1_1_2(t);
    func_varrho_1_2_1(t)=varrho_1_2_1(t);
    func_varrho_1_2_2(t)=varrho_1_2_2(t);
    
 
    

    if      t*step>=s_time
        break;
    end
end

sum1=0;
for t=1:50000
    sum1(t+1)=sum1(t)+abs(z_1_1(t))*step;
end
sum2=0;
for t=1:50000
    sum2(t+1)=sum2(t)+abs(z_2_1(t))*step;
end  
sum3=0;
for t=1:50000
    sum3(t+1)=sum3(t)+t*step*abs(z_1_1(t))*step;
end
sum4=0;
for t=1:50000
    sum4(t+1)=sum4(t)+t*step*abs(z_2_1(t))*step; 
end


figure()
 subplot(2,1,1)
 plot(time,func_x_i1_j1,'b-',time,y_1_d,'r--',time,k_a_1,'g-',time,fu_k_a_1,'c-','LineWidth',1.5)
 legend('$y_{1}$','$y_{1,d}$','${k_{{c_{1,1}}}}$','$-{k_{{c_{1,1}}}}$');
 set(legend,'Interpreter','latex')
%  ylabel('$y_{1}$,$y_{1,d}$ and $ \pm {k_{{c_{1,1}}}}$','Interpreter','latex');
 set(gca,'FontSize',12)
 subplot(2,1,2)
 plot(time,func_x_i2_j1,'b-',time,y_2_d,'r--',time,k_a_2,'g-',time,fu_k_a_2,'c-','LineWidth',1.5)
 legend('$y_{2}$','$y_{2,d}$','${k_{{c_{2,1}}}}$','$-{k_{{c_{2,1}}}}$');
 set(legend,'Interpreter','latex')
 xlabel('Time(Sec)','Interpreter','latex')
 set(gca,'FontSize',12)
ylabel('Trajectory tracking ($y_{1}$ and $y_{2}$)','Interpreter','latex')

 figure()
 subplot(2,1,1)
 plot(time,func_x_i1_j2,'b-',time,k_c_1,'r--',time,fu_k_c_1,'k--','LineWidth',1.5)
 legend('$x_{1,2}$','${k_{{c_{1,2}}}}$','$-{k_{{c_{1,2}}}}$');
 set(legend,'Interpreter','latex')
%  ylabel('$x_{1,2}$ and $\pm{k_{{c_{1,2}}}}$','Interpreter','latex')
 set(gca,'FontSize',12)
 subplot(2,1,2)
 plot(time,func_x_i2_j2,'b-',time,k_c_2,'r--',time,fu_k_c_2,'k--','LineWidth',1.5)
 legend('$x_{2,2}$','${k_{{c_{2,2}}}}$','$-{k_{{c_{2,2}}}}$');
 set(legend,'Interpreter','latex')
 ylabel('System states ($x_{1,2}$ and $x_{2,2}$)','Interpreter','latex')
 xlabel('Time(Sec)','Interpreter','latex')
 set(gca,'FontSize',12) 


figure()
subplot(2,1,1)
plot(time,func_Theta_1_1_1,'c-',time,func_Theta_1_1_2,'b-',time,func_varrho_1_1_1,'g-',time,func_varrho_1_1_2,'r-','linewidth',1.5)
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1,'LineWidth',1);
legend_1=legend('$\hat{\Theta}_{1,1}$','$\hat{\Theta}_{1,2}$','$\hat{D}_{1,1}$','$\hat{D}_{1,2}$');
set(legend_1,'Interpreter','latex')
ylabel('$\hat{\Theta}_{1,1}$,$\hat{\Theta}_{1,2}$,$\hat{D}_{1,1}$,$\hat{D}_{1,2}$','Interpreter','latex')
set(gca,'XTick',0:60/12:60,'FontSize',12)
subplot(2,1,2)
plot(time,func_Theta_1_2_1,'c-',time,func_Theta_1_2_2,'b-',time,func_varrho_1_2_1,'g-',time,func_varrho_1_2_2,'r-','linewidth',1.5)
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1,'LineWidth',1);
xlabel('Time(Sec)','Interpreter','latex')
legend_1=legend('$\hat{\Theta}_{2,1}$','$\hat{\Theta}_{2,2}$','$\hat{D}_{2,1}$','$\hat{D}_{2,2}$');
set(legend_1,'Interpreter','latex')
ylabel('$\hat{\Theta}_{2,1}$,$\hat{\Theta}_{2,2}$,$\hat{D}_{2,1}$,$\hat{D}_{2,2}$','Interpreter','latex')
set(gca,'XTick',0:60/12:60,'FontSize',12)

figure()
plot(time,u_1_1_2,'b-',time,v_1_1_2,'r--','linewidth',1.5)
xlabel('Time(Sec)','Interpreter','latex')
legend_1=legend('$u_1$','$v_1$');
set(legend_1,'Interpreter','latex')
ylabel('Self-triggered input $u_1$ and control input $v_1$','Interpreter','latex')
set(gca,'FontSize',12)
figure()
plot(time,u_1_2_2,'b-',time,v_1_2_2,'r--','LineWidth',1.5)
xlabel('Time(Sec)','Interpreter','latex') 
legend_6=legend('$u_2$','$v_2$'); 
set(legend_6,'Interpreter','latex') 
ylabel('Self-triggered input $u_2$ and control input $v_2$','Interpreter','latex')
set(gca,'FontSize',12)


