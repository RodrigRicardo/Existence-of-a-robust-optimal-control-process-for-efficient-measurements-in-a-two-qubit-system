clear all
clc

% We are going to produce a sample of pure and mixed two-qubit states and the associated concurrence
% using the parametrization of the following paper:
% Mixed state parametrization and two-qubit entaglement
% Otto C. W. Kong and Hock King Ting
%  arXiv:2112.10011v2
% For checking look at file parametrization density matriz pure states.wxmx

% We are going to mix pure and mixed states in the sample file. This is the difference with
% Qubit_sample_two_qubit_pure_and_mixed_states_new_parametrization_rand.m

function c=c(theta,psi)
  c=cos(theta/2)*exp(-1i*psi/2);

endfunction

function s=s(theta,psi)
  s=sin(theta/2)*exp(1i*psi/2);

endfunction

function e_0=e_0(qp,qm,zeta)

  e_0=[qp;
        0;
        0;
       qm*exp(1i*zeta)];

endfunction

function e_1=e_1(qp,qm,zeta,c21,s21,c32,s32)

  e_1=[-c21*qm*exp(-1i*zeta);
         s21*c32;
         s21*s32;
         c21*qp
        ];
endfunction


function e_2=e_2(qp,qm,zeta,c21,s21,c32,s32,c0,s0)

    e_2=[c0*conj(s21)*exp(-1i*zeta)*qm;
         c0*conj(c21)*c32-s0*conj(s32);
         c0*conj(c21)*s32+s0*conj(c32);
         -c0*conj(s21)*qp];

endfunction


function e_3=e_3(qp,qm,zeta,c21,s21,c32,s32,c0,s0)

     e_3=[-s0*s21*exp(1i*zeta)*qm;
         -s0*c21*conj(c32)+c0*s32;
         -s0*c21*conj(s32)+conj(c0)*c32;
          s0*s21*qp];

endfunction


function [rho,concurrence]=rhoState()

  sigma_y=[0,-i;i,0];
  sigmayy=kron(sigma_y,sigma_y);

% Defining the parameters of the density matrix

  r=(pi/2)*rand(1);
  qp=cos(r);
  qm=sin(r);
  zeta=2*pi*rand(1);
  num_par=randi([1,4]);
  theta= pi*rand(1,3);
  psi=2*pi*rand(1,3);
  c21=c(theta(1),psi(1));
  s21=s(theta(1),psi(1));
  c32=c(theta(2),psi(2));
  s32=s(theta(2),psi(2));
  c0=c(theta(3),psi(3));
  s0=s(theta(3),psi(3));
  e0=kron(e_0(qp,qm,zeta),ctranspose(e_0(qp,qm,zeta)));
  e1=kron(e_1(qp,qm,zeta,c21,s21,c32,s32),ctranspose(e_1(qp,qm,zeta,c21,s21,c32,s32)));
  e2=kron(e_2(qp,qm,zeta,c21,s21,c32,s32,c0,s0),ctranspose(e_2(qp,qm,zeta,c21,s21,c32,s32,c0,s0)));
  e3=kron(e_3(qp,qm,zeta,c21,s21,c32,s32,c0,s0),ctranspose(e_3(qp,qm,zeta,c21,s21,c32,s32,c0,s0)));
  % Let us work with the mixing parameters
  % They have to organize in descending order. We have 1 (pure states), and mixed states
  mu=rand(1,4);
  mu=sort(mu,'descend');
  sum_s=0;
  for i=1:num_par
     sum_s=sum_s+mu(i);
  endfor
  mu=mu/sum_s;
  rho1=zeros(4,4);
  A={mu(1)*e0,mu(2)*e1,mu(3)*e2,mu(4)*e3};
  for i = 1:num_par
     rho1=rho1+A{i};
  endfor
 % Now the concurrence part

  lambda=rho1*sigmayy*conj(rho1)*sigmayy;
  candi=sort(real(eig(lambda)),'descend');
  % Find the eigenvalues , take only the real part and sort them
  % in descending order.

  for i = 1:num_par
       sol(i)=candi(i);
  endfor
  if min(sol)>=0 % checking for non-negative matrix
        sol=sqrt(sol);
        sol=sort(sol); # Put them in ascending order
        val=sol(numel(sol)); # Take the greatest
        for i = 2:num_par
          val=val-candi(i);
        endfor
        concurrence=max(0,val);
        rho=rho1;
  endif
  #sol=real(eig(rho1))

endfunction


sample=10;

for j=1:sample
  j
  [rho,concurrence]=rhoState();
  states{j,1}=rho;
  states{j,2}=concurrence;
 endfor


str=mat2str(sample);
namefile=strcat('states_mixed_pure_',str,'_');
time=datestr(now,'yyyy_mm_dd_HH_MM_SS');
filename=strcat(namefile,time,'.mat');
save(filename,'states');
