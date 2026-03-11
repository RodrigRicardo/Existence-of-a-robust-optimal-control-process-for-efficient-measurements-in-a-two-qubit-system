
clear all
clc


%%%%%%%%%%%%%%%%%% THE CODE BEGINS HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PARAMETERS OF THE QUANTUM NEURAL NETWORK


hbar = 1.0/(2*pi);   % h = 1.0
ihbarinv = 1.0/(1i*hbar);


%%%%%%%%%%%%%%%%%%%%%%%
% DEFINITION OF THE ELEMENTS OF THE HAMILTONIAN

% set up the hamiltonian
sig_xa = [ 0    0   1   0
           0    0   0   1
           1    0   0   0
           0    1   0   0 ];
%
sig_xb = [ 0    1   0   0
           1    0   0   0
           0    0   0   1
           0    0   1   0 ];
%
sig_za = [ 1    0    0    0
           0    1    0    0
           0    0   -1    0
           0    0    0   -1 ];
%
sig_zb = [ 1     0   0    0
           0    -1   0    0
           0     0   1    0
           0     0   0   -1 ];
%
sig_zazb = sig_za*sig_zb;

sig_xaxb=[ 0     0   0    1
           0     0   1    0
           0     1   0    0
           1     0   0    0 ];


% Definition and initialization of the values for the weights
%Ka =Kainit;
%Kb =Kbinit;
%epsa=epsainit;
%epsb=epsbinit;
%zeta=zetainit;

%%%%%%%

% Calculate here the derivatives of H with respect to the weights
dHdKa = sig_xa;
dHdKb = sig_xb;
dHdepsa = sig_za;
dHdepsb = sig_zb;
dHdzeta = sig_zazb;
dHdnu = sig_xaxb;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% READING THE TRAINING SET

%state_file='states_mixed_pure_1000_2024_11_05_09_33_54.mat'
state_file='states_mixed_pure_100_2024_11_05_09_17_30.mat'

load(state_file,'states');

[sample,lines]=size(states);



  %%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization of learning rate for quantum neural network (this defines the step size)


  input_file='parameters_qq_';
  % This is the name of the file to save the results

  %%%%%%%%
 % Iniatilization of the parameters of the Hamiltonian
  deltat = 1.6*hbar;  %3*hbar;  % nanosec   % 10
  %% Learning constant for the quantum network

  iterations=1500;
  learning_const_qq=0.01;
  nlayers_min=4;
  % Minimum number of controls

%% USING THE GRADIENT METHOD
%% GUESSING ALL THE CONTROL VARIABLES &&
%% INITIALIZATION WITH RANDOM VALUES FOR ALL LAYERS IN THE QUANTUM NETWORK

  %%  Max and min values for the weights of the neural network, these need to be positive.



error_sample=zeros(1,sample);
concurrences_cal=zeros(1,sample);
concurrences=zeros(1,sample);

[rownum,colnum]=size(states{1,1});


  % Loop over number of training pairs
for jcount=1:sample

jcount
concurrences(jcount)=states{jcount,2};

##   rhoEvol(:,:,1)=states{jcount,1};
##   minval=-2;
##   maxval=2;
##   Ka = minval+(maxval-minval)*rand(nlayers-1,sample);
##   Kb = minval+(maxval-minval)*rand(nlayers-1,sample);
##   epsa = minval+(maxval-minval)*rand(nlayers-1,sample);
##   epsb= minval+(maxval-minval)*rand(nlayers-1,sample);
##   zeta= minval+(maxval-minval)*rand(nlayers-1,sample);
##   nu=minval+(maxval-minval)*rand(nlayers-1,sample);


iter=1;
nlayers=nlayers_min;

while iter<= iterations

    if iter==1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%

          % To keep the values over the sample loop we have the following vectors:

        rhoEvol=zeros(rownum,colnum,nlayers);
        H=zeros(rownum,colnum,nlayers-1);
        gamlam=zeros(rownum,colnum,nlayers);
        propagator=zeros(rownum,colnum,nlayers-1);
        invpropagator=zeros(rownum,colnum,nlayers-1);

        rhoEvol(:,:,1)=states{jcount,1};
        minval=-2;
        maxval=2;
        for j=1:nlayers-1
          Ka{j,jcount} = minval+(maxval-minval)*rand(1);
          Kb{j,jcount} = minval+(maxval-minval)*rand(1);
          epsa{j,jcount} =minval+(maxval-minval)*rand(1);
          epsb{j,jcount}= minval+(maxval-minval)*rand(1);
          zeta{j,jcount}= minval+(maxval-minval)*rand(1);
          nu{j,jcount}=minval+(maxval-minval)*rand(1);
        endfor
    endif

% FORWARD PROPAGATION OF DYNAMICS nlayers

   for ipath=1:nlayers-1

       % Defintion of the Hamiltonian, initial time ipath-1 and final time ipath
      H(:,:,ipath)=Ka{ipath,jcount}*sig_xa+Kb{ipath,jcount}*sig_xb+epsa{ipath,jcount}*sig_za+epsb{ipath,jcount}*sig_zb+zeta{ipath,jcount}*sig_zazb+nu{ipath,jcount}*sig_xaxb;

       % Definition of propagators

      propagator(:,:,ipath)=expm((ihbarinv*deltat)*H(:,:,ipath));
      invprop(:,:,ipath)=ctranspose(propagator(:,:,ipath));

        % Generating rho evolution

      rhoEvol(:,:,ipath+1)=propagator(:,:,ipath)*rhoEvol(:,:,ipath)*invprop(:,:,ipath);

     endfor % end for loop over the number of layers for forward propagation

    % Loop over number of iterations
    %%%%  NOW WE NEED TO CALCULATE THE DUAL VARIABLES LAMBDA

     %%  FIRST THE FINAL CONDITION:
     % states{jcount,19}=real(trace(sig_zazb*rhoEvol(:,:,nlayers)));
     error_qnn=states{jcount,2}-real(trace(sig_zazb*rhoEvol(:,:,nlayers)));
     gamlam(:,:,nlayers) =-error_qnn*sig_zazb;

%  Add this conditon to check if the algorithm has gone crazy and restart as needed

     if (abs(error_qnn)>0.5)&&(iter>(iterations*0.5))  % (abs(error_qnn)>0.5)&&(iter>(iterations*0.5))
       iter=1

     else

      %%  NOW THE OTHER VARIABLES:
          for ipath=nlayers-1:-1:1
                gamlam(:,:,ipath)=invprop(:,:,ipath)*gamlam(:,:,ipath+1)*propagator(:,:,ipath);

          endfor  % end for loop over the number of layers for backward propagation for dual variables



       %%  CORRECTIONS TO THE CONTROL PARAMETERS




        %
%        %%% CALCULATING THE DERIVATIVES AT EACH TIME
%        %%%% Just for safety reasons we take the real part of the derivatives
%
           Kadev(1)=0;
           Kbdev(1)=0;
           epsadev(1)=0;
           epsbdev(1)=0;
           zetadev(1)=0;
           nudev(1)=0;


          for ipath=1:nlayers-1

                % Ka derivative
                % Calculate here the commutator

                Kadot = ihbarinv*(dHdKa*rhoEvol(:,:,ipath)-rhoEvol(:,:,ipath)*dHdKa);
                Kadev(ipath)=Ka{ipath,jcount}+real(trace(gamlam(:,:,ipath+1)*Kadot));
                Ka{ipath,jcount} =Ka{ipath,jcount}-learning_const_qq*Kadev(ipath);


                % Kb corrections

                Kbdot = ihbarinv*(dHdKb*rhoEvol(:,:,ipath)-rhoEvol(:,:,ipath)*dHdKb);
                Kbdev(ipath)=Kb{ipath,jcount}+real(trace(gamlam(:,:,ipath+1)*Kbdot));
                Kb{ipath,jcount}=Kb{ipath,jcount}-learning_const_qq*Kbdev(ipath);

                % epsa corrections

                epsadot = ihbarinv*(dHdepsa*rhoEvol(:,:,ipath)-rhoEvol(:,:,ipath)*dHdepsa);
                epsadev(ipath)=epsa{ipath,jcount}+real(trace(gamlam(:,:,ipath+1)*epsadot));
                epsa{ipath,jcount} =epsa{ipath,jcount}-learning_const_qq*epsadev(ipath);

                % epsb  corrections

                epsbdot = ihbarinv*(dHdepsb*rhoEvol(:,:,ipath)-rhoEvol(:,:,ipath)*dHdepsb);
                epsbdev(ipath)=epsb{ipath,jcount}+real(trace(gamlam(:,:,ipath+1)*epsbdot));
                epsb{ipath,jcount} =epsb{ipath,jcount}-learning_const_qq*epsbdev(ipath);;

                % zeta corrections

                zetadot = ihbarinv*(dHdzeta*rhoEvol(:,:,ipath)-rhoEvol(:,:,ipath)*dHdzeta);
                zetadev(ipath)=zeta{ipath,jcount}+real(trace(gamlam(:,:,ipath+1)*zetadot));
                zeta{ipath,jcount} =zeta{ipath,jcount}-learning_const_qq*zetadev(ipath);

                % nu corrections

                nudot = ihbarinv*(dHdnu*rhoEvol(:,:,ipath)-rhoEvol(:,:,ipath)*dHdnu);
                nudev(ipath)=nu{ipath,jcount}+real(trace(gamlam(:,:,ipath+1)*nudot));
                nu{ipath,jcount}=nu{ipath,jcount}-learning_const_qq*nudev(ipath);

          endfor %    end for loop over the number of layers

     iter=iter+1;

   endif %  end of condition to continue.

   error_rel= abs(states{jcount,2}-real(trace(sig_zazb*rhoEvol(:,:,nlayers))))/states{jcount,2};
   if (error_rel<tolerance_qq)
      concurrences_cal(jcount)=real(trace(sig_zazb*rhoEvol(:,:,nlayers)));
      %error_sample(jcount)=states{jcount,2}-abs(error_qnn-real(trace(sig_zazb*rhoEvol(:,:,nlayers))));
      error_sample(jcount)=error_rel;
      break
   endif

   if (iter>iterations)
     disp('Relative error too big')
     iter=1;  % restart the loop of iterations
     nlayers=nlayers+1; % add one more layer
   endif

 end %  end of while over iterations

      %  /states{jcount,2};


endfor % End loop for the training states

    %%%%%%%  Let us calculate the  RMS error %%%%%%%%%%%%
rmserror=sqrt(sum(error_sample.^2)/sample);


save_results_quantum_nH(concurrences_cal,error_sample,learning_const_qq,deltat,input_file,rmserror,Ka,Kb,epsa,epsb,zeta,nu,state_file);

% save(state_file,'states')









