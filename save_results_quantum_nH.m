  function save_results_quantum_nH(concurrences_cal,error_sample,learning_const_qq,totaltime,input_file,rmserror,Ka,Kb,epsa,epsb,zeta,nu,state_file)
    %colNames={'rmserror','Ka','Kb', 'epsa','epsb','zeta'};
    %T=table(Ka',Kb', epsa', epsb', zeta', 'VariableNames',colNames);
      time=datestr(now,'yyyy_mm_dd_HH_MM_SS');
      %namefile=strcat('parameters','_',num2str(totaltime),'_',num2str(nepoch),'_',num2str(nlayers),'_');
      %namefile='parameters_qn_';
      filename=strcat(input_file,time,'.mat');
      save(filename,'concurrences_cal','error_sample','learning_const_qq','totaltime','rmserror','Ka','Kb','epsa','epsb','zeta','nu','state_file')

    %writetable(T,filename)
    %finalresults=results(iepoch-1,:);





