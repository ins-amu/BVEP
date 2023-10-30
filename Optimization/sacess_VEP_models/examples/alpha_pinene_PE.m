%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TITLE: Thermal isomerization of alfa-pinene
% REF:   R.E. Fuguitt y J. E. Hawkins, 1947. Rate of thermal isomerization
% of ?-pinene in the liquid phase. J. A. C. S., 69:461
% 
% NOTE!!!: [] indicates that the corresponding input may be omitted,
%              default value will be assigned
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%============================
% RESULTS PATHS RELATED DATA
%============================
inputs.pathd.results_folder='alpha_pinene'; % Folder to keep results (in Results) for a given problem                       
inputs.pathd.short_name='pinene';                       % To identify figures and reports for a given problem 
                                                     % ADVISE: the user may introduce any names related to the problem at hand 
inputs.pathd.runident='r1';                         % [] Identifier required in order not to overwrite previous results
                                                     %    This may be modified from command line. 'run1'(default)
                                                     

%============================
% MODEL RELATED DATA
%============================
inputs.model.input_model_type='charmodelC';          % Model introduction: 'charmodelC'|'c_model'|'charmodelM'|'matlabmodel'|'sbmlmodel'|                     
                                                     %                     'blackboxmodel'|'blackboxcost                   
inputs.model.n_st=5;                                 % Number of states                                  
inputs.model.n_par=5;                                % Number of model parameters                                  
inputs.model.n_stimulus=0;                           % Number of inputs, stimuli or control variables   
inputs.model.names_type='standard';                    % [] Names given to states/pars/inputs: 'standard' (x1,x2,...p1,p2...,u1, u2,...) 
                                                     %                                       'custom'(default)
                              
inputs.model.eqns=...                                % Equations describing system dynamics. Time derivatives are regarded 'd'st_name''
            char('dx1=-(p1+p2)*x1',...  
                 'dx2= p1*x1',...
                 'dx3= p2*x1-(p3+p4)*x3+p5*x5',...
                 'dx4= p3*x3',...
                 'dx5= p4*x3-p5*x5');
   
 p1=5.93e-5;  p2=2.96e-5;  p3=2.05e-5;  p4=27.5e-5;  p5=4e-5;
             
 inputs.model.par=[p1 p2 p3 p4 p5];                   % Nominal value for the parameters, this allows to fix known parameters
                                                      % These values may be updated during optimization  
%==================================
% EXPERIMENTAL SCHEME RELATED DATA
%==================================
 inputs.exps.n_exp=1;                                 % Number of experiments                                                                  
 inputs.exps.n_obs{1}=4;                              % Number of observed quantities per experiment                         
 inputs.exps.obs{1}=char('y1=x1','y2=x2','y3=x3','y4=x5');
 inputs.exps.exp_y0{1}=[100 0 0 0 0];       % Initial conditions for each experiment       
 inputs.exps.t_f{1}=36420;                            % Experiments duration
 inputs.exps.n_s{1}=8;                                % Number of sampling times
 inputs.exps.t_s{1}=[1230 3060 4920 7800 10680 15030 22620 36420];                         % [] Sampling times, by default equidistant
                                                            
%==================================
% EXPERIMENTAL DATA RELATED INFO
%==================================
 inputs.exps.data_type='real';                         % Type of data: 'pseudo'|'pseudo_pos'|'real'             
 inputs.exps.noise_type='homo';
 inputs.exps.std_dev{1}=[0.001 0.001 0.001 0.001];
 
%Experimental data 1: 
 inputs.exps.exp_data{1}=[
	88.35 7.3   2.3  1.75
    76.4  15.6  4.5  2.8
    65.1  23.1  5.3  5.8
    50.4  32.9  6.0  9.3
    37.5  42.7  6.0  12.0
    25.9  49.1  5.9  17.0
    14.0  57.4  5.1  21.0
    4.5   63.1  3.8  25.7
   ];


  
 %==================================
 % UNKNOWNS RELATED DATA
 %==================================
 
 % GLOBAL UNKNOWNS (SAME VALUE FOR ALL EXPERIMENTS)
 
 inputs.PEsol.id_global_theta='all';               % 'all'|User selected 
 inputs.PEsol.global_theta_max=[1 1 1 1 1];        % Maximum allowed values for the paramters
 inputs.PEsol.global_theta_min= [0 0 0 0 0];       % Minimum allowed values for the parameters

 
 

 % % GLOBAL INITIAL CONDITIONS
 inputs.PEsol.id_global_theta_y0='x1';               % [] 'all'|User selected| 'none' (default)
 inputs.PEsol.global_theta_y0_max=[110];                % Maximum allowed values for the initial conditions
 inputs.PEsol.global_theta_y0_min=[90];                % Minimum allowed values for the initial conditions
 %inputs.PEsol.global_theta_y0_guess=[];              % [] Initial guess
% 
% % LOCAL UNKNOWNS (DIFFERENT VALUES FOR DIFFERENT EXPERIMENTS)
% 
% inputs.PEsol.id_local_theta{1}='none';                % [] 'all'|User selected| 'none' (default)
% % inputs.PEsol.local_theta_max{iexp}=[];              % Maximum allowed values for the paramters
% % inputs.PEsol.local_theta_min{iexp}=[];              % Minimum allowed values for the parameters
% % inputs.PEsol.local_theta_guess{iexp}=[];            % [] Initial guess
% inputs.PEsol.id_local_theta_y0{1}='none';             % [] 'all'|User selected| 'none' (default)
% % inputs.PEsol.local_theta_y0_max{iexp}=[];           % Maximum allowed values for the initial conditions
% % inputs.PEsol.local_theta_y0_min{iexp}=[];           % Minimum allowed values for the initial conditions
% % inputs.PEsol.local_theta_y0_guess{iexp}=[];         % [] Initial guess
 
 
 %==================================
 % COST FUNCTION RELATED DATA
 %==================================
          
 inputs.PEsol.PEcost_type='lsq';                       % 'lsq' (weighted least squares default) | 'llk' (log likelihood) | 'user_PEcost' 
 inputs.PEsol.lsq_type='Q_I';
% inputs.PEsol.llk_type='homo_var';                     % [] To be defined for llk function, 'homo' | 'homo_var' | 'hetero' 
 
 
 %==================================
 % NUMERICAL METHODS RELATED DATA
 %==================================
 %
 % SIMULATION
 inputs.ivpsol.ivpsolver='cvodes';                     % [] IVP solver: 'radau5'(default, fortran)|'rkf45'|'lsodes'|
 inputs.ivpsol.senssolver='cvodes';                    % [] Sensitivities solver: 'cvodes' (C)


 inputs.ivpsol.rtol=1.0D-7;                            % [] IVP solver integration tolerances
 inputs.ivpsol.atol=1.0D-7; 

inputs.model.input_model_type='charmodelC';              % Model type must be 'charmodelC' to use C
inputs.model.exe_type='fullC';                           % Generates gcc command
inputs.model.odes_file=fullfile(pwd,'amigoRHS.c');
inputs.nlpsol.cvodes_gradient=1;                         % To activate gradient evaluation
[inputs privstruct]=AMIGO_Prep(inputs);


%inputs.plotd.plotlevel='full';                        % [] Display of figures: 'full'|'medium'(default)|'min' |'noplot' 
%inputs.plotd.figsave=1;
% inputs.plotd.epssave=0;                              % [] Figures may be saved in .eps (1) or only in .fig format (0) (default)
% inputs.plotd.number_max_states=8;                    % [] Maximum number of states per figure
% inputs.plotd.number_max_obs=8;                       % [] Maximum number of observables per figure
% inputs.plotd.n_t_plot=100;                           % [] Number of times to be used for observables and states plots
% inputs.plotd.number_max_hist=8;                      % [] Maximum number of unknowns histograms per figure (multistart)
%inputs.plotd.nx_contour=100;                          % Number of points for plotting the contours x and y direction
%inputs.plotd.ny_contour=100;                          % ADVISE: >50
