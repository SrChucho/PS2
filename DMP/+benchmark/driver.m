%
% Status : main Dynare file
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

tic0 = tic;
% Define global variables.
global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info estimation_info ys0_ ex0_
options_ = [];
M_.fname = 'benchmark';
M_.dynare_version = '5.5';
oo_.dynare_version = '5.5';
options_.dynare_version = '5.5';
%
% Some global variables initialization
%
global_initialization;
M_.exo_names = cell(1,1);
M_.exo_names_tex = cell(1,1);
M_.exo_names_long = cell(1,1);
M_.exo_names(1) = {'ez'};
M_.exo_names_tex(1) = {'ez'};
M_.exo_names_long(1) = {'ez'};
M_.endo_names = cell(6,1);
M_.endo_names_tex = cell(6,1);
M_.endo_names_long = cell(6,1);
M_.endo_names(1) = {'S'};
M_.endo_names_tex(1) = {'S'};
M_.endo_names_long(1) = {'S'};
M_.endo_names(2) = {'lambdaw'};
M_.endo_names_tex(2) = {'lambdaw'};
M_.endo_names_long(2) = {'lambdaw'};
M_.endo_names(3) = {'lambdaf'};
M_.endo_names_tex(3) = {'lambdaf'};
M_.endo_names_long(3) = {'lambdaf'};
M_.endo_names(4) = {'kappa'};
M_.endo_names_tex(4) = {'kappa'};
M_.endo_names_long(4) = {'kappa'};
M_.endo_names(5) = {'z'};
M_.endo_names_tex(5) = {'z'};
M_.endo_names_long(5) = {'z'};
M_.endo_names(6) = {'theta'};
M_.endo_names_tex(6) = {'theta'};
M_.endo_names_long(6) = {'theta'};
M_.endo_partitions = struct();
M_.param_names = cell(8,1);
M_.param_names_tex = cell(8,1);
M_.param_names_long = cell(8,1);
M_.param_names(1) = {'beta'};
M_.param_names_tex(1) = {'beta'};
M_.param_names_long(1) = {'beta'};
M_.param_names(2) = {'sigma'};
M_.param_names_tex(2) = {'sigma'};
M_.param_names_long(2) = {'sigma'};
M_.param_names(3) = {'gamma'};
M_.param_names_tex(3) = {'gamma'};
M_.param_names_long(3) = {'gamma'};
M_.param_names(4) = {'b'};
M_.param_names_tex(4) = {'b'};
M_.param_names_long(4) = {'b'};
M_.param_names(5) = {'rhoz'};
M_.param_names_tex(5) = {'rhoz'};
M_.param_names_long(5) = {'rhoz'};
M_.param_names(6) = {'sz'};
M_.param_names_tex(6) = {'sz'};
M_.param_names_long(6) = {'sz'};
M_.param_names(7) = {'alpha'};
M_.param_names_tex(7) = {'alpha'};
M_.param_names_long(7) = {'alpha'};
M_.param_names(8) = {'B'};
M_.param_names_tex(8) = {'B'};
M_.param_names_long(8) = {'B'};
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 1;
M_.endo_nbr = 6;
M_.param_nbr = 8;
M_.orig_endo_nbr = 6;
M_.aux_vars = [];
M_ = setup_solvers(M_);
M_.Sigma_e = zeros(1, 1);
M_.Correlation_matrix = eye(1, 1);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = true;
M_.det_shocks = [];
M_.surprise_shocks = [];
M_.heteroskedastic_shocks.Qvalue_orig = [];
M_.heteroskedastic_shocks.Qscale_orig = [];
options_.linear = false;
options_.block = false;
options_.bytecode = false;
options_.use_dll = false;
M_.orig_eq_nbr = 6;
M_.eq_nbr = 6;
M_.ramsey_eq_nbr = 0;
M_.set_auxiliary_variables = exist(['./+' M_.fname '/set_auxiliary_variables.m'], 'file') == 2;
M_.epilogue_names = {};
M_.epilogue_var_list_ = {};
M_.orig_maximum_endo_lag = 1;
M_.orig_maximum_endo_lead = 1;
M_.orig_maximum_exo_lag = 0;
M_.orig_maximum_exo_lead = 0;
M_.orig_maximum_exo_det_lag = 0;
M_.orig_maximum_exo_det_lead = 0;
M_.orig_maximum_lag = 1;
M_.orig_maximum_lead = 1;
M_.orig_maximum_lag_with_diffs_expanded = 1;
M_.lead_lag_incidence = [
 0 2 8;
 0 3 0;
 0 4 0;
 0 5 0;
 1 6 0;
 0 7 0;]';
M_.nstatic = 4;
M_.nfwrd   = 1;
M_.npred   = 1;
M_.nboth   = 0;
M_.nsfwrd   = 1;
M_.nspred   = 1;
M_.ndynamic   = 2;
M_.dynamic_tmp_nbr = [6; 0; 0; 0; ];
M_.model_local_variables_dynamic_tt_idxs = {
};
M_.equations_tags = {
  1 , 'name' , 'S' ;
  2 , 'name' , '2' ;
  3 , 'name' , 'lambdaw' ;
  4 , 'name' , 'lambdaf' ;
  5 , 'name' , 'kappa' ;
  6 , 'name' , 'theta' ;
};
M_.mapping.S.eqidx = [1 3 5 6 ];
M_.mapping.lambdaw.eqidx = [1 3 4 ];
M_.mapping.lambdaf.eqidx = [4 ];
M_.mapping.kappa.eqidx = [3 5 6 ];
M_.mapping.z.eqidx = [1 2 ];
M_.mapping.theta.eqidx = [4 6 ];
M_.mapping.ez.eqidx = [2 ];
M_.static_and_dynamic_models_differ = false;
M_.has_external_function = false;
M_.state_var = [5 ];
M_.exo_names_orig_ord = [1:1];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(6, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(1, 1);
M_.params = NaN(8, 1);
M_.endo_trends = struct('deflator', cell(6, 1), 'log_deflator', cell(6, 1), 'growth_factor', cell(6, 1), 'log_growth_factor', cell(6, 1));
M_.NNZDerivatives = [18; -1; -1; ];
M_.static_tmp_nbr = [6; 0; 0; 0; ];
M_.model_local_variables_static_tt_idxs = {
};
load dynareinput; 
M_.params(1) = xpar(1);
beta = M_.params(1);
M_.params(2) = xpar(2);
sigma = M_.params(2);
M_.params(3) = xpar(3);
gamma = M_.params(3);
M_.params(4) = xpar(4);
b = M_.params(4);
M_.params(5) = xpar(5);
rhoz = M_.params(5);
M_.params(6) = xpar(6);
sz = M_.params(6);
M_.params(7) = xpar(7);
alpha = M_.params(7);
M_.params(8) = xpar(8);
B = M_.params(8);
%
% SHOCKS instructions
%
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = 1;
%
% INITVAL instructions
%
options_.initval_file = false;
oo_.steady_state(1) = xsteady(1);
oo_.steady_state(2) = xsteady(2);
oo_.steady_state(3) = xsteady(3);
oo_.steady_state(6) = xsteady(4);
oo_.steady_state(4) = xsteady(5);
oo_.steady_state(5) = 1;
if M_.exo_nbr > 0
	oo_.exo_simul = ones(M_.maximum_lag,1)*oo_.exo_steady_state';
end
if M_.exo_det_nbr > 0
	oo_.exo_det_simul = ones(M_.maximum_lag,1)*oo_.exo_det_steady_state';
end
resid; 
steady;
options_.nograph = true;
options_.order = 1;
var_list_ = {};
[info, oo_, options_, M_] = stoch_simul(M_, options_, oo_, var_list_);


oo_.time = toc(tic0);
disp(['Total computing time : ' dynsec2hms(oo_.time) ]);
if ~exist([M_.dname filesep 'Output'],'dir')
    mkdir(M_.dname,'Output');
end
save([M_.dname filesep 'Output' filesep 'benchmark_results.mat'], 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'benchmark_results.mat'], 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'benchmark_results.mat'], 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'benchmark_results.mat'], 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'benchmark_results.mat'], 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'benchmark_results.mat'], 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'benchmark_results.mat'], 'oo_recursive_', '-append');
end
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
