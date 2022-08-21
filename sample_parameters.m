% Fit the Pol II model to different groups of data and compute the pausing
% time

sample_names = {'xrn3','GEX1A','fcafpa','DMSO','Col','fcafpa_Col_up',...
    'xrn3_Col_up'};
n_samples = 7;
list_params_opt = cell(1,n_samples);
mat_pausing_time = zeros(n_samples,10);

% Read global data
intron_count = readtable('data/Gene_intronNum');
gene_bed_info = readtable("data/TAIR10.gene.bed","FileType","text");
gene_bed_info.gene_length = gene_bed_info.Var3 - gene_bed_info.Var2;
max_intron_count = max(intron_count{:,2});

% Set parameter bounds
params_lb = [1e-6 1e-4 1e3 1e-6 300];
params_ub = [1e2 1 1e6 1 350];

%Sample 50,000 sets of random parameters within the feasible region
n_rand_params = 50000; 
rand_matrix = lhsdesign(n_rand_params, 5);
list_params_kept = cell(n_samples, 10);
list_pausing_time = cell(n_samples, 10);
fobj = zeros(n_rand_params,1);
fobj_cutoff = 10000;
n_kept = zeros(n_samples, 10);

% Set parameters for DSA
DSA_options.start_beta = 1e-6;
DSA_options.end_beta = 1e6;
DSA_options.expected_fobj = 100;
DSA_options.mc_length = 500;
DSA_options.cooling_rate = 1.1;
DSA_options.visual_output = 0;

for n = 1:n_samples
    % Read coverage profile
    filename = strcat('data\',sample_names{n},'_ser2.PCG.IntronNum.Matrix');
    coverage_by_gene = readtable(filename,'FileType','text');
    nc = size(coverage_by_gene,2);
    % Process coverage profile
    average_coverage = zeros(max_intron_count+1,400);
    genes_with_coverage = coverage_by_gene{:,4};
    average_gene_length = zeros(max_intron_count+1,1);
    for i = 0:max_intron_count
        genes_subset = intron_count{intron_count{:,2} == i,1};
        coverage_subset = coverage_by_gene{ismember(genes_with_coverage,genes_subset),nc-399:nc};
        length_subset = gene_bed_info{ismember(gene_bed_info.Var4,genes_subset),'gene_length'};
        average_coverage(i+1,:) = mean(coverage_subset,'omitnan');
        average_gene_length(i+1,:) = mean(length_subset);
    end
    
    figure;
    for i = 1:10
        test_data = average_coverage(i,:);
        gene_length = average_gene_length(i);
        params_lb(5) = gene_length;
        params_ub(5) = gene_length*1.1;
        
        [params_opt, f_opt, ~, ~] = DSA(@(x)dynamics_residues(test_data,...
            gene_length, x), log10(params_lb(:)), log10(params_ub(:)), DSA_options);
        fobj_cutoff = f_opt*2;
        params_opt = params_opt';
        log_lb_local = max(params_opt - 0.2, log10(params_lb));
        log_ub_local = min(params_opt + 0.2, log10(params_ub));
                
        logparams_rand_samp = linear_scale(rand_matrix, log_lb_local, log_ub_local);
        for j = 1:n_rand_params
            residues = dynamics_residues(test_data, gene_length, logparams_rand_samp(j,:));
            fobj(j) = sumsqr(residues);
        end
        params_kept = logparams_rand_samp(fobj<fobj_cutoff,:);
        list_params_kept{n,i} = params_kept;
        n_kept(n,i) = sum(fobj<fobj_cutoff);
        P_sim = zeros(n_kept(n,i), 301);
        pausing_time = zeros(n_kept(n,i), 1);
        for j = 1:n_kept(n,i)
            params = 10.^params_kept(j,:);
            P_sim(j,:) = dynamics_complete(params, gene_length);
            pausing_time(j) = time_of_pausing(1,params(2),params(3),...
            params(5),gene_length);
        end
        list_pausing_time{n,i} = pausing_time;
        
        xcor_before_300 = (0:200)*gene_length/200;
        xcor_after_300 = gene_length + 10*(301:400) - 3000; 
        xcor_comb = [xcor_before_300 xcor_after_300];
        region_x = [xcor_comb xcor_comb(end:-1:1)];
        region_y = [min(P_sim) max(P_sim(:,end:-1:1))];
        
        
        subplot(2,5,i);
        fill(region_x, region_y, [0.7 0.7 0.7], 'EdgeColor', 'none');
        
        hold on;
        plot(xcor_comb, test_data(100:400));
        
        legend('Simulated','Actual');
        xlabel('Position[bp from TSS]');
        ylabel('Coverage');
        title(strcat(sample_names(n),sprintf(', #Exons = %d',i)));   
    end
end

%% Compare parameters across different exon numbers and experiments
% The same experiment, different exons
param_names = {'F','k','b','d','PAS'};
sample_names = {'xrn3','GEX1A','fcafpa','DMSO','Col','fcafpa-Col-up',...
    'xrn3-Col-up'};
p_anova_exon_numbers = ones(n_samples,6);
p_anova_experiments = ones(10,6);

figure;
for n = 1:n_samples
    subplot(2,4,n);
    %{
    for i = 1:5
        params_vec = [];
        params_label = [];
        for j = 1:10
            params_kept = list_params_kept{n,j};
            params_vec = [params_vec;params_kept(:,i)];
            params_label = [params_label;j*ones(n_kept(n,j),1)];
        end
        subplot(2,3,i);
        title(param_names{i});
        violinplot(10.^params_vec, params_label, 'ShowData', false);
        p_anova_exon_numbers(n,i) = anova1(10.^params_vec, params_label, 'off');
        ylabel(param_names{i});
        xlabel('#Exons');
        xtickangle(0);
        xlim([0 11]);
        box on;
    end
    %}
    pause_vec = [];
    pause_label = [];
    for j = 1:10
        pause_kept = list_pausing_time{n,j};
        pause_vec = [pause_vec; pause_kept];
        pause_label = [pause_label;j*ones(n_kept(n,j),1)];
    end
    %subplot(2,3,6);
    title(sample_names{n});
    violinplot(pause_vec, pause_label, 'ShowData', false);
    xlabel('#Exons');
    ylabel('Pause time');
    xtickangle(0);
    xlim([0 11]);
    box on;
    %p_anova_exon_numbers(n,6) = anova1(pause_vec, pause_label, 'off');
end
    
%% The same exon number, different experiments
figure;
for n = 1:10
    subplot(2,5,n);
    %{
    for i = 1:5
        params_vec = [];
        params_label = [];
        for j = 1:n_samples
            params_kept = list_params_kept{j,n};
            params_vec = [params_vec;params_kept(:,i)];
            params_label = [params_label;repmat(sample_names(j),n_kept(j,n),1)];
        end
        subplot(2,3,i);
        title(param_names{i});
        violinplot(10.^params_vec, params_label, 'ShowData', false);
        p_anova_experiments(n,i) = anova1(10.^params_vec, params_label, 'off');
        ylabel(param_names{i});
        xlabel('Experiment');
        xlim([0 n_samples+1]);
        box on;
    end
    %}
    pause_vec = [];
    pause_label = [];
    for j = 1:n_samples
        pause_kept = list_pausing_time{j,n};
        pause_vec = [pause_vec; pause_kept];
        pause_label = [pause_label;repmat(sample_names(j),n_kept(j,n),1)];
    end
    title(sprintf('#Exons = %d', n));
    violinplot(pause_vec, pause_label, 'ShowData', false);
    xlabel('Experiment');
    ylabel('Pause time');
    xlim([0 n_samples+1]);
    box on;
    %p_anova_experiments(n,6) = anova1(pause_vec, pause_label, 'off');
end


%% Plot parameters against number of exons

figure;
for i = 1:5
    subplot(2,3,i);
    hold on;
    for j = 1:n_samples
        params_opt_comb = list_params_opt{j};
        plot(1:10,10.^params_opt_comb(:,i), '-o');
    end
    legend(sample_names);
    title(param_names{i});
    xlabel('#Exons');
    ylabel(param_names{i});
    box on;
end
subplot(2,3,6);
hold on;
for i = 1:n_samples
    plot(1:10,mat_pausing_time(i,:),'-o');
end
legend(sample_names);
title('Pausing time');
xlabel('#Exons');
ylabel('Pausing time');
box on;

function dydx = dynamics_after_TES(x,y,params)
%Parameters: v0,k,b,d
v0 = 1;%params(1);
k = params(2);
b = params(3);
d = params(4);
TES = params(5);

v = v0*((x-TES)^2+k*b)/((x-TES)^2+b);
dvdx = 2*b*(1-k)*v0*(x-TES)/(((x-TES)^2+b)^2);
dydx = -(dvdx+d)*y/v;
end

function P = dynamics_complete(params, gene_length)
%Parameters: F, v0, k, b, d
F = params(1);
v0 = 1;
k = params(2);
b = params(3);
d = params(4);
TES = params(5);
%gene_length = params(6);
xcor_before_300 = (0:200)*gene_length/200;
xcor_after_300 = gene_length + 10*(301:400) - 3000; 
xcor_comb = [xcor_before_300 xcor_after_300];

ystart = F/(k*v0);
xcor_before_TES = xcor_comb(xcor_comb<TES);
xcor_after_TES = xcor_comb(xcor_comb>=TES);

if length(xcor_after_TES)>1
    [~,P_after_TES] = ode23s(@(x,y)dynamics_after_TES(x,y,params),xcor_after_TES,ystart);
elseif length(xcor_after_TES)==0
    P_after_TES = [];
else
    P_after_TES = ystart;
end
if length(xcor_before_TES)>0
    P_before_TES = steady_P_before_TES(xcor_before_TES,F,v0,k,b,TES);
else
    P_before_TES = [];
end
P = [P_before_TES(:);P_after_TES(:)];
if length(P) < 301 || (~isreal(sum(P)))
    P = 1e8*ones(301,1);
end
end

function r = dynamics_residues(data, gene_length, logparams)
params = 10.^logparams;
P = dynamics_complete(params, gene_length);
%length(P)
x = P(:);
y = data(100:400)';
%size(x)
%size(y)
if length(x) ~= length(y)
    size(x)
end
c = 1;%(x'*y)/(x'*x);
r = c*P(:) - data(100:400)';
end

function x_scale = linear_scale(x, lb, ub)
x_scale = x.*(ub - lb) + lb;
end
