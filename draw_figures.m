%% Compare pausing time across different exon numbers, col-0
n = 5; %Col-0 is the 5-th experiment
pause_vec = [];
pause_label = [];
for i = 1:10
    pause_vec = [pause_vec;list_pausing_time{n,i}];
    n_add = length(list_pausing_time{n,i});
    pause_label = [pause_label;i*ones(n_add,1)];
end
blue = [68 114 196]/255;
orange = [237 125 49]/255;
figure;
violinplot(pause_vec, pause_label, 'ShowData', false,...
    'Violincolor', blue);
xlim([0 11]);
xlabel('#Exons');
ylabel('Pausing time');
box on;
title('Pausing time for genes with different exon numbers');

%% Compare pausing time between two different conditions: col-0 vs xrn3
n1 = 5;
n2 = 1;
name1 = 'Col-0';
name2 = 'xrn3';
pause_vec = [];
pause_label = [];
labels = cell(1,20);
for i = 1:10
    p1_add = list_pausing_time{n1,i};
    p2_add = list_pausing_time{n2,i};
    nadd_1 = length(p1_add);
    nadd_2 = length(p2_add);
    pause_vec = [pause_vec;p1_add;p2_add];
    labels{2*i-1} = sprintf('%s, #exons=%d', name1, i);
    labels{2*i} = sprintf('%s, #exons=%d', name2, i);
    pause_label = [pause_label;(2*i-1)*ones(nadd_1,1);(2*i)*ones(nadd_2,1)];
end
figure;colororder([blue;orange]);
violinplot(pause_vec, pause_label, 'ShowData', false);
xlim([0 21]);
ylabel('Pausing time');
box on;
title('Pausing time, Col-0 vs xrn3');
xticks(1.5:2:19.5);
xticklabels(1:10);
xlabel('#Exons');

%% Compare pausing time between two different conditions: col-0 vs fcafpa
n1 = 5;
n2 = 3;
name1 = 'Col-0';
name2 = 'fcafpa';
pause_vec = [];
pause_label = [];
labels = cell(1,20);
purple = [204 102 255]/255;
p_ranksum = zeros(10,1);
for i = 1:10
    p1_add = list_pausing_time{n1,i};
    p2_add = list_pausing_time{n2,i};
    p_ranksum(i) = ranksum(p1_add, p2_add, "tail", "right");
    nadd_1 = length(p1_add);
    nadd_2 = length(p2_add);
    pause_vec = [pause_vec;p1_add;p2_add];
    labels{2*i-1} = sprintf('%s, #exons=%d', name1, i);
    labels{2*i} = sprintf('%s, #exons=%d', name2, i);
    pause_label = [pause_label;(2*i-1)*ones(nadd_1,1);(2*i)*ones(nadd_2,1)];
end
figure;colororder([blue;purple]);
violinplot(pause_vec, pause_label, 'ShowData', false);
xlim([0 21]);
ylabel('Pausing time');
box on;
title('Pausing time, Col-0 vs fcafpa');
xticks(1.5:2:19.5);
xticklabels(1:10);
xlabel('#Exons');

%% Compare pausing time between two different conditions: DMSO vs GEX1A
n1 = 4;
n2 = 2;
name1 = 'DMSO';
name2 = 'GEX1A';
pause_vec = [];
pause_label = [];
labels = cell(1,20);
grey = [0.5 0.5 0.5];
green = [112 173 71]/255;
for i = 1:10
    p1_add = list_pausing_time{n1,i};
    p2_add = list_pausing_time{n2,i};
    nadd_1 = length(p1_add);
    nadd_2 = length(p2_add);
    p_ranksum(i) = ranksum(p1_add, p2_add, "tail", "right");
    pause_vec = [pause_vec;p1_add;p2_add];
    labels{2*i-1} = sprintf('%s, #exons=%d', name1, i);
    labels{2*i} = sprintf('%s, #exons=%d', name2, i);
    pause_label = [pause_label;(2*i-1)*ones(nadd_1,1);(2*i)*ones(nadd_2,1)];
end
figure;colororder([grey;green]);
violinplot(pause_vec, pause_label, 'ShowData', false);
xlim([0 21]);
ylabel('Pausing time');
box on;
title('Pausing time, DMSO vs GEX1A');
xticks(1.5:2:19.5);
xticklabels(1:10);
xlabel('#Exons');