function plotgenedistribution()
clc;
clear all;
close all;
methy_data_dir = 'methy_data/';
fig_dir = 'mat_plot/';
subdirs = {'hyper';'hypo';'std'};

cancer_names={'BRCA';'COAD';'LIHC';'LUAD';'LUSC'};
stage_names = {'all_stage';'normal';'i';'ii';'iii';'iv';'x'};
gene_names = {'ACVR1B';'APC'};%;'ARID1A';'ARID1B';'ARID2';'ASXL1';'ATM';'ATRX';'AXIN1';'B2M';'BAP1';'BCOR';'BRCA1';'BRCA2';'CASP8';'CDC73';'CDH1';'CDKN2A';'CEBPA';'CIC';'CREBBP';'CYLD';'DAXX';'EP300';'FBXW7';'FUBP1';'GATA1';'GATA3';'HNF1A';'KDM5C';'KDM6A';'MAP3K1';'MEN1';'MLH1';'MSH2';'MSH6';'NCOR1';'NF1';'NF2';'NOTCH1';'NOTCH2';'NPM1';'PAX5';'PBRM1';'PHF6';'PIK3R1';'PRDM1';'PTCH1';'PTEN';'RB1';'RNF43';'RUNX1';'SETD2';'SMAD2';'SMAD4';'SMARCA4';'SMARCB1';'SOCS1';'SOX9';'STAG2';'STK11';'TET2';'TNFAIP3';'TRAF7';'TP53';'TSC1';'VHL';'WT1'};

size_gene_names = size(gene_names);
size_cancer_names = size(cancer_names);
size_stage = size(stage_names);

fig_counter = 0;
for i=1: size_gene_names(1)
    gene_name = char(gene_names(i));
    fig_counter = fig_counter + 1;
    fig=figure(fig_counter);
    clf();
    hold on;
    for j=1:size_cancer_names(1)
        cancer_name = char(cancer_names(j));
        % 第一张图, 每个基因一个图, 每个癌症一个子图, 画散点和箱线图
        data1_path= [methy_data_dir, gene_name, '_xy_', cancer_name,'.dat'];
        data1 = load(data1_path);
        subplot(5,3,3*j-2);
        plot(data1(:, 1), data1(:, 2),'b.');
        set(gca, 'YTick', [0 0.2 0.4 0.6 0.8 1.0]);
        set(gca, 'XTick', [0 1 2 3 4 5 6]);
        axis([0 6.5 0.0 1.0]);
        xlabel('stage');
        ylabel('methy');
    end
    fig_save_path = [fig_dir,gene_name,'.eps'];
    exportfig(fig, fig_save_path , 'FontMode', 'fixed', 'color', 'cmyk', 'width',8,'height',10,'FontSize', 12,'Resolution',300,'LineWidth',0.5);
    close all;
end

