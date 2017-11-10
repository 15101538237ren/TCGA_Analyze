function plotgenedistribution()
clc;
clear all;
close all;
methy_data_dir = 'methy_data/';
fig_dir = 'mat_plot/';
subdirs = {'hyper';'hypo';'std'};

cancer_names={'BRCA';'COAD';'LIHC';'LUAD';'LUSC'};
cancer_markers = {'r-s';'g-p';'b-o';'k-x';'c-*'};
mean_std_names = {'mean';'std'};
stage_names = {'normal';'i';'ii';'iii';'iv'};%;'x'
gene_names = {'APC';}%'GATA3';'PAX5';'WT1'};%'ACVR1B';'APC';'ARID1A';'ARID1B';'ARID2';'ASXL1';'ATM';'ATRX';'AXIN1';'B2M';'BAP1';'BCOR';'BRCA1';'BRCA2';'CASP8';'CDC73';'CDH1';'CDKN2A';'CEBPA';'CIC';'CREBBP';'CYLD';'DAXX';'EP300';'FBXW7';'FUBP1';'GATA1';'GATA3';'HNF1A';'KDM5C';'KDM6A';'MAP3K1';'MEN1';'MLH1';'MSH2';'MSH6';'NCOR1';'NF1';'NF2';'NOTCH1';'NOTCH2';'NPM1';'PAX5';'PBRM1';'PHF6';'PIK3R1';'PRDM1';'PTCH1';'PTEN';'RB1';'RNF43';'RUNX1';'SETD2';'SMAD2';'SMAD4';'SMARCA4';'SMARCB1';'SOCS1';'SOX9';'STAG2';'STK11';'TET2';'TNFAIP3';'TRAF7';'TP53';'TSC1';'VHL';'WT1'};

size_gene_names = size(gene_names);
size_cancer_names = size(cancer_names);
size_stage = size(stage_names);
size_ms_names = size(mean_std_names);

rows = 5;
cols = 7;
dx=0.01;
fig_counter = 0;

alw = 0.75;    % AxesLineWidth
fsz = 14;      % Fontsize
lw = 1.5;      % LineWidth
msz = 8;       % MarkerSize

for i=1: size_gene_names(1)
    gene_name = char(gene_names(i));
    fig_counter = fig_counter + 1;
    fig=figure(fig_counter);
    clf();
    hold on;
    % 第一张图, 每个基因一个图, 每个癌症一个子图, 画散点和箱线图
    for j=1:size_cancer_names(1)
        cancer_name = char(cancer_names(j));
        data1_path= [methy_data_dir, gene_name, '_xy_', cancer_name,'.dat'];
        data1 = load(data1_path);
        subplot(rows,cols,cols*(j-1)+1);
        plot(data1(:, 1), data1(:, 2),'b.');
        set(gca, 'YTick', [0 0.2 0.4 0.6 0.8 1.0]);
        set(gca, 'XTick', [1 2 3 4 5 6]);
        set(gca,'XTicklabel',{'n';'i';'ii';'iii';'iv';'x'});
        axis([0 6.5 0.0 1.0]);
        title(cancer_name);
    end
    hold on;
    for j=1:size_cancer_names(1)
        cancer_name = char(cancer_names(j));
        data2_path= [methy_data_dir, gene_name, '_all_stage_', cancer_name,'.dat'];
        data2 = load(data2_path);
        subplot(rows,cols,cols*(j-1)+2);
        [u,x]=hist(data2, 0:dx:1);
        plot(x,u/(dx*sum(u)));
        box on;
        set(gca, 'XTick', [0 0.2 0.4 0.6 0.8 1.0]);
        %axis([0 6.5 0.0 1.0]);
        title('Distribution');
    end
    
    for j=1:size_stage(1)
        hold on;
        stage_name = char(stage_names(j));
        for k=1:size_cancer_names(1)
            cancer_name = char(cancer_names(k));
            data3_path= [methy_data_dir, gene_name, '_', stage_name, '_', cancer_name,'.dat'];
            data3 = load(data3_path);
            subplot(rows,cols,cols*(k-1)+2+j);
            hold on;
            [u,x]=hist(data3,0:dx:1);
            bar(x,u/(dx*sum(u)));
            set(gca, 'XTick', [0 0.2 0.4 0.6 0.8 1.0]);
            xlim([0 1]);
            ylim([0 60]);
            title(stage_name);
        end
    end
    set(fig,'paperunits','centimeter')
    set(fig,'papersize',[80,50])
    set(fig,'paperposition',[0 0 80 50]);
    set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties
    fig_save_path = [fig_dir,gene_name,'.pdf'];
    print(fig,[gene_name,'.pdf'],'-dpdf','-opengl','-r300');
    close all;
    %exportfig(fig, fig_save_path , 'FontMode', 'fixed', 'color', 'cmyk', 'width',32,'height',10,'FontSize', 12,'Resolution',300,'LineWidth',0.5);
    %
end

