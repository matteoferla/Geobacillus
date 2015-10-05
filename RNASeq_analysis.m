%olocus=locus; %protect variables..
%osamples=samples;

%% complicated import
%truth table
%d=logical(descriptor);
%rows are aerobicity, richness and rapidity
%rpkm is a 12 x 4062 matrix
n_samples=12;
fr=1:n_samples; %change if less.
n_repl=2;
ri=[1 2; 3 4;5 6; 7 8; 9 10; 11 12];
n_cond=6;

%kill zeros.
samples=osamples; %unless subsetting
%nzi=geomean(orpkm,2)>0;
%rpkm=orpkm(nzi,fr);
%locus=olocus(nzi);
locus=olocus;
rpkm=cds;

%% scaling
%adapted from http://uk.mathworks.com/help/bioinfo/examples/identifying-differentially-expressed-genes-from-rna-seq-data.html
pseudo_ref_sample =geomean(rpkm,2); %rows
nzi = pseudo_ref_sample>0; % ignore genes with zero geometric mean
ratios = bsxfun(@rdivide, rpkm(nzi,:), pseudo_ref_sample(nzi)); %elementwise division
sizeFactors = median(ratios, 1);
base_rpkm = bsxfun(@rdivide,rpkm,sizeFactors);
figure;
subplot(2,1,1)
maboxplot(log2(rpkm),samples,'title','Raw Read Counts','orientation', 'horizontal');
subplot(2,1,2)
maboxplot(log2(base_rpkm),samples,'title','Size Factor Adjusted Read Counts','orientation', 'horizontal');
xlabel('log2 rpkm');

%checking scaling worked.
verif=zeros(n_samples,n_samples);
veriz=zeros(n_samples,n_samples);
for s=1:n_samples
    for t=1:n_samples
        verif(s,t)=mean(base_rpkm(nzi,s),1)/mean(base_rpkm(nzi,t),1)-sizeFactors(s)/sizeFactors(t);
        veriz(s,t)=mean(base_rpkm(:,s),1)/mean(base_rpkm(:,t),1)-sizeFactors(s)/sizeFactors(t);
    end
end
figure;
imagesc(verif);
colorbar
set(gca,'XTickLabel',samples)
set(gca,'YTickLabel',samples)
title('$$\frac{\sum_i k_{ij} / n}{\sum_i k_{ij"} / n} - \frac{s_j}{s_{j"}}$$','interpreter','latex');
set(gca,'XTick',1:n_samples)
set(gca,'YTick',1:n_samples)
%hold on
%plot(veriz,'.b');


%% Replicate quality
% unfortunately the replicates are not amazing

figure
for s=1:n_cond
    p=subplot(2,3,s)
    loglog([1,10^5],[1,10^5],'k', 'linewidth', 2); %epimirically determined.
    p.XTick=[1 10 100 1000 10000 100000];
    p.YTick=[1 10 100 1000 10000 100000];
    hold on
    loglog(base_rpkm(:,s*2-1),base_rpkm(:,s*2),'.b')
    xlim([1,10^5]);
    ylim([1,10^5]);
    title(cond(s))
    text(2,10^4,['\rho = ' num2str(corr(base_rpkm(:,s*2-1),base_rpkm(:,s*2)))])
    
end
suptitle('Replicate {\it vs.} replicate')

%% clustergram
% Log2 tranformations and the zscore of those isn't the most current way of
% doing things, but it is quick and cheerful.
% The data shows that 0.01 h-1 don't make much of a difference.

logmean=log2(base_rpkm+1); % log2(N+1) does not work well with rpkm as it is a ratio
norlog=zscore(logmean,0,2);
clustergram(norlog,'RowLabels',locus,'columnlabels',samples)

%% checking...
%[cond(1),cond(3)]
%mairplot(meanRS(nzi,1),meanRS(nzi,3),'Labels',label(nzi),'Factor',20)

%% Variance analysis
cond=cell(n_cond,1);
meanRS=zeros(size(rpkm,1),n_cond);
varRS=zeros(size(rpkm,1),n_cond);
noiseRS=zeros(size(rpkm,1),n_cond);
rawvarRS=zeros(size(rpkm,1),n_cond);
fitvarRS=zeros(size(rpkm,1),n_cond);
funRS=cell(n_cond,1);
for s=1:n_cond
    cond(s)=samples(s*2);
    meanRS(:,s)=mean(base_rpkm(:,ri(s,:)),2);
    varRS(:,s) = var(base_rpkm(:,ri(s,:)),0, 2);
    noiseRS(:,s) = meanRS(:,s) * mean(1./sizeFactors(ri(s,:)));
    funRS{s}=estimateNBVarFunc(meanRS(:,s),varRS(:,s),sizeFactors(ri(s,:)));
    rawvarRS(:,s)  =  funRS{s}(meanRS(:,s));
    fitvarRS(:,s)  =  rawvarRS(:,s) + noiseRS(:,s);
end
%plot
figure
for s=1:n_cond
    subplot(2,3,s)
    loglog(meanRS(:,s), varRS(:,s), '*')
    hold on
    loglog(meanRS(:,s), fitvarRS(:,s), '.r')
    ylabel('Base Variances')
    xlabel('Base Means')
    title(cond(s))
end
suptitle('Dependence of the Variance on the Mean for each sample pair')

%% ECDF
degrees_of_freedom = 1; %2-1
ratiovarRS = varRS ./ fitvarRS;
pchisq = chi2cdf(degrees_of_freedom * ratiovarRS, degrees_of_freedom);
count_levels = [0 5 10 25 50 100 200];
labels = {'0-5','6-10','11-25','26-50','51-100','101-200','> 201'};
figure
cm = jet(7);
grps=zeros(size(rpkm,1),n_cond);
tally=zeros(7,n_cond);
for s=1:n_cond
    grps(:,s) = sum(bsxfun(@ge,meanRS(:,s),count_levels),2); % stratification
    subplot(2,3,s)
    hold on
    for i = 1:7
        [Y1,X1] = ecdf(pchisq(grps(:,s)==i,s));
        plot(X1,Y1,'LineWidth',2,'color',cm(i,:))
    end
    plot([0,1],[0,1] ,'k', 'linewidth', 2)
    ax = gca;
    ax.Box = 'on';
    legend(labels,'Location','NorthWest')
    xlabel('\chi^2 probability of residual')
    ylabel('Empirical CDF')
    title(cond(s))
    tally(:,s)=accumarray(grps(:,s),1);
end
suptitle('Residuals ECDF plot for each sample pair')
hold off
% This makes me wonder...
figure
bar(tally)
title('Histogram of sizes (to be fixed)')
set(gca,'XTickLabel',cond)



%% Negative binomial
%For some reason computePVal is not liking me and gives NaN.
%Tried everything.
%0 * 0 due to 0=nbinpdf(k, r=(mu.^2)./(var-mu), p=mu ./var), which is the pdf of x given (mu.^2)./(var-mu) successes each with p= mu ./var
%Checking the values, for locus(1), 409 successes happen before 160 failures even though the chance of success is 0.14. the mu is 0.9k, while var is 6k.
%The spread is diabolical, but why is k smaller than  mu?
%k is the sum of the normalised counts rowwise per cond. eq 10 on Anders
%poolmean is q, eq 12
%eq 13
%I should give the R version a go.


poolmean=zeros(size(rpkm,1),n_cond,n_cond);
meank=zeros(size(rpkm,1),n_cond,n_cond);
vark=zeros(size(rpkm,1),n_cond,n_cond);
k=zeros(size(rpkm,1),n_cond);
for s=1:n_cond
    k(:,s)=sum(rpkm(:,ri(s,:)),2); %test statistic
    for t=1:n_cond
        %poolmean(:,s,t) = mean(rpkm(:,[s*2-1:s*2 t*2-1:t*2]),2);
        poolmean(:,s,t) = sum(base_rpkm(:,[ri(s,:) ri(t,:)]),2); %q in Anders and Huber.
        %meank(:,s,t) = poolmean(:,s,t) * sum(sizeFactors(s*2-1:s*2));
        meank(:,s,t)= poolmean(:,s,t) * sum(sizeFactors(ri(s,:)));
        vark(:,s,t) = meank(:,s,t) + funRS{s}(poolmean(:,s,t)) * sum(sizeFactors(ri(s,:)).^2);
    end
end

pvals=zeros(size(rpkm,1),n_cond,n_cond);
%nzp=zeros(sum(nzi),n_cond,n_cond);


for s=1:n_cond
    for t=1:n_cond
        pvals(:,s,t) = computePVal(k(:,s), meank(:,s,t), vark(:,s,t), k(:,t), meank(:,t,s), vark(:,t,s));
        %nzp(:,s,t) = computePVal(k(nzi,s), meank(nzi,s,t), vark(nzi,s,t), k(nzi,t), meank(nzi,t), vark(nzi,t));
    end
end

%[kA, muA, varA, kB, muB, varB]=[k(:,s), meank(:,s,t), vark(:,s,t), k(:,t), meank(:,t,s), vark(:,t,s)]

%false discovery rate (FDR) with the Benjamini-Hochberg procedure [7] using the mafdr

fold_change = mean_B ./ mean_A;

%% Silly way
p=ones(size(base_rpkm,1),n_cond);
h=ones(size(base_rpkm,1),n_cond);
for s=1:n_cond
    for t=1:n_cond
        for x=1:size(base_rpkm,1)
            [h(x),p(x)]=ttest2(base_rpkm(x,ri(s,:)),base_rpkm(x,ri(t)),'alpha',0.05/numel(locus));
        end
    end
end
%histogram(p,100)

%% Matlab R

figure
subplot(2,3,1);
maboxplot(log2(cds+1),samples,'title','Raw (shifted log)','orientation', 'horizontal');
subplot(2,3,2);
maboxplot(log2(orpkm+1),samples,'title','RPKM (shifted log)','orientation', 'horizontal');
subplot(2,3,3);
maboxplot(rld,samples,'title','regulirised log','orientation', 'horizontal');
subplot(2,3,4);
maboxplot(vsd,samples,'title','var-stabilised','orientation', 'horizontal');

