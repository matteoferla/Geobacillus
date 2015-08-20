%%Is there any trend in the density of genes?
%I am assuming the genes are too close to make much sense, so there might
%be more data in the spaces.
%data from length_table python function.

%imported CSV via GUI

%% Which loci overlap

inter=start(2:end,1)-stop(1:end-1,1);
inter_start=start(2:end,1);
len=stop-start;
len2=circshift(len,1);

q=inter<0;
[sqi,i]=sort(inter(q));
qlc=locus(q);
sqlc=qlc(i);
qln=len(q);
sqln=qln(i);
qln2=len2(q);
sqln2=qln2(i);
table(abs(sqi),sqln,sqln2, 'Rownames',sqlc,'VariableNames',{'overlap' 'length' 'next_length'})

%% Gene density

g_len=3859975; %genome size.
bp=zeros(g_len,1);
for i=(1:numel(start))
    bp(start(i):stop(i))=bp(start(i):stop(i))+1; %Matlab, why no autoincrementor?!
end

figure;
plot(1:g_len,bp,'.b','MarkerSize',0.1)
hold on;
w=50000;
wl=floor(length(bp)/w);
wbp = mean(reshape(bp(1:end-mod(end,w)),[w,wl]))';
plot(linspace(1,g_len,numel(wbp)),wbp,'-r')
axis([0 g_len 0 3.1])
ylabel('Coding density')
xlabel('position bp')
title('Coding density in G. thermoglucosidans')
hold off;

%% plotted Spaces between genes across the genome.


windowSize = 10;
b = (1/windowSize)*ones(1,windowSize);
a = 2;

windowSize2 = 100;
b2 = (1/windowSize2)*ones(1,windowSize2);

figure
hold on;
plot(inter_start/1000,inter/1000,'.r');
plot(inter_start/1000,filter(b,a, inter/1000),'-b');
plot(inter_start/1000,filter(b2,a, inter/1000),'-c');
n=10;
plot(inter_start/1000,polyval(polyfit(inter_start,inter,n),inter_start)/1000,'-k');
xlabel('position(kb)');
ylabel('intergenic space(kb)');
title('G. thermoglucosidius intergenic spaces');
legend('Raw','MM10','MM100','polyfit10');
hold off;

%% Plotted on circle
%well... kind of. Hack using polar plot.
figure

inter_plot=polar(inter_start/g_len*2*pi,1+inter/g_len*500,'.r');
hold on;
genome=polar(linspace(0,pi*2,360),ones(1,360),'-k');
%workaround for bug in polar linespec
set(genome,'LineWidth',3)

%theta and rho
inter_plot=polar(inter_start/g_len*2*pi,1+inter/g_len*500,'.r');
inter_plot2=polar(inter_start/g_len*2*pi,1+filter(b,a,inter)/g_len*500,'-b');
inter_plot3=polar(inter_start/g_len*2*pi,1+filter(b2,a,inter)/g_len*500,'-c');
plot(inter_start/1000,1+polyval(polyfit(inter_start,inter,n),inter_start)/g_len*500,'-k');
hold off;
%axis([-1.5,1.5,-1.5,1.5])
view([90 -90])
xlabel('position(kb)');
ylabel('intergenic space(kb)');
title('G. thermoglucosidius intergenic spaces');
legend('Raw','Genome','Raw','MM10','MM100','polyfit'); %polar workaround is messy.
