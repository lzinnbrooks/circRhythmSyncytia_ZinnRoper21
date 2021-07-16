%Generates color map for two compartment division of labor simulations
clear;
close all;

%load data file (for nonuniform or uniform compartments)
load mRNAcounts_equalComp.mat

maxNum = max(max(numTrans));
countMx = zeros(maxNum+2,maxNum+2);
numRows = size(numTrans,1);

for i=1:numRows
    xval = numTrans(i,2);
    yval = numTrans(i,1);

    countMx(xval+1,yval+1) = countMx(xval+1,yval+1)+1;
end

pcolor(0:maxNum+1, 0:maxNum+1, countMx);
shading flat;
cmap = colormap('gray(100)');
cmap = 1-cmap;
colormap(cmap);
colorbar;
hold on
plot(0:maxNum+1,0:maxNum+1,'r','LineWidth',2)
axis square;
xticks([0.5:5:0.5+5*floor(maxNum/5)])
xticklabels({'0','5','10','15','20','25','30'})
yticks([0.5:5:0.5+5*floor(maxNum/5)])
yticklabels({'0','5','10','15','20','25','30'})

set(gca,'FontSize',16)
xlabel('\# of mRNAs transcribed (comp.\ 1)','interpreter','Latex')
ylabel('\# of mRNAs transcribed (comp.\ 2)','interpreter','Latex')
title('1:1','interpreter','Latex')