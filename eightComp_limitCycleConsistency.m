%Limit cycle consistency in an 8-compartment syncytium (with and without protein sharing)
%See main_GillespieSim.m for fully commented code (most of the steps are the same here)
close all;
clear;

%load cycles2.mat;

%run for 5 circadian "days"
totTime = 220;
dt = 0.001;
tvec = 0:dt:totTime;
tsteps = length(tvec);

T_0 = 23.2; %autonomous period
T = 23.2; %free-running period (23.2 is default)
theta = T_0/T;

%all of these toggles are set to zero for these simulations
burst = 0; %toggle bursty translation
powerSpec = 0; %toggle computing power spectrum
randProt = 0; %toggle random protein redistribution 
timecourse = 0; %toggle plotting of mRNA and protein timecourses
nonUniform = 0; %toggle non-uniform compartments
spectrumPlot = 0; %toggle power spectrum plot
corrCalc = 0; %calculate correlation coefficient

nComps = 1; %number of compartments (use 1 for independent compartments, 8 for protein sharing case)
trials = 8; %number of trials (use 8 for independent compartment, 1 for protein sharing case)

compLength = 5; %average compartment length (mcm)
minLength = 1; %minimum compartment length (mcm)

if nonUniform == 0
    compBounds = 0:1/nComps:1; %compartment boundaries 
else
    lengthVec = [2;1];
    compBounds = [0;cumsum(lengthVec)/sum(lengthVec)];
end

alpha_0 = 18;%*[1; 10; 100; 1000; 10000]; %transcription rate (in dark) h^(-1)
V_n = 0.1; %nuclear volume (pL)
V_c = 2; %cytoplasmic volume (pL)
beta = theta*10; %translation rate h^(-1)
eta = 2*pi/22; %import/decay rates h^(-1)
gamma_m = theta*eta; %mRNA export rate h^(-1)
gamma_p = theta*eta; %protein import rate h^(-1)
delta_m = theta*eta; %degradation of mRNA (in cytoplasm) h^(-1)
delta_p = theta*eta; %degradation of protein (in nucleus) h^(-1)
K = 200; %equilibrium constant pL^(-1)
r = 5; %number of protein binding sites
temp = [];
Kv = 1/(K*V_n);

totProt = 0; %total number of cytoplasmic proteins

oscLeft = alpha_0/V_n*beta;
oscRight = gamma_m^2*K*4*5^5;

tic

%spectral analysis
spectrumMx{length(alpha_0),nComps} = [];
spectrumTot{length(alpha_0)} = [];

%measures of limit cycle quality
qFactor1 = zeros(length(alpha_0),nComps);
qFactorTot = zeros(length(alpha_0),1);

%smoothFn
tskip = 100;
tvec2 = tvec(1:tskip:round(0.9*tsteps));
tsteps2 = length(tvec2);
halfPt = floor(tsteps2/2);
smoothFn = zeros(length(tvec2),trials,nComps);

minPR = zeros(nComps,1);

for i=1:length(alpha_0)
    %i
    
    %compartment 1 data for each trial
    compTime{trials,nComps} = [];
    compData{trials,nComps} = [];
    
    for z=1:trials
        z
        
        %initialize mRNA and protein counts in each compartment
        MnData{nComps} = []; McData{nComps} = [];
        PcData{nComps} = []; PnData{nComps} = [];
        
        %initialize time vectors
        tMn{nComps} = []; tMc{nComps} = [];
        tPc{nComps} = []; tPn{nComps} = [];
        tTrans{nComps} = []; %times that transcription events occurred
        
        t = 0;
        
        for c=1:nComps
            %trough
            MnData{c} = 4;%round(10.9*V_n); 
            McData{c} = 1;%round(1.88*V_c);
            PcData{c} = 17;%round(129.21*V_c); 
            PnData{c} = 8;%round(1000*V_n);
            
            %peak
%             MnData{c} = 1;
%             McData{c} = 2;
%             PcData{c} = 50; 
%             PnData{c} = 60;
            
            %MnData{c} = MnStart; McData{c} = McStart;
            %PcData{c} = PcStart; PnData{c} = PnStart;
            
            tMn{c} = 0; tMc{c} = 0; tPc{c} = 0; tPn{c} = 0; 
        end
        
        while t<totTime   
            for c=1:nComps
                minPR(c) = min(PnData{c}(end),r);
                
                for j=0:minPR(c)
                    temp(j+1,c) = Kv^j*nchoosek(PnData{c}(end),j)*factorial(r)/factorial(r-j);
                end
                
                pAlpha(c) = alpha_0(i)*1/sum(temp(:,c));
                
                MnExport(c) = MnData{c}(end)*eta; %number of mRNAs exported from nucleus
                McDecay(c) = McData{c}(end)*eta; %number of mRNAs that degrade in cytoplasm
                
                if (burst == 1)
                    PcBorn(c) = beta;
                else
                    PcBorn(c) = McData{c}(end)*beta; %number of proteins translated
                end
                
                PcImport(c) = PcData{c}(end)*eta; %number of proteins imported into nucleus
                PnDecay(c) = PnData{c}(end)*eta; %number of proteins that degrade in nucleus
            end
            
            w1 = sum(pAlpha);
            w2 = sum(pAlpha+MnExport);
            w3 = sum(pAlpha+MnExport+McDecay);
            w4 = sum(pAlpha+MnExport+McDecay+PcBorn);
            w5 = sum(pAlpha+MnExport+McDecay+PcBorn+PcImport);
            w = sum(pAlpha+MnExport+McDecay+PcBorn+PcImport+PnDecay);
            
            dtau = -log(rand)/w;
            t = t+dtau;
            
            randNo = rand;
            
            alph = cumsum(pAlpha)/sum(pAlpha);
            mExp = cumsum(MnExport)/sum(MnExport);
            mDec = cumsum(McDecay)/sum(McDecay);
            pTrans = cumsum(PcBorn)/sum(PcBorn);
            pImp = cumsum(PcImport)/sum(PcImport);
            pDec = cumsum(PnDecay)/sum(PnDecay);
            
            %transcription
            if randNo<w1/w
                %which comp. does transcription occur in?
                compInd = find(alph>rand,1);
                
                MnData{compInd}(end+1) = MnData{compInd}(end)+1;
                tMn{compInd}(end+1) = t;
                tTrans{compInd}(end+1) = t;
                
            %mRNA export
            elseif randNo>w1/w && randNo<w2/w
                %which comp. does mRNA export occur in?
                compInd = find(mExp>rand,1);
                
                MnData{compInd}(end+1) = MnData{compInd}(end)-1;
                McData{compInd}(end+1) = McData{compInd}(end)+1;
                tMn{compInd}(end+1) = t;
                tMc{compInd}(end+1) = t;
                
            %mRNA decay
            elseif randNo>w2/w && randNo<w3/w
                %which comp. does mRNA decay occur in?
                compInd = find(mDec>rand,1);
                
                McData{compInd}(end+1) = McData{compInd}(end)-1;
                tMc{compInd}(end+1) = t;
                
            %protein translation
            elseif randNo>w3/w && randNo<w4/w
                %which comp. does translation occur in?
                compInd = find(pTrans>rand,1);
                
                if(burst == 1)
                    for c=1:nComps
                        if c==compInd
                            PcData{c}(end+1) = PcData{c}(end)+McData{c}(end);
                        else
                            PcData{c}(end+1) = PcData{c}(end);
                        end
                        
                        tPc{c}(end+1) = t;
                    end
                else
                    for c=1:nComps
                        if c==compInd
                            PcData{c}(end+1) = PcData{c}(end)+1;
                        else
                            PcData{c}(end+1) = PcData{c}(end);
                        end
                        
                        tPc{c}(end+1) = t;
                    end
                end
                
                for c=1:nComps
                    totProt = totProt+PcData{c}(end);
                end
                
                %random protein redistribution
                if randProt==1
                    randVec = rand(totProt,1);
                    
                    for c=1:nComps
                        PcData{c}(end) = sum(randVec>compBounds(c) & randVec<compBounds(c+1));
                    end
                %uniform protein redistribution
                else
                    evenProt = floor(totProt/nComps); extraProt = totProt-nComps*evenProt;
                    if extraProt>0
                        extraInd = randsample(nComps,extraProt);
                    else
                        extraInd = 0;
                    end
                    
                    for c=1:nComps
                        PcData{c}(end) = evenProt+max(c==extraInd);
                    end
                end
                
            %protein import
            elseif randNo>w4/w && randNo<w5/w
                %which comp. does import occur in?
                compInd = find(pImp>rand,1);
                
                for c=1:nComps
                    if c==compInd
                        PcData{c}(end+1) = PcData{c}(end)-1;
                        PnData{c}(end+1) = PnData{c}(end)+1;
                        tPn{c}(end+1) = t;
                    else
                        PcData{c}(end+1) = PcData{c}(end);
                    end
                    
                    tPc{c}(end+1) = t;
                end
                
                for c=1:nComps
                    totProt = totProt+PcData{c}(end);
                end
                
                %random protein redistribution
                if randProt==1
                    randVec = rand(totProt,1);
                    
                    for c=1:nComps
                        PcData{c}(end) = sum(randVec>compBounds(c) & randVec<compBounds(c+1));
                    end
                %uniform protein redistribution
                else
                    evenProt = floor(totProt/nComps); extraProt = totProt-nComps*evenProt;
                    if extraProt>0
                        extraInd = randsample(nComps,extraProt);
                    else
                        extraInd = 0;
                    end

                    for c=1:nComps
                        PcData{c}(end) = evenProt+max(c==extraInd);
                    end
                end
                
            %protein decay
            else
                %which comp. does protein decay occur in?
                compInd = find(pDec>rand,1);
                
                PnData{compInd}(end+1) = PnData{compInd}(end)-1;
                tPn{compInd}(end+1) = t;
            end
            
            temp = [];
            
            totProt = 0;
        end
        
        for c=1:nComps
            compTime{z,c} = tPn{c}; compData{z,c} = PnData{c};
        end
    end
    
    smoothFn(1,:,:) = PnData{1}(1);
    
    for c=1:nComps
        for j=2:length(tvec2)
            for z=1:trials
                closeInd = find(compTime{z,c}>j*tskip*dt,1);
                distAhead = compTime{z,c}(closeInd)-j*tskip*dt;
                distBehind = j*tskip*dt-compTime{z,c}(closeInd-1);
                slope = (compData{z,c}(closeInd)-compData{z,c}(closeInd-1))...
                    /(distAhead+distBehind);
                smoothFn(j,z,c) = slope*distBehind + compData{z,c}(closeInd-1);
            end
        end
    end
    
    if (nComps==1)
        comp1func = sum(smoothFn,2); %sum over trials
    else
        comp1func = sum(smoothFn,3); %sum over compartments
    end
    
    comp1smooth = smooth(comp1func);
    tSmooth = tvec2(1:5:end);
    comp1smooth2 = comp1smooth(1:5:end);
    %find peaks in nuclear protein expression (defines bounds of circadian days)
    [pks,locs] = findpeaks(comp1smooth2);
    actLocs = 0.5*(locs-1);
    
    rng = 30;
    newLocs = [];
    maxVal = [];
    
    for g = 1:(length(locs)-1)
        leftInd = max(1,locs(g)-rng);
        rightInd = min(length(comp1smooth2),locs(g)+rng);
        subVec = comp1smooth2(leftInd:rightInd);
        newLocs(g) = (leftInd-1)+round(mean(find(subVec == max(subVec))));
        maxVal(g) = max(subVec);
    end
    
    actNewLocs = 0.5*(newLocs-1);
    
    if powerSpec == 1
        
        y = zeros(length(smoothFn(halfPt:end,1,1)),trials);
        
        for c=1:nComps
            x = smoothFn(halfPt:end,:,c);   avgPn = mean(mean(x));
            y = y + smoothFn(halfPt:end,:,c);
            tNew = tvec2(halfPt:end);
            
            meanPSD = 0;
            
            for z_0 = 1:trials
                [f,power] = powerspectrum(x(:,z_0),tNew);
                df = f(2)-f(1);
                meanPSD = meanPSD+power;
            end
            
            meanPSD = meanPSD/trials;
            spectrumMx{i,c} = meanPSD;
            
            %calculate percentage of spectrum within 2 hrs of peak period
            [maxP,indP] = max(spectrumMx{i,c});
            peakPer = 1/f(indP); %peak period (h)
            lower = 1/(peakPer+2);  upper = 1/(peakPer-2); %upper and lower frequency bounds
            freqStep = (upper-lower)/20;
            freqInt = (lower:freqStep:upper)';
            yVal = interp1q(f',spectrumMx{i,c},freqInt); %linear interpolation between freq. bounds
            int = trapz(freqInt,yVal); %total area between frequency bounds
            
            %total area
            totInt = trapz(f,spectrumMx{i,c});
            
            %quality factor
            qFactor1(i,c) = int/totInt;
        end
        
        %compute quality factor for total protein
        meanPSD = 0;
        
        for z_0 = 1:trials
            [f,power] = powerspectrum(y(:,z_0),tNew);
            df = f(2)-f(1);
            meanPSD = meanPSD+power;
        end
        
        meanPSD = meanPSD/trials;
        spectrumTot{i} = meanPSD;
        
        %calculate percentage of spectrum within 2 hrs of peak period
        [maxP,indP] = max(spectrumTot{i});
        peakPer = 1/f(indP); %peak period (h)
        lower = 1/(peakPer+2);  upper = 1/(peakPer-2); %upper and lower frequency bounds
        freqStep = (upper-lower)/20;
        freqInt = (lower:freqStep:upper)';
        yVal = interp1q(f',spectrumTot{i},freqInt); %linear interpolation between freq. bounds
        int = trapz(freqInt,yVal); %total area between frequency bounds
        
        %total area
        totInt = trapz(f,spectrumTot{i});
        
        %quality factor
        qFactorTot(i) = int/totInt;
    end
end

%qVec = mean(qFactor1,2)  %mean quality factor for each \alpha

toc

if timecourse == 1
    figure(1)
    for j=1:nComps
        plot(tMn{j}, MnData{j}, 'LineWidth', 1)
        hold on
    end
    
    xlabel('time (h)')
    ylabel('number of mRNAs in nucleus')
    xlim([0 totTime])
    set(gca,'FontSize',16)
    
    figure(2)
    for j=1:nComps
        plot(tMc{j}, McData{j}, 'LineWidth', 1)
        hold on
    end
    
    xlabel('time (h)')
    ylabel('number of mRNAs in cytoplasm')
    xlim([0 totTime])
    set(gca,'FontSize',16)
    
    figure(3)
    for j=1:nComps
        plot(tPc{j}, PcData{j}, 'LineWidth', 1)
        hold on
    end
    
    xlabel('time (h)')
    ylabel('number of proteins in cytoplasm')
    xlim([0 totTime])
    set(gca,'FontSize',16)
    
    figure(4)
    for j=1:nComps
        plot(tPn{j}, PnData{j}, 'LineWidth', 1.5)
        hold on
    end
    
    %plot(actLocs,pks,'k*','MarkerSize',10)
    
    xlabel('Time ($h$)', 'Interpreter', 'Latex')
    ylabel('Number of proteins in nucleus $P_n$', 'Interpreter', 'Latex')
    %title('Deterministic vs. stochastic results - $\alpha = 180$', 'Interpreter', 'Latex')
    xlim([0 totTime])
    set(gca,'FontSize',16)
end

if spectrumPlot == 1 && powerSpec == 1
    for i=1:length(alpha_0)
        plot(f,spectrumMx{i,1},'LineWidth',2)
        hold on
    end
    xlim([0 0.1])
    
    xlabel('Frequency $\omega$ ($h^{-1}$)', 'Interpreter', 'Latex')
    ylabel('Power $P(\omega)$', 'Interpreter', 'Latex')
    set(gca,'FontSize',16)
    
    qFactor1
end

actNewLocs = unique(actNewLocs);
locDiff = diff(actNewLocs);
maxVal = comp1smooth2(2*actNewLocs+1);

%clean up peak data (remove false peaks)
while min(locDiff)<8
    for i=1:length(locDiff)
        if locDiff(i)<8%(actNewLocs(i+1)-actNewLocs(i)) < 5
            if maxVal(i+1) > maxVal(i)
                actNewLocs(i) = actNewLocs(i+1);
            elseif maxVal(i) > maxVal(i+1)
                actNewLocs(i+1) = actNewLocs(i);
            else
                actNewLocs(i) = (actNewLocs(i)+actNewLocs(i+1))/2;
            end
        end
    end
    
    actNewLocs = unique(actNewLocs);
    locDiff = diff(actNewLocs);
    maxVal = comp1smooth2(2*actNewLocs+1);
end

actNewLocs = unique(actNewLocs);
maxVal = comp1smooth2(2*actNewLocs+1);

periodVec = diff(actNewLocs);
maxPer = ceil(max(periodVec));
edges = 0:maxPer;

if (actNewLocs(1) < 5)
    actNewLocs(1) = [];
    maxVal(1) = [];
end

%plot total nuclear protein expression (first five peak times marked with asterixes)
figure(5)
plot(tvec2,comp1func,'LineWidth',2)
hold on
plot(actNewLocs,maxVal,'r*','MarkerSize',20)

%display peak times
sprintf('%.1f ',actNewLocs(1:6))