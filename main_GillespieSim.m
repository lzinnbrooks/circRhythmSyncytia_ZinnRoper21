%Run Gillespie algorithm for various number of nuclear compartments (main code)
close all;
clear;

totTime = 2200; %simulation time (h)
dt = 0.001; %time step (h)
tvec = 0:dt:totTime;
tsteps = length(tvec);

T_0 = 23.2; %autonomous period
T = 23.2; %free-running period (23.2 is default)
theta = T_0/T;

burst = 0; %toggle bursty translation (all mRNAs translate protein at once - not used in final writeup)
randProt = 0; %toggle random protein redistribution 
timecourse = 1; %toggle plotting of mRNA and protein timecourses
nonUniform = 0; %toggle non-uniform compartment sizes
spectrumPlot = 1; %toggle power spectrum plot

trials = 100; %number of simulations
nComps = 2; %number of compartments (we tested 1, 2, 4, 8, and 16)
compLength = 5; %average compartment length (mcm)
minLength = 1; %minimum compartment length (mcm)

if nonUniform == 0
    compBounds = 0:1/nComps:1; %compartment boundaries for uniform case
else
    lengthVec = minLength+exprnd(compLength-minLength,nComps,1); %exponentially distributed random comp. lengths
                                                                 %(not used in final writeup)
    compBounds = [0;cumsum(lengthVec)/sum(lengthVec)]; %compartment boundaries for nonuniform case
end

alpha_0 = 18; %*[1; 10; 100; 1000; 10000]; %transcription rate (in dark) h^(-1) - can also test a set of trans. rates
V_n = 0.1; %nuclear volume (pL)
V_c = 2; %cytoplasmic volume (pL)
beta = 10; %translation rate h^(-1)
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

%spectral analysis
spectrumMx{length(alpha_0),1} = [];
%quality factor (measure of limit cycle quality)
qFactor1 = zeros(length(alpha_0),1);

%smoothFn (smooth data to calculate power spectrum)
tskip = 100;
tvec2 = tvec(1:tskip:round(0.9*tsteps));
tsteps2 = length(tvec2);
halfPt = floor(tsteps2/2);
smoothFn = zeros(length(tvec2),trials);

%minimum protein count
minPR = zeros(nComps,1);

%loop through different transcription rates
for i=1:length(alpha_0)
    %i
    
    %compartment 1 data for each trial (used for power spectrum)
    comp1Time{trials} = [];
    comp1Data{trials} = [];
    
    for z=1:trials
        z
        
        %initialize mRNA and protein counts in each compartment
        MnData{nComps} = []; McData{nComps} = [];
        PcData{nComps} = []; PnData{nComps} = [];
        
        %initialize time vectors
        tMn{nComps} = []; tMc{nComps} = [];
        tPc{nComps} = []; tPn{nComps} = [];
        
        t = 0;
        
        for c=1:nComps
            %initial conditions from Wang and Peskin (2018)
            MnData{c} = round(10.9*V_n); McData{c} = round(1.88*V_c);
            PcData{c} = round(129.21*V_c); PnData{c} = round(1000*V_n);
            
            tMn{c} = 0; tMc{c} = 0; tPc{c} = 0; tPn{c} = 0;
        end
        
        while t<totTime   
            for c=1:nComps
                %compute probability of transcription (from Wang and Peskin, 2018)
                minPR(c) = min(PnData{c}(end),r);
                
                for j=0:minPR(c)
                    temp(j+1,c) = Kv^j*nchoosek(PnData{c}(end),j)*factorial(r)/factorial(r-j);
                end
                
                pAlpha(c) = alpha_0(i)*1/sum(temp(:,c));
                
                MnExport(c) = MnData{c}(end)*eta; %number of mRNAs exported from nucleus
                McDecay(c) = McData{c}(end)*eta; %number of mRNAs that degrade in cytoplasm
                
                %bursty translation (all protein translated at once)
                if (burst == 1)
                    PcBorn(c) = beta(i);
                else
                    PcBorn(c) = McData{c}(end)*beta(i); %number of proteins translated
                end
                
                PcImport(c) = PcData{c}(end)*eta; %number of proteins imported into nucleus
                PnDecay(c) = PnData{c}(end)*eta; %number of proteins that degrade in nucleus
            end
            
            %setup for Gillespie (reaction weights)
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
                
                %increment number of nuclear mRNAs in chosen compartment
                MnData{compInd}(end+1) = MnData{compInd}(end)+1;
                
                %record time that reaction occurred
                tMn{compInd}(end+1) = t;
                
            %mRNA export
            elseif randNo>w1/w && randNo<w2/w
                %which comp. does mRNA export occur in?
                compInd = find(mExp>rand,1);
                
                %decrement number of nuclear mRNAs in chosen compartment
                MnData{compInd}(end+1) = MnData{compInd}(end)-1;
                %increment number of cytoplasmic mRNAs in chosen compartment
                McData{compInd}(end+1) = McData{compInd}(end)+1;
                
                tMn{compInd}(end+1) = t;
                tMc{compInd}(end+1) = t;
                
            %mRNA decay
            elseif randNo>w2/w && randNo<w3/w
                %which comp. does mRNA decay occur in?
                compInd = find(mDec>rand,1);
                
                %decrement number of cytoplasmic mRNAs in chosen compartment
                McData{compInd}(end+1) = McData{compInd}(end)-1;
                tMc{compInd}(end+1) = t;
                
            %protein translation
            elseif randNo>w3/w && randNo<w4/w
                %which comp. does translation occur in?
                compInd = find(pTrans>rand,1);
                
                %bursty translation
                if(burst == 1)
                    for c=1:nComps
                        if c==compInd
                            %all mRNAs translate proteins at once
                            PcData{c}(end+1) = PcData{c}(end)+McData{c}(end);
                        else
                            PcData{c}(end+1) = PcData{c}(end);
                        end
                        
                        tPc{c}(end+1) = t;
                    end
                else
                    for c=1:nComps
                        if c==compInd
                            %increment number of cytoplasmic proteins in chosen compartment
                            PcData{c}(end+1) = PcData{c}(end)+1;
                        else
                            PcData{c}(end+1) = PcData{c}(end);
                        end
                        
                        tPc{c}(end+1) = t;
                    end
                end
                
                for c=1:nComps
                    %count total number of proteins in cell
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
                    evenProt = floor(totProt/nComps); 
                    extraProt = totProt-nComps*evenProt; %remaining protein when proteins cannot divide evenly
                    
                    if extraProt>0
                        extraInd = randsample(nComps,extraProt); %randomly redistribute extra protein
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
                        %decrement number of cytoplasmic proteins in chosen compartment
                        PcData{c}(end+1) = PcData{c}(end)-1;
                        %increment number of nuclear proteins in chosen compartment
                        PnData{c}(end+1) = PnData{c}(end)+1;
                        tPn{c}(end+1) = t;
                    else
                        PcData{c}(end+1) = PcData{c}(end);
                    end
                    
                    tPc{c}(end+1) = t;
                end
                
                for c=1:nComps
                    %update total protein count
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
                
                %decrement number of nuclear proteins in chosen compartment
                PnData{compInd}(end+1) = PnData{compInd}(end)-1;
                tPn{compInd}(end+1) = t;
            end
            
            temp = [];
            
            totProt = 0;
        end
        
        %save data from compartment 1
        if nonUniform == 0
            comp1Time{z} = tPn{1}; comp1Data{z} = PnData{1};
        else
            [m,indx] = max(lengthVec);
            comp1Time{z} = tPn{indx}; comp1Data{z} = PnData{indx};
        end
    end
    
    %generate smoothFn (uses equal time steps so that power spectrum can be calculated)
    for j=2:length(tvec2)
        for z=1:trials
            closeInd = find(comp1Time{z}>j*tskip*dt,1);
            distAhead = comp1Time{z}(closeInd)-j*tskip*dt;
            distBehind = j*tskip*dt-comp1Time{z}(closeInd-1);
            slope = (comp1Data{z}(closeInd)-comp1Data{z}(closeInd-1))/(distAhead+distBehind);
            smoothFn(j,z) = slope*distBehind + comp1Data{z}(closeInd-1);
        end
    end
    
    x = smoothFn(halfPt:end,:);   avgPn = mean(mean(x));
    tNew = tvec2(halfPt:end);
    
    meanPSD = 0;
    
    for z_0 = 1:trials
        %compute power spectrum for each simulation
        [f,power] = powerspectrum(x(:,z_0),tNew);
        df = f(2)-f(1);
        meanPSD = meanPSD+power;
    end
    
    %average power spectrum across all trials
    meanPSD = meanPSD/trials;
    spectrumMx{i,1} = meanPSD;
    
    %calculate percentage of spectrum within 2 hrs of peak period (quality factor)
    [maxP,indP] = max(spectrumMx{i,1});
    peakPer = 1/f(indP); %peak period (h)
    lower = 1/(peakPer+2);  upper = 1/(peakPer-2); %upper and lower frequency bounds
    freqStep = (upper-lower)/20;
    freqInt = (lower:freqStep:upper)';
    yVal = interp1q(f',spectrumMx{i,1},freqInt); %linear interpolation between freq. bounds
    int = trapz(freqInt,yVal); %total area between frequency bounds
    
    %total area
    totInt = trapz(f,spectrumMx{i,1});
    
    %quality factor
    qFactor1(i,1) = int/totInt;
end

%plot time courses for each state variable in each compartment
if timecourse == 1
    figure(1)
    
    subplot(2,2,1)
    for j=1:nComps
        plot(tMn{j}, MnData{j}, 'LineWidth', 1)
        hold on
    end
    
    xlabel('Time ($h$)', 'Interpreter', 'Latex')
    %ylabel('number of mRNAs in nucleus')
    ylabel('$m_n$', 'Interpreter', 'Latex')
    xlim([0 totTime])
    set(gca,'FontSize',16)
    
    subplot(2,2,2)
    for j=1:nComps
        plot(tMc{j}, McData{j}, 'LineWidth', 1)
        hold on
    end
    
    xlabel('Time ($h$)', 'Interpreter', 'Latex')
    ylabel('$m_c$', 'Interpreter', 'Latex')
    %ylabel('number of mRNAs in cytoplasm')
    xlim([0 totTime])
    set(gca,'FontSize',16)
    
    subplot(2,2,3)
    for j=1:nComps
        plot(tPc{j}, PcData{j}, 'LineWidth', 1)
        hold on
    end
    
    xlabel('Time ($h$)', 'Interpreter', 'Latex')
    %ylabel('number of proteins in cytoplasm')
    ylabel('$p_c$', 'Interpreter', 'Latex')
    xlim([0 totTime])
    set(gca,'FontSize',16)
    
    subplot(2,2,4)
    for j=1:nComps
        plot(tPn{j}, PnData{j}, 'LineWidth', 1.5)
        hold on
    end
    
    xlabel('Time ($h$)', 'Interpreter', 'Latex')
    %ylavel('number of proteins in nucleus')
    ylabel('$p_n$', 'Interpreter', 'Latex')
    xlim([0 totTime])
    set(gca,'FontSize',16)
    
    title = strcat('Time courses for mRNA and protein');
    sgtitle(title, 'Interpreter', 'Latex')
end


%plot power spectrum
if spectrumPlot == 1
    figure(2)
    
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