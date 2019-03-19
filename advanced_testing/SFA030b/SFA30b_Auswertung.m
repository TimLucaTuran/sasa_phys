%% load and plot measurement data from MF
% clear
load('SFA030b.mat')

Names = {'D100','D110','D121','D133'};
Pos = {'I','O'};
CO = {'k','r'};

% plot "amplitudes"
figure(101);clf
for n = 1:length(Names)
    subplot(2,2,n)
    for m = 1:2
        wl = Samples.([Names{n} Pos{m}]).Amp.wl*1e9;
        plot(wl, Samples.([Names{n} Pos{m}]).Amp.PX,['-' CO{m}],...
            wl, Samples.([Names{n} Pos{m}]).Amp.PY,[':' CO{m}],'Linewidth',1.5)
        hold all
    end
    plot( wl, Samples.FirstLayer.Amp.PX,['-b'],...
        wl, Samples.FirstLayer.Amp.PY,[':b'],'Linewidth',1.5)
    
    xlim([475 1200])
    ylim([0 1])
    title(Names{n})
    legend('both PX','both PY','1st PX','1st PY','2nd PX','2nd PY','Location','Southeast')
end

set(gcf,'PaperOrientation','landscape','PaperUnits','normalized','PaperPosition',[0 0 1 1])
% print('SFA030b-1 Amp.pdf','-dpdf')

% plot phases
figure(102)
for n = 1:length(Names)
    subplot(2,2,n)
    for m = 1:2
        wl = Samples.([Names{n} Pos{m}]).Arg.wl*1e9;
        plot(wl, Samples.([Names{n} Pos{m}]).Arg.PX,['-' CO{m}],...
            wl, Samples.([Names{n} Pos{m}]).Arg.PY,[':' CO{m}],'Linewidth',1.5)
        hold all
    end
    plot( wl, Samples.FirstLayer.Arg.PX,['-b'],...
        wl, Samples.FirstLayer.Arg.PY,[':b'],'Linewidth',1.5)
    
    xlim([475 1200])
    ylim([-0.7 1.2])
    title(Names{n})
    legend('both PX','both PY','1st PX','1st PY','2nd PX','2nd PY','Location','Southeast')
end

set(gcf,'PaperOrientation','landscape','PaperUnits','normalized','PaperPosition',[0 0 1 1])
% print('SFA030b-1 Arg.pdf','-dpdf','-painters')


%% wires

addpath ../../calc/utilities


% Smat of gold wires embedded in glass
load2       = 'Wire_single_1';
file_path_2 = [load2,'/',load2,'_Daten_gesamt.mat'];
clear SMAT_

load(file_path_2)
MAT2_ = SMAT_;

%% Pick a point in parameter space and plot with SFA30b-eguivalent

% S-matrix dimensions: 11 x 26 x 31 x 64 x 4 x 4
% Parameter space: Height_ x L_vec x W_vec x Lambda_

% simulation parameters
L_vec   = 170:195; % length in y-direction
W_vec   = 90:120; % width in x-direction
Height_ = 35:45; % metal layer thickness in nm
wavel2  = linspace(0.47,1.2,64)*1000;

% transmittance factor is equal to one since n_clad = n_subs = silica
INT_SMAT2 = abs(MAT2_).^2;

% choose point in parameter space
H_ind = 11;
L_ind = 1;
W_ind = 1;

% set plot vector
PT2Y = squeeze( INT_SMAT2(H_ind, L_ind, W_ind, :, 1, 1) );
PT2X = squeeze( INT_SMAT2(H_ind, L_ind, W_ind, :, 2, 2) );

figure(106)
clf
for n = 1:length(Names)
    subplot(2,2,n)
    
    wl = Samples.([Names{n} Pos{2}]).Amp.wl*1e9;
    p_hand = plot(wl, Samples.([Names{n} Pos{2}]).Amp.PX,'-b',...
        wl, Samples.([Names{n} Pos{2}]).Amp.PY,':b',...
        wavel2, PT2X,'-r',...
        wavel2, PT2Y,':r');
    hold all
    LabelStr.x = 'wavelength in nm';
    LabelStr.y = 'amplitude';
    omni_plot(p_hand,17,LabelStr)
    legend('measured AmpX','measured AmpY','simulated AmpX',...
        'simulated AmpY','Location','best')
    
    xlim([475 1200])
    ylim([0 1])
    title(Names{n})
end
hold off

%% plot measurements of specific fields (dosiges)

n = 4;
wl = Samples.([Names{n} Pos{2}]).Amp.wl*1e9;
freq_ = 300./wl * 1000;

LabelStr.x = 'frequency in THz';
LabelStr.y = 'amplitude';

figure(320)
clf
p_hand = plot(freq_, Samples.([Names{n} Pos{2}]).Amp.PX,...
    freq_, Samples.([Names{n} Pos{2}]).Amp.PY);

omni_plot(p_hand,17,LabelStr)
legend('wires AmpX','wires AmpY')

figure(321)
p_hand = plot(freq_, Samples.FirstLayer.Amp.PX, ...
    freq_, Samples.FirstLayer.Amp.PY);

omni_plot(p_hand,17,LabelStr)
legend('squares AmpX','squares AmpY','location','best')

figure(322)
clf
p_hand = plot(freq_, Samples.([Names{n} Pos{1}]).Amp.PX,...
    freq_, Samples.([Names{n} Pos{1}]).Amp.PY);

omni_plot(p_hand,17,LabelStr)
legend('stack AmpX','stack AmpY','location','best')


%% Export current figures to pdf and fig

expfigPath = '/Volumes/sperrhake/Research/Programming/MatLab/Calculation/calc/ExportFig';

addpath(expfigPath)

fighand  = gcf;
svname   = 'stack_amplitude_D133';
svpath   = ['/Users/jan/Documents/Figures/Simulation results/',...
    'SASA4Experiment/patch_wire_exp_rizzl/']; 
savepath = [svpath,svname];
saveas(fighand,[savepath,'.png']);
export_fig(savepath,'-pdf','-transparent',fighand);

disp('done');







