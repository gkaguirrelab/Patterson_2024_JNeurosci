% plot watson model filters

% adjustable model parameters
g = [1 2 3];
gs = [0.5 0.7 0.9];
t_c = [0.01 0.02 0.03];

% fixed indpendent model parameters
nc = 9;
ns = 10;

% Frequencies
w = [2 4 8 16 32 64];
wDelta = min(diff(log10(w)));
upScale = 10;
wFit = 10.^(log10(min(w))-wDelta+wDelta/upScale:wDelta/upScale:log10(max(w))+wDelta);

Lw=[0.5 1 2 3];

figure
for x=1:length(g)
    G = g(x);
    Gs = gs(2);
    tc = t_c(2);
    ts = tc.*1.33;
    
    Hc = (wFit.*1i*2*pi*tc+1).^-1*nc;
    Hs = Gs*((wFit.*1i*2*pi*ts+1).^-1*ns);
    y = real(G*(Hc-Hs));
    
    subplot(1,2,1)
    hold on
    plot(wFit,real(Hc),'-r','Linewidth',Lw(x))
    plot(wFit,real(Hs).*-1,'-b','Linewidth',Lw(x))
    ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.XScale='log';
    title('Varying Amplitude')
    ylabel('center and surround filters')
    xlabel('frequency [hz]')
    legend(num2str(G))

    subplot(1,2,2)
    hold on
    plot(wFit,y,'-m','Linewidth',Lw(x))
    ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.XScale='log';
    ylabel('center+surround')
    xlabel('frequency [hz]')
    legend(num2str(G))
end

figure
for x=1:length(gs)
    G = g(2);
    Gs = gs(x);
    tc = t_c(2);
    ts = tc.*1.33;
    
    Hc = (wFit.*1i*2*pi*tc+1).^-1*nc;
    Hs = Gs*((wFit.*1i*2*pi*ts+1).^-1*ns);
    y = real(G*(Hc-Hs));
    
    subplot(1,2,1)
    hold on
    plot(wFit,real(Hc),'-r','Linewidth',Lw(x))
    plot(wFit,real(Hs).*-1,'-b','Linewidth',Lw(x))
    ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.XScale='log';
    title('Varying Surround Amplitude')
    ylabel('center and surround filters')
    xlabel('frequency [hz]')

    subplot(1,2,2)
    hold on
    plot(wFit,y,'-m','Linewidth',Lw(x))
    ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.XScale='log';
    ylabel('center+surround')
    xlabel('frequency [hz]')
end

figure
for x=1:length(t_c)
    G = g(2);
    Gs = gs(2);
    tc = t_c(x);
    ts = tc.*1.33;
    
    Hc = (wFit.*1i*2*pi*tc+1).^-1*nc;
    Hs = Gs*((wFit.*1i*2*pi*ts+1).^-1*ns);
    y = real(G*(Hc-Hs));
    
    subplot(1,2,1)
    hold on
    plot(wFit,real(Hc),'-r','Linewidth',Lw(x))
    plot(wFit,real(Hs).*-1,'-b','Linewidth',Lw(x))
    ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.XScale='log';
    title('Varying Time Constant')
    ylabel('center and surround filters')
    xlabel('frequency [hz]')

    subplot(1,2,2)
    hold on
    plot(wFit,y,'-m','Linewidth',Lw(x))
    ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.XScale='log';
    ylabel('center+surround')
    xlabel('frequency [hz]')
end