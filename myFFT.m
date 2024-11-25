function [ampSpec, f] = myFFT(x, Fs, plotFlag)

N = length(x) ;

modes  = (0:N-1)' ;  % equivalent to the index k from DFT

% find the number of the highest mode (account for odd/even N)
N_fastest = ceil((N+1)/2) ;

% set mode number of negative modes (account for odd/even N)
modes(N_fastest+1 : end) = modes(N_fastest+1 : end) - N ;
f = modes * Fs / N ;  % full frequency vector

Y = fft(x);               % at last, perform the Fast Fourier Transform

ampSpec = abs(Y/N);  % calc amplitude of entire Y and normalize by 1/N


% rearrage spectrum data for plotting the negative and positive modes
ind_for_plot = [ N_fastest+1:N , 1:N_fastest] ;

f = f(ind_for_plot) ;

% flip spectrum info such that amp at zero freq is in the middle
ampSpec    = ampSpec(ind_for_plot) ;

if exist('plotFlag')
    if plotFlag
        
        tvec = (0:N-1)/Fs ;
        spectrum_scale = 'linear' ;
        
        mycolors = lines(2);
        
        fig1 = figure('position',[80    65    1180    540]) ;
        subplot(2,1,1) ;
        plot(tvec, x,'.-','color',mycolors(1,:)) ;
        title('Signal in the time domain');
        xlabel('Time [sec]') ;
        ylabel('$x[t_n]$') ;
        grid on ; box on ;
        axis tight
        YL = ylim ; ylim(YL*1.1) ;
        
        subplot(2,1,2) ;
        plot(f, ampSpec, '.-', 'color', mycolors(2,:)) ;
        title('Double sided amplitude spectrum');
        xlabel('Frequency, $f$ [Hz]') ;
        ylabel('$\frac{1}{N} |X^d(f_k)|$') ;
        grid on ; box on ;
        axis tight
        set(gca,'yscale',spectrum_scale) ;
        YL = ylim ; ylim(YL*1.1) ;        
    end
end

return

