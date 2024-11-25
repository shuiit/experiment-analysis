function [S, tvec, f] = mySpectrogram(x, Fs, W, plot_flag,db_flag, cutoff_db)
% assume length(W) is even
if mod(length(W),2)~=0
    error('Here we assume window length in even')
end

% divide x into blocks of size N=length(W)

N    = length(W) ;
Ntot = length(x) ;

t_block_start = 1: (N/2) : Ntot-N+1 ;
t_block_mid   = t_block_start + N/2 ;

% let ind_j = t_block_start(j)
% then the k'th block will be taken from:
% x( ind_j : ind_j + N - 1) ;

Nblocks = length(t_block_start) ;

dt = 1/Fs ; % time interval between adjacent samples

% allocate spectrogram matrix S
S = zeros(N, Nblocks) ;
tvec = t_block_mid * dt ;
tvec_tot = (0:Ntot-1)*dt ;

for j=1:Nblocks    
    ind_j = t_block_start(j) ;            % find where current block begins
    xblock = x( ind_j : ind_j + N - 1) ;  % cut current block from signal
    xblock = xblock .* W ;                % apply window
    [ampSpec, f] = myFFT(xblock, Fs) ;    % calculate FFT
    S(:,j) = ampSpec ;                    % store FFT of current window in S
end

if plot_flag
    % from here on we just plot ...

    mycolors = lines(2);
    fig = figure('position',[80    65    1180    680]) ;
    
    ax1 = subplot(4,1,1) ;
    plot(tvec_tot, x,'.-','color',mycolors(1,:)) ;
    title('Entire signal in the time domain');
    %xlabel('Time [sec]') ;
    ylabel('$x(t)$') ;
    grid on ; box on ;
    axis tight
    YL = ylim ; ylim(YL*1.1) ;
    p1 = get(gca,'position') ;
    
    ax2 = subplot(4,1,2:4);
    
    S_db = 20*log10(S) ; % convert to db
    
    if exist('cutoff_db','var')  % remove small entries in S (for better visualization)
        S_db(S_db<cutoff_db) = cutoff_db ;
    end
    
    if db_flag
        imagesc(tvec, f, S_db) ;
    else
        imagesc(tvec, f, S) ;
    end
    
    set(gca,'ydir','normal') ; % make vertical axis grow bottom up
    % title('Spectrogram in dB') ;
    cbar = colorbar ; %#ok<NASGU>
    p2 = get(gca,'position') ; %#ok<NASGU>
    xlabel('Time [sec]') ;
    ylabel('Frequency [Hz]') 
    
    p1(3) = 0.725 ; 
    set(ax1,'position',p1) ;
    
    linkaxes([ax1, ax2],'x') ;
end

return

