%% Multi-band filtering

function xf = multfilt(data,bandstopfrequencies,frequencyRes)

k = length(bandstopfrequencies)*4;

% set parameters

% generate cuts
count = 1;
for kk = 1:4:k
cuts = round(bandstopfrequencies(count),1);
fcuts(kk:kk+3) = [cuts-0.2 cuts-0.1 cuts+0.1 cuts+0.2];
count = count+1;
end

mags = [ones(1,length(bandstopfrequencies)); zeros(1,length(bandstopfrequencies))];
mags = [mags(:); 1]';

devs = [repmat(0.01,1,length(bandstopfrequencies)); repmat(0.05,1,length(bandstopfrequencies))];
devs = [devs(:); 0.01]';

% Kaiser window design
[n,Wn,beta,ftype] = kaiserord(fcuts,mags,devs,frequencyRes);          % Kaiser Window FIR Specification
n = n + rem(n,2);
hh = fir1(n,Wn,ftype,kaiser(n+1,beta),'noscale');   

% zero-phase digital filtering
xf = filtfilt(hh,1,data);

