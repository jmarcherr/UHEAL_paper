function [est,coeff] = NNfloorEstim(X,f,foi,fitF,feedback)
%% Neural Noise Floor Estimator (NNfloorEstim) - Sam Watson - sadaw@dtu.dk - 12/09/2024
% Estimate the noise floor at a particular frequency(s) by fitting a c/(f^m) (aka 1/f) curve to a spectrum.
% Curve is fit by converting curve to log domain and solving for linear
% gradient and y-intercept: log10(P) = log(c) - m.log10(f)

% As such this is suited to neural noise which has an inherant distribution
% proportional to 1/f, particulalry when needing an estimate of the noose
% floor at low frequencies.

%Outputs
% est = noise floor power estimate (linear input units) at the requested foi
% m = estimated noise floor gradient coefficent
% c = estimated noise floor scaling coeffcient

%Inputs
% X = linear units POWER spectrum to be analysed/fit (positive frequencies only)
% f = corresponding frequency vector
% foi = exact frequency(s) of interest (must match frequency bins)
% fitF = exact frequencies to use for curve fitting (must match frequency bins)
% feedback = true/false to turn on/off plotting for review

%% Start function

%check inputs
validateattributes(X,{'double'},{'nonnegative','row'},'','power spectrum',1)
validateattributes(feedback,{'logical'},{'scalar'},'','feedback logical',5)

%check foi and fitF
if ~ all(ismember(foi,f))
    error("Frequencies in 'foi' do not match exact frequency bins")
elseif ~all(ismember(fitF,f))
    error("Frequencies in 'fitF' do not match exact frequency bins")
elseif any(ismember(foi,fitF))
    %warning("Frequencies to be predicted are included in fit data!")
end

%% Start fitting
%isolate fitting data
fitX = X(ismember(f, fitF));

%log transform
logfitX = log10(fitX);
logfitF = log10(fitF);

%linear param fit
coeff = polyfit(logfitF,logfitX,1); %robust least squares fit aross fitting data range

%use fit to predict values at foi
estlog = polyval(coeff,log10(foi));
est = 10.^estlog; %back to linear


if feedback
    modelx = log10(fitF);
    modely = coeff(1)*modelx + coeff(2);

    figure()
    tiledlayout(1,2)
    nexttile
    plot(f,X)
    xlabel("Frequency Hz")
    ylabel("Power \muV^2")
    title("Input Power Spectrum (Linear)")
    xlim([min(fitF),max(fitF)])

    nexttile

    plot(logfitF,logfitX,'--x') %fit data
    hold on
    plot(modelx,modely) %model line
    plot(log10(foi),estlog,'ro')
    xlabel("log10(Freq)")
    ylabel("log10(Power)")
    legend(["Input Data";"Model";"foi Estim."])
    title("log-log noise floor model")
    xlim(log10([min(fitF),max(fitF)]))
    
end


