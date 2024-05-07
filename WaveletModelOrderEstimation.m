clc;
clear all;
% basic commands to clear the workspace

% ------------------------ synthetic signal with noise ----------------- %

% time samples - 251 samples with 12.5 Hz for 20 seconds
t=0:0.08:20;

% signal samples for Model order 4 
% data = 0.021*exp(-0.02*t).*cos(2*pi*0.2*t+ pi*0.15)+0.022*exp(-0.15*t).*cos(2*pi*0.30*t+pi*0.5)+0.06*exp(-1.2*t).*cos(2*pi*1.50*t+ 0*pi)+0.06*exp(-2.50*t).*cos(2*pi*1.80*t+ 0.5*pi) ;

% Model order 6 
data = 0.4*exp(-0.11*t).*cos(2*pi*0.25*t+ pi*0.167)+ 0.2*exp(-3.40*t).*cos(2*pi*2.90*t+ pi*0.833)+ 0.5*exp(-0.17*t).*cos(2*pi*0.4*t+ 0.333*pi)+ 0.09*exp(-0.22*t).*cos(2*pi*0.6*t+ 0.7*pi)+ 0.02*exp(-0.31*t).*cos(2*pi*1*t+ pi*0.6)+ 0.9*exp(-0.73*t).*cos(2*pi*1.8*t);

% Model Order 2
% data= 0.002*cos(2*pi*0.25*t + 1.5*pi).*exp(0.19*t) +0.02*cos(2*pi*0.73*t + 0.5*pi).*exp(-0.15*t);

% transpose the row vector of 251 signal samples
x=data';

% assuming eigenValues snr of 40 and converting into absolute noise
SNR_dB=35;
SNR_abs=10^(SNR_dB/10);

% sum of squares / length of the vector - (known as the mean square value of the signal)
Px=sum(x.^2)/length(x);

% noise variance 
noiseVariance=Px/SNR_abs;

% rms noise (standard deviation of the noise)
noiseStandardDeviation=sqrt(noiseVariance);

% randn function generates matrix with r x c dimensions
% here we are generating eigenValues row vector of length of x and add it to the
% signal
noise=noiseStandardDeviation.*randn(length(x),1);
x_noise=x+noise;

% plotting the signal 
% plot(t,x_noise);
% xlabel('Time');  % Label for x-axis
% ylabel('Signal amplitude');  % Label for autoCorrelationMatrix-axis
% title('Input Signal Plot');  % Title of the plot

% ------------------------ Denoising the hankel matrix ----------------- %

% constructing the hankel matrix - dimensions are 126 x 126
% ind_mat = hankel(1:126 , 126:251); 
% signalClusterIndex = x_noise(ind_mat);

% decompose the formed Hankel matrix to Sparse + Low-Rank Matrix
% [L, S] = DenoisingALOHA(signalClusterIndex);

% % Denoising the signal using wavelet denoising
waveletSignalDenoiser(x_noise);

x_noise = wdenoise(x_noise,2, ...
    Wavelet='coif5', ...
    DenoisingMethod='Bayes', ...
    ThresholdRule='Median', ...
    NoiseEstimate='LevelIndependent');
% 
% % Constructing the hankel matrix 
ind_mat = hankel(1:126 , 126:251); 
H = x_noise(ind_mat);

% ------------------------ Estimating the model order ----------------- %

% construct the autocorrelation matrix and find the eigenvalues
autoCorrelationMatrix = (H'*H)/126;
unorderedEigenValues = abs(eig(autoCorrelationMatrix));

% sort the eigen values in decending order
eigenValues = sort (unorderedEigenValues, 'descend');

tempOrderArray=[];
for j=1:10
    
    trace = sum(eigenValues);
    
    % eigenvalues - the mean of the eigenvalues
    meanEigenValues = eigenValues- mean(eigenValues);


    % medoid vector is eigenValues 1D vector in which all the elements range from 1
    % to k(2) and medoid(i) denotes that the ith point is part of the
    % medoid(i)th cluster.

    % centroid(i)th element is the ith cluster
    [Medoid,Centroid]=kmedoids(meanEigenValues,2);
    signalClusterIndex=max(Centroid); 
    index = find(Centroid==signalClusterIndex);

% ------------------ PAM clustering iterations ------------------ %

    % if the eigenvalue belong to the higher cluster append in signalSubspaceCluster
    % signalSubspaceCluster denotes the eigenvalues in the higher cluster
    signalSubspaceCluster=[];
    for i=1:length(meanEigenValues)
        if Medoid(i)==index
            signalSubspaceCluster=[signalSubspaceCluster;meanEigenValues(i)];
        end
    end

    % size of the eigenvalues in higher cluster
    tempModelOrder=length(signalSubspaceCluster);
    
    % remove the eigenvalues in the signal subspace and partition again
    if tempModelOrder ~= length(meanEigenValues) && tempModelOrder ~= 0
      eigenValues = eigenValues (tempModelOrder+1:end);
    else
      break;
    end
    
    % append the number of eigenvalues in the signal cluster to the temp array 
    tempOrderArray(end+1) = tempModelOrder;
    
end

% find the 
cumulativeOrderSum = cumsum(tempOrderArray);
eigenValues = sort (unorderedEigenValues, 'descend');

% ------------------------ validation condition ----------------

maxModelOrder = 0;

% for each iteration it will calculate the quadratic sum and the linear sum
% and compare it with the smallest value greater than the ratio
for i=1:10
    quadraticSum = 0;
    linearSum = 0;
       
    % find both the linear and quadratic sum of eigenvalues
    for j=(cumulativeOrderSum(i)+1):length(eigenValues)
        linearSum = linearSum + (eigenValues(j)); 
        quadraticSum = quadraticSum + (((eigenValues(j))).^2);  
    end
    
    % find the eigenvalue > 5 times (quadratic sum / linear sum) that is
    % the value we need for finding the model order.
    if (eigenValues(cumulativeOrderSum(i))) >= 5*(quadraticSum/linearSum)
          maxModelOrder = cumulativeOrderSum(i);
    end
end

% ----------------------------------------------------------------- %

% calculating the model order - for even M/2, for odd (M + 1)/2
        
if rem(max(maxModelOrder),2) == 0
    modelOrder = max(maxModelOrder)/2;
else
    modelOrder = (max(maxModelOrder)+1)/2;
end

% print the model order
X = ['Model order is ',num2str(modelOrder)];
disp(X)