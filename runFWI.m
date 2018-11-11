clear,clc
%% Setting up dimensions
dims.dy =     10; % [m]
dims.dx =     10; % [m]
dims.dt = 1.0e-3; % [s]

dims.ny = 201; % Cells in y-direction
dims.nx = 301; % Cells in x-direction
dims.nt = 801; % Amount of time steps

%% Model dimensions
dims.modely = 100:150;
dims.modelx = 100:200;
dims.my = length(dims.modely);
dims.mx = length(dims.modelx);

%% Source locations
sx = min(dims.modelx):max(dims.modelx);
sy = min(dims.modely)*ones(1,length(sx));
dims.srcPos = sy + dims.ny*sx;

%% Receiver locations
rx = min(dims.modelx):max(dims.modelx);
ry = min(dims.modely)*ones(1,length(rx));
dims.recPos = ry+dims.ny*rx;

%% Creating background model
bg = zeros(dims.ny,dims.nx,'single');
bg(:) = 2.0e3;         % [m/s] - Background
bg(115:end,:) = 2.3e3; % [m/s] - Layer

%% Begin iteration
model = bg;     % Starting model
dims.ds = 10;   % Grid point distance between sources
maxIter = 10;   % Maximum number of iterations per frequency
freqs = [4];  % Frequencies to use in inversion

errVec = zeros(1,maxIter*length(freqs));
stepVec = errVec;

it = 1; tic;
for f = freqs
    %% Generating ricker source signature wavelet 
    source = rickerWave(f,dims);
    %% Load true recording
    load (['trueRec_',num2str(f),'Hz.mat']);   
    for i = 1:maxIter
        %% Calculate gradient ## IMPLEMENT ##
        [gradient,err, chi] = calculateGradient(dims, trueRec, source, model,f,it);
            
        %% Taper gradient
        gradient = taperGradient(dims,gradient);

        %% Calculate step length ## IMPLEMENT ##
        [stepLength,err] = calculateStepLength(dims,gradient,err,model,source,trueRec,f,it);

        %% Update model
        model = model + stepLength*gradient;
        stepVec(it) = stepLength;
        errVec(it) = err; 
        it = it + 1;
        figure(1)
        imagesc(model(dims.modely,dims.modelx))
        figure(2)
        plot(errVec)
        figure(3)
        plot(stepVec)
        toc
    end
end