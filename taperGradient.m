function G = taperGradient(dims, dm1)
    %% Taper gradient near sources to avoid interference
    G = zeros(dims.ny, dims.nx);
    G(dims.modely,dims.modelx) = dm1;
    
    [~,n] = size(G);
    taper = sin(linspace(0,pi/2,16));
    taper = exp(linspace(0,1,16))-1;
    maxtaper = max(taper);
    taper = taper./maxtaper;
    taper = repmat(taper,n,1)';   
    G(100:115,:) = taper.*G(100:115,:);
    G(100:102) = 0;
    
    %% Normalise gradient
    scale = 1.0/max(abs(G(:)));
    G = scale*G; 
end