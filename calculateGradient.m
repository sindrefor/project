function [gradient, err, chi] = calculateGradient(dims, trueRec, source, bg,f,it) 
%% Depending on how you solve the problem you most probably need more input variables
    recording = zeros(dims.nt,length(dims.recPos),'single');
    gradient  = zeros(dims.my,dims.mx,'single');
    forwardField = zeros(dims.my,dims.mx,dims.nt,'single'); 
    adjointField = zeros(dims.my,dims.mx,dims.nt,'single');
    err = 0;
    S = zeros(dims.ny,dims.nx);
    S(:) = dims.dt^2.*bg(:).^2./dims.dx^2; 
    for s = 1:dims.ds:length(dims.srcPos)
        %% Run forward simulation on background model
            U0 = zeros(dims.ny,dims.nx); % for time t-1
            U1 = U0; % for time t
            u = U0; % for time t+1
        for t = 1:dims.nt
            
U0(:,:)=U1(:,:);
U1(:,:)=u(:,:);
U1(dims.srcPos(s))=U1(dims.srcPos(s))+source(t); % adds the driving force

for k=5:dims.ny-4 % updating inner nodes (?-point stencil)
    for j=5:dims.nx-4
        u(k,j)=U1(k,j)*(2+2*(-205/72)*S(k,j))-U0(k,j)+((U1(k,j+1)+U1(k,j-1))+(U1(k+1,j)+U1(k-1,j)))*S(k,j)*(8/5)+(U1(k,j-2)+U1(k,j+2)+U1(k-2,j)+U1(k+2,j))*S(k,j)*(-1/5)+(U1(k,j-3)+U1(k,j+3)+U1(k-3,j)+U1(k+3,j))*S(k,j)*(8/315)+(U1(k,j-4)+U1(k,j+4)+U1(k-4,j)+U1(k+4,j))*S(k,j)*(-1/560);
    end
end
            % Record traces
            recording(t,:) = u(dims.recPos);
            % Save forward field for use in correlation
            forwardField(:,:,t) = u(dims.modely,dims.modelx);
        end
        %% Calculate difference and error
        recordingabsmax = recording(:,:); recordingabsmax = max(abs(recordingabsmax(:)));
        tRecabsmax = trueRec(:,:,s); tRecabsmax = max(abs(tRecabsmax(:)));
        chi = recording./recordingabsmax - trueRec(:,:,s)./tRecabsmax;
        chi(:,s)=0;
        if s>3 && s<(dims.mx-3) && f==4, chi(150:245,:)=0;
        elseif f==6, chi(50:250,:)=0; 
        elseif f==8, chi(78:120,:)=0; 
        end
        err=err+sum(abs(chi(:)));
        %% Run adjoint simulation
            U0 = zeros(dims.ny,dims.nx); % for time t-1
            U1 = U0; % for time t
            u = U0; % for time t+1
        for t = 1:dims.nt
            % Solve wave equation using the difference (chi) as sources
            U0(:,:)=U1(:,:);
            U1(:,:)=u(:,:);
            U1(dims.srcPos)=U1(dims.srcPos)+chi(end+1-t,:); % adds the driving force
            for k=5:dims.ny-4 % updating inner nodes (?-point stencil)
                for j=5:dims.nx-4
                    u(k,j)=U1(k,j)*(2+2*(-205/72)*S(k,j))-U0(k,j)+((U1(k,j+1)+U1(k,j-1))+(U1(k+1,j)+U1(k-1,j)))*S(k,j)*(8/5)+(U1(k,j-2)+U1(k,j+2)+U1(k-2,j)+U1(k+2,j))*S(k,j)*(-1/5)+(U1(k,j-3)+U1(k,j+3)+U1(k-3,j)+U1(k+3,j))*S(k,j)*(8/315)+(U1(k,j-4)+U1(k,j+4)+U1(k-4,j)+U1(k+4,j))*S(k,j)*(-1/560);
                end
            end
            % Save adjoint field for use in correlation
            adjointField(:,:,dims.nt-t+1) = u(dims.modely,dims.modelx);
        end
        %% Correlate
        for t = 2:dims.nt-1
            % Calculate the time derivative of the displacement to
            % gradient.
            for i=1:dims.my
                for j=1:dims.mx
                    gradient(i,j) = gradient(i,j)+((forwardField(i,j,t+1)-forwardField(i,j,t-1))/(2*dims.dt)) * ((adjointField(i,j,t+1)-adjointField(i,j,t-1))/(2*dims.dt));
                end
            end
        end
    end
    err=err/length(1:dims.ds:length(dims.srcPos));
end