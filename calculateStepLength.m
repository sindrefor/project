function [stepLength,newErr] = calculateStepLength(dims,gradient,oldErr,bg, source,trueRec,f,it)
%% Depending on how you solve the problem you most probably need more input variables                      
    recording = zeros(dims.nt,length(dims.recPos),'single');
    stepLength = 256;
%     cheatStep = 2*[256 128 128 64 128 32 32 32 16 16 128 64 32 128 32 8 32 16 8 8 64 32 16 32 16 16 16 16 16 16];
%     if it<length(cheatStep), stepLength = cheatStep(it); end
    newErr = inf;
    while (newErr > oldErr)
        newErr = 0;
        stepLength = stepLength/2;
        bg0=bg+stepLength*gradient;
        S = zeros(dims.ny,dims.nx);
        S(:) = dims.dt^2.*bg0(:).^2./dims.dx^2;
        for s = 1:dims.ds:length(dims.srcPos)
                        U0 = zeros(dims.ny,dims.nx); % for time t-1
            U1 = U0; % for time t
            u = U0; % for time t+1
            for t = 1:dims.nt
                %  Solve wave equation using test model update
                U0(:,:)=U1(:,:);
U1(:,:)=u(:,:);
U1(dims.srcPos(s))=U1(dims.srcPos(s))+source(t); % adds the driving force

for k=5:dims.ny-4 % updating inner nodes (?-point stencil)
    for j=5:dims.nx-4
        u(k,j)=U1(k,j)*(2+2*(-205/72)*S(k,j))-U0(k,j)+((U1(k,j+1)+U1(k,j-1))+(U1(k+1,j)+U1(k-1,j)))*S(k,j)*(8/5)+(U1(k,j-2)+U1(k,j+2)+U1(k-2,j)+U1(k+2,j))*S(k,j)*(-1/5)+(U1(k,j-3)+U1(k,j+3)+U1(k-3,j)+U1(k+3,j))*S(k,j)*(8/315)+(U1(k,j-4)+U1(k,j+4)+U1(k-4,j)+U1(k+4,j))*S(k,j)*(-1/560);
    end
end
                %  Record traces
                recording(t,:) = u(dims.recPos);
            end
            %% Calculate new error and check against old

        recordingabsmax = recording(:,:); recordingabsmax = max(abs(recordingabsmax(:)));
        tRecabsmax = trueRec(:,:,s); tRecabsmax = max(abs(tRecabsmax(:)));
        chi = recording./recordingabsmax - trueRec(:,:,s)./tRecabsmax;
        chi(:,s)=0;
        if s>3 && s<(dims.mx-3) && f==4, chi(150:245,:)=0;
        elseif f==6, chi(50:250,:)=0; 
        elseif f==8, chi(78:120,:)=0; 
        end
        newErr=newErr+sum(abs(chi(:)));
        end
        newErr=newErr/length(1:dims.ds:length(dims.srcPos));
    end
end