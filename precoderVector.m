function[bBest] = precoderVector( H, sigma,  lambda  , Nsolver2, Nsolver, delta)
%%% calculate the beamformer using the method from "Solving the Max-Min SNR
%%% Optimization Problem of Array Antenna Signal Processing"
streams=length(lambda)+1;
[L,~]=size(H);
t=length(sigma)-1;
bBest=zeros(L,1); %%% to pick best beamformer between the Nsolver2 explorations
Nu=null(H(:,lambda)'); %%% calculate null to update beta iff no solution can be found
for n2=1:Nsolver2
    theta=2*pi*rand(1,t+1); %%% init theta
    for n=1:Nsolver   %%% iterations to refine the solution
       %%% calculate beamformer/precoder for the XOR
        cvx_begin
        cvx_quiet('true');
        variable X(L,1) complex
        minimize( norm(X) )
        subject to
        for k=1:t+1
            real( exp(-1j*theta(k)) * H(:,sigma(k))'*X)>=1   
        end
        %%% only when streams>1
        for mm=1:streams-1
            H(:,lambda(mm))'*X == 0
        end
        cvx_end
        if isnan(X)==0      %%% if a value has been reached
            bNew=X/norm(X); %% beamformer/precoder for the XOR
%           disp('NOT NaN')
        else
%       disp('Is NaN')      %%% randomize new beamformer and theta
            if streams==1
                bNew=rand(L,1)+1j*rand(L,1); %%% if a value cannot be reached
            else
                bNew=Nu*rand(L-streams+1,1);
            end
            bNew=bNew./norm(bNew);
            theta=2*pi*rand(1,t+1); %#ok<NASGU> %used in cvx
        end
        
        %%% compare between previous and new
        if min( ((H(:,sigma)'*bNew)').'.*(H(:,sigma)'*bNew) )...
                        >= min( ((H(:,sigma)'*bBest)').'.*(H(:,sigma)'*bBest))
        %disp('Updating beta')
            theta=angle( H(:,sigma)'*bNew );
            bBest=bNew;       
            
        else  %%% found a worse beta, thus updating theta
              %%% disp('Worse beta, thus updating theta')
             theta=2*pi*rand(1,t+1);
        end
                %%% all users have very close SNR thus almost optimal
        if max( ((H(:,sigma)'*bBest)').'.*(H(:,sigma)'*bBest))...
                        -min( ((H(:,sigma)'*bBest)').'.*(H(:,sigma)'*bBest))<delta
            %disp('Small error')
            break;
        end
    end    
end
end


