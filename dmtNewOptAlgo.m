    clc; clear variables

    %%% dmt with the method from "Solving the Max-Min SNR Optimization
    %%% Problem of Array Antenna Signal Processing"

    K=30; % number of users
    iter=5e1;
    Nsolver=10;
    Nsolver2=10;  %% pick various initial points to explore
    delta=1e-3;
    Kvec=1:K;
    Pdb=-100:1:200;              power=db2pow(Pdb);
    
    
for L=2
    for t=8
%     timeBA=zeros(L, length(Pdb));             % to use for best one
    timeOpt=zeros(L, length(Pdb));            % for optimal
%     timeGS=timeBA;                            % for Gram-Schmidt
%     timeAAver=timeBA;                         % for antenna average
%     timeRV=timeBA;                            % random

%     rateBA=zeros(L, length(Pdb));             % to use for best one
    rateOpt=zeros(L, length(Pdb));            % for optimal
%     rateGS=rateBA;                           % for Gram-Schmidt
%     rateAAver=rateBA;                        % for antenna average
%     rateRV=rateBA;                           % random
    bOptMulticast=zeros(L,1);
    H=randn(L,K)+1j*randn(L,K); %%% Channel matrix Complex Gaussian
    for streams=1:min(L, K-t-1)
        for s=1:iter %%% To average various user subsets
             bOpt=zeros(L,streams);  
             rOpt=zeros(streams,1);
%              rBA=rOpt;
%              rAAver=rOpt;
%              rRV=rOpt;
%              rGS=rOpt;
             sigma  =  sort(Kvec(randperm(K, t+1))); % select users associated to the XOR
             LAMBDA =  setdiff(Kvec, sigma);  % set of remaining users
             lambda =  LAMBDA(randperm( K-t-1, streams-1) ); % select precoded users
             phi    =  sigma(randperm( t+1, 1));    %% select precoded user from the
                                                    %%% XOR % useful only when streams>1
%              if mod(s,5)==0 %%% print useful information!!
                clc;   
                t
                L           %%% number of antennas
                streams             %%% number of streams
                s                   %%% sample number
%              end             
             %%% calculate the OPTIMAL beamforming/precoding vectors
             bOpt(:,1)=precoderVector( H, sigma,  lambda  , Nsolver2, Nsolver, delta);
             for i=1:streams-1
                 lambdaTemp=setdiff(lambda, lambda(i));
                 others=[phi, lambdaTemp];
                 bOpt(:,i+1)=precoderVector( H, lambda(i),  others  , Nsolver2, Nsolver, delta);
             end
             rOpt(1)=min( abs( H(:,sigma)'*bOpt(:,1))  )^2;
             for i=1:streams-1
                rOpt(i+1) = abs( H(:,lambda(i))'*bOpt(:,i+1))^2;
             end
             hOpt=powerAllocation( rOpt);
%              
%              %%% calculate the BEST ANTENNA beamforming/precoding vectors
%              rBA(1)=bestAntenna( H, sigma,  lambda );
%              for i=1:streams-1
%                 lambdaTemp=setdiff(lambda, lambda(i));
%                 others=[phi, lambdaTemp] ;
%                 rBA(i+1) = bestAntenna( H, lambda(i), others);
%              end
%              hBA=powerAllocation(rBA);
%              %%% calculate the ANTENNA AVERAGE beamforming/precoding vectors
%              rAAver(1)= sumAntennas( H, sigma,  lambda );
%              for i=1:streams-1
%                 lambdaTemp=setdiff(lambda, lambda(i));
%                 others=[phi, lambdaTemp] ;
%                 rAAver(i+1) = sumAntennas( H, lambda(i), others);
%              end
%              hAAver=powerAllocation(rAAver);
%              
%              %%% calculate RANDOM beamforming/precoding vectors
%              rRV(1)= randomAllocation( H, sigma,  lambda );
%              for i=1:streams-1
%                 lambdaTemp=setdiff(lambda, lambda(i));
%                 others=[phi, lambdaTemp] ;
%                 rRV(i+1) = randomAllocation( H, lambda(i), others);
%              end
%              hRV=powerAllocation(rRV);
%              %%% calculate Gram-Schmidt beamforming/precoding vectors
%              bOptNEW=precoderVector(H, sigma, [],  Nsolver2, Nsolver, delta );
%              bGS = gramSchmidt( H, lambda, bOptNEW );
%              rGS(1)=min( abs( H(:,sigma)'*bGS(:,1))  )^2;
%              for i=1:streams-1
%                 rGS(i+1) = abs(  H(:,lambda(i))'*bGS(:,i) )^2;
%              end
%              hGS=powerAllocation(rGS);
             %%% i am only interested in the ...
             %%% rate achieved by the XOR due to the addition of new (precoded) users.
%            timeOpt(streams,:)=timeOpt(streams,:)+ timeCalculation( t, 1, 1, rOpt(1), power, iter   );
             timeOpt(streams,:)=timeOpt(streams,:)+ timeCalculation( t, streams, hOpt, rOpt, power, iter);
%              timeAAver(streams,:)=timeAAver(streams,:)+timeCalculation( t, streams, hAAver, rAAver, power, iter   );
%              timeBA(streams,:)=timeBA(streams,:)+timeCalculation( t, streams, hBA, rBA, power, iter   );
%              timeRV(streams,:)=timeRV(streams,:)+timeCalculation( t, streams, hRV, rRV, power, iter   );
%              timeGS(streams,:)=timeGS(streams,:)+timeCalculation( t, streams, hGS, rGS, power, iter   );
        end
        % per-user rate => easier to compare with other settings
        rateOpt(streams,:)=1./timeOpt(streams,:)/K;
%         rateBA(streams,:)= 1./timeBA(streams,:)/K;
%         rateAAver(streams,:)= 1./timeAAver(streams,:)/K;
%         rateRV(streams,:)= 1./timeRV(streams,:)/K;
%         rateGS(streams,:)=1./timeGS(streams,:)/K;
    end
% %         hold all
% %         plot(Pdb, rateOpt)
% %         plot(Pdb, rateGS)
% %         plot(Pdb, rateRV)
% %         plot(Pdb, rateAAver)
% %         plot(Pdb, rateBA)

str=['K', num2str(K), 'L', num2str(L), 't', num2str(t), 'iter', num2str(iter)];
% str=['K', num2str(K), 'L', num2str(L), 't', num2str(t), 'iter', num2str(iter), 'onlyOptimal'];
   
save(str)
    end        

end

























