    clc; clear variables

    %%% dmt with the method from "Solving the Max-Min SNR Optimization
    %%% Problem of Array Antenna Signal Processing"

    K=30; % number of users
    iter=5e1; % number of realizations
    Nsolver=10; %% for the optimization problem
    Nsolver2=10;  %% pick various initial points to explore
    delta=1e-3; %% for the optimization problem
    Kvec=1:K;  
    Pdb=-100:1:200;              power=db2pow(Pdb);
    
    LL=1:4;     % select all the antenna values of interest
    tt=1:4;     % select all the cache redundancy values of interest
    
    
for L=LL
    for t=tt
    timeOpt=zeros(L, length(Pdb));            % for optimal

    rateOpt=zeros(L, length(Pdb));            % for optimal


    bOptMulticast=zeros(L,1);
    H=randn(L,K)+1j*randn(L,K); %%% Channel matrix Complex Gaussian
    for streams=1:min(L, K-t-1)
        for s=1:iter %%% To average various user subsets
             bOpt=zeros(L,streams);  
             rOpt=zeros(streams,1);


             sigma  =  sort(Kvec(randperm(K, t+1))); % select users associated to the XOR
             LAMBDA =  setdiff(Kvec, sigma);  % set of remaining users
             lambda =  LAMBDA(randperm( K-t-1, streams-1) ); % select precoded users
             phi    =  sigma(randperm( t+1, 1));    %% select precoded user from the
                                                    %%% XOR % useful only when streams>1
                %%% print useful information!!
                clc;   
                t
                L           %%% number of antennas
                streams             %%% number of streams
                s                   %%% sample number
              %%% print useful information!!
     
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

             timeOpt(streams,:)=timeOpt(streams,:)+ timeCalculation( t, streams, hOpt, rOpt, power, iter);
        end
        % per-user rate => easier to compare with other settings
        rateOpt(streams,:)=1./timeOpt(streams,:)/K;
    end
    
str=['K', num2str(K), 'L', num2str(L), 't', num2str(t), 'iter', num2str(iter)];
   
save(str)
    end        

end

























