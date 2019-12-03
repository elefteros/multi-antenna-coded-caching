function[hp] = powerAllocation( rates )
    
streams=length(rates);
if streams>1
    cvx_begin
        cvx_quiet('true');
        variable X(streams,1)
        maximize( rates'*X )
        subject to
            X'*X<=1
    cvx_end
    hp=X;                   
    hp=hp/norm(hp);
else
    hp=1;
end

end

