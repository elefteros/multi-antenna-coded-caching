function[time] = timeCalculation(t, streams, hp, rp, power, iter)

    S= (t+ streams) * log( 1+ power.*(hp'*rp)/streams)*iter;
    
    time=S.^(-1);
    
end
