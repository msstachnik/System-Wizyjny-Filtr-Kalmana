function [ Y ] = down_up_limit(X, min, max )
% [ Y ] = up_down_limit(X, min, max )
% function limit x to uper or don bound

if X < min
    Y = min;
elseif X > max
    Y = max;
else
    Y = X;
end
        


end

