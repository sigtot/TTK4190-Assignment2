function [saturated_value, is_saturated] = saturate(value, limit)
%saturate Saturates a given value by a given limit
%   Given a value and limit, returns a value which is in the 
%   range -limit to limit. Also returns a boolean value on
%   wether or not the value is saturated. 
if value > limit
    saturated_value = limit; 
    is_saturated = true; 
elseif value < -limit
    saturated_value = -limit; 
    is_saturated = true; 
else
    saturated_value = value; 
    is_saturated = false; 
end

end

