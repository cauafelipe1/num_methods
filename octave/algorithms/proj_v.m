function proj = proj_v(u, v)
% proj_v
% calculates the vectorial projection of vector u onto vector v.
%   p = proj(u, v) returns the vector representing the orthogonal projection 
%   of u onto v.

    % check for zero vector to prevent division by zero
    if (norm(v) < 1e-15)
        error("division by zero: cannot project onto the null vector");
    endif

    proj = (dot(u, v) / dot(v, v)) * v;
end