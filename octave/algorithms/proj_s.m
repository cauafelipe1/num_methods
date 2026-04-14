function proj = proj_s(u, v)
% proj_s
% calculates the scalar projection of vector u onto vector v.
%   proj = proj_s(u, v) returns the length of the orthogonal projection 
%   of u onto v.

    % check for zero vector to prevent division by zero
    if (norm(v) < 1e-15)
        error("division by zero: cannot project onto a zero vector");
    endif

    % scalar projection formula: (u . v) / ||v||
    proj = dot(u, v) / norm(v);
end