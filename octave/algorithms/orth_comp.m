function w = orth_comp(u, v)
% orthogonal complement
% calculates the orthogonal complement of vector u relative to v.
%   w = orth_comp(u, v) returns the vector component of u that is 
%   perpendicular to v.

    % formula: w = u - projv(u, v)
    w = u - proj_v(u, v);
endfunction