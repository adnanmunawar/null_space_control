kp = [2,2]';
kv = sqrt(kp);
vmax = 2.0;
dx = [1.1, 1.1]';
xyz = [5.0, 5.0]';
target_xyz = [0.0, 0.0]';
lamb = kp ./ kv
x_tilde = xyz - target_xyz
sat = vmax ./ (lamb .* abs(x_tilde))
scale = ones(2,1)
if any(sat < 1)
    [val, index] = min(sat)
    unclipped = kp .* x_tilde(index)
    clipped = kv .* vmax .* sign(x_tilde(index))
    scale = ones(2,1) .* clipped ./ unclipped
    scale(index) = 1
end
clipped = sat ./ scale
clipped(clipped > 1) = 1
clipped(clipped < 0) = 0
u_xyz = -kv .* (dx + clipped .* scale .* lamb .* x_tilde)