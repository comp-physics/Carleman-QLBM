function cal_feq(rho0, vx, vy, wm, vcx, vcy, a, b, c, d)
  return rho0 .* wm .* (a .+ b .* (vx*vcx .+ vy*vcy) .+ c.*(vx*vcx .+ vy*vcy).^2 .+ d.*(vx.^2 .+ vy.^2))
end
