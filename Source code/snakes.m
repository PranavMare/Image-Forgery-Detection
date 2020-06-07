function v = snakes (I, v)
  % INPUT: N by M image I, a contour v of n control points
  % OUTPUT: converged contour v of n control points

  E_image = generateImageEnergy (I);

  while not converged
    F_cont = weight.alpha * contourDerivative(v, 2);
    F_curv = weight.beta * contourDerivative(v, 4);
    F_image = interp2 (E_image, v(:,2), v(:,1));
    F_image_norm = weight.k * F_image ./ norm (F_image);
    F_con = inputForces();

    F_internal = F_cont + weight.external * F_curv;
    F_external = weight.external * (F_image + F_con);

    v = updateSnake(v, F_internal, F_external);

    checkConvergence ();
  end

end