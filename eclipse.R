Elipse centered at (x0, y0) with axes _a_ and _b_:
  
  theta <- seq(0, 2 * pi, length=(big enough))
  x <- x0 + a * cos(theta)
  y <- y0 + b * sin(theta)
  plot(x, y, type = "l")
  
  :-)

Or, if you want to rotate it alpha degrees from the x-axis:
  
  x <- x0 + a * cos(theta) * cos(alpha) - b * sin(theta) * sin(alpha)
y <- y0 + a * cos(theta) * sin(alpha) + b * sin(theta) * cos(alpha)
