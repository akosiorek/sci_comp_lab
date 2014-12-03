function [] = PlotResults(T, Nx, Ny, index, Ns, Dts)

  x = 0: 1/(Nx+1): 1;   % Set x and y
  y = 0: 1/(Ny+1): 1;
  subplot(Ns, Dts, index); %Assign subplot position
  set(gcf, 'renderer', 'zbuffer');
  mesh(x, y, T); %Surface plot  
  
  pbaspect([1 1 1])
  xlabel('x');
  ylabel('y');
  zlabel('T');
  
end
