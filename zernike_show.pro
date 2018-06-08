function zernike_show

  DIMMAT = 256;
  FIRSTZERNIKE = 1;
  NZERNIKE = 22;


  aa = POLZER(DIMMAT,DIMMAT-1,FIRSTZERNIKE,NZERNIKE)

  coeff =[0.000041,  0,  0,  0.00004839,  0,  0,  0,  0,  0,  0,  0.00000249,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  -0.00001405];
  coeff =[0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0.00000249,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  -0.00001405]; %piston,TT, defocus removed
  wf = make_array(DIMMAT,DIMMAT)
  for idx=0,NZERNIKE-1 do wf = wf+ coeff[idx]*aa[*,*,idx]

  tvsg, wf, 3,0
  rw, 800, 800, 1
  shade_surf, wf
    
  

stop&stop

end