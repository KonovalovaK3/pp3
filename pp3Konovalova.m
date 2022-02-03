clc; clear all; close all;
delta = 0.2; 
N = 256;
x = 8.0*pi*(1:N)'/N;
u = 15 - 15*tanh(0.5*(x-4*pi)) -15*tanh(0.5*(x-4*pi)).^2 + 15*tanh(0.5*(x-4*pi)).^3;
v = fft(u);
dt = 1.0 / 64.0;
k = [ 0:N/2-1, 0.0, -N/2+1:-1 ]' / 2.0;
L = k.^2 - k.^4 - 1i*delta*k.^3;
E = exp ( dt * L );
E2 = exp ( dt * L / 2.0 );
M = 16;
r = exp (1i * pi * ( (1:M) - 0.5 ) / M );
LR = dt * L(:,ones(M,1)) + r(ones(N,1),:);

  Q  = dt * real ( mean ( ...
    ( exp ( LR / 2.0 ) - 1.0 ) ./ LR, 2 ) );
  f1 = dt * real ( mean ( ...
     ( - 4.0 - LR + exp ( LR ) .* ( 4.0 - 3.0 * LR + LR.^2 )) ./ LR.^3, 2 ) );
  f2 = dt * real ( mean ( ...
    ( 2.0 + LR + exp ( LR ) .* ( - 2.0 + LR ) ) ./ LR.^3, 2 ) );
  f3 = dt * real ( mean ( ...
    ( - 4.0 - 3.0 * LR - LR.^2 + exp ( LR ) .* ( 4.0 - LR ) ) ./ LR.^3, 2 ) );

uu = u;
  tt = 0.0;
  tmax = 50.0;
  nmax = round ( tmax / dt );
  nplt = floor ( ( tmax / 100 ) / dt );
  g = - 0.5 * 1i * k;
  
  for j = 1 : nmax

    t = j * dt;

    Nv = g .* fft ( real ( ifft ( v ) ) .^2 );
    a = E2 .* v + Q .* Nv;
    Na = g .* fft ( real ( ifft ( a ) ) .^2 );
    b = E2 .* v + Q .* Na;
    Nb = g .* fft ( real ( ifft ( b ) ) .^2 );
    c = E2 .* a + Q .* ( 2.0 * Nb - Nv );
    Nc = g .* fft ( real ( ifft ( c ) ) .^2 );

    v = E .* v + Nv .* f1 + 2.0 * ( Na + Nb ) .* f2 + Nc .* f3;

    if ( mod ( j, nplt ) == 0 )
      u = real ( ifft ( v ) );
      uu = [ uu, u ];
      tt = [ tt, t ];
    end
    
  end

  uu = uu';
waterfall(x, tt, uu);
