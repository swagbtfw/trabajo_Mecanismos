format long;

%% RESOLUCIÓN CALCULOS PROYECTO MECANISMOS
% Datos: 
enunciado = imread("enunciado_proyecto.png");
imshow(enunciado)


%Calculamos las posiciones de los puntos, sabiendo que O2 = (0.4684, 0.0233) y que
%02=250º.
%Calulamos SEN(02) y COS(02) sabiendo que mide 0.06m
O2 = [0.4684, 0.0233]
ax =O2(1,1) + 0.06*cos(250*pi/180);
ay =O2(1,2) + 0.06*sin(250*pi/180);
pA = [ax,ay];

%Ahora podemos calcular los demás puntos creando circunferencias y buscando
%por la diferencia de distancias donde van a coincidir

%Calculamos B
syms Bx By

% Coordenadas conocidas
O4x = 0.73;  O4y = 0.023;
% Ecuaciones de circunferencias
eq1 = (Bx - O4x)^2 + (By - O4y)^2 == 0.125^2;
eq2 = (Bx - ax)^2  + (By - ay)^2  == 0.24^2;

% Resolver
sol = solve([eq1, eq2], [Bx, By], 'Real', true); %Condicion para que tenga soluciones reales.

eval(sol.Bx)
eval(sol.By)

%Se obtienen 2 soluciones para x y 2 para y, simplemente comprobamos los
%valores en otras aplicaciones como geogebra para encontrar cual es el
%valor mas optimo. En este caso : B= (0.640425728166680, 0.110186293793982);
pB = [ 0.640425728166680,  0.110186293793982];
% Calculamos C
syms Cx Cy

Bx  = pB(1); By = pB(2);

eqC1 = (Cx - O4x)^2 + (Cy - O4y)^2 == 0.25^2;
eqC2 = (Cx - Bx)^2  + (Cy - By)^2  == 0.125^2;

solC = solve([eqC1, eqC2], [Cx, Cy], 'Real', true); %Condición para que tenga soluciones reales.

eval(solC.Cx);
eval(solC.Cy);

% Selección de la rama correcta según geometría
pC = [double(solC.Cx(2)), double(solC.Cy(2))];  
%D
% Calculamos D
syms Dx Dy

O6x = 0.315; O6y = 0.185;
Cx  = pC(1); Cy = pC(2);

eqD1 = (Dx - O6x)^2 + (Dy - O6y)^2 == 0.43^2;
eqD2 = (Dx - Cx)^2  + (Dy - Cy)^2  == 0.60^2;

solD = solve([eqD1, eqD2], [Dx, Dy], 'Real', true); %Condición para que tenga soluciones reales.

eval(solD.Dx);
eval(solD.Dy);

pD = [double(solD.Dx(1)), double(solD.Dy(1))];  

DC = 0.6;
DE = 0.8879;
EC = sqrt(DE^2 + DC^2 - 2 * DE * DC * cos(6.63*pi/180))
%Con la constante de la ley de cosenos, podemos calcular lo demás : 
k= EC/sin(6.63*pi/180)

sinC = DE/k; asin(sinC)*180/pi 
sinD =DC/k; asin(sinD)*180/pi 
% Direccion absoluta de DC desde D
vDC = pC - pD;
thetaDC = atan2(vDC(2), vDC(1));

alpha = 6.63 * pi/180;
thetaE = thetaDC + alpha;

% Coordenadas de E+
pE = pD + DE * [cos(thetaE), sin(thetaE)];

pA
pB
pC
pD
pE
%% Graficar mecanismo
figure; hold on; axis equal; grid on;

% Dibujar puntos
plot(pA(1), pA(2), 'ro', 'MarkerSize', 8, 'DisplayName','A');
plot(pB(1), pB(2), 'go', 'MarkerSize', 8, 'DisplayName','B');
plot(pC(1), pC(2), 'bo', 'MarkerSize', 8, 'DisplayName','C');
plot(pD(1), pD(2), 'mo', 'MarkerSize', 8, 'DisplayName','D');
plot(pE(1), pE(2), 'ko', 'MarkerSize', 8, 'DisplayName','E');

% Dibujar pivotes fijos
plot(O2(1), O2(2), 'ks', 'MarkerSize', 8, 'MarkerFaceColor','y','DisplayName','O2');
plot(O4x, O4y, 'ks', 'MarkerSize', 8, 'MarkerFaceColor','y','DisplayName','O4');
plot(O6x, O6y, 'ks', 'MarkerSize', 8, 'MarkerFaceColor','y','DisplayName','O6');

% Dibujar barras
plot([O2(1) pA(1)], [O2(2) pA(2)], 'r-', 'LineWidth',2);
plot([pA(1) pB(1)], [pA(2) pB(2)], 'g-', 'LineWidth',2);
plot([pB(1) pC(1)], [pB(2) pC(2)], 'b-', 'LineWidth',2);
plot([O4x pB(1)], [O4y pB(2)], 'm-', 'LineWidth',2);
plot([O4x pC(1)], [O4y pC(2)], 'm-', 'LineWidth',2);
plot([pC(1) pD(1)], [pC(2) pD(2)], 'c-', 'LineWidth',2);
plot([O6x pD(1)], [O6y pD(2)], 'k-', 'LineWidth',2);
plot([pD(1) pE(1)], [pD(2) pE(2)], 'y-', 'LineWidth',2);
plot([pC(1) pE(1)], [pC(2) pE(2)], 'Color',[0.5 0 0.5],'LineWidth',2,'DisplayName','C-E'); % nueva barra

title('Mecanismo de seis barras - Configuración θ2=250°');
legend show;

%% CÁLCULO DE VELOCIDADES (CINEMÁTICA)
fprintf('\n--- RESULTADOS DE VELOCIDADES ---\n');

% 1. Datos de Entrada
w2 = 2; % rad/s (k)

% 2. Velocidad de A (vA = w2 x rA)
rA = pA - O2; % Vector O2->A
vA = [-w2 * rA(2), w2 * rA(1)]; % Producto cruz plano (-w*y, w*x)
fprintf('vA = [%.4f, %.4f] m/s\n', vA(1), vA(2));

% Resolver w3 y w4 para hallar vB)
rAB = pB - pA;
rO4B = pB - [O4x, O4y];

Matriz1 = [-rAB(2),  rO4B(2); 
            rAB(1), -rO4B(1)];
TerminoIndep1 = [-vA(1); -vA(2)];

sol_w34 = Matriz1 \ TerminoIndep1; % Resolver sistema Ax=B
w3 = sol_w34(1);
w4 = sol_w34(2);

fprintf('w3 (Biela) = %.4f rad/s\n', w3);
fprintf('w4  = %.4f rad/s\n', w4);

vB = [-w4 * rO4B(2), w4 * rO4B(1)];
fprintf('vB = [%.4f, %.4f] m/s\n', vB(1), vB(2));


% 4. Velocidad de C
rO4C = pC - [O4x, O4y];
vC = [-w4 * rO4C(2), w4 * rO4C(1)]; 
fprintf('vC = [%.4f, %.4f] m/s \n', vC(1), vC(2));


% (Resolver w5 y w6 para hallar vD)
rCD = pD - pC;
rO6D = pD - [O6x, O6y];
Matriz2 = [-rCD(2),  rO6D(2); 
            rCD(1), -rO6D(1)];
TerminoIndep2 = [-vC(1); -vC(2)];

sol_w56 = Matriz2 \ TerminoIndep2;
w5 = sol_w56(1);
w6 = sol_w56(2);

fprintf('w5 (Acoplador C-D) = %.4f rad/s\n', w5);
fprintf('w6 (Salida D-E) = %.4f rad/s\n', w6);

% Calcular vD
vD = [-w6 * rO6D(2), w6 * rO6D(1)];
fprintf('vD = [%.4f, %.4f] m/s\n', vD(1), vD(2));


% 6. PUNTO FINAL: Velocidad de E
rO6E = pE - [O6x, O6y];
vE = [-w6 * rO6E(2), w6 * rO6E(1)];

fprintf('vE = [%.4f, %.4f] m/s (Resultado Final)\n', vE(1), vE(2));
mag_vE = norm(vE);
fprintf('Magnitud |vE| = %.4f m/s\n', mag_vE);

%% CÁLCULO DE ACELERACIONES
fprintf('\n--- RESULTADOS DE ACELERACIONES ---\n');
alpha2 = 0; 

% Aceleración de A
aA_n = -w2^2 * rA;
aA_t = [-alpha2 * rA(2), alpha2 * rA(1)];
aA = aA_n + aA_t;

fprintf('aA (Normal)   : [%.4f, %.4f]\n', aA_n);
fprintf('aA (Tangencial): [%.4f, %.4f]\n', aA_t);
fprintf('aA (TOTAL)    : [%.4f, %.4f] m/s^2\n\n', aA);

% Aceleración de B
a_n_BA  = -w3^2 * rAB;
a_n_O4B = -w4^2 * rO4B;

% Resolvemos alphas 3 y 4
TerminoIndep_Acel1 = [(a_n_O4B(1) - aA(1) - a_n_BA(1));
                      (a_n_O4B(2) - aA(2) - a_n_BA(2))];
sol_alpha34 = Matriz1 \ TerminoIndep_Acel1;
alpha3 = sol_alpha34(1);
alpha4 = sol_alpha34(2);

aB_n = a_n_O4B;
aB_t = [-alpha4 * rO4B(2), alpha4 * rO4B(1)];
aB = aB_n + aB_t;

fprintf('aB (Normal)   : [%.4f, %.4f]\n', aB_n);
fprintf('aB (Tangencial): [%.4f, %.4f]\n', aB_t);
fprintf('aB (TOTAL)    : [%.4f, %.4f] m/s^2\n\n', aB);

% Aceleración de C
aC_n = -w4^2 * rO4C;
aC_t = [-alpha4 * rO4C(2), alpha4 * rO4C(1)];
aC = aC_n + aC_t;

fprintf('aC (Normal)   : [%.4f, %.4f]\n', aC_n);
fprintf('aC (Tangencial): [%.4f, %.4f]\n', aC_t);
fprintf('aC (TOTAL)    : [%.4f, %.4f] m/s^2\n\n', aC);

% Aceleración de D
a_n_DC  = -w5^2 * rCD;
a_n_O6D = -w6^2 * rO6D;

% Resolvemos alphas 5 y 6
TerminoIndep_Acel2 = [(a_n_O6D(1) - aC(1) - a_n_DC(1));
                      (a_n_O6D(2) - aC(2) - a_n_DC(2))];
sol_alpha56 = Matriz2 \ TerminoIndep_Acel2;
alpha5 = sol_alpha56(1);
alpha6 = sol_alpha56(2);

aD_n = a_n_O6D;
aD_t = [-alpha6 * rO6D(2), alpha6 * rO6D(1)];
aD = aD_n + aD_t;

fprintf('aD (Normal)   : [%.4f, %.4f]\n', aD_n);
fprintf('aD (Tangencial): [%.4f, %.4f]\n', aD_t);
fprintf('aD (TOTAL)    : [%.4f, %.4f] m/s^2\n\n', aD);

% Aceleración de E
aE_n = -w6^2 * rO6E;
aE_t = [-alpha6 * rO6E(2), alpha6 * rO6E(1)];
aE = aE_n + aE_t;

fprintf('aE (Normal)   : [%.4f, %.4f]\n', aE_n);
fprintf('aE (Tangencial): [%.4f, %.4f]\n', aE_t);
fprintf('aE (TOTAL)    : [%.4f, %.4f] m/s^2\n', aE);
fprintf('|aE| Magnitud : %.4f m/s^2\n', norm(aE));
