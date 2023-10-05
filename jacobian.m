syms phi(t) theta(t) beta(t) r(t) R(t) t;
syms s; % Generalised co-ordinate(s)

l_1 = 490.0; % [mm]
l_2 = 160.0; % [mm]
l_3 = 400.0;   % [mm] - can be any value within reason
l_4 = (l_3^2 - l_2^2)^0.5; % [mm] - pythagoras thereom

v = 964.0; % [mm/s] - velocity of the acuator
s(t) = v * t + 390; % [mm] - displacement of acuator

phi(t) = acosd((l_1^2 + l_2^2 - s^2) / (2 * l_1 * l_2));
theta(t) = acosd((s^2 + l_1^2 - l_2^2) / (2 * s * l_1));
beta(t) = acosd((s^2 + l_2^2 - l_1^2) / (2 * s * l_2));

A = theta(0);
B = acosd((l_2^2 + l_3^2 - l_4^2) / (2 * l_2 * l_3));
% Assumption - Not neccessary; makes calculations slightly more simple.
C = 90; 

% Rotates the final matrix to give values that make more sense.
rotation_angle = -90;
rotation_matrix = [cosd(rotation_angle), -sind(rotation_angle);
                   sind(rotation_angle), cosd(rotation_angle)];

% Assumption - acuator is parallel to the x-axis at t = 0.
r_1_horz = l_1 * [cosd(A); 
                  sind(A)];
r_1 = simplify(rotation_matrix * r_1_horz);

r_2_horz(t) = -l_2 * [cosd(phi(t) + A);
                      sind(phi(t) + A)];
r_2 = simplify(rotation_matrix * r_2_horz);

r_3_horz(t) = -l_3 * [cosd(phi(t) + A + B);
                      sind(phi(t) + A + B)];
r_3(t) = simplify(rotation_matrix * r_3_horz(t));

R_1_horz(t) = s * [cosd(-(theta(t) - A)); 
                      sind(-(theta(t) - A))];
R_1(t) = simplify(rotation_matrix * R_1_horz(t));

R_2_horz(t) = -l_4 * [cosd((360 - C - beta(t) - theta(t) + A)); 
                      sind((360 - C - beta(t) - theta(t) + A))];
R_2(t) = simplify(rotation_matrix * R_2_horz(t));

% tip of the acuator
r(t) = r_1 + r_2;  
R(t) = R_1(t);

vector_closure = r(t) - R(t);
q_set = t;
dq_dt = diff(q_set, t);
dJ_dt = simplify(diff(vector_closure, t))
dJ_dq = simplify(dJ_dt / dq_dt);

% vector closure equation r - R = 0
r = round(eval(r(0.05)), 4);
R = round(eval(R(0.05)), 4);    







