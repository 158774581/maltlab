% syms q1 q2 q1d q2d q1dd q2dd m1 m2 L1 L2 r1 r2 r3 R1 R2 R3 I1 I2 Jm1...
%     Jm2 Tc11 Tc12 Tc21 Tc22 B1 B2 G1 G2 real;
% L=[L L2];
% L(1) = Revolute('d', 0, 'a', L1, 'alpha', 0, ...
%     'I', [0, 0, I1, 0, 0, 0], ...
%     'r', [r1,  r2,  r3], ...
%     'm', m1, ...
%     'Jm', Jm1, ...
%     'G', G1, ...
%     'B', B1, ...
%     'Tc', [Tc21 Tc12], ...
%     'qlim', [-90 90]*deg );
% 
% L(2) = Revolute('d', 0, 'a', L2, 'alpha', 0, ...
%     'I', [0, 0, I2, 0, 0, 0], ...
%     'r', [R1,R2,R3], ...
%     'm', m2, ...
%     'Jm', Jm2, ...
%     'G', G2, ...
%     'B', B2, ...
%     'Tc', [Tc21 Tc22], ...
%     'qlim', [-160 160]*deg );



%-----------------trajectory------------------
qz = [0 pi/3 ]; % zero angles, L shaped pose
qr = [0 pi/6 ]; % ready pose, arm up
qm = [0 0 0];
t = [0:16/36:16-15/36]; 		% create time vector
t = [0:.056:2]; 		% create time vector
[q,qd,qdd] = jtraj(qz, qr, t); % compute joint coordinate trajectory


scara.inertia(q(1,1:2))*qdd(18,1:2)'
scara.coriolis(q(18,1:2),qd(18,1:2))*qd(18,1:2)'
tau(18,1:2)'-scara.inertia(q(1,1:2))*qdd(18,1:2)'-scara.coriolis(q(18,1:2),qd(18,1:2))*qd(18,1:2)'
plot(t, tau(:,1:2)); xlabel('Time (s)'); ylabel('Joint torque (Nm)');
figure(2)
plot(qdd);