clear;
%------------- prm for planar------------------------
q=[pi/6 pi/6]; qd=[0 0]; qdd=[0 0];fext={};
n=2;
m=[2.247,11.078];
L=[0.25 0.25];
centroid={[0.11652   0] [0.11661  0.0000]};
%-------------- 1 -----------------
%step 1 :get module
robot.NB = n;
robot.parent = [0:n-1];

for i = 1:n
  robot.jtype{i} = 'r';
  if i == 1
    robot.Xtree{i} = plnr( 0, [0 0]);
  else
    robot.Xtree{i} = plnr( 0, [L 0]);
  end
  robot.I{i} = mcI( m(i), centroid{i}, 1/12 );
end
model=robot;

%------------- 2 ------------------
%step 2 :new_euler
%tau=ID(model,q,qd,qdd)
a_grav = get_gravity(model);

for i = 1:model.NB
  [ XJ, S{i} ] = jcalc( model.jtype{i}, q(i) );
  vJ = S{i}*qd(i);
  Xup{i} = XJ * model.Xtree{i};
%   XJ
%   model.Xtree{i}
%    Xup{i}
  if model.parent(i) == 0
    v{i} = vJ;
    a{i} = Xup{i}*(-a_grav) + S{i}*qdd(i);
  else
    v{i} = Xup{i}*v{model.parent(i)} + vJ;
    a{i} = Xup{i}*a{model.parent(i)} + S{i}*qdd(i) + crm(v{i})*vJ;
  end
  f{i} = model.I{i}*a{i} + crf(v{i})*model.I{i}*v{i};
end
 S{i}
  f{i}
for i = model.NB:-1:1
  tau(i,1) = S{i}' * f{i};
  if model.parent(i) ~= 0
    f{model.parent(i)} = f{model.parent(i)} + Xup{i}'*f{i};
  end
end
tau