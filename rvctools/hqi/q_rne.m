% NEW_EULER RECURTION
%standart D_H in P.194 of robotics vision.




function tau = q_rne(robot, q, qd,qdd)
      z0 = [0;0;1];
      grav = robot.gravity;   % default gravity from the object.
      fext = zeros(6, 1);
      n = robot.n;

      Rb=t2r(robot.base)';
      w=zeros(3,1);%base not move.
      wd=zeros(3,1);
      vd=Rb*grav;%base have acc from gravity.
  
  for i=1:n
      link = robot.links(i);
      if link.RP == 'R'
          d = link.d;
      else
          d = q(i);
      end
      alpha = link.alpha;
      % O_{j-1} to O_j in {j}, negative inverse of link xform
      pstar = [link.a; d*sin(alpha); d*cos(alpha)];
      pstarm{i}= pstar;
      R{i}=t2r(robot.links(i).A(q(i)));%{i} in {i-1} standart D_H
  end
  %the forward recurtion
  for i=1:n
      r=robot.links(i).r;
      
      if robot.links(i).RP == 'R'  %rotation
      wd=R{i}.'*(wd+z0*qdd(i)+cross(w,z0*qd(i)));
      w=R{i}.'*(w+z0*qd(i));
      vd=R{i}.'*vd+cross(wd,pstarm{i})+cross(w,cross(w,pstarm{i})); %or   vd = cross(wd,pstar) +cross(w, cross(w,pstar)) +Rt*vd
      
      elseif robot.links(i).RP== 'P'  %prismatic
      w=R{i}'*w;
      wd=R{i}*wd;
      
      vd=R{i}.'*vd+cross(wd,pstarm{i})+cross(w,cross(w,pstarm{i}))+2*cross(w,R{i}.'*(z0*qd(i))); %ware
      else
          error('wrong link type');
      end
      vcd=cross(wd,r')+cross(w,cross(w,r'))+vd;
      F{i}=robot.links(i).m*vcd;
      N{i}=robot.links(i).I*wd+cross(w,robot.links(i).I*w); 
  end
      f=fext(1:3);
      nn=fext(4:6);
  
% the backward recursion
  for i=n:-1:1
      link=robot.links(i);
      r=link.r;
      if i==n
          Rm=eye(3,3);
      else
          Rm=R{i+1};
      end
      nn=Rm*(nn+cross(Rm'*pstarm{i},f))+cross(pstarm{i}+r',F{i})+N{i};
      f=Rm*f+F{i};
      
      Rm=R{i};
      if link.RP == 'R'
            % revolute
            t = nn.'*(Rm.'*z0) + ...
                link.G^2 * link.Jm*qdd(i) - ...
                 link.friction(qd(i));
            tau(i) = t;
        else
            % prismatic
            t = f.'*(Rm.'*z0) + ...
                link.G^2 * link.Jm*qdd(i) - ...
                link.friction(qd(i));
            tau(i) = t;
      end
  end

   