mesh = abaqus_reader('mesh-t3-r2-hchintak.inp');
p.faces = mesh.conn';
p.facecolor = '#DDDDFF';
p.vertices = mesh.x';
p.linewidth = 0.25;
p.edgealpha = 0.25;
clf;
patch(p);
axis equal;
order = height(mesh.conn)/3;
%%
K_e = calc_stiffness_e(mesh);
F_e = calc_force_e(mesh,3,force_conn);
U_e = calc_disp_e(K_e,F_e,nodes2);
F_h = calc_force_heat(mesh,U_e);
K_h = calc_stiffness_h(mesh);
U_h = calc_disp_h(K_h,F_h,nodes1,nodes2);
total = sum(F_h)
max_heat = max(s)
max_temp = max(U_h)
%% calculate disp heat
function [U_h] = calc_disp_h(K_h,F_h,nodes1,nodes2)
    bc = [nodes1 nodes2];
    K_h(bc,:) = 0;
    K_h(:,bc) = 0;
    K_h(bc,bc) = eye(length(bc));
    F_h(bc) = 0;
    U_h = zeros(length(F_h),1);
    U_h = K_h\F_h;
    assignin('base',"U_h",U_h)
end
%% calculate stiffness heat
function[K_h] = calc_stiffness_h(mesh)
    % electrical resistivity k = 1/ro = 1/(460*10^-6)
    % thermal resistivity k = 0.37
    k = 0.37;
    nne = size(mesh.conn,1);
    if nne == 3
        shape = @shapeT3;
        qpts = [1/6 2/3 1/6;
               1/6 1/6 2/3;
               1/6 1/6 1/6];
    elseif nne == 6
        shape = @shapeT6;
        qpts = [1/6 2/3 1/6;
               1/6 1/6 2/3;
               1/6 1/6 1/6];
    end
    K_h = sparse(length(mesh.x),length(mesh.x));
    for c = mesh.conn    
        xe = mesh.x(:,c);
        Ke= zeros(length(c));
        for q = qpts
            [N,dNdp] = shape([q(1),q(2)]);
            x = xe*N;
            J = xe*dNdp;
            B = dNdp/J;
            Ke = Ke + B * k*eye(2) * B'* det(J)*q(end);
        end
        K_h(c,c) = K_h(c,c) + Ke;
    end
    assignin('base','K_h',K_h)
end
%% calculate force heat
function [F_h] = calc_force_heat(mesh,phi)
    nne = size(mesh.conn,1);
    k = 1/(460*10^-6);
    if nne == 3
        shape = @shapeT3;
        qpts = [1/6 2/3 1/6;
               1/6 1/6 2/3;
               1/6 1/6 1/6];
    elseif nne == 6
        shape = @shapeT6;
        qpts = [1/6 2/3 1/6;
               1/6 1/6 2/3;
               1/6 1/6 1/6];
    end
    F_h = zeros(length(mesh.x),1);
    C_d = zeros(length(mesh.conn),1);
    J_h = zeros(length(mesh.conn),1);
    s = zeros(length(mesh.x),1);
    i = 1;
    for c = mesh.conn
        xe = mesh.x(:,c); 
        Ae=0;
        C = 0;
        Je = 0;
        for q = qpts
            [N,dNdp] = shape([q(1),q(2)]);
            J = xe*dNdp;
            B = dNdp/J;
            E = -B'*phi(c);
            Je = Je + k*E*det(J)*q(end);
            Ja = k*E;
            J_h(c) = E'*Ja*det(J)*q(end);
            Ae = Ae + det(J)*q(end);
            F_h(c) = F_h(c) + N*E'*Ja*det(J)*q(end);
            s(c) = dot(Ja,E);
        end
        C = sqrt(Je(1)^2+Je(2)^2);
        C_d(i) = C/Ae;
        i = i+1;
    end
    assignin("base","C_d",C_d)
    assignin("base","F_h",F_h)
    assignin("base","J_h",J_h)
    assignin("base","s",s)
end
%% calculate disp electric
function [U_e] = calc_disp_e(K_e,F_e,nodes2)
    bc = nodes2;
    K_e(bc,:) = 0;
    K_e(:,bc) = 0;
    K_e(bc,bc) = eye(length(bc));
    F(bc) = 0;
    U_e = zeros(length(F_e),1);
    U_e = K_e\F_e;
    assignin('base',"U_e",U_e)
end
%% calculate force electric
function [F_e] = calc_force_e(mesh,p,force_conn)
    nne = size(mesh.conn,1);
    if nne == 3 
        shape = @shapelin;
        qpts = [0;2];
    elseif nne == 6
        shape = @shapequad;
        qpts = [sqrt(1/3) -sqrt(1/3);
               1 1];
    end
    F_e = sparse(length(mesh.x),1);
    for c = force_conn
        xe = mesh.x(2,c);
        for q = qpts
            [N, dNdp] = shape(q(1));
            J = xe*dNdp;
            F_e(c) = F_e(c) + N * p * J*q(end);
        end
    end
    assignin('base', 'F_e', F_e);
end
%% calculate_stiffness_e
function[K_e] = calc_stiffness_e(mesh)
    % electrical resistivity k = 1/ro = 1/(460*10^-6)
    % thermal resistivity k = 0.37
    k = 1/(460*10^-6);
    nne = size(mesh.conn,1);
    if nne == 3
        shape = @shapeT3;
        qpts = [1/6 2/3 1/6;
               1/6 1/6 2/3;
               1/6 1/6 1/6];
    elseif nne == 6
        shape = @shapeT6;
        qpts = [1/6 2/3 1/6;
               1/6 1/6 2/3;
               1/6 1/6 1/6];
    end
    K_e = sparse(length(mesh.x),length(mesh.x));
    for c = mesh.conn    
        xe = mesh.x(:,c);
        Ke= zeros(length(c));
        for q = qpts
            [N,dNdp] = shape([q(1),q(2)]);
            x = xe*N;
            J = xe*dNdp;
            B = dNdp/J;
            Ke = Ke + B * k*eye(2) * B'* det(J)*q(end);
        end
        K_e(c,c) = K_e(c,c) + Ke;
    end
    assignin('base','K_e',K_e)
end

function [N, dNdp] = shapelin(p)
    N = [(1-p(1))/2;(1+p(1))/2];
    dNdp = [-0.5;0.5];
end 

function [N, dNdp] = shapequad(p)
    N = [0.5*p*(p-1.0); 1.0-p*p; 0.5*p*(p+1.0)];
    dNdp = [ p - 0.5; -2.0*p; p + 0.5];
end
function [N, dNdp] = shapeT3(p)
    N = [p(1);p(2);1-p(1)-p(2)];
    dNdp = [1 0;0 1;-1 -1];
end 

function [N, dNdp] = shapeT6(p) 
    N = [p(1)*(2*p(1)-1);
         p(2)*(2*p(2)-1);
         (1-p(1)-p(2))*(2*(1-p(1)-p(2))-1);
         4*p(1)*p(2);
         4*p(2)*(1-p(1)-p(2));
         4*p(1)*(1-p(1)-p(2))];
    dNdp = [4*p(1)-1 0;
            0 4*p(2)-1;
            4*p(1)+4*p(2)-3 4*p(2)+4*p(1)-3;
            4*p(2) 4*p(1);
            -4*p(2) -4*(2*p(2)+p(1)-1);
            -4*(2*p(1)+p(2)-1) -4*p(1)];
end 