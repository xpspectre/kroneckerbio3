function sol = integrateMFK(m, con, opts)

% Constants
nx = m.nx;
nu = m.nu; 
omega = m.Compartments(1).Values;
S = m.S;

A1=m.A1; A2=m.A2; A3=m.A3; A4=m.A4; A5=m.A5; A6=m.A6;
D1=m.D1; D2=m.D2; D3=m.D3; D4=m.D4; D5=m.D5; D6=m.D6;

T_1_nx  = sparse(Tmatrix([1 nx]));
T_nx_nx = sparse(Tmatrix([nx nx]));

aux1 = A3*(kron(T_1_nx,eye(nu)));
aux2 = A4*(kron(T_1_nx,eye(nu)));

Inx=sparse(eye(nx));
Inu=sparse(eye(nu));

% Construct system
[der, jac] = constructSystem();


% Initial conditions
if opts.UseModelICs
    ic = m.x0;
else
    ic = con.x0;
end


% Append the initial conditions of Covariance Matrix
ic = cat(1,ic,opts.V0(:));


% Input
if opts.UseModelInputs
    u = m.u;
else
    u = con.u;
end

% Integrate x over time
sol = accumulateOde(der, jac, 0, con.tF, ic, u, con.Discontinuities, 1:nx, opts.RelTol, opts.AbsTol);
sol.u = u;
sol.C1 = m.C1;
sol.C2 = m.C2;
sol.c  = m.c;

% End of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% The system for integrating x and dxdT %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [der, jac] = constructSystem()
        
%        Ip = speye(nT);

        VStart = nx+1;
        VEnd   = nx+nx*nx;
        
        f       = m.f;
        dfdx    = m.dfdx;
        dfdT    = @dfdTSub;
        d2fdx2  = m.d2fdx2;
        d2fdTdx = @d2fdTdxSub;
        
        der = @derivative;
        jac = @jacobian;
        
        % Derivative of [x; dxdT] with respect to time
        function val = derivative(t, joint, u)
            u = u(t);            
            x = joint(1:nx);
            V = reshape(joint(VStart:VEnd), nx,nx);
            
            % and devide by volume to turn into concentrations
            x=x/omega;
            u=u/omega;
            V=V/(omega^2);
            
            
            
            dXdt = f(t, x*omega, u)/omega + A2*V(:);
            
            M =   A1    +    A2 * ( kron(eye(nx), x) + kron(x,eye(nx))) ...
                +    A3*kron(T_1_nx,eye(nu))  * kron(eye(nx),u) ...
                +    A4*kron(eye(nu),T_1_nx) * kron(u,eye(nx)) ;
            
%           M =     A1  +      A2 *     ( kron(Inx, x) + kron(x,Inx))     +                 aux1*(kron(Inx,u))         +                     aux2*kron(u,Inu) ;          
            
            
            Lamb = diag( D1*x + D2*kron(x,x) + D3*kron(u,x) + D4*kron(x,u) + D5*kron(u,u) + D6*u ) ;
                       
            
            dVdt = M*V + V*M' + (1/omega)*(S*Lamb*S');
            
            % Combine and change back to amounts and not concetrations
            val = [dXdt*(omega); dVdt(:)*(omega^2)];
            
            
           
        end
        
        % Jacobian of [X; V] derivative
        function val = jacobian(t, joint, u)
            u = u(t);            
            x = joint(1:nx); % x_
            V = reshape(joint(VStart:VEnd), nx,nx);
            
            % and devide by volume to turn into concentrations
             x=x/omega;
             u=u/omega;
             V=V/(omega^2);
             
             
            M =  A1 + A2 * ( kron(eye(nx), x) + kron(x,eye(nx))) +  A3*(kron(T_1_nx,eye(nu))*(kron(eye(nx),u)))  + A4*(kron(T_1_nx,eye(nu))*kron(u,eye(nx))) ;
           
            % d(dv/dt)/dx
            dMdx_dA2 = kron(vec(eye(nx)), eye(nx)) + kron(T_nx_nx, eye(nx)) * (kron(eye(nx),vec(eye(nx))));
            dMdx_dA2 = reshape(dMdx_dA2, [nx*nx nx nx]);
            
            dMdx = [];
                for i=1:nx
                   dMdx = cat(3,dMdx,A2*dMdx_dA2(:,:,i));
                end
                
            dMVdx = [];
                for i=1:nx
                   dMVdx = cat(3,dMVdx,dMdx(:,:,i)*V);
                end
            
            dVMtdx =[];
                for i=1:nx
                   dVMtdx = cat(3,dVMtdx,V*dMdx(:,:,i)');
                end
                
            dvRdx = D1 + D2 * ( kron(eye(nx), x) + kron(x,eye(nx))) +  D3*(kron(T_1_nx,eye(nu))*kron(u,eye(nx)))   + D4*(kron(T_1_nx,eye(nu))*(kron(eye(nx),u))) ;
            
            
            dSRStdx =[];
                for i=1:nx
                   dSRStdx = cat(3,dSRStdx,S*diag(dvRdx(:,i))*S');
                end
            
           
            ddvdtdx =  dMVdx + dVMtdx + (1/omega)*dSRStdx;
            ddvdtdx = reshape(ddvdtdx, [nx*nx nx]);
                
              
              
             % d(dVdt)/dv
             
             dVdv_dM = [];
             for i=1:nx^2
                 aux=zeros(nx,nx);
                 aux(i) = 1;
                 
                 dVdv_dM = cat(3, dVdv_dM, aux);
             end
             
             dMVdv = [];
             for i=1:nx^2
                 aux=zeros(nx,nx);
                 aux(i) = 1;
                 
                 dMVdv = cat(3, dMVdv, M*dVdv_dM(:,:,i));
             end
             
             
             dVMtdv = [];
             for i=1:nx^2
                 aux=zeros(nx,nx);
                 aux(i) = 1;
                 
                 dVMtdv = cat(3, dVMtdv, dVdv_dM(:,:,i)*M');
             end
             
             
            ddVdtdv = dMVdv + dVMtdv;
            ddVdtdv = reshape(ddVdtdv, [nx*nx nx*nx]);
             
            % Combine
            val = [m.dfdx(t, x*omega, u),  A2*eye(nx^2)/omega;
                       ddvdtdx*(omega), ddVdtdv];
                     
       
        end
        
        
        
    end

end