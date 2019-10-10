classdef CDMD_ADMM < handle
     
    properties
        param               % parameters
        
        % input
        X,Y                 % states/observables
        r                   % #modes
        
        % optimization
        rho
        steps               % #admm_steps
        tolabs,tolrel       % tolerances
        
        % stats
        verbose
        EG,PRES,DRES
    end
    
    methods
        function [ cdmd ] = CDMD_ADMM(X,Y,varargin)
            param = inputParser;
            addOptional(param,'r',0);
            addOptional(param,'rho',10);
            addOptional(param,'steps',10000);
            addOptional(param,'tolabs',1e-8);
            addOptional(param,'tolrel',1e-6);
            addOptional(param,'verbose',0);
            
            parse(param,varargin{:});
            
            cdmd.r      = param.Results.r;
            cdmd.rho    = param.Results.rho;
            cdmd.steps  = param.Results.steps;
            cdmd.tolabs = param.Results.tolabs;
            cdmd.tolrel = param.Results.tolrel;
            cdmd.verbose= param.Results.verbose;
            
            cdmd.X      = X;
            cdmd.Y      = Y;
            
            % stats
            cdmd.EG = []; cdmd.PRES = []; cdmd.DRES = [];
        end
        
        function [EV,EL,A,B,U] = compCDMD( cdmd )
            
            % prerequisites
            ID = eye( cdmd.r );
            rh = cdmd.rho; rh_as = 5; rh_a = 2;
            
            oX = cdmd.X; oY = cdmd.Y; oZ = [oX oY(:,end)]; I = 1:cdmd.r; 
            
            [U,S,V] = svd(oZ,'econ');
            
            U = U(:,I); V = V(:,I); S = S(I,I); IS = diag(1.0 ./ diag(S));
            
            sX = U'*oX; sY = U'*oY;
                        
            % initialization for A,B,Q1,Q2,rho
            A = sY*pinv(sX);
            B = sX*pinv(sY);
            Q1 = zeros( size(A) );
            Q2 = zeros( size(A) );
            
            % iterate
            for st = 1:cdmd.steps
                                
                % solve for A
                Ao = A;
                M1 = rh*(B'*B);
                M2 = sX*sX' + rh*(B*B');
                M3 = sY*sX' + 2*rh*B' - rh*Q1*B' - rh*B'*Q2;
                A = sylvester(M1,M2,M3);
                
                % solve for B
                Bo = B;
                N1 = rh*(A'*A);
                N2 = sY*sY' + rh*(A*A');
                N3 = sX*sY' + 2*rh*A' - rh*A'*Q1 - rh*Q2*A';
                B = sylvester(N1,N2,N3);
                
                % update Q
                Q1 = Q1 + A*B - ID;
                Q2 = Q2 + B*A - ID;
                
                % stats
                if mod(st,10) == 0
                    
                    Ad = A*sX-sY; Bd = sX-B*sY;
                    eg = .5*norm(Ad,'fro')^2 + .5*norm(Bd,'fro')^2;
                    pres = .5*norm(A*B-ID,'fro') + .5*norm(B*A-ID,'fro');
                    dres = .5*norm(rh*[A-Ao; B-Bo],'fro');

                    if pres > rh_as * dres
                        rh = rh * rh_a;
                        Q1 = Q1 / rh_a;
                        Q2 = Q2 / rh_a;
                    elseif dres > rh_as * pres
                        rh = rh / rh_a;
                        Q1 = Q1 * rh_a;
                        Q2 = Q2 * rh_a;
                    end
                    
                    % stopping condition 
                    epspri  = sqrt(cdmd.r)*cdmd.tolabs + cdmd.tolrel*max(norm(A*B,'fro'),max(norm(B*A,'fro'),norm(ID,'fro')));
                    epsdual  = sqrt(size(Q1,1)+size(Q2,1))*cdmd.tolabs + cdmd.tolrel*norm(rh*[Q1; Q2],'fro'); 
                    
                    if mod(st,100) == 0 && cdmd.verbose == 1
                        fprintf('%d o %.08f | pr %.8f dr %.8f | rho %.8f | peps %.8f deps %.8f\n', ...
                                st, eg, pres, dres, rh, epspri, epsdual );
                    end
                    cdmd.EG = [cdmd.EG; eg]; cdmd.PRES = [cdmd.PRES; pres]; cdmd.DRES = [cdmd.DRES; dres];
                    
                    % stopping condition
                    if pres < epspri && dres < epsdual; break; end
                end
            end
            if (st == cdmd.steps); fprintf('\ndid not converge!\n'); end
            
            % eigendecompose the integrated system
            [ev,el] = eig(A);
            
            EV = oY*V(2:end,:)*IS*ev;
            EL = el;
            
        end
    end
    
    methods (Static)
    end
    
end