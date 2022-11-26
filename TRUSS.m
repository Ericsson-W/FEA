classdef TRUSS
    % CIVE95006: Computational Methods II
    % CLASS TRUSS
    % Ching Wai Wong
    % 01704926
    
    properties
        % Type of test
        type
        % failure mode, indicates how element failed 1 for plastic squash 2
        % for buckling
        mode
        % L        
        L % will be initial input
        
        % P
        P % will be intial input
        Pcrit % MNA critical load
        Pcr % GMNA critical load
        % nodal coordinates
        coords
        
        % elemental connectivity
        ELEMENTS
              
        % degree of freedoms
        dofs
        
        % number of nodes
        nodes
        
        % nodal force vector
        F
        
        % free dofs
        dofs_free
        
        % restrained dofs
        dofs_restrained
        
        % λ vector
        lamda_total
        dlamda
        lamda
        
        % Deflection vector
        U
        uR
        uF
        %displacement at E
        dU_E 
        % E
        E % will be initial input
        
        % R
        R % will be initial input
        
        % EA
        EA % calculated
        
        % sigma
        sigma % will be initial input
        
        % matrix
        K
        KFF
        KFR
        KRR
        KRF
        fF
        fR
        
        % Axial forces of current iteration
        F_axial
        % Axial forces of all iterations
        F_axialtotal
       
        % element indication vector
        indicator % shows which elements have failed
        pointer  % shows which element has just failed in the current iteration
        
        % increment counter
        counter
        % lamda of element if buckling happened
        bucklinglamda
    end
    
    methods
        %% Assembling Truss
        function obj=TRUSS(coords,ELEMENTS,dofs,dofs_free,dofs_restrained,E,R,P,L,sigma,lamda_total,type,counter)
            
            % Assign values to obj
            obj.coords=coords;
            obj.ELEMENTS=ELEMENTS;
            obj.dofs=dofs;
            obj.dofs_free=dofs_free;
            obj.dofs_restrained=dofs_restrained;
            obj.E=E;
            obj.R=R;
            obj.P=P;
            obj.L=L;
            obj.sigma=sigma;
            obj.lamda_total=lamda_total;
            obj.type=type;
            obj.counter=counter;
            obj.nodes = size(obj.coords,1);
            obj.F = zeros(2*obj.nodes,1); % initialising en empty a 20x1 column vector for convenience
            obj.indicator=ones(size(obj.ELEMENTS,1),1);
        end
        

        
        %% Assemble stiffness matrix        
        function obj= stiffness(obj)            
            elements = size(obj.ELEMENTS,1); 
            k=elements;
            obj.EA = pi*obj.R^2*obj.E;

            elements = size(obj.ELEMENTS,1); 
            % Constructing the global stiffness matrix
            obj.K = zeros(2*obj.nodes); % initialisaing en empty 20 x 20 matrix; 10 nodes at 2 dofs/node
            for EL = 1:elements % loop through all elements & build stiffness matrix
                if obj.indicator(EL)==1 % element still here
                    n1 = obj.ELEMENTS(EL,1); n2 = obj.ELEMENTS(EL,2); % identify element node numbers
                    x1 = obj.coords(n1,1); y1 = obj.coords(n1,2); % element node 1 - x,y coordinates
                    x2 = obj.coords(n2,1); y2 = obj.coords(n2,2); % element node 2 - x,y coordinates
                    dof11 = obj.dofs(n1,1); dof12 = obj.dofs(n1,2); % element node 1 - dofs
                    dof21 = obj.dofs(n2,1); dof22 = obj.dofs(n2,2); % element node 2 - dofs
                    alpha = atan2(y2-y1,x2-x1); % angle of inclination relative to the POSITIVE direction of the x axis
                    c = cos(alpha); c2 = c*c; s = sin(alpha); s2 = s*s; cs = c*s; % angle parameters
                    Le = sqrt( (x2 - x1)^2 + (y2 - y1)^2 ); % element length
                    ke = obj.EA/Le; % element axial stiffness


                    % Updating global stiffness matrix [K] coefficients 
                    % Row 1 - element dof11
                    obj.K(dof11,dof11) = obj.K(dof11,dof11) + ke*c2; % Col 1 - element dof11
                    obj.K(dof11,dof12) = obj.K(dof11,dof12) + ke*cs; % Col 2 - element dof12
                    obj.K(dof11,dof21) = obj.K(dof11,dof21) - ke*c2; % Col 3 - element dof21
                    obj.K(dof11,dof22) = obj.K(dof11,dof22) - ke*cs; % Col 4 - element dof22

                    % Row 2 - element dof12
                    obj.K(dof12,dof11) = obj.K(dof12,dof11) + ke*cs; % Col 1 - element dof11
                    obj.K(dof12,dof12) = obj.K(dof12,dof12) + ke*s2; % Col 2 - element dof12
                    obj.K(dof12,dof21) = obj.K(dof12,dof21) - ke*cs; % Col 3 - element dof21
                    obj.K(dof12,dof22) = obj.K(dof12,dof22) - ke*s2; % Col 4 - element dof22

                    % Row 3 - element dof21
                    obj.K(dof21,dof11) = obj.K(dof21,dof11) - ke*c2; % Col 1 - element dof11
                    obj.K(dof21,dof12) = obj.K(dof21,dof12) - ke*cs; % Col 2 - element dof12
                    obj.K(dof21,dof21) = obj.K(dof21,dof21) + ke*c2; % Col 3 - element dof21
                    obj.K(dof21,dof22) = obj.K(dof21,dof22) + ke*cs; % Col 4 - element dof22

                    % Row 4 - element dof22
                    obj.K(dof22,dof11) = obj.K(dof22,dof11) - ke*cs; % Col 1 - element dof11
                    obj.K(dof22,dof12) = obj.K(dof22,dof12) - ke*s2; % Col 2 - element dof12
                    obj.K(dof22,dof21) = obj.K(dof22,dof21) + ke*cs; % Col 3 - element dof21
                    obj.K(dof22,dof22) = obj.K(dof22,dof22) + ke*s2; % Col 4 - element dof22                
                else % element failed (indicator=0)
                    fprintf('we dont need to consider failed element %g in K matrix \n',EL);
                end
            end

            % Constructing the global nodal force vector 
            if obj.counter==1
                obj.F(14)=obj.F(14)-3*obj.P*obj.lamda_total(end);
                obj.F(16)=obj.F(16)-1.5*obj.P*obj.lamda_total(end);
                obj.F(18)=obj.F(18)-obj.P*obj.lamda_total(end);
            else
                % since in the Truss1.critical function we scaled all the
                % loads by lamda, we dont need to do this here again
            end
                
                      


            % Specification of submatrices directly
            obj.KRR = obj.K(obj.dofs_restrained,obj.dofs_restrained);
            obj.KRF = obj.K(obj.dofs_restrained,obj.dofs_free);
            obj.KFR = obj.K(obj.dofs_free,obj.dofs_restrained);
            obj.KFF = obj.K(obj.dofs_free,obj.dofs_free);
            obj.fF = obj.F(obj.dofs_free);


            % Solution for the unknown nodal dofs
            obj.uR = [0 0 0 0 0 0 0 0]'; 
            obj.uF = obj.KFF\(obj.fF - obj.KFR*obj.uR); 
            obj.U(obj.dofs_restrained)=obj.uR;
            obj.U(obj.dofs_free)=obj.uF;

            % Solution for the unknown reactions
            obj.fR = obj.KRF*obj.uF + obj.KRR*obj.uR; 
            obj.F(obj.dofs_restrained)=obj.fR;
            obj.F(obj.dofs_free)=obj.fF;            
            
            % store the displacement at joint E 
            obj.dU_E(obj.counter)=obj.uF(4);
        end
        %% plotting
        function obj= plot(obj) 
            elements = size(obj.ELEMENTS,1);
            new_coords = zeros(size(obj.coords)); 
            amp_coords = zeros(size(obj.coords)); 
            amp = 2; % amplification factor for plotting purposes only 
            for I = 1:size(obj.coords,1)
                for J = 1:size(obj.coords,2)    
                    amp_coords(I,J) = obj.coords(I,J) + obj.U(obj.dofs(I,J))*amp;
                    new_coords(I,J) = obj.coords(I,J) + obj.U(obj.dofs(I,J));
                end
            end

            % Plotting
            figure; hold all; grid on; tol = 1e-3;
            xmin = min(obj.coords(:,1)); xmax = max(obj.coords(:,1)); difx = xmax - xmin;
            ymin = min(obj.coords(:,2)); ymax = max(obj.coords(:,2)); dify = ymax - ymin; fac = 0.25;
            axis([xmin-difx*fac  xmax+difx*fac  ymin-dify*fac  ymax+dify*fac]); 
            if obj.counter==1 % first iteration
                for EL = 1:elements
                    n1 = obj.ELEMENTS(EL,1); n2 = obj.ELEMENTS(EL,2); % identify element node numbers

                    % Plotting original structure
                    x1 = obj.coords(n1,1); y1 = obj.coords(n1,2); % element node 1 - x,y original coordinates
                    x2 = obj.coords(n2,1); y2 = obj.coords(n2,2); % element node 2 - x,y original coordinates
                    Le = sqrt( (x2 - x1)^2 + (y2 - y1)^2 ); % element length
                    ke = obj.EA/Le; % element axial stiffness
                    alpha = atan2(y2-y1,x2-x1); % angle of inclination relative to the POSITIVE x axis direction
                    plot([x1,x2],[y1,y2],'Color',[0.5 0.5 0.5],'Linewidth',3); 

                    % Check on changes in member lengths and plotting amplified deformed structure
                    x1_amp = amp_coords(n1,1); y1_amp = amp_coords(n1,2); % element node 1 - x,y amplified deformed coordinates
                    x2_amp = amp_coords(n2,1); y2_amp = amp_coords(n2,2); % element node 2 - x,y amplified deformed coordinates     
                    x1_new = new_coords(n1,1); y1_new = new_coords(n1,2); % element node 1 - x,y actual deformed coordinates
                    x2_new = new_coords(n2,1); y2_new = new_coords(n2,2); % element node 2 - x,y actual deformed coordinates    
                    u1 = x1_new - x1; v1 = y1_new - y1; u2 = x2_new - x2; v2 = y2_new - y2; % reconstruction of element global dofs
                    
                    up1 = cos(alpha)*u1 + sin(alpha)*v1; up2 = cos(alpha)*u2 + sin(alpha)*v2; dup = up2 - up1; % reconstruction of element local dofs     
                    % calculate axial force
                    obj.F_axial(EL)=ke*dup;

                    if dup < -tol % element length has decreased - member in compression
                        plot([x1_amp,x2_amp],[y1_amp,y2_amp],'b','Linewidth',3); % blue colour
                    elseif dup > tol % element length as increased - member in tension
                        plot([x1_amp,x2_amp],[y1_amp,y2_amp],'r','Linewidth',3); % red colour
                    else % no change in element length
                        plot([x1_amp,x2_amp],[y1_amp,y2_amp],'k','Linewidth',3); % black colour
                    end
 
                    % Plotting nodes 
                    plot(x1,y1,'ko','Markersize',7,'MarkerFaceColor','w');  
                    plot(x2,y2,'ko','Markersize',7,'MarkerFaceColor','w');
                    plot(x1_amp,y1_amp,'ko','Markersize',7,'MarkerFaceColor','y');  
                    plot(x2_amp,y2_amp,'ko','Markersize',7,'MarkerFaceColor','y');
                end
            else
               for EL = 1:elements
                   if obj.indicator(EL)==1 %element still here
                        n1 = obj.ELEMENTS(EL,1); n2 = obj.ELEMENTS(EL,2); % identify element node numbers

                        % Plotting original structure
                        x1 = obj.coords(n1,1); y1 = obj.coords(n1,2); % element node 1 - x,y original coordinates
                        x2 = obj.coords(n2,1); y2 = obj.coords(n2,2); % element node 2 - x,y original coordinates
                        Le = sqrt( (x2 - x1)^2 + (y2 - y1)^2 ); % element length
                        ke = obj.EA/Le; % element axial stiffness
                        alpha = atan2(y2-y1,x2-x1); % angle of inclination relative to the POSITIVE x axis direction
                        plot([x1,x2],[y1,y2],'Color',[0.5 0.5 0.5],'Linewidth',3); 

                        % Check on changes in member lengths and plotting amplified deformed structure
                        x1_amp = amp_coords(n1,1); y1_amp = amp_coords(n1,2); % element node 1 - x,y amplified deformed coordinates
                        x2_amp = amp_coords(n2,1); y2_amp = amp_coords(n2,2); % element node 2 - x,y amplified deformed coordinates     
                        x1_new = new_coords(n1,1); y1_new = new_coords(n1,2); % element node 1 - x,y actual deformed coordinates
                        x2_new = new_coords(n2,1); y2_new = new_coords(n2,2); % element node 2 - x,y actual deformed coordinates    
                        u1 = x1_new - x1; v1 = y1_new - y1; u2 = x2_new - x2; v2 = y2_new - y2; % reconstruction of element global dofs
                        up1 = cos(alpha)*u1 + sin(alpha)*v1; up2 = cos(alpha)*u2 + sin(alpha)*v2; dup = up2 - up1; % reconstruction of element local dofs     

                        % calculate axial force 
                        obj.F_axial(EL)=ke*dup;
                        
                        if dup < -tol % element length has decreased - member in compression
                            plot([x1_amp,x2_amp],[y1_amp,y2_amp],'b','Linewidth',3); % blue colour
                        elseif dup > tol % element length as increased - member in tension
                            plot([x1_amp,x2_amp],[y1_amp,y2_amp],'r','Linewidth',3); % red colour
                        else % no change in element length
                            plot([x1_amp,x2_amp],[y1_amp,y2_amp],'k','Linewidth',3); % black colour
                        end

                        % Plotting nodes last!
                        plot(x1,y1,'ko','Markersize',7,'MarkerFaceColor','w');  
                        plot(x2,y2,'ko','Markersize',7,'MarkerFaceColor','w');
                        plot(x1_amp,y1_amp,'ko','Markersize',7,'MarkerFaceColor','y');  
                        plot(x2_amp,y2_amp,'ko','Markersize',7,'MarkerFaceColor','y');
                   else
                        obj.F_axial(EL)=0; % failed element
                        n1 = obj.ELEMENTS(EL,1); n2 = obj.ELEMENTS(EL,2); % identify element node numbers

                        % Plotting original structure
                        x1 = obj.coords(n1,1); y1 = obj.coords(n1,2); % element node 1 - x,y original coordinates
                        x2 = obj.coords(n2,1); y2 = obj.coords(n2,2); % element node 2 - x,y original coordinates
                        plot(x1,y1,'ko','Markersize',7,'MarkerFaceColor','w');  
                        
                        x1_amp = amp_coords(n1,1); y1_amp = amp_coords(n1,2); % element node 1 - x,y amplified deformed coordinates
                        x2_amp = amp_coords(n2,1); y2_amp = amp_coords(n2,2); % element node 2 - x,y amplified deformed coordinates                        
                        % Plot nodes
                        plot(x2,y2,'ko','Markersize',7,'MarkerFaceColor','w');
                        plot(x1_amp,y1_amp,'ko','Markersize',7,'MarkerFaceColor','y');  
                        plot(x2_amp,y2_amp,'ko','Markersize',7,'MarkerFaceColor','y');
                        
                        fprintf('no plot since element %g has failed \n',EL);
                   end
                end
            end
            xlabel('x coordinate','FontSize',16);
            ylabel('y coordinate','FontSize',16);


            % Printing computed dofs & reactions 
            for dof = 1:length(obj.uF)
%                 disp(['The value of dof ',num2str(obj.dofs_free(dof)),' is ',num2str(obj.uF(dof))]);    
            end
            disp(' ');
            for react = 1:length(obj.fR)
%                 disp(['The value of the reaction at dof ',num2str(obj.dofs_restrained(react)),' is ',num2str(obj.fR(react))]);
            end
            disp(' '); disp('Vertical equilibrium check:');
            disp(['Total vertical reactions = ',num2str(obj.fR(2) + obj.fR(4)+obj.fR(6)+obj.fR(8))]);
            disp(['Total applied vertical loads = ',num2str(-5.5*obj.P*sum(obj.lamda_total))]);
            if abs(obj.fR(2) + obj.fR(4) - obj.P) < 1e-6; disp('Ok.'); end
            disp(' '); disp('Horizontal equilibrium check:');
            disp(['Total horizontal reactions = ',num2str(obj.fR(1) + obj.fR(3)+obj.fR(5)+obj.fR(7))]);
            disp('Total applied horizontal loads = 0');
            if abs(obj.fR(1) + obj.fR(3)+obj.fR(5)+obj.fR(7)) < 1e-6; disp('Ok.'); end
        end
        
        
        %%
        function obj=critical(obj) % find critical criterion
            obj.Pcrit=pi*obj.R^2*obj.sigma; %MNA
            obj.Pcrit=obj.Pcrit/1000; %kN

        obj.pointer=ones(size(obj.ELEMENTS,1),1);
        if obj.counter==1 %first iteration
            obj.lamda_total(obj.counter)=abs(obj.Pcrit)/max(abs((obj.F_axial))); % find first lamda
            obj.F_axialtotal(:,obj.counter)=obj.F_axial'*obj.lamda_total(obj.counter); % scaling up
            
            elements = size(obj.ELEMENTS,1);
            for i=1:elements
                if abs(obj.F_axialtotal(i,obj.counter))>=abs(obj.Pcrit) %Plastic Squash
                    obj.indicator(i)=0;
                    fprintf('element %g has failed \n',i);
                    obj.pointer(i)=0;
                    obj.mode=1;
                else
                    
                end                
            end
        else % 2nd iteration or above
            % dlamda stores lamda values for each element
            elements = size(obj.ELEMENTS,1);
            for i=1:elements
               if obj.indicator(i)==1 % element still here     
                   % MNA
                   if (obj.Pcrit-obj.F_axialtotal(i,obj.counter-1))/obj.F_axial(i)<=0 % tension but lamda<0
                       %fprintf('dlamda cannot be smaller than 0 \n');
                   elseif (obj.Pcrit-obj.F_axialtotal(i,obj.counter-1))/obj.F_axial(i)>0 % element fail in tension
                       obj.dlamda(i)=(obj.Pcrit-obj.F_axialtotal(i,obj.counter-1))/obj.F_axial(i);
                   end    
                   if (-obj.Pcrit-obj.F_axialtotal(i,obj.counter-1))/obj.F_axial(i)<=0 % compression but lamda<0
                       %fprintf('dlamda cannot be smaller than 0 \n');
                   elseif (-obj.Pcrit-obj.F_axialtotal(i,obj.counter-1))/obj.F_axial(i)>0 % element fail in compression
                       obj.dlamda(i)=(-obj.Pcrit-obj.F_axialtotal(i,obj.counter-1))/obj.F_axial(i);
                   end                   
               else % element failed
                   obj.dlamda(i)=0;
               end
            end
            obj.mode=1; % failure by plastic squashing
            obj.lamda_total(obj.counter)=min(nonzeros(obj.dlamda));
            for i=1:elements  % scaling up the axial forces in remaining elements
                if obj.indicator(i)==1
                    obj.F_axialtotal(i,obj.counter)=obj.F_axial(i)'*obj.lamda_total(end)+obj.F_axialtotal(i,obj.counter-1);
                else
                end
            end
            % check which element failed
            for i=1:elements 
                if abs(obj.F_axialtotal(i,obj.counter))>=abs(obj.Pcrit)
                    obj.indicator(i)=0;
                    obj.pointer(i)=0;
                    fprintf('element %g failed \n',i);
                else
                end
            end

        end
        % calculate total lamda
        obj.lamda=sum(obj.lamda_total(1:end-1));
        end
        %%
        function obj=update(obj)
            % update nodal forces;
            elements = size(obj.ELEMENTS,1);             
            for EL = 1:elements % loop through all elements
                n1 = obj.ELEMENTS(EL,1); n2 = obj.ELEMENTS(EL,2); % identify element node numbers
                x1 = obj.coords(n1,1); y1 = obj.coords(n1,2); % element node 1 - x,y coordinates
                x2 = obj.coords(n2,1); y2 = obj.coords(n2,2); % element node 2 - x,y coordinates
                dof11 = obj.dofs(n1,1); dof12 = obj.dofs(n1,2); % element node 1 - dofs
                dof21 = obj.dofs(n2,1); dof22 = obj.dofs(n2,2); % element node 2 - dofs
                alpha = atan2(y2-y1,x2-x1); 
                Le = sqrt( (x2 - x1)^2 + (y2 - y1)^2 ); % element length
                
                               
                if obj.pointer(EL)==0  % current failed element
                    % GMNA
                    if obj.mode==2 % element failed by buckling
                        Pc=-pi^2*obj.E*(pi/4*obj.R^4)/Le^2;
                        fprintf('element %g failed by buckling',EL);
                    else % element failed by plastic squash
                        Pc=obj.Pcrit;
                    end
                    
                    if alpha>=0 && alpha<=0.5*pi %First Quandrant
                        if obj.F_axialtotal(EL,obj.counter)<0 || obj.mode==2 % Compression or buckling
                            obj.F(dof11)=obj.F(dof11)+(Pc)*abs(cos(alpha));
                            obj.F(dof12)=obj.F(dof12)+(Pc)*abs(sin(alpha));
                            obj.F(dof21)=obj.F(dof21)-(Pc)*abs(cos(alpha));
                            obj.F(dof22)=obj.F(dof22)-(Pc)*abs(sin(alpha));
                        else  % Tension
                            obj.F(dof11)=obj.F(dof11)-Pc*abs(cos(alpha));
                            obj.F(dof12)=obj.F(dof12)-Pc*abs(sin(alpha));
                            obj.F(dof21)=obj.F(dof21)+Pc*abs(cos(alpha));
                            obj.F(dof22)=obj.F(dof22)+Pc*abs(sin(alpha)); 
                        end
                    elseif alpha>0.5*pi && alpha<=pi %Second Quadrant
                        if obj.F_axialtotal(EL,obj.counter)<0 || obj.mode==2  % Compression or buckling
                            obj.F(dof11)=obj.F(dof11)-(Pc)*abs(cos(alpha));
                            obj.F(dof12)=obj.F(dof12)+(Pc)*abs(sin(alpha));
                            obj.F(dof21)=obj.F(dof21)+(Pc)*abs(cos(alpha));
                            obj.F(dof22)=obj.F(dof22)-(Pc)*abs(sin(alpha));
                        else % Tension
                            obj.F(dof11)=obj.F(dof11)+Pc*abs(cos(alpha));
                            obj.F(dof12)=obj.F(dof12)-Pc*abs(sin(alpha));
                            obj.F(dof21)=obj.F(dof21)-Pc*abs(cos(alpha));
                            obj.F(dof22)=obj.F(dof22)+Pc*abs(sin(alpha)); 
                        end      
                    elseif alpha>=-pi && alpha<-0.5*pi %Third Quadrant
                        if obj.F_axialtotal(EL,obj.counter)<0 || obj.mode==2  % Compression or buckling
                            obj.F(dof11)=obj.F(dof11)-(Pc)*abs(cos(alpha));
                            obj.F(dof12)=obj.F(dof12)-(Pc)*abs(sin(alpha));
                            obj.F(dof21)=obj.F(dof21)+(Pc)*abs(cos(alpha));
                            obj.F(dof22)=obj.F(dof22)+(Pc)*abs(sin(alpha));
                        else % Tension
                            obj.F(dof11)=obj.F(dof11)+Pc*abs(cos(alpha));
                            obj.F(dof12)=obj.F(dof12)+Pc*abs(sin(alpha));
                            obj.F(dof21)=obj.F(dof21)-Pc*abs(cos(alpha));
                            obj.F(dof22)=obj.F(dof22)-Pc*abs(sin(alpha)); 
                        end
                    elseif alpha>=-0.5*pi && alpha<0 % Forth Quadrant
                        if obj.F_axialtotal(EL,obj.counter)<0  || obj.mode==2 % Compression or buckling
                            obj.F(dof11)=obj.F(dof11)+(Pc)*abs(cos(alpha));
                            obj.F(dof12)=obj.F(dof12)-(Pc)*abs(sin(alpha));
                            obj.F(dof21)=obj.F(dof21)-(Pc)*abs(cos(alpha));
                            obj.F(dof22)=obj.F(dof22)+(Pc)*abs(sin(alpha));
                        else % Tension
                            obj.F(dof11)=obj.F(dof11)-Pc*abs(cos(alpha));
                            obj.F(dof12)=obj.F(dof12)+Pc*abs(sin(alpha));
                            obj.F(dof21)=obj.F(dof21)+Pc*abs(cos(alpha));
                            obj.F(dof22)=obj.F(dof22)-Pc*abs(sin(alpha)); 
                        end 
                    end                    

%                     fprintf('nodal forces updated for element %g \n',EL);
                else % this element has not failed
                    %fprintf('no nodal forces to update! \n');
                end
                    
                   
            end

        end
            
          
    
    function obj=GMNA(obj) % GMNA
            elements = size(obj.ELEMENTS,1);
            obj.pointer=ones(elements,1);
            new_coords = zeros(size(obj.coords)); 
            amp_coords = zeros(size(obj.coords)); 
            amp = 2; % amplification factor for plotting purposes only 
            for I = 1:size(obj.coords,1)
                for J = 1:size(obj.coords,2)    
                    amp_coords(I,J) = obj.coords(I,J) + obj.U(obj.dofs(I,J))*amp;
                    new_coords(I,J) = obj.coords(I,J) + obj.U(obj.dofs(I,J));
                end
            end

            % Plotting
            figure; hold all; grid on; tol = 1e-3;
            xmin = min(obj.coords(:,1)); xmax = max(obj.coords(:,1)); difx = xmax - xmin;
            ymin = min(obj.coords(:,2)); ymax = max(obj.coords(:,2)); dify = ymax - ymin; fac = 0.25;
            axis([xmin-difx*fac  xmax+difx*fac  ymin-dify*fac  ymax+dify*fac]); 

            % Plastic Squashing
            obj.Pcrit=pi*obj.R^2*obj.sigma; %MNA
            obj.Pcrit=obj.Pcrit/1000; %kN
            if obj.counter==1 % first iteration
                for EL = 1:elements
                    n1 = obj.ELEMENTS(EL,1); n2 = obj.ELEMENTS(EL,2); % identify element node numbers

                    % Plotting original structure
                    x1 = obj.coords(n1,1); y1 = obj.coords(n1,2); % element node 1 - x,y original coordinates
                    x2 = obj.coords(n2,1); y2 = obj.coords(n2,2); % element node 2 - x,y original coordinates
                    Le = sqrt( (x2 - x1)^2 + (y2 - y1)^2 ); % element length

                    
                    % Buckling
                    % different Pcr for each element depending on length
                    obj.Pcr=pi^2*obj.E*(pi/4*obj.R^4)/Le^2;
                    
                   
                    ke = obj.EA/Le; % element axial stiffness
                    alpha = atan2(y2-y1,x2-x1); % angle of inclination relative to the POSITIVE x axis direction
                    plot([x1,x2],[y1,y2],'Color',[0.5 0.5 0.5],'Linewidth',3); 

                    % Check on changes in member lengths and plotting amplified deformed structure
                    x1_amp = amp_coords(n1,1); y1_amp = amp_coords(n1,2); % element node 1 - x,y amplified deformed coordinates
                    x2_amp = amp_coords(n2,1); y2_amp = amp_coords(n2,2); % element node 2 - x,y amplified deformed coordinates     
                    x1_new = new_coords(n1,1); y1_new = new_coords(n1,2); % element node 1 - x,y actual deformed coordinates
                    x2_new = new_coords(n2,1); y2_new = new_coords(n2,2); % element node 2 - x,y actual deformed coordinates    
                    u1 = x1_new - x1; v1 = y1_new - y1; u2 = x2_new - x2; v2 = y2_new - y2; % reconstruction of element global dofs
                    up1 = cos(alpha)*u1 + sin(alpha)*v1; up2 = cos(alpha)*u2 + sin(alpha)*v2; dup = up2 - up1; % reconstruction of element local dofs     
                  
                    % calculate axial force
                    obj.F_axial(EL)=ke*dup;

                    if dup < -tol % element length has decreased - member in compression
                        plot([x1_amp,x2_amp],[y1_amp,y2_amp],'b','Linewidth',3); % blue colour
                    elseif dup > tol % element length as increased - member in tension
                        plot([x1_amp,x2_amp],[y1_amp,y2_amp],'r','Linewidth',3); % red colour
                    else % no change in element length
                        plot([x1_amp,x2_amp],[y1_amp,y2_amp],'k','Linewidth',3); % black colour
                    end
 
                    % Plotting nodes last!
                    plot(x1,y1,'ko','Markersize',7,'MarkerFaceColor','w');  
                    plot(x2,y2,'ko','Markersize',7,'MarkerFaceColor','w');
                    plot(x1_amp,y1_amp,'ko','Markersize',7,'MarkerFaceColor','y');  
                    plot(x2_amp,y2_amp,'ko','Markersize',7,'MarkerFaceColor','y');
                    
                    % find lamda using GMNA
                    lam(EL)=obj.Pcr/obj.F_axial(EL);
                end
                    % scale up axial forces
                    lam=nonzeros(lam);
                    obj.lamda_total(obj.counter)=min(abs(lam));
                    obj.F_axialtotal(:,obj.counter)=obj.F_axial'*obj.lamda_total(obj.counter);
                    % check which element failed
                    for EL=1:elements
                        n1 = obj.ELEMENTS(EL,1); n2 = obj.ELEMENTS(EL,2); % identify element node numbers
                        x1 = obj.coords(n1,1); y1 = obj.coords(n1,2); % element node 1 - x,y original coordinates
                        x2 = obj.coords(n2,1); y2 = obj.coords(n2,2); % element node 2 - x,y original coordinates
                        Le = sqrt( (x2 - x1)^2 + (y2 - y1)^2 ); % element length
                        % check which element failed under buckling
                        obj.Pcr=pi^2*obj.E*(pi/4*obj.R^4)/Le^2;
                        if abs(obj.F_axialtotal(EL,obj.counter))>=obj.Pcr
                            fprintf('element %g has failed \n',EL)
                            obj.indicator(EL)=0;
                            obj.mode=2; % failed by buckling
                            obj.pointer(EL)=0;
                        end    
                    end
            else % second iteration or above
               for EL = 1:elements
                   if obj.indicator(EL)==1 %element still here
                        n1 = obj.ELEMENTS(EL,1); n2 = obj.ELEMENTS(EL,2); % identify element node numbers

                        % Plotting original structure
                        x1 = obj.coords(n1,1); y1 = obj.coords(n1,2); % element node 1 - x,y original coordinates
                        x2 = obj.coords(n2,1); y2 = obj.coords(n2,2); % element node 2 - x,y original coordinates
                        Le = sqrt( (x2 - x1)^2 + (y2 - y1)^2 ); % element length
                        
                        % different values of Pcr for each element
                        % depending on length
                        obj.Pcr=pi^2*obj.E*(pi/4*obj.R^4)/Le^2; 
                        
                        ke = obj.EA/Le; % element axial stiffness
                        alpha = atan2(y2-y1,x2-x1); % angle of inclination relative to the POSITIVE x axis direction
                        plot([x1,x2],[y1,y2],'Color',[0.5 0.5 0.5],'Linewidth',3); 

                        % Check on changes in member lengths and plotting amplified deformed structure
                        x1_amp = amp_coords(n1,1); y1_amp = amp_coords(n1,2); % element node 1 - x,y amplified deformed coordinates
                        x2_amp = amp_coords(n2,1); y2_amp = amp_coords(n2,2); % element node 2 - x,y amplified deformed coordinates     
                        x1_new = new_coords(n1,1); y1_new = new_coords(n1,2); % element node 1 - x,y actual deformed coordinates
                        x2_new = new_coords(n2,1); y2_new = new_coords(n2,2); % element node 2 - x,y actual deformed coordinates    
                        u1 = x1_new - x1; v1 = y1_new - y1; u2 = x2_new - x2; v2 = y2_new - y2; % reconstruction of element global dofs
                        up1 = cos(alpha)*u1 + sin(alpha)*v1; up2 = cos(alpha)*u2 + sin(alpha)*v2; dup = up2 - up1; % reconstruction of element local dofs     
                        
                        % calculate axial force 
                        obj.F_axial(EL)=ke*dup;
                        
                        if dup < -tol % element length has decreased - member in compression
                            plot([x1_amp,x2_amp],[y1_amp,y2_amp],'b','Linewidth',3); % blue colour
                        elseif dup > tol % element length as increased - member in tension
                            plot([x1_amp,x2_amp],[y1_amp,y2_amp],'r','Linewidth',3); % red colour
                        else % no change in element length
                            plot([x1_amp,x2_amp],[y1_amp,y2_amp],'k','Linewidth',3); % black colour
                        end

                        % Plotting nodes last!
                        plot(x1,y1,'ko','Markersize',7,'MarkerFaceColor','w');  
                        plot(x2,y2,'ko','Markersize',7,'MarkerFaceColor','w');
                        plot(x1_amp,y1_amp,'ko','Markersize',7,'MarkerFaceColor','y');  
                        plot(x2_amp,y2_amp,'ko','Markersize',7,'MarkerFaceColor','y');
                        
                        % GMNA, element can fail under plastic yielding and
                        % buckling
                        % check for plastic yielding
                        if (obj.Pcrit-obj.F_axialtotal(EL,obj.counter-1))/obj.F_axial(EL)<=0 % tension but lamda<0
                            %fprintf('dlamda cannot be smaller than 0 \n');
                        elseif (obj.Pcrit-obj.F_axialtotal(EL,obj.counter-1))/obj.F_axial(EL)>0 % element fail in tension
                            obj.dlamda(EL)=(obj.Pcrit-obj.F_axialtotal(EL,obj.counter-1))/obj.F_axial(EL);
                        end    
                        if (-obj.Pcrit-obj.F_axialtotal(EL,obj.counter-1))/obj.F_axial(EL)<=0 % compression but lamda<0
                            %fprintf('dlamda cannot be smaller than 0 \n');
                        elseif (-obj.Pcrit-obj.F_axialtotal(EL,obj.counter-1))/obj.F_axial(EL)>0 % element fail in compression
                            obj.dlamda(EL)=(-obj.Pcrit-obj.F_axialtotal(EL,obj.counter-1))/obj.F_axial(EL);
                        end 
                        % check for buckling
                        if (-obj.Pcr-obj.F_axialtotal(EL,obj.counter-1))/obj.F_axial(EL)<=0 % buckling but lamda<0
                           %fprintf('dlamda cannot be smaller than 0 \n');
                        elseif (-obj.Pcr-obj.F_axialtotal(EL,obj.counter-1))/obj.F_axial(EL)>0
                            obj.bucklinglamda(EL)=(-obj.Pcr-obj.F_axialtotal(EL,obj.counter-1))/obj.F_axial(EL);

                            if obj.bucklinglamda(EL)<obj.dlamda(EL)% choose the smallest lamda values, comparing buckling lamda with plastic lamda
                               obj.dlamda(EL)=obj.bucklinglamda(EL);                              
                            else
                            end  
                        end
                   else
                        obj.F_axial(EL)=0; % failed element
                        n1 = obj.ELEMENTS(EL,1); n2 = obj.ELEMENTS(EL,2); % identify element node numbers

                        % Plotting original structure
                        x1 = obj.coords(n1,1); y1 = obj.coords(n1,2); % element node 1 - x,y original coordinates
                        x2 = obj.coords(n2,1); y2 = obj.coords(n2,2); % element node 2 - x,y original coordinates
                        plot(x1,y1,'ko','Markersize',7,'MarkerFaceColor','w');  
                        
                        x1_amp = amp_coords(n1,1); y1_amp = amp_coords(n1,2); % element node 1 - x,y amplified deformed coordinates
                        x2_amp = amp_coords(n2,1); y2_amp = amp_coords(n2,2); % element node 2 - x,y amplified deformed coordinates                        
                        % Plot nodes
                        plot(x2,y2,'ko','Markersize',7,'MarkerFaceColor','w');
                        plot(x1_amp,y1_amp,'ko','Markersize',7,'MarkerFaceColor','y');  
                        plot(x2_amp,y2_amp,'ko','Markersize',7,'MarkerFaceColor','y');
                        
%                         fprintf('no plot since element %g has failed \n',EL);
                        obj.dlamda(EL)=10000; %place holder otherwise the value 0 will be chosen for lamda
                        obj.bucklinglamda(EL)=10000;
                   end
               end
                    % scale up axial forces
                    obj.lamda_total(obj.counter)=min(abs(obj.dlamda)); % smallest lamda across elements
                    % calculate total lamda
                    obj.lamda=sum(obj.lamda_total(1:end-1));
                    % check which element failed
                    for EL=1:elements
                        n1 = obj.ELEMENTS(EL,1); n2 = obj.ELEMENTS(EL,2); % identify element node numbers

                        % Plotting original structure
                        x1 = obj.coords(n1,1); y1 = obj.coords(n1,2); % element node 1 - x,y original coordinates
                        x2 = obj.coords(n2,1); y2 = obj.coords(n2,2); % element node 2 - x,y original coordinates
                        Le = sqrt( (x2 - x1)^2 + (y2 - y1)^2 ); % element length
                        
                        % different Pcr depending on member length
                        obj.Pcr=pi^2*obj.E*(pi/4*obj.R^4)/Le^2; 
                        
                        if obj.indicator(EL)==1 % element still here
                            obj.F_axialtotal(EL,obj.counter)=obj.F_axial(EL)'*obj.lamda_total(end)+obj.F_axialtotal(EL,obj.counter-1);
                        else % element failed, no need to scale the axial force (equals to 0)
                        end
                        if abs(obj.F_axialtotal(EL,obj.counter))==obj.Pcr % buckling
                            fprintf('element %g has failed \n',EL)
                            obj.indicator(EL)=0;
                            obj.mode=2;
                            obj.pointer(EL)=0;
                            
                        elseif abs(obj.F_axialtotal(EL,obj.counter))>=obj.Pcrit % plastic yielding
                            fprintf('element %g has failed \n',EL)
                            obj.indicator(EL)=0;
                            obj.mode=1;
                            obj.pointer(EL)=0;
                            
                        end                       
                    end    
            end
            xlabel('x coordinate','FontSize',16);
            ylabel('y coordinate','FontSize',16);


            % Printing computed dofs & reactions 
            for dof = 1:length(obj.uF)
%                 disp(['The value of dof ',num2str(obj.dofs_free(dof)),' is ',num2str(obj.uF(dof))]);    
            end
            disp(' ');
            for react = 1:length(obj.fR)
%                 disp(['The value of the reaction at dof ',num2str(obj.dofs_restrained(react)),' is ',num2str(obj.fR(react))]);
            end
            disp(' '); disp('Vertical equilibrium check:');
            disp(['Total vertical reactions = ',num2str(obj.fR(2) + obj.fR(4)+obj.fR(6)+obj.fR(8))]);
            disp(['Total applied vertical loads = ',num2str(-5.5*obj.P*obj.lamda_total)]);
            if abs(obj.fR(2) + obj.fR(4) - obj.P) < 1e-6; disp('Ok.'); end
            disp(' '); disp('Horizontal equilibrium check:');
            disp(['Total horizontal reactions = ',num2str(obj.fR(1) + obj.fR(3)+obj.fR(5)+obj.fR(7))]);
            disp('Total applied horizontal loads = 0');
            if abs(obj.fR(1) + obj.fR(3)+obj.fR(5)+obj.fR(7)) < 1e-6; disp('Ok.'); end
        end



        end
end
    
        
        




      
        

