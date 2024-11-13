% function musclepath_poly(muscles,mydir,musclepolyfile)
clearvars
clear all
clc
osimfile = 'das3_mod_full.osim';
model = das3_readosim(osimfile);
mydir = 'polys_from_abd';
muscles = model.muscles;


for imus = 1:length(muscles)
    
    for EULorQ = 1:2

    musfilename = [mydir,'\path_',muscles{imus}.name,'.mat'];
    ma = load(musfilename);
    mus = muscles{imus};
    

    ndofs = length(mus.dof_indeces); % number of dofs spanned by this muscle
    order = 4; % polynomial order
    
    if EULorQ == 1
        jnts = ma.alljnts;
        jacobs = ma.allmomarms;
        current_poly = 'Euler';
    elseif EULorQ == 2
        jnts = ma.alljntsQA;
        jacobs = ma.quat_J;
        current_poly = 'Quaternion';
    end

    num_data = size(ma.alljnts,1);
    fprintf(1,'Muscle name:      %s\n',mus.name);
    %
    % count how many parameters the polynomial model for muscle length has
    npar = prod(1:ndofs+order)/prod(1:ndofs)/prod(1:order);
    % fprintf(1,'Number of DOFs:   %d\n',ndofs);
    % fprintf(1,'Polynomial order: %d\n',order);
    % fprintf(1,'Potential number of polynomial terms: %d\n',npar);
    tot_data = num_data*(ndofs+1);	% total number of data points
    A = zeros(tot_data, npar);      % allocate memory space for A
    b = zeros(tot_data, 1);         % allocate memory space for b
    
    % get angle values for the dofs the muscle crosses
    musdof_indeces = zeros(ndofs,1);
    for idof = 1:ndofs
        imusdof = mus.dof_indeces(idof);
        musdof_indeces(idof) = imusdof;
    end
    
    ang = (jnts(:,musdof_indeces) + 1e-6);	% protect against angle = 0.0
    
    maxmomdof = zeros(1,ndofs);
    for idof = 1:ndofs
        maxmomdof(idof) = max(abs(jacobs(:,idof)))*1000;
    end    
    maxall = max(maxmomdof);
    % this normalises all moment arms
    ml_weight = (maxmomdof)/maxall;
    % Stopping criterion: error less than 10% of maximum moment arm (in mm) 
    % for the muscle or 2mm, whichever is greater
    momarm_error = max(0.1*maxall,2); 
    % momarm_error = 0.01*maxall; 
    %
    for idof = 1:ndofs
        % read moment arm from allmomarms matrix
        % and angles from alljnts matrix
        b((idof-1)*num_data+1:idof*num_data) = -jacobs(:,idof)/ml_weight(idof); %
    
        % generate the npar polynomial terms, and for each term, add a column to A
        polylist = zeros(npar,ndofs);
        expon = zeros(num_data,ndofs);	% start with all exponents zero
        for ii=1:npar
            polylist(ii,:) = expon(1,:);
            A((idof-1)*num_data+1:idof*num_data,ii) = (expon(:,idof).*prod(ang.^expon,2)./ang(:,idof))/ml_weight(idof); % % contribution of this term to moment arm idof
            % generate the next set of exponents, but remain within model order
            k = 1;
            while (1)
                expon(:,k) = expon(:,k)+1;
                if (sum(expon(1,:)) > order && ii<npar)
                    expon(:,k)=0;
                    k = k+1;
                else
                    break;
                end
            end
        end     % done generating model terms
    end		% done reading all data for this muscle
    % <num_data> more rows for muscle length
    % read length from alllengths vector
    
    % and angles from alljnts matrix
    
    b(ndofs*num_data+1:(ndofs+1)*num_data) = ma.alllengths;
    % generate the npar polynomial terms, and for each term, add a column to A    
    %
    
    for ipar=1:npar
        help = repmat(polylist(ipar,:),size(ang,1),1);
        A(ndofs*num_data+1:(ndofs+1)*num_data,ipar) = (prod(ang.^help,2));
    end 
    
    % fprintf('Total number of data points: %d\n',tot_data);
    
    % now we have all data for this muscle stored in A and b
    % A = A(1:3174,end-11);
    % b = b(1:3174);
    % A = A(3175:end,1:10);
    % b = b(3175:end);
    
    % solve the full model with npar-1 terms
    p = A\b;		% compute coefficients of the best fitting model
    bpred = A*p;	% these are the moment arms predicted by this model
    res = (bpred-b);	% residuals
    resJac = res(ndofs*num_data);
    resLength = res((ndofs+1)*num_data);
    RMSfull = (sqrt(sum(res.^2)/tot_data))*1000;
    RMSfullLength = (sqrt(sum(resLength.^2)/length(resLength)))*1000;
    RMSfullJac = (sqrt(sum(resJac.^2)/length(resJac)))*1000;
    disp(['RMS fit error of the full model: ',num2str(RMSfull)]);
    disp(['maximum moment arm: ', num2str(maxall)]);
    disp(['Max error in res: ',num2str(max(abs(res)))])
    disp(['Sum error of J:' num2str(sum(abs(res(1:ndofs*num_data))))])
    disp(['Sum error of length: ', num2str(sum(abs(res(ndofs*num_data+1:(ndofs+1)*num_data))))])
    disp('---------------------------')

    mus_model{imus}.name = mus.name;
    mus_model{imus}.num_lparams = npar;
    mus_model{imus}.lparams = polylist; % size of lparams: num_lparams x num_dofs
    mus_model{imus}.lcoef = p;
    npar
    RMSbar(imus,EULorQ) = RMSfull;
    RMSbarLength(imus,EULorQ) = RMSfullLength;
    RMSbarJac(imus,EULorQ) = RMSfullJac;
    xbar{imus} = mus.name;
    
    % if ~mod(imus,1)
    % examine_momarms(mus_model{imus}, mus.dof_names, jacobs, ang);	
    % end

    clear ma A b ang res resJac resLength;
    end
end

%%
Sxbar = categorical(string(xbar));
start_ind = 20;
end_ind = 60;
space = 4;
Sxbar = Sxbar(start_ind:space:end_ind);
Sxbar = removecats(Sxbar);
% Sxbar = removecats(Sxbar);
% figure
% bar(Sxbar,RMSbar)
% title('RMS jacobian + lengths')
% legend('Eul','Quat')
% set(gca, 'YScale', 'log')
figure
tiledlayout(1,2)
fig = gcf;
nexttile
bar(Sxbar,RMSbarLength(start_ind:space:end_ind,:))
title({'RMS of lengths approximations','w.r.t OpenSim values'})
legend('Eul','Quat','Location','northwest')
ylabel('RMS in log scale [mm]','FontWeight','bold')
set(gca, 'YScale', 'log')

nexttile
bar(Sxbar,RMSbarJac(start_ind:space:end_ind,:))
title({'RMS of Jacobians approximations','w.r.t OpenSim values'})
legend('Eul','Quat','Location','northwest')
ylabel('RMS in log scale [mm]','FontWeight','bold')
set(gca, 'YScale', 'log')

sppi = get(fig,"ScreenPixelsPerInch");
exportgraphics(fig,'full_RMS.png','Resolution',1000);



function examine_momarms(musmodel, dof_names, moment_arms, sangles)
% plot momentarm-angle data

% choose a subset of "angles" that contains only 100 points
% a = 1;
% b = size(angles,1);
% r = a + (b-a).*rand(100,1);
% indeces = ceil(sort(r));
% sangles = angles(indeces,:);

% ...or use all angles
%indeces = 1:size(angles,1);
%sangles = angles;

% calculate moment arms from polynomial
pmoment_arms = zeros(length(sangles(:,1)),length(dof_names));
for iframe = 1:length(sangles(:,1))
    for i=1:musmodel.num_lparams

        % add this term's contribution to the muscle length 
        term = musmodel.lcoef(i);

        for j=1:length(dof_names)
            for k=1:musmodel.lparams(i,j)
                term = term * sangles(iframe,j); % this creates lcoeff(i) * product of all angles to the power lparams(i,j) 
            end
        end

        % first derivatives of length with respect to all q's
        for  k=1:length(dof_names)
            % derivative with respect to q_k is zero unless exponent is 1 or higher and q is not zero
            if ((musmodel.lparams(i,k) > 0) && (sangles(iframe,k)))	
                dterm = musmodel.lparams(i,k)*term/sangles(iframe,k);
                pmoment_arms(iframe,k) = pmoment_arms(iframe,k) + dterm;
            end
        end
    end
end

figure;
for idof=1:length(dof_names)
    subplot(length(dof_names),1,idof);
    plot(moment_arms(:,idof),'bx-'); hold on; plot(-pmoment_arms(:,idof),'ro-'); 	
    title([dof_names{idof}, ' momentarms for ',musmodel.name],'Interpreter', 'none'); 
end
legend('osim','poly');

end