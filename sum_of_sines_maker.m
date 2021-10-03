function sim_eqn = sum_of_sines_maker(coefficients,project_file,dType)
    % project_file is a boolean that determines whether or not the output function is for a sim or project file
    % The equations for these two different file types are quite different
    % dType: char: decomposition type ('sumSines','fourier')
    
    if nargin < 3
        dType = 'sumSines';
    end
    
    %Animatlab doesn't play well with small exponential numbers, this is pre-processing to remove especially small numbers
    coefficients(log10(cell2mat(coefficients(:,2)))<-4,2)={0};
    
    if project_file && strcmp(dType,'sumSines')
        % This section allows the user to create a sum of sines equation in the format used by the .aproj files from Animatlab
        % For automated simulations, most injections occur in the .asim file, instead, which uses a different format
        %syms a b c t
%         eqn = sym(zeros(length(coefficients)/3,1));
%         count = 0;
        sim_eqn = [];
            for i=1:3:length(coefficients)
%                 count = count + 1;
                a = coefficients{i,2};
                b = coefficients{i+1,2};
                c = coefficients{i+2,2};
%                 eqn(count,1) = vpa(a*sin(b*t+c),5);
                if a >= 0 && i ~=1
                    astring = ['+',num2str(a)];
                else
                    if a < 0
                        error('Animatlab can''t handle negative ''a'' values for some reason. Change cftool bounds on ''a'' and run again.')
                    else
                        astring = num2str(a);
                    end
                end
                if c > 0
                    cstring = ['+',num2str(c)];
                else
                    cstring = num2str(c);
                end
                sim_eqn = [sim_eqn,astring,'*sin(',num2str(b),'*t',cstring,')'];
            end
%             sim_eqn = vpa(sum(eqn(1:end)),5);
    elseif project_file && strcmp(dType,'fourier')
        % This section allows the user to create a sum of sines equation in the format used by the .aproj files from Animatlab
        % For automated simulations, most injections occur in the .asim file, instead, which uses a different format
        %syms a b c t
%         eqn = sym(zeros(length(coefficients)/3,1));
%         count = 0;
        sim_eqn = [];
        w = num2str(coefficients{end,2});
        if w < 0
            error('w is negative, this may cause problems')
        end
        a0 = coefficients{1,2};
        sim_eqn = num2str(a0);
        count = 1;
            for i=2:2:length(coefficients)-2
                a = coefficients{i,2};
                b = coefficients{i+1,2};
                if a >= 0
                    astring = ['+',num2str(a)];
                else
                    astring = num2str(a);
                end
                if b >= 0
                    bstring = ['+',num2str(b)];
                else
                    bstring = num2str(b);
                end
                sim_eqn = [sim_eqn,astring,'*cos(',num2str(count),'*t*',w,')',bstring,'*sin(',num2str(count),'*t*',w,')'];
                count = count + 1;
            end
    else
        % This section creates sum of sines equations formatted for .asim files
        % These equations are more compact than the .aproj counterparts and store sine algebra in a "tail" refered here as the plus_trail
        sim_eqn = [];
        plus_trail = [];
        
        for j=1:3:length(coefficients)
            a = coefficients{j,2};
            b = coefficients{j+1,2};
            c = coefficients{j+2,2};
            ignoreplustrail = 0;
            if j == 1
                if a < 0
                    astring = ['0,',strrep(num2str(abs(a)), '0.', '.'),','];
                else
                    astring = [strrep(num2str(abs(a)), '0.', '.'),','];
                    ignoreplustrail = 1;
                end
            else
                if a < 0
                    astring = ['0,*,',strrep(num2str(abs(a)), '0.', '.'),','];
                else
                    astring = ['*,',strrep(num2str(abs(a)), '0.', '.'),','];
                end
            end
            if b < 0
                bstring = ['0,',num2str(abs(b)),','];
                sinstring = '-,sin,';
            else
                bstring = [num2str(b),','];
                sinstring = 'sin,';
            end
            if c < 0
                cstring = [num2str(abs(c)),',-,'];
            else
                cstring = [num2str(c),',+,'];
            end
            sinString = [astring,bstring,'t,*,',cstring,sinstring];
            sim_eqn = [sim_eqn,sinString];
             if ~ignoreplustrail
                if a > 0
                    plus_trail = [plus_trail,',+'];
                else
                    plus_trail = [plus_trail,',-'];
                end
             end
        end
        if length(coefficients) == 3 && a > 0
            sim_eqn = [sim_eqn,'*',flip(plus_trail(2:end))];
        else
            sim_eqn = [sim_eqn,'*,',flip(plus_trail(2:end))];   
        end
    end   
end