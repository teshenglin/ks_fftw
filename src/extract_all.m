% %%% Extract branches from b.all:
% fid = fopen('b.tt1');
% 
% tline = fgets(fid);
% while ischar(tline)
%     while ischar(tline) & str2num(tline(1:4))==0
%         tline = fgets(fid);
%     end
%     BR=str2num(tline(1:4))
%     filename=['BR_',num2str(BR),'.dat'];
%     BR_file_ID = fopen(filename,'w');
%     while ischar(tline) & str2num(tline(1:4))~=0
%         fprintf(BR_file_ID,'%s',tline);
%         % add the line to the file
%         tline = fgets(fid);
%     end
%     fclose(BR_file_ID);
%     % save the file
%     % disp(tline);
%     % tline = fgets(fid);
% end
% 
% fclose(fid);

%%% Extract LABs from b.all:
fid = fopen('b.tt1');
filename=['LABS.dat'];
LABs_file_ID = fopen(filename,'w');
tline = fgets(fid);
while ischar(tline)
    while ischar(tline) & str2num(tline(1:4))==0
        tline = fgets(fid);
    end
    while ischar(tline) & str2num(tline(1:4))~=0
        if str2num(tline(15:19))~=0
            fprintf(LABs_file_ID,'%s',tline);
            % add the line to the file
        end
        tline = fgets(fid);
    end
    % disp(tline);
    % tline = fgets(fid);
end
fclose(LABs_file_ID);
fclose(fid);

%%% Extract solutions from s.all
fid = fopen('s.tt1');

tline = fgets(fid);
while ischar(tline)
    if tline(6)~=' '
        LAB=str2num(tline(19:24))
        filename=['sol_LAB_',num2str(LAB),'.dat'];
        LAB_file_ID = fopen(filename,'w');
        while ischar(tline)
            tline = fgets(fid);
            if tline(5)~=' '
                break
            end
            fprintf(LAB_file_ID,'%s',tline(1:(end-1)));
            % add the line to the file
        end
        fclose(LAB_file_ID);
        % save the file
    end
%    disp(tline);
    tline = fgets(fid);
end

fclose(fid);