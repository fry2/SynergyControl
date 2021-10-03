
sagvec = -obj.joint_obj{1}.uuw_joint(:,1)';

m = length(obj.joint_obj{1}.uuw_joint);
n = length(forces);
if m ~= n
    forcesBig = interp1(1:n,forces,linspace(1,n,m));
end

for ii = 1:38
    muscle = obj.musc_obj{ii};
    whichsegmentisfree = diff(cell2mat(muscle.pos_attachments(:,3)));
    seg2use = find(whichsegmentisfree,1,'first'); % This is not a stringent search function <<< edit later
    matt1 = muscle.pos_attachments{seg2use,4};
    matt2 = muscle.pos_attachments{seg2use+1,4};
    %tdeg = zeros(length(matt1),1);
    for jj = 1:length(matt1)
        mvec = matt2(jj,:)-matt1(jj,:);
        theta(jj,1) = atan2d(norm(cross(sagvec,mvec)),dot(sagvec,mvec));
        tdeg(jj,ii) = sin(atan2d(norm(cross(sagvec,mvec)),dot(sagvec,mvec))*(pi/180));
    end
    outforces(:,ii) = forcesBig(:,ii)./tdeg(:,ii);
end

figure
subplot(3,1,1)
plot(forces)
subplot(3,1,2)
plot(outforces)
subplot(3,1,3)
plot(outforces-forcesBig)