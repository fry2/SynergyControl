% Reads in motion data from Data/Abnormal Gait Processing and plots it over time for the sagittal and coronal views

% Run abnormalGaitProcessing first
abnormalGaitProcessing

dataContent = table2array(tbl(1:end,3:end));
startInd = find(sum(dataContent,2),1,'first');
endInd = find(sum(dataContent,2),1,'last');

xvals = dataContent(startInd:endInd,1:3:end);
yvals = dataContent(startInd:endInd,2:3:end);
zvals = dataContent(startInd:endInd,3:3:end);

xbounds = [min(nonzeros(xvals)) max(nonzeros(xvals))]+[-.5 .5];
ybounds = [min(nonzeros(yvals)) max(nonzeros(yvals))]+[-.5 .5];
zbounds = [min(nonzeros(zvals)) max(nonzeros(zvals))]+[-.5 .5];

figHandles = get(groot, 'Children');
if ~isempty(figHandles)
    priorFigures = contains({figHandles(:).Name},{'WalkFig'});
    close(figHandles(priorFigures))
end

fig = figure('name','WalkFig');
H = uicontrol('Style', 'PushButton', ...
              'String', 'Break', ...
              'Callback', 'delete(gcbf)');

a =diff(xbounds);
b =diff(ybounds);
c =diff(zbounds);
axeslengths = [a b c];
normedaxes = axeslengths/norm(axeslengths);
ii = 301;

while (ishandle(H)) || 1
        femur = [HipX(ii),HipY(ii),HipZ(ii);KneeX(ii),KneeY(ii),KneeZ(ii)];
        tibia = [KneeX(ii),KneeY(ii),KneeZ(ii);AnkleX(ii),AnkleY(ii),AnkleZ(ii)];
        subplot(2,1,1)
            plot3(MPelantemed1X(ii),MPelantemed1Y(ii),MPelantemed1Z(ii),'mo')
            hold on
            plot3(MPelantelat1X(ii),MPelantelat1Y(ii),MPelantelat1Z(ii),'mx')
            plot3(MPelpostmed1X(ii),MPelpostmed1Y(ii),MPelpostmed1Z(ii),'mo')
            plot3(MPelpostlat1X(ii),MPelpostlat1Y(ii),MPelpostlat1Z(ii),'mx')
            if ~any(sum(femur,2)==0) || ~any(sum(tibia,2)==0)
                plot3(femur(:,1),femur(:,2),femur(:,3),'k--','LineWidth',2)
                plot3(tibia(:,1),tibia(:,2),tibia(:,3),'k--','LineWidth',2)
            end
            plot3(HipX(ii),HipY(ii),HipZ(ii),'r*')
            plot3(MFemante1X(ii),MFemante1Y(ii),MFemante1Z(ii),'mo')
            plot3(MFempost1X(ii),MFempost1Y(ii),MFempost1Z(ii),'mx')
            plot3(KneeX(ii),KneeY(ii),KneeZ(ii),'g*')
            plot3(MTibante1X(ii),MTibante1Y(ii),MTibante1Z(ii),'mo')
            plot3(MTibpost1X(ii),MTibpost1Y(ii),MTibpost1Z(ii),'mx')
            plot3(AnkleX(ii),AnkleY(ii),AnkleZ(ii),'b*')
            grid on
            xlabel('X')
            ylabel('Y')
            zlabel('Z')
            xlim(xbounds)
            ylim(ybounds)
            zlim(zbounds)
            view([-180 -90])
    %         view([-180 (-90/130)*abs(ii-550)])
            pbaspect(normedaxes)
            hold off
        subplot(2,1,2)
            plot3(MPelantemed1X(ii),MPelantemed1Y(ii),MPelantemed1Z(ii),'mo')
            hold on
            plot3(MPelantelat1X(ii),MPelantelat1Y(ii),MPelantelat1Z(ii),'mx')
            plot3(MPelpostmed1X(ii),MPelpostmed1Y(ii),MPelpostmed1Z(ii),'mo')
            plot3(MPelpostlat1X(ii),MPelpostlat1Y(ii),MPelpostlat1Z(ii),'mx')
            if ~any(sum(femur,2)==0) || ~any(sum(tibia,2)==0)
                plot3(femur(:,1),femur(:,2),femur(:,3),'k--','LineWidth',2)
                plot3(tibia(:,1),tibia(:,2),tibia(:,3),'k--','LineWidth',2)
            end
            plot3(HipX(ii),HipY(ii),HipZ(ii),'r*')
            plot3(MFemante1X(ii),MFemante1Y(ii),MFemante1Z(ii),'mo')
            plot3(MFempost1X(ii),MFempost1Y(ii),MFempost1Z(ii),'mx')
            plot3(KneeX(ii),KneeY(ii),KneeZ(ii),'g*')
            plot3(MTibante1X(ii),MTibante1Y(ii),MTibante1Z(ii),'mo')
            plot3(MTibpost1X(ii),MTibpost1Y(ii),MTibpost1Z(ii),'mx')
            plot3(AnkleX(ii),AnkleY(ii),AnkleZ(ii),'b*')
            grid on
            xlabel('X')
            ylabel('Y')
            zlabel('Z')
            xlim(xbounds)
            ylim(ybounds)
            zlim(zbounds)
            view([-90 0])
            camroll(-90)
    %         view([-180 (-90/130)*abs(ii-550)])
            pbaspect(normedaxes)
            hold off
        drawnow
        pause(.0017)
        if ii == endInd
            ii = startInd;
        else
            ii = ii + 1;
        end
end