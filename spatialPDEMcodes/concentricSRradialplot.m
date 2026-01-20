[mySurface,parameters,types] = runTheModel();

%%

figure(1)
subplot(1,2,1)

detect = 0.0;

S = types.S;
S = S.*(S > detect);
R = types.R;
R = R.*(R > detect);

radialPlot(1,S,'g')
hold on
radialPlot(1,R,'r')
lighting GOURAUD
legend('S','R')
view([-60 44])
%zlim([detect 0.6])
axis tight

subplot(1,2,2)
radialPlot(1,S,'g')
hold on
radialPlot(1,R,'r')
lighting flat
zlim([detect 0.6])
view(2)

set(1,'position',[53         320        1387         509])
export_fig(['figures/RFP-GFP-radialPlot-',date,'.png'])

