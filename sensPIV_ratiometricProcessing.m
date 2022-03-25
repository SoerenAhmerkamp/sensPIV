clear oxynano
%%

%%
files = dir('*.bmp');
imgT = zeros(2048,2448,1);
saveR = zeros(2048,2448,1);
saveG = zeros(2048,2448,1);
c = 0;

extr = [];
extrStd = [];
extrMax = [];
extrMin = [];

for m = 10:100%length(files)
    c = c+1;
    img = imread(files(m).name);
    img = im2double(img);
    R = img(:,:,1);
    G = img(:,:,2);
    B = img(:,:,3);
    imgT = imgT+(R)./(G);
    saveR = saveR+R;
    saveG = saveG+G;

    rect = medfilt2(imgT(1000:1050,1000:1050),[4 4]);
    %nanmean(rect(:))
    extr(end+1) = nanmedian(rect(:))/c;
    extrStd(end+1) = nanstd(rect(:)/c);
    extrMax(end+1) = quantile(rect(:)/c,0.65);
    extrMin(end+1) = quantile(rect(:)/c,0.35);
     
    m
end
%
% hold on
% fill([1:91 92-(1:91)]-1,[extrMax fliplr(extrMin)],[0.8 0.8 0.8])
% plot(extr,'--')
% %hline(extr)
% xlim([0 91])
% figure
% fill([1:91 92-(1:91)]-1,[extrMax./extr fliplr(extrMin./extr)],[0.8 0.8 0.8])
% hline(1)
% xlabel('Imagenumber')
% figure
%
img = imgT / c;
saveR = saveR / c;
saveG = saveG / c;

%
subplot(1,2,1)
hold on
imagesc(oxyCalib(flipud(medfilt2(img,[4,4]))))
%caxis([0.3 0.6]), axis equal
%caxis([0.2 0.5]), axis equal
oxynano.oxygen.oxyRaw = medfilt2(saveR./saveG,[4,4]);

%% Requires oxygen calibration function oxycalib
[Fx,Fy]=gradient(oxynano.oxygen.oxyRaw,100);
oxynano.oxygen.dCdx = Fx;
oxynano.oxygen.dCdy = Fy;

colormap jet, axis tight, axis equal
%caxis([550 650]), colorbar
%imagesc(oxynano.oxygen.dCdy)

subplot(1,2,2)
plot(oxyCalib(nanmean(img(:,1550:1600),2)),'x')

%% Convert to BW
mkdir BW
files = dir('*.bmp');

img0 = imread(files(1).name);
img0 = im2double(img0);
img = img0;
o = 0; 
for n = 10:length(files) %length(files)
    o = o+1;
    imgT = imread(files(n).name);
    imgT = im2double(imgT); 
    img = imgT + img;
    
    highpsize = 25;
   % h = fspecial('gaussian',highpsize,highpsize);
    %in_roi=double(in_roi-(imfilter(in_roi,h,'replicate')));
	%imgT=(imgT(:,:,2)-(imfilter(imgT(:,:,2),h,'replicate')));
    %imwrite(imgT,sprintf('./BW/%04g.bmp',n))    
    imwrite(imgT(:,:,2),sprintf('./BW/%04g.bmp',n))
end
imwrite(img(:,:,2)/o,sprintf('./Background.bmp'))
imagesc(img(:,:,2)/o)

%% Coral Masking and Surface Detection
figure
cd BW
imgRef = im2double(imread('../Background.bmp'));
oxynano.imgRef = imgRef;
% this might need modification
%imgRef(1800:end,:) = 1;
%imgRef(1600:end,1200:end) = 1;
%imgRef(1:1200,:) = 0;

imgBW = imbinarize(imgRef,0.55);
mask = ~imgBW;
%imagesc(imgBW)
trace =  bwconncomp(imgBW);

st = regionprops(trace,  'PixelList');
xy = cat(1,st(1).PixelList);
xy(:,2) = abs(xy(:,2) - size(trace,1)) + 1;
x = xy(:,1);
y = xy(:,2);
%plot(xy(:,1),xy(:,2),'g.')

xS = unique(x);
yS = zeros(1,length(xS));

for n = 1:length(yS)
    yS(n) = min(y(x == xS(n))); 
end

oxynano.mask = mask;
mask = double(mask);
mask(mask == 0) = NaN;
oxynano.maskN = mask;
oxynano.surface2.x = xS;
oxynano.surface2.y = yS';

figure
hold on
imagesc(imgRef),shading flat
plot(xS,yS','r','LineWidth',2), axis tight
cd ..

%% Pathlines

cd BW
files = dir('*.bmp');
figure
I1_stream = zeros(size(imgRef));
counter = zeros(size(imgRef));

imgRef = imread('../Background.bmp');
imgRef = im2double(imgRef);

ths = 0.1;
for n = 10:length(files)
   n
   
   imgT = imread(files(n).name);
   I1 = (im2double(imgT) - imgRef)./imgRef;
   %I1 = (im2double(imgT));
   I1(I1<ths) = 0;
  % imagesc(I1)

   %imshow(I1), pause(0.1)   
   counter(I1>ths) = counter(I1>ths) + 1;
   I1_stream = I1_stream+I1;
   pause(0.1)
   
end
cd ..
%
figure
imagesc(medfilt2(I1_stream./counter,[4,4])), colormap gray, axis equal, axis tight,% caxis([0, 0.01])
axis equal
oxynano.streamlines = medfilt2(I1_stream./counter,[4,4]);
%ylim([2048-1300,2048-400])
%xlim([300, 2400])
%% Particle Tracking
% Uses lagrangian particle tracker
%cd D:\Temp\CoralFlow\Backup
figure
cd BW
[vtracks,ntracks] = PredictiveTracker('./*.bmp',10,500,'../Background.bmp',3,0,0);
%[vtracks,ntracks] = PredictiveTracker('./*.bmp',10,100,'../Background.bmp',15,0,0);

%
hold on
im = imread('../background.bmp');
imagesc(im), shading flat, colormap gray

for n = 1:length(vtracks)
    plot(vtracks(n).X,vtracks(n).Y,'-','Linewidth',2)
end

axis tight
cd ..
%%
%hold on

[u,v,x,y,t,tr] = Velocities(vtracks([vtracks(:).len] > 10),[1 200],1);

scale = v.^2*2;
scale2 = 10./abs(u);
scale2(scale2>10000) = 1;
%scale(u>10) = NaN;
%scale(scale>100) = NaN;
%scale(y<1000) = NaN;
%

oxynano.tracks.u = u;
oxynano.tracks.v = v;
oxynano.tracks.mag = sqrt(u.^2+v.^2);
oxynano.tracks.x = x;
oxynano.tracks.y = y;

figure
scatter(oxynano.tracks.x,-oxynano.tracks.y,[],oxynano.tracks.v,'.')


%% Smooth Scatter data
oxynano.tracks.vSmooth = scatstat2(x,y,v,25,@median); 
oxynano.tracks.uSmooth = scatstat2(x,y,u,25,@median); 


%% Grid Velocity

u = oxynano.tracks.uSmooth;
v = oxynano.tracks.vSmooth;
x = oxynano.tracks.x;
y = oxynano.tracks.y;

[xGr, yGr] = meshgrid(1:2448,1:2048);
U = griddata(x,y,u,xGr,yGr);
V = griddata(x,y,v,xGr,yGr);

U = medfilt2(U,[10,10]);
V = medfilt2(V,[10,10]);

oxynano.velocity.U = U;
oxynano.velocity.V = V;
oxynano.velocity.Mag = sqrt(U.^2+V.^2);
oxynano.velocity.x = xGr;
oxynano.velocity.y = yGr;

oxynano.oxygen.x = xGr;
oxynano.oxygen.y = yGr;

%%
imagesc(V.*oxynano.mask),axis equal
caxis([0 50])

%% Downsampling

figure
Us = imresize(U.*oxynano.mask, 1/100, 'bilinear');
Vs = imresize(V.*oxynano.mask, 1/100, 'bilinear');
xGrS = imresize(xGr, 1/100, 'bilinear');
yGrS = imresize(yGr, 1/100, 'bilinear');
%quiver(Us,Vs,0.5,'k')

oxynano.velocity.Us = Us;
oxynano.velocity.Vs = Vs;
oxynano.velocity.xS = xGrS;
oxynano.velocity.yS = yGrS;

hold on
q = quiver(xGrS,-yGrS,Us,-Vs,1,'k')
plot(oxynano.surface2.x,-oxynano.surface2.y,'r-')
q.ShowArrowHead = 'on';

axis tight, axis equal

%%

scale = (max(oxynano.tracks.mag(:))+1-oxynano.tracks.mag).^2;

figure
scatter(oxynano.tracks.x,-oxynano.tracks.y,scale/10,oxynano.tracks.v,'.')
caxis([-5 5])

figure
scatter(oxynano.tracks.x,-oxynano.tracks.y,[],scale,'.')
caxis([-5 5])


%figure
%scatter(oxynano.tracks.x,-oxynano.tracks.y,[],oxynano.tracks.v-oxynano.tracks.vSmooth+oxynano.tracks.u-oxynano.tracks.uSmooth,'.')
%caxis([-5 5])


%% Smooth Surface and Slope calculation
smoothedSurf =  smooth(oxynano.surface2.y,1000);
oxynano.surface_smooth2.x = xS;
oxynano.surface_smooth2.y = smoothedSurf;

oxynano.surface_slope = smoothedSurf(2:end)-smoothedSurf(1:end-1);
oxynano.surface_slope(end+1) = oxynano.surface_slope(end);

angle = atand(oxynano.surface_slope);

oxynano.angle_surface = smooth(angle,2000);

figure
plot(smooth(angle,2000));


%% Polyfit for surface

x = oxynano.surface2.x;
y = oxynano.surface2.y;
slope = zeros(1,length(x));

for n = 1:length(x)
    shift = n+1000;
    if shift <= length(x)
        p = polyfit(x(n:shift),y(n:shift),1); 
    else
        p = polyfit(x(n:end),y(n:end),1);         
    end
    slope(n) = p(1);
end

angle = atand(slope);

figure
plot(angle)
oxynano.anglePoly_surface = angle;


%% Angle of upper velocities

x = oxynano.velocity.xS;
y = oxynano.velocity.yS;
U = oxynano.velocity.Us;
V = oxynano.velocity.Vs;

angle = atand(V./U);

hold on
pcolor(x,-y,angle),shading flat
quiver(x,-y,Us,-Vs,0.8,'k')
plot(oxynano.surface2.x,-oxynano.surface2.y,'r-')
plot(oxynano.surface2.x,-oxynano.surface2.y+300,'k-')


%%
hold on
plot(nanmean(x(8:10,:)),nanmean(angle(8:10,:)));

%oxynano.angleVel_surface = improfile(x,-y,angle,oxynano.surface2.x,-oxynano.surface2.y+300);
oxynano.angleVel_surface = interp1(nanmean(x(8:10,:)),nanmean(angle(8:10,:)),oxynano.surface2.x);
plot(smooth(oxynano.angleVel_surface,200))
%plot(C)
%% Angle Corrections

uS = oxynano.tracks.uSmooth;
vS = oxynano.tracks.vSmooth;
xS = oxynano.tracks.x;
yS = oxynano.tracks.y;

for n = 1:length(vS)
        theta = oxynano.angleVel_surface(round(xS(n)));
        rot = [uS(n) vS(n)]*[1;1i].*exp(-1i*theta*pi/180);
        %u(n) = real(rot);
        oxynano.tracks.v_rot(n) = imag(rot);
        oxynano.tracks.u_rot(n) = real(rot);
end

mag = sqrt(oxynano.tracks.v_rot.^2+oxynano.tracks.u_rot.^2);

[k,dist] = dsearchn([oxynano.surface2.x,oxynano.surface2.y],[oxynano.tracks.x,oxynano.tracks.y]);
oxynano.tracks.dist = dist;

%
sc = 1000000./dist.^2;
sc(sc>100) = 100;

figure
hold on
scatter(oxynano.tracks.x,-oxynano.tracks.y,sc,-oxynano.tracks.v_rot,'.')
colormap jet
caxis([-5 5])
plot(oxynano.surface2.x,-oxynano.surface2.y,'r-')
plot(oxynano.surface2.x,-oxynano.surface2.y+100,'k-')
quiver(x,-y,Us,-Vs,0.8,'k')

%% Grid Velocity

u = oxynano.tracks.u_rot;
v = oxynano.tracks.v_rot;
x = oxynano.tracks.x;
y = oxynano.tracks.y;

[xGr, yGr] = meshgrid(1:2448,1:2048);
U = griddata(x,y,u,xGr,yGr);
V = griddata(x,y,v,xGr,yGr);

U = medfilt2(U,[10,10]);
V = medfilt2(V,[10,10]);

oxynano.velocity.U_rot = U;
oxynano.velocity.V_rot = V;

%% Extract indivdual profiles

numb = 1:100;

C = zeros(length(oxynano.surface2.x),length(numb));
DC = zeros(length(oxynano.surface2.x),length(numb));
VelPr = zeros(length(oxynano.surface2.x),length(numb));

hold on
%pcolor(oxynano.velocity.x,-oxynano.velocity.y,oxynano.oxyraw.*oxynano.maskN), shading flat
%plot(oxynano.surface2.x,-oxynano.surface2.y,'k-')
%caxis([1 1.5]), axis equal

for m = 1:length(numb)
    temp = interp2(oxynano.oxygen.x,-oxynano.oxygen.y,oxynano.oxygen.oxyRaw.*oxynano.maskN,oxynano.surface2.x,-oxynano.surface2.y+m);
    tempVel = interp2(oxynano.velocity.x,-oxynano.velocity.y,oxynano.velocity.V_rot.*oxynano.maskN,oxynano.surface2.x,-oxynano.surface2.y+m);
    
    %plot(oxynano.surface2.x,-oxynano.surface2.y+m)
    %plot(smooth(temp))
    %pause(0.1)
    %plot(tempVel)
    DC(:,m) = temp;
    VelPr(:,m) = tempVel;
    m
end

C = nanmean(DC,2);
%
VelPr_min = quantile(VelPr',0.05);
VelPr_max = quantile(VelPr',0.95);
VelPr = nanmean(VelPr,2);
slope = zeros(length(oxynano.surface2.x),1);
slopeErr = zeros(length(oxynano.surface2.x),1);

%
for m = 1:length(oxynano.surface2.x)
    [temp,S,MU] = polyfit(numb,DC(m,:),1);
    slope(m) = temp(1);
    [y,DELTA] = polyval(temp,numb,S,MU);
    slopeErr(m) = nanmean(DELTA);
    % Delta is standard error
end

plot(smooth(slope,20))

oxynano.profiles.oxygen_grad = smooth(slope,20);
oxynano.profiles.oxygen_grad_std = smooth(slopeErr,20);
oxynano.profiles.oxygen_surface = C;
oxynano.profiles.velocity_surface = VelPr;
oxynano.profiles.velocity_surface_min = VelPr_min;
oxynano.profiles.velocity_surface_max = VelPr_max;

%%
plot(VelPr)




