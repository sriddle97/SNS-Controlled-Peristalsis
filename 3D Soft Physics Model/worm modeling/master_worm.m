clc;
clear;

% Originally written by Clayton Jackson

% first we run this to make the points for the composites and sites
% make_helix_vec(num_helix, num_segs_bt_intersections, worm_length, worm_radius,freq)
% the freq as in sin(freq*x)  (just how tight the spiral is)

%originally had make_helix_vecs(12, 3, 0.42, 0.1,14) and in worm_skeleton had bend=15e8 
% now have: bend=5e8 to match polyethylene
%           tendon stiffnesses set to 19.9
%           muscle scale set to 350000
%           6 tendon springs per segment

helix_vecs = make_helix_vecs(12, 3, 0.56, 0.16,14);     %radius of 0.16m for Yifan robot
% this is in a bit of a weird format so lets fix that
%       LH - left handed (counterclockwise)
%       RH - right handed (clockwise)
dummy = helix_vecs("left_handed");
LH_helixes = dummy{1,1};
dummy = helix_vecs("right_handed");
RH_helixes = dummy{1,1};
dummy = helix_vecs("intersection_points");
intersection_pts = dummy{1,1};

%% Now lets put these into an xml file
%       worm_skeleton.xml just has two bodies
%             1 composite
%             1 ball
%       just uses these as a framework for the worm_new.xml

input_to_xml(helix_vecs,'worm_skeleton.xml','worm_new.xml')

%% Making the excludes and connects
%       the input_to_xml function only inputs the bodies and sites
%       for the excludes and connects:
%             run this section 
%             copy from the command window
%                 paste the excludes in the <contact> section of worm_new.xml
%                 paste the connects in the <equality> section 
make_excludes(helix_vecs)
fprintf('\n\n\n\n\n\n\n\n\')
make_connects(helix_vecs)

% for adding tendons: the sites that you attach the tendon to will all have
% the same second index and the first index will go up to num_segments/2. For example: 
% for a worm with triangular muscle ring, the tendon section of the
% xml looked like this: 
% <tendon>
%   <spatial name="ring0" >
%         <site site="site02" />
%         <site site="site12" />
%         <site site="site22" />
%         <site site="site02" />
%    </spatial>
%   <spatial name="ring1" >
%         <site site="site02" />
%         <site site="site12" />
%         <site site="site22" />
%         <site site="site02" />
%    </spatial>
%   <spatial name="ring0" >
%         <site site="site04" />
%         <site site="site14" />
%         <site site="site24" />
%         <site site="site04" />
%    </spatial>
%   <spatial name="ring2" >
%         <site site="site06" />
%         <site site="site16" />
%         <site site="site26" />
%         <site site="site06" />
%    </spatial>
% </tendon


% Then for muscles you just specify the tendon and you can add other properties or just use mujoco's default
% <actuator>
%     <muscle name="ring0_muscle" tendon="ring0" class="muscle" />
%     <muscle name="ring1_muscle" tendon="ring1" class="muscle" />
%     <muscle name="ring2_muscle" tendon="ring2" class="muscle" />
% </actuator>

%% Plot the worm and connection points
%   just in case you want to see it in here

% strand length for doing mass/length calc (tuning model params)
strand_length = 0;
for i = 2:length(LH_helixes{1,1})
    point_diffs = LH_helixes{1,1}(i,:)-LH_helixes{1,1}(i-1,:);
    part_length = sqrt(sum(point_diffs.^2));
    strand_length = strand_length+part_length;
end
disp(strand_length)
disp(strand_length*12)      %full length of tubing in robot (12 strands)

% rhombus side length for doing height or diameter plots (currently assumes
% 3 discretizations per rhmobus side length (hence 5-7))
% l0 for Yifan robot is 3.5-3.6in (0.089-0.092m)
l0 = 0;
for i = 5:7
    point_diffs = LH_helixes{1,1}(i,:)-LH_helixes{1,1}(i-1,:);
    part_length = sqrt(sum(point_diffs.^2));
    l0 = l0+part_length;
end
disp(l0)

figure(1)
clf('reset') 
hold on

for i = 1: length(LH_helixes)
    plot3(LH_helixes{1,i}(:,1),LH_helixes{1,i}(:,2),LH_helixes{1,i}(:,3),'b','Linewidth',1.2)
    plot3(RH_helixes{1,i}(:,1),RH_helixes{1,i}(:,2),RH_helixes{1,i}(:,3),'b','Linewidth',1.2)
    plot3(intersection_pts{1,i}(:,1),intersection_pts{1,i}(:,2),intersection_pts{1,i}(:,3),'rx','MarkerSize',10,'LineWidth',1.5)
end

axis equal