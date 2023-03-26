function helix_vecs = make_helix_vecs(num_helix, num_segs_bt_intersections,worm_length, worm_radius,freq)
    
    num_helix_each_dir = num_helix/2;
    offsets = linspace(0,2*pi,num_helix_each_dir+1);
    offsets = offsets(1:end-1);
       
    t_intersections = [];
    for i = 1:length(offsets)
        for ii = 1:length(offsets)
            syms x
            xsol = solve(cos(freq*x + offsets(i)) == sin(freq*x + offsets(ii)), x, 'IgnoreProperties',true, 'ReturnConditions', true);
            y = matlabFunction(xsol.x);
            for k=-100:100
                t_intersections(end+1) = y(k);
            end
        end
    end
    
    t_intersections = t_intersections(t_intersections>0);
    t_intersections = t_intersections(t_intersections<worm_length);
    t_intersections = unique(round(sort(t_intersections),8));
    t_intersections = [0,t_intersections,worm_length];
    
    segmented_t = [];
    for i = 2:length(t_intersections)
        segmented_t = [segmented_t,linspace(t_intersections(i-1),t_intersections(i),num_segs_bt_intersections+1)];
    end
    segmented_t = unique(segmented_t);

    LH_helix = {};
    RH_helix ={};
    intersect_points = {};
    for i = 1:num_helix_each_dir
        sines = round(worm_radius + worm_radius*sin(freq*segmented_t + offsets(i)),5);
        cosines = round(worm_radius + worm_radius*cos(freq*segmented_t + offsets(i)),5);
        LH_helix{i} = [segmented_t',sines',cosines'];
        RH_helix{i} = [segmented_t',cosines',sines'];
        
        sines_intersect = round(worm_radius + worm_radius*sin(freq*t_intersections(2:end-1) + offsets(i)),5);
        cosines_intersect = round(worm_radius + worm_radius*cos(freq*t_intersections(2:end-1) + offsets(i)),5);

        intersect_points{i} = [t_intersections(2:end-1)', sines_intersect',cosines_intersect'];
    end    

%     helix_vecs = struct('left_handed',LH_helix,'right_handed',RH_helix,'intersection_points',intersect_points,'dummy',dummy);
    names = ["left_handed" "right_handed" "intersection_points"];
    data = {LH_helix, RH_helix, intersect_points};
    helix_vecs = dictionary(names,data);
end