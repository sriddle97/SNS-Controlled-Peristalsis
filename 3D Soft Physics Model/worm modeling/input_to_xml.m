function input_to_xml(helix_vecs,xml_skeleton_path,final_xml_path)
    
    dummy = helix_vecs("left_handed");
    LH_helixes = dummy{1,1};
    dummy = helix_vecs("right_handed");
    RH_helixes = dummy{1,1};
    dummy = helix_vecs("intersection_points");
    intersection_pts = dummy{1,1};
    
    xml_file = xml2struct(xml_skeleton_path);
    
    composite_body_copy = xml_file.mujoco.worldbody.body{1,1};
    ball_body_copy = xml_file.mujoco.worldbody.body{1,2};
    
    for i = 1:2*length(LH_helixes)
        xml_file.mujoco.worldbody.body{1,i} = composite_body_copy;
    end
    for i = 1:length(intersection_pts)*length(intersection_pts{1})
        xml_file.mujoco.worldbody.body{1,i+2*length(LH_helixes)} = ball_body_copy;
    end
    
    for i = 1:length(RH_helixes)
        xml_file.mujoco.worldbody.body{1,i}.composite.Attributes.prefix = append('RH',num2str(i-1));
        xml_file.mujoco.worldbody.body{1,i}.composite.Attributes.vertex = num2str(reshape(RH_helixes{i}.',1,[]));
        xml_file.mujoco.worldbody.body{1,i}.composite.geom.Attributes.size = '0.002';
        xml_file.mujoco.worldbody.body{1,i}.composite.Attributes.initial = 'none';
    
        xml_file.mujoco.worldbody.body{1,i+length(RH_helixes)}.composite.Attributes.prefix = append('LH',num2str(i-1));
        xml_file.mujoco.worldbody.body{1,i+length(RH_helixes)}.composite.Attributes.vertex = num2str(reshape(LH_helixes{i}.',1,[]));
        xml_file.mujoco.worldbody.body{1,i+length(RH_helixes)}.composite.geom.Attributes.size = '0.002';
        xml_file.mujoco.worldbody.body{1,i+length(RH_helixes)}.composite.Attributes.initial = 'none';
    
    end
    
    
    for i = 0:length(intersection_pts)-1
        for ii = 1:length(intersection_pts{1})
            xml_file.mujoco.worldbody.body{1,ii+2*length(RH_helixes) + i*length(intersection_pts{1})}.Attributes.name = append('ball',num2str(i),num2str(ii-1));
            xml_file.mujoco.worldbody.body{1,ii+2*length(RH_helixes) + i*length(intersection_pts{1})}.Attributes.pos = num2str(intersection_pts{i+1}(ii,:));
            xml_file.mujoco.worldbody.body{1,ii+2*length(RH_helixes) + i*length(intersection_pts{1})}.site.Attributes.name = append('site',num2str(i),num2str(ii-1));
        end
    
    
        struct2xml(xml_file,final_xml_path)
    end

end

