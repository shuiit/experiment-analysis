classdef HullToCsv<handle
    %HullToCsv converts hull to csv file
    %   this class generates a csv file from Shull - a structure of the hull of a fly. 
    %   it gets the values from a user defined features, save it in a table and save it to csv. 
    
    properties
        name_of_features_wings
        name_of_features_body
        name_of_general_features
        name_of_video_features
        name_of_features_wings_angle
        name_of_features_body_angle
        name_of_features_body_coords
        hull
        frame_table
        time_table
        path_to_save
        rotmat
    end
    
    methods
        function obj = HullToCsv(hull,path_to_save)
            
            obj.name_of_features_wings = {'span','chord','nrml'};
            obj.name_of_features_body = {'strkPlan','X','Y','Z'};
            obj.name_of_general_features = {'frames'};
            obj.name_of_video_features = {'timeframe'};
            obj.name_of_features_wings_angle = {'phi','theta','psi'};
            obj.name_of_features_body_angle = {'pitch','yaw','roll'};
            obj.name_of_features_body_coords = {'CM_real'}
            obj.hull = hull;
            obj.frame_table = array2table(obj.hull.video.timeframe',"VariableNames",{'time'});
            obj.time_table = array2table(obj.hull.frames',"VariableNames",{'frames'});
            obj.path_to_save = path_to_save;
            writematrix(obj.hull.rotmat_EWtoL,[obj.path_to_save,'_ew_to_lab_rotmat.csv']) 

            
        end
        function [table_to_keep,feat_for_header] = get_features_from_structure(obj,struct,name_of_features,feature_name)
            % Get features from structure. concatenate them and generate a
            % table
            array_to_keep = [];
            for feat_idx =  1:1:length(name_of_features)
                feature = struct.(name_of_features{feat_idx});
                if size(feature,2) > size(feature,1)
                    feature = feature';
                end
                array_to_keep = horzcat(array_to_keep,feature);
            end
            
            feat_for_header = {strcat(name_of_features,'_',feature_name)};
            if size(feature,2) > 1
                axes_name = {'x','y','z'};
                for ax_idx = 1:1:length(name_of_features)
                    feat_for_header{ax_idx} =strcat(name_of_features{ax_idx},'_',axes_name,'_',feature_name);
                end
            end
            table_to_keep = array2table(array_to_keep,"VariableNames",[feat_for_header{:}]);
        end
        
        function angles_to_csv(obj)
            % Generate table for the angles variables and save as a csv
            [table_rw_angle] = obj.get_features_from_structure(obj.hull.rightwing.angles,obj.name_of_features_wings_angle,'rw');
            [table_lw_angle] = obj.get_features_from_structure(obj.hull.leftwing.angles,obj.name_of_features_wings_angle,'lw');
            [table_body_angle] = obj.get_features_from_structure(obj.hull.body.angles,obj.name_of_features_body_angle,'body');
            [table_body_coords] = obj.get_features_from_structure(obj.hull.body.coords,obj.name_of_features_body_coords,'body');

            vectors_table_for_csv = [obj.frame_table,obj.time_table,table_rw_angle,table_lw_angle,table_body_angle,table_body_coords];
            writetable(vectors_table_for_csv,[obj.path_to_save,'_angles_cm.csv'])
        end
        
        function vectors_to_csv(obj)
            % Generate table for the vectors variables and save as a csv
            
            [table_rw] = obj.get_features_from_structure(obj.hull.rightwing.vectors,obj.name_of_features_wings,'rw');
            [table_lw] = obj.get_features_from_structure(obj.hull.leftwing.vectors,obj.name_of_features_wings,'lw');
            [table_body] = obj.get_features_from_structure(obj.hull.body.vectors,obj.name_of_features_body,'body');

            vectors_table_for_csv = [obj.frame_table,obj.time_table,table_rw,table_lw,table_body];
            writetable(vectors_table_for_csv,[obj.path_to_save,'_vectors.csv'])
            
        end
        
        
    end
end

