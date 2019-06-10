function [del] = setProj(obj,proj,projtype)
    %del = setProj(obj,proj,projtype)
    % kjr generic function to parse projected space options
    % returns del flag to delete overlapping elements when plotting
    % global models.

    % process bounds of mesh
    lon_mi = min(obj.p(:,1)); lon_ma = max(obj.p(:,1));
    lat_mi = min(obj.p(:,2)); lat_ma = max(obj.p(:,2));
    lat_mea = mean(obj.p(:,2)); lon_mea = mean(obj.p(:,1));
    % some defaults
    rad = 100; rot = 15;
    del = 0 ;

    % process the projtype as a varargin type
    I = find(strcmp(projtype,'long'));
    if ~isempty(I)
        if size(projtype(I+1)) == 1
            lon_mea = projtype{I+1};
        else
            lon_mi = projtype{I+1}(1); lon_ma = projtype{I+1}(2);
        end
    end
    I = find(strcmp(projtype,'lat'));
    if ~isempty(I)
        if size(projtype(I+1)) == 1
            lat_mea = projtype{I+1};
        else
            lat_mi = projtype{I+1}(1); lat_ma = projtype{I+1}(2);
        end
    end
    I = find(strcmp(projtype,'rad'));
    if ~isempty(I)
        rad = projtype{I+1};
    end
    I = find(strcmp(projtype,'rot'));
    if ~isempty(I)
        rot = projtype{I+1};
    end
    if ~ischar(projtype)
       projtype = projtype{1};
    end
    if proj == 0
        % normal geographic coordinates
        m_proj('equi','lat',[lat_mi lat_ma],'long',[lon_mi lon_ma]) ;
        del = 1;
    else
        if startsWith(projtype,'ste','IgnoreCase',true)
            % Special treatment of Stereographic projection
            if lat_ma < 0
                % center Antarctica
                m_proj(projtype,'lat',-90,...
                      'long',0.5*(lon_mi+lon_ma),...
                      'radius',lat_ma+90);
            else
                % center Arctic
                lat_mi = max(-88.0001,lat_mi);
                m_proj(projtype,'lat',90,...
                      'long',0.5*(lon_mi+lon_ma),...
                      'radius',90-lat_mi,'rot',180);
            end
            m_proj('get') ;
        elseif startsWith(projtype,'ort','IgnoreCase',true) || ...
               startsWith(projtype,'gno','IgnoreCase',true) || ...
               startsWith(projtype,'azi','IgnoreCase',true) || ...
               startsWith(projtype,'sat','IgnoreCase',true)
            % Azimuthal type projections;
            m_proj(projtype,'lat',lat_mea,'long',lon_mea,...
                   'radius',rad,'rot',rot);
            m_proj('get') ;
            del = 0;
        elseif startsWith(projtype,'obl','IgnoreCase',true)
            % Oblique Mercator projection
            asp = (lon_ma-lon_mi)/(lat_ma - lat_mi);
            dir = 'hor';
            if asp > 1 
                asp = 1./asp; dir = 'ver';
            end
            m_proj(projtype,'lon',[lon_mi lon_ma],...
                            'lat',[lat_mi lat_ma],...
                            'aspect',asp,'dir',dir) ;
            m_proj('get') ;
            del = 1;
        else
            % Cylindrical, Conic or Global type projections
            del = 1;
            m_proj(projtype,'lon',[lon_mi lon_ma],...
                            'lat',[lat_mi lat_ma]) ;
            m_proj('get') ;
        end
    end
end