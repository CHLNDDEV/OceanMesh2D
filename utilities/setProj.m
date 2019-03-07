        % kjr generic function to parse projected space options
        % returns del flag to delete overlapping elements when plotting
        % global models.
        function [del] = setProj(obj,proj,projtype)
            lon_mi = min(obj.p(:,1)); lon_ma = max(obj.p(:,1));
            lat_mi = min(obj.p(:,2)); lat_ma = max(obj.p(:,2));
            lat_mea = mean(obj.p(:,2));lon_mea = mean(obj.p(:,1));
            
            del = 0 ;
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
                    % Azimuthal type projections
                    m_proj(projtype,'lat',lat_mea,'long',lon_mea,...
                           'radius',50);
                    m_proj('get') ;
                    del = 1;
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