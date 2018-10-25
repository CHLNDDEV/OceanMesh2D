        % kjr generic function to parse projected space options
        % returns del flag to delete overlapping elements when plotting
        % global models.
        function [del] = setProj(obj,proj,projtype)
            lon_mi = min(obj.p(:,1)); lon_ma = max(obj.p(:,1));
            lat_mi = min(obj.p(:,2)); lat_ma = max(obj.p(:,2));
            lat_mea = mean(obj.p(:,2));lon_mea = mean(obj.p(:,1));
            
            del = 0 ;
            if proj==0
                % normal geographic coordinates
                m_proj('equidist','lat',[lat_mi lat_ma],'long',[lon_mi lon_ma]) ;
            else
                switch projtype
                    case('stereo')
                        del = 1 ;
                        if lat_ma < 0
                            % center Antarctica
                            m_proj('stereo','lat',-90,...
                                'long',0.5*(lon_mi+lon_ma),...
                                'radius',lat_ma+90);
                        else
                            % center Arctic
                            lat_mi = max(-88.0001,lat_mi);
                            m_proj('stereo','lat',90,...
                                'long',0.5*(lon_mi+lon_ma),...
                                'radius',90-lat_mi);
                        end
                        m_proj('get') ;
                        del = 1 ;
                    case('trans')
                        disp('INFO: Default projected space...') ;
                        m_proj('Trans','lon',[lon_mi lon_ma],'lat',[lat_mi lat_ma]) ;
                        m_proj('get') ;
                    case('robinson')
                        m_proj('robinson')
                        m_proj('get') ;
                        del = 1 ; 
                    case('miller')
                        m_proj('miller','lat',[lat_mi lat_ma],'lon',[lon_mi lon_ma]);
                        m_proj('get') ;
                        del = 1 ; 
                    case('ortho')
                        m_proj('ortho','lat',[lat_mea],'long',[lon_mea]);
                        m_proj('get') ;
                    case('lambert')
                        m_proj('lambert','lon',[lon_mi lon_ma],'lat',[lat_mi lat_ma]);
                        m_proj('get') ;
                    otherwise
                        fprintf(1, [ ...
                            ' Unrecognized projected space. Available options are  \n', ...
                            ' Sterographic (stereo), Transverse Mercator (trans), Miller (miller) \n', ...
                            ' Lambet Conformal (lambert), Orthographic (ortho)\n', ...
                            ] ) ;
                end
            end
        end