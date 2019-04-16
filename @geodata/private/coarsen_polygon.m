        function polyobj = coarsen_polygon(polyobj,iboubox)
            % Coarsen a polygon using a brute-force vertex reordering algorithm.
            idx = find(isnan(polyobj(:,1)));
            C = mat2cell(polyobj,diff([0;idx]),2);
            new = cell(length(C),1);
            % for each segment
            for iii = 1 : length(C)
                j = 0 ; k = 0;
                [in] = inpoly(C{iii},iboubox(1:end-1,:));
                % Initialise for speed
                new{iii} = zeros(length(C{iii}),2);
                while j < length(C{iii})-1
                    j = j + 1 ;
                    if in(j) % point is in domain keep it
                        k = k + 1;
                        new{iii}(k,:) = C{iii}(j,:);
                    else % pt is out of domain
                        bd = min([j+200,length(in)-1]) ;
                        exte = min(200,bd - j);
                        if sum(in(j:bd))==0 % if next hundred points are out, then we can decimate
                            k = k + 1 ;
                            new{iii}(k,:) = C{iii}(j,:);
                            k = k + 1 ;
                            new{iii}(k,:) = C{iii}(j+exte,:);
                            j = j + exte ;
                        else % otherwise keep
                            k = k + 1 ;
                            new{iii}(k,:) = C{iii}(j,:);
                        end
                        
                    end
                end
                new{iii}(k+1,:) = [NaN NaN] ;
                new{iii}(k+2:end,:) = [];
            end
            polyobj = cell2mat(new);
        end