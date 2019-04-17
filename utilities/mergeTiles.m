 function m3 = mergeTiles(tiles) 
            % Merge mesh tiles together 
            numTiles = length(tiles) ; 
            
            for i = 1 : numTiles-1
                
                if i == 1
                   dmy=load(tiles{i}) ;
                   m1 = dmy.m ; 
                else
                   m1 = m3 ;
                end
                dmy=load(tiles{i+1}) ; 
                m2 = dmy.m ; 

                
                dt = delaunayTriangulation(m1.p) ;
                dt.Points(end+(1:length(m2)),:) = m2.p ;
                
                m3 = msh() ;
                m3.p = dt.Points ;
                m3.t = dt.ConnectivityList ;
                
                bc = baryc(m3) ;
                
                bou1 = m1.getBoundaryOfMesh ;
                ee1   = Get_poly_edges(bou1) ;
                in1  = inpoly(bc,bou1,ee1) ;
                
                bou2 = m2.getBoundaryOfMesh ;
                ee2   = Get_poly_edges(bou2) ;
                in2  = inpoly(bc,bou2,ee2) ;
                
                rm = ~in1 & ~in2 ;
                
                m3.t(find(rm),:) = [ ] ;
                
                [m3.p,m3.t]=fixmesh(m3.p,m3.t) ;
                
                m3 = renum(m3) ;
            end
            
        end