    function [BC,mask] = BCmatch(data,Intersections,iShape,jShape)

        rho = data.rho;
        flux = data.flux;
    
        BC = zeros(size(rho));
        
        Iij = Intersections(iShape,jShape);
        Iji = Intersections(jShape,iShape);

        if(~isempty(Iij.PtsMask))
            if(~Iij.Flip)
                BC(Iij.PtsMask) = rho(Iij.PtsMask) - rho(Iji.PtsMask);
                BC(Iji.PtsMask) = (Iij.Normal + Iji.Normal)*flux;
            else
                BC(Iij.PtsMask) = rho(Iij.PtsMask) - flipud(rho(Iji.PtsMask));
                BC(Iji.PtsMask) = flipud(Iij.Normal*flux) + Iji.Normal*flux;
                
            end
        end

        mask = (Iij.PtsMask | Iji.PtsMask);

    end