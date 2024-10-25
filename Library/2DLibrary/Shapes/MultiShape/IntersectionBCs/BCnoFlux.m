    function [BC,mask] = BCnoFlux(data,Intersections,iShape,jShape)

        rho = data.rho;
        flux = data.flux;

        BC = zeros(size(rho));
        
        Iij = Intersections(iShape,jShape);
        Iji = Intersections(jShape,iShape);

        if(~isempty(Iij.PtsMask))

            BC(Iij.PtsMask) = Iij.Normal * flux;
            BC(Iji.PtsMask) = Iji.Normal * flux;

        end
        
        mask = (Iij.PtsMask | Iji.PtsMask);

    end