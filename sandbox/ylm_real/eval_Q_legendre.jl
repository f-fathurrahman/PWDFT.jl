function eval_Q_legendre( lmax::Int64, cost::Float64, sent::Float64 )

    Q = zeros(lmax+1,lmax+1)

    fpi = 4*pi

    for l = 0:lmax
        
        c = sqrt( (2*l+1) / fpi )
        
        if l == 0
            
            Q[1,1] = 1.0

        else if l == 1
            
            Q[2,1] = cost
            Q[2,2] = -sent/sqrt(2.0)

        else
            for m = 0:l-2
                Q[l+1,m+1] = cost *
                             (2*l-1)/sqrt(l*l-m*m) * Q[l,m+1] -
                             sqrt( (l-1)*(l-1)-m*m ) / sqrt(l*l-m*m) * Q[l-1,m+1]
            end
            Q[l+1,l] = cost * sqrt(2*l-1) * Q[l,l]
            Q[l+1,l+1] = -sqrt( 2*l-1 ) / sqrt(2*l)*sent * Q[l,l]
            
        end

    end


end
