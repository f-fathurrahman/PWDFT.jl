"""
Taken from fft_support file in QE
"""


function allowed( nr::Int64 )
# find if the fft dimension is a good one
# a "bad one" is either not implemented (as on IBM with ESSL)
# or implemented but with awful performances (most other cases)

    factors = [2, 3, 5, 7, 11]

    # find the factors of the fft dimension
    mr  = nr
    pwr = zeros(Int64,5)
    do_factor_loop = true

    for i = 1:5
        if !do_factor_loop
            break
        end
        fac = factors[i]
        maxpwr = round(log(Float64(mr))/log(Float64(fac))) + 1
        for p = 1:maxpwr
            if ( mr == 1 )
                do_factor_loop = false
                break
            end
            if mr%fac == 0
                mr = mr/fac
                pwr[i] = pwr[i] + 1
            end
        end
    end

    if nr != ( mr * 2^pwr[1]* 3^pwr[2] * 5^pwr[3] * 7^pwr[4] * 11^pwr[5] )
        println("Error in fft_support.jl: alllowed")
        exit()
    end


    if mr != 1
        # fft dimension contains factors > 11 : no good in any case
        return false
    else
        # fftw and all other cases: no factors 7 and 11
        return ( pwr[4] == 0 ) & ( pwr[5] == 0 )
    end

end


function good_fft_order( nr::Int64 )
#
#    This function find a "good" fft order value greater or equal to "nr"
#
#    nr  (input) tentative order n of a fft
#
#    Output: the same if n is a good number
#         the closest higher number that is good
#         an fft order is not good if not implemented (as on IBM with ESSL)
#         or implemented but with awful performances (most other cases)
#
    nfftx = 2049

    new_nr = nr

    while !allowed(new_nr) & ( new_nr <= nfftx )
        new_nr = new_nr + 1
    end

    if new_nr > nfftx
        println("Error in good_fft_order: fft_order too large")
        exit()
    end
    return new_nr

end
