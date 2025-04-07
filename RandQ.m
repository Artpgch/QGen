function  [nM, flM] = RandQ(qN, spL, dZ, sM, syBw, smBw, minBw, maxBw, minFl, maxFl, spFl)

rM = Q1(qN, spL); % creating a random matrix
rM = Q2(dZ, qN, rM); % zeros on the diagonal
rM = Q3(sM, qN, rM); % symmetry of the matrix
nM = Q4(qN, smBw, minBw, rM, maxBw); % bandwidth matrix
nM = Q5(qN, syBw, nM); % symmetry of the bandwidth matrix
flM = Q6(qN, minFl, maxFl, spFl); % flow matrix
flM = Q7(dZ, qN, flM); % zeros on the diagonal


function rM = Q1(qN, spL)
    rM = ones(qN);
    for n = 1:spL
            rM = rM.*randi([0, 1], qN);
    end
end

function rM = Q2(dZ, qN, rM) % zeros on the diagonal
    if dZ
        for m = 1:qN
          rM(m,m) = 0;
        end
    end
end

function rM = Q3(sM, qN, rM) % symmetry of the matrix
    if sM
        for p = 1 : qN
            for q = p : qN
                rM(q, p) = rM(p, q); 
            end
        end
    end
end

    function nM = Q4(qN, smBw, minBw, rM, maxBw) % bandwidth matrix

    if smBw
        nM = minBw.*rM;
    else
        nM = (minBw + (maxBw-minBw).*rand(qN)).*rM;
    end

end

    function nM = Q5(qN, syBw, nM) % symmetry of the bandwidth matrix

     if syBw
        for p = 1 : qN
            for q = p : qN
                nM(q, p) = nM(p, q); 
            end
        end
    end

end

function flM = Q6(qN, minFl, maxFl, spFl) % flow matrix
    flM = (minFl + (maxFl - minFl).*rand(qN));
    for n = 1:spFl
            flM = flM.*randi([0, 1], qN);
    end
end

function flM = Q7(dZ, qN, flM) % zeros on the diagonal
    if dZ
        for m = 1:qN
           flM(m,m) = 0;
        end
    end
end

end
