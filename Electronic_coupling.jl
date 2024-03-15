using LinearAlgebra
using Printf
elapsed_time = @elapsed begin

function find_nbasis(fchk)
    nbasis=0
    open(fchk, "r") do file
        seekstart(file)
        for line in eachline(file)
            words = split(strip(line))
            if length(words)<6
                continue
            end
            if words[3]=="basis" && words[4]=="functions"
                nbasis = parse(Int64, words[6])
            end
        end
    end
    return nbasis
end

function find_moene(fchk)
    moene=Float64[]
    nbasis = 0.0
    open(fchk, "r") do file
        seekstart(file)
        for line in eachline(file)
            words = split(strip(line))
            if length(words) < 6
                continue
            end
            if words[1]=="Alpha" && words[2]=="Orbital" && words[3]=="Energies" && words[4]=="R" && words[5]=="N="
                nbasis = parse(Float64, words[6])
                break
            end
        end
        num_lines = Int(ceil(nbasis/5.0))
        for i in collect(1:num_lines)
            line = readline(file)
            words = split(strip(line))
            append!(moene, parse(Float64, ss) for ss in words)
        end
    end
    return moene
end

function find_mocoeff(fchk)
    mocoeff = Float64[]
    ncoeff = 0.0
    open(fchk, "r") do file
        seekstart(file)
        for line in eachline(file)
            words = split(strip(line))
            if length(words) < 6
                continue
            end
            if words[1]=="Alpha" && words[2]=="MO" && words[3]=="coefficients" && words[4]=="R" && words[5]=="N="
                ncoeff = parse(Float64, words[6])
                break
            end
        end
        num_lines = Int(ceil(ncoeff/5.0))
        for i in collect(1:num_lines)
            line = readline(file)
            words = split(strip(line))
            append!(mocoeff, parse(Float64, ss) for ss in words)
        end
    end
    moc = reshape(mocoeff, Int(sqrt(ncoeff)),  Int(sqrt(ncoeff)))
    return moc
end

function find_overlap(log, nbasis)
    overlap=zeros(nbasis, nbasis)
    open(log, "r") do file
        for line in eachline(file)
            words = split(strip(line))
            if length(words)<3
                continue
            end
            if words[1]=="***" && words[2]=="Overlap" && words[3]=="***"
                break
            end
        end
        for i in collect(0:(Int(floor(nbasis/5.0))-1))
            line = readline(file)
            line = readline(file)
            overlap[i*5+1, i*5+1]=parse(Float64, replace(split(line)[2], "D" => "e"))
            line = readline(file)
            overlap[i*5+2, i*5+1]=parse(Float64, replace(split(line)[2], "D" => "e"))
            overlap[i*5+2, i*5+2]=parse(Float64, replace(split(line)[3], "D" => "e"))
            line = readline(file)
            overlap[i*5+3, i*5+1]=parse(Float64, replace(split(line)[2], "D" => "e"))
            overlap[i*5+3, i*5+2]=parse(Float64, replace(split(line)[3], "D" => "e"))
            overlap[i*5+3, i*5+3]=parse(Float64, replace(split(line)[4], "D" => "e"))
            line = readline(file)
            overlap[i*5+4, i*5+1]=parse(Float64, replace(split(line)[2], "D" => "e"))
            overlap[i*5+4, i*5+2]=parse(Float64, replace(split(line)[3], "D" => "e"))
            overlap[i*5+4, i*5+3]=parse(Float64, replace(split(line)[4], "D" => "e"))
            overlap[i*5+4, i*5+4]=parse(Float64, replace(split(line)[5], "D" => "e"))
            for j in (i*5+5):nbasis
                line = readline(file)
                overlap[j, i*5+1]=parse(Float64, replace(split(line)[2], "D" => "e"))
                overlap[j, i*5+2]=parse(Float64, replace(split(line)[3], "D" => "e"))
                overlap[j, i*5+3]=parse(Float64, replace(split(line)[4], "D" => "e"))
                overlap[j, i*5+4]=parse(Float64, replace(split(line)[5], "D" => "e"))
                overlap[j, i*5+5]=parse(Float64, replace(split(line)[6], "D" => "e"))
            end
        end
        line=readline(file)
        k = nbasis % 5
        line=readline(file)
        if k==1
            overlap[nbasis, nbasis]=parse(Float64, replace(split(line)[2], "D" => "e"))
        elseif k==2
            overlap[nbasis-1, nbasis-1]=parse(Float64, replace(split(line)[2], "D" => "e"))
            line=readline(file)
            overlap[nbasis, nbasis-1]=parse(Float64, replace(split(line)[2], "D" => "e"))
            overlap[nbasis, nbasis]=parse(Float64, replace(split(line)[3], "D" => "e"))
        elseif k==3 
            overlap[nbasis-2, nbasis-2]=parse(Float64, replace(split(line)[2], "D" => "e"))
            line=readline(file)
            overlap[nbasis-1, nbasis-2]=parse(Float64, replace(split(line)[2], "D" => "e"))
            overlap[nbasis-1, nbasis-1]=parse(Float64, replace(split(line)[3], "D" => "e"))
            line=readline(file)
            overlap[nbasis, nbasis-2]=parse(Float64, replace(split(line)[2], "D" => "e"))
            overlap[nbasis, nbasis-1]=parse(Float64, replace(split(line)[3], "D" => "e"))
            overlap[nbasis, nbasis]=parse(Float64, replace(split(line)[4], "D" => "e"))
        else k==4
            overlap[nbasis-3, nbasis-3]=parse(Float64, replace(split(line)[2], "D" => "e"))
            line=readline(file)
            overlap[nbasis-2, nbasis-3]=parse(Float64, replace(split(line)[2], "D" => "e"))
            overlap[nbasis-2, nbasis-2]=parse(Float64, replace(split(line)[3], "D" => "e"))
            line=readline(file)
            overlap[nbasis-1, nbasis-3]=parse(Float64, replace(split(line)[2], "D" => "e"))
            overlap[nbasis-1, nbasis-2]=parse(Float64, replace(split(line)[3], "D" => "e"))
            overlap[nbasis-1, nbasis-1]=parse(Float64, replace(split(line)[4], "D" => "e"))
            line=readline(file)
            overlap[nbasis, nbasis-3]=parse(Float64, replace(split(line)[2], "D" => "e"))
            overlap[nbasis, nbasis-2]=parse(Float64, replace(split(line)[3], "D" => "e"))
            overlap[nbasis, nbasis-1]=parse(Float64, replace(split(line)[4], "D" => "e"))
            overlap[nbasis, nbasis]=parse(Float64, replace(split(line)[5], "D" => "e"))
        end
    end
    for i in 1:nbasis
        for j in 1:nbasis
            overlap[i, j] = overlap[j, i]
        end
    end
    return overlap
end


monomer1fchk =  "1" * ".fchk"
monomer1log = "1" * ".log"
monomer2fchk =  "801" * ".fchk"
monomer2log = "801" * ".log"
dimerfchk = "dimer.fchk"
dimerlog = "dimer.log"

NBasis = find_nbasis(dimerfchk)
MOene = find_moene(dimerfchk)
MOcoeff = find_mocoeff(dimerfchk)
Overlap = find_overlap(dimerlog, NBasis)
Fock = Overlap * MOcoeff * diagm(MOene) * inv(MOcoeff)  # S * C * e * C^-1

NBasis1 = find_nbasis(monomer1fchk)
MOene1 = find_moene(monomer1fchk)
MOcoeff1 = find_mocoeff(monomer1fchk)

NBasis2 = find_nbasis(monomer2fchk)
MOene2 = find_moene(monomer2fchk)
MOcoeff2 = find_mocoeff(monomer2fchk)

cfrag = zeros(NBasis, NBasis)
cfrag[1:NBasis1, 1:NBasis1] = MOcoeff1
cfrag[NBasis1+1:end, NBasis1+1:end] = MOcoeff2


ihomo1, ilumo1 = 73, 73
ihomo2, ilumo2 = 73, 73
@printf("%5s%5s%10s%10s%10s%10s%10s%10s%10s%13s\n", "MO1", "MO2", "e1 eV", "e2 eV", "J meV","S", "ee1 eV", "ee2 eV", "Je meV", "Je kcal/mol")

for i in ihomo1:ilumo1
    for j in (ihomo2+NBasis1):(ilumo2+NBasis1)
        e1 = dot(cfrag[:, i], Fock*cfrag[:, i])*27.21138505
        e2 = dot(cfrag[:, j], Fock*cfrag[:, j] )*27.21138505
        J12 = dot(cfrag[:, i], Fock*cfrag[:, j])*27.21138505
        S12 = dot(cfrag[:, i], Overlap*cfrag[:, j])
        ee1 = 0.5*((e1+e2)-2.0*J12*S12+(e1-e2)*sqrt(1.0-S12^2))/(1-S12^2)
        ee2 = 0.5*((e1+e2)-2.0*J12*S12-(e1-e2)*sqrt(1.0-S12^2))/(1-S12^2)
        Je12 = (J12 - 0.5 * (e1 + e2) * S12) / (1 - S12^2)
        @printf("%5d%5d%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%13.3f\n", i, j-NBasis1, e1, e2, J12 * 1000, S12,   ee1, ee2, Je12 * 1000, Je12 * 23.05)
    end
end


end
println("Elapsed time: $elapsed_time seconds")