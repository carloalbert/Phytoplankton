function mle(data::AbstractVector{T}) where T <:AbstractFloat
    xmin = data[1]
    acc = zero(T)
    xlast = convert(T, Inf)
    ncount = 0
    for x in data
        if xlast == x
            continue
        end
        xlast = x
        ncount += 1
        acc += log(x / xmin)
    end
    ahat = 1 + ncount / acc
    stderr = (ahat - 1) / sqrt(ncount)
    return (ahat, stderr)
end

using Plots
#using MaximumLikelihoodPower

# Allocate memory for population:

P = Array{Float64, 1}(10^6);

# Initialize variables:

N = 1;
P[1] = 1;
i_max = 1;            # index of largest particle
x_max = P[i_max];     # trait of largest particle

omega_1 = 1;          # sets time unit (keep this fixed)
omega_2 = 1.1;
u = 1;                # sets trait unit (keep this fixed)
v = 10^3;
c = 10.0^(-5);

rng = MersenneTwister(1234);

# population updates:

for i = 1:10000

    # calculate times of earliest birth, death and "hit the wall" events:

    t_birth = randexp(rng, Float32)/(N*omega_2);
    t_death = randexp(rng, Float32)/(N^2*c);
    t_htw   = omega_1^(-1) * log((v+u)/(P[i_max]+u));

    # find out which event takes place earliest and at what time:

    event = indmin([t_birth,t_death,t_htw])
    t_star = [t_birth,t_death,t_htw][event]

    # let all particles grow until event happens:

    P[1:N] = exp(omega_1 * t_star)*P[1:N] + u*(exp(omega_1 * t_star) - 1)

    # execute event:

    if event == 1     # birth
        index = rand(1:N)
        P[index] = P[index]/2
        P[N+1] = P[index]
        N = N+1
    elseif event == 2 # death
        index = rand(1:N)
        P[index] = P[N]
        N = N-1
        else          # "hit the wall"
        index = i_max
        P[i_max] = v/2
        P[N+1] = v/2
        N = N+1
    end

    # find new largest particle, if necessary:

    if index == i_max
        i_max = indmax(P[1:N])
    end
end

# plot results:

histogram(log10(P[1:N]),yaxis = (:log10))

PP = sort(P[1:N]);
i=1;
while PP[i]<10 i=i+1 end
PP = PP[i:N];
mle(PP)
