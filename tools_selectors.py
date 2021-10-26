from numpy import *


class ActionSelector:
    """ ActionSelector class
    A selector is determined by the catalog size N, time horizon T, maximum reward Rmax (optional) and eta (optional), correlated_arms (optional).

    Methods
    ----------
        feedback(r)
        Update action selector
        arm()
        Get a random state
        arm_corr()
        Get a correlated random state
        fractional_state()
        Get the distribution over the catalog

    """

    def __init__(self, N, eta=None, correlated_arms=False):
        self.__correlated_arms = correlated_arms
        # Optimal fixed learning rate
        self.__N = N
        self.__eta = eta
        # Fractional state (distribution over the catalog).
        self.__distribution = ones(N) / N
        # Randomly initialized integral state in {0,1,...,N-1}
        self.__state = self.__sample()

    def feedback(self, r):
        r = array(r)
        delta = 1e-7
        self.__distribution_next = self.__distribution * exp(self.__eta * r) + delta
        self.__distribution_next /= linalg.norm(self.__distribution_next, 1)

        if self.__correlated_arms:
            self.__update_corr_state()
        self.__distribution = self.__distribution_next
        # print(self.__distribution)

    def fractional_state(self):
        return self.__distribution

    def arm(self):
        # print('Here')
        # print(self.__distribution.max())
        if self.__correlated_arms:
            return self.__state
        else:
            return self.__sample()

    def __sample(self):
        return random.choice(arange(self.__N), p=self.__distribution)

    def __update_corr_state(self):
        state = self.__state
        distribution = self.__distribution
        distribution_next = self.__distribution_next
        if distribution_next[state] - distribution[state] >= 0:
            pass
        else:
            flow = self.__obtain_flow(distribution_next - distribution)
            for key in flow:
                i, j = key
                delta = flow[key]
                state, pi, pj = self.__elementary_jump(state, distribution[i],
                                                       distribution[j], i, j,
                                                       delta)
                distribution[i] = pi
                distribution[j] = pj
            self.__state = state

    def __elementary_jump(self, state, pi, pj, i, j, delta):
        if state == i:
            c = random.choice([0, 1], p=[(pi - delta) / pi, delta / pi])
            if c == 0:
                return i, pi - delta, pj + delta
            else:
                return j, pi - delta, pj + delta
        else:
            return state, pi - delta, pj + delta

    def __obtain_flow(self, deltas):
        I = where(deltas < 0)[0]
        J = where(deltas > 1e-10)[0]
        deltas = abs(deltas)
        flow = {}
        while I.size > 0 and J.size > 0:
            i, j = I[0], J[0]
            c = argmin([deltas[i], deltas[j]])
            flow[(i, j)] = deltas[i] if c == 0 else deltas[j]
            if c == 0:
                I = delete(I, 0)
                deltas[j] -= deltas[i]
                deltas[i] -= deltas[i]
            else:
                J = delete(J, 0)
                deltas[i] -= deltas[j]
                deltas[j] -= deltas[j]
        return flow


if __name__ == "__main__":
    # Testing on a single cache with a single slot, (reward is 1 if it's in the cache otherwise 0)
    N = 100
    T = 1_000
    s = .6


    def zipf_distribution(s, N):
        c = sum((1 / arange(1, N + 1) ** s))
        return arange(1, N + 1) ** (-s) / c


    Requests = random.choice(arange(N), p=zipf_distribution(s, N), size=T)
    selector = ActionSelector(N, T)
    states = []
    for q in Requests:
        r = zeros(N)
        r[q] = 1
        selector.feedback(r)
        states.append(selector.arm())
    print(('State at T: ', states[-1]))
