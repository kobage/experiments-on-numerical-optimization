import numpy as np

def build_profile(lst,maxX, fig_num):
    if len(lst) <= 1:
        print("Error! We need 2 or 3 files!")
        return

    a = np.loadtxt(lst[0])
    a1 = np.loadtxt(lst[1])
    a = np.concatenate(([a], [a1]), axis=0)
    if len(lst) == 3:
        a1 = np.loadtxt(lst[2])
        a = np.concatenate((a, [a1]), axis=0)

    b = np.copy(a)

    def max3(a, b, c):
        return np.minimum(a, np.minimum(b, c))

    columns = b.shape[1]
    for i in range(columns):
        if len(lst) == 2:
            m = np.minimum(b[0, i], b[1, i])
            b[0, i] /= m;
            b[1, i] /= m
        else:
            m = max3(b[0, i], b[1, i], b[2, i])
            b[0, i] /= m;
            b[1, i] /= m;
            b[2, i] /= m;

    def ro(b, t):
        return (np.count_nonzero(b <= t) / columns)

    text = 'Fig. ' + str(fig_num) +'.  Performance profile on [1,' + str(50) + '], based on CPU time'
    import matplotlib.pyplot as plt

    x = np.linspace(1, maxX, 3 * maxX)

    y = np.array([])
    z = np.array([])
    u = np.array([])

    for t in x:
        y = np.append(y, [ro(b[0], t)])
        z = np.append(z, [ro(b[1], t)])
        if len(lst) == 3:
            u = np.append(u, [ro(b[2], t)])

    plt.plot(x, y, linestyle='dashed', color="black", label=lst[0].replace('.txt', ''))
    plt.plot(x, z, linestyle='dotted', color="black", label=lst[1].replace('.txt', ''))
    if len(lst) ==3:
        plt.plot(x, u, linestyle='dashdot', color="black", label=lst[2].replace('.txt', ''))

    plt.xlim(1, maxX)
    plt.ylim(0, 1)
    leg = plt.legend();
    plt.xlabel(text)
    plt.show()


def build_profile_ax(ax,lst,maxX, fig_num):
    if len(lst) <= 1:
        print("Error! We need 2 or 3 files!")
        return

    a = np.loadtxt(lst[0])
    a1 = np.loadtxt(lst[1])
    a = np.concatenate(([a], [a1]), axis=0)
    if len(lst) == 3:
        a1 = np.loadtxt(lst[2])
        a = np.concatenate((a, [a1]), axis=0)

    b = np.copy(a)

    def max3(a, b, c):
        return np.minimum(a, np.minimum(b, c))

    columns = b.shape[1]
    for i in range(columns):
        if len(lst) == 2:
            m = np.minimum(b[0, i], b[1, i])
            b[0, i] /= m;
            b[1, i] /= m
        else:
            m = max3(b[0, i], b[1, i], b[2, i])
            b[0, i] /= m;
            b[1, i] /= m;
            b[2, i] /= m;

    def ro(b, t):
        return (np.count_nonzero(b <= t) / columns)

    x = np.linspace(1, maxX, 3 * maxX)

    y = np.array([])
    z = np.array([])
    u = np.array([])

    for t in x:
        y = np.append(y, [ro(b[0], t)])
        z = np.append(z, [ro(b[1], t)])
        if len(lst) == 3:
            u = np.append(u, [ro(b[2], t)])

    ax.plot(x, y, linestyle='dashed', color="black", label=lst[0].replace('.txt', ''))
    ax.plot(x, z, linestyle='dotted', color="black", label=lst[1].replace('.txt', ''))
    if len(lst) ==3:
        ax.plot(x, u, linestyle='dotted', color="black", label=lst[2].replace('.txt', ''))

    ax.set_xlim(1, maxX)
    ax.set_ylim(0, 1)
    ax.legend(loc='lower right');
    ax.set_xlabel('Fig. ' + str(fig_num) +'. Performance profile on [1,' + str(maxX) + '], based on CPU time')

