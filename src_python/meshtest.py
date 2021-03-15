import numpy as np
import matplotlib.pyplot as plt

# specify params


def meshgenA(nzf,
             niw,
             nrw,
             iwMinMax=[0.1, 6.5],
             rwMinMax=[0.001, 0.85],
             plti=False):

    n = nzf * niw * (nrw + 1)

    sensitivity = 0.19  # 0 means no movement, 1 means max distance is init_dist

    # compute grid shape based on number of points
    num_y = np.int32(nzf * niw)
    num_x = np.int32(nrw + 1)

    # create regularly spaced neurons
    x = np.linspace(iwMinMax[0], iwMinMax[1], num_x, dtype=np.float32)
    y = np.linspace(rwMinMax[0], rwMinMax[1], num_y, dtype=np.float32)
    coords = np.stack(np.meshgrid(x, y), -1).reshape(-1, 2)

    # compute spacing
    init_dist = np.min((x[1] - x[0], y[1] - y[0]))
    min_dist = init_dist * (1 - sensitivity)

    print(init_dist, min_dist)
    assert init_dist >= min_dist

    # perturb points
    max_movement = (init_dist - min_dist) / 2
    noise = np.random.uniform(
        low=-max_movement, high=max_movement, size=(len(coords), 2))
    coords += noise

    # plot
    if plti:
        plt.figure(figsize=(10, 10))
        plt.scatter(coords[:, 0], coords[:, 1], s=3)
        plt.show()

    return coords


def meshgenB(nzf,
             niw,
             nrw,
             mindi,
             iwMinMax=[0.1, 6.5],
             rwMinMax=[0.001, 0.85],
             plti=False):

    nziw = nzf * niw
    nzrw = nzf * niw * nrw

    # compute grid shape based on number of points
    iwgrd = np.random.uniform(iwMinMax[0], iwMinMax[1], nzf * niw)

    rwgrd = [np.random.uniform(rwMinMax[0], rwMinMax[1])]

    iter = 0

    while ((len(rwgrd) < nziw * nzrw) & (iter < 10**5)):

        trl = np.random.uniform(rwMinMax[0], rwMinMax[1])
        dist = [np.linalg.norm(trl - elem) for elem in rwgrd]
        if np.min(dist) > mindi:
            rwgrd.append(trl)
        iter += 1

    tmp = np.array(rwgrd).reshape((-1, nziw))
    grd = np.append(tmp, iwgrd).reshape((-1, nziw))

    ycord = np.reshape(grd[:-1, :], (1, -1))
    xcord = np.repeat(iwgrd, nrw, axis=0)

    coords = []
    for iw in range(len(iwgrd)):
        for rw in range(nrw):
            coords.append([iwgrd[iw], grd[iw, rw]])
    coord = np.array(coords)
    if plti:
        plt.figure(figsize=(10, 10))
        plt.scatter(coord[:, 0], coord[:, 1], s=3)
        plt.show()

    return tmp.reshape((-1, 1)), iwgrd
