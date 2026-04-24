import math
import numpy as np
from copy import deepcopy
from gdsCAD import *
from circuit_elements.circle_points import circle


def resonator(L, L_sq, n, a, w):

    s = L / L_sq * w

    pd = 2 * np.pi / n
    spd = np.sin(pd / 2)
    cpd = np.cos(pd / 2)

    ra = a / (2 * spd)
    rb = (s + n * a * (1 / (2 * spd) - 0.5)) / (n * (1 + spd))

    start_point = [(rb * cpd, rb * spd)]

    points = []
    points.append(start_point[0])

    p = pd / 2 + pd * np.arange(n + 1)

    for i in range(n):

        if i % 2 == 0:
            points.append((ra * np.cos(p[i]), ra * np.sin(p[i])))
            points.append((ra * np.cos(p[i + 1]), ra * np.sin(p[i + 1])))
        else:
            points.append((rb * np.cos(p[i]), rb * np.sin(p[i])))
            points.append((rb * np.cos(p[i + 1]), rb * np.sin(p[i + 1])))

    print(s)
    print(calc_length(points))
    points = round_points(points)
    print(calc_length(points))

    superinductor = core.Path(points, w)

    elist = core.Elements()
    elist.add(superinductor)

    return elist


def flatronium(L, Ld, L_sq, n, a, w_ebeam):

    w = 0.05
    s = L / L_sq * w
    sd = Ld / L_sq * w

    po = 30.0 / 180.0 * np.pi
    pd = (2 * np.pi - po) / (n - 1)
    spd = np.sin(pd / 2)

    ra = a / (2 * spd)
    rb = 2 * sd / (n * (1 - spd)) + ra

    #b = 4 * sd * spd / (n * (1 - spd)) + a

    rc = (s + (ra + rb) * (1 - spd) * n / 2) / (n + (n - 2) * spd)

    # circle
    alpha = 0.25
    beta = 1.0

    start_point = [((rc + beta) * np.cos(po / 2) + alpha, (rc + beta) * np.sin(po / 2) - alpha)]
    start_point = round_points(start_point)

    points = []
    points.append(start_point[0])
    points.append((start_point[0][0] - alpha, start_point[0][1] + alpha))

    p = po / 2 + pd * np.arange(n)

    for i in range(0, n - 1):

        if i % 2 == 0:
            if i < n // 2:
                points.append((ra * np.cos(p[i]), ra * np.sin(p[i])))
                points.append((ra * np.cos(p[i + 1]), ra * np.sin(p[i + 1])))
            else:
                points.append((rb * np.cos(p[i]), rb * np.sin(p[i])))
                points.append((rb * np.cos(p[i + 1]), rb * np.sin(p[i + 1])))
        else:
            points.append((rc * np.cos(p[i]), rc * np.sin(p[i])))
            points.append((rc * np.cos(p[i + 1]), rc * np.sin(p[i + 1])))

    points.append((start_point[0][0] - alpha, - start_point[0][1] - alpha))
    points.append((start_point[0][0], - start_point[0][1]))

    print(calc_length(points))
    points = round_points(points)
    print(calc_length(points))

    superinductor = core.Path(points, w_ebeam)

    elist = core.Elements()
    elist.add(superinductor)

    # pads

    pad_width = 0.5
    pad_length = 0.8


    x = start_point[0][0]
    y = start_point[0][1]

    points = [(x + pad_width / 2, y + pad_width / 2), (x + pad_width / 2 + pad_length, y + pad_width / 2 - pad_length),
              (x - pad_width / 2 + pad_length, y - pad_width / 2 - pad_length), (x - pad_width / 2, y - pad_width / 2)]
    pad = core.Boundary(points)
    elist.add(pad)
    elist.add(deepcopy(pad).reflect('x'))

    ###########

    start_point = [(- rc * np.cos(pd / 2), rc * np.sin(pd / 2))]
    start_point = round_points(start_point)

    x = start_point[0][0]
    y = start_point[0][1]

    pad_width = 0.8
    pad_length1 = 2 * (y + w_ebeam / 2)
    pad_length2 = pad_length1 + pad_width
    pad_length3 = pad_length1 + 2 * pad_width

    points = [(x - w_ebeam / 2, y + w_ebeam / 2), (x - w_ebeam / 2 - pad_length1, y + w_ebeam / 2),
              (x - w_ebeam / 2 - pad_length2, pad_width / 2), (x - w_ebeam / 2 - pad_length3, pad_width / 2),
              (x - w_ebeam / 2 - pad_length3, - pad_width / 2), (x - w_ebeam / 2 - pad_length2, - pad_width / 2),
              (x - w_ebeam / 2 - pad_length1, - y - w_ebeam / 2), (x - w_ebeam / 2, - y - w_ebeam / 2)]
    pad = core.Boundary(points)
    elist.add(pad)

    return elist

########Flatronium+Junction:Mahya#######
def flatronium_junction(L, Ld, L_sq, n, a, w_ebeam):

    w = 0.05
    s = L / L_sq * w
    sd = Ld / L_sq * w
    if w_ebeam == -1:
        w_ebeam = w



    #pd = (2 * np.pi - po) / (n - 1)
    #spd = np.sin(pd / 2)

    #ra = a / (2 * spd)
    #rb = 2 * sd / (n * (1 - spd)) + ra

    #b = 4 * sd * spd / (n * (1 - spd)) + a

    alpha = 0.25
    beta = 1.0


    po = 30.0 / 180.0 * np.pi
    pd = 2 * np.pi / n
    p = po / 2 + pd * np.arange(n)
    spd = np.sin(pd / 2)
    cpd = np.cos(pd / 2)


    ra = a / (2 * spd)
    rb = calc_outer_radius(s, a, n)
    rc = (s + (ra + rb) * (1 - spd) * n / 2) / (n + (n - 2) * spd)


    start_point = [((rc + beta) * np.cos(po / 2) + alpha, (rc + beta) * np.sin(po / 2) - alpha)]
    start_point = round_points(start_point)
    points = []



    circle_version = True

    # circle version
    if circle_version:

        # inner circle 1
        p = pd * np.arange(n)
        rc_center_inner1 = ra * (cpd + spd ** 2 / cpd)
        rc_inner1 = ra * spd / cpd
        angle_inner = np.pi / 2 - pd / 2

        # inner circle 2
        p = pd * np.arange(n)
        rc_center_inner2 = rb * (cpd + spd ** 2 / cpd)
        rc_inner2 = rb * spd / cpd
        angle_inner = np.pi / 2 - pd / 2

        # outer circle
        p = pd * np.arange(n)
        rc_center_outer = rc * (cpd + spd ** 2 / cpd)
        rc_outer = rc * spd / cpd
        angle_outer = np.pi / 2 + pd / 2

        #points.append((rc_center_outer + rc_outer, 0.0))
        #c = circle(rc_center_outer, 0, rc_outer, 0.02, angle_outer)
        #points += c.tolist()

        for i in range(6, n-5):

            if i % 2 == 0:
                if i < n // 2:
                    c = circle(rc_center_inner1 * np.cos(p[i]),
                                   rc_center_inner1 * np.sin(p[i]), rc_inner1, p[i] - np.pi + angle_inner,
                                   p[i] - np.pi - angle_inner)
                    points += c.tolist()

                else:

                    c = circle(rc_center_inner2 * np.cos(p[i]),
                                   rc_center_inner2 * np.sin(p[i]), rc_inner2, p[i] - np.pi + angle_inner,
                                   p[i] - np.pi - angle_inner)
                    points += c.tolist()
            else:

                    c = circle(rc_center_outer * np.cos(p[i]),
                               rc_center_outer * np.sin(p[i]), rc_outer, p[i] - angle_outer, p[i] + angle_outer)
                    points += c.tolist()


        #c = circle(rc_center_outer, 0, rc_outer, - angle_outer, -0.02)
        #points += c.tolist()
        #points.append((rc_center_outer + rc_outer, 0.0))

    # no circle version
    if not circle_version:
        p = pd / 2 + pd * np.arange(n)
        points.append((rb * cpd, 0.0))
        for i in range(n):

            if i % 2 == 0:
                points.append((rb * np.cos(p[i]), rb * np.sin(p[i])))
                points.append((ra * np.cos(p[i]), ra * np.sin(p[i])))
            else:
                points.append((ra * np.cos(p[i]), ra * np.sin(p[i])))
                points.append((rb * np.cos(p[i]), rb * np.sin(p[i])))
        points.append((rb * cpd, 0.0))



    #points.append((start_point[0][0] - alpha, - start_point[0][1] - alpha))
    #points.append((start_point[0][0], - start_point[0][1]))

    #points.append(start_point[0])
    #points.append((start_point[0][0] - alpha, start_point[0][1] + alpha))

    print(calc_length(points))
    points = round_points(points)
    print(calc_length(points))

    superinductor = core.Path(points, w_ebeam)

    elist = core.Elements()
    elist.add(superinductor)

    # pads

    pad_width = 0.5
    pad_length = 2



    x = start_point[0][0]
    y = start_point[0][1]

    points = [(x + pad_width / 2, y + pad_width / 2), (x + pad_width / 2 + pad_length, y + pad_width / 2 - pad_length),
              (x - pad_width / 2 + pad_length , y - pad_width / 2 - pad_length), (x - pad_width / 2, y - pad_width / 2)]
    pad1 = core.Boundary(points)
    elist.add(pad1)

    points = [(x + pad_width / 2 + pad_length - pad_width , y + pad_width / 2), (x + pad_width / 2 + pad_length, y + pad_width / 2 - pad_width),
              (x - pad_width / 2 + pad_width , y - pad_width / 2 - pad_length), (x - pad_width / 2, y - pad_width / 2 - pad_length + pad_width)]
    pad2 = core.Boundary(points)
    elist.add(pad2)
    elist.add(deepcopy(pad1).reflect('x'))
    elist.add(deepcopy(pad2).reflect('x'))

    x1 = x
    y1 = y - pad_length
    pad_width1 = 0.15
    pad_length1 = 2

    points = [(x1 + pad_width1 / 2, y1 - pad_width1 / 2), (x1 + pad_width1 / 2 - pad_length1, y1 - pad_width1 / 2 - pad_length1),
              (x1 - pad_width1 / 2 - pad_length1 , y1 + pad_width1 / 2 - pad_length1), (x1 - pad_width1 / 2, y1 + pad_width1 / 2)]
    pad3 = core.Boundary(points)
    elist.add(pad3)
    elist.add(deepcopy(pad3).reflect('x'))

    ###########

    start_point = [(- rc * np.cos(pd / 2), rc * np.sin(pd / 2))]
    start_point = round_points(start_point)

    x = start_point[0][0]
    y = start_point[0][1]

    pad_width = 0.8
    pad_length1 = 2 * (y + w_ebeam / 2)
    pad_length2 = pad_length1 + pad_width
    pad_length3 = pad_length1 + 2 * pad_width

    points = [(x - w_ebeam / 2, y + w_ebeam / 2), (x - w_ebeam / 2 - pad_length1, y + w_ebeam / 2),
              (x - w_ebeam / 2 - pad_length2, pad_width / 2), (x - w_ebeam / 2 - pad_length3, pad_width / 2),
              (x - w_ebeam / 2 - pad_length3, - pad_width / 2), (x - w_ebeam / 2 - pad_length2, - pad_width / 2),
              (x - w_ebeam / 2 - pad_length1, - y - w_ebeam / 2), (x - w_ebeam / 2, - y - w_ebeam / 2)]
    pad = core.Boundary(points)
    elist.add(pad)

    return elist


def flatronium_junction_backup(L, Ld, L_sq, n, a, w_ebeam):

    w = 0.05
    Ld = 300
    s = L / L_sq * w
    sd = Ld / L_sq * w

    po = 30.0 / 180.0 * np.pi
    pd = (2 * np.pi - po) / (n -1)
    spd = np.sin(pd / 2)

    ra = a / (2 * spd)
    rb = 2 * sd / (n * (1 - spd)) + ra

    b = 4 * sd * spd / (n * (1 - spd)) + a

    rc = (s + (ra + rb) * (1 - spd) * n / 2) / (n + (n - 2) * spd)

    # circle
    alpha = 0.25
    beta = 1.0

    start_point = [((rc + beta) * np.cos(po / 2) + alpha, (rc + beta) * np.sin(po / 2) - alpha)]
    start_point = round_points(start_point)

    points = []
    points.append(start_point[0])
    points.append((start_point[0][0] - alpha, start_point[0][1] + alpha))

    p = po / 2 + pd * np.arange(n)

    for i in range(0, n - 1):

        if i % 2 == 0:
            if i < n // 2:
                points.append((ra * np.cos(p[i]), ra * np.sin(p[i])))
                points.append((ra * np.cos(p[i + 1]), ra * np.sin(p[i + 1])))
            else:
                points.append((rb * np.cos(p[i]), rb * np.sin(p[i])))
                points.append((rb * np.cos(p[i + 1]), rb * np.sin(p[i + 1])))
        else:
            points.append((rc * np.cos(p[i]), rc * np.sin(p[i])))
            points.append((rc * np.cos(p[i + 1]), rc * np.sin(p[i + 1])))

    points.append((start_point[0][0] - alpha, - start_point[0][1] - alpha))
    points.append((start_point[0][0], - start_point[0][1]))

    print(calc_length(points))
    points = round_points(points)
    print(calc_length(points))

    superinductor = core.Path(points, w_ebeam)

    elist = core.Elements()
    elist.add(superinductor)

    # pads

    pad_width = 0.5
    pad_length = 2



    x = start_point[0][0]
    y = start_point[0][1]

    points = [(x + pad_width / 2, y + pad_width / 2), (x + pad_width / 2 + pad_length, y + pad_width / 2 - pad_length),
              (x - pad_width / 2 + pad_length , y - pad_width / 2 - pad_length), (x - pad_width / 2, y - pad_width / 2)]
    pad1 = core.Boundary(points)
    elist.add(pad1)

    points = [(x + pad_width / 2 + pad_length - pad_width , y + pad_width / 2), (x + pad_width / 2 + pad_length, y + pad_width / 2 - pad_width),
              (x - pad_width / 2 + pad_width , y - pad_width / 2 - pad_length), (x - pad_width / 2, y - pad_width / 2 - pad_length + pad_width)]
    pad2 = core.Boundary(points)
    elist.add(pad2)
    elist.add(deepcopy(pad1).reflect('x'))
    elist.add(deepcopy(pad2).reflect('x'))

    x1 = x
    y1 = y - pad_length
    pad_width1 = 0.15
    pad_length1 = 2

    points = [(x1 + pad_width1 / 2, y1 - pad_width1 / 2), (x1 + pad_width1 / 2 - pad_length1, y1 - pad_width1 / 2 - pad_length1),
              (x1 - pad_width1 / 2 - pad_length1 , y1 + pad_width1 / 2 - pad_length1), (x1 - pad_width1 / 2, y1 + pad_width1 / 2)]
    pad3 = core.Boundary(points)
    elist.add(pad3)
    elist.add(deepcopy(pad3).reflect('x'))

    ###########

    start_point = [(- rc * np.cos(pd / 2), rc * np.sin(pd / 2))]
    start_point = round_points(start_point)

    x = start_point[0][0]
    y = start_point[0][1]

    pad_width = 0.8
    pad_length1 = 2 * (y + w_ebeam / 2)
    pad_length2 = pad_length1 + pad_width
    pad_length3 = pad_length1 + 2 * pad_width

    points = [(x - w_ebeam / 2, y + w_ebeam / 2), (x - w_ebeam / 2 - pad_length1, y + w_ebeam / 2),
              (x - w_ebeam / 2 - pad_length2, pad_width / 2), (x - w_ebeam / 2 - pad_length3, pad_width / 2),
              (x - w_ebeam / 2 - pad_length3, - pad_width / 2), (x - w_ebeam / 2 - pad_length2, - pad_width / 2),
              (x - w_ebeam / 2 - pad_length1, - y - w_ebeam / 2), (x - w_ebeam / 2, - y - w_ebeam / 2)]
    pad = core.Boundary(points)
    elist.add(pad)

    return elist


def round_points(points):

    for i, point in enumerate(points):

        points[i] = (np.round(point[0], 3), np.round(point[1], 3))

    return points

def calc_length(points):

    L = 0

    for i in range(len(points) - 1):

        L += np.sqrt((points[i][0] - points[i + 1][0]) ** 2 + (points[i][1] - points[i + 1][1]) ** 2)

    return L


def calc_outer_radius(s, a, n):
    pd = 2 * np.pi / n
    spd = np.sin(pd / 2)
    rb = (s + n * a * (1 / (2 * spd) - 0.5)) / (n * (1 + spd))

    return rb


def create_resolution_test(w_ebeam):

    elist = core.Elements()

    line = ((0.0, -3.0), (0.0, 3.0))
    line = core.Path(line, w_ebeam)

    N = 5
    distances = np.linspace(2 * w_ebeam, N * w_ebeam, 2 * (N - 2) + 1)

    for k in range(distances.shape[0]):
        for i in range(51):

            elist.add(deepcopy(line).translate(((k - (N - 2)) * 25 + (i - 25) * distances[k], 0)))


    return elist
