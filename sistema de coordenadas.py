def coordenadas():
    u, v = var('y, z')
    u_ = (u, -2, 2)
    v_ = (v, -2, 2)

    #xz = implicit_plot3d(x == 0, (-2, 2), (-2,2), (-2,2), opacity=0.3)
    yz = parametric_plot((u, v, 0), u_, v_, opacity = 0.1, color='red')
    xy = parametric_plot((u, 0, v), u_, v_, opacity=0.1, color='yellow')
    xz = parametric_plot((0, u, v), u_, v_, opacity=0.1)

    x = parametric_plot((u, 0, 0), u_,  opacity = 1, color='red')
    y = parametric_plot((0, u, 0), u_,  opacity=1, color='yellow')
    z = parametric_plot((0, 0, u), u_,  opacity=1)

    return x + y + z #eixos
    #return xz +xy + yz # planos
