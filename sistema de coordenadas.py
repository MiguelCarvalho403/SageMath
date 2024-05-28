def coordenadas(ax_size = 2, plane = False, opacity = 0.1):
    u, v = var('y, z')
    u_ = (u, -ax_size, ax_size) # dimensão dos eixos / planos
    v_ = (v, -ax_size, ax_size) # dimensão dos eixos / planos

    x = parametric_plot((u, 0, 0), u_,  opacity=opacity, color='red')
    y = parametric_plot((0, u, 0), u_,  opacity=opacity, color='yellow')
    z = parametric_plot((0, 0, u), u_,  opacity=opacity) # Blue
    
    if plane == True:
        
        #xz = implicit_plot3d(x == 0, (-2, 2), (-2,2), (-2,2), opacity=0.3)
        yz = parametric_plot((u, v, 0), u_, v_, opacity=opacity, color='red')
        xy = parametric_plot((u, 0, v), u_, v_, opacity=opacity, color='yellow')
        xz = parametric_plot((0, u, v), u_, v_, opacity=opacity)

        return xz +xy + yz + x + y + z # planos

    return x + y + z #eixos
