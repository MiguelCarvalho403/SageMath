def coordenadas(ax_size = (-2,2), plane = False, opacity = 0.1):
    u, v = var('y, z')
    size = (ax_size[1] - ax_size[0])
    u_ = (u, ax_size[0], ax_size[1]) # dimensão dos eixos / planos
    v_ = (v, ax_size[0], ax_size[1]) # dimensão dos eixos / planos

    x = parametric_plot((u, 0, 0), u_,  opacity=1, color='red')
    x += vector([ax_size[1],0,0]).plot(color = 'red')
    labels = text3d('x',1.1*vector([ax_size[1],0,0]),fontsize= str(50*size)+'%')
    
    y = parametric_plot((0, u, 0), u_,  opacity=1, color='yellow')
    y += vector([0,ax_size[1],0]).plot(color = 'yellow')
    labels += text3d('y',1.1*vector([0,ax_size[1],0]),fontsize= str(50*size)+'%')
    
    z = parametric_plot((0, 0, u), u_,  opacity=1) # Blue
    z += vector([0,0,ax_size[1]]).plot(color = 'blue')
    labels += text3d('z',1.1*vector([0,0,ax_size[1]]),fontsize= str(50*size)+'%')
    
    if plane == True:
        
        cont = ax_size[0]
        xy = parametric_plot((u, cont, 0), u_, opacity=opacity, color='black')
        xy += parametric_plot((cont, u, 0), u_, opacity=opacity, color='black')
        while cont != ax_size[1]:

            cont += size/10
            xy += parametric_plot((u, cont, 0), u_, opacity=opacity, color='black')
            xy += parametric_plot((cont, u, 0), u_, opacity=opacity, color='black')
        return xy + x + y + z + labels # planos

    return x + y + z + labels #eixos
    
