def coordenadas(x_size = (-2,2), y_size = (-2,2), plane = True, opacity = 0.7, grad_x = 1, grad_y = 1, **kwargs):
    # grad: graduação a qual a malha é subdividida, por padrão 1 unidade
    # x_size: intervalo que a malha se extenderá em x
    # y_size: intervalo que a malha se extenderá em y
    
    u, v = var('y, z') # variaveis simbólicas
    interval_x = (x_size[1] - x_size[0]) # Tamanho do intervalo x
    interval_y = (y_size[1] - y_size[0]) # Tamanho do intervalo y
    u_ = (u, x_size[0], x_size[1]) # dimensão dos eixos / planos
    v_ = (v, y_size[0], y_size[1]) # dimensão dos eixos / planos

    x = parametric_plot((u, 0, 0), u_,  opacity=1, color='red') # vetor de referencia em x
    x += vector([x_size[1],0,0]).plot(color = 'red')
    labels = text3d('x',1.1*vector([x_size[1],0,0]),fontsize= str(250)+'%')
    
    y = parametric_plot((0, v, 0), v_,  opacity=1, color='yellow') # vetor de referencia em y
    y += vector([0,y_size[1],0]).plot(color = 'yellow')
    labels += text3d('y',1.1*vector([0,y_size[1],0]),fontsize= str(250)+'%')
    
    z = parametric_plot((0, 0, u), u_,  opacity=1, frame = False) # vetor de referencia em z
    z += vector([0,0,x_size[1]]).plot(color = 'blue')
    labels += text3d('z',1.1*vector([0,0,x_size[1]]),fontsize= str(250)+'%')

    razao_x = interval_x / grad_x
    razao_y = interval_y / grad_y
    
    if plane == True:
        
        cont_x = x_size[0]
        cont_y = y_size[0]
        
        malha_xy = parametric_plot((u, cont_y, 0), u_, opacity=opacity, color='black', **kwargs)
        malha_xy += parametric_plot((cont_x, v, 0), v_, opacity=opacity, color='black')
        
        while cont_x != x_size[1]:

            cont_x += interval_x/razao_x
            malha_xy += parametric_plot((cont_x, v, 0), v_, opacity=opacity, color='black')
            
        while cont_y != y_size[1]:

            cont_y += interval_y/razao_y
            malha_xy += parametric_plot((u, cont_y, 0), u_, opacity=opacity, color='black')

            
        return malha_xy + x + y + z + labels # planos

    return x + y + z + labels #eixos
