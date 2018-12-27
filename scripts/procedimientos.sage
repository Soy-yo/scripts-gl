#\s
# Módulo para el control de la impresión de los procedimientos por pantalla.
#
# Para usuarios: activar los procedimientos con mostrar_procedimiento(True) si se quieren ver las operaciones que hacen
# las funciones por dentro (si el desarrollador ha indicado algo). Para dejar de mostrar, mostrar_procedimiento(False).
#
# Para desarrolladores: usar paso(texto) cuando la función haga algún tipo de operación que se consifere relevante.
#
# Autor: Pablo Sanz Sanz
#

_mostrar_procedimiento = [False]

#\f
# Función para desarrolladores. Muestra por pantalla el texto especificado si la variable _mostrar_procedimiento está
# activada (no hace nada en caso contrario). Se utiliza para mostrar los pasos que realiza una función.
#
# Por ejemplo, \\
# paso("Calculamos el determinante de la matriz M:") \\
# paso("det", M, "=", M.det())
#
# Parámetros \\
# *objetos: cualquiera - objetos a mostrar
#
def paso(*objetos):
    if _mostrar_procedimiento[-1]:
        s = ""
        for obj in objetos:
            s = s + str(latex(obj))
        show(LatexExpr(s))

#\f
# Activa o desactiva que se muestren los pasos que ejecutan las funciones.
#
# Parámetros \\
# b: booleano - True para activarlos, False para desactivarlos (True si no se especifica nada)
#
def mostrar_procedimiento(b = True):
    global _mostrar_procedimiento
    _mostrar_procedimiento = [b]

#\f
# Para desarroladores. Desactiva temporalmente los procedimientos. Se usa cuando vamos a llamar a otras funciones que
# sabemos que muestran procedimientos, pero no queremos que aparezcan por pantalla porque son cálculos poco relevantes
# para el cometido de nuestra función (por ejemplo, si queremos comprobar que realmente el punto que nos han pasado
# pertenece a un subespacio, pero no queremos que se muestre cómo lo hace porque el usuario da por hecho que su entrada
# es correcta).
#
# Uso \\
# _no_pasos() \\
# funciones que no queremos que muestren sus procedimientos \\
# _no_pasos(False) # IMPORTANTE NO OLVIDAR DESACTIVARLO AL FINAL O CAUSARÁ PROBLEMAS
#
# Parámetros \\
# b: booleano - True para desactivarlos, False para recuperarlos de nuevo (por defecto True)
def _no_pasos(b = True):
    global _mostrar_procedimiento
    if b:
        _mostrar_procedimiento.append(False)
    else:
        _mostrar_procedimiento.pop()
