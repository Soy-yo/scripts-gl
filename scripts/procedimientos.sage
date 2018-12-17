#\s
# Módulo para el control de la impresión de los procedimientos por pantalla.
#
# Para usuarios: activar los procedimientos con mostrar_procedimiento(True) si se quieren ver las operaciones que hacen
# las funciones por dentro (si el desarrollador ha indicado algo). Para dejar de mostrar, mostrar_procedimiento(False).
#
# Para desarrolladores: usar pasos(texto) cuando la función haga algún tipo de operación que se consifere relevante.
#

_mostrar_procedimiento = False

#\f
# Función para desarrolladores. Muestra por pantalla el texto especificado si la variable _mostrar_procedimiento está
# activada (no hace nada en caso contrario). Se utiliza para mostrar los pasos que realiza una función.
#
# Por ejemplo, \\
# paso("Calculamos el determinante de la matriz M:") \\
# paso("det" + str(M) + " = " + str(M.det()))
#
# Parámetros \\
# texto: string - texto a mostrar
#
def paso(texto):
    if _mostrar_procedimiento:
        print texto
        
#\f
# Activa o desactiva que se muestren los pasos que ejecutan las funciones.
#
# Parámetros \\
# b: booleano - True para activarlos, False para desactivarlos
#
def mostrar_procedimiento(b):
    global _mostrar_procedimiento
    _mostrar_procedimiento = b