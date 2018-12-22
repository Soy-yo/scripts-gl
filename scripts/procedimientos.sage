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

_mostrar_procedimiento = False
_forzar_ocultar = False

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
    if _mostrar_procedimiento and not _forzar_ocultar:
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
    global _forzar_ocultar
    _mostrar_procedimiento = b
    _forzar_ocultar = False
      
#\f
# Desactiva temporalmente que se muestren los pasos que ejecutan las funciones. Sirve para desarrolladores, si llaman a
# funciones auxiliares que generarán un procedimiento, pero que no es relevante para el cálculo que se está efectuando.
#
def ocultar_procedimiento():
    global _forzar_ocultar
    _forzar_ocultar = True
    
#\f
# Activa de nuevo que se muestren los pasos que ejecutan las funciones si estaban activados antes de llamar a
# ocultar_procedimiento().
#
def reanudar_procedimiento():
    global _forzar_ocultar
    _forzar_ocultar = False
