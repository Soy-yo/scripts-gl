load("espacios.sage")
#\s
# Este archivo proporciona diferentes funciones sobre razones dobles y homografías de la recta.
#
# Autor: Pablo Sanz Sanz
#

#\f
# Calcula la razón doble ({P0, P1, P2, P3}) de los puntos alineados dados.
#
# Implementación \\
# Utiliza que la razón doble es invariante por proyecciones y la fórmula de la 
# razón doble de los determinantes. Se elige una proyección que siga siendo una recta.
#
# Parámetros \\
# p0: vector(n>1) - primer punto de la razón doble \\
# p1: vector(n>1) - segundo punto de la razón doble \\
# p2: vector(n>1) - tercer punto de la razón doble \\
# p3: vector(n>1) - cuarto punto de la razón doble
#
def razon_doble_puntos(p0, p1, p2, p3):
    assert len(p0) > 1, "Para coordenadas de la recta usar 'razon_doble'"
    assert subespacio(p0, p1, p2, p3).dim() == 1, "Los puntos deben estar alineados"
    fin = len(p0) == 2
    i = 0
    j = 1
    q0 = p0
    q1 = p1
    q2 = p2
    q3 = p3
    while not fin:
        # Proyectamos en las coordenadas i y j
        q0 = vector([p0[i], p0[j]])
        q1 = vector([p1[i], p1[j]])
        q2 = vector([p2[i], p2[j]])
        q3 = vector([p3[i], p3[j]])
        # Si nos sigue quedando una recta (y no un punto) ya podemos calcular la razón doble
        fin = subespacio(q0, q1, q2, q3).dim() == 1
        j = j + 1
        # Creo que ni debería cumplirse esta condición
        if j == len(p0):
            i = i + 1
            j = i + 1
    det0 = matrix([[q0[0], q2[0]], [q0[1], q2[1]]]).det()
    det1 = matrix([[q1[0], q2[0]], [q1[1], q2[1]]]).det()
    det2 = matrix([[q0[0], q3[0]], [q0[1], q3[1]]]).det()
    det3 = matrix([[q1[0], q3[0]], [q1[1], q3[1]]]).det()
    if det1 == 0 or det2 == 0:
        return Infinity
    return (det0 / det1) * (det3 / det2)
    
#\f
# Calcula la razón doble ({theta0, theta1, theta2, theta3}) de los puntos de la recta dados
# (con parámetro no homogéneo).
#
# Implementación \\
# Utiliza la razón doble definida para puntos tras homogeneizar las coordenadas.
#
# Parámetros \\
# theta0: complejo - primer punto de la razón doble \\
# theta1: complejo - segundo punto de la razón doble \\
# theta2: complejo - tercer punto de la razón doble \\
# theta3: complejo - cuarto punto de la razón doble
#
def razon_doble(theta0, theta1, theta2, theta3):
    p0 = vector([1, 0]) if theta0 == Infinity else vector([theta0, 1])
    p1 = vector([1, 0]) if theta1 == Infinity else vector([theta1, 1])
    p2 = vector([1, 0]) if theta2 == Infinity else vector([theta2, 1])
    p3 = vector([1, 0]) if theta3 == Infinity else vector([theta3, 1])
    return razon_doble_puntos(p0, p1, p2, p3)

#\f
# Determina si los puntos dados forman una cuaterna armónica, esto es, {P0, P1, P2, P3} = -1.
# 
# Implementación \\
# Devuelve {P0, P1, P2, P3} == -1 utilizando las funciones anteriores.
#
# Parámetros \\
# p0: vector(n>1) - primer punto de la razón doble \\
# p1: vector(n>1) - segundo punto de la razón doble \\
# p2: vector(n>1) - tercer punto de la razón doble \\
# p3: vector(n>1) - cuarto punto de la razón doble
#
def es_cuaterna_armonica_puntos(p0, p1, p2, p3):
    return razon_doble_puntos(p0, p1, p2, p3) == -1
        
#\f
# Determina si los puntos dados forman una cuaterna armónica, esto es,
# {theta0, theta1, theta2, theta3} = -1.
# 
# Implementación \\
# Devuelve {theta0, theta1, theta2, theta3} == -1 utilizando las funciones anteriores.
#
# Parámetros \\
# theta0: complejo - primer punto de la razón doble \\
# theta1: complejo - segundo punto de la razón doble \\
# theta2: complejo - tercer punto de la razón doble \\
# theta3: complejo - cuarto punto de la razón doble
#
def es_cuaterna_armonica(theta0, theta1, theta2, theta3):
    return razon_doble(theta0, theta1, theta2, theta3) == -1
    
#\f
# Devuelve el conjugado armónico de theta2 respecto de theta0 y theta1.
# 
# Implementación \\
# Resuelve la ecuación {theta0, theta1, theta2, theta} == -1 en theta utilizando las funciones anteriores.
#
# Parámetros \\
# theta0: complejo - primer punto de la razón doble \\
# theta1: complejo - segundo punto de la razón doble \\
# theta2: complejo - tercer punto de la razón doble \\
# theta: variable - variable que se quiere que sea el conjugado armónico
#
def conjugado_armonico(theta0, theta1, theta2, theta):
    return solve([es_cuaterna_armonica(theta0, theta1, theta2, theta)], theta)
    
    
    
    
    
    
    
