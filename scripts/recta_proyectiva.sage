#\s
# Este archivo proporciona diferentes funciones sobre razones dobles y homografías de la recta.
#
# Autor: Pablo Sanz Sanz
#

load("espacios.sage")
load("procedimientos.sage")

# Funciones globales

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
    ocultar_procedimiento()
    assert subespacio(p0, p1, p2, p3).dim() == 1, "Los puntos deben estar alineados"
    reanudar_procedimiento()
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
        ocultar_procedimiento()
        fin = q0 != 0 and q1 != 0 and q2 != 0 and q3 != 0 and subespacio(q0, q1, q2, q3).dim() == 1
        reanudar_procedimiento()
        j = j + 1
        # Creo que ni debería cumplirse esta condición
        if j == len(p0):
            i = i + 1
            j = i + 1
    if len(p0) > 2:
        paso("Como la razón doble es invariante por proyecciones podemos usar los siguientes vectores:")
        paso(q0, ", ", q1, ", ", q2, ", ", q3)
    m0 = matrix([[q0[0], q2[0]], [q0[1], q2[1]]])
    m1 = matrix([[q1[0], q2[0]], [q1[1], q2[1]]])
    m2 = matrix([[q0[0], q3[0]], [q0[1], q3[1]]])
    m3 = matrix([[q1[0], q3[0]], [q1[1], q3[1]]])
    det0 = m0.det()
    det1 = m1.det()
    det2 = m2.det()
    det3 = m3.det()
    paso("Calculamos los determinantes:")
    paso("det", m0, " = ", det0)
    paso("det", m1, " = ", det1)
    paso("det", m2, " = ", det2)
    paso("det", m3, " = ", det3)
    paso("Calculamos la razón doble: {P0, P1; P2, P3} = (", det0, "/", det1, ") : (", det2, "/", det3, ")")
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
# theta0: complejo/Infinity - primer punto de la razón doble \\
# theta1: complejo/Infinity - segundo punto de la razón doble \\
# theta2: complejo/Infinity - tercer punto de la razón doble \\
# theta3: complejo/Infinity - cuarto punto de la razón doble
#
def razon_doble(theta0, theta1, theta2, theta3):
    p0 = vector([1, 0]) if theta0 == Infinity else vector([theta0, 1])
    p1 = vector([1, 0]) if theta1 == Infinity else vector([theta1, 1])
    p2 = vector([1, 0]) if theta2 == Infinity else vector([theta2, 1])
    p3 = vector([1, 0]) if theta3 == Infinity else vector([theta3, 1])
    paso("Homogeneizando las coordenadas quedan los vectores:")
    paso(p0, ", ", p1, ", ", p2, ", ", p3)
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
# theta0: complejo/Infinity - primer punto de la razón doble \\
# theta1: complejo/Infinity - segundo punto de la razón doble \\
# theta2: complejo/Infinity - tercer punto de la razón doble \\
# theta3: complejo/Infinity - cuarto punto de la razón doble
#
def es_cuaterna_armonica(theta0, theta1, theta2, theta3):
    return razon_doble(theta0, theta1, theta2, theta3) == -1

#\f
# Devuelve el conjugado armónico de p2 respecto de p0 y p1.
# 
# Implementación \\
# Como la razón doble es la coordenada del cuarto punto respecto a la referencia {p0, p1; p2}, obtiene la recta que
# forman p0, p1 y p2 y devuelve el punto correspondiente a la coordenada theta == -1.
#
# Parámetros \\
# p0: vector(n>1) - primer punto de la razón doble \\
# p1: vector(n>1) - segundo punto de la razón doble \\
# p2: vector(n>1) - tercer punto de la razón doble
#
def conjugado_armonico_puntos(p0, p1, p2):
    paso("Montamos la recta con la referencia {P0=", p0, ", P1=", p1, "; P2=", p2, "}")
    paso("Así, P0: theta = Infinity, P1: theta = 0, P2: theta = 1 y sustituimos en -1")
    return recta_proyectiva(p0, p1, p2)[-1]
    
#\f
# Devuelve el conjugado armónico de theta2 respecto de theta0 y theta1.
# 
# Implementación \\
# Resuelve la ecuación {theta0, theta1, theta2, theta} == -1 en theta utilizando las funciones anteriores.
#
# Parámetros \\
# theta0: complejo/Infinity - primer punto de la razón doble \\
# theta1: complejo/Infinity - segundo punto de la razón doble \\
# theta2: complejo/Infinity - tercer punto de la razón doble
#
def conjugado_armonico(theta0, theta1, theta2):
    theta = var('theta')
    paso("Resolvemos {", theta0, ", ", theta1, ", ", theta2, ", ", theta, "} = -1")
    res = solve([es_cuaterna_armonica(theta0, theta1, theta2, theta)], theta)
    # No ha habido solución
    if len(res) == 0:
        return Infinity
    return res[0].rhs()
    
# Clases

#\c
# Clase que representa una recta proyectiva para cualquier dimension ambiente. Sirve para cálculos de razones dobles
# o para homografías entre rectas.
#
class recta_proyectiva:
    
    #\i
    # Inicializa la recta dada la referencia {p0, p1; p2} que se quiere utilizar.
    #
    # Implementación \\
    # Obtiene la matriz asociada a la referencia dada (esto es, multiplicando el vector asociado a p0 y p1 por un coeficiente
    # de forma que a*p0 + b*p1 == p2) y construye el punto b*p1 + theta*a*p0.
    #
    # Parámetros \\
    # p0: vector(n) - punto cuya coordenada será theta == Infinity
    # p1: vector(n) - punto cuya coordenada será theta == 0
    # p2: vector(n) - punto cuya coordenada será theta == 1
    #
    def __init__(self, p0, p1, p2):
        ocultar_procedimiento()
        self._subespacio = subespacio(p0, p1, p2)
        reanudar_procedimiento()
        assert self._subespacio.dim() == 1, "Los puntos deben estar alineados"
        # Forzamos p2: theta == 1
        var('a b')
        paso("Forzamos P2: theta = 1:")
        paso(a, "*", p0, " + ", b, "*", p1, " = ", a * p0 + b * p1, " = ", p2)
        coef = solve((a * p0 + b * p1 - p2).list(), a, b)
        paso(coef[0])
        self._theta = var('theta', latex_name = '\theta')
        self._infinito = p0
        self._valores = (b * p1 + self._theta * a * p0).substitute(coef[0])
        paso("Sustituimos y combinamos con theta: ", self._valores)
        
    # Métodos accedentes
    
    #\m
    # Devuelve el punto con coordenada theta.
    #
    # Uso: r[theta_0] (theta_0 es la coordenada del punto y r es una recta(_proyectiva)).
    #
    # Parámetros \\
    # theta: complejo/Infinity - coordenada del punto de la recta que se quiere obtener
    #
    def __getitem__(self, theta):
        if theta == Infinity:
            return self._infinito
        return self._valores.substitute(self._theta == theta)
    
    #\m
    # Devuelve esta recta proyectiva como objeto del tipo subespacio.
    def subespacio(self):
        return self._subespacio;
    
    # Otros métodos
    
    #\m
    # Devuelve la coordenada no homogénea asociada al punto p dado.
    #
    # Implementación \\
    # Resuelve la ecuación recta(theta) == alpha * p para theta.
    #
    # Parámetros \\
    # p: vector(n) - punto de la recta del que se quiere conocer su coordenada
    #
    def coordenada(self, p):
        assert p in self._subespacio, "El punto debe pertenecer a la recta"
        var('alpha')
        inf = len(solve([self._infinito[i] == alpha * p[i] for i in range(len(p))], alpha)) > 0
        if inf:
            return Infinity
        return solve([self._valores[i] == alpha * p[i] for i in range(len(p))], alpha, self._theta)[0][1].rhs()
    
    def __repr__(self):
        return "<Recta proyectiva " + str(self._valores) + ">"
        
#\c
# Clase que representa una homografía de una recta en sí misma. En general, se crearán usando funciones auxiliares
# y no con su constructor.
#
class homografia_recta:
    
    #\i
    # Construye una homografía sobre la recta dada cuya matriz sea la dada. En general, este constructor no se
    # usará directamente.
    #
    # Parámetros \\
    # matriz: matriz(2, 2) - matriz asociada a la homografía \\
    # recta: recta sobre la que actúa la homografía (por defecto una recta con referencia {Infinity, 0; 1})
    #
    def __init__(self, matriz, recta = recta_proyectiva(vector([1, 0]), vector([0, 1]), vector([1, 1]))):
        assert matriz.nrows() == 2 and matriz.ncols() == 2, "La matriz debe ser 2x2"
        assert matriz.det() != 0, "Para que sea una homografía el rango de su matriz asociada debe ser máximo"
        self._matriz = matriz
        self._recta = recta
        self._theta1 = var('theta', latex_name = '\theta')
        self._theta2 = var('theta2', latex_name = '\theta_2')
        self._expresion_mobius = theta2 == (matriz[0][0] * theta + matriz[0][1]) / (matriz[1][0] * theta + matriz[1][1])
        # Multiplicamos por el denominador y restamos el numerador
        self._expresion = (theta2 * self._expresion_mobius.rhs().denominator() - self._expresion_mobius.rhs().numerator()).expand() == 0
        
    # Métodos accedentes
    
    #\m
    # Devuelve la matriz asociada a esta homografía.
    def matriz_asociada(self):
        return self._matriz
        
    #\m
    # Devuelve la recta sobre la que se aplica esta homografía.
    def recta(self):
        return self._recta
    
    #\m
    # Devuelve la homografía expresada como una transformación de M
    def expresion_mobius(self):
        return self._expresion_mobius
    
    #\m
    # Devuelve la ecuación en theta y theta' de esta homografía.
    def ecuacion(self):
        return self._expresion
    
    #\m
    # Calcula la imagen mediante esta homografía del punto dado.
    #
    # Uso: r(x) ó r(v) (donde r es una recta_proyectiva, x un complejo/Infinity y v un vector).
    #
    # Implementación \\
    # Dado un parámetro no homogéneo sustituye en la expresión de Möbius y devuelve el resultado. \\
    # Dado un vector del espacio, primero calcula su coordenada en la recta y luego vuelve a llamar a este método.
    #
    # Parámetros \\
    # x: complejo/Infinity/vector(n) - punto del que se quiere calcular su imagen
    #
    def __call__(self, x):
        # Estamos recibiendo un parámetro no homogéneo, simplemente sustituimos
        if x == Infinity or len(x.list()) == 1:
            return self.__sustituir(x)
        # Asumimos que es un vector, si es del mismo espacio calculamos su coordenada en la recta primero
        assert len(x) == self._recta.subespacio().dimension_ambiente() + 1, \
                "El punto debe pertenecer al mismo espacio ambiente"
        assert x in self._recta.subespacio(), "El punto debe pertenecer a la recta"
        coord = self._recta.coordenada(x)
        paso("La coordenada de ", x, " en la recta es: ", coord, ". Sustituimos en la expresión y recuperamos el punto")
        return self._recta[self(coord)]
            
    # Métodos auxiliares
            
    def __sustituir(self, x):
        # Para infinito se pierden las componentes sin theta
        if x == Infinity:
            return self._matriz[0][0] / self._matriz[1][0]
        # En cualquier otro caso se sustituye en la expresión de Möbius
        return self._expresion_mobius.substitute(self._theta1 == x).rhs()
        
    def __repr__(self):
        return "<Homografía " + str(self._expresion) + " de " + str(self._recta) + ">"
