#\s
# Este archivo proporciona diferentes funciones y clases sobre razones dobles y homografías de la recta.
#
# Suponemos que los puntos expresados con coordenada inhomogénea theta se obtiene mediante x/y.
#
# Autor: Pablo Sanz Sanz
#

# Funciones globales

#\f
# Función auxiliar que determina si el parámetro recibido es un punto (tipo vector) o un parámetro
# inhomogéneo de una recta proyectiva.
#
# Parámetros \\
# x: complejo/Infinity/variable/vector(n>1) - elemento a considerar
#
def es_parametro(x):
    # Infinity no tiene método .list(), se ve apartte
    # La lista tendrá un solo elemento a menos que sea una variable (que la trata como un polinomio)
    return x == Infinity or x == -Infinity or len(x.list()) == 1 or type(x) == sage.symbolic.expression.Expression
    
#\f
# Calcula la razón doble {A, B; C, D}. Es válido tanto para puntos como para parámetros inhomogéneos
# de una recta proyectiva (ver razon_doble_puntos y razon_doble_theta, respectivamente).
#
# Parámetros \\
# a: complejo/Infinity/vector(n>1) - primer punto de la razón doble \\
# b: complejo/Infinity/vector(n>1) - segundo punto de la razón doble \\
# c: complejo/Infinity/vector(n>1) - tercer punto de la razón doble \\
# d: complejo/Infinity/vector(n>1) - cuarto punto de la razón doble
#
def razon_doble(a, b, c, d):
    if es_parametro(a):
        assert es_parametro(b) and es_parametro(c) and es_parametro(d), \
                "Todos los puntos deben ser del mismo tipo"
        return razon_doble_theta(a, b, c, d)
    return razon_doble_puntos(a, b, c, d)

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
    _no_pasos()
    assert subespacio(p0, p1, p2, p3).dim() == 1, "Los puntos deben estar alineados"
    _no_pasos(False)
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
        _no_pasos()
        fin = q0 != 0 and q1 != 0 and q2 != 0 and q3 != 0 and subespacio(q0, q1, q2, q3).dim() == 1
        _no_pasos(False)
        j = j + 1
        # Creo que ni debería cumplirse esta condición
        if j == len(p0):
            i = i + 1
            j = i + 1
    if len(p0) > 2:
        paso("Como la razon doble es invariante por proyecciones podemos usar los siguientes vectores:")
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
    paso("Calculamos la razon doble: {P0, P1; P2, P3} = (", det0, "/", det1, ") : (", det2, "/", det3, ")")
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
def razon_doble_theta(theta0, theta1, theta2, theta3):
    p0 = vector([1, 0]) if theta0 == Infinity or theta == -Infinity else vector([theta0, 1])
    p1 = vector([1, 0]) if theta1 == Infinity or theta == -Infinity else vector([theta1, 1])
    p2 = vector([1, 0]) if theta2 == Infinity or theta == -Infinity else vector([theta2, 1])
    p3 = vector([1, 0]) if theta3 == Infinity or theta == -Infinity else vector([theta3, 1])
    paso("Homogeneizando las coordenadas quedan los vectores:")
    paso(p0, ", ", p1, ", ", p2, ", ", p3)
    return razon_doble_puntos(p0, p1, p2, p3)

#\f
# Determina si los puntos dados forman una cuaterna armónica, esto es, {A, B; C, D} = -1. Es válido tanto
# para puntos como para parámetros inhomogéneos de una recta proyectiva (ver es_cuaterna_armonica_puntos y
# es_cuaterna_armonica_theta, respectivamente).
#
# Parámetros \\
# a: complejo/Infinity/vector(n>1) - primer punto de la razón doble \\
# b: complejo/Infinity/vector(n>1) - segundo punto de la razón doble \\
# c: complejo/Infinity/vector(n>1) - tercer punto de la razón doble \\
# d: complejo/Infinity/vector(n>1) - cuarto punto de la razón doble
#
def es_cuaterna_armonica(a, b, c, d):
    if es_parametro(a):
        assert es_parametro(b) and es_parametro(c) and es_parametro(d), \
                "Todos los puntos deben ser del mismo tipo"
        return es_cuaterna_armonica_theta(a, b, c, d)
    return es_cuaterna_armonica_puntos(a, b, c, d)
    
    
#\f
# Determina si los puntos dados forman una cuaterna armónica, esto es, {P0, P1, P2, P3} = -1.
# 
# Implementación \\
# Devuelve {P0, P1; P2, P3} == -1 utilizando las funciones anteriores.
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
# Devuelve {theta0, theta1; theta2, theta3} == -1 utilizando las funciones anteriores.
#
# Parámetros \\
# theta0: complejo/Infinity - primer punto de la razón doble \\
# theta1: complejo/Infinity - segundo punto de la razón doble \\
# theta2: complejo/Infinity - tercer punto de la razón doble \\
# theta3: complejo/Infinity - cuarto punto de la razón doble
#
def es_cuaterna_armonica_theta(theta0, theta1, theta2, theta3):
    return razon_doble(theta0, theta1, theta2, theta3) == -1

#\f
# Devuelve el conjugado armónico de c respecto de a y b. Es válido tanto para puntos como para parámetros inhomogéneos
# de una recta proyectiva (ver conjugado_armonico_puntos y conjugado_armonico_theta, respectivamente).
#
# Parámetros \\
# a: complejo/Infinity/vector(n>1) - primer punto de la razón doble \\
# b: complejo/Infinity/vector(n>1) - segundo punto de la razón doble \\
# c: complejo/Infinity/vector(n>1) - tercer punto de la razón doble \\
#
def conjugado_armonico(a, b, c):
    if es_parametro(a):
        assert es_parametro(b) and es_parametro(c), "Todos los puntos deben ser del mismo tipo"
        return conjugado_armonico_theta(a, b, c)
    return conjugado_armonico_puntos(a, b, c)

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
    paso("Asi, P0: theta = Infinity, P1: theta = 0, P2: theta = 1 y sustituimos en -1")
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
def conjugado_armonico_theta(theta0, theta1, theta2):
    theta = var('theta', latex_name = r'theta')
    paso("Resolvemos {", theta0, ", ", theta1, ", ", theta2, ", ", theta, "} = -1")
    res = solve([es_cuaterna_armonica_theta(theta0, theta1, theta2, theta)], theta)
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
    # p0: vector(n) - punto cuya coordenada será theta == Infinity \\
    # p1: vector(n) - punto cuya coordenada será theta == 0 \\
    # p2: vector(n) - punto cuya coordenada será theta == 1
    #
    def __init__(self, p0, p1, p2):
        _no_pasos()
        self._subespacio = subespacio(p0, p1, p2)
        _no_pasos(False)
        assert self._subespacio.dim() == 1, "Los puntos deben estar alineados"
        # Forzamos p2: theta == 1
        var('a b')
        paso("Forzamos P2: theta = 1:")
        paso(a, "*", p0, " + ", b, "*", p1, " = ", a * p0 + b * p1, " = ", p2)
        coef = solve((a * p0 + b * p1 - p2).list(), a, b)
        paso(coef[0])
        self._theta = var('theta', latex_name = r'theta')
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
        if theta == Infinity or theta == -Infinity:
            return self._infinito
        return self._valores.substitute(self._theta == theta)
    
    #\m
    # Devuelve la referencia de esta recta proyectiva en forma de lista [Infinity, 0, 1].
    def referencia(self):
        return [self[Infinity], self[0], self[1]]
    
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
        assert p in self, "El punto debe pertenecer a la recta"
        var('alpha')
        inf = len(solve([self._infinito[i] == alpha * p[i] for i in range(len(p))], alpha)) > 0
        if inf:
            return Infinity
        return solve([self._valores[i] == alpha * p[i] for i in range(len(p))], alpha, self._theta)[0][1].rhs()
    
    # \m
    # Operador ==. Determina si dos rectas proyectivas son iguales.
    #
    # Uso: r == s (r y s son rectas proyectivas).
    #
    # Implementación \\
    # Comrpueba que el subespacio sea el mismo.
    #
    # Parámetros \\
    # otra: recta_proyectiva - recta a comprobar la igualdad
    #
    def __eq__(self, otra):
        return self.subespacio() == otra.subespacio()
        
    # \m
    # Operador in. Determina si un punto está contenido en este subespacio o no.
    #
    # Uso: P in r (P es un punto y r una recta proyectiva).
    #
    # Implementación \\
    # Comprueba que el punto pertenezca al subespacio asociado.
    #
    # Parámetros \\
    # punto: vector(n) - punto que comprobar si pertenece a la recta
    #
    def __contains__(self, punto):
        return punto in self._subespacio
    
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
        assert matriz.det() != 0, "Para que sea una homografia el rango de su matriz asociada debe ser maximo"
        self._aplicacion = aplicacion_proyectiva(matriz)
        self._matriz = matriz
        self._recta = recta
        self._theta1 = var('theta', latex_name = r'theta')
        self._theta2 = var('theta2', latex_name = r'theta_2')
        self._expresion_mobius = theta2 == (matriz[0][0] * theta + matriz[0][1]) / (matriz[1][0] * theta + matriz[1][1])
        # Multiplicamos por el denominador y restamos el numerador
        self._expresion = (theta2 * self._expresion_mobius.rhs().denominator() - self._expresion_mobius.rhs().numerator()).expand() == 0
        
    # Métodos accedentes
    
    #\m
    # Devuelve esta homografía como una aplicacion_proyectiva.
    def aplicacion(self):
        return self._aplicacion
    
    #\m
    # Devuelve la matriz asociada a esta homografía.
    def matriz_asociada(self):
        return self._matriz
        
    #\m
    # Devuelve la recta sobre la que se aplica esta homografía.
    def recta(self):
        return self._recta
        
    #\m
    # Devuelve la referencia de la recta en la que está expresada esta homografía.
    def referencia(self):
        return self._recta.referencia()
    
    #\m
    # Devuelve la homografía expresada como una transformación de M
    def expresion_mobius(self):
        return self._expresion_mobius
    
    #\m
    # Devuelve la ecuación en theta y theta' de esta homografía.
    def ecuacion(self):
        return self._expresion
        
    # Otros métodos
    
    #\m
    # Calcula la imagen mediante esta homografía del punto dado.
    #
    # Uso: h(x) ó h(p) (donde r es una homografia_recta, x un complejo/Infinity y p un punto).
    #
    # Implementación \\
    # Dado un parámetro no homogéneo sustituye en la expresión de Möbius y devuelve el resultado. \\
    # Dado un punto del espacio, primero calcula su coordenada en la recta y luego vuelve a llamar a este método.
    #
    # Parámetros \\
    # x: complejo/Infinity/vector(n) - punto del que se quiere calcular su imagen
    #
    def __call__(self, x):
        # Estamos recibiendo un parámetro no homogéneo, simplemente sustituimos
        if es_parametro(x):
            return self.__sustituir(x)
        # Asumimos que es un vector, si es del mismo espacio calculamos su coordenada en la recta primero
        assert x in self._recta.subespacio(), "El punto debe pertenecer a la recta"
        coord = self._recta.coordenada(x)
        paso("La coordenada de ", x, " en la recta es: ", coord, ". Sustituimos en la expresion y recuperamos el punto")
        return self._recta[self(coord)]
    
    #\m
    # Devuelve los autovalores de esta homografía.
    #
    # Implementación \\
    # Ver autovalores() de aplicacion_proyectiva (aplicacion_proyectiva.sage).
    #
    def autovalores(self):
        return self._aplicacion.autovalores()
    
    #\m
    # Devuelve los parámetros de los puntos fijos de esta homografía. Devuelve una lista vacía en caso de que sea la identidad.
    #
    # Implementación \\
    # Calcula los autovectores y deshomogeniza.
    #
    def puntos_fijos_theta(self):
        paso("Podemos obtener los puntos fijos directamente de:", self._expresion.substitute(self._theta1 == self._theta2), \
                ", pero vamos a usar los autovectores porque es más sencillo para el algoritmo")
        if self.es_identidad():
            return []
        return map(lambda x: Infinity if x[1] == 0 else x[0] / x[1], self._aplicacion.puntos_fijos())
        
    #\m
    # Devuelve los puntos fijos de esta homografía. Devuelve una lista vacía en caso de que sea la identidad.
    #
    # Implementación \\
    # Calcula los autovectores de la matriz asociada.
    #
    def puntos_fijos(self):
        fijos = self.puntos_fijos_theta()
        paso("Sustituimos los puntos fijos ", fijos, "en la recta")
        return map(lambda theta: self._recta[theta], fijos)
    
    #\m
    # Determina si esta homografía es la identidad.
    #
    # Implementación \\
    # Sustituye theta' por theta y resta un lado de la expresión de Möbius al otro para comparar con 0.
    #
    def es_identidad(self):
        return bool(self._expresion_mobius.lhs().substitute(self._theta2 == self._theta1) - self._expresion_mobius.rhs() == 0)
    
    #\m
    # Determina si esta homografía es una elación, esto es si tiene un punto fijo doble.
    #
    # Implementación \\
    # Calcula los puntos fijos y compara.
    #
    def es_elacion(self):
        return len(self.puntos_fijos()) == 1
    
    #\m
    # Determina si esta homografía es una involución, esto es si self^2 == self.
    #
    # Implementación \\
    # Comprueba que la traza sea nula.
    #    
    def es_involucion(self):
        return self._matriz.trace() == 0
    
    #\m
    # Calcula el módulo de esta homografía, esto es el cociente de sus autovalores (en cualquier orden).
    #
    def modulo(self):
        paso("El módulo se calcula mediante el cociente de los autovalores:")
        autovalores = self.autovalores()
        paso("Autovalores:", autovalores)
        # Autovalor doble: la división siempre da 1
        if len(autovalores) == 1:
            return 1
        return autovalores[0] / autovalores[1]
        
    #\m
    # Devuelve esta misma homografía en una referencia adecuada de forma que su ecuación quede de la forma más
    # simple posible. Para las involuciones no devuelve la forma general theta*theta' = 1. Para ello usar simplificar_involucion().
    #
    # Implementación \\
    # En función de los puntos fijos devuelve: \\
    # infinitos -> ella misma (identidad) \\
    # 1         -> matriz (1, 1; 0, 1) con referencia el punto fijo y dos arbitrarios, multiplicada por lambda = delta/beta \\
    # 2         -> matriz (-beta/gamma, 0; 0, 1) con referencia los puntos fijos y uno arbitrario
    #
    def simplificar(self):
        paso("Calculamos los puntos fijos")
        fijos = self.puntos_fijos()
        # Identidad: ya simplificada
        if len(fijos) == 0:
            return self
        ref = self._recta.referencia()
        # Elación: matriz de Jordan fija y referencia el punto fijo y dos arbitrarios
        if len(fijos) == 1:
            # Cogemos dos de los puntos de la referencia que no sean el fijo para la nueva referencia
            i = 0 if ref[0] != fijos[0] else 2
            j = 1 if ref[1] != fijos[0] else 2
            lamb = self._matriz[0][1] / -self._matriz[0][0]
            paso("Usamos el punto fijo como primer punto de la nueva referencia y dejamos dos de los otros dos:")
            paso("Cambiamos de referencia theta' = lambda*theta, con lambda = ", \
                self._matriz[0][1], "/", -self._matriz[0][0], "=", lamb)
            # Calculamos la nueva referencia obteniendo
            nueva0 = self._recta[lamb * self._recta.coordenada(fijos[0])]
            nueva1 = self._recta[lamb * self._recta.coordenada(ref[i])]
            nueva2 = self._recta[lamb * self._recta.coordenada(ref[j])]
            paso("R={", nueva0, ",", nueva1, ";", nueva2, "}")
            return homografia_recta(matrix([[1, 1], [0, 1]]), recta_proyectiva(nueva0, nueva1, nueva2))
        # Dos puntos fijos distintos, podría ser involución
        # Aseguramos que el tercer punto de la referencia no coincida con los fijos
        i = 0 if ref[0] != fijos[0] and ref[0] != fijos[1] else \
            1 if ref[1] != fijos[0] and ref[1] != fijos[1] else \
            2
        _no_pasos()
        modulo = self.modulo()
        # Se debe cumplir que {M, N; P, P'} = modulo (con M, N fijos, P arbitrario)
        # Si no, hay que escogerlos al revés
        j = 0 if razon_doble(fijos[0], fijos[1], ref[i], self(ref[i])) == modulo else 1
        _no_pasos(False)
        paso("Usamos los puntos fijos como primeros puntos (en el orden correcto) de la nueva referencia y mantenemos el otro:")
        paso("R={", fijos[j], ",", fijos[(j + 1) % 2], ";", ref[i], "}")
        return homografia_recta(matrix([[modulo, 0], [0, 1]]), recta_proyectiva(fijos[j], fijos[(j + 1) % 2], ref[i]))
    
    #\m
    # Devuelve esta misma involución en una referencia adecuada de forma que su ecuación quede theta*theta' = 1.
    #
    # Implementación \\
    # Asumiendo que esta homografía es una involución, los dos puntos fijos separan armónicamente cualquier par de la involución.
    # Por tanto, basta con elegir cualquier punto y su imagen como Infinity y 0 y un punto fijo como 1. De esta forma, el otro
    # punto fijo tendrá coordenada -1 y con ello obtenemos la ecuación deseada.
    #
    def simplificar_involucion(self):
        assert self.es_involucion(), "Para simplificar como involucion la homografia debe ser una involucion"
        paso("Calculamos los puntos fijos")
        fijos = self.puntos_fijos()
        ref = self._recta.referencia()
        i = 0 if ref[0] != fijos[0] and ref[0] != fijos[1] else \
            1 if ref[1] != fijos[0] and ref[1] != fijos[1] else \
            2
        paso("Usamos un punto arbitrario y su imagen como primeros puntos de la referencia y uno de los fijos como unidad")
        paso("R={", ref[i], ",", self(ref[i]), ",", fijos[0], "}")
        return homografia_recta(matrix([[0, 1], [1, 0]]), recta_proyectiva(ref[i], self(ref[i]), fijos[0]))
    
    #\m
    # Operador *. Devuelve la composición de las homografías ((self o otra)(theta) = self(otra(theta))).
    #
    # Uso: h * k (h y k son homografías de una misma recta).
    #
    # Implementación \\
    # Devuelve una nueva homografía cuya matriz sea el producto de las de ambas.
    #
    # Parámetros \\
    # otra: homografia_recta - homografía con la que componer
    #
    def __mul__(self, otra):
        assert self._recta == otra._recta, "Las homografias deben ser de la misma recta"
        return homografia_recta(self._matriz * otra._matriz, self._recta)
        
    #\m
    # Operador ^ (ó **). Devuelve el resultado de componer una homografía con sí misma n veces (^-1 devuelve la inversa).
    #
    # Uso h^n (ó h**n) (h es una homografía de la recta y n un entero).
    #
    # Implementación \\
    # Devuelve una nueva homografía cuya matriz es la de esta elevada a n.
    #
    # Parámetros \\
    # n: entero - exponente al que elevar
    #
    def __pow__(self, n):
        return homografia_recta(self._matriz^n, self._recta)
            
    # Métodos auxiliares
            
    def __sustituir(self, x):
        # Para infinito se pierden las componentes sin theta
        if x == Infinity or x == -Infinity:
            if self._matriz[1][0] == 0:
                return Infinity
            else:
                return self._matriz[0][0] / self._matriz[1][0]
        # En cualquier otro caso se sustituye en la expresión de Möbius
        return self._expresion_mobius.substitute(self._theta1 == x).rhs()
        
    def __repr__(self):
        return "<Homografia " + str(self._expresion) + " de " + str(self._recta) + ">"
